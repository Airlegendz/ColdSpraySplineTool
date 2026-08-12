"""
Parametric cold-spray nozzle geometry generator.

Produces a 2D axisymmetric wall profile (x, r) for a converging-diverging
nozzle: fixed convergent section -> circular fillet arc at the throat ->
cubic Bezier divergent section -> optional straight barrel section(s).

COORDINATE CONVENTION (see README.md for full detail):
  - x is the axial coordinate, in meters, increasing in the flow direction.
  - x = 0 is defined at the nozzle inlet (start of the fixed convergent
    section).
  - r is the wall radius, in meters, measured from the centerline.
  - The profile is a single wall trace (upper half); mirror about r = 0 to
    get the full axisymmetric cross-section (done by --plot).

DESIGN CHOICE (flagged for confirmation): the divergent section is
parameterized as a single cubic Bezier curve through 4 control points
(throat point, control point 1, control point 2, exit point). A cubic
Bezier was chosen over a B-spline because it has a fixed, small number of
control points that map 1:1 onto the stated design parameters, with no
extra knot-vector or degree choices to make. This does NOT interpolate
control_point_1/2 exactly in the strict sense of "the curve passes through
them" for an unconstrained Bezier -- for a *cubic* Bezier with exactly 4
control points, the curve *does* pass through all 4 control points only
at the endpoints (P0, P3); interior points (P1, P2) act as tangent/shape
handles, not as points the curve necessarily touches. Given the task's
requirement that the curve pass through control_point_1 and control_point_2
as well, we instead build a clamped cubic Bezier *spline* (piecewise cubic,
C1-continuous) through all four points: throat, cp1, cp2, exit. This is
documented here as the load-bearing design decision.
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Optional

logger = logging.getLogger("nozzle_geometry")

Point = Tuple[float, float]


# ----------------------------------------------------------------------------
# Defaults that are NOT established by the task and MUST be confirmed before
# a full CFD sweep is run. They are exposed as configurable parameters (see
# GeometryConfig below), not hardcoded into the algorithm.
# ----------------------------------------------------------------------------
DEFAULT_INLET_RADIUS = 0.010          # m, fixed convergent-section inlet radius. UNCONFIRMED.
DEFAULT_CONVERGENT_LENGTH = 0.030     # m, axial length of fixed convergent section. UNCONFIRMED.
# Single authoritative code-fallback default for each threshold (used only
# when a caller builds GeometryConfig directly without specifying the
# field -- e.g. ad hoc scripts or tests). config_example.yaml's `fixed`
# section sets the same values explicitly and, for anything driven through
# sweep.py, THAT value always wins: sweep.py forwards every key present in
# the YAML config as an explicit GeometryConfig kwarg, which overrides the
# dataclass field default below. These two are kept numerically identical
# on purpose so "the code default" and "the example sweep default" are the
# same UNCONFIRMED number until real values are validated.
DEFAULT_MAX_FILLET_FRACTION = 0.6     # throat_fillet_radius <= this * throat_radius. UNCONFIRMED.
DEFAULT_MAX_CURVATURE = 800.0         # 1/m, cap on |d2r/dx2| along the divergent curve. UNCONFIRMED.
DEFAULT_DIVERGENT_SAMPLES = 200       # samples along the divergent spline (>=50 required)
DEFAULT_FILLET_SAMPLES = 30
DEFAULT_CONVERGENT_SAMPLES = 40
DEFAULT_BARREL_SAMPLES_PER_SECTION = 10


@dataclass
class GeometryConfig:
    throat_radius: float
    throat_fillet_radius: float
    control_point_1_radius: float
    control_point_1_position: float
    control_point_2_radius: float
    control_point_2_position: float
    exit_radius: float
    exit_position: float  # axial distance from the throat to the exit point, meters. Independently swept -- see README.
    barrel_length: float
    barrel_count: int
    barrel_positions: List[float] = field(default_factory=list)

    # Fixed / convergent-section geometry (not part of the swept design
    # space, but exposed so they aren't hidden magic numbers).
    inlet_radius: float = DEFAULT_INLET_RADIUS
    convergent_length: float = DEFAULT_CONVERGENT_LENGTH

    # Constraint thresholds -- configurable, UNCONFIRMED defaults.
    max_fillet_fraction: float = DEFAULT_MAX_FILLET_FRACTION
    max_curvature: float = DEFAULT_MAX_CURVATURE

    # Sampling density.
    divergent_samples: int = DEFAULT_DIVERGENT_SAMPLES
    fillet_samples: int = DEFAULT_FILLET_SAMPLES
    convergent_samples: int = DEFAULT_CONVERGENT_SAMPLES
    barrel_samples_per_section: int = DEFAULT_BARREL_SAMPLES_PER_SECTION

    def to_dict(self) -> dict:
        return asdict(self)


class GeometryValidationError(Exception):
    """Raised (or logged) when a sampled parameter set produces an invalid geometry."""


@dataclass
class GeometryResult:
    points: List[Point]
    valid: bool
    rejection_reason: Optional[str]
    config: GeometryConfig


# ----------------------------------------------------------------------------
# Convergent section (fixed profile, not parameterized)
# ----------------------------------------------------------------------------
def _convergent_profile(cfg: GeometryConfig, throat_x: float) -> List[Point]:
    """
    Smooth fixed contraction from inlet_radius down to throat_radius, using a
    cosine (smoothstep-like) blend so slope is zero at the inlet and matches
    the fillet's tangent direction near the throat.
    """
    n = cfg.convergent_samples
    pts: List[Point] = []
    x0 = throat_x - cfg.convergent_length
    for i in range(n):
        t = i / (n - 1)
        # Raised-cosine contraction: smooth, monotonically decreasing radius.
        s = 0.5 * (1 - math.cos(math.pi * t))
        r = cfg.inlet_radius - (cfg.inlet_radius - cfg.throat_radius) * s
        x = x0 + t * cfg.convergent_length
        pts.append((x, r))
    return pts


def _convergent_slope_at_throat(cfg: GeometryConfig) -> float:
    """dr/dx of the convergent profile evaluated at its end (x = throat_x)."""
    # r(t) = inlet_r - (inlet_r - throat_r) * 0.5*(1-cos(pi t)), t in [0,1], x = x0 + t*L
    # dr/dt at t=1 = -(inlet_r - throat_r) * 0.5*pi*sin(pi*1) = 0
    return 0.0


# ----------------------------------------------------------------------------
# Throat fillet: circular arc tangent to the (flat, zero-slope) convergent
# wall at the throat and tangent to the divergent curve's start.
# ----------------------------------------------------------------------------
def _fillet_arc(cfg: GeometryConfig, throat_x: float) -> Tuple[List[Point], Point, float]:
    """
    Builds a circular arc of radius throat_fillet_radius that is tangent to
    the horizontal convergent wall (slope 0) at the throat centerline point
    and blends into the divergent section. Returns (points, end_point,
    end_slope) where end_point/end_slope feed the divergent Bezier spline's
    first tangent constraint.

    Since the convergent wall arrives at slope 0 (a straight, constant-radius
    approach as t->1), the fillet is modeled as a quarter-ish arc that curves
    the wall from flat (slope 0) into the divergent section's initial slope,
    center offset outward by throat_fillet_radius from the throat point.
    """
    rf = cfg.throat_fillet_radius
    rt = cfg.throat_radius
    if rf <= 0:
        return [(throat_x, rt)], (throat_x, rt), 0.0

    # Estimate the target initial slope of the divergent section as the
    # chord slope from throat to control_point_1 (a reasonable tangent
    # target that keeps the fillet small and local to the throat).
    dx = cfg.control_point_1_position
    dr = cfg.control_point_1_radius - rt
    target_slope = dr / dx if dx > 0 else 0.0
    # Clamp to non-negative: the fillet blends from the flat (slope-0)
    # convergent wall into the diverging wall, so it must not curve
    # backward in x or r. A negative target here (control_point_1_radius
    # < throat_radius) is a downstream parameterization error that the
    # monotonicity check will flag on its own merits -- the fillet itself
    # should never regress.
    target_slope = max(target_slope, 0.0)
    theta_end = math.atan(target_slope)  # arc sweep angle from flat (0) to target_slope

    n = cfg.fillet_samples
    pts: List[Point] = []
    # Arc center is above the throat point by rf (center of curvature for a
    # wall that curves radius-outward as x increases).
    cx, cr = throat_x, rt + rf
    for i in range(n):
        theta = theta_end * (i / (n - 1)) if n > 1 else theta_end
        x = cx + rf * math.sin(theta)
        r = cr - rf * math.cos(theta)
        pts.append((x, r))
    end_point = pts[-1]
    end_slope = math.tan(theta_end)
    return pts, end_point, end_slope


# ----------------------------------------------------------------------------
# Divergent section: piecewise-cubic Hermite spline (converted to Bezier for
# evaluation) through throat_end -> cp1 -> cp2 -> exit, C1-continuous.
#
# Tangents at the two interior knots (cp1, cp2) are estimated Catmull-Rom
# style from their actual neighboring knots. The two end tangents are
# *clamped* rather than extrapolated from a phantom point, so the curve
# is smoothly tangent to its neighbors on both sides:
#   - at the throat, the tangent equals the fillet arc's actual exit slope
#     (so throat fillet -> divergent spline is C1-continuous)
#   - at the exit, the tangent is forced horizontal (slope 0), matching the
#     constant-radius barrel section (or a free jet exit) that follows, so
#     divergent spline -> barrel is C1-continuous instead of kinking.
# ----------------------------------------------------------------------------
def _hermite_to_bezier(p0, p1, m0, m1):
    """Convert one Hermite segment (p0->p1, tangent vectors m0/m1) to a cubic Bezier."""
    b0 = p0
    b1 = (p0[0] + m0[0] / 3.0, p0[1] + m0[1] / 3.0)
    b2 = (p1[0] - m1[0] / 3.0, p1[1] - m1[1] / 3.0)
    b3 = p1
    return b0, b1, b2, b3


def _bezier_point(b0, b1, b2, b3, t):
    mt = 1 - t
    x = mt**3 * b0[0] + 3 * mt**2 * t * b1[0] + 3 * mt * t**2 * b2[0] + t**3 * b3[0]
    r = mt**3 * b0[1] + 3 * mt**2 * t * b1[1] + 3 * mt * t**2 * b2[1] + t**3 * b3[1]
    return x, r


def _divergent_spline(cfg: GeometryConfig, throat_point: Point, throat_slope: float) -> Tuple[List[Point], Point]:
    tx, tr = throat_point
    cp1 = (tx + cfg.control_point_1_position, cfg.control_point_1_radius)
    cp2 = (tx + cfg.control_point_2_position, cfg.control_point_2_radius)
    # exit_position is an independent, swept input (axial distance from the
    # throat) -- NOT derived from cp1/cp2 spacing. Divergent length was
    # found to be the single most influential parameter for particle
    # velocity (Badali et al.), so it must be freely explorable by the
    # sampler rather than locked to a fixed multiple of the cp1-cp2 span.
    exit_point = (tx + cfg.exit_position, cfg.exit_radius)

    knots = [throat_point, cp1, cp2, exit_point]

    # Tangent vectors at each knot, scaled to the local chord length so
    # segment-to-segment speed stays reasonable (standard Catmull-Rom scaling).
    m_throat = (1.0, throat_slope)
    chord0 = knots[1][0] - knots[0][0]
    m_throat = (chord0, throat_slope * chord0)

    m_cp1 = (knots[2][0] - knots[0][0], knots[2][1] - knots[0][1])
    m_cp1 = tuple(c / 2.0 for c in m_cp1)

    m_cp2 = (knots[3][0] - knots[1][0], knots[3][1] - knots[1][1])
    m_cp2 = tuple(c / 2.0 for c in m_cp2)

    chord_last = knots[3][0] - knots[2][0]
    if cfg.barrel_count > 0:
        # Clamp the exit tangent horizontal so the spline blends smoothly
        # (C1) into the constant-radius barrel that follows.
        m_exit = (chord_last, 0.0)
    else:
        # No barrel: use the last chord's own slope as the end tangent (a
        # "cardinal spline" boundary condition). This lets the curve open
        # freely at the exit and, notably, makes the spline degenerate
        # exactly to a straight line when all four knots are colinear.
        m_exit = (chord_last, knots[3][1] - knots[2][1])

    tangents = [m_throat, m_cp1, m_cp2, m_exit]

    n_segments = len(knots) - 1
    samples_per_segment = max(cfg.divergent_samples // n_segments, 12)

    pts: List[Point] = []
    for seg in range(n_segments):
        p0, p1 = knots[seg], knots[seg + 1]
        m0, m1 = tangents[seg], tangents[seg + 1]
        b0, b1, b2, b3 = _hermite_to_bezier(p0, p1, m0, m1)
        last = samples_per_segment
        for i in range(last):
            t = i / (last - 1) if last > 1 else 0.0
            if seg > 0 and i == 0:
                continue  # avoid duplicating the shared knot between segments
            pts.append(_bezier_point(b0, b1, b2, b3, t))
    return pts, exit_point


# ----------------------------------------------------------------------------
# Barrel section(s): straight, constant-radius runs inserted at given axial
# positions (relative to the exit of the divergent section).
# ----------------------------------------------------------------------------
def _barrel_sections(cfg: GeometryConfig, exit_point: Point) -> List[Point]:
    if cfg.barrel_count == 0:
        return []
    pts: List[Point] = []
    ex, er = exit_point
    for pos in sorted(cfg.barrel_positions):
        x0 = ex + pos
        n = cfg.barrel_samples_per_section
        for i in range(n):
            x = x0 + (i / (n - 1)) * cfg.barrel_length
            pts.append((x, er))
    return pts


# ----------------------------------------------------------------------------
# Constraint checks
# ----------------------------------------------------------------------------
# Every rejection reason is prefixed with one of these stable category tags
# so callers (e.g. sweep.py's rejection-rate breakdown) can bucket by
# constraint type instead of by the full, value-specific message text.
REASON_POSITIVE_LENGTHS = "positive_lengths"
REASON_FILLET_BOUND = "fillet_bound"
REASON_BARREL_OVERLAP = "barrel_overlap"
REASON_MONOTONICITY = "monotonicity"
REASON_CURVATURE = "curvature"


def _check_positive_lengths(cfg: GeometryConfig) -> Optional[str]:
    if not (0 < cfg.control_point_1_position < cfg.control_point_2_position < cfg.exit_position):
        return (
            f"[{REASON_POSITIVE_LENGTHS}] non-increasing divergent-section positions: "
            f"cp1_pos={cfg.control_point_1_position}, cp2_pos={cfg.control_point_2_position}, "
            f"exit_position={cfg.exit_position}"
        )
    if cfg.barrel_count > 0 and cfg.barrel_length <= 0:
        return f"[{REASON_POSITIVE_LENGTHS}] non-positive barrel_length={cfg.barrel_length}"
    if cfg.convergent_length <= 0:
        return f"[{REASON_POSITIVE_LENGTHS}] non-positive convergent_length={cfg.convergent_length}"
    return None


def _check_fillet_bound(cfg: GeometryConfig) -> Optional[str]:
    if cfg.throat_fillet_radius < 0:
        return f"[{REASON_FILLET_BOUND}] negative throat_fillet_radius={cfg.throat_fillet_radius}"
    limit = cfg.max_fillet_fraction * cfg.throat_radius
    if cfg.throat_fillet_radius > limit:
        return (
            f"[{REASON_FILLET_BOUND}] throat_fillet_radius={cfg.throat_fillet_radius} exceeds "
            f"max_fillet_fraction*throat_radius={limit} "
            f"(max_fillet_fraction={cfg.max_fillet_fraction}, UNCONFIRMED default)"
        )
    return None


def _check_barrel_overlap(cfg: GeometryConfig, exit_x: float) -> Optional[str]:
    if cfg.barrel_count <= 1:
        return None
    if len(cfg.barrel_positions) != cfg.barrel_count:
        return (
            f"[{REASON_BARREL_OVERLAP}] barrel_positions length {len(cfg.barrel_positions)} != "
            f"barrel_count {cfg.barrel_count}"
        )
    intervals = sorted((p, p + cfg.barrel_length) for p in cfg.barrel_positions)
    for (s0, e0), (s1, e1) in zip(intervals, intervals[1:]):
        if s1 < e0:
            return f"[{REASON_BARREL_OVERLAP}] overlapping barrel sections: [{s0},{e0}] and [{s1},{e1}]"
    return None


def _check_monotonic_and_curvature(points: List[Point], throat_x: float, cfg: GeometryConfig) -> Optional[str]:
    # Only check from the throat onward (radius may legitimately decrease
    # through the convergent section before the throat).
    seg = [(x, r) for x, r in points if x >= throat_x - 1e-9]
    seg.sort(key=lambda p: p[0])

    prev_r = None
    for x, r in seg:
        if prev_r is not None and r < prev_r - 1e-9:
            return f"[{REASON_MONOTONICITY}] non-monotonic radius at x={x:.6g}: r={r:.6g} < prev_r={prev_r:.6g}"
        prev_r = r

    # Second-derivative curvature proxy via finite differences on non-uniform
    # samples (the divergent spline sampling is not perfectly even in x).
    xs = [p[0] for p in seg]
    rs = [p[1] for p in seg]
    for i in range(1, len(seg) - 1):
        h1 = xs[i] - xs[i - 1]
        h2 = xs[i + 1] - xs[i]
        if h1 <= 0 or h2 <= 0:
            continue
        d2r = 2 * (h1 * rs[i + 1] - (h1 + h2) * rs[i] + h2 * rs[i - 1]) / (h1 * h2 * (h1 + h2))
        if abs(d2r) > cfg.max_curvature:
            return (
                f"[{REASON_CURVATURE}] curvature exceeds max_curvature at x={xs[i]:.6g}: "
                f"|d2r/dx2|={abs(d2r):.6g} > {cfg.max_curvature} "
                f"(UNCONFIRMED default threshold)"
            )
    return None


# ----------------------------------------------------------------------------
# Top-level generator
# ----------------------------------------------------------------------------
def generate_geometry(cfg: GeometryConfig) -> GeometryResult:
    """
    Builds the full wall profile and runs all constraint checks. Invalid
    geometries are still returned (with valid=False and a rejection_reason)
    rather than raised, so callers (e.g. the sweep driver) can log and skip
    them without a try/except per-sample.
    """
    reason = _check_positive_lengths(cfg)
    if reason:
        logger.warning("Rejected geometry: %s", reason)
        return GeometryResult(points=[], valid=False, rejection_reason=reason, config=cfg)

    reason = _check_fillet_bound(cfg)
    if reason:
        logger.warning("Rejected geometry: %s", reason)
        return GeometryResult(points=[], valid=False, rejection_reason=reason, config=cfg)

    throat_x = cfg.convergent_length  # inlet at x=0

    convergent_pts = _convergent_profile(cfg, throat_x)
    fillet_pts, throat_end_point, fillet_end_slope = _fillet_arc(cfg, throat_x)
    # Drop the duplicate junction point (convergent's last sample coincides
    # with the fillet's first sample) to avoid a zero-length segment.
    if convergent_pts and fillet_pts and convergent_pts[-1] == fillet_pts[0]:
        convergent_pts = convergent_pts[:-1]
    divergent_pts, exit_point = _divergent_spline(cfg, throat_end_point, fillet_end_slope)
    barrel_pts = _barrel_sections(cfg, exit_point)

    reason = _check_barrel_overlap(cfg, exit_point[0])
    if reason:
        logger.warning("Rejected geometry: %s", reason)
        return GeometryResult(points=[], valid=False, rejection_reason=reason, config=cfg)

    all_pts = convergent_pts + fillet_pts + divergent_pts + barrel_pts

    reason = _check_monotonic_and_curvature(all_pts, throat_x, cfg)
    if reason:
        logger.warning("Rejected geometry: %s", reason)
        return GeometryResult(points=all_pts, valid=False, rejection_reason=reason, config=cfg)

    return GeometryResult(points=all_pts, valid=True, rejection_reason=None, config=cfg)


# ----------------------------------------------------------------------------
# Regression test: degenerate case -> near-straight cone
# ----------------------------------------------------------------------------
def _regression_test_straight_cone(plot: bool = False) -> bool:
    """
    throat_fillet_radius ~ 0 and control points placed on the straight line
    from throat to exit should collapse the spline to a near-straight cone.
    Returns True if the max deviation from the ideal straight line is small.
    """
    throat_r = 0.003
    exit_r = 0.010

    cp1_pos = 0.008
    cp2_pos = 0.016
    exit_pos = 0.040  # independent exit_position, chosen beyond cp2_pos

    # Place cp1/cp2 radii exactly on the straight line throat->exit.
    slope = (exit_r - throat_r) / exit_pos
    cp1_r = throat_r + slope * cp1_pos
    cp2_r = throat_r + slope * cp2_pos

    cfg = GeometryConfig(
        throat_radius=throat_r,
        throat_fillet_radius=1e-6,
        control_point_1_radius=cp1_r,
        control_point_1_position=cp1_pos,
        control_point_2_radius=cp2_r,
        control_point_2_position=cp2_pos,
        exit_radius=exit_r,
        exit_position=exit_pos,
        barrel_length=0.0,
        barrel_count=0,
        barrel_positions=[],
        max_curvature=1e9,  # a straight line still has ~0 curvature; keep generous
    )
    result = generate_geometry(cfg)
    throat_x = cfg.convergent_length
    max_dev = 0.0
    for x, r in result.points:
        if x < throat_x - 1e-9:
            continue
        dx = x - throat_x
        r_ideal = throat_r + slope * dx
        max_dev = max(max_dev, abs(r - r_ideal))

    tol = 0.02 * (exit_r - throat_r)  # 2% of the total radius change
    passed = result.valid and max_dev < tol
    logger.info(
        "Regression test (straight-cone degeneracy): valid=%s max_dev=%.6g tol=%.6g -> %s",
        result.valid, max_dev, tol, "PASS" if passed else "FAIL",
    )

    if plot:
        from plotting import plot_profile
        plot_profile(result.points, "regression_straight_cone.png",
                     title="Regression: degenerate spline vs straight cone")

    return passed


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Nozzle geometry regression test.")
    parser.add_argument("--plot", action="store_true", help="Render the regression-test profile to a PNG.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    ok = _regression_test_straight_cone(plot=args.plot)
    print("Regression test:", "PASS" if ok else "FAIL")
