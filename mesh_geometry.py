"""
Gmsh-based 2D axisymmetric meshing pipeline for cold-spray nozzle geometry
CSVs produced by nozzle_geometry.py / sweep.py.

Consumes a single geometry_XXXX.csv (the wall (x, r) polyline, meters,
x=0 at the inlet, ordered inlet -> end of last barrel section -- see
README.md's "Coordinate convention" section) and builds a closed 2D
axisymmetric flow domain:

    wall (spline through the CSV points)   -- top
    outlet (vertical line at x = x_max)    -- right
    axis (r = 0, the centerline)           -- bottom
    inlet (vertical line at x = 0)         -- left

meshed with a structured (transfinite, recombined-to-quad) grid, radially
clustered toward the wall for boundary-layer resolution, and exported to
Gmsh's native .msh format for Fluent import.

Physical groups (see MESHING.md for the exact Fluent import steps):
    "wall"   (1D) -- the nozzle wall spline
    "axis"   (1D) -- the centerline; Fluent must be told this is the
                     axisymmetric Axis boundary type, not a regular wall
    "inlet"  (1D) -- x = 0 edge
    "outlet" (1D) -- x = x_max edge
    "fluid"  (2D) -- the domain interior
"""

from __future__ import annotations

import csv
import json
import logging
import math
from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple

import gmsh

logger = logging.getLogger("mesh_geometry")

Point = Tuple[float, float]

# Points closer than this (in x AND r, meters) are treated as duplicates and
# collapsed -- the generator emits exact-duplicate junction points at
# section boundaries (e.g. throat, divergent-to-barrel), which would
# otherwise produce zero-length Gmsh curve segments.
DEDUP_TOL = 1e-9


# ----------------------------------------------------------------------------
# Mesh configuration -- UNCONFIRMED defaults, same status as
# nozzle_geometry.py's max_curvature/max_fillet_fraction: exposed as
# configurable parameters with a logged default, not validated engineering
# values. The right values depend on the actual flow Reynolds number and
# haven't been established. See MESHING.md.
# ----------------------------------------------------------------------------
DEFAULT_N_AXIAL = 150          # transfinite node count along wall/axis (flow direction). UNCONFIRMED.
DEFAULT_N_RADIAL = 40           # transfinite node count along inlet/outlet (radial direction). UNCONFIRMED.
DEFAULT_FIRST_CELL_HEIGHT = 5e-6  # m, target wall-adjacent cell height for BL clustering. UNCONFIRMED.
DEFAULT_BL_GROWTH_RATIO_FALLBACK = 1.2  # used only if first_cell_height can't be solved for. UNCONFIRMED.
DEFAULT_MESH_FORMAT = "msh2"    # "msh2" (MSH2, older/pickier Fluent versions) or "msh4". UNCONFIRMED
                                 # which your specific Fluent install expects -- see MESHING.md.
DEFAULT_QUALITY_THRESHOLD = 0.02  # SOFT/informational threshold on minSICN. UNCONFIRMED, and read the
                                   # note below before treating this as "the pass/fail bar":
                                   #
                                   # Elements with quality <= 0 are inverted/degenerate -- always a real
                                   # meshing failure, checked separately (n_inverted_elements) and always
                                   # fatal regardless of this threshold.
                                   #
                                   # Elements with 0 < quality < quality_threshold are just WARNED about,
                                   # not failed, because boundary-layer meshes (which this pipeline
                                   # deliberately builds -- see first_cell_height) intentionally produce
                                   # very thin, high-aspect-ratio cells at the wall. minSICN penalizes
                                   # aspect ratio as well as skew, so a tight, otherwise-correct BL first
                                   # cell can legitimately score quite low here (observed ~0.003 on the
                                   # regression fixture at first_cell_height=5e-6m) without being wrong.
                                   # This default is set low enough to avoid flagging that expected case;
                                   # it still needs validating against your actual Reynolds number and
                                   # solver's own mesh-quality expectations (Fluent reports orthogonal
                                   # quality / aspect ratio separately, which is the more meaningful check
                                   # for BL cells specifically -- see MESHING.md).


@dataclass
class MeshConfig:
    n_axial: int = DEFAULT_N_AXIAL
    n_radial: int = DEFAULT_N_RADIAL
    first_cell_height: float = DEFAULT_FIRST_CELL_HEIGHT
    bl_growth_ratio_fallback: float = DEFAULT_BL_GROWTH_RATIO_FALLBACK
    mesh_format: str = DEFAULT_MESH_FORMAT
    quality_threshold: float = DEFAULT_QUALITY_THRESHOLD
    check_quality: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class MeshResult:
    geometry_id: str
    success: bool
    msh_path: Optional[str]
    png_path: Optional[str]
    n_nodes: int
    n_elements: int
    min_quality: Optional[float]
    n_inverted_elements: int
    quality_flagged: bool
    error: Optional[str]


# ----------------------------------------------------------------------------
# CSV I/O
# ----------------------------------------------------------------------------
def read_geometry_csv(csv_path: str) -> List[Point]:
    """
    Reads a geometry_XXXX.csv (x, r columns, meters) and dedups consecutive
    near-identical points (see DEDUP_TOL) so the wall spline has no
    zero-length segments at section junctions.
    """
    points: List[Point] = []
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        if header != ["x", "r"]:
            raise ValueError(f"{csv_path}: expected header ['x','r'], got {header}")
        for row in reader:
            x, r = float(row[0]), float(row[1])
            if points and abs(x - points[-1][0]) < DEDUP_TOL and abs(r - points[-1][1]) < DEDUP_TOL:
                continue
            points.append((x, r))
    if len(points) < 2:
        raise ValueError(f"{csv_path}: fewer than 2 usable points after dedup")
    return points


# ----------------------------------------------------------------------------
# Boundary-layer grading: solve for the geometric-progression ratio that
# achieves a target first-cell height over n_radial intervals spanning a
# known radial length, so first_cell_height is a real configurable
# parameter rather than just a documented aspiration.
# ----------------------------------------------------------------------------
def _solve_growth_ratio(length: float, n_intervals: int, first_cell_height: float,
                         fallback_ratio: float) -> float:
    """
    Solves length = first_cell_height * (ratio**n - 1) / (ratio - 1) for
    ratio > 1 via bisection (the sum of a geometric series of n cells whose
    first term is first_cell_height). Falls back to fallback_ratio if the
    target first_cell_height is infeasible (e.g. requests a first cell
    larger than a uniform split would give, or n_intervals <= 1).
    """
    if n_intervals <= 1 or first_cell_height <= 0 or first_cell_height >= length:
        return fallback_ratio

    def total_length(ratio: float) -> float:
        if abs(ratio - 1.0) < 1e-12:
            return first_cell_height * n_intervals
        return first_cell_height * (ratio ** n_intervals - 1) / (ratio - 1)

    # Uniform spacing (ratio=1) gives the minimum possible total length for
    # a fixed first_cell_height and n_intervals as ratio increases from 1;
    # if even that exceeds the target length, first_cell_height is too
    # large for this n_radial and curve length -- fall back rather than
    # returning a nonsensical (ratio < 1, growing the wrong direction) fit.
    if total_length(1.0) >= length:
        logger.warning(
            "first_cell_height=%.3g infeasible for n_radial=%d over length=%.3g "
            "(uniform spacing alone already exceeds target length) -- falling back to ratio=%.3g",
            first_cell_height, n_intervals, length, fallback_ratio,
        )
        return fallback_ratio

    lo, hi = 1.0 + 1e-9, 10.0
    while total_length(hi) < length and hi < 1e6:
        hi *= 2
    for _ in range(100):
        mid = (lo + hi) / 2
        if total_length(mid) < length:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ----------------------------------------------------------------------------
# Domain construction + meshing
# ----------------------------------------------------------------------------
def _resample_uniform_in_x(points: List[Point], n: int) -> List[Point]:
    """
    Resamples a monotonic-in-x polyline to exactly n points uniformly
    spaced in x (linear interpolation), preserving the exact endpoints.

    This exists to fix a real mesh-skew defect: Gmsh's transfinite curve
    meshing places nodes uniform in *arc length* by default. The axis
    curve is a straight line, so its arc-length parameterization is
    already uniform in x -- but the wall is a curved spline, where uniform
    arc-length spacing is NOT uniform in x (steep sections get
    disproportionately more nodes per unit x). Since transfinite surface
    meshing connects same-index nodes between opposite curves (wall<->axis)
    with straight "radial" lines, that parameterization mismatch visibly
    tilts/skews those lines wherever the wall curves -- confirmed visually
    in the throat-region zoom of the mesh preview PNG before this fix.
    Building the wall spline's underlying geometry from x-uniform points
    (matching the axis's natural spacing) instead of the raw CSV samples
    keeps wall and axis node indices aligned in x, eliminating that skew.
    """
    x0, x1 = points[0][0], points[-1][0]
    if n < 2 or x1 <= x0:
        return points
    step = (x1 - x0) / (n - 1)
    orig_xs = [p[0] for p in points]
    orig_rs = [p[1] for p in points]

    resampled: List[Point] = []
    j = 0
    for i in range(n):
        x = x0 + i * step
        while j < len(orig_xs) - 2 and orig_xs[j + 1] < x:
            j += 1
        x_lo, x_hi = orig_xs[j], orig_xs[j + 1]
        r_lo, r_hi = orig_rs[j], orig_rs[j + 1]
        t = (x - x_lo) / (x_hi - x_lo) if x_hi > x_lo else 0.0
        resampled.append((x, r_lo + t * (r_hi - r_lo)))
    resampled[0] = points[0]
    resampled[-1] = points[-1]
    return resampled


@dataclass
class _Domain:
    surfaces: List[int]
    wall_curves: List[int]
    axis_curves: List[int]
    inlet_curve: int
    outlet_curve: int
    first_wall_point: int


def _build_domain(points: List[Point], cfg: MeshConfig) -> _Domain:
    """
    Builds the closed 2D domain and applies transfinite meshing, split into
    TWO transfinite blocks at the throat (the minimum-radius point on the
    wall), rather than one single 4-sided block spanning the whole nozzle.

    This split is load-bearing, not cosmetic: a single-block transfinite
    (bilinear/Coons) surface blends its radial node distribution between
    only the inlet and outlet boundary conditions. For a converging-
    diverging "necked" domain like this nozzle, that blend does NOT locally
    adapt to the much smaller wall radius at intermediate (throat) stations
    -- confirmed by direct comparison: even with boundary-layer clustering
    disabled entirely (uniform radial spacing), a single-block mesh still
    showed pronounced radial-line shear through the throat region (visible
    in the mesh preview PNG's zoomed panel). Splitting into two blocks at
    the throat, with a THIRD radial curve (the "interface") graded from the
    throat's own local wall radius, resolved it directly: minSICN improved
    ~19x on the geometry used to diagnose this (0.00018 -> 0.0035), with
    the radial mesh lines visibly near-vertical through the throat
    afterward. Both blocks share the same interface curve and throat point
    (not two coincident-but-separate points), so the mesh is conformal
    across the split -- Fluent sees one continuous fluid zone, not two
    disconnected ones.
    """
    throat_idx = min(range(len(points)), key=lambda i: points[i][1])
    # Guard against a degenerate split if the minimum happens to sit at
    # either endpoint (shouldn't occur for this generator's geometry, but
    # a single-block domain is still valid and better than a crash).
    if throat_idx <= 0 or throat_idx >= len(points) - 1:
        return _build_single_block_domain(points, cfg)

    pts_a = points[:throat_idx + 1]
    pts_b = points[throat_idx:]
    x_throat, r_throat = points[throat_idx]

    # Split n_axial between the two blocks proportional to each segment's
    # share of the original CSV point count (a reasonable proxy for its
    # share of the nozzle's arc length / geometric complexity).
    n_a = max(int(round(cfg.n_axial * len(pts_a) / len(points))), 4)
    n_b = max(cfg.n_axial - n_a + 1, 4)

    wa = _resample_uniform_in_x(pts_a, n_a)
    wb = _resample_uniform_in_x(pts_b, n_b)

    x0, r0 = wa[0]
    x1, r1 = wb[-1]

    pa = [gmsh.model.geo.addPoint(x, r, 0) for x, r in wa]
    # Share the throat point tag between blocks (wb[0] == wa[-1]
    # geometrically) rather than creating a second, coincident-but-separate
    # point -- otherwise the two blocks aren't conformal at the interface.
    pb = [pa[-1]] + [gmsh.model.geo.addPoint(x, r, 0) for x, r in wb[1:]]

    p_axis_start = gmsh.model.geo.addPoint(x0, 0.0, 0)
    p_axis_throat = gmsh.model.geo.addPoint(x_throat, 0.0, 0)
    p_axis_end = gmsh.model.geo.addPoint(x1, 0.0, 0)

    c_wall_a = gmsh.model.geo.addSpline(pa)
    c_wall_b = gmsh.model.geo.addSpline(pb)
    c_axis_a = gmsh.model.geo.addLine(p_axis_start, p_axis_throat)
    c_axis_b = gmsh.model.geo.addLine(p_axis_throat, p_axis_end)
    c_inlet = gmsh.model.geo.addLine(pa[0], p_axis_start)
    c_interface = gmsh.model.geo.addLine(pa[-1], p_axis_throat)
    c_outlet = gmsh.model.geo.addLine(pb[-1], p_axis_end)

    loop_a = gmsh.model.geo.addCurveLoop([c_wall_a, c_interface, -c_axis_a, -c_inlet])
    surf_a = gmsh.model.geo.addPlaneSurface([loop_a])
    loop_b = gmsh.model.geo.addCurveLoop([c_wall_b, c_outlet, -c_axis_b, -c_interface])
    surf_b = gmsh.model.geo.addPlaneSurface([loop_b])

    gmsh.model.geo.mesh.setTransfiniteCurve(c_wall_a, n_a)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_axis_a, n_a)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_wall_b, n_b)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_axis_b, n_b)

    # Radial direction: inlet/interface/outlet, each graded from ITS OWN
    # local wall radius -- this is what the single-block version couldn't
    # do (it only had inlet and outlet to blend between). setTransfiniteCurve
    # numbers nodes starting at the curve's first point -- all three start
    # at the wall end, so coef > 1 grows spacing away from the wall,
    # clustering small cells at the wall as intended.
    inlet_ratio = _solve_growth_ratio(r0, cfg.n_radial - 1, cfg.first_cell_height, cfg.bl_growth_ratio_fallback)
    throat_ratio = _solve_growth_ratio(r_throat, cfg.n_radial - 1, cfg.first_cell_height, cfg.bl_growth_ratio_fallback)
    outlet_ratio = _solve_growth_ratio(r1, cfg.n_radial - 1, cfg.first_cell_height, cfg.bl_growth_ratio_fallback)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_inlet, cfg.n_radial, meshType="Progression", coef=inlet_ratio)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_interface, cfg.n_radial, meshType="Progression", coef=throat_ratio)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_outlet, cfg.n_radial, meshType="Progression", coef=outlet_ratio)

    gmsh.model.geo.mesh.setTransfiniteSurface(surf_a, cornerTags=[pa[0], pa[-1], p_axis_throat, p_axis_start])
    gmsh.model.geo.mesh.setTransfiniteSurface(surf_b, cornerTags=[pb[0], pb[-1], p_axis_end, p_axis_throat])
    gmsh.model.geo.mesh.setRecombine(2, surf_a)
    gmsh.model.geo.mesh.setRecombine(2, surf_b)

    return _Domain(
        surfaces=[surf_a, surf_b], wall_curves=[c_wall_a, c_wall_b], axis_curves=[c_axis_a, c_axis_b],
        inlet_curve=c_inlet, outlet_curve=c_outlet, first_wall_point=pa[0],
    )


def _build_single_block_domain(points: List[Point], cfg: MeshConfig) -> _Domain:
    """Fallback single-block domain, used only if the throat-split guard in
    _build_domain can't find a valid interior split point."""
    x0, r0 = points[0]
    x1, r1 = points[-1]

    wall_points = _resample_uniform_in_x(points, cfg.n_axial)
    wall_point_tags = [gmsh.model.geo.addPoint(x, r, 0) for x, r in wall_points]
    p_axis_start = gmsh.model.geo.addPoint(x0, 0.0, 0)
    p_axis_end = gmsh.model.geo.addPoint(x1, 0.0, 0)

    c_wall = gmsh.model.geo.addSpline(wall_point_tags)
    c_axis = gmsh.model.geo.addLine(p_axis_start, p_axis_end)
    c_inlet = gmsh.model.geo.addLine(wall_point_tags[0], p_axis_start)
    c_outlet = gmsh.model.geo.addLine(wall_point_tags[-1], p_axis_end)

    loop = gmsh.model.geo.addCurveLoop([c_wall, c_outlet, -c_axis, -c_inlet])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.mesh.setTransfiniteCurve(c_wall, cfg.n_axial)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_axis, cfg.n_axial)

    inlet_ratio = _solve_growth_ratio(r0, cfg.n_radial - 1, cfg.first_cell_height, cfg.bl_growth_ratio_fallback)
    outlet_ratio = _solve_growth_ratio(r1, cfg.n_radial - 1, cfg.first_cell_height, cfg.bl_growth_ratio_fallback)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_inlet, cfg.n_radial, meshType="Progression", coef=inlet_ratio)
    gmsh.model.geo.mesh.setTransfiniteCurve(c_outlet, cfg.n_radial, meshType="Progression", coef=outlet_ratio)

    gmsh.model.geo.mesh.setTransfiniteSurface(
        surface, cornerTags=[wall_point_tags[0], wall_point_tags[-1], p_axis_end, p_axis_start]
    )
    gmsh.model.geo.mesh.setRecombine(2, surface)

    return _Domain(
        surfaces=[surface], wall_curves=[c_wall], axis_curves=[c_axis],
        inlet_curve=c_inlet, outlet_curve=c_outlet, first_wall_point=wall_point_tags[0],
    )


def _plot_mesh_png(out_path: str, wall_curves: List[int], title: str = "") -> None:
    """
    Renders a two-panel PNG of the generated quad mesh via matplotlib,
    extracting node/element data directly from the active gmsh model.
    gmsh's own gmsh.write(*.png) requires a graphical (OpenGL) context
    that isn't available in headless batch runs, so this is a
    self-contained alternative using the same node/quality data the
    pipeline already has.

    Left panel: full-domain overview, equal aspect ratio. For a typical
    nozzle (tens of cm long, mm-scale radius) this necessarily looks
    "squashed" -- that's the true physical proportions, not a rendering
    bug -- so it alone can't show cell-level skew or boundary-layer
    clustering, which is the actual purpose of this spot-check image.

    Right panel: zoomed into the throat region (found from the wall
    curve's minimum-radius node) at non-equal aspect, close enough to see
    individual cells -- this is what actually lets a skewed cell near a
    tight fillet be caught visually before it's fed to Fluent.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coord_by_tag = {tag: (node_coords[3 * i], node_coords[3 * i + 1]) for i, tag in enumerate(node_tags)}

    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    polys = []
    for etype, tags, node_lists in zip(elem_types, elem_tags, elem_node_tags):
        _, _, _, num_nodes_per_elem, *_ = gmsh.model.mesh.getElementProperties(etype)
        n = len(tags)
        for i in range(n):
            elem_nodes = node_lists[i * num_nodes_per_elem:(i + 1) * num_nodes_per_elem]
            polys.append([coord_by_tag[t] for t in elem_nodes])

    wall_pts = []
    for c_wall in wall_curves:
        wall_node_tags, wall_coords, _ = gmsh.model.mesh.getNodes(1, c_wall, includeBoundary=True)
        wall_pts.extend((wall_coords[3 * i], wall_coords[3 * i + 1]) for i in range(len(wall_node_tags)))
    throat_x, throat_r = min(wall_pts, key=lambda p: p[1])

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(14, 5))

    ax_full.add_collection(PolyCollection(polys, facecolor="none", edgecolor="steelblue", linewidth=0.3))
    ax_full.autoscale()
    ax_full.set_aspect("equal", adjustable="datalim")
    ax_full.set_xlabel("x (m)")
    ax_full.set_ylabel("r (m)")
    ax_full.set_title("full domain (equal aspect)")

    zoom_half_width = max(throat_r * 8, 1e-4)
    ax_zoom.add_collection(PolyCollection(
        [p for p in polys], facecolor="none", edgecolor="steelblue", linewidth=0.6
    ))
    ax_zoom.set_xlim(throat_x - zoom_half_width, throat_x + zoom_half_width)
    ax_zoom.set_ylim(0, throat_r * 3)
    ax_zoom.set_xlabel("x (m)")
    ax_zoom.set_ylabel("r (m)")
    ax_zoom.set_title("throat region zoom (non-equal aspect, for cell-quality spot-check)")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def mesh_geometry(csv_path: str, out_dir: str, geometry_id: str, cfg: MeshConfig,
                   json_path: Optional[str] = None) -> MeshResult:
    """
    Full single-geometry pipeline: read CSV -> build domain -> transfinite
    quad mesh -> physical groups -> export .msh (+ optional quality check
    and PNG preview). Initializes/finalizes its own gmsh session so it's
    safe to call in a loop without cross-geometry state leaking.

    Never raises for expected failure modes (bad geometry, degenerate
    mesh) -- returns a MeshResult with success=False and .error set, so
    mesh_sweep.py can log and continue rather than aborting the batch.
    """
    import os
    os.makedirs(out_dir, exist_ok=True)
    msh_path = os.path.join(out_dir, f"{geometry_id}.{ 'msh' }")
    png_path = os.path.join(out_dir, f"{geometry_id}.png")

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(geometry_id)

        points = read_geometry_csv(csv_path)
        domain = _build_domain(points, cfg)
        gmsh.model.geo.synchronize()

        gmsh.model.addPhysicalGroup(1, domain.wall_curves, name="wall")
        gmsh.model.addPhysicalGroup(1, domain.axis_curves, name="axis")
        gmsh.model.addPhysicalGroup(1, [domain.inlet_curve], name="inlet")
        gmsh.model.addPhysicalGroup(1, [domain.outlet_curve], name="outlet")
        gmsh.model.addPhysicalGroup(2, domain.surfaces, name="fluid")

        gmsh.model.mesh.generate(2)

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        _, elem_tags, _ = gmsh.model.mesh.getElements(2)
        n_elements = sum(len(t) for t in elem_tags)

        min_quality = None
        n_inverted = 0
        quality_flagged = False
        if cfg.check_quality and len(elem_tags) > 0:
            all_tags = [int(t) for tags in elem_tags for t in tags]
            qualities = gmsh.model.mesh.getElementQualities(all_tags, "minSICN")
            if len(qualities) > 0:
                min_quality = float(min(qualities))
                n_inverted = sum(1 for q in qualities if q <= 0)

            if n_inverted > 0:
                # Always fatal: an inverted/degenerate element means the mesh
                # is structurally broken, not merely a thin (valid)
                # boundary-layer cell -- distinct from the soft threshold
                # below. Fail the geometry rather than writing an unusable mesh.
                raise RuntimeError(
                    f"{n_inverted} inverted/degenerate element(s) (minSICN<=0); "
                    f"mesh is not usable"
                )

            if min_quality is not None and min_quality < cfg.quality_threshold:
                quality_flagged = True
                logger.warning(
                    "%s: mesh quality below soft threshold (informational -- expected for "
                    "tight boundary-layer cells, see DEFAULT_QUALITY_THRESHOLD note): "
                    "min minSICN=%.4g < %.4g",
                    geometry_id, min_quality, cfg.quality_threshold,
                )

        if cfg.mesh_format == "msh2":
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        elif cfg.mesh_format == "msh4":
            gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        else:
            raise ValueError(f"Unknown mesh_format {cfg.mesh_format!r}, expected 'msh2' or 'msh4'")
        gmsh.write(msh_path)

        if json_path and os.path.exists(json_path):
            dest_json = os.path.join(out_dir, f"{geometry_id}.json")
            if os.path.abspath(json_path) != os.path.abspath(dest_json):
                import shutil
                shutil.copy(json_path, dest_json)

        _plot_mesh_png(png_path, domain.wall_curves, title=f"{geometry_id} mesh ({n_elements} elements)")

        return MeshResult(
            geometry_id=geometry_id, success=True, msh_path=msh_path, png_path=png_path,
            n_nodes=len(node_tags), n_elements=n_elements,
            min_quality=min_quality, n_inverted_elements=n_inverted,
            quality_flagged=quality_flagged, error=None,
        )
    except Exception as e:
        logger.warning("%s: meshing failed: %s", geometry_id, e)
        return MeshResult(
            geometry_id=geometry_id, success=False, msh_path=None, png_path=None,
            n_nodes=0, n_elements=0, min_quality=None, n_inverted_elements=0,
            quality_flagged=False, error=str(e),
        )
    finally:
        gmsh.finalize()


# ----------------------------------------------------------------------------
# Regression test: mesh the generator's own straight-cone degeneracy case
# and confirm the pipeline works end-to-end before trusting it on real data.
# ----------------------------------------------------------------------------
def _regression_test_straight_cone_meshes_cleanly(tmp_dir: str = "/tmp/mesh_geometry_regression") -> bool:
    import os
    os.makedirs(tmp_dir, exist_ok=True)

    from nozzle_geometry import GeometryConfig, generate_geometry

    throat_r = 0.003
    exit_r = 0.010
    cp1_pos, cp2_pos, exit_pos = 0.008, 0.016, 0.040
    slope = (exit_r - throat_r) / exit_pos
    cfg = GeometryConfig(
        throat_radius=throat_r,
        throat_fillet_radius=1e-6,
        control_point_1_radius=throat_r + slope * cp1_pos,
        control_point_1_position=cp1_pos,
        control_point_2_radius=throat_r + slope * cp2_pos,
        control_point_2_position=cp2_pos,
        exit_radius=exit_r,
        exit_position=exit_pos,
        barrel_length=0.005,
        barrel_count=1,
        barrel_positions=[0.0],
        max_curvature=1e9,
    )
    geom_result = generate_geometry(cfg)
    if not geom_result.valid:
        logger.error("Regression test: fixture geometry itself is invalid: %s", geom_result.rejection_reason)
        return False

    csv_path = os.path.join(tmp_dir, "geometry_regtest.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["x", "r"])
        writer.writerows(geom_result.points)

    json_path = os.path.join(tmp_dir, "geometry_regtest.json")
    with open(json_path, "w") as f:
        json.dump(cfg.to_dict(), f, indent=2)

    mesh_cfg = MeshConfig(n_axial=80, n_radial=20, check_quality=True)
    result = mesh_geometry(csv_path, tmp_dir, "geometry_regtest", mesh_cfg, json_path=json_path)

    checks = {
        "mesh_generated": result.success,
        "msh_file_written": bool(result.msh_path) and os.path.exists(result.msh_path or ""),
        "png_written": bool(result.png_path) and os.path.exists(result.png_path or ""),
        "json_copied": os.path.exists(os.path.join(tmp_dir, "geometry_regtest.json")),
        "nonzero_nodes_elements": result.n_nodes > 0 and result.n_elements > 0,
        "quality_checked": result.min_quality is not None,
        # The real pass bar is "no inverted/degenerate elements" -- see
        # DEFAULT_QUALITY_THRESHOLD's note on why a tight, correctly-built
        # boundary-layer mesh can legitimately score below the soft
        # quality_threshold on minSICN without being broken.
        "no_inverted_elements": result.n_inverted_elements == 0,
    }

    physical_groups_ok = False
    if result.success:
        gmsh.initialize()
        try:
            gmsh.open(result.msh_path)
            names = set()
            for dim, tag in gmsh.model.getPhysicalGroups():
                names.add(gmsh.model.getPhysicalName(dim, tag))
            expected = {"wall", "axis", "inlet", "outlet", "fluid"}
            physical_groups_ok = expected.issubset(names)
            if not physical_groups_ok:
                logger.error("Regression test: physical groups %s missing from expected %s", names, expected)
        finally:
            gmsh.finalize()
    checks["physical_groups_named_and_present"] = physical_groups_ok

    passed = all(checks.values())
    for name, ok in checks.items():
        logger.info("  [%s] %s", "PASS" if ok else "FAIL", name)
    logger.info(
        "Regression test (straight-cone meshes cleanly): min_quality=%s -> %s",
        result.min_quality, "PASS" if passed else "FAIL",
    )
    return passed


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Single-geometry Gmsh meshing pipeline / regression test.")
    parser.add_argument("--csv", help="Path to a geometry_XXXX.csv to mesh (omit to just run the regression test).")
    parser.add_argument("--json", help="Path to the matching geometry_XXXX.json (copied into --out-dir).")
    parser.add_argument("--out-dir", default="mesh_output", help="Output directory for .msh/.png.")
    parser.add_argument("--n-axial", type=int, default=DEFAULT_N_AXIAL)
    parser.add_argument("--n-radial", type=int, default=DEFAULT_N_RADIAL)
    parser.add_argument("--first-cell-height", type=float, default=DEFAULT_FIRST_CELL_HEIGHT)
    parser.add_argument("--mesh-format", choices=["msh2", "msh4"], default=DEFAULT_MESH_FORMAT)
    parser.add_argument("--quality-threshold", type=float, default=DEFAULT_QUALITY_THRESHOLD)
    parser.add_argument("--mesh-quality-check", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not args.csv:
        ok = _regression_test_straight_cone_meshes_cleanly()
        print("Regression test:", "PASS" if ok else "FAIL")
        raise SystemExit(0 if ok else 1)

    import os
    geometry_id = os.path.splitext(os.path.basename(args.csv))[0]
    mesh_cfg = MeshConfig(
        n_axial=args.n_axial, n_radial=args.n_radial, first_cell_height=args.first_cell_height,
        mesh_format=args.mesh_format, quality_threshold=args.quality_threshold,
        check_quality=args.mesh_quality_check,
    )
    result = mesh_geometry(args.csv, args.out_dir, geometry_id, mesh_cfg, json_path=args.json)
    print(result)
