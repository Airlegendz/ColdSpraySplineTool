# Cold Spray Nozzle Geometry Generator

Generates 2D axisymmetric wall profiles for converging-diverging cold-spray
nozzles from spline-based design parameters, for batch feeding into Gmsh (or
any other CFD preprocessor) via a Latin hypercube sweep.

## Files

- `nozzle_geometry.py` — core generator: convergent + fillet + divergent
  spline + barrel(s), constraint checks, and a hardcoded degeneracy
  regression test (`python3 nozzle_geometry.py --plot`).
- `sweep.py` — batch driver: reads a YAML bounds config, draws a Latin
  hypercube sample, generates + validates each geometry, writes per-geometry
  CSV/JSON and a master index CSV.
- `plotting.py` — shared matplotlib helper (mirrored wall plot to PNG).
- `config_example.yaml` — example parameter bounds.

## Coordinate convention

- Units: **meters** throughout (radii and axial positions).
- `x` is the axial coordinate, increasing in the flow direction.
- **`x = 0` is the nozzle inlet** — the start of the fixed convergent
  section.
- `r` is the wall radius from the centerline (the profile is the upper
  wall trace only; mirror about `r = 0` for the full cross-section, which
  is what `--plot` renders).
- Each output `geometry_XXXX.csv` is an ordered list of `(x, r)` points
  tracing the wall from inlet to the end of the last barrel section, dense
  enough (≥50 points, default ~200) across the divergent section for the
  mesher to reconstruct a smooth curve via spline/polyline import.

## Geometry sections

1. **Convergent** (fixed, not parameterized): a raised-cosine contraction
   from a fixed `inlet_radius` to `throat_radius` over `convergent_length`,
   with zero slope at both ends (both are configurable, not swept). Zero
   slope at the throat end is required, not incidental: the throat is by
   definition the point of minimum radius, and any smooth curve has zero
   slope at a true minimum. See "Throat transition smoothness" below.
2. **Throat fillet**: a circular arc of radius `throat_fillet_radius`,
   tangent to the flat convergent wall and blending into the divergent
   section's initial slope.
3. **Divergent**: a piecewise-cubic Hermite spline (evaluated via its Bezier
   form), C1-continuous, passing exactly through the throat point, control
   point 1, control point 2, and the exit point, in that order. **Design
   choice**: a true single cubic Bezier (4 control points) does not
   generally pass through its two interior control points — only through
   its endpoints — so a piecewise Hermite/Bezier spline was used instead to
   satisfy the requirement that the curve pass through `control_point_1`
   and `control_point_2` exactly. Tangents are *clamped* at both ends
   rather than extrapolated from a phantom point: the start tangent matches
   the throat fillet's actual exit slope (fillet → spline is C1), and the
   end tangent is forced horizontal when a barrel follows (spline → barrel
   is C1) or left as the natural last-chord slope when there is no barrel
   (so the spline degenerates exactly to a straight line when all four
   knots are colinear — see the regression test). See the module docstring
   in `nozzle_geometry.py` for the full rationale.

   `exit_position` (axial distance from the throat to the exit point) is an
   **independent, directly-swept parameter** — it is *not* derived from
   `control_point_2_position`. This was a deliberate fix: an earlier
   version derived it as `control_point_2_position + 1.5 * (cp2_position -
   cp1_position)`, which silently locked overall divergent length to the
   cp1↔cp2 spacing and removed the sweep's ability to explore divergent
   length independently — a problem because Badali et al. found divergent
   length to be the single most influential parameter for particle
   velocity. Bounds must satisfy `control_point_1_position <
   control_point_2_position < exit_position`, enforced by the
   `positive_lengths` constraint check.
4. **Barrel(s)**: `barrel_count` straight, constant-radius sections at the
   given `barrel_positions` (axial offsets from the divergent section's
   exit), each of length `barrel_length`. Skipped entirely if
   `barrel_count = 0`.

### Throat transition smoothness

Zoomed plots of the convergent → fillet → divergent-spline transition show
a visually flat-looking run of a few hundred microns to ~1mm around the
throat (reproduce with `nozzle_geometry.py`'s throat-zoom check). This is
**expected, C1/C2-continuous behavior, not a defect**: the throat is
defined as the wall's minimum radius, so `dr/dx = 0` there by construction,
and any smooth curve necessarily flattens as it approaches a true minimum
from both sides — independent of `throat_fillet_radius` (verified across
fillet radii from ~0.0008–0.0024m; the flat-looking span is essentially
unchanged). Forcing the convergent section to retain a nonzero slope at the
throat, as an alternative fix, was considered and rejected: it would mean
the wall is still narrowing as it enters the fillet, which the
`monotonicity` constraint (radius non-decreasing from throat to exit) would
then have to reject. No slope discontinuity exists at either the
convergent→fillet or fillet→divergent-spline joins — confirmed numerically
(finite-difference slope is continuous to within sampling noise, ~4e-4,
across both joins).

## Constraints enforced (reject, don't silently mesh)

Each rejection reason is tagged with a stable `[category]` prefix (e.g.
`[curvature] ...`) so `sweep.py` can report a clean breakdown by constraint
type rather than by raw, value-specific message text.

- **`[monotonicity]`**: wall radius must be non-decreasing from throat to
  exit.
- **`[curvature]`**: `|d²r/dx²|` is capped at `max_curvature` (finite-
  difference proxy) — a stand-in for flow-separation risk and
  manufacturability, not a value derived from CFD or experiment yet.
- **`[positive_lengths]`**: `control_point_1_position <
  control_point_2_position < exit_position`, positive `convergent_length`,
  positive `barrel_length` when `barrel_count > 0`.
- **`[fillet_bound]`**: `throat_fillet_radius <= max_fillet_fraction *
  throat_radius`.
- **`[fillet_extent]`**: the fillet arc's axial (x) extent must stay
  strictly less than `control_point_1_position` — otherwise the fillet
  geometrically overruns control point 1. This can happen for a steep
  target slope (`(control_point_1_radius - throat_radius) /
  control_point_1_position`, see the fillet-target-slope note below)
  combined with a small `control_point_1_position`, and was previously
  unconstrained.
- **`[barrel_overlap]`**: when `barrel_count > 1`, barrel intervals
  `[position, position + barrel_length]` must not overlap each other.

Every rejection is logged with the specific reason (which check failed and
the offending value), not dropped silently — this is what
`rejection_reason` in the index CSV and the sweep log records.

### ⚠️ Values that need confirmation before a full sweep

`max_curvature` and `max_fillet_fraction` are implemented as configurable
parameters with **one authoritative default each**, not as validated
engineering values:

- `max_curvature = 800.0` (1/m)
- `max_fillet_fraction = 0.6`

Both defaults live in `nozzle_geometry.py` (`DEFAULT_MAX_CURVATURE`,
`DEFAULT_MAX_FILLET_FRACTION`) and are duplicated in `config_example.yaml`'s
`fixed` section with the same values, so there is no drift between "the
code default" and "the example sweep default" — they're intentionally the
same UNCONFIRMED number. **Override behavior**: when driving generation
through `sweep.py`, any threshold set in the YAML config always wins —
`sweep.py` forwards every key present in the config as an explicit
`GeometryConfig` keyword argument, which overrides the dataclass field
default. The code-level default only applies if you construct
`GeometryConfig` directly (e.g. ad hoc scripts, tests) without specifying
the field.

Also unconfirmed: `inlet_radius`, `convergent_length` (fixed
convergent-section shape), and the **fillet target-slope heuristic** in
`_fillet_arc`/`_fillet_target_slope` — the fillet's sweep angle is set from
the straight-chord slope `(control_point_1_radius - throat_radius) /
control_point_1_position`, a real modeling decision (a different heuristic,
e.g. matching the divergent spline's actual solved tangent, would produce a
different fillet shape), not a validated one.

Confirm these against the actual CFD/experimental setup before committing
to a full production sweep. With the current `config_example.yaml` bounds
and the defaults above, a 200-sample sweep (`--seed 1`) sees a **38.5%
rejection rate** (123/200 valid), broken down as:

| category | count |
|---|---|
| `fillet_bound` | 38 |
| `monotonicity` | 26 |
| `curvature` | 13 |
| `fillet_extent` | 0 |

`fillet_extent` (the fillet-overrunning-control_point_1 check) catches
**zero** rejections at the current bounds — `throat_fillet_radius` ∈
`[0.0015, 0.0025]` and `control_point_1_position` ∈ `[0.004, 0.010]` keep
the fillet's x-extent well under `control_point_1_position` in practice, so
this is currently a theoretical edge case rather than a live problem. It's
still enforced because a future bounds change (smaller
`control_point_1_position` or larger fillet radii) could trigger it
silently otherwise.

This breakdown was re-measured *after* making `exit_position` an
independent parameter (see above) — the earlier ~50% figure included
artifacts from the old derived exit-position formula compressing the
divergent spline into a tighter-than-intended span, which spuriously
tripped the curvature cap. `fillet_bound` is now the largest bucket:
`throat_fillet_radius` bounds `[0.0015, 0.0025]` and `throat_radius` bounds
`[0.003, 0.005]` combined with `max_fillet_fraction = 0.6` reject any
sample where `throat_fillet_radius > 0.6 * throat_radius` — this is a
genuine bounds/threshold interaction, not an artifact, and is the next
thing to tighten once real fillet-fraction guidance is available.

## Usage

Regression test (degenerate near-straight-cone case):
```bash
python3 nozzle_geometry.py --plot
```

Batch sweep:
```bash
python3 sweep.py --config config_example.yaml --n-samples 200 --seed 42 --out-dir sweep_output --plot
```

Outputs in `sweep_output/`:
- `geometry_XXXX.csv` — `x, r` columns for each **valid** geometry.
- `geometry_XXXX.json` — exact parameter values used for that geometry.
- `geometry_XXXX.png` — mirrored wall-profile plot (only with `--plot`).
- `index.csv` — one row per sampled geometry (valid and rejected), with all
  parameters, `valid`, and `rejection_reason`.

The sweep log reports the overall rejection rate and a breakdown by failure
reason at the end of the run.

## Handing off to Gmsh / CFD preprocessing

Each `geometry_XXXX.csv` is a simple ordered `(x, r)` polyline of the upper
wall, in meters, with `x = 0` at the inlet. To build an axisymmetric CFD
mesh:

1. Import the CSV points as an ordered spline/polyline in Gmsh (e.g. via a
   `.geo` script reading the CSV and creating a `Spline` or `BSpline`
   through the points).
2. Revolve the 2D profile (or mesh the 2D half-profile directly for a
   quasi-2D axisymmetric solver) about the `x`-axis (`r = 0`) to form the
   3D/axisymmetric domain.
3. Close the domain with an inlet plane at `x = 0`, an outlet plane at the
   last barrel section's end, and the centerline (`r = 0`) as the axis of
   symmetry.

No units conversion is needed if your CFD preprocessor also expects meters;
otherwise scale the CSV columns accordingly before import.
