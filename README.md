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
   with zero slope at both ends (both are configurable, not swept).
2. **Throat fillet**: a circular arc of radius `throat_fillet_radius`,
   tangent to the flat convergent wall and blending into the divergent
   section's initial slope.
3. **Divergent**: a piecewise-cubic (Catmull-Rom-derived Bezier) spline,
   C1-continuous, passing exactly through the throat point, control point 1,
   control point 2, and the exit point, in that order. **Design choice**: a
   true single cubic Bezier (4 control points) does not generally pass
   through its two interior control points — only through its endpoints —
   so a clamped Catmull-Rom-to-Bezier spline was used instead to satisfy the
   requirement that the curve pass through `control_point_1` and
   `control_point_2` exactly. See the module docstring in
   `nozzle_geometry.py` for the full rationale. The exit axial position is
   derived as `control_point_2_position + 1.5 * (control_point_2_position -
   control_point_1_position)` — also a documented default, not a swept
   parameter.
4. **Barrel(s)**: `barrel_count` straight, constant-radius sections at the
   given `barrel_positions` (axial offsets from the divergent section's
   exit), each of length `barrel_length`. Skipped entirely if
   `barrel_count = 0`.

## Constraints enforced (reject, don't silently mesh)

- **Monotonicity**: wall radius must be non-decreasing from throat to exit.
- **Curvature cap**: `|d²r/dx²|` is capped at `max_curvature` (finite-
  difference proxy), a **configurable, currently UNCONFIRMED** threshold
  (default `800.0` 1/m in `config_example.yaml`, `500.0` in the code
  default) — a stand-in for flow-separation risk and manufacturability,
  not a value derived from CFD or experiment yet.
- **Positive lengths**: `control_point_1_position < control_point_2_position
  < exit_position`, positive `convergent_length`, positive `barrel_length`
  when `barrel_count > 0`.
- **Fillet radius bound**: `throat_fillet_radius <= max_fillet_fraction *
  throat_radius`, with `max_fillet_fraction` **configurable and
  UNCONFIRMED** (default `0.5`–`0.6` depending on config).
- **Barrel non-overlap**: when `barrel_count > 1`, barrel intervals
  `[position, position + barrel_length]` must not overlap each other.

Every rejection is logged with the specific reason (which check failed and
the offending value), not dropped silently — this is what
`rejection_reason` in the index CSV and the sweep log records.

### ⚠️ Values that need confirmation before a full sweep

The following are implemented as configurable parameters with a logged
default, **not** as validated engineering values:

- `max_curvature` (curvature/smoothness cap)
- `max_fillet_fraction` (fillet-radius-to-throat-radius cap)
- `inlet_radius`, `convergent_length` (fixed convergent-section shape)
- The divergent-section exit-position extension factor (`1.5x` the
  `cp1→cp2` spacing)

Confirm these against the actual CFD/experimental setup before committing
to a full production sweep — a first-pass `config_example.yaml` sweep at
these defaults sees a ~50% rejection rate, which is itself a signal the
bounds (or thresholds) need tightening, as flagged in the sweep log.

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
