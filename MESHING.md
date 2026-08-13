# Meshing pipeline: nozzle geometry CSVs -> Fluent-ready `.msh`

Consumes `geometry_XXXX.csv` files produced by `nozzle_geometry.py`/`sweep.py`
and produces a 2D axisymmetric structured (quad, transfinite) Gmsh mesh per
geometry, exported to Gmsh's native `.msh` format for Fluent import.

## Files

- `mesh_geometry.py` — core single-geometry pipeline: read CSV -> build
  domain -> transfinite quad mesh -> physical groups -> export `.msh` +
  mesh-preview `.png`, plus a regression test
  (`python3 mesh_geometry.py`).
- `mesh_sweep.py` — batch driver: globs `geometry_*.csv` in a `sweep.py`
  output directory, meshes each one, writes a `mesh_index.csv` and a
  success/failure/quality-flagged summary — analogous to `sweep.py`'s
  rejection-rate reporting, but for meshing, which is a **separate failure
  mode from geometry validity**: a geometry can pass every
  `nozzle_geometry.py` constraint check and still mesh poorly (or fail to
  mesh) at the given mesh-density settings.

## Requirements

```bash
pip install gmsh
```

## Running it

Regression test (meshes the generator's own straight-cone degeneracy case,
same fixture family as `nozzle_geometry.py`'s regression test, and confirms
the pipeline works end-to-end):
```bash
python3 mesh_geometry.py
```

Mesh one geometry manually:
```bash
python3 mesh_geometry.py --csv sweep_output_v1/geometry_0001.csv \
  --json sweep_output_v1/geometry_0001.json --out-dir mesh_output --mesh-quality-check
```

Batch-mesh an entire sweep output directory:
```bash
python3 mesh_sweep.py --sweep-dir sweep_output_v1 --out-dir mesh_output_v1 --mesh-quality-check
```

Outputs in `mesh_output_v1/`, one triple per successfully-meshed geometry:
- `geometry_XXXX.msh` — the mesh, Fluent-importable.
- `geometry_XXXX.json` — copied from the sweep output, so every mesh stays
  traceable back to its exact source parameters.
- `geometry_XXXX.png` — two-panel spot-check image: a full-domain overview
  (equal aspect — necessarily looks "squashed" for a long, thin nozzle;
  that's true physical proportions, not a bug) and a zoomed throat-region
  panel at non-equal aspect, which is what actually lets you see individual
  cell skew before trusting the mesh to Fluent.
- `mesh_index.csv` — one row per attempted geometry: success, quality
  metrics, node/element counts, errors.

The batch summary log reports total meshed / quality-flagged / failed, and
warns if the failure rate exceeds 20% (same pattern as `sweep.py`'s
rejection-rate warning).

## Domain construction: two transfinite blocks split at the throat

The 2D axisymmetric domain is the wall (top), the centerline `r=0` (bottom),
closed by an inlet edge at `x=0` and an outlet edge at the wall's final `x`.
**It's built as two transfinite blocks, split at the throat** (the wall's
minimum-radius point), not one single 4-sided block spanning the whole
nozzle. This split is load-bearing, not cosmetic:

A single-block transfinite (bilinear/Coons) surface blends its radial node
distribution only between the inlet and outlet boundary conditions. For a
converging-diverging "necked" domain like this nozzle, that blend does not
locally adapt to the much smaller wall radius at the throat. This was
diagnosed directly during development: even with boundary-layer clustering
disabled entirely (uniform radial spacing), a single-block mesh still showed
pronounced radial-line shear through the throat region, visible in the mesh
preview's zoomed panel. Splitting into two blocks at the throat — with a
third, shared "interface" curve graded from the throat's own local wall
radius — fixed it directly: minimum mesh quality (minSICN) improved roughly
19x on the geometry used to diagnose this (0.00018 -> 0.0035), with the
radial mesh lines visibly near-vertical through the throat afterward. Both
blocks share the same interface curve and throat point (not two
coincident-but-separate points), so the mesh is conformal across the split —
Fluent sees one continuous fluid zone.

Each wall segment's underlying spline geometry is also resampled to be
uniform in `x` (`_resample_uniform_in_x`) before building the curve, so its
transfinite node distribution stays aligned in `x` with the (straight, and
therefore naturally `x`-uniform) axis segment it's paired with.

## Physical groups (for Fluent boundary-condition assignment on import)

| Physical group | Type | Fluent boundary type |
|---|---|---|
| `wall` | 1D | Wall |
| `axis` | 1D | **Axis** (see import steps below — must be set manually) |
| `inlet` | 1D | Velocity inlet / pressure inlet / mass-flow inlet (your choice) |
| `outlet` | 1D | Pressure outlet / outflow (your choice) |
| `fluid` | 2D | (interior — the flow domain) |

## Exact Fluent import steps

1. `File -> Import -> Mesh...`, select the `.msh` file.
2. `General -> Solver -> 2D Space -> Axisymmetric`. This is required —
   without it, Fluent treats the domain as a planar 2D slice, not a
   revolved axisymmetric one, and the `axis` boundary has no special
   meaning.
3. `Boundary Conditions`: Fluent will have imported `wall`, `axis`, `inlet`,
   `outlet` as named zones from the physical groups above, but **`axis` is
   very likely still typed as `wall` on import** — Gmsh physical groups
   don't carry the "this is degenerate axisymmetric axis" boundary-type
   information Fluent needs. Select the `axis` zone and explicitly change
   its **Type** to **Axis**. Getting this wrong (leaving it as `wall`) will
   give you a solid centerline instead of a symmetry axis and visibly wrong
   results near `r=0`.
4. Assign `inlet`/`outlet` boundary condition types per your actual flow
   setup (velocity/pressure/mass-flow inlet, pressure outlet, etc. — not
   prescribed by this pipeline).

## Mesh format

Default export is **MSH2** (`Mesh.MshFileVersion = 2.2`), selectable via
`--mesh-format {msh2,msh4}`. **UNCONFIRMED which your specific Fluent
version expects** — older Fluent releases are known to be picky about MSH2
specifically vs. the newer MSH4 format; this default favors the more
broadly-compatible older format, but confirm against your actual Fluent
version before a full batch import, the same way `nozzle_geometry.py`
flags its own unconfirmed defaults.

## Configurable mesh-density parameters (all UNCONFIRMED defaults)

Same status as `nozzle_geometry.py`'s `max_curvature`/`max_fillet_fraction`:
exposed as configurable parameters with a logged default, not validated
engineering values. The right values depend on your actual flow Reynolds
number and haven't been established — confirm before a production meshing
run.

- `--n-axial` (default **3000**): transfinite node count along the
  wall/axis (flow direction), split between the two throat blocks
  proportional to each block's share of the CSV's point count (see
  `DEFAULT_N_AXIAL`'s code comment for why that split turned out fine —
  it's actually *more* generous to the fixed-length convergent block than
  a pure length-proportional split would be, not less). This default is
  much higher than an initial guess of 150 for a real, diagnosed reason —
  see "Quality-metric note" below.
- `--n-radial` (default 40): transfinite node count along inlet/interface/
  outlet (radial direction).
- `--first-cell-height` (default 5e-6 m): target wall-adjacent cell height
  for boundary-layer clustering. The pipeline solves for the geometric
  growth ratio that achieves this over `n_radial` cells (falls back to
  `--bl-growth-ratio-fallback`, default 1.2, if the target height is
  infeasible for the given curve length/cell count). Checked as a possible
  second contributor to the quality issue below (it's a fixed *absolute*
  value regardless of local throat scale) — for the actual 80-sample sweep
  this was tested against, `throat_radius` only varies ~1.27x across the
  batch, so `first_cell_height`'s relative scale barely varies either; not
  a live problem *for this batch's bounds*, but the same class of
  fragility as `n_axial` was, and worth revisiting if `throat_radius`
  bounds are ever widened substantially in a future sweep.
- `--quality-threshold` (default **0.05**, only used with
  `--mesh-quality-check`): **soft/informational** threshold on minSICN
  (Gmsh's signed inverted condition number, a scaled-Jacobian-like 0..1
  quality measure). Elements with quality <= 0 are inverted/degenerate and
  are **always** treated as a hard failure regardless of this threshold —
  a genuinely broken mesh. Elements scoring between 0 and this threshold
  are only **warned** about — see the quality-metric note below for what
  that's calibrated against now.

## `--mesh-quality-check`

Computes `minSICN` over all elements after meshing. Always fails the
geometry (marks it unmeshed, with an error) if any element is inverted or
degenerate (`quality <= 0`). Otherwise warns (does not fail) if the minimum
quality falls below `--quality-threshold`.

### Quality-metric note: core mesh vs. boundary-layer cells (corrected)

An earlier version of this document attributed low overall `minSICN`
readings (~0.003, observed with an initial `n_axial=150` default) entirely
to intentional boundary-layer cell thinness at the wall. **That was wrong**,
and is corrected here rather than left as-is.

Verified directly by peeling the mesh into structured layers outward from
the wall (by node adjacency) and tracking `minSICN` per layer: at
`n_axial=150`, quality climbed only gradually across *nearly the entire
radial extent* — even the layer immediately next to the axis (as far from
the wall as the mesh gets) only reached ~0.1, not the comfortably-high
value a "just the BL row" explanation would predict. The real cause: axial
cell width (domain length / `n_axial`) was far larger than the throat's
radial extent for a slender nozzle (e.g. ~1.8mm axial cells against a
~1.1mm-radius throat), so *nearly every cell in the mesh*, not just the
thin first row at the wall, had a poor (elongated) aspect ratio. This
affected every geometry in a real 80-sample sweep to varying degrees,
which is exactly why all of them showed similarly low minimums rather than
one obvious outlier.

Fixed by raising `n_axial` to 3000 (see that default's code comment for
the derivation — solved for the worst-case, i.e. smallest, `throat_radius`
in the batch, then checked it wasn't excessive for the largest-radius and
largest-`exit_position` geometries too). Re-verified with the same
layer-peeling method after the fix, on the geometry that had the batch's
worst overall reading:

| layer from wall | min quality | max quality | mean quality |
|---|---|---|---|
| 0 (wall-adjacent) | 0.055 | 0.122 | 0.078 |
| 3 | 0.077 | 0.179 | 0.109 |
| 8 | 0.133 | 0.343 | 0.192 |
| 12 | 0.207 | 0.560 | 0.300 |
| 20 | 0.425 | 0.998 | 0.651 |
| 38 (axis-adjacent) | 0.160 | 1.000 | 0.578 |

That's the shape a genuine "BL-row-only" quality dip should look like: a
fast, clean climb through the first ~10-20 layers, not a slow crawl across
the whole domain. Re-ran the full 72-geometry batch after the fix: **72/72
meshed, 0 failed, 0 inverted elements, quality range 0.055-0.089**, all
above the (now correspondingly raised) 0.05 soft threshold — versus
0.0027-0.0043 (all flagged) before the fix.

The remaining low points — the wall-adjacent row (thin by design) and a
narrow dip right at the axis-adjacent layer — are real and expected: the
first is deliberate BL clustering, working as intended; the second appears
tied to how the two throat blocks' independently-solved radial growth
ratios (each based on its own local wall radius) meet at the shared axis
point, and is a much smaller, localized effect worth a closer look later
but not the systemic problem the pre-fix numbers indicated. Fluent's own
orthogonal-quality / aspect-ratio diagnostics remain the more meaningful
check for boundary-layer cells specifically, and are still worth running
before a full solve.

## Known limitations / follow-ups

- The two-block throat split handles the single-throat convergent-divergent
  shape this generator produces. A geometry with multiple radius minima
  (not currently producible by `nozzle_geometry.py`) would need a more
  general multi-block split to get the same benefit.
- Mesh-quality assessment here is a single blended metric (minSICN); it's a
  useful automated gate for "is this mesh usable at all" (no inversions)
  and a rough informational signal, but isn't a substitute for checking
  Fluent's own mesh-quality diagnostics (orthogonal quality, aspect ratio,
  skewness) before a full solve.
- `first_cell_height`/`n_radial`/`n_axial` defaults are not derived from any
  actual Reynolds-number or y+ target — `n_axial` is now at least
  geometrically self-consistent (derived from matching axial and radial
  cell scales at the throat, see above), but none of the three have been
  validated against real flow physics. Confirm against your flow
  conditions before a production meshing run, same caveat as the geometry
  generator's own unconfirmed defaults.
- `first_cell_height` is a fixed absolute value, not scaled to local
  geometry (e.g. `throat_radius`). Checked and found not to be a live
  problem for the current sweep's narrow `throat_radius` range (~1.27x
  spread), but this is the same class of fragility that made the old
  `n_axial=150` default silently wrong for a wider geometry range — revisit
  if `throat_radius` bounds are ever widened substantially.
