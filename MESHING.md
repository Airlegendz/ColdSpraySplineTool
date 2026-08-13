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

- `--n-axial` (default 150): transfinite node count along the wall/axis
  (flow direction), split proportionally between the two throat blocks.
- `--n-radial` (default 40): transfinite node count along inlet/interface/
  outlet (radial direction).
- `--first-cell-height` (default 5e-6 m): target wall-adjacent cell height
  for boundary-layer clustering. The pipeline solves for the geometric
  growth ratio that achieves this over `n_radial` cells (falls back to
  `--bl-growth-ratio-fallback`, default 1.2, if the target height is
  infeasible for the given curve length/cell count).
- `--quality-threshold` (default 0.02, only used with
  `--mesh-quality-check`): **soft/informational** threshold on minSICN
  (Gmsh's signed inverted condition number, a scaled-Jacobian-like 0..1
  quality measure). Read this carefully before changing it: elements with
  quality <= 0 are inverted/degenerate and are **always** treated as a hard
  failure regardless of this threshold — a genuinely broken mesh. Elements
  scoring between 0 and this threshold are only **warned** about, because
  a correctly-built boundary-layer mesh legitimately produces very thin,
  high-aspect-ratio cells at the wall, and minSICN penalizes aspect ratio
  as well as skew — a tight, otherwise-correct first cell can score as low
  as ~0.003 (observed on the regression fixture) without being wrong. This
  default is set low enough to avoid flagging that expected case. Fluent
  reports orthogonal quality / aspect ratio as separate metrics, which are
  more meaningful for assessing boundary-layer cells specifically than a
  single blended quality number — worth checking there too before trusting
  a borderline mesh.

## `--mesh-quality-check`

Computes `minSICN` over all elements after meshing. Always fails the
geometry (marks it unmeshed, with an error) if any element is inverted or
degenerate (`quality <= 0`). Otherwise warns (does not fail) if the minimum
quality falls below `--quality-threshold` — see the note above on why that's
deliberately a soft check for this kind of mesh.

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
  actual Reynolds-number or y+ target — confirm against your flow
  conditions before a production meshing run, same caveat as the geometry
  generator's own unconfirmed defaults.
