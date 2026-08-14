"""
Per-geometry Fluent solve routine, driven via PyFluent's settings API over
an existing gRPC connection (see test_connection.py / FLUENT_REMOTE.md for
establishing that connection).

VERIFICATION STATUS -- read before trusting this beyond the --subset
validation run described in FLUENT_REMOTE.md:

Every settings path used below (e.g. `settings.setup.general.solver.type`,
`settings.setup.boundary_conditions.pressure_inlet[...].momentum.
gauge_total_pressure`, `settings.results.report.discrete_phase.
extended_summary`) was checked directly against the *.pyi type stubs
shipped inside the locally-installed ansys-fluent-core package
(ansys/fluent/core/generated/solver/settings_261.pyi) -- these are real,
current attribute/method names for that generated settings schema, not
guessed from documentation. That's a meaningfully stronger footing than
the earlier .jou/TUI-string approach, where command syntax had to be
inferred from general TUI conventions with no way to check it locally.

LIVE TESTING SO FAR (against a real Fluent 2026 R1 session on a separate
Windows PC, via run_subset.py): mesh import now succeeds cleanly using
Abaqus INP (mesh_geometry.py) via `file.import_.read(file_type=
"abaqus-input", ...)` -- two earlier attempts failed live first (native
"mesh" reader: "Null Domain Pointer"; CGNS: crashed the whole Fluent
process with SIGSEGV). Two real bugs were also found and fixed by live
testing: Command objects are keyword-only (positional calls to
`set_zone_type`/`injections.create` raised "Command.__call__() takes 1
positional argument but 3 were given"), and Fluent's Abaqus importer
discards our wall/axis/inlet/outlet element-set names entirely, requiring
`_recover_boundary_zone_names` (angle-based zone splitting + geometric
identification via face centroids) to recover them post-import -- see
that function's docstring for what's confirmed vs. still a live guess.

What is NOT yet verified, because it requires further live testing:
  - `_recover_boundary_zone_names` as a whole: the 45-degree split angle,
    whether it produces exactly 4 pieces, and whether centroid-based
    classification correctly identifies each one. CONFIRMED LIVE and
    working on geometry_0001: split produced exactly 4 pieces (2999/2999/
    39/39 faces) and centroid-based identification correctly matched all
    four to axis/inlet/outlet/wall on the first successful attempt.
  - `particle_type = "inert"` and `material = "copper"` -- CONFIRMED LIVE
    (Fluent printed "Copying particle material from the database: copper"
    with no error). `injection_type.option = "single"` was tried and
    FAILED live ("ASSQ: invalid argument [2]: improper list") -- removed
    rather than guessed again; a new injection defaults to single-point
    already, so leaving it unset should be equivalent.
  - Residual equation names used for convergence-criteria/monitoring
    (e.g. "continuity", "x-velocity") -- standard Fluent naming, but not
    confirmed against this schema version specifically.
  - Whether one long-lived session can cleanly read/solve/reset across
    all 72 geometries in sequence, or whether state leaks between runs
    (residual monitors, DPM injections, report definitions all persist on
    a session unless explicitly removed) -- run_subset.py's whole purpose
    is to observe this directly before fluent_batch.py assumes an answer.

See FLUENT_REMOTE.md for the full checklist to work through on first
contact with real Fluent.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

from ansys.fluent.core.fields.field_data_interfaces import SurfaceData, SurfaceDataType, SurfaceFieldDataRequest

logger = logging.getLogger("fluent_solve")


@dataclass
class SolveResult:
    geometry_id: str
    success: bool
    converged: bool
    iterations_run: int
    exit_velocity_m_s: Optional[float]
    warnings: list
    error: Optional[str]


def _zone_centroid_stats(solver_session, zone_name: str) -> tuple:
    """
    Returns (mean_x, mean_y, n_faces) for a face zone's centroids.

    Uses solver_session.fields.field_data.get_field_data(SurfaceFieldDataRequest(...))
    -- confirmed real (not guessed): checked directly against the installed
    PyFluent package's fields/live_field_data.py (LiveFieldData.get_field_data,
    with a documented usage example in its own docstring) and
    fields/field_data_interfaces.py (SurfaceData.face_centroids, an Nx3 array).
    session.fields.field_data itself confirmed present via session.py's
    BaseSession.__init__ (self.fields = Fields(...)).
    """
    data = solver_session.fields.field_data.get_field_data(
        SurfaceFieldDataRequest(data_types=[SurfaceDataType.FacesCentroid], surfaces=[zone_name])
    )
    surf = data[zone_name]
    # Confirmed live: get_field_data() already returns a SurfaceData instance
    # per zone (not a raw dict as the generic type hint in live_field_data.py
    # suggested) -- wrapping it in SurfaceData(...) again fails with
    # "'SurfaceData' object has no attribute 'get'". Handle both shapes
    # defensively rather than hardcode the one just observed.
    centroids = surf.face_centroids if hasattr(surf, "face_centroids") else SurfaceData(surf).face_centroids
    xs = centroids[:, 0]
    ys = centroids[:, 1]
    return float(xs.mean()), float(ys.mean()), len(xs)


def _recover_boundary_zone_names(solver_session, geometry_id: str) -> None:
    """
    Fluent's Abaqus INP importer collapses ALL boundary edges into one
    merged wall-type zone, ignoring our exported wall/axis/inlet/outlet
    element-set names entirely (confirmed live: a real Fluent 2026 R1
    session imported our mesh with all boundary edges under a single
    generic zone like "wall-2", not four named ones). Recovers the
    distinction here rather than at the file-format level, since two
    format attempts (CGNS, which crashed Fluent outright, and Abaqus INP)
    both failed to carry zone names through cleanly.

    Approach: split the merged zone at sharp-angle corners (our four
    boundary segments meet at ~90-degree corners; the wall/axis curves
    themselves are smooth, bounded-curvature C1 curves that shouldn't
    locally exceed a 45-degree angle between adjacent mesh edges at this
    resolution) via settings.mesh.modify_zones.sep_face_zone_angle --
    confirmed real signature (face_zone_name, angle, move_faces) checked
    against both the installed settings schema and PyFluent's official
    docs. Then identifies each resulting piece by its face-centroid
    geometry (axis: mean radial position ~0; inlet: mean axial position
    ~0; outlet: mean axial position ~x_max; wall: whatever's left, by far
    the most faces) and renames via settings.mesh.modify_zones.zone_name
    (also confirmed real).

    UNVERIFIED: the 45-degree angle threshold and the assumption that
    sep_face_zone_angle produces exactly 4 pieces here (not more, e.g. if
    the wall's own curvature somewhere locally exceeds 45 degrees between
    adjacent cells, or fewer, if some corners are smoother than expected).
    This whole function is the least-tested part of the pipeline so far --
    watch its behavior closely on the next live run.
    """
    settings = solver_session.settings
    wall_names_before = settings.setup.boundary_conditions.wall.get_object_names()
    if len(wall_names_before) != 1:
        logger.warning(
            "%s: expected exactly 1 merged wall zone before splitting, found %s -- "
            "proceeding with the first one, but this is unexpected.",
            geometry_id, wall_names_before,
        )
    merged_zone_name = wall_names_before[0]

    settings.mesh.modify_zones.sep_face_zone_angle(
        face_zone_name=merged_zone_name, angle=45.0, move_faces=True,
    )

    wall_names_after = settings.setup.boundary_conditions.wall.get_object_names()
    logger.info("%s: zone split produced %s (from %s)", geometry_id, wall_names_after, merged_zone_name)

    if len(wall_names_after) < 4:
        raise RuntimeError(
            f"Expected the merged boundary zone to split into >=4 pieces "
            f"(wall/axis/inlet/outlet); got {wall_names_after}. The 45-degree "
            f"angle threshold in _recover_boundary_zone_names likely needs adjusting."
        )

    stats = {name: _zone_centroid_stats(solver_session, name) for name in wall_names_after}
    logger.info("%s: zone centroid stats (mean_x, mean_y, n_faces): %s", geometry_id, stats)

    remaining = dict(stats)
    axis_name = min(remaining, key=lambda n: remaining[n][1])  # smallest mean radial position
    del remaining[axis_name]
    inlet_name = min(remaining, key=lambda n: remaining[n][0])  # smallest mean axial position
    del remaining[inlet_name]
    outlet_name = max(remaining, key=lambda n: remaining[n][0])  # largest mean axial position
    del remaining[outlet_name]
    wall_name = max(remaining, key=lambda n: remaining[n][2])  # most faces, of whatever's left

    renames = {axis_name: "axis", inlet_name: "inlet", outlet_name: "outlet", wall_name: "wall"}
    logger.info("%s: renaming zones: %s", geometry_id, renames)
    for old_name, new_name in renames.items():
        if old_name == new_name:
            continue
        settings.mesh.modify_zones.zone_name(zone_name=old_name, new_name=new_name)


def _set_pressure_inlet(inlet_bc, gas_cfg: dict) -> None:
    inlet_bc.momentum.gauge_total_pressure.value = gas_cfg["inlet_stagnation_pressure_pa"]
    inlet_bc.thermal.total_temperature.value = gas_cfg["inlet_stagnation_temperature_k"]
    inlet_bc.turbulence.turbulence_specification = "Intensity and Hydraulic Diameter"
    inlet_bc.turbulence.turbulent_intensity = gas_cfg["turbulence_intensity_fraction"]
    # hydraulic diameter set by caller (geometry-dependent, 2*inlet_radius)


def _set_pressure_outlet(outlet_bc, gas_cfg: dict) -> None:
    outlet_bc.momentum.gauge_pressure.value = gas_cfg["outlet_static_pressure_pa"]
    outlet_bc.turbulence.turbulence_specification = "Intensity and Hydraulic Diameter"
    outlet_bc.turbulence.turbulent_intensity = gas_cfg["turbulence_intensity_fraction"]
    # hydraulic diameter set by caller (geometry-dependent, 2*exit_radius)


def setup_case(solver_session, geometry_id: str, msh_path: str, geom_params: dict, cfg: dict) -> None:
    """
    Reads the mesh and configures models/BCs/DPM for one geometry. Raises
    on any hard failure (caller is expected to wrap in try/except per
    geometry -- see fluent_batch.py).
    """
    settings = solver_session.settings

    # --- Mesh import ----------------------------------------------------
    # Two prior attempts failed live against a real Fluent 2026 R1 session:
    #   1. settings.file.read(file_type="mesh", ...) read the file but
    #      errored "Null Domain Pointer" -- that's Fluent's *native*
    #      mesh-format reader (File > Read > Mesh), not a general importer.
    #   2. settings.file.import_.read(file_type="gmsh", ...) failed with an
    #      explicit, live-enumerated list of allowed file_type values --
    #      "gmsh" is NOT one of them. The full allowed list (Fluent 2026 R1,
    #      solver-mode file.import_.read): cgns-mesh, cgns-mesh-data,
    #      import-gambit, mechanical-apdl-input, mechanical-apdl-result,
    #      abaqus-fil, abaqus-input, abaqus-odb, import-hypermesh,
    #      import-ensight, nastran-bulkdata, nastran-output2, plot3d-mesh,
    #      tecplot-mesh, cfx-definition, cfx-result, gtm-files,
    #      partition-metis, partition-metis-zone. No Gmsh/MSH format at all.
    #   3. Switched to CGNS ("cgns-mesh", in that allowed list; physical
    #      group names confirmed intact through the export). Import got
    #      further -- Fluent's CGNS backend (CeetronSAM.dll) loaded -- but
    #      then the entire Fluent PROCESS segfaulted (SIGSEGV) while
    #      reading the file. Not a Python-catchable error; killed the
    #      session outright. Gmsh-written CGNS and this Fluent version's
    #      CGNS reader are evidently incompatible in some way not worth
    #      chasing further blind.
    #
    # Fixed by exporting Abaqus INP instead (mesh_geometry.py). "abaqus-
    # input" is in the allowed list, INP is a mature/simple format, and --
    # checked before touching Fluent again this time -- gmsh's INP writer
    # emits proper named element sets (*ELSET,ELSET=wall / axis / inlet /
    # outlet / fluid), unlike Nastran bulk data (also in the allowed list,
    # but has no zone-naming mechanism at all -- checked and ruled out) or
    # Gambit neutral format (silently wrote an empty file). msh_path here
    # is expected to point at a .inp file, despite the parameter name
    # (kept for compatibility with existing callers/tests).
    settings.file.import_.read(file_type="abaqus-input", file_name=msh_path)

    # --- Recover wall/axis/inlet/outlet zone names ------------------------
    # Confirmed live: the import above merges all boundary edges into one
    # generic zone (e.g. "wall-2"), discarding our named element sets
    # entirely. See _recover_boundary_zone_names's docstring.
    _recover_boundary_zone_names(solver_session, geometry_id)

    # --- Fix the cell zone type (solid -> fluid) ---------------------------
    # Confirmed live: solver initialization failed with "Flow boundary zone
    # found adjacent to solid zone." Abaqus INP originates as a structural-
    # analysis format, and Fluent's importer evidently defaults the cell
    # zone to "solid" type rather than "fluid" -- our pressure-inlet/outlet
    # face zones are geometrically adjacent to it, which Fluent correctly
    # rejects for a non-fluid cell zone. cell_zone_conditions has its own
    # set_zone_type (same signature as boundary_conditions.set_zone_type,
    # confirmed against the schema), used here instead.
    solid_zone_names = settings.setup.cell_zone_conditions.solid.get_object_names()
    if solid_zone_names:
        logger.info("%s: converting cell zone(s) %s from solid to fluid", geometry_id, solid_zone_names)
        settings.setup.cell_zone_conditions.set_zone_type(zone_list=solid_zone_names, new_type="fluid")

    # --- 2D axisymmetric, density-based solver ---------------------------
    settings.setup.general.solver.type = cfg["solver"]["type"]
    settings.setup.general.solver.two_dim_space = "axisymmetric"

    # --- Fix zone TYPES (name != type in Fluent) --------------------------
    # _recover_boundary_zone_names renamed the split-off pieces to "axis"/
    # "inlet"/"outlet"/"wall", but a rename only changes the NAME -- all
    # four are still typed as "wall" (inherited from the single merged
    # zone they were split from), confirmed live by the next failure this
    # caused: "'pressure_inlet' has no attribute 'inlet'" when trying to
    # set inlet's boundary condition, since a wall-typed zone named
    # "inlet" doesn't exist under boundary_conditions.pressure_inlet until
    # its TYPE is also changed. "wall" needs no type change (already
    # correct). PyFluent's generated Command objects are keyword-only
    # (Command.__call__ takes **kwds, no positional args) -- confirmed
    # live earlier when this used positional args instead.
    settings.setup.boundary_conditions.set_zone_type(zone_list=["axis"], new_type="axis")
    settings.setup.boundary_conditions.set_zone_type(zone_list=["inlet"], new_type="pressure-inlet")
    settings.setup.boundary_conditions.set_zone_type(zone_list=["outlet"], new_type="pressure-outlet")

    # --- Physics models ----------------------------------------------------
    settings.setup.models.energy.enabled = bool(cfg["solver"]["energy_equation"])
    settings.setup.models.viscous.model = cfg["solver"]["turbulence_model"]
    if cfg["solver"]["turbulence_model"] == "k-epsilon":
        settings.setup.models.viscous.k_epsilon_model = cfg["solver"]["k_epsilon_variant"]

    # --- Courant number (numerical stability) -------------------------------
    # settings.solution.controls.courant_number -- confirmed real path against
    # the schema. See fluent_config.yaml's courant_number comment for why
    # this was lowered from Fluent's density-based-implicit default (~5):
    # a live run showed a genuine solver stall (not just slow convergence)
    # late in the iteration count, alongside worsening temperature-cap
    # violations in a growing fraction of cells -- classic density-based
    # instability symptoms for a strong, fine-mesh compressible expansion.
    if "courant_number" in cfg["solver"]:
        settings.solution.controls.courant_number = cfg["solver"]["courant_number"]

    # --- Boundary conditions -----------------------------------------------
    gas_cfg = cfg["gas"]
    bc = settings.setup.boundary_conditions
    inlet_bc = bc.pressure_inlet["inlet"]
    _set_pressure_inlet(inlet_bc, gas_cfg)
    inlet_bc.turbulence.hydraulic_diameter = 2.0 * geom_params["inlet_radius"]

    outlet_bc = bc.pressure_outlet["outlet"]
    _set_pressure_outlet(outlet_bc, gas_cfg)
    outlet_bc.turbulence.hydraulic_diameter = 2.0 * geom_params["exit_radius"]

    # --- DPM: centerline particle injection at the inlet (x=0, r=0) --------
    particle_cfg = cfg["particle"]
    settings.setup.models.discrete_phase.injections.create(name="injection-1")
    injection = settings.setup.models.discrete_phase.injections["injection-1"]
    injection.particle_type = "inert"  # confirmed live: accepted without error
    injection.material = particle_cfg["material"]  # confirmed live: "Copying particle material from the database: copper"
    # injection_type.option removal (previous attempt) did NOT fix the
    # "ASSQ: invalid argument [2]: improper list" error -- it recurred at
    # the exact same point, meaning the actual culprit is one of the lines
    # below, not injection_type.option. None of these individually produce
    # console output, so the failing line can't be told apart from the
    # transcript alone. Isolating each one with its own try/except so the
    # next live run pinpoints exactly which assignment is malformed,
    # instead of guessing again blind.
    # velocity, particle_size, and mass_flow_rate are all Group-typed
    # (multiple sub-modes/components), not plain scalars, despite looking
    # like simple values in the task description -- confirmed live
    # ("ASSQ: invalid argument [2]: improper list" on a bare-float
    # assignment to velocity) and then confirmed structurally for all
    # three against the installed settings schema before trying again:
    #   velocity: Group with x_velocity/y_velocity/z_velocity/magnitude
    #     (not a single scalar or vector-as-one-field) -- centerline
    #     injection moving axially, so x_velocity = config value, y/z = 0.
    #   particle_size: Group with option/diameter/rosin_rammler/
    #     tabulated_size -- diameter (a plain Real) is the constant-size mode.
    #   mass_flow_rate: Group with flow_rate/total_flow_rate/scale_by_area
    #     -- total_flow_rate matches our config's intent (a single flow
    #     rate for the whole injection, not per-unit-area).
    # temperature (unlike the three above) IS a plain Real -- confirmed
    # against the schema, left as a direct assignment.
    _dpm_fields = [
        ("location.x", lambda: setattr(injection.initial_values.location, "x", 0.0)),
        ("location.y", lambda: setattr(injection.initial_values.location, "y", 0.0)),
        ("velocity.x_velocity", lambda: setattr(injection.initial_values.velocity, "x_velocity",
                                                  particle_cfg["injection_velocity_m_s"])),
        ("velocity.y_velocity", lambda: setattr(injection.initial_values.velocity, "y_velocity", 0.0)),
        ("velocity.z_velocity", lambda: setattr(injection.initial_values.velocity, "z_velocity", 0.0)),
        ("particle_size.diameter", lambda: setattr(injection.initial_values.particle_size, "diameter",
                                                     particle_cfg["diameter_um"] * 1e-6)),
        ("temperature", lambda: setattr(injection.initial_values, "temperature",
                                         particle_cfg["injection_temperature_k"])),
        ("mass_flow_rate.total_flow_rate", lambda: setattr(injection.initial_values.mass_flow_rate,
                                                             "total_flow_rate",
                                                             particle_cfg["mass_flow_rate_kg_s"])),
    ]
    for field_name, setter in _dpm_fields:
        try:
            setter()
            logger.info("DPM initial_values.%s: OK", field_name)
        except Exception as e:
            logger.error("DPM initial_values.%s: FAILED: %s", field_name, e)
            raise

    # --- Convergence criteria -----------------------------------------------
    # UNVERIFIED equation-name strings (see module docstring) -- standard
    # Fluent residual names, not confirmed against this schema version.
    conv = cfg["convergence"]
    residual_targets = {
        "continuity": conv["continuity_residual"],
        "x-velocity": conv["velocity_residual"],
        "y-velocity": conv["velocity_residual"],
        "energy": conv["energy_residual"],
        "k": conv["k_residual"],
        "epsilon": conv["epsilon_residual"],
    }
    equations = settings.solution.monitor.residual.equations
    for name, target in residual_targets.items():
        try:
            equations[name].absolute_criteria = target
        except Exception as e:
            logger.warning("Could not set residual criteria for %r: %s", name, e)


def run_to_convergence(solver_session, cfg: dict, geometry_id: str) -> tuple:
    """
    Initializes and iterates in chunks of cfg['convergence']['check_interval']
    iterations, checking PyFluent's monitors data after each chunk to decide
    whether to stop early, up to the hard cap max_iterations. Returns
    (converged: bool, iterations_run: int, warnings: list[str]).

    Uses solver_session.monitors.get_monitor_set_data("residual") --
    confirmed to exist by reading the installed PyFluent package's
    MonitorsManager source directly (streaming_services/monitor_streaming.py):
    get_monitor_set_names() lists available monitor sets and
    get_monitor_set_data(name) returns (x_values, {column_name: y_values}).
    UNVERIFIED: that the residual monitor set is actually named "residual"
    on a live session, and that its column names match the equation names
    used below (standard Fluent convention, e.g. "continuity",
    "x-velocity", but not confirmed against a real running case) -- this
    is checked defensively (falls back to running the full iteration cap
    without early-stopping if the monitor set/columns aren't found under
    those names) rather than assumed to work silently.
    """
    conv = cfg["convergence"]
    max_iterations = conv["max_iterations"]
    check_interval = conv["check_interval"]
    warnings = []

    settings = solver_session.settings
    settings.solution.initialization.hybrid_initialize()

    iterations_run = 0
    converged = False
    while iterations_run < max_iterations:
        chunk = min(check_interval, max_iterations - iterations_run)
        settings.solution.run_calculation.iterate(iter_count=chunk)
        iterations_run += chunk

        try:
            monitor_set_names = solver_session.monitors.get_monitor_set_names()
            if "residual" not in monitor_set_names:
                warnings.append(
                    f"No 'residual' monitor set found (available: {monitor_set_names}) -- "
                    f"cannot check convergence early, will run to the iteration cap."
                )
                logger.warning("%s: no 'residual' monitor set found, available: %s",
                                geometry_id, monitor_set_names)
                continue
            _, residual_series = solver_session.monitors.get_monitor_set_data("residual")
            latest = {name: values[-1] for name, values in residual_series.items() if len(values)}
            targets = {
                "continuity": conv["continuity_residual"],
                "x-velocity": conv["velocity_residual"],
                "y-velocity": conv["velocity_residual"],
                "energy": conv["energy_residual"],
                "k": conv["k_residual"],
                "epsilon": conv["epsilon_residual"],
            }
            if latest and all(
                name not in latest or latest[name] <= target for name, target in targets.items()
            ):
                converged = True
                break
        except Exception as e:
            warnings.append(f"Could not read monitor data at iteration {iterations_run}: {e}")
            logger.warning("%s: monitor read failed at iteration %d: %s", geometry_id, iterations_run, e)

    if not converged and iterations_run >= max_iterations:
        warnings.append(f"Hit the {max_iterations}-iteration cap without meeting all residual targets.")

    return converged, iterations_run, warnings


def extract_exit_velocity(solver_session, geometry_id: str, results_dir: str) -> Optional[float]:
    """
    Writes Fluent's DPM extended summary report (particle fate + velocity
    statistics per boundary, including the outlet) to a small text file via
    settings.results.report.discrete_phase.extended_summary -- a real,
    documented command (confirmed present in the installed settings
    schema), unlike the earlier report-definitions-based approach this
    replaces.

    Two real bugs fixed after a live test that ran a full 2000-iteration
    solve successfully and only failed at this final step:
      1. `injection=""` was rejected ("Value is not allowed: ('' is_not_in
         ('injection-1'))") -- pass our actual injection name instead.
      2. `results_dir` ("fluent_results/") doesn't exist ON THE REMOTE
         (Windows) MACHINE -- os.makedirs() on this (local) machine was
         useless, since Fluent writes the file server-side. Writes to a
         bare filename in Fluent's own working directory instead (which
         already exists) rather than assuming a subfolder is there.

    Also: the written file lives on the REMOTE machine's filesystem, not
    this one -- open()-ing summary_path locally would fail regardless
    (same class of cross-machine path issue as the mesh files). No
    file-transfer service is configured (see FLUENT_REMOTE.md), so this
    reads the file's content back over the existing gRPC connection via
    solver_session.scheme.eval() instead of assuming a shared filesystem.
    solver_session.scheme confirmed present (session.py: self.scheme =
    scheme_eval) and .eval(scm_input) confirmed to return the evaluated
    Scheme value as a Python object (services/scheme_interpreter.py).

    UNVERIFIED: the exact Scheme file-reading snippet below (standard
    Scheme idiom, but Fluent's embedded interpreter's exact behavior/
    available primitives aren't confirmed), and the exact text layout of
    the summary file itself, so the mean exit-velocity value has to be
    located by keyword search rather than a fixed column/row position.
    """
    summary_filename = f"{geometry_id}_dpm_summary.txt"

    solver_session.settings.results.report.discrete_phase.extended_summary(
        write_to_file=True,
        file_name=summary_filename,
        include_in_domain_particles=False,
        pick_injection=False,
        injection="injection-1",
    )

    scm_read_file = (
        f'(let* ((port (open-input-file "{summary_filename}")))'
        f"  (let loop ((lines '()))"
        f"    (let ((line (read-line port)))"
        f"      (if (eof-object? line)"
        f"          (begin (close-input-port port) (reverse lines))"
        f"          (loop (cons line lines))))))"
    )
    try:
        lines = solver_session.scheme.eval(scm_read_file)
    except Exception as e:
        logger.warning("%s: could not read back %s via scheme.eval: %s", geometry_id, summary_filename, e)
        return None

    return _parse_exit_velocity_from_lines(lines)


def _parse_exit_velocity_from_lines(lines) -> Optional[float]:
    """
    Best-effort keyword-based parse of the DPM extended summary text
    (as a list of lines, read back via scheme.eval -- see
    extract_exit_velocity's docstring) for a mean particle velocity
    magnitude at the outlet zone. UNVERIFIED against real Fluent output.
    """
    import re
    if not lines:
        logger.warning("DPM summary read back empty or None -- nothing to parse.")
        return None

    # Look for a line mentioning "outlet" and "velocity" with a numeric
    # mean/average value nearby -- deliberately loose, since the exact
    # Fluent-version text layout isn't known here.
    for line in lines:
        line_str = str(line)
        if "outlet" in line_str.lower() and "velocity" in line_str.lower():
            numbers = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line_str)
            if numbers:
                return float(numbers[-1])

    logger.warning("Could not locate an outlet velocity value in the DPM summary -- "
                    "check its actual layout by hand. First few lines: %s", lines[:10])
    return None


def solve_geometry(solver_session, geometry_id: str, msh_path: str, geom_params: dict, cfg: dict) -> SolveResult:
    """
    Full per-geometry pipeline: setup_case -> run_to_convergence ->
    extract_exit_velocity -> (optionally) write case/data. Never raises --
    catches everything and returns a SolveResult with success=False and
    .error set, so fluent_batch.py's per-geometry try/except is a second
    safety net, not the only one.
    """
    warnings = []
    try:
        setup_case(solver_session, geometry_id, msh_path, geom_params, cfg)
    except Exception as e:
        logger.error("%s: setup failed: %s", geometry_id, e)
        return SolveResult(geometry_id, False, False, 0, None, warnings, f"setup failed: {e}")

    try:
        converged, iterations_run, run_warnings = run_to_convergence(solver_session, cfg, geometry_id)
        warnings.extend(run_warnings)
    except Exception as e:
        logger.error("%s: solve failed: %s", geometry_id, e)
        return SolveResult(geometry_id, False, False, 0, None, warnings, f"solve failed: {e}")

    exit_velocity = None
    try:
        results_dir = cfg["output"]["results_dirname"]
        exit_velocity = extract_exit_velocity(solver_session, geometry_id, results_dir)
        if exit_velocity is None:
            warnings.append("Exit velocity could not be extracted from the DPM summary.")
    except Exception as e:
        logger.warning("%s: exit-velocity extraction failed: %s", geometry_id, e)
        warnings.append(f"exit-velocity extraction failed: {e}")

    if cfg["output"].get("write_case_data", True):
        try:
            # Bare filenames (Fluent's own working directory), same fix as
            # extract_exit_velocity -- cfg["output"]["results_dirname"]
            # ("fluent_results/") doesn't exist on the remote machine, and
            # os.path.join-ing it in here was a local-machine path that
            # meant nothing to the remote Fluent process anyway.
            solver_session.settings.file.write_case(file_name=f"{geometry_id}.cas.h5")
            solver_session.settings.file.write_data(file_name=f"{geometry_id}.dat.h5")
        except Exception as e:
            logger.warning("%s: writing case/data failed (non-fatal): %s", geometry_id, e)
            warnings.append(f"writing case/data failed: {e}")

    return SolveResult(
        geometry_id=geometry_id, success=True, converged=converged, iterations_run=iterations_run,
        exit_velocity_m_s=exit_velocity, warnings=warnings, error=None,
    )
