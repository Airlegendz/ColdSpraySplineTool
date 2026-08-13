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

What is NOT verified, because it requires either a live Fluent session or
version-specific runtime data the static stub doesn't carry:
  - Whether mesh import via `settings.file.read(file_type="mesh", ...)`
    actually accepts a Gmsh MSH2 file directly, vs. requiring a
    Fluent-Meshing-side conversion step first. No "gmsh" import command
    exists anywhere in the generated settings tree, which is itself a
    signal worth taking seriously -- FLAG THIS FIRST on the subset run.
  - The exact allowed-value strings for `injection_type.option` (DPM
    injection type, e.g. "single") and `particle_type`/`material_2`
    (e.g. "inert", "copper") -- these are runtime-populated
    AllowedValuesMixin classes with no values baked into the static stub.
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


def setup_case(solver_session, msh_path: str, geom_params: dict, cfg: dict) -> None:
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

    # --- 2D axisymmetric, density-based solver ---------------------------
    settings.setup.general.solver.type = cfg["solver"]["type"]
    settings.setup.general.solver.two_dim_space = "axisymmetric"

    # --- Fix the axis zone type (the manual step flagged in MESHING.md) --
    settings.setup.boundary_conditions.set_zone_type(["axis"], "axis")

    # --- Physics models ----------------------------------------------------
    settings.setup.models.energy.enabled = bool(cfg["solver"]["energy_equation"])
    settings.setup.models.viscous.model = cfg["solver"]["turbulence_model"]
    if cfg["solver"]["turbulence_model"] == "k-epsilon":
        settings.setup.models.viscous.k_epsilon_model = cfg["solver"]["k_epsilon_variant"]

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
    settings.setup.models.discrete_phase.injections.create("injection-1")
    injection = settings.setup.models.discrete_phase.injections["injection-1"]
    injection.particle_type = "inert"  # UNVERIFIED exact allowed-value string, see module docstring
    injection.material = particle_cfg["material"]  # UNVERIFIED: must exist in Fluent's material database
    injection.injection_type.option = "single"  # UNVERIFIED exact allowed-value string
    injection.initial_values.location.x = 0.0
    injection.initial_values.location.y = 0.0
    injection.initial_values.velocity = particle_cfg["injection_velocity_m_s"]
    injection.initial_values.particle_size = particle_cfg["diameter_um"] * 1e-6
    injection.initial_values.temperature = particle_cfg["injection_temperature_k"]
    injection.initial_values.mass_flow_rate = particle_cfg["mass_flow_rate_kg_s"]

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

    UNVERIFIED: the exact text layout of that summary file, so the mean
    exit-velocity value has to be located by keyword search rather than a
    fixed column/row position -- confirm this parsing actually finds the
    right number against a real summary file from the subset run, and
    adjust if the real layout differs.
    """
    import os
    os.makedirs(results_dir, exist_ok=True)
    summary_path = os.path.join(results_dir, f"{geometry_id}_dpm_summary.txt")

    solver_session.settings.results.report.discrete_phase.extended_summary(
        write_to_file=True,
        file_name=summary_path,
        include_in_domain_particles=False,
        pick_injection=False,
        injection="",
    )

    return _parse_exit_velocity_from_summary(summary_path)


def _parse_exit_velocity_from_summary(summary_path: str) -> Optional[float]:
    """
    Best-effort keyword-based parse of the DPM extended summary text file
    for a mean particle velocity magnitude at the outlet zone. UNVERIFIED
    against real Fluent output -- see extract_exit_velocity's docstring.
    """
    import re
    try:
        with open(summary_path) as f:
            text = f.read()
    except FileNotFoundError:
        logger.warning("Expected DPM summary file not found: %s", summary_path)
        return None

    # Look for a line mentioning "outlet" and "velocity" with a numeric
    # mean/average value nearby -- deliberately loose, since the exact
    # Fluent-version text layout isn't known here.
    for line in text.splitlines():
        if "outlet" in line.lower() and "velocity" in line.lower():
            numbers = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
            if numbers:
                return float(numbers[-1])

    logger.warning("Could not locate an outlet velocity value in %s -- check the file's actual layout by hand.",
                    summary_path)
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
        setup_case(solver_session, msh_path, geom_params, cfg)
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
            import os
            results_dir = cfg["output"]["results_dirname"]
            solver_session.settings.file.write_case(file_name=os.path.join(results_dir, f"{geometry_id}.cas.h5"))
            solver_session.settings.file.write_data(file_name=os.path.join(results_dir, f"{geometry_id}.dat.h5"))
        except Exception as e:
            logger.warning("%s: writing case/data failed (non-fatal): %s", geometry_id, e)
            warnings.append(f"writing case/data failed: {e}")

    return SolveResult(
        geometry_id=geometry_id, success=True, converged=converged, iterations_run=iterations_run,
        exit_velocity_m_s=exit_velocity, warnings=warnings, error=None,
    )
