# Remote Fluent batch solve via PyFluent (gRPC)

Drives the Fluent DPM cold-spray batch solve for all 72 nozzle meshes
**from this machine**, over a network connection to Fluent running on a
separate Windows PC. There is no file-transfer step and no hand-written
`.jou` script -- the solve loop, boundary conditions, and result
extraction are all real Python code here (`fluent_solve.py`), talking to
the remote Fluent session through PyFluent's `connect_to_fluent` gRPC
client. The Windows PC's only job is to have Fluent running with its gRPC
server exposed and reachable on the network.

## What has and hasn't been verified

**Verified** (checked directly against the `ansys-fluent-core` package
installed in this environment -- real source/type-stub inspection, not
documentation-from-memory):
- `pyfluent.connect_to_fluent`'s exact signature (`ip`, `port`, `address`,
  `server_info_file_name`, `password`, `allow_remote_host`, etc.).
- Every PyFluent settings path used in `fluent_solve.py` (mesh reading,
  2D axisymmetric mode, solver type, the `axis` zone-type fix, viscous/
  energy models, pressure-inlet/outlet boundary condition structure, DPM
  injection structure, residual monitor structure, run-calculation/
  initialization commands, case/data writing) -- confirmed against the
  generated settings schema shipped in the package
  (`ansys/fluent/core/generated/solver/settings_261.pyi`).
- The DPM result-extraction command,
  `settings.results.report.discrete_phase.extended_summary(...)` --
  confirmed to exist with that exact signature. This replaced an earlier,
  much shakier plan (a `report-definitions`-based "dpm-average" report
  type) once static inspection showed no such report type actually exists
  in the schema.
- The monitors API used for early-stopping,
  `session.monitors.get_monitor_set_names()` /
  `get_monitor_set_data(name)` -- confirmed against
  `streaming_services/monitor_streaming.py`'s actual `MonitorsManager`
  class (an initial guess at a method called `get_monitor_set_data_all`
  was wrong and was corrected after checking the real source).

**NOT verified** (no live Fluent session was available while building
this -- these need first-contact validation on your Windows machine,
starting with `run_subset.py` on 2-3 geometries, before trusting the full
batch):
- **Mesh import**: `settings.file.read(file_type="mesh", file_name=...)`
  is used to read each `geometry_XXXX.msh` (a Gmsh MSH2 file -- see
  `MESHING.md`). No dedicated "gmsh" import command exists anywhere in
  the generated settings schema, which is itself a signal worth taking
  seriously. **This is the single most likely thing to fail first** --
  if it does, the fallback is importing via Fluent Meshing mode and
  re-exporting a native Fluent case before the solver reads it (a real
  extra step, not yet implemented here, in case it's needed).
- Exact allowed-value strings for `injection_type.option` (guessed:
  `"single"`), `particle_type` (guessed: `"inert"`), and whether
  `"copper"` matches an exact material name in Fluent's material
  database on the target machine. These are runtime-populated
  `AllowedValuesMixin` classes with no values baked into the static
  stub -- genuinely can't be checked without a live session.
- Residual equation names (`"continuity"`, `"x-velocity"`, `"y-velocity"`,
  `"energy"`, `"k"`, `"epsilon"`) used both for setting convergence
  criteria and for reading back monitor data -- standard Fluent naming,
  not confirmed against this schema version specifically.
- Whether the residual monitor set is actually named `"residual"`.
- The exact text layout of the DPM extended-summary file, so
  `fluent_solve._parse_exit_velocity_from_summary` uses a loose
  keyword-based search (a line containing both "outlet" and "velocity")
  rather than a fixed row/column position. Confirm this actually finds
  the right number against a real summary file from the subset run.
- Whether one long-lived session reading/solving/resetting across all 72
  geometries in sequence is reliable, or whether state leaks between runs
  (residual monitors, DPM injections, and report definitions all persist
  on a session unless explicitly cleared). `run_subset.py` exists
  specifically to observe this before `fluent_batch.py` assumes an
  answer -- watch for anything from geometry 2 onward behaving
  differently than geometry 1 in the subset run.

Every one of these is also flagged inline in `fluent_solve.py`'s
docstrings, next to the code it applies to.

## 1. Start Fluent's gRPC server on the Windows PC

Three ways to do this -- pick whichever fits how Fluent is normally
launched there:

**A. From the Fluent Launcher GUI** (simplest if starting fresh): when
launching Fluent, there's typically a "Start python server" / gRPC server
option in the launcher's advanced settings, or you can start Fluent
normally and then run the TUI command below once it's open.

**B. Via the TUI command**, from within an already-running Fluent
session (Console/TUI panel):
```
server/start-server
```
This prints connection info (port and password) directly in the console.

**C. Via the command line**, launching Fluent with the server started
automatically and its connection info written to a file:
```
fluent 2ddp -sifile=server_info.txt
```
(`-sifile=` is the flag that tells Fluent to write a server-info file;
`2ddp` matches the 2D double-precision solver this project's meshes need
-- see `MESHING.md`.)

Whichever method you use, Fluent will report (or write to the server-info
file) three things: an **IP/port** and a **password**. You need at least
one of:
- the **server-info file** itself (copy it to this machine), or
- the **IP, port, and password** individually.

**TLS vs. insecure mode**: the Fluent Launcher's General Options tab has
an "Allow Remote gRPC Host" checkbox (must be checked for any of this to
work) alongside a "gRPC Certificates Folder" field and a "gRPC Insecure
Mode" checkbox. Fluent's remote gRPC defaults to TLS, which means
matching certificates would need to be generated and copied to this
machine too. For a private, trusted two-machine setup, the simpler option
is: leave the certificates folder blank and check **"gRPC Insecure Mode"**
-- then pass `--insecure` to every script here (`test_connection.py`,
`run_subset.py`, `fluent_batch.py`) to match on this side.

## 2. Confirm connectivity with `test_connection.py`

From this machine, with the server-info file copied over (or the IP/port/
password on hand):

```bash
pip install -r requirements.txt   # installs ansys-fluent-core if not already present

python3 test_connection.py --server-info-file path/to/server_info.txt
# or
python3 test_connection.py --ip 192.168.1.50 --port 12345 --password abc123
```

This does nothing except connect, print the reported Fluent version, and
disconnect. **Do not proceed past this step until it succeeds.** If it
fails:
- Check the Windows PC's firewall allows inbound connections on the gRPC
  port Fluent printed.
- Confirm both machines are on the same network / can reach each other
  (a simple `ping <windows-ip>` from this machine is a reasonable first
  check).
- Server-info files and passwords are per-launch -- if Fluent was
  restarted on the Windows side, get a fresh one.
- `allow_remote_host=True` is passed automatically by every script here
  (PyFluent defaults to *not* allowing non-localhost connections
  otherwise) -- if you hit a "remote host not allowed"-style error even
  though these scripts should already be setting that, that's worth
  reporting back, since it would mean something about this project's
  assumption here was wrong.

Once this connects cleanly, note down whichever connection details worked
-- `run_subset.py` and `fluent_batch.py` take the identical arguments.

## 3. Validate on a small subset

```bash
python3 run_subset.py --server-info-file path/to/server_info.txt --n 2
```

This solves 2 geometries (raise `--n` to 3 if you want) end-to-end --
mesh import, boundary conditions, DPM injection, iterate-to-convergence,
exit-velocity extraction -- and prints a summary. Go through the checklist
it prints at the end:

1. Did the mesh import succeed for every geometry? (Watch for `setup
   failed` errors -- this is where the unverified `gmsh` import command
   would show up as broken, if it is.)
2. Did at least one case converge within the iteration cap
   (`max_iterations` in `fluent_config.yaml`, default 2000)?
3. Is the exit velocity a physically plausible number -- not `None`, not
   `0`, not absurdly large? A `None` here means
   `_parse_exit_velocity_from_summary`'s keyword search didn't find what
   it expected in the real summary file layout; open the summary file by
   hand (written to `fluent_config.yaml`'s `output.results_dirname`,
   `<geometry_id>_dpm_summary.txt`) and adjust the parser if needed.
4. Do any warnings mention unrecognized residual/monitor names? If so,
   `fluent_solve.py`'s `run_to_convergence`/`setup_case` need their
   equation-name strings adjusted to match what your Fluent version
   actually calls them.

**Only move on to the full batch once this looks right.** If something's
broken, fix it in `fluent_solve.py` and re-run `run_subset.py` -- iterate
here, not on the full 72-geometry run.

## 4. Run the full batch

```bash
python3 fluent_batch.py --server-info-file path/to/server_info.txt
```

Connects once, solves every `geometry_XXXX.msh` in `mesh_output_v1/` in
sequence, and writes results **incrementally** to
`fluent_batch_results.csv` (one row per geometry, appended as each one
finishes) so a crash partway through doesn't lose everything already
computed. Each geometry is wrapped in its own `try/except`, so one failure
(a bad mesh, a solver crash, a lost connection mid-geometry) is logged and
the batch continues rather than stopping. Progress is logged as
`[N/72] geometry_XXXX: ...` with running converged/failed counts, since a
full batch may take a long time (72 geometries x up to 2000 iterations
each) -- there's no built-in time estimate here since solve time per
geometry will vary and hasn't been observed yet.

If the connection drops partway through, `fluent_batch_results.csv` still
has every row completed up to that point -- you can inspect what's there,
fix whatever caused the drop, and (for now) re-run from scratch, since
there's no resume-from-checkpoint logic here yet.

## 5. Build the surrogate training set

Once `fluent_batch_results.csv` exists (fully or partially):

```bash
python3 build_training_set.py --fluent-results fluent_batch_results.csv \
    --sweep-dir sweep_output_v1 --out-csv training_data.csv
```

Joins each result row back to its geometry's parameters
(`sweep_output_v1/geometry_XXXX.json`) and writes one master CSV: geometry
parameters + `exit_velocity_m_s` + `convergence_status`
(`converged`/`ran_not_converged`/`failed`) + `iterations_run` + `warnings`
+ `error`. Prints a convergence-rate breakdown in the same style used
elsewhere in this project (the geometry-rejection and mesh-quality
reports). Non-converged and failed runs are kept in the CSV, not silently
dropped -- filter on `convergence_status == "converged"` before training
the GPR surrogate, or handle the others deliberately.

## Config: `fluent_config.yaml`

Process conditions, particle properties, solver settings, and convergence
targets -- same status as `nozzle_geometry.py`'s `max_curvature`/
`max_fillet_fraction` and `mesh_geometry.py`'s `n_axial`/
`first_cell_height`: every numeric value is a configurable parameter with
a logged placeholder default, not a validated engineering value. The
paper's exact published process conditions (Badali et al.) weren't
available while building this pipeline, so these are typical nitrogen
cold-spray literature ballparks standing in for the real numbers --
individually flagged in the YAML file's comments. Confirm against the
actual paper / experimental setup before trusting results beyond a
pipeline smoke test.

## Network / firewall notes

Nothing beyond `allow_remote_host=True` (already set in every script
here) was found to be necessary from PyFluent's side while building this
-- but that's a documentation-level claim, not something confirmed
against a real cross-machine connection, since no live Fluent instance
was available. The realistic failure modes to check first if
`test_connection.py` doesn't connect are ordinary network ones (Windows
Firewall blocking the gRPC port, both machines not actually being on the
same reachable network/VPN, a stale port/password from an old Fluent
launch) rather than anything PyFluent-specific.
