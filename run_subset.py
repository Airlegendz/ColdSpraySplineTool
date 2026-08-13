"""
Validation step: connects to the remote Fluent session ONCE and solves
just a handful of geometries (2-3 by default), reporting results for
manual review before committing to the full 72-geometry batch via
fluent_batch.py.

Run this after test_connection.py has confirmed basic connectivity. This
step is where the real boundary-condition/DPM/mesh-import setup gets its
first live test -- everything in fluent_solve.py up to this point has only
been checked against PyFluent's static settings schema, not run.

IMPORTANT: Fluent's file.read() command needs a path Fluent itself (on the
Windows PC) can see -- it does NOT tunnel file contents through the gRPC
solver connection the way session.upload() would (and upload() is a no-op
unless a separate file-transfer service/server is configured, which this
project does not set up). So the .msh files this script points Fluent at
must already exist somewhere on the WINDOWS machine's filesystem --
--mesh-dir (below) is just used locally to enumerate which geometries
exist and load their .json parameters; --remote-mesh-dir is the Windows-
side path where the matching .msh files need to have been copied to
first (network share, USB, cloud sync, whatever's convenient).

Usage:
    python3 run_subset.py --server-info-file path/to/server_info.txt \\
        --remote-mesh-dir "C:\\fluent_meshes" --n 2
    python3 run_subset.py --ip 192.168.1.50 --port 12345 --password abc123 \\
        --remote-mesh-dir "C:\\fluent_meshes" --n 3
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os

import yaml

from fluent_solve import solve_geometry

logger = logging.getLogger("run_subset")


def connect(args):
    import ansys.fluent.core as pyfluent

    if args.server_info_file:
        return pyfluent.connect_to_fluent(
            server_info_file_name=args.server_info_file,
            password=args.password,
            allow_remote_host=True,
            insecure_mode=args.insecure,
        )
    return pyfluent.connect_to_fluent(
        ip=args.ip, port=args.port, password=args.password, allow_remote_host=True,
        insecure_mode=args.insecure,
    )


def main():
    parser = argparse.ArgumentParser(description="Solve a small subset of geometries to validate the setup.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--server-info-file")
    group.add_argument("--ip")
    parser.add_argument("--port", type=int)
    parser.add_argument("--password")
    parser.add_argument("--insecure", action="store_true",
                         help="Match Fluent Launcher's 'gRPC Insecure Mode' checkbox.")
    parser.add_argument("--mesh-dir", default="mesh_output_v1",
                         help="Local directory used to enumerate geometries and load their .json params.")
    parser.add_argument("--remote-mesh-dir", required=True,
                         help="Windows-side directory where the matching .msh files have already been copied "
                              "to (e.g. 'C:\\\\fluent_meshes'). Fluent reads from here, not from --mesh-dir.")
    parser.add_argument("--config", default="fluent_config.yaml")
    parser.add_argument("--n", type=int, default=2, help="Number of geometries to solve (2-3 recommended).")
    args = parser.parse_args()

    if args.ip and args.port is None:
        parser.error("--port is required when using --ip")

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    msh_paths = sorted(glob.glob(os.path.join(args.mesh_dir, "geometry_*.msh")))[:args.n]
    if not msh_paths:
        logger.error("No geometry_*.msh files found in %s", args.mesh_dir)
        return

    logger.info("Connecting to Fluent...")
    session = connect(args)
    logger.info("Connected. Fluent version: %s", session.get_fluent_version())

    results = []
    try:
        for msh_path in msh_paths:
            geometry_id = os.path.splitext(os.path.basename(msh_path))[0]
            json_path = os.path.join(args.mesh_dir, f"{geometry_id}.json")
            with open(json_path) as f:
                geom_params = json.load(f)

            remote_msh_path = f"{args.remote_mesh_dir.rstrip(chr(92)).rstrip('/')}\\{geometry_id}.msh"

            logger.info("=" * 60)
            logger.info("Solving %s (remote path: %s) ...", geometry_id, remote_msh_path)
            result = solve_geometry(session, geometry_id, remote_msh_path, geom_params, cfg)
            results.append(result)

            logger.info(
                "%s: success=%s converged=%s iterations=%d exit_velocity=%s warnings=%s error=%s",
                geometry_id, result.success, result.converged, result.iterations_run,
                result.exit_velocity_m_s, result.warnings, result.error,
            )
    finally:
        session.exit()
        logger.info("Session closed.")

    print()
    print("=" * 60)
    print("SUBSET VALIDATION SUMMARY")
    print("=" * 60)
    for r in results:
        status = "OK" if r.success and r.converged else ("RAN BUT DID NOT CONVERGE" if r.success else "FAILED")
        print(f"  {r.geometry_id}: {status}  exit_velocity={r.exit_velocity_m_s}  "
              f"iterations={r.iterations_run}  error={r.error}")
    print()
    print("Before running the full batch, check:")
    print("  1. Did the mesh import succeed for every geometry? (see 'setup failed' errors above)")
    print("  2. Did at least one case converge within the iteration cap?")
    print("  3. Is the exit velocity a physically plausible number (not None, not 0, not absurdly large)?")
    print("  4. Do the warnings list anything about unrecognized residual/monitor names -- if so, "
          "fluent_solve.py's convergence-checking needs its equation names adjusted.")
    print("See FLUENT_REMOTE.md's validation checklist for more detail.")


if __name__ == "__main__":
    main()
