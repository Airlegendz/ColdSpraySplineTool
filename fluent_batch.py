"""
Full 72-geometry Fluent batch driver. Connects to the remote Fluent
session ONCE, loops over every geometry_XXXX.msh in the mesh output
directory, and solves each one via fluent_solve.solve_geometry -- with a
try/except around each geometry so one failure doesn't kill the whole
batch, and incremental CSV writes so a crash partway through doesn't lose
already-computed results.

Do not run this before run_subset.py has validated the setup on a small
number of geometries -- see FLUENT_REMOTE.md.

Usage:
    python3 fluent_batch.py --server-info-file path/to/server_info.txt
    python3 fluent_batch.py --ip 192.168.1.50 --port 12345 --password abc123
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import logging
import os
import time

import yaml

from fluent_solve import solve_geometry

logger = logging.getLogger("fluent_batch")

CSV_FIELDNAMES = [
    "geometry_id", "success", "converged", "iterations_run",
    "exit_velocity_m_s", "warnings", "error",
]


def connect(args):
    import ansys.fluent.core as pyfluent

    if args.server_info_file:
        return pyfluent.connect_to_fluent(
            server_info_file_name=args.server_info_file,
            password=args.password,
            allow_remote_host=True,
        )
    return pyfluent.connect_to_fluent(
        ip=args.ip, port=args.port, password=args.password, allow_remote_host=True,
    )


def write_row(csv_path: str, result, write_header: bool) -> None:
    mode = "w" if write_header else "a"
    with open(csv_path, mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow({
            "geometry_id": result.geometry_id,
            "success": result.success,
            "converged": result.converged,
            "iterations_run": result.iterations_run,
            "exit_velocity_m_s": result.exit_velocity_m_s,
            "warnings": " | ".join(result.warnings),
            "error": result.error or "",
        })


def main():
    parser = argparse.ArgumentParser(description="Run the full Fluent batch over all geometries.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--server-info-file")
    group.add_argument("--ip")
    parser.add_argument("--port", type=int)
    parser.add_argument("--password")
    parser.add_argument("--mesh-dir", default="mesh_output_v1")
    parser.add_argument("--config", default="fluent_config.yaml")
    parser.add_argument("--out-csv", default="fluent_batch_results.csv")
    args = parser.parse_args()

    if args.ip and args.port is None:
        parser.error("--port is required when using --ip")

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    msh_paths = sorted(glob.glob(os.path.join(args.mesh_dir, "geometry_*.msh")))
    if not msh_paths:
        logger.error("No geometry_*.msh files found in %s", args.mesh_dir)
        return

    total = len(msh_paths)
    logger.info("Connecting to Fluent...")
    session = connect(args)
    logger.info("Connected. Fluent version: %s. %d geometries to solve.", session.get_fluent_version(), total)

    n_converged = 0
    n_ran_not_converged = 0
    n_failed = 0
    write_header = True
    t_start = time.time()

    try:
        for i, msh_path in enumerate(msh_paths, start=1):
            geometry_id = os.path.splitext(os.path.basename(msh_path))[0]
            json_path = os.path.join(args.mesh_dir, f"{geometry_id}.json")

            logger.info("[%d/%d] %s ...", i, total, geometry_id)

            if not os.path.exists(json_path):
                logger.warning("[%d/%d] %s: no matching .json, skipping", i, total, geometry_id)
                continue

            try:
                with open(json_path) as f:
                    geom_params = json.load(f)
                result = solve_geometry(session, geometry_id, os.path.abspath(msh_path), geom_params, cfg)
            except Exception as e:
                # Belt-and-suspenders: solve_geometry already catches its own
                # errors internally, but a truly unexpected exception (e.g.
                # a lost connection) shouldn't take the whole batch down either.
                logger.error("[%d/%d] %s: unexpected exception, logging and continuing: %s",
                             i, total, geometry_id, e)
                from fluent_solve import SolveResult
                result = SolveResult(geometry_id, False, False, 0, None, [], f"unexpected exception: {e}")

            write_row(args.out_csv, result, write_header)
            write_header = False

            if result.success and result.converged:
                n_converged += 1
            elif result.success:
                n_ran_not_converged += 1
            else:
                n_failed += 1

            elapsed_min = (time.time() - t_start) / 60.0
            logger.info(
                "[%d/%d] %s: %s  (running totals: %d converged, %d ran-but-not-converged, %d failed; %.1f min elapsed)",
                i, total, geometry_id,
                "converged" if result.converged else ("ran, no convergence" if result.success else "FAILED"),
                n_converged, n_ran_not_converged, n_failed, elapsed_min,
            )
    finally:
        session.exit()
        logger.info("Session closed.")

    logger.info(
        "Batch complete: %d/%d converged, %d ran but did not converge, %d failed. Results in %s",
        n_converged, total, n_ran_not_converged, n_failed, args.out_csv,
    )


if __name__ == "__main__":
    main()
