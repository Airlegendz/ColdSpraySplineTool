"""
Batch driver for the Gmsh meshing pipeline: globs geometry_*.csv files from
a nozzle_geometry.py/sweep.py output directory, meshes each one via
mesh_geometry.mesh_geometry(), and reports a success/failure/quality-flagged
summary -- analogous to sweep.py's rejection-rate reporting, but for
meshing, which is a distinct failure mode from geometry validity: a
geometry can pass every nozzle_geometry.py constraint check and still
produce a poor-quality or failed mesh (e.g. a very tight fillet radius
straining the transfinite radial grading), so this is logged separately
rather than reused from the generator's index.csv.
"""

from __future__ import annotations

import argparse
import csv
import glob
import logging
import os
from typing import List

from mesh_geometry import MeshConfig, MeshResult, mesh_geometry, DEFAULT_N_AXIAL, DEFAULT_N_RADIAL, \
    DEFAULT_FIRST_CELL_HEIGHT, DEFAULT_MESH_FORMAT, DEFAULT_QUALITY_THRESHOLD

logger = logging.getLogger("mesh_sweep")


def find_geometry_csvs(sweep_dir: str) -> List[str]:
    pattern = os.path.join(sweep_dir, "geometry_*.csv")
    return sorted(glob.glob(pattern))


def run_mesh_sweep(sweep_dir: str, out_dir: str, cfg: MeshConfig) -> List[MeshResult]:
    os.makedirs(out_dir, exist_ok=True)
    csv_paths = find_geometry_csvs(sweep_dir)
    if not csv_paths:
        logger.warning("No geometry_*.csv files found in %s", sweep_dir)
        return []

    results: List[MeshResult] = []
    for csv_path in csv_paths:
        geometry_id = os.path.splitext(os.path.basename(csv_path))[0]
        json_path = os.path.join(sweep_dir, f"{geometry_id}.json")
        json_path = json_path if os.path.exists(json_path) else None

        result = mesh_geometry(csv_path, out_dir, geometry_id, cfg, json_path=json_path)
        results.append(result)

        if result.success and not result.quality_flagged:
            logger.info("%s: OK (%d nodes, %d elements)", geometry_id, result.n_nodes, result.n_elements)
        elif result.success and result.quality_flagged:
            logger.warning(
                "%s: meshed but QUALITY FLAGGED (min minSICN=%.4g < %.4g)",
                geometry_id, result.min_quality, cfg.quality_threshold,
            )
        else:
            logger.warning("%s: MESHING FAILED: %s", geometry_id, result.error)

    return results


def write_mesh_index(results: List[MeshResult], out_dir: str) -> None:
    index_path = os.path.join(out_dir, "mesh_index.csv")
    if not results:
        return
    with open(index_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "geometry_id", "success", "quality_flagged", "n_nodes", "n_elements",
            "min_quality", "n_inverted_elements", "error", "msh_path", "png_path",
        ])
        writer.writeheader()
        for r in results:
            writer.writerow({
                "geometry_id": r.geometry_id, "success": r.success, "quality_flagged": r.quality_flagged,
                "n_nodes": r.n_nodes, "n_elements": r.n_elements, "min_quality": r.min_quality,
                "n_inverted_elements": r.n_inverted_elements, "error": r.error or "",
                "msh_path": r.msh_path or "", "png_path": r.png_path or "",
            })


def report_summary(results: List[MeshResult], cfg: MeshConfig) -> None:
    total = len(results)
    succeeded = [r for r in results if r.success]
    failed = [r for r in results if not r.success]
    quality_flagged = [r for r in succeeded if r.quality_flagged]
    clean = [r for r in succeeded if not r.quality_flagged]

    logger.info(
        "Mesh sweep complete: %d/%d meshed (%d clean, %d quality-flagged), %d failed",
        len(succeeded), total, len(clean), len(quality_flagged), len(failed),
    )
    if failed:
        logger.info("Failed geometries:")
        for r in failed:
            logger.info("  %s: %s", r.geometry_id, r.error)
    if quality_flagged:
        logger.info("Quality-flagged geometries (min minSICN < %.4g, meshed but worth a closer look):",
                     cfg.quality_threshold)
        for r in quality_flagged:
            logger.info("  %s: min minSICN=%.4g", r.geometry_id, r.min_quality)

    if total > 0:
        failure_rate = len(failed) / total
        if failure_rate > 0.2:
            logger.warning(
                "Mesh failure rate > 20%% -- mesh_sweep.py's n_axial/n_radial/first_cell_height "
                "config likely needs tuning for this geometry set, not just retrying."
            )


def main():
    parser = argparse.ArgumentParser(
        description="Batch-mesh all geometry_*.csv files from a sweep.py output directory."
    )
    parser.add_argument("--sweep-dir", required=True, help="Directory containing geometry_*.csv (+ .json) files.")
    parser.add_argument("--out-dir", default="mesh_output", help="Output directory for .msh/.png/mesh_index.csv.")
    parser.add_argument("--n-axial", type=int, default=DEFAULT_N_AXIAL,
                         help="Transfinite node count along wall/axis (flow direction). UNCONFIRMED default.")
    parser.add_argument("--n-radial", type=int, default=DEFAULT_N_RADIAL,
                         help="Transfinite node count along inlet/outlet (radial direction). UNCONFIRMED default.")
    parser.add_argument("--first-cell-height", type=float, default=DEFAULT_FIRST_CELL_HEIGHT,
                         help="Target wall-adjacent cell height (m) for boundary-layer clustering. "
                              "UNCONFIRMED default -- depends on flow Reynolds number.")
    parser.add_argument("--mesh-format", choices=["msh2", "msh4"], default=DEFAULT_MESH_FORMAT,
                         help="Gmsh export format version. UNCONFIRMED which your Fluent install expects.")
    parser.add_argument("--quality-threshold", type=float, default=DEFAULT_QUALITY_THRESHOLD,
                         help="Soft minSICN threshold for the quality-flagged warning. UNCONFIRMED default.")
    parser.add_argument("--mesh-quality-check", action="store_true",
                         help="Compute mesh quality metrics after meshing and flag geometries below threshold.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    cfg = MeshConfig(
        n_axial=args.n_axial, n_radial=args.n_radial, first_cell_height=args.first_cell_height,
        mesh_format=args.mesh_format, quality_threshold=args.quality_threshold,
        check_quality=args.mesh_quality_check,
    )

    results = run_mesh_sweep(args.sweep_dir, args.out_dir, cfg)
    write_mesh_index(results, args.out_dir)
    report_summary(results, cfg)


if __name__ == "__main__":
    main()
