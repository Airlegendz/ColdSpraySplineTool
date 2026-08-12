"""
Batch driver for the cold-spray nozzle geometry sweep.

Reads parameter bounds from a YAML config, draws a Latin hypercube sample,
generates + validates each geometry via nozzle_geometry.generate_geometry,
and writes one CSV + JSON per valid geometry plus a master index CSV.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
from typing import Dict, List

import numpy as np
import yaml

from nozzle_geometry import GeometryConfig, generate_geometry

logger = logging.getLogger("sweep")


def latin_hypercube(n_samples: int, n_dims: int, seed: int) -> np.ndarray:
    """Basic Latin hypercube sample in [0, 1)^n_dims. Uses scipy if available,
    otherwise falls back to a manual LHS implementation."""
    try:
        from scipy.stats import qmc
        sampler = qmc.LatinHypercube(d=n_dims, seed=seed)
        return sampler.random(n=n_samples)
    except ImportError:
        rng = np.random.default_rng(seed)
        result = np.zeros((n_samples, n_dims))
        for d in range(n_dims):
            perm = rng.permutation(n_samples)
            jitter = rng.random(n_samples)
            result[:, d] = (perm + jitter) / n_samples
        return result


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def build_param_sets(cfg: dict, n_samples: int, seed: int) -> List[Dict]:
    bounds: Dict[str, List[float]] = cfg["bounds"]
    fixed: Dict = cfg.get("fixed", {})
    param_names = list(bounds.keys())
    n_dims = len(param_names)

    unit_samples = latin_hypercube(n_samples, n_dims, seed)

    param_sets = []
    for row in unit_samples:
        params = dict(fixed)
        for name, u in zip(param_names, row):
            lo, hi = bounds[name]
            params[name] = float(lo + u * (hi - lo))
        param_sets.append(params)
    return param_sets


def params_to_geometry_config(params: Dict) -> GeometryConfig:
    barrel_count = int(params.get("barrel_count", 0))
    barrel_positions = [
        params[f"barrel_position_{i}"]
        for i in range(barrel_count)
        if f"barrel_position_{i}" in params
    ]

    known_fields = set(GeometryConfig.__dataclass_fields__.keys())
    kwargs = {k: v for k, v in params.items() if k in known_fields}
    kwargs["barrel_count"] = barrel_count
    kwargs["barrel_positions"] = barrel_positions
    return GeometryConfig(**kwargs)


def run_sweep(config_path: str, n_samples: int, seed: int, out_dir: str, make_plots: bool) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    os.makedirs(out_dir, exist_ok=True)

    cfg = load_config(config_path)
    param_sets = build_param_sets(cfg, n_samples, seed)

    index_rows = []
    n_valid = 0
    n_rejected = 0
    rejection_counts: Dict[str, int] = {}

    for i, params in enumerate(param_sets):
        geom_id = f"geometry_{i + 1:04d}"
        geom_cfg = params_to_geometry_config(params)
        result = generate_geometry(geom_cfg)

        row = {"geometry_id": geom_id, "valid": result.valid, "rejection_reason": result.rejection_reason or ""}
        row.update(geom_cfg.to_dict())
        index_rows.append(row)

        if result.valid:
            n_valid += 1
            csv_path = os.path.join(out_dir, f"{geom_id}.csv")
            with open(csv_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["x", "r"])
                writer.writerows(result.points)

            json_path = os.path.join(out_dir, f"{geom_id}.json")
            with open(json_path, "w") as f:
                json.dump(geom_cfg.to_dict(), f, indent=2)

            if make_plots:
                from plotting import plot_profile
                plot_profile(result.points, os.path.join(out_dir, f"{geom_id}.png"), title=geom_id)
        else:
            n_rejected += 1
            # Rejection reasons are tagged with a stable "[category]" prefix
            # by nozzle_geometry.py (e.g. "[curvature] ..."); bucket on that
            # tag rather than the full message so the breakdown groups by
            # constraint type instead of by value-specific text.
            reason = result.rejection_reason or "[unknown] no reason recorded"
            if reason.startswith("[") and "]" in reason:
                reason_key = reason[1:reason.index("]")]
            else:
                reason_key = "unknown"
            rejection_counts[reason_key] = rejection_counts.get(reason_key, 0) + 1
            logger.info("%s REJECTED: %s", geom_id, result.rejection_reason)

    index_path = os.path.join(out_dir, "index.csv")
    if index_rows:
        fieldnames = list(index_rows[0].keys())
        with open(index_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(index_rows)

    total = len(param_sets)
    rejection_rate = n_rejected / total if total else 0.0
    logger.info("Sweep complete: %d/%d valid, %d rejected (%.1f%% rejection rate)",
                n_valid, total, n_rejected, rejection_rate * 100)
    logger.info("Rejection breakdown: %s", rejection_counts)
    if rejection_rate > 0.5:
        logger.warning(
            "Rejection rate > 50%% -- parameter bounds in %s likely need tightening, "
            "not just re-sampling.", config_path
        )


def main():
    parser = argparse.ArgumentParser(description="Latin hypercube sweep of cold-spray nozzle geometries.")
    parser.add_argument("--config", default="config_example.yaml", help="Path to YAML bounds config.")
    parser.add_argument("--n-samples", type=int, required=True, help="Number of geometries to sample.")
    parser.add_argument("--seed", type=int, default=0, help="RNG seed for reproducibility.")
    parser.add_argument("--out-dir", default="sweep_output", help="Output directory.")
    parser.add_argument("--plot", action="store_true", help="Render a PNG for each valid geometry.")
    args = parser.parse_args()

    run_sweep(args.config, args.n_samples, args.seed, args.out_dir, args.plot)


if __name__ == "__main__":
    main()
