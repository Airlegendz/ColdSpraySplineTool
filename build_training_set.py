"""
Joins fluent_batch.py's results CSV back to each geometry's parameters
(geometry_XXXX.json from the sweep output) to build the GPR surrogate's
training data: one master CSV of geometry parameters + particle exit
velocity + convergence status.

Non-converged and failed runs are flagged, not silently included or
dropped -- the resulting CSV has a `convergence_status` column so the
surrogate-training step can decide how to handle them (typically:
train only on `converged` rows, but keep the rest visible for
diagnosis rather than discarding them at this stage).

Usage:
    python3 build_training_set.py --fluent-results fluent_batch_results.csv \\
        --sweep-dir sweep_output_v1 --out-csv training_data.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os

logger = logging.getLogger("build_training_set")


def classify_status(row: dict) -> str:
    if row["success"] != "True":
        return "failed"
    if row["converged"] != "True":
        return "ran_not_converged"
    return "converged"


def build(fluent_results_csv: str, sweep_dir: str, out_csv: str) -> None:
    with open(fluent_results_csv, newline="") as f:
        fluent_rows = list(csv.DictReader(f))

    if not fluent_rows:
        logger.warning("%s is empty -- nothing to join.", fluent_results_csv)
        return

    joined_rows = []
    n_missing_params = 0
    for row in fluent_rows:
        geometry_id = row["geometry_id"]
        json_path = os.path.join(sweep_dir, f"{geometry_id}.json")
        if not os.path.exists(json_path):
            n_missing_params += 1
            logger.warning("%s: no matching %s -- excluded from training set (not silently kept with blank params).",
                            geometry_id, json_path)
            continue

        with open(json_path) as f:
            geom_params = json.load(f)

        status = classify_status(row)
        joined = {"geometry_id": geometry_id}
        joined.update(geom_params)
        joined["exit_velocity_m_s"] = row["exit_velocity_m_s"]
        joined["iterations_run"] = row["iterations_run"]
        joined["convergence_status"] = status
        joined["warnings"] = row["warnings"]
        joined["error"] = row["error"]
        joined_rows.append(joined)

    if not joined_rows:
        logger.error("No rows could be joined -- check that --sweep-dir matches the geometries in --fluent-results.")
        return

    fieldnames = list(joined_rows[0].keys())
    # barrel_positions is a list in the JSON -- stringify so it's a single CSV cell.
    for row in joined_rows:
        if isinstance(row.get("barrel_positions"), list):
            row["barrel_positions"] = ";".join(str(v) for v in row["barrel_positions"])

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(joined_rows)

    n_total = len(joined_rows)
    n_converged = sum(1 for r in joined_rows if r["convergence_status"] == "converged")
    n_ran_not_converged = sum(1 for r in joined_rows if r["convergence_status"] == "ran_not_converged")
    n_failed = sum(1 for r in joined_rows if r["convergence_status"] == "failed")

    logger.info("Wrote %s: %d rows joined (%d skipped for missing params)", out_csv, n_total, n_missing_params)
    logger.info(
        "Convergence breakdown: %d converged (%.1f%%), %d ran but did not converge (%.1f%%), %d failed (%.1f%%)",
        n_converged, 100 * n_converged / n_total,
        n_ran_not_converged, 100 * n_ran_not_converged / n_total,
        n_failed, 100 * n_failed / n_total,
    )
    if n_converged < n_total:
        logger.warning(
            "%d/%d rows are NOT clean converged solutions -- the GPR surrogate should train on "
            "convergence_status == 'converged' rows only, or explicitly account for the others, "
            "not treat every row as equally trustworthy.",
            n_total - n_converged, n_total,
        )


def main():
    parser = argparse.ArgumentParser(description="Build the GPR surrogate training CSV from Fluent batch results.")
    parser.add_argument("--fluent-results", default="fluent_batch_results.csv")
    parser.add_argument("--sweep-dir", default="sweep_output_v1", help="Directory with geometry_XXXX.json files.")
    parser.add_argument("--out-csv", default="training_data.csv")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    build(args.fluent_results, args.sweep_dir, args.out_csv)


if __name__ == "__main__":
    main()
