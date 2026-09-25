"""Paired bootstrap intervals for descriptive P2C seed-level differences."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase2/p2c"
N_RESAMPLES = 30000


def main() -> None:
    with (OUT / "seed_metrics.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    lookup = {(int(row["seed"]), int(row["budget"]), row["policy"]): row for row in rows}
    result = {"method": "percentile paired-seed bootstrap of the mean difference",
              "resamples": N_RESAMPLES, "confidence_level": .95, "results": []}
    metrics = ("fleet_rmse_m", "worst_node_rmse_m", "relative_shape_rmse_m")
    for budget in (180, 360, 720):
        for policy_index, policy in enumerate(("geometry", "information")):
            for reference_index, reference in enumerate(("cyclic", "random")):
                for metric_index, metric in enumerate(metrics):
                    differences = np.array([
                        float(lookup[seed, budget, policy][metric]) -
                        float(lookup[seed, budget, reference][metric])
                        for seed in range(20)
                    ])
                    rng = np.random.default_rng(20260925 + 1000 * budget +
                                                100 * policy_index +
                                                10 * reference_index + metric_index)
                    means = np.mean(rng.choice(differences, (N_RESAMPLES, 20), replace=True),
                                    axis=1)
                    lo, hi = np.quantile(means, [.025, .975])
                    result["results"].append({
                        "budget": budget, "policy": policy, "reference": reference,
                        "metric": metric, "mean_difference_m": float(np.mean(differences)),
                        "interval_m": [float(lo), float(hi)],
                        "paired_wins": int(np.count_nonzero(differences < 0)),
                    })
    (OUT / "paired_uncertainty.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
