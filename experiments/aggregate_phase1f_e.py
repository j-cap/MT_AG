from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

from mt_ag.paper_pf import generate_paper_trajectory
from mt_ag.simulation import auxiliary_trajectory
from run_phase1f_e import aggregate_method, compact_geometry, distribution, observability_summary

ROOT = Path(__file__).resolve().parents[1]


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1f_e.yaml").read_text())
    shard_dir = ROOT / "results/phase1/p1f_e/shards"
    shard_files = sorted(shard_dir.glob("*.json"))
    if not shard_files:
        raise RuntimeError("No P1F-E shard files found")

    by_geometry = {name: [] for name in cfg["geometry_cases"]}
    observed = {}
    for path in shard_files:
        payload = json.loads(path.read_text())
        name = payload["geometry"]
        by_geometry[name].extend(payload["runs"])
        observed.setdefault(name, payload["observability"])
        if observed[name] != payload["observability"]:
            raise RuntimeError(f"Observability payload mismatch across {name} shards")

    expected = cfg["comparison"]["n_seeds"]
    for name, runs in by_geometry.items():
        seeds = sorted(run["seed"] for run in runs)
        if len(runs) != expected or seeds != list(range(expected)):
            raise RuntimeError(
                f"Expected seeds 0..{expected - 1} for {name}, got {seeds}"
            )

    trajectory = generate_paper_trajectory(**cfg["simulation"])
    moving = auxiliary_trajectory(trajectory.t, "moving")
    moving_observability = observability_summary(trajectory, moving, cfg)

    raw_geometry = {}
    compact = {}
    for name, runs in by_geometry.items():
        n_particles = cfg["comparison"]["n_particles"]
        result = {
            "n": len(runs),
            "n_particles": n_particles,
            "observability": observed[name],
            "initial_support_count": distribution([r["initial_support_count"] for r in runs]),
            "initial_support_positive_fraction": float(
                np.mean([r["initial_support_count"] > 0 for r in runs])
            ),
            "pf": aggregate_method(runs, "pf"),
            "aacopf": aggregate_method(runs, "aacopf"),
            "runs": runs,
        }
        raw_geometry[name] = result
        compact[name] = compact_geometry(result)

    aco_cfg = cfg["frozen_aacopf"]
    summary = {
        "frozen_aacopf": {
            "alpha": aco_cfg["alpha"],
            "beta": aco_cfg["beta"],
            "c_lambda": aco_cfg["c_lambda"],
        },
        "negative_controls": compact,
        "moving_observability_reference": moving_observability,
        "acceptance": {
            "criterion": (
                "Frozen AACOPF must not show systematic full-pose recovery in either "
                "locally rank-deficient geometry."
            ),
            "stationary_exact_symmetry_note": (
                "A fixed single auxiliary has an exact global rotational gauge in the "
                "paper-state problem."
            ),
        },
    }
    raw = {
        "design": {
            "purpose": "P1F-E rank-deficient geometry negative controls",
            "geometries": cfg["geometry_cases"],
            "n_seeds": expected,
            "n_particles": cfg["comparison"]["n_particles"],
            "parallel_shards": len(shard_files),
        },
        "negative_controls": raw_geometry,
        "moving_observability_reference": moving_observability,
        "summary": summary,
    }

    out_dir = ROOT / "results/phase1/p1f_e"
    (out_dir / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
