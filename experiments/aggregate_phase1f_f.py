from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]


def distribution(values):
    x = np.asarray(values, dtype=float)
    if len(x) == 0:
        return None
    return {
        "mean": float(np.mean(x)),
        "median": float(np.median(x)),
        "std": float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def index_runs(runs):
    return {int(run["seed"]): run for run in runs}


def degradation(clean_runs, stressed_runs, method):
    clean = index_runs(clean_runs)
    stressed = index_runs(stressed_runs)
    common = sorted(set(clean) & set(stressed))
    if len(common) != len(clean) or len(common) != len(stressed):
        raise RuntimeError("Clean/stressed seed sets do not match")

    clean_lock = np.array([clean[s][method]["correct_lock"] for s in common], dtype=bool)
    stressed_lock = np.array([stressed[s][method]["correct_lock"] for s in common], dtype=bool)
    clean_pose = np.array([clean[s][method]["pose_success"] for s in common], dtype=bool)
    stressed_pose = np.array([stressed[s][method]["pose_success"] for s in common], dtype=bool)
    pos_delta = np.array(
        [
            stressed[s][method]["position_rmse_m"]
            - clean[s][method]["position_rmse_m"]
            for s in common
        ],
        dtype=float,
    )
    yaw_delta = np.array(
        [
            stressed[s][method]["yaw_rmse_deg"] - clean[s][method]["yaw_rmse_deg"]
            for s in common
        ],
        dtype=float,
    )
    return {
        "n": len(common),
        "correct_lock_fraction_drop": float(np.mean(clean_lock) - np.mean(stressed_lock)),
        "clean_correct_locks_lost": int(np.sum(clean_lock & ~stressed_lock)),
        "stressed_only_correct_locks": int(np.sum(~clean_lock & stressed_lock)),
        "pose_success_fraction_drop": float(np.mean(clean_pose) - np.mean(stressed_pose)),
        "clean_pose_successes_lost": int(np.sum(clean_pose & ~stressed_pose)),
        "position_rmse_increase_m": distribution(pos_delta),
        "yaw_rmse_increase_deg": distribution(yaw_delta),
    }


def compact(aggregate):
    return {
        "n": aggregate["n"],
        "corruption_count_mean": aggregate["corruption_count"]["mean"],
        "max_abs_additional_error_mean_m": aggregate["max_abs_additional_error_m"]["mean"],
        "pf_correct_lock_fraction": aggregate["pf"]["correct_lock_fraction"],
        "aacopf_correct_lock_fraction": aggregate["aacopf"]["correct_lock_fraction"],
        "pf_wrong_lock_fraction": aggregate["pf"]["wrong_lock_fraction"],
        "aacopf_wrong_lock_fraction": aggregate["aacopf"]["wrong_lock_fraction"],
        "pf_pose_success_fraction": aggregate["pf"]["pose_success_fraction"],
        "aacopf_pose_success_fraction": aggregate["aacopf"]["pose_success_fraction"],
        "pf_position_rmse_mean_m": aggregate["pf"]["position_rmse_m"]["mean"],
        "aacopf_position_rmse_mean_m": aggregate["aacopf"]["position_rmse_m"]["mean"],
        "pf_yaw_rmse_mean_deg": aggregate["pf"]["yaw_rmse_deg"]["mean"],
        "aacopf_yaw_rmse_mean_deg": aggregate["aacopf"]["yaw_rmse_deg"]["mean"],
        "pf_recovery_fraction": aggregate["pf"]["recovery_fraction"],
        "aacopf_recovery_fraction": aggregate["aacopf"]["recovery_fraction"],
        "aacopf_catastrophic_collapse_fraction": aggregate["aacopf"][
            "catastrophic_collapse_fraction"
        ],
        "aacopf_dominant_clone_fraction": aggregate["aacopf"]["dominant_clone_fraction"],
        "paired_lock_outcomes": aggregate["paired_lock_outcomes"],
    }


def make_markdown(summary, cfg):
    intro = (
        "Both methods keep the nominal Gaussian likelihood with "
        f"sigma={cfg['uwb']['assumed_likelihood_sigma_m']:.2f} m. Only the measurement "
        "process changes, so this is not a robust-likelihood comparison."
    )
    lines = [
        "# P1F-F summary — non-Gaussian / outlier UWB stress",
        "",
        intro,
        "",
    ]
    for scenario in cfg["initial_mode_scenarios"]:
        header = (
            "| condition | PF lock | AACOPF lock | PF wrong | AACOPF wrong | "
            "PF pos RMSE [m] | AACOPF pos RMSE [m] | AACOPF collapse |"
        )
        lines.extend(
            [
                f"## Initial cloud: {scenario}",
                "",
                header,
                "|---|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for condition in cfg["ranging_conditions"]:
            value = summary["results"][scenario][condition]
            lines.append(
                f"| {condition} | {value['pf_correct_lock_fraction']:.2f} | "
                f"{value['aacopf_correct_lock_fraction']:.2f} | "
                f"{value['pf_wrong_lock_fraction']:.2f} | "
                f"{value['aacopf_wrong_lock_fraction']:.2f} | "
                f"{value['pf_position_rmse_mean_m']:.3f} | "
                f"{value['aacopf_position_rmse_mean_m']:.3f} | "
                f"{value['aacopf_catastrophic_collapse_fraction']:.2f} |"
            )
        lines.extend(["", "### Degradation relative to clean", ""])
        lines.append(
            "| condition | PF lock drop | AACOPF lock drop | PF clean locks lost | "
            "AACOPF clean locks lost |"
        )
        lines.append("|---|---:|---:|---:|---:|")
        for condition in cfg["ranging_conditions"]:
            if condition == "clean":
                continue
            value = summary["degradation_vs_clean"][scenario][condition]
            lines.append(
                f"| {condition} | {value['pf']['correct_lock_fraction_drop']:.2f} | "
                f"{value['aacopf']['correct_lock_fraction_drop']:.2f} | "
                f"{value['pf']['clean_correct_locks_lost']} | "
                f"{value['aacopf']['clean_correct_locks_lost']} |"
            )
        lines.append("")
    guardrail = (
        "P1F-F tests whether the frozen particle-management rules are inherently tolerant "
        "to corrupted ranges. Any robustness advantage must be visible without changing the "
        "Gaussian likelihood, gating measurements, or retuning AACOPF parameters. If both "
        "methods degrade strongly, the correct conclusion is that explicit robust measurement "
        "handling is needed rather than that AACOPF is a robust-ranging method."
    )
    lines.extend(["## Interpretation guardrail", "", guardrail, ""])
    return "\n".join(lines)


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1f_f.yaml").read_text())
    shard_dir = ROOT / "results/phase1/p1f_f/shards"
    files = sorted(shard_dir.glob("*.json"))
    expected_count = len(cfg["ranging_conditions"]) * len(cfg["initial_mode_scenarios"])
    if len(files) != expected_count:
        raise RuntimeError(f"Expected {expected_count} shard files, found {len(files)}")

    payloads = {}
    for path in files:
        payload = json.loads(path.read_text())
        condition = payload["design"]["condition"]
        scenario = payload["design"]["initial_scenario"]
        key = (scenario, condition)
        if key in payloads:
            raise RuntimeError(f"Duplicate P1F-F shard for {key}")
        payloads[key] = payload

    expected_seeds = list(
        range(
            cfg["validation"]["seed_start"],
            cfg["validation"]["seed_start"] + cfg["validation"]["n_seeds"],
        )
    )
    for key, payload in payloads.items():
        seeds = sorted(int(run["seed"]) for run in payload["runs"])
        if seeds != expected_seeds:
            raise RuntimeError(f"Seed mismatch for {key}: {seeds}")

    results = {scenario: {} for scenario in cfg["initial_mode_scenarios"]}
    degradation_summary = {scenario: {} for scenario in cfg["initial_mode_scenarios"]}
    raw = {scenario: {} for scenario in cfg["initial_mode_scenarios"]}
    for scenario in cfg["initial_mode_scenarios"]:
        clean_runs = payloads[(scenario, "clean")]["runs"]
        for condition in cfg["ranging_conditions"]:
            payload = payloads[(scenario, condition)]
            results[scenario][condition] = compact(payload["aggregate"])
            raw[scenario][condition] = payload
            if condition != "clean":
                degradation_summary[scenario][condition] = {
                    "pf": degradation(clean_runs, payload["runs"], "pf"),
                    "aacopf": degradation(clean_runs, payload["runs"], "aacopf"),
                }

    summary = {
        "frozen_aacopf": {
            "alpha": cfg["frozen_aacopf"]["alpha"],
            "beta": cfg["frozen_aacopf"]["beta"],
            "c_lambda": cfg["frozen_aacopf"]["c_lambda"],
        },
        "claim_boundary": cfg["claim_boundary"]["description"],
        "results": results,
        "degradation_vs_clean": degradation_summary,
    }
    out_dir = ROOT / "results/phase1/p1f_f"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    (out_dir / "raw_results.json").write_text(json.dumps(raw, indent=2))
    (out_dir / "summary.md").write_text(make_markdown(summary, cfg))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
