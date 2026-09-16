from __future__ import annotations

import itertools
import json
import time
from pathlib import Path

import numpy as np
import yaml

from mt_ag.aacopf import literal_all_pairs_aco_transition, run_literal_small_aacopf
from mt_ag.geometry import wrap_angle
from mt_ag.paper_pf import (
    generate_paper_trajectory,
    initialize_random_annulus_particles,
    run_paper_bootstrap_pf,
)
from mt_ag.scalable_aacopf import (
    bounded_candidate_aco_transition,
    run_bounded_candidate_aacopf,
)
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory

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


def terminal_hold_time(mask, dt, hold_time_s):
    mask = np.asarray(mask, dtype=bool)
    hold_steps = max(1, round(hold_time_s / dt))
    false_indices = np.flatnonzero(~mask)
    start = 0 if len(false_indices) == 0 else int(false_indices[-1] + 1)
    if len(mask) - start < hold_steps:
        return None
    return float(start * dt)


def controlled_particles(trajectory, auxiliary, cloud_cfg, correct_fraction, rng):
    n_particles = cloud_cfg["n_particles"]
    n_correct = round(correct_fraction * n_particles)
    n_wrong = n_particles - n_correct
    position_std = cloud_cfg["position_std_m"]
    yaw_std = np.deg2rad(cloud_cfg["yaw_std_deg"])
    truth_position = trajectory.state[0, :2]
    truth_yaw = trajectory.state[0, 2]
    vector = truth_position - auxiliary[0]
    range0 = np.linalg.norm(vector)
    truth_bearing = np.arctan2(vector[1], vector[0])
    wrong_bearing = truth_bearing + np.deg2rad(cloud_cfg["wrong_bearing_offset_deg"])
    wrong_center = auxiliary[0] + range0 * np.array(
        [np.cos(wrong_bearing), np.sin(wrong_bearing)]
    )
    wrong_yaw = wrap_angle(truth_yaw + np.deg2rad(cloud_cfg["wrong_yaw_offset_deg"]))
    correct = np.column_stack(
        [
            rng.normal(truth_position[0], position_std, n_correct),
            rng.normal(truth_position[1], position_std, n_correct),
            wrap_angle(rng.normal(truth_yaw, yaw_std, n_correct)),
        ]
    )
    wrong = np.column_stack(
        [
            rng.normal(wrong_center[0], position_std, n_wrong),
            rng.normal(wrong_center[1], position_std, n_wrong),
            wrap_angle(rng.normal(wrong_yaw, yaw_std, n_wrong)),
        ]
    )
    particles = np.vstack([correct, wrong])
    rng.shuffle(particles, axis=0)
    return particles


def controlled_datasets(cfg, split_name):
    trajectory = generate_paper_trajectory(**cfg["simulation"])
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    split = cfg[split_name]
    datasets = []
    for scenario_index, (scenario, scenario_cfg) in enumerate(
        cfg["controlled_scenarios"].items()
    ):
        for offset in range(split["n_seeds"]):
            seed = split["seed_start"] + offset
            range_rng = np.random.default_rng(
                cfg["randomness"]["range_seed_offset"]
                + 10000 * scenario_index
                + seed
            )
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2],
                auxiliary,
                range_rng,
                cfg["uwb"]["sigma_range_m"],
            )
            init_rng = np.random.default_rng(
                cfg["randomness"]["initializer_seed_offset"]
                + 10000 * scenario_index
                + seed
            )
            particles = controlled_particles(
                trajectory,
                auxiliary,
                cfg["controlled_cloud"],
                scenario_cfg["correct_fraction"],
                init_rng,
            )
            datasets.append(
                {
                    "scenario": scenario,
                    "seed": seed,
                    "ranges": ranges,
                    "particles": particles,
                }
            )
    return trajectory, auxiliary, datasets


def mode_and_pose_metrics(result, truth, cloud_cfg, dt, signal, runtime_s):
    position_error = np.linalg.norm(result.estimate[:, :2] - truth[:, :2], axis=1)
    yaw_error = np.abs(wrap_angle(result.estimate[:, 2] - truth[:, 2]))
    correct_time = terminal_hold_time(
        signal >= cloud_cfg["concentration_fraction"], dt, cloud_cfg["hold_time_s"]
    )
    wrong_time = terminal_hold_time(
        signal <= cloud_cfg["wrong_lock_fraction"], dt, cloud_cfg["hold_time_s"]
    )
    pose_time = terminal_hold_time(
        (position_error < cloud_cfg["correct_mode_position_threshold_m"])
        & (yaw_error < np.deg2rad(cloud_cfg["correct_mode_yaw_threshold_deg"])),
        dt,
        cloud_cfg["hold_time_s"],
    )
    metrics = {
        "correct_lock": correct_time is not None,
        "correct_lock_time_s": correct_time,
        "wrong_lock": wrong_time is not None,
        "pose_success": pose_time is not None,
        "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "final_correct_mode_signal": float(signal[-1]),
        "runtime_s": float(runtime_s),
    }
    if hasattr(result, "unique_parent_fraction"):
        update_slice = slice(1, None)
        min_parent = float(np.min(result.unique_parent_fraction[update_slice]))
        max_mult = int(np.max(result.max_destination_multiplicity[update_slice]))
        n_particles = len(result.particles_final)
        metrics.update(
            {
                "minimum_unique_parent_fraction": min_parent,
                "catastrophic_collapse": bool(
                    min_parent < cloud_cfg["catastrophic_unique_parent_fraction"]
                ),
                "maximum_destination_multiplicity_fraction": float(
                    max_mult / n_particles
                ),
                "dominant_clone": bool(
                    max_mult >= cloud_cfg["dominant_clone_fraction"] * n_particles
                ),
                "mean_moved_fraction": float(
                    np.mean(result.moved_fraction[update_slice])
                ),
                "mean_transition_runtime_s": float(
                    np.mean(result.transition_runtime_s[update_slice])
                ),
            }
        )
    return metrics


def evaluate_bounded_dataset(dataset, setting, cfg, trajectory, auxiliary):
    aco_cfg = cfg["frozen_aacopf"]
    cloud_cfg = cfg["controlled_cloud"]
    start = time.perf_counter()
    result = run_bounded_candidate_aacopf(
        trajectory.increments,
        dataset["ranges"],
        auxiliary,
        dataset["particles"].copy(),
        sigma_uwb_m=cfg["uwb"]["sigma_range_m"],
        propagation_convention=cfg["paper_model"]["propagation_convention"],
        candidate_count=setting["candidate_count"],
        max_move_fraction=setting["max_move_fraction"],
        destination_capacity_fraction=setting["destination_capacity_fraction"],
        alpha=aco_cfg["alpha"],
        beta=aco_cfg["beta"],
        c_lambda=aco_cfg["c_lambda"],
        epsilon_distance=aco_cfg["epsilon_distance"],
        epsilon_weight=aco_cfg["epsilon_weight"],
        truth_state=trajectory.state,
        correct_mode_position_threshold_m=cloud_cfg[
            "correct_mode_position_threshold_m"
        ],
        correct_mode_yaw_threshold_rad=np.deg2rad(
            cloud_cfg["correct_mode_yaw_threshold_deg"]
        ),
    )
    runtime = time.perf_counter() - start
    return mode_and_pose_metrics(
        result,
        trajectory.state,
        cloud_cfg,
        trajectory.dt,
        result.correct_mode_fraction_post_aco,
        runtime,
    )


def aggregate_controlled_runs(runs):
    return {
        "n": len(runs),
        "correct_lock_fraction": float(np.mean([r["correct_lock"] for r in runs])),
        "wrong_lock_fraction": float(np.mean([r["wrong_lock"] for r in runs])),
        "pose_success_fraction": float(np.mean([r["pose_success"] for r in runs])),
        "catastrophic_collapse_fraction": float(
            np.mean([r.get("catastrophic_collapse", False) for r in runs])
        ),
        "dominant_clone_fraction": float(
            np.mean([r.get("dominant_clone", False) for r in runs])
        ),
        "position_rmse_m": distribution([r["position_rmse_m"] for r in runs]),
        "yaw_rmse_deg": distribution([r["yaw_rmse_deg"] for r in runs]),
        "runtime_s": distribution([r["runtime_s"] for r in runs]),
        "minimum_unique_parent_fraction": distribution(
            [r["minimum_unique_parent_fraction"] for r in runs]
        )
        if "minimum_unique_parent_fraction" in runs[0]
        else None,
        "mean_moved_fraction": distribution([r["mean_moved_fraction"] for r in runs])
        if "mean_moved_fraction" in runs[0]
        else None,
    }


def development_sweep(cfg):
    trajectory, auxiliary, datasets = controlled_datasets(cfg, "development")
    entries = []
    for candidate_count, max_move, destination_capacity in itertools.product(
        cfg["development"]["candidate_counts"],
        cfg["development"]["max_move_fractions"],
        cfg["development"]["destination_capacity_fractions"],
    ):
        setting = {
            "candidate_count": candidate_count,
            "max_move_fraction": max_move,
            "destination_capacity_fraction": destination_capacity,
        }
        runs = []
        for dataset in datasets:
            metrics = evaluate_bounded_dataset(
                dataset, setting, cfg, trajectory, auxiliary
            )
            runs.append(
                {
                    "scenario": dataset["scenario"],
                    "seed": dataset["seed"],
                    **metrics,
                }
            )
        entries.append(
            {
                "setting": setting,
                "aggregate": aggregate_controlled_runs(runs),
                "runs": runs,
            }
        )

    success_floor = cfg["development"]["success_floor"]
    wrong_ceiling = cfg["development"]["wrong_lock_ceiling"]

    def rank_key(entry):
        agg = entry["aggregate"]
        eligible = (
            agg["correct_lock_fraction"] >= success_floor
            and agg["wrong_lock_fraction"] <= wrong_ceiling
        )
        setting = entry["setting"]
        if eligible:
            return (
                0,
                agg["catastrophic_collapse_fraction"],
                agg["dominant_clone_fraction"],
                setting["candidate_count"],
                setting["max_move_fraction"],
                setting["destination_capacity_fraction"],
                -agg["correct_lock_fraction"],
            )
        return (
            1,
            -agg["correct_lock_fraction"],
            agg["wrong_lock_fraction"],
            agg["catastrophic_collapse_fraction"],
            setting["candidate_count"],
        )

    entries.sort(key=rank_key)
    selected = entries[0]
    selected["eligible"] = (
        selected["aggregate"]["correct_lock_fraction"] >= success_floor
        and selected["aggregate"]["wrong_lock_fraction"] <= wrong_ceiling
    )
    return {
        "selected": selected,
        "ranked_settings": entries,
    }


def controlled_validation(cfg, selected_setting):
    trajectory, auxiliary, datasets = controlled_datasets(
        cfg, "controlled_validation"
    )
    methods = {
        "candidate_only": {
            **selected_setting,
            "max_move_fraction": 1.0,
            "destination_capacity_fraction": 1.0,
        },
        "move_guard_only": {
            **selected_setting,
            "destination_capacity_fraction": 1.0,
        },
        "destination_guard_only": {
            **selected_setting,
            "max_move_fraction": 1.0,
        },
        "guarded": dict(selected_setting),
    }
    by_method = {name: [] for name in ["pf", "literal", *methods]}
    raw = []
    cloud_cfg = cfg["controlled_cloud"]
    aco_cfg = cfg["frozen_aacopf"]
    for dataset in datasets:
        pf_rng = np.random.default_rng(
            cfg["randomness"]["pf_seed_offset"] + dataset["seed"]
        )
        start = time.perf_counter()
        pf = run_paper_bootstrap_pf(
            trajectory.increments,
            dataset["ranges"],
            auxiliary,
            dataset["particles"].copy(),
            pf_rng,
            sigma_uwb_m=cfg["uwb"]["sigma_range_m"],
            propagation_convention=cfg["paper_model"]["propagation_convention"],
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=cloud_cfg[
                "correct_mode_position_threshold_m"
            ],
            correct_mode_yaw_threshold_rad=np.deg2rad(
                cloud_cfg["correct_mode_yaw_threshold_deg"]
            ),
        )
        pf_runtime = time.perf_counter() - start
        pf_metrics = mode_and_pose_metrics(
            pf,
            trajectory.state,
            cloud_cfg,
            trajectory.dt,
            pf.correct_mode_mass,
            pf_runtime,
        )
        by_method["pf"].append(pf_metrics)

        start = time.perf_counter()
        literal = run_literal_small_aacopf(
            trajectory.increments,
            dataset["ranges"],
            auxiliary,
            dataset["particles"].copy(),
            sigma_uwb_m=cfg["uwb"]["sigma_range_m"],
            propagation_convention=cfg["paper_model"]["propagation_convention"],
            alpha=aco_cfg["alpha"],
            beta=aco_cfg["beta"],
            c_lambda=aco_cfg["c_lambda"],
            epsilon_distance=aco_cfg["epsilon_distance"],
            epsilon_weight=aco_cfg["epsilon_weight"],
            truth_state=trajectory.state,
            correct_mode_position_threshold_m=cloud_cfg[
                "correct_mode_position_threshold_m"
            ],
            correct_mode_yaw_threshold_rad=np.deg2rad(
                cloud_cfg["correct_mode_yaw_threshold_deg"]
            ),
        )
        literal_runtime = time.perf_counter() - start
        literal_metrics = mode_and_pose_metrics(
            literal,
            trajectory.state,
            cloud_cfg,
            trajectory.dt,
            literal.correct_mode_fraction_post_aco,
            literal_runtime,
        )
        by_method["literal"].append(literal_metrics)

        method_metrics = {"pf": pf_metrics, "literal": literal_metrics}
        for name, setting in methods.items():
            metrics = evaluate_bounded_dataset(
                dataset, setting, cfg, trajectory, auxiliary
            )
            by_method[name].append(metrics)
            method_metrics[name] = metrics
        raw.append(
            {
                "scenario": dataset["scenario"],
                "seed": dataset["seed"],
                "methods": method_metrics,
            }
        )

    summary = {
        name: aggregate_controlled_runs(runs) for name, runs in by_method.items()
    }
    return {"summary": summary, "runs": raw}


def initial_support_count(particles, truth, position_threshold, yaw_threshold_rad):
    pos = np.linalg.norm(particles[:, :2] - truth[:2], axis=1)
    yaw = np.abs(wrap_angle(particles[:, 2] - truth[2]))
    return int(np.sum((pos < position_threshold) & (yaw < yaw_threshold_rad)))


def global_pose_metrics(estimate, truth, cfg, dt):
    position_error = np.linalg.norm(estimate[:, :2] - truth[:, :2], axis=1)
    yaw_error = np.abs(wrap_angle(estimate[:, 2] - truth[:, 2]))
    pose_ok = (
        position_error < cfg["convergence_position_threshold_m"]
    ) & (yaw_error < np.deg2rad(cfg["convergence_yaw_threshold_deg"]))
    convergence_time = terminal_hold_time(
        pose_ok, dt, cfg["convergence_hold_time_s"]
    )
    late_steps = max(1, round(cfg["late_window_s"] / dt))
    return {
        "pose_success": convergence_time is not None,
        "convergence_time_s": convergence_time,
        "position_rmse_m": float(np.sqrt(np.mean(position_error**2))),
        "yaw_rmse_deg": float(np.rad2deg(np.sqrt(np.mean(yaw_error**2)))),
        "late_position_rmse_m": float(
            np.sqrt(np.mean(position_error[-late_steps:] ** 2))
        ),
        "late_yaw_rmse_deg": float(
            np.rad2deg(np.sqrt(np.mean(yaw_error[-late_steps:] ** 2)))
        ),
        "final_position_error_m": float(position_error[-1]),
        "final_yaw_error_deg": float(np.rad2deg(yaw_error[-1])),
    }


def random_annulus_validation(cfg, selected_setting):
    global_cfg = cfg["random_annulus_validation"]
    trajectory = generate_paper_trajectory(
        dt=cfg["simulation"]["dt"], duration=global_cfg["duration"]
    )
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    sigma = cfg["uwb"]["sigma_range_m"]
    delta_d = global_cfg["initialization_delta_d_sigma_factor"] * sigma
    yaw_threshold = np.deg2rad(global_cfg["convergence_yaw_threshold_deg"])
    aco_cfg = cfg["frozen_aacopf"]
    raw = {str(budget): [] for budget in global_cfg["particle_budgets"]}

    for budget in global_cfg["particle_budgets"]:
        for seed in range(
            global_cfg["seed_start"],
            global_cfg["seed_start"] + global_cfg["n_seeds"],
        ):
            range_rng = np.random.default_rng(
                cfg["randomness"]["range_seed_offset"] + seed
            )
            ranges = generate_uwb_ranges(
                trajectory.state[:, :2], auxiliary, range_rng, sigma
            )
            init_rng = np.random.default_rng(
                cfg["randomness"]["initializer_seed_offset"] + seed
            )
            particles = initialize_random_annulus_particles(
                ranges[0], auxiliary[0], budget, init_rng, delta_d
            )
            support = initial_support_count(
                particles,
                trajectory.state[0],
                global_cfg["convergence_position_threshold_m"],
                yaw_threshold,
            )

            pf_rng = np.random.default_rng(
                cfg["randomness"]["pf_seed_offset"] + 100000 * budget + seed
            )
            start = time.perf_counter()
            pf = run_paper_bootstrap_pf(
                trajectory.increments,
                ranges,
                auxiliary,
                particles.copy(),
                pf_rng,
                sigma_uwb_m=sigma,
                resample_fraction=global_cfg["pf_resample_fraction"],
                propagation_convention=cfg["paper_model"]["propagation_convention"],
                truth_state=trajectory.state,
                correct_mode_position_threshold_m=global_cfg[
                    "convergence_position_threshold_m"
                ],
                correct_mode_yaw_threshold_rad=yaw_threshold,
            )
            pf_runtime = time.perf_counter() - start
            pf_metrics = global_pose_metrics(
                pf.estimate, trajectory.state, global_cfg, trajectory.dt
            )
            pf_metrics["runtime_s"] = float(pf_runtime)

            start = time.perf_counter()
            guarded = run_bounded_candidate_aacopf(
                trajectory.increments,
                ranges,
                auxiliary,
                particles.copy(),
                sigma_uwb_m=sigma,
                propagation_convention=cfg["paper_model"]["propagation_convention"],
                **selected_setting,
                alpha=aco_cfg["alpha"],
                beta=aco_cfg["beta"],
                c_lambda=aco_cfg["c_lambda"],
                epsilon_distance=aco_cfg["epsilon_distance"],
                epsilon_weight=aco_cfg["epsilon_weight"],
                truth_state=trajectory.state,
                correct_mode_position_threshold_m=global_cfg[
                    "convergence_position_threshold_m"
                ],
                correct_mode_yaw_threshold_rad=yaw_threshold,
            )
            guarded_runtime = time.perf_counter() - start
            guarded_metrics = global_pose_metrics(
                guarded.estimate, trajectory.state, global_cfg, trajectory.dt
            )
            guarded_metrics.update(
                {
                    "runtime_s": float(guarded_runtime),
                    "minimum_unique_parent_fraction": float(
                        np.min(guarded.unique_parent_fraction[1:])
                    ),
                    "maximum_destination_multiplicity_fraction": float(
                        np.max(guarded.max_destination_multiplicity[1:]) / budget
                    ),
                }
            )

            record = {
                "seed": seed,
                "initial_support_count": support,
                "pf": pf_metrics,
                "guarded": guarded_metrics,
            }
            if budget == global_cfg["candidate_only_at_budget"]:
                candidate_only_setting = {
                    **selected_setting,
                    "max_move_fraction": 1.0,
                    "destination_capacity_fraction": 1.0,
                }
                start = time.perf_counter()
                candidate_only = run_bounded_candidate_aacopf(
                    trajectory.increments,
                    ranges,
                    auxiliary,
                    particles.copy(),
                    sigma_uwb_m=sigma,
                    propagation_convention=cfg["paper_model"]["propagation_convention"],
                    **candidate_only_setting,
                    alpha=aco_cfg["alpha"],
                    beta=aco_cfg["beta"],
                    c_lambda=aco_cfg["c_lambda"],
                    epsilon_distance=aco_cfg["epsilon_distance"],
                    epsilon_weight=aco_cfg["epsilon_weight"],
                    truth_state=trajectory.state,
                    correct_mode_position_threshold_m=global_cfg[
                        "convergence_position_threshold_m"
                    ],
                    correct_mode_yaw_threshold_rad=yaw_threshold,
                )
                candidate_runtime = time.perf_counter() - start
                candidate_metrics = global_pose_metrics(
                    candidate_only.estimate, trajectory.state, global_cfg, trajectory.dt
                )
                candidate_metrics.update(
                    {
                        "runtime_s": float(candidate_runtime),
                        "minimum_unique_parent_fraction": float(
                            np.min(candidate_only.unique_parent_fraction[1:])
                        ),
                        "maximum_destination_multiplicity_fraction": float(
                            np.max(candidate_only.max_destination_multiplicity[1:])
                            / budget
                        ),
                    }
                )
                record["candidate_only"] = candidate_metrics
            raw[str(budget)].append(record)

    summary = {}
    for budget_text, runs in raw.items():
        methods = ["pf", "guarded"]
        if "candidate_only" in runs[0]:
            methods.append("candidate_only")
        method_summary = {}
        for method in methods:
            values = [r[method] for r in runs]
            method_summary[method] = {
                "pose_success_fraction": float(
                    np.mean([v["pose_success"] for v in values])
                ),
                "pose_success_count": int(np.sum([v["pose_success"] for v in values])),
                "late_position_rmse_m": distribution(
                    [v["late_position_rmse_m"] for v in values]
                ),
                "late_yaw_rmse_deg": distribution(
                    [v["late_yaw_rmse_deg"] for v in values]
                ),
                "runtime_s": distribution([v["runtime_s"] for v in values]),
            }
            if method != "pf":
                method_summary[method]["minimum_unique_parent_fraction"] = distribution(
                    [v["minimum_unique_parent_fraction"] for v in values]
                )
                method_summary[method][
                    "maximum_destination_multiplicity_fraction"
                ] = distribution(
                    [v["maximum_destination_multiplicity_fraction"] for v in values]
                )
        summary[budget_text] = {
            "n": len(runs),
            "initial_support_count": distribution(
                [r["initial_support_count"] for r in runs]
            ),
            "initial_support_positive_fraction": float(
                np.mean([r["initial_support_count"] > 0 for r in runs])
            ),
            "methods": method_summary,
        }
    return {"summary": summary, "runs": raw}


def runtime_scaling(cfg, selected_setting):
    rng = np.random.default_rng(701001)
    bounded = []
    for budget in cfg["runtime_scaling"]["particle_budgets"]:
        particles = np.column_stack(
            [
                rng.normal(size=budget),
                rng.normal(size=budget),
                rng.uniform(-np.pi, np.pi, size=budget),
            ]
        )
        weights = rng.lognormal(mean=0.0, sigma=1.0, size=budget)
        runtimes = []
        score_counts = []
        for _ in range(cfg["runtime_scaling"]["repeats"]):
            _, _, _, diag = bounded_candidate_aco_transition(
                particles,
                weights,
                **selected_setting,
                alpha=cfg["frozen_aacopf"]["alpha"],
                beta=cfg["frozen_aacopf"]["beta"],
                c_lambda=cfg["frozen_aacopf"]["c_lambda"],
                epsilon_distance=cfg["frozen_aacopf"]["epsilon_distance"],
                epsilon_weight=cfg["frozen_aacopf"]["epsilon_weight"],
            )
            runtimes.append(diag.runtime_s)
            score_counts.append(diag.candidate_score_count)
        bounded.append(
            {
                "n_particles": budget,
                "runtime_s": distribution(runtimes),
                "candidate_score_count": distribution(score_counts),
                "dense_pair_count": budget * budget,
                "score_fraction_of_dense_pairs": float(
                    np.mean(score_counts) / (budget * budget)
                ),
            }
        )

    literal = []
    for budget in cfg["runtime_scaling"]["literal_reference_budgets"]:
        particles = np.column_stack(
            [
                rng.normal(size=budget),
                rng.normal(size=budget),
                rng.uniform(-np.pi, np.pi, size=budget),
            ]
        )
        weights = rng.lognormal(mean=0.0, sigma=1.0, size=budget)
        runtimes = []
        for _ in range(cfg["runtime_scaling"]["repeats"]):
            _, _, _, diag = literal_all_pairs_aco_transition(
                particles,
                weights,
                alpha=cfg["frozen_aacopf"]["alpha"],
                beta=cfg["frozen_aacopf"]["beta"],
                c_lambda=cfg["frozen_aacopf"]["c_lambda"],
                epsilon_distance=cfg["frozen_aacopf"]["epsilon_distance"],
                epsilon_weight=cfg["frozen_aacopf"]["epsilon_weight"],
            )
            runtimes.append(diag.runtime_s)
        literal.append(
            {
                "n_particles": budget,
                "runtime_s": distribution(runtimes),
            }
        )

    def slope(entries):
        x = np.log([e["n_particles"] for e in entries])
        y = np.log([max(e["runtime_s"]["median"], 1e-12) for e in entries])
        return float(np.polyfit(x, y, 1)[0])

    return {
        "bounded": bounded,
        "literal": literal,
        "bounded_loglog_runtime_slope": slope(bounded),
        "literal_loglog_runtime_slope": slope(literal),
    }


def make_markdown(summary):
    selected = summary["development"]["selected"]["setting"]
    controlled = summary["controlled_validation"]["summary"]
    global_summary = summary["random_annulus_validation"]["summary"]
    lines = [
        "# P1F-G summary — scalable bounded-candidate AACOPF",
        "",
        "Selected development setting:",
        "",
        f"- candidate_count: {selected['candidate_count']}",
        f"- max_move_fraction: {selected['max_move_fraction']}",
        f"- destination_capacity_fraction: {selected['destination_capacity_fraction']}",
        "",
        "## Controlled held-out validation",
        "",
        "| method | correct lock | wrong lock | collapse | dominant clone | runtime [s] |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for method, value in controlled.items():
        lines.append(
            f"| {method} | {value['correct_lock_fraction']:.3f} | "
            f"{value['wrong_lock_fraction']:.3f} | "
            f"{value['catastrophic_collapse_fraction']:.3f} | "
            f"{value['dominant_clone_fraction']:.3f} | "
            f"{value['runtime_s']['mean']:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Random-annulus global validation",
            "",
            "| Np | method | pose success | late pos RMSE [m] | late yaw RMSE [deg] | runtime [s] |",
            "|---:|---|---:|---:|---:|---:|",
        ]
    )
    for budget, value in global_summary.items():
        for method, metrics in value["methods"].items():
            lines.append(
                f"| {budget} | {method} | {metrics['pose_success_fraction']:.3f} | "
                f"{metrics['late_position_rmse_m']['mean']:.3f} | "
                f"{metrics['late_yaw_rmse_deg']['mean']:.3f} | "
                f"{metrics['runtime_s']['mean']:.3f} |"
            )
    lines.extend(
        [
            "",
            "## Runtime scaling",
            "",
            f"Bounded empirical log-log runtime slope: {summary['runtime_scaling']['bounded_loglog_runtime_slope']:.3f}",
            "",
            f"Literal empirical log-log runtime slope: {summary['runtime_scaling']['literal_loglog_runtime_slope']:.3f}",
            "",
        ]
    )
    return "\n".join(lines)


def main():
    cfg = yaml.safe_load((ROOT / "configs/phase1f_g.yaml").read_text())
    development = development_sweep(cfg)
    selected_setting = development["selected"]["setting"]
    controlled = controlled_validation(cfg, selected_setting)
    global_validation = random_annulus_validation(cfg, selected_setting)
    scaling = runtime_scaling(cfg, selected_setting)
    summary = {
        "design": {
            "purpose": "P1F-G scalable bounded-candidate AACOPF adaptation",
            "frozen_aacopf": cfg["frozen_aacopf"],
            "selection_policy": cfg["development"]["selection_policy"],
        },
        "development": development,
        "controlled_validation": controlled,
        "random_annulus_validation": global_validation,
        "runtime_scaling": scaling,
    }
    out_dir = ROOT / "results/phase1/p1f_g"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "raw_results.json").write_text(json.dumps(summary, indent=2))
    compact = {
        "design": summary["design"],
        "development": {
            "selected": development["selected"],
            "n_settings": len(development["ranked_settings"]),
        },
        "controlled_validation": {"summary": controlled["summary"]},
        "random_annulus_validation": {"summary": global_validation["summary"]},
        "runtime_scaling": scaling,
    }
    (out_dir / "summary.json").write_text(json.dumps(compact, indent=2))
    (out_dir / "summary.md").write_text(make_markdown(compact))
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
