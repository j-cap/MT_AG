from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "report" / "figures" / "phase1"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.size": 8.5,
        "axes.labelsize": 8.5,
        "axes.titlesize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def _panel_label(ax, label):
    ax.text(
        0.01,
        0.99,
        label,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontweight="bold",
    )


def _load_json(path):
    return json.loads((ROOT / path).read_text())


def fig_p1a():
    df = pd.read_csv(ROOT / "results" / "phase1" / "trajectory_report.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax = axes[0]
    ax.plot(df.x_true, df.y_true, linewidth=1.8, label="Truth")
    ax.plot(df.x_dr, df.y_dr, linestyle="--", linewidth=1.3, label="Noisy IMU DR")
    ax.plot(df.x_pf, df.y_pf, linestyle="-.", linewidth=1.3, label="PF + UWB")
    ax.scatter([df.x_true.iloc[0]], [df.y_true.iloc[0]], s=18, marker="o", label="Start")
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    _panel_label(ax, "(a)")

    err_dr = np.hypot(df.x_dr - df.x_true, df.y_dr - df.y_true)
    err_pf = np.hypot(df.x_pf - df.x_true, df.y_pf - df.y_true)
    ax = axes[1]
    ax.plot(df.t, err_dr, linestyle="--", linewidth=1.3, label="Noisy IMU DR")
    ax.plot(df.t, err_pf, linewidth=1.3, label="PF + UWB")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position error [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / "p1a_baseline_trajectory.pdf")
    plt.close(fig)


def fig_p1b():
    summary = _load_json("results/phase1/p1b/summary.json")
    init_order = ["small_offset", "moderate_offset", "large_offset", "boundary_offset"]
    init_labels = ["Small", "Moderate", "Large", "Boundary"]
    init_pf = [summary["initial_state"][key]["pf_rmse_mean_m"] for key in init_order]
    init_dr = [summary["initial_state"][key]["dr_rmse_mean_m"] for key in init_order]

    particle = summary["particle_count"]
    counts = np.array(sorted(int(key) for key in particle))
    pf_rmse = np.array([particle[str(n)]["pf_rmse_mean_m"] for n in counts])
    runtime = np.array([particle[str(n)]["runtime_mean_s"] for n in counts])

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax = axes[0]
    x = np.arange(len(init_order))
    width = 0.36
    ax.bar(x - width / 2, init_dr, width, label="Dead reckoning")
    ax.bar(x + width / 2, init_pf, width, label="PF + UWB")
    ax.set_xticks(x, init_labels, rotation=15)
    ax.set_ylabel("Mean position RMSE [m]")
    ax.set_xlabel("Initial-state uncertainty")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, "(a)")

    ax = axes[1]
    ax.plot(counts, pf_rmse, marker="o", linewidth=1.3, label="PF RMSE")
    ax.set_xscale("log")
    ax.set_xlabel("Particle count $N_p$")
    ax.set_ylabel("Mean position RMSE [m]")
    ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(counts, runtime, marker="s", linestyle="--", linewidth=1.1, label="Runtime")
    ax2.set_ylabel("Mean runtime [s]")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="center right")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / "p1b_known_pose_robustness.pdf")
    plt.close(fig)


def fig_p1c():
    summary = _load_json("results/phase1/p1c/summary.json")
    grid = summary["grid_resolution"]
    g = pd.DataFrame(
        [
            {
                "n_particles": value["n_particles"],
                "success": value["pose_success_fraction"],
                "rmse": value["late_position_rmse_mean_m"],
            }
            for _, value in sorted(grid.items(), key=lambda kv: kv[1]["n_particles"])
        ]
    )
    res = summary["resampling_sensitivity"]
    r = pd.DataFrame(
        [
            {
                "gamma": float(key),
                "success": value["pose_success_fraction"],
                "rmse": value["late_position_rmse_mean_m"],
            }
            for key, value in sorted(res.items(), key=lambda kv: float(kv[0]))
        ]
    )

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax = axes[0]
    ax.plot(g.n_particles, g.rmse, marker="o", linewidth=1.3, label="Late position RMSE")
    ax.set_xscale("log")
    ax.set_xlabel("Particle count $N_p$")
    ax.set_ylabel("Late position RMSE [m]")
    ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(
        g.n_particles,
        g.success,
        marker="s",
        linestyle="--",
        linewidth=1.1,
        label="Pose success",
    )
    ax2.set_ylabel("Pose success fraction")
    ax2.set_ylim(-0.02, 0.24)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="upper right")
    _panel_label(ax, "(a)")

    ax = axes[1]
    ax.plot(r.gamma, r.rmse, marker="o", linewidth=1.3, label="Late position RMSE")
    ax.set_xlabel("Resampling threshold $\\gamma$")
    ax.set_ylabel("Late position RMSE [m]")
    ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(r.gamma, r.success, marker="s", linestyle="--", linewidth=1.1, label="Pose success")
    ax2.set_ylabel("Pose success fraction")
    ax2.set_ylim(-0.02, 0.24)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="lower right")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.3)
    fig.savefig(FIG / "p1c_global_search_sensitivity.pdf")
    plt.close(fig)


def fig_p1d():
    df = pd.read_csv(ROOT / "results" / "phase1" / "p1d" / "diagnostics_report.csv")
    order = ["stationary", "constant_bearing", "moving"]
    labels = {
        "stationary": "Stationary",
        "constant_bearing": "Constant bearing",
        "moving": "Moving",
    }
    styles = ["--", ":", "-"]
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax = axes[0]
    for geom, ls in zip(order, styles):
        data = df[df.geometry == geom]
        ax.semilogy(
            data.time_s,
            np.maximum(data.obs_sigma_ratio, 1e-18),
            linestyle=ls,
            linewidth=1.4,
            label=labels[geom],
        )
    ax.axhline(1e-4, linestyle="-.", linewidth=0.9, label="$10^{-4}$ diagnostic level")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("$\\eta_{\\mathcal{O}}=\\sigma_{\\min}/\\sigma_{\\max}$")
    ax.set_ylim(1e-18, 2e-3)
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(loc="lower right")
    _panel_label(ax, "(a)")

    ax = axes[1]
    for geom, ls in zip(order, styles):
        data = df[df.geometry == geom]
        ax.plot(
            data.time_s,
            data.position_error_mean_m,
            linestyle=ls,
            linewidth=1.4,
            label=labels[geom],
        )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mean position error [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / "p1d_geometry_observability.pdf")
    plt.close(fig)


def _flow_box(ax, xy, width, height, title, lines):
    box = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.02",
        linewidth=1.0,
        fill=False,
    )
    ax.add_patch(box)
    x, y = xy
    ax.text(x + width / 2, y + height * 0.72, title, ha="center", va="center", fontweight="bold")
    ax.text(x + width / 2, y + height * 0.36, lines, ha="center", va="center", fontsize=7.5)


def fig_p1e():
    fig, ax = plt.subplots(figsize=(6.7, 2.7))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    boxes = [
        (0.02, 0.24, 0.2, 0.5, "Source paper", "$[x,y,\\phi]$ state\n$\\Delta L,\\Delta\\phi$ input\nAACOPF concept"),
        (0.28, 0.24, 0.2, 0.5, "Audit", "Propagation conflict\nlikelihood ambiguity\n$\\alpha,\\beta,\\lambda$ missing\nrejection rule unclear"),
        (0.54, 0.24, 0.2, 0.5, "Repository contract", "Gaussian likelihood\nhigher-weight candidates\nsynchronous copy\nuniform reset"),
        (0.80, 0.24, 0.18, 0.5, "Evaluation path", "literal-small\ntune/freeze\nmatched tests\nscalable adaptation"),
    ]
    for x, y, width, height, title, lines in boxes:
        _flow_box(ax, (x, y), width, height, title, lines)
    for left, right in [(0.22, 0.28), (0.48, 0.54), (0.74, 0.80)]:
        ax.annotate(
            "",
            xy=(right, 0.49),
            xytext=(left, 0.49),
            arrowprops={"arrowstyle": "->", "linewidth": 1.1},
        )
    ax.text(
        0.5,
        0.08,
        "P1E fixes a transparent qualitative reproduction, not a claim of exact numerical replication.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(FIG / "p1e_reproduction_architecture.pdf")
    plt.close(fig)


def fig_p1fa():
    summary = _load_json("results/phase1/p1f_a/summary.json")
    particle = summary["particle_count"]
    counts = np.array(sorted(int(key) for key in particle))
    success = np.array([particle[str(n)]["pose_success_fraction"] for n in counts])
    rmse = np.array([particle[str(n)]["late_position_rmse_mean_m"] for n in counts])
    support = summary["initial_support"]["random_10000"]

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax = axes[0]
    ax.plot(counts, rmse, marker="o", linewidth=1.3, label="Late position RMSE")
    ax.set_xscale("log")
    ax.set_xlabel("Particle count $N_p$")
    ax.set_ylabel("Late position RMSE [m]")
    ax.grid(True, alpha=0.25)
    ax2 = ax.twinx()
    ax2.plot(counts, success, marker="s", linestyle="--", linewidth=1.1, label="Pose success")
    ax2.set_ylabel("Pose success fraction")
    ax2.set_ylim(0, 0.5)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="upper center")
    _panel_label(ax, "(a)")

    ax = axes[1]
    labels = ["Initial support", "Terminal success"]
    values = [
        support["fraction_with_at_least_one_correct_mode_particle"],
        summary["core_random"]["pose_success_fraction"],
    ]
    ax.bar(labels, values)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction of runs")
    ax.grid(True, axis="y", alpha=0.25)
    ax.text(
        0.5,
        0.52,
        f"Mean correct-region particles at $N_p=10^4$: {support['correct_mode_count_mean']:.1f}",
        ha="center",
        transform=ax.transAxes,
        fontsize=7.5,
    )
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.4)
    fig.savefig(FIG / "p1fa_global_pf_scaling.pdf")
    plt.close(fig)


def fig_p1fb():
    summary = _load_json("results/phase1/p1f_b/summary.json")
    runtime_rows = summary["runtime_scaling"]["rows"]
    counts = np.array([row["n_particles"] for row in runtime_rows])
    runtime_ms = np.array([1000 * row["runtime_median_s"] for row in runtime_rows])

    correct = summary["representative_trace_events"]["seed_0_correct_lock"]
    wrong = summary["representative_trace_events"]["seed_4_wrong_lock"]
    event_labels = ["Correct-lock event", "Wrong-lock event"]
    moved = [correct["moved_fraction"], wrong["moved_fraction_at_lock"]]
    unique = [correct["unique_parent_fraction"], wrong["unique_parent_fraction_at_lock"]]

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax = axes[0]
    ax.loglog(counts, runtime_ms, marker="o", linewidth=1.3)
    ax.set_xlabel("Particle count $N_p$")
    ax.set_ylabel("Median transition runtime [ms]")
    ax.grid(True, which="both", alpha=0.25)
    slope = summary["runtime_scaling"]["empirical_loglog_runtime_slope"]
    ax.text(0.06, 0.86, f"Empirical slope: {slope:.2f}", transform=ax.transAxes)
    _panel_label(ax, "(a)")

    ax = axes[1]
    x = np.arange(2)
    width = 0.35
    ax.bar(x - width / 2, moved, width, label="Moved fraction")
    ax.bar(x + width / 2, unique, width, label="Unique-parent fraction")
    ax.set_xticks(x, event_labels, rotation=12)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Fraction of particle cloud")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    ax.text(
        1,
        0.16,
        f"max copies = {wrong['max_destination_multiplicity_at_lock']}",
        ha="center",
        fontsize=7.5,
    )
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.4)
    fig.savefig(FIG / "p1fb_literal_mechanism.pdf")
    plt.close(fig)


def fig_p1fc():
    df = pd.read_csv(ROOT / "results" / "phase1" / "p1f_c" / "tuning_report.csv")
    selected = df.iloc[df["rank"].argmin()]
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))

    ax = axes[0]
    ax.scatter(
        df.catastrophic_collapse_fraction,
        df.correct_lock_fraction,
        s=22,
        alpha=0.55,
        label="Sweep settings",
    )
    ax.scatter(
        [selected.catastrophic_collapse_fraction],
        [selected.correct_lock_fraction],
        marker="*",
        s=120,
        label="Selected $(0.5,0,2)$",
    )
    ax.axhline(0.75, linestyle="--", linewidth=0.8, label="Development success floor")
    ax.set_xlabel("Catastrophic-collapse fraction")
    ax.set_ylabel("Correct-lock fraction")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower left")
    _panel_label(ax, "(a)")

    ax = axes[1]
    ax.scatter(
        df.wrong_lock_fraction,
        df.correct_lock_fraction,
        s=22,
        alpha=0.55,
        label="Sweep settings",
    )
    ax.scatter(
        [selected.wrong_lock_fraction],
        [selected.correct_lock_fraction],
        marker="*",
        s=120,
        label="Selected $(0.5,0,2)$",
    )
    ax.axhline(0.75, linestyle="--", linewidth=0.8)
    ax.set_xlabel("Wrong-lock fraction")
    ax.set_ylabel("Correct-lock fraction")
    ax.set_xlim(-0.03, max(0.55, df.wrong_lock_fraction.max() + 0.03))
    ax.set_ylim(-0.03, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower left")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / "p1fc_tuning_tradeoff.pdf")
    plt.close(fig)


def _matched_panel(ax, x, y, xlabel, ylabel, label):
    maxv = max(float(np.max(x)), float(np.max(y))) * 1.05
    ax.plot([0, maxv], [0, maxv], linestyle="--", linewidth=0.9, label="Equal error")
    ax.scatter(x, y, s=26, alpha=0.75)
    ax.set_xlim(0, maxv)
    ax.set_ylim(0, maxv)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, label)


def fig_p1fd():
    df = pd.read_csv(ROOT / "results" / "phase1" / "p1f_d" / "paired_report_1200.csv")
    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.85))
    _matched_panel(
        axes[0],
        df.pf_late_position_rmse_m.to_numpy(),
        df.aacopf_late_position_rmse_m.to_numpy(),
        "PF late position RMSE [m]",
        "AACOPF late position RMSE [m]",
        "(a)",
    )
    _matched_panel(
        axes[1],
        df.pf_late_yaw_rmse_deg.to_numpy(),
        df.aacopf_late_yaw_rmse_deg.to_numpy(),
        "PF late yaw RMSE [deg]",
        "AACOPF late yaw RMSE [deg]",
        "(b)",
    )
    pos_better = np.mean(df.aacopf_late_position_rmse_m < df.pf_late_position_rmse_m)
    yaw_better = np.mean(df.aacopf_late_yaw_rmse_deg < df.pf_late_yaw_rmse_deg)
    note = (
        f"AACOPF lower error in {pos_better:.0%} of position pairs and "
        f"{yaw_better:.0%} of yaw pairs; strict pose success PF 2/20, AACOPF 1/20."
    )
    fig.text(0.5, -0.015, note, ha="center", fontsize=7.5)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / "p1fd_matched_pf_aacopf.pdf")
    plt.close(fig)


def fig_p1fe():
    summary = _load_json("results/phase1/p1f_e/summary.json")
    geometries = ["stationary", "constant_bearing"]
    labels = ["Stationary", "Constant bearing"]
    pf_success = [summary[g]["pf_pose_success_fraction"] for g in geometries]
    aco_success = [summary[g]["aacopf_pose_success_fraction"] for g in geometries]
    pf_pos = [summary[g]["pf_late_position_rmse_mean_m"] for g in geometries]
    aco_pos = [summary[g]["aacopf_late_position_rmse_mean_m"] for g in geometries]

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.75))
    x = np.arange(2)
    width = 0.36
    ax = axes[0]
    ax.bar(x - width / 2, pf_success, width, label="PF")
    ax.bar(x + width / 2, aco_success, width, label="Frozen AACOPF")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 0.32)
    ax.set_ylabel("Strict pose-success fraction")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, "(a)")

    ax = axes[1]
    ax.bar(x - width / 2, pf_pos, width, label="PF")
    ax.bar(x + width / 2, aco_pos, width, label="Frozen AACOPF")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Mean late position RMSE [m]")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    false_fit = summary["stationary"]["failed_pf_late_map_range_fit_rmse_median_m"]
    false_pos = summary["stationary"]["failed_pf_late_position_rmse_mean_m"]
    ax.text(
        0.04,
        0.88,
        f"Stationary failed PF:\nrange-fit RMSE {false_fit:.3f} m\nposition RMSE {false_pos:.1f} m",
        transform=ax.transAxes,
        va="top",
        fontsize=7.4,
    )
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.4)
    fig.savefig(FIG / "p1fe_negative_controls.pdf")
    plt.close(fig)


def fig_p1ff():
    summary = _load_json("results/phase1/p1f_f/summary.json")
    conditions = [
        "clean",
        "sparse_impulses",
        "mixture_moderate",
        "mixture_severe",
        "nlos_burst_moderate",
        "nlos_burst_severe",
    ]
    labels = ["Clean", "Impulses", "Mix. mod.", "Mix. sev.", "NLOS mod.", "NLOS sev."]
    cells = summary["cells"]

    def avg(metric, condition, prefix):
        return np.mean(
            [
                cells["balanced"][condition][f"{prefix}_{metric}"],
                cells["minority_correct"][condition][f"{prefix}_{metric}"],
            ]
        )

    pf_correct = [avg("lock", cond, "pf") for cond in conditions]
    aco_correct = [avg("lock", cond, "aacopf") for cond in conditions]
    pf_wrong = [avg("wrong", cond, "pf") for cond in conditions]
    aco_wrong = [avg("wrong", cond, "aacopf") for cond in conditions]

    fig, axes = plt.subplots(1, 2, figsize=(6.7, 2.85))
    x = np.arange(len(conditions))
    width = 0.36
    ax = axes[0]
    ax.bar(x - width / 2, pf_correct, width, label="PF")
    ax.bar(x + width / 2, aco_correct, width, label="Frozen AACOPF")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Correct-mode lock fraction")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="lower left")
    _panel_label(ax, "(a)")

    ax = axes[1]
    ax.bar(x - width / 2, pf_wrong, width, label="PF")
    ax.bar(x + width / 2, aco_wrong, width, label="Frozen AACOPF")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylim(0, 0.45)
    ax.set_ylabel("Wrong-mode lock fraction")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.4)
    fig.savefig(FIG / "p1ff_ranging_stress.pdf")
    plt.close(fig)


def fig_p1fg():
    summary = _load_json("results/phase1/p1f_g/summary.json")
    controlled = summary["controlled_validation"]
    methods = [
        ("PF", "pf"),
        ("Literal", "literal"),
        ("Cand. only", "candidate_only"),
        ("Move guard", "move_guard_only"),
        ("Dest. guard", "destination_guard_only"),
        ("Guarded", "guarded"),
    ]
    correct = [controlled[key]["correct_lock_fraction"] for _, key in methods]
    wrong = [controlled[key]["wrong_lock_fraction"] for _, key in methods]

    global_res = summary["random_annulus_validation"]
    budgets = np.array([1200, 5000, 10000])
    pf_pos = np.array([global_res[str(n)]["pf"]["late_position_rmse_mean_m"] for n in budgets])
    guard_pos = np.array(
        [global_res[str(n)]["guarded"]["late_position_rmse_mean_m"] for n in budgets]
    )

    runtime = summary["runtime_scaling"]
    fig, axes = plt.subplots(1, 3, figsize=(6.7, 2.55))

    ax = axes[0]
    x = np.arange(len(methods))
    width = 0.36
    ax.bar(x - width / 2, correct, width, label="Correct")
    ax.bar(x + width / 2, wrong, width, label="Wrong")
    ax.set_xticks(x, [label for label, _ in methods], rotation=55, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Lock fraction")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    _panel_label(ax, "(a)")

    ax = axes[1]
    ax.plot(budgets, pf_pos, marker="o", linewidth=1.2, label="PF")
    ax.plot(budgets, guard_pos, marker="s", linestyle="--", linewidth=1.2, label="Guarded")
    ax.set_xscale("log")
    ax.set_xlabel("Particle count $N_p$")
    ax.set_ylabel("Late position RMSE [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    _panel_label(ax, "(b)")

    ax = axes[2]
    slopes = [runtime["literal_loglog_slope"], runtime["bounded_loglog_slope"]]
    ax.bar(["Literal", "Bounded"], slopes)
    ax.set_ylabel("Empirical runtime slope")
    ax.set_ylim(0, 2.2)
    ax.grid(True, axis="y", alpha=0.25)
    ax.text(
        0.5,
        0.12,
        "At $N_p=40000$:\n0.020% of dense pairs",
        transform=ax.transAxes,
        ha="center",
        fontsize=7.4,
    )
    _panel_label(ax, "(c)")
    fig.tight_layout(w_pad=1.0)
    fig.savefig(FIG / "p1fg_scalable_guarded.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_p1a()
    fig_p1b()
    fig_p1c()
    fig_p1d()
    fig_p1e()
    fig_p1fa()
    fig_p1fb()
    fig_p1fc()
    fig_p1fd()
    fig_p1fe()
    fig_p1ff()
    fig_p1fg()
    print("Generated report figures in", FIG)
