from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "report" / "figures" / "phase1"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
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
})


def _panel_label(ax, label):
    ax.text(0.01, 0.99, label, transform=ax.transAxes, va="top", ha="left", fontweight="bold")


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


def fig_p1c():
    summary = json.loads((ROOT / "results" / "phase1" / "p1c" / "summary.json").read_text())
    grid = summary["grid_resolution"]
    g = pd.DataFrame([
        {
            "n_particles": v["n_particles"],
            "success": v["pose_success_fraction"],
            "rmse": v["late_position_rmse_mean_m"],
        }
        for _, v in sorted(grid.items(), key=lambda kv: kv[1]["n_particles"])
    ])
    res = summary["resampling_sensitivity"]
    r = pd.DataFrame([
        {
            "gamma": float(k),
            "success": v["pose_success_fraction"],
            "rmse": v["late_position_rmse_mean_m"],
        }
        for k, v in sorted(res.items(), key=lambda kv: float(kv[0]))
    ])

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
    ax.legend(lines, [x.get_label() for x in lines], loc="upper right")
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
    ax.legend(lines, [x.get_label() for x in lines], loc="lower right")
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
        d = df[df.geometry == geom]
        ax.semilogy(
            d.time_s,
            np.maximum(d.obs_sigma_ratio, 1e-18),
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
        d = df[df.geometry == geom]
        ax.plot(d.time_s, d.position_error_mean_m, linestyle=ls, linewidth=1.4, label=labels[geom])
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Mean position error [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left")
    _panel_label(ax, "(b)")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(FIG / "p1d_geometry_observability.pdf")
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


if __name__ == "__main__":
    fig_p1a()
    fig_p1c()
    fig_p1d()
    fig_p1fc()
    fig_p1fd()
    print("Generated report figures in", FIG)
