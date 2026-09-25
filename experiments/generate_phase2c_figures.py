"""Regenerate P2C report figures from committed compact campaign data."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phase2/p2c"
OUT = ROOT / "report/figures/phase2"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "legend.fontsize": 7, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight",
})
STYLES = {
    "cyclic": ("Cyclic", "#555555", "o"),
    "random": ("Random", "#9467bd", "s"),
    "geometry": ("Geometry", "#1f77b4", "^"),
    "information": ("Information", "#2ca02c", "D"),
}


def _rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _finish(fig: plt.Figure, filename: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / filename)
    plt.close(fig)


def budget_tradeoff(summary: dict) -> None:
    budgets = summary["protocol"]["budgets"]
    entries = {(row["budget"], row["policy"]): row for row in summary["results"]}
    with (ROOT / "results/phase2/p2b/summary.json").open() as handle:
        reference = json.load(handle)
    by_count = {row["uwb_exchanges"]: row for row in reference["rates"]}
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.65))
    columns = [("fleet_rmse_m", "Fleet RMSE [m]"),
               ("worst_node_rmse_m", "Worst-node RMSE [m]"),
               ("relative_shape_rmse_m", "Shape RMSE [m]")]
    for ax, (metric, ylabel) in zip(axes, columns):
        for policy, (label, color, marker) in STYLES.items():
            group = [entries[budget, policy]["metrics"][metric] for budget in budgets]
            means = [row["mean"] for row in group]
            ax.plot(budgets, means, marker=marker, color=color, linewidth=1.1,
                    markersize=3.5, label=label)
        ax.plot(budgets, [by_count[b]["metrics"][metric]["mean"] for b in budgets],
                linestyle="--", marker="x", color="#b95030", linewidth=.9,
                markersize=4, label="P2B all-link")
        ax.set_xticks(budgets)
        ax.set_xlabel("UWB exchanges per run")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=.25)
    axes[0].legend(loc="best", ncol=1)
    for index, ax in enumerate(axes):
        ax.set_title(f"({chr(97 + index)})")
    _finish(fig, "p2c_budget_tradeoff.pdf")


def paired_differences(seed_rows: list[dict]) -> None:
    budgets = [180, 360, 720]
    policies = ["random", "geometry", "information"]
    lookup = {(int(row["seed"]), int(row["budget"]), row["policy"]): row
              for row in seed_rows}
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.65))
    keys = [("fleet_rmse_m", "Fleet RMSE difference [m]"),
            ("worst_node_rmse_m", "Worst-node difference [m]"),
            ("relative_shape_rmse_m", "Shape RMSE difference [m]")]
    for ax, (key, ylabel) in zip(axes, keys):
        ax.axhline(0, linestyle="--", color="black", linewidth=.8)
        for offset, policy in zip([-.16, 0, .16], policies):
            label, color, marker = STYLES[policy]
            differences = [np.array([
                float(lookup[seed, b, policy][key]) -
                float(lookup[seed, b, "cyclic"][key]) for seed in range(20)
            ]) for b in budgets]
            means = [np.mean(value) for value in differences]
            sds = [np.std(value, ddof=1) for value in differences]
            x = np.arange(len(budgets)) + offset
            ax.errorbar(x, means, yerr=sds, fmt=marker, markersize=3.8,
                        capsize=2.2, color=color, label=label)
        ax.set_xticks(range(len(budgets)), [str(b) for b in budgets])
        ax.set_xlabel("UWB exchanges per run")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=.22)
    axes[-1].legend(loc="best")
    for index, ax in enumerate(axes):
        ax.set_title(f"({chr(97 + index)})")
    _finish(fig, "p2c_paired_differences.pdf")


def link_usage(rows: list[dict]) -> None:
    pairs = ["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"]
    policies = list(STYLES)
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.55), sharey=True)
    for ax, budget in zip(axes, [180, 360, 720]):
        matrix = np.array([[np.mean([
            int(row["exchanges"]) / budget for row in rows
            if int(row["budget"]) == budget and row["policy"] == policy
            and row["pair"] == pair
        ]) for pair in pairs] for policy in policies])
        view = ax.imshow(matrix, aspect="auto", vmin=.10, vmax=.23, cmap="YlGnBu")
        ax.set_title(f"{budget} exchanges")
        ax.set_xticks(range(6), pairs, rotation=45, ha="right")
        ax.set_yticks(range(4), [STYLES[p][0] for p in policies])
        ax.set_xlabel("Pair")
    fig.colorbar(view, ax=axes, fraction=.025, pad=.035,
                 ticks=[.10, 1 / 6, .23],
                 label="Mean fraction of exchanges")
    fig.subplots_adjust(bottom=.22, right=.89, wspace=.14)
    fig.savefig(OUT / "p2c_link_usage.pdf", bbox_inches="tight")
    plt.close(fig)


def time_profiles(rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(6.8, 4.1), sharex=True)
    for ax, budget in zip(axes, [180, 720]):
        for policy, (label, color, _) in STYLES.items():
            subset = [row for row in rows if int(row["budget"]) == budget
                      and row["policy"] == policy]
            times = np.array([float(row["time_s"]) for row in subset])
            mean = np.array([float(row["mean_shape_rms_error_m"]) for row in subset])
            std = np.array([float(row["std_shape_rms_error_m"]) for row in subset])
            ax.plot(times, mean, color=color, linewidth=1.05, label=label)
            ax.fill_between(times, np.maximum(0, mean - std), mean + std,
                            color=color, alpha=.1, linewidth=0)
        ax.set_title(f"{budget} pairwise exchanges over 60 s")
        ax.set_ylabel("Shape RMS error [m]")
        ax.grid(True, alpha=.25)
    axes[0].legend(loc="upper right", ncol=4)
    axes[-1].set_xlabel("Time [s]")
    _finish(fig, "p2c_time_profiles.pdf")


def main() -> None:
    summary = json.loads((DATA / "summary.json").read_text())
    budget_tradeoff(summary)
    paired_differences(_rows(DATA / "seed_metrics.csv"))
    link_usage(_rows(DATA / "link_counts.csv"))
    time_profiles(_rows(DATA / "time_profile.csv"))


if __name__ == "__main__":
    main()
