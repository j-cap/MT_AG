# ruff: noqa: E402
"""P2E scientific figures from committed compact results."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/phase2/p2e'
OUT = ROOT / 'report/figures/phase2'
OUT.mkdir(parents=True, exist_ok=True)
s = json.loads((DATA / 'summary.json').read_text())
lookup = {(r['condition'], r['mode'], r['budget']): r for r in s['results']}
conditions = [c['name'] for c in s['protocol']['conditions']]
labels = {'reference': 'Reference', 'noise_half': 'White noise x0.5', 'noise_x4': 'White noise x4', 'bias_half': 'Bias drift x0.5', 'bias_x10': 'Bias drift x10', 'combined_poor': 'White x4 + bias x10'}
colors = {'reference': '#333333', 'noise_half': '#3673ad', 'noise_x4': '#3673ad', 'bias_half': '#2d955d', 'bias_x10': '#2d955d', 'combined_poor': '#c76424'}
styles = {'noise_half': '--', 'bias_half': '--'}
budgets = [b for b in s['protocol']['budgets'] if b > 0]
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})

def save(fig, name):
    temp = OUT / f'p2e_{name}.tmp.pdf'
    fig.savefig(temp, bbox_inches='tight')
    if not temp.read_bytes().rstrip().endswith(b'%%EOF'):
        raise RuntimeError(f'Incomplete PDF {temp}')
    temp.replace(OUT / f'p2e_{name}.pdf')
    fig.savefig(OUT / f'p2e_{name}.png', dpi=160, bbox_inches='tight')
    plt.close(fig)
fig, axs = plt.subplots(2, 2, figsize=(10.4, 6.2))
for ax, key, title in zip(axs.flat, ['relative_shape_rmse_m', 'centroid_rmse_m', 'fleet_rmse_m', 'worst_node_rmse_m'], ['(a) Relative shape', '(b) Fleet centroid', '(c) Absolute fleet', '(d) Worst node']):
    for c in conditions:
        ax.plot(budgets, [lookup[c, 'frozen_q', b]['metrics'][key]['mean'] for b in budgets], marker='o', ms=4, ls=styles.get(c, '-'), color=colors[c], label=labels[c])
    if key == 'relative_shape_rmse_m':
        ax.axhline(s['protocol']['shape_target_m'], color='#777777', ls=':', lw=1)
    ax.set_xscale('log')
    ax.set_xticks(budgets, [str(b) for b in budgets])
    ax.minorticks_off()
    ax.set_xlabel('Attempted UWB exchanges (log scale)')
    ax.set_ylabel('Mean whole-run RMSE [m]')
    ax.set_title(title)
    ax.grid(alpha=0.2)
fig.legend(*axs[0, 0].get_legend_handles_labels(), loc='lower center', ncol=3, fontsize=8)
fig.tight_layout(rect=(0, 0.08, 1, 1))
save(fig, 'budget_curves')
comp = {(r['condition'], r['mode'], r['budget'], r['reference']): r for r in s['comparisons']}
fig, axs = plt.subplots(1, 3, figsize=(10.4, 3.5))
for ax, c in zip(axs, s['protocol']['rescaled_q_conditions']):
    vals = [comp[c, 'white_rescaled_q', b, 'frozen_q_same_condition']['metrics']['relative_shape_rmse_m'] for b in budgets]
    mean = np.array([v['mean'] for v in vals])
    lo = np.array([v['ci95'][0] for v in vals])
    hi = np.array([v['ci95'][1] for v in vals])
    ax.errorbar(range(len(budgets)), mean, yerr=[mean - lo, hi - mean], fmt='o', color=colors[c], capsize=3)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_xticks(range(len(budgets)), [str(b) for b in budgets])
    ax.set_xlabel('Attempted exchanges')
    ax.set_title(labels[c])
    ax.grid(axis='y', alpha=0.2)
axs[0].set_ylabel('Rescaled-Q minus frozen-Q shape RMSE [m]')
fig.tight_layout()
save(fig, 'q_adjustment')
profiles = list(csv.DictReader((DATA / 'time_profile.csv').open()))
fig, axs = plt.subplots(2, 2, figsize=(10.4, 6.2), sharex=True)
for row, key in enumerate(['shape_rms_error_m', 'centroid_error_m']):
    for col, budget in enumerate([360, 3600]):
        ax = axs[row, col]
        for c in ['reference', 'noise_x4', 'bias_x10', 'combined_poor']:
            group = [r for r in profiles if r['condition'] == c and r['mode'] == 'frozen_q' and (int(r['budget']) == budget)]
            t = np.array([float(r['time_s']) for r in group])
            m = np.array([float(r['mean_' + key]) for r in group])
            std = np.array([float(r['std_' + key]) for r in group])
            ax.plot(t, m, color=colors[c], label=labels[c])
            ax.fill_between(t, np.maximum(0, m - std), m + std, color=colors[c], alpha=0.1)
        ax.set_title(f"{('Shape' if row == 0 else 'Centroid')}, {budget} exchanges")
        ax.set_ylabel('Instantaneous error [m]')
        ax.grid(alpha=0.2)
        if row == 1:
            ax.set_xlabel('Time [s]')
fig.legend(*axs[0, 0].get_legend_handles_labels(), loc='lower center', ncol=4, fontsize=8)
fig.tight_layout(rect=(0, 0.05, 1, 1))
save(fig, 'time_profiles')
print('Generated three P2E figures')
