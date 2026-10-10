"""Regenerate P2F figures from the audited compact design catalogue."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'report/figures/phase2'
OUT.mkdir(parents=True, exist_ok=True)
ROWS = list(csv.DictReader((ROOT / 'results/phase2/p2f/design_summary.csv').open()))
plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42, 'axes.spines.top': False,
                     'axes.spines.right': False})


def select(**filters):
    return sorted([r for r in ROWS if all(r[k] == str(v) for k, v in filters.items())],
                  key=lambda r: int(r['budget']))


def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.png', dpi=160, bbox_inches='tight')
    plt.close(fig)


fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
base = select(sensor_condition='reference', q_treatment='frozen_q', range_condition='clean', update='naive')
for ax, metric, title in zip(axes, ['relative_shape_rmse_m', 'fleet_rmse_m', 'worst_node_rmse_m'],
                           ['Relative shape', 'Absolute fleet', 'Worst node']):
    for schedule in ['uniform', 'cyclic', 'random', 'geometry', 'information']:
        points = [r for r in base if r['schedule'] in [schedule, 'all_pairs'] and int(r['budget']) > 0]
        x = [int(r['budget']) for r in points]
        y = [float(r[f'mean_{metric}']) for r in points]
        label = 'uniform (all-link rounds)' if schedule == 'uniform' else schedule
        ax.plot(x, y, 'o-', ms=3, label=label)
    front = [r for r in base if r[f'pareto_mean_{metric}'] == 'True' and int(r['budget']) > 0]
    ax.scatter([int(r['budget']) for r in front], [float(r[f'mean_{metric}']) for r in front],
               s=65, facecolors='none', edgecolors='black', linewidths=1, label='empirical frontier')
    ax.set(xscale='log', xlabel='Attempted exchanges / 60 s', ylabel='Mean RMSE (m)', title=title)
    ax.grid(alpha=.2)
axes[0].axhline(.1, color='grey', ls=':', lw=1)
axes[0].legend(fontsize=7)
fig.tight_layout()
save(fig, 'p2f_clean_frontiers')

fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
conditions = [('reference', 'frozen_q', 'reference'), ('noise_half', 'frozen_q', 'half white noise'),
              ('noise_half', 'white_rescaled_q', 'half noise, rescaled Q'),
              ('bias_half', 'frozen_q', 'half bias drift'), ('noise_x4', 'frozen_q', '4x white noise'),
              ('bias_x10', 'frozen_q', '10x bias drift'), ('combined_poor', 'frozen_q', 'both degraded')]
for sensor, q, label in conditions:
    points = [r for r in select(sensor_condition=sensor, q_treatment=q, range_condition='clean', update='naive')
              if r['schedule'] in ['cyclic', 'all_pairs']]
    x = [int(r['budget']) for r in points]
    for ax, stat in zip(axes, ['mean', 'p95']):
        ax.plot(x, [float(r[f'{stat}_relative_shape_rmse_m']) for r in points], 'o-', ms=3, label=label)
for ax, title in zip(axes, ['Mean across 20 seeds', 'Empirical 95th percentile across 20 seeds']):
    ax.axhline(.1, color='black', ls=':', label='10 cm target')
    ax.set(xscale='log', xlabel='Attempted exchanges / 60 s', ylabel='Shape RMSE (m)', title=title)
    ax.grid(alpha=.2)
axes[1].legend(fontsize=7)
fig.tight_layout()
save(fig, 'p2f_sensor_targets')

fig, ax = plt.subplots(figsize=(10, 4))
cases = ['clean', 'loss30', 'loss50', 'network_outage', 'node_outage', 'outliers', 'nlos']
x = np.arange(len(cases))
for offset, schedule, update, label in [(-.24, 'cyclic', 'naive', 'cyclic ordinary'),
                                      (-.08, 'information', 'naive', 'information ordinary'),
                                      (.08, 'cyclic', 'gated', 'cyclic gated'),
                                      (.24, 'information', 'gated', 'information gated')]:
    xx, yy, err = [], [], []
    for i, case in enumerate(cases):
        points = select(sensor_condition='reference', q_treatment='frozen_q', range_condition=case,
                        update=update, schedule=schedule, budget=360)
        if points:
            xx.append(i + offset)
            yy.append(float(points[0]['mean_relative_shape_rmse_m']))
            err.append(float(points[0]['std_relative_shape_rmse_m']))
    ax.errorbar(xx, yy, yerr=err, fmt='o', capsize=3, ms=4, label=label)
ax.axhline(.1, color='grey', ls=':')
ax.set(xticks=x, xticklabels=['clean', '30% loss', '50% loss', 'network\noutage', 'node\noutage', 'outliers', 'NLOS bias'],
       ylabel='Shape RMSE (m)', title='Reliability: 360 attempted exchanges only; mean ± sample SD')
ax.legend(fontsize=8, ncol=2)
ax.grid(axis='y', alpha=.2)
fig.tight_layout()
save(fig, 'p2f_reliability_slice')
