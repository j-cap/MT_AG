"""P2D: reception and range-quality stresses at a fixed attempted budget.

Run with PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python experiments/run_phase2d.py.
Seed checkpoints are fingerprinted against config and core implementation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
from experiments.run_phase2a import _fleet_metrics, _initial_covariance, _range_data
from experiments.run_phase2b import _time_errors, _write_csv, uniform_schedule
from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet
from mt_ag.imu import simulate_imu_measurements
from mt_ag.joint_ekf import run_joint_ekf
from mt_ag.link_selection import PredictiveLinkSelector, baseline_mask, eligible_mask, exchange_capacity


def failure_fields(t, n_nodes, seed, cfg, base_seed):
    """Draw potential failures for ALL links/times, independent of selection."""
    rng = np.random.default_rng(base_seed + 10000 * seed + 9501)
    shape = (len(t), n_nodes, n_nodes)
    loss = np.zeros(shape)
    bias = np.zeros(shape)
    for i, j in all_pairs(n_nodes):
        loss[:, i, j] = rng.random(len(t))
        loss[:, j, i] = loss[:, i, j]
    # Independent streams prevent the loss draw count changing outlier realizations.
    rng = np.random.default_rng(base_seed + 10000 * seed + 9601)
    for i, j in all_pairs(n_nodes):
        active = rng.random(len(t)) < cfg['outlier_probability']
        bias[:, i, j] = active * rng.uniform(cfg['outlier_bias_min_m'], cfg['outlier_bias_max_m'], len(t))
        bias[:, j, i] = bias[:, i, j]
    return loss, bias


def stress(case, t, loss, outliers, cfg):
    received = np.ones(loss.shape, dtype=bool)
    bias = np.zeros_like(loss)
    window = (t >= cfg['outage_start_s']) & (t < cfg['outage_end_s'])
    if case.startswith('loss'):
        received = loss >= int(case[4:]) / 100
    elif case == 'network_outage':
        received[window] = False
    elif case == 'node_outage':
        received[window, 0, :] = False
        received[window, :, 0] = False
    elif case == 'outliers':
        bias = outliers
    elif case == 'nlos':
        active = (t >= cfg['nlos_start_s']) & (t < cfg['nlos_end_s'])
        i, j = cfg['nlos_pair']
        bias[active, i, j] = bias[active, j, i] = cfg['nlos_bias_m']
    elif case != 'clean':
        raise ValueError(case)
    return received, bias


def run_seed(seed, cfg, base):
    sim, ic, uc = base['simulation'], base['imu'], base['uwb']
    init, ec, campaign = base['initialization'], base['ekf'], base['campaign']
    fleet = generate_four_node_fleet(dt=sim['dt_s'], duration=sim['duration_s'])
    if fleet.n_nodes != 4 or len(fleet.t) != 6001:
        raise ValueError('P2D requires the frozen P2A benchmark')
    bs = campaign['seed_offset']
    measured = np.zeros_like(fleet.ideal_imu)
    for node in range(4):
        measured[:, node] = simulate_imu_measurements(
            fleet.ideal_imu[:, node], fleet.dt,
            np.random.default_rng(bs + 10000 * seed + 100 * node + 1),
            **{key: value for key, value in ic.items() if key.startswith('sigma_')}).measured
    ranges = _range_data(fleet.state, uc['sigma_range_m'], np.random.default_rng(bs + 10000 * seed + 9001))
    initial = fleet.state[0].copy()
    initial[:, 2:4] += np.random.default_rng(bs + 10000 * seed + 8001).normal(
        0, init['velocity_error_std_mps'], (4, 2))
    args = dict(imu_measurements=measured, initial_state=initial,
                initial_covariance=_initial_covariance(init, 4), dt=fleet.dt,
                sigma_range_m=uc['sigma_range_m'], **ec)
    capacity = exchange_capacity(len(fleet.t), round(1 / (uc['max_rate_hz'] * fleet.dt)), cfg['budget'], 6)
    masks = {'cyclic': baseline_mask(capacity, 4, 'cyclic'),
             'information': eligible_mask(capacity, 4),
             'p2b_uniform': uniform_schedule(len(fleet.t), 4, fleet.dt, 1.0)}
    fields = failure_fields(fleet.t, 4, seed, cfg, bs)
    sample = np.arange(0, len(fleet.t), round(1 / (cfg['profile_sample_rate_hz'] * fleet.dt)))
    rows, profiles, links = [], [], []
    for case in cfg['cases']:
        reception, bias = stress(case, fleet.t, *fields, cfg)
        for mode in (['naive', 'gated'] if case in cfg['gated_cases'] else ['naive']):
            for policy in cfg['policies']:
                selector = (PredictiveLinkSelector('information', capacity, uc['sigma_range_m'], fleet.dt, 4)
                            if policy == 'information' else None)
                result = run_joint_ekf(**args, pairwise_ranges=ranges + bias, range_mask=masks[policy],
                                       reception_mask=reception, range_selector=selector,
                                       innovation_gate_nis=cfg['innovation_gate_nis'] if mode == 'gated' else None)
                if result.range_attempt_count != cfg['budget']:
                    raise RuntimeError('Attempt budget violated')
                if result.range_received_count != result.range_update_count + result.range_rejected_count:
                    raise RuntimeError('Reception accounting violated')
                m = _fleet_metrics(result.state, fleet.state)
                row = dict(seed=seed, case=case, mode=mode, policy=policy,
                           attempted=result.range_attempt_count, received=result.range_received_count,
                           accepted=result.range_update_count, rejected=result.range_rejected_count,
                           fleet_rmse_m=m['fleet_position_rmse_m'],
                           worst_node_rmse_m=m['worst_node_position_rmse_m'],
                           centroid_rmse_m=m['centroid_position_rmse_m'],
                           relative_shape_rmse_m=m['relative_shape_rmse_m'],
                           yaw_rmse_deg=m['fleet_yaw_rmse_deg'])
                errors = _time_errors(result.state, fleet.state)
                shape = errors['shape_rms_error_m']
                for name, lo, hi in [('pre', 10, 20), ('outage', 20, 30), ('recovery', 30, 40), ('late', 40, 60.001)]:
                    window = (fleet.t >= lo) & (fleet.t < hi)
                    row[f'{name}_shape_rmse_m'] = float(np.sqrt(np.mean(shape[window]**2)))
                row['shape_at_30s_m'] = float(shape[3000])
                row['shape_at_40s_m'] = float(shape[4000])
                rows.append(row)
                profiles.append(dict(case=case, mode=mode, policy=policy,
                                     **{key: values[sample].tolist() for key, values in errors.items()}))
                selections = ([(t, i, j) for t, i, j, _ in selector.records] if selector else
                              [(int(t), i, j) for t in np.flatnonzero(capacity if policy == 'cyclic' else np.any(masks[policy], axis=(1, 2)))
                               for i, j in all_pairs(4) if masks[policy][t, i, j]])
                for i, j in all_pairs(4):
                    chosen = [(t, a, b) for t, a, b in selections if (a, b) == (i, j)]
                    links.append(dict(seed=seed, case=case, mode=mode, policy=policy, pair=f'{i+1}-{j+1}',
                                      attempted=len(chosen), received=sum(int(reception[t, i, j]) for t, _, _ in chosen),
                                      corrupted_received=sum(int(reception[t, i, j] and bias[t, i, j] != 0) for t, _, _ in chosen)))
    # Clean runs must reproduce committed P2C and P2B endpoints for every seed.
    for policy, path in [('cyclic', 'p2c'), ('information', 'p2c'), ('p2b_uniform', 'p2b')]:
        with (ROOT / f'results/phase2/{path}/seed_metrics.csv').open() as f:
            ref = next(r for r in csv.DictReader(f) if int(r['seed']) == seed and
                       (r.get('policy') == policy and int(r['budget']) == cfg['budget'] if path == 'p2c'
                        else int(r['uwb_exchanges']) == cfg['budget']))
        row = next(r for r in rows if r['policy'] == policy and r['case'] == 'clean' and r['mode'] == 'naive')
        for key in ['fleet_rmse_m', 'worst_node_rmse_m', 'centroid_rmse_m', 'relative_shape_rmse_m']:
            if not np.isclose(row[key], float(ref[key]), atol=1e-10, rtol=1e-10):
                raise RuntimeError(f'Clean endpoint mismatch: {seed} {policy} {key}')
    return dict(rows=rows, profiles=profiles, links=links, times=fleet.t[sample].tolist())


def aggregate(runs, cfg, out):
    rows = [row for run in runs for row in run['rows']]
    _write_csv(out / 'seed_metrics.csv', rows)
    _write_csv(out / 'link_counts.csv', [row for run in runs for row in run['links']])
    groups = sorted({(r['case'], r['mode'], r['policy']) for r in rows})
    summary = dict(protocol=cfg, n_seeds=len(runs), results=[], comparisons=[])
    profile_rows = []
    keys = [k for k in rows[0] if k not in ('seed', 'case', 'mode', 'policy')]
    for case, mode, policy in groups:
        group = [r for r in rows if (r['case'], r['mode'], r['policy']) == (case, mode, policy)]
        summary['results'].append(dict(case=case, mode=mode, policy=policy, metrics={
            k: dict(mean=float(np.mean([r[k] for r in group])), std=float(np.std([r[k] for r in group], ddof=1))) for k in keys}))
        histories = [next(p for p in run['profiles'] if (p['case'], p['mode'], p['policy']) == (case, mode, policy)) for run in runs]
        for idx, t in enumerate(runs[0]['times']):
            record = dict(case=case, mode=mode, policy=policy, time_s=t)
            for key in ('fleet_rms_error_m', 'centroid_error_m', 'shape_rms_error_m'):
                values = [h[key][idx] for h in histories]
                record[f'mean_{key}'] = float(np.mean(values))
                record[f'std_{key}'] = float(np.std(values, ddof=1))
            profile_rows.append(record)
    rng = np.random.default_rng(cfg['bootstrap_seed'])
    boot = rng.integers(0, len(runs), (cfg['bootstrap_resamples'], len(runs)))
    lookup = {(r['seed'], r['case'], r['mode'], r['policy']): r for r in rows}
    for case, mode, policy in groups:
        references = [('same_case_cyclic', case, mode, 'cyclic'), ('clean_same_policy', 'clean', mode, policy)]
        if mode == 'gated':
            references.append(('same_case_naive', case, 'naive', policy))
        for name, rc, rm, rp in references:
            if not all((s, rc, rm, rp) in lookup for s in range(len(runs))):
                continue
            item = dict(case=case, mode=mode, policy=policy, reference=name, metrics={})
            for key in ('relative_shape_rmse_m', 'fleet_rmse_m', 'worst_node_rmse_m', 'recovery_shape_rmse_m', 'late_shape_rmse_m'):
                diff = np.array([lookup[s, case, mode, policy][key] - lookup[s, rc, rm, rp][key] for s in range(len(runs))])
                ci = np.quantile(diff[boot].mean(axis=1), [.025, .975])
                item['metrics'][key] = dict(mean=float(diff.mean()), std=float(diff.std(ddof=1)),
                                           wins=int(np.sum(diff < 0)), ci95=ci.tolist())
            summary['comparisons'].append(item)
    _write_csv(out / 'time_profile.csv', profile_rows, significant_digits=7)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    cfg = yaml.safe_load((ROOT / 'configs/phase2d.yaml').read_text())
    base = yaml.safe_load((ROOT / cfg['base_config']).read_text())
    fingerprint = hashlib.sha256(json.dumps([cfg, base], sort_keys=True).encode())
    for path in ['experiments/run_phase2d.py', 'src/mt_ag/joint_ekf.py', 'src/mt_ag/link_selection.py']:
        fingerprint.update((ROOT / path).read_bytes())
    digest = fingerprint.hexdigest()
    out = ROOT / 'results/phase2/p2d'
    checkpoint = out / '_checkpoints'
    checkpoint.mkdir(parents=True, exist_ok=True)
    runs = {}
    n = base['campaign']['n_seeds']
    for seed in range(n):
        path = checkpoint / f'{seed}.json'
        if path.exists():
            saved = json.loads(path.read_text())
            if saved['fingerprint'] == digest:
                runs[seed] = saved['run']
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_seed, s, cfg, base): s for s in range(n) if s not in runs}
        for future in as_completed(futures):
            seed = futures[future]
            run = future.result()
            tmp = checkpoint / f'{seed}.tmp'
            tmp.write_text(json.dumps(dict(fingerprint=digest, run=run)))
            os.replace(tmp, checkpoint / f'{seed}.json')
            runs[seed] = run
            print(f'P2D completed seed {seed+1}; {len(runs)}/{n} saved', flush=True)
    aggregate([runs[s] for s in range(n)], cfg, out)
    print(f'P2D complete: {sum(len(r["rows"]) for r in runs.values())} runs', flush=True)

if __name__ == '__main__':
    main()
