# ruff: noqa: E402
"""P2E: sensor sensitivity and a separately labeled process-noise adjustment."""
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
from experiments.run_phase2b import _time_errors, _write_csv
from mt_ag.fleet_simulation import generate_four_node_fleet
from mt_ag.imu import simulate_imu_measurements
from mt_ag.joint_ekf import run_joint_ekf
from mt_ag.link_selection import baseline_mask, exchange_capacity

def rescale_imu(ideal, measured, bias, white_scale, bias_scale):
    """Scale the same latent white-noise and random-walk realizations.

    Baseline remains byte-identical; improved sensors retain initialization error.
    Bias starts at zero in P2A; scaling does not introduce an initial offset.
    """
    if white_scale < 0 or bias_scale < 0:
        raise ValueError('IMU amplitude scales must be nonnegative')
    if white_scale == bias_scale == 1:
        return measured.copy()
    noise = measured - ideal - bias
    return ideal + bias_scale * bias + white_scale * noise

def process_settings(base_ekf, mode, white_scale):
    if mode not in ('frozen_q', 'white_rescaled_q'):
        raise ValueError(mode)
    factor = white_scale if mode == 'white_rescaled_q' else 1.0
    return {k: float(v) * factor for k, v in base_ekf.items()}

def run_seed(seed, cfg, base):
    sim, ic, uc = (base['simulation'], base['imu'], base['uwb'])
    init, ec, campaign = (base['initialization'], base['ekf'], base['campaign'])
    fleet = generate_four_node_fleet(dt=sim['dt_s'], duration=sim['duration_s'])
    if fleet.n_nodes != 4 or len(fleet.t) != 6001:
        raise ValueError('P2E requires frozen four-node trajectories')
    bs = campaign['seed_offset']
    measured = np.zeros_like(fleet.ideal_imu)
    bias = np.zeros_like(measured)
    for node in range(4):
        data = simulate_imu_measurements(fleet.ideal_imu[:, node], fleet.dt, np.random.default_rng(bs + 10000 * seed + 100 * node + 1), **{k: v for k, v in ic.items() if k.startswith('sigma_')})
        measured[:, node], bias[:, node] = (data.measured, data.bias)
    ranges = _range_data(fleet.state, uc['sigma_range_m'], np.random.default_rng(bs + 10000 * seed + 9001))
    initial = fleet.state[0].copy()
    initial[:, 2:4] += np.random.default_rng(bs + 10000 * seed + 8001).normal(0, init['velocity_error_std_mps'], (4, 2))
    p0 = _initial_covariance(init, 4)
    masks = {b: baseline_mask(exchange_capacity(len(fleet.t), 10, b, 6), 4, 'cyclic') for b in cfg['budgets']}
    sample = np.arange(0, len(fleet.t), round(1 / (cfg['profile_sample_rate_hz'] * fleet.dt)))
    rows, profiles = ([], [])
    for condition in cfg['conditions']:
        name, w, b = (condition['name'], condition['white_scale'], condition['bias_scale'])
        imu = rescale_imu(fleet.ideal_imu, measured, bias, w, b)
        modes = ['frozen_q'] + (['white_rescaled_q'] if name in cfg['rescaled_q_conditions'] else [])
        zero_history = None
        for mode in modes:
            for budget in cfg['budgets']:
                result = run_joint_ekf(imu_measurements=imu, pairwise_ranges=ranges, range_mask=masks[budget], initial_state=initial, initial_covariance=p0, dt=fleet.dt, sigma_range_m=uc['sigma_range_m'], **process_settings(ec, mode, w))
                if (result.range_attempt_count, result.range_received_count, result.range_update_count) != (budget, budget, budget):
                    raise RuntimeError('P2E clean-range budget accounting failed')
                if budget == 0:
                    if zero_history is not None and (not np.array_equal(result.state, zero_history)):
                        raise RuntimeError('Q adjustment changed zero-UWB state propagation')
                    zero_history = result.state.copy()
                m = _fleet_metrics(result.state, fleet.state)
                if not np.isclose(m['fleet_position_rmse_m'] ** 2, m['centroid_position_rmse_m'] ** 2 + m['relative_shape_rmse_m'] ** 2):
                    raise RuntimeError('Whole-run error decomposition failed')
                row = dict(seed=seed, condition=name, mode=mode, budget=budget, communication_ratio=budget / 3600, white_scale=w, bias_scale=b, attempted=result.range_attempt_count, received=result.range_received_count, accepted=result.range_update_count, fleet_rmse_m=m['fleet_position_rmse_m'], worst_node_rmse_m=m['worst_node_position_rmse_m'], centroid_rmse_m=m['centroid_position_rmse_m'], relative_shape_rmse_m=m['relative_shape_rmse_m'], yaw_rmse_deg=m['fleet_yaw_rmse_deg'])
                errors = _time_errors(result.state, fleet.state)
                for node, value in enumerate(m['node_position_rmse_m'], 1):
                    row[f'node{node}_rmse_m'] = value
                for key, values in errors.items():
                    row['late_' + key.replace('error', 'rmse')] = float(np.sqrt(np.mean(values[4000:] ** 2)))
                rows.append(row)
                profiles.append(dict(condition=name, mode=mode, budget=budget, **{k: v[sample].tolist() for k, v in errors.items()}))
    for path, budgets in [('p2c', [180, 360, 720]), ('p2b', [0, 3600])]:
        with (ROOT / f'results/phase2/{path}/seed_metrics.csv').open() as f:
            refs = list(csv.DictReader(f))
        for budget in budgets:
            ref = next((r for r in refs if int(r['seed']) == seed and (r.get('policy') == 'cyclic' and int(r['budget']) == budget if path == 'p2c' else int(r['uwb_exchanges']) == budget)))
            row = next((r for r in rows if r['condition'] == 'reference' and r['budget'] == budget))
            for key in ['fleet_rmse_m', 'worst_node_rmse_m', 'centroid_rmse_m', 'relative_shape_rmse_m', 'yaw_rmse_deg']:
                if not np.isclose(row[key], float(ref[key]), atol=1e-10, rtol=1e-10):
                    raise RuntimeError(f'Inherited endpoint mismatch: {seed}, {budget}, {key}')
    return dict(rows=rows, profiles=profiles, times=fleet.t[sample].tolist())

def aggregate(runs, cfg, base, digest, out):
    rows = [r for run in runs for r in run['rows']]
    _write_csv(out / 'seed_metrics.csv', rows)
    groups = sorted({(r['condition'], r['mode'], r['budget']) for r in rows})
    summary = dict(protocol=cfg, base_config_snapshot=base, n_seeds=len(runs), n_runs=len(rows), implementation_fingerprint=digest, numpy_version=np.__version__, results=[], comparisons=[], target_budgets=[])
    metric_keys = ['fleet_rmse_m', 'worst_node_rmse_m', 'centroid_rmse_m', 'relative_shape_rmse_m', 'yaw_rmse_deg', 'late_shape_rms_rmse_m', 'late_fleet_rms_rmse_m', 'late_centroid_rmse_m']
    profile_rows = []
    lookup = {(r['seed'], r['condition'], r['mode'], r['budget']): r for r in rows}
    rng = np.random.default_rng(cfg['bootstrap_seed'])
    boot = rng.integers(0, len(runs), (cfg['bootstrap_resamples'], len(runs)))
    for condition, mode, budget in groups:
        group = [lookup[s, condition, mode, budget] for s in range(len(runs))]
        metrics = {}
        for key in metric_keys:
            vals = np.array([r[key] for r in group])
            metrics[key] = dict(mean=float(vals.mean()), std=float(vals.std(ddof=1)), p95=float(np.quantile(vals, 0.95)))
        summary['results'].append(dict(condition=condition, mode=mode, budget=budget, metrics=metrics, seeds_meeting_shape_target=int(sum((r['relative_shape_rmse_m'] <= cfg['shape_target_m'] for r in group)))))
        histories = [next((p for p in run['profiles'] if (p['condition'], p['mode'], p['budget']) == (condition, mode, budget))) for run in runs]
        for index, t in enumerate(runs[0]['times']):
            record = dict(condition=condition, mode=mode, budget=budget, time_s=t)
            for key in ('fleet_rms_error_m', 'centroid_error_m', 'shape_rms_error_m'):
                vals = [h[key][index] for h in histories]
                record[f'mean_{key}'] = float(np.mean(vals))
                record[f'std_{key}'] = float(np.std(vals, ddof=1))
            profile_rows.append(record)
        references = [('reference_same_budget', 'reference', 'frozen_q', budget), ('full_same_condition', condition, mode, 3600)]
        if mode == 'white_rescaled_q':
            references.append(('frozen_q_same_condition', condition, 'frozen_q', budget))
        for label, rc, rm, rb in references:
            item = dict(condition=condition, mode=mode, budget=budget, reference=label, metrics={})
            for key in ['relative_shape_rmse_m', 'fleet_rmse_m', 'centroid_rmse_m', 'worst_node_rmse_m']:
                diff = np.array([lookup[s, condition, mode, budget][key] - lookup[s, rc, rm, rb][key] for s in range(len(runs))])
                item['metrics'][key] = dict(mean=float(diff.mean()), std=float(diff.std(ddof=1)), wins=int(np.sum(diff < 0)), ci95=np.quantile(diff[boot].mean(axis=1), [0.025, 0.975]).tolist())
            summary['comparisons'].append(item)
    for condition, mode in sorted({(c, m) for c, m, b in groups}):
        group = [r for r in summary['results'] if (r['condition'], r['mode']) == (condition, mode)]
        meet = [r['budget'] for r in group if r['metrics']['relative_shape_rmse_m']['mean'] <= cfg['shape_target_m']]
        tail_meet = [r['budget'] for r in group if r['metrics']['relative_shape_rmse_m']['p95'] <= cfg['shape_target_m']]
        summary['target_budgets'].append(dict(condition=condition, mode=mode, smallest_tested_mean_target_budget=min(meet) if meet else None, smallest_tested_p95_target_budget=min(tail_meet) if tail_meet else None))
    _write_csv(out / 'time_profile.csv', profile_rows, significant_digits=7)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    lines = ['# P2E IMU sensitivity', '', '| Condition | Q treatment | Exchanges | Shape RMSE [m], mean +/- SD | Fleet mean [m] | Centroid mean [m] |', '|---|---|---:|---:|---:|---:|']
    for r in summary['results']:
        m = r['metrics']
        v = m['relative_shape_rmse_m']
        lines.append(f"| {r['condition']} | {r['mode']} | {r['budget']} | {v['mean']:.4f} +/- {v['std']:.4f} | {m['fleet_rmse_m']['mean']:.4f} | {m['centroid_rmse_m']['mean']:.4f} |")
    (out / 'summary.md').write_text('\n'.join(lines) + '\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    cfg = yaml.safe_load((ROOT / 'configs/phase2e.yaml').read_text())
    base = yaml.safe_load((ROOT / cfg['base_config']).read_text())
    fingerprint = hashlib.sha256(json.dumps([cfg, base], sort_keys=True).encode())
    for path in ['experiments/run_phase2e.py', 'src/mt_ag/imu.py', 'src/mt_ag/joint_ekf.py', 'src/mt_ag/link_selection.py', 'src/mt_ag/fleet_simulation.py', 'experiments/run_phase2a.py', 'experiments/run_phase2b.py']:
        fingerprint.update((ROOT / path).read_bytes())
    digest = fingerprint.hexdigest()
    out = ROOT / 'results/phase2/p2e'
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
            temp = checkpoint / f'{seed}.tmp'
            temp.write_text(json.dumps(dict(fingerprint=digest, run=run)))
            os.replace(temp, checkpoint / f'{seed}.json')
            runs[seed] = run
            print(f'P2E seed {seed + 1} completed; {len(runs)}/{n} saved', flush=True)
    aggregate([runs[s] for s in range(n)], cfg, base, digest, out)
    print(f"P2E complete: {sum((len(run['rows']) for run in runs.values()))} runs", flush=True)
if __name__ == '__main__':
    main()
