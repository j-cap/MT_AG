# ruff: noqa: E402
"""P2G: matched continuing global-position information, three-node scoring."""
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

METRICS = ['fleet_rmse_m', 'worst_node_rmse_m', 'centroid_rmse_m',
           'relative_shape_rmse_m', 'yaw_rmse_deg', 'late_fleet_rmse_m', 'reference_rmse_m']


def reference_cases(fleet, cfg, seed, offset):
    node = cfg['reference_node_zero_based']
    # Independent stream; identical noisy reports are shared by both noisy cases.
    sigma = cfg['reference_noise_sigma_m']
    measured = fleet.state[:, node, :2] + np.random.default_rng(
        offset + 10000 * seed + 9501).normal(0, sigma, (len(fleet.t), 2))
    stride = round(1 / (cfg['reference_rate_hz'] * fleet.dt))
    mask = np.zeros(len(fleet.t), bool)
    mask[stride::stride] = True
    exact = np.ones(len(fleet.t), bool)
    exact[0] = False
    start, end = cfg['reference_outage_s']
    intermittent = mask & ~((fleet.t >= start) & (fleet.t < end))
    return {'unanchored': (None, None, 0),
            'exact_100hz': (fleet.state[:, node, :2], exact, 0),
            'noisy_1hz': (measured, mask, sigma),
            'noisy_1hz_outage': (measured, intermittent, sigma)}


def run_seed(seed, cfg, base):
    fleet = generate_four_node_fleet(dt=base['simulation']['dt_s'], duration=base['simulation']['duration_s'])
    if fleet.n_nodes != 4 or len(fleet.t) != 6001:
        raise ValueError('P2G requires frozen P2A four-node benchmark')
    offset = base['campaign']['seed_offset']
    measured = np.zeros_like(fleet.ideal_imu)
    for node in range(4):
        data = simulate_imu_measurements(fleet.ideal_imu[:, node], fleet.dt,
            np.random.default_rng(offset + 10000 * seed + 100 * node + 1),
            **{k: v for k, v in base['imu'].items() if k.startswith('sigma_')})
        measured[:, node] = data.measured
    ranges = _range_data(fleet.state, base['uwb']['sigma_range_m'],
                        np.random.default_rng(offset + 10000 * seed + 9001))
    initial = fleet.state[0].copy()
    initial[:, 2:4] += np.random.default_rng(offset + 10000 * seed + 8001).normal(
        0, base['initialization']['velocity_error_std_mps'], (4, 2))
    p0 = _initial_covariance(base['initialization'], 4)
    masks = {b: baseline_mask(exchange_capacity(len(fleet.t), 10, b, 6), 4, 'cyclic')
             for b in cfg['budgets']}
    scored = cfg['scored_nodes_zero_based']
    sample = np.arange(0, len(fleet.t), round(1 / (cfg['profile_sample_rate_hz'] * fleet.dt)))
    rows, profiles, representative, audits = [], [], [], []
    baseline_zero = None
    with (ROOT / 'results/phase2/p2e/seed_metrics.csv').open() as handle:
        inherited = {int(r['budget']): r for r in csv.DictReader(handle)
                     if int(r['seed']) == seed and r['condition'] == 'reference' and r['mode'] == 'frozen_q'}
    for case, (positions, refmask, sigma) in reference_cases(fleet, cfg, seed, offset).items():
        if case not in cfg['cases']:
            continue
        for budget in cfg['budgets']:
            kwargs = {} if positions is None else dict(position_reference=positions,
                position_reference_mask=refmask, position_reference_node=cfg['reference_node_zero_based'],
                sigma_position_reference_m=sigma)
            result = run_joint_ekf(measured, ranges, masks[budget], initial, p0, fleet.dt,
                sigma_range_m=base['uwb']['sigma_range_m'], **base['ekf'], **kwargs)
            if result.range_attempt_count != budget or result.range_update_count != budget:
                raise RuntimeError('Clean UWB budget mismatch')
            if result.position_reference_update_count != (0 if refmask is None else int(refmask.sum())):
                raise RuntimeError('Reference count mismatch')
            if case == 'exact_100hz' and not np.array_equal(result.state[:, 0, :2], fleet.state[:, 0, :2]):
                raise RuntimeError('Exact reference position constraint violated')
            if case == 'unanchored':
                full = _fleet_metrics(result.state, fleet.state)
                for key, raw in [('fleet_rmse_m', 'fleet_position_rmse_m'),
                                 ('centroid_rmse_m', 'centroid_position_rmse_m'),
                                 ('worst_node_rmse_m', 'worst_node_position_rmse_m'),
                                 ('relative_shape_rmse_m', 'relative_shape_rmse_m'),
                                 ('yaw_rmse_deg', 'fleet_yaw_rmse_deg')]:
                    if not np.isclose(full[raw], float(inherited[budget][key]), atol=1e-10, rtol=1e-10):
                        raise RuntimeError(f'Inherited endpoint mismatch {seed} {budget} {key}')
                audits.append(dict(seed=seed, budget=budget, inherited_source='p2e_reference_frozen_q', matched=True))
            state, truth = result.state[:, scored], fleet.state[:, scored]
            if budget == 0:
                if case == 'unanchored':
                    baseline_zero = state.copy()
                elif not np.array_equal(state, baseline_zero):
                    raise RuntimeError('Reference leaked into disconnected non-reference nodes')
            metric = _fleet_metrics(state, truth)
            errors = _time_errors(state, truth)
            if not np.isclose(metric['fleet_position_rmse_m']**2,
                metric['centroid_position_rmse_m']**2 + metric['relative_shape_rmse_m']**2):
                raise RuntimeError('Three-node decomposition mismatch')
            ref_error = np.linalg.norm(result.state[:, 0, :2] - fleet.state[:, 0, :2], axis=1)
            row = dict(seed=seed, case=case, budget=budget, communication_ratio=budget/3600,
                attempted=result.range_attempt_count, received=result.range_received_count,
                accepted=result.range_update_count, reference_updates=result.position_reference_update_count,
                fleet_rmse_m=metric['fleet_position_rmse_m'], worst_node_rmse_m=metric['worst_node_position_rmse_m'],
                centroid_rmse_m=metric['centroid_position_rmse_m'], relative_shape_rmse_m=metric['relative_shape_rmse_m'],
                yaw_rmse_deg=metric['fleet_yaw_rmse_deg'],
                late_fleet_rmse_m=float(np.sqrt(np.mean(errors['fleet_rms_error_m'][4000:]**2))),
                reference_rmse_m=float(np.sqrt(np.mean(ref_error**2))))
            for node, value in zip(scored, metric['node_position_rmse_m']):
                row[f'node{node+1}_rmse_m'] = value
            rows.append(row)
            profiles.append(dict(case=case, budget=budget,
                **{k: v[sample].tolist() for k, v in errors.items()},
                reference_error_m=ref_error[sample].tolist()))
            if seed == base['campaign']['representative_seed'] and budget in [0, 360, 3600]:
                for index in sample:
                    for node in range(4):
                        representative.append(dict(case=case,budget=budget,time_s=float(fleet.t[index]),node=node+1,
                            true_x_m=float(fleet.state[index,node,0]),true_y_m=float(fleet.state[index,node,1]),
                            estimate_x_m=float(result.state[index,node,0]),estimate_y_m=float(result.state[index,node,1])))
    return dict(rows=rows, profiles=profiles, representative=representative, audits=audits,
                times=fleet.t[sample].tolist())


def aggregate(runs, cfg, base, digest, out):
    rows = [r for run in runs for r in run['rows']]
    _write_csv(out/'seed_metrics.csv', rows)
    _write_csv(out/'endpoint_audit.csv', [r for run in runs for r in run['audits']])
    _write_csv(out/'representative.csv', runs[0]['representative'], significant_digits=8)
    lookup = {(r['seed'],r['case'],r['budget']):r for r in rows}
    summary = dict(protocol=cfg, base_config_snapshot=base, n_seeds=len(runs), n_runs=len(rows),
                   implementation_fingerprint=digest, numpy_version=np.__version__,
                   results=[], comparisons=[], target_budgets=[])
    boot = np.random.default_rng(cfg['bootstrap_seed']).integers(0,len(runs),
                       (cfg['bootstrap_resamples'],len(runs)))
    profile_rows, comparisons = [], []
    for case in cfg['cases']:
        for budget in cfg['budgets']:
            group = [lookup[s,case,budget] for s in range(len(runs))]
            metrics = {}
            for key in METRICS:
                a=np.array([r[key] for r in group]);metrics[key]=dict(mean=float(a.mean()),
                    std=float(a.std(ddof=1)),p95=float(np.quantile(a,.95)))
            summary['results'].append(dict(case=case,budget=budget,
                reference_updates=group[0]['reference_updates'], metrics=metrics,
                seeds_meeting_absolute_target=int(sum(r['fleet_rmse_m']<=cfg['absolute_target_m'] for r in group))))
            histories = [next(p for p in run['profiles'] if (p['case'],p['budget'])==(case,budget)) for run in runs]
            for i,t in enumerate(runs[0]['times']):
                row=dict(case=case,budget=budget,time_s=t)
                for key in ['fleet_rms_error_m','centroid_error_m','shape_rms_error_m','reference_error_m']:
                    vals=[p[key][i] for p in histories];row['mean_'+key]=float(np.mean(vals));row['std_'+key]=float(np.std(vals,ddof=1))
                profile_rows.append(row)
            if case!='unanchored':
                refs=[('unanchored','unanchored')]
                if case=='noisy_1hz_outage':refs.append(('noisy_1hz','noisy_1hz'))
                for label,ref in refs:
                    item=dict(case=case,budget=budget,reference=label,metrics={})
                    flat=dict(case=case,budget=budget,reference=label)
                    for key in METRICS[:-1]:
                        diff=np.array([lookup[s,case,budget][key]-lookup[s,ref,budget][key] for s in range(len(runs))])
                        ci=np.quantile(diff[boot].mean(axis=1),[.025,.975])
                        item['metrics'][key]=dict(mean=float(diff.mean()),std=float(diff.std(ddof=1)),
                                                 wins=int((diff<0).sum()),ci95=ci.tolist())
                        flat.update({f'mean_{key}':float(diff.mean()),f'ci95_low_{key}':float(ci[0]),
                                     f'ci95_high_{key}':float(ci[1]),f'wins_{key}':int((diff<0).sum())})
                    comparisons.append(flat);summary['comparisons'].append(item)
        points=[r for r in summary['results'] if r['case']==case]
        targets={}
        for stat in ['mean','p95']:
            passed=[r['budget'] for r in points if r['metrics']['fleet_rmse_m'][stat]<=cfg['absolute_target_m']]
            targets[f'smallest_tested_{stat}_absolute_target_budget']=min(passed) if passed else None
        summary['target_budgets'].append(dict(case=case,**targets))
    _write_csv(out/'time_profile.csv',profile_rows,significant_digits=8)
    _write_csv(out/'paired_comparisons.csv',comparisons)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# P2G continuing global reference','',
           'All headline metrics score nodes 2–4. Reference costs are separate. No global alignment is applied.', '',
           '| Case | UWB attempts | Reference updates | Absolute fleet mean ± SD (m) | Centroid mean (m) | Shape mean (m) | Worst mean (m) |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for r in summary['results']:
        m=r['metrics'];v=m['fleet_rmse_m'];lines.append(f"| {r['case']} | {r['budget']} | {r['reference_updates']} | {v['mean']:.4f} ± {v['std']:.4f} | {m['centroid_rmse_m']['mean']:.4f} | {m['relative_shape_rmse_m']['mean']:.4f} | {m['worst_node_rmse_m']['mean']:.4f} |")
    (out/'summary.md').write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--seeds',type=int,default=None,help='Smoke check only; full campaign uses base seed count')
    args=parser.parse_args();cfg=yaml.safe_load((ROOT/'configs/phase2g.yaml').read_text())
    base=yaml.safe_load((ROOT/cfg['base_config']).read_text());n=args.seeds or base['campaign']['n_seeds']
    fingerprint=hashlib.sha256(json.dumps([cfg,base],sort_keys=True).encode())
    for path in ['experiments/run_phase2g.py','experiments/run_phase2a.py','experiments/run_phase2b.py',
                 'src/mt_ag/joint_ekf.py','src/mt_ag/imu.py','src/mt_ag/link_selection.py','src/mt_ag/fleet_simulation.py']:
        fingerprint.update((ROOT/path).read_bytes())
    digest=fingerprint.hexdigest();out=ROOT/'results/phase2/p2g';checkpoint=out/'_checkpoints'
    checkpoint.mkdir(parents=True,exist_ok=True);runs={}
    for seed in range(n):
        p=checkpoint/f'{seed}.json'
        if p.exists():
            saved=json.loads(p.read_text())
            if saved['fingerprint']==digest:runs[seed]=saved['run']
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run_seed,s,cfg,base):s for s in range(n) if s not in runs}
        for future in as_completed(futures):
            seed=futures[future];run=future.result();temp=checkpoint/f'{seed}.tmp'
            temp.write_text(json.dumps(dict(fingerprint=digest,run=run)));os.replace(temp,checkpoint/f'{seed}.json')
            runs[seed]=run;print(f'P2G seed {seed+1} completed; {len(runs)}/{n} saved',flush=True)
    if n==base['campaign']['n_seeds']:
        aggregate([runs[s] for s in range(n)],cfg,base,digest,out)
        print(f'P2G complete: {n*len(cfg["cases"])*len(cfg["budgets"])} runs',flush=True)
    else:print('Smoke seeds saved; aggregate withheld until full campaign',flush=True)

if __name__=='__main__':main()
