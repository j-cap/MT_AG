"""P2F: provenance-preserving consolidation; no simulation or estimator tuning."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
METRICS = ('fleet_rmse_m', 'worst_node_rmse_m', 'centroid_rmse_m',
           'relative_shape_rmse_m', 'yaw_rmse_deg')
CONDITION_KEYS = ('sensor_condition', 'q_treatment', 'range_condition', 'update')
DESIGN_KEYS = CONDITION_KEYS + ('schedule', 'timing', 'budget')


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def pareto_mask(costs, errors, tolerance=1e-12):
    """Two-objective minimization. Exact ties remain on the empirical frontier."""
    costs, errors = np.asarray(costs), np.asarray(errors)
    if len(costs) != len(errors) or not np.all(np.isfinite(errors)):
        raise ValueError('Invalid frontier inputs')
    keep = np.ones(len(costs), dtype=bool)
    for index, (cost, error) in enumerate(zip(costs, errors)):
        weak = (costs <= cost) & (errors <= error + tolerance)
        strict = (costs < cost) | (errors < error - tolerance)
        keep[index] = not np.any(weak & strict)
    return keep


def target_decision(points, target, statistic):
    """Smallest tested passing budget; selects lowest error at tied budgets."""
    eligible = [p for p in points if p[f'{statistic}_relative_shape_rmse_m'] <= target]
    return min(eligible, key=lambda p: (p['budget'], p[f'{statistic}_relative_shape_rmse_m'], p['design_id'])) if eligible else None


def canonical_row(source, raw):
    sensor, q, case, update = 'reference', 'frozen_q', 'clean', 'naive'
    if source == 'p2b':
        budget = int(raw['uwb_exchanges'])
        schedule, timing = 'uniform', 'uniform_rounds'
    else:
        budget = int(raw.get('budget', raw.get('attempted')))
        if source == 'p2c':
            schedule = raw['policy']
        elif source == 'p2d':
            schedule = 'uniform' if raw['policy'] == 'p2b_uniform' else raw['policy']
            case, update = raw['case'], raw['mode']
        elif source == 'p2e':
            sensor, q = raw['condition'], raw['mode']
            schedule = 'cyclic'
        else:
            raise ValueError(source)
        timing = 'uniform_rounds' if schedule == 'uniform' else 'p2c_capacity'
    if budget == 0:
        schedule, timing = 'none', 'none'
    elif budget == 3600:
        schedule, timing = 'all_pairs', 'full_10hz'
    row = dict(seed=int(raw['seed']), sensor_condition=sensor, q_treatment=q,
               range_condition=case, update=update, schedule=schedule, timing=timing,
               budget=budget, communication_ratio=budget/3600,
               attempted=int(raw.get('attempted', budget)),
               received=int(raw.get('received', budget)), accepted=int(raw.get('accepted', budget)))
    row.update({k: float(raw[k]) for k in METRICS})
    return row


def consolidate(raw_sources, n_seeds):
    records, provenance, audits = {}, {}, []
    for source in ['p2b', 'p2c', 'p2d', 'p2e']:
        for raw in raw_sources[source]:
            row = canonical_row(source, raw)
            key = tuple(row[k] for k in DESIGN_KEYS) + (row['seed'],)
            if key in records:
                old = records[key]
                for metric in METRICS + ('attempted', 'received', 'accepted'):
                    if not np.isclose(row[metric], old[metric], atol=1e-10, rtol=1e-10):
                        raise ValueError(f'Duplicate endpoint mismatch: {source} {key} {metric}')
                audits.append(dict(source=source, seed=row['seed'], kind='duplicate_endpoint',
                                   reference_sources=','.join(provenance[key])))
                provenance[key].append(source)
            else:
                records[key], provenance[key] = row, [source]
    for raw in raw_sources['p2a']:
        for prefix, budget in [('imu', 0), ('full', 3600)]:
            key = ('reference', 'frozen_q', 'clean', 'naive',
                   'none' if budget == 0 else 'all_pairs',
                   'none' if budget == 0 else 'full_10hz', budget, int(raw['seed']))
            for metric in METRICS:
                if not np.isclose(float(raw[f'{prefix}_{metric}']), records[key][metric], atol=1e-10, rtol=1e-10):
                    raise ValueError(f'P2A mismatch {prefix} seed {raw["seed"]} {metric}')
            provenance[key].append('p2a')
            audits.append(dict(source='p2a', seed=int(raw['seed']), kind='baseline_endpoint', reference_sources='p2b'))
    groups = {}
    for key, row in records.items():
        design = key[:-1]
        groups.setdefault(design, []).append(row)
        if not 0 <= row['accepted'] <= row['received'] <= row['attempted'] == row['budget']:
            raise ValueError('Communication accounting mismatch')
        if not np.isclose(row['fleet_rmse_m']**2, row['centroid_rmse_m']**2 + row['relative_shape_rmse_m']**2, rtol=1e-9, atol=1e-10):
            raise ValueError('Squared error decomposition mismatch')
        row['sources'] = ','.join(provenance[key])
    for group in groups.values():
        if sorted(r['seed'] for r in group) != list(range(n_seeds)):
            raise ValueError('Missing or duplicated matched seed')
    return sorted(records.values(), key=lambda r: tuple(r[k] for k in DESIGN_KEYS) + (r['seed'],)), groups, audits


def main():
    cfg = yaml.safe_load((ROOT/'configs/phase2f.yaml').read_text())
    raw_sources, source_hashes = {}, {}
    for source in cfg['sources']:
        path = ROOT/f'results/phase2/{source}/seed_metrics.csv'
        raw_sources[source] = list(csv.DictReader(path.open()))
        source_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    rows, groups, audits = consolidate(raw_sources, cfg['n_seeds'])
    n = cfg['n_seeds']
    points = []
    point_lookup = {}
    for index, (design, group) in enumerate(sorted(groups.items())):
        group = sorted(group, key=lambda r: r['seed'])
        point = dict(design_id=f'd{index:03d}', **dict(zip(DESIGN_KEYS, design)),
                     n_seeds=n, sources=','.join(sorted({s for r in group for s in r['sources'].split(',')})))
        for k in METRICS + ('received', 'accepted'):
            values = np.array([r[k] for r in group])
            point[f'mean_{k}'] = float(values.mean())
            point[f'std_{k}'] = float(values.std(ddof=1))
            point[f'p95_{k}'] = float(np.quantile(values, .95))
        point['seeds_meeting_shape_target'] = sum(r['relative_shape_rmse_m'] <= cfg['shape_target_m'] for r in group)
        points.append(point)
        point_lookup[design] = point
        for r in group:
            r['design_id'] = point['design_id']
    for condition in sorted({tuple(p[k] for k in CONDITION_KEYS) for p in points}):
        group = [p for p in points if tuple(p[k] for k in CONDITION_KEYS) == condition]
        for metric in ['relative_shape_rmse_m', 'fleet_rmse_m', 'worst_node_rmse_m']:
            keep = pareto_mask([p['budget'] for p in group], [p[f'mean_{metric}'] for p in group], cfg['pareto_numeric_tolerance'])
            for p, value in zip(group, keep):
                p[f'pareto_mean_{metric}'] = bool(value)
        keep = pareto_mask([p['budget'] for p in group], [p['p95_relative_shape_rmse_m'] for p in group], cfg['pareto_numeric_tolerance'])
        for p, value in zip(group, keep):
            p['pareto_p95_relative_shape_rmse_m'] = bool(value)
    targets = []
    # Clean sensor frontiers are comparable by Q/sensor; reliability cases are
    # separately labeled single-budget evidence, not inferred stress curves.
    for condition in sorted({tuple(p[k] for k in CONDITION_KEYS) for p in points}):
        group = [p for p in points if tuple(p[k] for k in CONDITION_KEYS) == condition]
        for schedule_filter in (['all_available_designs', 'cyclic_reference', 'uniform_reference'] if condition == ('reference', 'frozen_q', 'clean', 'naive') else ['all_available_designs']):
            allowed = {'cyclic_reference': {'none', 'cyclic', 'all_pairs'},
                       'uniform_reference': {'none', 'uniform', 'all_pairs'}}
            selected = group if schedule_filter == 'all_available_designs' else [p for p in group if p['schedule'] in allowed[schedule_filter]]
            target = dict(zip(CONDITION_KEYS, condition))
            target.update(schedule_scope=schedule_filter, tested_budgets=','.join(str(b) for b in sorted({p['budget'] for p in selected})),
                          evidence='single_budget' if len({p['budget'] for p in selected}) == 1 else 'budget_sweep')
            for statistic in ['mean', 'p95']:
                best = target_decision(selected, cfg['shape_target_m'], statistic)
                target[f'{statistic}_budget'] = best['budget'] if best else ''
                target[f'{statistic}_schedule'] = best['schedule'] if best else ''
                target[f'{statistic}_design_id'] = best['design_id'] if best else ''
                target[f'{statistic}_passing_seeds'] = best['seeds_meeting_shape_target'] if best else ''
            targets.append(target)
    boot = np.random.default_rng(cfg['bootstrap_seed']).integers(0, n, (cfg['bootstrap_resamples'], n))
    paired = []
    def compare(a, b, label, timing_matched):
        first = sorted(groups[a], key=lambda r: r['seed'])
        second = sorted(groups[b], key=lambda r: r['seed'])
        item = dict(comparison=label, design_id=point_lookup[a]['design_id'], reference_id=point_lookup[b]['design_id'],
                    equal_attempts=a[-1] == b[-1], timing_matched=timing_matched)
        for metric in ['relative_shape_rmse_m', 'fleet_rmse_m', 'worst_node_rmse_m']:
            diff = np.array([x[metric]-y[metric] for x, y in zip(first, second)])
            lo, hi = np.quantile(diff[boot].mean(axis=1), [.025,.975])
            item[f'mean_{metric}'] = float(diff.mean())
            item[f'std_{metric}'] = float(diff.std(ddof=1))
            item[f'ci95_low_{metric}'], item[f'ci95_high_{metric}'] = float(lo), float(hi)
            item[f'wins_{metric}'] = int(np.sum(diff < 0))
        paired.append(item)
    base = ('reference','frozen_q','clean','naive')
    for budget in [180,360,720]:
        cyclic = base + ('cyclic','p2c_capacity',budget)
        for policy in ['random','geometry','information']:
            compare(base + (policy,'p2c_capacity',budget), cyclic, f'{policy}_minus_cyclic', True)
        compare(base + ('uniform','uniform_rounds',budget), cyclic, 'uniform_minus_cyclic', False)
    for case in ['clean','loss30','loss50','network_outage','node_outage','outliers','nlos']:
        for update in ['naive','gated']:
            cyclic = ('reference','frozen_q',case,update,'cyclic','p2c_capacity',360)
            info = ('reference','frozen_q',case,update,'information','p2c_capacity',360)
            if cyclic in groups and info in groups:
                compare(info, cyclic, f'{case}_{update}_information_minus_cyclic', True)
        naive = ('reference','frozen_q',case,'naive','cyclic','p2c_capacity',360)
        gated = ('reference','frozen_q',case,'gated','cyclic','p2c_capacity',360)
        if gated in groups:
            compare(gated, naive, f'{case}_cyclic_gate_minus_naive', True)
    out = ROOT/'results/phase2/p2f'
    out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'seed_catalog.csv', rows)
    write_csv(out/'design_summary.csv', points)
    write_csv(out/'target_decisions.csv', targets)
    write_csv(out/'paired_comparisons.csv', paired)
    write_csv(out/'endpoint_audit.csv', audits)
    summary = dict(protocol=cfg, source_sha256=source_hashes, unique_designs=len(points), unique_seed_results=len(rows),
                   source_run_records=40+sum(len(raw_sources[x]) for x in ['p2b','p2c','p2d','p2e']),
                   endpoint_checks=len(audits), points=points, targets=targets, paired=paired)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# P2F conditional accuracy-communication consolidation','',
           f'{len(points)} unique designs, {len(rows)} unique seed records, {len(audits)} matched endpoint checks.',
           'No new simulation. Empirical Pareto labels minimize attempted cost and mean error within fixed sensor/Q/range/update conditions.',
           'Bootstrap intervals are descriptive and unadjusted. Single-budget reliability evidence is not a communication curve.','',
           '| Sensor | Q | Range case | Update | Scope | Tested budgets | Mean 10 cm target | p95 10 cm target |',
           '|---|---|---|---|---|---|---|---|']
    for r in targets:
        lines.append(f"| {r['sensor_condition']} | {r['q_treatment']} | {r['range_condition']} | {r['update']} | {r['schedule_scope']} | {r['tested_budgets']} | {r['mean_budget']} {r['mean_schedule']} | {r['p95_budget']} {r['p95_schedule']} |")
    (out/'summary.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:summary[k] for k in ['unique_designs','unique_seed_results','source_run_records','endpoint_checks']},indent=2))

if __name__ == '__main__':
    main()
