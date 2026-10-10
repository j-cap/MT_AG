import importlib.util
from pathlib import Path
import numpy as np

spec = importlib.util.spec_from_file_location('phase2f', Path(__file__).resolve().parents[1]/'experiments/analyze_phase2f.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_frontier_handles_lower_cost_dominance_and_exact_ties():
    result = module.pareto_mask([1,2,3,4,2], [3,2,2.5,1,2])
    assert np.array_equal(result, [True,True,False,True,True])


def test_target_uses_smallest_tested_budget_without_interpolation():
    points=[dict(design_id='a',budget=180,mean_relative_shape_rmse_m=.11),
            dict(design_id='b',budget=720,mean_relative_shape_rmse_m=.07),
            dict(design_id='c',budget=360,mean_relative_shape_rmse_m=.09)]
    assert module.target_decision(points,.10,'mean')['budget']==360
    assert module.target_decision(points,.06,'mean') is None


def test_conflicting_repeated_endpoint_is_rejected():
    def raw(source,shape):
        values=dict(seed='0',budget='360',policy='cyclic',fleet_rmse_m='1',
                    worst_node_rmse_m='1',centroid_rmse_m='.9',relative_shape_rmse_m=str(shape),yaw_rmse_deg='.1')
        if source=='p2d':values.update(case='clean',mode='naive',attempted='360')
        return values
    sources=dict(p2a=[],p2b=[],p2c=[raw('p2c',.1)],p2d=[raw('p2d',.2)],p2e=[])
    try:
        module.consolidate(sources,1)
    except ValueError as error:
        assert 'Duplicate endpoint mismatch' in str(error)
    else:
        raise AssertionError('Conflicting experiments were silently pooled')
