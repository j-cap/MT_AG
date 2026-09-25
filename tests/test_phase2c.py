from __future__ import annotations

import numpy as np

from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet
from mt_ag.joint_ekf import run_joint_ekf
from mt_ag.link_selection import (
    PredictiveLinkSelector,
    baseline_mask,
    eligible_mask,
    exchange_capacity,
)


def test_matched_budget_and_tick_capacity():
    for budget in (180, 360, 720):
        capacity = exchange_capacity(6001, 10, budget, 6)
        assert capacity.sum() == budget
        for policy in ("cyclic", "random"):
            rng = np.random.default_rng(123) if policy == "random" else None
            mask = baseline_mask(capacity, 4, policy, rng)
            assert np.array_equal(np.count_nonzero(mask, axis=(1, 2)) // 2, capacity)
            assert np.array_equal(mask, mask.transpose(0, 2, 1))
        assert max(capacity) <= 2


def test_online_selection_is_pre_measurement_and_has_exact_cost():
    fleet = generate_four_node_fleet(dt=.01, duration=2.0)
    n_time, n_nodes, _ = fleet.state.shape
    cap = exchange_capacity(n_time, 10, 4, len(all_pairs(n_nodes)))
    mask = eligible_mask(cap, n_nodes)
    ranges = np.zeros((n_time, n_nodes, n_nodes))
    for i, j in all_pairs(n_nodes):
        true_distance = np.linalg.norm(fleet.state[:, i, :2] - fleet.state[:, j, :2], axis=1)
        ranges[:, i, j] = ranges[:, j, i] = true_distance
    kwargs = dict(imu_measurements=fleet.ideal_imu, range_mask=mask,
                  initial_state=fleet.state[0], initial_covariance=np.eye(20) * .0025,
                  dt=.01, sigma_range_m=.12,
                  sigma_accel_process_mps2=.03, sigma_gyro_process_rps=.0015)
    for policy in ("geometry", "information"):
        selector_a = PredictiveLinkSelector(policy, cap, .12, .01, n_nodes)
        reference = run_joint_ekf(pairwise_ranges=ranges, range_selector=selector_a, **kwargs)
        assert reference.range_update_count == 4
        assert len(selector_a.records) == 4
        # Alter the first range observation: the first pair choice must already be fixed.
        changed = ranges.copy()
        changed[50] += 1.0
        selector_b = PredictiveLinkSelector(policy, cap, .12, .01, n_nodes)
        run_joint_ekf(pairwise_ranges=changed, range_selector=selector_b, **kwargs)
        assert selector_a.records[0][:3] == selector_b.records[0][:3]
