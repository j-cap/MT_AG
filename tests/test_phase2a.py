from __future__ import annotations

import numpy as np

from mt_ag.fleet_simulation import all_pairs, generate_four_node_fleet, pairwise_distances
from mt_ag.imu import mechanize_planar
from mt_ag.joint_ekf import run_joint_ekf


def test_four_node_paths_are_smooth_and_separated():
    fleet = generate_four_node_fleet(dt=0.01, duration=60.0)
    assert fleet.state.shape == (6001, 4, 5)
    assert fleet.ideal_imu.shape == (6000, 4, 3)

    speed = np.linalg.norm(fleet.state[:, :, 2:4], axis=2)
    assert np.min(speed) > 0.15
    assert np.max(speed) < 1.0

    distances = pairwise_distances(fleet.state[:, :, :2])
    assert len(distances) == 6
    assert min(float(np.min(values)) for values in distances.values()) > 1.5


def test_trajectory_first_imu_reconstructs_truth():
    fleet = generate_four_node_fleet(dt=0.01, duration=10.0)
    for node in range(4):
        reconstructed = mechanize_planar(
            fleet.state[0, node],
            fleet.ideal_imu[:, node],
            fleet.dt,
        )
        position_error = np.linalg.norm(
            reconstructed[:, :2] - fleet.state[:, node, :2],
            axis=1,
        )
        assert float(np.max(position_error)) < 2.0e-5


def test_joint_ekf_exact_inputs_track_truth():
    fleet = generate_four_node_fleet(dt=0.01, duration=2.0)
    n_time, n_nodes, _ = fleet.state.shape
    ranges = np.zeros((n_time, n_nodes, n_nodes), dtype=float)
    mask = np.zeros_like(ranges, dtype=bool)

    for i, j in all_pairs(n_nodes):
        distance = np.linalg.norm(
            fleet.state[:, i, :2] - fleet.state[:, j, :2],
            axis=1,
        )
        ranges[:, i, j] = distance
        ranges[:, j, i] = distance
        mask[10::10, i, j] = True
        mask[10::10, j, i] = True

    covariance = np.eye(5 * n_nodes) * 1.0e-8
    result = run_joint_ekf(
        fleet.ideal_imu,
        ranges,
        mask,
        fleet.state[0],
        covariance,
        fleet.dt,
        sigma_range_m=1.0e-4,
        sigma_accel_process_mps2=1.0e-5,
        sigma_gyro_process_rps=1.0e-6,
    )

    error = np.linalg.norm(result.state[:, :, :2] - fleet.state[:, :, :2], axis=2)
    assert float(np.max(error)) < 5.0e-4
    assert result.range_update_count == 20 * len(all_pairs(n_nodes))
