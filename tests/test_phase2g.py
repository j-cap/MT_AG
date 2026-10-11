"""Reference conditioning, isolation and counting are scientific requirements."""
import numpy as np
from mt_ag.joint_ekf import _position_reference_update, run_joint_ekf


def test_reference_conditions_joint_cross_covariance_without_truth_velocity():
    p = np.eye(10)
    p[0, 2] = p[2, 0] = .2
    p[0, 5] = p[5, 0] = .3
    state, post = _position_reference_update(np.zeros(10), p, np.array([2., 1.]), 0, 0.)
    assert np.allclose(state[[0, 1, 2, 5]], [2, 1, .4, .6])
    assert state[3] == state[4] == 0
    assert np.allclose(post[:2], 0)
    assert np.linalg.eigvalsh(post).min() >= -1e-12
    assert np.allclose(post[2, 2], .96)


def test_noisy_reference_matches_linear_gaussian_conditioning():
    p = np.eye(10)
    s, post = _position_reference_update(np.zeros(10), p, np.ones(2), 0, .5)
    assert np.allclose(s[:2], .8)
    assert np.allclose(post[:2, :2], np.eye(2)*.2)


def test_reference_without_uwb_leaves_other_nodes_byte_identical():
    args = dict(imu_measurements=np.zeros((3, 2, 3)),
                pairwise_ranges=np.zeros((4, 2, 2)), range_mask=np.zeros((4, 2, 2), bool),
                initial_state=np.zeros((2, 5)), initial_covariance=np.eye(10), dt=.01,
                sigma_range_m=.12, sigma_accel_process_mps2=.03, sigma_gyro_process_rps=.0015)
    base = run_joint_ekf(**args)
    ref = run_joint_ekf(**args, position_reference=np.ones((4, 2)),
                        position_reference_mask=np.array([False, True, False, True]))
    assert np.array_equal(base.state[:, 1], ref.state[:, 1])
    assert ref.position_reference_update_count == 2
    assert ref.range_attempt_count == ref.range_update_count == 0
    assert np.array_equal(ref.state[[1, 3], 0, :2], np.ones((2, 2)))


def test_noisy_outage_reuses_reports_and_excludes_only_outage_times():
    from experiments.run_phase2g import reference_cases
    from mt_ag.fleet_simulation import generate_four_node_fleet
    import yaml
    from pathlib import Path
    cfg = yaml.safe_load(Path('configs/phase2g.yaml').read_text())
    fleet = generate_four_node_fleet(dt=.01, duration=60)
    cases = reference_cases(fleet, cfg, 0, 700000)
    noisy, mask, _ = cases['noisy_1hz']
    intermittent, outage_mask, _ = cases['noisy_1hz_outage']
    assert np.array_equal(noisy, intermittent)
    assert mask.sum() == 60 and outage_mask.sum() == 50
    assert not outage_mask[2000] and outage_mask[3000]
    assert cases['exact_100hz'][1].sum() == 6000
