"""Failure accounting and measurement-blind selection are protocol requirements."""
import numpy as np
from mt_ag.joint_ekf import run_joint_ekf


def _args():
    return dict(imu_measurements=np.zeros((2, 2, 3)),
                pairwise_ranges=np.ones((3, 2, 2)),
                range_mask=np.ones((3, 2, 2), dtype=bool),
                initial_state=np.array([[0., 0, 0, 0, 0], [1., 0, 0, 0, 0]]),
                initial_covariance=np.eye(10)*.01, dt=.01,
                sigma_range_m=.12, sigma_accel_process_mps2=.03,
                sigma_gyro_process_rps=.0015)


def test_lost_attempts_consume_budget_and_do_not_change_eligibility():
    args = _args()
    seen = []
    def selector(t, state, covariance, eligible):
        seen.append(eligible)
        return eligible
    result = run_joint_ekf(**args, range_selector=selector,
                           reception_mask=np.zeros((3, 2, 2), dtype=bool))
    assert (result.range_attempt_count, result.range_received_count, result.range_update_count) == (2, 0, 0)
    assert seen == [((0, 1),), ((0, 1),)]


def test_gate_rejection_preserves_prediction_and_counts_received():
    args = _args()
    predicted = run_joint_ekf(**args, reception_mask=np.zeros((3, 2, 2), dtype=bool))
    args['pairwise_ranges'] *= 20
    result = run_joint_ekf(**args, innovation_gate_nis=9)
    assert (result.range_attempt_count, result.range_received_count,
            result.range_update_count, result.range_rejected_count) == (2, 2, 0, 2)
    assert np.array_equal(result.state, predicted.state)
    assert np.array_equal(result.covariance, predicted.covariance)


def test_clean_receptions_are_accepted_and_nonfinite_ranges_are_rejected():
    args = _args()
    result = run_joint_ekf(**args, innovation_gate_nis=9)
    assert (result.range_attempt_count, result.range_received_count, result.range_update_count) == (2, 2, 2)
    args['pairwise_ranges'][1, 0, 1] = np.nan
    result = run_joint_ekf(**args)
    assert (result.range_received_count, result.range_update_count, result.range_rejected_count) == (2, 1, 1)
