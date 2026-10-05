# ruff: noqa: E402
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from experiments.run_phase2e import rescale_imu, process_settings
from mt_ag.joint_ekf import run_joint_ekf

def test_shared_latent_noise_and_bias_are_scaled_separately():
    ideal = np.arange(18).reshape(6, 3).astype(float)
    bias = np.full_like(ideal, 0.04)
    noise = np.linspace(-0.02, 0.02, 18).reshape(6, 3)
    measured = ideal + bias + noise
    assert np.array_equal(rescale_imu(ideal, measured, bias, 1, 1), measured)
    assert np.allclose(rescale_imu(ideal, measured, bias, 4, 1), ideal + bias + 4 * noise)
    assert np.allclose(rescale_imu(ideal, measured, bias, 1, 10), ideal + 10 * bias + noise)
    assert np.array_equal(rescale_imu(ideal, measured, bias, 0, 0), ideal)

def test_process_adjustment_changes_covariance_without_changing_dead_reckoning():
    base = dict(sigma_accel_process_mps2=0.03, sigma_gyro_process_rps=0.0015)
    assert process_settings(base, 'frozen_q', 4) == base
    adjusted = process_settings(base, 'white_rescaled_q', 4)
    assert adjusted['sigma_accel_process_mps2'] == 0.12
    args = dict(imu_measurements=np.ones((3, 2, 3)) * 0.01, pairwise_ranges=np.zeros((4, 2, 2)), range_mask=np.zeros((4, 2, 2), dtype=bool), initial_state=np.array([[0.0, 0, 0, 0, 0], [1.0, 0, 0, 0, 0]]), initial_covariance=np.eye(10) * 0.01, dt=0.01, sigma_range_m=0.12)
    a = run_joint_ekf(**args, **base)
    b = run_joint_ekf(**args, **adjusted)
    assert np.array_equal(a.state, b.state)
    assert not np.array_equal(a.covariance, b.covariance)
