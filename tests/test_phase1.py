import numpy as np

from mt_ag.geometry import wrap_angle
from mt_ag.imu import mechanize_planar, simulate_imu_measurements
from mt_ag.observability import observability_history_planar_imu, relative_geometry_metrics
from mt_ag.particle_filter import systematic_resample
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory


def test_angle_wrap():
    assert np.isclose(wrap_angle(np.pi), -np.pi)
    assert np.isclose(wrap_angle(3 * np.pi), -np.pi)


def test_ideal_imu_reconstructs_truth():
    trajectory = generate_curved_trajectory(dt=0.1, duration=10.0)
    reconstruction = mechanize_planar(trajectory.state[0], trajectory.ideal_imu, trajectory.dt)
    assert np.max(np.abs(reconstruction - trajectory.state)) < 1e-12


def test_noisy_imu_differs_from_ideal_measurement():
    trajectory = generate_curved_trajectory(dt=0.1, duration=2.0)
    rng = np.random.default_rng(1)
    data = simulate_imu_measurements(trajectory.ideal_imu, trajectory.dt, rng)
    assert not np.allclose(data.measured, trajectory.ideal_imu)


def test_ideal_uwb_is_euclidean_range():
    p = np.array([[0.0, 0.0], [3.0, 4.0]])
    a = np.array([[0.0, 0.0], [0.0, 0.0]])
    assert np.allclose(generate_uwb_ranges(p, a), [0.0, 5.0])


def test_systematic_resampling_preserves_particle_count():
    rng = np.random.default_rng(1)
    weights = np.array([0.05, 0.15, 0.8])
    indices = systematic_resample(weights, rng)
    assert len(indices) == len(weights)
    assert np.all((indices >= 0) & (indices < len(weights)))


def test_constant_bearing_auxiliary_keeps_line_of_sight_direction():
    trajectory = generate_curved_trajectory(dt=0.1, duration=5.0)
    auxiliary = auxiliary_trajectory(
        trajectory.t,
        "constant_bearing",
        target_positions=trajectory.state[:, :2],
    )
    metrics = relative_geometry_metrics(trajectory.state[:, :2], auxiliary)
    assert metrics["bearing_span_deg"] < 1e-10
    assert metrics["range_span_m"] > 0.0


def test_moving_geometry_is_better_conditioned_than_constant_bearing():
    trajectory = generate_curved_trajectory(dt=0.1, duration=20.0)
    constant_bearing = auxiliary_trajectory(
        trajectory.t,
        "constant_bearing",
        target_positions=trajectory.state[:, :2],
    )
    moving = auxiliary_trajectory(trajectory.t, "moving")
    degenerate = observability_history_planar_imu(
        trajectory.state,
        trajectory.ideal_imu,
        constant_bearing,
        trajectory.dt,
    )
    informative = observability_history_planar_imu(
        trajectory.state,
        trajectory.ideal_imu,
        moving,
        trajectory.dt,
    )
    assert informative["sigma_ratio"][-1] > degenerate["sigma_ratio"][-1]
