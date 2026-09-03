import numpy as np

from mt_ag.geometry import wrap_angle
from mt_ag.initialization import initialize_range_conditioned_particles
from mt_ag.particle_filter import run_imu_bootstrap_pf_from_particles
from mt_ag.sensors import generate_uwb_ranges
from mt_ag.simulation import auxiliary_trajectory, generate_curved_trajectory


def test_range_conditioned_initializer_covers_ring_and_yaw():
    rng = np.random.default_rng(4)
    aux = np.array([8.0, -4.0])
    particles = initialize_range_conditioned_particles(
        5.0,
        aux,
        20,
        16,
        rng,
        radial_std_m=0.0,
        yaw_mode="uniform",
        velocity_mode="aligned_fixed_speed",
        speed_mean_mps=0.75,
    )
    assert particles.shape == (320, 5)
    radius = np.linalg.norm(particles[:, :2] - aux, axis=1)
    assert np.allclose(radius, 5.0)
    speed = np.linalg.norm(particles[:, 2:4], axis=1)
    assert np.allclose(speed, 0.75)
    velocity_heading = np.arctan2(particles[:, 3], particles[:, 2])
    assert np.max(np.abs(wrap_angle(velocity_heading - particles[:, 4]))) < 1e-12
    assert np.ptp(particles[:, 4]) > 5.0


def test_range_conditioned_pf_does_not_reuse_z0():
    trajectory = generate_curved_trajectory(dt=0.1, duration=0.5)
    auxiliary = auxiliary_trajectory(trajectory.t, "moving")
    ranges = generate_uwb_ranges(trajectory.state[:, :2], auxiliary)
    rng_init = np.random.default_rng(5)
    particles = initialize_range_conditioned_particles(
        ranges[0],
        auxiliary[0],
        10,
        8,
        rng_init,
        radial_std_m=0.0,
        yaw_mode="uniform",
        velocity_mode="aligned_fixed_speed",
        speed_mean_mps=0.75,
    )
    rng_pf = np.random.default_rng(6)
    result = run_imu_bootstrap_pf_from_particles(
        trajectory.ideal_imu,
        ranges,
        auxiliary,
        particles,
        trajectory.dt,
        rng_pf,
        sigma_process_accel_mps2=0.0,
        sigma_process_gyro_rps=0.0,
        sigma_uwb_m=0.12,
        initial_particles_conditioned_on_z0=True,
    )
    assert np.isclose(result.neff[0], len(particles))
    assert result.estimate.shape == trajectory.state.shape
    assert result.position_spread[0] > 0.0
