from __future__ import annotations

import numpy as np

from .geometry import wrap_angle


def initialize_range_conditioned_particles(
    range_measurement_m,
    auxiliary_position,
    n_bearing,
    n_yaw,
    rng,
    *,
    radial_std_m=0.12,
    yaw_mode="uniform",
    known_yaw_rad=0.0,
    yaw_std_rad=0.0,
    velocity_mode="aligned_fixed_speed",
    speed_mean_mps=0.75,
    speed_std_mps=0.0,
    speed_min_mps=0.0,
    speed_max_mps=1.5,
):
    """Construct particles conditioned on the first UWB range measurement.

    The unknown position is represented by a ring around the known auxiliary
    node. Bearing hypotheses and, when requested, yaw hypotheses are laid out
    on uniform grids so that global angular coverage does not depend on a
    lucky random draw. Radial uncertainty is sampled from a Gaussian centered
    on the measured range.

    Parameters
    ----------
    range_measurement_m:
        First UWB range measurement z_0. The returned particle cloud is
        already conditioned on this measurement; a filter using these
        particles must therefore not apply z_0 a second time.
    auxiliary_position:
        Known [x, y] position of the auxiliary node at the first sample.
    n_bearing, n_yaw:
        Number of angular hypotheses around the range ring and in yaw.
    yaw_mode:
        ``uniform`` for unknown yaw or ``known`` for a narrow/known yaw.
    velocity_mode:
        ``aligned_fixed_speed``: velocity is aligned with each yaw and has
        speed ``speed_mean_mps``.
        ``aligned_gaussian_speed``: aligned with yaw with Gaussian speed.
        ``free_velocity``: speed and direction are sampled independently of
        yaw, representing a much weaker prior.
    """
    if n_bearing <= 0 or n_yaw <= 0:
        raise ValueError("n_bearing and n_yaw must be positive")
    if range_measurement_m <= 0.0:
        raise ValueError("range_measurement_m must be positive")
    if radial_std_m < 0.0:
        raise ValueError("radial_std_m must be non-negative")
    if speed_min_mps < 0.0 or speed_max_mps <= speed_min_mps:
        raise ValueError("invalid speed bounds")

    aux = np.asarray(auxiliary_position, dtype=float)
    if aux.shape != (2,):
        raise ValueError("auxiliary_position must have shape (2,)")

    bearing_offset = rng.uniform(-np.pi / n_bearing, np.pi / n_bearing)
    bearings = (
        np.linspace(-np.pi, np.pi, n_bearing, endpoint=False) + bearing_offset
    )

    if yaw_mode == "uniform":
        yaw_offset = rng.uniform(-np.pi / n_yaw, np.pi / n_yaw)
        yaws = np.linspace(-np.pi, np.pi, n_yaw, endpoint=False) + yaw_offset
    elif yaw_mode == "known":
        if n_yaw != 1:
            raise ValueError("yaw_mode='known' requires n_yaw=1")
        yaws = np.array([known_yaw_rad], dtype=float)
    elif yaw_mode == "gaussian":
        if n_yaw == 1:
            yaws = np.array([known_yaw_rad], dtype=float)
        else:
            quantiles = (np.arange(n_yaw) + 0.5) / n_yaw
            # A symmetric deterministic approximation avoids a scipy
            # dependency. The exact normal quantiles are not important here;
            # this mode is intended only for narrow known-yaw priors.
            centered = 2.0 * quantiles - 1.0
            yaws = known_yaw_rad + 2.5 * yaw_std_rad * centered
    else:
        raise ValueError(f"unknown yaw_mode: {yaw_mode}")

    bearing_grid, yaw_grid = np.meshgrid(bearings, yaws, indexing="ij")
    bearing_flat = bearing_grid.ravel()
    yaw_flat = wrap_angle(yaw_grid.ravel())
    n_particles = bearing_flat.size

    if radial_std_m == 0.0:
        radius = np.full(n_particles, float(range_measurement_m))
    else:
        radius = range_measurement_m + rng.normal(0.0, radial_std_m, n_particles)
        radius = np.maximum(radius, 1e-6)

    particles = np.zeros((n_particles, 5), dtype=float)
    particles[:, 0] = aux[0] + radius * np.cos(bearing_flat)
    particles[:, 1] = aux[1] + radius * np.sin(bearing_flat)
    particles[:, 4] = yaw_flat

    if velocity_mode == "aligned_fixed_speed":
        speed = np.full(n_particles, float(speed_mean_mps))
        velocity_heading = yaw_flat
    elif velocity_mode == "aligned_gaussian_speed":
        speed = rng.normal(speed_mean_mps, speed_std_mps, n_particles)
        speed = np.clip(speed, speed_min_mps, speed_max_mps)
        velocity_heading = yaw_flat
    elif velocity_mode == "free_velocity":
        speed = rng.uniform(speed_min_mps, speed_max_mps, n_particles)
        velocity_heading = rng.uniform(-np.pi, np.pi, n_particles)
    else:
        raise ValueError(f"unknown velocity_mode: {velocity_mode}")

    particles[:, 2] = speed * np.cos(velocity_heading)
    particles[:, 3] = speed * np.sin(velocity_heading)
    return particles
