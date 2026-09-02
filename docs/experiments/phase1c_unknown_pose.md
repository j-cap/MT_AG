# Phase P1C — Unknown initial pose

## Question

Can the IMU-driven bootstrap particle filter recover the target's global position and yaw when the initial position is only constrained by the first UWB range and the initial yaw is unknown?

P1C deliberately keeps the P1A/P1B trajectory, moving auxiliary trajectory, IMU model, Gaussian UWB noise, and 10 Hz ranging schedule fixed. The only conceptual change is the initialization problem.

## First-range geometry

At sample 0 the known auxiliary position is `p_A,0` and the first measured UWB range is `z_0`. With one range measurement, the target position is not a point. It is approximately restricted to the circle

`||p_0 - p_A,0|| = z_0`.

The initializer therefore distributes position hypotheses around this range ring. A radial perturbation with standard deviation equal to the UWB range-noise standard deviation represents uncertainty in `z_0`.

## Structured global pose coverage

Bearing hypotheses are placed uniformly around `[-pi, pi)`. For the unknown-yaw experiment, yaw hypotheses are also placed uniformly around `[-pi, pi)`. The Cartesian product gives

`N_p = N_bearing * N_yaw`

particles. This structured grid avoids accidentally missing an angular region because of random initialization and mirrors the basic position-heading hypothesis structure used by Han et al.

## Important first-measurement convention

The initial particle cloud is already conditioned on `z_0`. The filter therefore does **not** apply the `z_0` likelihood again. The first subsequent Bayesian UWB update uses `z_1`. Applying `z_0` twice would artificially sharpen the initial radial distribution.

## Initial velocity question

Unlike the three-state Han formulation, the thesis IMU mechanization contains global velocity states. P1C therefore makes the initial velocity prior explicit.

Three cases are tested:

1. **Aligned fixed speed:** the initial speed magnitude is known (`0.75 m/s`) and velocity direction is aligned with the yaw hypothesis. This isolates unknown position/yaw and is the cleanest analogue of the Han formulation, in which travelled-distance increments are already available from the DR front-end.
2. **Aligned uncertain speed:** velocity remains aligned with yaw, but speed is sampled around `0.75 m/s` with `0.20 m/s` standard deviation.
3. **Free velocity:** speed is sampled over `[0, 1.5] m/s` and its direction is independent of the yaw hypothesis. This is a deliberately weak prior and tests how important the kinematic alignment assumption is.

A known-yaw control experiment is also included to separate the range-ring position problem from the additional unknown-yaw problem.

## Convergence criteria

P1C is a global convergence experiment, so whole-trajectory RMSE alone is insufficient. The following diagnostic convergence criteria are used:

- position convergence: position error below `1.0 m` for at least `5 s` continuously;
- pose convergence: position error below `1.0 m` and absolute yaw error below `10 deg` for at least `5 s` continuously;
- late-window performance: RMSE over the final `20 s` of the trajectory.

These thresholds are diagnostic research thresholds, not final application requirements.

## Experiment groups

- **P1C-A:** unknown position on the range ring, known yaw and known aligned speed;
- **P1C-B:** unknown position and yaw, known aligned speed;
- **P1C-C:** sensitivity to the initial velocity prior;
- **P1C-D:** structured grid resolution / particle-count sensitivity.

## Expected interpretation

A successful P1C result would show that the moving auxiliary and target motion can collapse a globally distributed position/yaw hypothesis set to the correct mode when the velocity prior is sufficiently informative. If the free-velocity case fails while the aligned-speed cases succeed, the thesis must explicitly acknowledge that raw planar IMU measurements do not by themselves provide initial velocity and that an initialization condition, motion constraint, or additional source of speed information is required.
