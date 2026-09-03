# Phase P1D — Geometry and observability

## Question

How does the relative geometry between the target and one globally known auxiliary node affect local observability, global-pose convergence, and particle-cloud evolution in the P1C unknown-pose problem?

P1D keeps the target trajectory, IMU/UWB noise models, 10 Hz range schedule, global first-range initialization, aligned fixed initial speed, and bootstrap PF fixed. Only the auxiliary trajectory is changed.

## Hypothesis

A UWB range measurement is not equally informative under all relative motions. Geometry that changes the target-to-auxiliary line-of-sight direction should improve the conditioning of the finite-horizon observability matrix and make the initial position/yaw ambiguity easier to resolve. Stationary or deliberately constant-bearing geometries should remain poorly conditioned.

This follows the qualitative observability argument in Han et al. (Sensors 2020, 20(2):467), but the thesis evaluates the five-state IMU-driven model `[p_x,p_y,v_x,v_y,psi]` rather than copying the paper's three-state `Delta L, Delta psi` model.

## Controlled geometries

All auxiliary trajectories start from `[8,-4] m`, so the initial true range is identical.

1. **Stationary:** auxiliary fixed at `[8,-4] m` in the navigation frame.
2. **Constant bearing:** auxiliary motion is constructed from the target truth so that the navigation-frame line-of-sight direction is constant while the range increases. This is a synthetic degeneracy test, not a proposed physical controller.
3. **Moving:** the informative moving trajectory already used in P1A--P1C.

The constant-bearing construction isolates the information carried by changing relative direction: range values continue to change, but the line-of-sight direction does not.

## Observability diagnostic

For the linearized five-state IMU propagation and scalar UWB range observation, the cumulative finite-horizon observability matrix is evaluated along the true trajectory. P1D records its singular values and

`eta_O(k) = sigma_min(O_0:k) / sigma_max(O_0:k)`.

`eta_O` is a scale-normalized conditioning diagnostic. It is not a complete nonlinear observability proof and the configured `1e-3` threshold is only used to report a comparable diagnostic time.

The study also records relative-range span, relative-bearing span, and total relative-bearing variation.

## Filtering experiments

Two layers are used.

### Realistic layer

For each geometry, 20 matched stochastic runs use the P1C sensor and PF uncertainty with a `100 x 100 = 10,000` particle global position/yaw grid. Sensor-noise samples, initialization randomness, and PF process-noise samples are matched across geometries for each seed.

Metrics include terminal position/pose convergence, late-window position and yaw RMSE, final error, effective sample size, position spread, yaw resultant, and runtime. The time histories of mean position error, position spread, yaw resultant, and normalized effective sample size are stored together with the observability history.

### Idealized layer

For each geometry, three runs use `200 x 200 = 40,000` particles, exact IMU, exact UWB, zero radial uncertainty, zero PF process perturbation, and no resampling. This layer asks whether the geometry itself permits the global hypothesis set to collapse under ideal information.

## Particle-cloud diagnostics

For one fixed representative seed, weighted particle snapshots are stored at 0, 10, 30, and 60 s for each geometry. These snapshots are diagnostic only; aggregate conclusions are based on all matched runs.

## Decision rule

P1D is complete when the geometry cases show a clear, documented relationship between relative motion, finite-horizon observability conditioning, and localization behavior, or when the experiment demonstrates that the chosen observability diagnostic does not explain the observed PF behavior. Either outcome is scientifically useful and determines how the AACOPF reproduction should be interpreted.
