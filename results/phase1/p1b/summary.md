# Phase P1B — Robustness of the known-pose baseline

## Question

Is the P1A result robust across stochastic realizations and reasonable uncertainty in the initial state, IMU quality, particle count, and PF process-noise assumptions?

The trajectory, moving auxiliary node, Gaussian UWB model, and 10 Hz ranging schedule remain unchanged from P1A so that P1B validates the existing estimator rather than mixing in Phase-2 communication effects.

## Main result

Across 20 independent stochastic seeds with the P1A configuration:

| Metric | Noisy IMU DR | PF + UWB |
|---|---:|---:|
| Position RMSE, mean ± std | 1.491 ± 0.635 m | **0.375 ± 0.196 m** |
| PF RMSE 95th percentile across seeds | — | 0.636 m |
| Worst PF RMSE across seeds | — | 1.066 m |
| Heading RMSE, mean | **0.102 deg** | 0.387 deg |

The PF has lower position RMSE than DR in **20/20 seeds**. Comparing aggregate mean RMSE values gives a position-RMSE reduction of approximately **74.9%**. As in P1A, this does not imply an improvement in yaw: the direct gyro integration remains better for heading in this known-pose test.

## Initial-state sensitivity

A random error is applied to the initial estimate and the PF prior is centered on that erroneous estimate. Values below are the standard deviations used for `[p_x,p_y,v_x,v_y,psi]`.

| Initial-error scale | PF RMSE mean | PF RMSE q95 | PF better than DR |
|---|---:|---:|---:|
| [0.25 m, 0.25 m, 0.06 m/s, 0.06 m/s, 1 deg] | 0.409 m | 0.627 m | 10/10 |
| [0.5 m, 0.5 m, 0.12 m/s, 0.12 m/s, 2 deg] | 0.668 m | 1.244 m | 10/10 |
| [1 m, 1 m, 0.25 m/s, 0.25 m/s, 5 deg] | 8.790 m | 35.963 m | 8/10 |
| [2 m, 2 m, 0.5 m/s, 0.5 m/s, 10 deg] | 14.215 m | 36.615 m | 8/10 |

The important result is the sharp change between the moderate and large-offset cases. The current local known-pose PF is robust to moderate initialization error but can converge to a wrong mode when the initial state becomes too uncertain. This is not treated as a defect to tune away in P1B; it motivates the explicit unknown-pose initialization problem in P1C.

## IMU-quality sensitivity

White-noise standard deviations and bias-random-walk intensities are jointly scaled relative to P1A.

| IMU scale | DR RMSE mean | PF RMSE mean | PF RMSE q95 |
|---:|---:|---:|---:|
| 0.5x | 0.736 m | 0.292 m | 0.479 m |
| 1x | 1.473 m | 0.415 m | 0.832 m |
| 2x | 2.946 m | 0.837 m | 2.090 m |
| 4x | 5.890 m | 3.741 m | 8.723 m |

UWB correction remains useful as IMU quality degrades, but it cannot fully compensate arbitrarily poor inertial propagation. At 4x the baseline IMU error the result becomes highly variable.

## Unestimated constant-bias stress test

The current PF does not estimate accelerometer or gyroscope biases. Synthetic constant initial biases were therefore added as a stress test. These values are **not hardware-calibrated** and should not be interpreted as expected DW3000/IMU performance.

| Bias case | Accel. bias [m/s²] | Gyro bias [rad/s] | DR RMSE | PF RMSE |
|---|---|---:|---:|---:|
| zero | [0, 0] | 0 | 1.473 m | 0.415 m |
| small | [0.002, -0.002] | 0.0002 | 1.862 m | 0.782 m |
| moderate | [0.01, -0.01] | 0.001 | 7.361 m | 4.551 m |
| large | [0.03, -0.03] | 0.003 | 21.937 m | 9.873 m |

The PF still improves position in all tested seeds, but persistent unmodeled bias produces large residual error. This is the strongest P1B evidence that bias handling may need to be revisited after real IMU characterization. We should not augment the state yet without knowing whether the hardware biases are large enough to justify the added state dimension.

## Particle-count / runtime trade-off

All particle-count configurations use the same 20 sensor realizations.

| Particles | PF RMSE mean | PF RMSE q95 | Mean runtime/run |
|---:|---:|---:|---:|
| 500 | 0.618 m | 1.913 m | 0.078 s |
| 1000 | 0.485 m | 0.821 m | 0.130 s |
| 2500 | 0.411 m | 0.799 m | 0.276 s |
| 5000 | 0.375 m | 0.636 m | 0.524 s |
| 10000 | 0.365 m | 0.680 m | 0.963 s |

There is clear diminishing return above 2500–5000 particles. Increasing from 5000 to 10000 almost doubles runtime while changing mean RMSE from 0.375 m to 0.365 m. The conservative decision is to retain **5000 particles** for the remainder of Phase 1 until the harder unknown-pose problem is understood. A 2500-particle configuration is reasonable for exploratory sweeps.

## PF process-noise sensitivity

The particle acceleration and yaw-rate perturbations are jointly scaled relative to the P1A values.

| PF process-noise scale | PF RMSE mean | PF RMSE q95 | PF better than DR |
|---:|---:|---:|---:|
| 0.25x | 1.009 m | 2.269 m | 7/10 |
| 0.5x | 0.656 m | 2.038 m | 9/10 |
| 1x | 0.415 m | 0.832 m | 10/10 |
| 2x | 0.420 m | 0.791 m | 10/10 |
| 4x | 0.512 m | 0.814 m | 10/10 |

The filter is particularly sensitive to process noise that is too small: under-dispersed particles cannot adequately follow the inertial uncertainty. The original P1A process noise lies in a broad useful region; 1x–2x gives similar position performance. There is therefore no evidence that P1A depended on a narrowly tuned process-noise value.

## Decisions

1. **P1A is statistically supported for position:** PF+UWB improves position RMSE over DR in all 20 baseline seeds.
2. **Do not claim yaw improvement:** the scalar UWB range does not improve yaw in this known-pose configuration.
3. **Keep 5000 particles as the conservative Phase-1 default.** The gain above this is negligible for P1B, but P1C may be more demanding.
4. **Keep the current PF process noise.** The result is robust around the selected value and clearly worse when the particle dynamics are under-dispersed.
5. **Do not add bias states yet.** Instead, record the strong bias sensitivity and use hardware characterization to decide whether explicit bias estimation is necessary.
6. **Proceed to P1C.** The catastrophic variability at large initial-state uncertainty is exactly the boundary that P1C should address through a principled unknown-pose initialization rather than further local tuning.
