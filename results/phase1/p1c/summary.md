# Phase P1C — Aggregate results

## Core conclusion

The realistic unknown-pose bootstrap PF is not robust: with 10,000 globally distributed position/yaw hypotheses, only **2/20** runs reach terminal position and pose convergence. An idealized 40,000-hypothesis sanity check with exact IMU/UWB information converges in **3/3** runs. This indicates that the selected moving-auxiliary scenario can resolve the global pose ambiguity in principle, while the conventional PF struggles to preserve/discriminate the competing modes under realistic uncertainty.

## Main experiment

| Case | Particles | Pose success | Late position RMSE | Final position error | Late yaw RMSE |
|---|---:|---:|---:|---:|---:|
| Known-yaw control | 1,000 | 8/10 | 0.733 m | 0.811 m | 0.132 deg |
| Unknown pose, aligned fixed speed | 10,000 | 2/20 | 12.616 m | 8.784 m | 29.810 deg |
| Idealized unknown pose | 40,000 | 3/3 | 0.147 m | 0.146 m | 0.266 deg |

## Initial-velocity prior

| Prior | Pose success | Late position RMSE | Late yaw RMSE |
|---|---:|---:|---:|
| Aligned fixed speed | 1/10 | 16.237 m | 41.302 deg |
| Aligned uncertain speed | 0/10 | 16.393 m | 37.924 deg |
| Free velocity | 0/10 | 35.449 m | 81.236 deg |

The free-velocity case is particularly poor, confirming that raw accelerometer/yaw-rate measurements do not provide the initial global velocity. The initial velocity assumption must remain explicit.

## Grid / particle-count sensitivity

| Grid | Particles | Pose success | Late position RMSE | Runtime |
|---|---:|---:|---:|---:|
| 50 x 50 | 2,500 | 0/10 | 18.272 m | 0.439 s |
| 75 x 75 | 5,625 | 0/10 | 14.230 m | 0.906 s |
| 100 x 100 | 10,000 | 1/10 | 16.237 m | 1.558 s |
| 150 x 150 | 22,500 | 0/10 | 12.120 m | 3.409 s |

More particles reduce some errors but do not make convergence reliable.

## Resampling sensitivity

| Resampling fraction | Pose success | Late position RMSE | Final particle spread |
|---:|---:|---:|---:|
| 0.00 | 0/10 | 13.896 m | 0.0045 m |
| 0.05 | 2/10 | 17.032 m | 1.716 m |
| 0.10 | 1/10 | 16.203 m | 1.982 m |
| 0.25 | 2/10 | 17.961 m | 2.461 m |
| 0.50 | 1/10 | 16.237 m | 1.563 m |

Disabling resampling does not solve the problem. The no-resampling case can become extremely concentrated while remaining wrong, so wrong-mode collapse cannot be attributed only to conventional resampling.

## Idealized sanity check

Design: 3 seeds, 200 bearing hypotheses, 200 yaw hypotheses, 40,000 particles, exact IMU, exact UWB, zero radial uncertainty, zero process noise, no resampling, and a 0.02 m likelihood width.

| Seed | Pose convergence | Late position RMSE | Final position error | Late yaw RMSE |
|---:|---:|---:|---:|---:|
| 0 | 11.7 s | 0.053 m | 0.042 m | 0.141 deg |
| 1 | 5.8 s | 0.204 m | 0.207 m | 0.367 deg |
| 2 | 13.2 s | 0.183 m | 0.189 m | 0.290 deg |

The test is an identifiability diagnostic, not a practical performance claim.

## Decision

P1C is complete. Proceed to P1D with the aligned fixed-speed prior and explicitly study the relationship between auxiliary/target geometry, local observability, particle-cloud evolution, and convergence. P1C should also be retained as motivation for the later AACOPF audit and comparison.
