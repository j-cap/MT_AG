# Phase P1D — Geometry and observability results

## Core conclusion

The three controlled auxiliary geometries have different observable state dimensions. The stationary single auxiliary remains rank 4, the synthetic constant-bearing case rank 3, and only the informative moving auxiliary reaches full rank 5 for the five-state IMU-driven model. This hierarchy is reflected in global-pose localization: the two rank-deficient cases never converge to the full pose, while the moving case is the only geometry that succeeds. Full local observability is nevertheless not sufficient for robust conventional bootstrap-PF localization under realistic uncertainty.

## Geometry and observability

| Geometry | LOS bearing span | Rank | Final sigma_min / sigma_max | First ratio >= 1e-4 |
|---|---:|---:|---:|---:|
| Stationary | 81.99 deg | 4 | 3.716e-17 | never |
| Constant bearing | ~0 deg | 3 | 1.740e-17 | never |
| Moving | 102.59 deg | 5 | 6.569e-4 | 26.6 s |

Final singular values:

- stationary: `[931.969, 232.845, 49.461, 7.046, 3.463e-14]`
- constant bearing: `[861.384, 155.442, 10.683, 8.158e-14, 1.499e-14]`
- moving: `[916.485, 320.809, 45.580, 8.047, 0.602]`

The stationary case demonstrates that line-of-sight bearing change alone is not sufficient: a fixed single anchor still leaves the global rotation symmetry unresolved when initial global yaw is unknown.

## Realistic filtering — 20 matched runs per geometry

| Geometry | Pose success | Late position RMSE | Final position error | Late yaw RMSE | Final position spread |
|---|---:|---:|---:|---:|---:|
| Stationary | 0/20 | 33.874 m | 34.781 m | 85.574 deg | 31.717 m |
| Constant bearing | 0/20 | 3.975 m | 4.510 m | 1.992 deg | 4.594 m |
| Moving | 5/20 | 10.507 m | 9.148 m | 23.838 deg | 1.050 m |

For the five successful moving runs, median terminal pose convergence is 49.8 s.

## Idealized filtering — 3 runs per geometry

Exact IMU and UWB data, zero radial/process uncertainty, no resampling, `200 x 200 = 40,000` global hypotheses.

| Geometry | Pose success | Late position RMSE | Final position error | Late yaw RMSE |
|---|---:|---:|---:|---:|
| Stationary | 0/3 | 33.151 m | 33.838 m | 91.595 deg |
| Constant bearing | 0/3 | 4.897 m | 7.821 m | 3.204 deg |
| Moving | 2/3 | 2.893 m | 5.180 m | 8.951 deg |

For the successful moving idealized runs:
- seed 0: pose convergence 12.0 s, late position RMSE 0.069 m, final error 0.053 m, late yaw RMSE 0.184 deg;
- seed 2: pose convergence 5.7 s, late position RMSE 0.190 m, final error 0.192 m, late yaw RMSE 0.341 deg.

The failed idealized moving run shows that a finite point-particle global grid remains discretization-sensitive even when the geometry is informative.

## Interpretation

1. Geometry is a prerequisite, not merely an accuracy modifier. Rank-deficient geometries do not recover the full global pose even with ideal measurements.
2. The stationary case retains one global rotational degree of freedom despite substantial apparent bearing change caused by target motion.
3. Constant-bearing motion can recover yaw well while leaving position unresolved; sub-state accuracy and posterior concentration are not sufficient evidence of full-state observability.
4. The moving case becomes full rank and is the only case with successful global-pose recovery, but realistic bootstrap-PF success remains only 25%.
5. Therefore AACOPF or any other particle-management method should be evaluated only after distinguishing geometric non-observability from posterior-approximation failure. An algorithm cannot create information absent from the geometry.

## Reproducibility

- configuration: `configs/phase1d.yaml`
- experiment: `experiments/run_phase1d.py`
- observability implementation: `src/mt_ag/observability.py`
- geometry implementation: `src/mt_ag/simulation.py`
- detailed scientific note: `docs/experiments/phase1d_geometry_observability.md`
- full time histories and particle snapshots: GitHub Actions artifact `p1d-results` from workflow run `33715346727`

## Decision

P1D is complete. Proceed to P1E: audit the Han et al. AACOPF equations, transition mechanism, thresholds, and underspecified implementation choices before implementing AACOPF.
