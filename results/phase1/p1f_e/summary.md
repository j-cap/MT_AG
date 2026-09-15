# Phase P1F-E — rank-deficient geometry negative controls

## Decision

P1F-E passes its negative-control purpose. With the P1F-C tuple frozen at

`(alpha, beta, c_lambda) = (0.5, 0, 2)`, AACOPF achieves **0/20 strict full-pose recoveries** in both the stationary and constant-bearing geometries. It therefore does not appear to manufacture a robust global-pose solution when the three-state finite-horizon range-observability diagnostic is rank deficient.

The experiment also sharpens the interpretation of the two geometries. The stationary single-auxiliary case has an exact global rotational gauge and is genuinely non-unique. The constant-bearing construction is only a local/first-order rank deficiency; conventional PF still recovers the true pose in 5/20 runs, so higher-order nonlinear range information cannot be ruled out.

## Three-state observability check

| Geometry | Final rank | Final sigma_min / sigma_max | Interpretation |
|---|---:|---:|---|
| Stationary auxiliary | 2/3 | 5.97e-17 | exact global rotational gauge |
| Constant bearing | 2/3 | ~1.42e-16 | local / first-order degeneracy |
| Informative moving | 3/3 | ~5.53e-2 | positive full-rank reference |

The stationary final singular values are approximately `(175.95, 10.13, 1.05e-14)`. A representative constant-bearing SVD is `(380.17, 10.28, 5.42e-14)`. Tiny variation in the last singular value across independent runners is ordinary SVD floating-point roundoff and does not alter the numerical rank.

## Matched filtering results

Each geometry uses 20 seeds and 1,200 particles. The random-annulus initialization contains 2.7 diagnostic correct-region particles on average; 19/20 seeds have at least one. Thus the experiment retains the sparse-support regime exposed by P1F-D rather than the deliberately well-populated P1F-C tuning clouds.

| Geometry | Method | Strict pose success | Mean late position RMSE | Mean late yaw RMSE | Median final position spread |
|---|---|---:|---:|---:|---:|
| Stationary | PF | 1/20 | 32.77 m | 85.46 deg | ~3.1e-13 m |
| Stationary | AACOPF | 0/20 | 37.08 m | 109.12 deg | 0.162 m |
| Constant bearing | PF | 5/20 | 6.59 m | 3.74 deg | ~1.9e-13 m |
| Constant bearing | AACOPF | 0/20 | 11.04 m | 13.26 deg | ~1.3e-13 m |

The one stationary PF success is not evidence against the symmetry argument: with a finite random cloud the filter can happen to select the true member of the rotationally equivalent family. It does so in only one seed and there is no information in the range sequence that makes that global orientation uniquely preferable.

## The stationary case exposes false confidence particularly clearly

For the 19 stationary PF runs that do **not** converge to the true pose:

- mean late position RMSE: **34.49 m**;
- mean late yaw RMSE: **89.94 deg**;
- median late MAP range-fit RMSE: only **0.052 m**.

At the same time, the median final PF position spread is essentially zero and the median yaw resultant is 1.0. Therefore the particle cloud can be extremely concentrated, and its selected trajectory can fit the UWB ranges substantially better than the measurement-noise scale, while the global pose remains grossly wrong. This is direct empirical evidence for

`measurement consistency / posterior concentration != global-pose correctness`.

## AACOPF behavior

The frozen ACO transition does not repair either negative control. Final diagnostic correct-mode mass/fraction is zero in all 20 AACOPF runs for both geometries. Aggressive genealogical behavior remains present:

| Geometry | Catastrophic ancestry collapse | Dominant-clone event |
|---|---:|---:|
| Stationary | 4/20 = 20% | 18/20 = 90% |
| Constant bearing | 4/20 = 20% | 17/20 = 85% |

Thus AACOPF can still become highly concentrated through copying without acquiring missing geometric information. This reinforces the P1F-B/P1F-D conclusion that particle concentration is not itself evidence of correct global localization.

## Constant-bearing nuance

Conventional PF succeeds in 5/20 constant-bearing runs, whereas frozen AACOPF succeeds in 0/20. This does **not** contradict the reported rank-two local observability matrix. The SVD is a first-order finite-horizon diagnostic, not a complete nonlinear observability theorem; the nonlinear range function may contain higher-order information that distinguishes some hypotheses over a long trajectory. Accordingly, the constant-bearing result should be described as locally degenerate/weakly informative, not as globally impossible.

The experiment therefore supports two distinct claims:

1. a fixed single auxiliary has a genuine global rotational ambiguity that no particle-management rule can remove without additional information;
2. a rank-deficient local linearization identifies a geometrically weak case, but nonlinear recovery must still be assessed empirically rather than ruled out categorically.

## Reproducibility

The 40 runs per method were executed as eight matched five-seed GitHub Actions shards from workflow run `34964754324`. All eight filtering shards completed successfully. The aggregate job in that run initially rejected tiny cross-runner differences in the last SVD singular value because it compared floating-point observability payloads for exact JSON equality; this was an aggregation-only issue, not an experiment failure. The aggregator has been corrected to recompute the deterministic observability diagnostics centrally.

The canonical sequential runner is `experiments/run_phase1f_e.py`; the parallel shard runner reproduces the same seed/geometry definitions while reducing wall-clock time.
