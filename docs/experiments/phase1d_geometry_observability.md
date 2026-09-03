# Phase P1D — Geometry and observability

## Question

How does the relative geometry between the target and one globally known auxiliary node affect local observability, global-pose convergence, and particle-cloud evolution in the P1C unknown-pose problem?

P1D keeps the target trajectory, IMU/UWB noise models, 10 Hz range schedule, global first-range initialization, aligned fixed initial speed, and bootstrap PF fixed. Only the auxiliary trajectory is changed.

## Hypothesis

A UWB range is not equally informative under all relative motions. Geometry that removes the symmetry of the one-range problem should increase the rank and conditioning of the finite-horizon observability matrix and make global position/yaw recovery possible. Stationary or deliberately degenerate relative motion should retain unobservable directions.

This follows the qualitative observability argument in Han et al. (Sensors 2020, 20(2):467), but the thesis evaluates the five-state IMU-driven model `[p_x,p_y,v_x,v_y,psi]` rather than the paper's three-state `Delta L, Delta psi` interface.

## Controlled geometries

All auxiliary trajectories start from `[8,-4] m`, so the first true UWB range is identical.

1. **Stationary:** auxiliary fixed at `[8,-4] m` in the navigation frame.
2. **Constant bearing:** auxiliary motion is constructed from the target truth so that the navigation-frame target-to-auxiliary line-of-sight direction remains constant while the range increases. This is a synthetic degeneracy test, not a proposed physical controller.
3. **Moving:** the informative moving trajectory already used in P1A--P1C.

The constant-bearing construction separates changing scalar range from changing measurement direction. Range changes by 4.8 m, but the line-of-sight bearing changes only at numerical round-off level.

## Observability diagnostic

For the linearized five-state IMU propagation and scalar UWB range observation, the cumulative finite-horizon observability matrix `O_0:k` is evaluated along the true trajectory. The implementation computes the singular values directly from `O_0:k`; it deliberately does not infer them from `O^T O`, because forming the Gram matrix squares the condition number and can create misleading numerical rank in nearly singular cases.

P1D records

`eta_O(k) = sigma_min(O_0:k) / sigma_max(O_0:k)`.

`eta_O` is a scale-normalized local conditioning diagnostic. It is not a complete nonlinear observability proof and the diagnostic threshold `eta_O >= 1e-4` has no application-level meaning. It is used only to compare how quickly useful conditioning accumulates across the controlled cases.

The study also records range span, line-of-sight bearing span, and total bearing variation.

## Filtering experiments

### Realistic layer

For each geometry, 20 matched stochastic runs use the P1C sensor/PF uncertainty and a `100 x 100 = 10,000` particle global position/yaw grid. Sensor-noise samples, initialization randomness, and PF process-noise samples are matched across geometries for each seed.

### Idealized layer

For each geometry, three runs use `200 x 200 = 40,000` particles, exact IMU, exact UWB, zero radial uncertainty, zero PF process perturbation, and no resampling. This layer distinguishes geometric ambiguity from realistic sensor/filtering difficulty.

Terminal convergence uses the P1C rule: position error below 1 m and absolute yaw error below 10 deg continuously until the end, with at least 5 s remaining. Late-window RMSE uses the final 20 s.

## Result 1 — the three geometries have different observable dimensions

| Geometry | Bearing span | Final singular values of `O` | Rank | Final `eta_O` |
|---|---:|---|---:|---:|
| Stationary | 81.99 deg | `[931.97, 232.84, 49.46, 7.046, 3.46e-14]` | **4** | `3.72e-17` |
| Constant bearing | ~0 deg | `[861.38, 155.44, 10.68, 8.16e-14, 1.50e-14]` | **3** | `1.74e-17` |
| Moving | 102.59 deg | `[916.48, 320.81, 45.58, 8.047, 0.602]` | **5** | `6.57e-4` |

Only the moving auxiliary produces a full-rank five-state local observability matrix. Its conditioning ratio first reaches the diagnostic value `1e-4` at approximately **26.6 s**.

The stationary result is especially important. Although the target-to-anchor bearing changes by about 82 deg as the target moves, one state direction remains unobservable. A fixed single range reference retains a global rotational symmetry when the target's initial global yaw is unknown. Hence raw bearing variation alone is not a sufficient observability criterion.

The constant-bearing construction is even more degenerate: only three independent local state directions are retained. This confirms that changing range values alone are not enough when the range Jacobian repeatedly points in essentially the same direction.

## Result 2 — filtering follows the observability hierarchy

| Geometry | Realistic pose success | Realistic late pos. RMSE | Realistic late yaw RMSE | Ideal pose success |
|---|---:|---:|---:|---:|
| Stationary | **0/20** | 33.874 m | 85.574 deg | **0/3** |
| Constant bearing | **0/20** | 3.975 m | 1.992 deg | **0/3** |
| Moving | **5/20** | 10.507 m | 23.838 deg | **2/3** |

The two rank-deficient geometries never achieve full terminal pose convergence, even in the idealized layer. The full-rank moving geometry is the only case that succeeds at all.

This is the central P1D result: **geometry is a prerequisite for global pose recovery rather than merely a secondary accuracy factor.**

## Result 3 — observability is necessary, but not sufficient for robust bootstrap-PF localization

The moving case is locally full rank, yet the realistic bootstrap PF succeeds in only 5/20 runs. Among successful runs, median terminal pose convergence is approximately **49.8 s**. The deterministic conditioning ratio reaches `1e-4` already at 26.6 s, so useful geometry accumulates well before robust stochastic convergence.

The idealized moving case succeeds in 2/3 runs, with successful convergence at approximately 12.0 s and 5.7 s. The remaining idealized run converges to an incorrect discrete hypothesis. Thus even a dense 40,000-point global particle grid can remain sensitive to angular discretization and likelihood concentration.

The resulting hierarchy is therefore

`geometric observability -> possibility of correct global recovery`,

but not

`geometric observability -> guaranteed bootstrap-PF convergence`.

This distinction is important for the later AACOPF evaluation: a particle-management method may improve approximation of an informative multimodal posterior, but it cannot create information in a geometry that is structurally unobservable.

## Result 4 — sub-state confidence can be misleading

The constant-bearing case provides a useful counterexample. In realistic runs the late yaw RMSE is only about **1.99 deg** and the mean final yaw resultant is about **0.996**, yet full pose convergence is 0/20 and mean late position RMSE remains about **3.97 m**. In the idealized runs, the point-particle approximation can collapse almost completely to one hypothesis even though the geometry is rank deficient.

Therefore:

- accurate yaw does not imply globally correct position;
- a concentrated particle cloud does not prove observability;
- posterior concentration should not be used as a standalone correctness metric.

The stationary case shows the opposite failure mode: the unresolved global rotational symmetry remains visibly diffuse, with about 31.7 m mean final position spread and a mean yaw resultant of only 0.314 in realistic runs.

## Particle-cloud diagnostics

Representative weighted particle snapshots at 0, 10, 30, and 60 s for each geometry were generated by the experiment workflow. The full snapshot set and the time-resolved diagnostic CSV are retained as GitHub Actions artifacts rather than committed to normal Git history because they are multi-megabyte diagnostic outputs. The aggregate conclusions above use all matched runs, not the representative snapshot alone.

## Numerical and methodological limitations

- The SVD analysis is a local linearization along the true trajectory; it is not a global nonlinear observability proof.
- The constant-bearing trajectory is deliberately constructed using ground truth and exists only as a diagnostic degeneracy.
- The moving case uses one specific informative trajectory, so P1D establishes a controlled existence/comparison result rather than a universal geometry threshold.
- The idealized PF still uses a finite structured angular grid; one of three moving runs fails because point-particle coverage and likelihood concentration remain discretization-sensitive.
- The realistic PF result therefore mixes observability, sensor uncertainty, particle approximation, and process-model mismatch after geometry has made the state locally observable.

## Decision

P1D is complete. Its main conclusion is:

> The relative UWB geometry determines which global state directions are observable. A fixed single auxiliary leaves one rotational degree of freedom, constant-bearing relative motion is more degenerate, and the informative moving auxiliary is the only tested geometry that yields full five-state local observability and any global-pose convergence. Full observability is nevertheless not sufficient for robust conventional bootstrap-PF localization under realistic uncertainty.

Proceed to **P1E — the Han et al. AACOPF reproduction audit**. The audit must explicitly separate two questions:

1. which failures are caused by insufficient geometry and therefore cannot be repaired algorithmically; and
2. whether AACOPF's particle-transition mechanism improves global-mode preservation and convergence in geometries that are already informative.
