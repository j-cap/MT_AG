# Phase P1F-E — rank-deficient geometry negative controls

## Purpose

P1F-E tests a deliberately negative hypothesis: the frozen AACOPF transition must **not** be interpreted as creating information that is absent or only weakly present in the ranging geometry. The same repository-selected tuple from P1F-C,

`(alpha, beta, c_lambda) = (0.5, 0, 2)`,

is applied without retuning to the stationary and constant-bearing geometries.

The experiment is not designed to improve localization performance. Its purpose is to separate algorithmic particle management from geometric information content.

## Important model-specific check

P1D established rank deficiency for the five-state IMU-driven model. P1F-E uses the separate paper-level state `[x,y,phi]`, so the P1D rank result is **not transferred by assumption**. P1F-E therefore adds a three-state finite-horizon local range-observability diagnostic using the exact paper-level DR propagation Jacobian.

For every sample, the cumulative range Jacobian is evaluated by direct SVD. The stationary and constant-bearing cases must remain locally rank deficient for `[x,y,phi]`; the informative moving geometry is computed as a positive observability reference but is not rerun as a filter experiment because P1F-D already provides the matched moving-geometry comparison.

The two negative controls do not have identical theoretical status. With a stationary single auxiliary, a global rotation of the entire target trajectory and yaw about the fixed auxiliary leaves all ranges and DR increments invariant, giving an exact rotational gauge freedom. The constant-bearing construction is a first-order/local degeneracy: the range Jacobian lacks sensitivity in one state direction, but higher-order nonlinear range information may still exist. P1F-E therefore treats its SVD rank as a local conditioning diagnostic rather than a complete nonlinear impossibility result.

## Filter design

- deterministic 60 s paper-state trajectory;
- exact `[Delta L, Delta phi]` increments;
- Equation-(3) `pre_turn` propagation;
- Gaussian UWB noise, `sigma=0.12 m`;
- random first-range annulus + random yaw initialization;
- 1,200 particles, matching the primary P1F-D budget;
- 20 seeds;
- conventional PF and frozen literal AACOPF evaluated from the same initial cloud;
- identical additive UWB-noise realization across the two negative-control geometries for each seed;
- same initial cloud across both geometries because the initial auxiliary position/range is identical.

## Primary acceptance criterion

Neither conventional PF nor frozen AACOPF should show **systematic** full-pose recovery in either locally rank-deficient geometry. For the stationary case, systematic unique global-pose recovery would contradict the exact rotational symmetry and would trigger an implementation/metric audit. For constant-bearing motion, occasional recovery is not mathematically excluded by the local-rank test, but robust recovery comparable to informative moving geometry would indicate that the first-order diagnostic alone is not capturing the relevant nonlinear information.

## Diagnostics

For each geometry and method P1F-E records:

- strict terminal full-pose convergence;
- late position and yaw RMSE;
- final position spread and yaw resultant;
- final correct-mode mass/fraction;
- MAP-trajectory range-fit RMSE;
- conventional-PF resampling diversity;
- AACOPF unique-parent fraction, catastrophic collapse, dominant cloning and moved fraction;
- runtime.

The MAP range-fit diagnostic is important because a rank-deficient estimator can fit the available UWB measurements while remaining globally wrong. Low range residual is therefore not sufficient evidence for global-pose correctness.

## Claim boundary

P1F-E is an empirical negative control backed by a local finite-horizon rank diagnostic. It is not a general nonlinear observability theorem. The stationary case additionally has an exact rotational gauge symmetry; the constant-bearing case is interpreted more narrowly as local/first-order degeneracy.
