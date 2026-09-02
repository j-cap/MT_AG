# Phase P1C — Unknown initial pose

## Question

Can the IMU-driven bootstrap particle filter recover the target's global position and yaw when the initial position is only constrained by the first UWB range and the initial yaw is unknown?

P1C deliberately keeps the P1A/P1B target trajectory, moving auxiliary trajectory, IMU model, Gaussian UWB noise, and 10 Hz ranging schedule fixed. The only conceptual change is the initialization problem.

## First-range geometry

At sample 0 the known auxiliary position is `p_A,0` and the first measured UWB range is `z_0`. A single range measurement does not determine a unique target position; it constrains the target approximately to

`||p_0 - p_A,0|| = z_0`.

The initializer therefore distributes position hypotheses around this range ring. A radial perturbation with standard deviation equal to the UWB range-noise standard deviation represents uncertainty in `z_0`.

## Structured global pose coverage

Bearing hypotheses are placed uniformly around `[-pi, pi)`. For unknown yaw, yaw hypotheses are also placed uniformly around `[-pi, pi)`. The Cartesian product gives

`N_p = N_bearing * N_yaw`

particles. A small random common angular offset prevents the structured grid from being tied to a privileged angular origin while preserving uniform coverage.

For particle `(n,m)`, the position hypothesis is

`p_0^(n,m) = p_A,0 + r_(n,m) [cos(theta_n), sin(theta_n)]^T`,

with

`r_(n,m) = z_0 + epsilon_r`, `epsilon_r ~ N(0, sigma_r^2)`.

For the core unknown-yaw experiment,

`psi_m = -pi + 2*pi*m/N_yaw`.

## Important first-measurement convention

The initial particle cloud is already conditioned on `z_0`. The filter therefore does **not** apply the `z_0` likelihood again. The first subsequent Bayesian UWB update uses `z_1`. Applying `z_0` twice would artificially sharpen the initial radial distribution.

## Initial velocity question

Unlike the three-state Han formulation, the thesis IMU mechanization contains global velocity states. P1C therefore makes the initial velocity prior explicit.

Three cases are tested:

1. **Aligned fixed speed:** the initial speed magnitude is known (`0.75 m/s`) and velocity direction is aligned with the yaw hypothesis. This isolates unknown position/yaw and is the cleanest analogue of the Han formulation, in which travelled-distance increments are already available from the DR front-end.
2. **Aligned uncertain speed:** velocity remains aligned with yaw, but speed is sampled around `0.75 m/s` with `0.20 m/s` standard deviation.
3. **Free velocity:** speed is sampled over `[0, 1.5] m/s` and its direction is independent of yaw. This deliberately weak prior tests how important the kinematic alignment assumption is.

A known-yaw control experiment separates the range-ring position problem from the additional unknown-yaw problem.

## Convergence criteria

Whole-trajectory RMSE is not sufficient for a global-localization experiment because the estimator may spend an initial period resolving competing modes. P1C therefore uses terminal convergence:

- **position convergence:** position error is below `1.0 m` continuously from `t_conv` to the end of the experiment, with at least `5 s` remaining;
- **pose convergence:** the position condition holds and absolute yaw error is below `10 deg` continuously from `t_conv` to the end, again with at least `5 s` remaining;
- **late-window performance:** RMSE over the final `20 s` of the 60 s trajectory.

These are diagnostic research thresholds, not final application requirements.

## Experiment groups

- **P1C-A:** unknown position on the range ring, known yaw and known aligned speed;
- **P1C-B:** unknown position and yaw, known aligned speed;
- **P1C-C:** sensitivity to the initial velocity prior;
- **P1C-D:** structured grid resolution / particle-count sensitivity;
- **P1C-E:** resampling-threshold sensitivity to test whether premature resampling explains wrong-mode collapse;
- **P1C-F:** idealized identifiability sanity check with exact IMU/UWB information, 40,000 hypotheses, no process noise, no radial uncertainty, and no resampling.

## Results

### P1C-A — known-yaw control

With 1000 bearing hypotheses and known yaw, 8/10 runs reach terminal position and pose convergence. Mean late-window position RMSE is `0.733 m`, mean final position error is `0.811 m`, and mean late yaw RMSE is `0.132 deg`.

This establishes that the first-range ring itself is not the main difficulty. When yaw is known, the moving auxiliary and subsequent UWB sequence resolve the unknown position in most runs.

### P1C-B — core unknown-pose case

With `N_bearing=N_yaw=100`, i.e. 10,000 global pose hypotheses, only 2/20 runs (`10%`) reach terminal position and pose convergence. Mean whole-trajectory position RMSE is `10.728 m`, mean late-window RMSE is `12.616 m`, mean final position error is `8.784 m`, and mean late yaw RMSE is `29.810 deg`.

The final yaw resultant is nevertheless close to one on average (`0.997`), meaning that many failed runs become highly concentrated in yaw. The filter therefore frequently becomes **confident in a wrong global mode**, rather than merely remaining diffuse.

### P1C-C — velocity prior

Using the first ten matched seeds:

| Velocity prior | Pose success | Late position RMSE | Late yaw RMSE |
|---|---:|---:|---:|
| aligned fixed speed | 1/10 | 16.237 m | 41.302 deg |
| aligned uncertain speed | 0/10 | 16.393 m | 37.924 deg |
| free velocity | 0/10 | 35.449 m | 81.236 deg |

The free-velocity case is clearly the most difficult. Raw planar IMU acceleration and yaw rate do not determine the initial global velocity. P1C therefore confirms that an initial velocity condition, motion constraint, or additional source of speed information must be stated explicitly in the thesis.

### P1C-D — particle/grid resolution

| Angular grid | Particles | Pose success | Late position RMSE | Mean runtime |
|---:|---:|---:|---:|---:|
| 50 x 50 | 2,500 | 0/10 | 18.272 m | 0.439 s |
| 75 x 75 | 5,625 | 0/10 | 14.230 m | 0.906 s |
| 100 x 100 | 10,000 | 1/10 | 16.237 m | 1.558 s |
| 150 x 150 | 22,500 | 0/10 | 12.120 m | 3.409 s |

Higher resolution reduces some error metrics but does not make global convergence reliable. The problem is therefore not explained by insufficient particle count alone.

### P1C-E — resampling sensitivity

The resampling fraction was varied from no resampling (`gamma=0`) to the P1A/P1B value (`gamma=0.5`). Pose success remained between 0% and 20% across the tested values. With no resampling, mean late position RMSE was still `13.896 m`, and the final particle spread collapsed to only `0.0045 m` while the yaw resultant approached one.

Therefore premature conventional resampling is **not sufficient to explain** the failures. Even without resampling, likelihood accumulation can concentrate the posterior on an incorrect mode.

### P1C-F — idealized identifiability sanity check

The final diagnostic removes the main stochastic approximations:

- exact IMU inputs;
- exact geometric UWB ranges;
- zero initial radial uncertainty;
- zero PF process perturbation;
- no resampling;
- `N_bearing=N_yaw=200`, i.e. 40,000 global pose hypotheses;
- a small `0.02 m` likelihood width to evaluate exact-range consistency numerically.

All 3/3 runs converge to the correct terminal pose. The convergence times are 5.6–13.2 s. Aggregate results are:

| Metric | Idealized result |
|---|---:|
| Position convergence | 3/3 |
| Pose convergence | 3/3 |
| Mean whole-trajectory position RMSE | 1.815 m |
| Mean late-window position RMSE | 0.147 m |
| Mean final position error | 0.146 m |
| Mean late yaw RMSE | 0.266 deg |

This sanity check is intentionally not a statistical performance claim. Its role is to distinguish a fundamental model/geometry failure from a practical filtering failure. Under ideal information and sufficiently dense global coverage, the selected moving-auxiliary scenario resolves the initial position/yaw ambiguity.

## Interpretation

The P1C evidence supports a specific conclusion:

> The selected moving-auxiliary scenario is capable of resolving the initial global pose under ideal information, but the conventional bootstrap PF is not robust enough to preserve and discriminate the competing global pose hypotheses under the realistic P1C sensor uncertainty.

This interpretation is stronger than either of the simplistic alternatives. We should **not** conclude that the geometry is fundamentally unobservable, because the idealized check converges. We should also **not** conclude that more particles or merely disabling resampling solves the realistic problem, because the dedicated sweeps do not support that claim.

The sharp degradation from known yaw to unknown yaw identifies the global yaw ambiguity as the dominant difficulty in the current setup. The velocity-prior study further shows that the IMU-driven formulation introduces an initialization question that is absent from Han et al.'s higher-level `Delta L, Delta psi` interface.

## Decisions

1. **P1C is complete as a diagnostic phase.** It has established the boundary of the conventional bootstrap PF rather than producing a robust global-localization baseline.
2. **Do not tune P1C indefinitely.** Particle count and resampling threshold are not the main explanation for the failure.
3. **Keep the aligned fixed-speed prior for the next geometry study** so that P1D isolates observability/geometry rather than simultaneously weakening the velocity prior.
4. **Proceed to P1D geometry/observability.** The next task is to compare stationary, degenerate, and informative moving auxiliary trajectories and connect convergence behavior to observability diagnostics.
5. **Retain P1C as motivation for the AACOPF audit.** A central later question is whether the AACOPF transition mechanism preserves useful global modes better than the conventional PF under matched conditions.
6. **Treat the idealized test only as an identifiability sanity check.** It uses three runs, exact measurements, no process uncertainty, and a known aligned initial speed; it does not establish practical localization performance.
