# Phase P1F-A — paper-state conventional PF results

## Core conclusion

The audited three-state Han-model baseline is internally consistent, but conventional bootstrap-PF global localization remains fragile even when the DR increments are exact and the only sensor error is Gaussian UWB range noise. This establishes a deliberately clean reference for the AACOPF mechanism: the challenge is not caused by the five-state IMU mechanization, unknown initial velocity, or a mismatch between simulated and estimator dynamics.

## Model-consistency check

The deterministic truth generated from

`state = [x, y, phi]`, `input = [Delta L, Delta phi]`

is reconstructed with maximum absolute numerical error **0.0** when propagated with the Equation-(3) pre-turn convention. This verifies the paper-state implementation before stochastic filtering is considered.

## Primary random-annulus baseline

Twenty stochastic runs use:

- 10,000 particles;
- random first-range annulus and random yaw initialization;
- `Delta d = 3 sigma_UWB = 0.36 m`;
- exact DR increments;
- Gaussian UWB noise with `sigma_UWB = 0.12 m`;
- informative moving auxiliary;
- Equation-(3) pre-turn propagation;
- systematic resampling for `N_eff < 0.5 N_p`.

Results:

| Metric | Result |
|---|---:|
| Terminal pose convergence | **4/20 = 20%** |
| Median convergence time among successful runs | **9.45 s** |
| Mean late position RMSE | 8.414 m |
| Mean late yaw RMSE | 22.615 deg |
| Mean final initial-yaw-lineage error | 21.177 deg |
| Mean normalized effective sample size | 0.821 |
| Mean resampling events | 14.05 |
| Mean runtime | 1.59 s/run |

The mean errors are dominated by wrong-mode runs and should not be interpreted as the accuracy conditional on successful localization. Successful runs have sub-metre/sub-degree terminal behaviour, while failed runs often become concentrated on an incorrect global mode.

## Initial-support diagnostic

Failure is **not simply caused by the random initializer missing the correct region**. In all 20 primary runs the initial 10,000-particle cloud contains particles satisfying the diagnostic initial region `position error < 1 m` and `|yaw error| < 10 deg`.

The initial correct-region count is:

- mean: **19.75 particles**;
- median: 18.5;
- minimum: 13;
- maximum: 27.

Thus the correct region occupies only about **0.20%** of the initial global cloud, but it is represented in every tested realization. Successful and failed runs have similar initial correct-region counts. This supports the interpretation that conventional likelihood updates and resampling can eliminate or fail to select a represented correct mode.

This diagnostic uses the same thresholds as the later convergence definition. It does not prove that every initially accepted particle would satisfy the terminal criterion if propagated forever; it establishes only that zero initial support is not an adequate explanation for the observed failures.

## Propagation-convention sensitivity

The paper contains the P1E-audited inconsistency between Equation (3) and Algorithm 1. Ten matched seeds were therefore evaluated with both conventions while the truth follows Equation (3).

| Convention | Pose success | Mean late position RMSE | Mean late yaw RMSE |
|---|---:|---:|---:|
| Equation (3), pre-turn | 1/10 | 10.698 m | 28.925 deg |
| Algorithm 1, post-turn | 1/10 | 14.344 m | 40.930 deg |

The success count is identical in this small matched set, but the post-turn interpretation increases mean late position RMSE by about 3.65 m and mean late yaw RMSE by about 12.0 deg. The Equation-(3) convention therefore remains the primary reproduction choice. The post-turn case is retained only as a documented source-ambiguity sensitivity.

## Structured-initialization control

A `100 x 100 = 10,000` structured position-bearing/yaw grid was evaluated for ten matched seeds. It achieves **0/10** terminal pose convergences, with mean late position RMSE 15.282 m and mean late yaw RMSE 47.16 deg.

This is not because the structured grid lacks the diagnostic correct region: every structured run contains 17--20 particles satisfying the initial 1 m / 10 deg region. Therefore deterministic angular coverage alone does not resolve the subsequent mode-selection problem. This control is not presented as the Han et al. initializer because the paper describes random particle generation.

## Particle-count study

Twenty stochastic runs were evaluated at each particle count.

| Particles | Pose success | Mean late position RMSE | Mean late yaw RMSE | Mean runtime |
|---:|---:|---:|---:|---:|
| 2,500 | 1/20 (5%) | 12.537 m | 35.993 deg | 0.43 s |
| 5,000 | 1/20 (5%) | 17.198 m | 53.563 deg | 0.79 s |
| 10,000 | 4/20 (20%) | 8.414 m | 22.615 deg | 1.64 s |
| 20,000 | 2/20 (10%) | 13.620 m | 40.191 deg | 3.23 s |
| 40,000 | 9/20 (45%) | 3.532 m | 9.507 deg | 6.53 s |

Every tested particle count has at least one initial correct-region particle in **20/20** realizations. The average correct-region count grows from about 4 particles at 2,500 particles to about 74 at 40,000 particles.

Increasing the population to 40,000 clearly improves the aggregate result, but the relationship is not monotonic at intermediate populations and even 40,000 particles fail in 11/20 runs. The particle populations are independent random realizations for each configured size rather than nested prefixes, so individual-seed success is not a paired monotonic quantity. The defensible conclusion is therefore that more particles improve representational density but **particle count alone does not make the conventional global PF reliable**.

The cost increase is approximately linear for the conventional PF: 40,000 particles take about four times the runtime of 10,000 particles. This contrasts with the literal all-pairs AACOPF transition, whose candidate scoring is quadratic in particle count.

## Interpretation for P1F-B

P1F-A separates several mechanisms that were confounded earlier in Phase 1:

1. **The paper-state dynamics are not the problem.** Exact-model reconstruction is exact.
2. **The five-state IMU adaptation is not required to obtain global-mode failures.** They remain in the simpler three-state paper model with exact DR increments.
3. **Initial support alone is not sufficient.** Correct-region particles are present in every tested primary run, but conventional PF success is only 20%.
4. **More particles help, but do not solve the problem.** Even 40,000 particles reach only 45% success under the present campaign.
5. **The propagation ambiguity matters quantitatively.** Equation (3) is the better-supported primary convention and the post-turn variant is a model-mismatch sensitivity.
6. **Structured angular coverage does not automatically preserve the correct mode.** The failure mechanism lies after initialization as well as in representation density.

This is the failure mode in which an ACO-inspired particle-management mechanism can be meaningfully tested. P1F-B should therefore first verify the literal ACO transition on a computationally small, controlled cloud and track correct-mode mass, diversity, moved-particle fraction, destination multiplicity, and runtime. Any claim that AACOPF improves localization must use exactly matched initial particles, ranges, DR increments, and random seeds against this conventional-PF reference.

## Reproducibility

- configuration: `configs/phase1f_a.yaml`
- implementation: `src/mt_ag/paper_pf.py`
- experiment: `experiments/run_phase1f_a.py`
- initial-support diagnostic: `experiments/analyze_phase1f_a_initial_support.py`
- tests: `tests/test_paper_pf.py`
- full raw results and support diagnostics: GitHub Actions artifact `p1f-a-results`, workflow run `33726102884`, artifact `9882181349`
- aggregate machine-readable result: `results/phase1/p1f_a/summary.json`

## Decision

P1F-A is complete. Retain **10,000 particles** as the practical conventional-PF reference for the later full global comparison because it is materially cheaper than 40,000 while preserving the difficult wrong-mode regime that AACOPF is intended to address. The 40,000-particle case remains an upper-population conventional reference, not the default.

Proceed to **P1F-B: literal-small AACOPF transition mechanism**. The purpose of P1F-B is mechanism verification, not yet a headline localization-performance comparison.
