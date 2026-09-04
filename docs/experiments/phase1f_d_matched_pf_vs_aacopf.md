# Phase P1F-D — Matched conventional PF versus frozen literal AACOPF

## Purpose

P1F-D is the first direct performance comparison between the conventional paper-state PF from P1F-A and the audited literal-small AACOPF after its parameters were selected and frozen in P1F-C.

The research question is:

> Under an informative moving-auxiliary geometry and identical global initial hypotheses, does the frozen AACOPF improve global-mode selection and pose convergence relative to conventional systematic-resampling PF at the same particle budget?

The frozen AACOPF tuple is

`alpha = 0.5`, `beta = 0`, `c_lambda = 2`.

These values are repository-selected P1F-C values and are not attributed to Han et al.

## Claim boundary

The literal all-pairs ACO transition is computationally quadratic. P1F-B showed that the paper-implied population of 160,000 particles is infeasible for this literal interpretation. P1F-D therefore does **not** compare AACOPF at the 10,000-particle practical PF reference used in P1F-A.

Instead, P1F-D performs equal-budget matched comparisons at manageable populations:

- 400 particles;
- 800 particles;
- 1,200 particles.

The 1,200-particle case is the primary literal-small comparison. P1F-A's 10,000- and 40,000-particle conventional-PF results remain contextual upper-population references, not matched AACOPF comparisons.

## Matched experiment design

For each particle budget and random seed, the two methods receive exactly the same:

- 60 s paper-level truth trajectory;
- Equation-(3) pre-turn propagation convention;
- exact `Delta L` and `Delta phi` increments;
- informative moving globally known auxiliary trajectory;
- Gaussian UWB realization with `sigma_UWB = 0.12 m`;
- first UWB measurement `z0`;
- random first-range annulus particle cloud;
- random initial yaw hypotheses.

The first-range annulus half-width remains the P1E/P1F repository choice

`Delta d = 3 sigma_UWB = 0.36 m`.

The first measurement is used only to construct the initial cloud; likelihood updates begin at the next sample to avoid double-counting `z0`. Twenty unseen seeds are evaluated at each particle budget. These random streams are separate from the P1F-C development and validation streams.

## Methods

### Conventional PF

The conventional method is the P1F-A paper-state bootstrap PF with systematic resampling triggered when

`N_eff < 0.5 Np`.

There is no artificial process noise in `Delta L` or `Delta phi` in this experiment, so both methods receive the same exact DR propagation. The purpose is to compare particle replacement/mode management rather than process-noise handling.

### Frozen literal AACOPF

AACOPF uses the P1E/P1F-B audited transition with the P1F-C frozen parameters

`alpha = 0.5`, `beta = 0`, `c_lambda = 2`.

The posterior estimate is recorded before ACO movement, and the post-transition weights are reset uniformly. No conventional systematic resampling is mixed into this path.

## Success and error metrics

Pose convergence uses the Phase-1 paper-state criterion:

- position error below 1 m;
- absolute wrapped yaw error below 10 deg;
- both conditions hold continuously to the end of the run;
- at least 5 s remain after the convergence time.

Late-window RMSE is computed over the final 20 s. The primary comparison is the matched success table at each budget: both methods converge, AACOPF only, conventional PF only, or neither. Paired late-position and late-yaw RMSE differences are retained as secondary metrics because unconditional mean RMSE can be dominated by wrong-mode runs.

## Initial-support diagnostic

A finite global particle cloud cannot recover a mode that is never represented when neither method contains rejuvenation/process noise. P1F-D therefore counts the initial particles satisfying position error below 1 m and yaw error below 10 deg. Results are stratified by zero particles, one or two particles, three or more particles, and any positive support.

## Diversity diagnostics

For the conventional PF, P1F-D records the number of systematic-resampling events and the minimum unique-particle fraction after resampling. For AACOPF, it retains minimum unique-parent fraction, catastrophic ancestry collapse (`<0.1`), dominant cloning (one destination receives at least 50% of the population), maximum destination multiplicity, moved fraction, and runtime.

## Results

### Aggregate matched comparison

| Np | Initial support >0 | PF success | AACOPF success | PF late position RMSE | AACOPF late position RMSE | PF late yaw RMSE | AACOPF late yaw RMSE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 400 | 50% | 0/20 | 0/20 | 17.678 m | 16.894 m | 56.45 deg | 50.32 deg |
| 800 | 95% | 1/20 | 0/20 | 21.371 m | 19.080 m | 64.89 deg | 57.33 deg |
| 1,200 | 95% | 2/20 | 1/20 | 19.945 m | 12.782 m | 66.95 deg | 33.99 deg |

The strict pose-convergence result is therefore negative with respect to the primary hypothesis: at no tested budget does the frozen literal AACOPF achieve a higher convergence fraction than conventional PF.

At the primary 1,200-particle budget, the matched outcome counts are

- both methods succeed: 0;
- AACOPF only succeeds: 1;
- conventional PF only succeeds: 2;
- neither succeeds: 17.

The successful realizations are therefore disjoint rather than showing a consistent AACOPF rescue of conventional-PF failures.

### Secondary paired error result

The convergence result should not hide a secondary effect. At 1,200 particles, AACOPF has lower late position RMSE than PF in 65% of matched runs and lower late yaw RMSE in 70%. The mean paired differences, AACOPF minus PF, are -7.162 m in late position RMSE and -32.96 deg in late yaw RMSE. Aggregate mean late yaw RMSE decreases from 66.95 deg to 33.99 deg, while mean late position RMSE decreases from 19.95 m to 12.78 m.

This means the frozen ACO operator often selects a global hypothesis closer to the truth even though it rarely enters and retains the strict 1 m / 10 deg terminal convergence region. This is an error-shaping effect, not evidence of a higher global-localization success rate.

### Initial-support result and transfer from P1F-C

The finite random annulus cloud is extremely sparse in the diagnostic correct region:

- 400 particles: mean support 0.80, positive in 10/20 runs;
- 800 particles: mean support 2.05, positive in 19/20 runs;
- 1,200 particles: mean support 2.35, positive in 19/20 runs, range 0--5.

This is fundamentally different from the P1F-C tuning environment. There, the correct mode initially contained 200 particles in the balanced 50/50 case and 80 particles even in the minority-correct 20/80 case. P1F-D therefore shows that a tuple that performs extremely well when the useful mode already owns substantial population does not automatically transfer to a global random-annulus cloud in which only a few correct-region particles exist.

At 1,200 particles, conditioning on positive initial support gives 2/19 PF successes and 1/19 AACOPF success. Restricting to the nine runs with at least three correct-region particles gives 1/9 success for each method. Thus positive support is necessary but a handful of particles is still insufficient for reliable survival and mode selection.

### Genealogical behavior and runtime

| Np | AACOPF catastrophic collapse | AACOPF dominant clone | PF runtime/run | AACOPF runtime/run |
|---:|---:|---:|---:|---:|
| 400 | 5% | 85% | 0.130 s | 6.19 s |
| 800 | 30% | 80% | 0.190 s | 21.83 s |
| 1,200 | 10% | 85% | 0.248 s | 47.08 s |

The P1F-C freeze reduces the incidence of the most severe ancestry collapse compared with the provisional P1F-B tuple, but dominant cloning remains common. At 1,200 particles the literal AACOPF is approximately 190 times slower than the conventional PF in the current Python implementation. This is consistent with the dense quadratic transition structure already established in P1F-B.

## Interpretation

P1F-D provides an important correction to the optimistic controlled result from P1F-C. The frozen AACOPF does not provide a robust global-convergence advantage when the correct mode is represented by only a few initial particles. The ACO mechanism can amplify an existing mode, but without process noise or rejuvenation it cannot create missing state support, and aggressive cloning can still eliminate rare useful hypotheses.

The appropriate conclusion is therefore not that AACOPF is ineffective. At the primary budget it shows a meaningful tendency to reduce the severity of wrong-mode errors. Rather, the literal mechanism by itself does not solve the sparse-support global-initialization problem.

The frozen tuple must **not** be retuned using the P1F-D seeds. Doing so would invalidate the held-out comparison. Any change intended to improve sparse-mode preservation belongs to a new explicitly labelled experiment or to the later scalable adaptation.

## P1F-D decision

P1F-D is complete. The primary hypothesis of improved strict pose-convergence fraction is **not supported** by the matched equal-budget experiment. At 1,200 particles, conventional PF succeeds in 2/20 runs and AACOPF in 1/20. A secondary error-reduction effect is present, with AACOPF reducing mean late position and yaw RMSE and outperforming PF on those metrics in a majority of matched runs.

Proceed to **P1F-E** as planned. The stationary and constant-bearing rank-deficient geometries should now be used as negative controls with the same frozen AACOPF settings. Their purpose is not to optimize performance but to verify that the particle-transition mechanism does not appear to overcome genuine missing observability. P1F-F then addresses non-Gaussian/outlier ranging, while P1F-G must address both the quadratic computational burden and the sparse-support/diversity problem revealed here.

## Reproducibility

- workflow run: `33881931650`;
- artifact: `p1f-d-results`, ID `9941206807`;
- configuration: `configs/phase1f_d.yaml`;
- runner: `experiments/run_phase1f_d.py`;
- aggregate evidence: `results/phase1/p1f_d/summary.json` and `summary.md`;
- full per-run raw results remain in the workflow artifact.
