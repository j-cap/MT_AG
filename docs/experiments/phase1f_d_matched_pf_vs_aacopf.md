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

The 1,200-particle case is the primary literal-small comparison. This is the largest population already characterized in the P1F-B runtime study and is large enough that the random global initialization often contains at least one particle in the diagnostic correct region. The lower budgets deliberately expose the interaction between initial support and particle-management behavior.

P1F-A's 10,000- and 40,000-particle conventional-PF results remain contextual upper-population references, not matched AACOPF comparisons.

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

The first measurement is used only to construct the initial cloud; likelihood updates begin at the next sample to avoid double-counting `z0`.

Twenty unseen seeds are evaluated at each particle budget. These random streams are separate from the P1F-C development and validation streams.

## Methods

### Conventional PF

The conventional method is the P1F-A paper-state bootstrap PF with systematic resampling triggered when

`N_eff < 0.5 Np`.

There is no artificial process noise in `Delta L` or `Delta phi` in this experiment, so both methods receive the same exact DR propagation. The purpose is to compare particle replacement/mode management rather than process-noise handling.

### Frozen literal AACOPF

AACOPF uses the P1E/P1F-B audited transition with the P1F-C frozen parameters:

`alpha = 0.5`, `beta = 0`, `c_lambda = 2`.

The posterior estimate is recorded before ACO movement, and the post-transition weights are reset uniformly. No conventional systematic resampling is mixed into this path.

## Success and error metrics

Pose convergence uses the Phase-1 paper-state criterion:

- position error below 1 m;
- absolute wrapped yaw error below 10 deg;
- both conditions hold continuously to the end of the run;
- at least 5 s remain after the convergence time.

Late-window RMSE is computed over the final 20 s.

The primary comparison is the matched success table at each budget:

- both methods converge;
- AACOPF only converges;
- conventional PF only converges;
- neither converges.

Paired late-position and late-yaw RMSE differences are retained as secondary metrics because unconditional mean RMSE can be dominated by wrong-mode runs.

## Initial-support diagnostic

A finite global particle cloud cannot recover a mode that is never represented when neither method contains rejuvenation/process noise. P1F-D therefore explicitly counts the initial particles satisfying

- position error < 1 m;
- yaw error < 10 deg.

Results are stratified by initial correct-region support:

- zero particles;
- one or two particles;
- three or more particles;
- any positive support.

This diagnostic prevents failures caused by missing initial support from being incorrectly attributed to the resampling mechanism.

## Diversity diagnostics

For the conventional PF, P1F-D records the number of systematic-resampling events and the minimum unique-particle fraction after resampling.

For AACOPF, P1F-D retains the P1F-B/P1F-C genealogy diagnostics:

- minimum unique-parent fraction;
- catastrophic ancestry collapse: unique-parent fraction < 0.1 at any update;
- dominant clone: one destination receives at least 50% of the population;
- maximum destination multiplicity;
- mean moved-particle fraction;
- total literal-transition runtime.

The comparison therefore distinguishes localization success from genealogical health.

## Interpretation rules fixed before results

1. AACOPF is considered beneficial only if it improves matched convergence without hiding a severe support imbalance.
2. Zero-support runs are reported separately; neither method can be expected to reconstruct an absent global mode in this deterministic-propagation experiment.
3. Higher AACOPF convergence does not by itself imply deployment readiness if ancestry collapse remains frequent.
4. Runtime is interpreted as the cost of the literal reproduction only. P1F-G remains responsible for scalable adaptation.
5. The P1F-C parameters are not retuned after seeing P1F-D results.

## Reproducibility

- configuration: `configs/phase1f_d.yaml`
- runner: `experiments/run_phase1f_d.py`
- full generated results: GitHub Actions artifact
- aggregate summary: `results/phase1/p1f_d/summary.json` and `summary.md` after the campaign completes
