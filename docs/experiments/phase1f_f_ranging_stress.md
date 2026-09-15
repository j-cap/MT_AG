# Phase P1F-F — non-Gaussian / outlier UWB stress

## Purpose

P1F-F tests the robustness motivation of the paper-level PF/AACOPF comparison under corrupted ranging. The scientific question is deliberately narrow:

> With the likelihood model and all AACOPF parameters frozen, does the literal AACOPF particle-management mechanism provide inherent resilience to non-Gaussian or NLOS-like UWB errors relative to conventional resampling?

This is **not** a robust-likelihood experiment. Both conventional PF and AACOPF continue to evaluate every measurement with the same nominal Gaussian likelihood, `sigma_UWB = 0.12 m`. No innovation gate, Huber loss, mixture likelihood, NLOS detector, rejection rule, or parameter retuning is introduced.

## Why the controlled bimodal setting is used

P1F-D showed that random-annulus global initialization at manageable literal-AACOPF populations is dominated by sparse correct-mode support: at `Np=1200` only 0--5 particles lie in the diagnostic correct region. That setting is therefore poorly suited for isolating measurement-corruption robustness, because clean-data failures already dominate.

P1F-F instead reuses the P1F-C controlled two-mode state cloud with `Np=400` and an informative moving auxiliary. This produces a meaningful clean reference in which the correct global mode is represented before corrupted measurements are introduced. Two initial clouds are retained:

- `balanced`: 50% correct mode, 50% wrong mode;
- `minority_correct`: 20% correct mode, 80% wrong mode.

The wrong mode is offset by 90 degrees in bearing and yaw, as in P1F-C. P1F-F uses new validation seeds (`200--219`), not the P1F-C tuning/validation seeds.

## Frozen methods

Conventional PF uses systematic resampling at `N_eff < 0.5 Np`.

AACOPF uses the frozen P1F-C tuple

`(alpha, beta, c_lambda) = (0.5, 0, 2)`

with no retuning. Recall that `beta=0` removes the inverse-distance factor from the selected literal-small transition, so the frozen transition is driven by posterior-weight differences and the normalized movement threshold.

## Ranging stress conditions

Every condition starts from the same nominal Gaussian range-noise realization for a given seed. Corruption is then added on top, so stressed runs are paired directly with their clean counterpart.

1. **clean** — nominal Gaussian noise only.
2. **sparse_impulses** — two isolated gross errors of `1.5 m`, with random sign.
3. **mixture_moderate** — 8% of samples receive an additional zero-mean Gaussian gross-error term with `sigma=0.6 m`.
4. **mixture_severe** — the same corruption mask/standardized gross errors as the moderate mixture case, scaled to `sigma=1.2 m`.
5. **nlos_burst_moderate** — one contiguous `1.0 s` interval receives a `+0.6 m` positive bias plus `0.08 m` extra jitter.
6. **nlos_burst_severe** — the same random burst-start construction with a `2.0 s` duration, `+1.2 m` bias, and `0.12 m` extra jitter.

Corruption is not applied at the first sample. The initial controlled cloud is therefore identical across ranging conditions and does not inherit an outlier through first-range initialization.

## Metrics

The primary robustness outcome is correct-mode lock: the correct-mode posterior signal must remain at or above 90% to the end of the run for at least 1 s. Wrong-mode lock uses the complementary 10% threshold. P1F-F also records strict pose success, position/yaw RMSE, the minimum correct-mode signal after the first corrupted sample, and whether/when the mode recovers above 90% after the final corruption.

For conventional PF the mode signal is the pre-resampling correct-mode posterior mass. For AACOPF the lock signal is the post-transition correct-mode particle fraction; the pre-ACO correct-mode mass is retained separately.

The existing diversity diagnostics are retained:

- PF resampling count and minimum unique fraction;
- AACOPF minimum unique-parent fraction;
- catastrophic ancestry collapse below 10%;
- dominant cloning when one destination receives at least 50% of the cloud;
- moved fraction and transition runtime.

## Results

### Balanced initial cloud

| condition | PF correct lock | AACOPF correct lock | PF wrong lock | AACOPF wrong lock | PF position RMSE [m] | AACOPF position RMSE [m] | AACOPF catastrophic collapse |
|---|---:|---:|---:|---:|---:|---:|---:|
| clean | 1.00 | 0.95 | 0.00 | 0.00 | 1.715 | 2.251 | 0.05 |
| sparse impulses | 1.00 | 0.80 | 0.00 | 0.15 | 1.817 | 2.740 | 0.50 |
| mixture moderate | 0.90 | 0.70 | 0.10 | 0.25 | 2.422 | 2.649 | 0.75 |
| mixture severe | 0.95 | 0.65 | 0.05 | 0.35 | 2.095 | 2.954 | 0.90 |
| NLOS burst moderate | 0.75 | 0.60 | 0.25 | 0.35 | 3.031 | 3.819 | 0.50 |
| NLOS burst severe | 0.75 | 0.60 | 0.25 | 0.40 | 3.071 | 4.443 | 0.45 |

Relative to the clean matched seeds, the correct-lock fraction drops by `0.00/0.15` for PF/AACOPF under sparse impulses, `0.10/0.25` for the moderate mixture, `0.05/0.30` for the severe mixture, and `0.25/0.35` for both NLOS-burst severities. Thus the ACO transition does not preserve the clean-data advantage when measurement evidence is corrupted.

### Minority-correct initial cloud

| condition | PF correct lock | AACOPF correct lock | PF wrong lock | AACOPF wrong lock | PF position RMSE [m] | AACOPF position RMSE [m] | AACOPF catastrophic collapse |
|---|---:|---:|---:|---:|---:|---:|---:|
| clean | 1.00 | 0.95 | 0.00 | 0.00 | 2.745 | 3.508 | 0.15 |
| sparse impulses | 0.95 | 0.75 | 0.05 | 0.20 | 2.819 | 4.159 | 0.45 |
| mixture moderate | 0.90 | 0.80 | 0.10 | 0.15 | 3.354 | 4.147 | 0.55 |
| mixture severe | 0.85 | 0.80 | 0.10 | 0.20 | 3.007 | 4.550 | 0.70 |
| NLOS burst moderate | 0.80 | 0.70 | 0.20 | 0.25 | 3.894 | 4.945 | 0.30 |
| NLOS burst severe | 0.75 | 0.65 | 0.25 | 0.25 | 4.169 | 5.205 | 0.45 |

The minority-correct case gives the same qualitative result. Relative to clean, AACOPF loses more correct-lock realizations for every stress family. The severe-mixture correct-lock drop is equal (`0.15`) for PF and AACOPF, but AACOPF still starts from a lower clean correct-lock fraction and has higher wrong-lock frequency and larger aggregate errors.

### Pooled stressed-run result

Across the 200 stressed runs (five corruption conditions, two initial clouds, 20 seeds per cell), conventional PF obtains 172/200 correct-mode locks (`86.0%`) and 27/200 wrong locks (`13.5%`). Frozen AACOPF obtains only 141/200 correct-mode locks (`70.5%`) and 51/200 wrong locks (`25.5%`). AACOPF has a lower correct-lock fraction than PF in **all ten stressed cells**.

By corruption family, pooled correct-lock rates are:

- sparse impulses: PF `97.5%`, AACOPF `77.5%`;
- mixture contamination: PF `90.0%`, AACOPF `73.75%`;
- NLOS bursts: PF `76.25%`, AACOPF `63.75%`.

The corresponding recovery fractions after the final corruption are `97.5%/82.5%`, `55.0%/42.5%`, and `71.25%/65.0%` for PF/AACOPF, respectively.

### Ancestry collapse and wrong-mode amplification

Catastrophic AACOPF ancestry collapse occurs in 111/200 stressed runs (`55.5%`). Among these collapse runs, 40/111 (`36.0%`) end in wrong-mode lock. Without catastrophic collapse, wrong-mode lock occurs in 11/89 (`12.4%`). This is a descriptive association, not a causal proof, but it is consistent with the mechanism identified in P1F-B/P1F-C: corrupted likelihood evidence can temporarily favor an incorrect hypothesis, after which aggressive copying makes that selection difficult to reverse.

The sparse-impulse case is especially informative because only two samples are corrupted. In the balanced cloud, conventional PF retains `100%` correct lock and `0%` wrong lock, while AACOPF falls to `80%` correct lock and `15%` wrong lock. In the minority-correct cloud, the corresponding values are `95%/5%` for PF versus `75%/20%` for AACOPF. This is direct evidence that the cloning operator can amplify a short-lived outlier-induced weight mistake.

## Interpretation and decision

P1F-F does **not** support an inherent robustness advantage of the frozen literal AACOPF. Under a deliberately misspecified Gaussian likelihood, AACOPF is systematically less robust than conventional PF in this controlled campaign: it loses more clean successes, wrong-locks more often, and generally has larger position and yaw RMSE.

The result also sharpens the interpretation of the Han-style particle transition. The transition operates on posterior weights; it does not distinguish trustworthy measurement information from a transient outlier-induced likelihood peak. Consequently, stronger exploitation of the currently highest-weight mode can reduce robustness when the likelihood itself is wrong.

The correct transfer decision is therefore:

1. do **not** attribute UWB outlier/NLOS robustness to the audited AACOPF transition;
2. if later thesis phases require robust ranging, introduce that property explicitly through a separately identified likelihood, innovation gate, NLOS/reliability model, or measurement-selection rule;
3. keep P1F-G focused on scalable candidate selection plus diversity/support preservation rather than trying to retune `alpha`, `beta`, or `c_lambda` to repair measurement-model mismatch.

The non-monotonic finite-sample difference between moderate and severe mixture results (for example, balanced PF lock `0.90` versus `0.95`) is not interpreted as a severity benefit. Both variants share the corruption pattern but only 20 seeds are evaluated, so cell-level ordering can remain stochastic. The robust conclusion is the systematic PF-versus-AACOPF comparison across all matched stress cells.

## Reproducibility

Configuration: `configs/phase1f_f.yaml`.

Runner: `experiments/run_phase1f_f.py`.

Aggregator: `experiments/aggregate_phase1f_f.py`.

Workflow run: `35000315361` (`P1F-F Ranging Stress`).

Aggregate artifact: `p1f-f-results`, artifact ID `10409367838`.

Version-controlled compact evidence: `results/phase1/p1f_f/summary.json` and `results/phase1/p1f_f/summary.md`. The full per-run evidence remains in the workflow artifact.
