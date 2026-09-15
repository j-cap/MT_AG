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

## Interpretation rule

A robustness claim requires stress degradation to be systematically smaller for AACOPF than for conventional PF under paired seeds and frozen settings. If AACOPF instead loses the correct mode as often as or more often than PF, the conclusion is that copying high-weight particles is **not** a substitute for explicit robust measurement handling and may amplify an outlier-induced weight error.
