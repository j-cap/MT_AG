# Phase P1F-C — Safety-aware AACOPF parameter tuning

## Purpose

P1F-B showed that the literal-small audited ACO transition can rapidly amplify the correct global mode, but can also catastrophically collapse particle ancestry onto a wrong mode. P1F-C therefore tunes the three under-specified AACOPF parameters `alpha`, `beta`, and the normalized threshold factor `c_lambda` before any headline PF-versus-AACOPF comparison.

This is repository tuning, not recovery of parameters reported by Han et al. The source paper does not provide numerical values for `alpha`, `beta`, or `lambda`.

## Parameterization

For source particle `i` and strictly higher-weight candidate `j`, the P1E/P1F-B interpretation uses

`score_ij = (w_j - w_i + eps_w)^alpha * (1 / (distance_ij + eps_d))^beta`

with normalized transition probability over the candidate set. A move is accepted when

`max_j P_ij > c_lambda / K_i`,

where `K_i` is the number of higher-weight candidates. Thus `c_lambda` controls selectivity relative to the uniform candidate probability. Values at or below one are expected to be aggressive; larger values require increasingly concentrated destination evidence.

The full development sweep is:

- `alpha in {0.5, 1, 2}`;
- `beta in {0, 0.5, 1, 2}`;
- `c_lambda in {0.5, 1, 2, 4, 8}`.

This gives 60 candidate settings.

## Matched development cases

Each setting is tested on the same noisy UWB realizations and the same initial particle clouds. The paper-level DR increments remain exact so that P1F-C isolates the ACO transition rather than process-model uncertainty.

Two controlled bimodal clouds are used:

1. **balanced:** 50% of the initial particles are near the correct mode and 50% near a wrong mode;
2. **minority-correct:** only 20% are near the correct mode, providing a harder mode-preservation/amplification test.

Both use 400 particles, the same 90-degree wrong position-bearing and yaw offsets as P1F-B, Gaussian UWB noise with `sigma = 0.12 m`, and an 8 s moving-auxiliary trajectory. Eight development seeds per scenario give 16 matched runs for each parameter tuple.

## Safety diagnostics

P1F-C does not rank settings by RMSE alone. For every run it records:

- terminal correct-mode concentration: at least 90% of particles in the diagnostic correct region for the final 1 s;
- terminal wrong-mode lock: at most 10% correct particles for the final 1 s;
- catastrophic ancestry collapse: unique-parent fraction below 10% at any update;
- dominant-clone event: one destination receives at least 50% of the complete population at an update;
- final correct-mode particle fraction;
- minimum and mean unique-parent fraction;
- maximum destination multiplicity;
- mean moved fraction;
- position/yaw RMSE and ACO transition runtime.

These diagnostics follow directly from the P1F-B finding that uniform post-transition weights can hide severe ancestry collapse from effective-sample-size metrics.

## Selection rule

The selection is deliberately safety-aware and fixed before looking at validation seeds.

A setting is development-eligible if its correct-mode concentration fraction is at least 0.75 over the combined matched development runs. Eligible settings are ranked lexicographically by:

1. lower wrong-mode-lock fraction;
2. lower catastrophic-collapse fraction;
3. lower dominant-clone fraction;
4. higher correct-mode concentration fraction;
5. higher final correct-mode fraction;
6. higher minimum unique-parent fraction.

If no setting reaches the development success floor, the best available setting is reported but is not described as an acceptable frozen configuration.

## Held-out validation and freeze criterion

The selected tuple is then evaluated on 20 unseen seeds per scenario, i.e. 40 validation runs. It is accepted for P1F-D only if all three conditions hold:

- correct-mode concentration fraction >= 0.80;
- wrong-mode-lock fraction <= 0.05;
- catastrophic-collapse fraction <= 0.25.

Failure of this criterion is a valid P1F-C outcome. In that case P1F-D should not silently use a tuned literal AACOPF setting; the mechanism or threshold parameterization must first be reconsidered.

## Claim boundary

P1F-C selects an explicit repository configuration for a later matched comparison. It does not claim that the selected values are those used by Han et al., and it does not optimize against the P1F-D evaluation seeds.
