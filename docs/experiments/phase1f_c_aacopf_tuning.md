# Phase P1F-C — Safety-aware AACOPF parameter tuning

## Purpose

P1F-B showed that the literal-small audited ACO transition can rapidly amplify the correct global mode, but can also catastrophically collapse particle ancestry onto a wrong mode. P1F-C therefore tunes the three under-specified AACOPF parameters `alpha`, `beta`, and the normalized threshold factor `c_lambda` before any headline PF-versus-AACOPF comparison.

This is repository tuning, not recovery of parameters reported by Han et al. The source paper does not provide numerical values for `alpha`, `beta`, or `lambda`.

## Parameterization

For source particle `i` and strictly higher-weight candidate `j`, the P1E/P1F-B interpretation uses

`score_ij = (w_j - w_i + eps_w)^alpha * (1 / (distance_ij + eps_d))^beta`

with normalized transition probability over the candidate set. A move is accepted when

`max_j P_ij > c_lambda / K_i`,

where `K_i` is the number of higher-weight candidates. Thus `c_lambda` controls selectivity relative to the uniform candidate probability. Values at or below one are aggressive in the tested regime; larger values require increasingly concentrated destination evidence.

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

## Development results

The safety-aware ranking selects

`alpha = 0.5`, `beta = 0`, `c_lambda = 2`.

Across the 16 development runs this tuple gives:

| Metric | Result |
|---|---:|
| Correct-mode terminal lock | 15/16 = 93.75% |
| Wrong-mode terminal lock | 0/16 = 0% |
| Catastrophic ancestry collapse | 1/16 = 6.25% |
| Dominant-clone event | 10/16 = 62.5% |
| Mean final correct-mode fraction | 0.993125 |
| Mean minimum unique-parent fraction | 0.402344 |
| Mean moved-particle fraction | 0.039879 |

The next two safety-ranked settings are `(alpha,beta,c_lambda)=(1,0,4)` and `(2,0,8)`. Both also have `beta=0`. In contrast, the P1F-B mechanism-check point `(1,1,0.5)` gives only 75% correct lock, 25% wrong lock, catastrophic ancestry collapse in 100% of development runs, dominant-clone events in 100%, and a mean minimum unique-parent fraction of only 0.0316.

Two qualitative effects are therefore clear in this controlled study. First, increasing the normalized threshold above the aggressive P1F-B value substantially reduces premature cloning. Second, the best safety-ranked settings remove the inverse-distance contribution (`beta=0`). This does not establish that spatial distance is generally harmful; it establishes only that, for this controlled bimodal UWB experiment, the distance heuristic is unnecessary for reliable mode selection and can contribute to locally attractive cloning decisions.

## Held-out validation and freeze criterion

The selected tuple is frozen before validation and then evaluated on 20 unseen seeds per scenario, i.e. 40 validation runs. It is accepted for P1F-D only if all three conditions hold:

- correct-mode concentration fraction >= 0.80;
- wrong-mode-lock fraction <= 0.05;
- catastrophic-collapse fraction <= 0.25.

The held-out result is:

| Metric | All validation runs | Balanced | Minority-correct |
|---|---:|---:|---:|
| Runs | 40 | 20 | 20 |
| Correct-mode terminal lock | **100%** | **100%** | **100%** |
| Wrong-mode terminal lock | **0%** | **0%** | **0%** |
| Catastrophic ancestry collapse | **22.5%** | 25% | 20% |
| Mean final correct-mode fraction | 0.999625 | 0.999375 | 0.999875 |

Across all 40 validation runs, the mean minimum unique-parent fraction is 0.3414 and the mean moved-particle fraction is 0.0392. Mean full-horizon position and yaw RMSE are 2.846 m and 21.812 deg, respectively. These RMSE values include the deliberately ambiguous initial interval and are not used as final localization-performance claims.

The frozen tuple therefore **passes** the predeclared P1F-C acceptance criterion:

`freeze_accepted = true`.

## Interpretation of the residual diversity risk

Passing the freeze criterion does not mean that the literal transition has become benign. A dominant-clone event still occurs in 85% of held-out runs, and the catastrophic-collapse diagnostic is triggered in 9/40 runs. The difference relative to P1F-B is that these events no longer produce a terminal wrong-mode lock in the held-out controlled experiment.

The correct interpretation is therefore narrow:

- the tuple is sufficiently reliable to be frozen for an unbiased P1F-D comparison;
- the tuning substantially reduces the wrong-lock failure observed with the provisional P1F-B setting;
- the literal ACO mechanism remains genealogically aggressive;
- P1F-G is still required before making any diversity, scalability, or embedded-deployment claim.

The frequent dominant-clone diagnostic does not violate the freeze criterion because a maximum dominant-clone frequency was deliberately not part of the predeclared acceptance rule. It remains an important warning metric that must be carried into P1F-D.

## Claim boundary

P1F-C selects an explicit repository configuration for a later matched comparison. It does not claim that the selected values are those used by Han et al., and it does not optimize against the P1F-D evaluation seeds. The selected values are

`alpha = 0.5`, `beta = 0`, `c_lambda = 2`.

They must be used unchanged in P1F-D unless a separately documented methodological reason requires reopening P1F-C.

## Reproducibility

- configuration: `configs/phase1f_c.yaml`;
- runner: `experiments/run_phase1f_c.py`;
- aggregate result: `results/phase1/p1f_c/summary.json` and `summary.md`;
- workflow run: `33803059007`;
- full generated result artifact: `p1f-c-results`, artifact ID `9912216810`.

## Decision

P1F-C is complete. Proceed to **P1F-D**, using the frozen tuple `(alpha,beta,c_lambda)=(0.5,0,2)` for the literal-small AACOPF. P1F-D must compare conventional PF and AACOPF under matched trajectories, UWB realizations, initial particle clouds, and declared particle budgets, while retaining correct-mode survival, ancestry-collapse, destination-multiplicity, and runtime diagnostics alongside localization error.
