# P1F-F summary — non-Gaussian / outlier UWB stress

Both methods keep the nominal Gaussian likelihood with sigma=0.12 m. Only the measurement process changes, so this is not a robust-likelihood comparison.

## Initial cloud: balanced

| condition | PF lock | AACOPF lock | PF wrong | AACOPF wrong | PF pos RMSE [m] | AACOPF pos RMSE [m] | AACOPF collapse |
|---|---:|---:|---:|---:|---:|---:|---:|
| clean | 1.00 | 0.95 | 0.00 | 0.00 | 1.715 | 2.251 | 0.05 |
| sparse_impulses | 1.00 | 0.80 | 0.00 | 0.15 | 1.817 | 2.740 | 0.50 |
| mixture_moderate | 0.90 | 0.70 | 0.10 | 0.25 | 2.422 | 2.649 | 0.75 |
| mixture_severe | 0.95 | 0.65 | 0.05 | 0.35 | 2.095 | 2.954 | 0.90 |
| nlos_burst_moderate | 0.75 | 0.60 | 0.25 | 0.35 | 3.031 | 3.819 | 0.50 |
| nlos_burst_severe | 0.75 | 0.60 | 0.25 | 0.40 | 3.071 | 4.443 | 0.45 |

### Degradation relative to clean

| condition | PF lock drop | AACOPF lock drop | PF clean locks lost | AACOPF clean locks lost |
|---|---:|---:|---:|---:|
| sparse_impulses | 0.00 | 0.15 | 0 | 3 |
| mixture_moderate | 0.10 | 0.25 | 2 | 6 |
| mixture_severe | 0.05 | 0.30 | 1 | 7 |
| nlos_burst_moderate | 0.25 | 0.35 | 5 | 7 |
| nlos_burst_severe | 0.25 | 0.35 | 5 | 8 |

## Initial cloud: minority_correct

| condition | PF lock | AACOPF lock | PF wrong | AACOPF wrong | PF pos RMSE [m] | AACOPF pos RMSE [m] | AACOPF collapse |
|---|---:|---:|---:|---:|---:|---:|---:|
| clean | 1.00 | 0.95 | 0.00 | 0.00 | 2.745 | 3.508 | 0.15 |
| sparse_impulses | 0.95 | 0.75 | 0.05 | 0.20 | 2.819 | 4.159 | 0.45 |
| mixture_moderate | 0.90 | 0.80 | 0.10 | 0.15 | 3.354 | 4.147 | 0.55 |
| mixture_severe | 0.85 | 0.80 | 0.10 | 0.20 | 3.007 | 4.550 | 0.70 |
| nlos_burst_moderate | 0.80 | 0.70 | 0.20 | 0.25 | 3.894 | 4.945 | 0.30 |
| nlos_burst_severe | 0.75 | 0.65 | 0.25 | 0.25 | 4.169 | 5.205 | 0.45 |

### Degradation relative to clean

| condition | PF lock drop | AACOPF lock drop | PF clean locks lost | AACOPF clean locks lost |
|---|---:|---:|---:|---:|
| sparse_impulses | 0.05 | 0.20 | 1 | 4 |
| mixture_moderate | 0.10 | 0.15 | 2 | 4 |
| mixture_severe | 0.15 | 0.15 | 3 | 4 |
| nlos_burst_moderate | 0.20 | 0.25 | 4 | 6 |
| nlos_burst_severe | 0.25 | 0.30 | 5 | 7 |

## Pooled stressed-run diagnostic

Across the 200 stressed runs (five corruption conditions, two initial-cloud scenarios, 20 seeds each), conventional PF achieves 172/200 correct-mode locks (86.0%) and 27/200 wrong locks (13.5%). Frozen AACOPF achieves 141/200 correct-mode locks (70.5%) and 51/200 wrong locks (25.5%). In every stressed cell AACOPF has a correct-lock fraction no higher than PF; it is strictly lower in all ten cells.

AACOPF catastrophic ancestry collapse occurs in 111/200 stressed runs. Wrong-mode lock occurs in 40/111 (36.0%) of collapse runs versus 11/89 (12.4%) without catastrophic collapse. This is a descriptive association rather than a causal proof, but it is consistent with the mechanism identified in P1F-B/P1F-C: a corrupted likelihood can make an incorrect mode temporarily attractive and the copying transition can then make that selection difficult to reverse.

## Interpretation

The frozen literal AACOPF does **not** provide inherent robustness to non-Gaussian or NLOS-like ranging when both methods continue to use the nominal Gaussian likelihood. The clearest result is the two-impulse test: conventional PF retains 100%/95% correct lock in the balanced/minority-correct clouds, while AACOPF falls to 80%/75% and introduces 15%/20% wrong locks.

Mixture contamination and positive NLOS-like bursts degrade both methods, but the same pattern persists. AACOPF generally loses more clean successes, exhibits more wrong-mode locking, and has larger position/yaw RMSE. Therefore particle copying is not a substitute for explicit robust measurement handling. If robustness to NLOS/outliers is required later in the thesis, it should be introduced through a separately identified likelihood/gating/reliability mechanism rather than attributed to the audited AACOPF transition itself.
