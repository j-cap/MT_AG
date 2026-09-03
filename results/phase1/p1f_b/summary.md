# P1F-B result summary

## Main finding

The literal-small ACO transition is implemented and behaves exactly as the P1E interpretation intends, but the first controlled campaign exposes a central risk: **the transition can collapse particle diversity extremely aggressively and can lock onto either the correct or the wrong mode depending on the early likelihood realization**.

This is not yet a PF-vs-AACOPF performance result. It is the mechanism result P1F-B was designed to obtain.

## Controlled bimodal campaign

Configuration: 10 seeds, 400 particles, 12 s, exact paper-level DR increments, moving known auxiliary, Gaussian UWB with `sigma = 0.12 m`, `alpha = beta = 1`, provisional `c_lambda = 0.5`.

The initial cloud contains an equal correct/wrong split. The wrong mode lies 90 deg around the same first-range circle and has a 90 deg yaw offset, so the first range alone cannot distinguish the two modes.

Aggregate results:

| Metric | Result |
|---|---:|
| Runs that ultimately concentrate >=90% on the correct mode | 9/10 |
| Median concentration time among those runs | 0.60 s |
| Mean final correct-mode particle fraction | 0.90 |
| Mean final pre-ACO correct-mode posterior mass | 0.90 |
| Mean moved-particle fraction across all updates | 0.0442 |
| Mean unique-parent fraction across all updates | 0.9595 |
| Mean maximum destination multiplicity per run | 300.7 particles |
| Mean transition runtime at N=400 | 5.89 ms/update |

The low average moved fraction is misleading because movement is highly episodic. At informative early updates the transition can move almost the entire cloud.

### Representative correct lock: seed 0

At `t = 0.5 s`, 94.25% of particles move and only 6% of pre-transition parents remain represented. The pre-ACO posterior already assigns 79.7% mass to the correct region, and after the transition 66.25% of particles are in that region. At the next update (`t = 0.6 s`) the cloud reaches 100% correct-mode particles.

### Representative wrong lock: seed 4

At `t = 0.3 s`, 98.5% of particles move and the unique-parent fraction falls to 4%. Although the pre-ACO posterior still assigns 76.4% mass to the correct region, the transition leaves only 46.25% of particles there. By `t = 0.5 s`, another 92.75% of particles move, only 7.5% of parents survive, one destination is copied by 302 particles, and the correct-mode fraction becomes exactly zero. The run then remains locked to the wrong mode.

This shows that the ACO transition is not merely a gentle anti-starvation correction under the current provisional threshold. It is a **mode-selection operator with potentially irreversible cloning events**.

## Consequence for diagnostics

Pre-ACO effective sample size is not sufficient to diagnose this behavior. Because the audited implementation resets weights uniformly after each ACO transition, the filter can have a high effective sample size while the actual particle cloud has already lost most of its distinct ancestry. P1F-C and P1F-D must therefore retain explicit diversity diagnostics such as unique-parent fraction, destination multiplicity, and correct-mode survival.

## Literal all-pairs computational scaling

Median single-transition runtimes:

| Particles | Dense pairs | Higher-weight candidate scores | Median runtime |
|---:|---:|---:|---:|
| 50 | 2,500 | 1,225 | 0.315 ms |
| 100 | 10,000 | 4,950 | 0.633 ms |
| 200 | 40,000 | 19,900 | 2.08 ms |
| 400 | 160,000 | 79,800 | 7.03 ms |
| 800 | 640,000 | 319,600 | 33.05 ms |
| 1,200 | 1,440,000 | 719,400 | 73.32 ms |

The measured log-log runtime slope over this range is approximately 1.76. The structural implementation remains dense all-pairs and therefore has quadratic storage/work scaling by construction.

For context, the `N_p = 160,000` population implied by Han et al.'s `N=M=400` would contain `2.56e10` dense pair entries. A single float64 pairwise matrix alone would require about **204.8 GB**. A naive `N^2` extrapolation from the 1,200-particle timing gives roughly **22 minutes per transition**, but that timing extrapolation is not operationally meaningful because the dense memory requirement is already prohibitive and the current implementation uses several pairwise arrays.

## P1F-B decision

P1F-B is successful as a mechanism-validation experiment:

- deterministic tests establish the intended candidate, score, threshold, synchronous-copy, and lineage semantics;
- the integrated literal-small filter executes and exposes the expected particle-management behavior;
- the controlled experiment shows that ACO can rapidly resolve a competing mode, but also that premature wrong-mode collapse is possible;
- the dense all-pairs interpretation is computationally unsuitable for the paper's nominal 160,000-particle population and certainly unsuitable for the intended embedded target.

**Do not proceed directly to the headline PF-vs-AACOPF comparison.** The next experiment should be P1F-C: tune and stress-test `alpha`, `beta`, and especially normalized `c_lambda` on development seeds, with explicit penalties for wrong-mode lock and catastrophic diversity collapse. Only then should a setting be frozen for P1F-D.

Provenance: workflow run `33751243293`, artifact `9891624402`, commit `c24f751891a6e57fe6e75d86c63340854141066a`.
