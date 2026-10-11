# P2G: continuing position reference for one moving node

P2G is a separate global-reference extension to the completed infrastructure-free P2A–P2F sequence. It asks how much a continuing global position observation of node 1 benefits the other three nodes through UWB.

## Frozen protocol

- Same four P2A paths, 60 s, 100 Hz IMU, clean Gaussian UWB (sigma 0.12 m), reference IMU amplitudes and frozen EKF Q.
- Same seeds 0–19 (offset700000), initial velocity errors, IMU latent realizations, six links, and cyclic opportunities at attempted budgets 0,180,360,720,3600.
- Node1 is still moving and retains its IMU, velocity and yaw estimates. Only its position is observed. No other true position, velocity or yaw enters the estimator after initialization.
- Every headline metric scores nodes2–4, including the unanchored control. Centroid and relative shape are recomputed on these three nodes. No Procrustes or global alignment is applied. Four-node control metrics are used only to audit inherited P2E endpoints.
- `unanchored`: no continuing reference; known initial poses as before.
- `exact_100hz`: exact node1 position at every integration step after initialization, 6000 conditioning updates. This is an ideal information scenario, not a deployable cost-free device or a theorem guaranteeing an error lower bound.
- `noisy_1hz`: independent zero-mean Gaussian position noise, sigma0.10m **per axis**, covariance0.01I square metres; reports at1,…,60s, 60 updates.
- `noisy_1hz_outage`: identical noisy report realizations, unavailable in[20,30)s, 50 updates. Reference filtering resumes at30s. UWB and IMU remain uninterrupted.
- The ideal-versus-noisy comparison changes both accuracy and report rate; it does not isolate either factor.

At each step: IMU prediction → external position conditioning when available → scheduled UWB updates. The position observation conditions the complete joint covariance with a Joseph update. A zero-noise constraint sets its position exactly and zeroes corresponding covariance rows/columns; velocity/yaw are corrected only through covariance. Reference observation precedes UWB to expose the current reference position at coincident opportunities. The external source and its distribution/delivery network are not modeled.

Cost is a **vector**: attempted undirected UWB exchanges and reference updates. Do not add these unequal transactions as a measured energy total. Ideal/reference curves are conditional comparisons at fixed reference service. Reference update counts exclude known initial poses. All UWB attempts are received and accepted; there is no range gating in P2G.

## Execution and verification

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python experiments/run_phase2g.py --workers 4`, then `python experiments/generate_phase2g_figures.py` from the repository root. On Windows cmd set the two thread variables separately before running Python. Twenty per-seed atomic checkpoints allow resume; an implementation/config fingerprint rejects stale checkpoints. `--seeds 1` runs a smoke seed and withholds headline aggregation until the full campaign is complete. Checkpoints are ignored; compact CSV/JSON and figure PDFs are committed.

400 simulations:4 cases ×5 budgets ×20 seeds. The 100 four-node unanchored endpoints reproduce P2E at1e-10 tolerance. With zero UWB, all reference conditions leave nodes2–4 byte-identical to the unanchored control, preventing unintended direct truth access. Every exact position constraint, UWB count, reference count, and three-node squared error decomposition is checked. Unit checks cover analytical noisy/exact conditioning, cross-node propagation, covariance positive semidefiniteness, and report/outage masks. Existing Phase2 tests verify unchanged unanchored behavior.

Whole-run metrics use all6001 time indices, sample SD over20 seeds, and empirical linearly interpolated p95 of **seed trajectory RMSE**. Mean instantaneous error profiles are sampled at10Hz, with mean±one sample SD (negative lower envelopes clipped at0). Paired bootstrap95% intervals use30000 fixed-seed resamples and remain descriptive/unadjusted. The 0.10m target now applies to **absolute three-node fleet RMSE**, not P2F's four-node relative-shape criterion. Minimums mean smallest tested budget; no interpolation across the large720→3600 gap.

## Findings

At360 attempts, absolute fleet mean±SD is1.0613±0.4488m unanchored,0.1183±0.0243m exact,0.1313±0.0261m noisy,0.1319±0.0259m noisy with outage. Every reference case improves absolute RMSE for every seed at every positive tested UWB budget. The exact-reference paired reduction is0.9430m,95% interval[0.7614,1.1409]m (sign reversed from case-minus-control). Shape changes at360 are around millimetres and uncertain; centroid reduction explains the large absolute gain.

Ideal reference reaches mean absolute10cm at720 (0.0984m;10/20 seeds pass), but empirical p95 only at3600 (0.0713m;20/20 pass). Noisy reports reach the mean at3600 (0.0851m;17/20 pass), but p95 is0.1048m. Outage gives0.0868m and16/20 pass at3600, p950.1096m. Unanchored never reaches the absolute target.

At360, outage minus uninterrupted noisy mean is+0.00059m,paired95% interval[-0.00130,0.00257]m; at3600 it is+0.00169m,[-0.00049,0.00392]m. This short, specified outage has little effect on whole-run RMSE; it does not establish robustness to longer or differently placed outages. During the sampled[20,30)s window at360, mean reference error rises from0.0539m to0.0797m, while mean instantaneous three-node fleet error changes from0.1306m to0.1319m. These window values average10Hz mean error, not trajectory RMSE. Reference trajectories recover after reports resume.

The extension shows why a continuing global observation changes achievable absolute localization. It does not eliminate geometry-dependent uncertainty, initial velocity transients, or yaw uncertainty, and it does not prove general observability or validate a real external reference. P2A–P2G are complete. The next scientific decision is the intended application: infrastructure-free relative geometry versus globally referenced absolute localization, followed by a deployment-informed reference model/target or a distinct reliability-aware scheduling experiment.
