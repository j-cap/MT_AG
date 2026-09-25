# Phase P2C — matched-budget pairwise link selection

## Scientific question

At a fixed number of pairwise UWB exchanges, does choosing which two mobile nodes range improve relative fleet geometry, worst-node error, or absolute fleet error compared with an inexpensive cyclic or random choice? P2A and P2B already established that common translation remains effectively unobservable; consequently a small change in fleet RMSE need not mean that link selection has no effect on relative localization.

P2C keeps the P2A trajectory-first four-node simulator, 20 matched seeds, IMU and UWB noise draws, known initial poses, velocity perturbations, six candidate links, and joint EKF settings unchanged. Within each seed the four policies see identical measured IMUs and pre-generated potential range values. Range realizations are sampled for every possible link and opportunity before masking, so policy decisions cannot change the measurement random-number stream.

## Cost and timing control

The three budgets match P2B's 0.5, 1, and 2 Hz all-link conditions: **180, 360, 720** undirected exchanges over 60 s. There are 600 common 10 Hz decision opportunities. If the budget is \(B\), opportunity \(s=1,\ldots,600\) has capacity

\[
c_s=\left\lfloor \frac{sB}{600}\right\rfloor-
    \left\lfloor \frac{(s-1)B}{600}\right\rfloor.
\]

Every P2C policy uses the **same \(c_s\) at every opportunity**, not merely the same overall exchange count. At 180 and 360 exchanges, active ticks allow one link; at 720, 480 ticks allow one link and 120 allow two. Exactly \(B\) links are used and no pair is repeated within a tick. There are never more than six exchanges per tick. This makes the comparison among P2C policies about **pair selection**. P2B's matching-budget all-six-link rounds occur at different times, so differences from P2B combine timing and pair-choice effects and are labeled contextual rather than causal link-selection effects.

## Frozen link policies

All policies select among the six undirected pairs. The range is applied to both uncertain endpoint states by the same joint EKF.

1. **Cyclic:** walk through the six pair indices in their canonical order, resuming where the preceding opportunity ended. This gives nearly identical per-pair usage.
2. **Random:** draw the allowed number of distinct pairs uniformly at each opportunity using an independent scheduler seed. Its draw is independent of sensor-noise and initial-state random streams.
3. **Geometry:** use only the EKF-predicted positions to construct the two-dimensional, translation-invariant range Jacobian for each candidate. Maintain a rolling \(8\times8\) geometry Gramian, with 2 s exponential memory and a predeclared unit diagonal regularizer. Greedily add the candidate with the largest \(h^\top(G+I)^{-1}h\) after scaling the range row by \(1/\sigma_{\rm UWB}\); add its row outer product to the temporary Gramian before choosing a second link. This selects geometric directions that complement recently measured directions, independently of the EKF covariance.
4. **Information:** form the predicted EKF range Jacobian \(H_{ij}\) and covariance \(P^-\). For candidate pair \((i,j)\), let \(v=P^-H_{ij}^{\top}\) and \(S=H_{ij}P^-H_{ij}^{\top}+\sigma_{\rm UWB}^2\). Remove the per-axis mean over the four nodes from the position entries of \(v\), yielding \(v_{\rm rel}\). Greedily maximize \(\|v_{\rm rel}\|^2/(4S)\), the predicted reduction of average **relative-position covariance trace**. When a second link is allowed, update a temporary covariance by \(P\leftarrow P-vv^\top/S\) before scoring the remaining pairs.

Both adaptive choices are made **after IMU prediction but before any UWB update at that tick**. Their input is the predicted estimate, covariance, and eligible-pair set. They cannot inspect the current noisy range, simulation truth, or future sensor samples. The information score targets the observable relative state rather than the common-translation error mode.

The geometry memory and regularizer are fixed in `configs/phase2c.yaml` before the campaign. There is no hyperparameter sweep or held-out trajectory in P2C; all reported comparisons are for this one predeclared heuristic and benchmark.

## Evaluation and reproduction

Compute the 60 s fleet, worst-node, centroid, and centered relative-shape position RMSE for every policy and seed, plus yaw RMSE as a diagnostic. Use paired **within-seed differences** to compare each policy with cyclic and random on their identical temporal schedule. Report means and sample standard deviations across seeds, the number of paired wins, per-link usage, and time-resolved mean ± sample-standard-deviation errors. Compare P2B at equal total exchanges only with the timing caveat above. Count radio effort in individual pairwise ranges; packet sizes, contention, and measured energy are outside this experiment.

From the repository root:

```bash
PYTHONPATH=src python experiments/run_phase2c.py
PYTHONPATH=src python experiments/analyze_phase2c.py
PYTHONPATH=src python experiments/generate_phase2c_figures.py
```

The runner validates exact per-tick and total budgets and writes compact metrics and link-use records to `results/phase2/p2c/`. The analysis script uses 30,000 paired-seed bootstrap resamples with fixed random seeds to write descriptive percentile intervals to `paired_uncertainty.json`. The figure script produces the detailed report PDFs in `report/figures/phase2/`.

## Results

All 20 seeds and three budgets completed without schedule violations. Mean whole-run errors are:

| Budget | Policy | Fleet RMSE [m] | Worst-node RMSE [m] | Shape RMSE [m] |
|---:|---|---:|---:|---:|
| 180 | Cyclic | 1.0639 | 1.0990 | 0.1095 |
| 180 | Random | 1.0639 | 1.1058 | 0.1174 |
| 180 | Geometry | 1.0635 | 1.0910 | 0.1066 |
| 180 | Information | 1.0637 | 1.0889 | 0.1068 |
| 360 | Cyclic | 1.0610 | 1.0855 | 0.0824 |
| 360 | Random | 1.0621 | 1.0960 | 0.0914 |
| 360 | Geometry | 1.0607 | 1.0856 | 0.0853 |
| 360 | Information | 1.0603 | 1.0831 | 0.0789 |
| 720 | Cyclic | 1.0594 | 1.0755 | 0.0687 |
| 720 | Random | 1.0593 | 1.0752 | 0.0697 |
| 720 | Geometry | 1.0598 | 1.0801 | 0.0673 |
| 720 | Information | 1.0593 | 1.0738 | 0.0682 |

The largest information-policy shape advantage appears at 360 exchanges: 0.0789 m compared with 0.0824 m for cyclic (4.3% lower) and 0.0914 m for random (13.7% lower). On the matched 20 seeds it beats cyclic in 13 and random in 17. The paired information-minus-cyclic difference is −0.00356 ± 0.01174 m (mean ± sample SD), with a descriptive 95% paired-bootstrap interval of [−0.00859, 0.00144] m. Against random, the difference is −0.01250 ± 0.01677 m and the interval is [−0.01953, −0.00517] m. The cyclic comparison interval crosses zero, so P2C does not establish a general information-over-cyclic benefit on this fixed benchmark. These intervals are conditional on the frozen paths and unadjusted for inspecting multiple policies and budgets.

At 180 exchanges both adaptive policies have somewhat lower mean shape and worst-node error than cyclic, but their paired shape advantages are less consistent. At 720 exchanges the policies largely converge. Fleet-centroid RMSE stays around 1.056 m and absolute fleet RMSE changes little within any budget. The strongest defensible conclusion is **modest, budget-dependent benefit for intelligent pair selection relative to random**, with cyclic remaining a demanding low-complexity baseline.

The P2B all-link schedules at matching total budgets have different **per-tick exchange counts**. Any P2C--P2B difference mixes link selection and timing and cannot be interpreted as a pure effect of choosing a pair. P2D should stress packet loss, outages, and corrupted ranges on matched P2B/P2C schedules and count attempted and successful exchanges separately. The detailed report is `report/sections/03q_phase2c.tex`.
