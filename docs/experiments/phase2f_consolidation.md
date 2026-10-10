# P2F: conditional accuracy–communication consolidation

P2F analyzes committed P2A–P2E results; it adds no simulations or estimator tuning.

Run `python experiments/analyze_phase2f.py`, then `python experiments/generate_phase2f_figures.py` from the repository root. The configuration fixes the source commit, target, numerical tolerance, and bootstrap settings. Source SHA-256 hashes are retained in `summary.json`.

The audited input contains 1960 source records, 200 repeated endpoint checks, and 1760 distinct seed records across 88 designs. Repeated endpoints are checked rather than pooled as extra replications. Every design must contain seeds 0–19. Counts and the squared fleet/centroid/shape decomposition are checked. The five compact output CSVs retain seed-level provenance, design statistics, target decisions, paired comparisons, and endpoint audits.

Pareto labels minimize attempted exchanges and an error statistic within fixed sensor/Q/range/update conditions. Exact ties remain. Schedule and timing stay explicit; uniform all-link rounds versus P2C selection compares entire designs, while P2C cyclic/selection isolates pair choice at matched opportunities. Separate labels are calculated for mean shape, mean fleet, mean worst node, and empirical p95 shape. These are tested-design frontiers without uncertainty-based dominance or interpolation.

The 0.10 m target is descriptive. Report the smallest tested passing budget separately for mean and empirical p95, together with passing seed counts. Linear-interpolated percentiles use 20 seed-level trajectory RMSE values, not instantaneous samples. Paired differences use 30000 fixed-seed bootstrap resamples, without multiple-comparison adjustment. Any schedule selected after observing these same runs is an empirical ranking, not validated optimality.

Reference cyclic crosses the mean target at 360 attempts (18/20 seeds) and p95 at 720. Uniform rounds need the same tested budgets. Lower white noise with rescaled Q crosses the mean at 180, but only 11/20 seeds pass and p95 requires 360. Degraded sensors fail the mean target even at full communication. Reliability has only budget360; passing there does not establish its minimum budget. No reliability condition meets the empirical p95 target. Full-rate absolute error remains dominated by centroid drift. Hardware energy, general deployment guarantees, and commercial IMU grades are outside this evidence.

The next optional experiment is P2G, a continuing global position reference for one existing moving node, scoring the same three non-reference nodes at matched UWB budgets. Noisy/intermittent references must propagate uncertainty and count updates separately. Reliability-aware scheduling and bias-state estimation are separate follow-ups.
