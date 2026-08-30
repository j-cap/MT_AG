# MT_AG — Collaborative UWB–IMU Localization

Shared research repository for the master thesis project on **energy-efficient collaborative UWB–IMU localization**.

The repository is intended to be a communication layer between student and supervisor: code, experiment configurations, results, interpretation, and the evolving LaTeX report remain connected and auditable.

## Collaboration workflow

- **Research code and experiments** are developed on short-lived feature/experiment branches and reviewed through pull requests.
- **The LaTeX research report lives directly on `main`** so that the Overleaf GitHub integration can always access the current report.
- Mature experiment results are transferred into the report on `main` after review. The report is the persistent scientific record/backlog and the starting point for the final thesis.

## Current status

Phase 1 bootstrap is implemented: deterministic 2-D trajectory generation, ideal/noisy dead reckoning, Gaussian UWB ranging, a conventional bootstrap particle filter, a local observability diagnostic, tests, and a reproducible baseline experiment.

The current PF experiment is intentionally a **known-initial-pose debugging baseline**. It is not yet the full Han et al. AACOPF reproduction.

## Quick start

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows: .venv\Scripts\activate
pip install -e ".[dev]"
pytest -q
python experiments/run_phase1.py
```

The experiment writes metrics/data to `results/phase1/`.

## Structure

- `src/mt_ag/` — reusable simulation/estimation code
- `experiments/` — reproducible experiment entry points
- `configs/` — experiment configurations
- `tests/` — unit/integration tests
- `docs/experiments/` — research notebook and interpretation
- `results/` — generated metrics and report data
- `report/` — Overleaf-compatible LaTeX research report/thesis backlog; canonical version lives on `main`

## Collaboration rules

1. `main` should stay working.
2. One conceptual code/experiment change per branch/PR (`feat/...`, `exp/...`, `fix/...`).
3. The LaTeX report is maintained on `main` for continuous Overleaf access.
4. Every stochastic experiment has an explicit seed and committed config.
5. Reusable mathematics/code belongs in `src/`, not notebooks.
6. An experiment is finished only when its result **and interpretation** are recorded.
7. Mature findings migrate from Markdown experiment notes into the LaTeX report on `main`.
8. Do not claim paper reproduction until unknown-pose initialization and AACOPF have been audited and implemented.

## LaTeX / Overleaf

The report entry point is:

```text
report/main.tex
```

When importing/syncing this repository in Overleaf, use the `main` branch and select `report/main.tex` as the main document.

## Reference baseline

Y. Han, C. Wei, R. Li, J. Wang, and H. Yu, “A Novel Cooperative Localization Method Based on IMU and UWB,” *Sensors*, 20(2):467, 2020. DOI: 10.3390/s20020467.
