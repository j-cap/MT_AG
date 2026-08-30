# MT_AG — Collaborative UWB–IMU Localization

Shared research repository for the master thesis project on **energy-efficient collaborative UWB–IMU localization**.

The repository is intended to be a communication layer between student and supervisor: code, experiment configurations, results, interpretation, and the evolving LaTeX report live together and should remain mutually consistent.

## Current status

Phase 1 bootstrap is implemented: deterministic 2-D trajectory generation, ideal/noisy dead reckoning, Gaussian UWB ranging, a conventional bootstrap particle filter, a local observability diagnostic, tests, a reproducible baseline experiment, and a LaTeX report skeleton.

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
- `report/` — Overleaf-compatible LaTeX report/thesis backlog

## Collaboration rules

1. `main` should stay working.
2. One conceptual change per branch/PR (`feat/...`, `exp/...`, `fix/...`).
3. Every stochastic experiment has an explicit seed and committed config.
4. Reusable mathematics/code belongs in `src/`, not notebooks.
5. An experiment is finished only when its result **and interpretation** are recorded.
6. Mature findings migrate from Markdown experiment notes into the LaTeX report.
7. Do not claim paper reproduction until unknown-pose initialization and AACOPF have been audited and implemented.

## Reference baseline

Y. Han, C. Wei, R. Li, J. Wang, and H. Yu, “A Novel Cooperative Localization Method Based on IMU and UWB,” *Sensors*, 20(2):467, 2020. DOI: 10.3390/s20020467.
