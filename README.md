# MT_AG

Research code and report for the master thesis on energy-efficient collaborative UWB–IMU localization.

## Collaboration workflow

This repository is used as the communication layer between student and supervisor.

- **Research code and experiments** are developed on short-lived feature/experiment branches and reviewed through pull requests.
- **The LaTeX research report lives directly on `main`** so that the Overleaf GitHub integration can always access the current report.
- Mature experiment results are transferred into the report on `main` after they have been reviewed. The report is therefore the persistent scientific record/backlog and the starting point for the final thesis.

The current Phase 1 implementation is developed in PR #1 (`phase1/bootstrap`).

## LaTeX / Overleaf

The report entry point is:

```text
report/main.tex
```

When importing/syncing this repository in Overleaf, use the `main` branch and select `report/main.tex` as the main document.
