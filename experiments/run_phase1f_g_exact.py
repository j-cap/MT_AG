"""Run P1F-G with the P1F-D exact-DR conventional-PF comparator.

The main P1F-G experiment module imports ``run_paper_bootstrap_pf`` into its
module namespace. This wrapper supplies the exact paper-level DR increment
settings used in P1F-D before delegating to the experiment entry point. It is
kept explicit so the held-out PF comparator cannot silently fall back to the
non-zero process-noise defaults of ``run_paper_bootstrap_pf``.
"""

from __future__ import annotations

import run_phase1f_g as experiment

_original_pf = experiment.run_paper_bootstrap_pf


def _exact_increment_pf(*args, **kwargs):
    cfg = experiment.yaml.safe_load(
        (experiment.ROOT / "configs/phase1f_g.yaml").read_text()
    )
    kwargs.setdefault("sigma_delta_l_m", cfg["paper_model"]["sigma_delta_l_m"])
    kwargs.setdefault(
        "sigma_delta_phi_rad", cfg["paper_model"]["sigma_delta_phi_rad"]
    )
    return _original_pf(*args, **kwargs)


experiment.run_paper_bootstrap_pf = _exact_increment_pf


if __name__ == "__main__":
    experiment.main()
