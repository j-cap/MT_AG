from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"


def _tex_files():
    yield REPORT / "main.tex"
    yield REPORT / "summary_main.tex"
    yield from sorted((REPORT / "sections").glob("*.tex"))
    yield from sorted((REPORT / "summary_sections").glob("*.tex"))


def test_report_inputs_exist():
    pattern = re.compile(r"\\input\{([^}]+)\}")
    for entry in [REPORT / "main.tex", REPORT / "summary_main.tex"]:
        text = entry.read_text()
        for match in pattern.finditer(text):
            target = REPORT / f"{match.group(1)}.tex"
            assert target.exists(), f"{entry.name} references missing input {target}"


def test_report_graphics_exist():
    pattern = re.compile(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}")
    for tex in _tex_files():
        text = tex.read_text()
        for match in pattern.finditer(text):
            target = REPORT / match.group(1)
            assert target.exists(), f"{tex} references missing figure {target}"


def test_condensed_report_reuses_shared_figures():
    summary_text = "\n".join(
        path.read_text() for path in sorted((REPORT / "summary_sections").glob("*.tex"))
    )
    assert "summary_figures" not in summary_text
    assert "figures/phase1/" in summary_text


def test_detailed_report_no_longer_duplicates_visual_summary():
    main_text = (REPORT / "main.tex").read_text()
    assert "03j_phase1_visual_summary" not in main_text
