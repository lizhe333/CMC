"""Generate the Question 2 paper tables and figure from accepted result files.

This script is deliberately a presentation-only step.  It reads the already
validated paper CSVs (and validation JSON), checks their contract, writes the
two LaTeX tables, and copies the accepted Q2 PDFs into ``paper``.
It never calls the Q2 solver or recomputes a numerical trajectory.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path


EXPECTED_TIMES_H = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
EXPECTED_RADII_CM = (0.0, 0.5, 1.0, 1.5, 2.0)
EXPECTED_HEADERS = ("time_s", "r_0_cm", "r_0.5_cm", "r_1_cm", "r_1.5_cm", "r_2_cm")


def _read_paper_csv(path: Path) -> tuple[list[float], list[list[float]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows or tuple(rows[0]) != EXPECTED_HEADERS:
        raise ValueError(f"unexpected paper CSV header: {path}")
    if len(rows) != 1 + len(EXPECTED_TIMES_H):
        raise ValueError(f"expected six paper rows in {path}, got {len(rows) - 1}")
    times: list[float] = []
    values: list[list[float]] = []
    for row in rows[1:]:
        if len(row) != len(EXPECTED_HEADERS):
            raise ValueError(f"unexpected column count in {path}: {row}")
        times.append(float(row[0]) / 3600.0)
        values.append([float(item) for item in row[1:]])
    for actual, expected in zip(times, EXPECTED_TIMES_H):
        if abs(actual - expected) > 1e-10:
            raise ValueError(f"unexpected paper time {actual} h in {path}")
    return times, values


def _latex_table(
    *,
    caption: str,
    label: str,
    values: list[list[float]],
) -> str:
    lines = [
        r"\begin{table}[H]",
        r"  \centering",
        r"  \fontsize{10.5pt}{12.6pt}\selectfont",
        r"  \renewcommand{\arraystretch}{1.25}",
        r"  \captionsetup{font=small,labelfont=bf,textfont=bf}",
        r"  \caption{" + caption + r"}",
        r"  \label{" + label + r"}",
        r"  \begin{tabularx}{0.92\textwidth}{|>{\centering\arraybackslash}p{0.23\textwidth}|*{5}{>{\centering\arraybackslash}X|}}",
        r"    \hline",
        r"    \multirow{2}{*}{时间/h} & \multicolumn{5}{|c|}{到药材中心的距离/cm} \\",
        r"    \cline{2-6}",
        r"     & 0 & 0.5 & 1 & 1.5 & 2 \\",
        r"    \hline",
    ]
    for hour, row in zip(EXPECTED_TIMES_H, values):
        fields = " & ".join(f"{item:.4f}" for item in row)
        lines.append(f"    {hour:.1f} & {fields} \\\\")
        lines.append(r"    \hline")
    lines.extend(
        [
            r"  \end{tabularx}",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def generate(project_root: Path) -> dict[str, str]:
    results = project_root / "results" / "q2"
    paper = project_root / "paper"
    validation_path = results / "q2_validation.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") != "NUMERICAL_CHECKS_PASSED":
        raise RuntimeError(f"Q2 validation is not accepted: {validation.get('status')}")

    temperature_csv = results / "q2_temperature_paper.csv"
    moisture_csv = results / "q2_moisture_paper.csv"
    times_t, temperature = _read_paper_csv(temperature_csv)
    times_c, moisture = _read_paper_csv(moisture_csv)
    if times_t != times_c:
        raise ValueError("temperature and moisture paper times differ")
    if any(len(row) != len(EXPECTED_RADII_CM) for row in temperature + moisture):
        raise ValueError("paper tables must contain five radial values per time")

    table_dir = paper / "tables"
    figure_dir = paper / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    temperature_path = table_dir / "q2_temperature.tex"
    moisture_path = table_dir / "q2_moisture.tex"
    temperature_path.write_text(
        _latex_table(
            caption=r"3小时内药材的温度",
            label="tab:q2-temperature",
            values=temperature,
        ),
        encoding="utf-8",
    )
    moisture_path.write_text(
        _latex_table(
            caption=r"3小时内药材的水分浓度",
            label="tab:q2-moisture",
            values=moisture,
        ),
        encoding="utf-8",
    )

    source_figure = results / "q2_fields.pdf"
    target_figure = figure_dir / "q2_fields.pdf"
    shutil.copy2(source_figure, target_figure)
    source_spacetime = results / "q2_spacetime.pdf"
    target_spacetime = figure_dir / "q2_spacetime.pdf"
    if not source_spacetime.is_file():
        raise FileNotFoundError(f"missing accepted spacetime figure: {source_spacetime}")
    shutil.copy2(source_spacetime, target_spacetime)
    return {
        "temperature_table": str(temperature_path),
        "moisture_table": str(moisture_path),
        "figure": str(target_figure),
        "spacetime_figure": str(target_spacetime),
        "source_validation": str(validation_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="CUMCM_2026_A project root (default: this script's parent project)",
    )
    args = parser.parse_args()
    print(json.dumps(generate(args.project_root.resolve()), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
