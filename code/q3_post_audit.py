"""Finish Q3 auxiliary audits without repeating verified high-grid solves.

Use this only after q3_solver.py has produced a validation file whose core
high-grid gates passed and whose only failed gates are Radau and sensitivity.
The script reruns those moderate-grid checks, revalidates the candidate
workbook against the saved tight-grid fields, and publishes only if every gate
is true.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np

import q3_solver as q3


CORE_GATES = (
    "event",
    "platform",
    "monotonicity",
    "center_grid_crosscheck",
    "strict_event",
    "event_root",
    "workbook",
    "mass_balance",
    "convergence",
    "q2_same_grid_compatibility",
)


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=q3.DEFAULT_OUTPUT)
    parser.add_argument("--publish", action="store_true")
    args = parser.parse_args()

    output = args.output.resolve()
    validation_path = output / "q3_validation.json"
    validation = read_json(validation_path)
    gate = dict(validation.get("gate_summary", {}))
    failed_core = [name for name in CORE_GATES if not bool(gate.get(name, False))]
    if failed_core:
        raise RuntimeError(f"Cannot resume auxiliary audits; failed core gates: {failed_core}")
    if validation.get("convergence", {}).get("status") != "PASSED":
        raise RuntimeError("Saved high-grid convergence evidence is not PASSED.")

    fields_path = output / "q3_fields.npz"
    candidate = output / "result3_candidate.xlsx"
    if not fields_path.exists() or not candidate.exists():
        raise FileNotFoundError("Saved tight-grid fields or candidate workbook is missing.")

    validation["status"] = "AUXILIARY_AUDIT_RUNNING"
    q3._write_json(validation_path, validation)

    boundary, _ = q3.q2.read_boundary(q3.DEFAULT_BOUNDARY, end=q3.PLATEAU_START_S)
    diagnostic = q3.platform_diagnostic(boundary)
    if not diagnostic.get("pass", False):
        raise RuntimeError("Platform diagnostic no longer passes against the authoritative attachment.")

    temporal = validation["convergence"]["temporal"]
    base = temporal["base"]
    n_intervals = int(validation["n_intervals"])
    radau = q3.run_radau_audit(
        boundary,
        diagnostic,
        n=n_intervals,
        max_end_s=q3.DEFAULT_MAX_END_S,
        rtol=float(base["rtol"]),
        atol=float(base["atol"]),
        max_step_before_s=float(base["max_step_before_s"]),
        max_step_after_s=float(base["max_step_after_s"]),
    )

    event_time = float(validation["event"]["time_s"])
    sensitivity_rows = q3.run_sensitivities(
        boundary,
        diagnostic,
        baseline_event_time_s=event_time,
        n=n_intervals,
        max_end_s=q3.DEFAULT_MAX_END_S,
        rtol=float(base["rtol"]),
        atol=float(base["atol"]),
        max_step_before_s=float(base["max_step_before_s"]),
        max_step_after_s=float(base["max_step_after_s"]),
        output=output,
    )
    sensitivity_pass = bool(sensitivity_rows) and all(
        row.get("status") == "EVENT_LOCATED" for row in sensitivity_rows
    )
    sensitivity = {
        "status": "PASSED" if sensitivity_pass else "REVIEW_REQUIRED",
        "rows": sensitivity_rows,
    }

    with np.load(fields_path, allow_pickle=False) as payload:
        times = payload["time_s"].astype(int)
        moisture = payload["moisture_kg_per_kg"].copy()
    minute_end = int(math.floor(event_time / 60.0) * 60)
    mask = (times >= 60) & (times <= minute_end)
    workbook_report = q3.validate_result3_workbook(
        candidate,
        times[mask],
        moisture[mask],
        template=q3.DEFAULT_TEMPLATE,
    )

    validation["independent_radau"] = radau
    validation["sensitivity"] = sensitivity
    validation["workbook_validation"] = workbook_report
    gate["radau"] = radau.get("status") == "PASSED"
    gate["sensitivity"] = sensitivity_pass
    gate["workbook"] = bool(workbook_report.get("four_decimal_match", False))
    validation["gate_summary"] = gate
    validation["status"] = "NUMERICAL_CHECKS_PASSED" if all(bool(value) for value in gate.values()) else "REVIEW_REQUIRED"

    if validation["status"] == "NUMERICAL_CHECKS_PASSED" and args.publish:
        shutil.copy2(candidate, q3.OFFICIAL_RESULT3)
        validation["official_workbook_validation"] = q3.validate_result3_workbook(
            q3.OFFICIAL_RESULT3,
            times[mask],
            moisture[mask],
            template=q3.DEFAULT_TEMPLATE,
        )
    else:
        validation["official_workbook_validation"] = {
            "status": "NOT_OVERWRITTEN",
            "reason": "requires all checks and --publish",
        }
    q3._write_json(validation_path, validation)
    print(
        json.dumps(
            {
                "status": validation["status"],
                "radau": radau,
                "sensitivity_status": sensitivity["status"],
                "official_workbook": str(q3.OFFICIAL_RESULT3),
                "published": bool(validation["status"] == "NUMERICAL_CHECKS_PASSED" and args.publish),
            },
            ensure_ascii=False,
        )
    )
    if validation["status"] != "NUMERICAL_CHECKS_PASSED":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
