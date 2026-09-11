"""Independent read-back audit for a completed formal Q4 run.

This script does not integrate the model and does not rewrite solver output.
It checks the immutable-looking evidence produced by ``run_q4_formal.py``:
case summaries, output schemas, endpoint certificates, block NPZ files,
cross-block state continuity, and the formal comparison statuses.  The result
is written beside the run as ``post_audit.json``.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


REQUIRED_BLOCK_KEYS = {
    "accepted_time_s",
    "accepted_temperature_C",
    "accepted_moisture_kg_per_kg",
    "accepted_radius_m",
    "accepted_g_N",
    "midpoint_time_s",
    "midpoint_temperature_C",
    "midpoint_moisture_kg_per_kg",
    "midpoint_radius_m",
    "midpoint_g_N",
    "output_time_s",
    "output_g_N",
}


def safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(safe(value), ensure_ascii=False, indent=2), encoding="utf-8")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def latest_run(base: Path) -> Path:
    candidates = [p for p in base.glob("q4_formal_*") if p.is_dir() and (p / "q4_formal_validation.json").exists()]
    if not candidates:
        raise FileNotFoundError(f"No completed formal run below {base}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def audit_blocks(case_dir: Path) -> dict[str, Any]:
    manifest_path = case_dir / "q4_blocks_manifest.json"
    if not manifest_path.exists():
        return {"status": "FAIL", "reason": "missing_q4_blocks_manifest"}
    manifest = read_json(manifest_path)
    records = list(manifest.get("segment_records", []))
    if not records:
        return {"status": "FAIL", "reason": "empty_segment_records"}
    block_rows: list[dict[str, Any]] = []
    previous: dict[str, Any] | None = None
    max_temperature_jump = 0.0
    max_moisture_jump = 0.0
    max_radius_jump = 0.0
    # ``integrate_segment`` numbers blocks from zero separately in the two
    # environment segments.  Chronological start time, rather than the local
    # block index, is therefore the only valid global ordering key.
    ordered_records = sorted(records, key=lambda r: (float(r.get("start_s", math.inf)), str(r.get("segment", "")), int(r.get("block_index", -1))))
    for record in ordered_records:
        evidence_name = record.get("evidence_file")
        if not evidence_name:
            return {"status": "FAIL", "reason": "block_record_has_no_evidence_file", "record": record}
        evidence_path = Path(evidence_name)
        if not evidence_path.exists():
            return {"status": "FAIL", "reason": "missing_block_evidence", "path": str(evidence_path)}
        try:
            with np.load(evidence_path, allow_pickle=False) as block:
                keys = set(block.files)
                missing = sorted(REQUIRED_BLOCK_KEYS - keys)
                if missing:
                    return {"status": "FAIL", "reason": "block_schema_missing_keys", "missing": missing, "path": str(evidence_path)}
                accepted_t = np.asarray(block["accepted_time_s"], dtype=float)
                accepted_T = np.asarray(block["accepted_temperature_C"], dtype=float)
                accepted_C = np.asarray(block["accepted_moisture_kg_per_kg"], dtype=float)
                accepted_R = np.asarray(block["accepted_radius_m"], dtype=float)
                accepted_g = np.asarray(block["accepted_g_N"], dtype=float)
                midpoint_t = np.asarray(block["midpoint_time_s"], dtype=float)
                midpoint_T = np.asarray(block["midpoint_temperature_C"], dtype=float)
                midpoint_C = np.asarray(block["midpoint_moisture_kg_per_kg"], dtype=float)
                midpoint_R = np.asarray(block["midpoint_radius_m"], dtype=float)
                midpoint_g = np.asarray(block["midpoint_g_N"], dtype=float)
                output_t = np.asarray(block["output_time_s"], dtype=float)
                output_g = np.asarray(block["output_g_N"], dtype=float)
        except Exception as exc:  # preserve the exact read-back failure
            return {"status": "FAIL", "reason": "block_npz_read_failed", "path": str(evidence_path), "error": repr(exc)}
        arrays = (accepted_t, accepted_T, accepted_C, accepted_R, accepted_g, midpoint_t, midpoint_T, midpoint_C, midpoint_R, midpoint_g, output_t, output_g)
        if any(not np.isfinite(array).all() for array in arrays):
            return {"status": "FAIL", "reason": "nonfinite_block_evidence", "path": str(evidence_path)}
        if len(accepted_t) == 0 or accepted_T.shape[0] != len(accepted_t) or accepted_C.shape[0] != len(accepted_t) or len(accepted_R) != len(accepted_t) or len(accepted_g) != len(accepted_t):
            return {"status": "FAIL", "reason": "accepted_block_shape_mismatch", "path": str(evidence_path)}
        if len(midpoint_t) != len(midpoint_T) or len(midpoint_t) != len(midpoint_C) or len(midpoint_t) != len(midpoint_R) or len(midpoint_t) != len(midpoint_g):
            return {"status": "FAIL", "reason": "midpoint_block_shape_mismatch", "path": str(evidence_path)}
        if len(output_t) != len(output_g):
            return {"status": "FAIL", "reason": "output_monitor_shape_mismatch", "path": str(evidence_path)}
        if np.any(np.diff(accepted_t) <= 0.0):
            return {"status": "FAIL", "reason": "accepted_times_not_strictly_increasing", "path": str(evidence_path)}
        row = {
            "block_index": int(record.get("block_index", -1)),
            "segment": record.get("segment"),
            "path": str(evidence_path),
            "start_s": float(accepted_t[0]),
            "end_s": float(accepted_t[-1]),
            "accepted_states": int(len(accepted_t)),
            "midpoint_states": int(len(midpoint_t)),
            "output_monitor_states": int(len(output_t)),
        }
        if previous is not None:
            dT = float(np.max(np.abs(previous["temperature"][-1] - accepted_T[0])))
            dC = float(np.max(np.abs(previous["moisture"][-1] - accepted_C[0])))
            dR = abs(float(previous["radius"][-1]) - float(accepted_R[0]))
            dt = abs(float(previous["time"][-1]) - float(accepted_t[0]))
            row.update({"shared_time_gap_s": dt, "shared_temperature_max_abs": dT, "shared_moisture_max_abs": dC, "shared_radius_abs_m": dR})
            max_temperature_jump = max(max_temperature_jump, dT)
            max_moisture_jump = max(max_moisture_jump, dC)
            max_radius_jump = max(max_radius_jump, dR)
            if dt > 1.0e-7 or dT > 1.0e-9 or dC > 1.0e-12 or dR > 1.0e-14:
                return {"status": "FAIL", "reason": "adjacent_block_state_discontinuity", "row": row}
        block_rows.append(row)
        previous = {"time": accepted_t, "temperature": accepted_T, "moisture": accepted_C, "radius": accepted_R}
    switch_rows = [row for row in block_rows if abs(float(row["start_s"]) - 14400.0) <= 1.0e-7 or abs(float(row["end_s"]) - 14400.0) <= 1.0e-7]
    return {
        "status": "PASS",
        "block_count": len(block_rows),
        "blocks": block_rows,
        "max_adjacent_temperature_abs": max_temperature_jump,
        "max_adjacent_moisture_abs": max_moisture_jump,
        "max_adjacent_radius_abs_m": max_radius_jump,
        "switch_14400_rows": switch_rows,
    }


def audit_case(case_dir: Path) -> dict[str, Any]:
    summary_path = case_dir / "q4_case_summary.json"
    fields_path = case_dir / "q4_fields.npz"
    event_path = case_dir / "q4_event.npz"
    if not summary_path.exists() or not fields_path.exists():
        return {"status": "FAIL", "reason": "case_summary_or_fields_missing", "case_dir": str(case_dir)}
    summary = read_json(summary_path)
    record_path = case_dir / "formal_case_record.json"
    formal_record = read_json(record_path) if record_path.exists() else {}
    with np.load(fields_path, allow_pickle=False) as fields:
        arrays = {key: np.asarray(fields[key]) for key in fields.files}
    required_fields = {"time_s", "positions_cm", "temperature_C", "moisture_kg_per_kg", "surface_moisture_kg_per_kg", "surface_temperature_C", "radius_m", "mean_moisture", "domain_mask"}
    missing = sorted(required_fields - set(arrays))
    if missing:
        return {"status": "FAIL", "reason": "field_schema_missing_keys", "missing": missing}
    times = arrays["time_s"]
    schema = (
        summary.get("output_schema")
        or formal_record.get("certificate", {}).get("output_schema")
        or summary.get("run_config", {}).get("output_schema")
    )
    if len(times) == 0 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0.0):
        return {"status": "FAIL", "reason": "invalid_output_time_sequence"}
    event = summary.get("event")
    if event is None or event.get("certification_status") != "PASS":
        return {"status": "FAIL", "reason": "event_not_certified", "event": event}
    full = event.get("full_run_sequence_audit") or {}
    event_ok = bool(
        full.get("status") == "CHECKED"
        and full.get("prior_positive")
        and int(full.get("downward_sign_changes", 0)) == 1
        and int(full.get("upward_sign_changes", 0)) == 0
        and full.get("no_upward_recrossing")
        and float(event.get("bracket_width_s", math.inf)) <= 1.0
        and abs(float(event.get("root_residual_kg_per_kg", math.inf))) <= 1.0e-9
        and event.get("strict_before_pass")
        and event.get("strict_after_pass")
    )
    block_audit = audit_blocks(case_dir)
    status = "PASS" if event_ok and block_audit.get("status") == "PASS" and (schema is None or schema.get("pass", False)) else "FAIL"
    return {
        "status": status,
        "case_dir": str(case_dir),
        "event_ok": event_ok,
        "event_time_s": float(event["time_s"]),
        "event_certification_status": event.get("certification_status"),
        "full_run_sequence_audit": full,
        "output_schema": schema,
        "field_shapes": {key: list(value.shape) for key, value in arrays.items()},
        "block_audit": block_audit,
        "event_file_present": event_path.exists(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=None, help="formal run directory; omitted means the latest completed run")
    parser.add_argument("--base", type=Path, default=None, help="results/q4/formal directory used for latest-run discovery")
    args = parser.parse_args(argv)
    base = (args.base or Path(__file__).resolve().parents[1] / "results" / "q4" / "formal").resolve()
    run_root = (args.run_root.resolve() if args.run_root else latest_run(base))
    failures: list[str] = []
    validation_path = run_root / "q4_formal_validation.json"
    config_path = run_root / "run_config.json"
    if not validation_path.exists() or not config_path.exists():
        failures.append("formal_validation_or_run_config_missing")
        validation: dict[str, Any] = {}
        config: dict[str, Any] = {}
    else:
        validation = read_json(validation_path)
        config = read_json(config_path)
        if validation.get("status") != "PASS":
            failures.append("formal_validation_status_not_PASS")
        if config.get("status") != "PASS":
            failures.append("run_config_status_not_PASS")
    case_audits: dict[str, Any] = {}
    case_locations: dict[str, Path] = {}
    if validation:
        cases = validation.get("cases") or {}
        for case in ("A", "B", "C", "D"):
            record = cases.get(case) or {}
            location = Path(record.get("output_directory", "")) if record.get("output_directory") else run_root / "cases" / f"case_{case}"
            case_locations[case] = location
            audit = audit_case(location)
            case_audits[case] = audit
            if audit.get("status") != "PASS":
                failures.append(f"case_{case}_readback_failed")
    space = validation.get("space") or {}
    time_check = validation.get("time") or {}
    if space.get("status") != "PASS" or not all(item.get("pass", False) for item in space.get("comparisons", [])):
        failures.append("space_convergence_not_PASS")
    if time_check.get("status") != "PASS" or not (time_check.get("comparison") or {}).get("pass", False):
        failures.append("time_convergence_not_PASS")
    if (validation.get("factorial") or {}).get("status") != "PASS":
        failures.append("factorial_not_PASS")
    for section, label in (("q3_regression", "q3_regression"), ("sensitivity", "sensitivity"), ("radau", "radau"), ("independent_flux", "independent_flux")):
        if validation.get(section) is None or validation[section].get("status") != "PASS":
            failures.append(f"{label}_not_PASS")
    report = {
        "status": "PASS" if not failures else "FAIL",
        "run_id": run_root.name,
        "run_root": str(run_root),
        "audited_at": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(timespec="seconds"),
        "failures": failures,
        "case_audits": case_audits,
        "space_status": space.get("status"),
        "time_status": time_check.get("status"),
        "factorial_status": (validation.get("factorial") or {}).get("status"),
        "q3_regression_status": (validation.get("q3_regression") or {}).get("status"),
        "sensitivity_status": (validation.get("sensitivity") or {}).get("status"),
        "radau_status": (validation.get("radau") or {}).get("status"),
        "independent_flux_status": (validation.get("independent_flux") or {}).get("status"),
        "formal_validation_source": str(validation_path),
        "run_config_source": str(config_path),
    }
    write_json(run_root / "post_audit.json", report)
    print(json.dumps(safe(report), ensure_ascii=False))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
