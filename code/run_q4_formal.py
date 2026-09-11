"""Run and audit the formal Question 4 numerical calculations.

The script is deliberately separate from the solver and from the paper and
workbook stages.  A run claims a fresh directory below ``results/q4/formal``
and records all numerical settings and source modification times before any
long integration starts.  Each completed case is saved with the native Q4
solver's unrounded fields, endpoint certificate, flux history and block
evidence.  A run can be resumed explicitly after an interrupted process; an
ordinary invocation never reuses or overwrites a previous run directory.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import platform
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve()
ROOT = HERE.parents[1]
CODE = ROOT / "code"
RESULTS_ROOT = ROOT / "results" / "q4" / "formal"
Q4_MODELING = ROOT / "reports" / "questions" / "Q4_MODELING.md"
Q4_SOLVER = CODE / "q4_solver.py"
Q4_CONTRACTS = CODE / "test_q4_contracts.py"

if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
import q4_solver as q4  # noqa: E402


SPACE_GRIDS = (640, 1280, 2560, 5120, 10240, 20480)
MAX_END_S = q4.MAX_END_S
COMMON_OUTPUT_POSITIONS_CM = np.arange(21, dtype=float) * 0.1
TABLE6_POSITIONS_CM = np.arange(5, dtype=float) * 0.5
FORMAL_THRESHOLDS = {
    "event_abs_s": 30.0,
    "event_relative": 1.0e-4,
    "temperature_field_C": 5.0e-5,
    "moisture_field_kg_per_kg": 5.0e-5,
    "table6_four_decimal": True,
}
REGRESSION_THRESHOLDS = {
    "temperature_C": 1.0e-6,
    "moisture_kg_per_kg": 1.0e-7,
    "event_abs_s": 1.0,
    "event_relative": 1.0e-6,
}
EVENT_THRESHOLDS = {
    "bracket_width_s": 1.0,
    "root_residual_kg_per_kg": 1.0e-9,
    "strict_witness_offset_s": 0.5,
}


def _safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_safe(value), ensure_ascii=False, indent=2), encoding="utf-8")


def write_progress(run_root: Path, stage: str, **details: Any) -> None:
    """Write a small resumable progress marker after each completed stage."""
    write_json(run_root / "progress.json", {
        "run_id": run_root.name,
        "stage": stage,
        "updated_at": now_utc(),
        **details,
    })


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def source_stamp(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path.resolve()),
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "size_bytes": int(stat.st_size),
    }


def new_run_id() -> str:
    return f"q4_formal_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:10]}"


def claim_run(run_id: str | None, resume: bool) -> tuple[str, Path, bool]:
    requested = run_id or new_run_id()
    candidate = (RESULTS_ROOT / requested).resolve()
    if candidate.parent != RESULTS_ROOT.resolve():
        raise ValueError("run_id must remain below results/q4/formal")
    if candidate.exists():
        if not resume:
            raise FileExistsError(f"Formal run directory already exists: {candidate}")
        marker = candidate / "run_config.json"
        if not marker.exists():
            raise FileExistsError("--resume requires an existing formal run_config.json")
        return requested, candidate, True
    # Use the solver's guarded claim routine so the same run id cannot be
    # silently reused by another invocation.
    claimed = q4._claim_evidence_directory(RESULTS_ROOT, requested)
    return requested, claimed, False


def settings(kind: str) -> dict[str, float]:
    if kind == "base":
        return {
            "rtol": q4.BASE_RTL,
            "atol": q4.BASE_ATOL,
            "max_step_before_s": q4.BASE_MAX_STEP_BEFORE_S,
            "max_step_after_s": q4.BASE_MAX_STEP_AFTER_S,
        }
    if kind == "tight":
        return {
            "rtol": q4.TIGHT_RTL,
            "atol": q4.TIGHT_ATOL,
            "max_step_before_s": q4.TIGHT_MAX_STEP_BEFORE_S,
            "max_step_after_s": q4.TIGHT_MAX_STEP_AFTER_S,
        }
    raise ValueError(f"unknown settings kind: {kind}")


def fixed_radius(max_end_s: float) -> tuple[q4.RadiusModel, dict[str, Any]]:
    model = q4.RadiusModel(
        np.array([0.0, max(float(max_end_s), q4.SWITCH_S + 1.0)]),
        np.array([q4.RADIUS_INITIAL_M, q4.RADIUS_INITIAL_M]),
        interpolation="linear",
        post_mode="plateau",
    )
    metadata = {
        "source": "constant initial radius R0",
        "geometry": "fixed_R0",
        "interpolation": "not_applicable",
        "post_mode": "not_applicable",
        "start_radius_m": q4.RADIUS_INITIAL_M,
        "end_radius_m": q4.RADIUS_INITIAL_M,
    }
    return model, metadata


def run_config_for_case(
    *,
    run_id: str,
    stage: str,
    case: str,
    n: int,
    model: str,
    radius_metadata: dict[str, Any],
    radius_interpolation: str,
    radius_post_mode: str,
    cfg: dict[str, float],
    h: float = q4.H_BASE,
    hm: float = q4.HM_BASE,
    plateau_temperature_C: float = q4.PLATEAU_T_C,
    plateau_moisture_kg_per_kg: float = q4.PLATEAU_C,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "run_id": run_id,
        "case": case,
        "backend": "q4_native",
        "property_model": model,
        "geometry": "fixed_R0" if case in {"A", "C"} else "attachment2_radius",
        "radius_interpolation": radius_interpolation,
        "radius_post_mode": radius_post_mode,
        "radius_metadata": radius_metadata,
        "n_intervals": int(n),
        "max_end_s": float(MAX_END_S),
        "environment": {
            "measured_interval_s": [0.0, q4.SWITCH_S],
            "post_switch": {
                "temperature_C": float(plateau_temperature_C),
                "moisture_kg_per_kg": float(plateau_moisture_kg_per_kg),
            },
            "description": "Attachment 1 piecewise linear through 4 h, then Q3 nominal platform unless this is an explicit sensitivity.",
        },
        "exchange": {"h_W_m2K": float(h), "hm_m_s": float(hm)},
        "time_settings": dict(cfg),
        "output": {
            "time_set": "60 s from 60 plus one unique exact endpoint",
            "positions_cm": COMMON_OUTPUT_POSITIONS_CM.tolist(),
            "table6_positions_cm": TABLE6_POSITIONS_CM.tolist(),
        },
        "termination": {
            "criterion": "first full-grid maximum moisture downward crossing",
            "threshold_kg_per_kg": q4.THRESHOLD_C,
            "max_window_s": MAX_END_S,
        },
    }


def result_certificate(result: dict[str, Any]) -> dict[str, Any]:
    event = result.get("event")
    summary = q4.case_summary(result)
    if event is None:
        return {
            "event_present": False,
            "event_time_s": None,
            "event_certification_status": None,
            "full_audit": result.get("g_evidence_summary", {}).get("full_run_sequence_audit"),
            "output_schema": result.get("output_schema"),
            "case_summary": summary,
        }
    full = event.get("full_run_sequence_audit", {})
    return {
        "event_present": True,
        "event_time_s": float(event["time_s"]),
        "event_certification_status": event.get("certification_status"),
        "bracket_s": event.get("bracket_s"),
        "bracket_width_s": float(event.get("bracket_width_s", math.inf)),
        "root_residual_kg_per_kg": float(event.get("root_residual_kg_per_kg", math.inf)),
        "strict_before_time_s": event.get("strict_before_time_s"),
        "strict_after_time_s": event.get("strict_after_time_s"),
        "strict_before_g_N": event.get("strict_before_g_N"),
        "strict_after_g_N": event.get("strict_after_g_N"),
        "strict_before_pass": bool(event.get("strict_before_pass", False)),
        "strict_after_pass": bool(event.get("strict_after_pass", False)),
        "refinement_difference_s": event.get("refinement_difference_s"),
        "full_audit": full,
        "full_prior_positive": bool(event.get("full_prior_positive", False)),
        "full_downward_sign_changes": int(event.get("full_downward_sign_changes", 0)),
        "full_upward_sign_changes": int(event.get("full_upward_sign_changes", 0)),
        "full_no_upward_recrossing": bool(event.get("full_no_upward_recrossing", False)),
        "monitor_recovery_attempts": int(event.get("monitor_recovery_attempts", 0)),
        "output_schema": result.get("output_schema"),
        "case_summary": summary,
    }


def save_case(
    run_root: Path,
    *,
    result: dict[str, Any],
    output_dir: Path,
    evidence_base: Path,
    run_id: str,
    boundary_metadata: dict[str, Any],
    radius_metadata: dict[str, Any],
    stage: str,
    case: str,
    n: int,
    model: str,
    radius_interpolation: str,
    radius_post_mode: str,
    cfg: dict[str, float],
    h: float = q4.H_BASE,
    hm: float = q4.HM_BASE,
    plateau_temperature_C: float = q4.PLATEAU_T_C,
    plateau_moisture_kg_per_kg: float = q4.PLATEAU_C,
) -> dict[str, Any]:
    metadata = run_config_for_case(
        run_id=run_id,
        stage=stage,
        case=case,
        n=n,
        model=model,
        radius_metadata=radius_metadata,
        radius_interpolation=radius_interpolation,
        radius_post_mode=radius_post_mode,
        cfg=cfg,
        h=h,
        hm=hm,
        plateau_temperature_C=plateau_temperature_C,
        plateau_moisture_kg_per_kg=plateau_moisture_kg_per_kg,
    )
    result["formal_stage_config"] = metadata
    q4.save_result_files(
        output_dir,
        result,
        boundary_metadata=boundary_metadata,
        radius_metadata=radius_metadata,
        validation={"formal_stage": stage, "case": case, "run_id": run_id},
    )
    certificate = result_certificate(result)
    record = {
        "run_id": run_id,
        "stage": stage,
        "case": case,
        "n_intervals": int(n),
        "model": model,
        "radius_interpolation": radius_interpolation,
        "radius_post_mode": radius_post_mode,
        "h_W_m2K": float(h),
        "hm_m_s": float(hm),
        "plateau_temperature_C": float(plateau_temperature_C),
        "plateau_moisture_kg_per_kg": float(plateau_moisture_kg_per_kg),
        "settings": cfg,
        "output_directory": str(output_dir.resolve()),
        "evidence_base": str(evidence_base.resolve()),
        "certificate": certificate,
        "runtime_s": float(result.get("runtime_s", math.nan)),
        "saved_at": now_utc(),
    }
    write_json(output_dir / "formal_case_record.json", record)
    return record


def load_fields(path: Path) -> dict[str, np.ndarray]:
    with np.load(path / "q4_fields.npz", allow_pickle=False) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def load_event(path: Path) -> dict[str, Any] | None:
    summary_path = path / "q4_case_summary.json"
    if not summary_path.exists():
        return None
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    return summary.get("event")


def _unique_position_index(positions: np.ndarray, target_cm: float) -> int | None:
    """Return the unique output-grid match, or None for interpolation."""
    positions = np.asarray(positions, dtype=float).reshape(-1)
    matches = np.flatnonzero(np.isclose(positions, float(target_cm), rtol=0.0, atol=1.0e-10))
    if len(matches) > 1:
        raise ValueError(f"Output positions are not uniquely matched at {target_cm} cm.")
    return int(matches[0]) if len(matches) == 1 else None


def _map_position_row(
    values: np.ndarray,
    positions_cm: np.ndarray,
    targets_cm: np.ndarray,
    *,
    valid_mask: np.ndarray | None = None,
    radius_m: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Map one physical-position row, interpolating only inside the domain."""
    values = np.asarray(values, dtype=float).reshape(-1)
    positions_cm = np.asarray(positions_cm, dtype=float).reshape(-1)
    targets_cm = np.asarray(targets_cm, dtype=float).reshape(-1)
    if len(values) != len(positions_cm) or not np.isfinite(positions_cm).all():
        raise ValueError("Position mapping requires a finite value/grid pair.")
    if len(positions_cm) > 1 and not np.all(np.diff(positions_cm) > 0.0):
        raise ValueError("Output positions must be strictly increasing for physical interpolation.")
    if valid_mask is None:
        valid = np.isfinite(values)
    else:
        valid = np.asarray(valid_mask, dtype=bool).reshape(-1) & np.isfinite(values)
    mapped = np.full(len(targets_cm), np.nan, dtype=float)
    mapped_mask = np.zeros(len(targets_cm), dtype=bool)
    radius_tol_cm = 0.0
    if radius_m is not None:
        radius_tol_cm = 100.0 * q4.radius_mask_tolerance(float(radius_m))
    for column, target_cm in enumerate(targets_cm):
        index = _unique_position_index(positions_cm, float(target_cm))
        if index is not None:
            if valid[index]:
                mapped[column] = values[index]
                mapped_mask[column] = True
            continue
        if radius_m is not None and float(target_cm) > 100.0 * float(radius_m) + radius_tol_cm:
            continue
        valid_positions = positions_cm[valid]
        valid_values = values[valid]
        if len(valid_positions) < 2:
            continue
        if float(target_cm) < valid_positions[0] - 1.0e-10 or float(target_cm) > valid_positions[-1] + 1.0e-10:
            continue
        mapped[column] = float(np.interp(float(target_cm), valid_positions, valid_values))
        mapped_mask[column] = True
    return mapped, mapped_mask


def _map_table6_columns(fields: dict[str, np.ndarray], row_indices: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map temperature and moisture to table-6 physical positions together."""
    positions = np.asarray(fields["positions_cm"], dtype=float).reshape(-1)
    temperature = np.asarray(fields["temperature_C"], dtype=float)
    moisture = np.asarray(fields["moisture_kg_per_kg"], dtype=float)
    domain = np.asarray(fields["domain_mask"], dtype=bool)
    radii = np.asarray(fields["radius_m"], dtype=float).reshape(-1)
    row_indices = np.asarray(row_indices, dtype=int).reshape(-1)
    if temperature.shape != moisture.shape or temperature.shape != domain.shape:
        raise ValueError("Table-6 source fields have inconsistent shapes.")
    if temperature.shape[0] != len(radii):
        raise ValueError("Table-6 radii and source rows have inconsistent lengths.")
    mapped_temperature = np.full((len(row_indices), len(TABLE6_POSITIONS_CM)), np.nan, dtype=float)
    mapped_moisture = np.full_like(mapped_temperature, np.nan)
    mapped_domain = np.zeros_like(mapped_temperature, dtype=bool)
    for output_row, source_row in enumerate(row_indices):
        common_valid = domain[source_row] & np.isfinite(temperature[source_row]) & np.isfinite(moisture[source_row])
        mapped_temperature[output_row], temperature_mask = _map_position_row(
            temperature[source_row], positions, TABLE6_POSITIONS_CM,
            valid_mask=common_valid, radius_m=float(radii[source_row]),
        )
        mapped_moisture[output_row], moisture_mask = _map_position_row(
            moisture[source_row], positions, TABLE6_POSITIONS_CM,
            valid_mask=common_valid, radius_m=float(radii[source_row]),
        )
        mapped_domain[output_row] = temperature_mask & moisture_mask
        mapped_temperature[output_row, ~mapped_domain[output_row]] = np.nan
        mapped_moisture[output_row, ~mapped_domain[output_row]] = np.nan
    return mapped_temperature, mapped_moisture, mapped_domain


def _four_decimal_equal(left: Any, right: Any) -> bool:
    """Compare one table-6 value or a row vector after four-decimal rounding."""
    left_array = np.asarray(left, dtype=float)
    right_array = np.asarray(right, dtype=float)
    if left_array.ndim != right_array.ndim or left_array.shape != right_array.shape:
        return False
    return bool(np.array_equal(np.round(left_array, 4), np.round(right_array, 4), equal_nan=True))


def table6_view(fields: dict[str, np.ndarray], endpoint_s: float | None = None) -> dict[str, Any]:
    times = np.asarray(fields["time_s"], dtype=float)
    endpoint = float(times[-1] if endpoint_s is None else endpoint_s)
    wanted = [float(6 * 3600 * i) for i in range(1, int(math.floor(endpoint / (6 * 3600) + 1e-12)) + 1)]
    if not wanted or abs(wanted[-1] - endpoint) > 1e-8:
        wanted.append(endpoint)
    selected: list[int] = []
    selected_times: list[float] = []
    for value in wanted:
        match = np.where(np.isclose(times, value, rtol=0.0, atol=1.0e-7))[0]
        if len(match) > 1:
            raise ValueError(f"Output time {value} s is not uniquely represented.")
        if len(match):
            selected.append(int(match[0]))
            selected_times.append(float(times[match[0]]))
    row_indices = np.asarray(selected, dtype=int)
    T, C, mask = _map_table6_columns(fields, row_indices)
    surface_C = np.asarray(fields["surface_moisture_kg_per_kg"], dtype=float)[selected]
    surface_T = np.asarray(fields["surface_temperature_C"], dtype=float)[selected]
    radius = np.asarray(fields["radius_m"], dtype=float)[selected]
    return {
        "time_s": np.asarray(selected_times, dtype=float),
        "positions_cm": TABLE6_POSITIONS_CM,
        "temperature_C": T,
        "moisture_kg_per_kg": C,
        "surface_temperature_C": surface_T,
        "surface_moisture_kg_per_kg": surface_C,
        "radius_m": radius,
        "domain_mask": mask,
    }


def table6_position_regression() -> dict[str, Any]:
    """Directly check table-6 columns, physical interpolation, and masking."""
    targets = TABLE6_POSITIONS_CM.copy()
    full_positions = COMMON_OUTPUT_POSITIONS_CM.copy()
    full_radii = np.array([0.02000, 0.01250], dtype=float)
    full_temperature = 20.0 + full_positions[None, :]
    full_moisture = 0.20 + 0.01 * full_positions[None, :]
    full_radius_tolerance_m = np.maximum(q4.RADIUS_MASK_ABS_M, q4.RADIUS_MASK_REL * np.abs(full_radii))
    full_domain = full_positions[None, :] <= (100.0 * full_radii[:, None] + 100.0 * full_radius_tolerance_m[:, None])
    full_fields = {
        "positions_cm": full_positions,
        "temperature_C": np.repeat(full_temperature, len(full_radii), axis=0),
        "moisture_kg_per_kg": np.repeat(full_moisture, len(full_radii), axis=0),
        "domain_mask": full_domain,
        "radius_m": full_radii,
    }
    full_T, full_C, full_mask = _map_table6_columns(full_fields, np.array([0, 1], dtype=int))
    direct_expected_T = 20.0 + targets
    direct_expected_C = 0.20 + 0.01 * targets
    direct_pass = bool(
        np.allclose(full_T[0], direct_expected_T, rtol=0.0, atol=1.0e-14)
        and np.allclose(full_C[0], direct_expected_C, rtol=0.0, atol=1.0e-14)
        and np.all(full_mask[0])
        and np.array_equal(full_mask[1], np.array([True, True, True, False, False], dtype=bool))
        and np.isnan(full_T[1, 3:]).all()
        and np.isnan(full_C[1, 3:]).all()
    )

    sparse_positions = np.arange(11, dtype=float) * 0.2
    sparse_radius = np.array([0.01250], dtype=float)
    sparse_fields = {
        "positions_cm": sparse_positions,
        "temperature_C": (10.0 + 3.0 * sparse_positions)[None, :],
        "moisture_kg_per_kg": (0.30 + 0.02 * sparse_positions)[None, :],
        "domain_mask": (sparse_positions[None, :] <= 100.0 * sparse_radius[:, None]),
        "radius_m": sparse_radius,
    }
    sparse_T, sparse_C, sparse_mask = _map_table6_columns(sparse_fields, np.array([0], dtype=int))
    sparse_expected_T = 10.0 + 3.0 * targets
    sparse_expected_C = 0.30 + 0.02 * targets
    interpolation_pass = bool(
        np.allclose(sparse_T[0, :3], sparse_expected_T[:3], rtol=0.0, atol=1.0e-14)
        and np.allclose(sparse_C[0, :3], sparse_expected_C[:3], rtol=0.0, atol=1.0e-14)
        and np.array_equal(sparse_mask[0], np.array([True, True, True, False, False], dtype=bool))
        and np.isnan(sparse_T[0, 3:]).all()
        and np.isnan(sparse_C[0, 3:]).all()
    )
    duplicate_rejected = False
    try:
        _unique_position_index(np.array([0.5, 0.50000000005], dtype=float), 0.5)
    except ValueError:
        duplicate_rejected = True
    scalar_rounding_pass = _four_decimal_equal(0.12344, 0.123441)
    single_array_rounding_pass = _four_decimal_equal(np.array([0.12344]), np.array([0.123441]))
    multi_array_rounding_pass = _four_decimal_equal(
        np.array([0.12344, 0.56784]), np.array([0.123441, 0.567841])
    )
    scalar_array_shape_guard = not _four_decimal_equal(0.12344, np.array([0.12344]))
    passed = bool(
        direct_pass
        and interpolation_pass
        and duplicate_rejected
        and scalar_rounding_pass
        and single_array_rounding_pass
        and multi_array_rounding_pass
        and scalar_array_shape_guard
    )
    return {
        "pass": passed,
        "targets_cm": targets,
        "direct_grid_targets_pass": direct_pass,
        "physical_interpolation_pass": interpolation_pass,
        "duplicate_position_rejected": duplicate_rejected,
        "scalar_rounding_pass": scalar_rounding_pass,
        "single_element_array_rounding_pass": single_array_rounding_pass,
        "multi_element_array_rounding_pass": multi_array_rounding_pass,
        "scalar_array_shape_guard": scalar_array_shape_guard,
        "domain_mask_full_radius": full_mask[0],
        "domain_mask_shrunk_radius": full_mask[1],
        "domain_mask_interpolated": sparse_mask[0],
    }


def compare_records(coarse: dict[str, Any], fine: dict[str, Any], coarse_dir: Path, fine_dir: Path) -> dict[str, Any]:
    a = load_fields(coarse_dir)
    b = load_fields(fine_dir)
    ta = np.asarray(a["time_s"], dtype=float)
    tb = np.asarray(b["time_s"], dtype=float)
    common = np.intersect1d(ta, tb)
    ia = {float(t): i for i, t in enumerate(ta)}
    ib = {float(t): i for i, t in enumerate(tb)}
    if len(common):
        aidx = np.asarray([ia[float(t)] for t in common], dtype=int)
        bidx = np.asarray([ib[float(t)] for t in common], dtype=int)
        Tdiff = np.abs(np.asarray(a["temperature_C"])[aidx] - np.asarray(b["temperature_C"])[bidx])
        Cdiff = np.abs(np.asarray(a["moisture_kg_per_kg"])[aidx] - np.asarray(b["moisture_kg_per_kg"])[bidx])
        masks = np.asarray(a["domain_mask"], dtype=bool)[aidx] & np.asarray(b["domain_mask"], dtype=bool)[bidx]
        Tdiff = np.where(masks, Tdiff, np.nan)
        Cdiff = np.where(masks, Cdiff, np.nan)
        max_T = float(np.nanmax(Tdiff)) if np.isfinite(Tdiff).any() else math.inf
        max_C = float(np.nanmax(Cdiff)) if np.isfinite(Cdiff).any() else math.inf
    else:
        max_T = max_C = math.inf
    ea = load_event(coarse_dir)
    eb = load_event(fine_dir)
    ta_event = None if ea is None or ea.get("time_s") is None else float(ea["time_s"])
    tb_event = None if eb is None or eb.get("time_s") is None else float(eb["time_s"])
    event_diff = math.inf if ta_event is None or tb_event is None else abs(ta_event - tb_event)
    event_relative = math.inf if ta_event is None or ta_event == 0.0 else event_diff / abs(ta_event)
    va = table6_view(a, ta_event)
    vb = table6_view(b, tb_event)
    table6_common = np.intersect1d(va["time_s"], vb["time_s"])
    table6_stable = True
    table6_rows = 0
    if len(table6_common):
        va_i = {float(t): i for i, t in enumerate(va["time_s"])}
        vb_i = {float(t): i for i, t in enumerate(vb["time_s"])}
        for t in table6_common:
            i, j = va_i[float(t)], vb_i[float(t)]
            for key in ("moisture_kg_per_kg", "surface_moisture_kg_per_kg"):
                x = np.asarray(va[key])[i]
                y = np.asarray(vb[key])[j]
                if key == "moisture_kg_per_kg":
                    mask = np.asarray(va["domain_mask"])[i] & np.asarray(vb["domain_mask"])[j]
                    x, y = x[mask], y[mask]
                if not _four_decimal_equal(x, y):
                    table6_stable = False
            table6_rows += 1
    endpoint_round_stable = None
    if ea is not None and eb is not None:
        # At their own exact endpoints the two grids need not share time, so
        # compare the root's mapped table-6 values rather than an interpolated
        # rounded output row.
        try:
            with np.load(coarse_dir / "q4_event.npz", allow_pickle=False) as xa, np.load(fine_dir / "q4_event.npz", allow_pickle=False) as xb:
                ca = np.asarray(xa["root_mapped_C"], dtype=float)
                cb = np.asarray(xb["root_mapped_C"], dtype=float)
                ca_mask = np.asarray(xa["root_domain_mask"], dtype=bool)
                cb_mask = np.asarray(xb["root_domain_mask"], dtype=bool)
                ca6, ca6_mask = _map_position_row(
                    ca, np.asarray(a["positions_cm"], dtype=float), TABLE6_POSITIONS_CM,
                    valid_mask=ca_mask, radius_m=float(xa["root_radius_m"]),
                )
                cb6, cb6_mask = _map_position_row(
                    cb, np.asarray(b["positions_cm"], dtype=float), TABLE6_POSITIONS_CM,
                    valid_mask=cb_mask, radius_m=float(xb["root_radius_m"]),
                )
                common_mask = ca6_mask & cb6_mask
                table6_endpoint_pass = bool(
                    np.array_equal(ca6_mask, cb6_mask)
                    and np.array_equal(np.round(ca6[common_mask], 4), np.round(cb6[common_mask], 4))
                )
                ca_surface = float(np.asarray(xa["moisture_kg_per_kg"], dtype=float)[-1])
                cb_surface = float(np.asarray(xb["moisture_kg_per_kg"], dtype=float)[-1])
                endpoint_round_stable = bool(table6_endpoint_pass and np.round(ca_surface, 4) == np.round(cb_surface, 4))
        except Exception:
            endpoint_round_stable = False
    pass_flag = bool(
        len(common) > 0
        and event_diff <= FORMAL_THRESHOLDS["event_abs_s"]
        and event_relative <= FORMAL_THRESHOLDS["event_relative"]
        and max_T <= FORMAL_THRESHOLDS["temperature_field_C"]
        and max_C <= FORMAL_THRESHOLDS["moisture_field_kg_per_kg"]
        and table6_stable
        and (endpoint_round_stable is not False)
    )
    return {
        "coarse_n": int(coarse.get("n_intervals", -1)),
        "fine_n": int(fine.get("n_intervals", -1)),
        "coarse_event_time_s": ta_event,
        "fine_event_time_s": tb_event,
        "event_abs_difference_s": event_diff,
        "event_relative_difference": event_relative,
        "common_rows": int(len(common)),
        "common_end_s": float(common[-1]) if len(common) else None,
        "temperature_max_abs_difference_C": max_T,
        "moisture_max_abs_difference_kg_per_kg": max_C,
        "table6_common_rows": table6_rows,
        "table6_four_decimal_stable_common_rows": table6_stable,
        "table6_four_decimal_stable_endpoints": endpoint_round_stable,
        "thresholds": FORMAL_THRESHOLDS,
        "pass": pass_flag,
    }


def event_and_audit_pass(record: dict[str, Any]) -> bool:
    cert = record.get("certificate", {})
    full = cert.get("full_audit") or {}
    return bool(
        cert.get("event_present")
        and cert.get("event_certification_status") == "PASS"
        and float(cert.get("bracket_width_s", math.inf)) <= EVENT_THRESHOLDS["bracket_width_s"]
        and abs(float(cert.get("root_residual_kg_per_kg", math.inf))) <= EVENT_THRESHOLDS["root_residual_kg_per_kg"]
        and cert.get("strict_before_pass")
        and cert.get("strict_after_pass")
        and full.get("status") == "CHECKED"
        and bool(full.get("prior_positive"))
        and int(full.get("downward_sign_changes", 0)) == 1
        and int(full.get("upward_sign_changes", 0)) == 0
        and bool(full.get("no_upward_recrossing"))
    )


def record_path(stage_dir: Path, name: str) -> Path:
    return stage_dir / name / "formal_case_record.json"


def run_one_case(
    *,
    run_root: Path,
    out_dir: Path,
    evidence_base: Path,
    run_id: str,
    boundary: np.ndarray,
    boundary_metadata: dict[str, Any],
    radius: q4.RadiusModel,
    radius_metadata: dict[str, Any],
    n: int,
    model: str,
    case: str,
    stage: str,
    cfg: dict[str, float],
    radius_interpolation: str,
    radius_post_mode: str,
    h: float = q4.H_BASE,
    hm: float = q4.HM_BASE,
    plateau_temperature_C: float = q4.PLATEAU_T_C,
    plateau_moisture_kg_per_kg: float = q4.PLATEAU_C,
) -> dict[str, Any]:
    marker = out_dir / "formal_case_record.json"
    if marker.exists():
        return json.loads(marker.read_text(encoding="utf-8"))
    out_dir.mkdir(parents=True, exist_ok=False)
    run = q4.run_case(
        boundary,
        radius,
        n=n,
        model=model,
        max_end_s=MAX_END_S,
        rtol=cfg["rtol"],
        atol=cfg["atol"],
        max_step_before_s=cfg["max_step_before_s"],
        max_step_after_s=cfg["max_step_after_s"],
        chunk_before_s=q4.DEFAULT_CHUNK_S,
        chunk_after_s=q4.DEFAULT_CHUNK_S,
        h=h,
        hm=hm,
        plateau_temperature_C=plateau_temperature_C,
        plateau_moisture_kg_per_kg=plateau_moisture_kg_per_kg,
        case_name=case,
        geometry_label="fixed_R0" if case in {"A", "C"} else "attachment2_radius",
        radius_metadata=radius_metadata,
        evidence_dir=evidence_base,
        run_id=run_id,
    )
    return save_case(
        run_root,
        result=run,
        output_dir=out_dir,
        evidence_base=evidence_base,
        run_id=run_id,
        boundary_metadata=boundary_metadata,
        radius_metadata=radius_metadata,
        stage=stage,
        case=case,
        n=n,
        model=model,
        radius_interpolation=radius_interpolation,
        radius_post_mode=radius_post_mode,
        cfg=cfg,
        h=h,
        hm=hm,
        plateau_temperature_C=plateau_temperature_C,
        plateau_moisture_kg_per_kg=plateau_moisture_kg_per_kg,
    )


def run_d_space(
    run_root: Path,
    boundary: np.ndarray,
    boundary_metadata: dict[str, Any],
    radius: q4.RadiusModel,
    radius_metadata: dict[str, Any],
    *,
    resume: bool,
) -> tuple[int, dict[str, Any], list[dict[str, Any]]]:
    stage = run_root / "space"
    evidence = stage / "evidence"
    stage.mkdir(parents=True, exist_ok=True)
    cfg = settings("tight")
    records: list[dict[str, Any]] = []
    previous: dict[str, Any] | None = None
    previous_dir: Path | None = None
    comparisons: list[dict[str, Any]] = []
    selected_final: int | None = None
    for n in SPACE_GRIDS:
        name = f"D_n{n:05d}"
        out = stage / name
        run_id = f"{name}_tight"
        record = run_one_case(
            run_root=run_root,
            out_dir=out,
            evidence_base=evidence,
            run_id=run_id,
            boundary=boundary,
            boundary_metadata=boundary_metadata,
            radius=radius,
            radius_metadata=radius_metadata,
            n=n,
            model="appendix4",
            case="D",
            stage="space_convergence",
            cfg=cfg,
            radius_interpolation=getattr(radius, "interpolation", "pchip"),
            radius_post_mode=getattr(radius, "post_mode", "plateau"),
        )
        records.append(record)
        if not event_and_audit_pass(record):
            write_json(stage / "space_failure.json", {"n_intervals": n, "record": record, "reason": "event_or_full_audit_failed"})
            raise RuntimeError(f"Formal D space case failed event/full audit at N={n}")
        if previous is not None and previous_dir is not None:
            comp = compare_records(previous, record, previous_dir, out)
            comparisons.append(comp)
            write_json(stage / f"compare_{int(previous['n_intervals'])}_{n}.json", comp)
            if n >= 2560 and comp["pass"]:
                selected_final = n
        previous = record
        previous_dir = out
        # The first three grids are mandatory.  Continue through 10240 as the
        # planned candidate sequence; only the final pair decides whether the
        # 20480 fallback is required.
        if n == 10240:
            selected_final = n if comparisons and comparisons[-1]["pass"] else None
            if selected_final is not None:
                break
    if selected_final is None:
        # The last comparison failed the publication thresholds.  N=20480 is
        # the prescribed final attempt; do not lower the thresholds.
        n = 20480
        if not records or records[-1]["n_intervals"] != n:
            name = f"D_n{n:05d}"
            out = stage / name
            record = run_one_case(
                run_root=run_root,
                out_dir=out,
                evidence_base=evidence,
                run_id=f"{name}_tight",
                boundary=boundary,
                boundary_metadata=boundary_metadata,
                radius=radius,
                radius_metadata=radius_metadata,
                n=n,
                model="appendix4",
                case="D",
                stage="space_convergence",
                cfg=cfg,
                radius_interpolation=getattr(radius, "interpolation", "pchip"),
                radius_post_mode=getattr(radius, "post_mode", "plateau"),
            )
            records.append(record)
            if not event_and_audit_pass(record):
                write_json(stage / "space_failure.json", {"n_intervals": n, "record": record, "reason": "event_or_full_audit_failed"})
                raise RuntimeError("Formal D space case failed event/full audit at N=20480")
            if previous is not None and previous_dir is not None:
                comp = compare_records(previous, record, previous_dir, out)
                comparisons.append(comp)
                write_json(stage / f"compare_{int(previous['n_intervals'])}_{n}.json", comp)
        if not comparisons or not comparisons[-1]["pass"]:
            write_json(stage / "space_failure.json", {"records": records, "comparisons": comparisons, "reason": "20480_did_not_meet_formal_thresholds"})
            raise RuntimeError("D space convergence did not meet formal thresholds through N=20480")
        selected_final = n
    summary = {
        "status": "PASS",
        "selected_final_n": int(selected_final),
        "mandatory_grids": [int(r["n_intervals"]) for r in records],
        "comparisons": comparisons,
        "settings": cfg,
        "thresholds": FORMAL_THRESHOLDS,
        "endpoint_time_s": float(records[-1]["certificate"]["event_time_s"]),
    }
    write_json(stage / "space_convergence.json", summary)
    return int(selected_final), summary, records


def run_time_convergence(
    run_root: Path,
    boundary: np.ndarray,
    boundary_metadata: dict[str, Any],
    radius: q4.RadiusModel,
    radius_metadata: dict[str, Any],
    final_n: int,
    *,
    space_records: list[dict[str, Any]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    stage = run_root / "time"
    evidence = stage / "evidence"
    stage.mkdir(parents=True, exist_ok=True)
    base = run_one_case(
        run_root=run_root,
        out_dir=stage / f"D_n{final_n:05d}_base",
        evidence_base=evidence,
        run_id=f"D_n{final_n:05d}_base",
        boundary=boundary,
        boundary_metadata=boundary_metadata,
        radius=radius,
        radius_metadata=radius_metadata,
        n=final_n,
        model="appendix4",
        case="D",
        stage="time_base",
        cfg=settings("base"),
        radius_interpolation=getattr(radius, "interpolation", "pchip"),
        radius_post_mode=getattr(radius, "post_mode", "plateau"),
    )
    # The tight solve at the selected final space grid is already the space
    # record and is reused exactly; no duplicate long integration is made.
    tight = next(r for r in space_records if int(r["n_intervals"]) == int(final_n))
    base_dir = stage / f"D_n{final_n:05d}_base"
    tight_dir = run_root / "space" / f"D_n{final_n:05d}"
    comparison = compare_records(base, tight, base_dir, tight_dir)
    comparison["base_settings"] = settings("base")
    comparison["tight_settings"] = settings("tight")
    comparison["selected_final_n"] = int(final_n)
    comparison["pass"] = bool(
        comparison["event_abs_difference_s"] <= FORMAL_THRESHOLDS["event_abs_s"]
        and comparison["event_relative_difference"] <= FORMAL_THRESHOLDS["event_relative"]
        and comparison["temperature_max_abs_difference_C"] <= FORMAL_THRESHOLDS["temperature_field_C"]
        and comparison["moisture_max_abs_difference_kg_per_kg"] <= FORMAL_THRESHOLDS["moisture_field_kg_per_kg"]
        and comparison["table6_four_decimal_stable_common_rows"]
        and comparison["table6_four_decimal_stable_endpoints"] is not False
    )
    write_json(stage / "base_tight_comparison.json", comparison)
    if not comparison["pass"]:
        raise RuntimeError("Formal base/tight time convergence failed at selected grid")
    summary = {"status": "PASS", "comparison": comparison, "base": base, "tight": tight}
    write_json(stage / "time_convergence.json", summary)
    return base, summary


def run_ad_cases(
    run_root: Path,
    boundary: np.ndarray,
    boundary_metadata: dict[str, Any],
    radius: q4.RadiusModel,
    radius_metadata: dict[str, Any],
    fixed: q4.RadiusModel,
    fixed_metadata: dict[str, Any],
    final_n: int,
    *,
    space_records: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    stage = run_root / "cases"
    evidence = stage / "evidence"
    stage.mkdir(parents=True, exist_ok=True)
    cfg = settings("tight")
    results: dict[str, dict[str, Any]] = {}
    # Reuse the certified final D space result as the D factorial case.
    d_record = next(r for r in space_records if int(r["n_intervals"]) == int(final_n))
    results["D"] = d_record
    for case, model, geom, rad, meta, interp, post in (
        ("A", "appendix3", "fixed_R0", fixed, fixed_metadata, "not_applicable", "not_applicable"),
        ("B", "appendix3", "attachment2_radius", radius, radius_metadata, getattr(radius, "interpolation", "pchip"), getattr(radius, "post_mode", "plateau")),
        ("C", "appendix4", "fixed_R0", fixed, fixed_metadata, "not_applicable", "not_applicable"),
    ):
        out = stage / f"case_{case}"
        record = run_one_case(
            run_root=run_root,
            out_dir=out,
            evidence_base=evidence,
            run_id=f"case_{case}_tight_n{final_n:05d}",
            boundary=boundary,
            boundary_metadata=boundary_metadata,
            radius=rad,
            radius_metadata=meta,
            n=final_n,
            model=model,
            case=case,
            stage="factorial_tight",
            cfg=cfg,
            radius_interpolation=interp,
            radius_post_mode=post,
        )
        results[case] = record
    # Copy a small index rather than duplicating the D binary fields.
    write_json(stage / "case_index.json", {"final_n": final_n, "cases": results})
    if not all(event_and_audit_pass(row) for row in results.values()):
        raise RuntimeError("At least one A-D factorial case failed event/full audit")
    times = {case: float(row["certificate"]["event_time_s"]) for case, row in results.items()}
    effects = {
        "event_time_s": times,
        "delta_property_s": times["C"] - times["A"],
        "delta_shrink_3_s": times["B"] - times["A"],
        "delta_shrink_4_s": times["D"] - times["C"],
        "delta_interaction_s": (times["D"] - times["C"]) - (times["B"] - times["A"]),
        "status": "PASS",
    }
    write_json(stage / "factorial_differences.json", effects)
    return results


def q3_regression(
    run_root: Path,
    boundary: np.ndarray,
    final_n: int,
    fixed: q4.RadiusModel,
    *,
    q4_a_record: dict[str, Any],
) -> dict[str, Any]:
    """Compare final A with the existing Q3 solver as an external reference."""
    import q3_solver as q3  # noqa: WPS433,E402

    stage = run_root / "audit"
    stage.mkdir(parents=True, exist_ok=True)
    diagnostic = q3.platform_diagnostic(boundary)
    if not diagnostic.get("pass", False):
        raise RuntimeError("Q3 platform diagnostic did not pass")
    cfg = settings("tight")
    started = time.perf_counter()
    q3_result = q3.run_case(
        boundary,
        diagnostic,
        n=final_n,
        max_end_s=MAX_END_S,
        rtol=cfg["rtol"],
        atol=cfg["atol"],
        max_step_before_s=cfg["max_step_before_s"],
        max_step_after_s=cfg["max_step_after_s"],
        chunk_before_s=int(q4.DEFAULT_CHUNK_S),
        chunk_after_s=int(q4.DEFAULT_CHUNK_S),
    )
    q4_fields = load_fields(run_root / "cases" / "case_A")
    q3_times = np.asarray(q3_result["times_s"], dtype=float)
    q4_times = np.asarray(q4_fields["time_s"], dtype=float)
    common = np.intersect1d(q3_times, q4_times)
    q3_i = {float(t): i for i, t in enumerate(q3_times)}
    q4_i = {float(t): i for i, t in enumerate(q4_times)}
    q3_T = np.asarray([q3_result["temperature_C"][q3_i[float(t)]] for t in common], dtype=float)
    q4_T = np.asarray([q4_fields["temperature_C"][q4_i[float(t)]] for t in common], dtype=float)
    q3_C = np.asarray([q3_result["moisture_kg_per_kg"][q3_i[float(t)]] for t in common], dtype=float)
    q4_C = np.asarray([q4_fields["moisture_kg_per_kg"][q4_i[float(t)]] for t in common], dtype=float)
    max_T = float(np.max(np.abs(q3_T - q4_T))) if len(common) else math.inf
    max_C = float(np.max(np.abs(q3_C - q4_C))) if len(common) else math.inf
    q3_event = q3_result.get("event")
    q4_event = load_event(run_root / "cases" / "case_A")
    event_difference = math.inf if q3_event is None or q4_event is None else abs(float(q3_event["time_s"]) - float(q4_event["time_s"]))
    event_relative = math.inf if q3_event is None else event_difference / abs(float(q3_event["time_s"]))
    # Q3 stores a full state at its root; Q4's event NPZ stores the same.
    endpoint_T = endpoint_C = math.inf
    if q3_event is not None and q4_event is not None:
        with np.load(run_root / "cases" / "case_A" / "q4_event.npz", allow_pickle=False) as data:
            q4_event_T = np.asarray(data["temperature_C"], dtype=float)
            q4_event_C = np.asarray(data["moisture_kg_per_kg"], dtype=float)
        endpoint_T = float(np.max(np.abs(np.asarray(q3_event["state"])[0::2] - q4_event_T)))
        endpoint_C = float(np.max(np.abs(np.asarray(q3_event["state"])[1::2] - q4_event_C)))
    report = {
        "status": "PASS" if max_T <= REGRESSION_THRESHOLDS["temperature_C"] and max_C <= REGRESSION_THRESHOLDS["moisture_kg_per_kg"] and event_difference <= REGRESSION_THRESHOLDS["event_abs_s"] and event_relative <= REGRESSION_THRESHOLDS["event_relative"] and endpoint_T <= REGRESSION_THRESHOLDS["temperature_C"] and endpoint_C <= REGRESSION_THRESHOLDS["moisture_kg_per_kg"] else "FAIL",
        "reference": "external read-only q3_solver",
        "n_intervals": int(final_n),
        "common_rows": int(len(common)),
        "common_end_s": float(common[-1]) if len(common) else None,
        "max_temperature_abs_difference_C": max_T,
        "max_moisture_abs_difference_kg_per_kg": max_C,
        "q3_event_time_s": None if q3_event is None else float(q3_event["time_s"]),
        "q4_event_time_s": None if q4_event is None else float(q4_event["time_s"]),
        "event_abs_difference_s": event_difference,
        "event_relative_difference": event_relative,
        "endpoint_full_state_temperature_abs_difference_C": endpoint_T,
        "endpoint_full_state_moisture_abs_difference_kg_per_kg": endpoint_C,
        "thresholds": REGRESSION_THRESHOLDS,
        "q4_case_record": q4_a_record,
        "runtime_s": float(time.perf_counter() - started),
    }
    write_json(stage / "q3_fixed_radius_regression.json", report)
    if report["status"] != "PASS":
        raise RuntimeError("Formal A/Q3 regression failed")
    return report


def independent_flux_audit(run_root: Path, boundary: np.ndarray, radius: q4.RadiusModel) -> dict[str, Any]:
    """Reassemble every control volume through explicit loops."""
    stage = run_root / "audit"
    stage.mkdir(parents=True, exist_ok=True)
    n = 1000
    dxi, nodes, faces, weights = q4.radial_geometry(n)
    measured = q4.measured_environment(boundary)
    rhs_log, _jac, _geometry = q4.make_log_rhs(measured, radius, n, model="appendix4")
    solution = solve_ivp(
        rhs_log,
        (10800.0, 12600.0),
        q4._to_log_state(q4.initial_state(n), n, model="appendix4"),
        method="BDF",
        rtol=q4.TIGHT_RTL,
        atol=q4.TIGHT_ATOL,
        max_step=q4.TIGHT_MAX_STEP_BEFORE_S,
        dense_output=True,
    )
    # The short interval starts from the initial state only for an isolated
    # flux assembly test; the main RHS and independent RHS are evaluated on
    # identical nonuniform states.  A second continuation from 0 is used for
    # the physical audit below.
    first = q4.integrate_segment(
        q4.initial_state(n), n, 0.0, 10800.0, measured, radius,
        model="appendix4", rtol=q4.TIGHT_RTL, atol=q4.TIGHT_ATOL,
        max_step=q4.TIGHT_MAX_STEP_BEFORE_S, chunk_seconds=q4.DEFAULT_CHUNK_S,
        output_times=np.empty(0), detect_event=False,
    )
    continuation = solve_ivp(
        rhs_log,
        (10800.0, 12600.0),
        q4._to_log_state(first["state_final"], n, model="appendix4"),
        method="BDF",
        jac_sparsity=q4.jacobian_sparsity(n),
        rtol=q4.TIGHT_RTL,
        atol=q4.TIGHT_ATOL,
        max_step=q4.TIGHT_MAX_STEP_BEFORE_S,
        dense_output=True,
    )
    # Use q4's public physical RHS only as the comparison target.  The
    # independent implementation below intentionally contains no vectorized
    # face slices from the main solver.
    main_rhs, _ = q4.make_rhs(measured, radius, n, model="appendix4")

    def independent(state: np.ndarray, t: float) -> np.ndarray:
        T = np.empty(n + 1, dtype=float)
        C = np.empty(n + 1, dtype=float)
        for i in range(n + 1):
            T[i] = float(state[2 * i])
            C[i] = float(state[2 * i + 1])
        _rho, _cp, k, D = q4.material_properties(C, T, model="appendix4")
        T_env, C_env = measured(float(t))
        R = float(radius(float(t)))
        heat_flux = np.zeros(n, dtype=float)
        moisture_flux = np.zeros(n, dtype=float)
        for face_index in range(n):
            face = float(faces[face_index])
            k_face = 0.5 * (float(k[face_index]) + float(k[face_index + 1]))
            D_face = 0.5 * (float(D[face_index]) + float(D[face_index + 1]))
            heat_flux[face_index] = face * k_face * (T[face_index + 1] - T[face_index]) / dxi
            moisture_flux[face_index] = face * D_face * (C[face_index + 1] - C[face_index]) / dxi
        heat_surface = -R * q4.H_BASE * (T[-1] - T_env)
        moisture_surface = -R * q4.HM_BASE * (C[-1] - C_env)
        result = np.empty(2 * (n + 1), dtype=float)
        for i in range(n + 1):
            if i == 0:
                hdiv = heat_flux[0] / float(weights[0])
                mdiv = moisture_flux[0] / float(weights[0])
            elif i == n:
                hdiv = (heat_surface - heat_flux[n - 1]) / float(weights[n])
                mdiv = (moisture_surface - moisture_flux[n - 1]) / float(weights[n])
            else:
                hdiv = (heat_flux[i] - heat_flux[i - 1]) / float(weights[i])
                mdiv = (moisture_flux[i] - moisture_flux[i - 1]) / float(weights[i])
            result[2 * i] = hdiv / (R * R * float(_rho[i]) * float(_cp[i]))
            result[2 * i + 1] = mdiv / (R * R)
        return result

    samples = np.asarray([10800.0 + 300.0 * i for i in range(7)], dtype=float)
    max_T = 0.0
    max_C = 0.0
    for t in samples:
        z = np.asarray(continuation.sol(float(t))).reshape(-1)
        physical, _ = q4._physical_from_log(z, n)
        main = main_rhs(float(t), physical)
        alt = independent(physical, float(t))
        max_T = max(max_T, float(np.max(np.abs(main[0::2] - alt[0::2]))))
        max_C = max(max_C, float(np.max(np.abs(main[1::2] - alt[1::2]))))
    report = {
        "status": "PASS" if max_T <= 1.0e-10 and max_C <= 1.0e-12 else "FAIL",
        "implementation": "explicit per-control-face and per-control-volume loop",
        "n_intervals": n,
        "interval_s": [10800.0, 12600.0],
        "sample_times_s": samples,
        "max_temperature_abs_difference_C_per_s": max_T,
        "max_moisture_abs_difference_kg_per_kg_s": max_C,
        "thresholds": {"temperature_C_per_s": 1.0e-10, "moisture_kg_per_kg_s": 1.0e-12},
        "isolated_short_solution_success": bool(solution.success),
        "continuation_success": bool(continuation.success),
    }
    write_json(stage / "independent_flux_audit.json", report)
    if report["status"] != "PASS":
        raise RuntimeError("Independent control-volume flux audit failed")
    return report


def radau_audit(run_root: Path, boundary: np.ndarray, radius: q4.RadiusModel) -> dict[str, Any]:
    stage = run_root / "audit"
    stage.mkdir(parents=True, exist_ok=True)
    n = 1000
    first = q4.integrate_segment(
        q4.initial_state(n), n, 0.0, q4.SWITCH_S, q4.measured_environment(boundary), radius,
        model="appendix4", rtol=q4.TIGHT_RTL, atol=q4.TIGHT_ATOL,
        max_step=q4.TIGHT_MAX_STEP_BEFORE_S, chunk_seconds=q4.DEFAULT_CHUNK_S,
        output_times=np.empty(0), detect_event=False,
    )
    plateau = q4.plateau_environment()
    rhs, jac, _ = q4.make_log_rhs(plateau, radius, n, model="appendix4")
    initial = q4._to_log_state(first["state_final"], n, model="appendix4")
    common = {
        "rtol": q4.TIGHT_RTL,
        "atol": q4.TIGHT_ATOL,
        "max_step": q4.TIGHT_MAX_STEP_AFTER_S,
        "dense_output": True,
        "jac_sparsity": q4.jacobian_sparsity(n),
    }
    bdf = solve_ivp(rhs, (q4.SWITCH_S, 18000.0), initial, method="BDF", jac=jac, **common)
    radau = solve_ivp(rhs, (q4.SWITCH_S, 18000.0), initial, method="Radau", jac=jac, **common)
    sample = np.asarray([14400.0, 15300.0, 16200.0, 17100.0, 18000.0], dtype=float)
    bdf_state = q4._physical_from_transformed_matrix(bdf.sol(sample), n)
    radau_state = q4._physical_from_transformed_matrix(radau.sol(sample), n)
    dT = float(np.max(np.abs(bdf_state[0::2, :] - radau_state[0::2, :])))
    dC = float(np.max(np.abs(bdf_state[1::2, :] - radau_state[1::2, :])))
    report = {
        "status": "PASS" if bdf.success and radau.success and dT <= 5.0e-7 and dC <= 5.0e-8 else "FAIL",
        "n_intervals": n,
        "interval_s": [q4.SWITCH_S, 18000.0],
        "sample_times_s": sample,
        "bdf_success": bool(bdf.success),
        "radau_success": bool(radau.success),
        "max_temperature_abs_difference_C": dT,
        "max_moisture_abs_difference_kg_per_kg": dC,
        "thresholds": {"temperature_C": 5.0e-7, "moisture_kg_per_kg": 5.0e-8},
        "bdf_nfev": int(bdf.nfev),
        "radau_nfev": int(radau.nfev),
    }
    write_json(stage / "radau_audit.json", report)
    if report["status"] != "PASS":
        raise RuntimeError("BDF/Radau representative interval failed")
    return report


def run_sensitivity(
    run_root: Path,
    boundary: np.ndarray,
    boundary_metadata: dict[str, Any],
    radius: q4.RadiusModel,
    radius_metadata: dict[str, Any],
    final_n: int,
    baseline_time: float,
) -> dict[str, Any]:
    stage = run_root / "sensitivity"
    evidence = stage / "evidence"
    stage.mkdir(parents=True, exist_ok=True)
    cfg = settings("tight")
    cases: list[tuple[str, q4.RadiusModel, dict[str, Any], float, float, float, float, str, str]] = []
    for name, h in (("h_minus20", 0.8 * q4.H_BASE), ("h_plus20", 1.2 * q4.H_BASE)):
        cases.append((name, radius, radius_metadata, h, q4.HM_BASE, q4.PLATEAU_T_C, q4.PLATEAU_C, getattr(radius, "interpolation", "pchip"), getattr(radius, "post_mode", "plateau")))
    for name, hm in (("hm_minus20", 0.8 * q4.HM_BASE), ("hm_plus20", 1.2 * q4.HM_BASE)):
        cases.append((name, radius, radius_metadata, q4.H_BASE, hm, q4.PLATEAU_T_C, q4.PLATEAU_C, getattr(radius, "interpolation", "pchip"), getattr(radius, "post_mode", "plateau")))
    linear, linear_metadata = q4.read_radius(q4.DEFAULT_RADIUS, interpolation="linear", post_mode="plateau")
    cases.append(("radius_linear", linear, linear_metadata, q4.H_BASE, q4.HM_BASE, q4.PLATEAU_T_C, q4.PLATEAU_C, "linear", "plateau"))
    # Q3's two long-term boundary scenarios are retained as explicit
    # sensitivities; the nominal 50 C/0.05 platform remains the model used by
    # A-D and by the selected endpoint.
    import q3_solver as q3  # noqa: WPS433,E402
    diagnostic = q3.platform_diagnostic(boundary)
    scenarios = [
        ("q3_long_final_hour_mean", float(diagnostic["temperature_C"]["mean"]), float(diagnostic["moisture_kg_per_kg"]["mean"])),
        ("q3_long_last_observation", float(boundary[-1, 1]), float(boundary[-1, 2])),
    ]
    for name, temp, moisture in scenarios:
        cases.append((name, radius, radius_metadata, q4.H_BASE, q4.HM_BASE, temp, moisture, getattr(radius, "interpolation", "pchip"), getattr(radius, "post_mode", "plateau")))
    records: list[dict[str, Any]] = []
    for name, rad, meta, h, hm, temp, moisture, interp, post in cases:
        record = run_one_case(
            run_root=run_root,
            out_dir=stage / name,
            evidence_base=evidence,
            run_id=f"{name}_n{final_n:05d}",
            boundary=boundary,
            boundary_metadata=boundary_metadata,
            radius=rad,
            radius_metadata=meta,
            n=final_n,
            model="appendix4",
            case="D",
            stage="sensitivity",
            cfg=cfg,
            radius_interpolation=interp,
            radius_post_mode=post,
            h=h,
            hm=hm,
            plateau_temperature_C=temp,
            plateau_moisture_kg_per_kg=moisture,
        )
        if not event_and_audit_pass(record):
            raise RuntimeError(f"Sensitivity case failed event/full audit: {name}")
        event_time = float(record["certificate"]["event_time_s"])
        record["delta_time_s_vs_nominal_D"] = event_time - baseline_time
        records.append(record)
    output = {
        "status": "PASS",
        "n_intervals": int(final_n),
        "baseline_event_time_s": float(baseline_time),
        "settings": cfg,
        "cases": records,
        "interpretation": "Each row is compared with the nominal D result at the same grid and time settings; boundary scenarios are not refits.",
    }
    write_json(stage / "sensitivity.json", output)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=None, help="fresh run id; existing ids require --resume")
    parser.add_argument("--resume", action="store_true", help="continue an existing run without replacing completed case directories")
    parser.add_argument("--skip-q3-regression", action="store_true")
    parser.add_argument("--skip-sensitivity", action="store_true")
    parser.add_argument("--table6-regression", action="store_true", help="run the short physical-position table-6 mapping contract")
    args = parser.parse_args(argv)
    table6_regression = table6_position_regression()
    if args.table6_regression:
        print(json.dumps(_safe(table6_regression), ensure_ascii=False))
        return 0 if table6_regression.get("pass") is True else 1
    if table6_regression.get("pass") is not True:
        raise RuntimeError(f"Table-6 position regression failed: {table6_regression}")
    run_id, run_root, resumed = claim_run(args.run_id, args.resume)
    boundary, boundary_metadata = q4.read_boundary()
    radius, radius_metadata = q4.read_radius(interpolation="pchip", post_mode="plateau")
    fixed, fixed_metadata = fixed_radius(MAX_END_S)
    contract_path = ROOT / "results" / "q4" / "dev" / "q4_dev_validation.json"
    if not contract_path.exists():
        raise FileNotFoundError(f"Required formal contract evidence is missing: {contract_path}")
    contract_evidence = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract_evidence.get("run_mode") != "formal_gate" or contract_evidence.get("pass") is not True:
        raise RuntimeError("The required q4 test_q4_contracts.py --formal gate is not PASS")
    if not resumed:
        write_json(run_root / "run_config.json", {
            "run_id": run_id,
            "run_root": str(run_root),
            "started_at": now_utc(),
            "host": platform.node(),
            "python": sys.version,
            "command": " ".join(sys.argv),
            "source": {"solver": source_stamp(Q4_SOLVER), "contracts": source_stamp(Q4_CONTRACTS), "modeling": source_stamp(Q4_MODELING)},
            "inputs": {"boundary": boundary_metadata, "radius": radius_metadata},
            "space_grids": SPACE_GRIDS,
            "max_end_s": MAX_END_S,
            "formal_thresholds": FORMAL_THRESHOLDS,
            "event_thresholds": EVENT_THRESHOLDS,
            "regression_thresholds": REGRESSION_THRESHOLDS,
            "table6_position_regression": table6_regression,
            "status": "RUNNING",
            "resume_supported": True,
            "no_workbook_or_paper_writes": True,
        })
    try:
        final_n, space_summary, space_records = run_d_space(run_root, boundary, boundary_metadata, radius, radius_metadata, resume=args.resume)
        write_progress(run_root, "space_convergence", selected_final_n=final_n, grids=space_summary.get("mandatory_grids"))
        base_record, time_summary = run_time_convergence(run_root, boundary, boundary_metadata, radius, radius_metadata, final_n, space_records=space_records)
        write_progress(run_root, "time_convergence", selected_final_n=final_n, pass_status=time_summary.get("comparison", {}).get("pass"))
        cases = run_ad_cases(run_root, boundary, boundary_metadata, radius, radius_metadata, fixed, fixed_metadata, final_n, space_records=space_records)
        write_progress(run_root, "factorial_cases", cases=list(cases))
        factorial = json.loads((run_root / "cases" / "factorial_differences.json").read_text(encoding="utf-8"))
        q3_report = None
        if not args.skip_q3_regression:
            q3_report = q3_regression(run_root, boundary, final_n, fixed, q4_a_record=cases["A"])
            write_progress(run_root, "q3_regression", status=q3_report.get("status"))
        d_time = float(cases["D"]["certificate"]["event_time_s"])
        sensitivity = None if args.skip_sensitivity else run_sensitivity(run_root, boundary, boundary_metadata, radius, radius_metadata, final_n, d_time)
        write_progress(run_root, "sensitivity", status=None if sensitivity is None else sensitivity.get("status"))         
        radau = radau_audit(run_root, boundary, radius)
        flux = independent_flux_audit(run_root, boundary, radius)
        write_progress(run_root, "audits", radau=radau.get("status"), independent_flux=flux.get("status"))
        validation = {
            "status": "PASS",
            "run_id": run_id,
            "run_root": str(run_root),
            "finished_at": now_utc(),
            "formal_contract": {"path": str(ROOT / "results" / "q4" / "dev" / "q4_dev_validation.json"), "status": "PASS", "required_before_formal": True},
            "space": space_summary,
            "time": time_summary,
            "factorial": factorial,
            "cases": cases,
            "q3_regression": q3_report,
            "sensitivity": sensitivity,
            "radau": radau,
            "independent_flux": flux,
            "final_grid_n": int(final_n),
            "D_event_time_s": d_time,
            "D_event_time_h": d_time / 3600.0,
            "slow1pct_required": bool(d_time > q4.RADIUS_MEASURED_END_S),
            "artifacts_pending": ["q4_post_audit.py", "generate_q4_paper_assets.py"],
        }
        write_json(run_root / "q4_formal_validation.json", validation)
        # The marker is intentionally written last so a partially completed
        # run is never mistaken for a publishable formal result.
        write_json(run_root / "run_status.json", {"status": "PASS", "run_id": run_id, "finished_at": validation["finished_at"], "final_grid_n": final_n})
        write_json(run_root / "run_config.json", {
            **json.loads((run_root / "run_config.json").read_text(encoding="utf-8")),
            "status": "PASS",
            "finished_at": validation["finished_at"],
            "final_grid_n": int(final_n),
            "D_event_time_s": d_time,
        })
        print(json.dumps(_safe({"status": "PASS", "run_id": run_id, "final_grid_n": final_n, "D_event_time_s": d_time}), ensure_ascii=False))
        return 0
    except Exception as exc:
        failure = {"status": "FAIL", "run_id": run_id, "failed_at": now_utc(), "error": repr(exc)}
        write_json(run_root / "run_status.json", failure)
        print(json.dumps(_safe(failure), ensure_ascii=False), file=sys.stderr)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
