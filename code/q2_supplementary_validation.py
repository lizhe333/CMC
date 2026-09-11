"""Supplementary Q2 validation: constitutive degradation and flux balance.

The module intentionally contains computation and artifact serialization only;
plotting is provided by ``plot_q2_supplementary_validation.py``.  The two experiments
write under ``results/q2_verification`` and never modify the formal Q2 result,
Excel workbook, or paper sources.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import q1_solver as q1
import q2_solver as q2


ROOT = CODE_DIR.parent
DEFAULT_BOUNDARY = ROOT.parent / "附件1.xlsx"
DEFAULT_OUTPUT = ROOT / "results" / "q2_verification"
DEGRADATION_GRIDS = (320, 640, 1280)
DEGRADATION_END_S = 1800
BALANCE_N = 10240
BALANCE_END_S = 10800
DEGRADATION_Q1_BASELINE_RTOL = 1e-11
DEGRADATION_Q1_REFINED_RTOL = 1e-12
DEGRADATION_Q1_INTERNAL_MAX_STEP = 5.0  # q1.solve_bdf's fixed internal setting
DEGRADATION_Q2_BASELINE_RTOL = 1e-11
DEGRADATION_Q2_BASELINE_ATOL = 1e-13
DEGRADATION_Q2_BASELINE_MAX_STEP = 2.5
DEGRADATION_Q2_REFINED_RTOL = 1e-12
DEGRADATION_Q2_REFINED_ATOL = 1e-14
DEGRADATION_Q2_REFINED_MAX_STEP = 1.25
BALANCE_Q2_BASELINE_RTOL = 1e-10
BALANCE_Q2_BASELINE_ATOL = 1e-12
BALANCE_Q2_BASELINE_MAX_STEP = 2.5
BALANCE_Q2_REFINED_RTOL = 1e-11
BALANCE_Q2_REFINED_ATOL = 1e-13
BALANCE_Q2_REFINED_MAX_STEP = 1.25
BALANCE_CANDIDATE_ORDERS = (4, 8, 16, 32)
BALANCE_QUADRATURE_THRESHOLD = 1e-9
FIELD_THRESHOLD_T_C = 5e-5
FIELD_THRESHOLD_C = 5e-5
DEGRADATION_TIME_THRESHOLD_T_C = 1e-7
DEGRADATION_TIME_THRESHOLD_C = 1e-8
BALANCE_THRESHOLD_C = 1e-7
BALANCE_TIME_PLATFORM_THRESHOLD = 1e-9
MEAN_CROSSCHECK_THRESHOLD = 1e-10
TOTAL_VOLUME = q2.RADIUS_M * q2.RADIUS_M / 2.0


def _stable_reproduction_commands() -> tuple[str, ...]:
    """Return the complete four-command validation/plotting recipe."""
    python = Path(sys.executable).resolve().as_posix()
    return (
        f"& '{python}' ./code/q2_supplementary_validation.py degradation",
        f"& '{python}' ./code/q2_supplementary_validation.py balance",
        f"& '{python}' ./code/q2_supplementary_validation.py all",
        f"& '{python}' ./code/plot_q2_supplementary_validation.py all",
    )


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_safe(payload), ensure_ascii=False, indent=2), encoding="utf-8"
    )


def refresh_combined_summary(output: Path) -> dict[str, Any]:
    """Combine completed experiment JSON files without invoking either solver."""
    combined: dict[str, Any] = {}
    for key, filename in (
        ("degradation", "q2_degradation.json"),
        ("balance", "q2_balance.json"),
    ):
        path = output / filename
        if path.exists():
            combined[key] = json.loads(path.read_text(encoding="utf-8"))
    _write_json(output / "q2_supplementary_validation.json", combined)
    return combined


def _environment_payload(boundary_path: Path, output: Path) -> dict[str, Any]:
    return {
        "created_by": str(Path(__file__).resolve()),
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "boundary": str(boundary_path.resolve()),
        "output": str(output.resolve()),
        "source_snapshot": str(
            (output / "source_snapshot" / "q2_solver.py").resolve()
        ),
        "formal_results_untouched": True,
        "paper_untouched": True,
    }


def _location(diff: np.ndarray, times: np.ndarray, radius_cm: np.ndarray) -> dict[str, Any]:
    index = np.unravel_index(int(np.argmax(diff)), diff.shape)
    return {
        "max_abs": float(np.max(diff)),
        "time_s": int(times[index[0]]),
        "radius_cm": float(radius_cm[index[1]]),
        "time_index": int(index[0]),
        "radius_index": int(index[1]),
    }


def _field_comparison(
    q1_temperature: np.ndarray,
    q2_temperature: np.ndarray,
    q1_moisture: np.ndarray,
    q2_moisture: np.ndarray,
    times: np.ndarray,
    radius_cm: np.ndarray,
    paper_times: np.ndarray | None = None,
    paper_radius_cm: np.ndarray | None = None,
    temperature_threshold: float = FIELD_THRESHOLD_T_C,
    moisture_threshold: float = FIELD_THRESHOLD_C,
) -> dict[str, Any]:
    temperature_difference = np.abs(q1_temperature - q2_temperature)
    moisture_difference = np.abs(q1_moisture - q2_moisture)
    if paper_times is None:
        paper_times = np.array([1800], dtype=int) if times[-1] < 10800 else q2.PAPER_SECONDS
    if paper_radius_cm is None:
        paper_radius_cm = np.array([0.0, 2.0]) if times[-1] < 10800 else np.arange(5) * 0.5
    time_mask = np.isin(times, paper_times)
    radius_mask = np.isin(radius_cm, paper_radius_cm)
    if not time_mask.any() or not radius_mask.any():
        raise ValueError("Paper comparison points are absent from the supplied fields.")
    paper_temperature = temperature_difference[np.ix_(time_mask, radius_mask)]
    paper_moisture = moisture_difference[np.ix_(time_mask, radius_mask)]
    return {
        "temperature": _location(temperature_difference, times, radius_cm),
        "moisture": _location(moisture_difference, times, radius_cm),
        "temperature_mean_abs": float(np.mean(temperature_difference)),
        "moisture_mean_abs": float(np.mean(moisture_difference)),
        "all_output_rounded_equal": bool(
            np.array_equal(np.round(q1_temperature, 4), np.round(q2_temperature, 4))
            and np.array_equal(np.round(q1_moisture, 4), np.round(q2_moisture, 4))
        ),
        "paper_times_s": np.asarray(paper_times, dtype=int),
        "paper_radius_cm": np.asarray(paper_radius_cm, dtype=float),
        "paper_temperature_max_abs_C": float(np.max(paper_temperature)),
        "paper_moisture_max_abs_kg_per_kg": float(np.max(paper_moisture)),
        "paper_rounded_equal": bool(
            np.array_equal(
                np.round(q1_temperature[np.ix_(time_mask, radius_mask)], 4),
                np.round(q2_temperature[np.ix_(time_mask, radius_mask)], 4),
            )
            and np.array_equal(
                np.round(q1_moisture[np.ix_(time_mask, radius_mask)], 4),
                np.round(q2_moisture[np.ix_(time_mask, radius_mask)], 4),
            )
        ),
        "temperature_threshold_C": float(temperature_threshold),
        "moisture_threshold_kg_per_kg": float(moisture_threshold),
        "temperature_threshold_pass": bool(np.max(temperature_difference) < temperature_threshold),
        "moisture_threshold_pass": bool(np.max(moisture_difference) < moisture_threshold),
        "pass": bool(
            np.max(temperature_difference) < temperature_threshold
            and np.max(moisture_difference) < moisture_threshold
        ),
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty supplementary CSV.")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(_json_safe(rows))


def _append_next_degradation_grid(
    planned_grids: list[int], current_grid: int, comparison_pass: bool
) -> list[int]:
    """Append at most one allowed finer grid after a failed finest attempt."""
    updated = list(planned_grids)
    if comparison_pass or current_grid >= 10240 or current_grid != updated[-1]:
        return updated
    next_grid = current_grid * 2
    if next_grid <= 10240 and next_grid not in updated:
        updated.append(next_grid)
    return updated


def _degradation_gate(
    grid_records: list[dict[str, Any]], q1_time_pass: bool, q2_time_pass: bool
) -> bool:
    """Gate only the finest attempted cross-program comparison and time checks."""
    if not grid_records:
        return False
    return bool(grid_records[-1]["comparison"]["pass"] and q1_time_pass and q2_time_pass)


def _balance_time_gate(
    baseline_selected: dict[str, Any],
    refined_selected: dict[str, Any],
) -> dict[str, Any]:
    """Apply platform exception or strict reduction beyond quadrature floors."""
    baseline_integer = float(baseline_selected["integer_seconds"]["max_abs_B"])
    baseline_accepted = float(baseline_selected["accepted_nodes"]["max_abs_B"])
    refined_integer = float(refined_selected["integer_seconds"]["max_abs_B"])
    refined_accepted = float(refined_selected["accepted_nodes"]["max_abs_B"])
    baseline_joint = max(baseline_integer, baseline_accepted)
    refined_joint = max(refined_integer, refined_accepted)
    baseline_floor = float(baseline_selected["quadrature_floor_max_abs_B"] or 0.0)
    refined_floor = float(refined_selected["quadrature_floor_max_abs_B"] or 0.0)
    floor_uncertainty = baseline_floor + refined_floor
    platform_pass = bool(
        max(baseline_joint, refined_joint) <= BALANCE_TIME_PLATFORM_THRESHOLD
    )
    strict_drop = baseline_joint - refined_joint
    reduction_pass = bool(strict_drop > floor_uncertainty)
    return {
        "platform_threshold": BALANCE_TIME_PLATFORM_THRESHOLD,
        "baseline_integer_max_abs_B": baseline_integer,
        "baseline_accepted_max_abs_B": baseline_accepted,
        "baseline_joint_max_abs_B": baseline_joint,
        "refined_integer_max_abs_B": refined_integer,
        "refined_accepted_max_abs_B": refined_accepted,
        "refined_joint_max_abs_B": refined_joint,
        "baseline_quadrature_floor_max_abs_B": baseline_floor,
        "refined_quadrature_floor_max_abs_B": refined_floor,
        "floor_uncertainty_sum": floor_uncertainty,
        "strict_drop": strict_drop,
        "strict_drop_required": floor_uncertainty,
        "platform_pass": platform_pass,
        "strict_reduction_pass": reduction_pass,
        "pass": bool(platform_pass or reduction_pass),
        "formula": (
            "platform: max(B_base_integer,B_base_accepted,B_ref_integer,B_ref_accepted) "
            "<=1e-9; otherwise require B_base_joint-B_ref_joint > "
            "floor_base+floor_ref"
        ),
    }


def _balance_gate(
    baseline_selection: dict[str, Any],
    refined_selection: dict[str, Any],
    baseline_selected: dict[str, Any],
    refined_selected: dict[str, Any],
    mean_crosscheck_pass: bool,
    field_time_pass: bool,
    time_balance_pass: bool,
) -> bool:
    """Pure final balance gate used by the report and contract tests."""
    return bool(
        baseline_selection["pass"]
        and refined_selection["pass"]
        and baseline_selected["integer_seconds"]["pass"]
        and baseline_selected["accepted_nodes"]["pass"]
        and refined_selected["integer_seconds"]["pass"]
        and refined_selected["accepted_nodes"]["pass"]
        and mean_crosscheck_pass
        and field_time_pass
        and time_balance_pass
    )


def _annotate_q1_checks(checks: dict[str, Any]) -> dict[str, Any]:
    """Make the scope of Q1's same spatial discretization label explicit."""
    annotated = dict(checks)
    annotated["method_scope"] = (
        "The same-spatial-discretization phrase describes only Q1's internal "
        "heat/moisture grid; it does not assert cross-program boundary or "
        "thermal-surface discretization identity."
    )
    annotated["cross_program_discretization_scope"] = (
        "Q1 and Q2 retain their original spatial and boundary discretizations; "
        "their thermal surface treatments differ."
    )
    return annotated


def run_degradation(
    boundary_path: Path,
    output: Path,
    grids: tuple[int, ...] = DEGRADATION_GRIDS,
    end: int = DEGRADATION_END_S,
) -> dict[str, Any]:
    """Compare Q1's BDF with Q2's q1 constitutive degradation.

    Each implementation retains its own original spatial and boundary
    discretisation, thermal surface treatment, and BDF assembly.  The
    comparison is an independent reproduction check, not a claim that either
    constitutive model is physically superior.
    """
    boundary_q1, q1_sheet = q1.read_boundary(boundary_path, end=end)
    boundary_q2, q2_sheet = q2.read_boundary(boundary_path, end=end)
    np.testing.assert_allclose(boundary_q1, boundary_q2, atol=0.0, rtol=0.0)
    times = np.arange(end + 1, dtype=int)
    radius_cm = np.arange(q2.OUTPUT_NODES, dtype=float) / 10.0
    rows: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    planned_grids = list(dict.fromkeys(int(grid) for grid in grids))
    if not planned_grids:
        raise ValueError("At least one degradation grid is required.")
    baseline_cache: dict[int, dict[str, Any]] = {}
    started = time.perf_counter()
    grid_index = 0
    while grid_index < len(planned_grids):
        n = planned_grids[grid_index]
        if n < 20 or n % 20:
            raise ValueError(f"Invalid degradation grid N={n}.")
        q1_started = time.perf_counter()
        q1_t, q1_c, q1_balance, q1_checks = q1.solve_bdf(
            boundary_q1, n=n, rtol=DEGRADATION_Q1_BASELINE_RTOL, end=end
        )
        q1_checks = _annotate_q1_checks(q1_checks)
        q1_runtime = time.perf_counter() - q1_started
        q2_started = time.perf_counter()
        q2_result = q2.integrate_case(
            boundary_q2,
            n=n,
            end=end,
            rtol=DEGRADATION_Q2_BASELINE_RTOL,
            atol=DEGRADATION_Q2_BASELINE_ATOL,
            max_step=DEGRADATION_Q2_BASELINE_MAX_STEP,
            property_model="q1",
        )
        q2_runtime = time.perf_counter() - q2_started
        comparison = _field_comparison(
            q1_t,
            q2_result["temperature_C"],
            q1_c,
            q2_result["moisture_kg_per_kg"],
            times,
            radius_cm,
        )
        baseline_cache[n] = {
            "q1_temperature_C": q1_t,
            "q1_moisture_kg_per_kg": q1_c,
            "q1_balance": q1_balance,
            "q1_checks": q1_checks,
            "q2_result": q2_result,
        }
        record = {
            "n_intervals": n,
            "dr_m": q2.RADIUS_M / n,
            "end_s": end,
            "q1_settings": {
                "solver": "q1.solve_bdf",
                "rtol": DEGRADATION_Q1_BASELINE_RTOL,
                "atol": DEGRADATION_Q1_BASELINE_RTOL * 0.01,
                "max_step_s": DEGRADATION_Q1_INTERNAL_MAX_STEP,
                "sheet": q1_sheet,
            },
            "q2_degraded_settings": {
                "solver": "q2.integrate_case(property_model='q1')",
                "rtol": DEGRADATION_Q2_BASELINE_RTOL,
                "atol": DEGRADATION_Q2_BASELINE_ATOL,
                "max_step_s": DEGRADATION_Q2_BASELINE_MAX_STEP,
                "property_model": "q1",
                "sheet": q2_sheet,
            },
            "q1_runtime_s": q1_runtime,
            "q2_runtime_s": q2_runtime,
            "q1_checks": q1_checks,
            "q2_integrator": {
                key: q2_result[key] for key in ("nfev", "njev", "nlu", "positivity_retries")
            },
            "comparison": comparison,
        }
        records.append(record)
        rows.append(
            {
                "case": f"N={n}",
                "n_intervals": n,
                "temperature_max_abs_C": comparison["temperature"]["max_abs"],
                "temperature_time_s": comparison["temperature"]["time_s"],
                "temperature_radius_cm": comparison["temperature"]["radius_cm"],
                "moisture_max_abs_kg_per_kg": comparison["moisture"]["max_abs"],
                "moisture_time_s": comparison["moisture"]["time_s"],
                "moisture_radius_cm": comparison["moisture"]["radius_cm"],
                "all_output_rounded_equal": comparison["all_output_rounded_equal"],
                "paper_rounded_equal": comparison["paper_rounded_equal"],
                "q1_runtime_s": q1_runtime,
                "q2_runtime_s": q2_runtime,
            }
        )
        np.savez_compressed(
            output / f"q2_degradation_N{n}.npz",
            time_s=times,
            radius_cm=radius_cm,
            q1_temperature_C=q1_t,
            q1_moisture_kg_per_kg=q1_c,
            q2_degraded_temperature_C=q2_result["temperature_C"],
            q2_degraded_moisture_kg_per_kg=q2_result["moisture_kg_per_kg"],
            temperature_abs_difference_C=np.abs(q1_t - q2_result["temperature_C"]),
            moisture_abs_difference_kg_per_kg=np.abs(
                q1_c - q2_result["moisture_kg_per_kg"]
            ),
            q1_balance=q1_balance,
            q2_mass_balance_residual=q2_result["mass_balance_residual"],
        )
        # Coarse grids are retained as convergence evidence.  Only the
        # finest attempted grid gates the degradation audit; if it fails,
        # append one next grid and stop as soon as that finer comparison passes.
        planned_grids = _append_next_degradation_grid(
            planned_grids, n, comparison["pass"]
        )
        grid_index += 1

    final_grid = planned_grids[-1]
    final_baseline = baseline_cache[final_grid]
    q1_t = final_baseline["q1_temperature_C"]
    q1_c = final_baseline["q1_moisture_kg_per_kg"]
    baseline = final_baseline["q2_result"]
    q1_refined_t, q1_refined_c, _, q1_refined_checks = q1.solve_bdf(
        boundary_q1, n=final_grid, rtol=DEGRADATION_Q1_REFINED_RTOL, end=end
    )
    q1_refined_checks = _annotate_q1_checks(q1_refined_checks)
    refined = q2.integrate_case(
        boundary_q2,
        n=final_grid,
        end=end,
        rtol=DEGRADATION_Q2_REFINED_RTOL,
        atol=DEGRADATION_Q2_REFINED_ATOL,
        max_step=DEGRADATION_Q2_REFINED_MAX_STEP,
        property_model="q1",
    )
    q1_time_comparison = _field_comparison(
        q1_t,
        q1_refined_t,
        q1_c,
        q1_refined_c,
        times,
        radius_cm,
        temperature_threshold=DEGRADATION_TIME_THRESHOLD_T_C,
        moisture_threshold=DEGRADATION_TIME_THRESHOLD_C,
    )
    q2_time_comparison = _field_comparison(
        baseline["temperature_C"],
        refined["temperature_C"],
        baseline["moisture_kg_per_kg"],
        refined["moisture_kg_per_kg"],
        times,
        radius_cm,
        temperature_threshold=DEGRADATION_TIME_THRESHOLD_T_C,
        moisture_threshold=DEGRADATION_TIME_THRESHOLD_C,
    )
    records.append(
        {
            "case": "time_refinement",
            "n_intervals": final_grid,
            "q1": {
                "baseline_settings": {
                    "solver": "q1.solve_bdf",
                    "rtol": DEGRADATION_Q1_BASELINE_RTOL,
                    "atol": DEGRADATION_Q1_BASELINE_RTOL * 0.01,
                    "max_step_s": DEGRADATION_Q1_INTERNAL_MAX_STEP,
                    "property_model": "q1",
                },
                "refined_settings": {
                    "solver": "q1.solve_bdf",
                    "rtol": DEGRADATION_Q1_REFINED_RTOL,
                    "atol": DEGRADATION_Q1_REFINED_RTOL * 0.01,
                    "max_step_s": DEGRADATION_Q1_INTERNAL_MAX_STEP,
                    "property_model": "q1",
                },
                "comparison": q1_time_comparison,
                "refined_checks": q1_refined_checks,
            },
            "q2": {
                "baseline_settings": {
                    "rtol": DEGRADATION_Q2_BASELINE_RTOL,
                    "atol": DEGRADATION_Q2_BASELINE_ATOL,
                    "max_step_s": DEGRADATION_Q2_BASELINE_MAX_STEP,
                    "property_model": "q1",
                },
                "refined_settings": {
                    "rtol": DEGRADATION_Q2_REFINED_RTOL,
                    "atol": DEGRADATION_Q2_REFINED_ATOL,
                    "max_step_s": DEGRADATION_Q2_REFINED_MAX_STEP,
                    "property_model": "q1",
                },
                "comparison": q2_time_comparison,
            },
        }
    )
    rows.append(
        {
            "case": "time_refinement",
            "n_intervals": final_grid,
            "q1_temperature_max_abs_C": q1_time_comparison["temperature"]["max_abs"],
            "q1_moisture_max_abs_kg_per_kg": q1_time_comparison["moisture"]["max_abs"],
            "q1_pass": q1_time_comparison["pass"],
            "q2_temperature_max_abs_C": q2_time_comparison["temperature"]["max_abs"],
            "q2_moisture_max_abs_kg_per_kg": q2_time_comparison["moisture"]["max_abs"],
            "q2_pass": q2_time_comparison["pass"],
            "all_output_rounded_equal_q1": q1_time_comparison["all_output_rounded_equal"],
            "all_output_rounded_equal_q2": q2_time_comparison["all_output_rounded_equal"],
            "q1_runtime_s": np.nan,
            "q2_runtime_s": float(refined["runtime_s"]),
        }
    )
    for record in records:
        if "n_intervals" in record and record.get("case") != "time_refinement":
            record["finest_attempted_grid"] = bool(record["n_intervals"] == final_grid)
    final_grid_record = next(
        record for record in records if record.get("n_intervals") == final_grid
    )
    _write_csv(output / "q2_degradation.csv", rows)
    degradation_pass = _degradation_gate(
        [final_grid_record], q1_time_comparison["pass"], q2_time_comparison["pass"]
    )
    payload = {
        "experiment": "q2_degradation_reproduction",
        "status": "PASS" if degradation_pass else "FAIL",
        "boundary": str(boundary_path.resolve()),
        "grids": planned_grids,
        "end_s": end,
        "records": records,
        "runtime_s": time.perf_counter() - started,
        "gate_thresholds": {
            "temperature_C": FIELD_THRESHOLD_T_C,
            "moisture_kg_per_kg": FIELD_THRESHOLD_C,
            "time_temperature_C": DEGRADATION_TIME_THRESHOLD_T_C,
            "time_moisture_kg_per_kg": DEGRADATION_TIME_THRESHOLD_C,
            "paper_times": [1800],
            "paper_radius_cm": [0.0, 2.0],
        },
        "source_snapshot": str((output / "source_snapshot" / "q2_solver.py").resolve()),
        "comparison_scope": {
            "q1_method_label": "BDF with analytic sparse Jacobian, same spatial discretization",
            "q1_method_scope": (
                "The same-spatial-discretization phrase describes only Q1's "
                "internal heat/moisture grid and does not assert cross-program "
                "boundary or surface discretization identity."
            ),
            "cross_program_scope": (
                "Q1 and Q2 retain their original spatial and boundary "
                "discretizations; their thermal surface treatments differ."
            ),
        },
        "interpretation": (
            "Q1 solve_bdf and q2.integrate_case(property_model='q1') are compared "
            "as a constitutive degradation/reproduction audit; this is not a model "
            "accuracy claim."
        ),
    }
    _write_json(output / "q2_degradation.json", payload)
    return payload


class _BalanceAccumulator:
    """Consume observer segments and integrate signed surface fluxes."""

    def __init__(
        self,
        boundary: np.ndarray,
        end: int,
        orders: tuple[int, ...] = BALANCE_CANDIDATE_ORDERS,
    ) -> None:
        self.boundary = np.asarray(boundary, dtype=float)
        self.end = int(end)
        self.orders = tuple(int(order) for order in orders)
        self.mean_moisture = np.full(self.end + 1, np.nan, dtype=float)
        self.loss_by_order = {
            order: np.full(self.end + 1, np.nan, dtype=float) for order in self.orders
        }
        self.cumulative_loss = {order: 0.0 for order in self.orders}
        self.initial_mean: float | None = None
        self.segment_records: list[dict[str, Any]] = []
        self.accepted_records: list[dict[str, np.ndarray]] = []

    @staticmethod
    def _integer_times(start: float, end: float) -> np.ndarray:
        first = int(np.ceil(start - 1e-10))
        last = int(np.floor(end + 1e-10))
        return np.arange(first, last + 1, dtype=int)

    def _mean_at_nodes(
        self,
        dense_log_solution: Any,
        node_times: np.ndarray,
        weights: np.ndarray,
    ) -> np.ndarray:
        if node_times.size == 0:
            return np.empty(0, dtype=float)
        means = np.empty(node_times.size, dtype=float)
        for chunk_start in range(0, node_times.size, 32):
            chunk = node_times[chunk_start : chunk_start + 32]
            transformed = np.asarray(dense_log_solution(chunk.astype(float)), dtype=float)
            with np.errstate(over="raise", invalid="raise", under="ignore"):
                physical_moisture = np.exp(transformed[1::2, ...])
            if not np.isfinite(physical_moisture).all() or np.any(physical_moisture <= 0.0):
                raise q2.Q2DomainError("Supplementary log(C) restoration failed without clipping.")
            if physical_moisture.ndim == 1:
                physical_moisture = physical_moisture[:, None]
            means[chunk_start : chunk_start + chunk.size] = (
                weights @ physical_moisture / TOTAL_VOLUME
            )
        return means

    def _mean_at_integer_times(
        self,
        dense_log_solution: Any,
        integer_times: np.ndarray,
        weights: np.ndarray,
    ) -> None:
        means = self._mean_at_nodes(dense_log_solution, integer_times, weights)
        if integer_times.size:
            self.mean_moisture[integer_times] = means

    def _union_nodes(self, payload: dict[str, Any]) -> np.ndarray:
        start = float(payload["start_s"])
        end = float(payload["end_s"])
        integer_times = self._integer_times(start, end).astype(float)
        union = np.concatenate(
            [
                np.asarray(payload["accepted_step_nodes"], dtype=float),
                np.asarray(payload["boundary_nodes"], dtype=float),
                integer_times,
                np.array([start, end], dtype=float),
            ]
        )
        union = np.unique(np.round(union, decimals=12))
        union = union[(union >= start - 1e-9) & (union <= end + 1e-9)]
        if union.size < 2 or np.any(np.diff(union) <= 0.0):
            raise ArithmeticError("Observer segment union nodes are not increasing.")
        return union

    def consume(self, payload: dict[str, Any]) -> None:
        start = float(payload["start_s"])
        end = float(payload["end_s"])
        if end <= start:
            raise ArithmeticError("Observer segment has non-positive duration.")
        weights = np.asarray(payload["control_volume_weights_m2"], dtype=float)
        radius_nodes = np.asarray(payload["radius_nodes_m"], dtype=float)
        if weights.shape != radius_nodes.shape:
            raise ArithmeticError("Observer radius and control-volume weights disagree.")
        if self.initial_mean is None:
            first_log = np.asarray(payload["dense_log_solution"](start), dtype=float)
            with np.errstate(over="raise", invalid="raise", under="ignore"):
                first_moisture = np.exp(first_log[1::2])
            self.initial_mean = float(weights @ first_moisture / TOTAL_VOLUME)

        integer_times = self._integer_times(start, end)
        self._mean_at_integer_times(payload["dense_log_solution"], integer_times, weights)
        nodes = self._union_nodes(payload)
        local_integrals = {order: 0.0 for order in self.orders}
        running_local = {order: 0.0 for order in self.orders}
        loss_at_nodes = {order: np.zeros(nodes.size, dtype=float) for order in self.orders}
        gauss_cache = {order: np.polynomial.legendre.leggauss(order) for order in self.orders}
        surface_log_solution = payload["moisture_surface_log_solution"]
        hm = float(payload["hm_m_s"])
        start_index = int(round(start))
        if abs(start - start_index) <= 1e-9 and 0 <= start_index <= self.end:
            for order in self.orders:
                self.loss_by_order[order][start_index] = self.cumulative_loss[order]
                loss_at_nodes[order][0] = self.cumulative_loss[order]
        for interval_index, (left, right) in enumerate(zip(nodes[:-1], nodes[1:])):
            half = 0.5 * (right - left)
            center = 0.5 * (right + left)
            for order in self.orders:
                abscissas, gauss_weights = gauss_cache[order]
                quadrature_times = center + half * abscissas
                with np.errstate(over="raise", invalid="raise", under="ignore"):
                    surface_c = np.exp(surface_log_solution(quadrature_times))
                if not np.isfinite(surface_c).all() or np.any(surface_c <= 0.0):
                    raise q2.Q2DomainError(
                        "Supplementary surface exp(log(C)) failed without clipping."
                    )
                environment_c = np.interp(quadrature_times, self.boundary[:, 0], self.boundary[:, 2])
                flux = q2.RADIUS_M * hm * (surface_c - environment_c)
                contribution = float(half * np.dot(gauss_weights, flux))
                local_integrals[order] += contribution
                running_local[order] += contribution
                right_index = int(round(right))
                if abs(right - right_index) <= 1e-9 and 0 <= right_index <= self.end:
                    self.loss_by_order[order][right_index] = (
                        self.cumulative_loss[order] + running_local[order]
                    )
                loss_at_nodes[order][interval_index + 1] = (
                    self.cumulative_loss[order] + running_local[order]
                )

        accepted_nodes = np.unique(
            np.round(np.asarray(payload["accepted_step_nodes"], dtype=float), decimals=12)
        )
        accepted_indices = np.searchsorted(nodes, accepted_nodes)
        if np.any(accepted_indices >= nodes.size) or not np.allclose(
            nodes[accepted_indices], accepted_nodes, atol=1e-10, rtol=0.0
        ):
            raise ArithmeticError("Accepted observer nodes were lost from the union.")
        accepted_means = self._mean_at_nodes(
            payload["dense_log_solution"], accepted_nodes, weights
        )
        accepted_record: dict[str, np.ndarray] = {
            "time_s": accepted_nodes,
            "mean_moisture": accepted_means,
        }
        for order in self.orders:
            accepted_loss = loss_at_nodes[order][accepted_indices]
            accepted_record[f"cumulative_loss_gl{order}"] = accepted_loss
            accepted_record[f"B_gl{order}_signed"] = (
                accepted_means + accepted_loss / TOTAL_VOLUME - self.initial_mean
            )
            accepted_record[f"B_gl{order}_abs"] = np.abs(
                accepted_record[f"B_gl{order}_signed"]
            )
        self.accepted_records.append(accepted_record)

        for order in self.orders:
            self.cumulative_loss[order] += local_integrals[order]
            final_index = int(round(end))
            if 0 <= final_index <= self.end:
                self.loss_by_order[order][final_index] = self.cumulative_loss[order]
        self.segment_records.append(
            {
                "start_s": start,
                "end_s": end,
                "accepted_step_count": int(np.asarray(payload["accepted_step_nodes"]).size),
                "boundary_node_count": int(np.asarray(payload["boundary_nodes"]).size),
                "union_node_count": int(nodes.size),
                "union_interval_count": int(nodes.size - 1),
                "local_integral_by_order": local_integrals,
                "accepted_step_trapezoid_weight_sum_s": float(
                    np.sum(payload["accepted_step_trapezoid_weights"])
                ),
            }
        )

    def finalize(self) -> dict[str, np.ndarray]:
        if self.initial_mean is None or not np.isfinite(self.mean_moisture).all():
            raise ArithmeticError("Supplementary balance did not populate all integer means.")
        for order in self.orders:
            values = self.loss_by_order[order]
            if not np.isfinite(values).all():
                raise ArithmeticError(f"Gauss-Legendre order {order} did not populate all losses.")
        result: dict[str, np.ndarray] = {
            "mean_moisture": self.mean_moisture.copy(),
        }
        for order in self.orders:
            result[f"cumulative_loss_gl{order}"] = self.loss_by_order[order].copy()
            result[f"B_gl{order}_signed"] = (
                self.mean_moisture + self.loss_by_order[order] / TOTAL_VOLUME - self.initial_mean
            )
            result[f"B_gl{order}_abs"] = np.abs(result[f"B_gl{order}_signed"])
        if not self.accepted_records:
            raise ArithmeticError("No accepted observer nodes were recorded.")
        accepted_times = np.concatenate([record["time_s"] for record in self.accepted_records])
        accepted_order = np.argsort(accepted_times, kind="stable")
        sorted_times = accepted_times[accepted_order]
        keep_last = np.r_[
            np.round(sorted_times[:-1], decimals=12)
            != np.round(sorted_times[1:], decimals=12),
            True,
        ]
        result["accepted_time_s"] = sorted_times[keep_last]
        for key in ("mean_moisture",) + tuple(
            f"{prefix}{order}{suffix}"
            for order in self.orders
            for prefix, suffix in (
                ("cumulative_loss_gl", ""),
                ("B_gl", "_signed"),
                ("B_gl", "_abs"),
            )
        ):
            values = np.concatenate([record[key] for record in self.accepted_records])
            result[f"accepted_{key}"] = values[accepted_order][keep_last]
        return result


def _integrate_with_balance_observer(
    boundary: np.ndarray,
    n: int,
    end: int,
    rtol: float,
    atol: float,
    max_step: float,
) -> tuple[dict[str, Any], _BalanceAccumulator, dict[str, np.ndarray], q2.SegmentObserver]:
    accumulator = _BalanceAccumulator(boundary, end, BALANCE_CANDIDATE_ORDERS)
    observer = q2.SegmentObserver(retain=False, callback=accumulator.consume)
    result = q2.integrate_case(
        boundary,
        n=n,
        end=end,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        property_model="q2",
        segment_observer=observer,
    )
    return result, accumulator, accumulator.finalize(), observer


def _select_quadrature_orders(balance: dict[str, np.ndarray]) -> dict[str, Any]:
    comparisons: list[dict[str, Any]] = []
    selected: list[int] | None = None
    for left, right in zip(BALANCE_CANDIDATE_ORDERS[:-1], BALANCE_CANDIDATE_ORDERS[1:]):
        integer_difference = float(
            np.max(np.abs(balance[f"B_gl{right}_signed"] - balance[f"B_gl{left}_signed"]))
        )
        accepted_difference = float(
            np.max(
                np.abs(
                    balance[f"accepted_B_gl{right}_signed"]
                    - balance[f"accepted_B_gl{left}_signed"]
                )
            )
        )
        difference = max(integer_difference, accepted_difference)
        item = {
            "lower_order": left,
            "higher_order": right,
            "integer_max_abs_signed_B_difference": integer_difference,
            "accepted_max_abs_signed_B_difference": accepted_difference,
            "joint_max_abs_signed_B_difference": difference,
            "max_abs_signed_B_difference": difference,
            "threshold": BALANCE_QUADRATURE_THRESHOLD,
            "pass": bool(difference <= BALANCE_QUADRATURE_THRESHOLD),
        }
        comparisons.append(item)
        if selected is None and item["pass"]:
            selected = [left, right]
    return {
        "comparisons": comparisons,
        "selected_orders": selected,
        "pass": selected is not None,
    }


def _accepted_summary(balance: dict[str, np.ndarray]) -> dict[str, Any]:
    times = balance["accepted_time_s"]
    summary: dict[str, Any] = {"accepted_node_count": int(times.size)}
    for order in BALANCE_CANDIDATE_ORDERS:
        values = balance[f"accepted_B_gl{order}_abs"]
        location = int(np.argmax(values))
        summary[f"gl{order}"] = {
            "max_abs_B": float(values[location]),
            "time_s": float(times[location]),
            "mean_moisture": float(balance["accepted_mean_moisture"][location]),
        }
    return summary


def _selected_balance_summary(
    balance: dict[str, np.ndarray], selection: dict[str, Any]
) -> dict[str, Any]:
    """Summarize selected-order integer and accepted-node residuals."""
    selected_orders = selection.get("selected_orders")
    selected_order = (
        int(selected_orders[1])
        if selected_orders is not None
        else int(BALANCE_CANDIDATE_ORDERS[-1])
    )
    integer_abs = np.asarray(balance[f"B_gl{selected_order}_abs"], dtype=float)
    accepted_abs = np.asarray(balance[f"accepted_B_gl{selected_order}_abs"], dtype=float)
    integer_index = int(np.argmax(integer_abs))
    accepted_index = int(np.argmax(accepted_abs))
    floor = None
    if selected_orders is not None:
        floor = next(
            item["max_abs_signed_B_difference"]
            for item in selection["comparisons"]
            if item["lower_order"] == selected_orders[0]
            and item["higher_order"] == selected_orders[1]
        )
    return {
        "selected_order": selected_order,
        "quadrature_floor_max_abs_B": floor,
        "integer_seconds": {
            "max_abs_B": float(integer_abs[integer_index]),
            "time_s": integer_index,
            "threshold": BALANCE_THRESHOLD_C,
            "pass": bool(integer_abs[integer_index] <= BALANCE_THRESHOLD_C),
        },
        "accepted_nodes": {
            "max_abs_B": float(accepted_abs[accepted_index]),
            "time_s": float(balance["accepted_time_s"][accepted_index]),
            "threshold": BALANCE_THRESHOLD_C,
            "pass": bool(accepted_abs[accepted_index] <= BALANCE_THRESHOLD_C),
        },
    }


def _checkpoint_summary(
    balance: dict[str, np.ndarray],
    trapezoid_signed: np.ndarray,
    selected_order: int,
    times: np.ndarray,
) -> list[dict[str, Any]]:
    rows = []
    for hour in (0.5, 1.0, 2.0, 3.0):
        second = int(round(hour * 3600.0))
        index = int(np.where(times == second)[0][0])
        rows.append(
            {
                "time_h": hour,
                "time_s": second,
                "mean_moisture": float(balance["mean_moisture"][index]),
                "B_gl_selected_signed": float(
                    balance[f"B_gl{selected_order}_signed"][index]
                ),
                "B_gl_selected_abs": float(balance[f"B_gl{selected_order}_abs"][index]),
                "B_trapezoid_signed": float(trapezoid_signed[index]),
            }
        )
    return rows


def run_balance(
    boundary_path: Path,
    output: Path,
    n: int = BALANCE_N,
    end: int = BALANCE_END_S,
) -> dict[str, Any]:
    """Run formal Q2 settings and independently integrate the surface flux."""
    boundary, boundary_sheet = q2.read_boundary(boundary_path, end=end)
    started = time.perf_counter()
    baseline, baseline_accumulator, balance, observer = _integrate_with_balance_observer(
        boundary,
        n,
        end,
        BALANCE_Q2_BASELINE_RTOL,
        BALANCE_Q2_BASELINE_ATOL,
        BALANCE_Q2_BASELINE_MAX_STEP,
    )
    times = baseline["times_s"].astype(int)
    outward_rate = q2.RADIUS_M * baseline["hm_m_s"] * (
        baseline["moisture_kg_per_kg"][:, -1] - baseline["environment_moisture_kg_per_kg"]
    )
    trapezoid_loss = np.zeros(end + 1, dtype=float)
    trapezoid_loss[1:] = np.cumsum(0.5 * (outward_rate[:-1] + outward_rate[1:]))
    balance["cumulative_loss_trapezoid"] = trapezoid_loss
    balance["B_trapezoid_signed"] = (
        baseline["mean_moisture"]
        + trapezoid_loss / TOTAL_VOLUME
        - baseline["initial_mean_moisture"]
    )
    balance["B_trapezoid_abs"] = np.abs(balance["B_trapezoid_signed"])
    balance["q2_mass_balance_residual"] = baseline["mass_balance_residual"]
    independent_mean_difference = balance["mean_moisture"] - baseline["mean_moisture"]
    balance["independent_mean_difference"] = independent_mean_difference

    quadrature_selection = _select_quadrature_orders(balance)
    selected_order = (
        quadrature_selection["selected_orders"][1]
        if quadrature_selection["selected_orders"] is not None
        else BALANCE_CANDIDATE_ORDERS[-1]
    )
    mean_crosscheck = {
        "max_abs": float(np.max(np.abs(independent_mean_difference))),
        "time_s": int(times[int(np.argmax(np.abs(independent_mean_difference)))]),
        "threshold": MEAN_CROSSCHECK_THRESHOLD,
        "pass": bool(np.max(np.abs(independent_mean_difference)) <= MEAN_CROSSCHECK_THRESHOLD),
    }

    # The time-refined state and its independent GL balance are both retained;
    # this keeps the time gate on B itself rather than only on sampled fields.
    refined, refined_accumulator, refined_balance, refined_observer = _integrate_with_balance_observer(
        boundary,
        n,
        end,
        BALANCE_Q2_REFINED_RTOL,
        BALANCE_Q2_REFINED_ATOL,
        BALANCE_Q2_REFINED_MAX_STEP,
    )
    refined_selection = _select_quadrature_orders(refined_balance)
    refined_order = (
        refined_selection["selected_orders"][1]
        if refined_selection["selected_orders"] is not None
        else BALANCE_CANDIDATE_ORDERS[-1]
    )
    baseline_selected_summary = _selected_balance_summary(balance, quadrature_selection)
    refined_selected_summary = _selected_balance_summary(refined_balance, refined_selection)
    field_time_difference = _field_comparison(
        baseline["temperature_C"],
        refined["temperature_C"],
        baseline["moisture_kg_per_kg"],
        refined["moisture_kg_per_kg"],
        times,
        baseline["radius_output_cm"],
    )
    time_balance_gate = _balance_time_gate(
        baseline_selected_summary, refined_selected_summary
    )
    balance_time_difference = {
        "selected_order": selected_order,
        "refined_selected_order": refined_order,
        "criterion": time_balance_gate,
        "refined_quadrature_selection": refined_selection,
    }

    np.savez_compressed(output / "q2_balance.npz", time_s=times, **balance)
    np.savez_compressed(
        output / "q2_balance_time_refined.npz",
        time_s=times,
        **{f"refined_{key}": value for key, value in refined_balance.items()},
        baseline_temperature_C=baseline["temperature_C"],
        baseline_moisture_kg_per_kg=baseline["moisture_kg_per_kg"],
        refined_temperature_C=refined["temperature_C"],
        refined_moisture_kg_per_kg=refined["moisture_kg_per_kg"],
    )

    csv_rows: list[dict[str, Any]] = []
    for index, second in enumerate(times):
        row: dict[str, Any] = {
            "time_s": int(second),
            "mean_moisture": float(balance["mean_moisture"][index]),
            "cumulative_loss_trapezoid": float(trapezoid_loss[index]),
            "B_trapezoid_signed": float(balance["B_trapezoid_signed"][index]),
            "B_trapezoid_abs": float(balance["B_trapezoid_abs"][index]),
            "q2_mass_balance_residual": float(balance["q2_mass_balance_residual"][index]),
        }
        for order in BALANCE_CANDIDATE_ORDERS:
            row[f"cumulative_loss_gl{order}"] = float(balance[f"cumulative_loss_gl{order}"][index])
            row[f"B_gl{order}_signed"] = float(balance[f"B_gl{order}_signed"][index])
            row[f"B_gl{order}_abs"] = float(balance[f"B_gl{order}_abs"][index])
        csv_rows.append(row)
    _write_csv(output / "q2_balance.csv", csv_rows)

    accepted_rows: list[dict[str, Any]] = []
    accepted_count = balance["accepted_time_s"].size
    for index in range(accepted_count):
        row = {
            "time_s": float(balance["accepted_time_s"][index]),
            "mean_moisture": float(balance["accepted_mean_moisture"][index]),
        }
        for order in BALANCE_CANDIDATE_ORDERS:
            row[f"cumulative_loss_gl{order}"] = float(
                balance[f"accepted_cumulative_loss_gl{order}"][index]
            )
            row[f"B_gl{order}_signed"] = float(balance[f"accepted_B_gl{order}_signed"][index])
            row[f"B_gl{order}_abs"] = float(balance[f"accepted_B_gl{order}_abs"][index])
        accepted_rows.append(row)
    _write_csv(output / "q2_balance_accepted_nodes.csv", accepted_rows)

    gl_vs_trapezoid = {
        f"gl{order}_vs_trapezoid_max_abs_B": float(
            np.max(np.abs(balance[f"B_gl{order}_signed"] - balance["B_trapezoid_signed"]))
        )
        for order in BALANCE_CANDIDATE_ORDERS
    }
    balance_gate = _balance_gate(
        quadrature_selection,
        refined_selection,
        baseline_selected_summary,
        refined_selected_summary,
        mean_crosscheck["pass"],
        field_time_difference["pass"],
        time_balance_gate["pass"],
    )
    checkpoint_rows = _checkpoint_summary(
        balance, balance["B_trapezoid_signed"], selected_order, times
    )
    summary = {
        "experiment": "q2_integral_balance",
        "status": "PASS" if balance_gate else "FAIL",
        "boundary": str(boundary_path.resolve()),
        "boundary_sheet": boundary_sheet,
        "n_intervals": n,
        "end_s": end,
        "settings": {
            "formal": {
                "rtol": BALANCE_Q2_BASELINE_RTOL,
                "atol": BALANCE_Q2_BASELINE_ATOL,
                "max_step_s": BALANCE_Q2_BASELINE_MAX_STEP,
                "property_model": "q2",
            },
            "time_refined": {
                "rtol": BALANCE_Q2_REFINED_RTOL,
                "atol": BALANCE_Q2_REFINED_ATOL,
                "max_step_s": BALANCE_Q2_REFINED_MAX_STEP,
                "property_model": "q2",
            },
            "gauss_legendre_candidate_orders": list(BALANCE_CANDIDATE_ORDERS),
            "gauss_legendre_selected_orders": quadrature_selection["selected_orders"],
            "node_union": ["accepted_step_nodes", "original_boundary_nodes", "integer_seconds"],
        },
        "observer_segments": observer.total_segments,
        "refined_observer_segments": refined_observer.total_segments,
        "segment_records": baseline_accumulator.segment_records,
        "refined_segment_records": refined_accumulator.segment_records,
        "quadrature_selection": quadrature_selection,
        "refined_quadrature_selection": refined_selection,
        "balance_summary": {
            "initial_mean_moisture": baseline_accumulator.initial_mean,
            "max_abs_B_integer_seconds": {
                f"gl{order}": float(np.max(balance[f"B_gl{order}_abs"]))
                for order in BALANCE_CANDIDATE_ORDERS
            },
            "accepted_nodes": _accepted_summary(balance),
            "selected_order": selected_order,
            "selected_order_residuals": baseline_selected_summary,
            "refined_selected_order_residuals": refined_selected_summary,
            "refined_accepted_nodes": _accepted_summary(refined_balance),
            "max_abs_B_trapezoid": float(np.max(balance["B_trapezoid_abs"])),
            "max_abs_q2_original_residual": float(
                np.max(np.abs(balance["q2_mass_balance_residual"]))
            ),
            "gl_vs_original_trapezoid": gl_vs_trapezoid,
            "mean_crosscheck": mean_crosscheck,
            "checkpoints": checkpoint_rows,
        },
        "time_refinement": {
            "field_comparison": field_time_difference,
            "balance_comparison": balance_time_difference,
            "criterion": time_balance_gate,
        },
        "runtime_s": time.perf_counter() - started,
        "source_snapshot": str((output / "source_snapshot" / "q2_solver.py").resolve()),
        "interpretation": (
            "B is a signed moisture balance residual: mean(C) plus integrated "
            "outward surface loss divided by the cylindrical volume minus the "
            "initial mean. Gauss--Legendre uses the accepted-segment dense log "
            "trajectory and independently restores C=exp(z); no clipping is used."
        ),
    }
    _write_json(output / "q2_balance.json", summary)
    _write_json(
        output / "q2_balance_segments.json",
        {
            "segments": baseline_accumulator.segment_records,
            "refined_segments": refined_accumulator.segment_records,
        },
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("degradation", "balance", "all"))
    parser.add_argument("--boundary", type=Path, default=DEFAULT_BOUNDARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    output = args.output.resolve()
    boundary = args.boundary.resolve()
    output.mkdir(parents=True, exist_ok=True)
    _write_json(output / "environment.json", _environment_payload(boundary, output))
    (output / "reproduction_command.txt").write_text(
        "\n".join(_stable_reproduction_commands()) + "\n", encoding="utf-8"
    )

    results: dict[str, Any] = {}
    if args.mode in ("degradation", "all"):
        results["degradation"] = run_degradation(boundary, output)
    if args.mode in ("balance", "all"):
        results["balance"] = run_balance(boundary, output)
    refresh_combined_summary(output)
    print(json.dumps(_json_safe(results), ensure_ascii=False))


if __name__ == "__main__":
    main()
