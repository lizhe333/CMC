"""Named development contracts for the Question 4 moving-radius solver.

The script is intentionally an audit harness rather than a generic smoke
test.  Each check below targets a specified modelling, discretisation,
regression, or output risk.  It writes only development evidence under
``results/q4/dev``; it never touches the official workbook or paper files.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp

CODE_DIR = Path(__file__).resolve().parent
PROJECT = CODE_DIR.parent
DEV = PROJECT / "results" / "q4" / "dev"
sys.path.insert(0, str(CODE_DIR))

import q4_solver as q4
import q3_solver as q3
import q2_solver as q2


RESULTS: dict[str, Any] = {
    "script": str(Path(__file__).resolve()),
    "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    "contracts": {},
}
SOURCE_PATH = Path(q4.__file__).resolve()
RESULTS["run_id"] = f"q4_dev_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:10]}"
RESULTS["source_path"] = str(SOURCE_PATH)
RESULTS["source_mtime_utc"] = datetime.fromtimestamp(SOURCE_PATH.stat().st_mtime, timezone.utc).isoformat()
RESULTS["test_source_mtime_utc"] = datetime.fromtimestamp(Path(__file__).resolve().stat().st_mtime, timezone.utc).isoformat()
FAILURES: list[str] = []
ARTIFACTS: dict[str, np.ndarray] = {}


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


def check(name: str, passed: bool, **metrics: Any) -> None:
    record = {"pass": bool(passed), **_safe(metrics)}
    RESULTS["contracts"][name] = record
    if not passed:
        FAILURES.append(name)
    status = "PASS" if passed else "FAIL"
    print(f"[{status}] {name}")
    if metrics:
        print(json.dumps(_safe(metrics), ensure_ascii=False, sort_keys=True))


def _assert_no_exception(name: str, fn: Callable[[], dict[str, Any]]) -> None:
    try:
        metrics = fn()
    except Exception as exc:  # noqa: BLE001 - preserve failure evidence
        check(name, False, exception=repr(exc))
        return
    check(name, True, **metrics)


def _constant_radius(radius_m: float = q4.RADIUS_INITIAL_M) -> q4.RadiusModel:
    return q4.RadiusModel(
        np.array([0.0, 14.0 * 24.0 * 3600.0]),
        np.array([radius_m, radius_m]),
        interpolation="linear",
        post_mode="plateau",
    )


def _state_fields(t: float, xi: np.ndarray) -> tuple[np.ndarray, ...]:
    """Smooth manufactured fields and their derivatives in (xi,t)."""
    xi = np.asarray(xi, dtype=float)
    e_t = np.exp(-float(t) / 7200.0)
    e_c = np.exp(-float(t) / 5400.0)
    e_s = np.exp(-float(t) / 3600.0)
    one_s = 1.0 - e_s
    T = 32.0 + 2.0 * e_t * (1.0 - xi**2) + 0.5 * one_s * xi**4
    C = 0.40 + 0.35 * e_c * (1.0 - xi**2) + 0.08 * one_s * xi**4
    T_t = -2.0 * e_t / 7200.0 * (1.0 - xi**2) + 0.5 * e_s / 3600.0 * xi**4
    C_t = -0.35 * e_c / 5400.0 * (1.0 - xi**2) + 0.08 * e_s / 3600.0 * xi**4
    T_x = -4.0 * e_t * xi + 2.0 * one_s * xi**3
    C_x = -0.70 * e_c * xi + 0.32 * one_s * xi**3
    T_xx = -4.0 * e_t + 6.0 * one_s * xi**2
    C_xx = -0.70 * e_c + 0.96 * one_s * xi**2
    return T, C, T_t, C_t, T_x, C_x, T_xx, C_xx


def _manufactured_environment(radius: Callable[[float], float], h: float, hm: float) -> Callable[[float], tuple[float, float]]:
    def evaluate(t: float) -> tuple[float, float]:
        xi = np.array([1.0])
        T, C, _Tt, _Ct, Tx, Cx, _Txx, _Cxx = _state_fields(float(t), xi)
        rho, cp, k, D = q4.material_properties(C, T, model="appendix4")
        R = float(radius(float(t)))
        return (
            float(T[0] + k[0] * Tx[0] / (h * R)),
            float(C[0] + D[0] * Cx[0] / (hm * R)),
        )

    return evaluate


def _manufactured_source(radius: Callable[[float], float], h: float, hm: float) -> Callable[..., tuple[np.ndarray, np.ndarray]]:
    def source(t: float, radius_m: float, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        T, C, T_t, C_t, T_x, C_x, T_xx, C_xx = _state_fields(float(t), nodes)
        rho, cp, k, D = q4.material_properties(C, T, model="appendix4")
        rho_C, cp_C, k_C, D_C, D_T = q4.property_derivatives(C, T, model="appendix4")
        heat_operator = np.empty_like(nodes)
        moisture_operator = np.empty_like(nodes)
        positive = nodes > 0.0
        heat_operator[positive] = (
            k[positive] * T_x[positive] / nodes[positive]
            + k_C[positive] * C_x[positive] * T_x[positive]
            + k[positive] * T_xx[positive]
        )
        moisture_operator[positive] = (
            D[positive] * C_x[positive] / nodes[positive]
            + (D_C[positive] * C_x[positive] + D_T[positive] * T_x[positive]) * C_x[positive]
            + D[positive] * C_xx[positive]
        )
        heat_operator[~positive] = 2.0 * k[~positive] * T_xx[~positive]
        moisture_operator[~positive] = 2.0 * D[~positive] * C_xx[~positive]
        return (
            rho * cp * T_t - heat_operator / (float(radius_m) ** 2),
            C_t - moisture_operator / (float(radius_m) ** 2),
        )

    return source


def _manufactured_solution(n: int, radius: Callable[[float], float], *, omit_r2: bool = False, omit_r1: bool = False, end_s: float = 14400.0, max_step: float = 5.0) -> tuple[Any, np.ndarray, np.ndarray]:
    _dxi, nodes, _faces, _weights = q4.radial_geometry(n)
    h, hm = q4.H_BASE, q4.HM_BASE
    env = _manufactured_environment(radius, h, hm)
    source = _manufactured_source(radius, h, hm)
    T0, C0, *_ = _state_fields(0.0, nodes)
    initial = q4.interleaved_state(T0, C0)
    transformed = q4._to_log_state(initial, n, model="appendix4")
    rhs, jac, geometry = q4.make_log_rhs(
        env,
        radius,
        n,
        model="appendix4",
        h=h,
        hm=hm,
        source=source,
        omit_r2=omit_r2,
        omit_r1=omit_r1,
    )
    solution = solve_ivp(
        rhs,
        (0.0, float(end_s)),
        transformed,
        method="BDF",
        jac=None,
        jac_sparsity=q4.jacobian_sparsity(n),
        rtol=2.0e-9,
        atol=2.0e-11,
        max_step=float(max_step),
        dense_output=True,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution, nodes, geometry[3]


def contract_input_data(boundary: np.ndarray, radius: q4.RadiusModel, radius_metadata: dict[str, Any]) -> None:
    _assert_no_exception(
        "input_radius_and_environment",
        lambda: {
            "radius_record_count": radius_metadata["record_count"],
            "radius_start_m": radius_metadata["start_radius_m"],
            "radius_end_m": radius_metadata["end_radius_m"],
            "radius_monotone": bool(np.all(np.diff(radius(radius.times_s)) <= 1.0e-12)),
            "derivative_finite": bool(np.isfinite(radius.derivative(np.linspace(0.0, radius.end_s, 401))).all()),
            "boundary_rows": int(len(boundary)),
            "boundary_end_s": float(boundary[-1, 0]),
        },
    )
    check(
        "input_radius_contract",
        radius_metadata["record_count"] == 145
        and np.isclose(radius_metadata["start_radius_m"], 0.020)
        and np.isclose(radius_metadata["end_radius_m"], 0.01198)
        and radius_metadata["interpolation"] == "pchip"
        and radius_metadata["outside_evaluation"].startswith("explicit"),
        record_count=radius_metadata["record_count"],
        start_radius_m=radius_metadata["start_radius_m"],
        end_radius_m=radius_metadata["end_radius_m"],
    )
    plateau_value = float(radius(radius.end_s + 3600.0))
    check(
        "radius_no_extrapolation_and_plateau",
        np.isclose(plateau_value, radius.final_radius_m, rtol=0.0, atol=1.0e-14)
        and np.isclose(radius.derivative(radius.end_s + 3600.0), 0.0, rtol=0.0, atol=1.0e-20),
        post_radius_m=plateau_value,
    )
    try:
        radius(-1.0)
    except ValueError:
        before_guard = True
    else:
        before_guard = False
    check("radius_preinterval_guard", before_guard)
    slow = q4.RadiusModel(radius.times_s, radius.radii_m, interpolation="pchip", post_mode="slow1pct")
    check(
        "radius_slow1pct_scenario",
        slow(slow.end_s + 10.0 * 86400.0) < radius.final_radius_m
        and slow(slow.end_s + 10.0 * 86400.0) > 0.98 * radius.final_radius_m,
        slow_radius_m=float(slow(slow.end_s + 10.0 * 86400.0)),
        slow_derivative_m_s=float(slow.derivative(slow.end_s + 10.0 * 86400.0)),
    )
    check(
        "environment_measured_range_guard",
        np.allclose(q4.measured_environment(boundary)(0.0), boundary[0, 1:])
        and np.allclose(q4.measured_environment(boundary)(q4.SWITCH_S), boundary[-1, 1:]),
    )


def contract_geometry_and_purity() -> None:
    radius = lambda t: 0.02 * (1.0 - 0.1 * float(t) / 100.0)
    environment = lambda _t: (50.0, 0.05)
    n = 40
    state = q4.initial_state(n)
    rhs, geometry = q4.make_rhs(environment, radius, n, model="appendix4", h=0.0, hm=0.0)
    direct = rhs(0.0, state)
    solution = solve_ivp(rhs, (0.0, 100.0), state, method="BDF", rtol=1.0e-10, atol=1.0e-12, max_step=5.0)
    C_change = float(np.max(np.abs(solution.y[1::2, :] - q4.INITIAL_C)))
    check(
        "pure_shrink_uniform_moisture",
        solution.success and np.max(np.abs(direct[1::2])) <= 1.0e-14 and C_change <= 1.0e-11,
        solver_success=bool(solution.success),
        max_initial_rhs_C=float(np.max(np.abs(direct[1::2]))),
        max_C_change=C_change,
    )
    _dxi, nodes, _faces, weights = q4.radial_geometry(n)
    T = 40.0 + 4.0 * nodes**2 + 0.2 * np.sin(3.0 * nodes)
    C = 0.20 + 0.20 * nodes**2 + 0.02 * np.sin(2.0 * nodes)
    state = q4.interleaved_state(T, C)
    rhs, _ = q4.make_rhs(environment, radius, n, model="appendix4", h=0.0, hm=0.0)
    derivative = rhs(17.0, state)
    rho, cp, _k, _D = q4.material_properties(C, T, model="appendix4")
    moisture_sum = float(np.dot(weights, derivative[1::2]))
    heat_sum = float(np.dot(weights, rho * cp * derivative[0::2]))
    check(
        "internal_flux_telescoping",
        abs(moisture_sum) <= 1.0e-12 and abs(heat_sum) <= 1.0e-9,
        moisture_weighted_sum=moisture_sum,
        heat_capacity_weighted_sum=heat_sum,
    )


def contract_jacobian(boundary: np.ndarray, radius: q4.RadiusModel) -> None:
    measured = q4.measured_environment(boundary)
    plateau = q4.plateau_environment()
    n = 20
    _dxi, nodes, _faces, _weights = q4.radial_geometry(n)
    errors: list[float] = []
    block_errors: list[dict[str, float]] = []
    cross_blocks: list[float] = []
    for t in (7200.0, 10800.0, 14400.0, 18000.0):
        env = measured if t <= q4.SWITCH_S else plateau
        state = q4.interleaved_state(
            42.0 + 3.0 * nodes**2 + 0.2 * np.sin(4.0 * nodes + t / 5000.0),
            0.18 + 0.25 * nodes**2 + 0.025 * np.cos(3.0 * nodes + t / 7000.0),
        )
        transformed = q4._to_log_state(state, n, model="appendix4")
        rhs, jac, _ = q4.make_log_rhs(env, radius, n, model="appendix4")
        analytic = jac(t, transformed).toarray()
        numerical = np.empty_like(analytic)
        for column in range(len(transformed)):
            step = 1.0e-6 * max(1.0, abs(float(transformed[column])))
            plus = transformed.copy()
            minus = transformed.copy()
            plus[column] += step
            minus[column] -= step
            numerical[:, column] = (rhs(t, plus) - rhs(t, minus)) / (2.0 * step)
        denominator = 1.0 + np.abs(numerical)
        relative = float(np.max(np.abs(analytic - numerical) / denominator))
        errors.append(relative)
        block_errors.append({
            "TT": float(np.max(np.abs(analytic[0::2, 0::2] - numerical[0::2, 0::2]) / (1.0 + np.abs(numerical[0::2, 0::2])))),
            "TC": float(np.max(np.abs(analytic[0::2, 1::2] - numerical[0::2, 1::2]) / (1.0 + np.abs(numerical[0::2, 1::2])))),
            "CT": float(np.max(np.abs(analytic[1::2, 0::2] - numerical[1::2, 0::2]) / (1.0 + np.abs(numerical[1::2, 0::2])))),
            "CC": float(np.max(np.abs(analytic[1::2, 1::2] - numerical[1::2, 1::2]) / (1.0 + np.abs(numerical[1::2, 1::2])))),
        })
        cross_blocks.append(float(max(np.max(np.abs(analytic[0::2, 1::2])), np.max(np.abs(analytic[1::2, 0::2])))))
    check(
        "analytic_log_jacobian_central_difference",
        max(errors) <= 3.0e-5 and min(cross_blocks) > 0.0,
        times_s=[7200.0, 10800.0, 14400.0, 18000.0],
        max_relative_error=max(errors),
        block_relative_error=block_errors,
        cross_block_max_abs=cross_blocks,
    )


def contract_manufactured_solution(radius: q4.RadiusModel, *, formal: bool = False) -> None:
    """Run the Q4_MODELING 9.1 contract or a separately labelled precheck."""
    formal_grids = [640, 1280, 2560]
    formal_times = np.array([3600.0, 7200.0, 10800.0, 14400.0])
    quick_grids = [40, 80, 160]
    quick_times = np.array([300.0, 600.0, 900.0, 1200.0])
    ns = formal_grids if formal else quick_grids
    sample_times = formal_times if formal else quick_times
    duration = 14400.0 if formal else 1200.0
    errors_T: list[float] = []
    errors_C: list[float] = []
    errors_2_T: list[float] = []
    errors_2_C: list[float] = []
    wrong_r2_T: list[float] = []
    wrong_r2_C: list[float] = []
    wrong_r1_T: list[float] = []
    wrong_r1_C: list[float] = []
    for n in ns:
        solution, nodes, weights = _manufactured_solution(n, radius, end_s=duration)
        physical = q4._physical_from_transformed_matrix(solution.sol(sample_times), n)
        exact_T = np.stack([_state_fields(t, nodes)[0] for t in sample_times], axis=0)
        exact_C = np.stack([_state_fields(t, nodes)[1] for t in sample_times], axis=0)
        error_T = np.abs(physical[0::2, :].T - exact_T)
        error_C = np.abs(physical[1::2, :].T - exact_C)
        errors_T.append(float(np.max(error_T)))
        errors_C.append(float(np.max(error_C)))
        errors_2_T.append(float(np.max(np.sqrt(2.0 * np.sum(weights[None, :] * error_T**2, axis=1)))))
        errors_2_C.append(float(np.max(np.sqrt(2.0 * np.sum(weights[None, :] * error_C**2, axis=1)))))
        for fault_name, fault_kwargs, target_T, target_C in (
            ("omit_R_squared_scale", {"omit_r2": True}, wrong_r2_T, wrong_r2_C),
            ("omit_R_boundary_scale", {"omit_r1": True}, wrong_r1_T, wrong_r1_C),
        ):
            faulty, faulty_nodes, faulty_weights = _manufactured_solution(n, radius, end_s=duration, **fault_kwargs)
            faulty_physical = q4._physical_from_transformed_matrix(faulty.sol(sample_times), n)
            fault_exact_T = np.stack([_state_fields(t, faulty_nodes)[0] for t in sample_times], axis=0)
            fault_exact_C = np.stack([_state_fields(t, faulty_nodes)[1] for t in sample_times], axis=0)
            fault_error_T = np.abs(faulty_physical[0::2, :].T - fault_exact_T)
            fault_error_C = np.abs(faulty_physical[1::2, :].T - fault_exact_C)
            target_T.append(float(np.max(fault_error_T)))
            target_C.append(float(np.max(fault_error_C)))
    order_T = [math.log(errors_T[i] / errors_T[i + 1], 2.0) for i in range(len(ns) - 1)]
    order_C = [math.log(errors_C[i] / errors_C[i + 1], 2.0) for i in range(len(ns) - 1)]
    order_2_T = [math.log(errors_2_T[i] / errors_2_T[i + 1], 2.0) for i in range(len(ns) - 1)]
    order_2_C = [math.log(errors_2_C[i] / errors_2_C[i + 1], 2.0) for i in range(len(ns) - 1)]
    label = "manufactured_solution_formal_gate" if formal else "manufactured_solution_quick_precheck"
    check(
        label,
        min(order_T) >= (1.5 if formal else 1.25) and min(order_C) >= (1.5 if formal else 1.25)
        and min(order_2_T) >= (1.5 if formal else 1.25) and min(order_2_C) >= (1.5 if formal else 1.25),
        contract_kind="formal_gate" if formal else "quick_precheck",
        grids=ns,
        sample_times_s=sample_times,
        duration_s=duration,
        max_error_T_C=errors_T,
        max_error_C=errors_C,
        max_error_2_T_C=errors_2_T,
        max_error_2_C=errors_2_C,
        observed_order_T=order_T,
        observed_order_C=order_C,
        observed_order_2_T=order_2_T,
        observed_order_2_C=order_2_C,
        formal_specification={"grids": formal_grids, "sample_times_s": formal_times.tolist(), "expected_order": 2.0, "last_two_order_threshold": 1.5},
    )
    def injection_detected(wrong_T: list[float], wrong_C: list[float]) -> bool:
        wrong_order_T = [math.log(wrong_T[i] / wrong_T[i + 1], 2.0) for i in range(len(wrong_T) - 1)]
        wrong_order_C = [math.log(wrong_C[i] / wrong_C[i + 1], 2.0) for i in range(len(wrong_C) - 1)]
        return (
            (wrong_order_T and min(wrong_order_T) < 1.0)
            or (wrong_order_C and min(wrong_order_C) < 1.0)
            or wrong_T[-1] >= 5.0 * errors_T[-1]
            or wrong_C[-1] >= 5.0 * errors_C[-1]
        )
    check(
        label + "_omit_R_squared_injection",
        injection_detected(wrong_r2_T, wrong_r2_C),
        fault="omit_R_squared_scale",
        wrong_error_T_C=wrong_r2_T,
        wrong_error_C=wrong_r2_C,
        correct_fine_error_T_C=errors_T[-1],
        correct_fine_error_C=errors_C[-1],
    )
    check(
        label + "_omit_R_boundary_injection",
        injection_detected(wrong_r1_T, wrong_r1_C),
        fault="omit_R_boundary_scale",
        wrong_error_T_C=wrong_r1_T,
        wrong_error_C=wrong_r1_C,
        correct_fine_error_T_C=errors_T[-1],
        correct_fine_error_C=errors_C[-1],
    )
    if not formal:
        RESULTS.setdefault("formal_contracts_not_run", []).append({"name": "manufactured_solution", "required_grids": formal_grids, "required_times_s": formal_times.tolist(), "reason": "development run intentionally limited to quick precheck"})
    ARTIFACTS["manufactured_grid_N"] = np.asarray(ns, dtype=int)
    ARTIFACTS["manufactured_error_T_C"] = np.asarray(errors_T, dtype=float)
    ARTIFACTS["manufactured_error_C"] = np.asarray(errors_C, dtype=float)
    ARTIFACTS["manufactured_error_2_T_C"] = np.asarray(errors_2_T, dtype=float)
    ARTIFACTS["manufactured_error_2_C"] = np.asarray(errors_2_C, dtype=float)
    ARTIFACTS["manufactured_wrong_r2_error_T_C"] = np.asarray(wrong_r2_T, dtype=float)
    ARTIFACTS["manufactured_wrong_r2_error_C"] = np.asarray(wrong_r2_C, dtype=float)
    ARTIFACTS["manufactured_wrong_r1_error_T_C"] = np.asarray(wrong_r1_T, dtype=float)
    ARTIFACTS["manufactured_wrong_r1_error_C"] = np.asarray(wrong_r1_C, dtype=float)


def contract_event_mapping_schema() -> None:
    n = 8
    _dxi, nodes, _faces, weights = q4.radial_geometry(n)

    def dense(times: np.ndarray) -> np.ndarray:
        times = np.asarray(times, dtype=float).reshape(-1)
        states = []
        for t in times:
            C = np.full(n + 1, 0.20 - 5.0e-5 * t)
            states.append(q4._to_log_state(q4.interleaved_state(np.full(n + 1, 40.0), C), n))
        return np.stack(states, axis=1)

    root, state, metadata = q4.locate_crossing(dense, 0.0, 2000.0, n)
    strict_before = float(np.max(np.exp(dense(np.array([root - 0.5]))[1::2, 0])))
    strict_after = float(np.max(np.exp(dense(np.array([root + 0.5]))[1::2, 0])))
    time_rows = q4.formal_time_set(123.4)
    check(
        "event_root_bracket_and_strict_sign",
        abs(root - 1000.0) <= 1.0e-8
        and metadata["bracket_width_s"] <= 1.0
        and abs(metadata["root_residual_kg_per_kg"]) <= 1.0e-9
        and strict_before > q4.THRESHOLD_C
        and strict_after < q4.THRESHOLD_C,
        root_time_s=root,
        bracket_s=metadata["bracket_s"],
        bracket_width_s=metadata["bracket_width_s"],
        residual=metadata["root_residual_kg_per_kg"],
        strict_before_C=strict_before,
        strict_after_C=strict_after,
    )
    mapped = q4.map_fixed_positions(
        nodes,
        np.full((1, n + 1), 0.4),
        np.array([0.013]),
        positions_cm=np.array([1.2, 1.3, 1.4]),
    )
    check(
        "fixed_position_domain_mask_and_surface_column",
        np.isfinite(mapped["moisture"][0, 0])
        and np.isfinite(mapped["moisture"][0, 1])
        and np.isnan(mapped["moisture"][0, 2])
        and np.isclose(mapped["surface_moisture"][0], 0.4),
        mapped_row=mapped["moisture"][0],
        surface=float(mapped["surface_moisture"][0]),
    )
    check(
        "formal_time_schema",
        np.array_equal(q4.formal_time_set(120.0), np.array([60.0, 120.0]))
        and np.array_equal(time_rows, np.array([60.0, 120.0, 123.4]))
        and 0.0 not in time_rows
        and len(np.unique(time_rows)) == len(time_rows),
        rows=time_rows,
    )


def contract_arbitrary_grid_position_mapping() -> None:
    checks: list[dict[str, Any]] = []
    for n in (4, 20, 40):
        _dxi, nodes, _faces, weights = q4.radial_geometry(n)
        T = 30.0 + 8.0 * nodes**2
        C = 0.2 + 0.3 * nodes**2
        transformed = q4._to_log_state(q4.interleaved_state(T, C), n, model="appendix4")[:, None]
        positions = np.array([0.0, 0.7, 1.3, 1.9, 2.0])
        mapped_C, mapped_T, surface_C, surface_T, _mean, _physical, mask = q4._sample_output(
            transformed,
            n,
            nodes,
            weights,
            np.array([0.013]),
            positions,
        )
        checks.append({
            "n": n,
            "T_shape": list(mapped_T.shape),
            "C_shape": list(mapped_C.shape),
            "mask": mask[0].tolist(),
            "surface_C": float(surface_C[0]),
            "surface_T": float(surface_T[0]),
        })
    passed = all(item["T_shape"] == [1, 5] and item["C_shape"] == [1, 5] and item["mask"] == [True, True, True, False, False] for item in checks)
    check(
        "arbitrary_grid_temperature_moisture_mapping",
        passed,
        grids=checks,
        contract="same physical position interpolation and synchronized mask for T/C",
    )
    _dxi, nodes, _faces, weights = q4.radial_geometry(20)
    noncentral = q4.interleaved_state(40.0 + nodes, 0.20 + 0.4 * nodes**2)
    audit = q4._audit_state(q4._to_log_state(noncentral, 20), 20, weights)
    check(
        "max_minus_center_audit_sign",
        audit["max_minus_center_kg_per_kg"] > 0.0,
        max_minus_center_kg_per_kg=audit["max_minus_center_kg_per_kg"],
    )


def contract_event_failure_and_multicross_guards() -> None:
    n = 4
    radius = lambda _t: q4.RADIUS_INITIAL_M
    environment = lambda _t: (50.0, 0.05)
    initial = q4.interleaved_state(np.full(n + 1, 30.0), np.full(n + 1, 0.1512))

    def source_ok(t: float, _radius: float, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return np.zeros_like(nodes), np.full_like(nodes, -1.0e-3)

    result = q4.integrate_segment(
        initial,
        n,
        0.0,
        3.0,
        environment,
        radius,
        model="appendix4",
        rtol=1.0e-9,
        atol=1.0e-11,
        max_step=0.5,
        chunk_seconds=1.0,
        output_times=np.array([1.0, 2.0]),
        detect_event=True,
        h=0.0,
        hm=0.0,
        source=source_ok,
    )
    event = result.get("event")
    cross_block_pass = (
        event is not None
        and event["before_state_error"] is None
        and abs(float(event["strict_before_time_s"]) - (float(event["time_s"]) - 0.5)) <= 1.0e-12
        and event["strict_before_pass"]
        and event["strict_after_pass"]
    )
    check(
        "event_cross_block_strict_witness",
        cross_block_pass,
        event=None if event is None else {"time_s": event["time_s"], "before_state_error": event["before_state_error"], "certification_status": event["certification_status"]},
    )
    root_time = float(event["time_s"]) if event is not None else math.nan
    root_row_count = int(np.sum(np.isclose(np.asarray(result["times_s"], dtype=float), root_time, rtol=0.0, atol=1.0e-10))) if event is not None else 0
    check(
        "event_root_appended_to_unified_output",
        event is not None
        and root_row_count == 1
        and np.isclose(float(result["times_s"][-1]), root_time, rtol=0.0, atol=1.0e-10),
        root_time_s=root_time,
        root_row_count=root_row_count,
        output_times_s=np.asarray(result["times_s"], dtype=float),
    )

    original_solve_ivp = q4.solve_ivp

    def fail_post_root_solver(*args: Any, **kwargs: Any) -> Any:
        t_span = args[1] if len(args) > 1 else kwargs.get("t_span")
        if t_span is not None and float(t_span[0]) > 1.0 and float(t_span[1]) - float(t_span[0]) < 0.75:
            raise RuntimeError("intentional post-root witness failure")
        return original_solve_ivp(*args, **kwargs)

    q4.solve_ivp = fail_post_root_solver  # type: ignore[assignment]
    try:
        failed = q4.integrate_segment(
            initial,
            n,
            0.0,
            3.0,
            environment,
            radius,
            model="appendix4",
            rtol=1.0e-9,
            atol=1.0e-11,
            max_step=0.5,
            chunk_seconds=1.0,
            output_times=np.array([1.0, 2.0]),
            detect_event=True,
            h=0.0,
            hm=0.0,
            source=source_ok,
        )
    finally:
        q4.solve_ivp = original_solve_ivp  # type: ignore[assignment]
    failed_event = failed.get("event")
    check(
        "event_post_root_failure_blocks_certificate",
        failed_event is not None
        and failed_event.get("certification_status") == "FAIL"
        and failed_event.get("after_state_error") is not None,
        certification_status=None if failed_event is None else failed_event.get("certification_status"),
        after_state_error=None if failed_event is None else failed_event.get("after_state_error"),
    )
    failed_output = DEV / "failed_event_save_probe" / str(RESULTS["run_id"])
    try:
        q4.save_result_files(failed_output, failed)
    except RuntimeError as exc:
        failed_save_blocked = "Event certificate failed" in str(exc)
        failed_save_message = str(exc)
    else:
        failed_save_blocked = False
        failed_save_message = None
    check(
        "failed_event_save_redirects_and_blocks",
        failed_save_blocked
        and (failed_output / "failed_event_evidence" / "q4_fields.npz").exists()
        and not (failed_output / "q4_fields.npz").exists(),
        redirected_path=str(failed_output / "failed_event_evidence"),
        exception=failed_save_message,
    )
    full_audit_failed_result = dict(result)
    full_audit_failed_result["event"] = dict(result["event"])
    full_audit_failed_result["g_evidence_summary"] = dict(result.get("g_evidence_summary", {}))
    full_audit_failed_result["g_evidence_summary"]["full_run_sequence_audit"] = {
        "status": "CHECKED",
        "prior_positive": True,
        "downward_sign_changes": 1,
        "upward_sign_changes": 1,
        "no_upward_recrossing": False,
    }
    full_audit_output = DEV / "failed_full_audit_save_probe" / str(RESULTS["run_id"])
    try:
        q4.save_result_files(full_audit_output, full_audit_failed_result)
    except RuntimeError as exc:
        full_audit_save_blocked = "Event certificate failed" in str(exc)
        full_audit_save_message = str(exc)
    else:
        full_audit_save_blocked = False
        full_audit_save_message = None
    check(
        "full_sequence_failure_redirects_and_blocks",
        full_audit_save_blocked
        and (full_audit_output / "failed_event_evidence" / "q4_fields.npz").exists()
        and not (full_audit_output / "q4_fields.npz").exists(),
        redirected_path=str(full_audit_output / "failed_event_evidence"),
        exception=full_audit_save_message,
    )

    sequence = q4._event_sequence_audit(np.array([0.0, 1.0, 2.0]), np.array([0.1, -0.1, 0.1]), event_time=1.0)
    check(
        "event_multicross_sample_guard",
        sequence["upward_sign_changes"] == 1 and not sequence["no_upward_recrossing"],
        sequence=sequence,
    )

    # Integral-level counterexample: with a smooth source the threshold is
    # crossed downward and upward between the endpoints of one coarse solve.
    # The accepted-step/midpoint sequence must reject that coarse evidence,
    # while an event solve with a halved step locates the first crossing.
    center = 1.2
    curvature = 1.0e-3
    initial_multi = q4.interleaved_state(
        np.full(n + 1, 30.0),
        np.full(n + 1, q4.THRESHOLD_C - 2.5e-4 + curvature * center**2),
    )

    def source_multi(t: float, _radius: float, nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return np.zeros_like(nodes), np.full_like(nodes, 2.0 * curvature * (float(t) - center))

    coarse = q4.integrate_segment(
        initial_multi,
        n,
        0.0,
        2.0,
        environment,
        radius,
        model="appendix4",
        rtol=1.0e-7,
        atol=1.0e-9,
        max_step=2.0,
        chunk_seconds=2.0,
        output_times=np.array([1.0, 2.0]),
        detect_event=False,
        h=0.0,
        hm=0.0,
        source=source_multi,
    )
    coarse_times = np.concatenate([coarse["accepted_monitor_times_s"], coarse["midpoint_monitor_times_s"]])
    coarse_g = np.concatenate([coarse["accepted_monitor_g_N"], coarse["midpoint_monitor_g_N"]])
    coarse_sequence = q4._event_sequence_audit(coarse_times, coarse_g, event_time=2.0)
    coarse_rejected = (
        coarse_sequence.get("downward_sign_changes", 0) >= 1
        and coarse_sequence.get("upward_sign_changes", 0) >= 1
        and not coarse_sequence.get("no_upward_recrossing", True)
    )
    coarse_event = q4.integrate_segment(
        initial_multi,
        n,
        0.0,
        2.0,
        environment,
        radius,
        model="appendix4",
        rtol=1.0e-7,
        atol=1.0e-9,
        max_step=2.0,
        chunk_seconds=2.0,
        output_times=np.array([1.0, 2.0]),
        detect_event=True,
        h=0.0,
        hm=0.0,
        source=source_multi,
    )
    fine_event = q4.integrate_segment(
        initial_multi,
        n,
        0.0,
        2.0,
        environment,
        radius,
        model="appendix4",
        rtol=1.0e-7,
        atol=1.0e-9,
        max_step=1.0,
        chunk_seconds=2.0,
        output_times=np.array([1.0, 2.0]),
        detect_event=True,
        h=0.0,
        hm=0.0,
        source=source_multi,
    )
    original_solve_ivp = q4.solve_ivp
    callback_state = {"hidden": False}

    def hide_first_event_callback(*args: Any, **kwargs: Any) -> Any:
        if not callback_state["hidden"] and kwargs.get("events") is not None:
            callback_state["hidden"] = True
            kwargs = dict(kwargs)
            kwargs["events"] = None
        return original_solve_ivp(*args, **kwargs)

    q4.solve_ivp = hide_first_event_callback  # type: ignore[assignment]
    try:
        recovered = q4.integrate_segment(
            initial_multi,
            n,
            0.0,
            2.0,
            environment,
            radius,
            model="appendix4",
            rtol=1.0e-7,
            atol=1.0e-9,
            max_step=2.0,
            chunk_seconds=2.0,
            output_times=np.array([2.0]),
            detect_event=True,
            h=0.0,
            hm=0.0,
            source=source_multi,
        )
    finally:
        q4.solve_ivp = original_solve_ivp  # type: ignore[assignment]
    recovered_event = recovered.get("event")
    check(
        "production_missed_callback_recovery",
        recovered_event is not None
        and recovered_event.get("certification_status") == "PASS"
        and int(recovered_event.get("monitor_recovery_attempts", 0)) >= 1
        and abs(float(recovered_event["time_s"]) - (center - 0.5)) <= 1.0e-4
        and recovered_event.get("refinement_pass", False)
        and recovered_event.get("sequence_audit", {}).get("prior_positive", False)
        and int(recovered_event.get("sequence_audit", {}).get("downward_sign_changes", 0)) == 1
        and int(recovered_event.get("sequence_audit", {}).get("upward_sign_changes", 0)) == 0,
        event=None if recovered_event is None else {
            "time_s": recovered_event.get("time_s"),
            "certification_status": recovered_event.get("certification_status"),
            "monitor_recovery_attempts": recovered_event.get("monitor_recovery_attempts"),
            "sequence_audit": recovered_event.get("sequence_audit"),
        },
    )
    coarse_root = coarse_event.get("event")
    fine_root = fine_event.get("event")
    refinement_consistent = (
        coarse_root is not None
        and fine_root is not None
        and abs(float(coarse_root["time_s"]) - float(fine_root["time_s"])) <= 1.0e-4
        and abs(float(coarse_root["time_s"]) - (center - 0.5)) <= 1.0e-4
    )
    check(
        "integral_level_multicross_step_refinement_guard",
        coarse_rejected and refinement_consistent,
        coarse_sequence=coarse_sequence,
        coarse_event_time_s=None if coarse_root is None else float(coarse_root["time_s"]),
        fine_event_time_s=None if fine_root is None else float(fine_root["time_s"]),
        expected_first_crossing_s=center - 0.5,
        coarse_evidence_rejected=coarse_rejected,
    )


def contract_floating_positive_event_endpoint() -> None:
    """Recover a callback endpoint with round-off-positive monitor g.

    This is a regression for the N=1280 formal run: solve_ivp returned the
    accepted event endpoint with a tiny positive residual.  The production
    locator may accept only that bounded floating-point case, must retain a
    <=1 s bracket, and must still reject a genuinely positive endpoint.
    """
    n = 4
    base_log = q4._to_log_state(q4.initial_state(n), n)
    endpoint_g = 5.0e-16

    def dense_near_zero(times: np.ndarray) -> np.ndarray:
        values = np.asarray(times, dtype=float).reshape(-1)
        states = np.repeat(base_log[:, None], len(values), axis=1)
        moisture = 0.15 + endpoint_g + 1.0e-3 * (1.0 - values)
        states[1::2, :] = np.log(moisture[None, :])
        return states

    root, _state, metadata = q4.locate_crossing(dense_near_zero, 0.0, 1.0, n)
    normalized_g = q4._normalize_event_monitor_g(endpoint_g)
    sequence = q4._event_sequence_audit(
        np.array([0.0, root], dtype=float),
        np.array([1.0e-3 + endpoint_g, normalized_g], dtype=float),
        event_time=root,
    )

    large_positive_g = 2.0e-9

    def dense_large_positive(times: np.ndarray) -> np.ndarray:
        values = np.asarray(times, dtype=float).reshape(-1)
        states = np.repeat(base_log[:, None], len(values), axis=1)
        moisture = 0.15 + large_positive_g + 1.0e-3 * (1.0 - values)
        states[1::2, :] = np.log(moisture[None, :])
        return states

    large_positive_rejected = False
    try:
        q4.locate_crossing(dense_large_positive, 0.0, 1.0, n)
    except ValueError:
        large_positive_rejected = True

    check(
        "floating_positive_event_endpoint_recovery",
        abs(root - 1.0) <= 1.0e-12
        and metadata["bracket_width_s"] <= 1.0
        and abs(metadata["root_residual_kg_per_kg"] - endpoint_g) <= 1.0e-15
        and normalized_g == 0.0
        and sequence["downward_sign_changes"] == 1
        and sequence["prior_positive"]
        and large_positive_rejected,
        root_time_s=root,
        bracket_width_s=metadata["bracket_width_s"],
        raw_endpoint_g=metadata["root_residual_kg_per_kg"],
        normalized_endpoint_g=normalized_g,
        downward_sign_changes=sequence["downward_sign_changes"],
        large_positive_g=large_positive_g,
        large_positive_rejected=large_positive_rejected,
    )


def contract_output_evidence_and_manifest(boundary: np.ndarray, radius: q4.RadiusModel) -> None:
    """Read back one short two-segment run and its per-block full-state evidence."""
    artifact_run_id = f"artifact_{uuid.uuid4().hex[:8]}"
    output = DEV / "artifact_contract_probe" / artifact_run_id
    evidence_base = DEV / "evidence_runs"
    result = q4.run_case(
        boundary,
        radius,
        n=4,
        model="appendix4",
        max_end_s=15000.0,
        rtol=1.0e-7,
        atol=1.0e-9,
        max_step_before_s=20.0,
        max_step_after_s=120.0,
        chunk_before_s=600.0,
        chunk_after_s=600.0,
        evidence_dir=evidence_base,
        run_id=artifact_run_id,
    )
    q4.save_result_files(output, result)
    fields_path = output / "q4_fields.npz"
    manifest_path = output / "q4_blocks_manifest.json"
    with np.load(fields_path) as fields:
        field_keys = set(fields.files)
        time_s = np.asarray(fields["time_s"], dtype=float)
        temp = np.asarray(fields["temperature_C"], dtype=float)
        moist = np.asarray(fields["moisture_kg_per_kg"], dtype=float)
        masks = np.asarray(fields["domain_mask"], dtype=bool)
        field_shape_pass = (
            temp.shape == moist.shape == masks.shape
            and temp.shape[0] == len(time_s)
            and np.array_equal(masks, np.isfinite(temp) & np.isfinite(moist))
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = manifest.get("segment_records", [])
    evidence_files = [Path(record["evidence_file"]) for record in records]
    evidence_exists = bool(records) and all(path.exists() for path in evidence_files)
    block_schema_pass = False
    block_continuity_pass = False
    block_state_max_difference = math.inf
    block_time_max_gap = math.inf
    switch_state_difference = math.inf
    switch_radius_difference = math.inf
    if evidence_exists:
        with np.load(evidence_files[0]) as block:
            required = {
                "accepted_time_s", "accepted_temperature_C", "accepted_moisture_kg_per_kg",
                "accepted_radius_m", "accepted_g_N", "midpoint_time_s", "midpoint_temperature_C",
                "midpoint_moisture_kg_per_kg", "midpoint_radius_m", "midpoint_g_N",
                "output_time_s", "output_g_N",
            }
            block_schema_pass = required.issubset(block.files) and (
                block["accepted_temperature_C"].shape[0] == len(block["accepted_time_s"])
                and block["accepted_moisture_kg_per_kg"].shape == block["accepted_temperature_C"].shape
                and len(block["accepted_radius_m"]) == len(block["accepted_time_s"])
                and len(block["accepted_g_N"]) == len(block["accepted_time_s"])
                and block["midpoint_moisture_kg_per_kg"].shape == block["midpoint_temperature_C"].shape
                and len(block["midpoint_radius_m"]) == len(block["midpoint_time_s"])
                and len(block["midpoint_g_N"]) == len(block["midpoint_time_s"])
                and len(block["output_g_N"]) == len(block["output_time_s"])
            )
        ordered = list(records)
        block_continuity_pass = all(
            np.isfinite(float(record["start_s"]))
            and float(record["actual_end_s"]) >= float(record["start_s"])
            and int(record["accepted_states"]) >= 2
            and int(record["midpoint_states"]) == int(record["accepted_states"]) - 1
            for record in ordered
        )
        state_differences: list[float] = []
        time_gaps: list[float] = []
        for previous_record, current_record in zip(ordered[:-1], ordered[1:]):
            with np.load(Path(previous_record["evidence_file"])) as previous_block, np.load(Path(current_record["evidence_file"])) as current_block:
                state_differences.append(max(
                    float(np.max(np.abs(previous_block["accepted_temperature_C"][-1] - current_block["accepted_temperature_C"][0]))),
                    float(np.max(np.abs(previous_block["accepted_moisture_kg_per_kg"][-1] - current_block["accepted_moisture_kg_per_kg"][0]))),
                ))
                time_gaps.append(abs(float(previous_block["accepted_time_s"][-1]) - float(current_block["accepted_time_s"][0])))
                if previous_record["segment"] != current_record["segment"]:
                    switch_state_difference = state_differences[-1]
                    switch_radius_difference = abs(float(previous_block["accepted_radius_m"][-1]) - float(current_block["accepted_radius_m"][0]))
        block_state_max_difference = max(state_differences) if state_differences else math.inf
        block_time_max_gap = max(time_gaps) if time_gaps else math.inf
        block_continuity_pass = (
            block_continuity_pass
            and block_state_max_difference <= 1.0e-12
            and block_time_max_gap <= 1.0e-10
            and switch_state_difference <= 1.0e-12
            and switch_radius_difference <= 1.0e-14
        )
    run_config = manifest.get("run_config", {})
    g_summary = result.get("g_evidence_summary", {})
    g_segments = g_summary.get("segments", {})
    g_summary_pass = (
        set(g_segments) == {"measured", "plateau"}
        and g_summary.get("accepted_count", 0) > 0
        and g_summary.get("midpoint_count", 0) > 0
        and g_summary.get("output_count", 0) > 0
        and g_segments["measured"].get("accepted_count", 0) > 0
        and g_segments["plateau"].get("accepted_count", 0) > 0
        and g_segments["measured"].get("output_count", 0) > 0
        and g_segments["plateau"].get("output_count", 0) > 0
        and g_summary.get("full_run_sequence_audit", {}).get("status") == "CHECKED"
    )
    try:
        q4._claim_evidence_directory(evidence_base, artifact_run_id)
    except FileExistsError as exc:
        unique_run_guard = True
        unique_run_guard_message = str(exc)
    else:
        unique_run_guard = False
        unique_run_guard_message = None
    check(
        "unified_fields_and_block_manifest",
        result.get("output_schema", {}).get("pass", False)
        and fields_path.exists()
        and manifest_path.exists()
        and {"time_s", "temperature_C", "moisture_kg_per_kg", "surface_moisture_kg_per_kg", "surface_temperature_C", "radius_m", "mean_moisture", "domain_mask"}.issubset(field_keys)
        and field_shape_pass
        and np.array_equal(time_s, q4.formal_time_set(15000.0))
        and evidence_exists
        and block_schema_pass
        and block_continuity_pass
        and g_summary_pass
        and run_config.get("run_id") == artifact_run_id
        and run_config.get("segments", {}).get("measured", {}).get("max_step_s") == 20.0
        and run_config.get("segments", {}).get("plateau", {}).get("max_step_s") == 120.0,
        output_directory=str(output),
        output_schema=result.get("output_schema"),
        field_shape=list(temp.shape),
        block_count=len(records),
        block_schema_pass=block_schema_pass,
        block_continuity_pass=block_continuity_pass,
        block_state_max_difference=block_state_max_difference,
        block_time_max_gap_s=block_time_max_gap,
        switch_state_difference=switch_state_difference,
        switch_radius_difference_m=switch_radius_difference,
        g_evidence_summary=g_summary,
        unique_run_guard=unique_run_guard,
        unique_run_guard_message=unique_run_guard_message,
        run_config_segments=run_config.get("segments"),
    )
    check(
        "explicit_evidence_directory_records_files",
        result.get("run_config", {}).get("evidence_directory") is not None
        and all(record.get("evidence_file") is not None for record in records),
        evidence_directory=run_config.get("evidence_directory"),
        record_evidence_files=[record.get("evidence_file") for record in records[:2]],
    )


def contract_independent_flux(boundary: np.ndarray, radius: q4.RadiusModel) -> None:
    n = 20
    dxi, nodes, faces, weights = q4.radial_geometry(n)
    measured = q4.measured_environment(boundary)
    first = q4.integrate_segment(
        q4.initial_state(n),
        n,
        0.0,
        10800.0,
        measured,
        radius,
        model="appendix4",
        rtol=q4.BASE_RTL,
        atol=q4.BASE_ATOL,
        max_step=q4.BASE_MAX_STEP_BEFORE_S,
        chunk_seconds=q4.DEFAULT_CHUNK_S,
        output_times=np.empty(0),
        detect_event=False,
    )
    rhs_log, jac_log, _ = q4.make_log_rhs(measured, radius, n, model="appendix4")
    trajectory = solve_ivp(
        rhs_log,
        (10800.0, 12600.0),
        q4._to_log_state(first["state_final"], n, model="appendix4"),
        method="BDF",
        jac=jac_log,
        jac_sparsity=q4.jacobian_sparsity(n),
        rtol=q4.BASE_RTL,
        atol=q4.BASE_ATOL,
        max_step=q4.BASE_MAX_STEP_BEFORE_S,
        dense_output=True,
    )
    if not trajectory.success:
        raise RuntimeError(trajectory.message)
    sample_times = np.array([10800.0 + 300.0 * m for m in range(7)])
    physical = q4._physical_from_transformed_matrix(trajectory.sol(sample_times), n)
    rhs, _geometry = q4.make_rhs(measured, radius, n, model="appendix4")
    errors_T: list[float] = []
    errors_C: list[float] = []
    for column, t in enumerate(sample_times):
        state = physical[:, column]
        temperature_nodes = np.empty(n + 1, dtype=float)
        moisture_nodes = np.empty(n + 1, dtype=float)
        for node_index in range(n + 1):
            temperature_nodes[node_index] = float(state[2 * node_index])
            moisture_nodes[node_index] = float(state[2 * node_index + 1])
        rho, cp, k, D = q4.material_properties(moisture_nodes, temperature_nodes, model="appendix4")
        env_T, env_C = measured(float(t))
        R = float(radius(float(t)))
        heat_flux = np.empty(n, dtype=float)
        moisture_flux = np.empty(n, dtype=float)
        for face_index in range(n):
            k_face = 0.5 * (float(k[face_index]) + float(k[face_index + 1]))
            D_face = 0.5 * (float(D[face_index]) + float(D[face_index + 1]))
            heat_flux[face_index] = float(faces[face_index]) * k_face * (temperature_nodes[face_index + 1] - temperature_nodes[face_index]) / dxi
            moisture_flux[face_index] = float(faces[face_index]) * D_face * (moisture_nodes[face_index + 1] - moisture_nodes[face_index]) / dxi
        heat_surface = -R * q4.H_BASE * (temperature_nodes[n] - env_T)
        moisture_surface = -R * q4.HM_BASE * (moisture_nodes[n] - env_C)
        heat_div = np.empty(n + 1, dtype=float)
        moisture_div = np.empty(n + 1, dtype=float)
        for control_index in range(n + 1):
            if control_index == 0:
                heat_div[control_index] = heat_flux[0] / float(weights[0])
                moisture_div[control_index] = moisture_flux[0] / float(weights[0])
            elif control_index == n:
                heat_div[control_index] = (heat_surface - heat_flux[n - 1]) / float(weights[n])
                moisture_div[control_index] = (moisture_surface - moisture_flux[n - 1]) / float(weights[n])
            else:
                heat_div[control_index] = (heat_flux[control_index] - heat_flux[control_index - 1]) / float(weights[control_index])
                moisture_div[control_index] = (moisture_flux[control_index] - moisture_flux[control_index - 1]) / float(weights[control_index])
        independent = np.empty_like(state)
        for node_index in range(n + 1):
            independent[2 * node_index] = heat_div[node_index] / (R * R * float(rho[node_index]) * float(cp[node_index]))
            independent[2 * node_index + 1] = moisture_div[node_index] / (R * R)
        difference = rhs(float(t), state) - independent
        errors_T.append(max(abs(float(difference[2 * node_index])) for node_index in range(n + 1)))
        errors_C.append(max(abs(float(difference[2 * node_index + 1])) for node_index in range(n + 1)))
    check(
        "independent_surface_and_internal_flux",
        max(errors_T) <= 1.0e-10 and max(errors_C) <= 1.0e-12,
        sample_times_s=sample_times,
        max_temperature_abs_error=errors_T,
        max_moisture_abs_error=errors_C,
        implementation="separate control-face loop on main-solution states",
        contract_interval_s=[10800.0, 12600.0],
    )


def contract_q3_regression(boundary: np.ndarray) -> None:
    fixed = _constant_radius()
    measured = q4.measured_environment(boundary)
    n_rhs = 20
    _dxi, nodes, _faces, _weights = q4.radial_geometry(n_rhs)
    native_state = q4.interleaved_state(
        42.0 + 3.0 * nodes**2 + 0.2 * np.sin(4.0 * nodes),
        0.18 + 0.25 * nodes**2 + 0.025 * np.cos(3.0 * nodes),
    )
    native_q4, _ = q4.make_rhs(measured, fixed, n_rhs, model="appendix3")
    native_q2, _ = q2.make_rhs(boundary, n_rhs)
    native_difference = native_q4(7200.0, native_state) - native_q2(7200.0, native_state)
    check(
        "q3_fixed_radius_native_rhs_equivalence",
        np.max(np.abs(native_difference)) <= 1.0e-12,
        max_absolute_difference=float(np.max(np.abs(native_difference))),
        threshold=1.0e-12,
        purpose="independent Q4 native finite-volume RHS remains equivalent to Q2 before compatibility integration",
    )
    diagnostic = q3.platform_diagnostic(boundary)
    tight_config = {
        "rtol": q4.TIGHT_RTL,
        "atol": q4.TIGHT_ATOL,
        "max_step_before_s": q4.TIGHT_MAX_STEP_BEFORE_S,
        "max_step_after_s": q4.TIGHT_MAX_STEP_AFTER_S,
    }
    base_config = {
        "rtol": q4.BASE_RTL,
        "atol": q4.BASE_ATOL,
        "max_step_before_s": q4.BASE_MAX_STEP_BEFORE_S,
        "max_step_after_s": q4.BASE_MAX_STEP_AFTER_S,
    }

    def run_pair(config: dict[str, float]) -> tuple[dict[str, Any], dict[str, Any]]:
        q3_result = q3.run_case(
            boundary,
            diagnostic,
            n=20,
            max_end_s=240000,
            rtol=config["rtol"],
            atol=config["atol"],
            max_step_before_s=config["max_step_before_s"],
            max_step_after_s=config["max_step_after_s"],
            chunk_before_s=int(q4.DEFAULT_CHUNK_S),
            chunk_after_s=int(q4.DEFAULT_CHUNK_S),
        )
        q4_result = q4.run_case(
            boundary,
            fixed,
            n=20,
            model="appendix3",
            max_end_s=240000,
            rtol=config["rtol"],
            atol=config["atol"],
            max_step_before_s=config["max_step_before_s"],
            max_step_after_s=config["max_step_after_s"],
            chunk_before_s=q4.DEFAULT_CHUNK_S,
            chunk_after_s=q4.DEFAULT_CHUNK_S,
        )
        return q3_result, q4_result

    def compare_pair(q3_result: dict[str, Any], q4_result: dict[str, Any]) -> dict[str, Any]:
        q3_times = np.asarray(q3_result["times_s"], dtype=float)
        q4_times = np.asarray(q4_result["times_s"], dtype=float)
        common = np.intersect1d(q3_times, q4_times)
        q3_index = {float(t): i for i, t in enumerate(q3_times)}
        q4_index = {float(t): i for i, t in enumerate(q4_times)}
        q3_T = np.asarray([q3_result["temperature_C"][q3_index[float(t)]] for t in common])
        q4_T = np.asarray([q4_result["temperature_C"][q4_index[float(t)]] for t in common])
        q3_C = np.asarray([q3_result["moisture_kg_per_kg"][q3_index[float(t)]] for t in common])
        q4_C = np.asarray([q4_result["moisture_kg_per_kg"][q4_index[float(t)]] for t in common])
        max_T = float(np.max(np.abs(q3_T - q4_T))) if len(common) else math.inf
        max_C = float(np.max(np.abs(q3_C - q4_C))) if len(common) else math.inf
        q3_event = q3_result.get("event")
        q4_event = q4_result.get("event")
        event_difference = math.inf if q3_event is None or q4_event is None else abs(float(q3_event["time_s"]) - float(q4_event["time_s"]))
        endpoint_T = endpoint_C = endpoint_mapped_C = endpoint_mapped_T = math.inf
        if q3_event is not None and q4_event is not None:
            endpoint_T = float(np.max(np.abs(np.asarray(q3_event["state"])[0::2] - np.asarray(q4_event["state"])[0::2])))
            endpoint_C = float(np.max(np.abs(np.asarray(q3_event["state"])[1::2] - np.asarray(q4_event["state"])[1::2])))
            q3_mapped = np.asarray(q3_event["state"][1::2], dtype=float)
            endpoint_mapped_C = float(np.max(np.abs(q3_mapped - np.asarray(q4_event.get("root_mapped_C", []), dtype=float))))
            q3_mapped_T = np.asarray(q3_event["state"][0::2], dtype=float)
            endpoint_mapped_T = float(np.max(np.abs(q3_mapped_T - np.asarray(q4_event.get("root_mapped_T_C", []), dtype=float))))
        return {
            "common_rows": int(len(common)),
            "common_end_s": float(common[-1]) if len(common) else None,
            "max_temperature_abs_difference_C": max_T,
            "max_moisture_abs_difference": max_C,
            "q3_event_time_s": None if q3_event is None else float(q3_event["time_s"]),
            "q4_event_time_s": None if q4_event is None else float(q4_event["time_s"]),
            "event_time_difference_s": event_difference,
            "endpoint_full_state_temperature_abs_difference_C": endpoint_T,
            "endpoint_full_state_moisture_abs_difference": endpoint_C,
            "endpoint_mapped_temperature_abs_difference_C": endpoint_mapped_T,
            "endpoint_mapped_moisture_abs_difference": endpoint_mapped_C,
            "q4_event_certification": None if q4_event is None else q4_event.get("certification_status"),
            "q4_full_run_sequence_audit": None if q4_event is None else q4_event.get("full_run_sequence_audit"),
            "q4_full_prior_positive": None if q4_event is None else q4_event.get("full_prior_positive"),
            "q4_full_downward_sign_changes": None if q4_event is None else q4_event.get("full_downward_sign_changes"),
            "q4_full_upward_sign_changes": None if q4_event is None else q4_event.get("full_upward_sign_changes"),
            "q4_full_no_upward_recrossing": None if q4_event is None else q4_event.get("full_no_upward_recrossing"),
            "output_schema_pass": bool(q4_result.get("output_schema", {}).get("pass", False)),
        }

    tight_q3, tight_q4 = run_pair(tight_config)
    tight_metrics = compare_pair(tight_q3, tight_q4)
    base_q3, base_q4 = run_pair(base_config)
    base_diagnostic = compare_pair(base_q3, base_q4)
    event_relative = (
        tight_metrics["event_time_difference_s"] / tight_metrics["q3_event_time_s"]
        if tight_metrics["q3_event_time_s"] is not None
        else math.inf
    )
    check(
        "q3_fixed_radius_appendix3_regression",
        tight_metrics["common_rows"] > 1000
        and tight_metrics["max_temperature_abs_difference_C"] <= 1.0e-6
        and tight_metrics["max_moisture_abs_difference"] <= 1.0e-7
        and tight_metrics["event_time_difference_s"] <= 1.0
        and event_relative <= 1.0e-6
        and tight_metrics["endpoint_full_state_temperature_abs_difference_C"] <= 1.0e-6
        and tight_metrics["endpoint_full_state_moisture_abs_difference"] <= 1.0e-7
        and tight_metrics["endpoint_mapped_moisture_abs_difference"] <= 1.0e-7
        and tight_metrics["endpoint_mapped_temperature_abs_difference_C"] <= 1.0e-6
        and tight_metrics["output_schema_pass"]
        and tight_metrics["q4_event_certification"] == "PASS"
        and tight_metrics["q4_full_run_sequence_audit"].get("status") == "CHECKED"
        and tight_metrics["q4_full_prior_positive"]
        and tight_metrics["q4_full_downward_sign_changes"] == 1
        and tight_metrics["q4_full_upward_sign_changes"] == 0
        and tight_metrics["q4_full_no_upward_recrossing"],
        # The event occurs after the 4-hour split in this external-reference
        # comparison, so its certificate must expose the merged two-segment
        # monitor audit rather than only the local plateau sequence.
        comparison_config=tight_config,
        tight_comparison=tight_metrics,
        tight_event_relative_difference=event_relative,
        base_diagnostic=base_diagnostic,
        hard_thresholds={"T_C": 1.0e-6, "C": 1.0e-7, "event_s": 1.0, "event_relative": 1.0e-6},
    )
    check(
        "no_implicit_evidence_without_directory",
        tight_q4.get("run_config", {}).get("evidence_directory") is None
        and all(record.get("evidence_file") is None for record in tight_q4.get("segment_records", [])),
        evidence_directory=tight_q4.get("run_config", {}).get("evidence_directory"),
        record_evidence_files=[record.get("evidence_file") for record in tight_q4.get("segment_records", [])[:2]],
    )


def contract_split_flux_balance() -> None:
    """The 4-hour discontinuity must retain one-sided fluxes."""
    initial_mean = 0.50
    first_times = np.array([0.0, q4.SWITCH_S])
    second_times = np.array([q4.SWITCH_S, 18000.0])
    first_rate = np.array([1.0e-5, 1.0e-5])
    second_rate = np.array([2.0e-5, 2.0e-5])
    output_times = np.array([60.0, q4.SWITCH_S, 15000.0, 18000.0])
    expected_mean = np.array([
        initial_mean - 60.0 * 1.0e-5,
        initial_mean - q4.SWITCH_S * 1.0e-5,
        initial_mean - q4.SWITCH_S * 1.0e-5 - (15000.0 - q4.SWITCH_S) * 2.0e-5,
        initial_mean - q4.SWITCH_S * 1.0e-5 - (18000.0 - q4.SWITCH_S) * 2.0e-5,
    ])
    result = {
        "first_segment": {
            "end_s_requested": q4.SWITCH_S,
            "flux_times_s": first_times,
            "flux_rate_mean_per_s": first_rate,
        },
        "second_segment": {
            "end_s_requested": 18000.0,
            "flux_times_s": second_times,
            "flux_rate_mean_per_s": second_rate,
        },
        "times_s": output_times,
        "mean_moisture": expected_mean,
    }
    balance = q4._effective_balance(result, initial_mean)
    check(
        "four_hour_one_sided_flux_split",
        balance.get("status") == "CHECKED"
        and float(balance["max_abs_residual"]) <= 1.0e-14
        and np.isclose(balance["split_flux_one_sided"]["left_rate_mean_per_s"], 1.0e-5)
        and np.isclose(balance["split_flux_one_sided"]["right_rate_mean_per_s"], 2.0e-5)
        and np.isclose(balance["split_flux_one_sided"]["jump_mean_per_s"], 1.0e-5),
        max_abs_residual=balance.get("max_abs_residual"),
        split_flux_one_sided=balance.get("split_flux_one_sided"),
        segment_checks=balance.get("segment_checks"),
    )
def contract_radau_interval(boundary: np.ndarray, radius: q4.RadiusModel) -> None:
    n = 20
    measured = q4.measured_environment(boundary)
    first = q4.integrate_segment(
        q4.initial_state(n),
        n,
        0.0,
        q4.SWITCH_S,
        measured,
        radius,
        model="appendix4",
        rtol=q4.BASE_RTL,
        atol=q4.BASE_ATOL,
        max_step=q4.BASE_MAX_STEP_BEFORE_S,
        chunk_seconds=q4.DEFAULT_CHUNK_S,
        output_times=np.empty(0),
        detect_event=False,
    )
    plateau = q4.plateau_environment()
    rhs, jac, _geometry = q4.make_log_rhs(plateau, radius, n, model="appendix4")
    initial = q4._to_log_state(first["state_final"], n, model="appendix4")
    kwargs = dict(
        rtol=q4.BASE_RTL,
        atol=q4.BASE_ATOL,
        max_step=60.0,
        dense_output=True,
        jac=jac,
        jac_sparsity=q4.jacobian_sparsity(n),
    )
    bdf = solve_ivp(rhs, (q4.SWITCH_S, 18000.0), initial, method="BDF", **kwargs)
    radau = solve_ivp(rhs, (q4.SWITCH_S, 18000.0), initial, method="Radau", **kwargs)
    sample = np.array([14400.0, 15300.0, 16200.0, 17100.0, 18000.0])
    bdf_state = q4._physical_from_transformed_matrix(bdf.sol(sample), n)
    radau_state = q4._physical_from_transformed_matrix(radau.sol(sample), n)
    diff_T = float(np.max(np.abs(bdf_state[0::2, :] - radau_state[0::2, :])))
    diff_C = float(np.max(np.abs(bdf_state[1::2, :] - radau_state[1::2, :])))
    check(
        "bdf_radau_representative_interval",
        bdf.success and radau.success and diff_T <= 5.0e-7 and diff_C <= 5.0e-8,
        interval_s=[q4.SWITCH_S, 18000.0],
        bdf_success=bool(bdf.success),
        radau_success=bool(radau.success),
        max_temperature_difference_C=diff_T,
        max_moisture_difference=diff_C,
        thresholds={"T_C": 5.0e-7, "C": 5.0e-8},
    )


def write_evidence() -> None:
    DEV.mkdir(parents=True, exist_ok=True)
    RESULTS["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    RESULTS["pass"] = not FAILURES
    RESULTS["failures"] = FAILURES
    RESULTS["artifacts"] = sorted(ARTIFACTS)
    q4.write_json(DEV / "q4_dev_validation.json", RESULTS)
    if ARTIFACTS:
        np.savez_compressed(DEV / "q4_manufactured_evidence.npz", **ARTIFACTS)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run named Q4 development contracts.")
    parser.add_argument("--formal", action="store_true", help="run the Q4_MODELING 9.1 high-grid gate; omitted for repair probes")
    args = parser.parse_args(argv)
    RESULTS["run_mode"] = "formal_gate" if args.formal else "repair_probe_with_quick_precheck"
    try:
        boundary, boundary_metadata = q4.read_boundary()
        radius, radius_metadata = q4.read_radius()
        contract_input_data(boundary, radius, radius_metadata)
        contract_geometry_and_purity()
        contract_jacobian(boundary, radius)
        contract_manufactured_solution(radius, formal=args.formal)
        contract_event_mapping_schema()
        contract_arbitrary_grid_position_mapping()
        contract_event_failure_and_multicross_guards()
        contract_floating_positive_event_endpoint()
        contract_output_evidence_and_manifest(boundary, radius)
        contract_independent_flux(boundary, radius)
        contract_q3_regression(boundary)
        contract_split_flux_balance()
        contract_radau_interval(boundary, radius)
        RESULTS["inputs"] = {"boundary": boundary_metadata, "radius": radius_metadata}
    except Exception as exc:  # noqa: BLE001 - persist the exact blocking cause
        RESULTS["fatal_exception"] = repr(exc)
        if "fatal_exception" not in FAILURES:
            FAILURES.append("fatal_exception")
    finally:
        write_evidence()
    print(json.dumps(_safe({"pass": not FAILURES, "failures": FAILURES}), ensure_ascii=False))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    raise SystemExit(main())
