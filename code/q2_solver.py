"""Question 2 variable-property coupled heat and moisture solver.

The solver is deliberately independent of q1_solver.py.  It uses a fixed-radius
axisymmetric finite-volume discretisation and a coupled BDF integration.  The
public functions are also used by test_q2_contracts.py for named contracts.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import shutil
import sys
import time
from copy import copy as copy_style
from pathlib import Path
from typing import Any

import numpy as np
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix, diags


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BOUNDARY = ROOT.parent / "附件1.xlsx"
DEFAULT_TEMPLATE = ROOT.parent / "A题" / "A题" / "附件" / "附件3" / "result2.xlsx"
DEFAULT_OUTPUT = ROOT / "results" / "q2"
OFFICIAL_RESULT2 = DEFAULT_TEMPLATE
TEMPLATE_BACKUP = ROOT / "source" / "result2_template_original.xlsx"

RADIUS_M = 0.02
LENGTH_M = 0.25
INITIAL_T_C = 28.0
INITIAL_C = 2.55
H_BASE = 25.0
HM_BASE = 8e-7
PAPER_SECONDS = np.array([1800, 3600, 5400, 7200, 9000, 10800], dtype=int)
PAPER_NODE_INDEX = np.array([0, 5, 10, 15, 20], dtype=int)
# Radial profiles in the four-panel figure must cover the same six paper times
# as Tables 3--4, including 9000 s (2.5 h).
PROFILE_SECONDS = np.array([1800, 3600, 5400, 7200, 9000, 10800], dtype=int)
OUTPUT_NODES = 21
DEFAULT_END_S = 10800
CONVERGENCE_THRESHOLD_T_C = 5e-5
CONVERGENCE_THRESHOLD_C = 5e-5
DEFAULT_CHUNK_S = 600
# The physical moisture field is C>0.  BDF is integrated in the reversible
# coordinate z=log(C), then converted back with C=exp(z) for every reported
# value.  This prevents an otherwise harmless Newton trial from crossing the
# singular boundary of Appendix 3.  No physical output is clipped.
MOISTURE_PARAMETERIZATION = "log(C) for BDF internal state; reported C=exp(log(C))"


class Q2DomainError(ValueError):
    """The state is outside the domain of the empirical properties."""


SUPPORTED_PROPERTY_MODELS = ("q2", "q1")


def _validate_property_model(property_model: str) -> str:
    """Validate and normalize the constitutive model selector.

    ``q2`` is the historical/default model.  ``q1`` is deliberately exposed
    only as a constitutive degradation path for independent reproducibility
    checks; it does not alter the published Q2 default.
    """
    if not isinstance(property_model, str) or property_model not in SUPPORTED_PROPERTY_MODELS:
        raise ValueError(
            f"property_model must be one of {SUPPORTED_PROPERTY_MODELS}, got {property_model!r}."
        )
    return property_model


class SegmentObserver:
    """Read-only BDF segment observer used by supplementary diagnostics.

    The observer is opt-in.  A callback can consume each segment immediately,
    which avoids retaining large dense-output arrays on fine grids.  With
    ``retain=True`` the compact metadata and dense callables are retained for
    small-grid contract tests.  Nothing returned here is fed back to the
    integrator or used to modify an accepted state.
    """

    def __init__(self, *, retain: bool = True, callback: Any | None = None) -> None:
        self.retain = bool(retain)
        self.callback = callback
        self.segments: list[dict[str, Any]] = []
        self.total_segments = 0

    def record(self, payload: dict[str, Any]) -> None:
        self.total_segments += 1
        if self.callback is not None:
            self.callback(payload)
        if self.retain:
            self.segments.append(payload)


def radial_geometry(n: int) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    """Return dr, node radii, face radii, and truncated radial volumes."""
    if int(n) != n or n < 20 or n % 20 != 0:
        raise ValueError("n must be an integer multiple of 20 and at least 20.")
    n = int(n)
    dr = RADIUS_M / n
    nodes = np.arange(n + 1, dtype=float) * dr
    faces = (np.arange(n, dtype=float) + 0.5) * dr
    edge = np.arange(n + 1, dtype=float)
    left = np.maximum(0.0, (edge - 0.5) * dr)
    right = np.minimum(RADIUS_M, (edge + 0.5) * dr)
    volumes = (right * right - left * left) / 2.0
    if not np.isclose(volumes.sum(), RADIUS_M * RADIUS_M / 2.0, rtol=0.0, atol=1e-15):
        raise ArithmeticError("Cylindrical control-volume weights do not close.")
    return dr, nodes, faces, volumes


def material_properties(
    C: np.ndarray | float,
    T_C: np.ndarray | float,
    property_model: str = "q2",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate one of the two constitutive models with strict domain checks.

    The default ``q2`` branch is unchanged from the formal solver.  ``q1``
    reproduces Q1's constant thermal properties and concentration-dependent
    diffusivity, allowing a solver-level degradation/reproduction audit.
    """
    property_model = _validate_property_model(property_model)
    c = np.asarray(C, dtype=float)
    t = np.asarray(T_C, dtype=float)
    if not np.isfinite(c).all() or not np.isfinite(t).all():
        raise Q2DomainError("Nonfinite temperature or moisture state.")
    if np.any(c <= 0.0):
        raise Q2DomainError("Appendix 3 requires strictly positive moisture C.")
    T_K = t + 273.15
    if np.any(T_K <= 0.0):
        raise Q2DomainError("Temperature must satisfy T + 273.15 > 0 K.")
    if property_model == "q1":
        rho = np.full_like(c, 820.0, dtype=float)
        cp = np.full_like(c, 2600.0, dtype=float)
        k = np.full_like(c, 0.36, dtype=float)
        D = 7.0e-9 * np.exp(-0.89 / c)
    else:
        rho = 650.0 + 128.0 * c
        cp = 1450.0 + 2736.0 * c / (c + 1.0)
        k = 0.21 + 0.38 * c / (c + 1.0)
        D = 2.4e-3 * np.exp(-0.45 / c) * np.exp(-3850.0 / T_K)
    if not all(np.isfinite(v).all() for v in (rho, cp, k, D)):
        raise Q2DomainError("Appendix 3 produced a nonfinite property.")
    return rho, cp, k, D


def property_derivatives(
    C: np.ndarray | float,
    T_C: np.ndarray | float,
    property_model: str = "q2",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return local derivatives used by audits and optional analytic Jacobians."""
    property_model = _validate_property_model(property_model)
    rho, cp, k, D = material_properties(C, T_C, property_model=property_model)
    c = np.asarray(C, dtype=float)
    t_k = np.asarray(T_C, dtype=float) + 273.15
    if property_model == "q1":
        rho_c = np.zeros_like(c, dtype=float)
        cp_c = np.zeros_like(c, dtype=float)
        k_c = np.zeros_like(c, dtype=float)
        D_c = D * 0.89 / c**2
        D_t = np.zeros_like(c, dtype=float)
    else:
        rho_c = np.full_like(c, 128.0, dtype=float)
        cp_c = 2736.0 / (c + 1.0) ** 2
        k_c = 0.38 / (c + 1.0) ** 2
        D_c = D * 0.45 / c**2
        D_t = D * 3850.0 / t_k**2
    return rho_c, cp_c, k_c, D_c, D_t


def read_boundary(
    path: str | Path,
    sheet: str | None = None,
    end: int | float = DEFAULT_END_S,
) -> tuple[np.ndarray, str]:
    """Read Attachment 1 and reject missing coverage or implicit extrapolation."""
    path = Path(path).expanduser().resolve()
    requested_end = float(end)
    if not path.exists():
        raise FileNotFoundError(f"Boundary file does not exist: {path}")
    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb[sheet] if sheet else wb.worksheets[0]
        rows = list(ws.values)
        header_index = None
        for idx, row in enumerate(rows[:20]):
            if len(row) < 3:
                continue
            normalized = [str(v).replace(" ", "") if v is not None else "" for v in row[:3]]
            if all(token in normalized[j] for j, token in enumerate(("时间", "温度", "水分浓度"))):
                header_index = idx
                break
        if header_index is None:
            raise ValueError("Expected columns 时间, 温度, 水分浓度.")
        records: list[list[float]] = []
        for row_number, row in enumerate(rows[header_index + 1 :], header_index + 2):
            if all(v is None for v in row):
                continue
            if len(row) < 3 or any(v is None or isinstance(v, bool) for v in row[:3]):
                raise ValueError(f"Invalid boundary row {row_number}.")
            try:
                values = [float(v) for v in row[:3]]
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid boundary row {row_number}.") from exc
            records.append(values)
        boundary = np.asarray(records, dtype=float)
        if boundary.ndim != 2 or len(boundary) < 2 or not np.isfinite(boundary).all():
            raise ValueError("At least two finite boundary rows are required.")
        if np.any(np.diff(boundary[:, 0]) <= 0.0):
            raise ValueError("Boundary times must be strictly increasing.")
        if boundary[0, 0] > 0.0 or boundary[-1, 0] < requested_end:
            raise ValueError(
                f"Boundary covers {boundary[0, 0]}--{boundary[-1, 0]} s, "
                f"but requested interval ends at {requested_end} s; no extrapolation."
            )
        if np.any(boundary[:, 1] <= -273.15) or np.any(boundary[:, 2] < 0.0):
            raise ValueError("Boundary contains an invalid temperature or concentration.")
        return boundary, ws.title
    finally:
        wb.close()


def environment_at(t: float, boundary: np.ndarray) -> tuple[float, float]:
    """Piecewise-linear environment evaluation with an explicit range guard."""
    if t < boundary[0, 0] - 1e-8 or t > boundary[-1, 0] + 1e-8:
        raise ValueError("Requested time is outside the supplied boundary; no extrapolation.")
    t_checked = min(max(float(t), float(boundary[0, 0])), float(boundary[-1, 0]))
    return (
        float(np.interp(t_checked, boundary[:, 0], boundary[:, 1])),
        float(np.interp(t_checked, boundary[:, 0], boundary[:, 2])),
    )


def jacobian_sparsity(n: int):
    """Block-tridiagonal sparsity mask for interleaved ``(T_i,C_i)`` nodes.

    A neighbouring 2x2 node block reaches offsets +/-1, +/-2 and +/-3 in the
    interleaved scalar ordering.  Including +/-3 is essential for the
    cross-component dependency from node ``i`` to ``C_(i+1)`` (or ``T_(i-1)``).
    """
    size = 2 * (n + 1)
    offsets = list(range(-3, 4))
    diagonals = [np.ones(size - abs(offset), dtype=float) for offset in offsets]
    return diags(diagonals, offsets=offsets, shape=(size, size), format="csc")


def initial_state(
    n: int,
    temperature_C: float = INITIAL_T_C,
    moisture: float = INITIAL_C,
    property_model: str = "q2",
) -> np.ndarray:
    """Create an interleaved uniform state."""
    property_model = _validate_property_model(property_model)
    if moisture <= 0.0:
        raise Q2DomainError("Initial moisture must be positive.")
    _, _, _, _ = material_properties(
        np.array([moisture]), np.array([temperature_C]), property_model=property_model
    )
    y = np.empty(2 * (n + 1), dtype=float)
    y[0::2] = temperature_C
    y[1::2] = moisture
    return y


def _fluxes(
    temperature_C: np.ndarray,
    moisture: np.ndarray,
    n: int,
    dr: float,
    faces: np.ndarray,
    property_model: str = "q2",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rho, cp, k, D = material_properties(
        moisture, temperature_C, property_model=property_model
    )
    k_face = 0.5 * (k[:-1] + k[1:])
    D_face = 0.5 * (D[:-1] + D[1:])
    heat_flux = faces * k_face * np.diff(temperature_C) / dr
    moisture_flux = faces * D_face * np.diff(moisture) / dr
    return rho, cp, heat_flux, moisture_flux, D


def make_rhs(
    boundary: np.ndarray,
    n: int,
    h: float = H_BASE,
    hm: float = HM_BASE,
    property_model: str = "q2",
):
    """Build the coupled finite-volume right-hand side and geometry."""
    property_model = _validate_property_model(property_model)
    if h <= 0.0 or hm <= 0.0:
        raise ValueError("h and hm must be positive.")
    dr, nodes, faces, volumes = radial_geometry(n)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        if y.size != 2 * (n + 1):
            raise ValueError("Interleaved state has the wrong size.")
        temperature_C = np.asarray(y[0::2], dtype=float)
        moisture = np.asarray(y[1::2], dtype=float)
        T_environment, C_environment = environment_at(float(t), boundary)
        rho, cp, heat_flux, moisture_flux, _ = _fluxes(
            temperature_C, moisture, n, dr, faces, property_model=property_model
        )
        heat_divergence = np.empty(n + 1, dtype=float)
        moisture_divergence = np.empty(n + 1, dtype=float)
        heat_divergence[0] = heat_flux[0] / volumes[0]
        moisture_divergence[0] = moisture_flux[0] / volumes[0]
        if n > 1:
            heat_divergence[1:-1] = (
                heat_flux[1:] - heat_flux[:-1]
            ) / volumes[1:-1]
            moisture_divergence[1:-1] = (
                moisture_flux[1:] - moisture_flux[:-1]
            ) / volumes[1:-1]
        heat_surface = -RADIUS_M * h * (temperature_C[-1] - T_environment)
        moisture_surface = -RADIUS_M * hm * (moisture[-1] - C_environment)
        heat_divergence[-1] = (heat_surface - heat_flux[-1]) / volumes[-1]
        moisture_divergence[-1] = (moisture_surface - moisture_flux[-1]) / volumes[-1]
        result = np.empty_like(y, dtype=float)
        result[0::2] = heat_divergence / (rho * cp)
        result[1::2] = moisture_divergence
        return result

    return rhs, (dr, nodes, faces, volumes)


def assemble_rhs(
    y: np.ndarray,
    t: float,
    boundary: np.ndarray,
    n: int,
    h: float = H_BASE,
    hm: float = HM_BASE,
    property_model: str = "q2",
) -> np.ndarray:
    """Evaluate one RHS value, mainly for targeted contract tests."""
    rhs, _ = make_rhs(boundary, n, h=h, hm=hm, property_model=property_model)
    return rhs(t, y)


def make_log_rhs(
    boundary: np.ndarray,
    n: int,
    h: float = H_BASE,
    hm: float = HM_BASE,
    property_model: str = "q2",
):
    """Build the RHS used by BDF in ``(T, z=log(C))`` coordinates.

    The physical finite-volume equations remain in :func:`make_rhs`.  Only
    the integration coordinate is changed: if ``g_C`` is the physical water
    derivative, then ``z_t=g_C/C`` and ``C=exp(z)`` exactly.  This is a
    reversible positive parameterisation, not a concentration floor or a
    post-step correction.  The same five-diagonal block sparsity mask remains
    valid because each node's ``(T_i,z_i)`` still depends only on the node and
    its two neighbours; in particular the D(C,T) cross-dependencies are
    preserved.
    """
    property_model = _validate_property_model(property_model)
    physical_rhs, geometry = make_rhs(
        boundary, n, h=h, hm=hm, property_model=property_model
    )
    dr, nodes, faces, volumes = geometry

    def _physical_from_log(z_state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        z_state = np.asarray(z_state, dtype=float)
        if z_state.size != 2 * (n + 1) or not np.isfinite(z_state).all():
            raise Q2DomainError("Nonfinite transformed temperature or moisture state.")
        # Do not clip or floor exp(z): underflow to zero is an explicit domain
        # failure, while normal positive values are carried exactly into the
        # physical RHS.
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            try:
                moisture = np.exp(z_state[1::2])
            except FloatingPointError as exc:
                raise Q2DomainError("log(C) trial overflowed; no clipping was applied.") from exc
        if not np.isfinite(moisture).all() or np.any(moisture <= 0.0):
            raise Q2DomainError("log(C) trial did not map to a strictly positive C.")
        physical_state = np.empty_like(z_state)
        physical_state[0::2] = z_state[0::2]
        physical_state[1::2] = moisture
        return physical_state, moisture

    def rhs_log(t: float, z_state: np.ndarray) -> np.ndarray:
        physical_state, moisture = _physical_from_log(z_state)
        physical_derivative = physical_rhs(t, physical_state)
        transformed_derivative = np.empty_like(z_state)
        transformed_derivative[0::2] = physical_derivative[0::2]
        transformed_derivative[1::2] = physical_derivative[1::2] / moisture
        if not np.isfinite(transformed_derivative).all():
            raise Q2DomainError("Nonfinite derivative in log(C) coordinates.")
        return transformed_derivative

    def jac_log(t: float, z_state: np.ndarray):
        """Analytic sparse Jacobian of ``rhs_log`` in interleaved z coordinates."""
        physical_state, moisture = _physical_from_log(z_state)
        temperature_C = physical_state[0::2]
        rho, cp, k, D = material_properties(
            moisture, temperature_C, property_model=property_model
        )
        rho_c, cp_c, k_c, D_c, D_t = property_derivatives(
            moisture, temperature_C, property_model=property_model
        )
        T_environment, C_environment = environment_at(float(t), boundary)
        delta_T = np.diff(temperature_C)
        delta_C = np.diff(moisture)
        k_face = 0.5 * (k[:-1] + k[1:])
        D_face = 0.5 * (D[:-1] + D[1:])
        heat_flux = faces * k_face * delta_T / dr
        moisture_flux = faces * D_face * delta_C / dr
        heat_surface = -RADIUS_M * h * (temperature_C[-1] - T_environment)
        moisture_surface = -RADIUS_M * hm * (moisture[-1] - C_environment)

        # Divergence values used for the diagonal quotient-rule terms.
        heat_divergence = np.empty(n + 1, dtype=float)
        moisture_divergence = np.empty(n + 1, dtype=float)
        heat_divergence[0] = heat_flux[0] / volumes[0]
        moisture_divergence[0] = moisture_flux[0] / volumes[0]
        if n > 1:
            heat_divergence[1:-1] = (heat_flux[1:] - heat_flux[:-1]) / volumes[1:-1]
            moisture_divergence[1:-1] = (moisture_flux[1:] - moisture_flux[:-1]) / volumes[1:-1]
        heat_divergence[-1] = (heat_surface - heat_flux[-1]) / volumes[-1]
        moisture_divergence[-1] = (moisture_surface - moisture_flux[-1]) / volumes[-1]

        # Face derivatives are with respect to physical (T,C).  The C-column
        # is multiplied by C_j below for the z_j=log(C_j) coordinate.
        heat_T_left = -faces * k_face / dr
        heat_T_right = -heat_T_left
        heat_C_left = faces * (0.5 * k_c[:-1]) * delta_T / dr
        heat_C_right = faces * (0.5 * k_c[1:]) * delta_T / dr
        moisture_T_left = faces * (0.5 * D_t[:-1]) * delta_C / dr
        moisture_T_right = faces * (0.5 * D_t[1:]) * delta_C / dr
        moisture_C_left = faces * ((0.5 * D_c[:-1]) * delta_C - D_face) / dr
        moisture_C_right = faces * ((0.5 * D_c[1:]) * delta_C + D_face) / dr

        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []

        def add(row: int, col: int, value: float) -> None:
            if value != 0.0:
                rows.append(row)
                cols.append(col)
                data.append(float(value))

        for i in range(n + 1):
            # Heat divergence H_i and moisture divergence M_i are sums of
            # signed face terms.  Each term contributes to both T_j and C_j,
            # retaining the adjacent-node heat/moisture cross blocks.
            if i == 0:
                heat_terms = [
                    (0, heat_T_left[0], heat_C_left[0]),
                    (1, heat_T_right[0], heat_C_right[0]),
                ]
                moisture_terms = [
                    (0, moisture_T_left[0], moisture_C_left[0]),
                    (1, moisture_T_right[0], moisture_C_right[0]),
                ]
            elif i < n:
                heat_terms = [
                    (i, heat_T_left[i], heat_C_left[i]),
                    (i + 1, heat_T_right[i], heat_C_right[i]),
                    (i - 1, -heat_T_left[i - 1], -heat_C_left[i - 1]),
                    (i, -heat_T_right[i - 1], -heat_C_right[i - 1]),
                ]
                moisture_terms = [
                    (i, moisture_T_left[i], moisture_C_left[i]),
                    (i + 1, moisture_T_right[i], moisture_C_right[i]),
                    (i - 1, -moisture_T_left[i - 1], -moisture_C_left[i - 1]),
                    (i, -moisture_T_right[i - 1], -moisture_C_right[i - 1]),
                ]
            else:
                heat_terms = [
                    (n - 1, -heat_T_left[-1], -heat_C_left[-1]),
                    (n, -RADIUS_M * h - heat_T_right[-1], -heat_C_right[-1]),
                ]
                moisture_terms = [
                    (n - 1, -moisture_T_left[-1], -moisture_C_left[-1]),
                    (n, -moisture_T_right[-1], -RADIUS_M * hm - moisture_C_right[-1]),
                ]

            heat_capacity = rho[i] * cp[i] * volumes[i]
            heat_capacity_c = volumes[i] * (rho_c[i] * cp[i] + rho[i] * cp_c[i])
            for j, dH_dT, dH_dC in heat_terms:
                add(2 * i, 2 * j, dH_dT / heat_capacity)
                add(2 * i, 2 * j + 1, dH_dC * moisture[j] / heat_capacity)
            # Local heat-capacity quotient rule in the z_i column.
            add(
                2 * i,
                2 * i + 1,
                -heat_divergence[i]
                * volumes[i]
                * heat_capacity_c
                * moisture[i]
                / heat_capacity**2,
            )

            moisture_scale = volumes[i] * moisture[i]
            for j, dM_dT, dM_dC in moisture_terms:
                add(2 * i + 1, 2 * j, dM_dT / moisture_scale)
                add(2 * i + 1, 2 * j + 1, dM_dC * moisture[j] / moisture_scale)
            # Quotient rule for z_t=(C_t/C), whose local extra term is -z_t.
            add(2 * i + 1, 2 * i + 1, -moisture_divergence[i] / moisture[i])

        jacobian = coo_matrix(
            (np.asarray(data), (np.asarray(rows), np.asarray(cols))),
            shape=(2 * (n + 1), 2 * (n + 1)),
        ).tocsc()
        if not np.isfinite(jacobian.data).all():
            raise Q2DomainError("Nonfinite analytic Jacobian in log(C) coordinates.")
        return jacobian

    return rhs_log, jac_log, geometry


def _trapezoid_node_weights(nodes: np.ndarray) -> np.ndarray:
    """Return diagnostic trapezoid weights for an increasing node sequence."""
    nodes = np.asarray(nodes, dtype=float)
    if nodes.ndim != 1 or nodes.size == 0 or not np.isfinite(nodes).all():
        raise ValueError("Observer nodes must be a nonempty finite 1-D array.")
    if nodes.size == 1:
        return np.zeros(1, dtype=float)
    if np.any(np.diff(nodes) <= 0.0):
        raise ValueError("Observer nodes must be strictly increasing.")
    weights = np.empty_like(nodes)
    weights[0] = 0.5 * (nodes[1] - nodes[0])
    weights[-1] = 0.5 * (nodes[-1] - nodes[-2])
    weights[1:-1] = 0.5 * (nodes[2:] - nodes[:-2])
    return weights


def _make_dense_component_solution(dense_solution: Any, component: int):
    """Build a scalar dense evaluator without materialising all state rows.

    SciPy's BDF dense interpolants expose the polynomial coefficient table
    ``D``.  Using one row is important for the N=10240 balance diagnostic:
    surface flux quadrature needs only the surface concentration, not a large
    state matrix at every Gauss point.  A conservative full-state fallback is
    kept for alternative SciPy dense-output implementations.
    """
    try:
        ts = np.asarray(dense_solution.ts, dtype=float)
        interpolants = tuple(dense_solution.interpolants)
        ascending = bool(ts[-1] >= ts[0])
        ts_sorted = ts if ascending else ts[::-1]
        side = "left" if ascending else "right"
    except AttributeError:
        ts = ts_sorted = None
        interpolants = ()
        ascending = True
        side = "left"

    def one(value: float) -> float:
        value = float(value)
        if not interpolants:
            state = np.asarray(dense_solution(value), dtype=float)
            return float(state[component])
        index = int(np.searchsorted(ts_sorted, value, side=side))
        segment = min(max(index - 1, 0), len(interpolants) - 1)
        if not ascending:
            segment = len(interpolants) - 1 - segment
        interpolant = interpolants[segment]
        if all(hasattr(interpolant, name) for name in ("D", "t_shift", "denom")):
            x = (value - interpolant.t_shift) / interpolant.denom
            p = np.cumprod(x)
            return float(np.dot(interpolant.D[1:, component], p) + interpolant.D[0, component])
        state = np.asarray(interpolant(value), dtype=float)
        return float(state[component])

    def evaluate(query: float | np.ndarray) -> float | np.ndarray:
        values = np.asarray(query, dtype=float)
        if values.ndim == 0:
            return one(float(values))
        result = np.fromiter((one(value) for value in values.ravel()), dtype=float, count=values.size)
        return result.reshape(values.shape)

    return evaluate


def _notify_segment_observer(observer: Any, payload: dict[str, Any]) -> None:
    """Send a read-only segment payload to an opt-in observer/callback."""
    if observer is None:
        return
    if hasattr(observer, "record"):
        observer.record(payload)
    elif callable(observer):
        observer(payload)
    else:
        raise TypeError("segment_observer must provide record(payload) or be callable.")


def integrate_case(
    boundary: np.ndarray,
    n: int,
    end: int = DEFAULT_END_S,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    max_step: float = 5.0,
    h: float = H_BASE,
    hm: float = HM_BASE,
    initial: np.ndarray | None = None,
    start_time: int = 0,
    chunk_seconds: int = DEFAULT_CHUNK_S,
    property_model: str = "q2",
    segment_observer: Any | None = None,
) -> dict[str, Any]:
    """Integrate in chunks, retaining only 21-point history and final state.

    ``property_model`` and ``segment_observer`` are appended keyword
    arguments, preserving all historical positional calls.  The observer is
    strictly read-only and is invoked only after a segment succeeds.
    """
    property_model = _validate_property_model(property_model)
    if end < start_time or int(end) != end or int(start_time) != start_time:
        raise ValueError("end and start_time must be ordered integer seconds.")
    if end > boundary[-1, 0] or start_time < boundary[0, 0]:
        raise ValueError("Boundary does not cover requested integration interval.")
    if rtol <= 0.0 or atol <= 0.0 or max_step <= 0.0 or chunk_seconds <= 0:
        raise ValueError("Integration tolerances, max_step, and chunk_seconds must be positive.")
    n = int(n)
    # Keep the physical RHS available for audits/contracts, but integrate the
    # moisture component in its reversible log coordinate.  A finite
    # difference/Newton trial in raw C can otherwise cross C=0 at fine grids
    # even when every accepted physical state is valid.
    rhs_log, jac_log, geometry = make_log_rhs(
        boundary, n, h=h, hm=hm, property_model=property_model
    )
    dr, nodes, faces, volumes = geometry
    state = (
        initial_state(n, property_model=property_model)
        if initial is None
        else np.asarray(initial, dtype=float).copy()
    )
    if state.size != 2 * (n + 1):
        raise ValueError("Initial state has the wrong size for this grid.")
    material_properties(state[1::2], state[0::2], property_model=property_model)
    with np.errstate(divide="raise", invalid="raise"):
        try:
            transformed_state = state.copy()
            transformed_state[1::2] = np.log(state[1::2])
        except FloatingPointError as exc:
            raise Q2DomainError("Initial moisture cannot be represented as log(C).") from exc
    times = np.arange(int(start_time), int(end) + 1, dtype=int)
    output_count = len(times)
    sample_nodes = np.arange(0, n + 1, n // 20, dtype=int)
    if len(sample_nodes) != OUTPUT_NODES or sample_nodes[-1] != n:
        raise ArithmeticError("Output grid does not contain exactly 21 nodes.")
    temperature_history = np.empty((output_count, OUTPUT_NODES), dtype=float)
    moisture_history = np.empty((output_count, OUTPUT_NODES), dtype=float)
    mean_moisture = np.empty(output_count, dtype=float)
    nfev = njev = nlu = 0
    positivity_retries = 0
    smallest_step_used = float(max_step)
    started = time.perf_counter()
    cursor = float(start_time)
    out_index = 0
    initial_mean = float(np.dot(state[1::2], volumes) / (RADIUS_M * RADIUS_M / 2.0))
    while cursor < float(end) - 1e-12:
        segment_end = min(cursor + float(chunk_seconds), float(end))
        step_used = float(max_step)
        while True:
            try:
                solution = solve_ivp(
                    rhs_log,
                    (cursor, segment_end),
                    transformed_state,
                    method="BDF",
                    jac=jac_log,
                    jac_sparsity=jacobian_sparsity(n),
                    rtol=rtol,
                    atol=atol,
                    max_step=step_used,
                    dense_output=True,
                )
                break
            except Q2DomainError:
                # A Newton trial state can leave C>0 even when the accepted
                # physical state is positive.  Retry the same segment with a
                # shorter BDF step; accepted/output states remain strict and
                # are never clipped to a hidden concentration floor.
                if step_used <= max_step / 64.0:
                    raise
                step_used /= 2.0
                positivity_retries += 1
        smallest_step_used = min(smallest_step_used, step_used)
        if not solution.success:
            raise RuntimeError(f"BDF integration failed at {cursor}--{segment_end}: {solution.message}")
        nfev += solution.nfev
        njev += solution.njev or 0
        nlu += solution.nlu or 0
        integer_times = np.arange(int(round(cursor)), int(round(segment_end)) + 1, dtype=int)
        transformed_values = solution.sol(integer_times.astype(float))
        if not np.isfinite(transformed_values).all():
            raise ArithmeticError("Nonfinite transformed state returned by BDF.")
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            try:
                physical_values = transformed_values.copy()
                physical_values[1::2, :] = np.exp(transformed_values[1::2, :])
            except FloatingPointError as exc:
                raise Q2DomainError("BDF log(C) output overflowed; no clipping was applied.") from exc
        temperatures = physical_values[0::2, :].T
        moistures = physical_values[1::2, :].T
        if not np.isfinite(physical_values).all():
            raise ArithmeticError("Nonfinite physical state returned by BDF.")
        if np.any(moistures <= 0.0):
            raise Q2DomainError("BDF log(C) output was not strictly positive; no clipping was applied.")
        if np.any(temperatures + 273.15 <= 0.0):
            raise Q2DomainError("BDF produced a nonphysical absolute temperature.")
        first = int(integer_times[0] - int(start_time))
        last = int(integer_times[-1] - int(start_time)) + 1
        temperature_history[first:last, :] = temperatures[:, sample_nodes]
        moisture_history[first:last, :] = moistures[:, sample_nodes]
        mean_moisture[first:last] = moistures @ volumes / (RADIUS_M * RADIUS_M / 2.0)
        transformed_state = solution.y[:, -1].copy()
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            try:
                state = transformed_state.copy()
                state[1::2] = np.exp(transformed_state[1::2])
            except FloatingPointError as exc:
                raise Q2DomainError("BDF log(C) final state overflowed; no clipping was applied.") from exc
        if not np.isfinite(state).all() or np.any(state[1::2] <= 0.0):
            raise Q2DomainError("BDF log(C) final state was not strictly positive.")
        if segment_observer is not None:
            accepted_nodes = np.asarray(solution.t, dtype=float).copy()
            accepted_nodes.setflags(write=False)
            accepted_step_trapezoid_weights = _trapezoid_node_weights(accepted_nodes)
            accepted_step_trapezoid_weights.setflags(write=False)
            boundary_mask = (
                (boundary[:, 0] >= cursor - 1e-10)
                & (boundary[:, 0] <= segment_end + 1e-10)
            )
            boundary_nodes = np.asarray(boundary[boundary_mask, 0], dtype=float)
            boundary_nodes = np.unique(
                np.concatenate(([float(cursor), float(segment_end)], boundary_nodes))
            )
            boundary_nodes.setflags(write=False)
            control_volume_weights = np.asarray(volumes, dtype=float).copy()
            control_volume_weights.setflags(write=False)
            radius_nodes = np.asarray(nodes, dtype=float).copy()
            radius_nodes.setflags(write=False)
            dense_solution = solution.sol

            def dense_log_solution(
                query: float | np.ndarray,
                dense_solution: Any = dense_solution,
            ) -> np.ndarray:
                """Expose the accepted segment's unmodified (T, log C) path."""
                values = np.asarray(dense_solution(query), dtype=float)
                if not np.isfinite(values).all():
                    raise Q2DomainError("Observer dense log state is nonfinite.")
                return values

            def physical_dense_solution(
                query: float | np.ndarray,
                dense_solution: Any = dense_solution,
            ) -> np.ndarray:
                transformed = np.asarray(dense_solution(query), dtype=float)
                with np.errstate(over="raise", invalid="raise", under="ignore"):
                    try:
                        physical = transformed.copy()
                        if physical.ndim == 1:
                            physical[1::2] = np.exp(transformed[1::2])
                        else:
                            physical[1::2, ...] = np.exp(transformed[1::2, ...])
                    except FloatingPointError as exc:
                        raise Q2DomainError(
                            "Observer dense log(C) output overflowed; no clipping was applied."
                        ) from exc
                if not np.isfinite(physical).all() or np.any(physical[1::2, ...] <= 0.0):
                    raise Q2DomainError(
                        "Observer dense output was not strictly positive; no clipping was applied."
                    )
                return physical

            moisture_surface_log_solution = _make_dense_component_solution(
                dense_solution, 2 * n + 1
            )
            def moisture_surface_solution(
                query: float | np.ndarray,
                component_solution: Any = moisture_surface_log_solution,
            ) -> np.ndarray:
                with np.errstate(over="raise", invalid="raise", under="ignore"):
                    values = np.exp(component_solution(query))
                if not np.isfinite(values).all() or np.any(values <= 0.0):
                    raise Q2DomainError(
                        "Observer surface moisture was not strictly positive; no clipping was applied."
                    )
                return values

            payload = {
                "start_s": float(cursor),
                "end_s": float(segment_end),
                "accepted_step_nodes": accepted_nodes,
                "accepted_step_trapezoid_weights": accepted_step_trapezoid_weights,
                "boundary_nodes": boundary_nodes,
                "radius_nodes_m": radius_nodes,
                "control_volume_weights_m2": control_volume_weights,
                "h_W_m2K": float(h),
                "hm_m_s": float(hm),
                "property_model": property_model,
                "dense_log_solution": dense_log_solution,
                "dense_solution": physical_dense_solution,
                "temperature_solution": _make_dense_component_solution(dense_solution, 0),
                "moisture_surface_log_solution": moisture_surface_log_solution,
                "moisture_surface_solution": moisture_surface_solution,
                "nfev": int(solution.nfev),
                "njev": int(solution.njev or 0),
                "nlu": int(solution.nlu or 0),
            }
            _notify_segment_observer(segment_observer, payload)
        cursor = segment_end
        out_index = last
    if output_count and not np.isfinite(mean_moisture).all():
        raise ArithmeticError("Mean moisture history was not fully populated.")
    env_temp = np.interp(times.astype(float), boundary[:, 0], boundary[:, 1])
    env_moisture = np.interp(times.astype(float), boundary[:, 0], boundary[:, 2])
    outward_rate = RADIUS_M * hm * (moisture_history[:, -1] - env_moisture)
    cumulative_loss = np.zeros(output_count, dtype=float)
    if output_count > 1:
        cumulative_loss[1:] = np.cumsum(0.5 * (outward_rate[:-1] + outward_rate[1:]))
    residual = mean_moisture + cumulative_loss / (RADIUS_M * RADIUS_M / 2.0) - initial_mean
    result = {
        "n_intervals": n,
        "dr_m": dr,
        "times_s": times,
        "radius_nodes_m": nodes,
        "radius_output_cm": nodes[sample_nodes] * 100.0,
        "temperature_C": temperature_history,
        "moisture_kg_per_kg": moisture_history,
        "mean_moisture": mean_moisture,
        "cumulative_loss": cumulative_loss,
        "mass_balance_residual": residual,
        "environment_temperature_C": env_temp,
        "environment_moisture_kg_per_kg": env_moisture,
        "final_time_s": float(end),
        "temperature_final_C": state[0::2].copy(),
        "moisture_final_kg_per_kg": state[1::2].copy(),
        "control_volume_weights_m2": volumes,
        "runtime_s": time.perf_counter() - started,
        "nfev": nfev,
        "njev": njev,
        "nlu": nlu,
        "positivity_retries": positivity_retries,
        "smallest_max_step_used_s": smallest_step_used,
        "state_min_temperature_C": float(temperature_history.min()),
        "state_max_temperature_C": float(temperature_history.max()),
        "state_min_moisture_kg_per_kg": float(moisture_history.min()),
        "state_max_moisture_kg_per_kg": float(moisture_history.max()),
        "initial_mean_moisture": initial_mean,
        "h_W_m2K": float(h),
        "hm_m_s": float(hm),
        "rtol": float(rtol),
        "atol": float(atol),
        "max_step_s": float(max_step),
        "boundary_start_s": float(boundary[0, 0]),
        "boundary_end_s": float(boundary[-1, 0]),
        "boundary_interpolation": "piecewise linear, no extrapolation",
        "model": (
            "Q2 fixed-radius axisymmetric variable-property coupled FV+BDF"
            if property_model == "q2"
            else "Q1 constitutive degradation in fixed-radius coupled FV+BDF"
        ),
        "moisture_parameterization": MOISTURE_PARAMETERIZATION,
        "property_model": property_model,
    }
    return result


def compare_results(first: dict[str, Any], second: dict[str, Any]) -> dict[str, Any]:
    """Compare temperature and moisture independently, with locations."""
    if not np.array_equal(first["times_s"], second["times_s"]):
        raise ValueError("Cannot compare results with different output times.")
    dt = np.abs(first["temperature_C"] - second["temperature_C"])
    dc = np.abs(first["moisture_kg_per_kg"] - second["moisture_kg_per_kg"])
    paper_mask = np.isin(first["times_s"], PAPER_SECONDS)
    paper_nodes = np.isin(first["radius_output_cm"], np.arange(5) * 0.5)
    def details(diff: np.ndarray, threshold: float) -> dict[str, Any]:
        all_index = np.unravel_index(int(np.argmax(diff)), diff.shape)
        paper_diff = diff[np.ix_(paper_mask, paper_nodes)]
        paper_index = np.unravel_index(int(np.argmax(paper_diff)), paper_diff.shape)
        paper_times = first["times_s"][paper_mask]
        paper_radii = first["radius_output_cm"][paper_nodes]
        return {
            "max_all": float(diff.max()),
            "threshold": float(threshold),
            "all_time_s": int(first["times_s"][all_index[0]]),
            "all_radius_cm": float(first["radius_output_cm"][all_index[1]]),
            "max_paper": float(paper_diff.max()),
            "paper_time_s": int(paper_times[paper_index[0]]),
            "paper_radius_cm": float(paper_radii[paper_index[1]]),
            "same_rounded_paper": bool(
                np.array_equal(
                    np.round(first["temperature_C" if diff is dt else "moisture_kg_per_kg"][np.ix_(paper_mask, paper_nodes)], 4),
                    np.round(second["temperature_C" if diff is dt else "moisture_kg_per_kg"][np.ix_(paper_mask, paper_nodes)], 4),
                )
            ),
            "pass_all": bool(diff.max() < threshold),
        }
    return {"temperature": details(dt, CONVERGENCE_THRESHOLD_T_C),
            "moisture": details(dc, CONVERGENCE_THRESHOLD_C)}


def write_csv(path: Path, header: list[str], values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(values.tolist())


def save_result_files(
    output: Path,
    result: dict[str, Any],
    boundary: np.ndarray,
    boundary_sheet: str,
    validation: dict[str, Any],
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    times = result["times_s"]
    radius_cm = result["radius_output_cm"]
    np.savez_compressed(
        output / "q2_fields.npz",
        time_s=times,
        radius_cm=radius_cm,
        temperature_C=result["temperature_C"],
        moisture_kg_per_kg=result["moisture_kg_per_kg"],
        mean_moisture=result["mean_moisture"],
        cumulative_loss=result["cumulative_loss"],
        mass_balance_residual=result["mass_balance_residual"],
    )
    np.savez_compressed(
        output / "q2_restart.npz",
        final_time_s=np.array(result["final_time_s"]),
        radius_nodes_m=result["radius_nodes_m"],
        control_volume_weights_m2=result["control_volume_weights_m2"],
        temperature_final_C=result["temperature_final_C"],
        moisture_final_kg_per_kg=result["moisture_final_kg_per_kg"],
        n_intervals=np.array(result["n_intervals"]),
        radius_m=np.array(RADIUS_M),
        length_m=np.array(LENGTH_M),
        h_W_m2K=np.array(result["h_W_m2K"]),
        hm_m_s=np.array(result["hm_m_s"]),
        rtol=np.array(result["rtol"]),
        atol=np.array(result["atol"]),
        max_step_s=np.array(result["max_step_s"]),
        boundary_time_s=boundary[:, 0],
        boundary_temperature_C=boundary[:, 1],
        boundary_moisture_kg_per_kg=boundary[:, 2],
        boundary_sheet=np.array(boundary_sheet),
        boundary_interpolation=np.array("piecewise linear, no extrapolation"),
        moisture_parameterization=np.array(result["moisture_parameterization"]),
        model=np.array(result["model"]),
    )
    radii_header = [f"r_{value:.1f}_cm" for value in radius_cm]
    write_csv(
        output / "q2_temperature_full.csv",
        ["time_s"] + radii_header,
        np.column_stack([times, result["temperature_C"]]),
    )
    write_csv(
        output / "q2_moisture_full.csv",
        ["time_s"] + radii_header,
        np.column_stack([times, result["moisture_kg_per_kg"]]),
    )
    write_csv(
        output / "boundary_used.csv",
        ["time_s", "temperature_C", "moisture_kg_per_kg"],
        np.column_stack([times, result["environment_temperature_C"], result["environment_moisture_kg_per_kg"]]),
    )
    paper_mask = np.isin(times, PAPER_SECONDS)
    paper_radii = radius_cm[PAPER_NODE_INDEX]
    paper_header = ["time_s"] + [f"r_{value:g}_cm" for value in paper_radii]
    write_csv(
        output / "q2_temperature_paper.csv",
        paper_header,
        np.column_stack([times[paper_mask], result["temperature_C"][np.ix_(paper_mask, PAPER_NODE_INDEX)]]),
    )
    write_csv(
        output / "q2_moisture_paper.csv",
        paper_header,
        np.column_stack([times[paper_mask], result["moisture_kg_per_kg"][np.ix_(paper_mask, PAPER_NODE_INDEX)]]),
    )
    write_csv(
        output / "q2_mass_balance.csv",
        ["time_s", "mean_moisture", "cumulative_loss_factor", "effective_balance_residual"],
        np.column_stack([times, result["mean_moisture"], result["cumulative_loss"], result["mass_balance_residual"]]),
    )
    (output / "q2_validation.json").write_text(
        json.dumps(validation, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (output / "q2_parameters_environment.json").write_text(
        json.dumps(
            {
                "parameters": validation.get("parameters", {}),
                "environment": validation.get("environment", {}),
                "boundary_source": validation.get("boundary_source"),
                "boundary_sheet": validation.get("boundary_sheet"),
                "boundary_interpolation": validation.get("boundary_interpolation"),
                "model": validation.get("model"),
                "moisture_parameterization": validation.get("moisture_parameterization"),
                "requested_end_s": validation.get("requested_end_s"),
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )


def _copy_cell_style(source, target) -> None:
    if source.has_style:
        target._style = copy_style(source._style)
    if source.number_format:
        target.number_format = source.number_format
    if source.alignment:
        target.alignment = copy_style(source.alignment)


def write_workbook(
    destination: str | Path,
    temperature_C: np.ndarray,
    moisture: np.ndarray,
    template: str | Path,
    end: int,
) -> None:
    """Write a 1..end workbook from the original two-sheet template."""
    destination = Path(destination)
    template = Path(template)
    temperature_C = np.asarray(temperature_C, dtype=float)
    moisture = np.asarray(moisture, dtype=float)
    if temperature_C.shape != (end + 1, OUTPUT_NODES) or moisture.shape != temperature_C.shape:
        raise ValueError("Workbook fields must contain t=0..end and 21 radial points.")
    wb = load_workbook(template)
    expected = {"温度", "水分浓度"}
    if set(wb.sheetnames) != expected:
        wb.close()
        raise ValueError(f"Expected sheets {expected}, got {wb.sheetnames}.")
    for sheet_name, field in (("温度", temperature_C), ("水分浓度", moisture)):
        ws = wb[sheet_name]
        if ws.cell(1, 1).value is None:
            wb.close()
            raise ValueError("Template A1 header is missing.")
        template_max_column = ws.max_column
        if template_max_column < 2:
            wb.close()
            raise ValueError("Template must contain at least one distance column.")
        header_style_row = list(ws[1])
        data_style_row = list(ws[2]) if ws.max_row >= 2 else list(ws[1])
        template_data_column = min(template_max_column, OUTPUT_NODES + 1)
        template_data_letter = get_column_letter(template_data_column)
        template_data_width = ws.column_dimensions[template_data_letter].width
        template_data_hidden = ws.column_dimensions[template_data_letter].hidden
        template_data_outline = ws.column_dimensions[template_data_letter].outlineLevel
        template_data_collapsed = ws.column_dimensions[template_data_letter].collapsed

        # Clear values only.  Styles and existing A:F widths remain those of
        # the original template; new G:V columns are extended from its final
        # real data column (F), not from a hard-coded width or default style.
        for row in ws.iter_rows(min_row=2, max_row=max(ws.max_row, end + 1), min_col=1, max_col=OUTPUT_NODES + 1):
            for cell in row:
                cell.value = None

        # Preserve the template's A-column style for every newly materialized
        # data row, including rows beyond the five-row demonstration table.
        source_a_data = data_style_row[0]
        source_row_height = ws.row_dimensions[2].height if ws.max_row >= 2 else None
        for row_index in range(2, end + 2):
            target_a = ws.cell(row_index, 1)
            _copy_cell_style(source_a_data, target_a)
            if source_row_height is not None:
                ws.row_dimensions[row_index].height = source_row_height

        for j in range(OUTPUT_NODES):
            output_column = j + 2
            # B:D are real template distance columns.  E is the template's
            # ellipsis marker, so all actual E:V distance columns inherit the
            # final real distance style F.
            source_column = output_column if output_column <= 4 else template_data_column
            source_header = header_style_row[source_column - 1]
            source_data = data_style_row[source_column - 1]
            header = ws.cell(1, output_column)
            _copy_cell_style(source_header, header)
            header.value = j / 10.0
            for row_index, second in enumerate(range(1, end + 1), start=2):
                cell = ws.cell(row_index, output_column)
                _copy_cell_style(source_data, cell)
                cell.value = float(f"{field[second, j]:.4f}")
                # The template demonstration cells are General, while the
                # formal result contract requires four-decimal data values.
                # This format is applied uniformly to B:V, including G:V.
                cell.number_format = "0.0000"

            if output_column > template_max_column:
                column_letter = get_column_letter(output_column)
                dimension = ws.column_dimensions[column_letter]
                dimension.width = template_data_width
                dimension.hidden = template_data_hidden
                dimension.outlineLevel = template_data_outline
                dimension.collapsed = template_data_collapsed
        for row_index, second in enumerate(range(1, end + 1), start=2):
            ws.cell(row_index, 1).value = second
        ws.freeze_panes = "B2"
    wb.save(destination)
    wb.close()


def validate_workbook(
    path: str | Path,
    temperature_C: np.ndarray,
    moisture: np.ndarray,
    end: int,
    template: str | Path,
) -> dict[str, Any]:
    """Read every cell back and compare with the four-decimal source arrays."""
    path = Path(path)
    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        if set(wb.sheetnames) != {"温度", "水分浓度"}:
            raise AssertionError(f"Unexpected sheet names: {wb.sheetnames}")
        template_wb = load_workbook(template, read_only=True, data_only=True)
        try:
            expected_a1 = template_wb["温度"].cell(1, 1).value
        finally:
            template_wb.close()
        report: dict[str, Any] = {"path": str(path.resolve()), "sheets": {}}
        for sheet_name, field in (("温度", temperature_C), ("水分浓度", moisture)):
            ws = wb[sheet_name]
            if ws.max_row != end + 1 or ws.max_column != OUTPUT_NODES + 1:
                raise AssertionError(
                    f"{sheet_name} dimensions {ws.max_row}x{ws.max_column}, "
                    f"expected {end + 1}x{OUTPUT_NODES + 1}."
                )
            if ws.cell(1, 1).value != expected_a1:
                raise AssertionError("A1 was not preserved from the original template.")
            headers = [ws.cell(1, j + 2).value for j in range(OUTPUT_NODES)]
            np.testing.assert_allclose(np.asarray(headers, dtype=float), np.arange(21) / 10.0, atol=1e-12, rtol=0.0)
            expected = np.round(field[1:], 4)
            observed = np.empty_like(expected)
            for row_index, row in enumerate(ws.iter_rows(min_row=2, max_row=end + 1, max_col=OUTPUT_NODES + 1), start=0):
                if row[0].value != row_index + 1:
                    raise AssertionError(f"{sheet_name} time column mismatch at row {row_index + 2}.")
                observed[row_index] = np.asarray([cell.value for cell in row[1:]], dtype=float)
            np.testing.assert_array_equal(observed, expected)
            report["sheets"][sheet_name] = {
                "rows": int(ws.max_row - 1),
                "columns": int(ws.max_column - 1),
                "four_decimal_match": True,
            }
        return report
    finally:
        wb.close()


def validate_workbook_template_styles(
    path: str | Path,
    template: str | Path,
    end: int,
) -> dict[str, Any]:
    """Verify template styles and widths after extending the 21-point table.

    Existing A:F column widths are required to remain byte-for-value equal to
    the source template.  New G:V columns must inherit the final real template
    distance column (F) width and cell style.  Header cells use the template
    header style; data cells use the template data style plus the explicit
    four-decimal output format required by the result contract.
    """
    path = Path(path)
    template = Path(template)
    wb = load_workbook(path, read_only=False, data_only=False)
    source_wb = load_workbook(template, read_only=False, data_only=False)

    def style_components(cell) -> tuple[Any, ...]:
        # StyleProxy equality is workbook-local; copy the underlying style
        # records before comparing a generated workbook with its source.
        return tuple(copy_style(getattr(cell, name)) for name in ("font", "alignment", "fill", "border", "protection"))

    try:
        expected = {"温度", "水分浓度"}
        if set(wb.sheetnames) != expected or set(source_wb.sheetnames) != expected:
            raise AssertionError("Workbook/template sheet names do not match.")
        report: dict[str, Any] = {"path": str(path.resolve()), "template": str(template.resolve()), "sheets": {}}
        for sheet_name in ("温度", "水分浓度"):
            ws = wb[sheet_name]
            source = source_wb[sheet_name]
            if ws.max_row != end + 1 or ws.max_column != OUTPUT_NODES + 1:
                raise AssertionError(f"{sheet_name} dimensions are not the formal output dimensions.")
            template_max_column = source.max_column
            source_distance_column = min(template_max_column, OUTPUT_NODES + 1)
            for col in range(1, template_max_column + 1):
                letter = get_column_letter(col)
                if ws.column_dimensions[letter].width != source.column_dimensions[letter].width:
                    raise AssertionError(f"{sheet_name} existing column width changed at {letter}.")
            source_width = source.column_dimensions[get_column_letter(source_distance_column)].width
            for col in range(template_max_column + 1, OUTPUT_NODES + 2):
                letter = get_column_letter(col)
                if ws.column_dimensions[letter].width != source_width:
                    raise AssertionError(f"{sheet_name} new column width mismatch at {letter}.")

            # Check representative rows: header, first data row, and final
            # data row.  All output columns are checked, not only G:V.
            for output_col in range(2, OUTPUT_NODES + 2):
                source_col = output_col if output_col <= 4 else source_distance_column
                if style_components(ws.cell(1, output_col)) != style_components(source.cell(1, source_col)):
                    raise AssertionError(f"{sheet_name} header style mismatch at {get_column_letter(output_col)}.")
                if ws.cell(1, output_col).number_format != source.cell(1, source_col).number_format:
                    raise AssertionError(f"{sheet_name} header number format mismatch at {get_column_letter(output_col)}.")
                for row in (2, end + 1):
                    if style_components(ws.cell(row, output_col)) != style_components(source.cell(2, source_col)):
                        raise AssertionError(
                            f"{sheet_name} data style mismatch at {get_column_letter(output_col)}{row}."
                        )
                    if ws.cell(row, output_col).number_format != "0.0000":
                        raise AssertionError(
                            f"{sheet_name} data number format mismatch at {get_column_letter(output_col)}{row}."
                        )
            report["sheets"][sheet_name] = {
                "template_style_match": True,
                "existing_widths_preserved": True,
                "new_columns_width_source": get_column_letter(source_distance_column),
                "new_columns": "G:V",
            }
        return report
    finally:
        source_wb.close()
        wb.close()


def plot_results(output: Path, result: dict[str, Any]) -> None:
    """Create the Q2 four-panel figure directly from the saved result arrays."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    available = {font.name for font in font_manager.fontManager.ttflist}
    font_name = next(
        (name for name in ("Microsoft YaHei", "SimHei", "SimSun", "Noto Sans CJK SC") if name in available),
        None,
    )
    if font_name:
        plt.rcParams["font.family"] = font_name
    plt.rcParams["font.size"] = 9
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False
    times = result["times_s"]
    radius_cm = result["radius_output_cm"]
    minutes = times / 60.0
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4), constrained_layout=True)
    colors = ["#174A70", "#2584A6", "#77AABD", "#D88A41", "#963E3E", "#6A4C93"]
    profile_seconds = PROFILE_SECONDS.tolist()
    for second, color in zip(profile_seconds, colors):
        idx = int(np.where(times == second)[0][0])
        axes[0, 0].plot(radius_cm, result["temperature_C"][idx], color=color, label=f"{second / 3600:g} h")
        axes[1, 0].plot(radius_cm, result["moisture_kg_per_kg"][idx], color=color, label=f"{second / 3600:g} h")
    for index, color in zip(PAPER_NODE_INDEX, colors[:5]):
        axes[0, 1].plot(minutes, result["temperature_C"][:, index], color=color, label=f"{index / 10:g} cm")
        axes[1, 1].plot(minutes, result["moisture_kg_per_kg"][:, index], color=color, label=f"{index / 10:g} cm")
    axes[0, 0].set(xlabel="径向距离（cm）", ylabel="温度（°C）")
    axes[1, 0].set(xlabel="径向距离（cm）", ylabel="干基含水率（kg/kg）")
    axes[0, 1].set(xlabel="时间（min）", ylabel="温度（°C）")
    axes[1, 1].set(xlabel="时间（min）", ylabel="干基含水率（kg/kg）")
    for ax in axes.ravel():
        ax.grid(alpha=0.18)
        ax.legend(frameon=False, fontsize=7, ncol=2)
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / "q2_fields.pdf", bbox_inches="tight")
    fig.savefig(output / "q2_fields.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def sensitivity_runs(
    boundary: np.ndarray,
    result: dict[str, Any],
    output: Path,
    end: int,
) -> list[dict[str, Any]]:
    """Run the required one-factor +/-20 percent h and hm sensitivity cases."""
    cases = [
        ("baseline", result["h_W_m2K"], result["hm_m_s"], result),
        ("h_minus20", 0.8 * H_BASE, HM_BASE, None),
        ("h_plus20", 1.2 * H_BASE, HM_BASE, None),
        ("hm_minus20", H_BASE, 0.8 * HM_BASE, None),
        ("hm_plus20", H_BASE, 1.2 * HM_BASE, None),
    ]
    rows: list[dict[str, Any]] = []
    baseline_values = {
        "T_center_C": float(result["temperature_C"][-1, 0]),
        "T_surface_C": float(result["temperature_C"][-1, -1]),
        "C_center": float(result["moisture_kg_per_kg"][-1, 0]),
        "C_surface": float(result["moisture_kg_per_kg"][-1, -1]),
    }
    for name, h, hm, existing in cases:
        current = existing or integrate_case(
            boundary,
            n=result["n_intervals"],
            end=end,
            rtol=result["rtol"],
            atol=result["atol"],
            max_step=result["max_step_s"],
            h=h,
            hm=hm,
            chunk_seconds=DEFAULT_CHUNK_S,
        )
        values = {
            "case": name,
            "h_W_m2K": float(h),
            "hm_m_s": float(hm),
            "T_center_C": float(current["temperature_C"][-1, 0]),
            "T_surface_C": float(current["temperature_C"][-1, -1]),
            "C_center_kg_per_kg": float(current["moisture_kg_per_kg"][-1, 0]),
            "C_surface_kg_per_kg": float(current["moisture_kg_per_kg"][-1, -1]),
            "runtime_s": float(current["runtime_s"]),
        }
        values.update(
            {
                "delta_T_center_C": values["T_center_C"] - baseline_values["T_center_C"],
                "delta_T_surface_C": values["T_surface_C"] - baseline_values["T_surface_C"],
                "delta_C_center_kg_per_kg": values["C_center_kg_per_kg"] - baseline_values["C_center"],
                "delta_C_surface_kg_per_kg": values["C_surface_kg_per_kg"] - baseline_values["C_surface"],
            }
        )
        rows.append(values)
    headers = list(rows[0].keys())
    output.mkdir(parents=True, exist_ok=True)
    with (output / "q2_sensitivity.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def run_convergence(
    boundary: np.ndarray,
    output: Path,
    end: int,
    max_n: int,
    rtol: float,
    atol: float,
    max_step: float,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]], dict[str, Any]]:
    """Run space refinement and a tenfold time refinement."""
    if max_n < 20 or max_n % 20 != 0:
        raise ValueError("max_n must be a multiple of 20.")
    records: list[dict[str, Any]] = []
    diagnostics_dir = output / "diagnostics"
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    previous: dict[str, Any] | None = None
    final: dict[str, Any] | None = None
    status = "REFINEMENT_REQUIRED"
    n = 20
    while n <= max_n:
        started = time.perf_counter()
        try:
            current = integrate_case(
                boundary,
                n=n,
                end=end,
                rtol=rtol,
                atol=atol,
                max_step=max_step,
            )
        except Exception as exc:
            record = {"n_intervals": n, "status": "SOLVER_FAILED", "error": repr(exc)}
            records.append(record)
            (diagnostics_dir / f"q2_n{n}_failure.json").write_text(
                json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
            )
            break
        record: dict[str, Any] = {
            "n_intervals": n,
            "dr_m": float(current["dr_m"]),
            "runtime_s": float(time.perf_counter() - started),
            "integrator": {
                "nfev": int(current["nfev"]),
                "njev": int(current["njev"]),
                "nlu": int(current["nlu"]),
            },
        }
        if previous is not None:
            differences = compare_results(previous, current)
            record["difference_from_previous"] = differences
            record["temperature_threshold_pass"] = differences["temperature"]["pass_all"]
            record["moisture_threshold_pass"] = differences["moisture"]["pass_all"]
            record["paper_rounding_pass"] = (
                differences["temperature"]["same_rounded_paper"]
                and differences["moisture"]["same_rounded_paper"]
            )
        records.append(record)
        (diagnostics_dir / f"q2_n{n}.json").write_text(
            json.dumps(_json_safe(record), ensure_ascii=False, indent=2), encoding="utf-8"
        )
        # Persist the successful comparison before leaving the refinement
        # loop.  The old ordering broke here first, silently dropping the
        # passing N record and attaching time refinement to the prior grid.
        if previous is not None and (
            record["temperature_threshold_pass"]
            and record["moisture_threshold_pass"]
            and record["paper_rounding_pass"]
        ):
            status = "SPATIAL_CHECK_PASSED"
            final = current
            break
        previous = current
        final = current
        n *= 2
    if final is None:
        return None, records, {"status": status, "reason": "No successful grid run."}
    if status != "SPATIAL_CHECK_PASSED":
        return final, records, {
            "status": status,
            "reason": f"Did not pass separate e_T/e_C thresholds by n={max_n}.",
        }
    try:
        tight = integrate_case(
            boundary,
            n=final["n_intervals"],
            end=end,
            rtol=rtol / 10.0,
            atol=atol / 10.0,
            max_step=max_step / 2.0,
        )
        temporal_difference = compare_results(final, tight)
        time_pass = (
            temporal_difference["temperature"]["pass_all"]
            and temporal_difference["moisture"]["pass_all"]
            and temporal_difference["temperature"]["same_rounded_paper"]
            and temporal_difference["moisture"]["same_rounded_paper"]
        )
        records[-1]["time_refinement"] = temporal_difference
        records[-1]["time_refinement_pass"] = time_pass
        if time_pass:
            final = tight
            status = "NUMERICAL_CHECKS_PASSED"
        else:
            status = "REFINEMENT_REQUIRED"
    except Exception as exc:
        records[-1]["time_refinement_error"] = repr(exc)
        status = "REFINEMENT_REQUIRED"
    validation = {
        "status": status,
        "spatial_records": records,
        "final_n_intervals": int(final["n_intervals"]),
        "separate_thresholds": {
            "temperature_C": CONVERGENCE_THRESHOLD_T_C,
            "moisture_kg_per_kg": CONVERGENCE_THRESHOLD_C,
        },
        "time_refinement": records[-1].get("time_refinement"),
        "accuracy_statement": "Empirical grid/time differences, not rigorous error bounds.",
    }
    return final, records, validation


def prepare_template(template: Path) -> Path:
    """Preserve the original result2 template before any generated workbook exists."""
    if not template.exists():
        raise FileNotFoundError(f"Template file does not exist: {template}")
    TEMPLATE_BACKUP.parent.mkdir(parents=True, exist_ok=True)
    if not TEMPLATE_BACKUP.exists():
        shutil.copy2(template, TEMPLATE_BACKUP)
    return TEMPLATE_BACKUP


def build_validation(
    final: dict[str, Any],
    boundary: np.ndarray,
    boundary_sheet: str,
    convergence: dict[str, Any],
    workbook_report: dict[str, Any],
    sensitivity: list[dict[str, Any]],
) -> dict[str, Any]:
    result = dict(convergence)
    result.update(
        {
            "status": "NUMERICAL_CHECKS_PASSED",
            "boundary_source": str(DEFAULT_BOUNDARY.resolve()),
            "boundary_sheet": boundary_sheet,
            "boundary_samples": int(len(boundary)),
            "boundary_time_range_s": [float(boundary[0, 0]), float(boundary[-1, 0])],
            "requested_end_s": int(final["final_time_s"]),
            "formal_workbook_time_range_s": [1, int(final["final_time_s"])],
            "formal_workbook_interpretation": (
                "User-approved Q2 scope: official result2.xlsx contains t=1..10800 s. "
                "The wording 'entire drying process' may require extension after Q3 determines "
                "the endpoint and post-14400 s boundary."
            ),
            "parameters": {
                "radius_m": RADIUS_M,
                "length_m": LENGTH_M,
                "initial_temperature_C": INITIAL_T_C,
                "initial_moisture_kg_per_kg": INITIAL_C,
                "h_W_m2K": final["h_W_m2K"],
                "hm_m_s": final["hm_m_s"],
                "h_assumption": "Appendix 2 is labelled Q1; Q2 has no new h, so this is a missing-parameter continuation assumption.",
                "hm_assumption": "Appendix 2 is labelled Q1; Q2 has no new hm, so this is a missing-parameter continuation assumption.",
            },
            "model": final["model"],
            "moisture_parameterization": final["moisture_parameterization"],
            "historical_diagnostics": {
                "obsolete": True,
                "directory": "diagnostics/archive/historical_raw_c",
                "note": "Raw-C BDF failure evidence is historical only; current status is read from q2_validation.json and q2_n*.json.",
            },
            "boundary_interpolation": "piecewise linear, no extrapolation",
            "restart_file": "q2_restart.npz",
            "workbook_validation": workbook_report,
            "sensitivity": sensitivity,
            "environment": {
                "python": sys.version,
                "executable": sys.executable,
                "numpy": np.__version__,
                "scipy": __import__("scipy").__version__,
                "openpyxl": __import__("openpyxl").__version__,
                "platform": platform.platform(),
            },
        }
    )
    result["boundary_source"] = str(DEFAULT_BOUNDARY.resolve())
    return _json_safe(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundary", type=Path, default=DEFAULT_BOUNDARY)
    parser.add_argument("--sheet", default=None)
    parser.add_argument("--end", type=int, default=DEFAULT_END_S)
    parser.add_argument("--max-n", type=int, default=10240)
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-11)
    parser.add_argument("--max-step", type=float, default=5.0)
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-sensitivity", action="store_true")
    args = parser.parse_args()

    output = args.output.resolve()
    boundary_path = args.boundary.resolve()
    template_path = args.template.resolve()
    boundary, boundary_sheet = read_boundary(boundary_path, sheet=args.sheet, end=args.end)
    template_backup = prepare_template(template_path)
    final, records, convergence = run_convergence(
        boundary,
        output=output,
        end=args.end,
        max_n=args.max_n,
        rtol=args.rtol,
        atol=args.atol,
        max_step=args.max_step,
    )
    if final is None or convergence.get("status") != "NUMERICAL_CHECKS_PASSED":
        failure = {
            "status": convergence.get("status", "REFINEMENT_REQUIRED"),
            "boundary_source": str(boundary_path),
            "template_backup": str(template_backup),
            "requested_end_s": args.end,
            "convergence": convergence,
            "records": records,
        }
        (output / "diagnostics").mkdir(parents=True, exist_ok=True)
        (output / "diagnostics" / "q2_failure_summary.json").write_text(
            json.dumps(_json_safe(failure), ensure_ascii=False, indent=2), encoding="utf-8"
        )
        print(json.dumps(_json_safe(failure), ensure_ascii=False))
        raise SystemExit(2)

    candidate = output / "result2_candidate_3h.xlsx"
    write_workbook(
        candidate,
        final["temperature_C"],
        final["moisture_kg_per_kg"],
        template=template_backup,
        end=args.end,
    )
    workbook_report = validate_workbook(
        candidate,
        final["temperature_C"],
        final["moisture_kg_per_kg"],
        end=args.end,
        template=template_backup,
    )
    workbook_style_report = validate_workbook_template_styles(
        candidate,
        template_backup,
        end=args.end,
    )
    sensitivity = [] if args.skip_sensitivity else sensitivity_runs(boundary, final, output, args.end)
    validation = build_validation(
        final,
        boundary,
        boundary_sheet,
        convergence,
        workbook_report,
        sensitivity,
    )
    validation["workbook_style_validation"] = workbook_style_report
    save_result_files(output, final, boundary, boundary_sheet, validation)
    plot_results(output, final)
    official = OFFICIAL_RESULT2.resolve()
    shutil.copy2(candidate, official)
    official_report = validate_workbook(
        official,
        final["temperature_C"],
        final["moisture_kg_per_kg"],
        end=args.end,
        template=template_backup,
    )
    official_style_report = validate_workbook_template_styles(
        official,
        template_backup,
        end=args.end,
    )
    validation["official_workbook"] = str(official)
    validation["official_workbook_validation"] = official_report
    validation["official_workbook_style_validation"] = official_style_report
    (output / "q2_validation.json").write_text(
        json.dumps(_json_safe(validation), ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    summary = {
        "status": validation["status"],
        "boundary": str(boundary_path),
        "template_backup": str(template_backup),
        "candidate": str(candidate.resolve()),
        "official": str(official),
        "n_intervals": final["n_intervals"],
        "runtime_s": final["runtime_s"],
        "T_3h_center_surface_C": final["temperature_C"][-1, [0, -1]].tolist(),
        "C_3h_center_surface": final["moisture_kg_per_kg"][-1, [0, -1]].tolist(),
        "max_effective_moisture_balance_residual": float(np.max(np.abs(final["mass_balance_residual"]))),
    }
    print(json.dumps(_json_safe(summary), ensure_ascii=False))


if __name__ == "__main__":
    main()
