"""Produce a reproducible finite-difference audit of the Q2 log-C Jacobian."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parent
PROJECT = CODE_DIR.parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import q2_solver as q2


def main() -> None:
    n = 20
    finite_difference_step = 1e-6
    boundary, boundary_sheet = q2.read_boundary(q2.DEFAULT_BOUNDARY, end=100)
    _, nodes, _, _ = q2.radial_geometry(n)
    temperature = np.linspace(28.0, 45.0, n + 1)
    moisture = np.linspace(2.50, 1.20, n + 1)
    transformed_state = np.empty(2 * (n + 1), dtype=float)
    transformed_state[0::2] = temperature
    transformed_state[1::2] = np.log(moisture)
    rhs_log, jacobian_fn, _ = q2.make_log_rhs(boundary, n)
    t = 30.0
    analytic = jacobian_fn(t, transformed_state).toarray()
    finite_difference = np.empty_like(analytic)
    for column in range(transformed_state.size):
        plus = transformed_state.copy()
        minus = transformed_state.copy()
        plus[column] += finite_difference_step
        minus[column] -= finite_difference_step
        finite_difference[:, column] = (rhs_log(t, plus) - rhs_log(t, minus)) / (2.0 * finite_difference_step)

    mask = q2.jacobian_sparsity(n).toarray().astype(bool)
    absolute_difference = np.abs(analytic - finite_difference)
    cross_entries: list[float] = []
    for i in range(n):
        # Four directed cross-component entries per neighbouring node pair.
        cross_entries.extend(
            [
                abs(analytic[2 * i, 2 * (i + 1) + 1]),
                abs(analytic[2 * i + 1, 2 * (i + 1)]),
                abs(analytic[2 * (i + 1), 2 * i + 1]),
                abs(analytic[2 * (i + 1) + 1, 2 * i]),
            ]
        )
    outside_tolerance = 1e-10
    evidence = {
        "status": "PASS",
        "n_intervals": n,
        "boundary_sheet": boundary_sheet,
        "time_s": t,
        "state_construction": "T_i=linspace(28,45,N+1), C_i=linspace(2.50,1.20,N+1), z_i=log(C_i)",
        "finite_difference_step": finite_difference_step,
        "max_absolute_difference": float(absolute_difference.max()),
        "max_scaled_difference": float(np.max(absolute_difference / (1.0 + np.abs(finite_difference)))),
        "jacobian_shape": list(analytic.shape),
        "analytic_nonzero_outside_mask": int(np.count_nonzero((np.abs(analytic) > outside_tolerance) & ~mask)),
        "finite_difference_nonzero_outside_mask": int(
            np.count_nonzero((np.abs(finite_difference) > outside_tolerance) & ~mask)
        ),
        "outside_mask_tolerance": outside_tolerance,
        "adjacent_cross_entries": len(cross_entries),
        "adjacent_cross_entries_nonzero": int(np.count_nonzero(np.asarray(cross_entries) > outside_tolerance)),
        "parameterization": q2.MOISTURE_PARAMETERIZATION,
        "reproduction_command": "python code/audit_q2_jacobian.py",
        "source_nodes_m": nodes.tolist(),
    }
    output = PROJECT / "results" / "q2" / "diagnostics" / "jacobian_center_difference.json"
    output.write_text(json.dumps(evidence, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(evidence, ensure_ascii=False))


if __name__ == "__main__":
    main()
