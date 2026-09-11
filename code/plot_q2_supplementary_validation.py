"""Plot-only companion for ``q2_supplementary_validation.py``.

This script reads supplementary CSV/NPZ/JSON artifacts and never invokes a
solver.  It produces candidate figures under ``results/q2_verification``;
the computation and plotting stages therefore remain independently rerunnable.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "q2_verification"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 10,
            "font.sans-serif": [
                "Microsoft YaHei",
                "SimHei",
                "Arial Unicode MS",
                "DejaVu Sans",
            ],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
        }
    )
    return plt


def plot_degradation(output: Path) -> list[str]:
    rows = _read_csv(output / "q2_degradation.csv")
    grid_rows = [row for row in rows if row.get("case", "").startswith("N=")]
    if not grid_rows:
        raise ValueError("No N-grid degradation rows found.")
    n = np.array([int(row["n_intervals"]) for row in grid_rows], dtype=float)
    t_error = np.array([float(row["temperature_max_abs_C"]) for row in grid_rows])
    c_error = np.array([float(row["moisture_max_abs_kg_per_kg"]) for row in grid_rows])
    plt = _matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), layout="constrained")
    axes[0].loglog(n, t_error, "o-", color="#1f5a8a", label=r"$T$")
    axes[0].axhline(5e-5, color="#999999", ls="--", lw=0.9, label="threshold")
    axes[0].set(xlabel="N", ylabel="max |Q1-Q2 degraded| (°C)")
    axes[1].loglog(n, c_error, "o-", color="#a35335", label=r"$C$")
    axes[1].axhline(5e-5, color="#999999", ls="--", lw=0.9, label="threshold")
    axes[1].set(xlabel="N", ylabel="max |Q1-Q2 degraded| (kg/kg)")
    for axis in axes:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False)
    paths = ["q2_degradation_comparison.png", "q2_degradation_comparison.pdf"]
    fig.savefig(output / paths[0], dpi=180, bbox_inches="tight")
    fig.savefig(output / paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def plot_balance(output: Path) -> list[str]:
    npz_path = output / "q2_balance.npz"
    refined_npz_path = output / "q2_balance_time_refined.npz"
    json_path = output / "q2_balance.json"
    with np.load(npz_path, allow_pickle=False) as payload:
        with np.load(refined_npz_path, allow_pickle=False) as refined_payload:
            time_s = np.asarray(payload["time_s"], dtype=float)
            with json_path.open("r", encoding="utf-8") as handle:
                metadata: dict[str, Any] = json.load(handle)
            selected = metadata["settings"].get("gauss_legendre_selected_orders") or [4, 8]
            selected_order = int(selected[1])
            refined_selected = (
                metadata.get("refined_quadrature_selection", {}).get("selected_orders")
                or selected
            )
            refined_order = int(refined_selected[1])
            baseline_trapezoid = np.abs(np.asarray(payload["B_trapezoid_signed"], dtype=float))
            baseline_gl = np.abs(
                np.asarray(payload[f"B_gl{selected_order}_signed"], dtype=float)
            )
            refined_gl = np.abs(
                np.asarray(refined_payload[f"refined_B_gl{refined_order}_signed"], dtype=float)
            )
    plt = _matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9), layout="constrained")
    for axis, mask in (
        (axes[0], time_s > 0.0),
        (axes[1], (time_s > 0.0) & (time_s <= 60.0)),
    ):
        axis.plot(
            time_s[mask] / 3600.0 if axis is axes[0] else time_s[mask],
            baseline_trapezoid[mask],
            color="#555555",
            lw=0.9,
            ls="--",
            label="逐秒梯形",
        )
        axis.plot(
            time_s[mask] / 3600.0 if axis is axes[0] else time_s[mask],
            baseline_gl[mask],
            color="#1f5a8a",
            lw=1.1,
            label=f"基准 GL{selected_order}",
        )
        axis.plot(
            time_s[mask] / 3600.0 if axis is axes[0] else time_s[mask],
            refined_gl[mask],
            color="#a35335",
            lw=1.0,
            label=f"时间加密 GL{refined_order}",
        )
        axis.set_yscale("log")
        axis.set_ylabel(r"积分平衡残差 $E_{\mathrm{bal}}=|B|$/(kg/kg)")
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=8)
    axes[0].set_xlabel("时间/h")
    axes[1].set_xlabel("时间/s")
    axes[1].set_xlim(0.0, 60.0)
    axes[0].set_title("全时段")
    axes[1].set_title("初始60s")
    paths = ["q2_balance_residuals.png", "q2_balance_residuals.pdf"]
    fig.savefig(output / paths[0], dpi=180, bbox_inches="tight")
    fig.savefig(output / paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("degradation", "balance", "all"))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    paths: list[str] = []
    if args.mode in ("degradation", "all"):
        paths.extend(plot_degradation(output))
    if args.mode in ("balance", "all"):
        paths.extend(plot_balance(output))
    print(json.dumps({"status": "PLOTTED", "files": paths}, ensure_ascii=False))


if __name__ == "__main__":
    main()
