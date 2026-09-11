"""Q3 supplementary endpoint, boundary, and diffusivity experiments.

This entry point writes only to ``results/q3/supplement``.  It never invokes
the publication path, never writes the official workbook, and does not modify
the accepted Q3 result files.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import time
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

import q2_solver as q2
import q3_solver as q3


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "q3" / "supplement"
BASE = {"rtol": 1.0e-9, "atol": 1.0e-11,
        "max_step_before_s": 5.0, "max_step_after_s": 60.0}
TIGHT = {"rtol": 1.0e-10, "atol": 1.0e-12,
         "max_step_before_s": 2.5, "max_step_after_s": 30.0}
MAX_END_S = 168 * 3600
CASES = {
    "nominal": (50.0, 0.050),
    "temperature_49": (49.0, 0.050),
    "temperature_51": (51.0, 0.050),
    "moisture_0045": (50.0, 0.045),
    "moisture_0055": (50.0, 0.055),
}


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(json_safe(value), ensure_ascii=False, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def merge_prefix_and_continuation(
    first: dict[str, Any], second: dict[str, Any], plateau_boundary: np.ndarray,
    h: float = q2.H_BASE, hm: float = q2.HM_BASE,
) -> dict[str, Any]:
    """Mirror q3.run_case's exact split merge while reusing one 0--4 h run."""
    mask = second["times_s"] > q3.PLATEAU_START_S
    merged = dict(second)
    merged["times_s"] = np.concatenate([first["times_s"], second["times_s"][mask]])
    merged["temperature_C"] = np.vstack([first["temperature_C"], second["temperature_C"][mask]])
    merged["moisture_kg_per_kg"] = np.vstack([first["moisture_kg_per_kg"], second["moisture_kg_per_kg"][mask]])
    merged["audit"] = np.vstack([first["audit"], second["audit"][mask]])
    merged["flux_times_s"] = np.concatenate([first["flux_times_s"], second["flux_times_s"]])
    merged["flux_rates_kg_m_s"] = np.concatenate([first["flux_rates_kg_m_s"], second["flux_rates_kg_m_s"]])
    merged["initial_mean_moisture"] = first["initial_mean_moisture"]
    fm, sm = first["monitor"], second["monitor"]
    merged["monitor"] = {
        "max_positive_neighbor_jump_kg_per_kg": max(fm["max_positive_neighbor_jump_kg_per_kg"], sm["max_positive_neighbor_jump_kg_per_kg"]),
        "max_center_gap_kg_per_kg": max(fm["max_center_gap_kg_per_kg"], sm["max_center_gap_kg_per_kg"]),
        "min_moisture_kg_per_kg": min(fm["min_moisture_kg_per_kg"], sm["min_moisture_kg_per_kg"]),
        "max_moisture_kg_per_kg": max(fm["max_moisture_kg_per_kg"], sm["max_moisture_kg_per_kg"]),
        "checked_accepted_states": fm["checked_accepted_states"] + sm["checked_accepted_states"],
        "checked_midpoint_states": fm["checked_midpoint_states"] + sm["checked_midpoint_states"],
        "checked_output_states": fm["checked_output_states"] + sm["checked_output_states"],
        "segment_records": fm["segment_records"] + sm["segment_records"],
    }
    merged["event"] = second["event"]
    merged["first_segment"] = first
    merged["second_segment"] = second
    merged["final_time_s"] = second["final_time_s"]
    merged["state_final"] = second["state_final"]
    merged["runtime_s"] = float(first["runtime_s"] + second["runtime_s"])
    merged["nfev"] = int(first["nfev"] + second["nfev"])
    merged["njev"] = int(first["njev"] + second["njev"])
    merged["nlu"] = int(first["nlu"] + second["nlu"])
    merged["positivity_retries"] = int(first["positivity_retries"] + second["positivity_retries"])
    merged["boundary_plateau"] = plateau_boundary
    merged["h_W_m2K"] = float(h)
    merged["hm_m_s"] = float(hm)
    return merged


def continue_case(
    first: dict[str, Any], n: int, temperature_C: float, moisture: float,
    settings: dict[str, float], diagnostic: dict[str, Any],
) -> tuple[dict[str, Any], float]:
    plateau = q3.make_plateau_boundary(
        diagnostic, MAX_END_S, temperature_C=temperature_C,
        moisture_kg_per_kg=moisture,
    )
    transmitted = np.asarray(first["state_final"], dtype=float).copy()
    continuity = float(np.max(np.abs(transmitted - first["state_final"])))
    second = q3.integrate_segment(
        transmitted, n, q3.PLATEAU_START_S, MAX_END_S, plateau,
        rtol=settings["rtol"], atol=settings["atol"],
        max_step=settings["max_step_after_s"], detect_event=True,
        event_kind="grid_max", h=q2.H_BASE, hm=q2.HM_BASE,
    )
    return merge_prefix_and_continuation(first, second, plateau), continuity


def summarize_case(
    name: str, n: int, temperature_C: float, moisture: float,
    settings_name: str, settings: dict[str, float], result: dict[str, Any],
    boundary: np.ndarray, q2_history: dict[str, np.ndarray], continuity: float,
) -> dict[str, Any]:
    event = q3._event_metadata(result)
    # ``_event_metadata`` is shared with the formal publication pipeline and
    # names its persisted event-state file.  Supplement cases intentionally
    # persist scalar audit evidence only, so do not retain that formal path.
    event.pop("state_file", None)
    event["state_storage"] = "not persisted; scalar event evidence is embedded in this case record"
    fields = q3.build_output_arrays(result, q2_history, boundary)
    residual = np.asarray(fields["mass_balance_residual"], dtype=float)
    relative_balance = float(np.max(np.abs(residual)) / result["initial_mean_moisture"])
    monitor = result["monitor"]
    located = event.get("status") == "EVENT_LOCATED"
    checks = {
        "event_located": located,
        "event_residual": located and abs(event["event_residual_C"]) <= q3.EVENT_CONCENTRATION_TOL,
        "event_bracket": located and event["grid_event_bracket_width_s"] <= q3.EVENT_TIME_TOL_S,
        "strict_before_after": located and event["strict_before_pass"] and event["strict_after_pass"],
        "center_grid_agreement": located and event["center_event_time_s"] is not None
            and abs(event["center_event_time_s"] - event["time_s"]) <= q3.EVENT_TIME_TOL_S,
        "center_controls": located and abs(event["center_minus_max_C"]) <= q3.MONOTONICITY_TOL_C,
        "radial_monotonicity": monitor["max_positive_neighbor_jump_kg_per_kg"] <= q3.MONOTONICITY_TOL_C,
        "positive_moisture": monitor["min_moisture_kg_per_kg"] > 0.0,
        "state_continuity": continuity == 0.0,
        "accepted_node_mass_balance": fields["mass_balance_method"] == "accepted_BDF_nodes_with_piecewise_boundary_flux"
            and relative_balance < q3.MASS_BALANCE_REL_TOL,
    }
    return {
        "case": name, "status": "PASSED" if all(checks.values()) else "REVIEW_REQUIRED",
        "n_intervals": n, "dr_m": q2.RADIUS_M / n,
        "plateau_temperature_C": temperature_C,
        "plateau_moisture_kg_per_kg": moisture,
        "h_W_m2K": q2.H_BASE, "hm_m_s": q2.HM_BASE,
        "settings_name": settings_name, "settings": settings,
        "observed_boundary_source": str(q3.DEFAULT_BOUNDARY.resolve()),
        "observed_boundary_end_s": q3.PLATEAU_START_S,
        "plateau_start_s": q3.PLATEAU_START_S,
        "event": event, "monitor": monitor,
        "mass_balance_method": fields["mass_balance_method"],
        "mass_balance_max_abs": float(np.max(np.abs(residual))),
        "mass_balance_max_relative": relative_balance,
        "state_continuity_max_abs": continuity,
        "checks": checks, "runtime_s_including_shared_prefix": result["runtime_s"],
    }


def run_grid_cases(n: int, settings_name: str, settings: dict[str, float],
                   selected_cases: list[str]) -> dict[str, dict[str, Any]]:
    boundary, _sheet = q2.read_boundary(q3.DEFAULT_BOUNDARY, end=q3.PLATEAU_START_S)
    diagnostic = q3.platform_diagnostic(boundary)
    if not diagnostic.get("pass", False):
        raise RuntimeError("The authoritative final-hour platform diagnostic no longer passes")
    q2_history = q3._load_q2_history()
    print(f"[{settings_name}] N={n}: integrating shared observed 0--4 h prefix", flush=True)
    first = q3.integrate_segment(
        q2.initial_state(n), n, 0, q3.PLATEAU_START_S, boundary,
        rtol=settings["rtol"], atol=settings["atol"],
        max_step=settings["max_step_before_s"], detect_event=False,
        h=q2.H_BASE, hm=q2.HM_BASE,
    )
    summaries: dict[str, dict[str, Any]] = {}
    for case_name in selected_cases:
        temperature_C, moisture = CASES[case_name]
        started = time.perf_counter()
        try:
            result, continuity = continue_case(first, n, temperature_C, moisture, settings, diagnostic)
            summary = summarize_case(case_name, n, temperature_C, moisture,
                                     settings_name, settings, result, boundary,
                                     q2_history, continuity)
            summary["case_runtime_s_excluding_shared_prefix"] = time.perf_counter() - started
        except Exception as exc:  # preserve a failed row instead of silently dropping it
            summary = {
                "case": case_name, "status": "FAILED", "error": repr(exc),
                "n_intervals": n, "plateau_temperature_C": temperature_C,
                "plateau_moisture_kg_per_kg": moisture,
                "settings_name": settings_name, "settings": settings,
            }
        summaries[case_name] = summary
        atomic_json(OUTPUT / "cases" / f"{settings_name}_n{n}_{case_name}.json", summary)
        event_h = summary.get("event", {}).get("time_h")
        print(f"[{settings_name}] N={n} {case_name}: {summary['status']}, t={event_h}", flush=True)
    return summaries


def convergence_summary(base: dict[int, dict[str, dict[str, Any]]]) -> dict[str, Any]:
    archived = json.loads((ROOT / "results/q3/q3_convergence.json").read_text(encoding="utf-8"))
    archived_5120 = float(archived["spatial"]["event"]["coarse_time_s"])
    archived_10240 = float(archived["spatial"]["event"]["fine_time_s"])
    rows = [
        {"n_intervals": 2560, "dr_m": q2.RADIUS_M / 2560,
         "time_s": base[2560]["nominal"]["event"]["time_s"], "source": "new_full_audit"},
        {"n_intervals": 5120, "dr_m": q2.RADIUS_M / 5120,
         "time_s": base[5120]["nominal"]["event"]["time_s"], "source": "new_full_audit",
         "archived_time_s": archived_5120},
        {"n_intervals": 10240, "dr_m": q2.RADIUS_M / 10240,
         "time_s": archived_10240, "source": "accepted_q3_convergence_base"},
    ]
    for row in rows:
        row["time_h"] = row["time_s"] / 3600.0
    comparisons = []
    for coarse, fine in zip(rows[:-1], rows[1:]):
        signed = fine["time_s"] - coarse["time_s"]
        comparisons.append({
            "coarse_n": coarse["n_intervals"], "fine_n": fine["n_intervals"],
            "signed_difference_s": signed, "absolute_difference_s": abs(signed),
            "relative_difference": abs(signed) / abs(fine["time_s"]),
            "pass": abs(signed) < q3.CONVERGENCE_ABS_TIME_S
                and abs(signed) / abs(fine["time_s"]) < q3.CONVERGENCE_REL_TIME,
        })
    archive_difference_s = abs(rows[1]["time_s"] - archived_5120)
    trend = {
        "both_pairs_pass": all(x["pass"] for x in comparisons),
        "same_signed_direction": comparisons[0]["signed_difference_s"] * comparisons[1]["signed_difference_s"] > 0,
        "later_increment_smaller": comparisons[1]["absolute_difference_s"] < comparisons[0]["absolute_difference_s"],
        "new_5120_matches_archive": archive_difference_s <= 1.0e-6,
    }
    result = {"settings_name": "base", "settings": BASE, "rows": rows,
              "comparisons": comparisons, "trend": trend,
              "new_5120_archive_absolute_difference_s": archive_difference_s,
              "pass": all(trend.values())}
    atomic_json(OUTPUT / "q3_supplement_grid_convergence.json", result)
    with (OUTPUT / "q3_supplement_grid_convergence.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["n_intervals", "dr_m", "time_s", "time_h", "source"])
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in writer.fieldnames})
    return result


def sensitivity_summary(base: dict[int, dict[str, dict[str, Any]]]) -> dict[str, Any]:
    rows = []
    for n in (2560, 5120):
        nominal = base[n]["nominal"]["event"]["time_s"]
        for name, (temperature_C, moisture) in CASES.items():
            item = base[n][name]
            time_s = item.get("event", {}).get("time_s")
            delta = None if time_s is None else time_s - nominal
            rows.append({
                "case": name, "n_intervals": n,
                "plateau_temperature_C": temperature_C,
                "plateau_moisture_kg_per_kg": moisture,
                "time_s": time_s, "time_h": None if time_s is None else time_s / 3600.0,
                "delta_time_s_vs_nominal_same_n": delta,
                "delta_percent_vs_nominal_same_n": None if delta is None else 100.0 * delta / nominal,
                "status": item["status"],
                "mass_balance_max_relative": item.get("mass_balance_max_relative"),
                "event_bracket_width_s": item.get("event", {}).get("grid_event_bracket_width_s"),
                "monotonicity_max_positive_neighbor_jump": item.get("monitor", {}).get("max_positive_neighbor_jump_kg_per_kg"),
            })
    grid_checks = []
    for name in CASES:
        a = base[2560][name]["event"]["time_s"]
        b = base[5120][name]["event"]["time_s"]
        grid_checks.append({"case": name, "absolute_difference_s": abs(b-a),
                            "relative_difference": abs(b-a)/abs(b),
                            "pass": abs(b-a) < q3.CONVERGENCE_ABS_TIME_S
                                and abs(b-a)/abs(b) < q3.CONVERGENCE_REL_TIME})
    candidates = [r for r in rows if r["n_intervals"] == 5120 and r["case"] != "nominal"]
    largest = max(candidates, key=lambda r: abs(r["delta_time_s_vs_nominal_same_n"]))["case"]
    result = {"settings_name": "base", "settings": BASE, "rows": rows,
              "grid_checks": grid_checks, "all_grid_checks_pass": all(x["pass"] for x in grid_checks),
              "largest_effect_case_n5120": largest}
    atomic_json(OUTPUT / "q3_supplement_boundary_sensitivity.json", result)
    with (OUTPUT / "q3_supplement_boundary_sensitivity.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)
    return result


def temporal_summary(base5120: dict[str, dict[str, Any]], tight: dict[str, dict[str, Any]], largest: str) -> dict[str, Any]:
    comparisons = {}
    for name in ("nominal", largest):
        a = base5120[name]["event"]["time_s"]
        b = tight[name]["event"]["time_s"]
        comparisons[name] = {"base_time_s": a, "tight_time_s": b,
                             "absolute_difference_s": abs(b-a),
                             "relative_difference": abs(b-a)/abs(b),
                             "pass": abs(b-a) < q3.CONVERGENCE_ABS_TIME_S
                                 and abs(b-a)/abs(b) < q3.CONVERGENCE_REL_TIME}
    delta_base = base5120[largest]["event"]["time_s"] - base5120["nominal"]["event"]["time_s"]
    delta_tight = tight[largest]["event"]["time_s"] - tight["nominal"]["event"]["time_s"]
    result = {"largest_effect_case": largest, "base_settings": BASE, "tight_settings": TIGHT,
              "comparisons": comparisons, "delta_base_s": delta_base,
              "delta_tight_s": delta_tight,
              "delta_change_after_refinement_s": delta_tight-delta_base,
              "pass": all(x["pass"] for x in comparisons.values())}
    atomic_json(OUTPUT / "q3_supplement_temporal_refinement.json", result)
    return result


def diffusion_analysis() -> dict[str, Any]:
    source = ROOT / "results/q3"
    with np.load(source / "q3_fields.npz", allow_pickle=False) as payload:
        times = payload["time_s"].astype(float)
        temperature = payload["temperature_C"][:, [0, 10, 20]].copy()
        moisture = payload["moisture_kg_per_kg"][:, [0, 10, 20]].copy()
    with np.load(source / "q3_event.npz", allow_pickle=False) as payload:
        end_time = float(payload["time_s"])
        indices = [0, 5120, 10240]
        end_temperature = payload["temperature_C"][indices]
        end_moisture = payload["moisture_kg_per_kg"][indices]
    times = np.r_[times, end_time]
    temperature = np.vstack([temperature, end_temperature])
    moisture = np.vstack([moisture, end_moisture])
    diffusivity = q2.material_properties(moisture, temperature)[3]
    labels = ["center", "half_radius", "surface"]
    with (OUTPUT / "q3_supplement_diffusivity.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "time_h", "D_center_m2_s", "D_half_radius_m2_s", "D_surface_m2_s"])
        for t, row in zip(times, diffusivity):
            writer.writerow([t, t/3600.0, *row])
    summary: dict[str, Any] = {"source_fields": "results/q3/q3_fields.npz",
        "source_event": "results/q3/q3_event.npz", "radii_cm": [0.0, 1.0, 2.0],
        "formula_source": "q2_solver.material_properties q2 branch", "series": {}}
    for j, label in enumerate(labels):
        peak = int(np.argmax(diffusivity[:, j]))
        summary["series"][label] = {"peak_time_s": times[peak], "peak_time_h": times[peak]/3600.0,
            "peak_D_m2_s": diffusivity[peak, j], "endpoint_D_m2_s": diffusivity[-1, j]}
    summary["endpoint_center_surface_ratio"] = diffusivity[-1, 0] / diffusivity[-1, 2]
    np.savez_compressed(OUTPUT / "q3_supplement_diffusivity.npz", time_s=times,
                        temperature_C=temperature, moisture_kg_per_kg=moisture,
                        diffusivity_m2_s=diffusivity, radius_cm=np.array([0.0, 1.0, 2.0]))
    atomic_json(OUTPUT / "q3_supplement_diffusivity.json", summary)
    return summary


def configure_plotting() -> None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    chinese = next((x for x in ("Microsoft YaHei", "SimHei", "SimSun") if x in available), "DejaVu Sans")
    plt.rcParams.update({"font.family": chinese, "font.size": 10, "axes.unicode_minus": False,
        "pdf.fonttype": 42, "figure.dpi": 150, "savefig.dpi": 300,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": .18, "grid.linewidth": .5})


def save_figure(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUTPUT / f"{name}.pdf", bbox_inches="tight", pad_inches=.08)
    fig.savefig(OUTPUT / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=.08)
    plt.close(fig)


def plot_sensitivity(summary: dict[str, Any]) -> None:
    configure_plotting()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0), layout="constrained")
    styles = {2560: ("--", "o", "#56B4E9"), 5120: ("-", "s", "#0072B2")}
    for n, (line, marker, color) in styles.items():
        rows = [r for r in summary["rows"] if r["n_intervals"] == n]
        by_name = {r["case"]: r for r in rows}
        tx = np.array([49, 50, 51], dtype=float)
        ty = np.array([by_name["temperature_49"]["time_h"], by_name["nominal"]["time_h"], by_name["temperature_51"]["time_h"]])
        mx = np.array([.045, .050, .055])
        my = np.array([by_name["moisture_0045"]["time_h"], by_name["nominal"]["time_h"], by_name["moisture_0055"]["time_h"]])
        axes[0].plot(tx, ty, line, marker=marker, color=color, label=f"N={n}")
        axes[1].plot(mx, my, line, marker=marker, color=color, label=f"N={n}")
    axes[0].set(xlabel="长期环境温度/℃", ylabel="临界烘干时间/h", title="(a) 温度扰动")
    axes[1].set(xlabel="长期环境水分浓度/(kg/kg)", ylabel="临界烘干时间/h", title="(b) 环境水分扰动")
    for ax in axes:
        ax.legend(frameon=False); ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    save_figure(fig, "q3_supplement_boundary_sensitivity")


def plot_diffusivity() -> None:
    configure_plotting()
    with np.load(OUTPUT / "q3_supplement_diffusivity.npz", allow_pickle=False) as payload:
        hours = payload["time_s"] / 3600.0
        values = payload["diffusivity_m2_s"]
    fig, ax = plt.subplots(figsize=(6.8, 4.1), layout="constrained")
    for j, (label, color, line) in enumerate([
        ("圆心", "#0072B2", "-"), (r"$r=0.5R_0$", "#009E73", "--"), ("表面", "#D55E00", "-.")]):
        ax.semilogy(hours, values[:, j], line, color=color, linewidth=1.7, label=label)
    ax.axvline(4, color="#666666", linestyle=":", linewidth=1.1, label="长期边界起点")
    ax.set(xlabel="时间/h", ylabel=r"有效扩散系数 $D$/(m$^2$/s)",
           xlim=(0, hours[-1]), title="不同径向位置的有效扩散系数")
    ax.legend(ncol=2, frameon=False)
    save_figure(fig, "q3_supplement_diffusivity")


def write_audit_bundle(grid: dict[str, Any], sensitivity: dict[str, Any],
                       temporal: dict[str, Any], diffusion: dict[str, Any],
                       base: dict[int, dict[str, dict[str, Any]]],
                       tight: dict[str, dict[str, Any]]) -> None:
    case_pass = all(x["status"] == "PASSED" for bundle in base.values() for x in bundle.values())
    tight_pass = all(x["status"] == "PASSED" for x in tight.values())
    status = "READY_FOR_INDEPENDENT_AUDIT" if (case_pass and tight_pass and grid["pass"]
        and sensitivity["all_grid_checks_pass"] and temporal["pass"]) else "REVIEW_REQUIRED"
    bundle = {"status": status, "created_by": "q3_supplement_experiments.py",
        "formal_results_modified": False, "paper_modified": False,
        "grid_convergence": grid, "boundary_sensitivity": sensitivity,
        "temporal_refinement": temporal, "diffusivity": diffusion,
        "limitations": [
            "Long-term boundary changes are specified engineering scenarios, not confidence intervals.",
            "The endpoint remains a conditional prediction under a constant post-4-hour environment.",
            "Local diffusivity curves support but do not uniquely prove a controlling resistance mechanism.",
            "Mass balance is closure of the effective moisture equation, not total physical mass validation.",
        ]}
    atomic_json(OUTPUT / "q3_supplement_audit_bundle.json", bundle)


def regenerate_summaries() -> dict[str, Any]:
    """Rebuild lightweight summaries from completed case records only."""
    def load_case(setting: str, n: int, name: str) -> dict[str, Any]:
        path = OUTPUT / "cases" / f"{setting}_n{n}_{name}.json"
        item = json.loads(path.read_text(encoding="utf-8"))
        event = item.get("event", {})
        event.pop("state_file", None)
        event["state_storage"] = "not persisted; scalar event evidence is embedded in this case record"
        atomic_json(path, item)
        return item

    base = {
        n: {
            name: load_case("base", n, name)
            for name in CASES
        }
        for n in (2560, 5120)
    }
    grid = convergence_summary(base)
    sensitivity = sensitivity_summary(base)
    largest = sensitivity["largest_effect_case_n5120"]
    tight = {
        name: load_case("tight", 5120, name)
        for name in ("nominal", largest)
    }
    temporal = temporal_summary(base[5120], tight, largest)
    diffusion = json.loads((OUTPUT / "q3_supplement_diffusivity.json").read_text(encoding="utf-8"))
    plot_sensitivity(sensitivity)
    write_audit_bundle(grid, sensitivity, temporal, diffusion, base, tight)
    return {"status": "SUMMARIES_REGENERATED", "largest_effect_case": largest}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["all", "diffusion", "summarize"], default="all")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "summarize":
        print(json.dumps(regenerate_summaries(), ensure_ascii=False)); return
    diffusion = diffusion_analysis()
    plot_diffusivity()
    if args.stage == "diffusion":
        print(json.dumps(json_safe(diffusion), ensure_ascii=False)); return
    base = {n: run_grid_cases(n, "base", BASE, list(CASES)) for n in (2560, 5120)}
    if not all(x["status"] == "PASSED" for bundle in base.values() for x in bundle.values()):
        raise RuntimeError("One or more base cases failed their evidence checks; see per-case JSON")
    grid = convergence_summary(base)
    sensitivity = sensitivity_summary(base)
    largest = sensitivity["largest_effect_case_n5120"]
    tight = run_grid_cases(5120, "tight", TIGHT, ["nominal", largest])
    temporal = temporal_summary(base[5120], tight, largest)
    plot_sensitivity(sensitivity)
    write_audit_bundle(grid, sensitivity, temporal, diffusion, base, tight)
    print(json.dumps({"status": "READY_FOR_INDEPENDENT_AUDIT", "largest_effect_case": largest,
                      "output": str(OUTPUT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
