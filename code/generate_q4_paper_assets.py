"""Generate unrounded CSV sources and publication-style Q4 figures.

The script consumes only a PASS formal run.  It intentionally stops when the
formal validation marker is absent or failed, so figures cannot be generated
from a partial or exploratory calculation.  It does not create or modify a
paper file or an Excel workbook.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


TABLE6_POSITIONS_CM = np.arange(9, dtype=float) * 0.5
COLOURS = ("#174A70", "#2584A6", "#77AABD", "#D88A41", "#963E3E", "#6A4C93", "#3B7A57", "#B06A8A", "#6B6E70")


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


def latest_run(base: Path) -> Path:
    runs = [p for p in base.glob("q4_formal_*") if p.is_dir() and (p / "q4_formal_validation.json").exists()]
    if not runs:
        raise FileNotFoundError(f"No formal Q4 run under {base}")
    return max(runs, key=lambda p: p.stat().st_mtime)


def load_fields(case_dir: Path) -> dict[str, np.ndarray]:
    with np.load(case_dir / "q4_fields.npz", allow_pickle=False) as data:
        fields = {key: np.asarray(data[key]) for key in data.files}
    required = {"time_s", "positions_cm", "temperature_C", "moisture_kg_per_kg", "surface_moisture_kg_per_kg", "surface_temperature_C", "radius_m", "mean_moisture", "domain_mask"}
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError(f"Final fields are missing keys: {missing}")
    return fields


def row_value(value: Any) -> Any:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return value
    return "" if not np.isfinite(x) else x


def position_indices(positions: np.ndarray, selected: np.ndarray) -> list[int]:
    return [int(np.where(np.isclose(positions, value, rtol=0.0, atol=1e-10))[0][0]) for value in selected]


def write_final_fields(out: Path, fields: dict[str, np.ndarray]) -> Path:
    path = out / "final_fields.csv"
    positions = np.asarray(fields["positions_cm"], dtype=float)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "time_s", "time_h", "radius_m", "radius_cm", "mean_moisture_kg_per_kg",
            "position_cm", "domain_in_material", "temperature_C", "moisture_kg_per_kg",
            "surface_temperature_C", "surface_moisture_kg_per_kg",
        ])
        times = np.asarray(fields["time_s"], dtype=float)
        radius = np.asarray(fields["radius_m"], dtype=float)
        T = np.asarray(fields["temperature_C"], dtype=float)
        C = np.asarray(fields["moisture_kg_per_kg"], dtype=float)
        masks = np.asarray(fields["domain_mask"], dtype=bool)
        mean = np.asarray(fields["mean_moisture"], dtype=float)
        surface_T = np.asarray(fields["surface_temperature_C"], dtype=float)
        surface_C = np.asarray(fields["surface_moisture_kg_per_kg"], dtype=float)
        for i, t in enumerate(times):
            for j, p in enumerate(positions):
                writer.writerow([
                    row_value(t), row_value(t / 3600.0), row_value(radius[i]), row_value(radius[i] * 100.0),
                    row_value(mean[i]), row_value(p), bool(masks[i, j]), row_value(T[i, j]), row_value(C[i, j]),
                    row_value(surface_T[i]), row_value(surface_C[i]),
                ])
    return path


def table6_arrays(fields: dict[str, np.ndarray], endpoint_s: float) -> dict[str, np.ndarray]:
    times = np.asarray(fields["time_s"], dtype=float)
    wanted = [float(6 * 3600 * i) for i in range(1, int(math.floor(endpoint_s / (6 * 3600) + 1e-12)) + 1)]
    if not wanted or abs(wanted[-1] - endpoint_s) > 1e-8:
        wanted.append(float(endpoint_s))
    selected: list[int] = []
    for t in wanted:
        found = np.where(np.isclose(times, t, rtol=0.0, atol=1e-7))[0]
        if len(found):
            selected.append(int(found[0]))
    pos = np.asarray(fields["positions_cm"], dtype=float)
    indices = position_indices(pos, TABLE6_POSITIONS_CM)
    return {
        "row_index": np.asarray(selected, dtype=int),
        "time_s": times[selected],
        "radius_cm": np.asarray(fields["radius_m"], dtype=float)[selected] * 100.0,
        "temperature_C": np.asarray(fields["temperature_C"], dtype=float)[selected][:, indices],
        "moisture_kg_per_kg": np.asarray(fields["moisture_kg_per_kg"], dtype=float)[selected][:, indices],
        "domain_mask": np.asarray(fields["domain_mask"], dtype=bool)[selected][:, indices],
        "surface_temperature_C": np.asarray(fields["surface_temperature_C"], dtype=float)[selected],
        "surface_moisture_kg_per_kg": np.asarray(fields["surface_moisture_kg_per_kg"], dtype=float)[selected],
    }


def write_table6(out: Path, fields: dict[str, np.ndarray], endpoint_s: float) -> Path:
    table = table6_arrays(fields, endpoint_s)
    path = out / "table6_source.csv"
    headers = ["time_s", "time_h", "radius_cm"]
    headers += [f"temperature_{p:.1f}cm_C" for p in TABLE6_POSITIONS_CM]
    headers += [f"moisture_{p:.1f}cm_kg_per_kg" for p in TABLE6_POSITIONS_CM]
    headers += [f"domain_{p:.1f}cm" for p in TABLE6_POSITIONS_CM]
    headers += ["surface_temperature_C", "surface_moisture_kg_per_kg"]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for i, t in enumerate(table["time_s"]):
            row: list[Any] = [row_value(t), row_value(t / 3600.0), row_value(table["radius_cm"][i])]
            row.extend(row_value(x) for x in table["temperature_C"][i])
            row.extend(row_value(x) if mask else "" for x, mask in zip(table["moisture_kg_per_kg"][i], table["domain_mask"][i]))
            row.extend(bool(x) for x in table["domain_mask"][i])
            row.extend([row_value(table["surface_temperature_C"][i]), row_value(table["surface_moisture_kg_per_kg"][i])])
            writer.writerow(row)
    return path


def accepted_max_history(case_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    manifest_path = case_dir / "q4_blocks_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    times: list[np.ndarray] = []
    values: list[np.ndarray] = []
    for item in manifest.get("segment_records", []):
        path = Path(item["evidence_file"])
        with np.load(path, allow_pickle=False) as data:
            times.append(np.asarray(data["accepted_time_s"], dtype=float))
            values.append(np.asarray(data["accepted_g_N"], dtype=float) + 0.15)
    if not times:
        raise ValueError("No accepted-step event evidence found")
    t = np.concatenate(times)
    c = np.concatenate(values)
    order = np.argsort(t, kind="mergesort")
    t, c = t[order], c[order]
    keep = np.concatenate(([True], np.diff(t) > 1e-9))
    return t[keep], c[keep]


def write_max_history(out: Path, case_dir: Path) -> Path:
    times, maximum = accepted_max_history(case_dir)
    path = out / "maximum_moisture.csv"
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "time_h", "maximum_moisture_kg_per_kg", "g_N"])
        for t, value in zip(times, maximum):
            writer.writerow([row_value(t), row_value(t / 3600.0), row_value(value), row_value(value - 0.15)])
    return path


def write_radius_distribution(out: Path, fields: dict[str, np.ndarray], endpoint_s: float) -> Path:
    table = table6_arrays(fields, endpoint_s)
    path = out / "radius_distribution.csv"
    positions = TABLE6_POSITIONS_CM
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time_s", "time_h", "radius_cm", "position_cm", "domain_in_material", "temperature_C", "moisture_kg_per_kg"])
        for i, t in enumerate(table["time_s"]):
            for j, p in enumerate(positions):
                writer.writerow([
                    row_value(t), row_value(t / 3600.0), row_value(table["radius_cm"][i]), row_value(p),
                    bool(table["domain_mask"][i, j]), row_value(table["temperature_C"][i, j]),
                    row_value(table["moisture_kg_per_kg"][i, j]) if table["domain_mask"][i, j] else "",
                ])
    return path


def plot_radius_profiles(out: Path, fields: dict[str, np.ndarray], endpoint_s: float) -> Path:
    table = table6_arrays(fields, endpoint_s)
    times_h = np.asarray(table["time_s"], dtype=float) / 3600.0
    radius_cm = np.asarray(table["radius_cm"], dtype=float)
    positions = TABLE6_POSITIONS_CM
    fig, (ax_r, ax_c, ax_t) = plt.subplots(1, 3, figsize=(11.0, 3.3), constrained_layout=True)
    ax_r.plot(times_h, radius_cm, color=COLOURS[0], linewidth=1.8)
    ax_r.set(xlabel="Drying time (h)", ylabel="Radius (cm)", title="Measured radius and extension")
    ax_r.grid(alpha=0.2)
    for i, time_h in enumerate(times_h):
        mask = np.asarray(table["domain_mask"])[i]
        ax_c.plot(positions[mask], np.asarray(table["moisture_kg_per_kg"])[i][mask], marker="o", markersize=2.2, linewidth=1.0, color=COLOURS[i % len(COLOURS)], label=f"{time_h:g} h")
    ax_c.axhline(0.15, color=COLOURS[4], linestyle="--", linewidth=0.8, label="Target 0.15")
    ax_c.set(xlabel="Physical position (cm)", ylabel="Moisture C (kg/kg)", title="Radial moisture profiles")
    ax_c.grid(alpha=0.2)
    ax_c.legend(frameon=False, fontsize=6.5, ncol=2)
    for i, time_h in enumerate(times_h):
        mask = np.asarray(table["domain_mask"])[i]
        ax_t.plot(positions[mask], np.asarray(table["temperature_C"])[i][mask], marker="o", markersize=2.2, linewidth=1.0, color=COLOURS[i % len(COLOURS)], label=f"{time_h:g} h")
    ax_t.set(xlabel="Physical position (cm)", ylabel="Temperature (°C)", title="Radial temperature profiles")
    ax_t.grid(alpha=0.2)
    path = out / "q4_radius_and_distribution.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_max_history(out: Path, case_dir: Path, endpoint_s: float) -> tuple[Path, Path, Path]:
    times, maximum = accepted_max_history(case_dir)
    time_h = times / 3600.0
    fig, (ax, axz) = plt.subplots(1, 2, figsize=(8.0, 3.3), constrained_layout=True)
    ax.plot(time_h, maximum, color=COLOURS[0], linewidth=1.3)
    ax.axhline(0.15, color=COLOURS[4], linestyle="--", linewidth=0.8)
    ax.axvline(endpoint_s / 3600.0, color=COLOURS[4], linestyle=":", linewidth=0.8)
    ax.set(xlabel="Drying time (h)", ylabel="Maximum moisture (kg/kg)", title="Full-time threshold history")
    ax.grid(alpha=0.2)
    zoom = times >= endpoint_s - 900.0
    axz.plot((times[zoom] - endpoint_s), maximum[zoom] - 0.15, color=COLOURS[0], linewidth=1.3)
    axz.axhline(0.0, color=COLOURS[4], linestyle="--", linewidth=0.8)
    axz.axvline(0.0, color=COLOURS[4], linestyle=":", linewidth=0.8)
    axz.set(xlabel="Time relative to endpoint (s)", ylabel="g_N (kg/kg)", title="Endpoint crossing (15 min)")
    axz.grid(alpha=0.2)
    combined = out / "q4_max_moisture_full_and_endpoint.png"
    fig.savefig(combined, dpi=220, bbox_inches="tight")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.3, 3.2), constrained_layout=True)
    ax.plot(time_h, maximum, color=COLOURS[0], linewidth=1.3)
    ax.axhline(0.15, color=COLOURS[4], linestyle="--", linewidth=0.8)
    ax.set(xlabel="Drying time (h)", ylabel="Maximum moisture (kg/kg)", title="Full-time maximum moisture")
    ax.grid(alpha=0.2)
    full = out / "q4_max_moisture_full.png"
    fig.savefig(full, dpi=220, bbox_inches="tight")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.3, 3.2), constrained_layout=True)
    ax.plot(times[zoom] - endpoint_s, maximum[zoom] - 0.15, color=COLOURS[0], linewidth=1.3)
    ax.axhline(0.0, color=COLOURS[4], linestyle="--", linewidth=0.8)
    ax.axvline(0.0, color=COLOURS[4], linestyle=":", linewidth=0.8)
    ax.set(xlabel="Time relative to endpoint (s)", ylabel="g_N (kg/kg)", title="Endpoint crossing detail")
    ax.grid(alpha=0.2)
    zoom_path = out / "q4_max_moisture_endpoint_zoom.png"
    fig.savefig(zoom_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return full, zoom_path, combined


def write_ad_comparison(out: Path, validation: dict[str, Any]) -> Path:
    factorial = validation.get("factorial") or {}
    times = factorial.get("event_time_s") or {}
    path = out / "ad_comparison.csv"
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["case", "property_model", "geometry", "event_time_s", "event_time_h"])
        metadata = {"A": ("appendix3", "fixed_R0"), "B": ("appendix3", "attachment2_radius"), "C": ("appendix4", "fixed_R0"), "D": ("appendix4", "attachment2_radius")}
        for case in ("A", "B", "C", "D"):
            value = times.get(case)
            writer.writerow([case, *metadata[case], row_value(value), row_value(None if value is None else float(value) / 3600.0)])
        writer.writerow(["delta_property_C_minus_A", "", "", row_value(factorial.get("delta_property_s")), ""])
        writer.writerow(["delta_shrink_3_B_minus_A", "", "", row_value(factorial.get("delta_shrink_3_s")), ""])
        writer.writerow(["delta_shrink_4_D_minus_C", "", "", row_value(factorial.get("delta_shrink_4_s")), ""])
        writer.writerow(["delta_interaction", "", "", row_value(factorial.get("delta_interaction_s")), ""])
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=None)
    args = parser.parse_args(argv)
    base = Path(__file__).resolve().parents[1] / "results" / "q4" / "formal"
    run_root = args.run_root.resolve() if args.run_root else latest_run(base)
    validation_path = run_root / "q4_formal_validation.json"
    if not validation_path.exists():
        raise FileNotFoundError("Formal validation marker is missing; no assets generated")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if validation.get("status") != "PASS":
        raise RuntimeError("Formal validation is not PASS; no assets generated")
    out = run_root / "assets"
    out.mkdir(parents=True, exist_ok=False) if not out.exists() else None
    endpoint_s = float(validation["D_event_time_s"])
    final_n = int(validation["final_grid_n"])
    case_dir = run_root / "space" / f"D_n{final_n:05d}"
    fields = load_fields(case_dir)
    generated = {
        "final_fields": str(write_final_fields(out, fields)),
        "table6_source": str(write_table6(out, fields, endpoint_s)),
        "radius_distribution": str(write_radius_distribution(out, fields, endpoint_s)),
        "maximum_moisture": str(write_max_history(out, case_dir)),
        "ad_comparison": str(write_ad_comparison(out, validation)),
    }
    generated["radius_and_distribution_figure"] = str(plot_radius_profiles(out, fields, endpoint_s))
    full, zoom, combined = plot_max_history(out, case_dir, endpoint_s)
    generated["maximum_full_figure"] = str(full)
    generated["maximum_endpoint_zoom_figure"] = str(zoom)
    generated["maximum_combined_figure"] = str(combined)
    manifest = {"status": "PASS", "run_id": run_root.name, "final_grid_n": final_n, "D_event_time_s": endpoint_s, "generated_at": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(timespec="seconds"), "files": generated}
    write_json(out / "asset_manifest.json", manifest)
    print(json.dumps(safe(manifest), ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
