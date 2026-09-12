"""Render the unified six-panel Q2 field figure from accepted solver output."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager, ticker
import numpy as np


PROFILE_SECONDS = (1800, 3600, 5400, 7200, 9000, 10800)
PAPER_RADII_CM = (0.0, 0.5, 1.0, 1.5, 2.0)
PAPER_NODE_INDEX = (0, 5, 10, 15, 20)


def _configure_style() -> None:
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next(
        (name for name in ("Microsoft YaHei", "SimHei", "SimSun") if name in available),
        "DejaVu Sans",
    )
    plt.rcParams.update(
        {
            "font.family": family,
            "font.size": 8.2,
            "axes.titlesize": 9.2,
            "axes.titleweight": "bold",
            "axes.labelsize": 8.2,
            "legend.fontsize": 6.8,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
            "axes.linewidth": 0.65,
            "lines.linewidth": 1.35,
        }
    )


def _style_2d_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(color="#9aa7b2", alpha=0.20, linewidth=0.45)
    ax.tick_params(length=2.5, width=0.6, pad=1.5)
    ax.set_axisbelow(True)


def _style_surface_axis(ax, *, zlim: tuple[float, float], zticks: list[float]) -> None:
    ax.set_xlim(0, 3)
    ax.set_ylim(0, 2)
    ax.set_zlim(*zlim)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([0, 1, 2])
    ax.set_zticks(zticks)
    ax.zaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _: f"{value:g}"))
    ax.tick_params(labelsize=6.4, pad=-1, length=1.5)
    ax.set_xlabel("时间/h", labelpad=1)
    ax.set_ylabel("径向位置/cm", labelpad=1)
    ax.set_proj_type("ortho")
    ax.set_box_aspect((1.18, 0.90, 0.72))
    ax.view_init(elev=27, azim=-125)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((0.965, 0.975, 0.982, 0.48))
        axis.pane.set_edgecolor("#d7dfe5")
        axis._axinfo["grid"].update(color="#d7dfe5", linewidth=0.38)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "results" / "q2"
    source_path = output / "q2_fields.npz"
    with np.load(source_path, allow_pickle=False) as fields:
        time_s = fields["time_s"].copy()
        radius_cm = fields["radius_cm"].copy()
        temperature = fields["temperature_C"].copy()
        moisture = fields["moisture_kg_per_kg"].copy()

    if temperature.shape != moisture.shape or temperature.shape != (10801, 21):
        raise ValueError(f"unexpected Q2 field shape: {temperature.shape}, {moisture.shape}")
    if not (np.isfinite(temperature).all() and np.isfinite(moisture).all()):
        raise ValueError("Q2 fields contain non-finite values")
    if not np.array_equal(time_s, np.arange(10801)):
        raise ValueError("Q2 time grid is not the required 0--10800 s grid")
    if not np.allclose(radius_cm, np.arange(21) / 10.0, rtol=0.0, atol=1e-14):
        raise ValueError("Q2 radius grid is not the required 0--2 cm output grid")

    _configure_style()
    profile_colors = ["#173F5F", "#0072B2", "#56B4E9", "#E69F00", "#D55E00", "#7A2533"]
    radial_colors = ["#173F5F", "#0072B2", "#56B4E9", "#E69F00", "#7A2533"]
    line_styles = ["-", "--", "-.", ":", (0, (5, 1.5)), (0, (3, 1, 1, 1))]
    heat_cmap = colors.LinearSegmentedColormap.from_list(
        "q2_heat", ["#fff4d5", "#F0E442", "#E69F00", "#D55E00"]
    )
    moisture_cmap = colors.LinearSegmentedColormap.from_list(
        "q2_moisture", ["#edf6fb", "#56B4E9", "#0072B2", "#173F5F"]
    )

    fig = plt.figure(figsize=(7.35, 4.95), facecolor="white")
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.0, 1.0, 1.13),
        left=0.072,
        right=0.988,
        bottom=0.095,
        top=0.955,
        wspace=0.34,
        hspace=0.39,
    )
    ax_tp = fig.add_subplot(grid[0, 0])
    ax_tt = fig.add_subplot(grid[0, 1])
    ax_ts = fig.add_subplot(grid[0, 2], projection="3d")
    ax_cp = fig.add_subplot(grid[1, 0])
    ax_ct = fig.add_subplot(grid[1, 1])
    ax_cs = fig.add_subplot(grid[1, 2], projection="3d")

    for second, color, style in zip(PROFILE_SECONDS, profile_colors, line_styles):
        label = f"{second / 3600:g} h"
        ax_tp.plot(radius_cm, temperature[second], color=color, linestyle=style, label=label)
        ax_cp.plot(radius_cm, moisture[second], color=color, linestyle=style, label=label)
    time_h = time_s / 3600.0
    for index, radial_position, color in zip(PAPER_NODE_INDEX, PAPER_RADII_CM, radial_colors):
        label = f"{radial_position:g} cm"
        ax_tt.plot(time_h, temperature[:, index], color=color, label=label)
        ax_ct.plot(time_h, moisture[:, index], color=color, label=label)

    for ax in (ax_tp, ax_tt, ax_cp, ax_ct):
        _style_2d_axis(ax)
    ax_tp.set(
        title="(a) 温度径向剖面",
        xlabel="径向位置/cm",
        ylabel="温度/°C",
        xlim=(0, 2),
        xticks=[0, 0.5, 1, 1.5, 2],
    )
    ax_tt.set(
        title="(b) 温度时间响应",
        xlabel="时间/h",
        ylabel="温度/°C",
        xlim=(0, 3),
        xticks=[0, 1, 2, 3],
    )
    ax_cp.set(
        title="(d) 含水率径向剖面",
        xlabel="径向位置/cm",
        ylabel="干基含水率/(kg/kg)",
        xlim=(0, 2),
        xticks=[0, 0.5, 1, 1.5, 2],
    )
    ax_ct.set(
        title="(e) 含水率时间响应",
        xlabel="时间/h",
        ylabel="干基含水率/(kg/kg)",
        xlim=(0, 3),
        xticks=[0, 1, 2, 3],
    )
    ax_tp.legend(loc="lower right", ncol=2, frameon=False, columnspacing=0.65, handlelength=1.7)
    ax_cp.legend(loc="lower left", ncol=2, frameon=False, columnspacing=0.65, handlelength=1.7)
    ax_tt.legend(loc="lower right", ncol=2, frameon=False, columnspacing=0.65, handlelength=1.7)
    ax_ct.legend(loc="lower left", ncol=2, frameon=False, columnspacing=0.65, handlelength=1.7)

    surface_indices = np.unique(np.r_[np.arange(61), np.arange(60, 10801, 60)])
    surface_time_h, surface_radius_cm = np.meshgrid(
        time_h[surface_indices], radius_cm, indexing="ij"
    )
    surface_specs = (
        (
            ax_ts,
            temperature,
            heat_cmap,
            (28.0, 51.0),
            [28, 34, 40, 46, 50],
            "(c) 温度时空曲面",
            "温度/°C",
        ),
        (
            ax_cs,
            moisture,
            moisture_cmap,
            (1.0, 2.6),
            [1.0, 1.4, 1.8, 2.2, 2.55],
            "(f) 含水率时空曲面",
            "含水率/(kg/kg)",
        ),
    )
    for ax, data, cmap, zlim, zticks, title, zlabel in surface_specs:
        ax.plot_surface(
            surface_time_h,
            surface_radius_cm,
            data[surface_indices],
            rstride=1,
            cstride=1,
            cmap=cmap,
            vmin=zlim[0],
            vmax=zlim[1],
            linewidth=0,
            antialiased=False,
            shade=False,
            rasterized=True,
        )
        for second in PROFILE_SECONDS:
            ax.plot(
                np.full(radius_cm.shape, second / 3600.0),
                radius_cm,
                data[second],
                color="#263c4a",
                alpha=0.48,
                linewidth=0.48,
            )
        _style_surface_axis(ax, zlim=zlim, zticks=zticks)
        ax.set_zlabel(zlabel, labelpad=0)
        ax.set_title(title, pad=2.5)

    pdf_path = output / "q2_combined_fields.pdf"
    png_path = output / "q2_combined_fields.png"
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)

    np.savez_compressed(
        output / "q2_combined_fields_data.npz",
        profile_time_s=np.asarray(PROFILE_SECONDS, dtype=int),
        time_s=time_s,
        radius_cm=radius_cm,
        temperature_C=temperature,
        moisture_kg_per_kg=moisture,
        surface_indices=surface_indices,
    )
    metadata = {
        "source": "q2_fields.npz",
        "layout": "2 rows x 3 columns: temperature/moisture by profile/time/surface",
        "profile_times_s": list(PROFILE_SECONDS),
        "representative_radii_cm": list(PAPER_RADII_CM),
        "surface_sampling": "1 s for first minute; 60 s afterwards; all 21 radial nodes",
        "smoothing": False,
        "extrapolation": False,
        "temperature_range_C": [float(temperature.min()), float(temperature.max())],
        "moisture_range_kg_per_kg": [float(moisture.min()), float(moisture.max())],
    }
    (output / "q2_combined_fields_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps({"pdf": str(pdf_path), "png": str(png_path), **metadata}, ensure_ascii=False))


if __name__ == "__main__":
    main()
