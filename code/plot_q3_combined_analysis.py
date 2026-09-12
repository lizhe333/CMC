"""Combine the verified Q3 sensitivity and diffusivity results into one paper figure.

This script only reads existing CSV evidence. It does not rerun the PDE solver or
modify any formal numerical result.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPPLEMENT = ROOT / "results" / "q3" / "supplement"
PAPER_FIGURES = ROOT / "paper" / "figures"


def configure_plotting() -> None:
    available = {font.name for font in font_manager.fontManager.ttflist}
    chinese = next(
        (name for name in ("Microsoft YaHei", "SimHei", "SimSun") if name in available),
        "DejaVu Sans",
    )
    plt.rcParams.update(
        {
            "font.family": chinese,
            "font.size": 10,
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.18,
            "grid.linewidth": 0.5,
        }
    )


def main() -> None:
    configure_plotting()
    sensitivity = pd.read_csv(
        SUPPLEMENT / "q3_supplement_boundary_sensitivity.csv", encoding="utf-8-sig"
    )
    diffusivity = pd.read_csv(
        SUPPLEMENT / "q3_supplement_diffusivity.csv", encoding="utf-8-sig"
    )

    fig = plt.figure(figsize=(10.4, 6.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1.0, 1.05))
    ax_temp = fig.add_subplot(grid[0, 0])
    ax_moist = fig.add_subplot(grid[0, 1])
    ax_diff = fig.add_subplot(grid[1, :])

    styles = {
        2560: ("--", "o", "#56B4E9"),
        5120: ("-", "s", "#0072B2"),
    }
    for intervals, (line, marker, color) in styles.items():
        rows = sensitivity[sensitivity["n_intervals"] == intervals]
        nominal = rows[rows["case"] == "nominal"].iloc[0]
        low_temp = rows[rows["case"] == "temperature_49"].iloc[0]
        high_temp = rows[rows["case"] == "temperature_51"].iloc[0]
        low_moist = rows[rows["case"] == "moisture_0045"].iloc[0]
        high_moist = rows[rows["case"] == "moisture_0055"].iloc[0]

        ax_temp.plot(
            [49, 50, 51],
            [low_temp.time_h, nominal.time_h, high_temp.time_h],
            line,
            marker=marker,
            color=color,
            linewidth=1.6,
            markersize=4,
            label=f"N={intervals}",
        )
        ax_moist.plot(
            [0.045, 0.050, 0.055],
            [low_moist.time_h, nominal.time_h, high_moist.time_h],
            line,
            marker=marker,
            color=color,
            linewidth=1.6,
            markersize=4,
            label=f"N={intervals}",
        )

    ax_temp.set(
        title="(a) 长期温度扰动",
        xlabel="长期环境温度/℃",
        ylabel="临界烘干时间/h",
        xticks=[49, 50, 51],
    )
    ax_moist.set(
        title="(b) 环境水分扰动",
        xlabel="长期环境水分浓度/(kg/kg)",
        ylabel="临界烘干时间/h",
        xticks=[0.045, 0.050, 0.055],
    )
    ax_temp.legend(frameon=False)
    ax_moist.legend(frameon=False)
    ax_temp.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax_moist.ticklabel_format(axis="y", style="plain", useOffset=False)

    hours = diffusivity["time_h"].to_numpy()
    for column, label, color, line in (
        ("D_center_m2_s", "圆心", "#0072B2", "-"),
        ("D_half_radius_m2_s", r"$r=0.5R_0$", "#009E73", "--"),
        ("D_surface_m2_s", "表面", "#D55E00", "-."),
    ):
        ax_diff.semilogy(
            hours,
            diffusivity[column].to_numpy(),
            line,
            color=color,
            linewidth=1.7,
            label=label,
        )
    ax_diff.axvline(
        4,
        color="#666666",
        linestyle=":",
        linewidth=1.1,
        label="长期边界起点",
    )
    ax_diff.set(
        title="(c) 不同径向位置的有效扩散系数",
        xlabel="时间/h",
        ylabel=r"有效扩散系数 $D$/(m$^2$/s)",
        xlim=(0, hours[-1]),
    )
    ax_diff.legend(ncol=4, frameon=False, loc="upper right")

    for target in (
        SUPPLEMENT / "q3_combined_analysis.pdf",
        PAPER_FIGURES / "q3_combined_analysis.pdf",
    ):
        fig.savefig(target, bbox_inches="tight", pad_inches=0.08)
    for target in (
        SUPPLEMENT / "q3_combined_analysis.png",
        PAPER_FIGURES / "q3_combined_analysis.png",
    ):
        fig.savefig(target, dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


if __name__ == "__main__":
    main()
