"""Plot Attachment 1 observations and the existing linear boundary interpolation.

Read-only input; does not import or run a PDE solver or modify result files.
"""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from openpyxl import load_workbook


def main():
    root = Path(__file__).resolve().parents[1]
    source = root.parent / "附件1.xlsx"
    wb = load_workbook(source, read_only=True, data_only=True)
    try:
        rows = list(wb.worksheets[0].values)
    finally:
        wb.close()
    assert tuple(str(v).strip() for v in rows[0][:3]) == ("时间", "温度", "水分浓度")
    records = [row[:3] for row in rows[1:] if any(v is not None for v in row)]
    data = np.asarray(records, dtype=float)
    assert data.shape == (241, 3) and np.isfinite(data).all()
    time = data[:, 0]
    assert time[0] == 0 and time[-1] == 14400
    assert np.all(np.diff(time) == 60)
    seconds = np.arange(0, 14401, dtype=float)
    curves = np.column_stack([np.interp(seconds, time, data[:, j]) for j in (1, 2)])
    assert np.array_equal(curves[time.astype(int)], data[:, 1:])

    # Verify that the new display uses exactly the saved Q1 environment.
    previous = np.genfromtxt(root / "results/q1/boundary_used.csv", delimiter=",", skip_header=1)
    assert np.array_equal(previous[:, 0], seconds[:len(previous)])
    np.testing.assert_allclose(curves[:len(previous)], previous[:, 1:], rtol=0, atol=1e-14)

    fonts = {f.name for f in font_manager.fontManager.ttflist}
    family = next(f for f in ("Microsoft YaHei", "SimHei", "SimSun") if f in fonts)
    plt.rcParams.update({"font.family": family, "font.size": 10,
                         "pdf.fonttype": 42, "axes.unicode_minus": False,
                         "axes.linewidth": .6})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.65), sharex=True)
    for j, ax in enumerate(axes):
        ax.axvspan(0, .5, color="#e8edf2", alpha=.8, linewidth=0)
        ax.axvline(.5, color="#777777", linewidth=.7, linestyle="--")
        ax.plot(seconds / 3600, curves[:, j], color="#0072B2", linewidth=1.0,
                label="分段线性插值", zorder=2)
        ax.scatter(time / 3600, data[:, j+1], s=3.5, facecolors="white",
                   edgecolors="#303030", linewidths=.45,
                   label="附件一采样值", zorder=3)
        ax.set_xlim(0, 4)
        ax.set_xticks(np.arange(0, 4.01, .5))
        ax.grid(axis="y", color="#e5e5e5", linewidth=.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=3, width=.6)
    axes[0].set_ylabel("温度（℃）")
    axes[0].set_ylim(26, 53)
    axes[0].set_yticks([28, 35, 42, 50])
    axes[0].text(.25, 51.6, "第一问", ha="center", va="top", fontsize=9, color="#555555")
    axes[0].text(.98, .08, "(a) 烘房温度", transform=axes[0].transAxes, fontsize=10, ha="right")
    axes[1].set_ylabel("水分浓度（kg/kg）")
    axes[1].set_ylim(.017, .054)
    axes[1].set_yticks([.02, .03, .04, .05])
    axes[1].text(.98, .08, "(b) 烘房水分浓度", transform=axes[1].transAxes, fontsize=10, ha="right")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc="upper center", ncol=2,
               frameon=False, bbox_to_anchor=(.5, 1.0), fontsize=9)
    fig.supxlabel("时间（h）", fontsize=10, y=.035)
    fig.subplots_adjust(left=.09, right=.99, top=.83, bottom=.20, wspace=.27)
    fig.savefig(root / "paper/figures/ambient_boundary.pdf", bbox_inches="tight", pad_inches=.04)
    fig.savefig(root / "figures/ambient_boundary.png", dpi=300, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)
    record = {"source": "../附件1.xlsx", "sample_count": len(time),
              "sampling_interval_s": 60, "range_s": [0, 14400],
              "interpolation": "piecewise linear (numpy.interp)",
              "plot_interval_s": 1, "all_original_points_shown": True,
              "layout": "side-by-side (1x2)",
              "smoothing": False, "extrapolation": False,
              "q1_boundary_max_absolute_difference": float(np.max(np.abs(curves[:len(previous)]-previous[:, 1:]))),
              "first_sample": data[0].tolist(), "last_sample": data[-1].tolist()}
    (root / "reports/AMBIENT_FIGURE_DATA.json").write_text(
        json.dumps(record, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(record, ensure_ascii=False))


if __name__ == "__main__":
    main()
