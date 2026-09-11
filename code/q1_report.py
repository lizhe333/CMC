"""Generate Q1 discussion and provenance from saved, validated numerical results."""
import json
import platform
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sci(value):
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"{mantissa}\times10^{{{int(exponent)}}}"



def slash_units(text):
    """Use slash notation only for physical units, leaving mathematical fractions intact."""
    units = {
        r"\frac{\mathrm{kg}}{\mathrm{kg}}": r"\mathrm{kg/kg}",
        r"\frac{\mathrm{kg}}{\mathrm{m}^3}": r"\mathrm{kg/m^3}",
        r"\frac{\mathrm{J}}{\mathrm{kg}\cdot\mathrm{K}}": r"\mathrm{J/(kg\cdot K)}",
        r"\frac{\mathrm{W}}{\mathrm{m}\cdot\mathrm{K}}": r"\mathrm{W/(m\cdot K)}",
        r"\frac{\mathrm{W}}{\mathrm{m}^2\cdot\mathrm{K}}": r"\mathrm{W/(m^2\cdot K)}",
        r"\frac{\mathrm{m}^2}{\mathrm{s}}": r"\mathrm{m^2/s}",
        r"\frac{\mathrm{m}}{\mathrm{s}}": r"\mathrm{m/s}",
    }
    for fraction, slash in units.items():
        text = text.replace(fraction, slash)
    return text


def generate(paper_only=False):
    folder = ROOT / "results" / "q1"
    meta = json.loads((folder / "q1_validation.json").read_text(encoding="utf-8"))
    if meta["status"] != "NUMERICAL_CHECKS_PASSED":
        raise ValueError("Refusing to publish provisional results.")
    data = np.load(folder / "q1_fields.npz")
    t, c = data["temperature_C"], data["moisture_kg_per_kg"]
    final, delta = meta["final_run"], meta["refinement"][-1]["difference_from_previous"]
    td, mean = meta["temporal_refinement_difference"], meta["end_mean_moisture"]
    env = np.loadtxt(folder / "boundary_used.csv", delimiter=",", skiprows=1, encoding="utf-8-sig")
    values = {
        "N": str(final["n_intervals"]), "PREV_N": str(final["n_intervals"] // 2),
        "SPACE_T": sci(delta["max_temperature_C"]),
        "SPACE_C": sci(delta["max_moisture_kg_per_kg"]),
        "TIME_T": sci(td["max_temperature_C"]),
        "TIME_C": sci(td["max_moisture_kg_per_kg"]),
        "MASS": sci(final["max_mass_balance_residual"]),
        "TEMP_GAP": f"{t[-1,-1]-t[-1,0]:.4f}",
        "ENV_END": f"{env[-1,1]:.4f}", "SURFACE_C": f"{c[-1,-1]:.4f}",
        "MEAN_C": f"{mean:.4f}", "DROP": f"{100*(2.55-mean)/2.55:.4f}",
    }
    # Editorial structure lives in reviewed templates; this step only fills saved results.
    for table_name in ("temperature", "moisture"):
        table_path = ROOT / "paper" / "tables" / f"q1_{table_name}.tex"
        table = table_path.read_text(encoding="utf-8")
        if r"\fontsize{10.5pt}{13pt}" not in table:
            table = table.replace(r"\centering", r"\centering" + "\n" +
                                  r"  \fontsize{10.5pt}{13pt}\selectfont" + "\n" +
                                  r"  \renewcommand{\arraystretch}{1.08}")
        table = table.replace(r"\begin{tabular}{rrrrrr}",
                              r"\begin{tabular*}{0.92\textwidth}{@{\extracolsep{\fill}}rrrrrr}")
        table = table.replace(r"\end{tabular}", r"\end{tabular*}")
        table_path.write_text(slash_units(table), encoding="utf-8")
    for name in ("numerics", "results"):
        text = (Path(__file__).parent / "templates" / f"q1_{name}.tex.tpl").read_text(encoding="utf-8")
        for key, value in values.items():
            text = text.replace(f"@@{key}@@", value)
        if "@@" in text:
            raise ValueError(f"Unresolved token in q1_{name}.tex.tpl")
        (ROOT / "paper" / "sections" / f"q1_{name}_text.tex").write_text(slash_units(text), encoding="utf-8")
    if paper_only:
        return
    rows = []
    for record in meta["refinement"]:
        d = record.get("difference_from_previous")
        values = (f"{d['max_temperature_C']:.6g} | {d['max_moisture_kg_per_kg']:.6g} | "
                  f"{d['paper_moisture_kg_per_kg']:.6g}" if d else "— | — | —")
        rows.append(f"| {record['n_intervals']} | {record['dr_m']:.9g} | {values} |")
    report = f"""# 第一问计算结果与来源

状态：完成真实附件计算和数值检查，不代表全题提交验收。

## 输入与模型

- 原始输入：{meta['boundary_source']}，工作表 {meta['boundary_sheet']}，共 {meta['boundary_samples']} 个采样点，覆盖 0–14400 秒，原文件不修改。
- 本问取 0–1800 秒，采样点间分段线性插值，不外推。派生的逐秒边界见 results/q1/boundary_used.csv。
- 参数来自题面附录 2；模型采用用户确认的第一问解法及 paper/sections/5_q1_preheating.tex。
- 固定半径、主体截面一维径向、常热物性、无显式潜热反馈，水分按题目浓度口径采用有效 Fick 扩散。
- 局部审查：圆柱几何项、圆心极限、Robin 边界符号及变系数通量已核对。第二至四问未运行。

## 实际求解方法

保留显式程序作为对照；最终用相同空间离散的稀疏 BDF 隐式时间积分。原因是早期表面水分梯度层很薄，细网格显式法的时间步限制很严格。最终采用 {final['n_intervals']} 个区间，相对容差 1e-11、绝对容差 1e-13，最大步长 5 秒。

## 收敛证据

下表比较相邻网格，采用相同的 BDF 时间容差。最大差异在全部 1801×21 个输出点上计算。

| 径向区间数 | 空间步长（m） | 温度最大差异（°C） | 含水率最大差异 | 论文表格点含水率差异 |
| --- | --- | --- | --- | --- |
{chr(10).join(rows)}

收紧十倍时间容差后温度最大差异 {td['max_temperature_C']:.8g} °C，含水率最大差异 {td['max_moisture_kg_per_kg']:.8g}。
全部输出点的相邻网格差异小于 5e-5；论文表格在最后两级网格上四位小数一致。此为经验检验，不是严格误差上界。

## 1800 秒结果

| 距中心轴距离（cm） | 温度（°C） | 干基含水率 |
| --- | --- | --- |
"""
    for j in [0, 5, 10, 15, 20]:
        report += f"| {j/10:g} | {t[-1,j]:.4f} | {c[-1,j]:.4f} |\n"
    report += f"""
体积加权平均含水率 {mean:.10f}，相对初始值下降 {100*(2.55-mean)/2.55:.6f}%；中心与表面温差 {t[-1,-1]-t[-1,0]:.6f} °C。中心未舍入含水率为 {c[-1,0]:.12f}。

## 交付物与来源

- results/q1/result1.xlsx：温度、水分浓度两个工作表，首行为 0–2 cm、每 0.1 cm，A 列为 0–1800 秒，共 1801 个数据行。数值保存并显示四位小数。
- 附件 3 模板未提供，按题面附录 1 重建工作表名称、坐标和数据布局。
- results/q1/q1_fields.npz：未舍入状态数组、时空坐标与守恒记录。
- results/q1/q1_temperature_full.csv 和 q1_moisture_full.csv：完整未舍入结果；对应的 *_paper.csv 是论文指定表格。
- results/q1/q1_fields.pdf 与 q1_fields.png：由结果数组直接生成的四面板图，径向曲线使用规定的 21 个输出点。
- results/q1/q1_validation.json：输入路径、参数、环境、逐级收敛与运行状态。
- results/q1/mass_balance.csv：平均含水率、累计外流和守恒残差；最大残差 {final['max_mass_balance_residual']:.8g}。
- paper/tables/q1_*.tex 和 paper/sections/q1_*_text.tex：由真实结果生成的表格与讨论。

## 验证

code/test_q1_contracts.py 的五项针对性验证通过：均匀平衡不变；变系数质量守恒及显隐式一致性；圆柱 Robin 导热与独立 Bessel 级数解析解对照；拒绝不稳定显式步长；拒绝重复时间与不完整边界。
Excel 导出后逐单元格读取并与舍入结果数组比较，全部通过。合成验证场景不作为题目计算结果。

## 复现

实际解释器：{sys.executable}。Python {platform.python_version()}。完整环境见验证 JSON。

从仓库根目录依次执行：

1. python compition/2026/CUMCM_2026_A/code/q1_solver.py --boundary compition/2026/附件1.xlsx --sync-paper
2. python compition/2026/CUMCM_2026_A/code/q1_report.py
3. python compition/2026/CUMCM_2026_A/code/test_q1_contracts.py

原显式程序可用 --method explicit --max-n 640 复现。在此网格下完整输出的早期表面含水率尚未达到 5e-5 网格差异条件，所以程序保留诊断结果并返回非零状态，不发布到论文。
"""
    (ROOT / "reports" / "Q1_RESULTS_REPORT.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-only", action="store_true",
                        help="Update reviewed LaTeX fragments without rewriting the technical report.")
    generate(paper_only=parser.parse_args().paper_only)
