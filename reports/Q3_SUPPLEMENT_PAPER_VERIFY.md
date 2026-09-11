# 第三问补充实验论文验收

## 结论

**PASS（仅针对第三问本次修改）**。该结论不表示问题四、摘要占位内容或整篇论文已经达到提交状态。

## 修改范围

- 原位重写 `paper/sections/7_q3_end_time.tex` 的结果分析和数值检验部分。
- 新增 `paper/tables/q3_convergence.tex`。
- 新增两幅论文图：长期边界敏感性和不同位置扩散系数变化。
- 删除第三问正文对最大含水率图、时空曲面和径向达标次序图的引用；原图文件保留为备查材料。
- 问题二正文、问题二结果表及表 7 与修改前逐字节一致。

## 数值一致性

| 检查项 | 结果 | 来源 |
|---|---|---|
| 基准终点 $57.4722\,\mathrm h$ | PASS | `results/q3/q3_validation.json` |
| 三级网格终点与相邻差异 | PASS | `results/q3/supplement/q3_supplement_grid_convergence.json` |
| 五种长期边界情景 | PASS | `results/q3/supplement/q3_supplement_boundary_sensitivity.json` |
| 终点扩散系数及比值 | PASS | `results/q3/supplement/q3_supplement_diffusivity.json` |
| 最大相对平衡残差 | PASS | `results/q3/supplement/q3_supplement_independent_audit.json` |

正文将 $57.4722\,\mathrm h$ 明确限定为 4 h 后环境保持 $50\,{}^\circ\mathrm C$、$0.05\,\mathrm{kg/kg}$ 时的条件预测；网格差异只用于证明数值稳定，不解释为实际预测误差。扩散系数差异写成与后期迁移受限相一致的机理证据，没有将表层写成唯一控制阻力。

## 图表与文字

- 第三问正文只保留两幅分析图及一张新增网格收敛表，另保留题目规定的表 7。
- 两幅图均可追溯到 `code/q3_supplement_experiments.py` 及同目录 CSV、JSON、NPZ 数据。
- 图题、坐标轴、单位和图例清晰，未发现裁切、重叠或误导性平滑拟合。
- 收敛表采用与问题二结果表一致的细网格样式，表头无异常换行。
- 正文未出现内部代理名称、临时目录、审计文件路径或版本管理信息。

## 编译与视觉检查

- 使用 `latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex` 编译成功，输出 29 页 A4 PDF。
- 日志无未定义引用、重复标签、Overfull、缺字或 LaTeX Error。仅保留第一问既有的两处 Underfull，本轮未改第一问。
- 已逐页检查第三问第 19--22 页：公式、表 7、收敛表、两幅图、图题、页码及问题四衔接均正常。
- `6verity` 的通用文本脚本仍会把章节内以 `figures/...` 编写的合法路径按章节目录解析，并将第一问二级 `\input` 识别为重复一级章节；这些是已知静态脚本误报。实际 XeLaTeX 编译、文件存在性检查和 PDF 视觉检查均已通过。

## 仍需处理的问题

- 长期径向近似的物理误差尚未量化。
- 问题四和摘要仍有既有占位内容，因此不能据本报告宣称全稿可提交。
