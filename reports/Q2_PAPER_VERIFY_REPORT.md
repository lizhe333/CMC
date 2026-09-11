# 问题二论文同步专项验证报告

## 结论

**Q2 范围 PASS；整篇论文仍为协作草稿。**

本轮只同步问题二正文、问题二表图、摘要中的问题二段落、模型准备中的跨问接口和模型评价中的问题二结论。问题一、问题三、问题四的实质模型与未完成占位不在本轮验收范围内。

## 资产生成与来源

资产由以下只读结果驱动脚本生成，脚本不导入求解器、不重跑 PDE：

```powershell
Set-Location D:\mathsmatical_modeling\compition\2026\CUMCM_2026_A
& 'C:\Users\32898\anaconda3\python.exe' .\code\generate_q2_paper_assets.py
```

脚本读取 `results/q2/q2_temperature_paper.csv`、`q2_moisture_paper.csv` 和 `q2_validation.json`，验证六个时刻、五个半径位置、30 个温度值和 30 个含水率值，然后写出：

- `paper/tables/q2_temperature.tex`
- `paper/tables/q2_moisture.tex`
- `paper/figures/q2_fields.pdf`

两张表均为 10.5 pt、约 `0.92\textwidth`，单位采用斜杠写法，数值统一四位小数。表和图直接使用同一份已验收的未舍入结果来源。

## Q2 内容门禁

| 检查项 | 结果 |
| --- | --- |
| 四个二级小节 | `问题分析`、`模型建立`、`模型求解`、`结果分析`，PASS |
| 初值与边界 | 明确从 $t=0$ 重算；附件 1 分段线性边界；无 $1800\,\mathrm{s}$ 人为切换，PASS |
| 缺参口径 | $h,h_m$ 明确为附录 2 的延续假设，并报告 ±20% 同变量相对敏感性，PASS |
| 物理边界 | 明确 $T+273.15$、固定半径、无潜热/收缩/额外源项，PASS |
| 正性参数化 | 明确 $z=\log C$ 仅用于 BDF 内部，输出仍为 $C$，不裁剪，PASS |
| 离散与 BDF | 圆柱控制体、保守界面通量、Robin 表面通量、联立 BDF 与交叉 Jacobian，PASS |
| 表格 | 0.5--3 h 六时刻 × 0--2 cm 五位置，共 60 个值，PASS |
| 图形 | 四联图含六个径向时刻和五个位置时间线，PASS |
| Q2 占位 | Q2 正文、Q2 两表、摘要 Q2 段和评价 Q2 项无 `resultblank`/`draftnote`，PASS |

正文采用真实结果解释温度梯度衰减、含水率中心--表面差异、$D(C,T)$ 的温度/含水率竞争作用和 $h,h_m$ 敏感性；未跨摄氏度与 kg/kg 比较影响强弱。前三小时工作簿的题意边界已作为问题三续算接口保留。

## 数值溯源核对

`results/q2/q2_validation.json` 状态为 `NUMERICAL_CHECKS_PASSED`。论文中使用的关键值均来自当前结果报告和论文 CSV：

- $3\,\mathrm{h}$ 中心/表面温度：$49.8495/49.9664\,{}^\circ\mathrm{C}$；
- $3\,\mathrm{h}$ 中心/表面含水率：$1.7662/1.0081\,\mathrm{kg/kg}$；
- $N=5120\to10240$ 的最大差：$e_T=3.426781\times10^{-6}\,{}^\circ\mathrm{C}$、$e_C=3.731918\times10^{-6}\,\mathrm{kg/kg}$；
- 时间收紧最大差：$3.751870\times10^{-6}\,{}^\circ\mathrm{C}$、$4.548273\times10^{-8}\,\mathrm{kg/kg}$；
- 60 个论文值在空间/时间加密后四位小数稳定。

## LaTeX 编译

编译入口为 `paper/main.tex`，运行：

```powershell
Set-Location D:\mathsmatical_modeling\compition\2026\CUMCM_2026_A\paper
latexmk -xelatex -interaction=nonstopmode main.tex
```

结果：`paper/main.pdf`，25 页，A4（595.28 × 841.89 pt），PDF 非空，交叉引用稳定。Q2 在第 13--17 页；第 17 页为 Q2 结果分析的前三段，随后进入问题三。Q2 对应日志未发现 undefined reference、Q2 图片缺失或 Q2 overfull hbox。Q2 五页已逐页视觉检查：标题、公式编号、两张表、图题、图例和页边界均无裁切或重叠。

上述三个审计修订后于 2026-09-11 重新运行 `latexmk -xelatex -interaction=nonstopmode main.tex`；编译完成且目标为最新，Q2 相关日志仍无 `undefined`、`Overfull`、`error` 或缺图记录。

## 6verity 检查说明

已运行 `MathModelAgent/skills/6verity/scripts/writing_check.sh`。该脚本按 section 文件目录解析 `figures/...`，而本 LaTeX 工程按 `main.tex` 目录解析同一相对路径，因此直接运行会把现有的 `figures/q*.pdf` 报为路径误报；实际 `latexmk` 已成功加载 Q2、Q3、Q1 图形。脚本同时提示未引用的备用模型图和无 citation marker，这些不属于 Q2 数值同步错误。

为排除上述相对路径差异对正文门禁的影响，曾在 `paper/sections/figures/` 建立临时只读副本后再次运行同一脚本；本次 `writing text gate passed`。复核完成后已删除临时副本。该次运行仍仅保留非 Q2 的结构性提示：Q1 辅助文本未被主入口直接纳入、假设节偏短、备用流程图未引用以及当前工程没有 citation marker。

Q2 专项补充检查确认：

1. Q2 正文和 Q2 表格不存在占位符或内部 `results/`、`reports/` 路径泄露；
2. `paper/figures/q2_fields.pdf` 存在且图题明确六个规定时刻和五个位置；
3. 表格行列数和四位小数由生成器按 CSV 合同验证；
4. PDF 文本层含表 5、表 6 的六个时刻及 Q2 关键结果。

## 非 Q2 遗留项

`main.tex` 仍保持草稿开关，摘要中的问题一、问题三、问题四以及模型评价中 Q1/Q3/Q4 仍有待定内容；这些未在本轮填充。问题三、问题四的长期延拓、终点和收缩对照应由相应小问完成后再更新，不应由 Q2 同步阶段代填。
