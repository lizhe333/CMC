# 问题二用户集中修订记录

日期：2026-09-11。授权：用户逐项提出意见后明确要求“开始修改”。范围为第二问论文正文、摘要对应段、两张规定表及表格生成器；第一问仅用作写作参照和数值对照来源。

## 八项修改与验收

| 用户要求 | 对应修改 | 验收对象 |
| --- | --- | --- |
| 运算分式与式（7.1）对齐 | 单位保留斜杠，数学除法使用分式；四条物性公式等号对齐 | 正文源码和最终 PDF |
| 指数统一写为 e 的幂 | 行内与行间统一指数形式 | Q2 源码无指数函数命令 |
| 重要内容黑体加粗 | 关键方法与结果判断使用局部黑体 | 最终 PDF |
| 全文论文语言优化 | 按 5writing 重组现象、机理、证据，移除交接式描述 | 正文和图后解释 |
| 模型求解参照 Q1 | 连贯介绍空间离散、时间推进和精度检验 | 与 q1_numerics_text.tex 的组织方式对照 |
| 保留题面规定表格 | 网格线、两层合并表头、原字段与顺序，填入四位小数 | 两表及生成器，逐值对照 CSV |
| 与 Q1 对比印证 | 同初边条件下对照 1800 s 的结果，解释本构差异与状态反馈 | q1_q2_comparison.json 与两问原始场 |
| 增加时空变化分析图 | 温度与含水率并排三维曲面，补充整体演化分析 | q2_spacetime.pdf 及采样数据、绘图脚本 |

## 新增结论的来源

比较脚本为 `code/compare_q1_q2.py`，输入为两问已验收的 NPZ 和逐秒边界 CSV，输出为 `results/q2/q1_q2_comparison.json`。脚本校验共同时间、半径、初始状态及边界一致，记录两问初始物性和共同区间差值。正式表值仍来自 `q2_temperature_paper.csv`、`q2_moisture_paper.csv`；原四联图来自 `q2_fields.pdf`，新增时空图来自 `q2_spacetime.pdf`。

比较不将本构差异与耦合反馈混同，不把模型间差值宣称为实验误差。最终发布的 Q2 时间积分设置依据 `q2_restart.npz`，为相对容差 $10^{-10}$、绝对容差 $10^{-12}$、最大步长 $2.5\,\mathrm s$。

新增时空图采用 `academic-plotting` 的数据绘图路径，由 `code/plot_q2_spacetime.py` 生成 PDF、300 dpi PNG、`q2_spacetime_data.npz` 和元数据 JSON。图中使用前一分钟逐秒、此后逐分钟的原始状态值和全部 21 个径向节点，共 $240\times21$ 个采样点；不拟合、不改变模型结果。曲面单元仅连接相邻采样点。

## 验收状态

本轮第二问修订验收通过。最终 PDF 为 `paper/main.pdf`（27 页），第二问位于第 14--18 页；两张规定表位于第 16 页，原四联图位于第 17 页，新增时空曲面图为第 18 页图 5。

`code/verify_q2_paper_revision.py` 独立检查了两表的 60 个值与原始 NPZ、结果 CSV 在四位小数下逐值一致；检查了原题全网格和合并表头；核对时空图数据与原数组抽样逐元素完全一致；确认论文嵌入的两张 PDF 与结果目录源图逐字节一致。记录见 `results/q2/q2_paper_revision_validation.json`。本轮不重算 PDE，不改官方 Excel。

### 编译与目视核验

- `latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex` 退出码为 0，所有构建目标已更新，引用稳定；无未定义引用、无 Overfull。两条 Underfull 均来自现有第一问源码第 89--91、124--128 行，本轮未处理第一问。
- 主 agent 实际渲染检查第 14--18 页及相邻的第 13、19 页、摘要第 1 页；公式（7.1）和初边值式等号对齐，表头与网格完整，公式、图表及注释未越出版心。通过局部 `Needspace` 修复初边值引导句、结果分析标题和下一问标题的孤立分页。
- 时空图曲面采用 300 dpi 栅格层、文字及坐标轴保留矢量，消除纯矢量曲面在 PDF 阅读器中形成的细纹伪影。原始状态值、曲面采样和颜色范围均未因此改变。
- 按 `5writing` 去除交接文档式表述，以现象、定量结果和物理机理组织正文；主体修订由现有 Luna/max 执行，独立科学与数值对照审计由现有 Sol/high 执行，主 agent 完成最终数值复核及版面修订。

### 6verity 范围与工具限制

本轮执行 `6verity` 的写作检查及第二问专项人工检查。通用 `writing_check.sh` 的返回码为 1，不能记为工具自动通过：其重复章节提示来自第一问第 1 行被注释的旧 `section`；其缺图提示将 `figures/...` 错按 `sections/` 解析，而实际 XeLaTeX 从论文入口目录解析，五个相关图片均真实存在且编译成功。其第一问嵌套 `input` 警告同样不是遗漏。未为绕过这些报告而修改第一问或检查脚本。

第二问范围的结构、符号、单位、公式、60 个规定表值、图源、结论证据、引用稳定性与页面检查均通过。全稿摘要中第一问、第四问和总验证段的既有占位符仍保留，故本记录不宣称整篇论文已可直接提交。先前 Q2 论文核验报告中的页码与三线表描述以本记录为准。

### 复现

在项目目录下依次运行（绘图环境需 NumPy、Matplotlib）：

```powershell
python code/compare_q1_q2.py
python code/plot_q2_spacetime.py
python code/generate_q2_paper_assets.py
python code/verify_q2_paper_revision.py
Set-Location paper
latexmk -xelatex -interaction=nonstopmode -halt-on-error main.tex
```
