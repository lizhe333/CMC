# 2026 CUMCM A 题 Overleaf 工程

这是一套按“药材的烘干问题”四问结构定制的中文 LaTeX 工程。模型准备和第一问已按确认稿定稿，包含真实计算结果与数值检验；第二至四问及其他未完成部分继续保留红色草稿标识。

## 上传与编译

1. 直接把项目根目录的 `CUMCM_2026_A_Overleaf.zip` 上传到 Overleaf。
2. 将 Main document 设为 `main.tex`。
3. 将 Compiler 设为 XeLaTeX。仓库同时提供 `latexmkrc`，上传后通常会自动选用 XeLaTeX。
4. 本地在 paper 目录运行 `latexmk -xelatex main.tex`，自动重复编译至目录和交叉引用稳定。若直接用 XeLaTeX，首次需要重复编译，直到日志不再提示引用未定义或编号变化；不要把第一遍仍含 `??` 的 PDF 当作最终输出。

工程不依赖本机字体，使用 TeX Live 自带的 Fandol 中文字体，适合 Overleaf 直接编译。

## 文件分工

- `main.tex`：摘要、全局排版、章节入口与草稿开关。
- `sections/`：按题目逻辑拆分的正文，四问分别独立编辑。
- `tables/`：题面要求的温度、含水率结果表骨架。
- `figures/`：论文图片与跨问题模型框架图。
- `references.tex`：只加入已核验的真实文献。

## 协作约定

- 三位队员尽量分别编辑不同的 `sections/*.tex`，减少 Overleaf 合并冲突。
- 新图使用 `fig_q1_*`、`fig_q2_*`、`fig_q3_*`、`fig_q4_*` 命名；新表使用相同前缀。
- 论文关键数值必须能追溯到 `../results/`、`../reports/RESULTS_REPORT.md` 或代码输出。
- 草稿阶段保持 `main.tex` 中 `\draftmodetrue`；提交前改为 `\draftmodefalse`，并全文搜索“待填”“待确认”“TODO”。

## 2026-09-10 定稿更新

- 本次更新模型准备与第一问；参数用大括号表达，不再单列参数表。
- sections/model_preparation.tex：共用物理关系及初边值条件。
- sections/5_q1_preheating.tex：第一问分析、局部条件和完整模型。
- sections/q1_numerics_text.tex、q1_results_text.tex：求解检验与结果分析。
- 本地 code/q1_report.py --paper-only 从已保存结果更新上述两段文字，不运行 PDE。文字结构由 code/templates/ 下两个模板维护，避免重新生成旧稿。
- Overleaf 中可直接协作编辑 sections 下的正文；回同步本地时，将两段生成正文的文字调整同时合并到 code/templates/。
- 包内不需要 Python、Excel 或本机绝对路径；main.tex 使用 Fandol 字体，XeLaTeX 即可编译。
- 整篇仍处于协作阶段，不能以关闭草稿显示代替完成后续小问。
