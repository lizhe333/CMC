# code

当前实现范围：第一问。第二至四问尚未实现。

## 第一问复现

在仓库根目录运行：

```powershell
python -m pip install -r compition/2026/CUMCM_2026_A/code/requirements-q1.txt
python compition/2026/CUMCM_2026_A/code/q1_solver.py --boundary compition/2026/附件1.xlsx --sync-paper
python compition/2026/CUMCM_2026_A/code/q1_report.py
python compition/2026/CUMCM_2026_A/code/test_q1_contracts.py
```

本机已用 `C:\Users\32898\anaconda3\python.exe` 完成运行，不需要重复安装已有依赖。

- `q1_solver.py`：读取原始附件、显式基准算法、稀疏 BDF 求解、空间收敛与时间容差检查、结果导出及论文表格同步。
- `q1_report.py`：由已验证的结果生成第一问讨论和来源报告。
- `test_q1_contracts.py`：独立解析解、质量守恒、平衡状态、步长和输入有效性检查；合成测试数据不是题目答案。
- 默认 BDF 求解逐级加密至全部规定输出点满足差异阈值，再收紧时间容差。`--method explicit --max-n 640` 可复现显式算法的加密过程。
- `--template <原始result1.xlsx路径>` 用于保留官方模板格式；当前按题面附录 1 重建布局。
- 正式结果在 `../results/q1/`；精度检查失败时保存在 `diagnostics/` 子目录，不覆盖已有正式结果。
- `--sync-paper` 只写第一问的两张表和图；章节讨论由报告脚本更新。

内部使用米、秒、摄氏温度及题目干基含水率口径。输出距离为厘米，结果保留四位小数，未舍入数组另存 NPZ 和 CSV。

第一问论文表达约定：数学分式统一用 `\frac{分子}{分母}`，热流率用 `\mathrm{heat}`，不用带点的 Q。

## 仅更新论文正文

已存在验证结果时，运行 `python compition/2026/CUMCM_2026_A/code/q1_report.py --paper-only`，只同步 LaTeX 检验、结果文字及两张表的版式，不运行 PDE、不改 Excel 或技术结果报告。文字模板位于 `code/templates/`，修改正文时应同步模板，避免再生成时覆盖。
