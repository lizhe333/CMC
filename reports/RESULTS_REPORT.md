# 计算结果

## 运行环境

问题二使用 `C:\Users\32898\anaconda3\python.exe`、NumPy 2.1.3、SciPy 1.15.3、openpyxl 3.1.5。每问的完整参数和环境记录在对应 `results/q*/` 目录的验证 JSON；本文件只做全局索引。

## 数据读取与预处理

附件 1 的环境数据读取第一张表，检查时间严格递增、有限且覆盖请求区间。问题二在 0--10800 s 内对温度和有效水分浓度作分段线性插值，禁止超出原始 14400 s 数据端点外推。

## 问题一结果摘要

- 详细结果、检验及论文同步：[Q1_RESULTS_REPORT.md](Q1_RESULTS_REPORT.md)。
- 正式提交结果：`compition/2026/A题/A题/附件/附件3/result1.xlsx`；未舍入数组：`results/q1/q1_fields.npz`。
- 1800 s 时圆心、表面温度分别为 33.5753 °C、36.7856 °C，含水率分别为 2.5500、1.5102。

## 问题二结果摘要

- 详细结果：[Q2_RESULTS_REPORT.md](Q2_RESULTS_REPORT.md)；逐问建模与推算：[questions/Q2_MODELING.md](questions/Q2_MODELING.md)。
- 状态：`NUMERICAL_CHECKS_PASSED`；最终网格 $N=10240$，全时段 21 点相邻网格误差为 $e_T=3.426781\times10^{-6}$ °C、$e_C=3.731918\times10^{-6}$ kg/kg，论文表 60 个值四位小数稳定。
- 3 h 基准结果：中心 $(T,C)=(49.8494955852\ ^\circ\mathrm C,1.7661771011)$，表面 $(49.9664129357\ ^\circ\mathrm C,1.0081231229)$。
- 结果目录：`results/q2/`；验证主文件：`q2_validation.json`；参数/环境：`q2_parameters_environment.json`；完整场：`q2_fields.npz` 与完整 CSV；续算状态：`q2_restart.npz`；图形：`q2_fields.png`、`q2_fields.pdf`；灵敏度：`q2_sensitivity.csv`。
- 四联图径向剖面包含 0.5、1、1.5、2、2.5、3 h 六个题面时刻；Jacobian 证据为 `results/q2/diagnostics/jacobian_center_difference.json`。
- 正式工作簿：`compition/2026/A题/A题/附件/附件3/result2.xlsx`，温度和水分浓度两表各 10800 行（$t=1,\ldots,10800$ s），逐格回读与四位小数来源数组一致。

### 问题一与问题二的共同区间对照

本轮新增时空曲面图为 `results/q2/q2_spacetime.pdf` 与 PNG，由 `code/plot_q2_spacetime.py` 读取已验收场生成；采样数据保存在 `q2_spacetime_data.npz`，采样规则见 `q2_spacetime_metadata.json`。原四联图仍保留。

2026-09-11 按用户要求补充。`code/compare_q1_q2.py` 只读取两问已验收的未舍入场和逐秒边界，生成 `results/q2/q1_q2_comparison.json`。两问在 $0\sim1800\,\mathrm{s}$ 的初始状态、输出坐标与环境边界一致。

$1800\,\mathrm{s}$ 时，问题二中心和表面温度较问题一分别低 $1.3861\,{}^\circ\mathrm C$、$1.3725\,{}^\circ\mathrm C$；表面含水率高 $0.1384\,\mathrm{kg/kg}$。问题一初始体积热容量为 $2.132\times10^6\,\mathrm{J/(m^3\cdot K)}$，问题二为 $3.334695\times10^6\,\mathrm{J/(m^3\cdot K)}$；对应初始热扩散率分别为 $1.688555\times10^{-7}$ 与 $1.448282\times10^{-7}\,\mathrm{m^2/s}$。两问采用不同附录物性关系，该比较说明本构选择及状态反馈影响预测，不是耦合项的单因素消融，也不是对实验误差的度量。

读取 `q2_restart.npz` 确认：实际发布结果采用时间加密后的 $\mathrm{rtol}=10^{-10}$、$\mathrm{atol}=10^{-12}$、最大步长 $2.5\,\mathrm s$；$10^{-9}$、$10^{-11}$、$5\,\mathrm s$ 为此前空间加密比较的基准设置。

## 问题三结果摘要

- 详细模型与推算：[questions/Q3_MODELING.md](questions/Q3_MODELING.md)；验证主文件：`results/q3/q3_validation.json`。
- 状态：`NUMERICAL_CHECKS_PASSED`。4 h 后基准边界采用题设稳定工况 $T_\infty=50\,^\circ\mathrm C$、$C_\infty=0.05\,\mathrm{kg/kg}$，末一小时统计只作为平台证据和敏感性对照。
- 临界烘干时间为 $206900.064938\,\mathrm s=57.4722403\,\mathrm h$。临界点由圆心控制；临界点后 $0.5\,\mathrm s$ 的全网格最大含水率为 $0.1499998544\,\mathrm{kg/kg}$。
- 5120/10240 区间的临界时间差为 $0.638465\,\mathrm s$，10240 区间时间精化差为 $0.012902\,\mathrm s$。BDF/Radau 时间积分器对照和有效含水率累计通量检查均通过。
- 正式工作簿：`compition/2026/A题/A题/附件/附件3/result3.xlsx`，含 3448 个数据时刻和 21 个半径位置，时间为 $60,120,\ldots,206880\,\mathrm s$；逐格回读、四位小数和模板格式均通过。
- 未舍入场、精确事件、图、敏感性和审计证据位于 `results/q3/`。论文表 5 的结束行读取精确事件状态，不取最近整分钟替代。

## 灵敏度分析

问题二固定最终网格，对缺参延续假设的 $h$ 和 $h_m$ 分别执行 ±20% 单因素试验。完整表格和相对基准差值见 `results/q2/q2_sensitivity.csv`。敏感性只在同一输出变量内部比较相对变化：$h_m\pm20\%$ 使中心/表面含水率约变化 $+5.21\%/-4.20\%$、$+16.80\%/-13.34\%$，而 $h\pm20\%$ 对应约 $+0.59\%/-0.37\%$、$+0.15\%/-0.14\%$；不跨温度与含水率单位比较强弱。

问题三的时长敏感性使用 1000 区间同网格名义基准。末一小时均值相对 $50/0.05$ 只使终点延后 $6.08\,\mathrm s$，末次观测值使终点提前 $1090.60\,\mathrm s$；$h\pm20\%$ 使终点变化 $-46.92$ 至 $+70.96\,\mathrm s$，$h_m\pm20\%$ 使终点变化 $-3745.78$ 至 $+6112.30\,\mathrm s$。这些是模型情景比较，不是统计置信区间。

## 约束与一致性校验

问题二的 9 项针对性合同测试全部通过；测试覆盖附录 3 物性、均匀平衡、制造解、内部通量闭合、热湿交叉依赖、非法域与禁止外推、9000 s 图形时刻、Excel 数值及 G:V 样式/列宽继承、重启形状。原始模板副本已存在时仅验证未覆盖；正式文件由候选工作簿通过逐格数值与样式回读后更新。

问题二 BDF 内部使用可逆 $z=\log C$ 变量以保持严格正性，输出仍为物理 $C=\exp z$，无裁剪。交错 Jacobian 稀疏掩码覆盖 ±1、±2、±3 标量偏移，保留相邻节点的热湿交叉块。

历史 raw-$C$ 失败证据已归档至 `results/q2/diagnostics/archive/historical_raw_c/`，现行验证入口为 `q2_validation.json` 与 `diagnostics/q2_n*.json`。

### 问题二补充验证索引（2026-09-11）

补充实验与正式结果、Excel 和论文隔离，源快照保存在
results/q2_verification/source_snapshot/q2_solver.py。退化实验和有效含水率模型
积分平衡实验均为 PASS；后者是有效干基含水率变量的离散积分平衡检查，不表述为
湿物料总质量的严格守恒。
退化对照中 Q1 和 Q2 保留各自原始空间、边界和热表面离散；Q1 记录中的
“same spatial discretization”只表示 Q1 内部热方程与水分方程使用共同网格，
不表示两个程序的跨边界或热表面离散相同。

退化实验运行约 59.265 s，最终网格为 $N=1280$，未需要扩展到 $N=2560$：

| $N$ | 最大温度差（°C） | 最大含水率差（kg/kg） |
|---:|---:|---:|
| 320 | $1.938227\times10^{-6}$ | $1.186906\times10^{-9}$ |
| 640 | $4.815333\times10^{-7}$ | $1.495956\times10^{-9}$ |
| 1280 | $1.319033\times10^{-7}$ | $1.815060\times10^{-9}$ |

Q1 时间精化的最大温度和含水率差分别为
$2.715879\times10^{-8}\,^\circ\mathrm C$ 和 $1.955803\times10^{-9}$ kg/kg；
Q2 时间精化对应 $3.168729\times10^{-8}\,^\circ\mathrm C$ 和
$4.644156\times10^{-10}$ kg/kg。论文规定时刻和半径的四位小数子集保持一致。

积分平衡实验使用 $N=10240$ 的正式设置
($\mathrm{rtol}=10^{-10}$、$\mathrm{atol}=10^{-12}$、最大步长 2.5 s)及
时间精化设置 ($10^{-11}$、$10^{-13}$、1.25 s)。Gauss--Legendre 候选阶为
4、8、16、32，基准和精化均选 GL4--GL8，报告阶为 GL8。核心门禁值如下：

| 指标 | 基准 | 时间精化 |
|---|---:|---:|
| 逐秒节点最大 $|B|$（时刻） | $1.088507\times10^{-11}$（10560 s） | $1.403766\times10^{-12}$（7560 s） |
| accepted-node 最大 $|B|$（时刻） | $1.117773\times10^{-11}$（10560.0356 s） | $1.490363\times10^{-12}$（5760.0190 s） |
| 求积 floor | $8.881784\times10^{-16}$ | $8.881784\times10^{-16}$ |

积分门限 $10^{-7}$ kg/kg 和时间门禁均通过；平台阈值为 $10^{-9}$，基准/精化联合
最大残差分别为 $1.117773\times10^{-11}$ 和 $1.490363\times10^{-12}$。
逐秒梯形对照的最大 $|B|$ 为 $5.037692\times10^{-7}$ kg/kg；独立 dense 平均值
交叉核对最大差为 $8.881784\times10^{-15}$（阈值 $10^{-10}$）。

| 时刻 | 平均含水率（kg/kg） | GL8 $|B|$（kg/kg） | 逐秒梯形 $B$（kg/kg） |
|---:|---:|---:|---:|
| 0.5 h | 2.2826406262 | $2.311928\times10^{-12}$ | $5.032109\times10^{-7}$ |
| 1 h | 2.0649642254 | $2.948308\times10^{-12}$ | $5.036658\times10^{-7}$ |
| 2 h | 1.6900860834 | $5.765166\times10^{-12}$ | $5.037271\times10^{-7}$ |
| 3 h | 1.3825270550 | $1.049516\times10^{-11}$ | $5.037661\times10^{-7}$ |

补充文件索引：results/q2_verification/q2_degradation.json、q2_degradation.csv、
三个网格 NPZ、q2_balance.json、q2_balance.csv、q2_balance.npz、
q2_balance_time_refined.npz、q2_balance_accepted_nodes.csv、
q2_balance_segments.json，以及候选图 q2_degradation_comparison.png/pdf 和
q2_balance_residuals.png/pdf。平衡图改用对数纵轴并保留全时段及初始 60 s
局部视图，$t=0$ 零值不人为替换。27 项补充、Q2 与 Q3 合同测试全部通过。

## 与建模报告的一致性说明

问题二从 $t=0$ 重新初始化并全时段使用附录 3，不拼接问题一的 1800 s 状态；半径固定，不引入收缩、潜热或额外交叉源项。$h,h_m$ 的基准值是附录 2 缺参下的延续假设，已在结果报告中标明并做灵敏度。

## 可复现运行方式

```powershell
Set-Location D:\mathsmatical_modeling\compition\2026\CUMCM_2026_A
Set-Location .\code
& 'C:\Users\32898\anaconda3\python.exe' -m unittest test_q2_contracts -v
Set-Location ..
& 'C:\Users\32898\anaconda3\python.exe' .\code\q2_solver.py --max-n 10240 --skip-sensitivity
& 'C:\Users\32898\anaconda3\python.exe' .\code\run_q2_sensitivity.py
```

灵敏度在基准通过后调用 `q2_solver.sensitivity_runs`，结果写入 `results/q2/q2_sensitivity.csv`。正式论文数值必须继续从上述 JSON/CSV/图形数据溯源；本阶段不自动修改论文正文。


### 第三问时空分析图补充

新增含水率曲面和后期阈值等值图、代表时刻剖面和局部首次达标时间图。基于正式未舍入结果的分钟插值，表面约 12.55 h 达标，圆心约 57.47 h，二者相差约 44.92 h。局部时间为展示性插值估计，完整数值和来源见 results/q3/q3_analysis_figures.json 及 reports/questions/Q3_MODELING.md 第 10 节。未重新求解。
