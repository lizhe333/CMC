# 问题四建模报告：移动边界下的热湿耦合模型

## 1. 问题目标与建模边界

问题四要求在问题二、问题三的热湿耦合框架上，进一步考虑药材因失水而产生的径向尺寸变化，并根据附录 4 的经验物性重新计算烘干时长。这里有两个同时发生的变化：一是药材的计算区域随半径变化，二是密度、比热容、导热系数和水分扩散系数整体切换为附录 4。因此，问题四的时长与问题三的差异必须通过控制变量分解，不能全部解释为收缩效应。

药材仍近似为长度 $L=0.25\,\mathrm{m}$ 的轴对称圆柱，只保留径向传热和传质。初始半径为 $R_0=0.02\,\mathrm{m}$，初始温度和干基含水率分别为

$$
T(r,0)=28\,^{\circ}\mathrm{C},
\qquad
C(r,0)=2.55\,\mathrm{kg/kg}.
$$

烘干结束时刻定义为药材内部含水率首次达到临界值的时刻。由于题目要求“各处”低于 $0.15\,\mathrm{kg/kg}$，判据取全域最大含水率，而不是径向平均含水率：

$$
g(t)=\max_{0\le r\le R(t)}C(r,t)-0.15.
$$

临界时刻取 $g(t)$ 首次向下穿过零的时刻；等号时刻作为临界根，其后严格小于阈值的状态作为题意中“低于”的满足证据。

## 2. 数据、物性与共同参数

### 2.1 输入数据

问题四使用以下已核对的数据口径。

| 数据 | 采用方式 | 单位与范围 |
| --- | --- | --- |
| 附件 1 | 在原始观测范围内对烘房温度和环境水分浓度分别作分段线性插值 | 时间为 $\mathrm{s}$，温度为 $^{\circ}\mathrm{C}$，水分浓度为 $\mathrm{kg/kg}$ |
| 附件 2 | 读取 145 条半径记录，使用同一个半径插值函数同时得到 $R(t)$ 和 $\dot R(t)$ | 时间为 $\mathrm{s}$，半径由 $\mathrm{cm}$ 转为 $\mathrm{m}$；覆盖 $0$--$259200\,\mathrm{s}$ |
| 附录 4 | 在每一个空间节点以当前 $C$ 和 $T$ 计算局部物性 | $T$ 的指数变量使用开尔文温度 |
| 附件 3 | 作为 `result4.xlsx` 的结果模板 | 时间列为 $\mathrm{s}$，距离列为 $\mathrm{cm}$ |

附件 2 的半径从 $2.000\,\mathrm{cm}$ 单调下降至 $1.198\,\mathrm{cm}$，末端记录保持 $1.198\,\mathrm{cm}$。半径的单调性、时间递增性和数据有限性在求解前检查；方程内部统一使用秒和米，避免将厘米数据直接代入 SI 制参数。

### 2.2 附录 4 的局部物性

对归一化坐标中的状态 $\widetilde C(\xi,t)$ 和摄氏温度 $\widetilde T(\xi,t)$，定义

$$
\begin{aligned}
\rho(\widetilde C)&=760+90\widetilde C,\\
c_p(\widetilde C)&=1850+2150\frac{\widetilde C}{\widetilde C+1},\\
k(\widetilde C)&=0.12+0.20\frac{\widetilde C}{\widetilde C+1},\\
D(\widetilde C,\widetilde T)&=4.2\times10^{-4}
\exp\left(-\frac{0.30}{\widetilde C}\right)
\exp\left(-\frac{3850}{T_K}\right),
\qquad T_K=\widetilde T+273.15.
\end{aligned}
$$

其中 $\rho$ 的单位为 $\mathrm{kg/m^3}$，$c_p$ 的单位为 $\mathrm{J/(kg\cdot K)}$，$k$ 的单位为 $\mathrm{W/(m\cdot K)}$，$D$ 的单位为 $\mathrm{m^2/s}$。程序只在 $\widetilde C>0$ 且 $T_K>0$ 的物性定义域中评价这些函数；内部对数变量只用于保持正性，不改变物理输出的含义。

题目给出的 $h=25\,\mathrm{W/(m^2\cdot K)}$ 和 $h_m=8\times10^{-7}\,\mathrm{m/s}$ 沿用问题二、问题三的参数。问题四没有另行给出交换系数，因此该取值是缺少新参数时的延续假设，需通过 $h$ 和 $h_m$ 的单因素扰动说明其影响。

### 2.3 物理近似及其范围

1. 药材长度和轴向形状保持不变，附件 2 的表面半径代表整个截面的半径；不凭半径数据外推轴向收缩。
2. 内部材料点作均匀径向收缩，内部运动规律由材料坐标假设给出，而不是附件 2 直接观测得到的事实。
3. $C$ 是干基含水率，即水质量与干物质质量之比，不是单位体积水质量；纯压缩本身不使 $C$ 按体积比例放大。
4. 不加入潜热反馈、独立的蒸发源项或未经题面支持的力学方程。温度对水分迁移的影响仅通过附录 4 中的 $D(C,T)$ 体现。
5. $\rho(C)$ 仅作为有效热学系数中的经验密度使用。本模型不同时要求它、给定半径变化和固定干物质质量满足严格的湿物料密度闭合关系。

这些近似把附件 2 的尺寸序列作为外生几何输入，因而本问的结论是所选有效热湿模型下的烘干时长，而不是包含收缩力学和完整焓守恒的三维物理模型。

## 3. 移动材料坐标与连续模型

### 3.1 坐标移动和材料速度

实际径向坐标为 $r$，随时间变化的材料区域为 $0\le r\le R(t)$。采用归一化材料坐标

$$
\xi=\frac{r}{R(t)},
\qquad
0\le\xi\le1,
$$

并把物理场拉回固定区间：

$$
\widetilde T(\xi,t)=T(R(t)\xi,t),
\qquad
\widetilde C(\xi,t)=C(R(t)\xi,t).
$$

均匀径向收缩意味着固定材料标签 $\xi$ 在时刻 $t$ 的位置为 $r(t)=R(t)\xi$，所以材料速度为

$$
v(r,t)=\frac{\dot R(t)}{R(t)}r.
$$

对任意随体标量场 $F(r,t)$，令 $\widetilde F(\xi,t)=F(R(t)\xi,t)$。链式法则给出

$$
\left.\frac{\partial F}{\partial t}\right|_r
=
\left.\frac{\partial\widetilde F}{\partial t}\right|_\xi
-\frac{\dot R(t)}{R(t)}\xi
\frac{\partial\widetilde F}{\partial\xi},
\qquad
\frac{\partial F}{\partial r}
=
\frac{1}{R(t)}\frac{\partial\widetilde F}{\partial\xi}.
$$

因此，材料导数满足

$$
\frac{\partial F}{\partial t}+v\frac{\partial F}{\partial r}
=
\frac{\partial\widetilde F}{\partial t}.
$$

这一步决定了收缩项的处理方式：若在物理坐标中将热湿演化写成含材料导数的有效方程，材料速度与网格速度恰好相同，$\dot R$ 产生的坐标项与输运项抵消。因而映射后的方程中不再另加一项形式上的“收缩对流项”。收缩仍通过 $R(t)^{-2}$ 的扩散尺度因子、$R(t)^{-1}$ 的表面边界尺度因子以及输出坐标映射进入模型。

### 3.2 变换后的热湿方程

在物理坐标中，采用与问题二、问题三相同的有效热湿方程，并将时间导数理解为材料导数：

$$
\rho(C)c_p(C)
\left(\frac{\partial T}{\partial t}+v\frac{\partial T}{\partial r}\right)
=
\frac{1}{r}\frac{\partial}{\partial r}
\left(rk(C)\frac{\partial T}{\partial r}\right),
$$

$$
\frac{\partial C}{\partial t}+v\frac{\partial C}{\partial r}
=
\frac{1}{r}\frac{\partial}{\partial r}
\left(rD(C,T_K)\frac{\partial C}{\partial r}\right).
$$

利用 $r=R(t)\xi$ 及上一节的材料导数关系，得到固定区间 $0<\xi<1$ 上的方程

$$
\rho(\widetilde C)c_p(\widetilde C)
\frac{\partial\widetilde T}{\partial t}
=
\frac{1}{R(t)^2\xi}
\frac{\partial}{\partial\xi}
\left(\xi k(\widetilde C)
\frac{\partial\widetilde T}{\partial\xi}\right),
$$

$$
\frac{\partial\widetilde C}{\partial t}
=
\frac{1}{R(t)^2\xi}
\frac{\partial}{\partial\xi}
\left(\xi D(\widetilde C,\widetilde T_K)
\frac{\partial\widetilde C}{\partial\xi}\right),
\qquad
\widetilde T_K=\widetilde T+273.15.
$$

式中的 $R(t)^{-2}$ 来自两次径向导数变换，不能只把问题三程序中的固定半径 $R_0$ 替换为 $R(t)$ 而省略坐标运动和空间算子推导。

### 3.3 初始条件与边界条件

圆心处使用轴对称极限，满足

$$
\left.\frac{\partial\widetilde T}{\partial\xi}\right|_{\xi=0}=0,
\qquad
\left.\frac{\partial\widetilde C}{\partial\xi}\right|_{\xi=0}=0.
$$

表面状态记为 $\widetilde T_s(t)=\widetilde T(1,t)$ 和 $\widetilde C_s(t)=\widetilde C(1,t)$。由物理坐标中的交换条件

$$
-k(C_s)T_r(R(t),t)=h\,[T_s(t)-T_\infty(t)],
$$

$$
-D(C_s,T_{s,K})C_r(R(t),t)=h_m\,[C_s(t)-C_\infty(t)],
$$

变换后得到

$$
-\frac{k(\widetilde C_s)}{R(t)}
\left.\frac{\partial\widetilde T}{\partial\xi}\right|_{\xi=1}
=h\,[\widetilde T_s(t)-T_\infty(t)],
$$

$$
-\frac{D(\widetilde C_s,\widetilde T_{s,K})}{R(t)}
\left.\frac{\partial\widetilde C}{\partial\xi}\right|_{\xi=1}
=h_m\,[\widetilde C_s(t)-C_\infty(t)].
$$

初始条件为

$$
\widetilde T(\xi,0)=28,
\qquad
\widetilde C(\xi,0)=2.55,
\qquad 0\le\xi\le1.
$$

边界右端的正负号由材料外法向约定。若表面温度低于环境温度，热交换使热量进入药材；若表面含水率高于环境有效含水率，传质通量指向环境。

## 4. 干基含水率口径与通量闭合

### 4.1 含水率的守恒含义

干基含水率 $C$ 的定义是

$$
C=\frac{\text{水质量}}{\text{干物质质量}}.
$$

因此，随体材料块在水和干物质均未交换的纯压缩过程中，$C$ 不因体积缩小而改变。若要建立真实水质量，需要先给出当前体积中的干物质密度 $\rho_d$，再写成

$$
M_w(t)=\int_{\Omega(t)}\rho_d(\boldsymbol{x},t)
C(\boldsymbol{x},t)\,\mathrm dV.
$$

直接计算 $\int_{\Omega(t)}C\,\mathrm dV$ 只能得到带体积量纲的积分，不能称为总水质量。若把经验式 $\rho(C)$ 解释为湿药材体积密度，还需另行令 $\rho_d=\rho/(1+C)$，并检验其与 $R(t)$、固定长度及干物质守恒是否相容。本问不作这一未经数据支持的解释。

### 4.2 所选近似下的干基平均量

在初始干物质沿半径均匀分布、材料标签由 $\xi$ 保持的近似下，用参考截面的干物质基准定义归一化平均含水率

$$
\overline C_d(t)=2\int_0^1\widetilde C(\xi,t)\,\xi\,\mathrm d\xi.
$$

将水分方程乘以 $2\xi$ 后积分，得到

$$
\begin{aligned}
\frac{\mathrm d\overline C_d}{\mathrm dt}
&=2\int_0^1\frac{\partial\widetilde C}{\partial t}\xi\,\mathrm d\xi\\
&=\frac{2}{R(t)^2}\int_0^1
\frac{\partial}{\partial\xi}
\left(\xi D\frac{\partial\widetilde C}{\partial\xi}\right)\mathrm d\xi\\
&=\frac{2}{R(t)^2}
\left[\xi D\frac{\partial\widetilde C}{\partial\xi}\right]_{0}^{1}\\
&=-\frac{2h_m}{R(t)}
\left[\widetilde C_s(t)-C_\infty(t)\right].
\end{aligned}
$$

最后一步使用圆心无通量条件及表面边界条件

$$
D_s\left.\frac{\partial\widetilde C}{\partial\xi}\right|_{\xi=1}
=-R(t)h_m[\widetilde C_s(t)-C_\infty(t)].
$$

上述等式是所选有效扩散方程的干基通量闭合检查：它检验内部通量在离散后是否只通过表面边界离开参考控制域。由于模型没有同时求解干物质密度输运、严格湿物料质量关系和收缩力学，它不是完整水质量守恒定律，也不能由此宣称总能量守恒。

## 5. 半径和环境函数

### 5.1 半径插值

在附件 2 的观测区间 $[0,259200]$ 内，对数据点 $(t_j,R_j)$ 在米制单位下构造保形三次插值

$$
R(t)=\operatorname{PCHIP}\{(t_j,R_j)\},
\qquad
\dot R(t)=\frac{\mathrm dR(t)}{\mathrm dt},
$$

并从同一插值函数求导。PCHIP 的优点是保持单调收缩数据的形状，避免普通高次多项式的区间过冲；实现时显式关闭区间外外推。分段线性插值作为几何敏感性对照，不能与 PCHIP 的半径和导数混用。

若临界时刻未超过 $72\,\mathrm h$，只使用实测区间内的 PCHIP。若临界时刻超过 $72\,\mathrm h$，主方案在 $t_{72}=259200\,\mathrm{s}$ 后保持

$$
R(t)=R_{72}=0.01198\,\mathrm m,
\qquad
\dot R(t)=0.
$$

这一平台是对附件 2 之外时间的延拓假设，不是新的观测。若终点落在该延拓区间，另以“24 小时时间尺度再渐近收缩 1%”作情景对照：

$$
R_{\mathrm{slow}}(t)=R_{72}
\left[1-0.01\left(1-\exp\left(-\frac{t-t_{72}}{86400}\right)\right)\right],
\qquad t\ge t_{72},
$$

其导数由同一表达式解析求得。该情景只用于量化延拓不确定性，不用于重新拟合附件 2。

### 5.2 烘房边界

在附件 1 的原始覆盖区间内，对 $T_\infty(t)$ 和 $C_\infty(t)$ 分别作分段线性插值，不允许插值器自动外推。附件 1 末端后的主情景沿用问题三已审阅的长期边界：从 $4\,\mathrm h$ 起取

$$
T_\infty(t)=50\,^{\circ}\mathrm C,
\qquad
C_\infty(t)=0.05\,\mathrm{kg/kg}.
$$

该平台由末段观测的稳定性支持，但仍属于缺少长期观测时的环境延拓，末段均值和末次观测值只作为边界敏感性情景，不替代名义平台。

## 6. 四组控制变量对照

为了分离物性变化与几何收缩的作用，四组计算从同一初态出发，使用相同的环境函数、$h$、$h_m$、终止阈值和数值设置。

| 组别 | 物性 | 几何 | 比较目的 |
| --- | --- | --- | --- |
| A | 附录 3 | 固定半径 $R_0$ | 复现问题三的固定半径基准 |
| B | 附录 3 | 附件 2 的 $R(t)$ | 观察附录 3 物性下的收缩影响 |
| C | 附录 4 | 固定半径 $R_0$ | 分离物性公式变化 |
| D | 附录 4 | 附件 2 的 $R(t)$ | 问题四正式模型 |

记四组首次达标时间为 $t_A,t_B,t_C,t_D$。报告以下差值：

$$
\Delta_{\mathrm{property}}=t_C-t_A,
\qquad
\Delta_{\mathrm{shrink},4}=t_D-t_C,
$$

以及收缩效应在两套物性之间的交互差

$$
\Delta_{\mathrm{interaction}}
=\left(t_D-t_C\right)-\left(t_B-t_A\right).
$$

直接比较 $t_D-t_A$ 只表示两项改变的综合效果。若某一差值小于当前时长的数值分辨能力，应报告为未能分辨，而不是将其解释为确定的物理效应。

## 7. 数值离散与计算接口

### 7.1 归一化圆柱有限体积离散

在 $[0,1]$ 上取 $N$ 个区间，节点为 $\xi_i=i\Delta\xi$，其中 $\Delta\xi=1/N$；控制面为 $\xi_{i+1/2}$。节点 $i$ 的无量纲圆柱权重取为

$$
w_i=\frac{\xi_{i,+}^2-\xi_{i,-}^2}{2},
$$

其中 $\xi_{i,-}=\max(0,(i-1/2)\Delta\xi)$、$\xi_{i,+}=\min(1,(i+1/2)\Delta\xi)$。于是

$$
\sum_{i=0}^{N}w_i=\frac12,
$$

圆心半控制体和表面半控制体均被保留。

节点间的热、湿通量因子分别为

$$
\Phi^T_{i+1/2}=\xi_{i+1/2}k_{i+1/2}
\frac{\widetilde T_{i+1}-\widetilde T_i}{\Delta\xi},
$$

$$
\Phi^C_{i+1/2}=\xi_{i+1/2}D_{i+1/2}
\frac{\widetilde C_{i+1}-\widetilde C_i}{\Delta\xi},
$$

其中 $k_{i+1/2}$ 和 $D_{i+1/2}$ 按相邻节点局部物性的算术平均取值，并始终保留在通量散度内。圆心通量为零；表面通量由 Robin 条件直接给出

$$
\Phi^T_{N+1/2}=-R(t)h\,[\widetilde T_N-T_\infty(t)],
\qquad
\Phi^C_{N+1/2}=-R(t)h_m\,[\widetilde C_N-C_\infty(t)].
$$

因此，热方程和水分方程的半离散形式分别为

$$
\rho_i c_{p,i}w_i\frac{\mathrm d\widetilde T_i}{\mathrm dt}
=\frac{\Phi^T_{i+1/2}-\Phi^T_{i-1/2}}{R(t)^2},
$$

$$
w_i\frac{\mathrm d\widetilde C_i}{\mathrm dt}
=\frac{\Phi^C_{i+1/2}-\Phi^C_{i-1/2}}{R(t)^2}.
$$

这里 $\Phi_{-1/2}=0$，而表面半控制体使用上式给出的交换通量。该写法同时保留了内部通量抵消、$R^{-2}$ 空间尺度和表面 $R^{-1}$ 边界效应。

### 7.2 联立时间推进和正性变量

温度直接以摄氏度推进，含水率用内部变量

$$
z_i=\log\widetilde C_i,
\qquad
\widetilde C_i=\exp(z_i)
$$

推进；水分右端在变量变换后满足

$$
\frac{\mathrm dz_i}{\mathrm dt}
=\frac{1}{\widetilde C_i}
\frac{\mathrm d\widetilde C_i}{\mathrm dt}.
$$

采用联立 BDF 隐式积分，记录每个分块积分的接受步状态、步中点状态、逐时半径和全网格重启状态。解析 Jacobian 覆盖相邻节点的热湿交叉块以及 $z=\log C$ 的商法则；$R(t)$ 是外生输入，不作为状态变量求导。不得以浓度裁剪、重置物理时间或把外生半径当作待估状态来掩盖迭代失败。

### 7.3 可编码的输入与输出

求解器输入包括物性版本、半径插值方式、环境延拓方式、$h$、$h_m$、初值、网格数、BDF 容差和最大步长。输出至少包括：

- 未舍入的全网格 $\widetilde T(\xi_i,t)$、$\widetilde C(\xi_i,t)$、$R(t)$ 及接受步时间；
- 终止根、根前后状态、终点包围区间、通量闭合量和网格/时间收敛记录；
- 四组对照的时长、场量和统一的运行配置；
- 结果表需要的物理距离重构数据及图表数据。

## 8. 终点定位与结果映射

### 8.1 全域事件判定

计算中在全部径向节点、所有接受步端点、步中点和规定输出时刻监测

$$
g_N(t)=\max_{0\le i\le N}\widetilde C_i(t)-0.15.
$$

只在确认全程的中心含水率等于离散全域最大值后，才允许用中心值作为快速事件函数；否则始终使用全网格最大值。求得候选跨越后，在同一个稠密数值解上连续求根，保存宽度不超过 $1\,\mathrm{s}$ 的包围区间 $[t_-,t_+]$、根残差

$$
|g_N(t_d)|\le10^{-9}\,\mathrm{kg/kg},
$$

以及严格穿越证据

$$
g_N(t_d-0.5\,\mathrm{s})>0,
\qquad
g_N(t_d+0.5\,\mathrm{s})<0.
$$

同时用接受步、中点和步长精化检查单个大步内是否可能发生多次阈值穿越。若在最长 $14$ 天计算窗内未找到事件，按“该时间窗未达标”记录，不能为了符合背景描述而调整模型参数。

### 8.2 固定物理距离的输出

结果表的固定距离列为

$$
r_j=0,0.1,0.2,\ldots,2.0\,\mathrm{cm}.
$$

给定时刻先用当前半径计算

$$
\xi_j(t)=\frac{r_j}{R(t)}.
$$

当 $0\le r_j\le R(t)$ 时，在未舍入的 $\widetilde C(\xi,t)$ 上按 $\xi$ 作分段线性插值；当 $r_j>R(t)$ 时，该位置已在药材区域外，工作簿相应单元留空，不填零、环境值或表面值。另设独立的“药材表面”列，始终取 $\widetilde C(1,t)$，不能用最近的固定距离列替代。

工作簿按 $60\,\mathrm{s}$ 采样，并在末尾追加唯一的精确终点行；若终点恰为整分钟则与该分钟行合并。论文表 6 给出每隔 $6\,\mathrm h$ 的状态和精确终点状态，固定距离仍按 $0.5\,\mathrm{cm}$ 展开，并保留实时半径对应的表面含水率。所有展示值最后才按题面要求保留四位小数。

## 9. 进入代码前的审计接口

下表区分已由题面或数据确定的内容、需要明确声明的假设和必须由程序验证的内容，供模型二次审稿与数值审计使用。

| 类别 | 内容 |
| --- | --- |
| 题面/数据确定 | 初始状态、长度、阈值、附件 1 与附件 2 的字段和单位、附录 4 四项公式、结果文件的采样距离和时间要求 |
| 建模假设 | 均匀径向收缩、长度保持、$C$ 的干基含义、$\rho(C)$ 仅作有效热学系数、无潜热、附件 1 末端后的环境平台、附件 2 末端后的半径平台 |
| 需要代码核验 | $R^{-2}$ 内部算子、$R^{-1}$ Robin 边界、圆心极限、表面半控制体、内部通量抵消、正性、解析 Jacobian、全域最大值事件和物理距离映射 |

至少执行以下针对性检查后，才将结果用于论文：

1. 取已知源项的移动域制造解，验证空间离散中的 $R^{-2}$、边界中的 $R^{-1}$ 以及网格加密收敛；另检验无交换且初态均匀时纯收缩不改变 $C$。
2. 固定半径、切换附录 3 并采用相同网格和时间设置，复现问题三的完整规定输出和终点，确认新增模块没有改变基础方程。
3. 在多个半径时刻和非均匀状态上用中心差分核对解析 Jacobian，覆盖热湿交叉块及 $z=\log C$ 的商法则；检查不裁剪浓度、不把 $R$ 作为求导状态。
4. 以 $640,1280,2560,\ldots$ 个区间逐级加密至最多 $20480$ 个区间，并在固定网格上收紧 BDF 容差、减小最大步长。终点差应同时满足绝对值不超过 $30\,\mathrm{s}$、相对值不超过 $10^{-4}$；共同物理距离输出的温度和含水率误差分别不超过 $5\times10^{-5}$。
5. 通过连续根、宽度不超过 $1\,\mathrm{s}$ 的包围区间、根残差和根前后严格符号共同确认终点，不能以最后一个整分钟采样点代替事件时刻。
6. 以独立通量实现和代表性区间的 Radau 结果作交叉核对；在 D 组中分别扰动 $h,h_m$ $\pm20\%$，比较 PCHIP 与分段线性半径，并检查问题三长期环境边界情景。

当前文档只确定模型和可编码接口，不预先写入问题四的数值结论。所有时长、表格数值和论文结论均应在上述检查闭合后，从同一组未舍入结果生成。
