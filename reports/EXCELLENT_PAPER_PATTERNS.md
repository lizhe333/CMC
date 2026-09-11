# 优秀论文 Pattern Cards

本文件由 `mathmodel-paper-style-learning` 维护，用于记录每篇优秀论文可复用的结构、摘要、推导、图表、排版和验收模式。

## 使用规则

- 没有初稿时，只做学习和沉淀，不直接改写 `paper/`。
- 有初稿时，先诊断差距，再按用户确认进入优化。
- 同题范文可学习结果组织和逐问表达，但不得复制原文、原图、代码、数值或答案。
- 同类题型论文只迁移模型包装、图表选择和推导呈现方式。
- 跨题优秀论文只提取通用结构、摘要节奏、版式和验收标准。
- 新增优秀论文时，在下方追加 pattern card，不覆盖旧卡片。

## Batch <YYYY-MM-DD>: <paper set>

### <paper id or title>

- Source file/folder:
- Contest/problem:
- Current-problem relation: same-problem | same-type | cross-problem
- A/B boundary:
- Evidence quality:
- Use scope:
- Do-not-use boundary:

#### Structure Map

#### Abstract Pattern

#### Modeling Narrative

#### Derivation Presentation

#### Algorithm Stack and Complexity Control

- Algorithm stack:
- Role of each algorithm:
- Main-line algorithm vs auxiliary checks:
- How the paper keeps the story compact:
- Evidence that the complexity is necessary:
- Risk if imitated blindly:

#### Figure and Table Habits

#### Layout and Typesetting

#### Reusable Assets

#### Adoption Decision


## Batch 2026-09-10: 三篇 2018A 高温作业服论文

本批按全文文本层检索、正文页面总览及关键页原尺寸渲染核对。页码均为 PDF 文件页码。本题与范文同属非稳态传递及后续参数扩展问题，属于 same-type；同为 A 题不是迁移依据。仅学习表达和排版，原文中的服装分层、换热参数反演、算法及数值不移入药材模型。

### P18-1 高温作业服设计

- Source: exercise/优秀论文/2018A：高温作业服设计.pdf。
- Contest/problem: 2018A 高温作业服；31 页，附录从第 25 页开始。
- Evidence: 全文文本层；正文第 4—14 页总览，重点查看第 9、12 页。
- Use scope: 本次主要结构参考；不转用服装的分层关系或辐射假设。

#### Structure Map
第 4 页集中假设、符号；第 4—7 页模型准备；第一问从第 7 页至第 14 页上部，依次为问题分析、模型建立、求解、结果及检验、扩展；第二问从第 14 页开始，第三问从第 18 页开始。模型准备先铺基础，再为具体问题组织完整模型，后续问题明确使用前问框架。

#### Abstract Pattern
第 1 页先述研究对象，逐问交代方法与结果，再交代检验。当前只借鉴“各问目标—方法—答案”的节奏；第二至四问无结果，不补写其摘要结论。

#### Modeling Narrative
第 7—8 页先交代本问已知数据、未知量和后续用途；第 8 页分述控制方程与定解条件，第 9 页以“模型综合”成组呈现。当前借鉴完整初边值系统，但基础推导前移后不再逐项重复。

#### Derivation Presentation
第 8 页控制方程前仅一段守恒解释；第 9 页完整模型按物理角色标明层次。第 10—11 页仍有差分展开和稳定条件，并非全篇均少推导。当前采用短物理解释与成套方程，不迁移差分节点细节。

#### Algorithm Stack and Complexity Control
非稳态正向求解为主线，最小二乘参数搜索服务第一问，后续厚度优化复用正向模型；辐射扩展用于检验假设。角色各有用途。当前 Q1 参数已给定，不增加拟合或寻优；已有 BDF 求解与收敛检验已足够。

#### Figure and Table Habits
第 6、9、10 页的二维分层与离散示意服务变量定位，适合精确矢量绘制；当前如需示意，只画自身圆柱截面及边界。第 12 页拟合图后紧接误差和小表，下方并排两图；第 13 页列代表数据、完整数据放支撑文件。迁移“图表—证据—解释”邻接关系。当前题目规定的两张结果表必须完整保留，不能照范文用省略行。

#### Layout and Typesetting
A4；正文主要为 12 pt 宋体，一级标题有 14 pt 黑体，正文可见三级标题；第 8 页连续正文行距约 24.5 pt，明显比另两篇疏朗。三线表，表题在上、图题在下，公式居中编号靠右。第 4 页末存在标题与内容分离，第 9 页末也有标题悬置；这些不应模仿。

#### Reusable Assets
四节式小问框架；成套初边值模型；小型检验表；图后立即解释；按实际任务递进承接。

#### Adoption Decision
Adopt now：作为结构和结果组织主参考。Reject：通用传热知识综述、未用到的边界类型、逐节点展开、标题悬置和无关物理假设。

### P18-2 高温作业专用服装设计

- Source: exercise/优秀论文/2018A：高温作业专用服装设计.pdf。
- Contest/problem: 2018A 高温作业服；43 页，附录从第 29 页开始。
- Evidence: 全文文本层；正文第 4—21 页总览，重点查看第 18 页。
- Use scope: 模型层次、表图布局与局部检验；不迁移其两层到四层建模过程。

#### Structure Map
第 4 页进入建模，第 4—20 页主要为第一问，第二问从第 20 页开始，第三问从第 24 页开始，第 28 页模型评价及参考文献。没有独立的全局模型准备章节，基础知识在第一问内部展开。此结构不符合当前共用基础前置的目标。

#### Abstract Pattern
第 1 页用逐问段落组织模型、求解和结果。可借鉴平行组织，不复制原有句式或数字。

#### Modeling Narrative
第 4 页用一段说明本问目标和建模路线；第 7—13 页由两层耦合扩展四层耦合，连接条件单独说明。可学“局部关系组织为完整系统”，不在当前药材模型虚构材料分层。

#### Derivation Presentation
第 5—7 页基础定律与积分推导，第 13—18 页大量差分、矩阵、追赶法细节。推导比本次目标长得多；不据其优秀范文身份照搬篇幅。

#### Algorithm Stack and Complexity Control
隐式离散与追赶法承担正向求解，枚举搜索服务参数与厚度问题，结果与实测对比作检验。当前仅学习明确各方法用途，不用追赶法名称替换实际 BDF。

#### Figure and Table Habits
第 18 页表格采用蓝色表头、浅色交替底纹和细分隔线；第 19 页四幅位置曲线组成一组，随后给出解释与场分布图。需要数据为真实表格、状态场和观测值。当前保留自身四联图与规定结果表，表头可用克制浅蓝；不照搬深色满版表头、流程图或原图。示意属于二维精确制图，结果图由数据生成。

#### Layout and Typesetting
A4，正文主要 12 pt 宋体，正文部分连续行距约 16.3 pt，视觉比 P18-1 密集；文字可见左边界约 90 pt，比 P18-1 内缩。步骤标记较多，数学系统有占据大半页的情况。当前不通过挤压公式和正文复制这种密度。

#### Reusable Assets
紧凑表格、同主题子图组合、结果后就近给证据、将步骤合并为功能清晰的短段落。

#### Adoption Decision
Adopt now：表图组织。Reject：第一问占据十余页的比例、重复耦合推导、大量算法框图。Keep as candidate：未来真实优化若复杂再用流程图，不为装饰增加。

### P18-3 基于非稳态导热的高温作业专用服装设计

- Source: exercise/优秀论文/2018A：基于非稳态导热的高温作业专用服装设计.pdf。
- Contest/problem: 2018A 高温作业服；28 页，附录从第 19 页开始。
- Evidence: 全文文本层；正文第 4—14 页总览，重点查看第 11 页。
- Use scope: “基础—具体模型—求解—解释”的段落关系与紧凑排版；不将其细杆模型当作药材径向模型。

#### Structure Map
第 4 页符号后进入模型建立与求解，先有物理背景，第 5—6 页导热推导；第 6—8 页第一问模型，第 8—12 页求解与结果；第二问从第 12 页开始。基础前置于第一问具体模型，但位于建模大章节内部。

#### Abstract Pattern
第 1 页逐问写出任务、传热模型或优化工具、结果及解释。当前取其逐问聚焦，未完成小问的成果仍留空。

#### Modeling Narrative
先以简化几何和物理机制确定方程，再区分几何、初始、衔接和边界条件，最后求解。当前借鉴定解条件的清楚分组；药材为径向圆柱几何，不能移用沿细杆轴向的一维方程。

#### Derivation Presentation
第 5—6 页采用细杆热量收支推导；第 8—10 页先讨论多种差分再给最终方法，叙述仍较长。当前只保留薄壳守恒三步，不复述候选算法历程，也不迁移其物理结论。

#### Algorithm Stack and Complexity Control
主线最终采用后向欧拉和三对角求解，后续用优选搜索厚度；算法选择有其算力和稳定性背景。当前方法由真实代码决定，BDF 不改名为后向欧拉。

#### Figure and Table Habits
第 11 页先场分布曲面，再给若干时刻的二维剖面，并随后解释时间和空间差异。借鉴两种观察角度及图后解释；当前已有径向剖面与时间曲线，无需为了外观添加三维曲面。早期物理照片不服务当前定量推导，舍弃。图源必须是当前 results，适合程序绘图。

#### Layout and Typesetting
正文主要 12 pt 宋体，部分说明 10.5 pt，连续正文行距常见约 15.6 pt；小标题比 P18-1 更弱，图片占据较大页面比例。公式用多层复合编号，未统一迁移，当前保持 LaTeX 章节公式自动编号。

#### Reusable Assets
段内标签代替过多小节；时间与空间结果相互解释；正文末段自然承接下一问。

#### Adoption Decision
Adopt now：段内层次及结果解释节奏。Reject：长篇算法概念、细杆推导形式、原文“均匀介质温度梯度相同”等未经本题验证的断言、冗长复合编号。
