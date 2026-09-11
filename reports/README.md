# 报告目录

四问共同入口：[全局建模与数值求解注意事项](GLOBAL_MODELING_NOTES.md)。各问建模、代码审计和结果发布前读取；具体执行证据写入逐问报告和验证文件。

阶段产物由对应 skill 生成或覆盖：

- `PROBLEM_UNDERSTANDING.md`：`1start-mathmodel` 的阶段0题意理解记录。
- `ANALYSIS_MODELING_REPORT.md`：`2analysis-modeling` 的候选模型、主线总纲、逐问报告索引和跨问题一致性检查。
- `MODEL_AUDIT_REPORT.md`：`2-5model-audit` 的模型二次审稿记录，检查正确性、可行性、新颖性和全面性。
- `SYMBOL_STYLE_GUIDE.md`：本题全局符号、命名风格和逐问继承规则。
- `EXCELLENT_PAPER_PATTERNS.md`：`mathmodel-paper-style-learning` 的单篇优秀论文 pattern card 与追加记录。
- `PAPER_STYLE_GUIDE.md`：从同题、同类题和跨题优秀论文沉淀出的本题写作、推导、排版和摘要风格指南。
- `PAPER_FIGURE_PLAN.md`：写作前的图表规划，记录每张图表的表达目的、数据来源要求和参考风格。
- `questions/`：逐问报告目录，每一问使用独立的 `Q*_MODELING.md`。
- `RESULTS_REPORT.md`：`3coding-visual` 的结果记录。
- `DRAWIO_REPORT.md`：`4drawio` 的非数据图记录。
- `VERIFY_REPORT.md`：`6verity` 的验收记录。

通用格式约束：

- Markdown 报告中的数学符号和公式必须使用 `$...$` 或 `$$...$$`。
- 不要用标记为 `text` 的代码块或反引号代码样式代替数学公式。
- 反引号只用于真实文件名、命令、代码变量或数据字段。
- 训练报告中的等式统一用普通等号 `=`，不要用 `:=`。
- 优先使用直观符号；涉及对时间的变化率时，优先写成 `$\frac{d\theta(t)}{dt}$` 这类形式，不默认使用点号导数等简写。
- 公式要写成能渲染的完整表达，关系式必须有清楚的等号或不等号。
- 每一问都必须写清推算过程，不能只给最终公式或最终结果。
- 阶段2确定主线后，必须先更新 `SYMBOL_STYLE_GUIDE.md`；后续每一问只新增少量局部符号，并继承全局符号。
