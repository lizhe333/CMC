# 方案

项目：CUMCM_2026_A
创建时间：2026-09-10 19:04:52

## 用户偏好

- 排版引擎：LaTeX
- 竞赛类型：国赛
- 论文语言：中文
- 子问题数量：4
- 工作模式：半自动辅助，阶段产出后人工审稿

## 推进原则

- 四问共同遵循 [全局建模与数值求解注意事项](reports/GLOBAL_MODELING_NOTES.md)。建模、实现和验收分别核对；注意事项不是已通过的验证记录。
- 先定总纲：确认题型、各问依赖关系、共用数据处理、共同假设、主线模型族和论文叙事线。
- 再逐问推进：每一问都完成“局部题意复述 -> 推算过程 -> 局部模型细化 -> 代码实现 -> 结果解读 -> 论文小节同步”。
- 每问独立成文：每一问都写入 `reports/questions/Q*_MODELING.md`，全局报告只保留总纲、索引和跨问题一致性。
- 模型二次审稿：每一问进入代码或论文前，使用 `2-5model-audit` 检查正确性、可行性、新颖性和全面性。
- 允许回修前问：后续问题带来更深理解时，可以回到前一问修正假设、指标或实现，但必须记录原因。

## 正式比赛目标

- 稳定完成四问，保持“固定区域预热模型 -> 变物性全过程模型 -> 达标时刻判定 -> 收缩移动边界模型”的递进主线。
- 优先保证方程、边界条件、单位、结果模板和数值来源一致，再讨论亮点扩展。
- 论文工程使用中文 XeLaTeX，保持可直接上传 Overleaf 协作。

## workflow

| step | skill | 产物 | 人工审稿点 |
| --- | --- | --- | --- |
| 0 | `1start-mathmodel` | `reports/PROBLEM_UNDERSTANDING.md`, `plan.md`, `todo.md` | 确认题意、子问题、附件和关键歧义 |
| 1 | `2analysis-modeling` | `reports/ANALYSIS_MODELING_REPORT.md` 初稿 | 建模前确认 2-3 套候选模型 |
| 2 | `2analysis-modeling` | `reports/ANALYSIS_MODELING_REPORT.md` 总纲版，`reports/questions/Q*_MODELING.md` 骨架 | 确认主线模型族、共用假设、各问依赖、通用变量、指标口径和逐问顺序 |
| 2.5 | `2-5model-audit` | `reports/MODEL_AUDIT_REPORT.md` | 二次确认模型正确、可行、新颖、全面 |
| 3 | `3coding-visual` | `code/`, `results/`, `figures/`, `reports/RESULTS_REPORT.md`, `reports/questions/Q*_MODELING.md` | 每一问写代码前确认局部推算过程、变量、目标函数、约束、评价指标和输出表 |
| 4 | `mathmodel-tikz-figures` / `4drawio` | `figures/tikz/*.tex`, `figures/*.pdf`, `reports/TIKZ_FIGURE_REPORT.md` 或 `reports/DRAWIO_REPORT.md` | 只制作确实服务模型表达的图 |
| 5 | `5writing` | `paper/` | 写论文前确认主模型、创新点、图表清单和章节结构 |
| 6 | `6verity` | `reports/VERIFY_REPORT.md`, `review.md`, `score.md`, `failure_cases.md` | 提交前只修验收报告指出的小问题，并复盘 workflow / skill 改进 |

## 数值和数据规则

- 不编造数据、参考文献或结果。
- 论文关键数值必须来自 `reports/RESULTS_REPORT.md`、`results/`、代码输出或图表数据。
- 每次修改模型或结果后，必须同步更新结果报告和论文。
- 每一问的推算过程、局部模型、结果解释优先写入对应的 `reports/questions/Q*_MODELING.md`。
- 问题四必须区分尺寸收缩与附录 4 物性公式变化的影响，不能把问题三、四的结果差全部解释为收缩效应。

