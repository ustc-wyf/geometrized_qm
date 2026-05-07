# 163-三切片guarded-plus-only初始包验证

日期：2026-05-07

## 一句话结论

`no-active-edge + auto_bad_zero` 不是只在 `tau=-3.5` 偶然有效：它在 `tau=-3.5,0,+3.5` 三个关键切片上都能生成可推进的 D/gBCD `plus-only` 初始加速度包；三者一步物质重构都没有负质量壳判别式，\(\rho^\tilde\) 拉回到 A 表象的一步加权偏差均在 `~3e-5` 到 `~1e-4` 量级。

## 本轮要回答的问题

上一轮只在左侧分离态 `tau=-3.5` 上验证了 guarded `plus-only` 初始包。为了判断该方法是否具有普适性，本轮把同一流程扩展到：

- `tau=-3.5`：左侧分离态，上一轮最难切片；
- `tau=0`：干涉中心；
- `tau=+3.5`：右侧分离态。

流程保持一致：

1. 用 `metric_active_dilation=0` 重新解 `plus-only` \(g_+\)。
2. 导出 `plus-only` 初始包。
3. 使用 `auto_bad_zero` 局部可容许性守卫。
4. 做一步 D 物质推进，检查质量壳判别式和测度偏差。

## 汇总输出

- 汇总 JSON：
  - `visualizations/equation_first_gbcd_guarded_three_tau_summary/summary.json`
- 线性坐标汇总图：
  - `visualizations/equation_first_gbcd_guarded_three_tau_summary/guarded_three_tau_summary_linear.png`
- 对数坐标汇总图：
  - `visualizations/equation_first_gbcd_guarded_three_tau_summary/guarded_three_tau_summary_log.png`

图中变量含义：

- `D equation residual`：把 guarded 初值包代回
  \[
  \tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2=\mathcal C_{\mu\nu}
  \]
  后得到的相对残差，加权平均越小越好。
- `one-step measure deviation`：从该初始包推进一步后，在 \(g_+\) 上重构出的 \(\sqrt{|\tilde g|}\rho^\tilde\) 与 A/KG 下一切片对应测度的相对 L1 偏差。这里 A 只作比较参考，不作为演化约束。
- `guard frozen points`：`auto_bad_zero` 冻结 \(\delta g_+\) 的支撑区格点数。它表示局部不可容许点数量，不是物理源项。

## 三切片结果

### `tau=-3.5`

输出：

- `visualizations/equation_first_gbcd_plus_initial_package_n384_taum3p5_core10_noactiveedge_guarded/`
- `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_taum3p5_core10_noactiveedge_guarded/`

结果：

- fit points：`2550`
- support points：`6034`
- D 方程 residual weighted mean：`0.13220`
- residual p95：`0.57969`
- frozen support points：`2`
- frozen core10 points：`2`
- 初始 \(\rho^\tilde\) 拉回偏差：`1.40e-17`
- 一步 \(g_+\) 测度偏差：`1.091e-4`
- 一步 \(g_+\) 加权测度偏差：`1.043e-4`
- 负判别式比例：`0`
- 判别式最小值：`2.75e-7`

解释：

- 这是三者中残差最高的切片，仍是后续求解器优化重点。
- 但它已经从“旧 active-edge 一步测度偏差 `~1.4e4`”改善到 `~1e-4`，且不再有负判别式。

### `tau=0`

输出：

- `visualizations/equation_first_gbcd_plus_initial_package_n384_tau0_core10_noactiveedge_guarded/`
- `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_tau0_core10_noactiveedge_guarded/`

结果：

- fit points：`1024`
- support points：`4242`
- D 方程 residual weighted mean：`6.60e-5`
- residual p95：`1.56e-4`
- frozen support points：`0`
- 初始 \(\rho^\tilde\) 拉回偏差：`1.37e-17`
- 一步 \(g_+\) 测度偏差：`4.994e-4`
- 一步 \(g_+\) 加权测度偏差：`2.999e-5`
- 负判别式比例：`0`
- 判别式最小值：`9.28e-8`

解释：

- 干涉中心是当前最干净的切片。
- D 方程残差已经是 `1e-4` 以下量级，且不需要任何守卫冻结点。
- 未加权测度偏差比加权偏差大，说明误差更多来自支撑区边缘，而不是主概率权重区。

### `tau=+3.5`

输出：

- `visualizations/equation_first_gbcd_plus_initial_package_n384_taup3p5_core10_noactiveedge_guarded/`
- `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_taup3p5_core10_noactiveedge_guarded/`

结果：

- fit points：`3589`
- support points：`6211`
- D 方程 residual weighted mean：`0.01324`
- residual p95：`0.05676`
- frozen support points：`2`
- frozen core10 points：`0`
- 初始 \(\rho^\tilde\) 拉回偏差：`1.32e-17`
- 一步 \(g_+\) 测度偏差：`2.464e-4`
- 一步 \(g_+\) 加权测度偏差：`2.884e-5`
- 负判别式比例：`0`
- 判别式最小值：`9.90e-8`

解释：

- 右侧分离态残差比左侧分离态低约一个数量级。
- 守卫冻结 2 个支撑点，但不在 core10 主核心区。
- 一步物质重构稳定。

## 物理判断

这轮结果支持以下判断：

- `plus-only` 作为 D 初始加速度求解器是可行的初步路线。
- A-derived 初始切片不需要通过大幅改 \(\rho^\tilde\) 或 \(g_0\) 来满足 D 方程；保持 \(\tilde g_0,\rho^\tilde,u_i\) 后，由 D 方程决定 \(g_+\) 更自然。
- `auto_bad_zero` 处理的是分支/近退化点的数值可容许性，而不是改变物理方程；每次都必须报告冻结点数量和 residual 代价。
- 当前真正难点不再是“是否能构造一个 D 初始包”，而是“如何把这个 guarded 初始加速度求解推广成多步全演化，并把守卫从后处理升级成强约束或 patch-boundary 条件”。

## 下一步

1. 把三切片 guarded package 作为当前标准初值验证集。
2. 实现多步推进原型：每一步都进行 plus-only metric solve、matter step、可容许性守卫、D residual 记录。
3. 将 `auto_bad_zero` 从后处理升级为解线性系统时的硬约束或活动集 active-set 规则。
4. 对 `tau=-3.5` 继续优化：检查为何左侧分离态 residual `0.132` 远高于 `tau=0` 和 `tau=+3.5`，特别是是否与近退化 patch 分布、辅助张量代表元或时间非对称有关。
