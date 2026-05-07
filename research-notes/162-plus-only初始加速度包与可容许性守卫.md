# 162-plus-only初始加速度包与可容许性守卫

日期：2026-05-07

## 一句话结论

`plus-only` 路线已经从“残差诊断”推进到一个可读取的 D 初始加速度包：在 `tau=-3.5,n=384,core10` 左侧分离态上，`no-active-edge` 重新求解加 `auto_bad_zero` 局部可容许性守卫后，D 方程完整回代残差约 `0.132`，一步物质重构无负质量壳判别式，\(\rho^\tilde\) 拉回到 A 表象后的相对偏差约 `1.1e-4`。

## 本轮要回答的问题

上一轮确认 `center_plus + source` 不健康，而 `plus-only` 是真正合理的初始加速度自由度。但充分迭代的旧 `plus-only` 解有一个新问题：

- D 方程残差很好：exact weighted mean `0.13917`。
- 但把一步物质状态重构到 \(g_+\) 上时，约 `1.9%` 支撑点出现负质量壳判别式，测度偏差被放大到 `~1.4e4`。

本轮问题是：

- 这是 D 支物理模型失败，还是 metric update 在边界/退化点上没有施加可容许性约束？
- 能否在不改 D 方程的前提下，得到一个既满足 D 方程、又能推进物质的一步初始包？

## 新增与修改脚本

- `kg_examples/export_gbcd_plus_initial_package.py`
  - 新增 `--plus-guard auto_bad_zero`。
  - 导出可被后续演化器读取的 \(g_-,g_0,g_+\)、\(\rho^\tilde\)、\(u_\mu\)、守恒密度、辅助张量 \(\mathcal C_{\mu\nu}\) 和 D 方程残差。
  - `auto_bad_zero` 的含义：先用 candidate \(g_+\) 做一步物质重构；凡是支撑区内出现负质量壳判别式或错误行列式分支的格点，把该点的 \(\delta g_+\) 冻结为 0，即回到 A 参考未来切片。它是局部可容许性/信赖域处理，不是修改 D 方程。

- `kg_examples/test_gbcd_plus_initial_package_one_step.py`
  - 读取导出的初始包。
  - 只推进一小步 D 物质变量。
  - 同时比较在 \(g_0\) 和 \(g_+\) 上重构出的测度与 A/KG 下一切片测度的差异。

- `kg_examples/scan_gbcd_plus_alpha_admissibility.py`
  - 扫描全局 \(\alpha\delta g_+\)。
  - 用来验证“整体阻尼”能否同时保留 D 残差改善和物质可容许性。

- `kg_examples/scan_gbcd_plus_local_guard.py`
  - 新增局部守卫扫描。
  - 测试 `fit_only`、`fit_erode`、`rho_taper`、`relative clamp`、`auto_bad_zero` 等规则。
  - 目的不是定义新物理，而是找出下一版 metric solve 里应施加哪类局部硬约束。

## 关键结果

### 1. 旧 active-edge 解的问题

旧充分迭代输出：

- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_forceedge_coeff_cg6000/`

结果：

- D 方程 exact residual weighted mean：`0.13917`
- 一步 \(g_+\) 重构的测度偏差：`14404.69`
- 负质量壳判别式比例：`0.01922`
- 最小 `det(g_+)`：`-3.74e15`

局部诊断显示：

- 坏点主要来自支撑区边缘/低密度带。
- 全局阻尼无效：即使 `alpha=1e-5`，仍有约 `1.24%` 支撑点出现负判别式。
- 原因不是平均更新太大，而是少数点的相对更新极端大；旧解的 \(|\delta g_+|/|g_+|\) 最大值达到 `7.58e6`。

解释：

- 旧 active-edge 解虽然在拟合区降低了 D 残差，但在拟合区外沿/支撑边缘产生了不可容许的未来度规。
- 这说明不能靠全局 \(\alpha\) 阻尼修复，必须在 metric solve 中加入局部可容许性约束或局部信赖域。

### 2. `no-active-edge` 重新求解消除了大多数异常

新求解输出：

- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_noactiveedge_coeff_cg6000/`

设置：

- `fit_region=core10`
- `metric_active_dilation=0`
- `metric_variable_slices=plus`
- `solver=cg`
- `cg_maxiter=6000`

结果：

- 更新前 residual weighted mean：`1.24344`
- 线性更新后 residual weighted mean：`0.12830`
- 完整 nonlinear 回代 residual weighted mean：`0.13192`
- \(|\delta g_+|/|g_+|\) p95：`4.68e-5`
- \(|\delta g_+|/|g_+|\) max：`0.00931`

解释：

- 不让 metric 变量外扩到 core10 边界外后，求解器仍能把 D 方程残差压低到 `~0.13`。
- 同时 \(g_+\) 更新幅度变得非常温和，没有旧 active-edge 解的百万级局部异常。

### 3. 剩余问题是孤立近退化点，而不是整体失败

对 `no-active-edge` 解做一步可容许性检查：

- 一步 \(g_+\) 测度偏差仍被放大到 `231.998`。
- 负判别式比例为 `1/6034 = 0.0001657`。

定位坏点：

- 坐标：\((x,z)\approx(-3.392\,\mu m,-6.013\,\mu m)\)。
- `rho/max(rho)`：`0.148`，所以不是低密度边缘。
- `X/m^2=-0.40894`，`Q/m^2=-1.40894`。
- 参考 \(\tilde g\) 本身接近退化：
  - `det(g_plus_base) ~ 1.24e-24`
  - `det(g_center_base) ~ -5.04e-24`
  - 一个特征值约 `1e-17`
- `delta g_+ / g_+` 在该点只有 `3.1e-8`。

解释：

- 这个坏点不是更新幅度过大导致，而是 A-derived \(\tilde g\) 参考初态在该格点已经落在近退化分支面上。
- 这和我们之前理论讨论一致：在 \(\sqrt{|\tilde g|}\rho^\tilde\) 或质量壳求根退化处，速度/分支选择本来就不适定。
- 数值上应把这种点当作局部边界/patch 守卫点处理，而不是让质量壳求根强行穿过去。

### 4. `auto_bad_zero` 给出当前最健康初始包

局部守卫扫描输出：

- `visualizations/equation_first_gbcd_plus_local_guard_scan_n384_taum3p5_core10_noactiveedge_cg6000/`

最佳可容许候选：

- 名称：`auto_bad_zero`
- 冻结点：支撑区 `2` 个点，均在 core10。
- D 方程 residual weighted mean：`0.13220`
- residual p95：`0.57969`
- 一步 \(g_+\) 测度偏差：`1.0914e-4`
- 负质量壳判别式比例：`0`
- 非正 `det(g_+)` 比例：`0`
- `disc_plus_min_support`：`2.75e-7`
- `det_plus_min_support`：`9.54e-25`

解释：

- 只冻结 2 个不可容许点，D 方程残差几乎不变：`0.13192 -> 0.13220`。
- 但一步物质重构从灾难性 `231.998` 降到 `1.1e-4`。
- 因此当前问题不是 D 支整体不可解，而是 metric update 需要在退化/分支点上加入局部可容许性守卫。

### 5. guarded 初值包与一步测试

导出的 guarded 初值包：

- `visualizations/equation_first_gbcd_plus_initial_package_n384_taum3p5_core10_noactiveedge_guarded/`

关键诊断：

- exact residual weighted mean：`0.13220`
- \(\rho^\tilde\) 初始拉回偏差：`1.4e-17`
- mass-shell defect p95：`2.24e-13`
- 负判别式比例：`0`
- guard 冻结点：`2`

一步测试输出：

- `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_taum3p5_core10_noactiveedge_guarded/`

关键诊断：

- 初始测度与 A 下一切片偏差：`3.50e-7`
- 一步后在 \(g_0\) 上重构的测度偏差：`5.98e-6`
- 一步后在 \(g_+\) 上重构的测度偏差：`1.091e-4`
- 加权一步 \(g_+\) 测度偏差：`1.043e-4`
- \(g_0\) 上负判别式比例：`0`
- \(g_+\) 上负判别式比例：`0`
- 守恒密度一步相对变化：`5.79e-7`
- \(u_x\) 一步相对变化：`1.11e-8`
- \(u_z\) 一步相对变化：`3.71e-9`

解释：

- 这是目前第一版可用的 D 初始加速度包。
- 它没有证明长时间全演化已经完成，但证明了：至少在左侧分离态这个最难切片上，可以同时满足较低 D 残差、初始物质接近 A、质量壳可解和一步物质推进稳定。

## 物理解释与数值解释的分界

物理层面：

- D 方程不是被改写了；仍然使用
  \[
  \tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2=\mathcal C_{\mu\nu}.
  \]
- `plus-only` 的含义是：\(\tilde g_0,\rho^\tilde,u_i\) 是初始切片，D 方程决定 \(g_+\)，即初始加速度。
- 近退化点本来就是变换表象下分支/速度定义不适定的候选点，应作为 patch 边界或局部守卫点处理。

数值层面：

- `metric_active_dilation=0` 是求解器设置，表示不让核心拟合区外的边缘自由度参与本轮 \(g_+\) 解。
- `auto_bad_zero` 是局部可容许性守卫，表示对质量壳判别式或行列式分支不合法的点冻结 \(\delta g_+\)。
- 这不是刻意降低 D 方程发散；它的代价已经由 residual 统计显示出来：residual 从 `0.13192` 轻微升到 `0.13220`。

## 当前结论

- `center_plus + source` 仍不应作为主线。
- `plus-only` 现在可以视作 D 初始加速度求解器的第一版核心。
- 对 `tau=-3.5`，旧 active-edge 的失败主要是数值求解器边界/退化点处理不当，不是 D 支方程立刻失败。
- 下一版全演化器必须把“可容许性守卫”变成强约束，而不是事后手工修补。

## 下一步

1. 把 `no-active-edge + auto_bad_zero` 封装为标准 D 初始数据生成流程。
2. 对 `tau=0,+3.5` 同样生成 guarded package，确认干涉中心和右侧分离态也满足一步可容许性。
3. 将后续每步 metric solve 改成“先解 D 方程，再投影/守卫不可容许点，再记录 residual 代价”的迭代流程。
4. 对近退化点建立更物理的 body-fitted / patch-boundary 处理，避免长期演化时反复依赖逐点冻结。
