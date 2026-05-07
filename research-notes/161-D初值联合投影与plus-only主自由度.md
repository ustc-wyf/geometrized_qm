# 161-D初值联合投影与plus-only主自由度

日期：2026-05-07

## 目的

用户要求实现真正的 D 初值求解器：在 D 支方程作为硬目标的前提下，使 \(\rho^\tilde\) 拉回到 \(g\) 表象后的 \(\rho\) 与 A 支初态的加权偏差尽量小，并检查后续 full tensor matching 是否可保持。

本轮先做第一版联合投影原型，目的是区分哪些自由度是真正有效的：

- `center_plus + source`：同时允许改中心切片 \(g_0\)、未来切片 \(g_+\)，并允许局部改变
  \[
  \eta=\delta\log(\sqrt{|\tilde g|}\tilde\rho).
  \]
- `plus-only`：保持中心初态 \(\tilde g_0,\rho^\tilde,u_i\) 不变，只通过 \(g_+\) 改变时间二阶导数，即由 D 场方程确定初始加速度。

## 新增脚本

- `kg_examples/solve_gbcd_joint_initial_projection_sparse.py`
  - 解线性化联合系统
    \[
    \delta \tilde G_{\mu\nu}
    -
    \delta(\tilde T_{\mu\nu}/M_P^2)
    =
    \mathcal C_{\mu\nu}
    -
    (\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2).
    \]
  - 可选源项变量 `source_log_measure` 表示 \(\eta\)。
  - 求解后重新用修正后的 \(\tilde g\) 计算 \(\sqrt{|\tilde g|}\)、\(\tilde\rho\)，并重新解质量壳得到 \(u_t\)。
  - 输出完整 nonlinear back-substitution 残差、\(\eta\) 裁剪比例、\(\rho\) pullback 偏差、质量壳判别式和行列式诊断。

- `kg_examples/scan_gbcd_joint_projection_alpha.py`
  - 对联合解做 trust-region 扫描：
    \[
    \tilde g\to \tilde g+\alpha\delta \tilde g,\qquad
    \eta\to\alpha\eta.
    \]
  - 每个 \(\alpha\) 都完整回代 D 张量残差，并同时检查 \(\rho\) 偏差、质量壳判别式和 \(\det\tilde g\)。

## 关键结果

### 1. `center_plus + source` 不健康

高分辨率左侧分离态：

- 输出目录：
  - `visualizations/equation_first_gbcd_joint_initial_projection_sparse_n384_taum3p5_core10_centerplus_rhoclose1e6_v2/`
- 参数：
  - `n=384`
  - `tau=-3.5`
  - `core10`
  - `metric_variable_slices=center_plus`
  - `source_variable_mode=log_measure`
  - `source_closeness_weight=1e6`
  - `eta_clip=0.5`

结果：

- 更新前 exact residual weighted mean：`1.24344`
- 线性层面 residual weighted mean：`0.62368`
- 完整 nonlinear 回代 residual weighted mean：`0.81379`
- \(\eta\) 未裁剪值的 weighted mean 约 `1.05e12`
- \(\eta\) 裁剪比例：`1.0`
- \(\rho\) pullback 加权相对 L1 偏差：`0.51344`
- 负质量壳判别式比例：`0.01882`
- fit 区最小 `det(gtilde)`：约 `-3.93e4`

解释：

- 源项 \(\tilde T_{\mu\nu}/M_P^2\) 被 \(M_P^2\) 强烈压制，用 \(\rho^\tilde\) 幅值去修几何残差极不经济。
- 求解器把 \(\eta\) 推到荒唐量级，说明这不是一个物理健康的自由度。
- 即使裁剪后，仍需要约 `51%` 的 \(\rho\) pullback 偏差，且 exact residual 仍高。
- 因此 `center_plus + source` 不能作为当前 D 初值求解器路线。

### 2. trust-region 扫描否定“小步 center_plus”补救

输出目录：

- `visualizations/equation_first_gbcd_joint_projection_alpha_scan_n384_taum3p5_centerplus_rhoclose1e6_v2/`

结果：

- 若要求 \(\rho\) pullback 偏差 \(\le 5\%\)，可用 \(\alpha\) 约不超过 `0.1`。
- `alpha=0.1` 时 exact residual weighted mean 约 `1.16685`，只比原始 `1.24344` 略好。
- `alpha=1.0` 的 exact residual 最低约 `0.81379`，但对应 \(\rho\) 偏差 `0.51344`，并伴随判别式/行列式问题。
- 在默认“\(\rho\) 偏差 ≤5%、负判别式比例=0、det>0”的约束下，没有可接受的 feasible alpha。

解释：

- `center_plus + source` 不是简单过冲问题。
- 把步长缩小到物理可接受范围后，D 残差改善太弱。
- 因此不应继续沿“改中心度规和密度源项”这条路线硬推进。

### 3. `plus-only` 是当前健康主自由度

输出目录：

- `visualizations/equation_first_gbcd_joint_initial_projection_sparse_n384_taum3p5_core10_plus_nosource_v2/`

参数：

- `metric_variable_slices=plus`
- `source_variable_mode=none`
- `lsqr-maxiter=1000`

结果：

- 更新前 exact residual weighted mean：`1.24344`
- 线性 residual weighted mean：`0.18885`
- 完整 nonlinear 回代 residual weighted mean：`0.21334`
- \(\rho\) pullback 偏差：`0.0`
- \(\eta=0\)，没有源项裁剪问题。
- 中心切片 \(\delta g_0=0\)。
- 负质量壳判别式比例：`0.0`

对照既有充分迭代版本：

- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_forceedge_coeff_cg6000/`
- `plus-only + CG6000` exact residual weighted mean：`0.13917`

解释：

- `plus-only` 不改变初始 \(\tilde g_0,\rho^\tilde,u_i\)，而是让 D 场方程决定下一切片 \(g_+\)，也就是决定初始加速度/二阶时间导数。
- 这正符合 Cauchy 问题的结构：初始数据不一定通过改中心切片来满足方程，而应由约束方程和演化方程共同给出允许的初始速度/加速度。
- 因此当前真正可行的 D 初值求解方向是：
  1. 保持 A-derived 初始物质数据作为近似初态；
  2. 用 gBCD/full tensor 方程解出 \(g_+\) 或等价的 \(\partial_t^2\tilde g\)；
  3. 后续演化时同步推进物质守恒、质量壳和 \(\mathcal C_{\mu\nu}\) 守恒。

## 当前结论

- A-derived 初态直接作为 D 初态并非严格满足 D 方程；但它可以作为近似初始切片。
- 不应通过大幅修改 \(\rho^\tilde\) 来硬配 D 方程，因为 \(\tilde T/M_P^2\) 太小，数值上会要求荒唐的 \(\delta\rho\)。
- 允许改中心切片 \(g_0\) 会把质量壳、\(\rho^\tilde\) 拉回、\(\tilde T_{\mu\nu}\)、辅助应力和连续性全部耦合起来；第一版测试显示这条路线非线性不自洽。
- 当前最稳路线是 plus-only / 初始加速度路线：把 D 场方程看作确定下一切片或二阶时间导数的方程，再接入物质演化。

## 下一步

1. 把 `plus-only` 解法提升为正式 D 初始加速度求解器，而不是叫“初值投影”。
2. 用充分迭代/预条件版本复现 `0.13917` 或更低残差，并保存可被演化器直接读取的 \(g_-,g_0,g_+\) 三切片。
3. 在该三切片基础上推进一小步 D 物质方程，检查：
   - \(\rho^\tilde\) 拉回偏差；
   - D 方程 exact residual；
   - \(\tilde\nabla^\mu\mathcal C_{\mu\nu}\)；
   - 质量壳判别式；
   - full tensor interface matching。
4. 若一步演化后约束漂移明显，再加入 hard/nullspace projection 或 augmented-Lagrangian，而不是调 \(\rho^\tilde\) 幅值。
