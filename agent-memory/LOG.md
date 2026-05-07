# 工作日志

## 长记忆状态

- 当前状态：ON
- 生效方式：本项目后续对话默认启用长记忆，直到用户显式关闭
- 最近显式确认：2026-04-27 跨机器接管时，用户要求“把它加进你的 skill 里，然后开启”

## 轮回信息

- 当前轮回：2
- 轮回说明：跨电脑复制后的恢复轮回

## 日志记录

### 2026-05-07 gBCD 投影型 Einstein-like 方程候选

- 用户再次明确纠偏：最终需要的是“显式可写、可一般化推广的方程或约束条件”，不是停留在数值计算器或逐点拟合。
- 已把当前 equation-first 主线整理为投影型场方程，并写入 `research-notes/164-gBCD投影型Einstein-like方程候选.md`。
- 核心候选：
  \[
  \Pi_E^\perp\left(\tilde G_{\mu\nu}-\frac{1}{M_P^2}\tilde T_{\mu\nu}\right)=0,
  \]
  其中 \(\Pi_E\) 投影到
  \[
  E^I_{\mu\nu}=\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}.
  \]
- 守恒/测地线约束：
  \[
  \tilde\nabla^\mu\Pi_E\left(\tilde G_{\mu\nu}-\frac{1}{M_P^2}\tilde T_{\mu\nu}\right)=0.
  \]
  这等价于要求有效修正张量 \(\mathcal C_{\mu\nu}\) 单独守恒，从而结合物质连续性方程给出 \(\tilde g\)-测地线运动。
- 物质方程沿用此前已证明等价的变换后 Madelung/KG 方程：
  \[
  \tilde g^{\mu\nu}\partial_\mu S\partial_\nu S=m^2,\qquad
  \tilde\nabla_\mu(\tilde\rho\,\tilde g^{\mu\nu}\partial_\nu S)=0.
  \]
- 还必须加分支条件：
  \[
  Q\to0\Rightarrow \Pi_E\left(\tilde G_{\mu\nu}-\frac{1}{M_P^2}\tilde T_{\mu\nu}\right)\to0,
  \]
  否则投影方程只限制残差方向，不保证回到普通 Einstein 方程。
- 当前未完成问题：
  - 投影内积 \(W^{\mu\nu\rho\sigma}\) 的物理选择；
  - Gram 矩阵近退化时的 patch/branch 处理；
  - \(Q\)-门控辅助泛函是否必要；
  - 该 equation-first 候选是否满足 Helmholtz/self-adjoint integrability、能否来自作用量。

### 2026-05-07 gBCD 投影方程闭合条件

- 已继续完成投影方程本身的理论闭合分析，并写入 `research-notes/165-gBCD投影方程闭合条件.md`。
- 关键结论：
  - 投影方程的本质是
    \[
    \mathcal R_{\mu\nu}\in\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}.
    \]
  - 只要投影内积 \(W\) 非退化，且 Gram 矩阵 \(H_{IJ}\) 非退化，\(\Pi_E^\perp\mathcal R=0\) 的在壳解集不依赖 \(W\)；\(W\) 只影响离壳残差和代表元规范。
  - 对内积 \(\langle X,Y\rangle=X_{\mu\nu}Y^{\mu\nu}\)，若
    \[
    a=u^2,\quad b=r^2,\quad c=u\cdot r,\quad \Delta=ab-c^2,
    \]
    则
    \[
    \det H=\frac{d-2}{2}\Delta^3.
    \]
  - 因此该四方向投影在 \(d>2\)、\(\Delta\neq0\) 的 patch 内非退化；严格 \(1+1d\) 中 \(A\tilde g_{\mu\nu}\) 不是独立方向，必须单独降维处理。
  - HJ 壳方程 \(u^2=m^2\) 与 \(u=dS\) 本身已经推出 \(\tilde g\)-测地线；守恒条件的角色是保证引力侧修正与 Bianchi identity、物质连续性方程相容。
  - \(Q\to0\) 分支可更明确写成 \(\lambda_I=\chi(\mathcal Q)\hat\lambda_I\)、\(\chi(0)=0\)、\(\hat\lambda_I\) regular。
  - 是否存在局域作用量仍需 Helmholtz/self-adjoint 检查；当前只确认主部没有立刻显示不可变分。
- 下一步理论任务更新为：先做 3+1d 方程计数/主部分析、\(\Delta=0\) patch 条件、\(Q\to0\) 分支实现和 Helmholtz 条件首轮检查，再回到数值例子验证。

### 2026-05-07 gBCD 投影方程 3+1d 计数与主部

- 已完成 3+1d 方程计数与主部首轮分析，并写入 `research-notes/166-gBCD投影方程的3加1计数与主部.md`。
- 关键结论：
  - 在 3+1d 中，对称二阶张量有 10 个分量，非退化 \(E=\{\tilde g,uu,rr,ur\}\) 子空间有 4 个方向；
  - 因此 \(\Pi_E^\perp\mathcal R=0\) 给出 \(10-4=6\) 个独立条件，和 metric theory 去掉 4 个坐标规范自由度后的度规自由度数量一致；
  - 这只是计数合理，不是 Cauchy 适定证明。
- 主部结论：
  \[
  \delta\mathcal E_{\mu\nu}^{\rm prin}
  =
  \Pi_E^\perp\delta\tilde G_{\mu\nu}^{\rm prin}[h].
  \]
  在 harmonic gauge 下近似为
  \[
  \delta\mathcal E_{\mu\nu}^{\rm prin}
  \simeq
  -\frac12\Pi_E^\perp\tilde\square\bar h_{\mu\nu}.
  \]
- 物理含义：
  - 最高阶传播结构仍来自 \(\tilde g\) 的 Einstein 主部；
  - \(E\) 方向残差被解释为量子几何化诱导的有效源，而不是被设为零；
  - \(Q\to0\) 分支条件 \(\lambda_I=\chi(\mathcal Q)\hat\lambda_I,\chi(0)=0\) 会把 \(E\) 方向也压回零，从而恢复完整 Einstein 方程。
- 下一步剩余理论任务：写出 gauge-fixed principal symbol 在 \(\Pi_E^\perp\) 子空间上的双曲性条件，并继续做 Helmholtz/self-adjoint 检查。

### 2026-05-07 gBCD 主符号缺口与标量闭合条件

- 已继续检查 harmonic gauge 后的 reduced principal symbol，并写入 `research-notes/167-gBCD投影方程的主符号缺口与标量闭合条件.md`。
- 关键结论：
  - 投影方程加 harmonic gauge 的主符号可写成
    \[
    \xi^2\Pi_E^\perp\bar h_{\mu\nu}=0,\qquad
    \xi^\mu\bar h_{\mu\nu}=0.
    \]
  - 若 \(\xi^2\neq0\)，第一式要求 \(\bar h_{\mu\nu}=\lambda_I E^I_{\mu\nu}\)，第二式变为
    \[
    M_{\nu I}(\xi)\lambda_I=0,\qquad
    M_{\nu I}=\xi^\mu E^I_{\mu\nu}.
    \]
  - 因为 \(M\) 的列都在 \(\mathrm{span}\{\xi,u,r\}\) 内，泛型 rank 为 3，所以还剩一个 \(E\)-方向主部模式。
  - 该未定模式显式为
    \[
    N^I=(0,-q^2,-p^2,2pq),\qquad p=\xi\cdot u,\quad q=\xi\cdot r,
    \]
    或
    \[
    \bar h^{(0)}_{\mu\nu}=-w_\mu w_\nu,\qquad
    w_\mu=q u_\mu-p r_\mu.
    \]
- 理论判断：
  - 投影方程本身不是完整 Cauchy 动力学方程；
  - 必须额外加入一个协变标量状态方程，并且它的主部必须满足 \(s_I N^I(\xi)\neq0\)；
  - 这解释了此前辅助场守恒主符号 rank=3、nullity=1 的根源。
- 首个自然候选闭合类：
  - 在 massive branch 中构造 \(u-r\) 二平面正定收缩张量 \(\mathsf q_E^{\mu\nu}\)；
  - 加入
    \[
    \mathsf q_E^{\mu\nu}(\Pi_E\mathcal R)_{\mu\nu}
    =
    \chi(\mathcal Q)\Theta,\qquad \chi(0)=0.
    \]
  - 因为对未定模式 \(N_{\mu\nu}=-w_\mu w_\nu\)，有
    \[
    \mathsf q_E^{\mu\nu}N_{\mu\nu}<0,
    \]
    所以它能看见并固定该主部缺口。
- 下一步：在高斯干涉三切片上检验 \(\mathsf q_E\)-trace 闭合是否比普通 trace 更合理，并做 Helmholtz/action 检查。

### 2026-05-07 \(\mathsf q_E\)-trace 标量闭合三切片检验

- 已修改 `kg_examples/fit_gbcd_trace_closure.py`，新增闭合：
  - `qe0`: \(\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=0\)
  - `qe_q1`: \(\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=\alpha Q/Q_0\)
  - `qe_q2`: \(\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=\alpha Q/Q_0+\beta(Q/Q_0)^2\)
- 已在 `full_resolution=96, fit_region=trusted, force_region=trusted` 下跑三组中心切片：
  - `tau=-3.5`: `visualizations/equation_first_gbcd_qe_trace_closure_n96_taum3p5_trusted/`
  - `tau=0`: `visualizations/equation_first_gbcd_qe_trace_closure_n96_tau0_trusted/`
  - `tau=+3.5`: `visualizations/equation_first_gbcd_qe_trace_closure_n96_taup3p5_trusted/`
- 关键结果已写入 `research-notes/168-qEtrace闭合三切片检验.md`：
  - 普通 trace=0 的中心 residual 加权均值分别为 `0.1185, 0.0264, 0.1103`；
  - 最好的 \(\mathsf q_E\)-trace 版本分别约为 `0.1810, 0.2907, 0.1941`；
  - 因此最小硬闭合 \(\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=0\) 及其简单 \(Q,Q^2\) RHS 版本不被当前数据支持；
  - 普通 trace=0 在这一轮中反而是最小代数闭合。
- 理论修正：
  - 普通 trace 对主符号未定模式的收缩为 \(-w^2\)，所以它只在 \(w^2=0\) 时漏掉该模式；
  - \(\mathsf q_E\)-trace 对未定模式是负定收缩，更安全但不等于更符合代数数据。
- 已补充检查普通 trace 的退化面：
  - 输出 `visualizations/equation_first_gbcd_trace_w2_degeneracy_n96_3tau/summary.json`；
  - \(\tau=-3.5\) 约 `12%~14%` trusted 点在采样主方向上接近 \(w^2=0\)；
  - \(\tau=0\) 约 `2.7%`；
  - \(\tau=+3.5\) 约 `3.7%~5.5%`；
  - 这些点高度对应 \(r^\perp\) 退化/近退化区域，也就是 \(\mathsf q_E\) 自身也需要换 patch 的地方。
- 下一步应改为：
  - 检查普通 trace=0 的潜在退化面 \(w^2=0\) 是否在高斯干涉可信区域出现；
  - 或把 \(\mathsf q_E\)-trace 放入辅助场正定范数，而不是作为硬代数状态方程。

### 2026-05-07 放宽 pure-gtilde 限制后的 u/r 引力作用量路线

- 用户提出新路线：仍采用强约束，但暂时放开“引力作用量必须完全由 `gtilde` 衍生”的限制，允许在 `gtilde` 表象下显含 \(u^\mu\) 与 \(r^\mu\)。
- 已完成理论层面第一轮判断，并写入 `research-notes/144-放宽纯gtilde限制后的ur引力作用量理论分析.md`。
- 用户随后提醒并已修正：当前 `gtilde` 物质作用量不是原 KG/Madelung 物质作用量的逐项拉回；我们证明的是物质方程等价。因此“拉回 \(S_{\rm EH}[g]\) + 新 `gtilde` 物质项”不是原完整总作用量的变量重写，不能自动推出等价 Einstein 方程。
- 核心结论：
  - 若 disformal 变换局域可逆，则“完整总作用量整体拉回”可作为严格等价但同义反复的基准；
  - 仅拉回 EH 引力项：
    \[
    S_{\rm grav}^{\rm pullback}[\tilde g,\rho,S]
    =
    \frac{M_P^2}{2}\int\sqrt{-G[\tilde g,u,r,\rho,\ldots]}\,
    R[G[\tilde g,u,r,\rho,\ldots]].
    \]
    再加上另行构造的 `gtilde` 物质项，是 hybrid action，必须重新变分检查。
  - 显含 \(u,r\) 的引力项会泛型改变 \(S,\rho\) 的变分方程；只有在完整总作用量整体拉回时，这些额外项才自动保证与原方程重组等价。
  - 最自然的低阶可检验 ansatz 是各向异性 Ricci 耦合
    \[
    \int\sqrt{|\tilde g|}\,[A^{\mu\nu}(\mathcal I)\tilde R_{\mu\nu}-2V(\mathcal I)]
    \]
    其中 \(A^{\mu\nu}=a_0\tilde g^{\mu\nu}+a_1\hat u^\mu\hat u^\nu+a_2\hat r^\mu\hat r^\nu+a_3\hat u^{(\mu}\hat r^{\nu)}\)。
- 当前下一步若继续执行：在 1550nm 高斯干涉三切片上实现包含完整 \(\nabla\nabla A^{\mu\nu}\) 导数项的 universal residual fit，检验 \(a_i(\mathcal I),V(\mathcal I)\) 是否能比 pure \(f(R)\)、\(f(R,I_2)\) 显著改善。

### 2026-05-07 从 `g` 表象 Einstein 方程直接变换的分解路线

- 用户提出更结构化的新思路：不先猜 action，而是从 A 支 \(M_P^2G_{\mu\nu}[g]=T^A_{\mu\nu}\) 直接做当前变换，把变换后的方程拆成 pure-\(\tilde g\) 几何部分与必须显含 \(u,r,\rho\) 的剩余部分。
- 已写入 `research-notes/145-从g表象Einstein方程变换得到gtilde方程的分解路线.md`。
- 核心公式：
  - 设 \(g_{\mu\nu}=G_{\mu\nu}[\tilde g,u,r,\rho,\ldots]\)，定义
    \(C^\alpha{}_{\mu\nu}=\Gamma^\alpha{}_{\mu\nu}[G]-\tilde\Gamma^\alpha{}_{\mu\nu}\)；
  - 则
    \[
    G_{\mu\nu}[G]=\tilde G_{\mu\nu}[\tilde g]+\mathcal H_{\mu\nu}[\tilde g,u,r,\rho,\ldots].
    \]
  - 变换后的方程可写成
    \[
    M_P^2\tilde G_{\mu\nu}
    =
    T^{(\tilde m)}_{\mu\nu}
    +
    \mathcal R_{\mu\nu},
    \]
    \[
    \mathcal R_{\mu\nu}
    =
    T^{A}_{\mu\nu}[G,\rho,S]
    -T^{(\tilde m)}_{\mu\nu}[\tilde g,\tilde\rho,S]
    -M_P^2\mathcal H_{\mu\nu}.
    \]
- 判断：
  - \(\mathcal R_{\mu\nu}\) 是当前 \(\tilde g\) 物质 action + EH[\(\tilde g\)] 缺掉的有效源项；
  - 若它能写成 pure-\(\tilde g\) 作用量的 metric variation，则 pure 几何路线仍可行；
  - 若它显式依赖 \(u,r,\rho\) 的张量方向，则直接给出需要加入的非纯几何引力 sector。
- 下一步代码任务：实现 \(C,\mathcal H,\mathcal R\) 诊断脚本，在 1550nm 三切片上做结构分解。
- 用户追问 \(g_{\mu\nu}=G_{\mu\nu}[\tilde g,u,r,\rho,\ldots]\) 的含义，已澄清并更新笔记：
  - \(G[\cdots]\) 不是已知简单公式，而是当前 disformal 变换在指定分支上的形式反映射；
  - 由于 \(\tilde g\) 依赖 \(Q=\Box_g\sqrt\rho/\sqrt\rho\)，反解 \(g\) 一般可能是微分反问题，不是局域代数函数；
  - 第一版 \(\mathcal R_{\mu\nu}\) 诊断不需要反解 \(G\)，直接用 A 参考中已知的 \(g\) 与正向生成的 \(\tilde g\) 计算。

### 2026-05-05 D 支算法状态核查

- 用户询问“算法都搞好了吗”，已重新按长记忆协议核查 `LOG/TASKS/DECISIONS/README` 与当前脚本。
- 当前结论：
  - 高分辨率局部窗口短窗算法已经跑通：大域 A 正频 KG 快照裁剪到 `[-9,9]um`、正确初始化 D 支、用 `rho_D_to_A=m^2 ntilde_D/|X_g[D]|` 比较、`--evolve-mask trusted` 下 `tau=0` 与 `tau=-3.5` 5 分钟诊断稳定。
  - 但正式全窗/生产算法尚未完成：`trusted` bulk + support/active 缓冲层仍是局部窗口边界处理，不是最终 conservative boundary/interface flux。
  - 代码检索未发现已接入 `stencil_mask`、正式 patch-boundary flux 或 conservative interface flux；下一步仍应先把该环节落实，再谈全窗准确生产。
- 对用户的答复口径：不能说“全部搞好”；只能说“短窗可信框架已搞好，完整 D 支全动力学算法还差正式边界/通量闭合与长窗效率验证”。

### 2026-05-05 D 支 fixed-buffer stencil 边界通量地基

- 已完成第一版边界/缓冲层算法升级：
  - `kg_examples/simulate_d_tridomain_full_dynamics.py`
    的 `rk4_patchwise_local_time_matter_step` 支持 `stencil_mask` 与 `write_mask` 分离；
  - `rk4_local_time_matter_step` 支持只在 `write_mask` 上写回 RK 子步，同时用 `stencil_mask` 上的固定边界/halo 数据计算中心差分通量；
  - 旧调用不传 `stencil_mask/write_mask` 时保持原行为。
- `kg_examples/simulate_d_local_window_from_a_snapshot.py` 新增：
  - `--boundary-stencil-mask write|support|active`；
  - 当前推荐短窗试验使用 `--evolve-mask trusted --boundary-stencil-mask active`；
  - summary/fields 记录 `boundary_stencil_count`、`fixed_buffer_count` 和底层 `step_fixed_buffer_cell_count`。
- 物理/数值含义：
  - `trusted` 是写回区 / bulk evolution region；
  - `support/active` 外层只作为固定 buffer/ghost data 给差分 stencil 提供边界通量信息；
  - 这不是阻尼、削峰、滤波或改源项，也不是把 halo 偷偷当 bulk 演化。
- 已通过检查：
  - `python3 -m py_compile kg_examples/simulate_d_tridomain_full_dynamics.py kg_examples/simulate_d_local_window_from_a_snapshot.py`
- `tau=0` fixed-buffer 20 秒 smoke：
  - 输出：`visualizations/d_local_window_snapshot_tau0_fixed_buffer_smoke20/`
  - `steps=20`，无 halving；
  - `rho_pullback_rel_l1_trusted≈2.82e-05`；
  - `disc_min_trusted≈5.31e-08`；
  - `fixed_buffer_count=4112`，最终 step `step_fixed_buffer_cell_count=4332`。
- `tau=0` fixed-buffer 120 秒 smoke：
  - 输出：`visualizations/d_local_window_snapshot_tau0_fixed_buffer_smoke120/`
  - `steps=112`，无 halving；
  - `rho_pullback_rel_l1_trusted≈1.36e-04`；
  - `measure_rel_l1_trusted≈9.36e-04`；
  - `disc_min_trusted≈1.40e-09`。
  - 与旧 trusted-only 120 秒基线相比：误差同量级，质量壳裕度略小。
- `tau=0` fixed-buffer 165 步公平对照：
  - 输出：`visualizations/d_local_window_snapshot_tau0_fixed_buffer_steps165/`
  - `rho_pullback_rel_l1_trusted≈2.04e-04`，旧 trusted-only 165 步基线约 `2.02e-04`；
  - `measure_rel_l1_trusted≈1.41e-03`，旧基线约 `1.39e-03`；
  - `disc_min_trusted≈3.50e-10`，旧基线约 `4.12e-10`。
  - 结论：fixed-buffer stencil 基本不改变短时误差量级，但打开了正式边界通量数据通道。
- `tau=0` fixed-buffer 5 分钟：
  - 输出：`visualizations/d_local_window_snapshot_tau0_fixed_buffer_5min/`
  - `steps=250`，无 halving；
  - `rho_pullback_rel_l1_trusted≈3.29e-04`；
  - `measure_rel_l1_trusted≈2.21e-03`；
  - `disc_min_trusted≈2.78e-11`，仍为正但已接近边界。
- 当前判断：
  - fixed-buffer stencil 是必要的算法地基，但不是最终 conservative finite-volume/interface matching；
  - 它未破坏短窗稳定，但也没有解决质量壳裕度随演化逼近零的问题；
  - 下一步应升级为真正动态 conservative/interface flux：边界 buffer 不应永远冻结，而应由界面 matching / 守恒通量决定其交换量。

### 2026-05-05 有限体积 face flux / Rusanov flux 首轮

- 已在 `kg_examples/simulate_d_tridomain_full_dynamics.py` 中新增：
  - `finite_volume_flux_divergence(...)`；
  - `rk4_local_time_matter_step(..., flux_divergence_mode=central|face|rusanov)`；
  - `face` 模式把守恒量 `n_tau` 的通量散度写成非周期有限体积面通量；
  - `rusanov` 模式在面通量上加入局部 Lax-Friedrichs/Rusanov penalty。
- 已在 `kg_examples/simulate_d_local_window_from_a_snapshot.py` 中新增：
  - `--boundary-flux-mode central|face|rusanov`；

### 2026-05-06 trusted-buffer matching flux 接入

- 用户要求“让 trusted-buffer 界面上的通量由匹配条件决定”，已落实第一版物质守恒通量匹配。
- 代码改动：
  - `kg_examples/simulate_d_tridomain_full_dynamics.py`
    - `finite_volume_flux_divergence(...)` 新增 `mode="matching"`；
    - 在 `trusted-buffer` 面上不再直接用普通 centered/Rusanov 面通量，而是构造一个单一法向守恒流通量；
    - 当前匹配条件是“一侧 trusted/write、一侧 stencil buffer 的法向 conserved-current flux 单值”，数值上用两侧 one-sided flux 的加权 least-squares 解；
    - 新参数 `matching_flux_weight=equal|abs_n`，默认/本轮使用 `abs_n`，即按局部 `|n_tau|` 加权，避免低密度 halo 主导面通量；
    - patchwise local-time 路径新增 `matching_write_mask`，用全局写回区识别物理 `trusted-buffer` 界面，避免把普通 patch/chart 边界误判为物理界面。
  - `kg_examples/simulate_d_local_window_from_a_snapshot.py`
    - `--boundary-flux-mode` choices 扩展为 `central|face|rusanov|matching`；
    - 新增 `--matching-flux-weight equal|abs_n`；
    - `summary.json` 中记录 `matching_flux_weight` 与 `step_face_flux_fv_matching_*` 诊断。
- 重要限定：
  - 这一步是物质连续性方程的 conserved-current interface flux matching；
  - 它不是 full tensor Einstein interface matching，也不是把张量 jump 条件作为强边界条件接入；
  - `rusanov_strength` 在 `matching` 模式下只是保留在参数记录里，不参与匹配通量计算。
- 验证：
  - `python3 -m py_compile kg_examples/simulate_d_tridomain_full_dynamics.py kg_examples/simulate_d_local_window_from_a_snapshot.py kg_examples/coordinate_matter_evolution.py` 通过；
  - `tau=0, matching, steps=20`：
    - 输出 `visualizations/d_local_window_tau0_matching_flux_steps20/`
    - `rho_pullback_rel_l1_trusted≈2.798e-05`，core `≈5.025e-06`
    - `measure_rel_l1_trusted≈1.838e-04`
    - `disc_min_trusted≈5.310e-08`
    - 最后一步匹配界面面数：x `1574`，z `596`
    - `step_face_flux_fv_matching_flux_jump_l1≈2.622e-02`
  - `tau=0, matching, steps=60`：
    - 输出 `visualizations/d_local_window_tau0_matching_flux_steps60/`
    - `rho_pullback_rel_l1_trusted≈7.433e-05`，core `≈1.508e-05`
    - `measure_rel_l1_trusted≈5.066e-04`
    - `disc_min_trusted≈8.504e-09`
    - 最后一步匹配界面面数仍为 x `1574`，z `596`
    - 无 halving、无负质量壳。
- 与无耗散 face / `rusanov_strength=0` 对照：
  - 20 步 face trusted L1 `≈2.823e-05`，matching `≈2.798e-05`；
  - 60 步 face trusted L1 `≈7.507e-05`，matching `≈7.433e-05`；
  - 质量壳裕度在同一步数下基本相同；
  - 初步看 matching 没有引入 Rusanov 式耗散失真，但需要继续做 `equal/abs_n`、`tau=±3.5` 和更长步数对照。

### 2026-05-06 matching vs face 三切片对照

- 按用户要求先做对照，比较 `boundary_flux_mode=matching` 与 `face`。
- 同参数：
  - `full_resolution=640`
  - `window_um=9`
  - `dt_old=2.5e-7`
  - `evolve_mask=trusted`
  - `boundary_stencil_mask=active`
  - `freeze_uncovered_chart`
  - `chart_boundary_halo=1`
  - `local_time_max_tilt=5`
  - `local_time_step=0.5`
  - `trusted_erosion=1`
  - `active_dilation=2`
- 对照结果：
  - `tau=-3.5`：
    - matching 60 步：trusted L1 `1.909190e-04`，core `1.724189e-04`，measure trusted `9.233887e-05`，disc `2.809625e-07`，matching faces `2106`
    - face 60 步：trusted L1 `1.909189e-04`，core `1.724022e-04`，measure trusted `9.237704e-05`，disc `2.809625e-07`
    - 结论：分离前/左侧窗口中二者几乎重合，matching 未破坏稳定性。
  - `tau=0`：
    - matching 60 步：trusted L1 `7.433061e-05`，core `1.507620e-05`，measure trusted `5.066300e-04`，disc `8.504252e-09`，matching faces `2170`
    - face 60 步：trusted L1 `7.506592e-05`，core `1.565409e-05`，measure trusted `5.069982e-04`，disc `8.504252e-09`
    - 结论：干涉中心 matching 略低于 face，但差异仍小；还不能仅凭短窗宣称物理优越。
  - `tau=+3.5`：
    - matching 因 300s wall time 截断在 45 步：trusted L1 `3.608342e-05`，core `3.235098e-05`，measure trusted `4.141179e-05`，disc `5.889474e-08`，matching faces `2364`
    - face 因 300s wall time 截断在 51 步：trusted L1 `4.095739e-05`，core `3.667532e-05`，measure trusted `4.701376e-05`，disc `5.889474e-08`
    - 结论：两者均未物理失败；因步数不同，只能说明同阶稳定，不能作严格优劣比较。
- 当前判断：
  - 第一版 matching 是可运行的物质界面通量闭合，不是“找不到符合规则的数值解”；
  - 它和 face 的短窗差别很小，说明 trusted-buffer 界面当前不是主导误差的唯一来源；
  - 真正的 full matching 缺口在于还没有把张量几何匹配方程与物质通量一起组成耦合界面求解器。
  - summary 记录 `boundary_flux_mode`，并说明 `rusanov` 是数值通量选项，不是新物理源项。
- 已通过：
  - `python3 -m py_compile kg_examples/simulate_d_tridomain_full_dynamics.py kg_examples/simulate_d_local_window_from_a_snapshot.py`
- `tau=0, face flux, 20s`：
  - 输出 `visualizations/d_local_window_snapshot_tau0_face_flux_smoke20/`
  - `steps=20`；
  - `rho_pullback_rel_l1_trusted≈2.82e-05`；
  - `disc_min_trusted≈5.31e-08`；
  - 结果几乎等同 fixed-buffer central，因为 centered face flux 在内部与中心差分代数接近。
- `tau=0, rusanov flux, 20s`：
  - 输出 `visualizations/d_local_window_snapshot_tau0_rusanov_flux_smoke20/`
  - `steps=19`；
  - `rho_pullback_rel_l1_trusted≈2.92e-05`；
  - `disc_min_trusted≈5.31e-08`。
- `tau=0, rusanov flux, 120s`：
  - 输出 `visualizations/d_local_window_snapshot_tau0_rusanov_flux_smoke120/`
  - `steps=98`；
  - `rho_pullback_rel_l1_trusted≈1.31e-04`；
  - `measure_rel_l1_trusted≈8.27e-04`；
  - `disc_min_trusted≈8.06e-09`。
- 对照：
  - fixed-buffer central 120s：`steps=112`，`rho L1 trusted≈1.36e-04`，`disc_min_trusted≈1.40e-09`；
  - Rusanov 120s：步数较少但质量壳裕度更宽，误差同量级。
- 当前判断：
  - face 模式完成了“显式面通量表达”，但不是新的稳定机制；
  - Rusanov 模式可能改善质量壳裕度，但含数值耗散，不能直接作为最终物理算法；
  - 下一步必须做 Rusanov penalty 强度/收敛检查，或改用更物理的界面 matching flux，避免用耗散掩盖真实结构。

### 2026-05-06 Rusanov penalty 强度扫描

- 已新增参数：
  - `--rusanov-strength`；
  - `0` 等于 centered face flux；
  - `1` 是标准 Rusanov；
  - 中间值用于检查数值耗散强度趋势。
- 已通过：
  - `python3 -m py_compile kg_examples/simulate_d_tridomain_full_dynamics.py kg_examples/simulate_d_local_window_from_a_snapshot.py`
- `tau=0` 固定 20 步扫描：
  - `strength=0`：trusted L1 `2.82e-05`，core L1 `5.22e-06`，disc `5.31e-08`，`n_tau_delta≈1.09e-07`
  - `0.25`：trusted L1 `2.84e-05`，core L1 `5.42e-06`，disc `5.31e-08`，`n_tau_delta≈1.55e-07`
  - `0.5`：trusted L1 `2.87e-05`，core L1 `5.77e-06`，disc `5.31e-08`，`n_tau_delta≈2.55e-07`
  - `1.0`：trusted L1 `3.06e-05`，core L1 `7.51e-06`，disc `5.31e-08`，`n_tau_delta≈4.74e-07`
- `tau=0` 固定 60 步扫描：
  - `strength=0`：trusted L1 `7.51e-05`，core L1 `1.57e-05`，disc `8.50e-09`
  - `0.25`：trusted L1 `7.56e-05`，core L1 `1.63e-05`，disc `8.50e-09`
  - `0.5`：trusted L1 `7.65e-05`，core L1 `1.73e-05`，disc `8.50e-09`
  - `1.0`：trusted L1 `8.22e-05`，core L1 `2.25e-05`，disc `8.50e-09`
- 近 100 步并行扫描（受 wall time 影响，各组 `91-92` 步）：
  - `0`：trusted L1 `1.12e-04`，disc `1.23e-08`
  - `0.5`：trusted L1 `1.13e-04`，disc `1.07e-08`
  - `1.0`：trusted L1 `1.22e-04`，disc `1.07e-08`
- 当前判断：
  - 在固定步数下，Rusanov 强度增加并未改善质量壳裕度；
  - 误差与 `n_tau` 单步变化随强度单调增大；
  - 因而 Rusanov 不应作为默认物理算法，只保留为稳定性/耗散诊断；
  - 当前守恒通量主线应回到 `face/strength=0` 的无耗散有限体积基线，并继续发展 interface/tensor matching flux。

### 2026-04-29

- 用户明确要求：不要再停留在同步规范的失稳诊断上，直接换更稳的规范
- 已据此把 `B/C` 两支的全动力学原型从“固定 `N=1, N^i=0` 的同步规范”扩展为可选规范：
  - `gauge_mode = synchronous`
  - `gauge_mode = 1plog`
- 当前 `1+log` 实现保持零 shift，仅把 lapse 作为新演化变量接入：
  `lapse_t = -2 lapse K`
- 已完成的实现层修改：
  - `kg_examples/simulate_ab_adm_synchronous.py`
    中的 `geometry_rhs` 已加入 lapse Hessian 项
  - `kg_examples/simulate_bc_from_a_initial_data.py`
    中的 `B/C` 状态已加入 `lapse`
  - 物质特征边界条件已改为用当前 `lapse` 计算特征速度，而不再默认 `N=1`
- 已通过的最小回归：
  - `localized + C(flat) + 1plog + steps=40`
    下 `C(flat)` 仍逐项严格退回 `B`
  - 该回归中 `lapse` 已轻微偏离 `1`，说明新规范确实接管切片而非仍是同步规范
- 已完成的短中时间测试：
  - `localized + ell=0.2 + elliptic_lambda=1e-4 + 1plog + steps=800`
    即 `t=0.2`
  - 非平凡 `C` 仍与 `B` 分开：
    `max|h_xx^(C)-h_xx^(B)| ~ 5.07e-4`
    `max|rho_C-rho_B| ~ 1.64e-4`
- 当前轻量长时间诊断脚本：
  `kg_examples/diagnose_bc_overflow.py --gauge-mode 1plog --output visualizations/bc_overflow_diagnosis_1plog --checkpoint-every 200`
- 到当前最近检查点 `t=0.55` 为止，`1+log` 仍稳定：
  - `B(flat)`：
    `beta_max ~ 7.89e-7`
    `lapse_min ~ 0.9999196`
  - 非平凡 `C`：
    `beta_min ~ -5.04e-4`
    `lapse_min ~ 0.9969153`
    `phi_minus_1_abs_max ~ 1.00e-3`
- 当前尚未得到 `1+log` 的首次失败时刻；长时间诊断仍在继续

### 2026-04-27

- 围绕当前新约定 `Q = + \Box \sqrt{\rho} / \sqrt{\rho}`，继续重审纯 rank-1 ansatz 下的 `m=0` 分支
- 明确修正了一处此前过强的说法：要把变换后连续性方程恢复到原方程，不必逐点强制
  `sqrt(-g~) ρ~ (C + D X) = sqrt(-g) ρ`；更弱、也更本质的方程级条件是
  `∂_μ[ sqrt(-g~) ρ~ (κ/X) u^μ ] = 0`
- 在 `m=0` 时，原 HJ 方程变成 `X = Q`；若再要求 `Q -> 0` 时恢复 `g~ -> g` 与 `ρ~ -> ρ`，则固定非零 `κ` 会立刻失败：
  因为 `X -> 0` 时 `CX + DX^2 = κ` 不可能与 `C -> 1, D -> 0` 同时成立
- 同时确认：`κ = 0` 也不能救回当前框架；这时 HJ 变成 `X(C + DX)=0`，在一般 `X ≠ 0` 区域迫使 `C + DX = 0`，使变换后连续性方程塌成平凡的 `0 = 0`
- 当前最稳判断：
  在现有物质作用量
  `I~ = ∫ sqrt(-g~) ρ~ (g~^(μν) u_μ u_ν - κ)`
  和纯 rank-1 ansatz
  `g~^(μν) = C g^(μν) + D u^μ u^ν`
  下，`m=0` 分支不存在一个正则的固定常数 `κ` 与正则函数 `C`
  能同时满足：
  1. 变换后 HJ 方程
  2. 非平凡连续性方程恢复
  3. `Q -> 0` 时 `g~ -> g`、`ρ~ -> ρ`
- 物理解释也已收紧：
  这不是 `C=1` 选错，而是当前“密度乘壳约束”的物质作用量在 null shell 上先天退化；
  后续若认真处理 `m=0`，应改用带辅助场/拉格朗日乘子的 null-particle 型作用量，而不是继续沿用当前 massive-like 写法
- 继续追问“若去掉渐进恢复约束会怎样”后，得到更细的分支判断：
  - 若仍要求 `g~` 是非退化 Lorentzian 度规，并且连续性方程保持非平凡，则必须取 `κ ≠ 0`
  - 此时 `C` 不再被方程恢复固定，只需满足 `C ≠ 0`，而 `D` 由 `D = (κ - C X)/X^2` 给出
  - 但为了让 `sqrt(-g~)` 与 `ρ~` 保持实数，需要 `sign(C) = sign(X/κ)`；因此它只能在 `X = Q` 不过零且符号固定的局域 patch 内成立
  - 若坚持 `κ = 0`，则 `C + DX = 0`，使 `det(g~^(..)) = 0`，`g~` 退化成零壳方向上的退化几何，同时连续性方程塌成平凡恒等式
- 当前新增理解：
  `m=0` 时，去掉渐进条件后并不是“问题消失”，而是把问题从“无经典恢复”改成了“只能局域 patch 化”或“接受退化几何”
- 再次校正：`κ ≠ 0` 才有非平凡连续性方程、以及 `X` 过零/变号导致的 patch 问题，其实对 `m≠0` 与 `m=0` 都成立；`m=0` 真正特有的点只是 `Q -> 0` 与 `X -> 0` 重合，使弱量子/经典极限本身就落在分支边界上
- 围绕体积元与 signature change 又做了一轮外部核对，当前最稳的几何规范结论如下：
  - 必须区分“体积 n-形式（volume form）”与“正测度/密度（positive measure/density）”
  - 对定向的伪黎曼流形，规范的体积 n-形式在正定向坐标中写作 `vol_g = sqrt(|det g|) dx^1 ∧ ... ∧ dx^n`
  - 若坐标变换反转定向，则同一个 `vol_g` 的坐标表达式会多一个负号；这不是 `sqrt(|g|)` 失效，而是 n-形式本来就会随取向翻号
  - 若只想得到“始终为实且非负”的积分测度，则应使用密度 `dμ_g = sqrt(|det g|) |d^n x|`，而不是带定向符号的 n-形式
  - 在标准 GR 中常写 `sqrt(-g) d^4x`，只是因为预先固定了 Lorentzian signature，于是 `det g < 0`，这其实是 `sqrt(|det g|)` 的简写
  - 一个真正的坐标变换不会改变 signature；若某个场依赖变换让 signature 改变，那得到的是一个新的度规场，不是单纯换坐标
  - 若度规平滑且非退化，则 signature 在连通区域内保持常数；因此 signature change 必然伴随某处的度规退化（`det g = 0`）或不连续
- 标准 Lorentzian GR 本身并不涵盖 signature-changing 点；标准处理是把两侧各自当作固定 signature 的区域，在中间退化超曲面上另加拼接/边界条件或改用一阶/密度化变量
- 用户提出一个新的混合逆度规候选：
  `D=Z, E=-Y, F=X`
  （当前记号下的 ansatz 为
  `g~^(μν)=C g^(μν)+D u^μu^ν+E(u^μr^ν+r^μu^ν)+F r^μr^ν`）
- 已核对该解的核心代数结构，记
  `Δ := XZ - Y^2`
  后，有
  `M^(μν):=Z u^μu^ν - Y(u^μr^ν+r^μu^ν) + X r^μr^ν`
  满足：
  `M^(μν)u_ν = Δ u^μ`
  `M^(μν)r_ν = Δ r^μ`
  且对所有同时正交于 `u,r` 的向量 `v_ν`，有 `M^(μν)v_ν = 0`
- 这说明该混合块本质上是在 `span{u,r}` 这个二维平面上乘以 `Δ`，在其正交补上为零；因此它可理解为“Bohm 平面”的投影型修正
- 在此结构下：
  - `g~^(μν)u_ν = (C+Δ) u^μ`
  - `g~^(μν)u_μu_ν = (C+Δ) X`
  - 四维下 `det(g~^(..))/det(g^(..)) = C^2 (C+Δ)^2`
- 若要求 `g~` 系 HJ 方程 `g~^(μν)u_μu_ν = m^2` 精确等价于 `g` 系原 HJ 方程 `X = m^2 + Q`，则唯一的 `C` 为
  `C = 1 - Q/X - Δ`
- 此时
  `C + Δ = 1 - Q/X`
  因而：
  - `g~^(μν)u_μu_ν = X - Q`
  - `g~^(μν)u_ν = (1 - Q/X) u^μ`
  所以只要定义密度变换
  `sqrt(|g~|) ρ~ (1 - Q/X) = sqrt(|g|) ρ`
  就能把 `g~` 系连续性方程精确拉回 `g` 系原方程
- 当前判断：
  这个解的核心想法是正确的，而且结构非常漂亮；但“变换行列式恒正”要更精确地表述为：
  `det(g~^(..))/det(g^(..)) = C^2 (1-Q/X)^2`
  因此它是 manifestly 非负的，并且在 `C ≠ 0` 且 `1-Q/X ≠ 0` 的 patch 上严格为正
- 关键保留：
  - 若 `m=0` 且在壳上 `X=Q`，则第二因子为零，行列式退化
  - 即使 `m≠0`，若 `C=0`（即 `1 - Q/X = Δ`），也会出现退化面
- 继续分析用户提出的“对 `D,E,F` 同时乘或除以 `Δ` 或别的函数以恢复经典极限”后，当前最稳的统一写法是：
  令
  `M^(μν) := Z u^μu^ν - Y(u^μr^ν+r^μu^ν) + X r^μr^ν`
  并取
  `g~^(μν) = C g^(μν) + λ M^(μν)`
  其中真正起作用的不是 `λ` 本身，而是
  `s := λ Δ`
- 因为 `M` 在 `span{u,r}` 上按 `Δ` 作用、在正交补上为零，所以：
  - `g~^(μν)u_ν = (C+s)u^μ`
  - `g~^(μν)u_μu_ν = (C+s)X`
  - 若要求 `g~` 系壳方程等价于 `g` 系原 HJ 方程并取 `κ=m^2`，则必须有
    `C+s = m^2/X = 1 - Q/X`
    即
    `C = m^2/X - s`
- 因而混合块是否在经典极限下消失，完全由 `s=λΔ` 控制：
  - 要恢复 `g~ -> g`，必须要求 `s -> 0` 且 `m^2/X -> 1`
  - 原始漂亮解 `λ=1` 对应 `s=Δ`，所以它只有在额外满足 `Δ -> 0` 时才自动恢复原度规；单靠 `Q -> 0` 并不保证这一点
  - 若把 `D,E,F` 同时乘上 `Δ`，则 `s=Δ^2`，这会更容易在有限 `Δ` 情况下压回零
  - 若把 `D,E,F` 同时除以 `Δ`，则 `s=1`，反而会破坏经典恢复；除非更一般地取 `λ=f/Δ`，并要求 `f(Q,...) -> 0`
- 当前最自然的重参数化是把混合块写成投影型：
  在 `Δ ≠ 0` 的 patch 上定义 `P := M/Δ`
  （作为 `u-r` 平面上的投影算子），然后取
  `g~^(μν) = (m^2/X - f) g^(μν) + f P^(μν)`
  其中 `f(Q,...) -> 0` 时自动有经典恢复
- 继续要求 `C=1` 后，投影型解族中的自由函数被唯一固定为
  `f = m^2/X - 1 = -Q/X`
  因此混合投影解可写成
  `g~^(μν) = g^(μν) - (Q/X) P^(μν)`
  或等价地
  `g~^(μν) = g^(μν) - [Q/(XΔ)] M^(μν)`（`Δ ≠ 0` 的 patch 上）
- 这一步的含义是：
  - `C=1` 与恢复原 HJ 方程一起，把 `f` 完全锁死，不再自由
  - 对 `m≠0`，因 `X=m^2+Q`，有 `f -> 0` 当 `Q -> 0`，故自动恢复原度规
  - 对 `m=0`，壳上 `X=Q`，于是 `f=-1`，经典恢复仍不存在
- 针对用户随后要求“先搁置 `X=0`，只讨论 `m=0` 且仍取 `κ≠0`”又做了一次修正核对：
  - 若要从 `\tilde g` 系 HJ 方程 `\tilde g^(μν)u_μu_ν = κ` 推回 `g` 系原 HJ 方程 `X = Q`，投影块缩放不能写成 `κ/X` 型，而必须写成 `κ/Q` 型
  - 在 `C=1`、`Δ≠0` 的 patch 上，正确选择是
    `λ = (κ/Q - 1)/Δ`
    因而
    `D = Z(κ/Q - 1)/Δ`
    `E = -Y(κ/Q - 1)/Δ`
    `F = X(κ/Q - 1)/Δ`
  - 这时
    `\tilde g^(μν)u_ν = (κ/Q) u^μ`
    `\tilde g^(μν)u_μu_ν = (κ/Q) X`
    所以变换后 HJ 方程等价于 `X=Q`
  - 连续性方程也可通过
    `sqrt(|g~|) ρ~ (κ/Q) = sqrt(|g|) ρ`
    精确拉回
  - 行列式比值变成
    `det(g~^(..))/det(g^(..)) = (κ/Q)^2`
    因而在 `Q≠0` 的 patch 上严格为正
  - 该分支的代价仍是：`Q -> 0` 时系数发散，故无经典恢复；同时还需额外避开 `Δ=0`
- 已把上述结果系统整理成新研究笔记：
  [062-混合投影型逆度规解的统一整理.md]
  该笔记并排总结了：
  - `m≠0, C=1` 的混合投影主分支
  - `m=0, C=1, κ≠0` 的局域 patch 分支
  并明确把 `X=0`、`Q=0`、`Δ=0` 留作下一步单独讨论
### 2026-04-28

- 用户明确纠正当前短 PPT 的主线：要统一的不是“`X=0` 不是灾难”，而是“存在一类类似 Bekenstein 的变换，在不改变因果序的前提下，把 `g` 表象下的方程保持不变，并在 `g~` 表象下给出测地线解释”
- 据此把短 deck `presentation-workspaces/x0-boundary-summary` 彻底改写成新的 7 页结构：
  1. 允许的几何化变换的总主张
  2. 三条最小成功标准：因果序保持、原方程等价、`g~` 表象下测地线化
  3. 投影型混合块作为当前最自然的变换骨架
  4. `g~` 系 HJ / 连续性方程如何精确拉回原 `g` 系
  5. `g~` 表象中的测地线解释
  6. `X=0` 只作为边界与图册注记，不再占据总叙事中心
  7. 当前统一结论与后续检查方向
- 新版短 deck 路径：
  - PPTX：`presentation-workspaces/x0-boundary-summary/output/output.pptx`
  - 源码：`presentation-workspaces/x0-boundary-summary/src/deck.mjs`
  - 预览：`presentation-workspaces/x0-boundary-summary/scratch/previews/slide-01.png` 至 `slide-07.png`
- 已重新导出并跑 QA：
  - `quality-report.json` 无 failures、无 warnings
  - 抽查封面、变换结构页、边界注记页，版式正常
- 当前统一叙事已经收紧为：
  - 先定义“什么样的变换才算成功”
  - 这个成功标准由“因果序保持 + 方程等价 + 测地线化”给出
  - `X=0`、`Q=0`、`Δ=0` 的讨论只能放在后半部分当作边界/图册检查，不能反过来主导总叙事

- 围绕“度规退化、变换发散、奇点该如何理解”查阅了标准广义相对论、Bekenstein/disformal、signature change 和 mimetic gravity 文献
- 当前最稳的文献分层如下：
  1. 坐标奇点（coordinate singularity）：
     只是不好的 chart；标准做法是换 chart，例如 Schwarzschild 地平线用 Eddington-Finkelstein / Kruskal 坐标穿过
  2. 场重定义奇点（field-redefinition singularity / non-invertible transformation）：
     若 disformal 变换失去可逆性，就不能再把它当作“等价改写”；继续取这一极限通常对应新理论，典型例子是 mimetic gravity
  3. 真曲率奇点（curvature singularity）：
     若曲率不变量发散或 geodesically incomplete，则不能靠坐标变换消掉；这通常表示经典 GR 失效，而不是坐标没选好
  4. signature change / 退化度规：
     标准 Lorentzian GR 默认连续、非退化、固定 signature 的度规；signature 变化必然伴随退化面或不连续，且该面上没有公认唯一的 Einstein 方程标准写法
- 对用户提出的三个“绕过”方案做了阶段判断：
  - `iε` 式在分母加小虚部：可作计算正则化或复度规路径积分启发，但不是标准经典 GR / 标准 Bekenstein 处理；若字面采用，会把实 Lorentzian 度规改成复度规
  - “真空场会阻止这种情况”：没有一般性定理支持；Schwarzschild 本身就是真空解，却同时包含地平线坐标奇点和 `r=0` 真奇点
  - “可通过某些变换绕过发散”：只对坐标奇点和某些纯变量奇点成立；若遇到真正的曲率奇点或 disformal 非可逆面，则至多改写成 patch 边界或新理论，不能仍宣称与原理论全局等价
- 对本项目当前分支的具体启示：
  - `X=0`、`Q=0`、`Δ=0` 应先被当作“变换图册（atlas）边界 / 非可逆面候选”来审查，而不是先当作 `g` 系本身的物理奇点
  - 若原 `g,ρ,S` 侧的局域不变量都有限，而只有变换系数爆炸，那么更可能是“坏变量/坏 patch”，不是原理论坏掉
  - `Δ=0` 时投影子 `P=M/Δ` 会坏，但未归一化混合块 `M` 不一定坏；后续讨论退化面时应优先用 `M` 或 `s=λΔ` 这类不显含 `1/Δ` 的变量
- 将这一轮判断正式整理成新文档：
  `research-notes/063-度规退化、不可逆变换与真奇点的文献分层.md`
- 新增一个对后续很关键的分层结论：
  - 若发散只出现在变换系数里，而原 `g,ρ,S` 及其不变量都保持有限，则应优先按“场重定义奇点 / atlas 边界”处理
  - 只有当原几何不变量或 geodesic completeness 也失效时，才应把它升级理解为真曲率奇点
- 与用户就这一点继续收敛：
  - 当前形成的暂时共识是：若原场量 `ρ,S,g` 及其构造出的原理论局域不变量都不发散，则 `X=0`、`Q=0`、`Δ=0` 这类只在变换中引起 `1/X`、`1/Q`、`1/Δ` 爆炸的点，应优先视为 atlas 边界
  - 下一步应把这个“atlas 边界判据”写成更严格的项目内局域判据，并分别检查 `X=0`、`Q=0`、`Δ=0` 的特殊性
- 用户进一步提出一个很值得保留的观点：
  - 某些退化面上的发散，可能对应“当前 disformal 变换退化成更低秩的变换”，因此应改用更低层的表示，例如从混合 disformal 退回到纯 rank-1，甚至纯 conformal
- 当前对此的阶段判断：
  - 该观点对 `Δ=0` 最有说服力，尤其在 `X≠0` 时，`Δ=0` 意味着 `r_μ` 与 `u_μ` 失去独立性，混合平面结构塌缩，当前 `u-r` 投影型图册很可能应退回更低层表示
  - 但对 `Q=0` 与 `X=0` 不能直接一概而论；是否退回共形，要分分支和壳约束再看
- 又新增一个重要技术点：
  - `∂_μ ρ` 在 `ρ=0` 处未必发散，甚至常常为 0；真正更容易出问题的是 `r_μ = ∂_μ sqrt(ρ)`
  - 对简单节点，`sqrt(ρ)` 往往只 Lipschitz 连续而不可微，因此当前混合 `u-r` ansatz 的定义域本身就天然应限制在 `ρ>0` 的 patch 上，或改用别的振幅变量
- 用户继续追问：当前混合/纯 ansatz 是否算 Bekenstein 变换，以及是否有因果性问题
- 当前判断：
  - 这些变换最准确的定位是：Bekenstein 型（Bekenstein-type）或更一般的 disformal 变换
  - 纯 rank-1 `g~^(μν)=C g^(μν)+D u^μu^ν` 更贴近标准单方向 disformal；混合 `u-r` 版本则已超出 Bekenstein 1992 原始“单个标量梯度”形式，更像多梯度/广义 disformal
  - 因果性要分三层：
    1. `g~` 自身是否仍是可逆 Lorentzian 度规
    2. `g~` 的锥与原 `g` 的锥如何比较
    3. `g~` 是否被当作所有物质共同耦合的物理度规
  - 对当前混合投影型 ansatz，局域上保持同一 Lorentzian 分支的最简单充分条件是 `C>0` 且 `C+s>0`；这里 `s=λΔ`
  - 因此：
    - `m≠0, C=1` 主分支上有 `C+s=m^2/X`，在 `X>0` patch 内因果结构局域上是良性的
    - `m=0, C=1, κ≠0` patch 上需额外要求 `κ/Q>0`；否则 `u-r` 平面的符号会翻转，可能改变 signature / 因果解释
  - 若只把 `g~` 当作该单个 KG 标量 Bohm/HJ 流的重写几何，而不把它立刻宣布为所有场共享的物理时空度规，那么局域锥的差异本身不等于观测性因果悖论；真正的双锥冲突只有在把 `g` 和 `g~` 同时当作不同物质的物理时空时才会升级成 serious causality issue
- 用户要求回到纯共形变换，单独检查它能否满足此前对非共形变换提出的要求
- 已完成统一复核并整理成：
  `research-notes/064-纯共形变换能否满足几何化要求.md`
- 当前结论分层如下：
  - 在“先变分后代回”的方程级口径下，纯共形变换并不自动失败；若允许 `ρ~` 真正一起变换，则它可以恢复原 HJ 方程、恢复原连续性方程，并给出 `g~` 下的测地线解释
  - 对 `m≠0` 主分支，取 `κ=m^2` 时有 `C=m^2/X=m^2/(m^2+Q)`，因此 `Q->0` 时还能恢复 `g~->g`
  - 但若把 `C=κ/X` 和对应的 `ρ~` 变换直接拉回旧变量作用量，则 integrand 恒等塌成零；因此纯共形不能同时充当一个非退化的 pulled-back 有效作用量
  - 此外，纯共形不会改变 null cone，所以它不能承担那些依赖“新因果锥”或“类时/类空重解释”的任务
- 用户提出另一组混合逆度规解：
  `D=Z, E=-aY, F=aX, C=aXZ`
- 已完成检查并整理成：
  `research-notes/065-混合ansatz新解D等于Z-E等于负aY-F等于aX.md`
- 当前结论：
  - 这组解同样会把 `g~^(μν)u_ν` 压回 `u^μ` 方向，具体为
    `g~^(μν)u_ν = (XZ + aΔ) u^μ`
  - 因此只要在 `XΔ≠0` 的 patch 上取
    `a = [κ - X^2 Z] / (XΔ)`
    并定义
    `sqrt(|g~|)ρ~(XZ+aΔ)=sqrt(|g|)ρ`
    就能同时恢复 `g` 系 HJ 方程和连续性方程
  - 四维下的一般行列式比值为
    `det(g~^(..))/det(g^(..)) = a^3 X^2 Z^2 (2XZ - Y^2)(XZ + aΔ)`
  - 与投影型主解相比，这组解的退化结构更复杂：
    - bare ansatz 自身就在 `X=0`、`Z=0`、`2XZ-Y^2=0`、`a=0` 等面退化
    - 为恢复方程而解出 `a` 后，`Δ=0` 也变成 matching chart 的显式发散面
    - 额外还有 `κ = X^2 Z` 导致 `a=0, C=0` 的退化面
- 当前总体评价：这是一个代数上成立的解，但从经典极限自然性和退化结构复杂度看，弱于现有 `C=1` 的投影型主解
- 用户继续追问：一般 `C,D,E,F` 下，`g~` 的行列式何时恒正
- 已完成总公式整理并单列为：
  `research-notes/066-一般混合逆度规ansatz的行列式恒正条件.md`
- 当前结论：
  - 四维里的一般公式是
    `det(g~^(..))/det(g^(..)) = C^2 [ C^2 + C(DX+2EY+FZ) + (DF-E^2)Δ ]`
  - 因此 strict positivity 的充要条件就是：
    `C ≠ 0`
    且
    `C^2 + C(DX+2EY+FZ) + (DF-E^2)Δ > 0`
  - pure conformal、rank-1、投影型混合和新解 `D=Z,E=-aY,F=aX,C=aXZ` 都可作为该总公式的特例来读
  - 同时再次明确：
    `det` 比值严格为正只是“不退化且同号”的第一层条件，不自动等于保持同一 Lorentzian signature 或无双锥因果问题
- 用户继续提出一个很关键的检验点：若把 `κ` 改成不等于 `m^2` 的固定常数，`X=0` 的 blow-up 会不会“移动”
- 当前核对结论：
  - 对 pure conformal、rank-1 以及 `C=1` 的投影型混合主解，只要 `κ` 仍是固定非零常数，核心奇异面都不会从 `X=0` 挪走

### 2026-04-29

- 用户追问：既然物理问题已经固定为平直背景上两束分离高斯波包的干涉，为什么还会出现“初始几何”问题
- 当前统一后的回答是：
  - 若只把 `(\rho,S)` 当作给定物质初值、并把背景几何固定为平直 `g`，当然可以直接取“平直时空 + 两束高斯波包”作为初始条件
  - 但一旦任务改成“三种引力作用量下的全动力学演化”，就不再是固定背景问题，而是 Einstein 型约束—演化问题
  - 在这种情况下，初始空间度规与外曲率不能任意给定；它们必须与给定的 `\rho,S` 一起满足 Hamilton 约束和动量约束
  - 因而“平直几何 + 非零束流物质”一般只是弱场近似启动，不是严格的全动力学约束解
- 这条区分已成为当前数值主线的基本口径：
  - 固定背景 / 近似启动：可直接取平直初值
  - 严格全动力学：必须单独处理初始几何
- 用户进一步给出当前参数区间下的工作判断：
  - 对现实单粒子或弱物质量级，`A` 支的能动量回授很弱，`g` 可近似视为平直或固定背景
  - 因此 `A` 支可先作为估算 `(\rho,S,g)` 初值的参考
  - 真正值得投入主要全动力学计算精力的是 `B/C` 支，因为 `g~` 不必接近固定背景
- 当前对此达成的方法论共识：
  - 不否认 `A` 支的严格全动力学定义
  - 但在当前现实弱场目标下，数值工作优先次序改为：`A` 近似参考，`B/C` 重点全动力学
- 按这条新顺序，新建：
  `kg_examples/simulate_bc_geometry_from_a_reference.py`
- 第一版实现中，先错误地把“参考 `A`”仍写成旧的 `(\rho,S)` 近似推进器；短时间测试显示参考能量密度会迅速抬升到 `O(10^4)`，不适合作为新的主参考
- 已立即修正：当前参考 `A` 支改为真正的平直背景复 Klein-Gordon 演化，再由 `\psi=\phi_1+i\phi_2` 重构 `(\rho,S)`
- 修正后，参考 `A` 支的最大能量密度在短时间内从 `867.97` 缓慢下降到 `~861`，不再出现先前的异常抬升
- 同时引入线性化弱场误差估计：
  `\Delta \omega = - \mathcal E_A / M_P^2`
  得到当前 `A` 支近似平直背景的误差量级约为 `4e-3`
- 在同一参考物质驱动下，`B/C` 支的几何偏离量级已达到 `10^{-1}`，明显大于 `A` 支平直近似误差
- 但当前 `B/C` 的几何时间推进仍有明显时间步敏感性：
  同一物理时长下，`max|h_xx-1|` 约从 `0.452` 变到 `0.700`
- 已把这一轮方法与结果整理成：
  `research-notes/082-A支平直参考驱动BC几何的第一轮结果.md`
  - 原因是恢复 HJ 和连续性方程所必需的关键组合始终是 `κ/X`（massive-like 分支）或 `κ/Q`（massless patch 分支）
  - 因而 `κ` 改变的首先是 blow-up 的整体尺度与某些“无修正面/额外退化面”的位置，而不是那个由原 null 极限 `X=0` 决定的主奇异面
  - 例如在 `C=1` 的投影型主解里，`s = κ/X - 1`；因此 `s=0` 的“无修正面”会从 `X=m^2` 挪到 `X=κ`，但 `X=0` 的极点位置不变
- 用户进一步提出更细的想法：是否能把 `X=0` “平移”成 `g~` 下的正长度，从而避免行列式 blow-up
- 当前数学化后的结论：
  - 在当前所有“连续性恢复干净”的变换族中，都有
    `g~^(μν)u_ν = A u^μ`
    从而
    `g~^(μν)u_μu_ν = A X`
  - 若想把原 `X=0` 映成 `g~` 下的有限正长度 `L^2 > 0`，就等价于要求
    `A X -> L^2`
    因而必有
    `A ~ L^2 / X`
    发散
  - 这说明：在当前这类保持流方向不变、从而能干净恢复连续性方程的变换族里，“null -> timelike finite shift”本身就强制 `1/X` 型奇异
  - 若进一步要求 determinant 不发散，可以通过让 transverse 因子 `C` 同时趋零来补偿；但这只会把问题改写成横向缩退或系数不光滑，而不是得到一个在 `X=0` 处真正 regular 的同类图册
  - 若真想避免这个结论，必须放弃当前关键结构之一：
    1. 不再要求 `g~^(μν)u_ν ∥ u^μ`
    2. 或不再要求 `g~^(μν)u_μu_ν` 在 `X=0` 处取非零固定值
    3. 或改用新的 null-shell / 辅助场图册
- 用户继续追问：能否让 `ρ~` 在 `X=0` 时等于 0，从而吸收 `u_μ` 模长为 0 的性质
- 已完成系统检查并整理成：
  `research-notes/067-在X等于0处让rho波浪吸收奇异性的可能性.md`
- 当前结论：
  - 连续性方程真正要求匹配的是电流
    `J~^μ = sqrt(|g~|)ρ~ g~^(μν)u_ν`
    而不是单独的 `g~^(μν)u_ν`
  - 但只要原电流非零，点态匹配自动推出 `g~^(μν)u_ν ∥ u^μ`；所以“先乘 `ρ~` 再平行”在非零电流区域并没有真正放宽方向条件
  - 一般公式是
    `ρ~ = ρ sqrt(R) / A`
    其中 `A` 由 `g~^(μν)u_ν = A u^μ` 定义，`R = det(g~^(..))/det(g^(..))`
  - 在 rank-1 `C=1` 分支里，`ρ~ ~ sqrt(X)`，确实会在 `X -> 0` 时趋于 0
  - 在投影型混合 `C=1` 主解里，`ρ~ = ρ`，并不会趋于 0；这里电流的 regularization 主要由体积元吸收
  - 在 pure conformal 里，`ρ~ ~ 1/X`，反而会发散
  - 因此，让 `ρ~ -> 0` 的确可能帮助“电流级 regularization”，但它不能消除把 `X=0` 映成非零壳长所必需的 `A ~ 1/X` 主奇异
- 用户进一步质疑：对混合型的分析里，不应先假定 `A ~ 1/X`，因为控制组合也许可以 `1/X^2` 甚至更高阶发散
- 已完成严格局域分析并整理成：
  `research-notes/068-混合ansatz在X等于0附近的严格局域分析.md`
- 当前严格结论：
  - 在一个 `X≠0` 的 punctured neighborhood 里，若
    1. `u^μ` 与 `r^μ` 线性无关
    2. 原电流非零
    3. 要求点态恢复原连续性方程
  - 则自动推出
    `EX + FY = 0`
    从而
    `g~^(μν)u_ν = A_u u^μ`
    其中
    `A_u = C + DX + EY`
  - 再要求变换后 HJ 壳长固定为 `κ ≠ 0`，则精确得到
    `A_u = κ/X`
  - 因而：
    - 单个系数 `C,D,E,F` 可以比 `1/X` 发散得更快
    - 但真正受两条方程共同约束的主组合不可能是 `1/X^2`，而只能是 `κ/X`
- 用户追问 Bekenstein 变换的原则，以及这是否会约束 `κ`
- 当前判断：
  - Bekenstein 原始路线的核心不是“必须把光锥保持完全不变”，而是：在弱等效原理、协变性和因果性约束下，物质度规与引力度规之间的关系应落在 disformal 这一类里
  - 纯 conformal 是这个框架下的特例；只有 pure conformal 才严格保持同一组 null cone
  - 一般 disformal 会改变光锥开张程度，因此不自动保持完全相同的因果序；真正需要的是保持 Lorentzian 性、时间定向，以及在所选物理解释下不产生不可接受的双锥冲突
  - 对当前几何化常数 `κ`，Bekenstein 型因果原则给出的首先是符号/锥结构约束，而不是数值必须等于 `m^2`
  - `κ = m^2` 真正来自另一层要求：在 `m≠0` 主分支上希望 `Q->0` 时恢复 `g~->g`，并把变换后壳约束规范化成原来的质量壳
- 用户进一步指出一个关键直觉：
  - 对投影型主解 `D=Z, E=-Y, F=X, C=1`，行列式比值是 `(m^2/X)^2`，因此 `X=0` 时的 blow-up 更像“把原 `g` 下变成 null 的方向映到 `g~` 下固定非零壳长”所需的奇异重标定，而不必自动理解成度规退化
- 当前修正判断：
  - 这条直觉的核心是对的：对当前这类满足 `g~^(μν)u_μu_ν = κ ≠ 0` 的变换族，若原壳量 `X -> 0`，则为了保持变换后壳长固定非零，系数如 `κ/X` 必然 blow up；这是“非可逆/奇异 disformal 极限”的典型形态
  - 这与文献里对 `A+BX=0` 一类 disformal singular surface 的处理一致：若还想把它当成等价场重定义，就必须避开该面；若把该极限当真，就往往进入 non-invertible / mimetic 分支
  - 但“行列式比值恒正”仍不足以推出“signature 一定完全不变”；要保证保持同一 Lorentzian 分支，仍需额外检查 `u-r` 平面的惯性结构，而不能只看 determinant

### 2026-04-21

- 初始化 `agent-memory/` 目录
- 创建核心文件：`README.md`、`TASKS.md`、`LOG.md`、`DECISIONS.md`
- 根据用户要求将长记忆状态设为 `ON`
- 项目目录当前未发现可识别的源码、文档或配置文件
- 已将长期记忆的开关语义并入 skill：支持“长记忆开/关”“agent memory on/off”“am on/off”等指令
- 用户给出课题“量子势几何化”的初步研究设想：以 KG 标量场为最小模型，尝试通过 Bekenstein disformal transform 将量子势吸收到新度规中
- 已根据用户描述整理出第一版树状研究计划
- 当前推进重点切换为：先严格建立最小理论设定与成功判据，再进入具体变换构造
- 已创建 `research-notes/` 作为独立于 `agent-memory/` 的人工查阅留档目录
- 已完成第一份阶段文档：`001-KG作用量的Madelung分解与量子势定位.md`
- 本轮得到关键结论：量子势在原始作用量中编码为振幅梯度项，在分部积分后的等价作用量中可显式写成 `\rho Q`
- 已完成第二份阶段文档：`002-量子势几何化的成功失败与可接受判据.md`
- 已把后续 ansatz 的检查标准固定为：协变性、可逆性、作用量成功、方程等价、轨迹测地线化、经典极限、解释增益
- 已完成第三份阶段文档：`003-disformal-ansatz候选与第一轮筛选.md`
- 当前筛选结果：纯 conformal 仅作基线；`u_\mu` 单向 disformal 与 `n_\mu` 单向 disformal 保留；混合双向 disformal 为主攻路线；显含二阶导数的度规暂缓
- 用户新增统一记号约定：`X=g^{\mu\nu}\partial_\mu S \partial_\nu S`，`Y=g^{\mu\nu}\partial_\mu S \partial_\nu \sqrt{\rho}`，`Z=g^{\mu\nu}\partial_\mu \sqrt{\rho} \partial_\nu \sqrt{\rho}`
- 已将候选变换中的振幅方向记号统一为 `r_\mu=\partial_\mu\sqrt{\rho}`
- 已完成第四份阶段文档：`004-四类度规变换的初步比较分析.md`
- 已根据逻辑需要修正计划：把共形、相位单向、振幅单向、混合四类放入统一比较框架，不再只分析混合路线
- 用户要求直接落实到计算，不要盲目否定任何方案
- 已完成第五份阶段文档：`005-四类变换代回作用量与HJ方程的直接计算.md`
- 当前计算结论：四类方案都没有被直接算掉；共形、相位单向、振幅单向都能在某个层面吸收部分结构，混合方案则是最小的同时含 `X,Y,Z` 的方案
- 用户进一步提出问题：是否能证明存在更广泛的一类变换都可以实现这一目的
- 已完成第六份阶段文档：`006-HJ层面的吸收Q求解与更广泛存在性.md`
- 当前新增结论：在 HJ 层面，吸收 `Q` 对应的不是孤立解，而是一个更广泛的变换族；但这还不是完整几何化的最终存在性证明
- 用户进一步强调完整筛选标准：作用量兼容、场方程等价、测地线、经典极限、可逆性、因果结构
- 用户指出一个关键观察：包含 `r_mu` 的方案未必总有解，类 II 在 `A=1` gauge 下则给出唯一解
- 已完成第七份阶段文档：`007-快速筛选：作用量兼容性与存在条件.md`
- 当前筛选结论：共形最简单解与类 II 的 `A=1` 最简单解都只达到 on-shell 作用量兼容；类 III 为判别式条件存在；类 IV 为 `Delta ≠ 0` 型存在
- 用户指出关键方法论修正：物理上不要求作用量逐项等价，只要求作用量导出的物理方程等价即可
- 已完成第八份阶段文档：`008-方法论修正：作用量等价不是必要条件.md`
- 当前推进标准已修正：后续优先检查变分方程、测地线与物理等价，不再把 off-shell 作用量逐项匹配当作硬门槛
- 已完成第九份阶段文档：`009-变分筛选：最简单共形与类II方案.md`
- 当前新增结论：最简单共形方案在当前几何作用量下破坏原连续性方程；类 II 的 `A=1` 最简单解会使作用量恒等为零
- 已识别一个结构性约束：若通过选系数让 HJ 条件在 integrand 中恒等成立，则几何作用量会退化，不能这样构造完整理论
- 用户提出关键修正：`ρ~` 也应随度规变换，并满足 `sqrt(-g)ρ = sqrt(-g~)ρ~`
- 已完成第十份阶段文档：`010-密度随度规变换后的重新检验.md`
- 修正后结论：类 II 的 `A=1` 最简单解仍然退化为零作用量；最简单共形方案的问题仍然在于连续性方程不对，只是原因定位到 `g~^(mu nu)` 而不是体积因子
- 用户继续指出关键方法错误：不能先把变换代回作用量做判断，而应先对新作用量变分，再代回变换检查
- 已完成第十一份阶段文档：`011-方法修正后：先变分再代回变换.md`
- 已撤回上一轮对类 II 最简单解的直接排除：按正确顺序计算后，它给出的是修改后的连续性方程，而不是零动力学
- 已完成第十二份阶段文档：`012-一般类II能否同时满足HJ与连续性方程.md`
- 当前新增结论：在当前最简几何作用量 ansatz 下，一般类 II 一旦满足原 HJ 方程，其连续性方程就自动变成加权流形式，因此泛型地不能再现原连续性方程
- 用户怀疑类 II 的 `A=1` 计算有误，并指出自己算到的 `B` 含 `Z`
- 已完成第十三份阶段文档：`013-类II中A=1时B是否必然含Z的复核.md`
- 当前复核结论：在纯类 II 且只由 HJ 方程决定 `B` 的最小问题里，`B` 不被强制含 `Z`；若 `Z` 出现，说明引入了额外条件或已离开纯类 II
- 用户重新明确当前要做的核心任务：把 `g~^(mu nu)=C g^(mu nu)+D u^mu u^nu` 直接代回拉氏量，要求把 `Q` 项消去
- 已完成第十四份阶段文档：`014-按逆度规ansatz直接消去拉氏量中的Q项.md`
- 当前新增结论：在 `sqrt(-g)ρ = sqrt(-g~)ρ~` 约定下，直接消去 `Q` 的核心条件是 `D X^2 + (C-1)X - Q = 0`；最简单解是 `C=1, D=Q/X^2`
- 已完成第十五份阶段文档：`015-最简单拉氏量消项解的后续检验.md`
- 当前新增结论：最简单消项解 `g~^(mu nu)=g^(mu nu)+(Q/X^2)u^mu u^nu` 完美重现 HJ 方程，并在 HJ 层面支持测地线解释，但连续性方程仍带权重 `(1+Q/X)`，因此守恒结构未完成
- 用户指出连续性方程检查还不够细，要求把 `X` 的具体表达式代进去再算
- 已完成第十六份阶段文档：`016-连续性方程中代入X=m2-Q后的精确比较.md`
- 当前新增结论：把 `X=m^2-Q` 代入后，新旧连续性方程重合的精确条件是 `J^mu ∂_mu Q = 0`
- 用户进一步追问：`J^mu ∂_mu Q = 0` 是否与新几何下的测地线方程等价
- 已完成第十七份阶段文档：`017-J^mu∂_muQ与新几何测地线方程是否等价.md`
- 当前新增结论：两者不等价；测地线方程由新 HJ 方程自动导出，而 `J^mu ∂_mu Q = 0` 是让新旧连续性方程重合的额外守恒条件
- 已完成第十八份阶段文档：`018-最小混合逆度规下的连续性方程no-go.md`
- 当前新增结论：在线性最小混合逆度规 ansatz 下，不可能同时实现“拉氏量精确消项”和“原连续性方程逐字重合”，除非 `Q=0`
- 用户指出守恒方程变了意味着波函数本身也应随视角变化，并要求去掉 `sqrt(-g)ρ` 保持不变的假设重新考虑
- 已完成第十九份阶段文档：`019-去掉sqrt-g-rho不变假设后的纯u-ansatz翻盘.md`
- 当前新增结论：放开 `ρ~` 变换后，只要满足 `sqrt(-g~)ρ~ = sqrt(-g)ρ / (C + D X)`，纯 `u^mu u^nu` ansatz 就足以同时重现原 HJ 方程和连续性方程
- 已完成第二十份阶段文档：`020-最简单解的可逆性签名因果结构与经典极限.md`
- 当前新增结论：最简单纯 ansatz 的物理适用域是 `X = m^2 - Q > 0`；在该区域内它同时满足可逆性、Lorentz 签名保持、实幅度与正确经典极限
- 用户提出关键判断：`X < 0` 应对应其他变换分支
- 已完成第二十一份阶段文档：`021-X正负对应的纯ansatz双分支结构.md`
- 当前新增结论：纯 ansatz 天然分成 `X>0` 的 timelike branch 与 `X<0` 的 spacelike branch；若不切换分支，就会误把 `X<0` 误判成几何失效
- 用户给出 Cursor 对本项目的对照性总结，要求检查对方或我方的错误
- 已对照 Cursor 的 `DERIV-0015`、`DERIV-0016` 与我方 `014/019/020/021`
- 当前对照结论：Cursor 的核心 rank-1 + 共形严格匹配思路没有明显代数错误；我方前期真正需要撤回的是基于“密度固定不变”得到的若干连续性失败判断；两边最大的差异是采用了不同的等价性定义（Cursor 保持 `(R,S)` 不变并做严格 action-level matching，我方后期允许 `ρ~` 真正变换）
- 已完成第二十二份阶段文档：`022-对Cursor对照分析的判断.md`
- 用户要求做阶段性总结，并继续研究与 Einstein 作用量结合的部分
- 已完成第二十三份阶段文档：`023-阶段性总结与Einstein作用量耦合框架.md`
- 当前新增结论：`S_EH[g]` 与 `S_EH[g~]` 都应保留为候选；在弱量子势、近似平直背景下，两者差异主项被强烈压低，而在“强引力 + 强量子势”区域差异才会放大
- 用户补充提醒：即使背景 `g` 小，若量子势 `Q` 很大，则 `g~` 仍可能很大，不能简单用“近似平直背景”推出差异小
- 已完成第二十四份阶段文档：`024-S_EH[g]与S_EH[g~]的弱展开与差异量级.md`
- 当前新增结论：两种 Einstein 作用量差异是否小，不只取决于背景引力是否弱，还取决于 `ε=Q/m^2` 及其导数是否小；弱引力但强/快变量子势也可能带来显著差异
- 用户进一步追问：既然 `Q` 可以很大，那么这种情况下如何恢复到实验预测
- 当前阶段判断：这不是自动成立的，反而是当前理论最关键的一致性检验；若在已观测弱场系统中 `Q/m^2` 大且变化快，却又不导致可观测偏差，则必须依赖额外机制（例如只在不可观测区域出现、仅常数大而导数小、普适耦合到同一几何、或某种筛选/屏蔽机制），否则理论会失败
- 用户要求先建立一个 1+1d 双高斯干涉自由 KG 基线模型，并要求模块化，便于后续放入 Schwarzschild / 引力波背景
- 已创建脚本：`kg_examples/kg_double_packet.py`
- 已完成第二十五份阶段文档：`025-1+1d自由KG双高斯波包基线模型.md` 与 `026-1+1d双高斯KG基线算例结果.md`
- 当前新增结论：在当前基线参数组下，主支撑区内 `Q/m^2` 已达 `10^-1 ~ 1` 量级，而超过 1 的大尖峰主要集中在低密度、接近节点的区域
- 用户澄清双高斯例子的真正用途：用来判断“四种一致性机制”哪种可行、哪种应否定，而不是停在抽象猜想
- 已完成第二十六份阶段文档：`027-用双高斯KG例子检验四种一致性机制.md`
- 当前新增结论：在这个具体模型上，机制 1（大 Q 根本不出现）与机制 2（Q 虽大但近常数）一般不成立；机制 3（普适耦合）暂不能否定；机制 4（屏蔽/筛选）是当前最值得继续深推的路线
- 用户提出新的思路：考虑 `F(Q)R~` 或 `ρ~R~` 形式的引力作用量，例如 `exp(-(κQ)^2)R~`、`1/(1+(κQ)^2)R~`
- 已完成第二十七份阶段文档：`028-Q依赖引力作用量的一般要求与例子分析.md`
- 当前新增结论：若以当前双高斯模型作为第一轮检验，最值得优先深推的是 `F(Q)=1/(1+(κQ)^2)` 这一类平滑、有界、正定的屏蔽型修正；纯 `ρ~R~` 不适合作为主引力项
- 已完成第二十八份阶段文档：`029-双高斯模型上测试FQ=1-1pluskQ2.md`
- 当前新增结论：在双高斯模型上，`F(Q)=1/(1+(κQ)^2)` 在主支撑区通常较温和，但节点附近的 `∇F`、`∇∇F` 会迅速增大；因此真正的危险量往往是导数项而不是 `F` 本身
- 用户指出自己真正要的是：在这个 1+1d 双高斯模型下直接比较两种总作用量预测出的差异
- 已完成数值脚本：`kg_examples/kg_compare_total_actions_1p1.py`
- 已完成第三十一份阶段文档：`031-1+1d双高斯模型下两种总作用量的直接比较结果.md`
- 当前新增结论：在 1+1d 平直背景里，`EH[g]` 与 `EH[g~]` 的局域几何密度差异可以很大，但纯 EH 项不会自动给出 3+1d 意义下的新的 bulk 动力学预测
- 用户指出当前 3+1d 数值计算导致后台卡顿，希望控制内存占用
- 当前诊断：卡顿主要来自高维广播数组和高分辨率二维切片的瞬时内存压力，而不是单步数学错误；后续 3+1d 计算应改为按线、按块、按需求值

### 2026-04-24

- 用户要求切到 3+1d 并多次校正几何：先后尝试了屏幕双缝模型、对射高斯包模型、准单色双缝束流模型、3+1d 自由 KG 精确 Fourier 对射模型
- 用户指出多次数值/建模问题：单束伪多峰、几何配置与描述不符、后台卡顿、文件未真正生成、近似束流模型不适合正式计算 Bohm 量
- 用户明确否定此前的 reduced closure 联合演化结果，要求从完整作用量严格推导动力学方程再数值积分；这条要求已接受
- 已完成文档：`032`、`033`、`034`、`035`、`036`、`037`、`038`、`039`、`040`、`041`、`042`
- 当前最新主线：A/B/C 三种完整作用量已统一到一个约束框架中；C 分支必须先做辅助场重写；后续不再直接 brute-force 3+1d，而是在 2+1d 一致约化框架下推进

### 2026-04-25

- 用户多次显式发出 `agent memory on!!!`
- 长记忆状态继续保持为 `ON`
- 当前线程继续视为长记忆模式下的工作线程，不重置项目上下文
- 用户再次显式发出 `agent memory on!`
- 长记忆状态再次确认维持为 `ON`
- 用户追问：`2+1d` 约化的动力学本质能否与 `3+1d` 一样，以及是否意味着 `1+1d` 空间不存在引力。
- 当前判断：`2+1d` 一致约化只能保留带一个 Killing 对称的 `3+1d` 子扇区，不等于完整 `3+1d`；纯 `1+1d` Einstein-Hilbert 项没有局域度规动力学，但从 `3+1d` 一致约化到 `1+1d` 并保留横向尺度/形状场时，会得到二维 dilaton/标量-张量型引力系统，仍可作为最小可控模型研究三种作用量的动力学差异。
- 用户最终决定放弃 `1+1d`，回到 `2+1d` 研究 `x-z` 平面中 `±45°` 两束干涉，并要求从约化后的完整动力学方程出发做全动力学数值模拟。
- 已完成新文档：`044-2+1d下三分支动力学方程与数值闭合.md`
- 已创建新脚本：`kg_examples/simulate_2p1_full_dynamics.py`
- 当前数值器实际跑通的三支是：`flat reference`、`EH[g~]`、`F(R~)R~`
- 当前结论：`EH[g~]` 分支在第一组参数下只对平直参考支产生很小偏离；`F(R~)R~` 分支会出现更强回授，且比 `EH[g~]` stiff 得多
- 当前明确剩余问题：严格意义下的 A 分支还缺一个明确的 `2+1d reduced closure` 来把 `g` 与 `g~` 关起来；在未补上这一点前，数值里的 `flat reference` 只是参考支，不是完整 A 分支
- 对照用户认可的 `visualizations/double_slit_beamlike_xz/alpha_0.5_xz_linear.png` 后确认：当前 `simulate_2p1_full_dynamics.py` 的初值几何又偏回了“两个局域团块”，没有保持此前认可的 `z=±x` 两束相交几何，因此这轮 `2+1d` 数值结果只能看成方程原型测试，不能直接当作目标几何下的物理解。
- 用户进一步指出：A 分支“就是不用 Bekenstein 变换的那一支”，因此不应在 A 分支里同时出现 `g` 与 `g~` 两套几何变量。这说明此前对三分支的统一写法里，A 分支的物理含义被我混淆了。
- 用户查看了当前 `current_geometry_check` 下的 `x-z` 截面 `rho` 图后，认可这批静态束流几何基本没问题，并追问为何先前“全动力学”结果里看不到干涉条纹。
- 当前判断：先前全动力学代码求解的是一个自由 Cauchy 初值问题，而用户认可的束流图对应的是持续相干入射/准单色束流图样；两者不是同一个问题，因此条纹不会自动在后者那样保留下来。
- 已创建新脚本：`kg_examples/simulate_2p1_boundary_driven.py`
- 当前边界驱动结果：A 分支（原始平直 KG）已经出现与用户直觉一致的 `±45°` 束交叉和中心条纹；B 分支保留类似条纹但回授后峰值更高；C 分支在温和参数下不再立刻爆炸，但中央增强明显更强，仍显示出较强 stiffness
- 已按用户要求把边界驱动改成“只从 `z<0` 半平面向 `z>0` 入射”，同时提高分辨率并把 `F(R~)R~` 的线宽收窄。
- 当前新现象：A 分支会沿 `±45°` 方向向上推进；B/C 分支在当前边界实现下却明显失去向上传播，显示出边界条件与几何化方程不相容的迹象。当前判断是不应把它直接解释成作用量的物理预言。
- 用户进一步指出：A 分支与 `current_geometry_check` 的解析束流图样明显不一致。当前诊断：这是因为 A 分支的边界源仍是手写近似入射波，而 `current_geometry_check` 是整域解析束流解；两者不是同一个边界值问题。
- 已新增脚本：`kg_examples/check_benchmark_residual.py`
- 残差检查结论：当前 `beamlike` benchmark 代回平直 KG/Helmholtz 后，在主支撑区的平均相对残差约为 `3e-2`，局部最大可达 `O(1)`；因此它不是 full KG 时间推进器应收敛到的精确目标解
- 当前更新判断：若要维持用户认可的几何图样并进行可信的三分支比较，最自然的下一步是改写成稳态边界值问题，而不是继续拿超曲型时间推进器逼近该 benchmark
- 已创建新脚本：`kg_examples/exact_kg_2p1_benchmark.py`
- 已生成真正的 `2+1d` free KG 精确 benchmark，输出在 `visualizations/exact_kg_2p1_benchmark/`
- 当前结论：精确 benchmark 的单束是沿 `±45°` 方向的局域高斯包在重叠时刻的截面，双束则给出局域干涉条纹；它与此前 `beamlike` 图样属于不同的模型层级
- 用户要求把传播方向（`±45°` 传播轴）的动量展宽压到 `0.01 k0` 量级。
- 已将 `exact_kg_2p1_benchmark.py` 改成各向异性动量谱：`sigma_k_parallel=0.1, sigma_k_perp=1.0`（对应 `sigma_parallel/k0=0.01, sigma_perp/k0=0.1`），并重生成精确 benchmark。
- 当前新结果：单束已经明显沿 `±45°` 方向拉长，双束在原点附近形成更接近用户直觉的 `X` 形交叉与中心条纹；但由于束在传播方向上更长，当前 `[-5,5]` 观察窗的边界不再足够干净。
- 已把精确 benchmark 观察窗扩大到 `[-10,10]`，重新生成 `visualizations/exact_kg_2p1_benchmark/`
- 已创建新脚本：`kg_examples/calibrate_exact_kg_solver_2p1.py`
- 当前 A 分支快速校准结果：在粗时间步 `dt≈0.203`、`61x61` 网格下，对精确 benchmark 的相对 `L1` 误差约为 `7.9%`，相对 `L∞` 误差约为 `96%`；说明当前 leapfrog + 直接 Dirichlet 边界的时间推进器还没有调到可信比较所需的精度
- 已进一步把校准改成“边界缓存 + 多个时间步长”的收敛测试，避免每一步都重算整张精确图
- 当前校准结论：`dt≈0.203 -> 0.0997 -> 0.0499` 时，相对 `L1` 误差从 `7.9% -> 8.5% -> 9.9%`，没有收敛迹象
- 直接结论：当前 `leapfrog + 全边界 Dirichlet` 的 A 分支数值器不可信，问题不在时间步不够小，而在边界实现/离散方法本身
- 已尝试把 A 分支改写为 Helmholtz 边界值问题（`solve_a_helmholtz_2p1.py`），结果表明当前精确 benchmark 由于是波包而非单频场，不适合直接用单一 Helmholtz 方程拟合
- 已新增 `calibrate_exact_kg_fft_2p1.py`，用 FFT 正频传播去校准 A 分支
- 当前 FFT 校准结果：对 `alpha=0.5`，中心条纹几乎与精确 benchmark 重合，误差主要集中在边缘/下半部；相对 `L1` 约 `0.8%`，相对 `L∞` 约 `17%`
- 已新增 `calibrate_exact_kg_ps_fd_absorbing.py`，把 A 分支改成“谱 Laplacian + 吸收层 + 精确边界”的混合校准器
- 当前最新结果：内核区 `[-10,10]^2` 上，相对 `L1 ≈ 0.22%`，相对 `L∞ ≈ 4.69%`
- 当前判断：A 分支终于有了一个足够可信的 PDE 数值基线，后续可以开始把同样的空间离散/边界思想平移到 B/C
- 已新增 `simulate_bc_psfd_imex.py`，把 B/C 改成 IMEX/半隐式时间推进，同时保持与 A 相同的谱空间离散、吸收层和 benchmark 边界准备
- IMEX 版结果：B/C 至少可以在默认参数下稳定跑完，但与 A 的偏离极大
- 弱耦合回归测试（大 `M_P`、小 `ell`）仍未出现 `B,C -> A` 的趋势
- 当前结论：B/C 的 IMEX 版虽然数值上比显式版稳，但离散实现仍然不可信；下一步应回到“B 在弱耦合极限必须退回 A”这一 sanity check 上定位离散 bug
- 已新增 `scan_bc_weak_coupling.py`，对 IMEX 版做系统弱耦合扫描
- 扫描结果：`M_P=300 -> 1000 -> 3000 -> 10000`、`ell=0.02 -> 0.01 -> 0.005 -> 0.001` 时，`B_vs_A_relL1` 与 `C_vs_A_relL1` 都没有单调减小，且量级远大于 1
- 当前最强判断：B/C 当前离散实现不是“还没进入弱耦合区”，而是“方程实现或边界耦合本身仍有 bug”
- 已新增 `debug_b_boundary_compatibility.py` 和研究笔记 `047-B分支平直极限下的边界不相容性.md`
- 关键 debug 结果：即使把 B 的几何完全冻平，只保留其平直物质方程，并继续喂 A benchmark 的边界 `ρ,S`，依然得到 `Bflat_vs_A_relL1 ≈ 2.7e3`
- 当前判断：B/C 与 A 的严重不一致，至少有一大块原因来自“同一组 A 边界数据对 B/C 物质系统并不相容”；因此下一步需要改成三支共享 incoming characteristic data，而不是直接共享整个 `ρ,S` 边界图样
- 已新增 `simulate_abc_cauchy_same_initial.py`，尝试用“同一初始切片”而非共享整段边界历史来比较 A/B/C
- 当前结果：A 的 Cauchy 数值器在这版脚本中尚未校准好，B/C 则迅速塌成中心竖直尖峰；说明问题已经不只是边界历史不相容，而是 B/C 的离散实现本身还没有通过基本的 Cauchy sanity check
- 当前最新判断：后续 debug 应优先围绕 B 的平直极限和无几何回授极限，逐步恢复各项；不要再直接把完整 B/C 方程一股脑推进到长时间
- 用户提出新的比较入口：若 C 和 B 相似，则可以对 C 取很小的 `F(x)` 线宽，并把它作为与 A 比较的入口
- 已新增 `scan_c_small_ell_local.py`，专门扫描 C 分支小 `ell` 极限并分别统计全域与低曲率代理区（`|phi-1|<1e-2`）相对 A 的偏离
- 当前扫描结果：在现有实现里，`ell` 不是单调越小越接近 A；反而在 `ell≈1e-2` 左右偏离最小，而更小的 `ell` 会再次变差
- 当前判断：这一步说明“用小 `ell` 的 C 作为比较入口”在原则上是值得尝试的，但当前 C 的离散实现仍然掺杂较强数值效应，不能直接把扫描曲线当物理结论
- 已补充输出每个 `ell` 对应的 `rho_C` 图和 `|C-A|` 图，便于人工直接查看各线宽下 C 路的空间分布
- 用户要求把 C 的 `rho` 显示压到和 A 同量级，避免中心超大峰值毁掉整张图的可读性
- 已为每个 `ell` 补充输出一套以 `A` 的 `rho_max` 为统一色标上限的 `rho_C` 图（文件名含 `_matched_to_A`）
- 已新增 `debug_b_phase_gradients.py`，比较 B 平直物质系统中两种相位梯度定义（`unwrap(S)` 差分 vs `j/ρ`）
- 当前结果：`j/ρ` 版本比直接对 `unwrap(S)` 求导稳定一些，但两者都远离 A 的后续演化
- 当前判断：这一步支持了一个更强的认识——问题不只是相位梯度数值实现，而是 B 的平直物质系统从同一初始切片出发本来就不应期待跟 A 保持相同后续演化
- 已对单束情形 `alpha=0` / `alpha=1` 做最小测试：即使没有干涉节点，`Bflat` 也会形成异常尖峰，`Bfull` 与 `Bflat` 几乎重合，而 `Cfull` 更早塌到更大峰值
- 当前判断：这说明问题不只是“干涉节点/相位分支跳跃”引起的病态；至少在当前实现里，B/C 物质更新模块本身就有严重问题
- 已新增笔记：`048-Bflat中的caustic迹象.md`
- 新增证据：在 `Bflat` 中把空间分辨率从 `81 -> 121 -> 161` 提高时，中心轴附近的峰值迅速增大（约 `1.1e3 -> 5.1e6 -> 2.8e7`），显示出 caustic/特征线交叉迹象
- 当前判断：Bflat 的中心尖峰不应再简单视为普通数值 bug；后续比较应区分平滑区与 caustic 区
- 用户提出新的比较入口：对 C 分支取很小的 `F(x)` 线宽，使其在大多数低量子势/低曲率区域与 A 接近，仅在少数高曲率区偏离。
- 当前判断：这条路是合理的，而且比强行要求 B 贴近 A 更自然；后续可把 C 当作“局域屏蔽修正”去与 A 比较，并重点看偏离是否局限在高曲率小区域。
- 用户进一步追问：既然实验上可以在初始时刻准备出相同的 `ρ,S,g`，为何还说对 B/C 不相容。
- 当前修正判断：`ρ,S,g` 在单个初始时刻完全可以取和 A 相同；真正不相容的是把 A 的整套边界/时间历史（尤其全边界驱动的 `ρ,S`）原样搬给 B/C。后续必须区分“相同初始切片”与“相同完整边界-时间数据”。
- 已新增 `simulate_bc_psfd_absorbing.py`，尝试把 B/C 平移到与 A 同样的“谱 Laplacian + 吸收层 + benchmark 边界”骨架
- 当前测试：无论是全边界 benchmark 喂法还是只保留底边注入，B/C 在显式推进下仍迅速进入 overflow/NaN
- 当前最新判断：B/C 的主要问题已收缩为“几何场方程在当前显式时间推进下过 stiff”，后续应转向更稳定的时间积分（例如隐式/半隐式）或先做准静态/约束求解，而不是继续在同一显式骨架上蛮调参数
- 已把精确 benchmark 的波包中心移到 `x0=z0=8`，使 `[-10,10]` 观察窗与新的入射位置匹配
- 已新增 `calibrate_exact_kg_fd_absorbing.py`，测试“大域 + 吸收层”的 A 分支 PDE 求解
- 当前 `fd+sponge` 校准结果：在内核区 `[-10,10]^2` 上，相对 `L1` 约 `3.0%`，相对 `L∞` 约 `78%`；主几何轮廓与中心条纹已经对上，但中心峰和若干斜向结构仍有明显偏差
- 当前判断：A 分支现在已有两个可用校准参考，`FFT` 在中心区更贴 benchmark，`fd+sponge` 更像真正 PDE 求解器；后续若要把方法平移到 B/C，应优先沿 `fd+sponge` 这条线继续压误差
- 已继续改进 `fd+sponge` A 分支校准器：把空间 Laplacian 的 bulk 升到四阶，同时保留边缘二阶回退
- 当前改进后结果：内核区 `[-10,10]^2` 的相对 `L1` 误差从约 `3.0%` 降到约 `1.17%`，相对 `L∞` 从约 `78%` 降到约 `25.6%`
- 当前判断：A 分支的 PDE 求解器终于达到“可以作为 B/C 平移模板”的阶段；下一步应把同样的大域 + 吸收层 + 四阶空间离散骨架移到 B/C
### 2026-04-26

- Added `kg_examples/debug_b_matter_solver.py` to isolate the flat B matter update on a single-beam exact KG slice.
- Found that the old matter update was the first real numerical failure point:
  - `current` spectral-current update gave `total_n_rel_drift ~ 1.38e3`
  - `localgrad` still drifted by `O(1)`
  - a new finite-volume Rusanov update reduced the drift to about `1.1e-2`
- Moved the conservative Rusanov matter flux into `kg_examples/simulate_bc_psfd_imex.py`.
- Re-ran `debug_b_flat_limit.py`:
  - the no-sponge weak-gravity test is now essentially exact again
  - the old failure mode is strongly tied to damping/sponge being applied to matter transport variables
- Re-ran `debug_c_decomposition.py` after the matter update fix:
  - `Bflat` no longer shows the fake ultra-large spike; it now stays on the same scale as `Bfull`
  - this lowers confidence in the earlier interpretation that the old `Bflat` spike was mainly a physical caustic
- Compared `rhs_branch_c(phi=1)` against `rhs_branch_b` and found a concrete bug in the scalar-tensor potential implementation for branch C.
- Fixed `U_of_phi` and `uprime_of_phi` in `kg_examples/simulate_bc_psfd_imex.py` so the flat limit satisfies `U(1)=0` and `U'(1)=0`.
- After the potential fix:
  - `Cphi1` now returns to the same scale as `Bfull`
  - the remaining runaway is localized to the true `phi` dynamics in `Cfull`, not the flat-limit reduction
- Current narrowed problem:
  - A solver is acceptable as benchmark
  - B matter update is substantially healthier
  - C still has a stiff `phi` dynamics problem that needs its own targeted debug
- Follow-up checks after the above fixes:
  - Re-running `debug_c_decomposition.py` now gives `Bflat ~ Bfull` and `Cphi1 ~ Bfull`; the old catastrophic `Bflat` spike is gone once the conservative matter update is used.
  - `Cfull` still blows up to large `tau,beta`, so the remaining problem is genuinely in the dynamic `phi -> (tau,beta)` coupling.
- New isolation tests show:
  - evolving `phi` alone on fixed matter only produces a small deviation (`|phi-1| ~ 1.9e-2`)
  - but evolving `tau,beta,phi` together on fixed matter still drives `tau,beta` to the clipping ceiling
- Current interpretation:
  - the remaining C-branch failure is no longer in the B-style matter transport
  - it is now localized to the reduced `C` geometry source terms / `phi`-coupling itself
- Additional parameter probes confirm the same direction:
  - increasing `M_P` barely changes the C blow-up
  - decreasing `ell` makes the C branch markedly worse
- This indicates the dominant remaining stiffness is not the matter source `mu`, but the `phi`-variable representation itself.
- Added note `research-notes/049-C分支剩余stiffness来自phi变量而非物质源.md`.
- Added `kg_examples/debug_b_single_beam_propagation.py` to test single-beam (`alpha=1`) propagation for `Bflat/Bfull` against the exact A benchmark.
- Single-beam result:
  - `Bflat` and `Bfull` are almost indistinguishable
  - halving `dt` changes the solution only at about `4.5e-3` relative level
  - both remain far from `A exact` (`rel_to_exact ~ 0.327`)
- Interpretation: in the single-beam case, the remaining mismatch now looks much more like a genuine dynamical difference of the B matter system than a still-broken time integrator.
- Added note `research-notes/050-B分支单束下更像动力学差异而非积分器异常.md`.
- Performed a post-fix representative C scan against A for `ell = 2e-2, 1e-2, 5e-3` in both the double-beam (`alpha=0.5`) and single-beam (`alpha=1`) cases.
- Trend after the C fixes:
  - `ell` made smaller does **not** move C toward A
  - in both alpha cases the smallest tested `ell` is actually farther from A than `ell = 2e-2`
- Saved this representative trend in `visualizations/c_after_fix_scan_summary.json`.
- Added a weak-gravity representative scan after the C fixes:
  - file: `visualizations/c_weak_gravity_representative.json`
  - for `alpha=0.5`, sending `M_P` from `300` to `1e6` barely changes `B_rel_to_A` and only marginally changes `C_rel_to_A`
  - making `ell` smaller at fixed huge `M_P` makes `C` farther from `A`, not closer
- Rendered the representative weak-gravity C comparison as images in `visualizations/c_weak_gravity_representative/` for both `alpha=0.5` and `alpha=1.0`, including `A`, `C`, `|C-A|`, and centerline plots.
- User raised two important objections:
  - the current `A_exact_alpha_1` benchmark still looks spotty / multi-lobed rather than like a single smooth beam
  - if the intended comparison is that switching off the gravitational part should make A/B/C coincide, then the present B/C flat-limit implementation is conceptually suspect
- Current interpretation updated:
  - the spotty single-beam A benchmark likely still contains spectral/truncation artifacts and is not yet the ideal visual benchmark
  - more importantly, if the user’s intended theory requires `gravity -> 0` to collapse A/B/C to the same action, then the present flat-limit B/C matter PDE cannot be the right flat-limit implementation
- Added note `research-notes/051-gravity关掉时BC不应变成自由HJ系统.md`
- New corrected interpretation:
  - the previous `Bflat/Bfull/Cfull` flat-limit solvers were not testing the intended theory
  - they implicitly dropped the geometrization constraints and kept only the bare HJ+continuity part
  - therefore their mismatch with A cannot be used as a physical conclusion
- Added `kg_examples/check_constrained_gravity_off_equivalence.py`
- Numerical action-level check:
  - reconstructing `g~` and `rho~` from the exact A solution via the rank-1 constraints gives machine-precision current identity agreement
  - this confirms that the correct gravity-off B/C limit is “A plus algebraic reconstruction”, not an independently evolved flat HJ system
- Added note `research-notes/052-作用量层面看gravity关掉时的正确BC极限.md`
- Added `kg_examples/simulate_bc_gravity_off_target.py`
- This script is not a substitute solver for B/C; it is the correct gravity-off regression target built from the trusted A dynamics plus per-step constraint reconstruction.
- First run for `alpha=0.5` completed:
  - `relative_l1_A_solver ~ 7.1e-3`
  - meaning the A backbone inside this target script is consistent with the current calibrated A solver to sub-1% level
- Added `kg_examples/simulate_abc_zero_gravity_regression.py`
- This is a strict zero-gravity regression harness, not the final B/C dynamical solver:
  - all three branches share the same A backbone in the gravity-off limit
  - the produced `rho` maps verify the required regression target visually
- Current alpha=0.5 output:
  - after correcting the benchmark parameters to match `exact_kg_2p1_benchmark`, `A_num_vs_exact_relL1 ~ 2.0e-2`
  - `B_vs_A_relL1 = 0`
  - `C_vs_A_relL1 = 0`
- Important fix:
  - the earlier `abc_zero_gravity_regression` images were misleading because the harness used a different benchmark parameter set than the approved `exact2p1_alpha_0.5_xz_linear.png`
  - this has now been corrected so the regression harness uses the same `[-10,10]`, `nx/nz`, and `nkx/nkz` configuration as the approved benchmark
- Added `kg_examples/simulate_bc_small_gravity_backbone.py`
- This new harness is a solver-debug continuation model:
  - trusted A backbone evolves `psi`
  - B/C geometry is driven on top of it with a tunable `lambda_grav`
  - `rho_B/rho_C` are reconstructed so `lambda_grav=0` gives exact A agreement
- Verified continuity:
  - `lambda_grav=0` gives `B_vs_A_relL1 = 0`, `C_vs_A_relL1 = 0`
  - `lambda_grav=1e-3` keeps both B and C within `~1e-8` of A
  - `lambda_grav=1e-2` still keeps both branches extremely close to A
- Output directories:
  - `visualizations/bc_small_gravity_backbone_zero`
  - `visualizations/bc_small_gravity_backbone_1e3`
  - `visualizations/bc_small_gravity_backbone_1e2`
### 2026-04-26 additional

- Cleaned up `kg_examples/simulate_bc_small_gravity_backbone.py` so the displayed `rho_A_inner` is exactly the approved `exact_kg_2p1_benchmark` crop, not a separate outer-box A image.
- Re-ran:
  - `visualizations/bc_small_gravity_backbone_zero`
  - `visualizations/bc_small_gravity_backbone_1e3`
- With the corrected displayed benchmark:
  - `lambda_grav=0` now gives machine-precision `B/C -> A` agreement on the same shown `rho`
  - `lambda_grav=1e-3` still stays extremely close to A (`~1e-9` level)
- Pushed the same backbone harness to larger `lambda_grav`:
  - `lambda=0.1`: `B_vs_A ~ 3.13e-7`, `C_vs_A ~ 9.38e-6`
  - `lambda=1.0`: `B_vs_A ~ 3.13e-6`, `C_vs_A ~ 1.24e-3`
  - `lambda=1.0` with half `dt`: `C_vs_A` drops to `~5.7e-5`
  - `lambda=10`: both `B` and `C` sit near `~3e-5`, with `C_phi_dev_max = 0`
- Interpretation:
  - `B` remains smooth and near-continuous across these runs
  - `C` is still numerically suspect: it is dt-sensitive at `lambda=1`, and at `lambda=10` the auxiliary field is being numerically pinned back to `phi=1`, making C effectively collapse toward B
- Added an even weaker-gravity representative run:
  - `visualizations/bc_small_gravity_backbone_1e4`
  - for `lambda_grav=1e-4`, both `B` and `C` remain visually indistinguishable from `A`
  - summary: `B_vs_A_relL1 ~ 3.13e-10`, `C_vs_A_relL1 ~ 3.18e-10`
- Added `kg_examples/debug_c_phi_pathology.py`
- New C-solver diagnosis:
  - at the flat limit, `src_phi` is positive everywhere
  - after a single forward step, the entire field moves to `phi > 1`
  - then clipping forces it back to `phi = 1`
  - therefore the current C solver is numerically pinning `phi` off, not honestly evolving it
- Added note `research-notes/053-phi变量在平直极限的病态性.md`
- Added note `research-notes/054-fR模型在平直极限下先天贴近EH分支.md`
- Symbolic expansion confirms:
  - `f(R)=R/sqrt(1+(\ell^2 R)^2) = R - (1/2)\ell^4 R^3 + O(R^5)`
  - `f''(0)=0`
- Interpretation:
  - one reason C remains close to B on the shared weak-gravity baseline is theoretical, not purely numerical: this chosen `f(R)` model differs from EH only at cubic order in curvature near flat space
- Applied a consistency fix to `rhs_branch_c` so the `phi` source also vanishes in the `phi=1` / flat-phi limit.
- Re-ran `debug_c_phi_pathology.py`:
  - `src_phi` is no longer positive everywhere
  - one-step trial no longer overshoots above `phi=1`
  - clipping is no longer actively driving the dynamics at flat limit
- Re-ran the shared-baseline continuation:
  - `lambda_grav=1e-3`: `C_vs_A ~ 1.09e-8`, `C_vs_B ~ 7.75e-9`
  - `lambda_grav=1.0`: `C_vs_A ~ 1.81e-3`, `C_vs_B ~ 1.81e-3`
- Interpretation:
  - after the flat-phi source fix, C is no longer numerically pinned to B
  - C now shows a genuine, solver-resolved separation from B at order-one gravity strength
- Extended convergence study on the corrected C solver:
  - at `lambda_grav=1`, dt-halving gives
    - `dt~0.025`: `C_vs_A ~ 1.81e-3`
    - `dt~0.0125`: `C_vs_A ~ 5.76e-5`
    - `dt~0.00625`: `C_vs_A ~ 1.34e-5`
  - so `lambda~1` now looks like a genuinely converging regime
- Strong-gravity continuation:
  - at `lambda_grav=10`, dt-halving gives
    - `dt~0.025`: `C_vs_A ~ 3.13e-1`
    - `dt~0.0125`: `C_vs_A ~ 9.24e-2`
    - `dt~0.00625`: `C_vs_A ~ 3.71e-4`
  - this is a dramatic improvement with shrinking dt, but still indicates the solver is under-resolved / unstable at strong gravity on coarse steps
- Saved a compact numerical summary to `visualizations/c_dt_convergence_scan/summary.json`
- Extended the strong-gravity dt scan one level further:
  - `lambda_grav=10`, `dt~0.003125`: `C_vs_A ~ 1.24e-4`, `C_vs_B ~ 9.30e-5`
- Interpretation:
  - even at `lambda~10`, the corrected C solver continues to converge downward as dt shrinks
  - the previous huge deviations at coarse dt were overwhelmingly numerical, not trustworthy physics
- Added C-geometry subcycling to `kg_examples/simulate_bc_small_gravity_backbone.py`
- Strong-gravity (`lambda=10`, outer `dt~0.025`) results with C substeps:
  - `substeps=1`: `C_vs_A ~ 3.13e-1`
  - `substeps=2`: `C_vs_A ~ 1.03e-1`
  - `substeps=4`: `C_vs_A ~ 3.15e-4`
  - `substeps=8`: `C_vs_A ~ 1.26e-4`
- This matches the dedicated fixed-small-step run (`dt~0.003125`) very closely.
- Saved compact summary to `visualizations/c_subcycling_scan/summary.json`
- Interpretation:
  - the strong-gravity C solver can now be stabilized efficiently by subcycling the C geometry subsystem, without shrinking the global backbone step for A/B
- Compiled the current ABC comparison baseline into:
  - `visualizations/abc_comparison_backbone_scan/summary.json`
  - `visualizations/abc_comparison_backbone_scan/BA_CA_vs_lambda.png`
  - `visualizations/abc_comparison_backbone_scan/CB_vs_lambda.png`
  - This summary uses:
  - direct continuation values in the weak regime
  - subcycled C values in the strong regime
  - and keeps the unified zero-gravity baseline fixed

### 2026-04-27

- Created a formal handoff document for moving the project to another machine / another agent:
  - `E:\量子势几何化_codex\agent-memory\HANDOFF-2026-04-27.md`
- The handoff summarizes:
  - current theory status
  - solver status
  - trusted outputs
  - current ABC comparison baseline
  - recommended next step for the next agent
- Performed dt-halving tests on the corrected shared-baseline continuation solver:
  - `lambda_grav=1.0`, `dt~0.025`: `C_vs_A ~ 1.81e-3`
  - `lambda_grav=1.0`, `dt~0.0125`: `C_vs_A ~ 5.76e-5`
  - `lambda_grav=1.0`, `dt~0.00625`: `C_vs_A ~ 1.34e-5`
- Interpretation:
  - at `lambda_grav ~ 1`, the corrected C solver is strongly dt-sensitive but converging downward as dt shrinks
  - much of the earlier apparent C/A separation at `lambda=1` was still numerical
- Stronger-gravity probe:
  - `lambda_grav=10`, `dt~0.0125` gives `C_vs_A ~ 9.24e-2`, `C_tau_max = 2.0`, `C_phi_dev_max = 0.8`
  - so C is not solved yet in the strong-gravity regime; the solver still breaks there

## 轮回 2：2026-04-27 跨机恢复

- 判定为 restart：当前 agent 无法从 live context 可靠知道自己继承了旧项目的第几个轮回，因此按重启式接管处理
- 读取并确认保留旧记忆文件：`README.md`、`TASKS.md`、`LOG.md`、`DECISIONS.md`
- 额外读取交接文件：`agent-memory/HANDOFF-2026-04-27.md`
- 确认本项目是从另一台电脑复制到当前机器；历史 Windows 路径 `E:\量子势几何化_codex` 与当前 macOS 路径 `/Users/wangyunfei/Desktop/量子势几何化_codex` 指向同一份项目副本
- 从旧 `LOG.md` 恢复到的最新显式长记忆状态为 `ON`
- 用户本轮再次明确要求“把它加进你的 skill 里，然后开启”；长记忆状态在当前机器上续确认为 `ON`
- 本轮接管策略：不覆盖任何旧记忆文件，只追加新的轮回、状态确认、路径映射与协议固化信息
- 已在当前机器创建自定义 skill 定义文件：`/Users/wangyunfei/.codex/skills/agent-memory/SKILL.md`
- 本轮完成的长期记忆维护：
  - 在 `README.md` 中补充长记忆协议、跨机器路径映射和本机 skill 入口
  - 在 `TASKS.md` 中补充跨机器接管进展与后续每轮 reload 要求
  - 在 `DECISIONS.md` 中补充“跨机器接管只追加新轮回”的抽象决策
- 按长记忆协议重新读取了全部核心 memory 文件、交接手册、所有研究笔记标题结构、全部 Python 源码函数结构，以及整个项目文件清单
- 对当前活跃理论/数值主线补读了：
  - `research-notes/052-作用量层面看gravity关掉时的正确BC极限.md`
  - `research-notes/053-phi变量在平直极限的病态性.md`
  - `research-notes/054-fR模型在平直极限下先天贴近EH分支.md`
  - `kg_examples/simulate_bc_small_gravity_backbone.py`
  - `kg_examples/simulate_bc_psfd_imex.py`
  - `kg_examples/simulate_abc_zero_gravity_regression.py`
  - `kg_examples/simulate_bc_gravity_off_target.py`
  - `kg_examples/exact_kg_2p1_benchmark.py`
  - `kg_examples/calibrate_exact_kg_ps_fd_absorbing.py`
  - `kg_examples/check_constrained_gravity_off_equivalence.py`
- 当前机器缺少 `matplotlib`；已在项目根目录新增本地运行依赖目录 `.agent-deps/`，并用它补跑后续汇总图与 continuation 扫描
- 已更新 `kg_examples/scan_abc_comparison_backbone.py`：
  - 把稀疏扫描改成更密的 `lambda_grav` 网格
  - 引入分段 `c_substeps`
  - 追加阈值摘要输出
  - 增加“已有 case 直接复用”的可恢复扫描机制
- 已完成共享 baseline 下更密的 continuation 扫描，并新增阶段文档：
  - `research-notes/055-共享baseline下A-B-C偏离的更密lambda扫描.md`
- 当前 dense scan 的关键结果：
  - `B_vs_A_relL1 >= 1e-6` 首次出现在 `lambda ~ 1`
  - `B_vs_A_relL1 >= 1e-5` 首次出现在 `lambda ~ 10`
  - `C_vs_B_relL1` 在 `lambda ~ 3` 稳健进入 `1e-5` 量级
  - 在当前 `lambda <= 10` 扫描内，`C_vs_B_relL1` 仍未超过 `1e-4`
  - `C_vs_A_relL1 >= 1e-4` 首次出现在 `lambda ~ 10`
- 额外运行代价记录：
  - 单独补跑 `lambda=3, c_substeps=8` 用时约 `257 s`
  - 这进一步说明 continuation 脚本必须支持 resume，而不该反复全扫
- 已新增同初值动态图样比较脚本：
  - `kg_examples/compare_abc_same_initial_evolution.py`
- 这一步采用的口径是：
  - A 用可信的谱 Laplacian + 吸收层 + 精确边界骨架做动态演化
  - B/C 在同一组 A 的 `rho,S` 历史上推进几何子系统
  - 再由几何场重构 `rho_B, rho_C`
- 第一版脚本最初把 A 对精确 benchmark 的检查放在每一个时间步上，导致 `lambda=1` 的首轮作业耗时明显过高
- 已中止那条低效旧作业，并把脚本改成：
  - 只在最终时刻做一次 `A_num_vs_exact_relL1` spot check
  - 保留完整的 `B/A/C` 动态误差曲线与 montage 出图
- 提速后的脚本已通过 `py_compile`
- 已完成两组同初值动态比较：
  - `visualizations/abc_same_initial_evolution_lambda1/`
  - `visualizations/abc_same_initial_evolution_lambda3/`
- `lambda=1, c_substeps=4` 的最终结果：
  - `B_vs_A_relL1 = 4.22e-6`
  - `C_vs_A_relL1 = 1.89e-5`
  - `C_vs_B_relL1 = 1.47e-5`
  - `A_num_vs_exact_relL1 = 1.77e-2`
  - 首次过阈时间：
  - `B_vs_A >= 1e-6` 于 `t/T ~ 0.341`
  - `C_vs_A >= 1e-5` 于 `t/T ~ 0.505`
  - `C_vs_B >= 1e-5` 于 `t/T ~ 0.576`
- `lambda=3, c_substeps=8` 的最终结果：
  - `B_vs_A_relL1 = 1.26e-5`
  - `C_vs_A_relL1 = 4.55e-5`
  - `C_vs_B_relL1 = 3.29e-5`
  - `A_num_vs_exact_relL1 = 1.77e-2`
  - 首次过阈时间：
  - `B_vs_A >= 1e-6` 于 `t/T ~ 0.174`
  - `C_vs_A >= 1e-5` 于 `t/T ~ 0.299`
  - `C_vs_B >= 1e-5` 于 `t/T ~ 0.365`
- 当前动态图样层面的判断：
  - 三支在 `t=0` 时几乎重合
  - 主要偏离是在双束重叠过程中逐步积累出来的
  - B 相对 A 的偏离更弱、更晚、更平滑
  - C 相对 B 的偏离更早出现，也更集中在干涉主结构附近
- 已新增阶段文档：
  - `research-notes/056-同初值下2+1d高斯波包干涉的ABC后续演化比较.md`
- 已新增只跑 A/B 的实验量提取脚本：
  - `kg_examples/analyze_ab_measurable_observables.py`
- 已生成输出：
  - `visualizations/ab_measurable_observables_lambda1/`
- 当前默认测量示例：
  - 物理条纹间距 `d_phys = 1 micron`
  - 读出尺度 `L_readout = 1 m`
  - 纳赫兹 GW 代表应变 `h = 2.4e-15`
- `lambda=1` 下，中心条纹的实验量提取结果：
  - `fringe_spacing_rel_diff ~ 0`
  - `effective_phase_shift ~ -5.58e-12 rad`
  - `equivalent_fringe_shift / spacing ~ -8.88e-13`
  - `visibility_rel_diff ~ 3.54e-8`
  - `central_peak_rel_intensity_shift ~ 9.78e-5`
  - `mean_peak_rel_intensity_shift ~ 9.29e-5`
  - `support_pointwise_rel_diff_max ~ 1.05e-4`
- 当前最关键的新判断：
  - `lambda=1` 的 A/B 差异几乎不表现为条纹间距或相位滑移
  - 最现实的可测入口是亮纹强度约 `1e-4` 的系统性重标定
- 对 Earth Schwarzschild 的当前弱场估计：
  - 使用 `kappa_E = GM/(c^2 R^3) ~ 1.72e-23 m^-2`
  - `L_readout = 1 m` 时，代表取向的读出调制量级约为 `1e-23`
  - 最大方向为径向，最小方向为 magic angle `~54.7 deg`
- 对已观测纳赫兹 GW 背景的当前弱场估计：
  - 采用 NANOGrav 15-year 代表应变 `h ~ 2.4e-15`
  - 最佳取向下的读出调制量级约为 `h/2 ~ 1.2e-15`
  - 最大方向为横向且沿 `+` 偏振主轴，沿传播方向响应为 `0`
- 已新增阶段文档：
  - `research-notes/057-lambda1下AB可测量实验量与弱曲背景估计.md`
- 用户进一步追问两点：
  - `lambda=1` 下当前世界的弱实验差异是否支持两种作用量尚未被现有实验排除
  - 为什么当前 `C` 没有比 `B` 更贴近 `A`
- 本轮澄清：
  - 当前 `lambda=1` 的弱差异确实支持一个温和判断：在当前受控 scaffold 与当前这类 observables 上，现有实验量级并没有直接把 A/B 或 A/C 这种作用量选择排掉
  - 但这还不是“理论已经完全安全”，因为还没扫完所有背景、observable 和 full curved-background / full backreacted 版本
- 更关键的模型澄清：
  - 当前代码里真正实现并跑起来的 `C`，是 `f(R~)=R~/sqrt(1+(\ell^2 R~)^2)` 分支
  - 它不是项目早期那条直接按 `Q` 做屏蔽的 `F(Q)R~` 思路
  - 因此不能把“早期希望通过屏蔽减小实验偏差”的直觉，直接投射成“当前这个 `f(R~)` C 必须比 B 更接近 A”
- 当前最准确的解释口径：
  - `C` 的设计目标首先是不要在现实弱背景里带来巨大的新偏离
  - 这不等价于要求 `C` 在每个参数点都比 `B` 更贴近 `A`
  - 当前数值里 `B` 比 `C` 更接近 `A`，只说明当前选定的 `f(R~)` 修正没有专门去抵消 `B-A`，而不是说明 `C` 的物理动机本身失败
- 用户现已明确收束主线：
  - 暂时挂起 `C` 与更一般作用量
  - 当前项目只讨论 `A/B`
  - 当前阶段的目标不是“把理论完全证明安全”，而是找到至少一个可测、非数值误差的 A/B observable，并把它与噪声基准比较
- 已据此把任务与抽象决策更新为：
  - 先做 `lambda=1` 下 A/B observable 自身的收敛性检查
  - 再围绕用户指定的 Earth surface vs L1、GW vs flat 场景，以及额外的 orientation/null-test 场景，组织实验差异与噪声比较
- 已新增三套脚本，继续补完此前明确缺的三层：
  - `kg_examples/scan_ab_observable_convergence.py`
  - `kg_examples/analyze_ab_time_observables.py`
  - `kg_examples/compare_ab_background_scenarios.py`
- 本机缺 `matplotlib`，已按用户许可通过国内镜像补装到当前 Python 用户环境：
  - `python3 -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple matplotlib`
- 已完成 `lambda=1` 下 A/B observable 的直接 dt/grid 收敛扫描：
  - 输出：`visualizations/ab_observable_convergence/`
  - refined 主 observable：
  - `central_peak_rel_intensity_shift = 6.761437198341409e-05`
  - `mean_peak_rel_intensity_shift = 6.561680454841884e-05`
  - `support_pointwise_rel_diff_max = 7.820840064539984e-05`
  - 数值地板：
  - `central peak floor ~ 3.57e-09`
  - `mean peak floor ~ 6.70e-08`
  - `support floor ~ 2.21e-08`
  - 结论：
  - 中心亮纹强度差高出数值地板约 `1.9e4`
  - 平均亮纹强度差高出数值地板约 `9.8e2`
  - support max 高出数值地板约 `3.5e3`
  - `visibility` 虽非零，但只高出空间离散地板约 `O(1e0)`
  - `phase shift` 与 `equivalent fringe shift` 的信号远低于空间离散地板，当前不能作为主 observable
- 已完成 `lambda=1` 下 A/B 时间序列与积分型 readout 分析：
  - 输出：`visualizations/ab_time_observables_lambda1/`
  - 首次过阈时间：
  - `B_vs_A_relL1 >= 1e-6` 于 `t/T ~ 0.3725`
  - `center_peak_rel_shift >= 1e-5` 于 `t/T ~ 0.4066`
  - `central_disk_rel_shift >= 1e-5` 于 `t/T ~ 0.4165`
  - `bright_region_rel_shift >= 1e-5` 于 `t/T ~ 0.1297`
  - 最终积分型 readout：
  - `center_point_rel_shift = 6.761080361107943e-05`
  - `central_disk_rel_shift = 6.683999408006137e-05`
  - `bright_region_rel_shift = 6.625028280014768e-05`
  - `interference_box_rel_shift = 6.557413714432494e-05`
  - 当前判断：
  - A/B 分裂不是初始时刻就有，而是在双束重叠过程中建立
  - 积分型亮区 readout 比单个中心峰更早开始响应
  - 多种 readout 一致落在 `6.6e-05 ~ 6.8e-05`，说明亮纹强度重标定不是孤立像素伪影
- 已把弱曲背景层升级成场景梯度与取向扫描：
  - 输出：`visualizations/ab_background_scenarios_lambda1/`
  - 当前脚本优先读取 refinement 后的平直基线，而不再沿用旧的粗口径 `~9.8e-05`
  - 纳入场景：
  - Earth surface
  - LEO 400 km
  - GEO
  - Sun-Earth L1
  - `1 tonne @ 10 cm, readout = 1 cm` 的近源质量调制场景
  - 纳赫兹 GW 背景
  - 当前所有这些都仍是 beam region 内的 vacuum tidal / readout modulation
  - 代表结果：
  - Earth surface 轴向调制 `~3.43e-23`
  - Sun-Earth L1 轴向调制 `~3.54e-30`
  - 近源质量场景轴向调制 `~1.49e-25`
  - GW 最佳方向调制 `~1.2e-15`
  - 把这些调制乘到当前 refined 主 observable 上后，背景诱导改变量仍远低于当前数值地板
  - 例如 Earth surface 轴向仅相当于 `~6.5e-19` 个数值地板，GW 也只有 `~2.3e-11`
- 已新增阶段文档：
  - `research-notes/058-lambda1下AB三层补完：收敛、时间序列与背景场景.md`
- 已把 `kg_examples/analyze_ab_measurable_observables.py` 的默认口径对齐到 refined 基线：
  - `dt_target=0.0125`
  - `nx=nz=181`
  - `nkx=nkz=81`
  - 并已重写 `visualizations/ab_measurable_observables_lambda1/summary.json`
- 用户接受了把原始强度型 observable 进一步改写成“两个实验设置之间的差分量”的方向，并指定优先计算：
  - `phase-flip differential`
  - `interference excess`
- 已新增脚本：
  - `kg_examples/analyze_ab_differential_experiment_observables.py`
- 已在三档时间步下完成最小核对：
  - `visualizations/ab_differential_experiment_observables_lambda1_dt025/`
  - `visualizations/ab_differential_experiment_observables_lambda1/`
  - `visualizations/ab_differential_experiment_observables_lambda1_dt00625/`
- 定义口径：
  - `D_phi = (N(phi=0)-N(phi=pi)) / (N(phi=0)+N(phi=pi))`
  - `E_int = (N_both-(N_beam1+N_beam2)) / (N_beam1+N_beam2)`
  - 真正看的理论分裂是 `Delta_AB = observable_B - observable_A`
- 两个简单读出窗都已比较：
  - `central disk`：半径 `1.5`
  - `interference box`：`|x|<=2, |z|<=2`
- refined (`dt=0.0125`) 下结果：
  - `phase-flip differential, central disk`：
  - `Delta_AB final = -2.304989379919853e-09`
  - `max |Delta_AB| = 2.934160958992238e-08` at `t/T ~ 0.876`
  - `phase-flip differential, interference box`：
  - `Delta_AB final = -5.3196539377513286e-08`
  - `max |Delta_AB| = 5.3196539377513286e-08` at `t/T = 1`
  - `interference excess, central disk`：
  - `Delta_AB final = -2.39155306669403e-10`
  - `max |Delta_AB| = 3.3662228959840945e-07` at `t/T ~ 0.541`
  - `interference excess, interference box`：
  - `Delta_AB final = -1.1791510912129155e-07`
  - `max |Delta_AB| = 1.2542084717459545e-06` at `t/T ~ 0.577`
- 三档时间步核对后，以上四类数值都稳定：
  - `D_phi` 在 `dt=0.025,0.0125,0.00625` 下的最终值和峰值基本重合
  - `E_int` 尤其是 `interference box` 下的峰值 `~1.25e-06` 也基本重合
- 当前新增判断：
  - 这两类双设置差分量确实给出稳定非零的 A/B 分裂
  - 但它们会把此前 `~6.7e-05` 的原始强度重标定大幅共模抵消掉
  - 因而：
  - `raw intensity-type channel` 提供最大的理论分裂量级
  - `differential channel` 提供更实验化、更抗共模噪声的读出方式
  - 在差分通道内部，`interference excess` 明显优于 `phase-flip differential`
- 已新增阶段文档：
  - `research-notes/059-lambda1下AB的两类实验差分量：相位翻转与干涉超额.md`
- 用户继续追问：若在地球表面把装置轴向旋转，实验差异本身能有多大
- 已进一步澄清当前 Earth-surface 取向实验的口径：
  - 在当前外部 Schwarzschild 真空 + 弱潮汐 scaffold 下，只对“相对地球径向的倾角”敏感
  - 若保持同一倾角、仅绕竖直方向改方位角，则当前模型预测为 `0`
- 已把“旋转产生的 observable 差值”压成数字：
  - 对 refined 原始主 observable `central_peak_rel_intensity_shift ~ 6.7614e-05`
  - `vertical(axial) -> horizontal(transverse)` 的差值约：
  - `3.48e-27`
  - `vertical(axial) -> magic angle` 的差值约：
  - `2.32e-27`
  - 对最佳差分实验通道 `interference excess @ interference box`
  - 峰值 `~1.256e-06` 时，`vertical -> horizontal` 的差值也只有：
  - `6.46e-29`
  - 若看最终时刻值，则只有：
  - `6.07e-30`
- 当前判断：
  - 地表旋转取向实验在现有真空弱场口径下几乎为零
  - 它更适合作为理论上的 null/consistency check，而不是现实中最值得优先追的信号源
- 用户同意进入“噪声与系统误差基准”这一步
- 已新增脚本：
  - `kg_examples/analyze_ab_noise_baselines.py`
- 已新增输出：
  - `visualizations/ab_noise_baselines_lambda1/`
- 已把四条实验通道放入统一基准：
  - `raw center peak`
  - `raw interference box`
  - `phase-flip differential @ interference box`
  - `interference excess @ interference box`
- 当前 shot-noise 计数需求：
  - `raw center peak`：
  - `5 sigma ~ 5.47e9 counts`
  - `raw interference box`：
  - `5 sigma ~ 5.81e9 counts`
  - `phase flip`：
  - `5 sigma ~ 8.83e15 total counts across two phase states`
  - `interference excess`：
  - `5 sigma ~ 1.19e13 denominator counts`
  - `5 sigma ~ 1.79e13 total counts across all three runs`
- 当前系统误差容忍度：
  - raw 通道若只给单项系统误差 `10%` 预算，则残余 gain/source drift 需压到 `~6e-06`
  - `phase flip` 的 pair imbalance 若只给 `10%` 预算，则需压到 `~5e-09`
  - `interference excess` 的 single-beam normalization 若只给 `10%` 预算，则需压到 `~2.5e-07`
- 当前新增判断：
  - 在最简单的 shot-noise + 漂移 / 归一化模型下，`raw intensity` 反而比 differential channel 更适合作为第一发现通道
  - `interference excess` 仍值得保留，但更适合作为 companion cross-check
  - `phase flip` 当前统计和系统误差代价都太高，不适合主打
- 已新增阶段文档：
  - `research-notes/060-lambda1下AB实验通道的噪声基准比较.md`
- 本轮结束时的最稳结论：
  - `lambda=1` 下 A/B 最扎实的非零差异是亮纹/热点强度约 `6.7e-05` 的系统性重标定
  - 这件事已经通过 dt/grid refinement、时间序列和多种积分型 readout 共同支撑
  - Earth / L1 / GW / 近源质量这些当前真实弱背景，在真空弱场口径下都不会把它进一步放大成新的主实验入口
  - 下一轮若继续推进，应转到“实验噪声与系统误差基准”而不是继续纠缠相位/条纹平移量

## 轮回 2 · 2026-04-27 17:38:19 CST · 演示稿交付

- 用户请求：
  - 先不要进入 Overleaf 或论文正文，而是先做一份能够“把所有细节讲明白”的 PPT，并且必须包含数值模拟配图
- 已完成演示稿工作流：
  - 使用 `Presentations` skill
  - 新建 deck workspace：
    - `presentation-workspaces/quantum-potential-geometrization-ab/`
  - 新建源码：
    - `presentation-workspaces/quantum-potential-geometrization-ab/src/deck.mjs`
  - 成功导出原生 `pptx`：
    - `presentation-workspaces/quantum-potential-geometrization-ab/output/output.pptx`
  - 成功渲染逐页 PNG 预览：
    - `presentation-workspaces/quantum-potential-geometrization-ab/scratch/previews/slide-01.png` ... `slide-14.png`
  - 已生成 contact sheet：
    - `presentation-workspaces/quantum-potential-geometrization-ab/scratch/previews/contact-sheet.png`
- 当前 deck 内容结构：
  - 共 `14` 页
  - 主线为：
    - 理论问题与 action-level ambiguity
    - 几何化骨架与 A/B/C 分支定位
    - shared-backbone 数值 scaffold
    - ABC 动力学分叉
    - `lambda=1` 下 A/B 原始 observable
    - 收敛与数值地板
    - 时间结构
    - 背景时空场景
    - 差分实验量
    - 噪声基准
    - 当前可投稿主张
    - C/continuation appendix
    - 文件地图 appendix
- 数值模拟配图已实际嵌入：
  - `abc_same_initial_evolution_lambda1/rho_montage.png`
  - `abc_same_initial_evolution_lambda1/relL1_vs_time.png`
  - `ab_measurable_observables_lambda1/centerline_AB_rel_diff.png`
  - `ab_observable_convergence/dt_scan_observables.png`
  - `ab_time_observables_lambda1/signed_intensity_shifts_vs_time.png`
  - `ab_background_scenarios_lambda1/scenario_max_fractional_modulation.png`
  - `ab_differential_experiment_observables_lambda1/interference_excess_interference_box.png`
  - `ab_noise_baselines_lambda1/required_counts_by_channel.png`
  - `abc_same_initial_evolution_lambda3/diff_montage.png`
- 质量检查结果：
  - 已运行 `check_presentation_quality.js`
  - 最终 `quality-report.json` 为无失败、无警告
  - `slide_count = 14`
  - `media_count = 10`
  - `zero_byte_media = []`
  - `placeholder_text = []`
- 本轮踩坑与修复：
  - 初版 `pptx` 曾出现 `zero-byte media`
  - 原因不是图不存在，而是 `artifact-tool` 的图片资产若只给 `path/uri`，导出时可能不真正内嵌数据
  - 已修复方案：
    - 先把图复制到当前 deck workspace 的 `scratch/assets/`
    - 再转成 `data URL` 传给 `image(...)`
    - 这样 `pptx` 与预览都能稳定拿到真实图片 bytes
  - 另一个修复：
    - `drawSlideToCtx` 预览渲染应喂 `slide model + presentation model`
    - 不应直接喂 `toProto()` 的扁平结果
- 当前状态：
  - PPT 已可直接打开、继续编辑、删减成组会版或答辩版
  - 用户下一步大概率会在“论文正文 / Overleaf 同步 / 汇报精简版”三者里继续推进

## 轮回 2 · 2026-04-27 18:52:34 CST · PPT 口径纠偏

- 用户检查当前 PPT 后，明确质疑其中“shared backbone / 数值 scaffold / full-backreact”等术语的物理含义，并指出主结果与其原始要求“全动力学演化且包含物质对背景时空的反作用”不一致。
- 已重新核对当前主比较脚本与研究笔记：
  - 当前 PPT 主结果实际来自 `kg_examples/compare_abc_same_initial_evolution.py`
  - 其方法口径是“共享主干比较”：A 分支先可信演化，再把 A 的 `rho,S` 历史喂给 B/C 的几何子系统，最后重构 B/C 的 `rho`
  - 这不是“全反作用自洽演化”
- 已确认：项目中确有更接近全动力学目标的原型与方程留档（如 `research-notes/044-2+1d下三分支动力学方程与数值闭合.md` 与 `kg_examples/simulate_2p1_full_dynamics.py`），但它们没有被本轮 PPT 作为主结果使用。
- 已形成新的展示约束：
  - 凡是不是“全反作用自洽演化”的结果，都必须在标题、图注、正文和口头说明里显式标出其方法口径
  - 图、表、PPT 中所有展示量必须先定义对象、窗口、时间切片与归一化，再给出数值或曲线
  - 正文统一采用中文，必要术语写成“中文（English）”

## 轮回 2 · 2026-04-27 19:09:37 CST · 回扣当前物质场变换

- 用户要求回到理论细节，先重新厘清“当前对物质场的变换具体是什么”。
- 已重新对齐当前主线笔记：
  - `research-notes/019-去掉sqrt-g-rho不变假设后的纯u-ansatz翻盘.md`
  - `research-notes/020-最简单解的可逆性签名因果结构与经典极限.md`
  - `research-notes/023-阶段性总结与Einstein作用量耦合框架.md`
  - `research-notes/040-三种完整作用量的统一约束形式.md`
- 当前主线结论：
  - 物质场不是按“标准复标量场线性重定义”去变，而是在 Madelung 变量下保持相位 `S` 不变，只让密度/幅度 `ρ -> ρ~` 非线性重标定
  - 一般变换由
    - `g~^(μν) = C g^(μν) + D u^μ u^ν`
    - `sqrt(-g~) ρ~ (C + D X) = sqrt(-g) ρ`
    共同决定
  - 最简单主解下
    - `C = 1`
    - `D = Q / X^2`
    - `sqrt(-g~) ρ~ = sqrt(-g) ρ X / m^2`
    - `ρ~ = ρ sqrt(X/m^2)`（在当前主时间样分支 `X>0` 上）

## 轮回 2 · 2026-04-27 19:09:37 CST · 核对 g 视角下的 Madelung 符号

- 用户要求回到最原始的 g 视角，逐步检查 Madelung 分解怎样给出 HJ 方程与连续性方程，并特别检查 `Q` 的符号与 `η` 的号规。
- 已按项目第一份正式笔记重新对齐：
  - `η_{μν} = diag(1,-1,-1,-1)`
  - `\Box = ∂_μ ∂^μ = ∂_t^2 - \nabla^2`
  - `Q = - (\Box \sqrt{ρ}) / \sqrt{ρ}`
  - 因此 HJ 方程写成 `g^{μν} ∂_μ S ∂_ν S = m^2 - Q`

## 轮回 2 · 2026-04-27 19:30:05 CST · 对比 Bohm 原文与相对论文献的 Q 记号

- 用户要求把当前 `Q` 的定义和其他文献、尤其是 Bohm 原始文献做仔细对比，并注意不同号规。
- 已核对文献：
  - Bohm 1952 原文 `A Suggested Interpretation of the Quantum Theory in Terms of "Hidden" Variables. I`
  - Chavoya-Aceves 2003 `A de Broglie-Bohm Like Model for Klein-Gordon Equation`
  - Jalalzadeh & Capistrano 2020 `Bohmian mechanics of Klein-Gordon equation via quantum metric and mass`
- 当前关键澄清：
  - Bohm 原文在非相对论 Schrödinger 情形里用的是能量维度的量子势 `U = -(ħ^2/2m) ∇^2R / R`，并不直接涉及 `η` 号规
  - 我们项目里的 `Q = - \Box R / R` 是相对论 KG 口径下的“质量平方修正”，量纲与 Bohm 原文的 `U` 不同
  - 在我们的 `η = diag(1,-1,-1,-1)` 约定下，非相对论极限满足 `Q_KG ≈ -2m U_Bohm`
  - Jalalzadeh 2020 采用的是相反号规 `(-,+,+,+)`，并把 `Q` 定义成无量纲量 `Q = -(1/m^2)\Box R / R`；把号规和量纲统一后，与当前项目的 HJ 结构等价

## 轮回 2 · 2026-04-27 19:44:22 CST · 切换 Q 的默认记号

- 用户明确要求：虽然旧定义内部自洽，但为了让当前项目的 `Q` 在非相对论近似下与 Bohm 原始量子势保持同向，应去掉原定义前面的负号。
- 已执行一次“约定迁移”而不是零散替换：
  - 新默认约定改为 `Q = + \Box sqrt(ρ) / sqrt(ρ)`
  - 因而当前原始 HJ 方程改写为 `X = m^2 + Q`
  - 分部积分后的作用量口径改写为 `L ~ X - m^2 - Q`
  - 最简单 rank-1 解改写为 `g~^(μν) = g^(μν) - (Q/X^2) u^μ u^ν`
  - 对应最简单系数改为 `D = -Q/X^2`
- 已更新的核心理论笔记：
  - `research-notes/001-KG作用量的Madelung分解与量子势定位.md`
  - `research-notes/014-按逆度规ansatz直接消去拉氏量中的Q项.md`
  - `research-notes/019-去掉sqrt-g-rho不变假设后的纯u-ansatz翻盘.md`
  - `research-notes/020-最简单解的可逆性签名因果结构与经典极限.md`
  - `research-notes/023-阶段性总结与Einstein作用量耦合框架.md`
- 已新增统一迁移说明：
  - `research-notes/061-Q约定改为正BoxR除以R.md`
- 已更新的活跃脚本：
  - `kg_examples/simulate_bc_gravity_off_target.py`
  - `kg_examples/check_constrained_gravity_off_equivalence.py`
  - `kg_examples/kg_double_packet.py`
  - `kg_examples/kg_double_packet_3d.py`
  - `kg_examples/kg_double_slit_beamlike_bohm.py`
  - `kg_examples/kg_double_packet_curvature.py`
- 已做基础验证：
  - `python3 -m py_compile` 已通过，覆盖上述全部改动脚本
  - `python3 kg_examples/check_constrained_gravity_off_equivalence.py` 可正常运行并输出结果，说明这次修改至少没有破坏当前活跃回归脚本的可执行性
- 当前约束：
  - 历史研究笔记仍有大量旧记号存在，不应当作“当前默认约定”直接使用
  - 阅读旧笔记时，一律先按 `Q_old = -Q` 做心里换元

## 轮回 2 · 2026-04-27 20:12:19 CST · 重新审视变换后 HJ 常数与质量为零极限

- 用户指出我上一次分析中的一个关键问题：不应当过早把变换后 HJ 常数 `\kappa` 解到密度变换关系的右边分母上；原始未化简关系应保留成 `\kappa \sqrt{-\tilde g}\tilde\rho = X \sqrt{-g}\rho`。
- 重新推导后确认：
  - 变换后一般 HJ 方程取 `\tilde g^{μν}u_μu_ν = \kappa` 时，
    - HJ 给出 `X(C+DX)=\kappa`
    - 连续性精确匹配给出 `\sqrt{-\tilde g}\tilde\rho (C+DX)=\sqrt{-g}\rho`
  - 两式联立的最稳形式确实是
    - `\kappa \sqrt{-\tilde g}\tilde\rho = X \sqrt{-g}\rho`
  - 只有在 `\kappa ≠ 0` 时，才进一步化成
    - `\sqrt{-\tilde g}\tilde\rho = X \sqrt{-g}\rho / \kappa`
- 因而需要修正之前的口径：
  - “`m^2` 可以换成任意固定非零常数”这点保留
  - 但对 `\kappa=0` 的讨论必须单独视为奇异极限，不能和非零常数并列处理
- 当前新的理解：
  - `\tilde g^{μν}u_μu_ν = 0` 是否可行，不是简单的“定义错了”，而是当前 rank-1 ansatz、同一相位 `S`、且要求有限 `\tilde\rho` 时会遇到专门的 massless/null-shell 难点
  - 若真要认真处理 `m=0` 或 `\kappa=0`，需要把“如何恢复非平凡连续性方程”作为单独子问题重新推导，而不是套用 `\kappa ≠ 0` 的那套解法

## 轮回 2 · 2026-04-27 20:57:35 CST · 纠正连续性恢复条件与 C 的自由度

- 用户进一步指出两点：
  - 我之前把“恢复连续性方程”错误地强化成了括号内流密度必须逐点相等
  - 还需要单独检查“强制 `C=1` 是否才是问题所在”
- 重新整理后确认：
  - 变换后连续性方程写成
    - `∂_μ [ sqrt(-g~) ρ~ (C+DX) u^μ ] = 0`
  - 用 HJ 关系 `X(C+DX)=κ` 后，只能推出
    - `∂_μ [ sqrt(-g~) ρ~ (κ/X) u^μ ] = 0`
  - 若取 `κ = m^2`，则这一步已经足够得到
    - `∂_μ [ sqrt(-g~) ρ~ (m^2/X) u^μ ] = 0`
  - 这与原连续性方程逐点同流只是一个更强的局域选择，而不是逻辑必要条件
- 对 `C` 的结论更新为：
  - 对任意非零 `m^2`，并要求 `Q→0` 时 `g~→g`、`ρ~→ρ`：
    - `C` 并不唯一
    - 只需 `C→1`
    - 然后取 `D = (m^2 - C X)/X^2`
    - 就能恢复变换后 HJ 方程
  - 在这种非零质量分支里，连续性方程中的系数 `m^2/X` 其实与 `C` 无关
- 对 `m=0` 的新结论：
  - 若还要求 `κ = m^2 = 0` 且 `Q→0` 时 `g~→g`、`ρ~→ρ`
  - 则当前 rank-1 ansatz 下不存在满足这些要求的正则 `C`
  - 因此真正的问题不在“强制 `C=1`”，而在“零质量 + null-shell + 当前物质作用量结构”的整体退化

## 轮回 2 · 2026-04-28 14:36:00 CST · 混合 ansatz 下 `ρ~ -> 0` 的严格判据

- 针对用户“混合型里也许不该先假定 `A ~ 1/X`，因为可能出现 `1/X^2` 或更高阶”的质疑，重新做了更严格的局域分析
- 在一般混合 ansatz
  `g~^(μν)=C g^(μν)+D u^μu^ν+E(u^μr^ν+r^μu^ν)+F r^μr^ν`
  下，先定义
  `A_u = C + DX + EY`
  `A_r = EX + FY`
  使得
  `g~^(μν)u_ν = A_u u^μ + A_r r^μ`
- 在 punctured neighborhood `U \\ Σ_X` 上，若同时要求：
  1. `u^μ` 与 `r^μ` 线性无关
  2. 原电流不为零
  3. 点态恢复连续性方程
  4. 变换后 HJ 壳长固定为非零常数 `κ`
  则会严格推出
  `A_r = 0`
  `A_u = κ/X`
- 因而新的收紧理解是：
  - 单个 `C,D,E,F` 当然可以比 `1/X` 发散得更快
  - 但真正受“恢复 HJ + 恢复连续性”共同控制的有效组合，不可能是 `1/X^2` 或别的阶数，只能精确等于 `κ/X`
- 在同一假设下，
  `ρ~ = ρ X sqrt(R) / κ`
  其中 `R = det(g~^(..))/det(g^(..))`
  因而：
  - 若 `R = o(1/X^2)`，则 `ρ~ -> 0`
  - 若 `R ~ const / X^2`，则 `ρ~` 有限非零
  - 若 `R >> 1/X^2`，则 `ρ~` 发散
- 进一步给出了一个 genuinely mixed 的显式例子：
  在 `Y ≠ 0` 的 generic patch 上取
  `C=1`
  `E=αY/Δ`
  `F=-αX/Δ`
  并由 `A_u=κ/X` 解出
  `D = κ/X^2 - 1/X - αY^2/(XΔ)`
  则
  `R = (κ/X)(1-α)`
  `ρ~ = ρ sqrt[(1-α)X/κ]`
  因而在 `(1-α)κ/X > 0` 的 patch 上，混合 ansatz 也能实现 `ρ~ -> 0`
- 当前阶段结论因此进一步收紧为：
  `ρ~ -> 0` 在 rank-1 和 genuinely mixed 两类里都可以出现，但这只是 current-level regularization；只要还坚持点态恢复连续性方程和固定非零 `κ`，metric/HJ 级主奇异 `A_u = κ/X` 仍然不会消失

## 轮回 2 · 2026-04-28 14:49:00 CST · 澄清“散度等价”与“点态电流匹配”的区别

- 当前继续整理一个关键概念层级：
  - 点态电流匹配（strong current matching）指的是
    `\tilde J^μ = J^μ`
  - 散度等价（divergence-level equivalence）指的是只要求
    `∂_μ \tilde J^μ = 0` 与 `∂_μ J^μ = 0` 作为方程等价，而不要求 `\tilde J^μ` 与 `J^μ` 逐点相同
- 其中
  `J^μ = sqrt(|g|) ρ u^μ`
  `\tilde J^μ = sqrt(|g~|) \tildeρ \tilde g^(μν) u_ν`
- 一个很自然的充分条件是：
  `\tilde J^μ = J^μ + K^μ`
  且
  `∂_μ K^μ ≡ 0`
  例如 `K^μ = ∂_ν A^[μν]` 这类由反对称超势构成的改进项
- 更弱的等价条件也可以写成：
  `∂_μ \tilde J^μ = Ω ∂_μ J^μ`
  其中 `Ω` 是处处有限非零的标量函数；这样两条连续性方程的零集相同
- 这层区分的意义是：
  先前关于混合 ansatz 的严格结论 `A_r=0`、`A_u=κ/X`，是在“点态电流匹配”这个更强条件下推出的
  若只要求散度等价，则 `\tilde J^μ` 允许带有不逐点等于 `J^μ` 的附加项，理论上有机会绕开该强约束

## 轮回 2 · 2026-04-28 15:02:00 CST · 混合 ansatz 在散度等价下的主方程

- 已把“只要求散度等价而不要求点态电流匹配”的情形系统推了一遍，并整理成：
  `research-notes/069-混合ansatz下散度等价的主方程与绕开kappa除以X的条件.md`
- 在一般混合 ansatz 下，记
  `A_u = C + DX + EY`
  `A_r = EX + FY`
  `s = sqrt(|g~|)ρ~`
  以及原电流 `J^μ = j u^μ`, `j = sqrt(|g|)ρ`
- 变换后 HJ 方程 `g~^(μν)u_μu_ν = κ` 给出
  `X A_u + Y A_r = κ`
  因而
  `A_u = (κ - Y A_r)/X`
- 变换后电流可写成
  `J~^μ = (sκ/X)u^μ + B W^μ`
  其中
  `W^μ = X r^μ - Y u^μ`
  `B = s A_r / X`
  且 `u_μ W^μ = 0`
- 若只要求散度等价，即
  `∂_μ(J~^μ - J^μ) = 0`
  则得到主方程
  `∂_μ{ [sκ/X - j]u^μ + B(Xr^μ - Yu^μ) } = 0`
- 这一步非常关键，因为它表明：
  - `A_u = κ/X` 不是 HJ 单独推出的
  - 它是 `HJ + 点态电流匹配` 共同推出的强结论
  - 一旦只要求散度等价，`A_r` 不再被迫为零，`A_u` 也不再被迫等于 `κ/X`
- 绕开 `1/X` 的条件也已明确：
  - 若 `Y ≠ 0`，可取
    `A_r = κ/Y + Xχ`
    则
    `A_u = -Yχ`
    从而 `A_u` 在 `X -> 0` 时保持有限
  - 但若 `Y = 0`，HJ 直接退化成 `X A_u = κ`，此时仍然无路可绕，必须有 `A_u = κ/X`
- 当前阶段结论：
  混合 ansatz 的“`κ/X` 主奇异”并非在所有等价概念下都不可避免；
  它是强等价（点态电流匹配）下的 no-go，而在弱等价（散度等价）下，`Y ≠ 0` patch 上出现了真正的逃逸窗口

## 轮回 2 · 2026-04-28 15:18:00 CST · 关于“ρ~ 为零时速度无定义”的分支化判断

- 用户提出一个重要的物理解读：
  若在 `X=0` 邻域有 `ρ~ -> 0`，则即使 `g~` 表象下的速度场发散，也不对应真实粒子流，因为该处本来就没有粒子密度
- 当前判断做了更细的分层：
  1. 这一物理解读本身是合理的
  2. 但它只适用于那些确实满足 `ρ~ -> 0` 的分支
  3. 它不能自动推广到所有混合 ansatz
- 已新增专门笔记：
  `research-notes/070-关于rho波浪为零时速度无定义的物理解读与混合投影反例.md`
- 关键复核结果：
  - rank-1 `C=1` 分支上，
    `ρ~ = ρ sqrt(|X/κ|) -> 0`
    因而“零密度边界、速度无需定义”的解释是自洽的
  - 对未缩放的 raw 混合投影块
    `D=Z, E=-Y, F=X`
    若要求 `g~^(μν)u_μu_ν = κ`，则有
    `C+Δ = κ/X`
    `R = C^2 (κ/X)^2`
    强连续性匹配给出
    `ρ~ = ρ |C|`
    因此 `ρ~ = 0` 并不是自动推出的；只有 `C -> 0` 时才成立
    若原变量正则使 `Δ` 有限，则
    `C = κ/X - Δ ~ κ/X`
    反而导致 `ρ~` 发散
  - 对当前项目主用的 `C=1` 投影型分支，
    `R = (κ/X)^2`
    强连续性匹配给出
    `ρ~ = ρ`
    因而这里也没有 `ρ~ -> 0`
- 因而当前收紧后的结论是：
  “`ρ~ -> 0` 可以为速度发散提供物理解读”这点保留；
  但用 `D=Z,E=-Y,F=X` 或当前 `C=1` 投影型主分支来证明 `ρ~ = 0` 并不成立，除非再额外选取会使 `C -> 0` 的特殊分支

## 轮回 2 · 2026-04-28 15:26:00 CST · 修正行列式比值记号歧义

- 用户指出一个真实的记号问题：若把 `R` 读成协变度规行列式比值，则确实应有
  `sqrt(|g~|)/sqrt(|g|) = sqrt(R)`
- 重新核对后确认：
  - 我之前在若干笔记里把 `R` 用成了逆度规行列式比值
    `R_inv = det(g~^(..))/det(g^(..))`
  - 在这个定义下，
    `sqrt(|g~|)/sqrt(|g|) = 1/sqrt(R_inv)`
  - 若改用协变度规行列式比值
    `R_cov = det(g~)/det(g) = 1/R_inv`
    则
    `sqrt(|g~|)/sqrt(|g|) = sqrt(R_cov)`
- 已把 `research-notes/067`、`068`、`070` 中相关记号统一改成 `R_inv`，避免再把协变和逆度规行列式比值混淆
- 关键复核结论：
  这个记号修正不会改变那些分支上 `ρ~` 的最终结论；例如 raw 混合投影块里
  `ρ~ = ρ sqrt(R_inv)/(C+Δ) = ρ/(sqrt(R_cov)(C+Δ)) = ρ|C|`
  两种记法完全一致

## 轮回 2 · 2026-04-28 15:37:00 CST · 提升到测度密度口径后的统一解释

- 用户进一步指出：即使某些分支上 `ρ~` 本身不趋于 0，只要
  `sqrt(|g~|)ρ~ -> 0`
  ，那么在 `g~` 表象里该点的“粒子数密度”仍趋于 0，因此速度无需定义
- 当前复核后确认：这条修正后的解释是成立的，而且比讨论裸 `ρ~` 更一般
- 在强连续性匹配且
  `g~^(μν)u_ν = A_u u^μ`
  `A_u = κ/X`
  的所有分支里，都有统一关系
  `sqrt(|g~|)ρ~ A_u = sqrt(|g|)ρ`
  因而
  `sqrt(|g~|)ρ~ = (|X|/|κ|) sqrt(|g|)ρ -> 0`
  （假设原表象 `sqrt(|g|)ρ` 有限）
- 这说明：
  - rank-1 `C=1` 分支：`ρ~ -> 0`，同时测度密度也趋零
  - `C=1` 投影型主分支：`ρ~ = ρ`，但测度密度仍趋零
  - raw 混合投影块：`ρ~ = ρ|C|`，但测度密度同样趋零
- 当前收紧后的统一解释：
  `X=0` 可被理解为 `g~` 表象中的零测度密度边界；该处速度样对象不必有良好定义，只要电流保持有限，方程仍可有意义
- 但也保留一个边界：
  这解决的是 current-level / interpretational 问题，并不自动证明 `g~` 图册在 `X=0` 处本身已光滑正则

## 轮回 2 · 2026-04-28 16:12:00 CST · 三种引力作用量的全动力学方程与当前数值可信层级

- 用户把讨论切回引力作用量部分，要求重新统一三种总作用量：
  - `A: sqrt(-g) R[g]`
  - `B: sqrt(-g~) R[g~]`
  - `C: sqrt(-g~) f(R~)` / `sqrt(-g~) F(R~) R~`
- 本轮已把“形式上的全动力学方程”和“当前真正可信的数值结果”明确分层，并整理成新笔记：
  `research-notes/072-三种引力作用量在2+1d下的全动力学方程与当前数值状态.md`
- 当前统一约化框架仍是：
  - `3+1d -> 2+1d`
  - `y` 方向 Killing
  - `g~` 表象采用
    `ds^2 = e^{2τ}(dt^2-dx^2-dz^2) - e^{2β}dy^2`
  - 共同物质方程写成
    `HJ: S_t^2 - S_x^2 - S_z^2 = m^2 e^{2τ}`
    与
    `n_t = -∂_x[(n/E)S_x] - ∂_z[(n/E)S_z]`
- 当前正式可写成显式 PDE 并可数值积分的，是：
  - B 支：`(τ, β, n, S)` 的二阶/一阶耦合系统
  - C 支：`(τ, β, Φ, n, S)` 的辅助场二阶/一阶耦合系统
- B 支当前显式方程为：
  - `□τ = -1/2[(∂τ)^2 + (∂β)^2 + μ]`
  - `□β = -∂β·∂τ - (∂β)^2 - μ/2`
  - `μ = m^2 e^{2τ}ρ/M_P^2`
- C 支当前显式方程为：
  - `□τ = (C - 2A - 2ΦB)/(6Φ)`
  - `□β = -B - 2□τ`
  - `□Φ = -A - 2Φ□τ`
  - 其中 `A,B,C` 为由 `τ,β,Φ,U(Φ),μ` 组成的标准辅助场组合
- A 支当前已经明确：
  - 形式上的总作用量和形式场方程是清楚的
  - 但在 `2+1d` 下一旦把引力项写在 `g` 而物质写在 `g~`，若不额外指定一个明确的 `g <-> g~` reduced closure，就没有唯一的闭合 PDE
  - 因此 A 支当前还不能和 B/C 一样在同一闭合层级上做“正式 full backreacted 数值器”
- 已重新核对当前最可信的三支数值比较结果：
  - 可信口径仍然是 shared-baseline 同初值比较
  - 即：A 支先可信演化，再把同一条物质主干历史喂给 B/C 几何子系统
- 当前最可信代表结果：
  - `lambda=1`：
    `B_vs_A_relL1 = 4.2155e-06`
    `C_vs_A_relL1 = 1.8925e-05`
    `C_vs_B_relL1 = 1.4709e-05`
  - `lambda=3`：
    `B_vs_A_relL1 = 1.2649e-05`
    `C_vs_A_relL1 = 4.5499e-05`
    `C_vs_B_relL1 = 3.2850e-05`
- 同时重新核对了旧的 “same-initial full Cauchy” 原型：
  `kg_examples/simulate_abc_cauchy_same_initial.py`
  当前输出爆炸到 `1e7 ~ 1e8` 量级，不能作为正式物理结果引用
- 当前结论已收紧为：
  1. 三种引力作用量的形式全方程已经清楚
  2. B/C 的 `2+1d` 显式闭合 PDE 已清楚
  3. A 的 `2+1d` full closure 仍缺一个额外 reduced closure，尚未进入真正可信的 full numerical branch
  4. 因此当前对外可正式报告的三支数值结果，仍应限定为 shared-baseline 口径，而不能写成“三支完全自洽全反作用都已完成”

## 轮回 2 · 2026-04-28 16:26:00 CST · 变换改动后旧数值脚本不能直接沿用

- 用户指出了一个关键问题：当前若采用新的混合变换/新的 `ρ~` 关系，就不能把旧 shared-baseline 数值结果直接拿来当新变换下的物理结果
- 已核对旧脚本 `kg_examples/compare_abc_same_initial_evolution.py` 的硬编码位置：
  - 在记录与作图层，直接使用
    `rho_b_outer = rho_a_outer * exp(-(beta_b + tau_b))`
    `rho_c_outer = rho_a_outer * exp(-(beta_c + tau_c))`
  - 在几何源项层，也直接把 `rho_a_outer` 喂给 `rhs_branch_b/c`
- 这说明旧脚本并不是“只差解释口径”，而是把旧分支的密度重构和旧共享主干耦合方式硬编码进了：
  1. B/C 几何源项
  2. B/C 后处理 `ρ~` 重构
  3. 后续 AB / ABC observable 提取
- 当前结论：
  - 旧几何推进器的波动更新、吸收层、边界处理这些数值骨架可以借鉴
  - 但只要变换换成新的 `D=Z,E=-Y,F=X` 一类混合分支，旧脚本就必须至少重写：
    - `ρ~` / `sqrt(|g~|)ρ~` 的重构
    - B/C 的有效源项输入
    - 相关 observables 的定义与比较口径
- 因而当前旧的 shared-baseline 结果只能保留为“旧变换下的方法学 benchmark”，不能直接当作“新变换下的数值结果”

## 轮回 2 · 2026-04-28 16:43:00 CST · 新混合变换 shared-backbone 源项替换的第一步审计

- 已按“先走一步”的要求，新增诊断脚本：
  `kg_examples/audit_mixed_transform_shared_backbone.py`
- 这一步暂时不重写整套几何推进，而是先审计：
  - 当前 exact `A` benchmark 上的 `Q`
  - 由 `X = m^2 + Q` 得到的 `X`
  - 以及新混合变换下最自然的第一替换对象
    `N~ = sqrt(|g~|)rho~ = (|X|/|kappa|) rho_A`
    （带符号的 `X/kappa` 仅作 branch diagnostic）
- 同步新增正式笔记：
  `research-notes/073-新混合变换下shared-baseline源项替换的第一步审计.md`
- 当前审计结果很明确：
  - 在 `t/T = 0, 0.25, 0.5, 0.75, 1` 五个代表时刻上，
    `X` 在主支撑区内部都已经大量过零
  - `X<0` 的支撑区占比大约在 `0.44 ~ 0.48` 之间
  - 因而不能直接把旧脚本里的 `rho_A` 换成 `(X/kappa) rho_A` 就当作新的 shared-baseline 源项；正密度口径下应先看 `|X/kappa| rho_A`
- 当前收紧后的结论：
  - 新变换下的第一真正 blocker 不是数值器本身，而是 source prescription
  - 在真正重跑 B/C 之前，必须先明确 `X<0` 区域如何进入源项：
    - signed
    - abs
    - positive-part
    - 或其它明确的 branch prescription

## 轮回 2 · 2026-04-29 03:49:55 CST · 正测度密度口径修正

- 用户指出：
  `sqrt(|g~|)ρ~` 的正根口径下，前面写成 `(X/kappa)` 是不严谨的；若把 `R` 统一理解成协变度规体积元关系，就必须显式保留绝对值/正根
- 已重新核对并确认：
  - `C^2(C+Δ)^2` 一直对应的是逆变度规行列式比值
  - 体积元 `sqrt(|g~|)` 则由协变度规行列式给出
  - 在强匹配且正密度的单个物理 patch 上，应有 `κ/X > 0`
  - 因而真正的正测度密度关系应写成
    `sqrt(|g~|)ρ~ = (|X|/|κ|) sqrt(|g|)ρ`
  - 带符号的 `(X/kappa)` 只能保留成 branch diagnostic，不能直接当作物理正密度
- 已同步修正：
  - `research-notes/071-测度密度趋零与X等于0处速度无定义的统一解释.md`
  - `research-notes/073-新混合变换下shared-baseline源项替换的第一步审计.md`
- 这一步把 shared-baseline 的 blocker 又收紧了一层：
  - 当前不仅仅是“要不要选 signed / abs prescription”
  - 而是：只要 `X` 在主支撑区内部变号，就不存在一个覆盖整个主支撑区的单一固定-`κ` 正密度强匹配 patch

## 轮回 2 · 2026-04-29 04:08:58 CST · 放弃旧对角数值器并建立一般 `2+1 ADM` 新方法

- 用户明确要求：不要再和旧数值计算混在一起，而要为当前新混合变换重新发展一套 honest 的全动力学方法
- 已重新核对现有 full-dynamics 脚本：
  - `kg_examples/simulate_2p1_full_dynamics.py`
  - `kg_examples/simulate_abc_cauchy_same_initial.py`
  - 以及 shared-baseline 脚本 `kg_examples/compare_abc_same_initial_evolution.py`
- 当前确认的关键结构事实：
  - 旧 full-dynamics 数值器都把 `g~` 的 `t-x-z` 部分写死成对角 conformal gauge
  - 但当前混合变换在 `±45°` 高斯波包初值上会 generically 生成 `t-x`、`t-z`、`x-z` 非对角分量
  - 因此旧的 `τ,β` 数值器不再对应当前模型，不能继续拿来宣称是在算新变换下的全动力学
- 为此已新增一份正式方法笔记：
  `research-notes/074-新混合变换下的一般2加1ADM全动力学方法.md`
- 这份笔记完成了新方法的第一步：
  - 定义新的 `2+1 ADM` 度规
    `ds^2 = N^2 dt^2 - h_ij (dx^i + N^i dt)(dx^j + N^j dt) - e^{2β} dy^2`
  - 在这一一般度规下，把共同物质作用量重写成 `(n,S)` 的守恒系统
  - 给出
    `S_t = N^i S_i - N E`
    `n_t = -∂_i [ n N^i + (N n/E) h^(ij) S_j ]`
    `ρ = n / (e^β sqrt(h) E)`
- 同步新增第一块可运行代码：
  `kg_examples/adm_matter_2p1.py`
- 已做最小自检：
  - 在零 shift、平直 `h_ij = δ_ij`、`β=0` 的测试下
  - `ρ` 可由 `n` 精确回收
  - `S_t` 退回正频支
  - `n_t` 在简单线性相位测试上表现正常
- 当前新的真正推进顺序已收紧为：
  1. 先把三支共同的物质推进器切到一般 `ADM`
  2. 再分别重推 B/C 的 `ADM` 引力演化方程
  3. 最后为 A 支写出一般 `ADM` 下的 `g↔g~` 约束闭合

## 轮回 2 · 2026-04-29 04:22:47 CST · 完成 B/C 两支的一般 `ADM` 重写

- 已新增正式笔记：
  `research-notes/075-B与C两支在一般2加1ADM度规下的引力动力学.md`
- 这一步完成的内容：
  - 对 `B` 支，引入三维 Einstein 形式度规 `q_ab = e^{2β} h~_ab`
  - 把约化作用量重写成
    `2+1 Einstein + β + Bohm 物质`
  - 写出共同的物质方程、`β` 方程以及三维 Einstein 方程
  - 再把它组织成标准 `ADM` 初值问题的“约束方程 + 演化方程”结构
- 对 `C` 支：
  - 在同一个 `q_ab` 上完成辅助场 `Φ` 的一般作用量重写
  - 识别出它是一个 `2+1` 标量-张量系统
  - 写出 `Φ` 方程、`β` 方程和带非最小耦合的 Einstein 方程
  - 并说明其 `ADM` 变量应是
    `h_ij, K_ij, β, Π_β, Φ, Π_Φ`
- 当前阶段的清晰结论：
  - 三支里 B/C 两支现在都已经进入一般 `ADM` 语言
  - 共同物质骨架也已独立出来
  - 因而当前唯一真正剩下的理论瓶颈就是：
    A 支如何把 `g ↔ g~` 的当前混合变换翻成一般 `ADM` 变量之间的代数闭合关系

## 轮回 2 · 2026-04-29 04:46:32 CST · 修正 A 支的基本作用量定义

- 用户指出：若 A 支引力作用量取 `\sqrt{-g}R[g]`，则物质部分也应统一用 `g` 写；否则会把“基本作用量的定义”和“事后几何化解释”混成两层。
- 已接受这一修正，并新增正式笔记：
  `research-notes/076-A支应作为g上Einstein-KG系统而非g与g波浪混合基本作用量.md`
- 当前统一口径：
  - A 支 = `g` 上的 Einstein-Klein-Gordon 系统，以 `(\rho,S)` 变量表达
  - `g~` 只是在求出 `(g,\rho,S)` 后，由当前混合变换事后重构出来的有效几何
  - 因此 A 支后续的全动力学任务，应改写为：把 `g` 上的 Einstein-KG 系统写成一般 `2+1 ADM` 约束—演化方程
- 这一步也意味着：此前把 A 支描述成“引力在 `g` 上、物质在 `g~` 上”的口径作废，不再使用

## 轮回 2 · 2026-04-29 05:18:11 CST · 统一切到四维 ADM 并启动 A/B 原型

- 新增正式笔记：
  - `research-notes/077-三支理论在y方向Killing对称下的统一四维ADM初值形式.md`
  - `research-notes/078-A与B两支在统一四维ADM下的最小同步规范原型.md`
- 新增代码：
  - `kg_examples/adm_ykilling_geometry.py`
  - `kg_examples/adm_sources_abc.py`
  - `kg_examples/simulate_ab_adm_synchronous.py`
- 这一轮的关键修正：
  - 不再把三支的数值方法建立在 `2+1` 共形重参数化上
  - 统一改用四维 `ADM` 初值问题，只保留 `y` 方向 Killing 对称
- 当前已完成的内容：
  - `A` 支：用复标量 `(\phi_1,\phi_2)` 与共轭变量 `(\Pi_1,\Pi_2)` 写成标准 Einstein-Klein-Gordon 初值系统
  - `B` 支：把 Bohm 型 `(\rho,S)` 物质的一阶守恒系统与 Einstein 几何源项直接接入统一 `ADM` 语言
  - `A/B` 两支在最小同步规范 `N=1, N^x=N^z=0` 下已经真正开始时间推进
- 已跑的最小数值验证：
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 3 --dt 0.001`
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 10 --dt 0.0005`
  - 两次运行都正常结束，短时间内 `A/B` 两支都保持有限
- 当前阶段结论：
  - 新的统一四维 `ADM` 方法已经不只是公式，而是已进入可运行原型
  - 下一步应补上约束监测，并把 `C` 支辅助场 `(\Phi,\Pi_\Phi)` 与有效源项接入同一原型

## 轮回 2 · 2026-04-29 05:39:27 CST · 将 C 支接入统一原型并加入最小约束监测

- `kg_examples/adm_ykilling_geometry.py` 已补充：
  - `scalar_hessian_ykilling`
  - 更方便的 `\Phi` 空间 Hessian / Laplacian 工具
- `kg_examples/adm_sources_abc.py` 已补充：
  - `c_branch_effective_sources`
  - 用于把 `C` 支写成“等效 Einstein 源 + 辅助场演化”的第一版实现
- `kg_examples/simulate_ab_adm_synchronous.py` 已扩展为：
  - 同时推进 `A/B/C` 三支
  - 同时记录最基本的 Hamilton 约束残差
- 已运行最小测试：
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 3 --dt 0.0005`
- 当前最关键结果：
  - `A/B/C` 三支现在都已在同一份同步规范原型里真正往前推进
  - `A/B` 的短时间几何量级和 Hamilton 残差相近
  - `C` 支虽然可运行，但 Hamilton 残差明显偏大，已单独记录在：
    `research-notes/079-C支已接入统一原型但约束残差仍显著偏大.md`
- 当前阶段判断：
  - 下一步应优先检查 `\Phi=1,\Pi_\Phi=0` 是否数值退回 `B` 支
  - 并继续校正 `\nabla_\mu\nabla_\nu\Phi` 的 `ADM` 投影和 Jordan 形式的约束实现

## 轮回 2 · 2026-04-29 06:19:11 CST · C 支平直 `\Phi=1` 退回 B 支的校准完成

- 已修正 `kg_examples/adm_sources_abc.py`：
  - 将 `U(\Phi)`、`U_\Phi(\Phi)` 改成与旧稳定实现一致、满足
    `U(1)=0`、`U_\Phi(1)=0`
  - 在有效空间应力中引入
    `\Box\Phi_{\mathrm{corr}} = \Box\Phi - T^{(B)}/(3M_P^2)`
    的平直 `\Phi=1` 扣除
- 已修正 `kg_examples/simulate_ab_adm_synchronous.py`：
  - `\Pi_\Phi` 演化只使用 `\Box\Phi_{\mathrm{corr}}`
  - 边界阻尼后强制保持 `\Phi=1`、`\Pi_\Phi=0` 边界值，避免数值上人为打破冻结支
- 已做逐项源项回归：
  - 在标准双束初值上比较 `B` 与 `C`，当 `\Phi=1,\Pi_\Phi=0` 时
    `energy`、`mom_x`、`mom_z`、`s_xx`、`s_xz`、`s_zz`、`s_yy`
    的差值最大值全部为 `0`
- 已重跑最小同步规范原型：
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 3 --dt 0.0005`
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 10 --dt 0.0005`
- 当前结果：
  - `C` 支的 `r3` 历史与 `B` 支在 `3` 步和 `10` 步测试中完全重合
  - `C` 支的 Hamilton 约束残差历史与 `B` 支完全重合
  - `C` 支的 `rho` 与 `phi` 也保持在冻结支上
- 已新增正式笔记：
  - `research-notes/080-C支在统一四维ADM原型下已严格退回B支.md`
- 当前阶段判断：
  - `C` 支的平直辅助场归一化和有效 Einstein 源平直扣除已经校准完成
  - 后续真正的剩余问题已收缩成：
    1. 如何为 `C` 支构造非平凡而一致的 `\Phi` 初值
    2. 如何为三支求解满足约束的共同初始几何

## 轮回 2 · 2026-04-29 07:01:42 CST · 完成共同初始几何的第一版构造

- 已新增代码：
  - `kg_examples/adm_initial_data_abc.py`
- 这份新模块完成了三件事：
  1. 从同一组 `(\rho,S,g)` 在时间对称几何假设下正确重构 `A` 支的 `(\phi_1,\phi_2,\Pi_1,\Pi_2)` 初值
  2. 在共形平坦 ansatz 下，用平均能量源 `(\mathcal E_A+\mathcal E_B)/2` 解共同空间度规
  3. 用平均动量源 `(\mathcal P^{(A)}+\mathcal P^{(B)})/2` 解最小 York 型共同外曲率
- 已做关键对比：
  - 在同一组 `\rho,S` 和平直几何下，`A/B` 的源差异很大：
    - `energy` 相对 `L1` 差约 `0.645`
    - `mom_x` 相对 `L1` 差约 `2.997`
    - `s_xx` 相对 `L1` 差约 `0.791`
    - `s_xz` 相对 `L1` 差约 `1.230`
  - 因而“三支完全同一组 `g` 初值”不能同时严格满足三支约束，只能在弱场下近似成立
- 已在共同初始几何上比较 branchwise 共形因子量级：
  - `A` 源：`omega_max ≈ 5.66e-3`
  - `B` 源：`omega_max ≈ 7.25e-3`
  - 平均源：`omega_max ≈ 6.45e-3`
- 已把这套共同初始几何接回
  `kg_examples/simulate_ab_adm_synchronous.py`
  ，并新增：
  - 正确的 `A` 支复标量动量初值
  - 平均源共同空间度规
  - 平均源共同外曲率
  - `C` 支的 `r3_proxy` 非平凡初值入口
- 已跑最小测试：
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 3 --dt 0.0005 --c-phi-mode flat`
  - `python3 kg_examples/simulate_ab_adm_synchronous.py --steps 3 --dt 0.0005 --c-phi-mode r3_proxy`
- 当前结果：
  - 共同初始几何下，三支 `t=0` 的 Hamilton 残差已降到 `~8e-2`
  - `r3_proxy` 模式下，
    `Phi_min ≈ 0.99999999992`
    ，说明在当前 `M_P=300, ell=0.02` 的弱场参数下，这条最小非平凡 `\Phi` 初值几乎仍是平直冻结支
- 已新增正式笔记：
  - `research-notes/081-共同初始几何的第一版构造与AB源差异.md`
- 当前阶段判断：
  - `C` 支的下一个问题不再是校准，而是“怎样选一个真正会让 `C` 与 `B` 分开的非平凡 `\Phi` 初值”
  - 共同初始几何的下一个问题则是“显式计算动量约束残差，判断当前平均源近似是否已够用”

## 轮回 2 · 2026-04-29 08:32:10 CST · 更正 A 参考用法并完成 B/C 完整耦合同步规范测试

- 已先把 `kg_examples/simulate_bc_geometry_from_a_reference.py` 的时间推进从二阶 Heun 升到四阶 Runge-Kutta，并分别重跑：
  - `10` 步，`dt=5e-4`
  - `20` 步，`dt=2.5e-4`
  - `40` 步，`dt=1.25e-4`
- 结果显示：
  - 把积分器从二阶升到四阶后，`B/C` 的时间步敏感性几乎不变
  - 同一物理时长下，`B` 支几何偏离仍从约 `0.452` 增长到 `0.700`、再到 `0.910`
  - 因而当前问题不在积分阶数
- 已进一步把 `A` 支线性化弱场 `omega` 作为 `B/C` 的初始空间度规估计接入同一脚本重跑
- 结果显示：
  - 使用 `omega` 初始几何后，`10` 步与 `20` 步结果和此前几乎重合
  - 因而当前问题也不主要来自“完全平直”的 `B/C` 初始空间度规
- 在继续分析后，已明确修正一个更根本的方法口径：
  - 用户要求的是“`A` 支先给定初始 `(\rho,S)`，然后 `B/C` 自己做完整耦合演化”
  - 因此前一轮“把整条 `A` 历史拿来外部驱动 `B/C`”的做法只能算诊断，不是正式口径
- 已新增新脚本：
  - `kg_examples/simulate_bc_from_a_initial_data.py`
- 这份新脚本做的是：
  - 用平直背景 `A` 支给出初始 `(\rho,S)` 与初始弱场 `omega`
  - 用该 `omega` 初始化 `B/C` 的空间度规
  - 之后让 `B/C` 用各自的 `(n,S,g~)` 完整耦合推进
  - 其中 `n = e^\beta \sqrt{h}\rho E` 作为真正的守恒密度变量
- 当前已完成最小测试：
  - `5` 步，`dt=5e-4`
  - `10` 步，`dt=5e-4`
  - `20` 步，`dt=2.5e-4`
- 当前结果：
  - 初始 Hamilton 残差仅约 `7.5e-2`
  - 但一进入演化，残差与几何偏离仍迅速增长：
    - `10` 步后 `max|h_xx-1| ≈ 0.452`，`max|H| ≈ 35.4`
    - `20` 步后 `max|h_xx-1| ≈ 0.700`，`max|H| ≈ 130`
  - 当前 `C` 支仍取 `Phi=1, Pi_Phi=0`，因此与 `B` 严格重合
- 已新增正式笔记：
  - `research-notes/083-B与C两支从A支初始数据出发的完整耦合同步规范测试.md`
- 当前阶段判断：
  - 现在已经可以排除两个较弱解释：
    1. 问题不是由二阶时间推进造成的
    2. 问题不主要来自把 `B/C` 初始空间度规取成平直
  - 当前真正下一步应收缩到：
    - 同步规范是否引入了主要约束增长
    - `B` 支几何演化的约束传播是否需要更合适的规范或阻尼写法

## 轮回 2 · 2026-04-29 09:11:42 CST · 定位到 B/C 原型的主问题是几何边界条件而非演化方程本身

- 已在
  `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_bc_from_a_initial_data.py`
  中补入：
  - `B` 支初始动量源对应的最小 York 型外曲率
  - Hamilton 约束与动量约束残差的并列监测
- 补跑结果显示：
  - 初始动量约束残差已降到 `~1e-3`
  - 但若仍使用旧的“对所有变量统一乘法阻尼”边界条件，`10/20` 步的巨大几何偏离几乎不变
- 进一步逐步打印 `h_{xx}` 与 `k_{xx}` 的量级后，确认先前的主要爆长来自边界条件实现错误：
  - 几何变量 `h_xx,h_xz,h_zz,beta,k_xx,k_xz,k_zz,k_beta` 被错误地和物质变量一起做乘法阻尼
  - 这会把本应回到背景值的边界几何强行压向 `0`
  - 从而伪造出巨大曲率与巨大约束残差
- 已改正边界条件：
  - 物质变量 `n,S` 仍做边界阻尼
  - 几何变量改为回到背景值：
    - `h_xx -> 1`
    - `h_xz -> 0`
    - `h_zz -> 1`
    - `beta -> 0`
    - `k_xx,k_xz,k_zz,k_beta -> 0`
  - `C` 支边界继续取 `Phi -> 1, Pi_Phi -> 0`
- 修正后数值结果出现明显转折：
  - `5` 步，`dt=5e-4`：
    - `max|h_xx-1| ≈ 7.76e-3`
    - `max|H| ≈ 2.59e-1`
    - `max|M| ≈ 2.62e-3`
  - `10` 步，`dt=5e-4`：
    - `max|h_xx-1| ≈ 7.76e-3`
    - `max|H| ≈ 2.58e-1`
    - `max|M| ≈ 5.20e-3`
  - `20` 步，`dt=2.5e-4`：
    - `max|h_xx-1| ≈ 7.76e-3`
    - `max|H| ≈ 2.58e-1`
    - `max|M| ≈ 5.44e-3`
  - 说明：
    1. 先前的“时间步越细越爆”并非方程本身问题，而主要是几何边界条件伪影
    2. 当前 `B` 支几何响应在弱场范围内稳定，且时间步敏感性大幅消失
- 已把 `C` 支的 `r3_proxy` 非平凡初值入口接回当前稳定版本，并跑了
  - `200` 步，`dt=2.5e-4`
  - 在当前 `ell=0.02` 下得到
    - `max|Phi-1| ≈ 2.6e-8`
    - `C` 相对 `B` 的几何差异仍极小
- 已新增正式笔记：
  - `research-notes/084-B与C两支同步规范原型的主要问题来自几何边界条件而非时间推进.md`
- 当前阶段判断：
  - `B/C` 原型已经从“明显伪不稳定”转入“弱场稳定演化”阶段
  - 下一步不该再优先查时间推进器，而应转向：
    1. 延长物理时间窗口
    2. 设计更强但自洽的 `C` 支初始 `Phi`
    3. 量化 `A` 支平直背景近似的系统误差

## 轮回 2 · 2026-04-29 09:38:21 CST · `C` 支 `ell` 扫描完成，默认弱场下确实几乎退回 `B`

- 已在
  `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_bc_from_a_initial_data.py`
  中加入统一的 `B-C` 终态差异摘要：
  - `max|h_xx^{(C)}-h_xx^{(B)}|`
  - `max|h_xz^{(C)}-h_xz^{(B)}|`
  - `max|h_zz^{(C)}-h_zz^{(B)}|`
  - `max|beta^{(C)}-beta^{(B)}|`
  - `max|rho^{(C)}-rho^{(B)}|`
  - 以及 Hamilton / 动量约束残差差距
- 在固定
  `steps=200, dt=2.5e-4, c_phi_mode=r3_proxy`
  下完成 `ell` 扫描：
  - `ell=0.02`：
    - `max|Phi-1| ≈ 2.61e-8`
    - `max|h_xx^C-h_xx^B| ≈ 2.81e-8`
  - `ell=0.05`：
    - `max|Phi-1| ≈ 1.45e-7`
    - `max|h_xx^C-h_xx^B| ≈ 1.60e-7`
  - `ell=0.1`：
    - `max|Phi-1| ≈ 5.46e-7`
    - `max|h_xx^C-h_xx^B| ≈ 6.09e-7`
  - `ell=0.2`：
    - `max|Phi-1| ≈ 1.74e-6`
    - `max|h_xx^C-h_xx^B| ≈ 2.03e-6`
  - `ell=0.5`：
    - `max|Phi-1| ≈ 7.70e-6`
    - `max|h_xx^C-h_xx^B| ≈ 7.04e-6`
  - `ell=1.0`：
    - `max|Phi-1| ≈ 2.00e-4`
    - `max|h_xx^C-h_xx^B| ≈ 1.21e-4`
- 当前结论：
  - 默认 `ell=0.02` 时，`C` 在当前双束弱场问题上数值上确实几乎退回 `B`
  - 要看到明显 `C-B` 分离，至少需要 `ell` 进入 `0.5~1.0` 量级，或换更强的 `Phi` 初值
- 又做了一步符号检查：
  - 若试图用纯局域代数条件
    `2u(Phi)-Phi u_Phi(Phi)=trace(T_B)/M_P^2`
    直接构造非平凡 `Phi` 初值，则对当前模型
    `0 < Phi <= 1` 时左边非正、右边正，不能成立
  - 因此下一步若要真正把 `C` 拉开，不应继续找纯局域代数 `Phi` 公式，而应改为求解带空间拉普拉斯的椭圆型初值方程
- 已新增正式笔记：
  - `research-notes/085-C支在稳定同步规范原型下的ell扫描与Phi初值含义.md`

## 轮回 2 · 2026-04-29 09:47:58 CST · 当前边界条件的严格性判断

- 已重新检查当前物质与几何边界条件的实现：
  - 物质变量 `n,S` 仍使用乘法海绵层 `edge_damp`
  - 几何变量与 `Phi` 使用背景常数边界值
- 当前判断已收紧为：
  1. 物质子系统存在更规范的特征边界条件：应按边界法向特征速度区分 inflow / outflow，而不是统一乘法阻尼
  2. 当前同步规范 `ADM` 几何系统并没有现成的“严格可证明良定”的边界条件；若要严格，应切换到调和规范（harmonic gauge）或广义调和 / Z4c / BSSN 这类已有约束保持边界条件的形式
  3. 因而当前 `edge_damp + 背景常数几何边界` 只能视为数值上可用的近似边界，不应称为严格边界条件
- 外部主参考：
  - Kreiss-Reula-Sarbach-Winicour 2007：harmonic Einstein 初边值问题的良定性
  - Seiler-Szilagyi-Pollney-Rezzolla 2008：harmonic 形式下的约束保持边界条件
  - Winicour 2009：Sommerfeld 型 Einstein 边界数据的几何表达

## 轮回 2 · 2026-04-29 10:06:12 CST · 物质特征边界条件已接入且未破坏 `B/C` 稳定性

- 已在
  `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_bc_from_a_initial_data.py`
  中移除物质变量 `n,S` 的乘法海绵层，改为特征入流/出流边界条件：
  - 左/右边界按 `v^x` 符号判断 inflow/outflow
  - 上/下边界按 `v^z` 符号判断 inflow/outflow
  - inflow 点直接施加同一高斯束准备下的解析参考 `rho_ref,S_ref`
  - outflow 点用相邻内点的单边外流更新
- 几何变量边界条件保持不变：
  - `h_xx -> 1, h_xz -> 0, h_zz -> 1, beta -> 0`
  - `k_xx, k_xz, k_zz, k_beta -> 0`
  - `Phi -> 1, Pi_Phi -> 0`
- 已完成两组 `200` 步、`dt=2.5e-4` 的长时间测试：
  1. `C` 平直辅助场（即严格退回 `B`）
  2. `C` 取 `r3_proxy` 初值且 `ell=0.2`
- 长时间结果：
  - `B` 支仍保持弱场稳定：
    - `max|h_xx-1| ≈ 7.74e-3`
    - `max|h_xz| ≈ 5.33e-6`
    - `max|beta| ≈ 2.77e-9`
  - Hamilton 约束残差从 `~2.59e-1` 缓慢下降到 `~2.12e-1`
  - 动量约束残差增长到 `~4.81e-2`，但未出现灾难性失控
  - `ell=0.2` 时 `C-B` 偏离与旧海绵层版本几乎一致：
    - `max|h_xx^C-h_xx^B| ≈ 2.03e-6`
    - `max|beta^C-beta^B| ≈ 1.03e-6`
    - `max|Phi-1| ≈ 1.74e-6`
    - `|ΔH| ≈ 1.50e-4`
    - `|ΔM| ≈ 3.15e-5`
- 当前判断：
  1. 物质边界改成特征条件后，`B/C` 的稳定性没有被破坏
  2. 先前观察到的 `C-B` 弱分离不是旧海绵层的伪影
  3. 下一步不应继续修补物质边界，而应转向为 `C` 支构造椭圆型 `Phi` 初值方程

## 轮回 2 · 2026-04-29 10:17:48 CST · 原始椭圆型 `Phi` 初值已验证为“方向正确但整幅过强”

- 已在
  `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_bc_from_a_initial_data.py`
  中加入新的 `c_phi_mode = elliptic`
- 该模式在固定初始几何与初始 `(\rho,S)` 下，求解初始切片上的静态椭圆方程：
  - 取 `Pi_Phi = 0`
  - 取 `∂_t Pi_Phi = 0`
  - 令 `Phi` 在边界回到 `1`
- 短时间探针（`steps=20, dt=2.5e-4, ell=0.2`）给出非常明确的失败信号：
  - 初始时刻 `max|Phi-1| ≈ 0.999999`
  - 也就是主支撑区内 `Phi` 被直接拉到接近 `0`
  - `C` 支初始 Hamilton 残差约 `2.07e9`
  - 初始动量残差约 `6.41e2`
  - 只推进极短时间，`R^(3)` 先冲到 `~1e7`，随后主变量迅速出现 `NaN`
- 当前解释：
  - 这不是“椭圆型思路错了”
  - 而是“完整椭圆型 `Phi` 解”与当前弱场同步规范近似不兼容
  - 根本原因是当前 `C` 支有效 Einstein 源中含 `1/Phi`，所以一旦 `Phi` 被拉到接近 `0`，几何源会被极度放大
- 当前结论：
  1. `r3_proxy` 太弱
  2. 原始完整椭圆型 `Phi` 初值又太强
  3. 下一步若继续沿这条线走，必须引入“椭圆型初值的幅度控制/连续延拓”或改用更温和的辅助变量

## 轮回 2 · 2026-04-29 10:34:40 CST · 椭圆型 `Phi` 初值的小参数延拓窗口已初步锁定

- 已在同一设置下
  `steps=20, dt=2.5e-4, ell=0.2`
  继续扫描
  `elliptic_lambda = 1e-4, 3e-4, 1e-3, 3e-3`
- 扫描结果：
  - `1e-4`：稳定，`max|Phi-1| ~ 1e-4`，`max|h_xx^C-h_xx^B| ~ 1.31e-6`，`Hamiltonian gap ~ 1.26e-1`
  - `3e-4`：稳定但约束代价明显上升，`max|h_xx^C-h_xx^B| ~ 3.88e-6`，`Hamiltonian gap ~ 5.19e-1`
  - `1e-3`：`C-B` 分离更明显，但 `Hamiltonian gap ~ 1.91`，已经偏大
  - `3e-3`：`Hamiltonian gap ~ 5.88`，当前同步规范弱场原型下已明显过强
- 当前阶段判断：
  - 椭圆型 `Phi` 初值的受控延拓确实存在一个小而可用的稳定窗口
  - 当前最合适的首选基准是 `elliptic_lambda = 1e-4`
  - `3e-4` 可保留为更激进的对照值
- 又补做了长时间检查（`steps=200, dt=2.5e-4, ell=0.2`）：
  - `elliptic_lambda = 1e-4`：
    - `max|h_xx^C-h_xx^B| ~ 1.16e-4`
    - `Hamiltonian gap ~ 1.27e-1`
    - `Momentum gap ~ 9.89e-3`
    - 说明该窗口不是“只能短时使用”的假窗口，长时间下 `C-B` 分离会自然积累到 `10^{-4}` 量级
  - `elliptic_lambda = 3e-4`：
    - `max|h_xx^C-h_xx^B| ~ 3.32e-4`
    - `Hamiltonian gap ~ 3.85e-1`
    - `Momentum gap ~ 8.35e-2`
    - 仍可推进，但更适合作为激进对照，不宜当默认基准

## 轮回 2 · 2026-04-29 13:xx CST · 固定平直时空下局域双高斯波包参考初值已建立

- 按用户最新要求，暂停当前 `B/C` 椭圆型 `Phi` 初值线，转而先建立“固定平直时空、真正空间局域的双高斯波包”参考演化
- 新建脚本：
  `kg_examples/simulate_flat_localized_crossing_packets.py`
- 当前初值参数：
  - `m=1`
  - `k0=6`
  - `sigma_parallel=1.6`
  - `sigma_perp=1.2`
  - `alpha=0.5`
  - `phi0=0`
  - 初始中心 `(-6,-6)` 与 `(6,-6)`
  - 两束分别沿 `+45^\circ` 与 `-45^\circ` 方向传播
  - 固定区域 `x,z ∈ [-20,20]`
- 方法：
  - 在平直背景上直接采用正频 Klein-Gordon 谱演化
  - 从复场重构 `rho=|psi|^2` 与 `S=arg psi`
- 当前结果：
  - 理论相遇时间 `t_meet ≈ 8.60`
  - 采样时刻：
    `t=0, 2.67, 5.33, 8.00, 10.67, 13.33, 16.00`
  - 边界相对密度最大值：
    `max_t rho_boundary_max/rho_max ≈ 1.30e-28`
  - 说明在整个 `t∈[0,16]` 时间窗内，边界已可视为真空平直背景
- 新增正式笔记：
  `research-notes/089-固定平直时空下局域双高斯波包在±45度方向相交的演化.md`
- 当前方法论更新：
  - 以后 `B/C` 的主要全动力学测试，应优先建立在这组真正局域的平直参考初值上，而不再使用旧的束状初值
- 已把 `simulate_bc_from_a_initial_data.py` 默认切换为 `initial_mode = localized`
- 新增局域波包边界条件实现：
  - 物质边界参考不再来自旧的束状解析解
  - 而是来自同一组局域双高斯波包的精确平直 Klein-Gordon 演化
- 新的局域波包 `B/C` 基线测试（`steps=40, dt=2.5e-4, c_phi_mode=flat`）结果：
  - `B` 支初始 Hamilton 残差 `~5.75e-4`
  - `B` 支初始动量残差 `~3.02e-6`
  - `40` 步后 `max|h_xx-1| ~ 1.08e-3`
  - 边界相对密度最大值 `~5.95e-10`
  - `C(flat)` 与 `B` 各分量差异严格为 `0`
- 新的局域波包 `C` 非平凡测试（`ell=0.2, elliptic_lambda=1e-4`）：
  - `40` 步后：
    - `max|Phi-1| ~ 9.62e-5`
    - `max|h_xx^C-h_xx^B| ~ 2.99e-6`
    - `ΔH ~ 3.17e-2`
    - `ΔM ~ 1.43e-3`
    - 边界相对密度最大值仍 `~5.95e-10`
  - `200` 步后：
    - `max|Phi-1| ~ 1.38e-5`
    - `max|h_xx^C-h_xx^B| ~ 6.84e-5`
    - `max|h_zz^C-h_zz^B| ~ 7.72e-5`
    - `max|beta^C-beta^B| ~ 3.82e-5`
    - `ΔH ~ 1.50e-2`
    - `ΔM ~ 5.60e-3`
    - 边界相对密度最大值仍仅 `~1.32e-9`
- 新增正式笔记：
  `research-notes/090-局域双高斯波包初值下BC两支的同步规范早期演化.md`

## 轮回 2 · 2026-04-29 14:xx CST · A/B 同初值比较与 B 支图示已补齐

- 为回应用户“先告诉我 AB 的差异有多大，B 的演化结果图示给我”，新增专门分析脚本：
  `kg_examples/analyze_ab_localized_reference_vs_b.py`
- 在实现该比较时，发现 `B` 支初值里守恒密度 `n` 的构造此前并未严格从同一组 `(rho,S,g)` 直接给出，而是通过了一个不一致的中间重构
- 已在
  `kg_examples/simulate_bc_from_a_initial_data.py`
  中修正为从同一组初始 `(\rho,S,g)` 直接定义
  `n = e^beta sqrt(h) rho E`
  以保证 `A/B` 真正从同一组初值启动
- 修正后，`A/B` 在 `t=0` 的密度差已经回到数值噪声量级：
  - `relative L1 ~ 3.54e-10`
  - `max abs diff ~ 1.00e-12`
- 在局域双高斯波包、`steps=200`、`dt=2.5e-4`、即 `t=0.05` 的早期时间窗内，`A/B` 的密度差为：
  - `final relative L1 ~ 1.57e-4`
  - `final max abs rho diff ~ 4.47e-5`
  - 在采样时刻上的峰值与最终值相同，说明当前短窗口内差异单调缓慢积累
- 同一时刻 `B` 支几何偏离与约束量级为：
  - `max|h_xx-1| ~ 1.08e-3`
  - `max|h_xz| ~ 1.56e-5`
  - `max|h_zz-1| ~ 1.08e-3`
  - `max|beta| ~ 6.93e-9`
  - `max Hamilton residual ~ 1.09e-4`
  - `max momentum residual ~ 1.87e-5`
- `B` 支图示与 `A/B` 差异图现已输出到：
  `visualizations/ab_localized_reference_vs_b/`
  其中包括：
  - `rho_b_montage.png`
  - `rho_b_log10_montage.png`
  - `rho_b_minus_a_montage.png`
  - `hxx_b_minus_1_montage.png`
  - `ab_difference_vs_time.png`

## 轮回 2 · 2026-04-29 14:xx CST · 当前脚本耗时基准

- 已实测固定平直参考脚本
  `kg_examples/simulate_flat_localized_crossing_packets.py`
  在当前机器上直接完成 `t∈[0,16]`、`7` 个采样帧的总耗时约为：
  - `real 3.58 s`
- 已实测当前 `B/C` 同步规范全动力学脚本
  `kg_examples/simulate_bc_from_a_initial_data.py --initial-mode localized --c-phi-mode flat --steps 200 --dt 2.5e-4`
  在当前机器上的总耗时约为：
  - `real 65.69 s`
- 由此估算：
  - 该 `B/C` 原型在 `dt=2.5e-4` 下若推进到 `t=16`，需要 `64000` 步
  - 若近似线性缩放，则总耗时约为
    `65.69 × (64000/200) ≈ 2.10e4 s ≈ 5.8 h`

## 轮回 2 · 2026-04-29 14:xx CST · 稀疏有限性检查与长时预检

- 为在不改方程的前提下降低长跑诊断开销，已把
  `kg_examples/simulate_bc_from_a_initial_data.py`
  中的有限性检查改成可配置频率：
  `finite_check_every`
  - `1` 表示每步都查
  - `100` 表示每 `100` 步查一次，并保留第 `1` 步检查
- 该修改不改变任何演化方程、边界条件或离散格式，只降低 `NaN/Inf` 诊断频率
- 实测当前默认非平凡 `C` 设置
  `ell=0.2, elliptic_lambda=1e-4`
  在
  `steps=200, dt=2.5e-4, finite_check_every=100`
  下耗时约：
  - `real 84.95 s`
  因而外推到
  `t in [0,16]`
  （即 `64000` 步）总耗时约：
  - `~7.6 h`
- 已完成两轮更长时间预检：
  1. `B(flat)` 到 `t=0.2`（`steps=800`）：
     - 未出现任何非有限值
     - `max|h_xx-1| ~ 1.08e-3`
     - `max|h_xz| ~ 6.45e-5`
     - `max Hamilton residual ~ 1.06e-4`
     - `max momentum residual ~ 2.15e-5`
     - 边界相对密度 `~4.68e-10`
  2. 非平凡 `C` 到 `t=0.2`（`ell=0.2, elliptic_lambda=1e-4`）：
     - 未出现任何非有限值
     - `max|Phi-1| ~ 3.88e-4`
     - `max|h_xx^C-h_xx^B| ~ 4.94e-4`
     - `max|rho_C-rho_B| ~ 2.03e-4`
     - `Hamiltonian gap ~ 3.57e-2`
     - `Momentum gap ~ 1.20e-2`
     - 边界相对密度 `~4.68e-10`
- 结论：当前局域双高斯波包初值、物质特征边界条件和默认非平凡 `C` 初值，在 `t=0.2` 的中等时间窗内没有发现提前崩溃或非有限值机制

## 轮回 2 · 2026-04-29 14:xx CST · `t∈[0,16]` 长跑已启动

- 已用持续 PTY 会话启动正式长时间计算：
  - 脚本：
    `kg_examples/simulate_bc_from_a_initial_data.py`
  - 参数：
    `--initial-mode localized`
    `--c-phi-mode elliptic`
    `--ell 0.2`
    `--elliptic-lambda 1e-4`
    `--steps 64000`
    `--dt 2.5e-4`
    `--finite-check-every 100`
  - 输出目录：
    `visualizations/bc_localized_elliptic1e4_t16_sparse100_live/`
- 由于该脚本只在运行结束时写出 `summary.json`，中途不会持续刷新文本日志；后续需通过会话轮询确认进度与结束状态

## 轮回 2 · 2026-04-29 15:xx CST · `t∈[0,16]` 长跑失败

- 对当前默认非平凡 `C` 设置
  `ell=0.2, elliptic_lambda=1e-4`
  已用 PTY 会话启动长跑：
  `steps=64000, dt=2.5e-4, finite_check_every=100`
- 长跑并未完成，而是在中途失败退出；失败信息显示首先出现的是几何和物质重构中的指数溢出：
  - `exp(beta)` 溢出
  - `exp(-2 beta)` 溢出
  - 随后 `rho` 重构与几何源项出现 `NaN/Inf`
  - 最终在有限性检查处触发
    `FloatingPointError: 步进前: 变量 h_xx 出现非有限值`
- 这说明：
  - 当前同步规范、当前 `Phi` 初值与当前边界条件，在 `t=0.2` 前的中等时间窗内仍可工作
  - 但把时间窗直接拉到 `t=16` 时，系统会进入明显非线性失稳区，当前版本不能直接完成全时长演化
- 当前最直接的后续修正方向应是：
  1. 先给长跑加步数/物理时间 checkpoint 日志，定位失稳发生的大致时间
  2. 再决定是缩短研究时间窗，还是修改规范/变量以压制 `beta` 驱动的指数爆炸

## 轮回 2 · 2026-04-29 15:xx CST · 溢出定位诊断脚本已建立并启动

- 新增诊断脚本：
  `kg_examples/diagnose_bc_overflow.py`
- 该脚本不改任何动力学方程，只做三件事：
  1. 每隔固定步数记录一次 `t`、`beta` 的最小/最大值以及 `h_xx,h_xz,h_zz` 的量级
  2. 开启 `numpy` 的 `over='raise', invalid='raise'`
  3. 在第一次指数溢出或 `NaN/Inf` 出现时，把当时的 `t`、异常类型、异常信息和当前几何量写入结果文件
- 诊断输出目录：
  `visualizations/bc_overflow_diagnosis/`
  其中：
  - `checkpoints.jsonl`：检查点序列
  - `result.json`：首次失败或完整完成的摘要
- 当前硬结论：
  - 已知系统在 `t=0.2` 之前稳定
  - 原先那次 `t∈[0,16]` 长跑的确发生了 `exp(beta)` / `exp(-2 beta)` 相关溢出
  - 但首次溢出的精确 `t` 与当时 `beta` 数值，需要等该专门诊断脚本跑到失败点后才能给出

## 轮回 2 · 2026-04-29 16:xx CST · `yy` 分量中可严格吸收的指数因子已改写

- 按用户要求，专门检查“不改变方法”的前提下，代码中哪些 `exp(±2 beta)` 只是实现层面被拆开了、其实应先代数约掉。
- 已确认并修正两类严格可吸收的组合：
  1. `exp(-2 beta) * R^(3)_{yy}`  
     其中 `R^(3)_{yy}` 原先在
     `adm_ykilling_geometry.py`
     里写成
     `-exp(2 beta) * (lap_beta + grad_beta_sq)`；
     现改为直接返回
     `r3_yy_reduced = -(lap_beta + grad_beta_sq)`，
     几何演化方程只使用该 reduced 量。
  2. `exp(-2 beta) * S_{yy}`  
     `A/B/C` 三支源项中原先都显式构造了带 `exp(2 beta)` 的 `s_yy`；
     现统一改为直接构造
     `s_yy_reduced := exp(-2 beta) * s_yy`，
     并在几何演化中只使用 reduced 量。
- 相关文件已更新：
  - `kg_examples/adm_ykilling_geometry.py`
  - `kg_examples/adm_sources_abc.py`
  - `kg_examples/adm_initial_data_abc.py`
  - `kg_examples/simulate_ab_adm_synchronous.py`
  - `kg_examples/simulate_bc_geometry_from_a_reference.py`
  - `kg_examples/simulate_bc_from_a_initial_data.py`
- 语法检查：
  `python3 -m py_compile`
  已通过。
- 改写后的两组 `200` 步短时复查也已通过：
  1. `B(flat)`：结果与改写前保持一致量级，说明这一步只是实现修正，不是改模型；
  2. 非平凡 `C(ell=0.2, elliptic_lambda=1e-4)`：早期 `C-B` 分离仍存在，未被此实现修正抹掉。
- 当前判断：
  - `yy` 分量中的一部分指数溢出风险确实是实现把可抵消因子拆开造成的；
  - 但物质守恒密度
    `n = exp(beta) sqrt(h) rho E`
    及从 `n` 恢复 `rho` 的那组 `exp(beta)` 不是同类，它们目前还不能在不改变变量定义的前提下直接约掉。

## 轮回 2 · 2026-04-29 19:xx CST · 改写后长时间溢出诊断的中期结果

- 已用新的独立输出目录重新启动诊断：
  `kg_examples/diagnose_bc_overflow.py --output visualizations/bc_overflow_diagnosis_postabsorb --checkpoint-every 200`
- 当前检查点已推进到：
  `t = 0.70`（第 `2800` 步）
- 到 `t=0.70` 为止，**仍未发生溢出或非有限值**。
- 当前最重要的量级：
  - `B(flat)`：
    - `beta_max ~ 1.24e-6`
    - `max|h_xx-1| ~ 1.11e-3`
    - `max|h_xz| ~ 2.47e-4`
  - 非平凡 `C`（`ell=0.2, elliptic_lambda=1e-4`）：
    - `beta_min ~ -6.56e-4`
    - `max|Phi-1| ~ 1.28e-3`
    - `max|h_xx-1| ~ 2.11e-3`
    - `max|h_xz| ~ 4.46e-4`
- 阶段性判断：
  - 把 `yy` 分量里严格可约掉的 `exp(±2 beta)` 先吸收之后，系统已经明显超过此前只确认到 `t=0.2` 的稳定窗口；
  - 当前 `beta` 量级仍远小于机器指数溢出阈值，说明早期到中期演化中，`yy` 分量的拆分实现确实是重要问题之一；
  - 后续若仍发生失稳，更可能来自：
    1. 物质守恒密度 `n` 与 `rho` 反演中的 `exp(beta)`；
    2. 或同步规范本身的长时间几何失稳。

## 轮回 2 · 2026-04-29 19:xx CST · 改写后诊断已推进到 `t=3.4`

- 新诊断输出目录：
  `visualizations/bc_overflow_diagnosis_postabsorb/`
- 到 `t=3.4`（第 `13600` 步）为止，仍未出现溢出或非有限值。
- 当前量级：
  - `B(flat)`：
    - `beta_max ~ 9.59e-6`
    - `max|h_xx-1| ~ 1.86e-3`
    - `max|h_xz| ~ 1.54e-3`
  - 非平凡 `C`（`ell=0.2, elliptic_lambda=1e-4`）：
    - `beta_min ~ -2.87e-3`
    - `max|h_xx-1| ~ 8.90e-3`
    - `max|h_xz| ~ 2.44e-3`
    - `max|Phi-1| ~ 5.74e-3`
- 阶段性结论进一步加强：
  - 先前的长时间失败不能再简单归因于 `yy` 分量里那批本可先约掉的指数因子；
  - 改写后，系统已经稳定超过原先只知会在 `t>0.2` 失稳的窗口很多；
- 真正的剩余风险更像是晚时刻的 `n/rho` 指数反演或同步规范长期失稳。

## 轮回 2 · 2026-04-29 20:xx CST · 改写后长跑首次失败点已定位

- 改写后诊断长跑在
  `visualizations/bc_overflow_diagnosis_postabsorb/`
  中已产生正式结果文件：
  `result.json`
- 首次失败发生在：
  - `step = 31198`
  - `t = 7.79925`
  - 下一个子步时间：
    `t_next = 7.7995`
  - 异常：
    `FloatingPointError: overflow encountered in exp`
- 失败时当前已接受状态的量级：
  - `B(flat)`：
    - `beta_max ~ 3.43e-5`
    - `max|h_xx-1| ~ 6.62e-1`
    - `max|h_xz| ~ 7.08`
    - `max|h_zz-1| ~ 4.90e1`
  - 非平凡 `C`：
    - `beta_min ~ -4.21e-3`
    - `max|Phi-1| ~ 8.47e-3`
    - 几何量级与 `B` 相同
- 重要判断：
  - 当前失败点处，**已接受状态中的 `beta` 仍然很小**，远未达到 `exp(beta)` 或 `exp(-2 beta)` 的机器溢出阈值；
  - 因此这次 `overflow encountered in exp` 更可能来自：
    1. `RK4` 的某个中间试探态（trial state）中 `beta` 突然暴涨；
    2. 或长时间几何失稳导致某个中间态先跨出正常范围，而不是最终已接受状态本身的 `beta` 直接指数爆炸；
  - 换句话说：改写 `yy` 分量里可严格吸收的指数因子以后，先前那种“显然由已接受态中的 `e^beta` 拆分实现导致”的溢出，已经不再是主问题。

## 轮回 2 · 2026-04-29 22:xx CST · `1+log` 切片长跑结果

- 已完成零 shift 的 `1+log` 长时间溢出诊断，输出目录：
  `visualizations/bc_overflow_diagnosis_1plog/`
- 该长跑未能到达 `t=16`；首次失败发生在：
  - `step = 11645`
  - `t = 2.911`
  - `t_next = 2.91125`
  - 异常：
    `FloatingPointError: overflow encountered in exp`
- 失败时已接受状态量级：
  - `B(flat)`：
    - `beta_max ~ 8.46e-6`
    - `lapse_min ~ 0.999532`
    - `max|h_xx-1| ~ 1.52e-3`
    - `max|h_xz| ~ 9.06e-4`
  - 非平凡 `C`：
    - `beta_min ~ -2.46e-3`
    - `lapse_min ~ 0.270`
    - `lapse_max ~ 9.21e2`
    - `max|h_xx-1| ~ 1.81e3`
    - `max|h_xz| ~ 1.01e3`
    - `max|h_zz-1| ~ 2.36e4`
    - `max|Phi-1| ~ 5.01e-3`
- 当前判断：
  - `1+log` 对 `B(flat)` 明显更稳，但没有解决当前非平凡 `C` 支的长期病态；
  - 本次失败的主危险量已从 `beta` 转移到 `lapse` 与空间度规分量本身；
  - 在“零 shift + 1+log”下，`C` 支几何先进入强非线性，再触发子步中的 `exp(...)` 溢出。

## 轮回 2 · 2026-04-30 00:xx CST · `B` 支单独长时间诊断中期结果

- 已将 `B` 支从 `C` 中单独抽出，使用：
  `kg_examples/diagnose_b_overflow.py --gauge-mode 1plog --steps 64000 --dt 2.5e-4`
- 当前输出目录：
  `visualizations/b_overflow_diagnosis_1plog_t16/`
- 到 `t=3.0`（`step=12000`）为止，`B` 支仍未出现溢出或非有限值。
- 当前量级：
  - `beta_max ~ 8.68e-6`
  - `lapse_min ~ 0.999531`
  - `lapse_max ~ 1.000232`
  - `max|h_xx-1| ~ 1.54e-3`
  - `max|h_xz| ~ 9.26e-4`
  - `max|h_zz-1| ~ 1.54e-3`
- 当前判断：
  - 这进一步支持“当前主要长期稳定性难点集中在非平凡 `C` 支，而不是 `B` 支”
  - 若 `B` 最终可长跑到 `t=16`，则后续数值方法升级应优先围绕 `C` 支的几何-辅助场耦合稳定性来做

## 轮回 2 · 2026-04-30 00:xx CST · `B` 支单独长跑继续推进到 `t=6.0`

- `B` 支单独 `1+log` 长跑仍在继续，当前诊断目录：
  `visualizations/b_overflow_diagnosis_1plog_t16/`
- 到 `t=6.0`（`step=24000`）为止，仍未出现溢出或非有限值。
- 当前量级：
  - `beta_max ~ 1.55e-5`
  - `lapse_min ~ 0.98795`
  - `lapse_max ~ 1.01173`
  - `max|h_xx-1| ~ 9.23e-2`
  - `max|h_xz| ~ 1.01e-2`
  - `max|h_zz-1| ~ 9.22e-2`
- 当前判断：
  - `B` 支虽然到 `t≈5.0` 后几何偏离开始明显加快，但到 `t=6.0` 仍未进入类似非平凡 `C` 的病态区；
  - 这进一步支持当前方法对 `B` 支比对非平凡 `C` 支稳得多。

## 轮回 2 · 2026-04-30 00:xx CST · `B` 支单独 `1+log` 长跑首次失败点

- `B` 支单独长跑目录：
  `visualizations/b_overflow_diagnosis_1plog_t16/`
- 首次失败发生在：
  - `step = 27380`
  - `t = 6.84475`
  - `t_next = 6.845`
  - 异常：
    `FloatingPointError: overflow encountered in exp`
- 失败前的最后稳定检查点在：
  - `step = 27000`
  - `t = 6.75`
  - `beta_max ~ 2.16e-5`
  - `lapse_min ~ 0.539`
  - `lapse_max ~ 1.375`
  - `max|h_xx-1| ~ 1.94`
  - `max|h_xz| ~ 2.74e-1`
  - `max|h_zz-1| ~ 1.94`
- 首次失败时已接受状态已经明显病态：
  - `beta_max ~ 7.02e3`
  - `lapse_min ~ -9.21e7`
  - `lapse_max ~ 2.68e5`
  - `max|h_xx| ~ 1.73e22`
  - `max|h_xz| ~ 2.09e22`
  - `max|h_zz| ~ 1.34e24`
- 当前判断：
  - 当前“零 shift + 1+log”并非只对 `C` 支失稳；`B` 支在更晚时刻也会进入强非线性并失稳
  - 但 `C` 支仍然更早失稳：`t_C ~ 2.911 < t_B ~ 6.84475`
  - 因此当前问题不是 `C` 支特例，而是现有几何演化形式在强场长期区总体不够稳；`C` 只是在更早暴露这个问题

## 轮回 2 · 2026-04-30 00:xx CST · `B` 支 `BSSN` 型强场原型

- 新增文件：
  - `kg_examples/simulate_b_bssn_driver.py`
  - `kg_examples/diagnose_b_bssn_overflow.py`
- 新原型采用：
  - 共形-无迹变量：`chi, gt_xx, gt_xz, gt_zz, k_trace, at_xx, at_xz, at_zz`
  - 规范：`1+log` 切片
  - 非零 `shift`：一阶 `Gamma` 驱动型条件
- 已完成 `t=0.2` 短时间回归：
  - 输出目录：`visualizations/b_bssn_driver_t02/`
  - 终态量级：
    - `chi_min ~ 0.999289`
    - `lapse_min ~ 0.999988`
    - `shift_abs_max ~ 1.16e-5`
    - `max|h_xx-1| ~ 1.08e-3`
    - `max|h_xz| ~ 6.49e-5`
    - `max|h_zz-1| ~ 1.08e-3`
  - 与旧 `B` 支弱场窗口结果一致，说明新变量与新规范在早期没有引入明显物理偏差
- 已启动带检查点的长时间诊断：
  - 目录：`visualizations/b_bssn_overflow_diagnosis_t16/`
  - 目标：检查是否能跨过旧方法的 `t_B ~ 6.84475` 失稳时刻

## 轮回 3 · 2026-04-30 03:xx CST · 直接由变换生成初始 `\tilde g_0` 的一致性检查

- 新增模块：
  - `kg_examples/mixed_tilde_initial_data.py`
- 该模块做的事：
  1. 从局域双高斯波包的平直背景正频解，直接计算 `\rho,S` 的精确导数；
  2. 计算 `Q,X,Y,Z,\Delta`；
  3. 不再用裸的 `Q/(X\Delta) M^{\mu\nu}`，而改写成
     `\tilde g^{\mu\nu}=g^{\mu\nu}-(Q/X)P^{\mu\nu}`，
     其中 `P^{\mu\nu}` 通过二维 Gram 矩阵的 Moore-Penrose 伪逆构造；
  4. 尝试直接从该 `\tilde g^{\mu\nu}` 提取 `lapse/shift/空间度规` 作为 `B/C` 的严格初始几何。
- 结果：
  - `\Delta \to 0` 不是当前初始数据的主问题；伪逆投影形式可以把这类可去奇性吸收掉。
  - 真正的新问题是：在当前局域高斯波包初值下，主支撑区里 `\tilde g^{tt}` 会变成非正。
  - 单包和双包都存在这一现象，不只是双包干涉造成的。
- 具体数值（主支撑区定义为 `\rho > 10^{-3}\rho_{\max}`）：
  - 单包：
    - `X_min ~ -2.8666`
    - `min \tilde g^{tt} ~ -8.1268`
    - `730 / 1714` 个主支撑区网格点满足 `\tilde g^{tt} <= 0`
  - 双包：
    - `X_min ~ -2.8666`
    - `min \tilde g^{tt} ~ -8.1268`
    - `1424 / 3246` 个主支撑区网格点满足 `\tilde g^{tt} <= 0`
- 这说明：
  - 之前的 `B/C` 初值结果之所以能跑，是因为当时并没有真正用变换生成初始 `\tilde g_0`。
  - 一旦按用户要求严格生成 `\tilde g_0`，当前问题就从“分母奇性”转移成了“初始 `\tilde g` 是否仍然是可接受的 Lorentz 度规/同一分支”。
- 2026-04-30: 在用户明确接受 `\tilde g` 可具有实验室坐标下的奇异/时间反向现象后，开始转向“固定实验室坐标下直接推进 `\tilde g_{\mu\nu}` 分量”的路线，不再把 `ADM` 友好性当作理论前提。
- 2026-04-30: 新增 `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/coordinate_matter_evolution.py`，把物质方程改写成固定坐标下的守恒密度 `N=\sqrt{|g|}J^t` 与相位 `S` 的分量形式，完全不依赖 `ADM` 变量。
- 2026-04-30: 在平直 `g` 下，直接用 `g` 回归会失败，因为精确 KG 相位满足的是带 `Q` 的 Bohm-HJ 方程，不是无修正质量壳；因此正确回归对象必须是由变换生成的 `\tilde g`。
- 2026-04-30: 对当前局域双高斯波包，用严格变换生成的 `\tilde g_0` 再把 HJ 方程写成关于实验室时间导数 `S_t` 的二次方程时，判别式在大量网格点为负。这表明问题已不只是 `ADM` 切片不友好，而是“固定实验室时间下的实相位演化”本身也需要重新处理。
- 2026-04-30：新增 `mixed_tilde_initial_data.py` 中的 `localized_direct_tilde_coordinate_initial(...)`，可直接从 `g_0=\eta,\rho_0,S_0` 严格生成坐标分量形式的初始 `\tilde g_0` 与 `\partial_t \tilde g_0`，不再经由 `ADM` 的 lapse/shift 提取。
- 2026-04-30：重新核对局域双高斯波包的严格初值后，主支撑区内 `\tilde g^{ab}u_a u_b=m^2` 误差约 `2e-9`，说明当前初值与变换在物理相关区域是自洽的。
- 2026-04-30：新增 `coordinate_matter_evolution.py` 中基于协变量 `(u_x,u_z)` 的固定实验室坐标守恒流写法，避免再用 `unwrap(arg psi)` 的相位梯度重构作为数值主变量。
- 2026-04-30：新增 `simulate_b_coordinate_frozen_metric.py` 作为最小原型。直接在整个网格上反演密度会因稀薄尾部 `J^t\to 0` 病态而失败；按用户认可的“零测度密度区不要求速度定义”解释，把尾部真空化后，冻结 `\tilde g_0` 的短时间原型可稳定推进到 `t=0.02`。
- 2026-04-30：`coordinate_matter_evolution.py` 已扩展为支持局域质量平方场；后续 `B` 支固定实验室坐标直接演化时，可直接把 Hamilton-Jacobi 质量壳从常数 `m^2` 替换成 `m^2 e^{-2\beta}`。
- 2026-04-30：开始实现 `simulate_b_direct_harmonic_minimal.py`，尝试在固定实验室坐标下直接推进 `\tilde g_{ab}`（最小广义调和型原型），不再经由 `ADM/BSSN` 切片。
- 2026-04-30：严格初值模块的首个实现问题已经排除：`3×3` 逆度规块在极低密度尾部会接近奇异，导致直接矩阵求逆失败；现已在 `mixed_tilde_initial_data.py` 与 `coordinate_matter_evolution.py` 中统一改用 Moore-Penrose 伪逆，并显式记录尾部/近奇异掩膜。
- 2026-04-30：新的最小直接分量原型仍然会在极短时间窗内迅速失稳，但定位结果表明：当前主问题不再是低密度尾部奇异，也不是局部 `6×6` 线性系统病态；在主支撑区内，这个 `6×6` 系统条件数接近 `1`，而爆炸源于这版固定坐标几何方程本身要求过大的 `g_{ab,tt}`。
- 2026-04-30：将活跃区缩到高密度核心（例如 `rho > 0.1 rho_max`）并把活跃区外的 `metric_t,beta_t,u_i` 全部真空化后，原型仍会在 `t=0.002` 内把 `metric_abs_max` 推到 `10^8` 量级，说明当前“逐点解 `g_{ab,tt}`”的最小广义调和闭合还过于生硬，不能直接作为正式强场积分器。
- 2026-04-30：已在 `simulate_b_direct_harmonic_minimal.py` 中真正加入固定实验室坐标下的约束阻尼项，并完成 `kappa=1,10` 的短时间对照测试。阻尼对 `metric_t` 与 `beta_t` 的爆炸有明显抑制，但不足以从根本上稳定当前最小原型；问题已不只是约束模增长，还包括这版逐点 `g_{ab,tt}` 闭合本身过于刚性。
- 2026-04-30：已直接升级到一阶广义调和约束阻尼原型 `simulate_b_direct_gh_damped_full.py`，把 `(g_{ab},Pi_{ab},Phi_{i ab})` 与 `(beta,Pi_beta,Phi_{i beta})` 作为完整状态并采用 `RK4` 推进。
- 2026-04-30：新的完整原型在 `kappa=10`、`dt=2.5e-4` 的短时间窗内仍会因 `exp(-beta)` 与几何失真迅速病态化；进一步把步长缩小到 `dt=6.25e-5` 依然失败，说明当前问题不是单纯时间步过大，而是固定坐标强场闭合与物质质量项耦合的结构性不稳。
- 2026-04-30：阶段性评估完成。当前数值困难不再主要来自实现细节，而是来自理论闭合仍不够明确：尤其是 `C` 支辅助场初值、固定实验室坐标下应使用的基本变量体系、以及是否只在主支撑区求解/如何与外部真空拼接，这些都还需要先做理论整理。
- 2026-04-30：在当前 `m\neq0, C=1` 混合投影分支下，已重新整理精确比值关系：逆度规行列式比值 `det(\tilde g^{..})/det(g^{..})=(m^2/X)^2`，协变体积元比值 `sqrt(|\tilde g|)/sqrt(|g|)=|X|/m^2`，测度密度比值 `sqrt(|\tilde g|)\tilde\rho / (sqrt(|g|)\rho)=|X|/m^2 -> 0` 当 `X->0`，而流 `\tilde J^\mu = sqrt(|\tilde g|)\tilde\rho\,\tilde g^{\mu\nu}u_\nu` 与原流 `J^\mu=sqrt(|g|)\rho u^\mu` 精确相等。结论是：`X->0` 时真正稳定的对象是流与测度密度，而不是裸的“速度样”量 `\tilde g^{\mu\nu}u_\nu`。
- 2026-04-30：用户明确要求：凡遇到理论定义、初值构造、规范选择、近似层级等存在不确定性的地方，必须先向用户确认，不能再由代理自行作主张引入近似或替代定义。此前“用 A 支弱场背景近似代替严格 `\tilde g_0`”已被用户判定为走弯路，后续必须避免重犯。
- 2026-04-30：用户明确要求持续维护 `agent-memory`，并在关键理论分歧、方法路线变更、阶段性结论出现后及时回写 `LOG.md` 与 `DECISIONS.md`。后续不能只在数值推进后更新，理论口径修正也必须同步记录。
- 2026-04-30：已完成 `3+1d` 下 A/B/C 三支引力作用量密度的渐进一致性第一轮分析。结论：`B≈A` 需要 `|Q|/m^2 << 1` 且 `Q/m^2` 与 Bohm 平面投影子 `Π` 缓变；`C≈B` 还需要 `|ℓ^2 R~| << 1` 且 `f_R` 的导数项小。只有这两层同时成立时，才可期待 `A≈B≈C`。`|Q|>>m^2` 或 `X->0` 一般不是三支渐进一致区。
- 2026-04-30：用户补充要求分析两条额外渐近线：1）`g->η` 时，对一般 `Q,X` 的 `A/B/C` 作用量密度结构；2）`m=0` 时三支的相应极限。当前结论：在 `g=η` 时，`A` 支引力密度直接为零，而 `B/C` 支完全由变换诱导的 `\tilde g` 曲率控制；在有效 `1+1d` 时有精确公式 `sqrt(|g~|)R~ = -Box_η ln(X/m^2)`（massive）或其 `Q/κ` 版本（massless patch）。`m=0` 时不能把 `m≠0` 的 `X/m^2` 公式平滑取极限，必须使用单独的 `κ≠0` patch，且 `Q->0` 一般是退化/奇性边界而非平滑弱场极限。
## 2026-04-30 15:05 +0800
- User proposed a new qualitative strategy: A branch should stay near flat space because backreaction is suppressed by large Planck mass; B branch likely will not reproduce A full dynamics when `\tilde R` fluctuates strongly; C branch may suppress large-`|\tilde R|` regions if `\ell^2` is chosen large, while keeping near-Einstein behavior for small `|\tilde R|`.
- Need to analyze this carefully in terms of action densities and field equations, not only heuristic curvature size. Important correction: large `\ell^2` narrows the Einstein-like linear window `|\ell^2 \tilde R| \ll 1`, so matching ordinary weak-field experiments additionally requires a hierarchy between ordinary curvatures and interference-region curvatures.

## 2026-04-30 15:12 +0800
- User further clarified motivation for taking very large `\ell^2`: in A-branch simulations the quantum potential spans a huge range, so the induced `\tilde R` is expected to fluctuate strongly; to keep C-branch dynamics close to pre-transformation flat-space quantum mechanics, user wants most regions to lie in the saturated regime where `F(\tilde R) -> 0` (equivalently `f_R -> 0`) so gravitational response is strongly suppressed.

## 2026-04-30 15:20 +0800
- Need to explain carefully that in branch C with very large `\ell^2`, the saturated regime `|\ell^2 \tilde R| \gg 1` is not an Einstein limit but a degenerate metric-f(R) regime: the `f_R \tilde R_{\mu\nu}` piece is suppressed, while a residual cosmological-constant-like term `-(1/2)f \tilde g_{\mu\nu} ~ -sgn(\tilde R)\tilde g_{\mu\nu}/(2\ell^2)` remains, plus derivative terms from `\nabla\nabla f_R`.

## 2026-04-30 15:27 +0800
- New theoretical direction from user: analyze branch C using a thin-layer / delta-potential analogy. Bulk saturated regions may be nearly decoupled, while narrow transition layers where `\tilde R` changes rapidly should be treated by integrating the field equations to obtain jump / boundary conditions. Key object for junction analysis is the trace equation `f_R \tilde R - 2f + 3\tilde\Box f_R = \tilde T/M_P^2`.

## 2026-04-30 16:03 +0800
- Ran a first coarse diagnostic on the A-branch localized crossing reference, reconstructing `\tilde g` and evaluating C-branch trace-equation terms for `\ell=100` at `t=0,8,16` on a `96x96` grid.
- Result: `f_R \tilde R` is strongly suppressed and `-2f` stays near the saturated constant `2/\ell^2 = 2e-4`, but `3\tilde\Box f_R` develops very large spikes. This supports a "bulk saturation + sharp transition-layer" picture rather than smooth suppression everywhere.

## 2026-04-30 16:10 +0800
- Localized the largest `3\tilde\Box f_R` regions on the same coarse A-reference data. Conclusion: they are not all concentrated in destructive-interference nodes.
  - `t=0`: packets are still separated, yet strong spikes already exist.
  - `t=8`: the support/core maximum sits near the overlap region, but at relative density ~0.2, not at a deep density node.
  - `t=16`: the support/core maximum lies inside a dense packet region (`rho/rho_max ~ 0.55`).
- Interpretation: the large derivative term tracks transition layers of `\tilde R` / `f_R`, not simply low-density interference cancellation.

## 2026-04-30 16:24 +0800
- 用户明确指出“核心区”一词未先定义，造成阅读障碍。已据此修正两处：
  - 在 `agent-memory/README.md` 中加入强规则：任何自创工作名词、窗口名、区域名、掩膜名都必须在第一次出现处立即给出数学定义或数值阈值。
  - 在 `research-notes/094-A支参考解上C支场方程各项的粗网格诊断.md` 中补充定义：
    `主支撑区 := {rho > 1e-3 rho_max}`，
    `核心区 := {rho > 1e-1 rho_max}`。

## 2026-04-30 16:31 +0800
- 明确了一个关于 `C` 支 metric `f(R)` 方程的结构性结论：若要求导数项
  `(\nabla_\mu\nabla_\nu - g_{\mu\nu}\Box) f_R`
  对任意场构型都整体为零，则在标准 metric 变分下唯一可能是
  `f_{RR}=0`，即
  `f(R)=aR+b`。
- 因而任何真正非线性的饱和型 `f(R)`（例如当前 `R/sqrt(1+\ell^4 R^2)`）都不可能在一般变曲率背景中全局消掉导数项；它至多只能在某些特殊解上（例如 `R=const`）或在体区近似中把导数项压小。

## 2026-04-30 16:39 +0800
- 用户提醒“不必严格为零，只要有 `1/\ell^n` 压制即可”。据此补充更细判断：
  - 对标准饱和族 `f_\ell(R)=\ell^{-2}\Phi(\ell^2 R)`，在饱和体区 `|\ell^2 R|\gg1` 中，导数项确实可按逆幂被压小；
  - 但在过渡层 `|\ell^2 R|=O(1)` 中，链式法则会给出 `f_{RR}=O(\ell^2)`、`f_{RRR}=O(\ell^4)` 的系数，因此不存在“对所有区域统一按 `1/\ell^n` 压小”的非线性饱和型 metric `f(R)`。
- 这说明当前更准确的图景是：bulk suppression 可以成立，但 transition-layer suppression 不能靠 `f(R)` 选型单独保证。
- 另补充：若改用 Palatini 变分（度规与联络独立），metric 方程里不会出现 `(\nabla_\mu\nabla_\nu-g_{\mu\nu}\Box)f_R` 这一类导数项；这提供了一条结构性去掉高阶导数的替代路线，但它已经是不同于当前 metric `C` 支的理论。

## 2026-04-30 16:52 +0800
- 已完成球对称径向约化下的联络差与原 Einstein-Hilbert 项重写。对
  `ds^2 = h_ab dx^a dx^b + R(x)^2 dΩ_2^2`
  且 `\tilde h_ab = (X/m^2) h_ab`、角向部分不变的扇区，
  联络差只剩两类非零分量：
  1. 基底二维共形联络差 `ΔΓ^a_{bc}`；
  2. 角向抬升项 `ΔΓ^a_{ij} = (1 - m^2/X) (D^a R / R) g_{ij}`。
- 由此得到封闭公式：
  `sqrt(-g) R[g] = sqrt(-g~) [ R[g~] + □_{h~} ln(X/m^2) + 2/R^2 (m^2/X - 1) ]`
  其中 `□_{h~}` 是新径向二维基底度规 `\tilde h_ab` 的 d'Alembertian。
- 这说明在球对称径向问题里，A/B 引力项之差可被压缩成“二维共形 Laplacian 项 + 角向曲率半径修正项”。

## 2026-05-01 00:08 +0800
- 进一步澄清了“1+1d 里 A/B 看起来相等”的真正含义：
  - 对纯内禀 `1+1d` Einstein-Hilbert 项，`sqrt(-g)R` 在共形变换下只差总散度，因此 A/B 作为纯二维 EH 引力项只在边界项意义下不同；这与二维 EH 的拓扑/平凡动力学性质一致。
  - 但对 `3+1d` 球对称径向约化，A 项重写后会多出角向曲率半径修正 `2/R^2 (m^2/X-1)`，因此不再与 `sqrt(-g~)R[g~]` 仅差一个二维散度项。
- 同时区分了“本征表象中的基本作用量”与“变换后重写的作用量”：
  - 在 `3+1d` 里，A 支若以 `g` 为独立变量、B 支若以 `g~` 为独立变量，各自的纯引力部分都给 Einstein 型场方程；
  - 但把 A 重写到 `g~` 表象后，它不再是纯 Einstein-Hilbert 形式，因而不能再说“它对 `g~` 变分仍只是 Einstein 方程”。

## 2026-05-01 00:17 +0800
- 回答用户“把 A/B 两个 Einstein 方程换到同一表象下是否等价”时，结论已收紧为：
  - 纯内禀 `1+1d`：由于 EH 只差边界项且二维 Einstein 张量恒零，A/B 在引力方程层面只是一种“平凡等价”；
  - `3+1d` 球对称径向约化：A 重写到 `g~` 表象后等于 `sqrt(-g~)[R[g~] + □_{h~}\ln(X/m^2) + 2/R^2(m^2/X-1)]`，其中除边界项外还保留角向曲率修正，因此 A 的推前方程 = B 的 Einstein 方程 + 额外张量项，二者一般不等价，只在 `X/m^2 -> 1` 或 `R -> ∞` 等极限下渐进一致。
- 2026-05-01 关键澄清：当用户问“A支和B支换到同一表象下是否等价”时，正确比较对象是“把各自以自由场 `g` 与 `\tilde g` 变分得到的场方程推到同一变量表象后，是否给出同一个解集”；不能只比较作用量表面形式。当前结论：一般 `3+1d` 下不等价；只有纯内禀 `1+1d` 的特殊情形下，引力作用量才只差边界项。
- 2026-05-01 再次修正：上条若只基于引力项仍然不充分。判断 A/B 在同一表象下是否等价，必须比较**完整作用量**推前后的场方程；物质项不同意味着对应应力张量也不同，可能抵消引力侧额外项。正确判据是在同一表象下检查 `\Delta T_{\mu\nu}=M_P^2 H_{\mu\nu}` 是否成立，其中 `H_{\mu\nu}` 是 A 支引力项推前后相对 B 支 Einstein 项的额外几何项。
- 2026-05-01 再次纠正比较规则：用户要求比较的是“先分别以各自自由场 `g` 与 `\tilde g` 变分得到本征 Einstein 方程，再把 B 支本征方程中的变量代换回 `g` 表象，与 A 支本征 `g` 表象 Einstein 方程比照”；不是把 B 支作用量拉回 `g` 后再对 `g` 变分。今后该问题只按“先变分、后代换”的规则处理。
- 2026-05-01 已按“先变分、后代换”规则算出 B 支本征 Einstein 方程右边代回 `g` 表象后的精确形式，并与 A 支右边做差：在 `X>0` branch 上，
  `T~^(B→g)_{μν} = 2ρ u_μ u_ν - (g_{μν} + (Q/m^2)Π_{μν})ρ(X-m^2-Q)`；
  在共同 HJ 壳 `X-m^2-Q=0` 上，右边差压缩成
  `ΔT_{μν}^{(B-A)} = -2 r_μ r_ν + g_{μν}(r^2 + ρQ) = -2∂_μ√ρ∂_ν√ρ + g_{μν}[(∂√ρ)^2 + √ρ□_g√ρ]`。
- 2026-05-01 已算出 B 支本征 Einstein 张量代回 `g` 表象后的左边精确结构：`ΔG_{μν}^{(B-A)} = G~_{μν}[g~(g,ρ,S)]-G_{μν}[g]` 精确由联络差 `ΔΓ` 的 `K_{μν}` 组合给出；在弱形变 `h_{μν}=(Q/m^2)Π_{μν}` 的一般 `3+1d` 一阶展开中，`ΔG_{μν}` 至少含 `∂∂Q`，即一般含 `√ρ` 的四阶导数，而右边差在共同 HJ 壳上只到二阶。当前判断：一般 `3+1d` 下 A/B 两条本征 Einstein 方程代到同一 `g` 表象后不应期待严格等价，只可能在 `|Q|/m^2` 小且缓变时渐近一致。
- 2026-05-01 补充平直背景判断：设 `g=η`，则背景曲率项 `R[g], R_{μν}[g]` 确实消失，但 `ΔG_{μν}^{(B-A)}` 仍保留由 `ΔΓ[Q,Π]` 诱导的项，不会自动为零。最简反例：取固定 `t-x` Bohm 平面 `Π_{μν}=diag(1,-1,0,0)`，令 `q(x)=Q(x)/m^2`，则一阶有 `δG_{11}=-2 q''(x)`；故只要 `q''≠0`，即使背景完全平直，左边差也不消失。
- 2026-05-01 新澄清：若把引力度规也做场依赖变换，则“物质项混入引力作用量、引力项混入物质作用量”只是表象分块变化，不自动破坏完整变分等价。真正决定是否出现“测地线方程”的，不是某一项被归到引力部分还是物质部分，而是**哪一个度规上物质保持最小耦合**。若某表象中物质对该度规最小耦合，则测试粒子/几何光学极限沿该度规测地线；若换到另一表象后出现非最小耦合，则同一物理运动在该表象里会写成“测地线 + 额外力”形式，而不是该表象度规的纯测地线。
- 2026-05-01 进一步澄清：若把 A 支完整作用量拉回 `\tilde g` 表象并把 `(\tilde g,\rho,S)` 当独立变量，则由于 `g(\tilde g,\rho,S)` 依赖 `u_\mu=\partial_\mu S` 与 `Q(\rho)`，EH 项会显式贡献到 `S` 方程，因此该表象下的 `S` 方程一般不再是“纯最小耦合于 `\tilde g` 的测地线/HJ 方程”。但这些额外项不是任意的；按泛函链式法则，它们本质上是 A 支本征 Einstein 方程通过 `g(\tilde g,\rho,S)` 的场依赖混入 `S` 方程。故：作为独立基本理论，拉回后的 A 支不是“B 支那种最小耦合到 `\tilde g` 的理论”；作为 A 支的重写，它仍与 A 支完整方程组等价。
- 2026-05-01 新关键观察：若暂时不局限于当前混合变换，而考虑四维纯共形变换 `\tilde g_{\mu\nu}=\Omega^2 g_{\mu\nu}`，则
  `sqrt(-\tilde g)\tilde R = sqrt(-g)\Omega^2 [R - 6\Box\ln\Omega - 6(\nabla\ln\Omega)^2]
  = sqrt(-g)[\Omega^2 R - 6\Omega\Box\Omega]`.
  因而取 `\Omega=\sqrt{\rho}` 时，
  `sqrt(-\tilde g)\tilde R = sqrt(-g)[\rho R - 6\rho Q]`，其中 `Q=\Box\sqrt{\rho}/\sqrt{\rho}`。
  这说明：通过合适的新变换，物质作用量中的 `\rho Q` 结构确实可以在引力作用量里被精确几何化；但代价是引力项同时变成 `\rho R` 型非最小耦合，而不是仍保持纯 Einstein-Hilbert。
- 2026-05-01 若用户坚持使用四维共形变换 `\tilde g_{\mu\nu}=(\rho/\rho_*)g_{\mu\nu}`（即 `\Omega=\sqrt{\rho/\rho_*}`，无量纲且与 `\sqrt{\rho}` 成正比），并把新标量定义为 `\tilde\rho:=\rho_*/\rho`、`\tilde S:=S`，则原 A 支完整拉氏量
  `sqrt(-g)[(M_P^2/2)R + \rho(X-m^2-Q)]`
  精确重写成
  `sqrt(-\tilde g)[ \rho_* \tilde X - \rho_* m^2 \tilde\rho + (M_P^2/2)\tilde\rho\,\tilde R + (\rho_* - 3 M_P^2 \tilde\rho)\tilde Q ]`,
  其中 `\tilde X=\tilde g^{\mu\nu}\partial_\mu \tilde S \partial_\nu \tilde S`，`\tilde Q=\tilde\Box\sqrt{\tilde\rho}/\sqrt{\tilde\rho}`。
  结论：该共形因子把动能项常数化得很漂亮，但不存在一个**常数**倍 `\Omega=k\sqrt{\rho}` 能把 `Q` 在全拉氏量里全局消去；因为 `\tilde Q` 的系数是 `\rho_* - 3 M_P^2 \tilde\rho`，依赖场值而不是常数。
- 2026-05-01 进一步澄清：若某个“完全由 `\tilde g` 衍生的引力作用量”还额外显含 `\rho`，即 `S_{\rm grav}=S_{\rm grav}[\tilde g,\rho]` 但不含 `S/u`，则 `S` 变分（连续性方程）仍不受引力项直接影响；然而 `\rho` 变分会得到额外项，从而一般改写质量壳/HJ 方程 `\tilde g^{\mu\nu}u_\mu u_\nu-m^2=0`。因此若目标是保留纯 `\tilde g`-测地线/纯 HJ 结构，引力作用量原则上应对 `\rho` 也独立，除非其 `\rho` 依赖只贡献边界项或一个可吸收到质量中的常数。
- 2026-05-01 用户当前目标进一步收缩为：寻找一个“完全只由 `\tilde g` 衍生出的引力作用量”，并要求它在 `g` 很小/接近平直时导出的场方程与 `g` 表象中的 Einstein 场方程渐近一致；优先关心这种作用量是否存在、若存在应满足什么条件，而不是继续扩展含 `\rho` 或 `S` 的几何作用量。
- 2026-05-01 已按用户当前最认可的“残差判据”执行第一轮 `B/C/D` 对比：把由平直时空量子力学波函数 `(\rho,S)` 生成的 `\tilde g[\rho,S]`，直接代入三组纯 `\tilde g` 候选的 native 场方程，检查残差是否小。
- 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/analyze_bcd_residuals_from_a_reference.py`。
- 本轮粗网格设置：`nx=nz=96`, `t=0,8,16`, `\ell=100`, `M_P=300`，并沿用两个显式定义的统计掩膜：
  - `主支撑区`: `rho > 1e-3 * rho_max`
  - `高密度区`: `rho > 1e-1 * rho_max`
- 结果文件：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/bcd_residuals_from_a_reference_ell100_coarse/summary.json`
- 第一轮结论：
  - `B` 支（纯 `\tilde R`）残差远大于 `C/D`，与当前平直量子力学参考明显不相容；
  - `C` 支已显著优于 `B`，但在 `t=16` 的高密度区仍出现较强残差；
  - `D` 支（`tanh` 饱和型）在代表性统计量 `p95` 上优于 `C`，当前最符合“即使 `\tilde R` 很大，也尽量不偏离平直时空量子力学”的目标；
  - `C/D` 两支仍保留极薄过渡层尖峰，说明剩余偏差主要集中在过渡层，而不是体区普遍失配。
- 2026-05-01 已对 `C/D` 过渡结构做 `192x192` 高分辨率局部放大。当前采用三类显式定义：
  - 严格过渡层：`0.5 <= |ell^2 R_tilde| <= 2.0`
  - 粗粒度可视化过渡带：`1 <= |ell^2 R_tilde| <= 100`（仅用于画图定位）
  - 热点掩膜：主支撑区 `rho > 1e-3 * rho_max` 内残差范数前 `1%`
- 高分辨率结论：
  - 即使升到 `192x192`，严格过渡层在主支撑区内仍几乎没有采样点，说明它本身比当前网格还薄，存在明显的分辨率不足问题。
  - `t=8` 时，主支撑区残差热点位于中央重叠区，而严格/粗粒度过渡带都不在主支撑区内，因此此时的不规则热点不能主要归因于“已分辨出的过渡层”。
  - `t=16` 时，粗粒度过渡带开始进入主支撑区；对 `C` 支，热点与粗粒度过渡带的重合明显增大；对 `D` 支，重合有所增加但仍较弱。
  - 当前更稳妥的总判断是：这些不规则图案不是单纯网格伪影，也不是严格过渡层已被完全解析后的清晰几何图像，而是“真实多尺度结构 + 严格薄层仍未解析”共同造成的。
- 2026-05-01 又对 `t=16` 做了 `256x256` 的超局部三窗口放大（左弧、中央竖条、右弧），专门拆开“整图杂乱感”。
- 超局部结论：
  - 左右弧线不是随机噪声，而是镜像较好的组织化结构；`C/D` 两支在左右弧窗口里都有显著的“热点落在粗粒度过渡带上”的重合（大约 0.53~0.67）。
  - 中央竖条的 `|ell^2 R_tilde|` 结构很强，但对残差贡献很弱：`C` 的中心窗口没有热点，`D` 只有极少数孤立热点，因此中央杂乱竖纹不是当前主要误差来源。
  - 在 `256x256` 下，主支撑区内终于出现了极少量严格过渡层采样点（左右弧中仅个别像素），说明严格层开始被触及，但仍然显著欠分辨。
  - `D` 支在左右弧窗口里的代表性残差 `p95` 远小于 `C` 支，说明 `D` 的不规则图样更多是极少数尖点，而 `C` 的弧线上残差更连续地沿过渡带分布。
- 2026-05-01 已按用户要求重新导出 `C/D` 线性坐标版图集，不再对 `rho`、`R_tilde`、残差取对数：
  - 高分辨率总图目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/cd_transition_zoom_highres_linear_ell100/`
  - 超局部图目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/cd_transition_superzoom_linear_t16_256/`
  - 新脚本：`kg_examples/plot_cd_transition_superzoom.py`；`kg_examples/plot_cd_transition_zoom_highres.py` 也已改为输出线性色标图。
- 2026-05-01 已开始把“`D` 支 + 标量-张量重写 + 参考解校正 + IMEX/隐式 + 局部加密”路线落成可运行原型，而不再继续扩展纯 metric 显式推进器。
- 新原型脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_reference_imex.py`
  - 采用 `D` 支的带符号曲率代理 `chi_ref := R_tilde_ref`，不再把 `phi=f_R` 当主变量。
  - 用当前局域双高斯平直量子力学基准生成 native `tilde g_ref[rho,S]`。
  - 计算 `D` 支的 `2+1d` 迹残差、局部 IMEX 校正种子 `delta_chi_seed`、以及局部加密掩膜。
- 第一轮 `96x96` smoke test 已跑通，输出目录：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_scalar_tensor_reference_imex_prototype_96/`
  摘要在：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_scalar_tensor_reference_imex_prototype_96/summary.json`
- 第一轮原型的工作性结论：
  - 原型本身能稳定完成 `t=0,8,16` 三个参考时刻，不再像之前的纯 metric 显式推进器那样一上来就炸。
  - 局部加密掩膜只覆盖主支撑区约 `1%` 到 `1.7%`，说明“最大局部刚性”确实来自极少数局部区域，而不是整个 bulk。
  - `delta_chi_seed` 在主支撑区的代表性尺度已降到 `O(1)`，不再出现 flat-limit 下 `phi` 表示导致的病态爆大。
  - `t=0` 和 `t=16` 的建议全局显式时间步仍被极少数尖点压到极小，这再次支持“下一步必须做局部隐式/局部加密，而不是继续全局显式缩步”的判断。
- 2026-05-01 已继续实现第二层原型：局部隐式 `chi` 校正器
  - 脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_local_implicit_corrector.py`
  - 第一版“冻结导数 + 统一刚性尺度”的粗糙校正器在 `96x96` 上测试后，证明会把尖点残差推坏，不能继续沿用。
  - 第二版已改成“加密点点态数值 Jacobian + 线搜索”的局部 Newton/IMEX 校正。
  - 对 `t=8` 的 `64x64` 单时刻验证：
    - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_local_implicit_corrector_t8_64/`
    - 结果表明：在本来就很接近解的 bulk 区域，校正器会自动收敛到 `delta_chi = 0`，不会无谓地扰动参考解。
  - `t=16` 的单时刻 Jacobian 校正仍明显更重，后续应优先在局部窗口而不是整块域上做同类验证。
- 2026-05-01 已新增时间切片扫描器：
  - 脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/scan_d_reference_times.py`
  - 作用：对 `D` 支参考解做全时间粗扫描，输出主支撑区 residual `p95/p99/abs_max`、迹残差、局部加密比例、建议步长与综合 severity score。
  - `64x64`、`t=0:2:16` 的首轮扫描结果表明：不需要每个时间切片都做隐式校正；真正值得重点处理的是两类时刻：
    1. `t=0,2,4` 这类早期代表性残差偏大的切片；
    2. `t=14,16` 这类 `p95` 仍小、但 `abs_max` 和局部加密比例重新抬头的后期尖峰切片。
- 2026-05-01 已在 `t=16` 上完成 `64x64` 的点态 Jacobian 局部隐式校正验证：
  - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_local_implicit_corrector_t16_64/`
  - 结果表明：该校正器对 bulk 的 `p95` 几乎不造成扰动，但能明显压低局部尖峰：
    - 主支撑区迹残差 `abs_max`：`108.0 -> 62.7`
    - 主支撑区迹残差 `abs_mean`：`0.090 -> 0.0497`
    - 局部加密掩膜内迹残差 `p95`：`16.36 -> 6.11`
  - 这说明“坏切片 + 坏区域”的局部隐式校正是有效的，而不是只会停在 `delta_chi = 0`。
- 2026-05-01 已新增左弧固定窗口的 `C/D` 分辨率对照：
  - 脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/scan_cd_resolution_left_arc.py`
  - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/cd_left_arc_resolution_t16/`
  - 目标：不再看整图杂乱感，而是在同一物理窗口 `-7.4<=x<=-2.2, 2.6<=z<=6.2` 上比较 `64/96/128/192` 的结构收敛。
  - 结果：
    - 粗粒度过渡带的重心在 `C/D` 两支下都稳定落在左弧一带，说明左弧不是纯网格噪声，而是真实物理结构区。
    - 严格过渡层 `0.5 <= |ell^2 R_tilde| <= 2.0` 在这些分辨率下仍几乎没有主支撑区采样点，说明真正薄层仍显著欠分辨。
    - 最极端残差尖点的 `abs_max` 与位置随分辨率剧烈跳动，尚未收敛，因此这些尖点目前仍应视为“欠分辨薄层 + 数值尖峰”的混合物，而不能直接当作确定物理效应。
- 2026-05-01 又做了 `ell=10` 与 `ell=1` 的左弧窗口对照，作为“尖峰性质”诊断，而不是纯收敛测试：
  - `ell=10` 时，严格过渡层在左弧窗口内已经有明显采样，且 `C/D` 的热点会更清楚地贴着严格过渡层出现，说明 `ell=100` 下那些极薄尖峰与真实过渡层有关。
  - `ell=1` 时，左弧窗口里的孤立热点大多消失，但窗口 `p95` 仍保持较大，说明“尖峰”转化成了一整片更平滑的过渡失配，而不是凭空消失。
  - 当前最合理的综合判断是：左弧是真实物理结构区；而 `ell=100` 下最夸张的针尖型 `abs_max` 主要来自“很薄的真实过渡层在粗网格上的欠分辨表现”。
- 2026-05-01 继续把 `D` 支 `ell=100, t=16` 左弧窗口加密到 `256` 与 `320`：
  - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_left_arc_resolution_t16_hi_ell100/`
  - 结论进一步收紧：
    - 左弧窗口主支撑区 `p95` 已经基本稳定在 `~6.4e-4`，说明 bulk 代表性残差已明显收敛；
    - 粗粒度过渡带重心稳定在 `(x,z) ~ (-5.8, 4.3)` 一带，说明物理结构位置也已稳定；
    - 但最极端 `abs_max` 仍在 `1.7e6 -> 1.7e5 -> 6.9e5` 之间跳动，最大点位置仍移动，因此这些最尖针尖爆点仍不能视为已收敛物理量。
- 2026-05-01 已按用户要求暂停继续数值推进，先从理论层面评估 C/D 作为 pure-`\tilde g` 作用量的可行性。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/106-CD支作为候选纯tildeg作用量的理论评估.md`
  - 结论：C/D 是明确的 metric `f(\tilde R)` 协变候选，具有弱曲率 EH 极限和强 `\tilde R` 饱和屏蔽机制，因此比 B 支更符合“普通量子实验不显著偏离平直 QM”的目标。
  - 但 C/D 仍只能称为有效候选，不能称为已成功理论；主要理论红旗是 `f_R -> 0` 的退化/强耦合风险、过渡层中的 `\nabla\nabla f_R` 薄层项、以及正曲率区 `f_{RR}<0` 与标准 metric `f(R)` 稳定性条件的张力。
  - 当前判断：B 支可作为本目标下的反例基本排除；C 支保留为幂律饱和对照；D 支优先推进，但下一步应先检查 scalar-tensor/`chi` Cauchy 问题、线性稳定性与过渡层积分边界条件，再进入全作用量演化比较。
- 2026-05-01 已继续补充用户要求的概念解释与 `ell` 参数标定逻辑。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/107-Cauchy主部scalaron稳定性与ell标定路线.md`
  - 已明确：`Cauchy 主部` 是场方程最高阶导数结构，决定 Cauchy 初值问题是否适定；在 C/D 的 scalar-tensor 形式中，主部含 `phi=f_R` 乘在 Einstein-like 度规波算子前，因此 `f_R -> 0` 同时意味着屏蔽和主部退化风险。
  - 已明确：`scalaron` 是 metric `f(R)` 中的额外几何标量自由度 `phi=f_R`；常见健康条件为 `f_R>0` 与 `f_RR>0`。C/D 均满足 `f_R>0`，但在正 `R~` 区有 `f_RR<0`，因此作为普适 metric `f(R)` 理论有稳定性红旗。
  - 已明确：`ell` 应视为唯一普适可调长度尺度，等价曲率阈值 `R_* = ell^-2`；必须存在窗口 `|R_ordinary| << ell^-2 << |R_tilde_quantum|`。若解释为退相干有效参数，则必须另给 `ell_eff` 的独立决定律，否则会失去预测力。
- 2026-05-01 已完成第一轮 `ell-window scan` 与更细理论澄清。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/scan_cd_ell_window.py`
  - 输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/cd_ell_window_scan_64/summary.json`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/108-CD支主部退化正曲率稳定性与ell窗口首扫.md`
  - 结论：正 `R~` 区 `f_RR<0` 对 metric `f(R)` 意味着 scalaron tachyonic 倾向，但不是自动瞬时发散；是否灾难取决于增长时间、区域寿命、模式激发和 effective theory 截断。
  - 已澄清：当 `f_R -> 0` 时，`f -> 1/ell^2` 项是零阶代数项，不能作为新的健康主部；完全饱和区更像“度规动力学退化/代数约束区”，需要与 EH-like 区和 transition layer 分区匹配。
  - 首扫显示：增大 `ell` 会显著降低 C/D bulk residual，但同时让主支撑区几乎全部进入 `f_R << 1` 的饱和/退化区；因此 `ell` 选择是 residual suppression 与 Cauchy 主部健康性的张力，而不是越大越好。
- 2026-05-01 已按用户追问修正上一轮完成度表述，并补充 `M_P^2` 与过渡层 jump 的更完整分析。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/109-关于jump公式完成度与MP2钉住退化区的澄清.md`
  - 已承认：此前所谓“完成过渡层 jump 公式”只是写出迹方程积分关系，不等于完成 jump 分析；完整分析还需数值检查 `[\partial_n f_R]`、`\int \partial_n^2 f_R dn` 与 `\int |\partial_n^2 f_R|dn`。
  - 已澄清 `M_P^2`：标准平直量子力学中 `g≈eta` 的原因是大 `M_P^2` 使物质反作用 `G~T/M_P^2` 很小，而不是引力作用量项本身小；C/D 中主部有效系数是 `M_P^2 f_R`，因此退化判断必须看 `M_P^2 f_R` 与物质源/导数项的相对大小。
  - 已理解并记录用户思路：饱和退化区可尝试不作为普通 Cauchy PDE 自主演化，而用类似平直时空量子力学的选解原则，通过大 `M_P`、参考背景、边界条件或最小曲率原则固定度规；这会把理论推进为 singular-limit / matched-boundary effective theory。
- 2026-05-01 已实现第一版 `M_P^2` 有效刚度与过渡层 jump 代理量数值诊断。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_cd_transition_jump_mp2.py`
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_transition_jump_mp2_t16_160/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/c_transition_jump_mp2_t16_160/summary.json`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/110-MP2有效刚度与过渡层jump首个数值诊断.md`
  - 左弧窗口 `t=16, 160x160, M_P=300` 的结果：`ell=10` 下 C/D 的 `M_P^2 f_R` p95 仍远大于物质源 `T` p95，说明几何仍可能被有效钉住；`ell=30` 是过渡，C 仍较硬而 D 已明显软化；`ell=100,300` 进入强退化/下溢区。
  - 粗过渡带中的 `M_P^2 derivative / T` 对 `ell=10,30,100` 可达 `10^6~10^8`，支持“过渡带导数项是旧全动力学 stiffness 主要来源”的判断。
  - 严格过渡层在 `ell=100,300` 下没有被 `160x160` 采样，因此不能判断 jump 安全；`ell=10,30` 的局部 signed/absolute proxy 显示有明显正负抵消，但绝对层强度不为零。
- 2026-05-01 已将用户提出的“极薄过渡层可否视为表面/线/边界”正式整理为下一阶段数值路线。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/111-将极薄过渡层视为界面条件的数值路线.md`
  - 结论：这是正确方向，并与当前计划一致；准确说，应从“局部加密解析薄层”升级为“bulk 区域方程 + interface jump/matching 条件”。
  - 关键要求：不能简单删掉过渡层，必须保留 `[\partial_n f_R]`、`\int\partial_n^2 f_R dn`、`\int|\partial_n^2 f_R|dn`、界面两侧度规/曲率匹配和物质流守恒条件。
  - 下一步建议：优先 D 支 `ell=30,100`，在左弧窗口构造 `y=ell^2 R_tilde` 水平集界面，沿法向短线积分并检查分辨率收敛；若收敛，再将 full evolution 改写为 bulk-interface 问题。
- 2026-05-01 已回答并落实“界面固定还是动态”的问题：界面由 `y(t,x,z)=ell^2 R_tilde(t,x,z)=±1` 定义，物理上是动态水平集；左弧窗口只是当前诊断窗口。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/integrate_d_interface_jumps.py`
  - 输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_interface_jumps_t16_res_scan/summary.json`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/112-D支动态界面法向短线积分首轮结果.md`
  - 方法：在单元边上插值提取 `y=±1` 界面点，沿 `grad(y)` 法向取短线，先插值 `y(s)` 再用解析 `f_R=sech^2(y)` 重建短线上的 `f_R`，计算 endpoint jump、signed integral 与 absolute layer strength。
  - 结果：`ell=30` 的 `240/320` 界面积分开始接近稳定；`ell=100` 仍明显不收敛，说明其界面太薄或需专门局部界面算法。
- 2026-05-01 已执行用户要求的 `ell=30` bulk-interface toy system 首版。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_bulk_interface_toy.py`
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240/d_bulk_interface_toy_ell30_n240_t16.png`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/113-D支bulk-interface-toy-system首版.md`
  - 当前切片 `ell=30,t=16,n=240` 下，主支撑区约 `99.846%` 为 saturated bulk，EH-like bulk 仅约 `0.154%`；动态界面邻域约占主支撑区 `12.2%`。
  - 当前 `y=±1` interface near 覆盖 derivative top 1% hotspot 的约 `40.7%`，说明界面源抓住了相当一部分 stiff 区，但仍需多 level set 或改进界面邻域来覆盖更多热点。
- 2026-05-01 本轮重启后已按长记忆协议核对上述 `ell=30` bulk-interface toy system：
  - 已确认脚本、`summary.json`、PNG 图和 `research-notes/113` 均存在，数值摘要与日志一致。
  - 本轮没有新增理论结论；交付重点是把已执行结果清楚反馈给用户，并明确下一步为多 level set interface 覆盖率测试。
- 2026-05-01 已补充解释当前 toy 图里的 `hotspot` 与“和平直时空量子力学相差多少”：
  - 在 `prototype_d_bulk_interface_toy.py` 中，`hotspot` 严格定义为主支撑区内原始 `derivative_norm` 最大的 `1%` 网格点，不是实验热点，也不是概率密度热点。
  - 提高 interface 对 hotspot 的覆盖率，目的是确认 `bulk-interface` 模型是否抓住旧全动力学数值 stiffness 的主要来源；不能把它直接理解成提高物理可观测量拟合度。
  - 当前图使用的 `rho,S` 本身来自平直时空量子力学参考解，所以物质分布差异为零；真正差异指标是把这组参考解代入 D 支 native 场方程后的 residual。
  - `ell=30,t=16,n=240` 主支撑区 residual 统计：median `4.46e-3`，p90 `8.41e-2`，p95 `8.81e-1`，p99 `9.34e2`，max `1.31e7`。说明大多数区域残差较小/中等，但少数薄层尖点极大。
- 2026-05-01 已继续执行并修正 D 支 `ell=30` 多界面测试。
  - 发现并修正实现口径：此前 `prototype_d_bulk_interface_toy.py` 复用了左弧窗口界面提取函数，界面短线积分只在 `-7.4<=x<=-2.2, 2.6<=z<=6.2` 的左弧窗口内做；本轮已为脚本新增 `--bbox`，默认 `global`。
  - 新增 `--levels` 参数；`--levels 0.5,1,2` 会追踪 `y=±0.5,±1,±2`。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/114-D支全域多界面测试与当前数值路线澄清.md`
  - 全域单界面输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240_global_single/summary.json`
    - `interface_lines=3031`
    - `interface_near/support=0.6176`
    - `hotspot_near_interface=1.0`
    - signed integral p95 `1.58e3`
    - absolute integral p95 `7.34e3`
  - 全域多界面输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240_global_multilevel/summary.json`
    - `interface_lines=9091`
    - `interface_near/support=0.6186`
    - `hotspot_near_interface=1.0`
    - signed integral p95 `1.68e3`
    - absolute integral p95 `7.30e3`
  - 结论：需要多少个 level 不是纯拓扑问题；全域 `|y|=1` 已覆盖 top 1% derivative hotspot，多界面主要用于校准层内源项，而不是提高 hotspot 覆盖率。当前瓶颈是普通网格上的 interface-near 覆盖了约 `62%` 主支撑区，下一步应改做 subcell/body-fitted interface toy evolution。
- 2026-05-01 已开始执行 subcell/body-fitted interface toy evolution。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_subcell_interface_evolution.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/115-D支subcell界面toy演化首轮.md`
  - 方法：用 `y=ell^2 R_tilde=±1` 的动态界面作为子网格线段，不再把界面膨胀为厚带；用 cloud-in-cell 沉积线源，解 `(1-L^2 Δ)δχ=S_Γ` 的 screened-Poisson toy 方程，再代回 D 支 `2+1d` 迹方程残差。
  - `96x96` smoke 输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_subcell_interface_toy_smoke_96_seed_local/summary.json`
    - support `930`
    - segments `649`
    - before abs_max `3.86e5`
    - after abs_max `3.22e5`
    - before p95 `0.492`
    - after p95 `0.893`
  - `128x128` 复核输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_subcell_interface_toy_ell30_t16_128_seed_local/summary.json`
    - support `1669`
    - segments `888`
    - before abs_max `3.62e5`
    - after abs_max `1.63e5`
    - before abs_mean `362.1`
    - after abs_mean `338.0`
    - before p95 `1.014`
    - after p95 `1.177`
  - 结论：subcell interface 线源确实能压低最极端尖峰，但当前 frozen-background 单一 `δχ` toy 会把误差扩散到普通 bulk，使 p95 上升；下一步需要“界面压尖峰 + bulk IMEX/Newton 降 p95 + saturated 延拓”三块耦合。
- 2026-05-02 用户指出数值方法名词过多，并担心是否用人为规则降低发散、偏离 C/D 物理作用量。
  - 已明确修正口径：物理模型只有 C/D 作用量及其 Euler-Lagrange 方程；“界面法”只有在作为该方程极薄过渡层的弱形式/jump 条件时才是物理合规的。
  - `subcell/body-fitted/interface/IMEX/Newton` 等词都只是求解同一方程的数值表示或算法，不应被当作新增物理规律。
  - 上一轮的 screened-Poisson `δχ` toy correction 只能作为诊断：它说明子网格界面源能压尖峰，但不能作为物理预测或最终演化器，因为它没有完整来自 C/D 作用量变分。
  - 后续推进必须先写清：哪些方程来自 C/D 作用量，哪些 jump 条件来自穿过薄层积分，哪些只是离散求解算法；任何 ad hoc damping/clipping/residual minimization 都不能进入最终物理结论。
- 2026-05-02 已按“物理锁定”口径继续推进 D 支界面弱形式诊断。
  - 发现命名问题：`prototype_d_reference_imex.py` 里旧 `y_ref` 实际是 `tanh(ell^2 R_tilde)`，而界面变量应是 `raw_y=ell^2 R_tilde`；已新增 `raw_y_ref` 并修正 `prototype_d_subcell_interface_evolution.py` 使用 raw interface variable。
  - 因该命名问题，上一轮 `subcell_interface_toy` 结果降级为算法 smoke test，不能作为正式物理弱形式结论。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_trace_interface_weak_form.py`
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_interface_spacetime_kinematics.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/116-D支迹方程界面弱形式与时空界面必要性.md`
  - 空间法向弱形式 `128x128` 结果：
    - `normal_jump p95=1.44e-3`
    - `bulk_integral p95=2.18e-5`
    - `weak_residual p95=1.39e-3`
    - `trace_integral_direct p95=8.70e-1`
    - `direct_minus_weak p95=8.63e-1`
  - 时空运动学 `128x128` 结果：
    - coordinate normal speed median `0.109`
    - coordinate normal speed p95 `1.90`
    - raw_y_t median `2.18e3`
    - raw_y_t p95 `7.51e8`
  - 结论：固定时间切片上的空间界面弱形式不够；必须把过渡层作为时空动态界面处理，并从完整 D 支场方程跨越时空薄管积分得到 jump/matching 条件。
- 2026-05-02 已把 D 支时空界面弱形式推进到“界面速度律”层面。
  - 已修正脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_spacetime_interface_jump.py`
    - 默认 `--half-width-y` 从 `2.0` 改为 `1.0`，避免 `raw_y=+1` 的局部窗口跨到 `raw_y=-1` 层。
    - 新增 `algebraic_integral_raw_y`、`jump_law_residual_scaled` 与 `target_normal_norm_same_signature`，把原来的弱残差改写成可反求界面速度的形式。
  - 新增脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_interface_speed_law.py`
    - 对每个 `raw_y=±1` 界面点，用 leading trace/scalaron jump law 反求目标 `tilde g^{ab}F_aF_b`，再解局部二次方程得到满足 jump 条件的 `F_t` 与坐标法向速度。
    - 该脚本是 moving-interface reduced simulator 的核心模块，但还不是完整 tensor matching 全动力学。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/117-D支时空界面速度律首轮.md`
  - 主运行 `ell=30,t=16,128x128`：
    - 输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_interface_speed_law_ell30_t16_128/summary.json`
    - 总界面点 `1361`，可解点 `849`，可解比例 `0.624`。
    - 参考速度 p95 `2.516`，jump-law 修正速度 p95 `1.437`。
    - 参考法向范数 p95 `1.650e24`，修正法向范数 p95 `5.663e-3`。
    - 修正后 scaled residual p95 `2.17e-19`。
  - 时间切片对照 `96x96`：
    - `t=0` 可解比例 `0.381`，修正速度 p95 `1.003`。
    - `t=8` 可解比例 `0.365`，修正速度 p95 `2.351`。
    - `t=16` 可解比例 `0.641`，修正速度 p95 `1.326`。
  - 解释：可解点上，D 支 leading trace jump law 会把界面法向范数从参考历史中的巨大值压到稳定的 `O(10^-3)` 量级；这支持“真实 D 支界面接近特征化/近零法向范数运动”的判断。
  - 未解点不能用 clipping/damping 处理；应进入完整张量 jump/matching 或允许界面空间形状共同调整。
- 2026-05-02 已按用户要求把过渡层界面改为连续线段图，并启动 moving-interface reduced prototype。
  - 已修改：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_trace_interface_weak_form.py`
    - `find_level_segments` 现在保留每条界面线段的端点 `x0,z0,x1,z1`。
  - 已修改：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_interface_speed_law.py`
    - 速度律图从散点改为连续 `LineCollection`；
    - 可解段用颜色表示速度，未解段用黑/灰色虚线显示；
    - 以后这类图可直接用来分析界面结构，而不是误看成随机点。
  - 连续线图输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_interface_speed_law_lines_ell30_t16_128_fix/d_interface_speed_law_ell30_n128_t16.png`
    - `summary.json` 中仍为总线段 `1361`、可解 `849`、未解 `512`、可解比例 `0.624`。
  - 新增脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_moving_interface_reduced.py`
    - 每步重新抽取 `raw_y=±1` 全域界面；
    - 用 leading trace/scalaron jump law 给可解段的 `F_t`；
    - 用 cloud-in-cell 与邻近平均只做数值速度延拓；
    - 未解段不做 clipping/damping，只标记为 unresolved。
  - moving-interface reduced prototype 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_moving_interface_reduced_ell30_t16_96_fix/summary.json`
    - step 0：`866` 段，可解 `555`，可解比例 `0.641`；
    - step 6：`958` 段，可解 `657`，可解比例 `0.686`；
    - 修正速度 p95 维持约 `1.3-1.4`，修正法向范数 p95 维持约 `5.66e-3`。
  - 重要限制：该原型冻结 `metric_inv/stress_trace/rho`，只推进 `raw_y`，不是完整自洽全动力学；直接显式推进原始 `raw_y` 会暴露强 stiffness，因为 `raw_y` 梯度极大时 \(F_t=-v_n|\nabla raw_y|\) 可非常大。
  - 新结论：下一版完整模拟器不能把原始 `raw_y` 当普通平滑场推进，应改用 signed-distance / body-fitted interface 表示界面位置，并在每步由 D 支 jump/matching 条件给速度。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/118-D支连续界面线图与moving-interface原型.md`
- 2026-05-02 已按用户要求补充线性 rho 图、过渡层放大图和 signed-distance 预览。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/plot_d_interface_linear_zooms.py`
  - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_interface_linear_zooms_ell30_t16_192/`
  - 代表图：
    - `d_interface_linear_global_ell30_n192_t16.png`
    - `d_interface_linear_upper_transition_fan_ell30_n192_t16.png`
    - `d_interface_linear_central_tangle_ell30_n192_t16.png`
    - `d_interface_linear_left_leaf_ell30_n192_t16.png`
    - `d_interface_linear_right_leaf_ell30_n192_t16.png`
    - `d_interface_linear_low_density_edge_ell30_n192_t16.png`
  - 线性图三栏含义：
    1. 线性 `rho` + 全部 `raw_y=±1` 等值线；
    2. 线性 `rho` + `rho>1e-3 rho_max` 主支撑区保留界面；
    3. `raw_y` 场 + 主支撑区保留界面。
  - 断裂来源首轮定量：
    - 全域全部等值线段 `47379`；
    - 主支撑区保留线段 `2287`；
    - 支撑阈值排除 `45092`；
    - 全部等值线段中点 `rho` 中位数 `2.27e-25`，说明绝大多数全域等值线在极低密度尾部；
    - 主支撑区保留线段中点 `rho` 中位数 `7.98e-3`。
  - 解释：很多“断裂”来自主支撑区阈值截断；另外 `raw_y=±1` 是标量水平集，可闭合、分叉或在局部极值附近生成/消失，本来不需要延伸到计算边界。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_signed_distance_interface.py`
  - signed-distance 输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_interface_ell30_t16_128/d_signed_distance_interface_ell30_n128_t16.png`
  - signed-distance 定义：界面上为 0，\(|raw_y|<1\) 一侧为负，\(|raw_y|\ge1\) 一侧为正，大小是到最近支撑区界面线段的欧氏距离；它是数值重初始化变量，不是新物理场。
  - signed-distance 首轮统计：界面段数 `1361`，近界面支撑区 `|\nabla d|` p95 `0.998`，说明它可作为下一版 moving-interface solver 的候选数值变量。
- 2026-05-02 已把 moving-interface 原型升级为 signed-distance 表示。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_signed_distance_moving_interface.py`
  - 变量：
    - `d_plus=0` 表示 `raw_y=+1` 分支，初始化符号 `sign(raw_y-1)`；
    - `d_minus=0` 表示 `raw_y=-1` 分支，初始化符号 `sign(raw_y+1)`。
  - 算法：每步抽取 `d_plus/d_minus` 零水平集，用 D 支 leading trace/scalaron jump law 反求界面法向速度，按 `d_t+v_n|grad d|=0` 推进，并从新零水平集重初始化 signed distance。
  - 小时间步主输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/signed_distance_interface_overlay_ell30_n96_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/frames/signed_distance_interface_step_020.png`
  - 结果 `ell=30,t=16,n=96,dt=0.01,steps=20`：
    - step 0：界面段数 `864`，可解 `515`，可解比例 `0.596`，总界面长度 `307.053`，速度 p95 `1.129`；
    - step 20：界面段数 `712`，可解 `445`，可解比例 `0.625`，总界面长度 `214.710`，速度 p95 `3.105`。
  - 解释：signed-distance 重初始化没有导致整体崩溃，且避免了直接推进原始 `raw_y` 的显式 stiffness；但 trace-only jump law 仍留下约三分之一未解段，下一步必须补 tensor jump/matching。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/119-D支signed-distance-moving-interface原型.md`
- 2026-05-02 已完成 D 支 tensor jump 首轮诊断。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_tensor_interface_jump.py`
  - 张量薄层条件：
    \(\Delta(\partial_q f_R)(-F_aF_b+\tilde g_{ab}\tilde g^{cd}F_cF_d)+\int_{\rm layer}(f_R\tilde R_{ab}-\frac12 f\tilde g_{ab}-\tilde T_{ab}/M_P^2)dq=0\)。
  - 数值处理：
    - 界面上先插值 `metric_cov` 再现场求逆，避免 `g_ab` 与 `g^ab` 分别插值破坏迹；
    - Ricci 无迹部分冻结，迹部强制按 `R(q)=q/ell^2` 通过层变化；
    - 先用 trace 条件解 `F_t`，再检查完整张量残差；
    - 新增 free-covector rank-one 兼容性检查。
  - 主输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_interface_jump_ell30_t16_128_direction/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_interface_jump_ell30_t16_128_direction/d_tensor_interface_jump_ell30_n128_t16.png`
  - `ell=30,t=16,n=128`：
    - 总界面段 `1361`；
    - trace 数值验证通过 `614`，unresolved `747`，solved fraction `0.451`；
    - trace 目标法向范数 median `5.646e-3`、p95 `5.667e-3`；
    - trace solve 后残差 median `4.66e-10`、p95 `4.32e-6`；
    - 固定参考界面法向后的 `tensor_residual_relative` median `4.37e5`、p95 `9.08e7`，无迹部分同量级；
    - free-covector rank-one mismatch median `0.018`、p95 `0.431`；
    - 目标自由法向与当前参考界面空间法向对齐 median `0.486`。
  - 解释：trace-only 速度律稳定但不足；目标张量很多地方可以近似由某个 covector \(F_a\) 表示，但它的空间方向通常不是当前 `raw_y=±1` 参考界面的法向。下一步应写 local tensor-matching interface solver，让界面法向/形状成为未知量。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/120-D支tensor-jump首轮诊断.md`
- 2026-05-02 已按用户纠正，新增 direct tensor interface 条件，不再把 trace 作为闭合条件。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_direct_tensor_interface.py`
  - 直接张量方程：
    \(H_{ab}=-F_aF_b+\tilde g_{ab}F^2\)，其中 \(H_{ab}=-A_{ab}/\Delta(\partial_q f_R)\)。
  - 解法：由 \(F^2=\frac12\mathrm{Tr}_{\tilde g}H\) 得 \(F_aF_b=\tilde g_{ab}F^2-H_{ab}\)，直接检查右侧是否为 rank-one 实协向量外积。
  - 主输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_v2/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_v2/d_direct_tensor_interface_ell30_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_strict001/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_strict001/d_direct_tensor_interface_ell30_n128_t16.png`
  - `ell=30,t=16,n=128`：
    - 总界面段 `1361`，实 covector candidate `1317`；
    - `10%` 数值验收 accepted `456`，fraction `0.335`，accepted residual median `0.0254`、p95 `0.0887`，accepted speed p95 `1.012`；
    - `1%` 数值验收 accepted `62`，fraction `0.0456`，accepted residual median `0.00775`、p95 `0.00972`，accepted speed p95 `1.002`；
    - accepted 段的目标/当前空间法向 alignment median 仍约 `0.41-0.42`。
  - 结论：trace 不是充分条件；直接 tensor condition 显示速度幅值不是主要难点，主要难点是界面形状/法向必须由 tensor condition 重构，而不是沿当前 `raw_y` 法向推进。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/121-D支direct-tensor-interface条件.md`
- 2026-05-02 已新增 direct tensor 目标法向重构预览，避免后续滑回 trace-only。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/plot_d_direct_tensor_normal_reconstruction.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/122-D支direct-tensor目标法向重构预览.md`
  - 物理口径：`accepted` 只由完整 leading tensor condition 决定，即 \(H_{ab}=-F_aF_b+\tilde g_{ab}F^2\) 的 rank-one 实协向量和 direct tensor residual；不先解 trace speed。
  - `10%` 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/d_direct_tensor_normal_reconstruction_ell30_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/summary.json`
    - accepted `456/1361`，fraction `0.335`，accepted residual median `0.0254`、p95 `0.0887`，accepted speed median `0.9957`、p95 `1.0124`，target/current normal alignment median `0.4186`。
  - `1%` 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/d_direct_tensor_normal_reconstruction_ell30_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/summary.json`
    - accepted `62/1361`，fraction `0.0456`，accepted residual median `0.00775`、p95 `0.00972`，accepted speed median `0.9990`、p95 `1.0023`，target/current normal alignment median `0.4114`。
  - 解释：direct tensor 给出的速度幅值温和，主要矛盾仍是当前 `raw_y=±1` 界面法向与 tensor matching 目标法向不一致；下一步必须做 tensor-driven signed-distance/body-fitted interface reconstruction。
- 2026-05-02 已完成 tensor-driven 界面重构原型和 direct tensor rejected 模式分类。
  - 新脚本：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_tensor_reconstructed_interface.py`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_direct_tensor_rejection_modes.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/123-D支tensor-driven界面重构原型与失败分类.md`
  - 重构物理口径：只用 direct tensor accepted 段的目标 \((F_x,F_z)\) 构造目标切线场；`raw_y=+1/-1` 分支分开；未通过 direct tensor 的段不由 trace speed、damping、clipping 或插值补成可解段。
  - `10%` 重构输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128/d_tensor_reconstructed_interface_ell30_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128/summary.json`
    - accepted `456/1361`，重构曲线 `29` 条，分支 `-1:10,+1:19`，总长度 `165.86`，accepted 覆盖率 `1.0`。
  - `1%` 重构输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/d_tensor_reconstructed_interface_ell30_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/summary.json`
    - accepted `62/1361`，重构曲线 `16` 条，分支 `-1:4,+1:12`，总长度 `44.15`，accepted 覆盖率 `1.0`。
  - rejection mode 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128/d_direct_tensor_rejection_modes_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128_strict001/d_direct_tensor_rejection_modes_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128_strict001/summary.json`
  - `10%` 失败分类：accepted `33.5%`；`tensor_residual` `27.8%`；`rank_tail+negative+tensor_residual` `33.1%`；`no_real_covector` `3.23%`。
  - 结论：accepted 区域本身能积分成曲线，瓶颈不是重构算法；主要 unresolved 分成两类：rank-one 近似可行但完整 tensor residual 超限的段，和 rank-one/negative/residual 同时失败的段。下一步要分别检查完整层内/两侧 bulk matching 与分辨率/ell 稳定性。
- 2026-05-02 已完成 direct tensor 失败原因复查。
  - 新脚本：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_direct_tensor_failure_features.py`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_direct_tensor_least_squares.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/124-D支direct-tensor失败原因复查.md`
  - 局部特征结论：
    - 按密度分箱，accepted 从 `rho/rho_max=1e-3..1e-2` 的 `0.217` 升到 `>0.3` 的 `0.601`；
    - 按 `|grad raw_y|` 分箱，`<1e3` 时 accepted `0.827`，`1e5..1e8` 时只有 `0.070`；
    - 按局部 \(\tilde g\) 条件数分箱，`cond<50` accepted 约 `0.49-0.52`，`cond>100` 约 `0.03` 或更低。
  - 稳定性结论：
    - 分辨率 `n=96,128,160` 下 accepted `0.254 -> 0.335 -> 0.398`，说明解析度有影响；
    - 但 `rank_tail+negative+tensor_residual` 仍为 `0.400 -> 0.331 -> 0.287`，未消失；
    - `ell=10` 与 `ell=30` 在 `n=128` 下比例几乎相同；
    - `half-width=0.5` 会明显变差，`half-width=1` 与 `2` 基本相同，说明要覆盖完整过渡层但不是越宽越好。
  - 最小二乘 \(F_a\) 复查：
    - 新方法直接最小化完整张量残差，不再只用 rank-one 特征向量投影；
    - `10%` 下 accepted 从 `456/1361=0.335` 提高到 `867/1361=0.637`；
    - `1%` 下 accepted 从 `62/1361=0.0456` 提高到 `333/1361=0.2447`；
    - `tensor_residual` 型在 `10%` 下 `327/378` 被修复，在 `1%` 下 `223/277` 被修复，说明它主要是求解/投影算法问题；
    - `rank_tail+negative+tensor_residual` 大多仍失败，`10%` 下仅 `78/453` 被修复，`1%` 下仅 `44/787` 被修复。
  - 修正判断：
    - 后续不能再把 `tensor_residual` 型直接当成物理失败，应升级 direct tensor solver 为局部最小二乘；
    - 真正顽固失败集中在低密度、极陡过渡层、度规病态和 rank/negative 同时失败的段，需补完整层内 profile、两侧 bulk matching 或判定局部闭合失败。
- 2026-05-02 已把 tensor-driven 界面重构器升级为局部 least-squares \(F_a\) 求解，并完成 `n=128/160` 分辨率复查。
  - 修改脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_tensor_reconstructed_interface.py`
  - 新增参数：`--solver rank_projection|least_squares`、`--max-iter`；`least_squares` 仍求同一个 direct tensor condition，只改变 \(F_a\) 的局部数值求解方式。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/125-D支least-squares界面重构与分辨率复查.md`
  - `ell=30,t=16,n=128`：
    - row `1361`，candidate `1317`；
    - accepted `867`，accepted fraction `0.6370`；
    - accepted residual median `0.01489`、p95 `0.07362`；
    - target speed median `0.99858`、p95 `1.02843`；
    - 重构曲线 `30` 条，accepted coverage `1.0`。
  - `ell=30,t=16,n=160`：
    - row `1847`，candidate `1769`；
    - accepted `1214`，accepted fraction `0.6573`；
    - accepted residual median `0.01410`、p95 `0.06934`；
    - target speed median `0.99893`、p95 `1.02342`；
    - 重构曲线 `40` 条，accepted coverage `1.0`。
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_128/d_tensor_reconstructed_interface_ell30_n128_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_160/d_tensor_reconstructed_interface_ell30_n160_t16.png`
  - 分辨率结论：
    - 分辨率确实需要提高，因为旧 rank-one 下 `n=96,128,160` 的 accepted `0.254 -> 0.335 -> 0.398`；
    - 但主改进来自求解器升级：`n=128` 从 rank-one 的 `0.335` 升到 least-squares 的 `0.637`；
    - `n=128 -> 160` 在 least-squares 下只从 `0.637` 升到 `0.657`，说明加密是必要收敛步骤，但不能替代层内 profile、两侧 bulk matching 或真正张量闭合。
- 2026-05-02 已完成 D 支 least-squares 局部分辨率收敛与 multistart 求解器复查。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_ls_local_resolution_convergence.py`
  - 修改脚本：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_direct_tensor_least_squares.py`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_tensor_reconstructed_interface.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/126-D支局部收敛与multistart复查.md`
  - 局部收敛输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_ls_local_resolution_convergence_ell30_t16_128_160_192/d_ls_local_resolution_convergence_ell30_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_ls_local_resolution_convergence_ell30_t16_128_160_192/summary.json`
  - 自动选出的失败窗口：
    - `W1: x∈[-4,0], z∈[4,8]`
    - `W2: x∈[0,4], z∈[4,8]`
    - `W3: x∈[-4,0], z∈[8,12]`
    - `W4: x∈[0,4], z∈[8,12]`
  - `n=128,160,192` 全局 accepted fraction 为 `0.6370 -> 0.6573 -> 0.6458`，不单调。
  - 局部 accepted fraction：
    - `W1: 0.569 -> 0.601 -> 0.607`
    - `W2: 0.483 -> 0.596 -> 0.519`
    - `W3: 0.423 -> 0.453 -> 0.597`
    - `W4: 0.440 -> 0.486 -> 0.465`
  - 未解段 residual p95 仍在 `~1` 附近，全局 `1.008 -> 1.040 -> 1.010`。
  - multistart 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_multistart_ell30_t16_160/d_tensor_reconstructed_interface_ell30_n160_t16.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_multistart_ell30_t16_160/summary.json`
  - multistart 在 `n=160` 下 accepted `1214/1847=0.6573 -> 1242/1847=0.6724`，只额外修复 `28` 段。
  - 精确对比：
    - 修复段来源：`no_candidate 18`、`rank_tail+negative+residual 6`、`residual 4`；
    - 仍失败 `605` 段，其中旧初始类型最多为 `rank_tail+negative+residual 451`；
    - 仍失败段 residual median 从约 `1.000005` 到 `0.999989`，没有实质下降。
  - 当前判断：剩余失败不应优先解释为普通分辨率不足或单初值局部最小；下一步应在 `W1-W4` 做层内 profile、张量分量分解、切向/外曲率项和两侧 bulk matching 检查。
- 2026-05-02 已按用户质疑检查“低密度区失败是否可能只是计算精度误差”。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/127-D支低密度区是否为计算精度误差.md`
  - 密度分箱使用 `n=192` least-squares rows：
    - `rho/rho_max in [1e-3,1e-2)`: accepted `0.546`
    - `[1e-2,1e-1)`: accepted `0.624`
    - `[1e-1,0.3)`: accepted `0.839`
    - `[0.3,1.1)`: accepted `0.695`
  - 高密度区仍有未解段：
    - `rho/rho_max >= 1e-2`: unresolved fraction `0.308`
    - `>= 1e-1`: unresolved fraction `0.233`
    - `>= 0.3`: unresolved fraction `0.305`
    - 这些高密度未解段 residual median 仍约为 `1`。
  - `rho_floor` 敏感性：
    - `1e-10,1e-12,1e-14` 下结果完全一致到当前输出精度；
    - accepted 都是 `1214/1847=0.657282`，unresolved p95 都是 `1.03959`。
  - `probe_dt` 敏感性：
    - `2.5e-4`: `1204/1843=0.653283`，unresolved p95 `1.03709`
    - `5e-4`: `1214/1847=0.657282`，unresolved p95 `1.03959`
    - `1e-3`: `1214/1847=0.657282`，unresolved p95 `1.03959`
  - 层内 `samples` 敏感性：
    - `41/81/121/161` 下 accepted 都是 `1214/1847=0.657282`，unresolved p95 约 `1.03959`，`161` 时 `1.03667`。
  - 判断：
    - 普通计算精度误差、rho floor、时间有限差分和层内积分采样都不是 residual≈1 的主因；
    - 但低密度区仍应标记为 `low-density-untrusted`，不能用作否定 D 支的强证据；
    - 后续主证据应来自 `support-trusted` 区内的张量 residual 分量分解和 bulk matching 检查。
- 2026-05-02 已按用户纠正，修正 A 参考残差与 D 支动态匹配的解释口径。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/128-A参考残差与D支动态匹配的关系.md`
  - 用户指出：当前 jump law/residual 检查使用的是 A 分支平直量子力学 \((\rho_A,S_A)\)，不应期待它严格满足 D 分支完整动力学；若完全满足，反而意味着 A/D 在该问题上几乎动力学等价。
  - 修正结论：
    - A 参考解上的 residual 只能解释为“距离 D 支约束/匹配流形有多远”；
    - unresolved fraction 不是 D 支最终演化失败率；
    - 真正判据是能否构造 D 支一致初值，并演化 \((\rho_D,S_D,\tilde g_D)\) 与动态界面，使观测量不偏离实验太多。
  - 后续路线：
    - 以 A 实验准备态为参考，做 D 支 constrained projection / boundary matching；
    - jump law 用作动态界面条件，决定界面法向、速度及生成/消失；
    - \(\rho,S\) 必须按 D 支物质方程演化，允许并预期偏离 A 分支；
  - residual 诊断保留为初值投影难度和 A/D 偏离源定位，不再作为静态否定判据。
- 2026-05-02 已实现并运行 D 支同初态 reduced 动态演化首轮。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_reduced_dynamic_same_initial.py`
  - 修改：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/coordinate_matter_evolution.py`，把已计算的 `flux_x/flux_z` 返回给守恒型通量。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/129-D支同初态reduced动态演化首轮.md`
  - 关键实现修正：`n_cons=sqrt(|g~|) j^t` 在负频/混合时间方向下可带符号，不能裁剪为非负；Rusanov 速度必须用 `flux/n_cons`，不是 `flux/abs(n_cons)`。
  - frozen baseline 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_frozen_n64_dt1e6_active0_valid74/d_reduced_dynamic_vs_a_frozen_initial_ell30_n64_t7.4e-05.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_frozen_n64_dt1e6_active0_valid74/summary.json`
  - frozen baseline 最后可靠时间片 `t=7.4e-05`：`rho_rel_l1_support=1.5836e-3`，`rho_max_abs_support=1.4297e-2`，`phase_grad_rel_l1_support=5.878e-5`，`disc_min_support=3.52e-8`；下一步 `t=7.5e-05` 判别式变负。
  - instant transform 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_instant_n64_dt1e6_active0_stopdisc/d_reduced_dynamic_vs_a_instant_transform_ell30_n64_t2e-06.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_instant_n64_dt1e6_active0_stopdisc/summary.json`
  - instant transform 在 step 2 即出现 `disc_min_support=-2.14e11` 和 \(10^{35}\) 量级密度病态；判断为用动态 \(\rho\) 有限差分重构 \(Q=\Box sqrt(rho)/sqrt(rho)\) 的 \(1/dt^2\) 放大问题，不应作为主路线。
- 2026-05-02 用户指出并已确认：上述 D 支 reduced 动态初态没有正确执行 \(\rho_A\to\tilde\rho\) 映射。
  - 代码实际使用的是 `n_used=sqrt(|g~|) rho_A g~^{tν}u_ν`；
  - 当前混合变换、\(\kappa=m^2\) 强匹配口径下，应使用 `sqrt(|g~|) rho~ = |X| rho_A / m^2` 的正测度密度关系，或保留带符号分支 `m^2 sqrt(|g~|)rho~ = X rho_A`；
  - 因此上一轮 frozen baseline 只能作为错误密度映射下的方法学诊断，不能作为 D 支同初态物理演化与 A 支的有效比较；
  - 研究笔记 `129-D支同初态reduced动态演化首轮.md` 已追加该更正。
- 2026-05-02 已按用户三域规则完成理论可模拟性分析和 rho 映射错误审计。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/130-D支三域规则可模拟性与rho映射错误审计.md`
  - 三域规则判断：
    - 弱曲率 bulk：\(|y|=|\ell^2\tilde R|\ll1\)，可用近 Einstein/近零曲率分支；
    - 饱和 bulk：\(|y|\gg1\)，\(\phi=f_R\to0\) 且导数项消失，但 \(f\to\pm1/\ell^2\)，仍有 \(-M_P^2f\tilde g_{\mu\nu}/2\) 残留，不能无条件写成 `0=0`；
    - 过渡层：\(|y|\sim1\)，应退化成 interface，用 direct tensor jump/body-fitted matching 参与动力学。
  - 类似错误审计：
    - `simulate_d_reduced_dynamic_same_initial.py` 初态 `n_cons` 错用 \(\rho_A\)；
    - D 支 residual/jump/tensor 诊断中所有直接把 A snapshot `rho` 传给 `stress_tensor_tilde` 的 \(\tilde T_{\mu\nu}\) 源项都需重算；
    - `support=rho>...` 应改为或至少同时报告 \(\tilde N=\sqrt{|\tilde g|}\tilde\rho\) 支撑区；
    - `rho_D-rho_A` 图和指标必须替换为明确的裸 \(\tilde\rho\)、测度密度、流或拉回 observable；
    - `instant_transform` 既数值病态，也概念上误把 D 表象变量当作原表象 \(\rho\) 重构 \(Q\)。
- 2026-05-02 已实现并运行修正 \(\rho_A\to\tilde\rho\) 后的 D 支三域动力学首轮。
  - 新脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/131-D支修正rho映射后三域动力学首轮.md`
  - 初态修正：
    - 不再使用 `rho_tilde=rho_A`；
    - 使用正测度密度关系 `sqrt(|g~|)rho~=|X|rho_A/m^2`；
    - 输出比较对象改为 `ntilde_measure=sqrt(|g~|)rho~` 与 `n_cons=sqrt(|g~|)rho~ g~^{tν}u_ν`。
  - 稳定正式输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/d_tridomain_dynamics_inert_tridomain_ell300_n64_t5.7e-05_fixed.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/fields_final.npz`
  - 参数：`n=64`、`dt=1e-6`、`ell=300`、`M_P=300`、`t_completed=5.7e-05`。
  - 结果：`measure_rel_l1_support=2.29e-2`、`n_cons_rel_l1_support=5.92e-5`、`disc_min_support=4.52e-7`、`mass_shell_defect_p95=3.95e-14`。
  - 同参数尝试到 `t=5.8e-05` 时质量壳判别式变负 `disc_min_support=-4.63e-8`，因此可信截止为 `t=5.7e-05`。
  - `ell=300` 下初始支撑区全为 saturated；快速初始检查显示 `ell=1` 时支撑区仍约 `81%` saturated，`ell>=10` 时基本全 saturated。
  - 当前几何闭合为 `inert_tridomain`：bulk metric 作为惯性代表保持，interface/jump 只作诊断；这不是最终完整张量几何求解器。
  - plateau 诊断：`M_P^2/ell^2=1`，`plateau_lhs_p95=5.56e-3`、`plateau_rhs_p95=3.48e-4`、`plateau_relative_p95≈0.999`；解释为 metric 型满秩项与 rank-one matter 张量结构不匹配，而不是普通数值误差。
- 2026-05-02 记录用户提出的新理论观察：若 \(\tilde g=\eta_{\tilde{}}\) 是 \(\tilde R\) 很小区域的解，拉回 \(g\) 表象一般不意味着 \(g=\eta\)。
  - 关键判断：\(\tilde g\leftrightarrow g\) 是依赖 \(\rho,S,Q,X\) 的 disformal 场重定义，不是坐标变换；
  - massive 分支中 \(\tilde g^{\mu\nu}=g^{\mu\nu}-Q M^{\mu\nu}/(X\Delta)\)，因此固定 \(\tilde g=\eta_{\tilde{}}\) 后，\(g\) 一般是 matter-dependent 且可有非零曲率的有效 metric；
  - 只有 \(Q\to0\) 或 \(|Q|/m^2\ll1\) 且导数温和时，才有 \(\tilde g=\eta_{\tilde{}}\Rightarrow g\simeq\eta\)；
  - 这提示：\(\tilde R\) 小不是量子效应消失，而可能是量子势被转移到 \(g\) 表象的逆变换关系中。
- 2026-05-02 已按用户要求开始把数值模拟切到 1550nm 物理标定。
  - 新模块：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/physical_units.py`
  - 修改脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/132-1550nm物理标定参数准备.md`
  - 新增 CLI：`--physical-optical`、`--wavelength-nm`、`--mass-over-omega`、`--ell-over-planck`、`--dt-old-units`、`--normalize-probability`。
  - 内部单位：自然单位 \(\hbar=c=1\)，能量 eV，长度/时间 eV\(^{-1}\)；输出记录 μm/fs 换算。
  - 对 `1550nm`、自洽 `m=omega0/10`：
    - \(k_0=0.799898054\,{\rm eV}\)
    - \(\omega_0=0.803927793\,{\rm eV}\)
    - \(m=0.080392779\,{\rm eV}\)
    - \(\sigma_\parallel=2.368\,\mu{\rm m}\)
    - \(\sigma_\perp=1.776\,\mu{\rm m}\)
    - 初始中心绝对坐标 \(8.881\,\mu{\rm m}\)
    - 域半宽 \(29.603\,\mu{\rm m}\)
    - 相遇时间 \(42.10\,{\rm fs}\)
    - \(M_P=1.220890128\times10^{28}\,{\rm eV}\)、\(l_P=1/M_P=8.190745236\times10^{-29}\,{\rm eV}^{-1}\)
  - 仍需用户确认：具体 \(\ell/l_P\) 与波函数/光强归一化，否则不能把物理单位结果解释为真实实验强度预测。
- 2026-05-02 用户确认：\(\ell/l_P\) 按推荐扫描，波函数按 KG 方程归一化。
  - 修改脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 新增 `--normalize-kg`，使用正频 KG norm
    \(N_{\rm KG}=\int i(\psi^*\dot\psi-\psi\dot\psi^*)dx dz=1\)；
    对 \(\exp(-i\omega t)\) 模式等于 \(\int 2\omega|\psi|^2 dx dz=1\)。
  - 短程扫描输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_scan_n32_steps10/physical_1550nm_kg_ell_scan_summary.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_scan_n32_steps10/scan_summary.json`
  - 扫描点 `ell/l_P=1e29,1e35,1e45,1e60`：
    - `1e29`：transition `0.524`、saturated `0.476`、plateau lhs p95 `2.18e-1`
    - `1e35`：全 saturated、plateau lhs p95 `1.87e-13`
    - `1e45`：全 saturated、plateau lhs p95 `1.87e-33`
    - `1e60`：全 saturated、plateau lhs p95 `1.87e-63`，KG source rhs p95 `2.13e-60`
  - 物理标定下 64x64 首轮选 `ell/l_P=1e60`：
    - 输出：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_ell1e60_n64_stable17/`
    - 图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_ell1e60_n64_stable17/d_tridomain_dynamics_inert_tridomain_ell8.19075e+31_n64_t0.000127516.png`
    - KG norm after `1.0`
    - 可信截止 `step=17`，`t=1.275e-4 eV^-1 = 8.39e-5 fs`
    - `measure_rel_l1_support=1.542e-2`
    - `n_cons_rel_l1_support=1.809e-5`
    - `disc_min_support=2.60e-9`
    - `saturated_fraction=1.0`
    - `plateau_lhs_p95=2.06e-61`，`plateau_rhs_p95=9.05e-60`
    - 下一步 `step=18` 判别式变负 `disc_min_support=-2.51e-10`，因此不作为物理结果。
- 2026-05-02 已完成 1550nm 物理标定 D 支 active/support 边界诊断。
  - 修改脚本：`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/133-D支1550nm物理标定active边界与trusted核心诊断.md`
  - active dilation 扫描输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_active_dilation_sweep_ell1e60_n64_dt2p5em7/active_dilation_sweep_decision_plot.png`
    - `active_dilation=0` 在 step 98 质量壳判别式穿零；
    - `active_dilation=1/2/4` 均完成 180 步，其中 `2/4` 指标几乎相同；
    - `active_dilation=6` step 14 非有限值，说明过度扩张会纳入低密度病态尾部。
  - support 停止条件长跑：
    - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_steps400_dt2p5em7/`
    - 在 `active_dilation=2`、原始 support 停止条件下，step 231 停止；
    - step 231：`disc_min_support=-1.535e-10`、`measure_rel_l1_support=4.675`；
    - 失败点为 grid `(42,28)`，坐标 `(9.251,-3.700) um`，是 support 边缘但不是 active 边缘；
    - 同一时刻一格内缩 trusted core 上 `disc_min=3.865e-6`、`measure_rel=8.30e-5`。
  - 新增 `--trusted-erosion` 与 `--stop-mask support|trusted`：
    - trusted 是原始 support 按 8 邻域向内腐蚀，用于区分阈值边缘伪故障与主体演化；
    - 原始 support 指标仍保留报告。
  - trusted 核心区 400 步正式结果：
    - 输出目录：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_trusted1_steps400_dt2p5em7/`
    - 图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_trusted1_steps400_dt2p5em7/trusted_vs_support_time_diagnostics.png`
    - 图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_trusted1_steps400_dt2p5em7/trusted_final_spatial_diagnostic.png`
    - 完成 400 步，`t=4.937e-4 fs`；
    - 原始 support：`measure_rel_l1=17.981`、`disc_min=-5.554e-8`；
    - trusted core：`measure_rel_l1=1.437e-4`、`n_cons_rel_l1=9.622e-5`、`disc_min=3.867e-6`。
  - 当前判断：早期失败主要来自 active/support 阈值边缘处理，不是 trusted 核心区 D 支主体演化立即失效；下一步需做 `trusted_erosion=1` 下的分辨率和 dt 收敛。
- 2026-05-03 已完成 D 支 1550nm 物理标定时空收敛检查。
  - 新笔记：`/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/134-D支1550nm物理标定时空收敛检查.md`
  - 综合图：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_convergence_spacetime_summary/spacetime_convergence_summary.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_convergence_spacetime_summary/spacetime_convergence_summary.json`
  - 半窗口 `t≈2.4686e-4 fs`：
    - `n64,dt_old=2.5e-7`：trusted measure `7.189e-5`，support measure `7.124e-3`，support disc `1.001e-8`；
    - `n64,dt_old=1.25e-7`：trusted measure `7.189e-5`，support measure `7.119e-3`，support disc `1.001e-8`；
    - 因此时间步减半几乎无影响；
    - `n96,dt_old=2.5e-7`：trusted measure `1.448e-5`，support measure `1.648e-5`，support disc `8.080e-3`；
    - `n128,dt_old=2.5e-7`：trusted measure `1.165e-5`，support measure `1.228e-5`，support disc `6.199e-3`。
  - 全窗口 `t≈4.9372e-4 fs`：
    - `n64`：trusted measure `1.437e-4`，support measure `17.981`，support disc `-5.554e-8`；
    - `n96`：trusted measure `2.895e-5`，support measure `3.295e-5`，support disc `8.080e-3`。
  - 当前判断：主要病灶是空间分辨率/支撑边缘解析度，而不是时间步误差；后续默认工作分辨率应提升到 `n=96`，除非专门做更高精度收敛才跑 `n=128`。
- 2026-05-03 已按用户提醒检查 support 边缘表达的物理地位，并增加从 A 支任意时刻生成 D 初态的能力。
  - 修改：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/mixed_tilde_initial_data.py`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 新增：
    - `localized_direct_tilde_coordinate_initial(..., initial_time=...)`
    - `simulate_d_tridomain_full_dynamics.py --initial-time`
    - `simulate_d_tridomain_full_dynamics.py --initial-time-old-units`
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/135-D支干涉时刻初态与support边缘表达说明.md`
  - 方法判断：
    - `trusted` 只影响诊断和停止条件，不改物理方程；
    - `support -> active` 定义实际数值演化窗口，属于数值区域/外边界处理，不是新物理；
    - 平滑权重/tapered support/body-fitted 边界若只改变低密度边界表示和积分区域，则是数值方法；若删除物质或改源项则会改变物理，不能这么做。
  - 干涉时刻检查：
    - `t_old=8` 对应 `initial_time_fs=39.4977507 fs`，接近旧比例相遇时间 `~8.53`；
    - `t_old=8,n=96,dt_old=2.5e-7,steps=200` 完成；
    - 输出：
      - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_interference_time_check/interference_time_check_diagnostics.png`
      - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_interference_time_check/interference_time_check_summary.json`
    - 对比 `t_old=0,n=96`：trusted measure 从 `1.448e-5` 升到 `1.099e-3`，trusted disc 从 `3.190e-2` 降到 `4.198e-7`；
    - `t_old=8,n=128` 也完成，但没有改善：trusted measure `1.562e-3`，trusted disc `6.064e-8`；
    - 输出：
      - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_interference_resolution_check/interference_resolution_diagnostics.png`
      - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_interference_resolution_check/interference_resolution_summary.json`
  - 当前判断：
    - \(t=0\) 的参数稳定性不能外推到干涉区；
    - 干涉区困难不是单纯空间分辨率不足；
    - 下一步最小必要检查是 `t_old=8,n=96,dt_old=1.25e-7`。
- 2026-05-03 已完成 `t_old=8` 干涉窗口的时间步复查，并据此调整 support 边缘路线优先级。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/136-D支told8时间步复查与support路线调整.md`
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_dt_resolution_decision/told8_dt_resolution_decision.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_dt_resolution_decision/told8_dt_resolution_decision_summary.json`
  - `t_old=8,n96,dt_old=2.5e-7`：
    - trusted measure `1.098689981e-3`
    - trusted `n_cons` `2.360176039e-4`
    - trusted disc `4.198042102e-7`
    - support measure `1.126289735e-3`
  - `t_old=8,n96,dt_old=1.25e-7`：
    - trusted measure `1.098689979e-3`
    - trusted `n_cons` `2.360176039e-4`
    - trusted disc `4.198042101e-7`
    - support measure `1.126299837e-3`
  - `t_old=8,n128,dt_old=2.5e-7`：
    - trusted measure `1.561562759e-3`
    - trusted disc `6.063554930e-8`
  - 结论：
    - 把 `dt` 减半几乎逐项不改变诊断；
    - `n128` 没有改善干涉窗口；
    - 因而该窗口问题不是时间步太大，也不是普通空间加密能直接修复；
    - `dt_old=2.5e-7` 可继续作为当前 `inert_tridomain` 诊断默认时间步。
  - 路线调整：
    - support 平滑/taper/body-fitted 低密度边界仍需保留，但只作为数值边界卫生与收敛检查；
    - 它不能改变 \(\tilde\rho\)、\(\tilde T_{\mu\nu}\)、质量壳方程、D 支作用量或物理边界条件；
    - 主要瓶颈应转向 `inert_tridomain` 几何闭合在干涉区的不足，以及后续 metric/scalaron/interface 自洽调整。
- 2026-05-03 已按用户排队指令完成固定参数判定，并推进 support 边缘表达为后处理诊断。
  - 新脚本：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/postprocess_d_support_edge.py`
  - 修改：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
    - 新增诊断型 taper 权重、support edge band 指标；
    - `active_dilation` 改为无环绕扩张；对当前 `t_old=8` 三组已有数据，后处理检查 `active_wrap_diff_count=0`，因此不改变这些结果。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/137-D支support边缘表达后处理与固定参数判定.md`
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/support_edge_taper_overview.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/support_edge_taper_summary.json`
    - taper 宽度对照目录：
      - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess_taper0p5/`
      - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess_taper2p0/`
  - 固定参数判定：
    - `n=96` 固定为当前工作分辨率；
    - `dt_old=2.5e-7` 固定为当前工作时间步；
    - 暂不引入动态参数。
  - support 边缘表达后处理结果：
    - 0.5 decade taper 下，`n96 dt` taper measure `1.380657e-3`，`n96 dt/2` `1.380667e-3`，`n128 dt` `1.598716e-3`；
    - 1.0 decade taper 下，`n96 dt` `2.059287e-3`，`n96 dt/2` `2.059297e-3`，`n128 dt` `1.621821e-3`；
    - 2.0 decade taper 下，`n96 dt` `2.355253e-3`，`n96 dt/2` `2.355260e-3`，`n128 dt` `1.652621e-3`。
  - 判断：
    - support edge band 的误差确实略高于 trusted core；
    - 但 trusted core 本身也在 `1e-3` 量级；
    - 因此 support 边缘表达会影响诊断数值，但不是干涉窗口全部偏离的主因；
    - 它应继续作为后处理/数值边界诊断，不作为改物理或改源项的手段。
- 2026-05-03 已按用户要求继续推进饱和区代表 metric 的选择原则。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/138-D支饱和区代表度规的受限延拓原则.md`
  - 当前判断：
    - 饱和区 metric 不是任意平均，而是接口约束下的受限延拓；
    - 过渡层 boundary 是约束的一部分，但只约束边界值/法向 jump，不足以唯一决定饱和区内部代表；
    - 当前 `inert_tridomain` 仍只是第一版可执行闭合，不是唯一物理解。
- 2026-05-03 已将饱和区代表 metric 准则落实到代码。
  - 修改：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 新增几何闭合：
    - `--geometry-closure restricted_extension`
  - 新增准则：
    - 只在 saturated bulk 内部尝试更新代表 metric；
    - weak / transition / interface 侧作为物理锚点；
    - 无物理锚点的 saturated 连通块不允许被 support 外人工边界牵引，保持原代表；
    - 有锚点时用离散 harmonic / 最小 metric-gradient 延拓；
    - 延拓后必须通过 Lorentz 签名、质量壳判别式和正 \(\tilde\rho\) 检查，否则回溯或退回原代表。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/139-D支饱和区代表度规准则落地.md`
  - 验证：
    - 无锚点保护的初版最小曲率延拓会被人工 support 边界牵引，导致质量壳判别式快速穿零；
    - 加入锚点保护后，`t_old=8, ell/l_P=1e60` 被识别为全 saturated orphan bulk，因此不改代表 metric；
    - `n=96, dt_old=2.5e-7, steps=20` 与旧 `inert_tridomain` 在同一步完全一致：
      `measure_rel_l1_support=1.1397749603909711e-4`，
      `measure_rel_l1_trusted=1.1488058838534157e-4`，
      `n_cons_rel_l1_support=2.339903038536604e-5`，
      `disc_min_trusted=4.0995731809267966e-7`。
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_n96_steps20/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_n96_steps20/d_tridomain_dynamics_restricted_extension_ell8.19075e+31_n96_t3.75048e-05.png`
- 2026-05-03 已按用户新假设加入低物质边界平直锚点。
  - 理论判断：
    - 远端低物质边界可以选为平直度规，因为此时 D 支与平直量子力学近似等价；
    - 但不能把内部干涉低密度节点也当作平直边界。
  - 修改：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
    - 新增 `--flat-anchor-boundary-layers`
    - 新增 `--flat-anchor-rho-frac`
    - 新增 `--flat-anchor-measure-frac`
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/140-D支低物质边界平直锚点假设.md`
  - 测试：
    - `t_old=8, ell/l_P=1e60, n=96, steps=20`
    - `flat_anchor_count=97`
    - `metric_extension_solve_count=260`
    - `orphan_component_count=0`
    - `lorentz_fallback_count=96`
    - `metric_extension_admissible_blend=0`
  - 解释：
    - 平直边界锚点已生效，使饱和区不再是 orphan；
    - 但由其诱导的 harmonic metric update 未通过 admissible 检查，所以最终未改变内部代表；
    - 当前物质诊断仍与无 flat anchor / inert 情况一致。
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_flat_anchor_n96_steps20/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_flat_anchor_n96_steps20/d_tridomain_dynamics_restricted_extension_ell8.19075e+31_n96_t3.75048e-05.png`
- 2026-05-03 已按用户要求尝试 ADM 分块延拓，并确认质量壳判别口径。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/141-D支ADM分块延拓首轮与质量壳逐点判别.md`
  - 质量壳判别：
    - 当前代码逐点计算 `discriminant=b^2-ac`；
    - 停止条件使用 support/trusted mask 上的逐点最小值 `np.nanmin`；
    - 因此不是对支撑区做平均后判断，一个坏点也可以触发停止。
  - 修改：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_tridomain_full_dynamics.py`
    - 新增 `--metric-extension-variable-mode covariant|adm`
    - `adm` 模式把 \(g_{\mu\nu}\) 分解成 log lapse、shift、空间度规 Cholesky 变量后再做受限 harmonic extension。
  - 测试：
    - `t_old=8, ell/l_P=1e60, n=96, steps=20`
    - `geometry_closure=restricted_extension`
    - `metric_extension_variable_mode=adm`
    - flat boundary anchors 打开
  - 关键结果：
    - `metric_extension_flat_anchor_count=97`
    - `metric_extension_solve_count=5`
    - `metric_extension_solve_fraction=0.012165450121654502`
    - `metric_extension_lorentz_fallback_count=0`
    - `metric_extension_adm_valid_fraction=0.25790754257907544`
    - `metric_extension_roughness_before_p95=368.9097786294534`
    - `metric_extension_roughness_after_p95=1.6304060060012092e+53`
    - `metric_extension_admissible_blend=0`
    - `measure_rel_l1_trusted=0.00011488058838534157`
    - `disc_min_trusted=4.0995731809267966e-07`
  - 判断：
    - ADM 分块避免了 Lorentz fallback，但在当前固定实验室时间切片上只有约 `25.8%` active 点可作实 ADM 分解；
    - full ADM 变量延拓候选仍被 admissible 检查拒绝，最终不改变 metric 代表；
    - 当前下一步应收缩更新自由度，例如只更新 lapse 或 lapse+shift，或把 direct tensor interface matching 作为额外物理锚点接入。
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_adm_flat_anchor_n96_steps20/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_adm_flat_anchor_n96_steps20/d_tridomain_dynamics_restricted_extension_ell8.19075e+31_n96_t3.75048e-05.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_adm_flat_anchor_n96_steps20/support_edge_tapered_diagnostics.png`
- 2026-05-03 继续推进 ADM 收缩闭合，修正签名验收并得到首个可接受的小步 metric 松弛候选。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/142-D支分支保持ADM小步松弛闭合.md`
  - 关键修正：
    - ADM 无效点不再用 floor 后的伪 ADM 变量重构，而是保留原 metric；
    - admissible 签名检查从强制 \(+--\) 改为保持该点原有惯性分支不变；
    - 新增 `--metric-extension-adm-update-fields all|lapse|lapse_shift`；
    - 新增 `--metric-extension-max-rel-change` 作为饱和区代表 metric 小步松弛上限。
  - 诊断发现：
    - 旧 `admissible_blend=0` 的主要原因不是候选破坏质量壳或密度；
    - trusted 区初态已有 `7/128` 个点不是 \(+--\)，而候选在这些点 `delta_norm=0`；
    - 因此旧验收把已有非标准签名分支误判成候选更新失败。
  - `t_old=8,n=96,steps=20,ell/l_P=1e60,max_rel_change=1e-3` 对照：
    - inert/orphan：trusted measure `1.1488058838534157e-4`，interface solved `0.3761904761904762`；
    - ADM lapse-only：trusted measure `2.753600326353254e-4`，support measure `2.622926385981251e-4`，trusted disc `4.0995731775961275e-7`，accepted blend `1.0`，interface solved `0.44976076555023925`；
    - ADM lapse+shift：trusted measure `6.635578387816041e-4`，support measure `6.252942691344944e-4`，trusted disc `4.0995731775961275e-7`，accepted blend `1.0`，interface solved `0.4312796208530806`。
  - 判断：
    - 当前最稳 working closure 是 `branch-preserving + ADM lapse-only + max_rel_change=1e-3`；
    - 它比 inert 偏离 A 更大，但仍在 `10^-4` 到 `10^-3` 量级，没有质量壳崩溃；
    - `lapse+shift` 偏离更大，暂不作为默认优先闭合。
  - 输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_branch_preserve_cap1e3_n96_steps20/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_branch_preserve_cap1e3_n96_steps20/d_tridomain_dynamics_restricted_extension_ell8.19075e+31_n96_t3.75048e-05.png`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_shift_branch_preserve_cap1e3_n96_steps20/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_shift_branch_preserve_cap1e3_n96_steps20/d_tridomain_dynamics_restricted_extension_ell8.19075e+31_n96_t3.75048e-05.png`
- 2026-05-03 已完成 ADM lapse-only 小步上限扫描。
  - 新笔记：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/research-notes/143-D支ADM-lapse小步上限扫描.md`
  - 汇总输出：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_cap_scan/summary.json`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_cap_scan/adm_lapse_cap_scan.png`
  - 固定设置：
    - `t_old=8,n=96,steps=20,ell/l_P=1e60`
    - `restricted_extension + adm + lapse`
  - 扫描结果：
    - inert/orphan：trusted measure `1.1488e-4`，interface solved `0.3762`
    - cap `1e-4`：trusted measure `1.3096e-4`，interface solved `0.4381`
    - cap `3e-4`：trusted measure `1.6315e-4`，interface solved `0.4429`
    - cap `1e-3`：trusted measure `2.7536e-4`，interface solved `0.4498`
    - cap `3e-3`：trusted measure `5.9205e-4`，interface solved `0.4498`
  - 质量壳最小判别式各组稳定在 `4.0996e-7`，`n_cons` 误差也几乎不变。
  - 判断：
    - 当前最佳折中为 `--metric-extension-max-rel-change 3e-4`；
    - `1e-4` 太接近 inert，但很安全；
    - `1e-3` 和 `3e-3` 更新更强，interface 收益饱和但 D-A transformed measure 偏离明显增大。

- 2026-05-03 已把 full tensor interface matching 从诊断原型接入 `simulate_d_tridomain_full_dynamics.py` 的闭合链。
  - 新增的主链能力：
    - 通过 `prototype_d_tensor_reconstructed_interface.least_squares_tensor_rows(...)` 直接求解 full tensor interface 段；
    - 将 accepted tensor rows raster 成 `tensor_anchor_mask`，作为 saturated-bulk 延拓的额外物理锚点；
    - 在 `restricted_extension` 中把 tensor anchor 与 weak/transition/interface 锚点、flat boundary anchors 一起并入受限延拓；
    - 新增 `--tensor-rel-tol`、`--tensor-max-iter`、`--tensor-max-seeds`、`--tensor-multistart`、`--tensor-anchor-band-factor`。
  - 烟雾测试：
    - `python3 kg_examples/simulate_d_tridomain_full_dynamics.py --output /tmp/d_full_tensor_smoke1 --geometry-closure restricted_extension --resolution 16 --steps 1 --interface-every 1 --no-tensor-multistart --metric-extension-iterations 1 --metric-extension-max-rel-change 0 --flat-anchor-boundary-layers 0 --metric-extension-variable-mode covariant --metric-extension-adm-update-fields all`
    - 运行成功并输出 summary；
    - `tensor_interface_solved_count=0` 时，张量锚点为空，但主循环与闭合仍可正常推进；
    - `tensor_interface_solved_count=2` 的 `t=0` 步中，`tensor_interface_anchor_count=9`，说明 accepted full-tensor 段已能生成实际锚点带。
  - 当前判断：
    - full tensor matching 已经从“单纯残差诊断”提升为“参与饱和区代表选择的物理锚点层”；
    - 这仍不是唯一性证明，但比之前只看 trace jump 更接近用户要的强闭合逻辑。

- 2026-05-03 20:32 +0800 按用户纠正重定物理时间窗，并完成 D 支三切片试运行。
  - 时间窗换算：
    - 旧 `t_old in [0,16]` 在 1550nm 物理标定下对应 `T_max=120.01529382382436 eV^-1=78.99550140570793 fs`；
    - 旧 `t_old=8` 对应 `T_max/2=60.00764691191218 eV^-1=39.49775070285396 fs`。
  - 低分辨率全窗粗试：
    - `n=32, dt_old=0.1, steps=160` 未跑通，`step=14` 产生非有限 D matter RK4；
    - `n=32, dt_old=0.01, steps=1600` 仍未跑通，`step=21` 产生非有限 D matter RK4；
    - 判断：粗步长全窗失败是时间步/粗网格导致的物质 RK4 爆炸，不是命令行或张量锚点接入的语法 bug。
  - 当前可信短窗三切片均使用：
    - `physical_optical, ell/l_P=1e60, n=96, dt_old=2.5e-7, steps=20, normalize_kg`；
    - `restricted_extension + adm + lapse + metric_extension_max_rel_change=3e-4`；
    - `active_dilation=2, trusted_erosion=1, stop_mask=trusted`。
  - `T_mid=T_max/2`，每步刷新 full tensor matching (`interface_every=1`)：
    - 运行成功，无 bug，用时 `374.739 s`；
    - `measure_rel_l1_trusted=1.1488058849761347e-4`，`n_cons_rel_l1_trusted=2.3620497275342517e-5`；
    - `disc_min_trusted=4.0995731775961275e-7`；
    - tensor accepted `56/406`，`tensor_interface_anchor_count=127`，`direct_tensor_residual_relative_p95=0.09616`。
  - `T=0` 稀疏刷新 (`interface_every=20`)：
    - 运行成功，用时 `140.043 s`；
    - `measure_rel_l1_trusted=1.9452343687563538e-2`，`n_cons_rel_l1_trusted=2.011768905579346e-6`；
    - `disc_min_trusted=0.03190047390472728`；
    - tensor accepted `482/840`，`tensor_interface_anchor_count=533`；
    - 判断：早期分离态 full tensor anchor 很强，明显牵引 transformed measure，不能默认“锚点越多越好”。
  - `T=T_max` 稀疏刷新 (`interface_every=20`)：
    - 运行成功，用时 `162.458 s`；
    - `measure_rel_l1_trusted=3.5790782983999524e-5`，`n_cons_rel_l1_trusted=7.917449529927763e-7`；
    - `disc_min_trusted=5.0941967250528464e-5`；
    - tensor accepted `560/1386`，`tensor_interface_anchor_count=779`。
  - 输出目录：
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_fullwindow_trial_n32_dtold0p1_tensor_anchor/`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_fullwindow_trial_n32_dtold0p01_tensor_anchor/`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_Tmid_tensor_anchor_n96_steps20/`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_T0_tensor_anchor_n96_steps20_sparse/`
    - `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_Tmax_tensor_anchor_n96_steps20_sparse/`

- 2026-05-03 用户明确修正优先级：当前先建立完整 D 支数值计算器，不再用 D-A transformed-measure 偏差作为 full tensor anchor 的主要否决标准。
  - 新口径：
    - D-A 偏差仍可作为后续实验可观测差异诊断；
    - 但当前阶段判断一个 D 支演化器是否成立，应优先看 D 支内部自洽：
      质量壳逐点有实根、\(\tilde\rho>0\)、守恒量稳定、D 支场方程/约束残差、full tensor interface matching、饱和区代表 metric 规则、以及界面/ bulk 的一致拼接。
  - 对此前 `T=0` full tensor anchor 使 D-A 偏差到 `~2%` 的解释也随之修正：
    - 这不再自动表示 full tensor anchor “失败”；
    - 它只表示 D 支若按该闭合演化，可能早期就与 A 支产生明显差异；
    - 是否接受该闭合，应改由 D 支自身方程残差和约束闭合来决定。

- 2026-05-04 继续按用户要求优化 D 支数值计算，目标是寻找 1 天内覆盖完整旧 `[0,16]` / 物理 `[0,T_max]` 的方法，并特别检查 16G MacBook 内存风险。
  - 系统检查：未发现遗留 Python 大算例进程；Codex/Chrome 等占用较多内存；后续测试均用 `/usr/bin/time -l` 记录峰值 RSS。
  - 已修改 `kg_examples/simulate_d_tridomain_full_dynamics.py`：
    - 新增 `--geometry-every`：分频刷新昂贵的 Ricci / `R_tilde` / 三分区几何缓存；
    - 新增 `--metric-extension-every`：分频执行 saturated-bulk metric 代表延拓；
    - 新增 `--diagnostics-every`：分频计算 A-reference 对比、plateau algebraic residual 和记录行；
    - 新增 `--no-render`：计时扫描时跳过绘图，减少内存和运行时间；
    - 修复 summary：`final` 现在记录真实终止步，`last_diagnostic_record` 记录最后一次采样诊断，避免把早停后的最后采样误当最终状态；
    - 新增实验性 `--matter-variable-mode phase`，但当前 naive 相位差分与谱导数初值不一致，`t_old=8` 第 0 步质量壳判别式即为负，暂不作为生产路线；
    - 新增实验性 `--matter-projection` / `--matter-projection-mask` / `--matter-positivity-limiter`，用于测试显式 RK 后的质量壳与正密度约束保持。
  - 性能结果：
    - `n=64, dt_old=2.5e-7, steps=20, geometry/extension/diagnostics every 20, no-render`：完成，用时 `3.82s`，最大 RSS `77MB`；
    - `n=96` 同设置：完成，用时 `7.67s`，最大 RSS `96MB`；
    - `dt_old=1e-5, n=96, geometry/extension/diagnostics every 100` 的理想速度若能稳定，完整旧 `[0,16]` 约 `1.6e6` 步，按早停前速度可在数小时内完成，算力/内存本身满足 1 天目标。
  - 稳定性结果：
    - `dt_old=1e-5,5e-6,2.5e-6,2.5e-7` 均在几乎同一物理时间 `t_old≈0.00137` 附近触碰质量壳/正密度边界；缩小步长只缩小负判别式幅度，不改变失效物理时间；
    - `inert_tridomain` 也在同一时间失败，说明不是 `restricted_extension` / ADM lapse-only 延拓导致；
    - 失败点位于 trusted/support 内但接近低密度支撑边缘，典型点 `x≈-25.0,z=0`，`rho_A≈3.96e-6`、`measure_A≈1.47e-3`、`disc≈-3.4e-10`；trusted 内只有 1 个点先越界，active 外层已有更多负密度/负判别式点；
    - 支撑阈值提高到 `3e-3` 或 `1e-2` 可推迟早停，但最终仍出现非有限 RK4 和巨大负密度，说明单靠 support 边界表达不能根治；
    - `trusted_erosion=2` 可移除首个坏点，但只保留 support 内约 `70%` transformed measure，不能作为完整物理演化域；
    - 正性/质量壳投影和 limiter 可把早停推迟到 `t_old≈0.0029~0.0042`，但随后仍卡在 `rho=0` 或 `disc=0` 的约束边界附近，并显著增加单步开销。
  - 当前判断：
    - 16G 内存不是优化版主限制；峰值 RSS 基本 `80~105MB`，远低于 16G；
    - 真正阻塞一日全窗的是物质变量选择：当前用实验室时间切片上的守恒密度 `n=sqrt(|gtilde|) rho_tilde j^t` 反解 `rho_tilde`，当局部 `j^t` 接近或穿过 0 时会病态，显式 RK4 不能保持 `rho_tilde>0` 与 `disc>=0`；
    - 因此当前脚本经分频优化后“算力上可进 1 天”，但“数学上仍不能严谨跑完整窗”。下一步应换成避免 `j^t=0` 病态的物质积分器：局部协动/特征时间切片、semi-Lagrangian/finite-volume 正性保持格式、或以非负 `rho_tilde` 加质量壳根作为主变量而不是从 `n_cons` 病态反解。

- 2026-05-04 继续推进 D 支物质积分器退化诊断。
  - 新增脚本：`kg_examples/diagnose_d_jt_degeneracy.py`。
    - 功能：对 A 参考切片生成的 `tilde g, rho_tilde, u_mu` 计算 `jt_unit = gtilde^{t nu} u_nu`、`n_cons=sqrt(|gtilde|)rho_tilde jt_unit`、质量壳判别式，并统计 support/trusted 内符号、零线、近零测度占比；
    - 输出图和 summary：
      - `visualizations/d_physical_1550nm_jt_degeneracy_n96/summary.json`
      - `visualizations/d_physical_1550nm_jt_degeneracy_n96/jt_degeneracy_summary.png`
      - `visualizations/d_physical_1550nm_jt_degeneracy_scan_n96/summary.json`
      - `visualizations/d_physical_1550nm_jt_degeneracy_scan_n96/jt_degeneracy_summary.png`
  - 三切片结果：
    - `t_old=0` support 内 `jt+ 158 / jt- 258`，trusted 内 `jt+58 / jt-176`，trusted 负测度占比 `0.961`；
    - `t_old=8` support 内 `jt+64 / jt-196`，trusted 内 `jt+43 / jt-85`，trusted sign-change edges `72`，`|jt|<1e-3` 的 trusted measure fraction `0.00713`；
    - `t_old=16` support 内 `jt+114 / jt-746`，trusted 内 `jt+57 / jt-581`，trusted 负测度占比 `0.976`。
  - 全时间扫描 `t_old=0..16`：
    - `jt` 符号混合从一开始就存在；
    - 近零问题在 `t_old≈8..13` 最严重，`trusted |jt| p01` 低至 `1.589e-4`，`|jt|<1e-3` 的 trusted measure fraction 最高约 `9.14%`；
    - 因此当前实验室时间 `t` 上的 `n_cons` 不是全局良好密度变量。
  - 尝试固定倾斜时间 covector `d tau = dt + a dx + b dz`：
    - 粗网格搜索存在候选如 `a≈3,b≈1.7`，可降低全窗加权近零比例，整体 major sign measure fraction 约 `0.935`；
    - 但没有单一 `(a,b)` 使所有时刻、所有 trusted/support 内同号且远离零；
    - 结论：单一倾斜时间只能作局部图册，不能作为完整全窗 Cauchy 时间。
  - 当前路线判断：下一步应实现多 patch / 特征推进框架，而不是继续在实验室时间 `n_cons` 变量上做全窗显式 RK4。

- 2026-05-04 继续推进 D 支全窗优化与相位钟改造。
  - 已把 `kg_examples/phase_clock_projection.py` 改成局部支持集加权重建 `project_to_weighted_gradient`，并把 `diagnose_d_phase_time_clock.py` 默认切到 `--projection-mode weighted`；
  - 这一步不改物理方程，只改相位钟诊断的数值重建方式，避免周期边界 wrap-around 污染局域波包；
  - 已给 `kg_examples/simulate_d_tridomain_full_dynamics.py` 增加 `--fast-profile`，默认关闭渲染、稀疏几何/诊断/界面刷新，并关闭 tensor multistart，专门用于一日全窗成本评估；
  - smoke test 结果：
    - `diagnose_d_phase_time_clock.py --projection-mode weighted` 在 `32x32` 上运行正常，`j^S` 仍以 `~1e-14` 相对误差锁住 `m^2`；
    - `simulate_d_tridomain_full_dynamics.py --fast-profile --resolution 32 --steps 20 --no-render` 用时约 `3.10 s`，峰值 RSS 约 `68 MB`；
    - `simulate_d_tridomain_full_dynamics.py --fast-profile --resolution 64 --steps 10 --no-render` 用时约 `7.15 s`，峰值 RSS 约 `80 MB`；
  - 当前判断：
    - 16G MacBook 的内存不是主要瓶颈；
    - 一日全窗的主要问题变成了步数规模和物质积分器病态，而不是几何缓存或渲染；
    - 这条路线的下一步仍应是局部时间 / 特征推进 / 正性保持积分器，而不是继续靠加内存或降图像质量硬撑。

- 2026-05-04 应用户要求复核“距离完整数值计算器还有多久”的状态。
  - 已确认当前代码层已经具备：
    - D 支三分区/界面原型；
    - 物理参数版本的 1550nm/KG 归一化输入；
    - full tensor anchor / interface matching 诊断入口；
    - 几何/诊断分频和 `--fast-profile` 性能扫描；
    - 相位钟的局部加权重建诊断。
  - 但当前还不能称为完整 D 支数值计算器，原因不是内存或绘图，而是实验室时间变量下的守恒密度会在 `j^t` 接近或穿过 0 时病态，导致 `rho_tilde` 反解、质量壳实根和正密度保持失败。
  - 当前最短闭环应直接实现一个替代物质推进器：不用 `n_cons=sqrt(|gtilde|)rho_tilde j^t` 作为主变量，而改用局部时间/特征推进/正性保持格式，并把 full tensor matching 作为强约束或强诊断接入。

- 2026-05-04 继续按用户要求实现 D 支局部时间/特征推进原型，避免实验室时间 `j^t≈0` 病态。
  - 已修改 `kg_examples/simulate_d_tridomain_full_dynamics.py`：
    - 新增 `--matter-variable-mode local_time`；
    - 新增 patchwise 局部时间图册 `--local-time-patchwise`，每个 patch 使用 `tau=t+a x+b z` 推进；
    - 主变量改为正的 transformed measure `measure_density=sqrt(|gtilde|)rho_tilde`，不再从实验室 `n_cons/(sqrt|g|j^t)` 反解密度；
    - 输出仍重构到实验室 `t` 切片，保存 `measure_density_D`、`n_cons_D`、`u_t/u_x/u_z` 等诊断；
    - 新增 `--local-time-stationary-candidates`：对每个网格点解析求使 `gtilde^{mu nu}tau_mu tau_nu` 最大的局部时间候选，避免粗 `(a,b)` 网格漏掉窄时间锥；
    - 新增 `--local-time-patch-bboxes`：每个 patch 只在带一层 halo 的局部包围盒上做 RK4，保持中心差分 stencil，不改变物理方程。
  - 重要诊断：
    - 旧 patchwise 粗候选在 `n=64,t_old=0` 有 `14.36%` trusted transformed measure 无法被 `tau_norm>1e-12` 覆盖；
    - 扫描显示不是 `max_tilt` 太小，而是合法时间锥很窄，`0.5` 粗步长跳过了解析最佳点；
    - 加入解析 stationary candidates 后，同一切片 trusted 覆盖率升到 `100%`。
  - 测试结果：
    - `python3 -m py_compile kg_examples/simulate_d_tridomain_full_dynamics.py` 通过；
    - `n=16` smoke：local_time patchwise + stationary candidates 通过；
    - `n=64,t_old=0,steps=20,ell/l_P=1e60,KG norm`：
      - 非 bbox 版用时约 `100.93s`，RSS 约 `116MB`；
      - bbox 版用时约 `17.68s`，RSS 约 `117MB`；
      - trusted `measure_rel_l1=2.65e-6`，`disc_min_trusted=3.86e-6`，`rho_tilde_min_trusted=1.39e-5`；
      - `local_time_admissible_fraction=1.0`，`uncovered_trusted_measure_fraction=0`，`patch_count=84`，`candidate_count=222`，`bbox_area_fraction≈0.0224`。
    - `n=32,t_old=8,steps=20` 干涉中心切片：
      - 非 bbox 版用时约 `50.91s`，RSS 约 `84.5MB`；
      - bbox 版用时约 `13.55s`，RSS 约 `87.4MB`；
      - trusted `measure_rel_l1=2.63e-6`，`disc_min_trusted=7.05e-3`，`rho_tilde_min_trusted=1.59e-5`；
      - `local_time_admissible_fraction=1.0`，`uncovered_trusted_measure_fraction=0`，`patch_count=139`，`candidate_count=319`，`bbox_area_fraction≈0.0997`。
  - 当前判断：
    - 这是第一个真正避开实验室 `j^t=0` 密度反解病态的 D 支物质推进原型；
    - 它仍不是完整生产级全动力学计算器，因为 patch 边界还只是 mask/recombine，尚未实现严格 conservative patch-boundary flux，也尚未接入 `restricted_extension`/full tensor metric 强匹配；
    - 但它已经把主瓶颈从“数学变量病态”推进到“性能与 patch 间保守通量/几何闭合”。

- 2026-05-04 继续执行 local_time 步长扫描；中途用户断网后复查确认没有残留进程，100 步长窗已完整结束，无需重跑。
  - 扫描文件：
    - `visualizations/d_local_time_dt_scan_t0_n32/scan_summary.json`
    - `visualizations/d_local_time_dt_scan_t8_n32/scan_summary.json`
    - `visualizations/d_local_time_dt_scan_t16_n32/scan_summary.json`
    - `visualizations/d_local_time_dt_verify_n64/verify_summary.json`
    - `visualizations/d_local_time_dt_refine_t8_n64/scan_summary.json`
    - `visualizations/d_local_time_dt_refine_t8_n64_tilt5/scan_summary.json`
    - 汇总索引：`visualizations/d_local_time_dt_scan_compact_summary.json`
  - `n32` 快扫：
    - `t_old=0` 到 `dt_old=1e-4` 完成，质量壳裕度约 `3.36e-3`；
    - `t_old=8` 到 `dt_old=2.5e-4` 完成，质量壳裕度约 `7.05e-3`；
    - `t_old=16` 到 `dt_old=5e-4` 完成，但裕度仅 `~2e-7` 量级；
    - 结论：`n32` 对大步长过于乐观，不能单独决定生产步长。
  - `n64` 复核：
    - `t_old=0,dt_old=1e-4` 完成，trusted `measure_rel≈1.06e-3`，`disc_min≈3.88e-6`；
    - `t_old=8,max_tilt=3,dt_old=1e-4` 4 步早停，coverage 约 `0.94`，说明默认 tilt 太窄；
    - `t_old=16,dt_old=1e-4` 虽完成但 transformed measure 爆到 `~1.49e4`，不可用；
    - `t_old=16,dt_old=5e-5` 完成且 trusted `measure_rel≈2.64e-3`，`disc_min≈3.30e-6`，较可信。
  - 对 `t_old=8,n64` 的图册诊断：
    - `max_tilt=3,step=0.5` 下约 `2.67%` trusted transformed measure 未覆盖；
    - 把 `max_tilt` 放到 `5` 后覆盖恢复为 `100%`；
    - 因此把 local_time 默认 `--local-time-max-tilt` 从 `3` 改为 `5`。
  - `t_old=8,n64,max_tilt=5` 细扫：
    - `dt_old=5e-6`：20 步完成，trusted `measure_rel≈1.01e-4`，`disc_min≈1.32e-8`；
    - `dt_old=1e-5`：20 步完成，trusted `measure_rel≈2.01e-4`，`disc_min≈1.33e-8`；
    - `dt_old=2.5e-5`：20 步完成，trusted `measure_rel≈5.00e-4`，`disc_min≈6.19e-9`；
    - `dt_old=5e-5`：20 步完成但 `disc_min≈3.34e-11`，太贴近质量壳边界，不作为安全步长。
  - 诊断修正：
    - `local_time_measure_rel_delta` 只是 `sqrt(|gtilde|)rho_tilde` 在倾斜时间推进中的变化，不是局部时间守恒密度；
    - 新增 `local_time_n_tau_rel_delta_sum`，对应真正的局部坐标守恒量 `sqrt(|gtilde|)rho_tilde j^tau` 的 patch-summed 相对变化；
    - 敏感例 `t_old=8,n64,max_tilt=5,dt_old=2.5e-5,steps=20` 中：`measure_delta≈0.113`，但 `n_tau_rel_delta_sum≈7.25e-6`，所以不能把前者误判成守恒失败。
  - 长窗检查：
    - `t_old=8,n64,max_tilt=5,dt_old=2.5e-5,steps=100` 完成；
    - trusted `measure_rel≈5.16e-3`，`disc_min≈1.60e-8`，`rho_min≈1.95e-6`，coverage `1.0`，`n_tau_rel_delta_sum≈1.10e-5`；
    - wall time `~830s`，user time `~194s`，RSS `~173MB`；内存仍不是问题，但 Python 原型速度离全窗生产还很远。
  - 当前判断：
    - 安全候选步长暂定 `dt_old≈1e-5`；
    - `dt_old=2.5e-5` 可作为激进探针，但全窗生产前必须先继续优化 patch atlas/patch grouping/局部通量；
    - 下一步应减少 patch 数或缓存候选，并实现真正 conservative patch-boundary flux，而不是继续只扫步长。

- 2026-05-05 继续按用户要求寻找 1 天内完整旧 `[0,16]` / 物理 `[0,T_max]` 的 D 支 local_time 演化路线，并复查断网后的执行状态。
  - 进程检查：未发现遗留 `simulate_d_tridomain_full_dynamics.py` Python 进程，之前 100 步运行已结束。
  - 已修改 `kg_examples/coordinate_matter_evolution.py`：
    - `conservative_density_from_rho_u`、`recover_rho_from_conservative_density_and_u`、`coordinate_matter_rhs_covector` 支持传入预计算的 `metric_inv_txz` / `det_cov_txz`；
    - 这只是复用同一度规的逆和行列式，不改变质量壳方程或守恒方程。
  - 已修改 `kg_examples/simulate_d_tridomain_full_dynamics.py`：
    - local_time patch RK4 使用解析坐标变换得到局部 `tau=t+a x+b z` 的协变/逆变度规，避免每个 patch 反复做 3x3 SVD/pinv；
    - 每步复用同一 `metric_inv/det`，避免 `current_from_measure_and_u` 与 chart 选择重复求逆；
    - `--fast-profile` 不再把用户显式设为 `0` 的 `--diagnostics-every`、`--interface-every`、`--metric-extension-every` 强行改回 8；因此可以测真正纯演化成本；
    - 新增实验开关 `--local-time-tile-size` 和 `--local-time-atlas-every`，但二者暂不作为生产默认。
  - 编译检查：`python3 -m py_compile kg_examples/coordinate_matter_evolution.py kg_examples/simulate_d_tridomain_full_dynamics.py` 通过。
  - 关键性能结果（`t_old=8,n64,dt_old=2.5e-5,ell/l_P=1e60,KG norm, local_time patchwise, stationary candidates`）：
    - 原 `q=0.1` 同口径约 `29.5s/20 steps`；
    - 缓存逆度规后约 `24.5s/20 steps`；
    - 解析局部坐标度规变换后约 `17.56s/20 steps`；
    - 关闭中间诊断/界面/几何刷新后的纯演化约 `8.07s/20 steps`；
    - 进一步复用每步逆度规后约 `7.50s/20 steps`；
    - 100 步纯演化为 `21.30s`，峰值 RSS 约 `151MB`，`disc_min_trusted≈1.602e-8`，coverage `1.0`，`local_time_n_tau_rel_delta_sum≈1.01e-5`。
  - 性能摘要保存：`visualizations/d_local_time_performance_summary.json`。
  - 全窗成本估计：
    - 以 100 步纯演化 `0.213s/step` 粗估，`dt_old=2.5e-5` 完整 `[0,16]` 约 `37.9h`；
    - 扣除端点几何开销后，物质步进约 `0.17s/step`，同一步长约 `30h`；
    - 若使用更安全的 `dt_old=1e-5`，仍约数天，不能宣称已进入 1 天。
  - 否决/保留的优化：
    - `--local-time-tile-size>1` 会漏覆盖或让 trusted measure 偏差爆炸；不作为生产路线；
    - 粗 `--local-time-step=1.0/1.25/2.0` 虽快但漏覆盖、负判别式或 measure 爆炸；不合格；
    - `--local-time-step=1.5` 搭配更细 stationary quantization 仍约 `0.95%` trusted measure 未覆盖；不合格；
    - `--local-time-atlas-every>1` 保持覆盖但不明显提速，且 `20` 步复用时 `n_tau` 守恒诊断恶化；不作为主优化；
    - 试写轻量 RHS 与原 RHS 逐项等价但无明显提速，已撤回以避免维护两套方程。
  - 当前判断：
    - 在 Python/NumPy patchwise 原型中，不牺牲局部时间覆盖的安全优化已基本耗尽；
    - 要进入 1 天全窗，下一步更可能需要把 patch RK4 小数组循环迁移到编译/向量化后端，或改成真正 conservative multi-patch finite-volume/特征通量形式，而不是继续粗化 chart 图册。

### 2026-05-05 并行与快速调试提速复查

- 按用户要求继续检查并行与此前几个提速方案是否还能进一步降低调试等待时间。
- 新增/修改：
  - `kg_examples/simulate_d_tridomain_full_dynamics.py` 增加 `--quick-diagnostics` 和 `--skip-fields-npz`，短窗调试时跳过端点曲率/interface/plateau/A-reference 重型诊断和压缩场输出；这只改变诊断/输出，不改变演化方程；
  - 同一文件把 `matplotlib` 改为真正渲染时才延迟导入；
  - 新增 `kg_examples/run_d_local_time_parallel_scan.py`，用于多进程并行跑独立短窗扫描，并汇总每个 run 的 `summary.json`、`run.log`、质量壳、正密度、local-time coverage 与 `n_tau` 守恒诊断。
- 编译检查通过：
  - `python3 -m py_compile kg_examples/simulate_d_tridomain_full_dynamics.py kg_examples/run_d_local_time_parallel_scan.py`
- 快速诊断基准：
  - `t_old=8,n64,dt_old=2.5e-5,steps=20` 常规 fast-profile 为 `real 7.44s`；
  - 加 `--quick-diagnostics --skip-fields-npz` 后为 `real 3.88s`；
  - 延迟导入 matplotlib 后为 `real 3.77s`；
  - 稳定性指标一致：`disc_min_trusted=6.186e-09`，coverage missed `0`，`local_time_n_tau_rel_delta_sum=6.667e-06`。
- 并行扫描基准：
  - 2 进程并行 `t_old=0,8`、各 `50` 步总耗时 `10.94s`；单独 `t_old=8` 50 步约 `10.88s`，说明可把另一个短窗检查几乎并入同一等待时间；
  - 4 进程并行 `t_old=0,4,8,12`、各 `30` 步总耗时 `14.98s`，但单任务明显变慢，说明 CPU/内存带宽开始竞争；
  - 当前建议默认 `--max-workers 2`，最多 `4`，不要默认开满 `8`。
- 图册复用复查：
  - `atlas_every=1,2,5,10,20` 的 50 步扫描结果保存到 `visualizations/d_local_time_atlas_reuse_benchmark/atlas_reuse_summary.json`；
  - 50 步中 `atlas_every=20` 可从 `12.18s` 降到 `10.03s`，短窗指标基本不变；
  - 但 100 步复查中 `atlas_every=20` 的 `local_time_n_tau_rel_delta_sum` 恶化到 `0.123`，不可作为生产默认；
  - 生产/最终图仍应使用 `atlas_every=1` 或做严格对照。
- profile 结论：
  - quick 模式后主要耗时为 `rk4_patchwise_local_time_matter_step` / `rk4_local_time_matter_step`；
  - 20 步 profile 中约 `2.82s/3.89s` 在 patch RK4，约 `0.84s` 在 chart 选择；
  - 单条完整长演化若要再大幅提速，下一步应走编译/JIT patch 内核、真正 conservative patch flux 或新的特征推进格式，而不是继续关诊断或粗化图册。

### 2026-05-05 旧窗口 `t=0,4,8,12` 各约 5 分钟完整端点诊断

- 按用户要求，对旧切片 `t_old=0,4,8,12` 各跑约 5 分钟墙钟时间的 D 支 local_time 短窗诊断。
- 运行口径：
  - `physical-optical`，`ell/l_P=1e60`，KG normalization；
  - `n64`，`dt_old=2.5e-5`，`local_time patchwise + stationary candidates`；
  - 不使用 `quick-diagnostics`，因此保留端点 A-reference transformed-measure 偏差；
  - 关闭作图和 `fields_final.npz`，只为减少输出开销，不改演化方程。
- 输出汇总：
  - `visualizations/d_old_slices_5min_full_diagnostics/five_min_diagnostics_summary.json`
- 四个任务全部 `returncode=0` 且 `stopped_reason=completed`，即程序层面均未崩溃。
- 结果摘要：
  - `t_old=0`：1200 步，演化旧时间 `0.03`，trusted A 偏差 `1.624%`，但 `local_time_n_tau_rel_delta_sum≈2.419` 且 coverage missed `0.4458%`，判作数值红旗；
  - `t_old=4`：800 步，演化旧时间 `0.02`，trusted A 偏差 `1.115%`，coverage `0` missed，`n_tau≈2.27e-6`，当前最正常；
  - `t_old=8`：750 步，演化旧时间 `0.01875`，质量壳仍正且 coverage `0` missed，但 trusted A 偏差爆到 `1.139e34`，support/tapered 偏差更大，判作严重失败/爆发；
  - `t_old=12`：600 步，演化旧时间 `0.015`，trusted A 偏差 `1.205%`，但 tapered 偏差 `34.47`、`n_tau≈0.0119`、coverage missed `1.43e-4`，判作边缘/尾部风险。
- 当前判断：
  - “能正常完成进程”不等于物理/数值可信；
  - `t_old=4` 是当前 5 分钟短窗最干净切片；
  - `t_old=0` 与 `t_old=12` 需要检查 patch coverage / support-edge / rho_floor 触发；
  - `t_old=8` 需要优先定位 transformed measure 爆发原因，尤其要检查是否为局部时间 patch 重组、低密度 floor、或干涉相消区放大导致。

### 2026-05-05 `t_old=8` 五分钟演化图像化

- 用户指出 `t=8` 是干涉条纹最显著附近，要求展示 A 支和 D 支结果图。
- 重跑 `t_old=8, steps=750`，这次保存 `fields_final.npz` 以便自定义作图：
  - 输出目录：`visualizations/d_old_t8_5min_rendered/`
  - `summary.json` 与 `fields_final.npz` 已生成；
  - 初次命令同时给了 `--render --fast-profile`，由于 `fast-profile` 自动关闭 render，内置图未生成；随后基于 `fields_final.npz` 手动生成 A/D 对比图。
- 图像文件：
  - `ad_measure_rho_linear_robust.png`：A/D transformed measure 与 rho 的线性 robust scale 对比；
  - `ad_measure_full_linear_outlier.png`：完整线性尺度，显示 D outlier 主导色标；
  - `ad_ratio_rawy_location_linear_clipped.png`：D/A ratio、`raw_y` 与 support/trusted mask 的位置诊断。
- 关键数值：
  - `D_measure_max≈2.085e55`，位置约 `x=-8.326um,z=-0.925um`；
  - 该最大点不在 support/trusted 内，但 support 内 `D_measure_p99.5≈4.113e32`，说明不是单个外部孤立点；
  - `A_measure_max≈8.51e-2`，`A_measure_p99.5_support≈5.83e-2`；
  - support 内 `D/A` 比值 p99 约 `8.05e33`。
- 当前判断：
  - `t_old=8` 的爆发不是单纯图像色标或单个外部点问题；
  - support 内也已经有极端 transformed-measure 爆发，需要定位其与干涉相消区、patch label 边界、`rho_floor`/低密度区和局部时间通量重组的关系。

### 2026-05-05 修正 D-A 比较口径：应拉回 A/g 表象比较 `rho`

- 用户指出此前直接比较 `sqrt(|gtilde|)rho_tilde` 与 A 支不符合他要看的物理量；正确口径应是把 D 支变量先拉回 A/g 表象得到 `rho_D_to_A`，再与 `rho_A` 比较。
- 已据此基于 `visualizations/d_old_t8_5min_rendered/fields_final.npz` 重算：
  - `rho_D_to_A = m^2 * (sqrt(|gtilde|)rho_tilde)_D / |X_g[D]|`
  - `X_g[D]=u_t^2-u_x^2-u_z^2`，平直 A/g 表象下 `sqrt(|g|)=1`。
- 输出：
  - `corrected_rho_A_vs_D_pulled_back_linear.png`
  - `corrected_rho_A_vs_D_pulled_back_full_linear.png`
  - `corrected_rho_pullback_stats.json`
- 关键数值：
  - trusted 相对 `L1`：`rho_rel_l1_trusted≈9.475e33`；
  - support 相对 `L1`：`≈4.122e40`；
  - support 内中位比值 `rho_D_to_A/rho_A≈1.0011`；
  - support 内 p99 比值 `≈8.068e33`；
  - `rho_D_to_A_max≈1.873e52`，最大点不在 support/trusted 内；
  - trusted 内 `|X_g[D]|` 未接近零：`min≈0.0227`，p01≈`0.0390`。
- 修正后的判断：
  - 需要撤回“直接用 transformed measure 作为物质分布比较”的表述；
  - 但在正确拉回 A/g 表象后，`t_old=8` 仍存在极端爆发；
  - 该爆发不是因为 `|X_g[D]|` 在 trusted 内接近 0 简单造成的；
  - 大多数 support 点可能仍接近 A（中位比值约 1），问题集中在少数高分位/局部区域。
- 应用户要求，补做 A 与 D->A 共用 A 色标范围的图：
  - `corrected_rho_same_A_colorscale_p995.png`：A 与 D->A 都使用 `rho_A` 在 support 上的 p99.5 作为线性色标上限；
  - `corrected_rho_same_A_colorscale_fullmax.png`：A 与 D->A 都使用 `rho_A` 全局最大值作为线性色标上限；
  - `same_A_colorscale_stats.json` 保存色标数值与 D 饱和区域统计。
- 同色标统计：
  - `A_full_max≈0.0015556`；
  - `A_p99.5_support≈0.0011127`；
  - D->A 超过 A p99.5 色标的 active 点数 `178`，约占 active 区 `17.87%`；
  - D->A 超过 A 全局最大色标的 active 点数 `176`，约占 active 区 `17.67%`。
- 用户指出同色标图中并看不出“主支撑区域几乎相同”；复核后确认此前用 support 中位数表达“多数点接近”过于误导。
- 新增核心区分层统计：
  - `rho_pullback_core_region_stats.json`
  - `corrected_rho_same_A_scale_with_core_contours.png`
- 修正判断：
  - support 中位比值确实约 `1.001`，但这不能代表主支撑整体正常；
  - support 中约 `21.1%` 点满足 ratio `>2`，约 `19.9%` 点满足 ratio `>10`，约 `17.7%` 点超过 A 全局最大色标；
  - trusted 中约 `19.9%` 点满足 ratio `>2`，约 `19.0%` 点满足 ratio `>10`，约 `16.9%` 点超过 A 全局最大色标；
  - 高密度核心也有明显失控：`rho_A>=0.1 max` 的区域中约 `38.2%` 点 ratio `>2`，约 `35.3%` 点 ratio `>10`；
  - 因此应撤回“主支撑区域内几乎相同”的表述；正确说法是“中位数接近，但支撑/核心内存在相当比例的严重局部失控”。

### 2026-05-05 A 支 `t=8` 看不到干涉条纹的原因复查

- 用户指出同色标图里 A 支看起来只有一个波包、没有干涉条纹。
- 已复查当前图的时间和分辨率：
  - 渲染图对应绝对旧时间 `t_old=8.01875`；
  - 当前 1550nm 缩放参数下理论相遇时间为 `t_meet_old≈8.5280`，因此当前图还没到最强干涉中心；
  - `n=64` 时空间步长约 `0.925um`，而波长为 `1.55um`，只有约 `1.68` 个网格点/波长，低于 Nyquist 稳妥显示条纹所需，肉眼很容易看成单个包络。
- 已生成 A 支分辨率参考图：
  - 输出目录：`visualizations/a_branch_t8_resolution_check/`
  - 比较 `n=64,128,256` 在 `t_old=8.0`、`t_old=8.01875`、`t_old=t_meet≈8.5280` 的 A 支 `rho`；
  - `n=128` 约 `3.35` 点/波长，`n=256` 约 `6.70` 点/波长，更适合肉眼看 1550nm 条纹。
- 修正判断：
  - 当前 `n=64,t_old≈8.02` 图不适合作为“是否有干涉条纹”的视觉证据；
  - 要展示光学干涉条纹，应至少用 `n=128`，最好 `n=256`，并取接近 `t_meet_old≈8.528` 的 A 支切片；
  - 但 D 支完整五分钟演化若直接升到 `n=256`，计算成本会显著上升，应先用 A 支高分辨率图解释视觉问题，再决定是否做 D 支高分辨率短窗。

### 2026-05-05 干涉区分辨率硬门槛修正

- 用户指出 `n=64` 对 1550nm 干涉问题分辨率太低，不能作为可信 D 支物理模拟；该纠正成立。
- 当前两束波包的干涉条纹周期不是 `lambda=1.55um`，而是由两束斜入射波矢差决定：
  - 两束波矢约为 `(±k0/sqrt(2), k0/sqrt(2))`；
  - `|Delta k|=sqrt(2) k0`；
  - 条纹周期 `lambda_fringe=lambda/sqrt(2)≈1.096um`；
  - 旧单位下 `period_old=2*pi/(sqrt(2)*6)≈0.7405`。
- 因此当前全域 `old_half_range=20,n=64` 的 `dx_old=0.625` 只有约 `1.18` 点/条纹，远不足以分辨干涉振荡。
- 新的分辨率计划：
  - 以 `dx_old=0.0625` 作为主目标，约 `11.85` 点/干涉条纹；
  - 若保留旧全域 `[-20,20]`，需 `n=640`；
  - 为节省成本，主计划改用更小空间窗口。
- 当前建议主窗口：
  - 物理/旧时间窗口优先选 `[5,12]`，包含分离、进入干涉、离开干涉；
  - 空间窗口先取 `old_half_range=12`，这样仍包含 `t=0` 初始包的主要尾部，适配当前 FFT 正频 A 参考构造；
  - 主分辨率 `n=384`，因为 `2*12/384=0.0625 old units`；
  - 收敛阶梯：`n=256` 粗探针，`n=384` 主跑，`n=512` 或缩小窗口后的等效更细网格作确认。
- 需要捡回/重审的低分辨率误判：
  - `t_old=8,n64` D 支爆发只能说明低分辨率原型在干涉区失控，不能作为 D 支物理失败证据；
  - local_time atlas/tile/coarse chart 的失败可能部分被低分辨率相位/导数 aliasing 放大，应在高分辨率窗口重测；
  - support-edge/tapered support、rho_floor、低密度节点和 patch 边界的诊断都可能被一格接近条纹周期的网格污染；
  - full tensor/interface matching 失败比例、metric representative 选择困难也需要在解析条纹后复核。

### 2026-05-05 A 支时间原点平移到干涉中心并生成确认图

- 用户要求把 `t=8` 即干涉中点移到新 `t=0`，并先展示 A 支“初始、干涉、结束”图确认。
- 修正说明：
  - 在当前物理参数下真正包心相遇时间为旧 `t_meet≈8.5280`，不是旧 `8.0`；
  - 新时间定义为 `tau_old = old_t - t_meet_old`，因此 `tau=0` 对应旧 `t≈8.5280`。
- 已生成 A 支高分辨率参考图：
  - 输出目录：`visualizations/a_branch_shifted_interference_center/`
  - 全域计算 `old_half_range=20,n=640`；
  - `dx_old=0.0625`，干涉条纹周期 `0.74048 old`，约 `11.85` 点/条纹；
  - 显示裁剪 `old x,z in [-8,8]`，对应约 `±11.84um`；
  - 展示时刻：`tau=-3.5,0,+3.5`，对应旧 `t≈5.028,8.528,12.028`。
- 生成文件：
  - `a_rho_shifted_tau_montage_same_colorscale.png`
  - `a_rho_interference_center_zoom.png`
  - `a_rho_interference_center_lineout_z0.png`
  - `summary.json`
- 等用户确认 A 支窗口/时刻后，再系统复盘此前低分辨率 `n64` 可能导致的误判，并决定 D 支高分辨率小窗口模拟参数。

### 2026-05-05 用户澄清窗口为物理 `[-9,9]um`

- 用户澄清所说 `[-9,9]` 是物理坐标 `um`，不是旧无量纲坐标。
- 已重新生成 A 支 `[-9,9]um` 窗口三时刻图：
  - 输出目录：`visualizations/a_branch_shifted_interference_center_um9_window/`
  - `a_rho_shifted_tau_montage_window_pm9um.png`
  - `a_rho_interference_center_window_pm9um.png`
  - `a_rho_interference_center_lineout_pm9um.png`
  - `summary.json`
- 换算：
  - `[-9,9]um` 对应旧坐标半宽约 `6.0805`；
  - 若维持 `dx_old=0.0625`，该窗口每方向需要约 `195` 个点；实际可取 `n=192` 稍粗或 `n=208/224` 留边。
- A 支窗口内总量比例：
  - `tau=-3.5`：`99.945%`；
  - `tau=0`：`99.99997%`；
  - `tau=+3.5`：`99.593%`。
- 当前判断：
  - 物理窗口 `[-9,9]um` 足够覆盖 A 支分离-干涉-分离主体，边界为弱尾部；
  - D 支下一轮高分辨率小窗口可以围绕该物理窗口设计，优先用同一新时间 `tau=old_t-t_meet`。

### 2026-05-05 D 支高分辨率局部窗口诊断

- 已新增脚本：
  `kg_examples/simulate_d_local_window_from_a_snapshot.py`
- 该脚本的关键修正：
  - 不再把小窗口本身当成周期 FFT 初始域；
  - 先在旧全域 `[-20,20]`、`n=640` 上构造 1550nm A 支正频 KG 精确快照；
  - 再裁剪到物理窗口 `[-9,9]um`；
  - 用裁剪后的 `rho,S,gtilde` 初始化 D 支局部时间推进。
- 当前新参数：
  - 物理波长 `1550nm`；
  - `m=omega0/10`；
  - KG norm 归一化为 1；
  - 裁剪窗口实际网格 `195 x 195`；
  - `dx≈0.09251um`，旧单位 `dx_old=0.0625`；
  - 干涉条纹周期约 `1.096um`，约 `11.85` 点/条纹。
- 已确认比较口径：
  - 不能直接拿 D 的 `sqrt(|gtilde|)rho_tilde` 和 A 的 `rho` 比；
  - 正确主比较为
    `rho_D_to_A = m^2 * (sqrt(|gtilde|)rho_tilde)_D / |X_g[D]|`，
    其中 `X_g[D]=u_t^2-u_x^2-u_z^2`。
  - 初态该拉回误差在 support 上为 `~2.6e-17`，说明初始化表象变换本身对齐。
- `dt_old=2.5e-5` 的 smoke test 结论：
  - 在干涉中心 `tau=0` 一步后即出现负质量壳判别式；
  - 这不是窗口或 KG 初始化错误，而是该步长对高分辨率干涉节点/近退化区域太大；
  - 后续 5 分钟诊断改用此前可信的 `dt_old=2.5e-7`。
- `tau=-3.5` 分离态 5 分钟诊断结果：
  - 输出目录：`visualizations/d_local_window_snapshot_tau_m3p5_5min_dt2p5e7/`
  - 跑满 wall time，`steps=167`；
  - 旧时间推进量 `4.175e-05`；
  - `rho_D_to_A` 相对 L1：
    support `5.25e-4`，trusted `5.32e-4`，core `4.85e-4`；
  - 质量壳判别式最小值：
    support/trusted `2.77e-7`，未变负；
  - 注意：该时刻 support/active 触到窗口边界，边缘仍可能受边界口径影响。
- `tau=0` 干涉中心诊断结果：
  - 输出目录：`visualizations/d_local_window_snapshot_tau0_5min_dt2p5e7/`
  - 未跑满 5 分钟，`steps=28` 后停止；
  - 停止原因：trusted 区质量壳判别式变为 `-6.22e-9`；
  - 放宽停止阈值到 `1e-8` 的补测只多走到 `steps=29`，随后 `disc_min_trusted≈-1.66e-8`，说明不是单纯停止阈值过严；
  - core 区 `rho_D_to_A` 相对 L1 仍只有 `~7e-6`；
  - 但 support/trusted 的 L1 分别约 `0.112`/`0.082`，主要被低密度、近 `X_g[D]=0`、节点/支撑边缘点主导。
- 当前判断：
  - 高分辨率后，分离态局部 D 推进在当前短窗内可稳定运行；
  - 干涉中心暴露出真正的局部刚性/近退化问题，不能再用低分辨率 `n=64` 结果定性，但也不能说完整 D 演化已经可用；
  - 下一步应把干涉中心的失败点按 `rho`、`X_g[D]`、判别式和局部时间 chart 覆盖率定位出来，再决定是收缩 evolution mask 到物理 core、换特征/正性保持格式，还是引入对近 `X_g=0` 分支边界的专门处理。

### 2026-05-05 D 支局部窗口边界/缓冲层修正

- 已在 `kg_examples/simulate_d_local_window_from_a_snapshot.py` 中补充：
  - `chart_boundary`、`step_active`、`evolve_base` 等最终掩码保存到 `fields_final.npz`；
  - 图中增加红色 `chart_boundary` 轮廓；
  - `summary.json` 中新增 `evolved_trusted`、`evolved_trusted_no_chart_halo`、outlier 质量权重、pointwise p50/p95/p99/max 等诊断；
  - 修正图中文字中自适应步长下 `advanced_old` 仍按 `step*dt_old` 计算的误导性显示。
- 诊断确认：
  - 原 `tau=0` 干涉中心巨大 support L1 主要来自 support 非 trusted 的低密度边缘孤立爆点；
  - 这些爆点不是主概率支撑区整体物理偏离，也不是 KG 初始化/表象拉回错误；
  - 若显式推进区仍取 `active`，即使冻结不可覆盖 chart 单元，support L1 仍可被 1-2 个边缘点抬到 `~0.57`。
- 已新增关键选项：
  - `--evolve-mask active|support|trusted`
  - 当前推荐短窗物理诊断用 `--evolve-mask trusted`，即只把 trusted 内部当作显式 bulk 推进区，低密度 support/active 外缘作为局部窗口边界/缓冲层。
  - 这不是阻尼、削峰或改物理源项；它只是局部窗口数值边界处理，避免把低密度边界/过渡层单元错误当作普通 bulk 点推进。
- 新规则下 `tau=0` 干涉中心 120 秒 smoke：
  - 输出目录：`visualizations/d_local_window_snapshot_tau0_trusted_evolve_smoke120/`
  - `steps=111`，无需 step halving，`disc_min_trusted≈1.69e-9>0`；
  - `rho_D_to_A` 相对 L1：support `1.34e-4`，trusted `1.35e-4`，core `2.83e-5`；
  - 旧规则同条件 support L1 为 `~0.57`，说明边缘爆点已被处理。
- 新规则下两组 5 分钟诊断：
  - 分离态 `tau=-3.5`：
    输出目录 `visualizations/d_local_window_snapshot_tau_m3p5_trusted_evolve_5min/`；
    `steps=109`，`disc_min_trusted≈2.80e-7`；
    `rho_D_to_A` 相对 L1：support `3.55e-4`，trusted `3.67e-4`，core `3.29e-4`；
    pointwise max 相对误差 `~4.05%`，无 `>10%` outlier；
    注意该时刻 support 触碰 `[-9,9]um` 窗口边界，所以该运行可作稳定性诊断，不能单独作为最终边界设置。
  - 干涉中心 `tau=0`：
    输出目录 `visualizations/d_local_window_snapshot_tau0_trusted_evolve_5min/`；
    `steps=165`，`disc_min_trusted≈4.12e-10>0`；
    `rho_D_to_A` 相对 L1：support `2.00e-4`，trusted `2.02e-4`，core `4.16e-5`；
    transformed measure 相对 L1：support `1.37e-3`，trusted `1.39e-3`；
    只有 1 个点相对误差 `>1`，其 `rho_A` 权重约 `2.18e-5`。
- 当前结论：
  - `[-9,9]um,n_full=640` 下，D 支局部窗口物质推进的主要技术阻塞已从“干涉中心不可运行”缩小为“如何给低密度/过渡边界提供正式边界条件”；
  - trusted bulk 内部短时演化已经稳定，且正确拉回到 A 表象后偏差保持在 `10^-4` 量级；
  - 下一步若继续全动力学，应把 `trusted`/边界缓冲处理升级为正式 conservative boundary/interface flux，而不是回到 active 全域显式推进。

### 2026-05-06 D 支局部窗口接入 full tensor interface 诊断

- 已按“不能把物质通量 matching 冒充 full tensor matching”的修正，更新：
  `kg_examples/simulate_d_local_window_from_a_snapshot.py`
- 新增能力：
  - `--tensor-interface-diagnostics`；
  - `--ell`、`--ell-over-planck`、`--mp`；
  - `--half-width-y`、`--samples`、`--tensor-rel-tol`、`--tensor-max-iter`、`--tensor-max-seeds`、`--tensor-multistart`、`--tensor-anchor-band-factor`；
  - 运行后保存 `tensor_interface_rows.jsonl`、`tensor_interface_diagnostic.png`、`tensor_anchor_mask` 和 `tensor_raw_y`。
- 物理/数值口径：
  - full tensor interface 诊断直接检查 `raw_y=ell^2 Rtilde=±1` 界面段是否满足 leading thin-layer tensor jump condition；
  - accepted 段表示局部最小二乘 \(F_a\) 的 full tensor residual 小于容差，可作为后续几何闭合的强锚点候选；
  - rejected 段不能被强行当作物理边界；它们需要界面形状、两侧 bulk metric 或更完整耦合规则来处理；
  - 当前局部窗口的 metric 仍固定为初始 A-derived metric，所以这是 strong diagnostic/anchor-export，不是完整 coupled metric evolution。
- 验证：
  - `python3 -m py_compile kg_examples/simulate_d_local_window_from_a_snapshot.py kg_examples/simulate_d_tridomain_full_dynamics.py kg_examples/prototype_d_tensor_reconstructed_interface.py` 通过。
- 高分辨率 `tau=0`、`[-9,9]um`、`n_full=640`、`dt_old=2.5e-7`、`ell/l_P=1e60` 诊断：
  - 0 步输出：
    `visualizations/d_local_window_tau0_full_tensor_diag_steps0/`
  - 20 步输出：
    `visualizations/d_local_window_tau0_full_tensor_diag_steps20/`
  - 两者 full tensor 统计一致：`segment_count=13926`，`candidate_count=13264`，`solved_count=3842`，`solved_fraction≈0.2759`，`tensor_anchor_count=6058`；
  - accepted residual：median `≈0.04855`，p95 `≈0.09013`，max `≈0.09998`，在当前 `tensor_rel_tol=0.1` 内；
  - 20 步物质推进仍稳定：`rho_pullback_rel_l1_trusted≈2.798e-05`，`measure_rel_l1_trusted≈1.838e-04`，`disc_min_trusted≈5.31e-08`。
- 关键结论：
  - 不是“最小二乘完全找不到解”，而是 full tensor condition 只接受部分界面段；
  - 下一步应把 accepted tensor rows/mask 接入几何代表选择或界面贴体闭合，而不是把所有 `raw_y=±1` 段统一强匹配。

### 2026-05-06 D 支初值约束投影首轮可行性检查

- 用户提出新要求：重新确定初值，使 \(\tilde\rho\) 拉回到 \(g\) 表象后的 \(\rho\) 与 A 支逐点加权偏离最小，同时 \(\tilde g\) 表象满足 D 支方程。
- 已新增脚本：
  `kg_examples/analyze_d_initial_projection_feasibility.py`
- 该脚本目前做最保守的首轮检查：
  - 固定当前 A-derived \(\tilde g\)、\(u_\mu\) 与 metric jets；
  - 对 \(\sqrt{|\tilde g|}\tilde\rho\) 做全局缩放 \(\alpha\)，检查拉回 \(\rho\) 目标；
  - 同时比较 full tensor jump 几何代数张量范数与物质源项 \(T/M_P^2\) 的量级。
- 运行：
  - 输入案例：`visualizations/d_local_window_tau0_full_tensor_diag_steps20/`
  - 输出目录：`visualizations/d_initial_projection_feasibility_tau0_steps20/`
  - 图：`d_initial_projection_feasibility.png`
  - summary：`summary.json`
- 结果：
  - 固定几何与 \(u_\mu\) 时，加权 L1/L2 的最佳全局缩放均为 `alpha=1.0`；
  - 此时 `rho_weighted_rel_l1≈2.77e-05`，`rho_weighted_rel_l2≈1.69e-03`；
  - full tensor accepted 仍为 `3842/13926≈27.6%`；
  - 界面行上 `||T/M_P^2||` 中位数约 `9.82e-61`；
  - 几何代数张量范数中位数约 `41.66`；
  - 比值 `||T/M_P^2|| / ||geometric tensor||` 中位数约 `1.04e-62`；
  - 若只靠缩放 \(\tilde\rho\) 让物质源影响 D tensor jump 到一阶量级，中位数需要放大约 `9.66e61` 倍。
- 结论：
  - 只调整 \(\tilde\rho\) 不能同时最小化拉回 \(\rho\) 偏离并修复 D 支 full tensor 方程；
  - 真正的 D 初值投影必须主要调整 \(\tilde g\) 的几何初值、metric jets、界面形状或饱和区 representative；
  - 不能把“重新确定初值”简化成重新缩放物质密度。

### 2026-05-06 D 支方程认证模式

- 用户追问：能否保证初态和后续演化以及数值模拟方法满足 D 支方程；若可以就实施。
- 明确答复并实施可复现认证：
  - 当前代码不能数学保证完整 D 支方程严格成立；
  - 可实施的是 residual certificate，对当前已实现的 D 支方程/约束逐项认证。
- 已新增脚本：
  `kg_examples/certify_d_branch_equations.py`
- 认证范围：
  - \(\tilde g^{\mu\nu}u_\mu u_\nu=m^2\) 质量壳；
  - \(\tilde\rho\ge0\)；
  - 拉回 \(\rho\) 加权偏差目标；
  - full tensor interface rows；
  - saturated-bulk 近似代数方程 \(-\frac12 f(\tilde R)\tilde g_{\mu\nu}\approx \tilde T_{\mu\nu}/M_P^2\)。
- 对当前最好的 `tau=0`、20 步案例运行：
  - 输入：`visualizations/d_local_window_tau0_full_tensor_diag_steps20/`
  - 输出：`visualizations/d_local_window_tau0_full_tensor_diag_steps20_certificate/d_branch_equation_certificate.json`
- 结果：
  - `certified=false`；
  - 通过项：
    - mass shell p95 `≈3.16e-16 < 1e-10`；
    - \(\tilde\rho>0\)；
    - pullback \(\rho\) weighted L1 `≈4.02e-06 < 1e-3`；
  - 失败项：
    - full tensor interface：`3842/13926≈27.6%`，未达到要求 `100%`；
    - saturated-bulk algebraic relative residual p95 `≈0.9999`，未达到 `1e-2`。
- 结论：
  - 当前局部窗口演化不能称为满足 D 支完整方程；
  - 若要“保证”，下一步必须先写真正的 D 支几何约束/边值求解器或每步隐式约束投影器，至少要解决 full tensor interface 和 saturated-bulk algebraic 两个失败项。

### 2026-05-06 D 支硬约束初值求解器首轮

- 用户要求：实施真正的 D 初值求解器；求解后再用此前数值方法测试 full tensor matching、D 方程残差和 \(\tilde\rho\) 拉回 \(\rho\) 与 A 的偏差演化。
- 已新增脚本：
  `kg_examples/solve_d_initial_data_constrained.py`
- 该脚本不是 penalty optimizer，而是 hard-constrained gate：
  - 在饱和 bulk 使用当前三域 D 支规则的必要方程
    \(-\frac12 M_P^2 f(\tilde R)\tilde g_{\mu\nu}=\tilde T_{\mu\nu}\)；
  - 同时要求物质质量壳成立；
  - 若非零物质处 RHS 为秩一而 \(\tilde g_{\mu\nu}\) 必须满秩，则返回 `status=infeasible`，拒绝继续演化假候选。
- 运行结果 1：
  - 输入：`visualizations/d_local_window_tau0_full_tensor_diag_steps0/`
  - 输出：`visualizations/d_constrained_initial_solver_tau0/`
  - 结果：`status=infeasible`
  - 支撑区 `11806/11806` 全部 saturated；
  - matter saturated count `11806`；
  - rank obstruction candidate count `11806`；
  - infeasible count `11806`；
  - infeasible rho mass fraction `1.0`；
  - full tensor interface solved fraction 仍为 `3842/13926≈0.2759`。
- 运行结果 2：
  - 输入：`visualizations/d_local_window_tau0_full_tensor_diag_steps20/`
  - 输出：`visualizations/d_constrained_initial_solver_tau0_steps20/`
  - 同样 `status=infeasible`；
  - 旧候选 20 步后依然是支撑质量 `100%` 触发饱和区秩障碍。
- 关键结论：
  - 在当前“精确饱和 bulk 代数方程 + on-shell dust-like 标量物质”的 D 支三域模型下，A-like 光学初值不是“最小二乘没调好”，而是必要条件层面不可行；
  - 因此本轮没有运行后续 D 演化测试，因为那会变成演化非 D 方程解。

### 2026-05-06 修正：plateau 代数秩障碍不是完整 D 方程不可行证明

- 用户指出：饱和区方程应只是渐进成立；如果没有左右都趋于 0，则直接使用代数 plateau 方程并不合法；完整方程中的导数项可能补足张量秩。
- 已修正：
  - `kg_examples/solve_d_initial_data_constrained.py` 默认不再把 strict plateau algebraic rank test 当成 hard rejection；
  - 新默认状态为 `requires_full_fR_equation_solve`；
  - 只有显式传 `--strict-plateau-algebraic-hard` 才会把该极限诊断作为硬门槛。
- 已新增完整方程诊断脚本：
  `kg_examples/diagnose_d_full_fr_equation_local_window.py`
- 该脚本直接计算：
  \[
  \phi \tilde R_{\mu\nu}-\frac12 f\tilde g_{\mu\nu}
  +(\tilde g_{\mu\nu}\tilde\Box-\tilde\nabla_\mu\tilde\nabla_\nu)\phi
  -\tilde T_{\mu\nu}/M_P^2
  \]
  其中 \(f=\tanh(\ell^2\tilde R)/\ell^2\)，\(\phi=f_R=\mathrm{sech}^2(\ell^2\tilde R)\)。
- 对 `tau=0`、1550nm、`ell/l_P=1e60`、`dt_old=2.5e-7` 局部窗口运行：
  - 输出：`visualizations/d_full_fr_equation_tau0/`
  - 主支撑区 `11806` 个点全部 saturated；
  - `|raw_y|` 最小约 `2.60e60`，p50 约 `3.41e64`；
  - \(\phi=f_R\) 在双精度下全为 `0`；
  - derivative term norm 在支撑区全为 `0`；
  - algebraic norm p50 约 `1.04e-62`；
  - RHS norm p50 约 `1.04e-60`；
  - relative residual p50 约 `0.992`，p95 约 `0.9999`。
- 更正后的结论：
  - 秩障碍只否定“平滑饱和 bulk 中把完整方程退化为纯代数 plateau 方程”的强近似；
  - 它不否定完整 D 支方程；
  - 但在当前 A-derived 光学初态主支撑区，\(\phi\) 及其导数项已经极端饱和/数值为零，所以完整方程在该候选切片上实际没有可见的导数补秩项。

### 2026-05-06 \(\kappa\) 是否能使 \(\tilde T^{\mu\nu}\) 满秩

- 用户追问：检查 \(\kappa\) 的取值能否使 \(\tilde T^{\mu\nu}\) 为满秩。
- 口径确认：
  - 此处按物质壳常数 \(\tilde g^{\mu\nu}u_\mu u_\nu=\kappa\) 理解 \(\kappa\)，不是早期 \(F(Q)\) 里的尺度参数。
- 理论结论：
  - 物质应力张量为
    \[
    \tilde T_{\mu\nu}
    =2\tilde\rho u_\mu u_\nu
    -\tilde g_{\mu\nu}\tilde\rho(\tilde X-\kappa).
    \]
  - 若执行 \(\tilde\rho\) 变分得到的壳方程 \(\tilde X=\kappa\)，第二项消失；
  - 因而 on-shell 有 \(\tilde T_{\mu\nu}=2\tilde\rho u_\mu u_\nu\)，秩为 1（或 \(\tilde\rho=0/u=0\) 时更低），升指标不改变秩；
  - 所以任何固定 \(\kappa\) 都不能让满足物质方程的 \(\tilde T^{\mu\nu}\) 满秩。
- 数值对照：
  - 当前 `tau=0` 支撑区中 \(\tilde X\) 中位数/当前 \(\kappa=m^2\) 约 `0.006462998960209225`；
  - on-shell \(\kappa=\tilde X\) 时，支撑区 `11806/11806` 均 rank 1；
  - 若故意 off-shell 取 `κ=0`、`κ=0.25κ0` 或 `κ=2κ0`，支撑区约 `11803/11806` 点变成 rank 3；
  - 但这些 off-shell 选择的壳缺陷 p50 分别约 `0.00646`、`0.00485`、`0.00646`，不是物质场方程解。
- 结论：
  - \(\kappa\) 可以通过违反壳方程人为引入满秩 metric 项；
  - 但在当前 rank-1/dust-like 物质作用量的自洽解上，\(\kappa\) 不能解决 \(\tilde T\) 秩不足问题；
  - 若想 RHS 满秩，必须改物质模型/耦合结构，例如加入压力/势能型项、多个独立相位梯度或不让 on-shell 拉氏量项消失。

### 2026-05-06 plateau 代数式的秩/行列式论证再次修正

- 用户指出：
  \[
  -\frac12 M_P^2 f\tilde g_{\mu\nu}=2\tilde\rho u_\mu u_\nu
  \]
  中，右边 rank 1 行列式为 0，但左边若 \(f\to0\)，行列式也趋于 0，不能仅凭行列式说矛盾。
- 修正：
  - 用户判断正确；行列式极限不是合适判据；
  - 对任意固定有限 \(f\neq0\) 且非退化 \(\tilde g\)，左边作为矩阵仍满秩，只是整体范数很小；
  - \(f\to0\) 时矩阵族 \(f\tilde g\) 的行列式和范数趋零，但 rank 在极限处不连续；
  - 因此严格矛盾只适用于“固定有限 \(f\neq0\) 且要求逐点精确等式”的口径；
  - 在 \(f\to0\)、\(\tilde T/M_P^2\to0\) 的渐近口径下，方程可退化成近似 \(0\simeq0\)，但此时它不再能唯一确定 \(\tilde g\) 的演化/代表。
- 更精确的判断标准：
  - 不能只看 determinant；
  - 要看分量投影/相对残差：RHS 只沿 \(u_\mu u_\nu\) 方向有分量，而 LHS 的 metric 型项在横向方向也有分量；
  - 若所有这些分量的绝对值都被 \(M_P^2/\ell^2\) 等尺度压到可忽略，饱和区可作为有效退化/选解区；
  - 若要求精确解或相对残差小，则 rank-one RHS 仍不能匹配 finite-\(f\) 的满秩 metric 型 LHS。

### 2026-05-07 metric \(f(R)\) 方程中尝试令 \(f_R\tilde R_{\mu\nu}\) 逼近 \(R_{\mu\nu}\)

- 用户澄清：讨论的是 metric \(f(R)\) 方程，不是手写 \(F(\tilde R)\tilde G_{\mu\nu}=T_{\mu\nu}\)。
- 用户思路：
  - \(\tilde T_{\mu\nu}\) 拉回 \(g\) 表象后与 \(T_{\mu\nu}\) 同方向/可相等；
  - 因此希望选 \(f(\tilde R)\)，使 \(f_R\tilde R_{\mu\nu}\) 在 \(g\) 表象展开的首阶逼近 \(R_{\mu\nu}\)，从而更接近 A 支 Einstein 方程。
- 当前口径校准：
  - 下指标 \(\tilde T_{\mu\nu}=2\tilde\rho u_\mu u_\nu\)，在当前投影分支和正测度关系下通常 \(\tilde\rho\simeq\rho\)，所以拉回到 \(g\) 表象下可写成 \(2\rho u_\mu u_\nu\)；
  - 若用 \(g\) 升指标，也得到 \(2\rho u^\mu u^\nu\)；
  - 若用 \(\tilde g\) 升指标，则会带 \((m^2/X)^2\) 因子，不能混淆。
- 理论判断：
  - 设 \(\tilde R_{\mu\nu}=R_{\mu\nu}+K_{\mu\nu}[u,r,Q]\)；
  - 想要 \(f_R(\tilde R)\tilde R_{\mu\nu}\approx R_{\mu\nu}\)，需要同一个标量 \(f_R\) 同时满足所有张量分量比例；
  - 这只有在 \(K_{\mu\nu}\) 与 \(R_{\mu\nu}\) 近似同方向，或 \(K_{\mu\nu}\) 已被 \(f_R\) 压到可忽略时才可能；
  - 因为 \(f_R\) 只是 \(\tilde R\) 的标量函数，不能单独消去 \(\tilde R_{\mu\nu}\) 中的 \(ur\)、\(rr\)、横向 traceless 等分量。
- 下一步自然诊断：
  - 计算 \(\tilde R_{\mu\nu}\) 或 \(\tilde G_{\mu\nu}\) 在 \(\{u_\mu u_\nu,u_{(\mu}r_{\nu)},r_\mu r_\nu,\perp\}\) 上的投影；
  - 检查是否存在标量 \(\phi=f_R\) 使 \(\phi\tilde R_{\mu\nu}\) 的非 \(uu\) 分量足够小，同时 \(uu\) 分量匹配 \(R_{\mu\nu}\) 或 \(T_{\mu\nu}/M_P^2\)。

### 2026-05-07 metric \(f(R)\) 方向匹配诊断脚本与三切片结果

- 用户提醒：不仅要看 \(f_R\tilde R_{\mu\nu}\)，还要留意
  \[
  f_R\tilde R_{\mu\nu}-\frac12 f\tilde g_{\mu\nu},
  \]
  metric \(f(R)\) 方程左边整体方向，以及可能的 \(f\) 形式。
- 新增脚本：
  `kg_examples/diagnose_metric_fr_direction_matching.py`
- 脚本功能：
  - 复用物理 1550nm、高分辨率局部窗口、A-derived \(\tilde g\) 几何；
  - 以 \(\tilde T_{\mu\nu}/M_P^2\) 为目标张量；
  - 分别诊断三层最优局部匹配：
    1. \(\phi\tilde R_{\mu\nu}\)；
    2. \(\phi\tilde R_{\mu\nu}-\frac12 f\tilde g_{\mu\nu}\)；
    3. \(\phi\tilde R_{\mu\nu}-\frac12 f\tilde g_{\mu\nu}+\phi_R(g\Box R-\nabla\nabla R)+\phi_{RR}(g(\nabla R)^2-\nabla R\nabla R)\)；
  - 同时检查反推出的 \(\phi,f,\phi_R,\phi_{RR}\) 是否近似为 \(\tilde R\) 的单值函数。
- 运行输出：
  - `visualizations/metric_fr_direction_tau_m3p5/`
  - `visualizations/metric_fr_direction_tau0/`
  - `visualizations/metric_fr_direction_tau_p3p5/`
- trusted 区主要结果，残差定义为 `||fit-target||/||target||`：
  - `tau=-3.5`：
    - 最优 \(\phi\tilde R_{\mu\nu}\) p50/p95 `0.122/0.975`；
    - 加 \(-1/2 f g\) 后 p50/p95 `0.0386/0.648`；
    - 完整局部 jets 后 p50/p95 `0.0459/0.489`；
    - 当前 D 函数完整左边 p50/p95 `0.9997/1.124`。
  - `tau=0`：
    - 最优 \(\phi\tilde R_{\mu\nu}\) p50/p95 `0.891/0.999`；
    - 加 \(-1/2 f g\) 后 p50/p95 `0.204/0.720`；
    - 完整局部 jets 后 p50/p95 `0.154/0.700`；
    - 当前 D 函数完整左边 p50/p95 `0.9993/1.068`。
  - `tau=+3.5`：
    - 最优 \(\phi\tilde R_{\mu\nu}\) p50/p95 `0.405/0.988`；
    - 加 \(-1/2 f g\) 后 p50/p95 `0.1225/0.734`；
    - 完整局部 jets 后 p50/p95 `0.141/0.627`；
    - 当前 D 函数完整左边 p50/p95 `0.9998/1.024`。
- 初步结论：
  - 用户提醒的 \(-\frac12 f\tilde g_{\mu\nu}\) 项非常重要；它明显改善了张量方向匹配；
  - 单靠 \(f_R\tilde R_{\mu\nu}\) 在干涉中心尤其差；
  - 但反推出的 \(\phi=f_R\) 和 \(f\) 在同一 \(\tilde R\) bin 内离散度很大，尚不支持“纯 metric \(f(\tilde R)\) 单值函数足够”的结论；
  - 当前 D 函数 \(\tanh(\ell^2R)/\ell^2\) 在这三个切片上仍基本不匹配目标张量。

### 2026-05-07 metric \(f(R)\) 单值性复查

- 用户确认并追问：“单值性就是说这个 \(F\) 是不是只是 \(\tilde R\) 的函数的意思嘛？”
- 回答口径：
  - 是的。若作用量是纯 metric \(f(\tilde R)\)，则同一个 \(\tilde R\) 值不能在不同空间点要求不同的 \(f(\tilde R)\)、\(f_R(\tilde R)\)、\(f_{RR}(\tilde R)\)；
  - 否则说明局部匹配系数还依赖 \(\rho,u_\mu,r_\mu\)、位置、干涉条纹或其它曲率不变量。
- 已升级脚本 `kg_examples/diagnose_metric_fr_direction_matching.py`：
  - 新增 `core_rho_1pct` 与 `core_rho_10pct` 统计，排查低密度/节点污染；
  - 新增 invariant trace/Ricci contraction 方法反推 \(\phi,f\)，避免只依赖坐标 Frobenius 最小二乘；
  - 单值性报告新增 `median_global_mad` 与 `p95_global_mad`，避免系数接近 0 时相对 MAD 被放大。
- 新输出目录：
  - `visualizations/metric_fr_single_value_tau_m3p5/`
  - `visualizations/metric_fr_single_value_tau0/`
  - `visualizations/metric_fr_single_value_tau_p3p5/`
- 关键结果：
  - 高密度核心区确实改善单值性，因此低密度/节点会污染判断；
  - 但核心区改善后仍未达到“清楚的单值函数”：
    - `tau=0, core_rho_10pct` 最好，Frobenius algebraic \(\phi,f\) 的 relative MAD median 约 `0.65/0.48`，invariant 约 `0.48/0.45`；
    - `tau=-3.5, core_rho_10pct` 仍差，Frobenius \(\phi,f\) median relative MAD 约 `4.49/5.66`，invariant 约 `1.69/5.21`；
    - `tau=+3.5, core_rho_10pct` 也差，Frobenius 约 `9.82/4.34`，invariant 约 `1.70/3.35`。
  - invariant trace/Ricci 反推的 full tensor residual 通常明显比 Frobenius local fit 差；例如 core_rho_10pct p50 residual：
    - `tau=-3.5`: Frobenius algebraic `0.0468` vs invariant `0.823`；
    - `tau=0`: `0.1849` vs `0.2678`；
    - `tau=+3.5`: `0.1227` vs `0.936`。
- 当前结论：
  - “纯 \(f(\tilde R)\)”没有被完全否定，因为干涉中心高密度核心区显示出一定单值趋势；
  - 但三切片整体不支持它作为普适形式；
  - 更可能需要依赖额外不变量，例如 \(\tilde R_{\mu\nu}\tilde R^{\mu\nu}\)、\(\tilde T\) 相关标量、或由 \(u,r\) 诱导的几何标量。

### 2026-05-07 保留 metric \(f(R)\) 二阶导数项后的普适函数拟合

- 用户追问：考虑上 metric 方程左边的二阶导数项后，再看是否可以是普适函数。
- 新增脚本：
  `kg_examples/fit_metric_fr_universal_function.py`
- 脚本做法：
  - 设单一普适函数
    \[
    f(R)=\sum_n c_n T_n(t(\operatorname{asinh}(R/R_0)))
    \]
    并自动计算同一个函数的 \(f_R,f_{RR},f_{RRR}\)；
  - 回代完整 metric \(f(R)\) 局部结构：
    \[
    f_R R_{\mu\nu}-\frac12 f g_{\mu\nu}
    +f_{RR}(g_{\mu\nu}\Box R-\nabla_\mu\nabla_\nu R)
    +f_{RRR}(g_{\mu\nu}(\nabla R)^2-\nabla_\mu R\nabla_\nu R) ;
    \]
  - 在 `tau=-3.5,0,+3.5` 三个切片的 `core10` 上同时拟合同一个 \(f(\tilde R)\)。
- 第一版未归一化最小二乘因列尺度跨越过大，不能直接作为最终判断；随后已加入列归一化。
- 归一化输出目录：
  `visualizations/metric_fr_universal_core10_scaled/`
- 关键结果：
  - degree 6/10/14/18 的加权相对残差分别约：
    `0.99993, 0.99991, 0.99984, 0.99983`；
  - 三个切片 core10 回代 `residual_to_rhs` p50 均约 `0.99997~0.99999`；
  - 回代 LHS/RHS p50 仅约 `3e-5~5e-5`，说明拟合出的普适 \(f(R)\) 几何侧基本没有追上物质目标；
  - 提高 degree 没有实质改善，且条件数快速恶化。
- 当前结论：
  - 这一步已经把 metric \(f(R)\) 的二阶导数项纳入，并强制所有导数来自同一个 \(f(\tilde R)\)；
  - 在当前 1550nm 高斯干涉三切片、高密度核心区上，普适纯 \(f(\tilde R)\) 不能重现逐点局部拟合的低残差；
  - 因此此前“逐点可拟合”主要来自每点独立选择 \(\phi,f,\phi_R,\phi_{RR}\) 的自由度，不是来自一个普适纯 \(f(\tilde R)\)。

### 2026-05-07 Planck 质量与 2+1d 截面量纲 caveat

- 用户怀疑：普适 \(f(\tilde R)\) 拟合中 LHS/RHS 量级追不上，是否丢了 \(\hbar,c\) 或 \(M_P\)。
- 核查结果：
  - `physical_units.py` 中 `HBAR_C_EV_M=1.973269804e-7`，1550nm 被换算为 `k0≈0.799898 eV`，这条 \(\hbar c\) 换算链条看起来正确；
  - 代码当前 `PLANCK_MASS_EV=1.220890128e28 eV` 是非约化 Planck mass；
  - 若场方程写作 \(M_P^2G_{\mu\nu}=T_{\mu\nu}\)，则通常应使用约化 Planck mass \(M_{P,\mathrm{red}}\approx2.4353e27 eV\)，差一个 \(8\pi\)；
  - 当前高斯算例是 2+1d 截面 KG 归一化，若解释为 3+1d 光束截面，需要给物质张量补横向 profile/厚度因子 \(1/L_y\)。
- 数值尺度：
  - `sigma_perp≈9.001 eV^-1≈1.776um`；
  - 用 \(L_y=\sigma_\perp\) 或 \(\sqrt{2\pi}\sigma_\perp\) 估计，横向厚度与约化 Planck mass 的组合修正大约是 O(1) 到 O(3)，不是 \(10^4\) 或 \(10^5\) 级。
- 对当前普适 \(f(\tilde R)\) 结论的影响：
  - 全局常数重标定 RHS 会让最佳 \(f\) 系数整体重标定，但不会改变“同一个 \(f(\tilde R)\) 是否能同时拟合三切片”的相对残差；
  - 当前 universal fit 的问题是列空间/函数自由度无法覆盖目标张量方向，表现为相对残差仍约 1；
  - 因此该 caveat 必须修正/显式报告，但目前不足以拯救纯 \(f(\tilde R)\)。

### 2026-05-07 \(f(\tilde R,\tilde R_{\mu\nu}\tilde R^{\mu\nu})\) 第一轮可行性检验

- 用户要求考虑
  \[
  \sqrt{-\tilde g}\, f(\tilde R,\tilde R_{\mu\nu}\tilde R^{\mu\nu})
  \]
  是否可行。
- 新增脚本：
  - `kg_examples/diagnose_metric_fr_ricci2_direction_matching.py`
  - `kg_examples/fit_metric_fr_ricci2_universal_algebraic.py`
- 记号：
  - \(I_2=\tilde R_{\mu\nu}\tilde R^{\mu\nu}\)；
  - \(Q_{\mu\nu}=\tilde R_{\mu\alpha}\tilde R^\alpha{}_\nu\)；
  - metric \(f(R,I_2)\) 场方程的代数主方向为
    \[
    f_R \tilde R_{\mu\nu}+2 f_{I_2} Q_{\mu\nu}-\frac12 f\tilde g_{\mu\nu}.
    \]
- 第一脚本做局部逐点诊断：
  - 同一物理 1550nm KG 高斯干涉局部窗口；
  - 三切片 `tau=-3.5,0,+3.5`；
  - 比较纯 \(f(R)\) 局部 \(\phi R_{\mu\nu}-\frac12 fg_{\mu\nu}\) 与加入 \(I_2\) 后局部 \(\phi R_{\mu\nu}+2\psi Q_{\mu\nu}-\frac12 fg_{\mu\nu}\)；
  - 同时做 \((R,I_2)\) 二维单值性 bin 诊断和 invariant trace/Ricci/Q contraction 反推。
- 高分辨率输出：
  - `visualizations/metric_fr_ricci2_direction_tau_m3p5/`
  - `visualizations/metric_fr_ricci2_direction_tau0/`
  - `visualizations/metric_fr_ricci2_direction_tau_p3p5/`
- 局部最佳残差 `||fit-target||/||target||` 结果：
  - `tau=-3.5` trusted：纯 \(R\) p50/p95 `0.0386/0.648`，\(R,I_2\) p50/p95 `0.0284/0.391`；
  - `tau=0` trusted：纯 \(R\) `0.204/0.720`，\(R,I_2\) `0.0980/0.546`；
  - `tau=+3.5` trusted：纯 \(R\) `0.1225/0.734`，\(R,I_2\) `0.0543/0.463`；
  - `core_rho_10pct` 也类似改善，例如 `tau=0` 从 `0.1849/0.540` 到 `0.0837/0.400`。
- 重要 caveat：
  - invariant contraction 反推并不乐观：三切片 trusted p50 通常 `0.39~0.57`，p95 `~3`；
  - 二维单值性仍差，局部 \(\phi,\psi,f\) 的相对 MAD 常为 O(1)-O(10) 甚至更大；
  - 因此“每点多一个张量方向就能拟合得更好”不等于“存在普适 \(f(R,I_2)\)”。
- 第二脚本做更硬的普适二维函数代数回代：
  - \(f(R,I_2)=\sum c_{nm}T_n(t_R)T_m(t_{I_2})\)；
  - 强制 \(f,f_R,f_{I_2}\) 来自同一个函数；
  - 只回代 metric \(f(R,I_2)\) 的代数部分，尚未包含 \(\nabla\nabla f_R\) 与 \(\nabla\nabla(f_{I_2}R_{\mu\nu})\) 等导数项。
- 高分辨率输出：
  - `visualizations/metric_fr_ricci2_universal_algebraic_core10/`
- 普适代数拟合结果：
  - degree `4:4,6:6,8:8,10:10` 的 weighted relative residual 约 `0.920,0.908,0.906,0.903`；
  - 比纯 \(f(R)\) universal 的 `~0.9998` 稍好，但仍远未通过；
  - degree 提高后条件数从 `1.6e4` 到 `5.1e8`，改善很慢且病态增长；
  - core10 三切片回代 p50 仍大约 `0.78~0.95`。
- 当前结论：
  - 加入 \(I_2\) 在“局部张量方向”上明显有帮助；
  - 但纯粹的普适代数 \(f(R,I_2)\) 还不够；
  - 真正的下一步若继续此路线，必须推导/实现完整 metric \(f(R,I_2)\) 变分的导数项，尤其是 \(\nabla\nabla f_R\) 和 \(\nabla\nabla(f_{I_2}R_{\mu\nu})\) 的 contribution。

### 2026-05-07 完整 metric \(f(R,I_2)\) 导数项回代检验

- 用户要求：把完整 metric 方程中的导数项也加上考虑。
- 文献/公式口径：
  - 对 \(Y=R_{\mu\nu}R^{\mu\nu}\)，metric \(f(R,Y)\) 方程左边采用
    \[
    f_R R_{\mu\nu}-\frac12 f g_{\mu\nu}
    +(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)f_R
    +2f_Y R_{\mu\alpha}R^\alpha{}_\nu
    +\Box(f_YR_{\mu\nu})
    +g_{\mu\nu}\nabla^\alpha\nabla^\beta(f_YR_{\alpha\beta})
    -2\nabla_\alpha\nabla_\beta(f_YR^\alpha{}_{(\mu}\delta^\beta{}_{\nu)}).
    \]
  - 该结构与已查文献一致；注意某些网页排版中 \(f_R,f_Y\) 定义有笔误，但方程结构清楚。
- 新增脚本：
  - `kg_examples/fit_metric_fr_ricci2_universal_full.py`
- 脚本做法：
  - 固定 A-derived \(\tilde g[\rho,S]\) 背景；
  - 使用五个时间层构造中央切片、前后切片的 \(R,Y,R_{\mu\nu}\)；
  - 对每个 Chebyshev 二维基函数显式计算完整 \(f(R,Y)\) 左边；
  - 强制同一个 \(f(R,Y)\) 的 \(f,f_R,f_Y\) 参与所有项；
  - 进行三切片共同最小二乘回代；
  - 这是 fixed-background 回代/拟合诊断，不是完整 D 支动力学求解器。
- 数值输出：
  - 三切片中分辨率：
    `visualizations/metric_fr_ricci2_universal_full_core10_n160/`
  - 低分辨率趋势：
    `visualizations/metric_fr_ricci2_universal_full_core10_n96/`
  - 单切片高阶过拟合检查：
    `visualizations/metric_fr_ricci2_universal_full_tau0_n160_deg4/`
- 关键结果：
  - 三切片 `n=160, core10`：
    - degree `1:1,2:2,3:3` 的 weighted relative residual 分别约 `0.9985,0.9959,0.9942`；
    - 三个切片 core10 p50 残差基本都在 `~0.999`；
    - 导数项非零，但没有把 LHS 拉到 RHS 方向/量级上。
  - 低分辨率 `n=96` 有很小改善，但 `n=160` 消失，故不能把低分辨率改善当成物理信号。
  - 单切片 `tau=0,n=160` 若只拟合自己：
    - weighted residual 可到 `0.763`；
    - core10 p50 可到约 `0.69~0.71`；
    - 但 trusted p95 爆到 `~90~150`，条件数从 `32` 增长到 `2.1e7`；
    - 说明这更像不稳定局部过拟合，而不是普适作用量。
- 当前结论：
  - 在当前 1550nm 高斯干涉 fixed-background 回代口径下，完整 metric \(f(R,I_2)\) 普适函数不支持三切片；
  - 这比“代数项不足”的结论更强；
  - 若继续几何作用量路线，应考虑加入更多不变量/方向，或改变目标口径，而不是继续单纯提高 \(f(R,I_2)\) 阶数。

### 2026-05-07 完整 \(f(R,I_2)\) 更高阶检查

- 用户要求：对更高阶数再检查。
- 已完成三切片 `n=160, core10` 完整 metric 方程 `degree=4:4`：
  - 输出目录：
    `visualizations/metric_fr_ricci2_universal_full_core10_n160_deg4/`
  - weighted residual `0.99088`；
  - 条件数约 `4.49e3`；
  - core10 p50 residual：
    - `tau=-3.5`: `0.999656`
    - `tau=0`: `0.998826`
    - `tau=+3.5`: `0.999850`
  - core10 p95 residual：
    - `tau=-3.5`: `1.174816`
    - `tau=0`: `0.999919`
    - `tau=+3.5`: `1.107145`
- 与 `degree=3:3` 对比：
  - weighted residual 从 `0.99422` 降到 `0.99088`；
  - 但三切片 core p50 仍几乎等于 1；
  - 因此没有实质改善。
- 曾启动 `degree=4:4,5:5` 联合扫描，但 `5:5` 在当前 Python 原型下超过 7 分钟仍未完成；
  - 由于 `4:4` 已显示无实质改善，停止该长进程以免继续占用 CPU；
  - 若后续必须检查 `5:5` 以上，应先向量化 `fit_metric_fr_ricci2_universal_full.py` 中张量二阶协变导数计算。

### 2026-05-07 单变量 \(f(\tilde R)\) 高阶复查

- 用户澄清：要求检查单变量 \(f(\tilde R)\) 的更高阶，而不是 \(f(R,I_2)\)。
- 已重新运行：
  - `kg_examples/fit_metric_fr_universal_function.py`
  - 输出目录：
    `visualizations/metric_fr_universal_core10_highdeg_22_34/`
  - 参数：三切片 `tau=-3.5,0,+3.5`，`core10` 共同拟合，完整 metric \(f(R)\) 导数项已包含。
- 与旧结果对比：
  - 旧 degree `18` weighted residual `0.999826568`；
  - 新 degree `22` weighted residual `0.999730961`；
  - 新 degree `26` weighted residual `0.999599407`；
  - 新 degree `30` weighted residual `0.999617940`；
  - 新 degree `34` weighted residual `0.999453948`。
- 三切片 core10 残差仍几乎等于 1：
  - degree `34`：
    - `tau=-3.5`: p50/p95 `0.999980/1.00010`
    - `tau=0`: p50/p95 `0.999968/0.999994`
    - `tau=+3.5`: p50/p95 `0.999994/1.00015`
- 条件数已到 \(10^{14}\) 量级：
  - degree `22`: `2.6e14`
  - degree `26`: `2.8e14`
  - degree `30`: `3.0e14`
  - degree `34`: `2.26e14`
- LHS/RHS 典型量级仍很小：
  - degree `34` 下 core10 的 LHS/RHS p50 约 `3e-5~4e-5`；
  - 说明问题不是少数点坏，而是普适 \(f(R)\) 几何侧整体没有追上 \(\tilde T_{\mu\nu}/M_P^2\) 的量级/方向。
- 当前结论：
  - 单变量 \(f(\tilde R)\) 高阶到 degree `34` 仍未改善；
  - 更高阶还可能继续病态化，不建议作为主线。

### 2026-05-07 \(\mathcal R_{\mu\nu}\) 是否 pure-\(\tilde g\) 的高斯干涉诊断

- 用户指出记号问题：不能再用大写 `G[...]` 表示反解度规，因为这会和 Einstein tensor \(G_{\mu\nu}\) 混淆。
- 已修正 `research-notes/145-从g表象Einstein方程变换得到gtilde方程的分解路线.md`：
  - 反解度规记号统一改为 \(\mathfrak g_{\mu\nu}[\tilde g,u,r,\rho,\ldots]\)；
  - 仅保留一处说明“避免用大写 \(G[\cdots]\)”。
- 新增脚本：
  - `kg_examples/diagnose_mathcal_r_pure_geometry.py`
  - 它在 1550nm 高斯干涉三切片上计算
    \[
    \mathcal R_{\rm need}/M_P^2=\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2
    \]
    并用已修正的密度变换
    \(\sqrt{|\tilde g|}\tilde\rho=|X|\rho_A/m^2\) 构造 \(\tilde T_{\mu\nu}\)。
- 验证：
  - `python3 -m py_compile kg_examples/diagnose_mathcal_r_pure_geometry.py` 通过。
- 输出：
  - smoke: `visualizations/mathcal_r_pure_geometry_smoke_n96/`
  - 三切片中分辨率: `visualizations/mathcal_r_pure_geometry_n160_3tau/`
  - 三切片高分辨率: `visualizations/mathcal_r_pure_geometry_n320_3tau/`
- 高分辨率 `n=320,crop=97x97,taus=-3.5,0,+3.5,core10` 关键结果：
  - 完整 \(\mathcal R_{\rm need}/M_P^2\) 用 pure geometry 基底拟合几乎精确；
    - 不显式放 \(\tilde G\)，只用 \(\tilde R_{\mu\nu}\) 与 \(-\frac12\tilde R\tilde g_{\mu\nu}\) 组合也得到 weighted residual `3.73e-12`；
    - 这是代数恒等式层面的平凡结果，因为目标本来几乎就是 \(\tilde G_{\mu\nu}\)。
  - 物质源 \(\tilde T_{\mu\nu}/M_P^2\) 的尺度被真实 \(M_P\) 压到约 `1e-60`；
    - `|Gtilde|` p50 在三切片为约 `169,418,66`；
    - `|Ttilde|/Mp^2` p50 在三切片为约 `1.22e-60,1.05e-60,1.74e-60`。
  - 非平凡判据：拟合 \(\tilde T_{\mu\nu}/M_P^2\) 时，pure geometry without EH 的 weighted residual `0.996`，三切片 p50 残差约 `0.99`；
  - 加入简单 matter directions \(\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\) 后，weighted residual 降到 `0.594`，但仍不是闭合。
- 当前判断：
  - 不能把“\(\mathcal R_{\rm need}\) 可 pure geometry 拟合”解读为 pure-\(\tilde g\) 作用量路线被救活；
  - 它主要反映 \(\tilde G_{\mu\nu}\) 的平凡几何主项和 \(M_P\) 对物质源的巨大压制；
  - 若要构造严格动力学，仍需研究非平凡 \(u,r,\rho\)-dependent sector，或者从 action-level 重新变分构造。

### 2026-05-07 \(\mathcal R\) 是否是 \(\tilde g\) 函数/泛函的留一切片诊断

- 用户指出：如果 \(\mathcal R_{\mu\nu}\) 是纯几何量，则 pure-\(\tilde g\) 路线仍有希望；要求查看 \(\mathcal R\) 如何由 \(\tilde g\) 产生，是否是 \(\tilde g\) 的函数。
- 已新增脚本：
  - `kg_examples/diagnose_mathcal_r_gtilde_function.py`
  - 方法：留一时间切片 kNN 预测，用两个 \(\tau\) 切片的局部 \(\tilde g\) 几何特征预测第三个切片的非平凡源 \(\tilde T_{\mu\nu}/M_P^2\)；
  - 不把完整 \(\mathcal R_{\rm need}\) 当目标，因为完整目标会被 \(\tilde G_{\mu\nu}\) 平凡主项支配。
- 验证：
  - `python3 -m py_compile kg_examples/diagnose_mathcal_r_gtilde_function.py` 通过。
- 输出：
  - `visualizations/mathcal_r_gtilde_function_n160_3tau/`
  - `visualizations/mathcal_r_gtilde_function_n320_3tau/`
- 高分辨率 `n=320,crop=97x97,core10,k=16` 结果：
  - curvature scalars 特征（\(\tilde R,\mathrm{tr}\tilde R^2,\mathrm{tr}\tilde R^3,\log|\det\tilde g|\)）预测 \(\tilde T/M_P^2\)：overall p50 `0.747`，weighted mean `0.924`；
  - local geometry tensors 特征（\(\tilde g,\tilde R_{\mu\nu},\tilde G_{\mu\nu},\tilde R^2_{\mu\nu},I_2\tilde g\) 等）：overall p50 `0.598`，weighted mean `0.863`；
  - 加入 \(u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\) 后没有稳定改善：overall p50 `0.635`，weighted mean `0.997`。
- 当前判断：
  - 完整 \(\mathcal R_{\rm need}\) 的 pure-\(\tilde g\) 性质确实说明“还有希望”，但这主要是 EH 主项的希望；
  - 非平凡细节不支持“简单局部 \(\tilde g\) 函数”假设；
  - 仍未排除高阶、非局域或 action-level 特殊 pure-\(\tilde g\) 泛函；
  - 下一步若继续 pure-\(\tilde g\) 路线，应检查 Helmholtz/self-adjoint integrability 与 \(\tilde\nabla\)-守恒，而不是仅继续低阶曲率拟合。

### 2026-05-07 更高阶局域 pure-\(\tilde g\) 路线

- 用户要求先考虑“更高阶的方程（纯 `gtilde` 衍生，局域）”。
- 已在 `research-notes/145-从g表象Einstein方程变换得到gtilde方程的分解路线.md` 中补入新节：
  - 明确区分高次 `f(tilde R)` 与真正的有限 jet 局域度规泛函；
  - 给出自然的局域高阶基底：\(\tilde R^2\)、\(\tilde R_{\mu\nu}\tilde R^{\mu\nu}\)、\(\tilde\nabla\tilde R\cdot\tilde\nabla\tilde R\)、\(\tilde R\tilde\Box\tilde R\)、\(\tilde\nabla\tilde R_{\mu\nu}\cdot\tilde\nabla\tilde R^{\mu\nu}\) 等；
  - 说明其方程是四阶或更高阶局域 PDE，需检查主部、适定性、Gauss-Bonnet/曲率恒等式与 2+1d 约化下的独立不变量。
- 当前策略更新：
  - 不再把“继续提高单变量 `f(tilde R)` 次数”当作主路；
  - 若继续 pure-\(\tilde g\) 局域路线，应该先做 quadratic curvature basis，再加导数不变量，最后才检查 action-level integrability。

### 2026-05-07 4 阶 quadratic local scan 直接结论

- 为回答用户“我们的例子下，大概要到多少阶才可能拟合得比较好”，已对标准 4 阶局域 quadratic action 做直接扫描：
  \[
  S_{\rm grav}^{(4)}=\int\sqrt{|\tilde g|}\left[\Lambda+\frac{M_P^2}{2}\tilde R+\alpha \tilde R^2+\beta \tilde R_{\mu\nu}\tilde R^{\mu\nu}\right].
  \]
- 扫描方式：
  - 用其 metric variation 张量基底去拟合 1550nm 高斯干涉三切片上的非平凡源 \(\tilde T_{\mu\nu}/M_P^2\)；
  - 使用 `n=160`、`core10`、\(\tau=-3.5,0,3.5\)。
- 结果：
  - 单片加权残差：
    - \(\tau=-3.5\): `0.993246`
    - \(\tau=0\): `0.773260`
    - \(\tau=3.5\): `0.997984`
  - 三切片总体加权残差：`0.998489`
- 结论：
  - 2 阶局域 pure-\(\tilde g\) 基底不够；
  - 4 阶 quadratic local action 也还是不够，而且整体几乎没有把问题解决掉；
  - 若继续坚持 pure-\(\tilde g\) 且局域，至少应考虑 6 阶及以上的导数不变量，或者承认仅靠局域 pure-\(\tilde g\) 仍不足以闭合当前高斯干涉的非平凡剩余。
- 2026-05-07 6 阶纯 \(f(\tilde R)\) 复查：
  - 用同一 1550nm 三切片、`n=320`、`core10` 做最轻量 6 阶 pure-\(f(\tilde R)\) sanity check；
  - 总体加权残差为 `0.998034319760273`，与 4 阶 quadratic local action 几乎同级；
  - 单片 core10 中位相对残差仍约 `0.994~0.999`；
  - 这确认了“单纯继续抬纯 \(f(\tilde R)\) 的多项式次数”并不是有效方向。
- 2026-05-07 真 6 阶 finite-jet local 原型：
  - 新建 `kg_examples/fit_metric_finitejet6_local.py`，采用最小局域基底 \(\{1,\tilde R,I_2,\tilde R^2,\tilde R I_2,\tilde R^3\}\)；
  - 在 `n=160`、`core10`、三切片 \(\tau=-3.5,0,3.5\) 上，三切片总体加权残差仍为 `0.9986743707235821`；
  - 单片 core10 中位残差仍约 `0.991~0.999`；
  - 说明真 6 阶局域 finite-jet action 的最小子空间也没有显著改善当前高斯干涉的非平凡源。
- 2026-05-07 方程-first 路线更新：
  - 明确区分“方程可解/可闭合”与“方程可由局域作用量导出”两件事；
  - 当前可先采用 equation-first 路线：先写协变、对称、\(\tilde\nabla\)-守恒、平直极限正确的 pure-\(\tilde g\)+matter 闭合方程；
  - 作用量重建留到后验，作为 Helmholtz/self-adjoint integrability 的附加筛选，而不是第一道门槛。

### 2026-05-07 equation-first 最小张量方向与测地线限制

- 用户要求先找符合要求的最小方程，同时留意测地线限制。
- 新建并运行 `kg_examples/fit_equation_first_tensor_couplings.py`：
  - 方程形式：
    \[
    \tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2=\mathcal C_{\mu\nu}.
    \]
  - 在 `n=160, core10, tau=-3.5,0,3.5` 上比较张量方向。
- 结果：
  - 只用 `uu`：overall p50 residual `0.396`，weighted mean `0.466`，不够；
  - `g+uu`：p50 `0.189`，weighted mean `0.265`；
  - `uu+rr`：p50 `0.223`，weighted mean `0.296`；
  - `uu+rr+ur`：p50 `0.01797`，weighted mean `0.0723`；
  - `g+uu+rr+ur`：p50 `0.00520`，weighted mean `0.0343`。
- 当前判断：
  - 最小非平凡候选可先取
    \[
    \mathcal C_{\mu\nu}=Buu+Crr+D\,u_{(\mu}r_{\nu)}
    \]
    因为它已显著压低 residual，且不含最容易产生横向力的 \(A\tilde g_{\mu\nu}\) 项；
  - 若允许 \(A\tilde g_{\mu\nu}\)，拟合更好，但必须同时满足
    \(h^\nu{}_\alpha \tilde\nabla^\mu\mathcal C_\mu{}^\alpha=0\)，让 \(h\nabla A\) 被其它项抵消；
  - 下一步应做 constrained fit：最小化 residual 同时惩罚/约束 transverse force，而不是只做张量拟合。

### 2026-05-07 最小 equation-first 候选的横向力诊断

- 新建并运行 `kg_examples/diagnose_equation_first_transverse_force.py`。
- 目的：比较各候选 \(\mathcal C_{\mu\nu}\) 的 algebraic residual 与
  \(h^\nu{}_\alpha\tilde\nabla^\mu\mathcal C_\mu{}^\alpha\) 横向力。
- 结果（`n=160, core10`，三切片）：
  - `uu`：algebraic residual p50 约 `0.44`，但 transverse fraction p50 约 `2.5e-5`，即测地线友好但拟合差；
  - `uu+rr+ur`：algebraic residual p50 约 `0.0187`，但 transverse fraction p50 约 `1.47`，即拟合好但无约束逐点系数会破坏测地线；
  - `g+uu+rr+ur`：algebraic residual p50 约 `0.00466`，但 transverse fraction p50 约 `12.5`，拟合更好但测地线风险更大。
- 当前判断：
  - \(A\tilde g_{\mu\nu}\) 不适合作为最小模型首选；
  - \(Buu+Crr+D\,ur\) 是最小张量结构，但必须把横向力条件作为闭合方程；
  - 更干净的下一版可先试强条件 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\)，它自动保证不向 \(\tilde T_{\mu\nu}\) 注入力。

### 2026-05-07 constrained equation-first 拟合：BCD 软约束、BCD 硬约束、gBCD 硬约束

- 新建并升级 `kg_examples/fit_equation_first_constrained_bcd.py`。
- 脚本现在支持：
  - 软约束：最小化代数残差并惩罚 \(h^\nu{}_\alpha\tilde\nabla^\mu\mathcal C_\mu{}^\alpha\)；
  - 硬约束：通过 KKT 系统把横向力行作为线性等式；
  - 可选 atoms：`uu,rr,ur` 或 `g,uu,rr,ur`。
- BCD 软约束结果：
  - `n=160, tau=0, fit_region=trusted, force_region=trusted` 输出：
    `visualizations/equation_first_constrained_bcd_n160_tau0_trustedfit/`；
  - 最佳软权重 `lambda=100` 下，中心切片代数残差仍小：
    weighted mean `0.01493`，p50 `0.00887`，p95 `0.04634`；
  - 横向力自然尺度残差
    \(|J_\perp|/(|C|/h)\) weighted mean `2.16e-4`，p50 `1.10e-5`，p95 `8.71e-4`；
  - 但相邻时间层代数残差升高到 weighted mean 约 `0.037` 和 `0.031`，说明系数场的时间闭合不是免费的。
- BCD 硬约束结果：
  - `n=96, tau=0, trusted` 输出：
    `visualizations/equation_first_constrained_bcd_hard_n96_tau0_trustedfit/`；
  - 硬约束把 force rows 压到 `norm_Fx=8.15e-8`；
  - 但中心切片代数残差灾难性升高：weighted mean `13.85`，p50 `0.344`，p95 `64.10`；
  - 结论：严格横向力为零时，仅 \(Buu+Crr+D\,ur\) 在 trusted 区域无法同时保持场方程近似成立。
- gBCD 硬约束结果：
  - 使用 `--atoms g,uu,rr,ur --hard-constraint`；
  - `n=96, tau=0, trusted` 输出：
    `visualizations/equation_first_constrained_gbcd_hard_n96_tau0_trustedfit_cli/`；
  - force residual 几乎为零：`force_residual_norm=1.21e-15`，
    \(|J_\perp|/(|C|/h)\) weighted mean `9.78e-15`；
  - 中心切片代数残差仍好：weighted mean `0.01177`，p50 `0.00549`，p95 `0.03887`；
  - 相邻时间层代数残差约 `0.089` 与 `0.082`，提示仍需真正的时间演化/闭合方程，而不是仅三层投影。
- 当前判断：
  - 三项 BCD 是好的软投影结构，但不是严格闭合的最小方程；
  - 加入 \(A\tilde g_{\mu\nu}\) 后，gBCD 是当前最小可继续研究的 equation-first 候选；
  - 下一步应围绕 gBCD 推导系数方程/守恒闭合，而不是回到纯 \(f(\tilde R)\) 或低阶 pure-\(\tilde g\) action。

### 2026-05-07 gBCD 三时刻 hard constraint、full 守恒与本构关系初检

- 已把 `kg_examples/fit_equation_first_constrained_bcd.py` 升级：
  - 新增 `--force-mode transverse/full`；
  - `transverse` 强制 \(h^\nu{}_\alpha\tilde\nabla^\mu\mathcal C_\mu{}^\alpha=0\)；
  - `full` 强制 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\)。
- 已新建研究笔记：
  - `research-notes/146-equation-first-gBCD闭合方程与守恒检验.md`
  - 写出 gBCD 守恒方程展开式，并说明它等价于各向异性有效应力张量的守恒闭合问题。
- transverse hard 三时刻检验：
  - 输出 `visualizations/equation_first_constrained_gbcd_hard_n96_3tau_summary/`；
  - 中心切片 weighted mean 代数残差：
    \(\tau=-3.5\): `0.02858`，
    \(\tau=0\): `0.01177`，
    \(\tau=+3.5\): `0.04607`；
  - 横向力自然尺度约束压到 \(10^{-15}\sim10^{-11}\)。
- full divergence-free hard 三时刻检验：
  - 输出 `visualizations/equation_first_constrained_gbcd_hard_full_n96_3tau_summary/`；
  - 中心切片 weighted mean 代数残差：
    \(\tau=-3.5\): `0.02859`，
    \(\tau=0\): `0.01179`，
    \(\tau=+3.5\): `0.04641`；
  - full 守恒也能被 KKT 约束压到接近数值零；
  - 相邻 probe 时间层残差升高到 `0.14~0.23`，比 transverse-only 更强地暴露时间闭合问题。
- 新建 `kg_examples/diagnose_gbcd_constitutive_coefficients.py`：
  - 用 full hard 三时刻的 \(A,B,C,D\) 系数场，检查是否可由局部标量函数预测；
  - 特征集：\(\log\rho,\tilde R,I_2,u^2,r^2,u\cdot r\)；
  - 输出 `visualizations/equation_first_gbcd_constitutive_full_n96_3tau/`。
- 本构关系初检结果：
  - 简单线性/二次多项式的 in-sample weighted relative error 通常已经 \(\gg1\)；
  - 留一时刻泛化更差，尤其 \(\tau=-3.5\)；
  - kNN 也不能稳定跨时刻预测；
  - 结论：不能把 \(A,B,C,D\) 轻易当作这些简单局部标量的无新自由度函数。需要更复杂本构、辅助场演化，或处理 KKT null-space/gauge 后重测。

### 2026-05-07 gBCD 系数产生机制：ridge 代表元、单势函数与辅助应力场

- 用户要求“探究 \(A,B,C,D\) 的产生机制”。
- 先检查 KKT 代表元问题：
  - 对 `tau=0, full hard, gBCD` 扫描 ridge/minimum-norm 正则 `1e-12` 到 `1e-2`；
  - `1e-8~1e-6` 基本不破坏中心代数残差，仍约 `1.18%`；
  - 部分系数极值下降，说明本构失败中确实有 KKT null-space/ill-conditioning 代表元选择的成分；
  - 用 `ridge=1e-6` 补跑 \(\tau=\pm3.5\)，输出：
    `visualizations/equation_first_gbcd_hard_full_n96_taum3p5_ridge_1em6/`
    和
    `visualizations/equation_first_gbcd_hard_full_n96_taup3p5_ridge_1em6/`。
- 用 ridge=1e-6 重新做本构函数诊断：
  - 输出 `visualizations/equation_first_gbcd_constitutive_full_n96_3tau_ridge_1em6/`；
  - 相比未正则解，`uu` 与 `ur` 的 kNN 跨时刻误差明显下降；
  - 但 `rr` 方向仍很差，简单局部函数机制仍未通过。
- 提出并检验“单个局域势函数”机制：
  - 若 \(L=L(\rho,U,V,W)\)，\(U=u^2,V=r^2,W=u\cdot r\)，则 metric variation 产生
    \(L\tilde g_{\mu\nu}-2L_Uuu-2L_Vrr-2L_Wur\)，形式上正是 gBCD；
  - 新建 `kg_examples/fit_gbcd_potential_constitutive.py`；
  - degree 2 输出 `visualizations/equation_first_gbcd_potential_matter_deg2_n96_3tau/`，最佳 global weighted residual `0.98984`；
  - degree 3 输出 `visualizations/equation_first_gbcd_potential_matter_deg3_n96_3tau/`，最佳 global weighted residual `0.98085`；
  - 结论：低阶 \(L(\rho,u^2,r^2,u\cdot r)\) 单势函数机制失败。
- 新建研究笔记：
  - `research-notes/147-gBCD系数ABCD的可能产生机制.md`。
- 当前判断：
  - 最合理机制是“辅助各向异性应力场”：
    \[
    \mathcal C_{\mu\nu}=\lambda_I E^I_{\mu\nu},\quad
    \lambda_I=(A,B,C,D),
    \]
    其中 \(\lambda_I\) 由场方程投影、full 守恒 PDE、边界/初值条件和一个代表元规范共同产生；
  - 这不是无新自由度的简单局部函数，也不是低阶单势函数直接变分。

### 2026-05-07 gBCD 代表元规范与 nullspace hard constraint

- 新建并运行 `kg_examples/fit_gbcd_auxiliary_gauge.py`：
  - 在同一 gBCD/full-conservation 方程下比较代表元规范；
  - 支持 `norm`、`space`、`time`、`space_time`；
  - 新增 `--hard-solver nullspace`，先求约束矩阵零空间，再在零空间内最小化代数残差和平滑泛函。
- 发现旧 KKT hard constraint 在病态切片上可能数值泄漏：
  - tau=±3.5 的原始 divergence 异常不是物理失败，而是 KKT 条件数导致；
  - nullspace 版把缩放约束残差压到 `~1e-16`。
- tau=0 空间平滑扫描：
  - `space_weight=1e-6` 会把中心代数 weighted mean 从约 `0.0119` 恶化到约 `0.383`；
  - `space_weight=1e-7` 仍恶化到约 `0.072`；
  - `1e-9~1e-10` 基本不损伤残差，但也没有稳定空间粗糙度；
  - 因此空间平滑暂不能当作物理状态方程。
- nullspace 正式三切片输出：
  - `visualizations/equation_first_gbcd_auxiliary_gauge_nullspace_n96_taum3p5_norm_time/`
  - `visualizations/equation_first_gbcd_auxiliary_gauge_nullspace_n96_tau0_norm_time/`
  - `visualizations/equation_first_gbcd_auxiliary_gauge_nullspace_n96_taup3p5_norm_time/`
- 正式结果：
  - norm 代表元中心代数 weighted mean：
    tau=-3.5 `0.02946`，tau=0 `0.01187`，tau=+3.5 `0.04647`；
  - time 代表元中心代数 weighted mean：
    tau=-3.5 `0.03252`，tau=0 `0.01199`，tau=+3.5 `0.04692`；
  - hard constraint max 均为 `~1e-16`；
  - time 代表元可温和降低时间层间系数跳动，但不是最终物理本构。
- 新建研究笔记：
  - `research-notes/148-gBCD辅助应力场代表元规范与nullspace硬约束.md`。
- 新增并运行 `kg_examples/diagnose_gbcd_conservation_principal_symbol.py`：
  - 检查 \(\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0\) 对 \(\lambda_I=(A,B,C,D)\) 的一阶主符号；
  - 输出 `visualizations/equation_first_gbcd_conservation_principal_symbol_n96_trusted/`；
  - 在三切片、`lab_t/lab_x/lab_z` 方向上，主符号 rank 全部为 `3`，nullity 全部为 `1`；
  - 结论：在当前 2+1 约化中，full conservation 是必要条件，但不能唯一演化四个 \(\lambda_I\)，必须额外给一个状态方程/规范条件。

### 2026-05-07 gBCD 状态方程：trace 与剪切闭合首轮

- 用户提出 \(Q\to0\) 时应回到 Einstein 方程；已判断这是必要边界条件，但不自动给出有限 \(Q\) 的闭合。
- 新建并运行 `kg_examples/fit_gbcd_trace_closure.py`：
  - 在 gBCD/full-conservation/nullspace hard constraint 上额外加入代数状态方程；
  - 候选包括 `zero`、`q1`、`q2`、`shear0`。
- trace closure 定义：
  - `zero`: \(\mathcal C^\mu{}_\mu=0\)；
  - `q1`: \(\mathcal C^\mu{}_\mu=\alpha Q/Q_0\)；
  - `q2`: \(\mathcal C^\mu{}_\mu=\alpha Q/Q_0+\beta(Q/Q_0)^2\)；
  - 当前 \(Q=X_{\rm flat}-m^2\)。
- tau=0 结果：
  - `zero`: 中心 weighted mean `0.02638`，p95 `0.07935`；
  - `q1`: 中心 weighted mean `2.01557`，p95 `5.86496`；
  - `q2`: 中心 weighted mean `0.78269`，p95 `2.20743`。
- `trace=0` 三切片：
  - tau=-3.5 中心 weighted mean `0.11846`；
  - tau=0 中心 weighted mean `0.02638`；
  - tau=+3.5 中心 weighted mean `0.11031`；
  - 跨时间明显比无 trace 的 `1%~5%` 差。
- `shear0` 即 \(D=0\) 三切片：
  - tau=-3.5 中心 weighted mean `0.10032`；
  - tau=0 中心 weighted mean `0.03105`；
  - tau=+3.5 中心 weighted mean `0.13447`；
  - reduced 主符号满秩但 p95 条件数可到 `1e9` 量级。
- 结论：
  - 简单代数状态方程能形式补足缺失自由度，但会显著恶化残差或病态化；
  - 只依赖 \(Q\) 的 trace closure 明确失败；
  - 下一步应转向辅助场作用量/变分型闭合，例如最小范数、最小梯度或带 \(M^{IJ},K^{IJ\mu\nu}\) 的辅助泛函。
- 新建研究笔记：
  - `research-notes/149-gBCD状态方程trace与剪切闭合首轮检验.md`。

### 2026-05-07 \(Q\to0\) 条件作为辅助场变分权重

- 用户澄清：不是要求用 \(F(Q)\) 替代，而是问 \(Q\to0\) 回到 Einstein 方程能否作为约束条件。
- 已把该条件改写为辅助场 action 的权重条件：
  \[
  S_{\lambda,\mathrm{mass}}
  =
  \frac12\int\sqrt{|\tilde g|}\,
  \mu(Q)M^{IJ}\lambda_I\lambda_J,
  \qquad
  \mu(Q)=1/(\epsilon+|Q|/Q_0)^p.
  \]
  这样 \(Q\to0\) 时 \(\lambda_I=(A,B,C,D)\) 代价增大，倾向 \(\mathcal C_{\mu\nu}\to0\)。
- 升级 `kg_examples/fit_gbcd_auxiliary_gauge.py`：
  - 新增 `--q-gated-norm`、`--q-gate-eps`、`--q-gate-power`、`--q-gate-max`；
  - 默认本轮测试用 `eps=0.03,power=2,max=1e4`。
- 三切片输出：
  - `visualizations/equation_first_gbcd_aux_qgated_n96_taum3p5_norm_time/`
  - `visualizations/equation_first_gbcd_aux_qgated_n96_tau0_norm_time/`
  - `visualizations/equation_first_gbcd_aux_qgated_n96_taup3p5_norm_time/`
- 与未门控 `time` 代表元相比：
  - tau=-3.5 中心残差 `0.03252 -> 0.03612`；
  - tau=0 中心残差 `0.01199 -> 0.01431`；
  - tau=+3.5 中心残差 `0.04692 -> 0.05029`；
  - 残差小幅变差，但空间粗糙度和系数范数整体下降，full conservation 仍保持数值精度。
- tau=0 低 \(|Q|\) 诊断：
  - 未门控解在 \(|Q|/Q_{95}<0.05\) 区域已满足 \(|\mathcal C|/|\mathcal C|_{95}\approx0.00456\)；
  - 门控后 \(|Q|/Q_{95}<0.1\) 区域 weighted mean 从 `0.03982` 降到 `0.03720`；
  - 说明 \(Q\to0\) 条件与当前最小范数/时间平滑代表元相容。
- 判断：
  - \(Q\to0\) 回 Einstein 适合作为辅助场 action 的边界/权重条件；
  - 不应简单写成 trace=\(F(Q)\)；
  - 下一步应推导带 \(\mu(Q)M^{IJ}\lambda_I\lambda_J+K^{IJ\alpha\beta}\nabla_\alpha\lambda_I\nabla_\beta\lambda_J\) 的连续辅助场 Euler-Lagrange 方程。
- 新建研究笔记：
  - `research-notes/150-Q趋零条件作为辅助场变分权重.md`。
- 已完成辅助场变分闭合推导：
  - 新建 `research-notes/151-gBCD辅助场变分闭合的Euler-Lagrange方程.md`；
  - 固定参考切片上的泛函为：
    残差项 + \(\mu(Q)M^{IJ}\lambda_I\lambda_J\) + \(K^{IJ\alpha\beta}\nabla_\alpha\lambda_I\nabla_\beta\lambda_J\) + conservation 乘子；
  - Euler-Lagrange 方程：
    \(H_{IJ}\lambda_J-b_I+\mu M_{IJ}\lambda_J-\nabla_\alpha(K_{IJ}^{\alpha\beta}\nabla_\beta\lambda_J)-E^I_{\mu\nu}\nabla^\mu\xi^\nu=0\)，
    配合 \(\nabla^\mu(\lambda_I E^I_{\mu\nu})=0\)。
- 已做 tau=0 的 \(Q\)-门控强度扫描：
  - 输出 `visualizations/equation_first_gbcd_aux_qgated_scan_n96_tau0_summary/`；
  - 扫描 \(\epsilon=0.01,0.03,0.1\)，\(p=1,2\)；
  - \(p=1\) 基本不伤残差但正则弱；
  - \(p=2\) 正则更明显但残差上升；
  - 当前温和折中点可先取 \(\epsilon=0.1,p=2\)，中心残差约 `0.01280`，p95 `0.04352`。

### 2026-05-07 无 Q 门控的最小辅助场 action 扫描

- 用户指出 3+1d 中 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\) 有 4 个分量；已修正理解：
  - 此前“少一条方程”只直接对应 `2+1 txz` 约化；
  - 更精确问题是当前 gBCD 张量壳的守恒主符号列落在 \(\mathrm{span}\{k,u,r\}\)，一般 rank 至多 3；
  - 因此仍需要辅助场 action 或其它机制选择剩余自由度。
- 按用户要求暂时不考虑 \(Q\to0\)，只扫描最小辅助场 action：
  - norm 项；
  - time gradient 项；
  - space gradient 项。
- tau=0 时间梯度扫描：
  - 输出 `visualizations/equation_first_gbcd_aux_noq_time_scan_n96_tau0_tw_*/`；
  - `time_weight=0 -> 1e-4` 时，中心残差只从 `0.01187` 到 `0.01270`；
  - 时间粗糙度 p95 从 `0.00499` 降到 `0.00151`；
  - 时间梯度项有效且温和。
- tau=0 空间梯度扫描：
  - 输出 `visualizations/equation_first_gbcd_aux_noq_space_scan_n96_tau0_sw_*/`；
  - `space_weight<=1e-10` 几乎无改善；
  - `space_weight=1e-8` 将中心残差推到 `0.01996`，且空间粗糙度/系数 p95 变差；
  - 当前简单邻点差分空间项暂不作为默认闭合项。
- 三切片时间梯度复查：
  - 输出 `visualizations/equation_first_gbcd_aux_noq_time_scan_n96_3tau_summary/`；
  - `time_weight=1e-5` 是当前较稳折中：
    tau=-3.5 残差 `0.03506`，
    tau=0 残差 `0.01241`，
    tau=+3.5 残差 `0.05017`；
  - 它显著降低 tau=+3.5 系数 p95：`1.214 -> 0.875`。
- 新建研究笔记：
  - `research-notes/152-无Q门控的最小辅助场action扫描.md`。

### 2026-05-07 gBCD 显式理论候选 v0

- 按用户提醒，明确最终目标是显式理论而非数值拟合器。
- 新建 `research-notes/153-gBCD显式理论候选v0.md`：
  - 场变量：\(\tilde g_{\mu\nu},\rho,S,\lambda_I=(A,B,C,D)\)；
  - 主场方程：
    \[
    \tilde G_{\mu\nu}=\tilde T_{\mu\nu}/M_P^2+\lambda_I E^I_{\mu\nu};
    \]
  - 守恒条件：
    \[
    \tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0;
    \]
  - 闭合原则：
    在主场方程和守恒约束解空间中，使
    \[
    \mathcal A_\lambda=\frac12\int\sqrt{|\tilde g|}
    [M^{IJ}\lambda_I\lambda_J+K_t^{IJ}(n^\mu\nabla_\mu\lambda_I)(n^\nu\nabla_\nu\lambda_J)]
    \]
    最小。
- v0 默认参数：
  - `norm_weight=1e-8,time_weight=1e-5,space_weight=0,q_gate=False`。
- 写出受限变分形式，乘子 \(\Xi^{\mu\nu}\) 强制主场方程，\(\xi^\nu\) 强制守恒；
  \(\lambda_I\) 的 EL 方程为：
  \[
  M^{IJ}\lambda_J-\nabla_\alpha(K^{IJ\alpha\beta}\nabla_\beta\lambda_J)
  -\Xi^{\mu\nu}E^I_{\mu\nu}-E^I_{\mu\nu}\nabla^\mu\xi^\nu=0.
  \]
- 三层离散时间诊断：
  - tau=0 干涉中心 \(\lambda\) 时间二阶差分 p95 很小，约 `0.002~0.007`；
  - tau=-3.5 的 `rr/ur` 二阶差分 p95 可到 `1.38/1.07`；
  - 后续若做时间推进，`rr/ur` 是主要风险自由度，可能需要隐式/半隐式推进或不同权重。

### 2026-05-07 gBCD v0 一步 lambda 推进诊断

- 新建脚本：
  - `kg_examples/test_gbcd_v0_one_step_lambda.py`。
- 目的：
  - 固定 A 支背景，只用已知 \(\lambda_-,\lambda_0\) 预测 \(\lambda_+\)，检验是否能替代逐切片三层投影。
- 输出：
  - `visualizations/equation_first_gbcd_v0_one_step_n96_taum3p5/`
  - `visualizations/equation_first_gbcd_v0_one_step_n96_tau0/`
  - `visualizations/equation_first_gbcd_v0_one_step_n96_taup3p5/`
- 结果：
  - hard conservation 可满足到机器精度，`equality_max_abs ~ 1e-16`；
  - 一步预测 \(\lambda_+\) 基本复现完整三层投影中的 \(\lambda_+\)，tau=0/tau=+3.5 的 p95 平均差异小于 `0.4%`；
  - 但下一切片代数场方程残差仍高：weighted mean 约 `0.137~0.170`，p95 约 `0.57~0.69`。
- 解释：
  - 这不是代码失败，而是固定 A 支背景、只推进 \(\lambda\) 的物理/方程兼容性不足；
  - 真正 D 支动力学必须允许 \(\rho,S,\tilde g\) 或至少 \(\rho,S\) 与 \(\lambda\) 一起变化。
- 新建研究笔记：
  - `research-notes/154-gBCD-v0一步lambda推进诊断.md`。

### 2026-05-07 gBCD 约束-演化分裂与 principal 投影

- 发现并修正一步推进理解：
  - 真正的二阶几何演化中，中心切片场方程含 \(\tilde g_+\) 的二阶时间差分；
  - 不应把“固定 A 背景上的下一切片代数残差”当成显式一步验收条件；
  - 应把场方程分成当前切片约束组合和可由 \(\tilde g_+\) 推进的演化组合。
- 新建脚本：
  - `kg_examples/diagnose_gbcd_metric_principal_correction.py`；
  - `kg_examples/fit_gbcd_principal_constraint_projection.py`。
- 主部诊断：
  - \(\delta\tilde G_{\mu\nu}/\delta\tilde g_+\) 的 time-time 主部逐点 rank=3；
  - \(\Delta t/\Delta x\simeq6e-7\)，所以该主部是显式时间推进最强局部项；
  - rank=3 意味着 3 个演化组合 + 3 个约束组合。
- 直接用旧 v0 \(\lambda\) 时：
  - 主部修正只能轻微降低残差；
  - tau=0 几乎不降，说明旧 v0 残差主要在约束方向。
- principal-constraint projection：
  - 新硬约束为 full conservation 加主部左零空间约束 \(N_A^{\mu\nu}(C-r_{\rm need})_{\mu\nu}=0\)；
  - 三组切片都跑通，hard equality max 约 `1e-13~1e-14`；
  - 投影后中心残差 weighted mean：
    - tau=-3.5: `3.57e-4`
    - tau=0: `4.92e-5`
    - tau=+3.5: `2.79e-4`
  - 所需相对 \(|\delta\tilde g_+|\) p95 约 `1e-11~1e-9`。
- 完整几何代回复查：
  - 将 pointwise principal \(\delta\tilde g_+\) 直接代回完整 `metric_jets_full -> geometry_data` 后，exact residual 仍为：
    tau=-3.5 `0.3407`，tau=0 `0.1529`，tau=+3.5 `0.03685` weighted mean；
  - 说明逐点 principal correction 会激发 mixed/空间导数项，不能作为最终一步演化器；
  - 简单邻点平滑会破坏主部方程，反而变差。
- 解释：
  - 这比 lambda-only 推进更接近真正几何演化算法；
  - 但下一步必须构造全局 full-linear \(\delta\tilde g_+\) solve，再接入物质方程推进 \(\rho,S\)。
- 新建研究笔记：
  - `research-notes/155-gBCD约束演化分裂与principal投影.md`。

### 2026-05-07 高分辨率口径与 full-linear metric update 首轮

- 用户要求把“高分辨率处理干涉条纹”和“全局 full-linear metric update”一起推进。
- 分辨率复查：
  - `n=96` 局部窗口只有 `29x29`；
  - \(\Delta x\simeq0.617um\)，每个 1550nm 波长只有 `2.51` 点；
  - 确认只能做低分辨率结构诊断，不能作为干涉/曲率可信证据。
- 新建脚本：
  - `kg_examples/render_a_branch_resolution_check.py`。
- `n=384` A 支图：
  - 输出 `visualizations/a_branch_resolution_check_n384_window9um/`；
  - 局部窗口 `117x117`；
  - \(\Delta x\simeq0.154um\)，每波长 `10.05` 点；
  - 这是下一档最低可用条纹分辨率。
- 新建脚本：
  - `kg_examples/solve_gbcd_full_linear_metric_update.py`；
  - 对当前 Einstein tensor 离散 stencil 关于 \(\tilde g_+\) 做 full-linear dense prototype；
  - 包括 local `tt`、mixed `tx/tz` 和 connection-quadratic linear terms。
- tau=0 ridge 扫描：
  - `ridge=1e-10` 线性残差低但非线性代回差，exact wmean `0.2298`；
  - `ridge=1e-6` 当前最好，exact wmean `0.00661`，p95 `0.0495`；
  - `ridge=1e-4` 接近但略差；
  - `ridge=1e-2` 正则过强。
- 三切片 `ridge=1e-6`：
  - tau=-3.5 exact wmean `0.0679`，p95 `0.861`；
  - tau=0 exact wmean `0.00661`，p95 `0.0495`；
  - tau=+3.5 exact wmean `0.0824`，p95 `0.703`。
- 判断：
  - full-linear global update 明显优于 pointwise principal correction；
  - 但当前 dense prototype 不能直接上 `n=384`；
  - 下一步要把 principal-constraint projection 与 full-linear metric update 改成 sparse/matrix-free。
- 新建研究笔记：
  - `research-notes/156-高分辨率口径与full-linear-metric-update首轮.md`。

### 2026-05-07 sparse/matrix-free 求解器首轮

- 新建 `kg_examples/solve_gbcd_full_linear_metric_update_sparse.py`：
  - 稀疏列存储 full-linear metric update；
  - 支持 CG normal equation 和 LSQR augmented least-squares。
- 算子一致性检查：
  - sparse forward 与 dense \(A\) 相对差 `1.75e-16`；
  - sparse transpose 与 dense \(A^T\) 相对差 `1.54e-16`；
  - 稀疏装配本身正确。
- 迭代求解差距：
  - dense `ridge=1e-6,tau=0` exact residual wmean `0.00661`；
  - sparse CG 2000 次 exact residual wmean `0.01156`；
  - sparse LSQR 1000 次 exact residual wmean `0.01413`；
  - 说明当前瓶颈是迭代求解/预条件，不是算子装配。
- 新建 `kg_examples/fit_gbcd_principal_constraint_projection_sparse.py`：
  - sparse LSQR 版 principal projection；
  - 当前用 penalty rows 近似 hard force/principal constraints，不是最终 hard/nullspace 版本。
- `n=96,core10` 对照：
  - dense hard principal + dense full-linear exact residual wmean `0.02548`；
  - sparse soft-principal + sparse full-linear exact residual wmean `0.07885`；
  - 差距主要来自 sparse principal 还不是 hard constraint，以及 metric update 迭代不够强。
- 新建研究笔记：
  - `research-notes/157-sparse-matrixfree求解器首轮.md`。

### 2026-05-07 n=384 core10 首个高分辨率 D 支几何更新结果

- 按用户要求继续推进数值求解器和物理主线。
- `n=384,tau=0,core10` sparse pipeline 跑通：
  - principal projection:
    `ncols=19716, rows=55789, nnz=171920`；
  - metric update:
    `fit_points=1024, variables=9858, rows=6144, nonzero_column_blocks=27620`。
- 结果：
  - before residual weighted mean `1.13493`；
  - linear after residual weighted mean `5.53e-5`；
  - exact nonlinear after residual weighted mean `0.00644`；
  - exact nonlinear p95 `0.04337`；
  - exact nonlinear max `0.28371`；
  - exact residual p99 `0.12366`。
- corrected \(\tilde g_+\) determinant in core10:
  - min `13.12`；
  - p50 `2294.99`；
  - p95 abs `3475.72`；
  - negative count `0`，positive count `1024`。
- 解释：
  - 在能解析干涉条纹的高分辨率核心区，D 支约束投影 + full-linear 几何更新没有被否定；
  - 干涉中心核心区已进入物理可讨论阶段。
- `n=384,tau=-3.5,core10` 快速尝试：
  - metric update CG `1500` iterations 后 normal residual 仍 `0.00296`；
  - exact residual weighted mean `0.40698`；
  - ridge `1e-4` 仍不好；
  - 判断为求解器未收敛/分离态条件更差，不作为物理失败证据。
- block-Jacobi 预条件：
  - n=96 core10 上 normal residual 可降低，但 exact nonlinear residual 变差；
  - 暂不采用。
- 新建研究笔记：
  - `research-notes/158-n384-core10首个高分辨率D支几何更新结果.md`。

### 2026-05-07 n=384 三切片 gBCD 求解器升级

- 用户要求“尽快解决数值求解器的问题，然后继续推进物理问题的主线”。
- 已将 `kg_examples/solve_gbcd_full_linear_metric_update_sparse.py` 的 full-linear metric update 稀疏 matvec 从 Python column-loop 升级为 `entry-arrays + numpy.bincount`：
  - 新增 `--operator-backend entry-arrays|python-column-loops`；
  - 默认 `entry-arrays`；
  - 方程与线性算子定义不变，只优化求解器实现。
- 已将 `kg_examples/fit_gbcd_principal_constraint_projection_sparse.py` 的 sparse LSQR matvec 同样升级为 entry-array 向量化。
- 已新建 `kg_examples/plot_gbcd_metric_update_diagnostics.py`：
  - 同图展示 A 支 `rho`、更新前残差、线性更新后残差、完整非线性代回残差、\(|\delta g_+|/|g_+|\)、`det(corrected g+)`；
  - 图中白线已明确定义为 fitting mask 边界。
- 三切片 `n=384,core10` 结果：
  - `tau=0`：exact nonlinear residual weighted mean `0.006437`，p95 `0.04337`；`core50` weighted mean `7.56e-5`；det 全正。
  - `tau=-3.5`：将 projection 从旧 fast 500 步升级到 3000 步后，projection system residual 从 `6.65e-2` 降到 `5.39e-3`；full-linear update 后 exact weighted mean `0.14577`，p95 `0.83019`；det 有 1 个近零/非正点。
  - `tau=+3.5`：projection system residual `6.86e-3`；full-linear update 后 exact weighted mean `0.03864`，p95 `0.10831`；det 全正。
- 阻尼/线搜索检查：
  - 对 `tau=-3.5` 的 \(\delta \tilde g_+\) 扫描 \(\alpha=0\ldots1.5\)；
  - `alpha=1.0` 最优，欠阻尼不能降低残差；
  - 因此当前左侧分离态问题不是 Newton 过冲。
- 关键输出：
  - `visualizations/equation_first_gbcd_sparse_pipeline_n384_tau0_core10_pen1e4/diagnostics/gbcd_metric_update_diagnostic_maps.png`
  - `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_iter3000_coeff_cg6000/diagnostics/gbcd_metric_update_diagnostic_maps.png`
  - `visualizations/equation_first_gbcd_sparse_pipeline_n384_taup3p5_core10_pen1e4_iter3000_coeff_cg6000/diagnostics/gbcd_metric_update_diagnostic_maps.png`
- 新建研究笔记：
  - `research-notes/159-n384三切片gBCD求解器升级与物理主线判断.md`。
- 当前判断：
  - 求解器已从“不能高分辨率调试”推进到“可在数分钟内完成三切片核心诊断”；
  - 干涉中心和右侧分离态支持继续 gBCD 物理主线；
  - 左侧分离态仍是主要卡点，下一步应做误差分层、hard/nullspace 稀疏投影、代表元规范检查，而不是继续盲目换 CG/LSQR 或加阻尼。

### 2026-05-07 tau=-3.5 残差来源与 center-plus 初值投影

- 新建 `kg_examples/diagnose_gbcd_metric_update_residual_sources.py`：
  - 按 `core50/core10_only`、mask 边缘层、\(|det|\)、\(|\delta g|/|g|\)、张量分量分层诊断完整非线性残差；
  - 图中 `boundary_layer_0` 定义为 fitting mask 内至少有一个 8-邻居在 mask 外的边缘格点。
- `tau=-3.5,matter4,force-mask-erosion=1` 残差来源：
  - 全 mask exact residual weighted mean `0.14577`；
  - `boundary_layer_0` weighted mean `0.3092`，是最差区域；
  - `core50` weighted mean `0.1124`，说明问题不只在边缘；
  - 张量分量贡献分散，不是单个分量爆炸。
- `force-mask-erosion=0` 测试：
  - force rows `6213 -> 7650`；
  - full-linear update exact residual weighted mean `0.14577 -> 0.13917`；
  - 边缘守恒约束有帮助但不根治。
- `time_weight=0` 测试：
  - full-linear update exact residual weighted mean 变为 `0.16319`；
  - `tau=-3.5` 问题不是简单时间方向代表元偏置。
- `matter6_normal` 扩展：
  - 已给 projection/update 增加 `--atom-family matter4|matter4_plus_normal|matter6_normal`，并保存/自动读取 `atom_family/atom_names`；
  - `matter6_normal` projection system residual `5.39e-3 -> 3.60e-4`；
  - 但 full-linear exact residual `0.15269`，`matter6+edge` 为 `0.15187`，均不如 `matter4+edge`；
  - 说明补足局域张量基底改善投影层，但不是非线性闭合主因。
- `center_plus` 初值投影：
  - `solve_gbcd_full_linear_metric_update_sparse.py` 新增 `--metric-variable-slices plus|center_plus|minus_center_plus`；
  - `center_plus` 同时修正 \(g_0,g_+\)，旧默认 `plus` 不变；
  - `matter4+edge+center_plus,ridge=1e-6` 线性残差 weighted mean `0.08504`，但 exact nonlinear `0.92042`；
  - 新建 `kg_examples/scan_gbcd_metric_update_alpha.py`，线搜索显示 \(\alpha=1\) 是该方向最优但仍 `0.92042`；
  - `ridge=1e-3` 线性残差降到 `0.05835`，exact nonlinear 仍 `0.90481`。
- 当前判断：
  - 只改几何的 `center_plus` 会得到线性好、非线性不自洽的更新；
  - 若允许 \(g_0\) 改变，必须同步重算/约束 \(\tilde T_{\mu\nu}\)、\(\mathcal C_{\mu\nu}\)、\(u,r,\rho^\tilde\) 拉回、质量壳和连续性；
  - 下一步应实现联合 D 初值投影，而不是继续只调单独 metric update。
- 新建研究笔记：
  - `research-notes/160-tau负3p5残差来源与center-plus初值投影检验.md`。

### 2026-05-07 D 初值联合投影与 plus-only 主自由度

- Git 状态：
  - 已本地初始化仓库并提交 `8cf55c0 Initial project snapshot`；
  - 本轮新增内容已提交 `11a3ff0 Add joint D initial projection diagnostics`；
  - 已配置远端 `https://github.com/ustc-wyf/geometrized_qm.git`；
  - 初次推送时终端到 GitHub 超时；网络恢复后 HTTPS push 报 `could not read Username for 'https://github.com'`，说明命令行 Git 缺少 GitHub token/credential；
  - SSH 测试显示本机有 `~/.ssh/id_ed25519`，但未绑定 GitHub：`Permission denied (publickey)`；
  - 因此远端 push 尚未完成，本地仓库干净，远端配置保留。
  - 用户绑定 SSH key 后，`ssh -T git@github.com` 已认证为 `ustc-wyf`；
  - 已将 `origin` 切到 `git@github.com:ustc-wyf/geometrized_qm.git`；
  - 已成功 `git push -u origin main`，当前 `main...origin/main` 干净同步。
- 新建 `kg_examples/solve_gbcd_joint_initial_projection_sparse.py`：
  - 同一稀疏线性系统中求 \(\delta\tilde g\) 与可选 \(\eta=\delta\log(\sqrt{|\tilde g|}\tilde\rho)\)；
  - 求解后重新计算 \(\sqrt{|\tilde g|}\)、\(\tilde\rho\)、质量壳 \(u_t\)、\(\tilde T_{\mu\nu}\)，并做完整 nonlinear back-substitution。
- 新建 `kg_examples/scan_gbcd_joint_projection_alpha.py`：
  - 对联合解做 trust-region \(\alpha\) 扫描；
  - 同时检查 D 张量残差、\(\rho\) pullback 偏差、质量壳判别式、\(\det\tilde g\)。
- `center_plus+source` 高分辨率测试：
  - 输出 `visualizations/equation_first_gbcd_joint_initial_projection_sparse_n384_taum3p5_core10_centerplus_rhoclose1e6_v2/`；
  - exact residual weighted mean `0.81379`；
  - 未裁剪 \(\eta\) weighted mean 约 `1.05e12`，裁剪比例 `1.0`；
  - \(\rho\) pullback 加权相对 L1 偏差 `0.51344`；
  - 负质量壳判别式比例 `0.01882`；
  - 判断：源项密度幅值自由度在 \(M_P^2\) 压制下不健康，不能用于真实 D 初态。
- `center_plus+source` trust-region 扫描：
  - 输出 `visualizations/equation_first_gbcd_joint_projection_alpha_scan_n384_taum3p5_centerplus_rhoclose1e6_v2/`；
  - 若限制 \(\rho\) 偏差 \(\le 5\%\)，可用 \(\alpha\lesssim0.1\)，exact residual 仍约 `1.16685`；
  - 没有满足默认可接受条件的 feasible alpha。
- `plus-only,no-source` 高分辨率测试：
  - 输出 `visualizations/equation_first_gbcd_joint_initial_projection_sparse_n384_taum3p5_core10_plus_nosource_v2/`；
  - `LSQR1000` exact residual weighted mean `0.21334`；
  - \(\rho\) pullback 偏差 `0`，负判别式比例 `0`，中心切片 \(\delta g_0=0\)；
  - 既有 `plus-only+CG6000` 对照 exact residual weighted mean `0.13917`。
- 当前判断：
  - 真正健康的 D 初值自由度是 `plus-only`/初始加速度，即保持初始 \(\tilde g_0,\rho^\tilde,u_i\) 不变，由 D 场方程确定下一切片 \(g_+\) 或 \(\partial_t^2\tilde g\)；
  - 不应继续用大幅 \(\delta\rho^\tilde\) 修几何残差；
  - 下一步应把 plus-only 解法包装为正式 D 初始加速度求解器，并接入短步物质演化与约束漂移检查。
- 新建研究笔记：
  - `research-notes/161-D初值联合投影与plus-only主自由度.md`。

### 2026-05-07 沟通偏好修正

- 用户反馈刚才的数值结果报告“看不懂，语言跳跃太大”。
- 后续报告要求：
  - 先给一句话结论；
  - 再说明本轮“在问什么问题”；
  - 明确定义新变量、新名词和图中每个量；
  - 区分“物理结论”和“数值算法诊断”；
  - 给出关键数字时说明其含义、好坏方向和可接受性；
  - 不要直接堆 `weighted_mean/p95/eta/alpha` 等术语而不解释；
  - 复杂结果用“为什么做、怎么做、看到什么、所以怎么办”的顺序。

### 2026-05-07 plus-only 初始加速度包与局部可容许性守卫

- 按用户“继续”要求，接续 D/gBCD 初值求解器主线。
- 新增/修改脚本：
  - `kg_examples/export_gbcd_plus_initial_package.py`：导出 plus-only 三切片初始包，新增 `--plus-guard auto_bad_zero`；
  - `kg_examples/test_gbcd_plus_initial_package_one_step.py`：读取初始包并做一步物质推进/重构诊断；
  - `kg_examples/scan_gbcd_plus_alpha_admissibility.py`：扫描全局 \(\alpha\delta g_+\) 的 D 残差和物质可容许性；
  - `kg_examples/scan_gbcd_plus_local_guard.py`：扫描局部 trust-region/guard 规则。
- 旧 active-edge `plus-only+CG6000` 诊断：
  - 输出 `visualizations/equation_first_gbcd_plus_initial_package_n384_taum3p5_core10_cg6000/`；
  - D 方程 exact residual weighted mean `0.13917`；
  - 一步 \(g_+\) 测度偏差 `14404.69`；
  - 负质量壳判别式比例 `0.01922`；
  - 全局阻尼扫描显示 `alpha=1e-5` 仍有负判别式，说明不能靠整体缩步解决。
- 旧 active-edge 坏点分层：
  - 坏点主要集中在支撑区边缘/低密度带；
  - 但 p95 更新很小而 max 更新极大，说明少数局部自由度不可容许。
- 新跑 `metric_active_dilation=0` 的 no-active-edge 重新求解：
  - 输出 `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_noactiveedge_coeff_cg6000/`；
  - D 方程 exact residual weighted mean `0.13192`；
  - \(|\delta g_+|/|g_+|\) p95 `4.68e-5`，max `0.00931`；
  - 消除了旧解的百万级局部异常。
- no-active-edge 剩余坏点定位：
  - 只剩 `1/6034` 个支撑点出现负判别式；
  - 坐标约 \((-3.392\,\mu m,-6.013\,\mu m)\)；
  - `rho/max(rho)=0.148`，不是低密度边缘；
  - 参考 \(\tilde g\) 本身近退化，`det(g_plus_base) ~ 1.24e-24`，一个特征值约 `1e-17`；
  - 该点 `delta g_+/g_+ ~ 3.1e-8`，因此问题是分支/退化点，不是更新幅度过大。
- `auto_bad_zero` 局部守卫扫描：
  - 输出 `visualizations/equation_first_gbcd_plus_local_guard_scan_n384_taum3p5_core10_noactiveedge_cg6000/`；
  - 冻结 2 个不可容许支撑点；
  - D 方程 exact residual weighted mean `0.13220`；
  - 一步 \(g_+\) 测度偏差 `1.09e-4`；
  - 负判别式比例和非正 det 比例均为 `0`。
- guarded 初值包：
  - 输出 `visualizations/equation_first_gbcd_plus_initial_package_n384_taum3p5_core10_noactiveedge_guarded/`；
  - 初始 \(\rho^\tilde\) 拉回偏差 `1.4e-17`；
  - mass-shell defect p95 `2.24e-13`。
- guarded 一步测试：
  - 输出 `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_taum3p5_core10_noactiveedge_guarded/`；
  - 一步 \(g_+\) 测度偏差 `1.091e-4`；
  - 负判别式比例 `0`；
  - 守恒密度一步相对变化 `5.79e-7`。
- 语法检查：
  - `python3 -m py_compile` 覆盖新增/修改脚本，通过。
- 新建研究笔记：
  - `research-notes/162-plus-only初始加速度包与可容许性守卫.md`。

### 2026-05-07 三切片 guarded plus-only 初始包验证

- 按计划把 `no-active-edge + auto_bad_zero` guarded 初始加速度流程扩展到三组关键切片：`tau=-3.5,0,+3.5`。
- `tau=0`：
  - 重新求解输出 `visualizations/equation_first_gbcd_sparse_pipeline_n384_tau0_core10_pen1e4_noactiveedge_coeff_cg6000/`；
  - D 方程 exact residual weighted mean `6.604e-5`；
  - `auto_bad_zero` 冻结点 `0`；
  - 一步 \(g_+\) 测度偏差 `4.994e-4`，加权偏差 `2.999e-5`；
  - 负判别式比例 `0`。
- `tau=+3.5`：
  - 重新求解输出 `visualizations/equation_first_gbcd_sparse_pipeline_n384_taup3p5_core10_pen1e4_noactiveedge_coeff_cg6000/`；
  - D 方程 exact residual weighted mean `0.01324`；
  - `auto_bad_zero` 冻结支撑点 `2`，core10 冻结点 `0`；
  - 一步 \(g_+\) 测度偏差 `2.464e-4`，加权偏差 `2.884e-5`；
  - 负判别式比例 `0`。
- 三切片汇总：
  - 输出 `visualizations/equation_first_gbcd_guarded_three_tau_summary/summary.json`；
  - 线性图 `visualizations/equation_first_gbcd_guarded_three_tau_summary/guarded_three_tau_summary_linear.png`；
  - 对数图 `visualizations/equation_first_gbcd_guarded_three_tau_summary/guarded_three_tau_summary_log.png`。
- 当前判断：
  - guarded plus-only 初始包不是只在单个切片偶然可用；
  - 三切片一步物质重构均无负判别式；
  - 左侧分离态 `tau=-3.5` 仍是残差最高的瓶颈，但已经可推进；
  - 下一步应实现多步推进原型，并把 `auto_bad_zero` 升级为求解器内的活动集/强约束。
- 新建研究笔记：
  - `research-notes/163-三切片guarded-plus-only初始包验证.md`。

### 2026-05-07 目标层级纠偏：数值器不是主目标

- 用户纠正：项目目标不是“做一个能演化的 D 支数值计算器”。
- 重新核对记忆与研究笔记后确认，当前真正主线是：
  - 在高斯波包干涉试验场上寻找能描述变换后几何的、类似 Einstein 方程的显式 equation-first 场方程；
  - action-first 的 pure-\(\tilde g\) 或简单 \(f(R)\)、\(f(R,I_2)\) 路线已在该试验场上严重受挫，因此转为 equation-first；
  - 数值计算器、初始包和一步推进只是检验候选方程是否可行、是否闭合、是否保持测地线解释的工具，不是最终目标。
- 当前候选方程仍应表述为 gBCD：
  \[
  \tilde G_{\mu\nu}
  =
  \tilde T_{\mu\nu}/M_P^2
  +A\tilde g_{\mu\nu}
  +B u_\mu u_\nu
  +C r_\mu r_\nu
  +D u_{(\mu}r_{\nu)} .
  \]
- 下一步应回到方程形式的一般化：
  - 明确 \(A,B,C,D\) 的本构/辅助场闭合；
  - 用三切片 guarded 结果说明哪些系数/约束结构被支持；
  - 不应把“多步演化器”说成目标，只能说成检验候选方程的工具。

### 2026-05-08 trace=0 最小投影方程组整理

- 使用 `agent-memory` 规则重载长期记忆，确认当前主线是寻找显式 equation-first Einstein-like 方程，而不是把 D 支数值器当作最终目标。
- 新建研究笔记：
  - `research-notes/169-trace0最小投影方程组与patch条件.md`。
- 将当前最小候选整理为：
  \[
  \tilde G_{\mu\nu}
  =
  \tilde T_{\mu\nu}/M_P^2
  +\mathcal C_{\mu\nu},
  \quad
  \mathcal C_{\mu\nu}\in
  \mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\},
  \]
  \[
  \tilde g^{\mu\nu}\mathcal C_{\mu\nu}=0,
  \qquad
  \tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
  \]
- 明确在 \(\tilde g\) 表象中变量应写为 \((\tilde\rho,S)\)，不能把 \(\rho\) 和 \(\tilde\rho\) 混用；与 A 支比较时必须先把 \(\tilde\rho\) 拉回 \(g\) 表象。
- 记录当前 patch 条件：
  - \(d>2\)；
  - \(\Delta=u^2r^2-(u\cdot r)^2\neq0\)；
  - trace=0 主符号还需避开或 patch 处理 \(w^2=0\)；
  - \(r^\perp\to0\) 或严格 \(1+1d\) 需要换低维 basis。
- 当前判断：
  - trace=0 是目前最小 algebraic closure；
  - \(\mathsf q_E\) 暂保留为主符号/辅助场正定范数候选；
  - 仍需做 Helmholtz/self-adjoint 检查和 \(Q\to0\Rightarrow\mathcal C\to0\) 分支选择分析。

### 2026-05-08 trace=0 投影方程 Helmholtz 检查

- 新建研究笔记：
  - `research-notes/170-trace0投影方程的Helmholtz检查.md`。
- 检查对象：
  \[
  \Pi_E^\perp\mathcal R=0,\qquad
  \tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0.
  \]
- 主要判断：
  - 把 6 条投影方程 + 1 条 trace 方程嵌为 tensor equation 时，需要选择 trace 代表方向 \(\Sigma_{\mu\nu}\)，这已经不像 Einstein 方程那样天然唯一；
  - 线性化有
    \[
    \delta F=M\,\delta\mathcal R+(\delta M)\mathcal R,
    \]
    其中 \(M[\tilde g,u,r]\) 含投影算子；
  - \((\delta M)\mathcal R\) 一般破坏 Frechet self-adjointness；
  - 纯 metric/matter 协变作用量在物质壳上有 4 个 Noether 恒等式，而投影 6 条 + trace 1 条看起来有 7 条 metric 约束，存在计数压力。
- 当前结论：
  - trace=0 投影方程可继续作为 equation-first 候选；
  - 不能声称它已经来自普通局域纯 metric 作用量；
  - action 化更合理的路线是辅助应力场 action 或乘子约束系统；
  - 下一步应写最小辅助 action 原型，并在 `n=384` 复查 trace=0 优势。

### 2026-05-08 最小辅助应力场 action 原型

- 新建研究笔记：
  - `research-notes/171-最小辅助应力场action原型.md`。
- 关键修正：
  - action 化时，真正进入 Einstein 方程右边的是
    \[
    \Theta_{\mu\nu}
    =
    -2/\sqrt{|\tilde g|}\,
    \delta S_{\rm aux}/\delta\tilde g^{\mu\nu},
    \]
    不是裸 \(\lambda_I E^I_{\mu\nu}\)。
  - 若 \(C_{\mu\nu}\) 是独立、无迹协变张量，源项
    \[
    -\frac{M_P^2}{2}\int\sqrt{|\tilde g|}\tilde g^{\mu\nu}C_{\mu\nu}
    \]
    可以产生 \(M_P^2 C_{\mu\nu}\) 作为 metric 源。
  - 但一旦直接令 \(C_{\mu\nu}=\lambda_I E^I_{\mu\nu}[\tilde g,u,r]\)，metric 变分会产生额外项，物质变分也通常会改变 Hamilton-Jacobi/连续性方程。
- 当前判断：
  - 裸 \(\lambda_I E^I\) action 不是正确目标；
  - 正确目标是寻找辅助 sector，使其有效能动量 \(\Theta_{\mu\nu}/M_P^2\) 满足投影、无迹、守恒和 \(Q\to0\) 分支；
  - 若要保持物质测地线，必须要求 \(\delta S_{\rm aux}/\delta S\) 和 \(\delta S_{\rm aux}/\delta\tilde\rho\) 在目标物质壳上消失，或引入独立 \(U,R\) shadow fields。

### 2026-05-08 shadow-field 分支最小方程组

- 新建研究笔记：
  - `research-notes/172-shadow-field分支的最小方程组.md`。
- 定义 shadow fields：
  \[
  \Phi,\sigma,\qquad
  U_\mu=\partial_\mu\Phi,\qquad
  R_\mu=\tilde\nabla_\mu\ln\sqrt{\sigma}.
  \]
- 目的：
  - 让辅助 action 依赖 \((\Phi,\sigma)\)，而不是直接依赖物质变量 \((S,\tilde\rho)\)；
  - 从而保持
    \[
    \delta S_{\rm aux}/\delta S=0,\qquad
    \delta S_{\rm aux}/\delta\tilde\rho=0.
    \]
- 关键判断：
  - 不能在 action 中用硬约束 \(U=u,R=r\)，否则乘子会重新改写物质变分；
  - 应让 shadow fields 满足同型动力学
    \[
    U^2=m^2,\qquad
    \tilde\nabla_\mu(\sigma U^\mu)=0,
    \]
    并通过初始/边界分支选择 \(\Phi=S+\mathrm{const},\sigma=\tilde\rho\)；
  - 若完整耦合系统适定且唯一，物理分支上可保持 \(U=u,R=r,Q_{\rm sh}=Q\)。
- 当前结论：
  - shadow-field 分支可以保护测地线解释；
  - 代价是引入额外自由度和物理分支选择；
  - 下一步要检查 shadow stress 是否可被 \(C_{\mu\nu}\) 状态方程吸收，以及分支传播是否成立。

### 2026-05-08 shadow stress 与物理分支传播检查

- 新建研究笔记：
  - `research-notes/173-shadow-stress与物理分支传播条件.md`。
- 检查结论：
  - shadow-field 分支只保证 \(S_{\rm aux}\) 不直接变分 \(S,\tilde\rho\)；
  - 但 \(S_{\rm state}\) 对 \(\tilde g,\Phi,\sigma\) 的变分仍会产生 shadow stress 和 shadow-source；
  - 线性乘子约束即使在约束成立时也会贡献 metric stress，并且 \(\Lambda_\perp(\Pi^\perp C)\) 会通过 \(\delta\Pi^\perp\) 给 \(\Phi,\sigma\) 加源；
  - \(\mu(Q_{\rm sh})C^2\) 会通过 \(\mu'(Q_{\rm sh})C^2\delta Q_{\rm sh}\) 改写 \(\sigma\) 方程，所以 \(Q\to0\Rightarrow C\to0\) 暂应作为分支/边界条件，而不是直接作为 variational 权重。
- 新条件：
  \[
  \mathcal E_\Phi^{\rm state}|_{\rm phys}=0,
  \qquad
  \mathcal E_\sigma^{\rm state}|_{\rm phys}=0.
  \]
  即状态项在物理分支上对 shadow fields “一阶静默”。
- 当前判断：
  - shadow-field 分支仍可行，但需要 stealth-state action；
  - 简单线性乘子或 \(Q_{\rm sh}\)-门控 action 都不够健康；
  - 下一步应构造具体 stealth-state action 候选，或承认 action 化暂时失败并退回 equation-first。

### 2026-05-08 stealth-state action 候选最小 no-go

- 新建研究笔记：
  - `research-notes/174-stealth-state-action候选的最小no-go.md`。
- 检查候选：
  \[
  S_{\rm aux}^{\rm trial}
  =
  S_{\rm source}
  +S_{\rm tr}
  +S_{\rm stealth}^{(2)}
  +S_{\rm sh}.
  \]
- 关键推导：
  - 线性源项
    \[
    S_{\rm source}
    =
    -\frac{M_P^2}{2}\int\sqrt{|\tilde g|}\,\tilde g^{\mu\nu}C_{\mu\nu}
    \]
    可以产生 metric 源 \(M_P^2C_{\mu\nu}\)；
  - 但若 \(C_{\mu\nu}\) 是自由变分场，\(C\)-变分给出 \(-M_P^2\tilde g^{\mu\nu}/2\)；
  - 用无迹乘子抵消该项会令线性源整体在壳上消失，从而 metric 源也消失；
  - 平方型 stealth-state 项在目标约束面上一阶变分为零，不能平衡这个线性源。
- 当前结论：
  - 最小 stealth-state action 不能生成非零 trace=0 投影源；
  - action 化若继续，必须寻找真实辅助场 \(\chi\)-sector 的有效应力；
  - 否则应把 trace=0 投影方程作为 equation-first 理论继续推进，并优先做 `n=384` 复检。

### 2026-05-08 n=384 trace=0 局部复检

- 新增轻量诊断脚本：
  - `kg_examples/diagnose_gbcd_trace_local_closure.py`。
- 背景：
  - 旧 `fit_gbcd_trace_closure.py` 在 `n=384` 下会构造很大的全局硬约束/nullspace 系统，进程被系统杀掉；
  - 因此本轮先做逐点局部代数 sanity check，不声称已经满足完整守恒或时间耦合。
- 新建研究笔记：
  - `research-notes/175-n384-trace0局部复检.md`。
- 运行输出：
  - `visualizations/equation_first_gbcd_trace_local_closure_n384_trusted_3tau/summary.json`
  - `visualizations/equation_first_gbcd_trace_local_closure_n384_trusted_3tau/local_trace_closure_summary.png`
  - `visualizations/equation_first_gbcd_trace_local_closure_n384_core10_3tau/summary.json`
  - `visualizations/equation_first_gbcd_trace_local_closure_n384_core10_3tau/local_trace_closure_summary.png`
- 关键数据：
  - `trusted` weighted mean residual：
    - tau=-3.5：trace0 `0.014764`，qe0 `0.012302`，但 qe0 跳过 `1071/4999` 点；
    - tau=0：trace0 `0.003195`，qe0 `0.019128`；
    - tau=+3.5：trace0 `0.010518`，qe0 `0.015883`。
  - `core10` weighted mean residual：
    - tau=-3.5：trace0 `0.008586`，qe0 `0.004855`，但 qe0 跳过 `146/2550` 点；
    - tau=0：trace0 `0.008569`，qe0 `0.132874`；
    - tau=+3.5：trace0 `0.010399`，qe0 `0.013445`。
- 当前判断：
  - 四个张量方向 \(\{\tilde g,uu,rr,ur\}\) 的 unconstrained 局部拟合残差仍很低，支持 gBCD 投影壳；
  - 普通 trace=0 在干涉中心和右侧分离态明显更稳；
  - \(\mathsf q_E\)-trace 只在左侧分离态局部略优，且有不可定义/退化跳点，不应作为主闭合；
  - 当前主候选仍是普通 trace=0 + full conservation + patch/branch 条件的 equation-first 方程。

### 2026-05-08 trace=0 局部守恒兼容性检查

- 新增脚本：
  - `kg_examples/diagnose_gbcd_trace_local_conservation.py`。
- 新建研究笔记：
  - `research-notes/176-trace0局部守恒兼容性.md`。
- 检查口径：
  - `n=384, core10, window=[-9,9] um`；
  - 逐点求局部 \(C_{\mu\nu}\)，再直接计算中心切片 \(\tilde\nabla^\mu C_{\mu\nu}\)；
  - 这不是 hard conservation solve，而是判断局部闭合是否自然守恒。
- trace0 结果：
  - tau=-3.5：代数残差 `0.004064`，full divergence scale `0.783756`；
  - tau=0：代数残差 `0.000598`，full divergence scale `0.000937`；
  - tau=+3.5：代数残差 `0.006150`，full divergence scale `0.254636`。
- unconstrained 基线：
  - tau=-3.5：代数残差 `0.000648`，full divergence scale `0.023113`；
  - tau=0：代数残差 `0.000069`，full divergence scale `0.000236`；
  - tau=+3.5：代数残差 `0.001647`，full divergence scale `0.007245`。
- 当前判断：
  - 干涉中点 trace0 局部代表元几乎自然守恒；
  - 分离态 trace0 局部代表元不自然守恒，且比 unconstrained 基线更差；
  - 因此 trace0 仍是主闭合，但必须和 full conservation 联立求解，不能先逐点 trace0 再事后检查守恒。

### 2026-05-08 trace0 稀疏全局守恒扫描

- 新增脚本：
  - `kg_examples/fit_gbcd_trace0_sparse_conservation.py`。
- 新建研究笔记：
  - `research-notes/177-trace0稀疏全局守恒扫描.md`。
- 方法：
  - 逐点对 \(\tilde g^{\mu\nu}C_{\mu\nu}=0\) 做 SVD nullspace 硬消元；
  - 在消元后三自由度上建立 sparse 线性系统；
  - 扫描完整守恒行的 penalty 权重 `force_weight`。
- 输出：
  - `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_taum3p5_core10_scan/`
  - `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_tau0_core10_scan/`
  - `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_taup3p5_core10_scan/`
- 关键数据：
  - tau=-3.5：`force_weight=0` 时 row_alg `0.004797`、full_div `0.783756`；`force_weight=10` 时 row_alg `0.005478`、full_div `0.007740`。
  - tau=0：`force_weight=0` 时 row_alg `0.000675`、full_div `0.000937`；`force_weight=10` 时 row_alg `0.000675`、full_div `0.000458`。
  - tau=+3.5：`force_weight=0` 时 row_alg `0.007307`、full_div `0.254636`；`force_weight=10` 时 row_alg `0.007698`、full_div `0.006502`。
- 当前判断：
  - trace0 与完整守恒在高分辨率三切片上兼容；
  - 分离态的大自然散度不是理论矛盾，而是逐点后处理没有联立守恒；
  - 下一步应把 penalty scan 升级成 hard conservation 或 saddle-point/KKT 矩阵自由求解器。

### 2026-05-08 trace0 硬守恒 KKT 与 ALM 原型

- 已扩展 `kg_examples/fit_gbcd_trace0_sparse_conservation.py`：
  - `--hard-project`：只投影到守恒面，不最小化代数残差；
  - `--hard-kkt`：矩阵自由 KKT constrained least-squares；
  - `--hard-alm`：增广拉格朗日/乘子法，把硬守恒问题转化为多轮正定 least-squares。
- 新建研究笔记：
  - `research-notes/178-trace0硬守恒KKT与ALM原型.md`。
- 关键运行输出：
  - `visualizations/equation_first_gbcd_trace0_hard_alm_n384_taum3p5_core10/`
  - `visualizations/equation_first_gbcd_trace0_hard_alm_n384_tau0_core10/`
  - `visualizations/equation_first_gbcd_trace0_hard_alm_n384_taup3p5_core10/`
- `hard-kkt` 判断：
  - 在 `n=96,tau=0` 上可把 full divergence 压到 `~3.6e-13`，代数残差保持 `~6.6e-4`；
  - 但在 `n=384,tau=-3.5` 上，当前无 MINRES 的 LSQR saddle-point 实现病态，不能把该失败解读成理论矛盾。
- `hard-alm` 三切片结果：
  - `tau=-3.5`：full divergence `0.007740 -> 0.001987`，中心 row algebraic `0.005478 -> 0.006009`；
  - `tau=0`：full divergence `0.000458 -> 0.0000239`，中心 row algebraic 基本不变 `~0.000675`；
  - `tau=+3.5`：full divergence `0.006502 -> 0.001940`，中心 row algebraic `0.007698 -> 0.008172`。
- 当前判断：
  - `trace0 + full conservation` 可以作为 equation-first 候选的硬约束结构继续保留；
  - ALM 只是数值实现，不改变物理方程；
  - 这一步不是 action 化成功，也不是完整演化理论完成，而是说明当前候选没有被高斯干涉三切片硬守恒检验否定。
