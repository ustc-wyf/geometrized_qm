# D 支 standalone 求解器首轮原型

日期：2026-05-08

## 0. 目的

接续 `186-D支harmonic闭合定理草案与当前例子可演化性.md`。

本轮开始构建真正从 D 初值出发的数值求解器，并后验比较它和 A 支 KG 演化的差异。

重要口径：

- A 支只用于初值包和后验比较；
- 演化过程中不从 A 支后续快照读取 \(\tilde g,\tilde\rho,u_\mu,C\)；
- 当前版本是最小 generalized-harmonic principal prototype，不是最终 full nonlinear implicit D solver。

新增脚本：

```bash
kg_examples/simulate_d_harmonic_standalone.py
```

## 1. 输入初值

使用现有 D 初值包：

```bash
visualizations/equation_first_gbcd_plus_initial_package_n384_tau0_core10_noactiveedge_guarded/gbcd_plus_initial_package.npz
```

该包包含：

- \(g_-,g_0,g_+\) 三张 \(\tilde g\) 时间片；
- \(\tilde\rho,u_x,u_z,u_t\)；
- \(C_{\mu\nu}\)；
- 初始 D 方程残差；
- support/core10/fit mask。

包本身诊断：

- 初始 D 方程 residual weighted mean 约 `6.6e-5`；
- 初始 mass-shell defect p95 约 `8.6e-16`；
- 初始 \(\rho\) pullback 与 A 支逐点一致到机器精度；
- 初始包里的 \(C\) 来自 `principal_constraint_projection`，不是最终 trace0+full-conservation 版本。

这个最后一点很重要：它解释了为什么不能把本轮结果直接当成最终 D 支预测。

## 2. 演化规则

脚本推进变量：

\[
\tilde\rho,\quad u_i,\quad \tilde g_{\mu\nu},\quad C_{\mu\nu}.
\]

物质推进：

\[
\tilde g^{\mu\nu}u_\mu u_\nu=m^2
\]

求 \(u_t\)，再用

\[
\partial_t(\sqrt{|\tilde g|}\tilde\rho u^t)
+\partial_i(\sqrt{|\tilde g|}\tilde\rho u^i)=0
\]

推进守恒密度。

几何推进：

- 使用 \(g_-,g_0,g_+\) 的三层历史；
- 修正了一个关键时间层 bug：\(g_+\) 是下一张中心切片，不能同时当作新中心和新未来片；
- 当前最小版本从
  \[
  (g_{n-1},g_n,g_{n+1})
  \]
  推到
  \[
  (g_n,g_{n+1},g_{n+2})
  \]
  时，用显式 harmonic principal extrapolation 得到 \(g_{n+2}\)。

\(C\)-sector 当前有三种测试口径：

1. `frozen`：冻结初始 \(C\)；
2. `local` 无 trace0：每步局部解 \(C\in\mathrm{span}\{g,uu,rr,ur\}\)；
3. `local` trace0：每步局部解 \(C\in E\) 且 trace0。

当前还没有接入 full conservation：

\[
\tilde\nabla^\mu C_{\mu\nu}=0.
\]

所以本轮仍不是最终 D 支求解器。

## 3. 运行结果

### 3.1 冻结 \(C\)

命令：

```bash
python3 kg_examples/simulate_d_harmonic_standalone.py \
  --package visualizations/equation_first_gbcd_plus_initial_package_n384_tau0_core10_noactiveedge_guarded/gbcd_plus_initial_package.npz \
  --output visualizations/d_harmonic_standalone_tau0_steps2_timefix \
  --tau 0 \
  --steps 2 \
  --evolve-region core10 \
  --metric-corrector none
```

结果：

- step 0 residual wmean core10：`6.72e-5`；
- step 1 residual wmean core10：`0.525`；
- step 2 residual wmean core10：`0.531`；
- step 2 \(\rho\) pullback L1 core10：`9.18e-5`；
- active 区 mass-shell discriminant 无负值。

解释：

冻结 \(C\) 后，物质短时仍接近 A，但 D 方程残差迅速变大。因此 \(C\)-sector 必须随物质和几何演化。

### 3.2 每步局部重解 \(C\)，不加 trace0

命令：

```bash
python3 kg_examples/simulate_d_harmonic_standalone.py \
  --package visualizations/equation_first_gbcd_plus_initial_package_n384_tau0_core10_noactiveedge_guarded/gbcd_plus_initial_package.npz \
  --output visualizations/d_harmonic_standalone_tau0_steps10_localC \
  --tau 0 \
  --steps 10 \
  --evolve-region core10 \
  --auxiliary-update local \
  --no-auxiliary-trace0 \
  --metric-corrector none
```

结果：

- step 0 residual wmean core10：`0.0231`；
- step 1 residual wmean core10：`0.0389`；
- step 10 residual wmean core10：`0.0390`；
- step 10 residual p95 core10：`0.194`；
- step 10 \(\rho\) pullback L1 core10：`4.55e-4`；
- step 10 \(\rho\) pullback weighted L1 core10：`2.15e-4`；
- active 区 mass-shell discriminant 无负值；
- \(|R^\tilde|\) 量级稳定，没有出现前一版时间层 bug 导致的爆炸。

解释：

局部更新 \(C\) 能把 D 方程残差从冻结 \(C\) 的 `~0.53` 降到 `~0.039`，说明 \(C\)-sector 更新是主瓶颈。

### 3.3 每步局部重解 \(C\)，加 trace0

命令：

```bash
python3 kg_examples/simulate_d_harmonic_standalone.py \
  --package visualizations/equation_first_gbcd_plus_initial_package_n384_tau0_core10_noactiveedge_guarded/gbcd_plus_initial_package.npz \
  --output visualizations/d_harmonic_standalone_tau0_steps2_localC_trace0 \
  --tau 0 \
  --steps 2 \
  --evolve-region core10 \
  --auxiliary-update local \
  --auxiliary-trace0 \
  --metric-corrector none
```

结果：

- step 0 residual wmean core10：`0.0879`；
- step 2 residual wmean core10：`0.0970`；
- step 2 residual p95 core10：`0.355`；
- step 2 \(\rho\) pullback L1 core10：`9.18e-5`。

解释：

trace0 局部闭合会增大代数残差，但仍远好于冻结 \(C\)。这和前面的高分辨率静态诊断一致：trace0 不能逐点后处理，必须和 full conservation 联立。

## 4. 本轮直接回答

当前最小 standalone D 原型已经可以从 D 初值推进，不再依赖 A 后续快照。

在 \(\tau=0\)、`n=384, core10`、10 个小步内：

\[
\rho_{\rm D\to A}
\]

与 A 支的差异仍小：

\[
\mathrm{L1}_{core10}\approx4.5\times10^{-4}.
\]

但这不是最终 D 支预测，因为完整 D 方程残差没有被压到初始水平。最好的当前原型是局部 \(C\) 更新，无 trace0，残差约：

\[
\mathrm{residual}_{core10}\approx3.9\times10^{-2}.
\]

加 trace0 后约：

\[
\mathrm{residual}_{core10}\approx9.7\times10^{-2}.
\]

冻结 \(C\) 则约：

\[
\mathrm{residual}_{core10}\approx5.3\times10^{-1}.
\]

因此，关键缺口已经定位为：

\[
\boxed{
\text{\(C\)-sector 必须每步按 trace0 + full conservation 全局更新，而不能冻结或逐点局部重解。}
}
\]

## 5. 下一步

下一步不应该继续拉长当前原型的时间窗。

应实现：

\[
\boxed{
\text{per-step global \(C\)-update: \(C\in E,\ \mathrm{tr}C=0,\ \tilde\nabla C=0\)}
}
\]

并把它接入 `simulate_d_harmonic_standalone.py`。

技术路线：

1. 每步根据当前 \((\tilde g,\tilde\rho,u)\) 重建 \(E\) basis；
2. 用 trace0 nullspace 消元；
3. 用当前三层时间数据构造 \(\tilde\nabla^\mu C_{\mu\nu}=0\) 的全局稀疏约束；
4. 用 ALM 或更稳的 saddle solve 求 \(C\)；
5. 再做 metric plus corrector；
6. 检查 D 残差是否能回到 `1e-3` 以下，同时观察 D-A 物质偏差。

