# 130. D 支三域规则可模拟性与 rho 映射错误审计

日期：2026-05-02

## 用户提出的三域规则

用户把当前 D 支模拟区域分成三类：

1. \(\tilde R\) 很小的区域；
2. \(f(\tilde R)\) 饱和区；
3. 过渡区。

该分类应使用

\[
y=\ell^2\tilde R
\]

作为无量纲判据。粗略地说：

- 弱曲率 bulk：\(|y|\ll1\)；
- 饱和 bulk：\(|y|\gg1\)；
- 过渡层/interface：\(|y|\sim1\)。

D 支候选函数为

\[
f(\tilde R)=\frac{\tanh(\ell^2\tilde R)}{\ell^2}
=\frac{\tanh y}{\ell^2},
\qquad
\phi=f_R=\operatorname{sech}^2 y.
\]

## 对三域规则的理论判断

该三域图像是合理的 singular-perturbation / matched-asymptotic 路线，可以作为数值模拟的组织原则，但需要补上精确定义和边界条件。

### 1. 弱曲率 bulk

当 \(|y|\ll1\) 时：

\[
f(\tilde R)\approx \tilde R,\qquad \phi\approx1,
\]

且若 \(y\) 在该区域内平滑，高阶导数项

\[
(\tilde g_{\mu\nu}\tilde\Box-\tilde\nabla_\mu\tilde\nabla_\nu)\phi
\]

近似为小量。场方程近似为普通 Einstein 方程：

\[
M_P^2\tilde G_{\mu\nu}\approx \tilde T_{\mu\nu}.
\]

因此在实验室弱引力、低密度边界条件下，可以把 \(\tilde g\) 选为接近平直的零曲率度规，或用线性化约束求一个近似零曲率解。但严格地说，只有当 \(\tilde T_{\mu\nu}/M_P^2\) 可忽略且边界/拓扑允许时，才能直接取零曲率。

### 2. 饱和 bulk

当 \(|y|\gg1\) 时：

\[
\phi\to0,\qquad \tilde\nabla\phi\to0.
\]

所以高阶导数项确实近似消失。但需要注意，\(f\) 本身并不是零，而是：

\[
f\to \pm \frac1{\ell^2}.
\]

因此场方程近似为：

\[
-\frac{M_P^2}{2}f\,\tilde g_{\mu\nu}
\approx
\tilde T_{\mu\nu}.
\]

只有当 \(M_P^2/\ell^2\) 在当前无量纲单位下也可忽略，或者该项被吸收到有效 cosmological/background 规则中，才可以把饱和 bulk 近似说成 \(0\approx0\)。这点后续必须显式检查，不能直接跳过。

如果正测度密度

\[
\tilde N=\sqrt{|\tilde g|}\tilde\rho
\]

在该区很小，或者 \(u^\mu\) 近零，则 \(\tilde T_{\mu\nu}\) 可近似小。但这个判断必须使用 \(\tilde\rho\) 或 \(\tilde N\)，不能使用原始 \(\rho_A\)。

### 3. 过渡区/interface

当 \(|y|\sim1\) 时，\(\phi=\operatorname{sech}^2y\) 快速变化，导数项不能当作 bulk 小量。若层宽远小于其他尺度，可以把过渡区压缩为界面 \(\Gamma\)，并对法向积分得到 jump/matching 条件。

这时模拟不应把该区当作普通网格点显式推进，而应作为 free-boundary/interface 问题：

- bulk 内分别解弱曲率/饱和近似方程；
- interface 上施加由 \((\tilde g_{\mu\nu}\tilde\Box-\tilde\nabla_\mu\tilde\nabla_\nu)\phi\) 积分得到的 jump law；
- interface 的位置、法向、速度、生成/消失由 matching 条件决定。

此前 trace-only jump law 只能作为诊断；主闭合应使用 direct tensor jump。

## 是否可以进入数值模拟

理论上可以，但不是直接把现有网格 PDE 硬推。应改成 bulk-interface 数值系统：

1. 状态变量使用 \(\tilde N=\sqrt{|\tilde g|}\tilde\rho\)、\(u_i\)、bulk \(\tilde g\) 或 scalaron/interface 变量；
2. 使用 \(y=\ell^2\tilde R\) 分类 bulk 与 interface；
3. 弱曲率 bulk 解近似 Einstein/线性化约束，必要时取近零曲率分支；
4. 饱和 bulk 解退化后的代数/边界约束，不能无条件写成 \(0=0\)，必须处理残留 \(-M_P^2 f\tilde g_{\mu\nu}/2\)；
5. 过渡层用 direct tensor jump/body-fitted interface matching；
6. 允许 interface 出现、合并、断裂、消失，这应由 level-set/body-fitted 几何和 jump solvability 决定；
7. 初值必须先做 \(\rho_A\to\tilde\rho_0\) 映射，并满足或投影到 D 支约束。

## rho 映射错误审计

用户指出上一轮 reduced 动态中没有处理 \(\rho\neq\tilde\rho\)。复查确认该错误存在，而且会影响多处解释。

当前混合变换、非零质量 \(\kappa=m^2\) 分支下，应使用：

\[
m^2\sqrt{|\tilde g|}\tilde\rho
=
X\sqrt{|g|}\rho
\]

或正测度密度口径：

\[
\sqrt{|\tilde g|}\tilde\rho
=
\frac{|X|}{m^2}\sqrt{|g|}\rho.
\]

在平直 \(g=\eta\) 初态中：

\[
\tilde\rho_0
=
\frac{|X_0|}{m^2}
\frac{\rho_A}{\sqrt{|\tilde g_0|}}.
\]

审计结果：

- `simulate_d_reduced_dynamic_same_initial.py` 中，初始 `n_cons` 错用 \(\rho_A\)，应改为 \(\tilde\rho_0\) 或 \(\tilde N_0\)。
- D 支 residual/jump/tensor 诊断中，凡把 A snapshot 的 `rho` 直接传给 `stress_tensor_tilde` 的结果，都错误使用了 \(\tilde T_{\mu\nu}\) 的密度源。
- `support = rho > frac*rho.max()` 这类掩膜使用的是 \(\rho_A\)，物理上应同时报告或改用 \(\tilde N\) 支撑区。
- 图中的 `rho_D-rho_A` 不是物理 observable，必须改成明确对象：裸 \(\tilde\rho\)、测度密度 \(\tilde N\)、流 \(J^\mu\)，或拉回到 \(g\) 表象的等效密度。
- `instant_transform` 不只是数值上病态，也概念上有问题：它把 D 演化变量当成原表象 \(\rho\) 去重构 \(Q\)，而正确 D 表象变量应是 \(\tilde\rho\) 或 \(\tilde N\)。

## 结论

用户三域规则在理论上可以支撑数值模拟，但需要把模拟器明确改造成：

\[
\text{weak bulk}
\;+\;
\text{saturated bulk}
\;+\;
\text{direct-tensor interface}
\]

的 bulk-interface 系统。

同时，之前所有涉及 D 支源项、密度比较、support mask 的数值结论必须在 \(\rho_A\to\tilde\rho\) 修正后重算；未重算前只能保留为几何/方法学诊断，不能作为物理结论。
