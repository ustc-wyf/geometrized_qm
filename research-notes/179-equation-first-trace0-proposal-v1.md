# equation-first trace0 几何方程候选 v1

日期：2026-05-08

## 一句话版本

当前最小可写出的候选理论不是一个已知 action 的 Euler 方程，而是一个 equation-first 的 Einstein-like 系统：

\[
\boxed{
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+C_{\mu\nu},
\qquad
C_{\mu\nu}\in E,
\qquad
C^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0.
}
\]

其中

\[
E=
\mathrm{span}
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

这组方程的物理含义是：变换后几何不必满足普通 Einstein 方程，而允许一个由相位流方向 \(u\) 与振幅梯度方向 \(r\) 张成的、无迹且守恒的有效几何应力 \(C_{\mu\nu}\)。它是当前高斯干涉算例没有否定、且能兼容隐变量测地线解释的最小显式候选。

## 1. 基本变量

在 \(\tilde g\) 表象中取基本变量：

\[
\tilde g_{\mu\nu},\qquad S,\qquad \tilde\rho.
\]

定义：

\[
u_\mu=\partial_\mu S,
\qquad
r_\mu=\tilde\nabla_\mu\ln\sqrt{\tilde\rho}.
\]

重要约定：

\[
\boxed{\rho\neq\tilde\rho.}
\]

和 A 支或实验室 \(g\) 表象比较时，必须把 \(\tilde\rho\) 通过变换关系拉回 \(\rho\)，不能直接把 \(\tilde\rho\) 当作原表象密度。

## 2. 物质方程与测地线

物质部分沿用已经证明等价的 \(\tilde g\) 表象 Madelung/KG 方程：

\[
\boxed{
\tilde g^{\mu\nu}u_\mu u_\nu=m^2.
}
\]

\[
\boxed{
\tilde\nabla_\mu(\tilde\rho\,\tilde g^{\mu\nu}u_\nu)=0.
}
\]

由于 \(u_\mu=\partial_\mu S\)，有

\[
\tilde\nabla_\mu u_\nu=\tilde\nabla_\nu u_\mu.
\]

对质量壳方程取梯度：

\[
\tilde\nabla_\alpha(u^\mu u_\mu)=0.
\]

得到：

\[
u^\mu\tilde\nabla_\mu u_\alpha=0.
\]

所以在 \(m\neq0\) 的 timelike branch 中，积分曲线是 \(\tilde g\) 的测地线。无质量极限需改用 null Hamilton-Jacobi 分支单独处理。

这里的关键点是：引力方程中的 \(C_{\mu\nu}\) 不能直接改写上述物质变分方程。当前 equation-first 版本把物质方程作为独立方程列入系统；若未来 action 化，则辅助 sector 必须避免对 \(S,\tilde\rho\) 产生额外变分源，或通过 shadow fields/补偿场实现。

## 3. 几何方程

定义 Einstein 残差：

\[
\mathcal R_{\mu\nu}
\equiv
\tilde G_{\mu\nu}
-
\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

候选场方程是：

\[
\boxed{
\mathcal R_{\mu\nu}=C_{\mu\nu}.
}
\]

并要求：

\[
C_{\mu\nu}
=
A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C_r r_\mu r_\nu
+D u_{(\mu}r_{\nu)}.
\]

为避免与张量 \(C_{\mu\nu}\) 混淆，这里把 \(r_\mu r_\nu\) 方向的系数记为 \(C_r\)。

闭合条件：

\[
\boxed{
\tilde g^{\mu\nu}C_{\mu\nu}=0.
}
\]

\[
\boxed{
\tilde\nabla^\mu C_{\mu\nu}=0.
}
\]

第一条是当前数值支持的最小标量状态方程；第二条来自 Bianchi 恒等式与物质守恒的相容性，不是任意数值规则。

## 4. 投影等价写法

在 \(E\) 非退化的 patch 内，可以定义投影 \(\Pi_E\)。方程组可写为：

\[
\boxed{
\Pi_E^\perp\mathcal R_{\mu\nu}=0.
}
\]

\[
\boxed{
\tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0.
}
\]

\[
\boxed{
\tilde\nabla^\mu\Pi_E\mathcal R_{\mu\nu}=0.
}
\]

这个写法更接近“方程结构”：第一式说 Einstein 残差只能落在 \(E\) 子空间；第二式选择 trace0 状态分支；第三式保证该有效修正张量守恒。

## 5. patch 与分支条件

投影良定义需要 Gram 矩阵非退化。对内积 \(\langle X,Y\rangle=X_{\mu\nu}Y^{\mu\nu}\)，此前得到：

\[
\det H
=
\frac{d-2}{2}
\left(
u^2r^2-(u\cdot r)^2
\right)^3.
\]

因此主 patch 条件是：

\[
\boxed{
d>2,
\qquad
\Delta\equiv u^2r^2-(u\cdot r)^2\neq0.
}
\]

还需单独处理：

- \(r^\perp=0\) 或 \(u,r\) 线性相关的退化区域；
- trace 主符号可能漏掉的 \(w^2=0\) 面；
- 严格 \(1+1d\)，因为 \(\tilde g_{\mu\nu}\) 方向不再独立；
- \(m=0\) null branch；
- \(Q\to0\) 的 Einstein 分支。

当前最自然的经典/弱量子分支条件是：

\[
\boxed{
Q\to0
\quad\Longrightarrow\quad
C_{\mu\nu}\to0.
}
\]

这不是由 trace0 自动推出的，必须作为分支、边界或正则性选择加入。

## 6. 当前数值支持

在 1550nm 高斯波包干涉三切片、`n=384`、`core10`、窗口 `[-9,9]um` 上，trace0 先逐点 nullspace 硬消元，再联立 full conservation。

增广拉格朗日 hard conservation 结果：

| tau | penalty full divergence | ALM full divergence | penalty algebraic | ALM algebraic |
|---:|---:|---:|---:|---:|
| -3.5 | \(7.740\times10^{-3}\) | \(1.987\times10^{-3}\) | \(5.478\times10^{-3}\) | \(6.009\times10^{-3}\) |
| 0 | \(4.576\times10^{-4}\) | \(2.388\times10^{-5}\) | \(6.747\times10^{-4}\) | \(6.747\times10^{-4}\) |
| +3.5 | \(6.502\times10^{-3}\) | \(1.940\times10^{-3}\) | \(7.698\times10^{-3}\) | \(8.172\times10^{-3}\) |

解释：

- `hard-project` 不可作为代表元选择，因为它可能只满足守恒而破坏代数场方程；
- `hard-kkt` 在小规模有效，但当前 LSQR saddle-point 版本在左侧分离态病态；
- `hard-alm` 最稳，说明 trace0 + full conservation 可以接近硬约束实现，且不需要明显牺牲代数壳。

这支持当前 proposal，但不构成证明。它只说明这个方程候选没有在当前高分辨率干涉试验场上被直接排除。

## 7. 当前不声称的内容

本 proposal 不声称：

- 已经找到对应的局域 action；
- 它是纯 \(\tilde g\) 几何 action 的变分方程；
- 它已经给出完整 Cauchy 适定性证明；
- 它已经解决所有 patch/branch 退化面；
- 它已经证明所有实验场景都与平直量子力学一致。

已有 Helmholtz 检查表明，裸 trace0 投影方程不宜直接声称来自普通 pure-metric 局域作用量。若要 action 化，更可能需要真实辅助场 \(\chi\)-sector 的有效应力，而不是把 \(C_{\mu\nu}\) 当成一个可自由变分的裸张量源。

## 8. 下一步

最重要的理论任务：

1. 写出 \(Q\to0\Rightarrow C\to0\) 的具体分支/边界/正则性条件。
2. 检查约束传播：若初始切片满足 \(\Pi_E^\perp\mathcal R=0\)、trace0 和 \(\tilde\nabla C=0\)，演化是否保持。
3. 在 patch 交界面写清换图规则，特别是 \(\Delta=0,r^\perp=0,w^2=0\)。
4. 若继续 action 化，放弃最小 stealth-state ansatz，转向真实辅助场 \(\chi\)-sector。
5. 扩展数值反例测试：单高斯、洛伦兹线形、不同质量/无质量分支，以及 \(Q\to0\) 极限。
