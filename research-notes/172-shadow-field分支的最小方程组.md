# 172-shadow-field分支的最小方程组

日期：2026-05-08

## 一句话结论

shadow-field 分支可以解决一个关键问题：它能让辅助几何 sector 使用类似

\[
u_\mu,\qquad r_\mu
\]

的方向，同时不让辅助作用量直接变分物质变量

\[
S,\qquad \tilde\rho.
\]

但它会引入新自由度，所以不能把它当作免费胜利。它真正可行的条件是：

\[
\boxed{
U_\mu,R_\mu
\text{ 不能通过 action 里的硬约束强行等于 }u_\mu,r_\mu；
}
\]

而应通过同型动力学、初始/边界分支选择、以及唯一性来保证在物理分支上

\[
\boxed{
U_\mu=u_\mu,
\qquad
R_\mu=r_\mu.
}
\]

这条路的地位是：

\[
\boxed{
\text{它是目前最有希望保护测地线物质方程的 action 化路线，}
\text{但代价是引入 shadow 初值分支。}
}
\]

## 1. 名词定义：什么是 shadow field

这里的 shadow field 是一组独立辅助场：

\[
\Phi,\qquad \sigma.
\]

它们生成两个辅助方向：

\[
U_\mu=\partial_\mu\Phi,
\]

\[
R_\mu=\tilde\nabla_\mu\ln\sqrt{\sigma}.
\]

本笔记默认：

\[
\sigma>0,
\]

并先讨论 massive branch：

\[
m\neq0.
\]

无质量分支需要把后面的 \(U^2=m^2\) 改写成 null 分支条件，不能机械照搬。

注意：

\[
\boxed{
\Phi,\sigma
\text{ 不是物质场 }S,\tilde\rho。
}
\]

物质场仍是：

\[
S,\qquad \tilde\rho,
\]

并且由物质 action 给出：

\[
u_\mu=\partial_\mu S,
\]

\[
r_\mu=\tilde\nabla_\mu\ln\sqrt{\tilde\rho}.
\]

shadow field 的作用只是给辅助几何 sector 提供一个张量子空间：

\[
E_{\rm sh}
=
\mathrm{span}
\{
\tilde g_{\mu\nu},
U_\mu U_\nu,
R_\mu R_\nu,
U_{(\mu}R_{\nu)}
\}.
\]

这样辅助 action 可以依赖 \(\Phi,\sigma\)，但不直接依赖 \(S,\tilde\rho\)。于是物质变分不会被辅助 sector 直接改写。

## 2. 为什么不能用硬约束 \(U=u,R=r\)

最直观的想法是往 action 里加：

\[
\int\sqrt{|\tilde g|}
\left[
A^\mu(U_\mu-u_\mu)
+B^\mu(R_\mu-r_\mu)
\right].
\]

这会强制：

\[
U_\mu=u_\mu,
\qquad
R_\mu=r_\mu.
\]

但这正是不能做的事。

因为对 \(S\) 变分时，第一项给出：

\[
\delta_S
\int\sqrt{|\tilde g|}[-A^\mu\partial_\mu S]
=
\int\sqrt{|\tilde g|}
(\tilde\nabla_\mu A^\mu)\delta S.
\]

于是 \(S\) 的方程会变成：

\[
\mathcal E_S^{(m)}
+\tilde\nabla_\mu A^\mu
=0.
\]

除非

\[
\tilde\nabla_\mu A^\mu=0
\]

自动成立，否则 Hamilton-Jacobi/连续性结构会被改写。

同理，\(R=r\) 的硬约束会对 \(\tilde\rho\) 变分产生源项。

所以：

\[
\boxed{
\text{硬约束 }U=u,\ R=r
\text{ 会把 shadow 分支退回到上一轮已经发现的问题。}
}
\]

如果乘子 \(A^\mu,B^\mu\) 在壳上全为 0，物质方程不被改写；但那时硬约束又失去 enforce 能力。因此这条路不稳。

## 3. 正确策略：同型动力学 + 分支初值

正确策略不是在 action 中硬连 \(U,R\) 和 \(u,r\)，而是让 shadow fields 满足和物质场同型的动力学：

\[
\boxed{
U^2=m^2,
}
\]

\[
\boxed{
\tilde\nabla_\mu(\sigma U^\mu)=0.
}
\]

物质场满足：

\[
\boxed{
u^2=m^2,
}
\]

\[
\boxed{
\tilde\nabla_\mu(\tilde\rho u^\mu)=0.
}
\]

然后选择物理分支的初值：

\[
\boxed{
\Phi|_{\Sigma}=S|_{\Sigma}+{\rm const},
\qquad
\partial_n\Phi|_{\Sigma}=\partial_n S|_{\Sigma},
\qquad
\sigma|_{\Sigma}=\tilde\rho|_{\Sigma}.
}
\]

若这两个系统在同一个 \(\tilde g\) 上适定且解唯一，则后续演化保持：

\[
\Phi=S+{\rm const},
\qquad
\sigma=\tilde\rho,
\]

因此：

\[
\boxed{
U_\mu=u_\mu,
\qquad
R_\mu=r_\mu.
}
\]

这是一种“分支选择”，不是 action 中的硬约束。它类似选择同一个解支，而不是在方程里额外加一个力把两者绑住。

## 4. 最小 action 原型

总作用量写成：

\[
S_{\rm total}
=
S_{\rm EH}[\tilde g]
+S_m[\tilde g,\tilde\rho,S]
+S_{\rm aux}[\tilde g,C,\Phi,\sigma,\ldots].
\]

物质 action 只含：

\[
\tilde g,\tilde\rho,S.
\]

辅助 action 不直接含：

\[
S,\qquad \tilde\rho.
\]

它只含：

\[
\tilde g,\quad C_{\mu\nu},\quad \Phi,\quad \sigma,
\quad
\text{乘子或状态变量}.
\]

最小结构可写为：

\[
S_{\rm aux}
=
S_{\rm source}[\tilde g,C]
+S_{\rm state}[C,\tilde g,U,R,Q_{\rm sh}]
+S_{\rm sh}[\tilde g,\Phi,\sigma].
\]

其中：

\[
S_{\rm source}
=
-\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\tilde g^{\mu\nu}C_{\mu\nu}.
\]

shadow 量子势定义为：

\[
Q_{\rm sh}
=
\frac{\tilde\square\sqrt{\sigma}}{\sqrt{\sigma}}.
\]

在物理分支上：

\[
\sigma=\tilde\rho
\quad\Longrightarrow\quad
Q_{\rm sh}=Q.
\]

## 5. 状态约束

辅助状态条件仍是：

\[
\boxed{
C_{\mu\nu}\in E_{\rm sh}.
}
\]

即：

\[
\Pi_{E_{\rm sh}}^\perp C_{\mu\nu}=0.
\]

无迹条件：

\[
\boxed{
C^\mu{}_\mu=0.
}
\]

守恒条件：

\[
\boxed{
\tilde\nabla^\mu C_{\mu\nu}=0.
}
\]

\(Q\to0\) 回 Einstein 的分支条件改写为：

\[
\boxed{
Q_{\rm sh}\to0
\quad\Longrightarrow\quad
C_{\mu\nu}\to0.
}
\]

在物理分支上 \(Q_{\rm sh}=Q\)，所以这就是原来需要的条件。

一个形式上的状态 action 是：

\[
S_{\rm state}
=
\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\left[
\alpha\,C^\mu{}_\mu
+\zeta^\nu\tilde\nabla^\mu C_{\mu\nu}
+\Lambda_\perp^{\mu\nu}
(\Pi_{E_{\rm sh}}^\perp C)_{\mu\nu}
+\mu(Q_{\rm sh})C_{\mu\nu}C^{\mu\nu}
\right].
\]

这里：

\[
\mu(Q_{\rm sh})\to\infty
\quad
\text{as}
\quad
Q_{\rm sh}\to0,
\]

用于选择 Einstein 分支。

## 6. shadow 动力学

为了让 shadow fields 可以和物质场选同一解支，最小 shadow 方程应复制物质 hydrodynamic 方程：

\[
U_\mu=\partial_\mu\Phi,
\]

\[
U^2=m^2,
\]

\[
\tilde\nabla_\mu(\sigma U^\mu)=0.
\]

可以用乘子形式写：

\[
S_{\rm sh}
=
\int\sqrt{|\tilde g|}
\left[
\eta(U^2-m^2)
+\beta\,\tilde\nabla_\mu(\sigma U^\mu)
\right].
\]

对 \(\eta\) 变分给出：

\[
U^2=m^2.
\]

对 \(\beta\) 变分给出：

\[
\tilde\nabla_\mu(\sigma U^\mu)=0.
\]

对 \(\Phi,\sigma\) 的变分则决定乘子 \(\eta,\beta\) 的传播。

这个 \(S_{\rm sh}\) 仍会贡献 metric stress，但这是辅助 sector 的一部分，应该被计入：

\[
\Theta_{\mu\nu}
=
-\frac{2}{\sqrt{|\tilde g|}}
\frac{\delta S_{\rm aux}}{\delta\tilde g^{\mu\nu}}.
\]

所以最终 metric 方程是：

\[
M_P^2\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}^{(m)}
+\Theta_{\mu\nu}.
\]

## 7. 这条路保护了什么

因为 \(S_{\rm aux}\) 不显含 \(S,\tilde\rho\)，所以：

\[
\frac{\delta S_{\rm aux}}{\delta S}=0,
\]

\[
\frac{\delta S_{\rm aux}}{\delta\tilde\rho}=0.
\]

于是物质方程仍由 \(S_m\) 单独给出：

\[
u^2=m^2,
\]

\[
\tilde\nabla_\mu(\tilde\rho u^\mu)=0.
\]

因此仍可推出：

\[
u^\mu\tilde\nabla_\mu u_\nu=0.
\]

也就是说：

\[
\boxed{
\text{shadow-field 分支可以保护隐变量粒子的测地线解释。}
}
\]

## 8. 它没有自动解决什么

这条路也有明显代价。

第一，它引入了新自由度：

\[
\Phi,\qquad \sigma.
\]

如果不选择物理分支，它们可以偏离物质场：

\[
U_\mu\neq u_\mu,
\qquad
R_\mu\neq r_\mu.
\]

这会导致辅助几何 sector 由 shadow matter 决定，而不是由真实物质波函数决定。

第二，分支选择现在是初始/边界条件：

\[
\Phi=S+{\rm const},
\qquad
\sigma=\tilde\rho.
\]

这不是 action 自动强制的结果。

第三，若 shadow equations 与物质 equations 在完整耦合 \(\tilde g\) 下不完全同型，等号不一定传播。

因此这条路的严格要求是：

\[
\boxed{
\text{必须证明物理分支 }(\Phi,\sigma)=(S,\tilde\rho)
\text{ 在完整耦合系统中是保持的。}
}
\]

## 9. 与原 equation-first 方程的关系

在物理分支上：

\[
U=u,
\qquad
R=r,
\qquad
Q_{\rm sh}=Q.
\]

于是：

\[
E_{\rm sh}=E.
\]

如果再有：

\[
\Theta_{\mu\nu}/M_P^2=C_{\mu\nu}
\]

或至少

\[
\Theta_{\mu\nu}/M_P^2
\in E,
\qquad
\Theta^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu\Theta_{\mu\nu}=0,
\]

则 metric 方程变为：

\[
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}^{(m)}
+\mathcal C_{\mu\nu},
\]

其中

\[
\mathcal C_{\mu\nu}\in E,
\qquad
\mathcal C^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
\]

这就回到了 trace=0 投影方程候选。

因此 shadow-field 分支的逻辑是：

\[
\boxed{
\text{用 action 保护物质方程；}
\quad
\text{用 shadow 分支复制 }u,r\text{ 的几何方向；}
\quad
\text{用辅助有效应力产生投影修正。}
}
\]

## 10. 当前判断

shadow-field 分支不是无代价的，但它比硬约束 \(U=u,R=r\) 更健康。

可以给出三层判断：

1. 它可以保护物质 HJ/连续性方程，因此保留测地线解释；
2. 它可以在物理分支上恢复原来的 \(u,r\) 投影结构；
3. 它引入额外初值分支，所以必须把“物理分支选择”作为理论的一部分。

因此当前结论是：

\[
\boxed{
\text{shadow-field 分支是 action 化的可行候选，}
\text{但它把难点转移到分支选择和额外自由度控制上。}
}
\]

## 11. 下一步

下一步有两个方向。

理论方向：

1. 检查 \(S_{\rm sh}\) 的完整 metric stress 是否仍落在 \(E_{\rm sh}\) 或可被 \(C_{\mu\nu}\) 状态方程吸收；
2. 检查 \((\Phi,\sigma)=(S,\tilde\rho)\) 的分支在完整耦合系统中是否传播；
3. 若传播成立，写出完整 shadow-auxiliary action v1。

数值方向：

1. 在 `n=384` 上复查 trace=0 投影方程；
2. 若 trace=0 仍优于 \(\mathsf q_E\)-trace，则在三切片上把 \(U,R\) 直接取为 A 支 \(u,r\) 的物理分支，先检验 shadow 版本的状态 residual；
3. 再考虑 off-branch 扰动，检查分支是否稳定。
