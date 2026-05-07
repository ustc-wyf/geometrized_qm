# 164-gBCD投影型Einstein-like方程候选

日期：2026-05-07

## 一句话结论

当前最值得推进的显式方程不是“逐点拟合 \(A,B,C,D\)”，而是一个投影型 Einstein-like 方程：

\[
\boxed{
\Pi^\perp_E
\left[
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right]
=0,
}
\]

并配合

\[
\boxed{
\tilde\nabla^\mu
\left[
\Pi_E
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)
\right]
=0,
}
\]

以及 \(Q\to0\) 时投影部分消失的分支条件。这里 \(\Pi_E\) 是投影到

\[
E^I_{\mu\nu}
=
\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}
\]

所张成的张量子空间上，\(\Pi^\perp_E=1-\Pi_E\)。

这就是 gBCD 的更一般、更几何化的写法。

## 1. 基本变量与定义

变换后几何变量：

\[
\tilde g_{\mu\nu}.
\]

物质变量：

\[
\rho,\qquad S.
\]

定义两个协变方向：

\[
u_\mu=\partial_\mu S,
\]

\[
r_\mu=\tilde\nabla_\mu\ln\sqrt{\rho}
=\frac{\tilde\nabla_\mu\sqrt{\rho}}{\sqrt{\rho}}.
\]

如果后续坚持使用原 \(g\) 表象的量子势，也可以把 \(r_\mu\) 换成由原变换拉回的振幅方向；但在 \(\tilde g\) 表象中，最自然的局域写法是上式。

张量基底：

\[
E^0_{\mu\nu}=\tilde g_{\mu\nu},
\]

\[
E^1_{\mu\nu}=u_\mu u_\nu,
\]

\[
E^2_{\mu\nu}=r_\mu r_\nu,
\]

\[
E^3_{\mu\nu}=u_{(\mu}r_{\nu)}
=\frac12(u_\mu r_\nu+r_\mu u_\nu).
\]

记

\[
\lambda_I=(A,B,C,D),
\]

\[
\mathcal C_{\mu\nu}
=
\lambda_I E^I_{\mu\nu}
=
A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C r_\mu r_\nu
+D u_{(\mu}r_{\nu)}.
\]

## 2. 从 gBCD 到投影型方程

原 gBCD 候选写成

\[
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu}.
\]

令 Einstein 残差为

\[
\mathcal R_{\mu\nu}
\equiv
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

则 gBCD 等价于要求

\[
\mathcal R_{\mu\nu}
\in
\mathrm{span}
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

也就是：

\[
\boxed{
\Pi^\perp_E \mathcal R_{\mu\nu}=0.
}
\]

这比“先拟合 \(A,B,C,D\)”更像一个真正场方程，因为它直接约束 \(\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2\) 的张量方向。

## 3. 投影算子的显式定义

给对称二阶张量空间选一个局域内积。最简单选择是

\[
\langle X,Y\rangle
=
X_{\mu\nu}Y_{\rho\sigma}
W^{\mu\nu\rho\sigma},
\]

\[
W^{\mu\nu\rho\sigma}
=
\frac12
(
\tilde g^{\mu\rho}\tilde g^{\nu\sigma}
+
\tilde g^{\mu\sigma}\tilde g^{\nu\rho}
).
\]

也可以改用 DeWitt 型超度规；这会改变代表元规范，但不改变“残差落在 \(E^I\) 子空间”的核心思想。

定义 Gram 矩阵：

\[
H_{IJ}
=
\langle E_I,E_J\rangle.
\]

若 \(H_{IJ}\) 非退化，则投影系数为

\[
\lambda_I
=
(H^{-1})_{IJ}
\langle E^J,\mathcal R\rangle.
\]

投影张量为

\[
(\Pi_E\mathcal R)_{\mu\nu}
=
\lambda_I E^I_{\mu\nu}.
\]

因此主方程可以完全写成

\[
\boxed{
\mathcal R_{\mu\nu}
-
E^I_{\mu\nu}
(H^{-1})_{IJ}
\langle E^J,\mathcal R\rangle
=0.
}
\]

这就是无显式自由 \(\lambda\) 的写法。

若 \(H_{IJ}\) 近退化，则需要 patch 分支或 Moore-Penrose 伪逆，并把近退化面作为表象图册边界处理。这正对应数值中出现的 `guard` 问题。

## 4. 守恒约束与测地线

为了让物质仍然给出 \(\tilde g\) 测地线，要求修正张量单独守恒：

\[
\boxed{
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
}
\]

投影写法中即

\[
\boxed{
\tilde\nabla^\mu
(\Pi_E\mathcal R)_{\mu\nu}
=0.
}
\]

展开后为

\[
\begin{aligned}
0={}&
\partial_\nu A
+u_\nu \tilde\nabla_\mu(Bu^\mu)
+B u^\mu \tilde\nabla_\mu u_\nu \\
&+r_\nu \tilde\nabla_\mu(Cr^\mu)
+C r^\mu \tilde\nabla_\mu r_\nu \\
&+\frac12 r_\nu \tilde\nabla_\mu(Du^\mu)
+\frac12D u^\mu \tilde\nabla_\mu r_\nu \\
&+\frac12 u_\nu \tilde\nabla_\mu(Dr^\mu)
+\frac12D r^\mu \tilde\nabla_\mu u_\nu .
\end{aligned}
\]

这不是额外随意规则，而是由 Bianchi identity 给出的自然条件：

\[
\tilde\nabla^\mu\tilde G_{\mu\nu}=0.
\]

若主方程成立且 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\)，则

\[
\tilde\nabla^\mu\tilde T_{\mu\nu}=0.
\]

对 Madelung dust-like 能动张量，例如

\[
\tilde T_{\mu\nu}=\alpha_m\tilde\rho\,u_\mu u_\nu
\]

其中 \(\alpha_m\) 只取决于物质作用量归一化。结合连续性方程

\[
\tilde\nabla_\mu(\tilde\rho u^\mu)=0,
\]

可得

\[
u^\mu\tilde\nabla_\mu u_\nu=0.
\]

因此隐变量粒子轨迹是 \(\tilde g\) 测地线。

## 5. 物质方程

物质部分仍使用前面已完成等价性分析得到的变换后物质方程。

标准 massive KG 归一化下：

\[
\boxed{
\tilde g^{\mu\nu}u_\mu u_\nu=m^2.
}
\]

若使用更一般壳参数 \(\kappa\)，则写作

\[
\tilde g^{\mu\nu}u_\mu u_\nu=\kappa,
\]

但要回到 massive KG，应取 \(\kappa=m^2\)。

连续性方程：

\[
\boxed{
\tilde\nabla_\mu(\tilde\rho u^\mu)=0.
}
\]

密度与体积元变换沿用此前约定，例如在当前 massive 分支中使用

\[
\sqrt{|\tilde g|}\tilde\rho
=
\text{由原 }(\rho,S,g)\text{ 和变换规则确定的正测度密度}.
\]

这部分不是 gBCD 新加的，而是物质部分等价性已经完成的结果。

## 6. \(Q\to0\) 回 Einstein 的分支条件

投影型方程本身只要求 Einstein 残差落在 \(E\) 子空间，它还不自动要求残差为零。

因此必须加一个分支/极限条件：

\[
\boxed{
Q\to0
\quad\Longrightarrow\quad
\lambda_I\to0,
}
\]

等价地：

\[
\boxed{
Q\to0
\quad\Longrightarrow\quad
\Pi_E\mathcal R_{\mu\nu}\to0.
}
\]

于是

\[
\tilde G_{\mu\nu}
\to
\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

在 \(\tilde g\to g\) 且 \(Q\to0\) 的经典/弱量子势极限下，这就回到普通 Einstein 方程。

若要把这个条件写成方程，而不是只写成边界条件，可以引入 \(Q\)-门控辅助泛函：

\[
\mathcal A_\lambda
=
\frac12
\int d^4x\sqrt{|\tilde g|}
\left[
\mu(Q)M^{IJ}\lambda_I\lambda_J
+K^{IJ\alpha\beta}
\tilde\nabla_\alpha\lambda_I
\tilde\nabla_\beta\lambda_J
\right],
\]

其中

\[
\mu(Q)\to+\infty
\quad(Q\to0).
\]

这样最小化会在 \(Q\to0\) 时选择 \(\lambda_I\to0\) 的 Einstein 分支。

一个协变的 \(Q\)-like 标量可取

\[
\mathcal Q_{\tilde g}
=
\frac{\tilde\square\sqrt{\rho}}{\sqrt{\rho}}
=
\tilde\nabla_\mu r^\mu+r_\mu r^\mu.
\]

在数值实现中可以用带 regulator 的

\[
\mu(\mathcal Q)
=
\mu_0
\left(
\frac{Q_0^2}{\mathcal Q^2+\epsilon^2Q_0^2}
\right)^p,
\]

理论极限则取 \(\epsilon\to0\)。

## 7. 最小显式方程组

当前最小可写理论候选为：

\[
\boxed{
\Pi^\perp_E
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)=0,
}
\]

\[
\boxed{
\tilde\nabla^\mu
\Pi_E
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)=0,
}
\]

\[
\boxed{
\tilde g^{\mu\nu}\partial_\mu S\partial_\nu S=m^2,
\qquad
\tilde\nabla_\mu(\tilde\rho\tilde g^{\mu\nu}\partial_\nu S)=0,
}
\]

\[
\boxed{
Q\to0\Rightarrow
\Pi_E
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)
\to0.
}
\]

这组方程的含义：

- 第一式：几何场方程的非 gBCD 张量方向必须消失；
- 第二式：gBCD 修正张量单独守恒，从而保证物质测地线解释；
- 第三式：物质部分仍是已证明等价的变换后 Madelung/KG 方程；
- 第四式：选择回到 Einstein 的经典/弱量子势分支。

## 8. 与数值结果的关系

三切片 guarded plus-only 数值结果支持的是：

- 对当前高斯干涉试验场，\(\mathcal R_{\mu\nu}\) 的主要非 Einstein 修正确实可以被 \(\tilde g_{\mu\nu},uu,rr,ur\) 这四个方向解释；
- 单独 BCD 不够，必须保留 \(A\tilde g_{\mu\nu}\)；
- 近退化 \(H_{IJ}\) 或近退化 \(\tilde g\) 点需要 patch/branch 处理；
- 数值中的 `auto_bad_zero` 应理解为投影方程在退化图册边界上的临时实现，而不是物理方程的一部分。

## 9. 仍未完成的地方

1. 投影内积 \(W^{\mu\nu\rho\sigma}\) 是否唯一？

不唯一。不同 \(W\) 对应不同代表元规范。需要通过 \(Q\to0\)、稳定性、残差和测地线要求筛选。

2. \(Q\)-门控辅助泛函是否必须？

若只把 \(Q\to0\) 写成边界条件，可以暂时不引入。但若要得到唯一分支和稳定演化，最好引入。

3. 是否能来自作用量？

还没有证明。这个投影型方程是 equation-first 候选。后续要检查 Helmholtz/self-adjoint integrability。

4. 是否已证明普适？

没有。当前只说明它是高斯波包干涉试验场上未被否定、且比 pure \(f(R)\) 路线更有希望的显式方程形式。
