# 165-gBCD投影方程闭合条件

日期：2026-05-07

## 一句话结论

上一份笔记中的投影型方程可以进一步写得更清楚：

\[
\boxed{
\mathcal R_{\mu\nu}
\equiv
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\in
\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}.
}
\]

这个“属于某个张量子空间”的条件本身不依赖投影内积 \(W\)。\(W\) 只影响离壳时怎样定义最小残差、怎样选择代表元、以及数值诊断的权重。

真正的结构性条件是：

1. \(u_\mu,r_\mu\) 张成的二平面不能退化；
2. 有效维度不能是严格 \(1+1\) 维；
3. \(Q\to0\) 时投影残差必须选择消失分支；
4. 若要从作用量推出，还必须通过 Helmholtz/self-adjoint 检查。

## 1. 两种等价写法

定义

\[
\mathcal R_{\mu\nu}
=
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

### 写法 A：投影/成员关系形式

\[
\boxed{
\Pi_E^\perp \mathcal R_{\mu\nu}=0.
}
\]

等价于：

\[
\mathcal R_{\mu\nu}=\lambda_I E^I_{\mu\nu},
\]

其中

\[
E^I_{\mu\nu}
=
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

这个形式最好，因为它直接是一个张量方向约束。

### 写法 B：辅助应力形式

\[
\boxed{
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu},
}
\]

\[
\boxed{
\mathcal C_{\mu\nu}
=
A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C r_\mu r_\nu
+D u_{(\mu}r_{\nu)}.
}
\]

若把 \(A,B,C,D\) 当作辅助场，则还要加守恒/分支/代表元条件；若把它们直接定义为投影系数，则在非退化区由 \(\mathcal R_{\mu\nu}\) 唯一决定。

## 2. \(W\) 是否改变物理方程？

给对称二阶张量空间选内积

\[
\langle X,Y\rangle
=
X_{\mu\nu}Y_{\rho\sigma}W^{\mu\nu\rho\sigma}.
\]

只要 \(W\) 非退化，且 \(E^I\) 的 Gram 矩阵在该点非退化，投影算子 \(\Pi_E\) 满足：

\[
\Pi_E^\perp \mathcal R=0
\quad\Longleftrightarrow\quad
\mathcal R\in \mathrm{span}\{E^I\}.
\]

因此：

- 在壳方程 \(\Pi_E^\perp\mathcal R=0\) 的解集不依赖 \(W\)；
- 离壳残差大小依赖 \(W\)；
- 离壳拟合出的 \(\lambda_I\) 依赖 \(W\)；
- 若最终用一个辅助场 action 或范数最小原则选择 \(\lambda_I\)，则 \(W\) 会成为代表元规范的一部分。

这意味着 \(W\) 的正确定位是“诊断/规范选择”，不是投影方程本身的物理自由参数。

## 3. Gram 矩阵和退化条件

取最简单的张量内积

\[
\langle X,Y\rangle=X_{\mu\nu}Y^{\mu\nu}.
\]

定义三个标量：

\[
a=u^2=\tilde g^{\mu\nu}u_\mu u_\nu,
\]

\[
b=r^2=\tilde g^{\mu\nu}r_\mu r_\nu,
\]

\[
c=u\cdot r=\tilde g^{\mu\nu}u_\mu r_\nu.
\]

再定义二平面的 Gram 行列式：

\[
\Delta=ab-c^2.
\]

在 \(d\) 维时，基底

\[
\{\tilde g,uu,rr,ur\}
\]

的 Gram 行列式为

\[
\boxed{
\det H=\frac{d-2}{2}\,\Delta^3.
}
\]

所以：

- 在 \(3+1d\) 中，\(\det H=\Delta^3\)；
- 在 \(2+1d\) 中，\(\det H=\frac12\Delta^3\)；
- 在严格 \(1+1d\) 中，\(\det H=0\) 恒成立。

这给出一个非常清楚的判断：

\[
\boxed{
\text{投影型 gBCD 在 }d>2\text{ 且 }\Delta\neq0\text{ 的 patch 内是非退化的。}
}
\]

这里的 \(\Delta=0\) 不只表示 \(u,r\) 线性相关。在 Lorentzian 度规里，它也可能表示 \(u,r\) 张成的二平面本身是 null-degenerate 的。

这解释了两个现象：

1. 在 \(3+1d\) 或 \(2+1d\) 数值实验里，\(\{\tilde g,uu,rr,ur\}\) 可以作为四个独立方向；
2. 在严格 \(1+1d\) 里，\(\tilde g_{\mu\nu}\) 本来就可由 \(uu,rr,ur\) 表示，所以 \(A\tilde g_{\mu\nu}\) 不是独立方向，必须单独降维处理。

## 4. 守恒条件到底是额外方程还是一致性条件？

若采用辅助应力形式，需要写：

\[
\boxed{
\tilde\nabla^\mu \mathcal C_{\mu\nu}=0.
}
\]

这保证

\[
\tilde\nabla^\mu \tilde T_{\mu\nu}=0,
\]

因为

\[
\tilde\nabla^\mu\tilde G_{\mu\nu}=0.
\]

但是如果采用投影/成员关系形式，并且物质方程已经包含：

\[
\tilde g^{\mu\nu}u_\mu u_\nu=m^2,
\]

\[
\tilde\nabla_\mu(\tilde\rho u^\mu)=0,
\]

那么对 dust-like

\[
\tilde T_{\mu\nu}=\alpha_m\tilde\rho u_\mu u_\nu
\]

已经有

\[
\tilde\nabla^\mu\tilde T_{\mu\nu}=0.
\]

此时

\[
\tilde\nabla^\mu\mathcal R_{\mu\nu}=0
\]

由 Bianchi identity 自动成立。

所以更准确的说法是：

- 如果 \(\mathcal C_{\mu\nu}\) 是独立辅助场，则 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\) 是必须加的闭合方程；
- 如果 \(\mathcal C_{\mu\nu}\) 被定义为 \(\Pi_E\mathcal R_{\mu\nu}\)，且物质方程已经成立，则守恒是相容性条件和数值诊断；
- 它不是任意附加规则。

## 5. 测地线条件的两个来源

第一种来源是 Hamilton-Jacobi 壳方程。

因为

\[
u_\mu=\partial_\mu S,
\]

所以

\[
\tilde\nabla_\mu u_\nu=\tilde\nabla_\nu u_\mu.
\]

若

\[
u^2=m^2
\]

为常数，则

\[
u^\mu\tilde\nabla_\mu u_\nu
=
u^\mu\tilde\nabla_\nu u_\mu
=
\frac12\tilde\nabla_\nu(u^2)
=0.
\]

所以 HJ 壳方程本身已经给出 \(\tilde g\)-测地线。

第二种来源是应力守恒。

\[
\tilde\nabla^\mu(\tilde\rho u_\mu u_\nu)=0
\]

加上

\[
\tilde\nabla_\mu(\tilde\rho u^\mu)=0
\]

也推出

\[
u^\mu\tilde\nabla_\mu u_\nu=0.
\]

因此我们对引力侧的要求不是“再制造测地线”，而是不能破坏物质方程与 Bianchi identity 的一致性。

## 6. \(Q\to0\) 分支条件的显式写法

投影方程只要求

\[
\mathcal R_{\mu\nu}\in\mathrm{span}\{E^I\}.
\]

它不自动要求

\[
\mathcal R_{\mu\nu}=0.
\]

因此必须额外选择 Einstein 分支：

\[
\boxed{
Q\to0
\quad\Longrightarrow\quad
\lambda_I\to0.
}
\]

更一般地，可以写成：

\[
\boxed{
\lambda_I=\chi(\mathcal Q)\,\hat\lambda_I,
\qquad
\chi(0)=0,
\qquad
\hat\lambda_I \text{ regular}.
}
\]

其中一个纯 \(\tilde g\) 表象中可用的量子势型标量是

\[
\mathcal Q
=
\frac{\tilde\square\sqrt\rho}{\sqrt\rho}
=
\tilde\nabla_\mu r^\mu+r_\mu r^\mu.
\]

若未来能把原表象 \(Q_g\) 严格拉回成 \(\tilde g,\rho,S\) 的函数，则可用精确拉回的 \(Q_g\)。在只用 \(\tilde g\) 表象的当前方程中，\(\mathcal Q\) 是最自然的局域替代。

这一步非常关键：没有它，投影方程允许在 \(Q=0\) 时仍有非 Einstein 的各向异性有效源。

## 7. 与作用量的关系

并不是所有协变方程都一定来自局域作用量。

对一个候选 metric equation

\[
\mathcal E_{\mu\nu}[\tilde g,\rho,S]=0,
\]

它来自作用量的必要条件是线性化算子满足 Helmholtz/self-adjoint 条件：

\[
\int\sqrt{|\tilde g|}\,
h^{\mu\nu}\,
\delta\mathcal E_{\mu\nu}[k]
=
\int\sqrt{|\tilde g|}\,
k^{\mu\nu}\,
\delta\mathcal E_{\mu\nu}[h]
+\text{boundary}.
\]

对当前候选

\[
\mathcal E_{\mu\nu}
=
\Pi_E^\perp
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right),
\]

最高阶主部来自 \(\Pi_E^\perp\delta\tilde G_{\mu\nu}\)。

由于 \(\Pi_E^\perp\) 是代数投影，如果选正交投影，主部有希望保持自伴。但完整方程还包含

\[
\delta\Pi_E^\perp
\]

带来的低阶项，所以 action-level integrability 仍未证明。

当前判断：

- 投影型方程作为 equation-first 候选是明确可写的；
- 它的主部不像明显病态的非自伴方程；
- 但是否存在局域作用量，需要下一步显式检查 Helmholtz 条件；
- 若失败，仍可作为非作用量的有效场方程研究，但理论规格会低于 action-first 理论。

## 8. 下一步应做什么

下一步不应回到“先跑多步数值器”。

更合理的顺序是：

1. 对 \(\Pi_E^\perp\mathcal R=0\) 做 3+1d 方程计数与主部分析；
2. 对 \(\Delta=0\) 和 \(d=2\) 退化面写 patch/branch 条件；
3. 检查 \(Q\to0\) 分支条件是否能由 \(\lambda_I=\chi(\mathcal Q)\hat\lambda_I\) 稳定实现；
4. 做 Helmholtz/self-adjoint 条件的第一轮符号分析；
5. 再回到高斯干涉数值例子，只检验这些明确方程条件，而不是继续盲目优化演化器。

