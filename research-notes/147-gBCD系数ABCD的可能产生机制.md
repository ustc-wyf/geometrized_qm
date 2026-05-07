# 147. gBCD 系数 \(A,B,C,D\) 的可能产生机制

日期：2026-05-07

## 1. 问题

当前 equation-first 候选为

\[
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu},
\]

\[
\mathcal C_{\mu\nu}
=
A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C r_\mu r_\nu
+D u_{(\mu}r_{\nu)}.
\]

前一轮已经说明：这个张量壳可以在三时刻 hard constraint 下同时保持较低代数残差与守恒/测地线条件。但它还不是完整理论，因为必须说明 \(A,B,C,D\) 如何产生。

本文把“产生机制”分成三类。

## 2. 机制 I：简单局部本构函数

最简单想法是

\[
A=A(q),\quad B=B(q),\quad C=C(q),\quad D=D(q),
\]

其中 \(q\) 是局部标量集合，例如

\[
q=(\log\rho,\tilde R,I_2,u^2,r^2,u\cdot r).
\]

已实现检查：

`kg_examples/diagnose_gbcd_constitutive_coefficients.py`

输出：

`visualizations/equation_first_gbcd_constitutive_full_n96_3tau/`

结果：

- 线性/二次多项式不能稳定拟合；
- 留一时刻泛化很差；
- kNN 也不能稳定跨时刻预测；
- ridge/minimum-norm 正则可以改善一部分系数的极值和 kNN 表现，但 \(rr\) 方向仍很差。

当前判断：

\[
\boxed{
\text{简单局部标量函数机制目前不成立。}
}
\]

但这不是最终 no-go，因为 KKT 解仍有代表元选择问题，且当前特征集不含 \(\nabla u,\nabla r,\nabla R\) 等导数特征。

## 3. 机制 II：单个局域势函数 \(L\) 的度规变分

更强的机制是假设存在一个局域标量势函数

\[
L=L(\rho,U,V,W),
\]

其中

\[
U=u^2,\qquad V=r^2,\qquad W=u\cdot r.
\]

若

\[
S_L=\int d^dx\sqrt{|\tilde g|}\,L,
\]

且 \(L\) 不含度规导数，则其度规变分给出

\[
T^{(L)}_{\mu\nu}
=
L\tilde g_{\mu\nu}
-2L_U u_\mu u_\nu
-2L_V r_\mu r_\nu
-2L_W u_{(\mu}r_{\nu)}.
\]

因此 gBCD 形式会自然出现，并且系数不独立：

\[
A=L,\qquad
B=-2L_U,\qquad
C=-2L_V,\qquad
D=-2L_W,
\]

差一个整体符号约定。

这带来可检验的 integrability 条件：

\[
\partial_V B=\partial_U C,\qquad
\partial_W B=\partial_U D,\qquad
\partial_W C=\partial_V D.
\]

已实现检查：

`kg_examples/fit_gbcd_potential_constitutive.py`

输出：

- `visualizations/equation_first_gbcd_potential_matter_deg2_n96_3tau/`
- `visualizations/equation_first_gbcd_potential_matter_deg3_n96_3tau/`

结果：

- degree 2：最佳 global weighted residual `0.98984`；
- degree 3：最佳 global weighted residual `0.98085`；
- 单片 weighted residual 仍普遍远大于自由 gBCD hard constraint 的 `1%~5%` 水平。

当前判断：

\[
\boxed{
\text{低阶单势函数 }L(\rho,u^2,r^2,u\cdot r)
\text{ 不能解释当前 gBCD 源。}
}
\]

如果未来继续这个机制，必须至少加入导数不变量或曲率耦合，例如

\[
L=L(\rho,U,V,W,\nabla u,\nabla r,\tilde R,\ldots),
\]

但此时度规变分会产生更多导数张量项，不再只是纯 gBCD。

## 4. 机制 III：辅助各向异性应力场

当前最有希望的机制是把 \(A,B,C,D\) 看成辅助场，类似各向异性流体/弹性介质中的能量密度、两个方向压力和剪切应力。

定义张量基底

\[
E^I_{\mu\nu}
=
\{\tilde g_{\mu\nu},\,
u_\mu u_\nu,\,
r_\mu r_\nu,\,
u_{(\mu}r_{\nu)}\}.
\]

则

\[
\mathcal C_{\mu\nu}=\lambda_I E^I_{\mu\nu},
\qquad
\lambda_I=(A,B,C,D).
\]

场方程可以写成

\[
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
-\lambda_I E^I_{\mu\nu}=0.
\]

若给定 \(\tilde g,u,r,\tilde T\)，则 \(\lambda_I\) 可以通过投影近似确定：

\[
\lambda_I
=
(H^{-1})_{IJ}
\langle E^J,\,
\tilde G-\tilde T/M_P^2\rangle,
\]

其中

\[
H_{IJ}=\langle E_I,E_J\rangle.
\]

但这只是代数投影。真正闭合还需要守恒：

\[
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0.
\]

这就是 \(A,B,C,D\) 的一阶 PDE：

\[
\begin{aligned}
0={}&
\partial_\nu A
+u_\nu \tilde\nabla_\mu(Bu^\mu)
+B u^\mu \tilde\nabla_\mu u_\nu \\
&+r_\nu \tilde\nabla_\mu(Cr^\mu)
+C r^\mu \tilde\nabla_\mu r_\nu \\
&+\frac12 r_\nu \tilde\nabla_\mu(Du^\mu)
+\frac12 D u^\mu \tilde\nabla_\mu r_\nu \\
&+\frac12 u_\nu \tilde\nabla_\mu(Dr^\mu)
+\frac12 D r^\mu \tilde\nabla_\mu u_\nu .
\end{aligned}
\]

这个机制的含义是：

\[
\boxed{
A,B,C,D
\text{ 不是由简单函数直接给出，而是由“场方程投影 + 守恒 PDE + 边界/初值条件”共同产生。}
}
\]

这与流体力学很像：仅有应力张量形式还不够，还需要状态方程和守恒律共同决定演化。

## 5. 当前数值证据

自由 gBCD hard constraint：

- transverse hard 三时刻中心残差约 `1%~5%`；
- full divergence-free hard 三时刻中心残差也约 `1%~5%`；
- 守恒约束可压到近数值零。

这说明机制 III 在 fixed A-reference 背景上是可行的。

但相邻 probe 时间层残差在 full hard 下升到约 `14%~23%`，说明“每个切片独立投影”不能替代真实时间演化。真正 D 支动力学中，\(\tilde g,u,r,A,B,C,D\) 必须一起调整。

## 6. 下一步

下一步不应直接开全动力学，而应先补上机制 III 的闭合选择：

1. 选择 \(A,B,C,D\) 的代表元规范。

可选：

- 最小范数；
- 最小空间梯度；
- 最小时间梯度；
- 使 \(A,B,C,D\) 尽量接近某个局部本构函数。

2. 写出辅助场演化系统。

例如把

\[
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0
\]

作为 \(\lambda_I\) 的一阶约束/演化方程，再配合一个额外状态方程，例如 trace 关系、最小范数条件或由初态投影继承的条件。

3. 检查作用量重建。

若机制 III 能在数值上产生稳定演化，再检查是否存在带辅助场的作用量

\[
S_{\rm aux}[\tilde g,u,r,\lambda_I]
\]

使得对 \(\lambda_I\) 和 \(\tilde g\) 变分后恢复上述方程。

当前最稳妥的阶段性判断是：

\[
\boxed{
\text{gBCD 的产生机制最像“辅助各向异性应力场”，}
\text{而不是简单局部势函数。}
}
\]

补充：代表元规范与 hard constraint 的数值实现已在
`research-notes/148-gBCD辅助应力场代表元规范与nullspace硬约束.md`
中单独整理。该轮确认：KKT hard constraint 会在病态切片上产生数值泄漏；改用 nullspace 投影后，
full conservation 可被压到 \(10^{-16}\) 量级，同时中心代数残差仍保持在 `1%~5%`。
