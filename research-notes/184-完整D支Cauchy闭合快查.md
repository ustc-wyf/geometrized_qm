# 完整 D 支 Cauchy 闭合快查

日期：2026-05-08

## 问题

用户追问：当前几何候选

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu},
\qquad
C\in E,\quad
\mathrm{tr}_{\tilde g}C=0,\quad
\tilde\nabla^\mu C_{\mu\nu}=0
\]

再加物质方程

\[
\tilde g^{\mu\nu}u_\mu u_\nu=m^2,
\qquad
\tilde\nabla_\mu(\tilde\rho\,\tilde g^{\mu\nu}u_\nu)=0
\]

是否已经足够独立给出 D 支演化？

## 快速结论

在主非退化 patch 中，答案应改为：

\[
\boxed{
\text{方程在主部/计数意义上已经接近闭合；缺的不是“再加一条物质方程”，而是 gauge、约束传播和 patch 适定性证明。}
}
\]

更具体：

- 若采用投影写法，\(\Pi_E^\perp\mathcal R=0\) 给出 6 个 metric 主方程；
- harmonic gauge 给出 4 个坐标规范条件；
- 这 10 条主部条件仍有 1 个 \(E\)-方向 null mode；
- trace0 正是补这个 null mode 的标量状态方程；
- 物质 HJ + 连续性方程给出 \(S,\tilde\rho\) 的一阶演化；
- 因此在 \(\Delta\neq0,w^2\neq0\)、质量壳非退化、时间切片非特征的区域，当前系统已经可以作为局部 Cauchy 演化候选。

所以此前说“还不能独立演化”应理解为“还没证明适定/还没写成规范固定的生产级演化器”，不是说方程明显缺项。

## 投影写法下的独立方程

定义

\[
\mathcal R_{\mu\nu}
=
\tilde G_{\mu\nu}
-\tilde T_{\mu\nu}/M_P^2.
\]

主方程可写为：

\[
\Pi_E^\perp\mathcal R_{\mu\nu}=0,
\]

\[
\mathrm{tr}_{\tilde g}\Pi_E\mathcal R=0.
\]

在非退化 patch 中，第一条是 6 个独立二阶 metric 方程，第二条是 1 个补充标量主部条件。

注意：

\[
\tilde\nabla^\mu \Pi_E\mathcal R_{\mu\nu}=0
\]

在精确投影壳和物质守恒成立时，不是新的独立高阶方程。因为

\[
\mathcal R=\Pi_E\mathcal R,
\qquad
\tilde\nabla^\mu\mathcal R_{\mu\nu}
=
-\tilde\nabla^\mu\tilde T_{\mu\nu}/M_P^2
=0.
\]

但在数值求解/初值投影中，它必须硬施加，否则逐点投影代表元不自动满足导数相容性。

## 主部闭合

harmonic gauge 下，投影方程主部为

\[
-\frac12\Pi_E^\perp\tilde\square\bar h_{\mu\nu}.
\]

此前 `167` 已算出，\(\Pi_E^\perp\mathcal R=0\) 加 harmonic gauge 后仍有一个 \(E\)-方向 null mode：

\[
\bar h^{(0)}_{\mu\nu}
=
-w_\mu w_\nu,
\qquad
w_\mu=(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu.
\]

普通 trace 对这个模式的收缩是

\[
\tilde g^{\mu\nu}\bar h^{(0)}_{\mu\nu}=-w^2.
\]

因此只要

\[
w^2\neq0,
\]

trace0 就能看见并固定该主部模式。

这说明：**trace0 不是装饰性条件，而是局部 Cauchy 闭合所需的最后一个标量主部条件。**

## 物质方程是否够

物质方程为：

\[
u_\mu=\partial_\mu S,
\qquad
\tilde g^{\mu\nu}u_\mu u_\nu=m^2,
\]

\[
\tilde\nabla_\mu(\tilde\rho\,u^\mu)=0.
\]

在给定 \(\tilde g\) 和非特征时间切片后：

- HJ 方程可解出 \(\partial_t S\)，条件是质量壳判别式非负并选定分支；
- 连续性方程可给出 \(\partial_t\tilde\rho\)，条件是流的时间分量不退化；
- \(r_\mu=\partial_\mu\ln\sqrt{\tilde\rho}\) 因而由 \(\tilde\rho\) 的演化确定；
- \(\tilde T_{\mu\nu}\) 由 \(\tilde\rho,u,\tilde g\) 给出。

因此物质部分在 massive branch、非特征 patch 内是够的。

## 辅助场写法下的 \(C\) 演化

若显式写

\[
C_{\mu\nu}=\lambda_I E^I_{\mu\nu},
\]

则 \(\lambda_I\) 有 4 个。

trace0 给出一个代数关系，剩 3 个独立系数。

守恒方程

\[
\tilde\nabla^\mu C_{\mu\nu}=0
\]

对 \(\lambda_I\) 的主符号是

\[
M_{\nu I}(\xi)=\xi^\mu E^I_{\mu\nu}.
\]

此前已算出 \(M\) 泛型 rank 为 3，正好给出 3 个一阶传播条件；trace0 去掉的那个方向正是原来的 null direction。

所以在 \(w^2\neq0\) 的 patch 内，辅助场写法也在主部意义上闭合。

## 仍未完成的部分

当前不能直接声称“已经有完整 D 支演化器”，原因不是方程数量明显不足，而是以下内容还没完成：

1. 需要明确采用 harmonic gauge、ADM lapse-shift gauge 或其他规范；
2. 需要写出约束方程，并证明 harmonic/Einstein 约束传播；
3. 需要证明 trace0 与投影方程在演化中相容，而不是只在初始切片成立；
4. 需要处理 \(\Delta=0,r^\perp=0,w^2=0\) 的 patch transition；
5. 需要处理 massless/null branch；
6. 需要把数值器从“三切片相容求解”改成“给初值后推进”的 gauge-fixed solver。

## 当前判断

可以在文章中更积极地写：

\[
\boxed{
\text{The proposed equation-first system is locally closed at the level of principal-symbol counting in the non-degenerate massive patch.}
}
\]

但还应同时写：

\[
\boxed{
\text{A full well-posedness theorem and a stable standalone D-branch evolution code remain open.}
}
\]

中文表述：

当前方程组已经不只是“残差拟合形式”，而是一个在主 patch 内具备局部闭合希望的 Cauchy 候选系统；下一步工作是规范固定和约束传播证明，而不是再盲目寻找新的场方程。

