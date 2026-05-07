# 170-trace0投影方程的Helmholtz检查

日期：2026-05-08

## 一句话结论

当前 trace=0 投影方程组作为 equation-first 候选仍然有价值，但它不太可能直接就是某个普通局域协变作用量对 \(\tilde g_{\mu\nu}\) 变分得到的 Euler-Lagrange 方程。

主要原因有三条：

1. 投影算子 \(\Pi_E\) 依赖 \(\tilde g,u,r\)，线性化时会出现 \(\delta\Pi_E\cdot\mathcal R\)；
2. 即使原始 Einstein 残差 \(\mathcal R_{\mu\nu}\) 来自作用量，投影后的 \(M\mathcal R\) 一般不再保持 Frechet 导数的形式自伴性；
3. 纯 metric 协变作用量在物质壳上有 4 个 Noether 恒等式，而“投影 6 条 + trace 1 条”看起来有 7 条独立 metric 约束，除非 trace 条件来自辅助场方程而不是 metric 方程。

因此当前最稳妥的判断是：

\[
\boxed{
\text{trace=0 投影方程应暂时保留为 equation-first 理论候选；}
}
\]

\[
\boxed{
\text{若要 action-first 化，大概率需要辅助场/乘子作用量，而不是普通纯 }\tilde g\text{ 局域作用量。}
}
\]

## 1. Helmholtz/self-adjoint 检查是什么意思

对一组微分方程

\[
F_A[q]=0,
\]

它是否来自某个作用量

\[
S[q]=\int L(q,\partial q,\ldots)
\]

不是自动成立的。

局域变分方程的必要条件是：方程的 Frechet 线性化

\[
D F
\]

必须在分部积分意义下形式自伴：

\[
\boxed{
D F = (D F)^\dagger.
}
\]

直观地说，若

\[
F_A=\frac{\delta S}{\delta q^A},
\]

那么

\[
\frac{\delta F_A}{\delta q^B}
=
\frac{\delta^2 S}{\delta q^B\delta q^A}
\]

必须是对称的。这就是“二阶变分可交换”的无限维版本。

Einstein 方程满足这个条件，因为

\[
\tilde G_{\mu\nu}-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\]

来自 Einstein-Hilbert 作用量加物质作用量。

问题是：我们把它投影以后，是否仍然满足这个条件。

## 2. 当前方程的自然 tensor 嵌入

当前最小候选为：

\[
\mathcal R_{\mu\nu}
\equiv
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu},
\]

\[
\Pi_E^\perp\mathcal R_{\mu\nu}=0,
\]

\[
\tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0.
\]

其中

\[
E=\mathrm{span}
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

为了做 Helmholtz 检查，需要先把“6 条投影方程 + 1 条 trace 方程”嵌成一个 tensor equation。一个自然写法是

\[
F_{\mu\nu}
=
\Pi_E^\perp\mathcal R_{\mu\nu}
+\Sigma_{\mu\nu}\,
\tilde g^{\alpha\beta}\Pi_E\mathcal R_{\alpha\beta},
\]

其中 \(\Sigma_{\mu\nu}\) 是 \(E\) 子空间内某个非退化代表方向，例如可先取 \(\Sigma_{\mu\nu}\propto\tilde g_{\mu\nu}\)。

只要 \(\Sigma_{\mu\nu}\) 不落入坏的零方向，

\[
F_{\mu\nu}=0
\]

就等价于

\[
\Pi_E^\perp\mathcal R=0,
\qquad
\mathrm{Tr}(\Pi_E\mathcal R)=0.
\]

但注意：\(\Sigma_{\mu\nu}\) 的选择本身不是原方程给定的。这已经说明 trace=0 投影方程组不像普通 Einstein 方程那样天然就是一个唯一的 metric Euler tensor。

## 3. 线性化结构

把上面的 tensor 嵌入记成

\[
F_{\mu\nu}
=
M_{\mu\nu}{}^{\alpha\beta}[\tilde g,u,r]\,
\mathcal R_{\alpha\beta}.
\]

其中 \(M\) 包含 \(\Pi_E^\perp\)、trace 代表方向 \(\Sigma\)、以及 \(\Pi_E\)。

线性化为：

\[
\delta F_{\mu\nu}
=
M_{\mu\nu}{}^{\alpha\beta}\,
\delta\mathcal R_{\alpha\beta}
+
\delta M_{\mu\nu}{}^{\alpha\beta}\,
\mathcal R_{\alpha\beta}.
\]

第一项是“投影后的 Einstein Hessian”。

第二项是新出现的投影变化项。

如果没有投影，即 \(M=1\)，则

\[
\delta F=\delta\mathcal R
\]

继承 Einstein-Hilbert 作用量的自伴性。

但现在 \(M\) 依赖场：

\[
M=M[\tilde g,u,r].
\]

其中投影矩阵包含：

\[
H^{-1}_{IJ},
E^I_{\mu\nu},
\qquad
\tilde g^{\mu\nu}.
\]

当它作用到 \(\mathcal R\) 上时，才产生收缩

\[
\langle E_I,\mathcal R\rangle.
\]

所以一般有：

\[
\boxed{
\delta M\neq0.
}
\]

这会产生

\[
\boxed{
(\delta\Pi_E)\mathcal R
}
\]

这样的项。它没有理由等于某个二阶变分的对称部分。

## 4. 有限维类比

这个问题可以用有限维类比看得很清楚。

设原方程来自势函数：

\[
R_i(q)=\frac{\partial V}{\partial q^i}.
\]

因此

\[
\partial_j R_i
\]

是对称矩阵。

现在定义投影方程：

\[
F_i=P_i{}^j(q)R_j.
\]

则

\[
\partial_k F_i
=
P_i{}^j\partial_kR_j
+
(\partial_kP_i{}^j)R_j.
\]

即使 \(\partial_kR_j\) 对称，新的 \(\partial_kF_i\) 也通常不对称，除非同时满足：

1. \(P\) 与 Hessian 特殊交换；
2. \((\partial P)R\) 的反对称部分消失。

这两个条件在当前投影方程中都没有自然保证。

对应到本项目：

\[
P\leftrightarrow M[\tilde g,u,r],
\qquad
R\leftrightarrow \mathcal R_{\mu\nu}.
\]

所以裸投影方程一般不会自动是 Euler-Lagrange 方程。

## 5. Einstein 极限为什么没有这个问题

若进入 Einstein 分支：

\[
Q\to0,
\qquad
\mathcal C_{\mu\nu}\to0,
\]

则

\[
\mathcal R_{\mu\nu}\to0.
\]

此时麻烦项

\[
(\delta M)\mathcal R
\]

也趋于 0。

所以在经典/低量子势极限，投影方程可以渐近恢复 Einstein-Hilbert 的变分结构。

这说明：

\[
\boxed{
\text{Helmholtz 问题主要发生在 } \mathcal C_{\mu\nu}\neq0 \text{ 的量子修正区。}
}
\]

这和项目目标并不矛盾。我们本来就希望在弱 \(Q\) 或实验可观测范围内接近普通 Einstein/平直量子力学，而在强量子几何区允许新结构出现。

## 6. Noether 恒等式的计数压力

如果存在普通局域协变作用量

\[
S[\tilde g,\tilde\rho,S],
\]

则 metric Euler 方程记为

\[
\mathcal E_{\mu\nu}=0.
\]

微分同胚不变性给出 Noether 恒等式：

\[
\tilde\nabla^\mu \mathcal E_{\mu\nu}
+
\mathcal E_S\partial_\nu S
+
\mathcal E_{\tilde\rho}\partial_\nu\tilde\rho
=0.
\]

在物质方程成立时，

\[
\mathcal E_S=0,
\qquad
\mathcal E_{\tilde\rho}=0,
\]

于是

\[
\tilde\nabla^\mu\mathcal E_{\mu\nu}=0
\]

成为 4 个恒等式。

因此普通 metric 方程在 \(3+1d\) 中通常只有

\[
10-4=6
\]

个独立组合。

但当前候选的 metric 部分是：

\[
\Pi_E^\perp\mathcal R=0
\]

给出 6 条，再加

\[
\tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0
\]

给出 1 条，总共像是 7 条独立 metric 条件。

这不是严格 no-go 定理，因为：

1. trace 条件可能与 matter 方程或 patch 条件在某些区域相关；
2. 可以引入辅助场，使 trace 条件来自辅助场变分而不是 metric 变分；
3. 非局域作用量或非标准约束理论可能绕开普通计数。

但它给出一个强烈信号：

\[
\boxed{
\text{trace=0 投影方程很难直接等同于普通纯 metric Euler 方程。}
}
\]

## 7. 守恒条件不是 Helmholtz 的替代品

我们已经要求：

\[
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
\]

这保证了和 Bianchi 恒等式、物质测地线解释兼容。

但守恒不等于变分性。

很多非变分方程也可以构造为守恒形式。Helmholtz 条件要求的是更强的二阶变分对称性：

\[
D F=(D F)^\dagger.
\]

所以不能因为 \(\mathcal C_{\mu\nu}\) 守恒，就直接认定它来自某个局域作用量。

## 8. 可能的 action 化路线

### 路线 A：保持 equation-first

最保守的路线是承认当前理论暂时不是 action-first，而是一个协变 equation-first 理论：

\[
\Pi_E^\perp\mathcal R=0,
\qquad
\mathrm{Tr}_g(\Pi_E\mathcal R)=0,
\qquad
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
\]

这条路线的优点是方程清楚，缺点是理论规格较低，需要额外证明 Cauchy 适定性、约束传播和实验一致性。

### 路线 B：辅助应力场 action

把

\[
\mathcal C_{\mu\nu}=\lambda_I E^I_{\mu\nu}
\]

视为某种辅助介质/各向异性应力张量。

然后写一个辅助场作用量，使它的 stress tensor 等于 \(\mathcal C_{\mu\nu}\)，并由辅助场方程给出：

\[
\mathcal C^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0,
\qquad
Q\to0\Rightarrow\mathcal C\to0.
\]

这条路线更有希望 action 化，但代价是必须引入新的自由度或约束场。

### 路线 C：乘子约束 action

写形式作用量：

\[
S_{\rm con}
=
\int\sqrt{|\tilde g|}
\left[
\Lambda^{\mu\nu}_{\perp}
(\Pi_E^\perp\mathcal R)_{\mu\nu}
+\psi\,\tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}
\right].
\]

对 \(\Lambda_\perp,\psi\) 变分，可以强制投影方程和 trace 方程。

但对 \(\tilde g_{\mu\nu}\) 变分会产生新的高阶方程，且 \(\Pi_E\) 显含 \(u,r\)，还会影响物质变分。

所以这不是“免费得到原方程”的办法，而是一个新的约束系统。

## 9. 当前判断

现在可以把结论分成三层：

第一层，最小方程是否仍可作为显式理论候选？

\[
\boxed{\text{可以。}}
\]

它是清楚的、协变的、可数值检验的 equation-first 候选。

第二层，它是否已经是一个普通 action-first 理论？

\[
\boxed{\text{不能这样说。}}
\]

投影算子的场依赖和 7 条 metric 条件的计数压力，都说明它很可能不满足普通 Helmholtz 条件。

第三层，它是否完全不能 action 化？

\[
\boxed{\text{还不能否定。}}
\]

更合理的 action 化方向是引入辅助场/乘子/介质变量，让 trace 条件和 \(Q\to0\) 分支条件来自辅助场方程，而不是要求裸投影后的 metric equation 自己就是 Euler-Lagrange 方程。

## 10. 下一步

下一步应做两件事：

1. 理论上：把辅助应力场路线写成最小 action 原型，明确哪些变量是自由场，哪些方程来自 metric variation，哪些来自辅助场 variation。
2. 数值上：在 `n=384` 可信条纹分辨率上复查 trace=0 投影方程，判断 `n=96` 的 trace=0 优势是否仍成立。

如果 `n=384` 仍支持 trace=0，那么主线应转向：

\[
\boxed{
\text{trace=0 投影方程}
+
\text{辅助应力场 action 化}
+
\text{patch/branch 规则}
}
\]

而不是回到纯 \(f(R)\) 或纯 metric 作用量扫描。
