# 145 从 `g` 表象 Einstein 方程变换得到 `gtilde` 方程的分解路线

## 1. 核心想法

用户提出的新路线是：

从 A 支 `g` 表象 Einstein 方程
\[
M_P^2 G_{\mu\nu}[g]=T^{A}_{\mu\nu}[g,\rho,S]
\]
出发，直接做当前 disformal 变换，把方程写到 `gtilde` 表象。

然后把变换后的方程拆成两类：

1. 只由 \(\tilde g_{\mu\nu}\) 及其曲率/协变导数构成的 pure-\(\tilde g\) 部分；
2. 无法单纯由 \(\tilde g\) 衍生出来、必须显含 \(u_\mu,r_\mu,\rho\) 或其导数的剩余部分。

这比先猜 action 更结构化，因为它直接给出“若要等价于 A 支，D/hybrid 支缺的项到底是什么”。

## 2. 联络差分解

在一个局域可逆 patch 上形式地写
\[
g_{\mu\nu}=\mathfrak g_{\mu\nu}[\tilde g,u,r,\rho,\ldots].
\]

这里的 \(\mathfrak g_{\mu\nu}[\cdots]\) 不是爱因斯坦张量，也不是已经找到的简单显式公式，而是当前 disformal 变换在指定分支上的反映射记号。后续避免用大写 \(G[\cdots]\) 表示这个对象，以免和 Einstein tensor \(G_{\mu\nu}\) 混淆。对当前混合变换而言，\(\tilde g\) 依赖 \(g\)、\(u,r\) 以及 \(Q=\Box_g\sqrt\rho/\sqrt\rho\)，所以反解 \(g\) 一般不是局域代数函数，而可能是分支依赖的微分反问题。

因此在数值诊断中不应先强行反解 \(\mathfrak g\)。更可靠的做法是：从 A 参考解或 A 参考快照中直接使用已知的 \(g\)，同时使用由正向变换得到的 \(\tilde g\)，然后在这对 \((g,\tilde g)\) 上计算 \(C,\mathcal H,\mathcal R\)。只有当我们想把方程完全改写为 native `gtilde` 理论时，才需要研究反分支 \(\mathfrak g[\tilde g,u,r,\rho,\ldots]\) 的存在性、唯一性和正则性。

用 \(\tilde\nabla\) 表示 \(\tilde g\) 的 Levi-Civita 导数，定义两个联络之差
\[
C^\alpha{}_{\mu\nu}
=
\Gamma^\alpha{}_{\mu\nu}[\mathfrak g]
-
\tilde\Gamma^\alpha{}_{\mu\nu}[\tilde g].
\]

它是一个张量，并且
\[
C^\alpha{}_{\mu\nu}
=
\frac12 \mathfrak g^{\alpha\beta}
\left(
\tilde\nabla_\mu \mathfrak g_{\nu\beta}
+\tilde\nabla_\nu \mathfrak g_{\mu\beta}
-\tilde\nabla_\beta \mathfrak g_{\mu\nu}
\right).
\]

于是 Ricci 张量满足
\[
R_{\mu\nu}[\mathfrak g]
=
\tilde R_{\mu\nu}
+\mathcal D_{\mu\nu}[C],
\]
其中
\[
\mathcal D_{\mu\nu}[C]
=
\tilde\nabla_\alpha C^\alpha{}_{\nu\mu}
-\tilde\nabla_\nu C^\alpha{}_{\alpha\mu}
+C^\alpha{}_{\alpha\lambda}C^\lambda{}_{\nu\mu}
-C^\alpha{}_{\nu\lambda}C^\lambda{}_{\alpha\mu}.
\]

因此 Einstein 张量可写为
\[
G_{\mu\nu}[\mathfrak g]
=
\tilde G_{\mu\nu}[\tilde g]
+\mathcal H_{\mu\nu}[\tilde g,u,r,\rho,\ldots],
\]
其中
\[
\mathcal H_{\mu\nu}
=
\mathcal D_{\mu\nu}
-\frac12\mathfrak g_{\mu\nu}\mathfrak g^{\alpha\beta}\mathcal D_{\alpha\beta}
-\frac12\mathfrak g_{\mu\nu}\mathfrak g^{\alpha\beta}\tilde R_{\alpha\beta}
+\frac12\tilde g_{\mu\nu}\tilde g^{\alpha\beta}\tilde R_{\alpha\beta}.
\]

这里 \(\tilde G_{\mu\nu}\) 是最小 pure-\(\tilde g\) 部分，\(\mathcal H_{\mu\nu}\) 是由变换引入的剩余几何张量。它一般含有 \(u,r,\rho\) 及其导数。

## 3. 变换后的 Einstein 方程

把上式代入 A 支方程：
\[
M_P^2\left(\tilde G_{\mu\nu}+\mathcal H_{\mu\nu}\right)
=
T^{A}_{\mu\nu}[\mathfrak g,\rho,S].
\]

于是得到精确的 `gtilde` 表象方程：
\[
M_P^2\tilde G_{\mu\nu}[\tilde g]
=
T^{A}_{\mu\nu}[\mathfrak g,\rho,S]
-M_P^2\mathcal H_{\mu\nu}[\tilde g,u,r,\rho,\ldots].
\]

如果我们希望右边写成当前 `gtilde` 物质作用量给出的源
\[
T^{(\tilde m)}_{\mu\nu}[\tilde g,\tilde\rho,S]
\]
加额外项，则定义
\[
\mathcal R_{\mu\nu}
:=
T^{A}_{\mu\nu}[\mathfrak g,\rho,S]
-T^{(\tilde m)}_{\mu\nu}[\tilde g,\tilde\rho,S]
-M_P^2\mathcal H_{\mu\nu}.
\]

于是
\[
M_P^2\tilde G_{\mu\nu}
=
T^{(\tilde m)}_{\mu\nu}
+\mathcal R_{\mu\nu}.
\]

\(\mathcal R_{\mu\nu}\) 就是“当前 \(\tilde g\) 物质 action + 普通 EH[\(\tilde g\)]”缺掉的有效源项。

## 4. 这个分解告诉我们什么

如果 \(\mathcal R_{\mu\nu}=0\)，则普通 EH[\(\tilde g\)] + 当前 `gtilde` 物质 action 已经等价于 A 支。

如果 \(\mathcal R_{\mu\nu}\neq0\)，但可写成某个 pure-\(\tilde g\) 局域作用量的 metric variation，例如 \(f(\tilde R)\)、\(f(\tilde R,\tilde R_{\mu\nu}\tilde R^{\mu\nu})\) 或更高曲率项，则仍可保留 pure-\(\tilde g\) 引力理论。

如果 \(\mathcal R_{\mu\nu}\) 的张量方向或单值性依赖 \(u,r,\rho\)，则 pure-\(\tilde g\) 作用量不够；此时 \(\mathcal R_{\mu\nu}\) 直接告诉我们需要加入怎样的 \(u,r\)-dependent gravitational sector。

自然候选就是
\[
\int\sqrt{|\tilde g|}
\left[
A^{\mu\nu}(\mathcal I)\tilde R_{\mu\nu}
-2V(\mathcal I)
\right],
\]
因为 \(\mathcal H_{\mu\nu}\) 中的主结构正来自 \(\mathfrak g_{\mu\nu}\neq\tilde g_{\mu\nu}\) 以及 \(C=\Gamma[\mathfrak g]-\tilde\Gamma\)，它们会产生各向异性 Ricci 耦合和 \(\tilde\nabla\tilde\nabla A\) 型项。

## 5. 和 action 的关系

这条路线首先是 field-equation 级别的分解。它不会犯“只拉回 EH 就自动等价”的错误，因为我们显式保留了
\[
T^{A}_{\mu\nu}-T^{(\tilde m)}_{\mu\nu}
\]
这一项。

若要构造严格等价的 hybrid action，一个形式上正确的定义是
\[
S_{\rm needed}
=
\left(S_{\rm EH}[\mathfrak g]+S^{A}_{m}[\mathfrak g,\rho,S]\right)
-S^{(\tilde m)}_{m}[\tilde g,\tilde\rho,S].
\]

然后
\[
S_{\rm needed}+S^{(\tilde m)}_m
=
S_A^{\rm pullback}.
\]

这当然等价，但 \(S_{\rm needed}\) 一般非常复杂，也会显含 \(u,r,\rho\)。实际研究目标是检查 \(S_{\rm needed}\) 是否可以被一个简单的低阶作用量近似，例如 \(A^{\mu\nu}\tilde R_{\mu\nu}-2V\)。

## 6. 重要 caveat

1. 这个分解的“pure-\(\tilde g\) 部分”不是唯一的，除非先指定允许的 pure 几何基底。最小选择是把 \(\tilde G_{\mu\nu}\) 放左边，其余都放进 \(\mathcal R_{\mu\nu}\)。
2. 若允许 \(f(\tilde R)\)、\(f(\tilde R,I_2)\) 等高曲率项，也可以把它们从 \(\mathcal R_{\mu\nu}\) 中继续吸收到左边。
3. 当前 1550nm 平直 QFT 快照本身不是精确 A 支 Einstein 解，而是 \(M_P\to\infty\) 或弱引力近似下的物质参考解。因此数值检验时要区分：
   - 理论恒等式：从完整 A 支 Einstein-KG 解出发；
   - 实验近似检验：从平直 QFT 快照出发，看 \(\mathcal R_{\mu\nu}\) 的大小和结构。

## 7. 下一步

下一步最直接的工作是实现一个诊断脚本。第一版不需要求解反问题 \(g=\mathfrak g[\tilde g,u,r,\rho]\)，而是直接使用 A 参考中已知的 \(g\)：

1. 输入 A 参考快照 \((\rho,S,g)\) 与由变换得到的 \(\tilde g\)；
2. 计算 \(C^\alpha{}_{\mu\nu}\)、\(\mathcal H_{\mu\nu}\)；
3. 计算
   \[
   \mathcal R_{\mu\nu}
   =
   T^{A}_{\mu\nu}
   -T^{(\tilde m)}_{\mu\nu}
   -M_P^2\mathcal H_{\mu\nu};
   \]
4. 检查 \(\mathcal R_{\mu\nu}\) 在三切片上能否由 pure-\(\tilde g\) 几何张量解释；
5. 若不能，拟合它在
   \[
   \tilde g_{\mu\nu},\quad
   u_\mu u_\nu,\quad
   r_\mu r_\nu,\quad
   u_{(\mu}r_{\nu)}
   \]
   以及 \(\tilde\nabla\tilde\nabla A^{\mu\nu}\) 诱导方向上的分解。

这将把“猜 action”问题变成“对已知剩余张量做结构分解”问题。

## 8. 1550nm 高斯干涉上的第一版 \(\mathcal R\) 结构诊断

新增脚本：

- `kg_examples/diagnose_mathcal_r_pure_geometry.py`

它在平直 A 参考高斯干涉三切片 \(\tau_{\rm old}=-3.5,0,3.5\) 上计算：
\[
\mathcal R_{\rm need}/M_P^2
=
\tilde G_{\mu\nu}
-
\tilde T_{\mu\nu}/M_P^2,
\]
其中 \(\tilde T_{\mu\nu}\) 使用已修正的密度变换
\[
\sqrt{|\tilde g|}\tilde\rho=|X|\rho_A/m^2.
\]

第一版高分辨率诊断输出：

- `visualizations/mathcal_r_pure_geometry_n320_3tau/summary.json`
- `visualizations/mathcal_r_pure_geometry_n320_3tau/mathcal_r_pure_geometry_maps.png`

结果分两层理解：

1. 对完整 \(\mathcal R_{\rm need}/M_P^2\)，pure-\(\tilde g\) 几何基底当然能近乎精确拟合。即使不用显式 \(\tilde G_{\mu\nu}\)，只用 \(\tilde R_{\mu\nu}-\frac12\tilde R\tilde g_{\mu\nu}\) 也会得到 weighted residual \(\sim 3.7\times10^{-12}\)。这不是新的物理成功，而是因为目标本身几乎就是 \(\tilde G_{\mu\nu}\)。
2. 真正有判别力的是非平凡物质源 \(\tilde T_{\mu\nu}/M_P^2\)。在同一三切片、core10 区域上，用固定系数 pure geometry 基底拟合它的 weighted residual \(\sim0.996\)，逐点 p50 残差约 \(0.99\)。加入简单 matter directions \(\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\) 后 weighted residual 降到 \(\sim0.594\)，但仍不是精确闭合。

因此当前不能把“\(\mathcal R_{\rm need}\) 可被 pure geometry 拟合”解读为 pure-\(\tilde g\) 引力 action 路线被救活。它只是说明在普朗克质量下 \(\tilde T/M_P^2\sim10^{-60}\) 被 \(\tilde G\) 的数值尺度完全压住；若目标是构造与平直量子场论一致的严格动力学，非平凡剩余仍需要更细的 \(u,r,\rho\)-dependent sector 或更严格的 action-level 变分构造。

## 9. \(\mathcal R\) 是否是 \(\tilde g\) 的函数/泛函

用户指出：如果 \(\mathcal R_{\mu\nu}\) 是 pure geometry，就说明 pure-\(\tilde g\) 路线仍有希望。这一点成立，但需要区分两种层级：

1. 点态或有限 jet 函数：
   \[
   \mathcal R_{\mu\nu}(x)=F_{\mu\nu}\bigl(\tilde g,\partial\tilde g,\ldots,\partial^N\tilde g\bigr)(x).
   \]
2. 作用量变分产生的几何张量：
   \[
   \mathcal R_{\mu\nu}
   =
   -\frac{2}{\sqrt{|\tilde g|}}
   \frac{\delta S_{\rm extra}[\tilde g]}{\delta \tilde g^{\mu\nu}}.
   \]

第二个条件更强。即使某个张量可以由 \(\tilde g\) 的局部几何特征拟合出来，也还要满足 diffeomorphism covariance、Bianchi 型恒等式 \(\tilde\nabla^\mu\mathcal R_{\mu\nu}=0\)（或与左侧其他项组合后恒等守恒）以及 Helmholtz/self-adjoint integrability，才一定来自 pure metric action。

新增函数性诊断脚本：

- `kg_examples/diagnose_mathcal_r_gtilde_function.py`

诊断方法：留一时间切片 kNN 预测。用两个时间切片上的局部 \(\tilde g\) 几何特征训练，预测第三个切片上的非平凡源
\[
\tilde T_{\mu\nu}/M_P^2.
\]
不使用完整 \(\mathcal R_{\rm need}\) 作目标，因为完整目标会被平凡 \(\tilde G_{\mu\nu}\) 主项支配。

高分辨率输出：

- `visualizations/mathcal_r_gtilde_function_n320_3tau/summary.json`
- `visualizations/mathcal_r_gtilde_function_n320_3tau/mathcal_r_gtilde_function_knn.png`

结果：

- 只用曲率标量特征 \(\tilde R,\mathrm{tr}\tilde R^2,\mathrm{tr}\tilde R^3,\det\tilde g\)：overall p50 residual \(\sim0.747\)，weighted mean \(\sim0.924\)。
- 用更丰富的局部几何张量特征 \(\tilde g,\tilde R_{\mu\nu},\tilde G_{\mu\nu},\tilde R^2_{\mu\nu},I_2\tilde g\)：overall p50 residual \(\sim0.598\)，weighted mean \(\sim0.863\)。
- 加入 \(u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\) 作为对照并未稳定改善，overall p50 \(\sim0.635\)，weighted mean \(\sim0.997\)。

当前解读：

1. 完整 \(\mathcal R_{\rm need}\) 的 pure-\(\tilde g\) 性质确实给出一个“有希望”的方向，但这个希望主要来自 \(\tilde G_{\mu\nu}\) 主项。
2. 对非平凡细节，当前数据不支持“它是简单的局部 \(\tilde g\) 函数”。至少用低阶局部曲率/张量 jet 不能稳定预测。
3. 仍未排除更高阶、非局域、或 action-level 特殊组合的 pure-\(\tilde g\) 泛函；但如果要继续 pure-\(\tilde g\) 路线，下一步应从 Helmholtz/integrability 和 \(\tilde\nabla\)-守恒条件入手，而不是继续只做低阶曲率拟合。

## 10. 更高阶但仍局域的 pure-\(\tilde g\) 路线

用户提出先考虑“更高阶的方程（纯 \(\tilde g\) 衍生，局域）”。这里要和前面已经扫过的高次 \(f(\tilde R)\) 区分开来：

1. 高次 \(f(\tilde R)\) 只是单标量曲率的高阶多项式或解析函数；
2. 真正的“更高阶局域 pure-\(\tilde g\)”应当是有限 jet 的局域度规泛函，允许独立曲率不变量与导数不变量。

一个自然的局域 action 基底可以写成
\[
S_{\rm grav}^{\rm loc}[\tilde g]
=
\int\sqrt{|\tilde g|}
\left[
+\Lambda
+\frac{M_P^2}{2}\tilde R
+\alpha_1 \tilde R^2
+\alpha_2 \tilde R_{\mu\nu}\tilde R^{\mu\nu}
+\alpha_3 \tilde\nabla_\mu \tilde R\,\tilde\nabla^\mu \tilde R
+\alpha_4 \tilde R\,\tilde\Box\tilde R
+\alpha_5 \tilde\nabla_\lambda \tilde R_{\mu\nu}\tilde\nabla^\lambda \tilde R^{\mu\nu}
+\cdots
\right].
\]

对应的场方程会是四阶或更高阶的局域 PDE；其中常见的主结构包括
\[
(\tilde g_{\mu\nu}\tilde\Box-\tilde\nabla_\mu\tilde\nabla_\nu)\tilde R,\qquad
\tilde\Box \tilde R_{\mu\nu},\qquad
\tilde\nabla_\mu\tilde\nabla_\nu(\tilde R_{\alpha\beta}\tilde R^{\alpha\beta}),
\]
以及更高导数的张量项。

这条路的优点是仍然完全 pure-\(\tilde g\) 且局域；缺点是：

- principal part 会明显升阶，适定性和初值要求更苛刻；
- 不能只看代数项，必须同时检查最高阶导数项的符号、退化和混合主部；
- 在 4D 里要记住 Gauss-Bonnet/曲率恒等式会减少独立基底；在当前 2+1d 约化情形下独立不变量还会更少。

因此，若继续 pure-\(\tilde g\) 局域路线，下一步最合理的顺序不是继续加高 \(f(\tilde R)\) 阶数，而是：

1. 先做 quadratic curvature basis 的完整线性/非线性诊断；
2. 再加导数不变量；
3. 最后才考虑 action-level integrability 是否成立。

这比继续把单变量 \(f(\tilde R)\) 的次数抬高更有信息量。

## 11. 4 阶 quadratic local scan 的直接结论

为了回答“我们的例子下，大概要到多少阶才可能拟合得比较好”，补做了一个直接的 4 阶局域 action 扫描。基底取标准 quadratic curvature 局域项：

\[
S_{\rm grav}^{(4)}
=
\int\sqrt{|\tilde g|}
\left[
+\Lambda
+\frac{M_P^2}{2}\tilde R
+\alpha \tilde R^2
+\beta \tilde R_{\mu\nu}\tilde R^{\mu\nu}
\right].
\]

然后用它的 metric variation 张量基底去拟合 1550nm 高斯干涉三切片上的非平凡源 \(\tilde T_{\mu\nu}/M_P^2\)。

结果（`n=160`，`core10`，三切片 \(\tau=-3.5,0,3.5\)）：

- 单片加权残差：
  - \(\tau=-3.5\): `0.993246`
  - \(\tau=0\): `0.773260`
  - \(\tau=3.5\): `0.997984`
- 三切片总体加权残差：`0.998489`

这说明：

1. 2 阶局域 pure-\(\tilde g\) 基底不够；
2. 4 阶 quadratic local action 也还是不够，而且整体几乎没有把问题解决掉；
3. 若继续坚持 pure-\(\tilde g\) 且局域，至少应考虑 6 阶及以上的导数不变量，或者承认仅靠局域 pure-\(\tilde g\) 仍不足以闭合当前高斯干涉的非平凡剩余。

所以“较好拟合”的经验下限，当前可以保守地说是：

\[
\boxed{\text{至少超过 4 阶；更像是 6 阶或更高，甚至可能需要非局域或 }u,r\text{-dependent 结构。}}
\]

补充：随后又做了一个最轻量的 6 阶纯 \(f(\tilde R)\) sanity check（`n=320`，`core10`，三切片 \(\tau=-3.5,0,3.5\)），总体加权残差仍为 `0.998034`，单片中位残差仍在 `0.994\sim0.999` 附近。也就是说，把纯 \(f(\tilde R)\) 的多项式次数从 4 提到 6，几乎没有带来实质改善。

补充 2：按用户澄清后，又改成了真正的 6 阶 finite-jet 局域 action 原型，最小基底取
\[
\{1,\tilde R,I_2,\tilde R^2,\tilde R I_2,\tilde R^3\},
\qquad
I_2=\tilde R_{\mu\nu}\tilde R^{\mu\nu}.
\]
在 `n=160`、`core10`、三切片 \(\tau=-3.5,0,3.5\) 上，三切片总体加权残差仍为 `0.9986743707235821`，单片 core10 中位残差仍约 `0.995\sim0.999`。这说明即使进入真正的 6 阶局域 finite-jet 子空间，当前高斯干涉例子的非平凡源仍没有被显著解释掉。
