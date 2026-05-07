# 173-shadow-stress与物理分支传播条件

日期：2026-05-08

## 一句话结论

shadow-field 分支只解决了第一层问题：辅助 sector 可以不直接变分物质变量

\[
S,\qquad \tilde\rho,
\]

从而不直接破坏物质 Hamilton-Jacobi 方程和连续性方程。

但它没有自动解决第二层问题：辅助 action 自己对

\[
\tilde g_{\mu\nu},\qquad \Phi,\qquad \sigma
\]

的变分会产生额外的 shadow stress 和 shadow-source。若这些项不消失、不落在可接受的 \(E_{\rm sh}\) 子空间，或不保持物理分支

\[
\Phi=S+\mathrm{const},\qquad \sigma=\tilde\rho,
\]

那么 action 化仍然失败。

当前最重要的必要条件是：

\[
\boxed{
\left.
\frac{\delta S_{\rm state}}{\delta\Phi}
\right|_{\rm phys}=0,
\qquad
\left.
\frac{\delta S_{\rm state}}{\delta\sigma}
\right|_{\rm phys}=0,
}
\]

以及

\[
\boxed{
\Theta_{\mu\nu}/M_P^2
\text{ 在目标分支上仍满足投影、无迹、守恒和 }Q\to0\text{ 条件。}
}
\]

这里“phys”表示物理分支：

\[
U=u,\qquad R=r,\qquad \sigma=\tilde\rho.
\]

## 1. 需要检查的两个问题

上一轮的 shadow-field 方案引入：

\[
\Phi,\qquad \sigma,
\]

\[
U_\mu=\partial_\mu\Phi,
\qquad
R_\mu=\tilde\nabla_\mu\ln\sqrt{\sigma}.
\]

辅助张量空间为：

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

由于 \(S_{\rm aux}\) 不直接依赖 \(S,\tilde\rho\)，有：

\[
\frac{\delta S_{\rm aux}}{\delta S}=0,
\qquad
\frac{\delta S_{\rm aux}}{\delta\tilde\rho}=0.
\]

这保护了物质方程。

但还必须检查：

1. shadow stress 问题：
   \[
   \Theta_{\mu\nu}
   =
   -\frac{2}{\sqrt{|\tilde g|}}
   \frac{\delta S_{\rm aux}}{\delta\tilde g^{\mu\nu}}
   \]
   是否仍能作为我们想要的 \(\mathcal C_{\mu\nu}\)。
2. 分支传播问题：
   \[
   \Phi=S+\mathrm{const},\qquad \sigma=\tilde\rho
   \]
   是否在完整耦合方程中保持。

这两个问题都不是自动成立的。

## 2. shadow stress 的分解

设

\[
S_{\rm aux}
=
S_{\rm source}
+S_{\rm state}
+S_{\rm sh}.
\]

其中

\[
S_{\rm source}
=
-\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\tilde g^{\mu\nu}C_{\mu\nu}.
\]

则总有效应力为：

\[
\Theta_{\mu\nu}
=
\Theta_{\mu\nu}^{\rm source}
+\Theta_{\mu\nu}^{\rm state}
+\Theta_{\mu\nu}^{\rm sh}.
\]

若 \(C_{\mu\nu}\) 暂时当作独立 covariant tensor，则

\[
\Theta_{\mu\nu}^{\rm source}
=
M_P^2
\left(
C_{\mu\nu}
-\frac12\tilde g_{\mu\nu}C^\alpha{}_\alpha
\right).
\]

在

\[
C^\alpha{}_\alpha=0
\]

时，

\[
\Theta_{\mu\nu}^{\rm source}=M_P^2 C_{\mu\nu}.
\]

但完整 metric 源是：

\[
\boxed{
\Theta_{\mu\nu}/M_P^2
=
C_{\mu\nu}
+\frac{1}{M_P^2}
\left(
\Theta_{\mu\nu}^{\rm state}
+\Theta_{\mu\nu}^{\rm sh}
\right).
}
\]

因此要回到 trace=0 投影方程，必须要求：

\[
\boxed{
\Theta_{\mu\nu}^{\rm state}
+\Theta_{\mu\nu}^{\rm sh}
\text{ 要么在目标分支上为 0，}
\text{ 要么也落入允许的无迹守恒 }E_{\rm sh}\text{ 修正内。}
}
\]

这个条件很强。

## 3. 线性乘子约束的问题

上一轮形式上写过：

\[
S_{\rm state}
=
\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\left[
\alpha C^\mu{}_\mu
+\zeta^\nu\tilde\nabla^\mu C_{\mu\nu}
+\Lambda_\perp^{\mu\nu}
(\Pi_{E_{\rm sh}}^\perp C)_{\mu\nu}
+\mu(Q_{\rm sh})C_{\mu\nu}C^{\mu\nu}
\right].
\]

这个写法的优点是约束清楚。

但它有两个严重副作用。

第一，线性乘子项即使在约束成立时，metric variation 也不一定消失。例如：

\[
\alpha C^\mu{}_\mu
=
\alpha \tilde g^{\mu\nu}C_{\mu\nu}.
\]

即使

\[
C^\mu{}_\mu=0,
\]

对 metric 变分仍有：

\[
\delta_g(\alpha \tilde g^{\mu\nu}C_{\mu\nu})
\supset
\alpha C_{\mu\nu}\delta \tilde g^{\mu\nu}.
\]

所以 \(\alpha\) 会贡献 stress。

第二，投影约束项会给 shadow fields 加源：

\[
\delta_{\Phi,\sigma}
\left[
\Lambda_\perp^{\mu\nu}
(\Pi_{E_{\rm sh}}^\perp C)_{\mu\nu}
\right]
=
\Lambda_\perp^{\mu\nu}
(\delta_{\Phi,\sigma}\Pi_{E_{\rm sh}}^\perp C)_{\mu\nu}.
\]

这项即使

\[
\Pi_{E_{\rm sh}}^\perp C=0
\]

也不必为零，因为变分打在投影算子上。

所以：

\[
\boxed{
\text{线性乘子能强制约束，但通常会产生额外 stress 和 shadow-source。}
}
\]

这会威胁物理分支传播。

## 4. \(Q_{\rm sh}\)-门控项的问题

门控项

\[
\mu(Q_{\rm sh})C_{\mu\nu}C^{\mu\nu}
\]

看起来很自然，因为它可以在

\[
Q_{\rm sh}\to0
\]

时把 \(C_{\mu\nu}\) 压向 0。

但如果它真的在 action 中，并且

\[
Q_{\rm sh}
=
\frac{\tilde\square\sqrt{\sigma}}{\sqrt{\sigma}},
\]

那么对 \(\sigma\) 的变分会产生：

\[
\delta_\sigma
\left[
\mu(Q_{\rm sh})C^2
\right]
=
\mu'(Q_{\rm sh})C^2\,\delta_\sigma Q_{\rm sh}
+\cdots.
\]

除非：

\[
C^2=0,
\]

或

\[
\mu'(Q_{\rm sh})=0,
\]

否则 shadow density 方程会被改写。

因此：

\[
\boxed{
\text{把 }Q_{\rm sh}\text{ 门控作为 action 权重，会破坏 }
\sigma=\tilde\rho
\text{ 分支传播。}
}
\]

更安全的做法是：

\[
\boxed{
Q\to0\Rightarrow C\to0
\text{ 先作为分支/边界/正则性条件，而不是 variational 权重。}
}
\]

如果以后一定要 action 化这个门控，就必须加入补偿项，使 \(\delta S_{\rm state}/\delta\sigma\) 在物理分支上消失。

## 5. 平方型状态项的优点和缺点

为了减少乘子带来的 shadow-source，可以考虑平方型状态项：

\[
S_{\rm state}^{(2)}
=
\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\left[
\kappa_\perp
\|\Pi_{E_{\rm sh}}^\perp C\|^2
+\kappa_{\rm tr}(C^\mu{}_\mu)^2
+\kappa_{\rm div}
\|\tilde\nabla^\mu C_{\mu\nu}\|^2
\right].
\]

它的优点是：在目标约束面

\[
\Pi_{E_{\rm sh}}^\perp C=0,
\qquad
C^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0
\]

上，\(S_{\rm state}^{(2)}\) 及其一阶变分通常为 0。

所以它比较容易满足：

\[
\left.
\frac{\delta S_{\rm state}^{(2)}}{\delta\Phi}
\right|_{\rm target}=0,
\qquad
\left.
\frac{\delta S_{\rm state}^{(2)}}{\delta\sigma}
\right|_{\rm target}=0.
\]

这有利于分支传播。

但缺点也很明显：

\[
\delta_C S_{\rm state}^{(2)}
\]

在目标约束面也趋于 0，因此它不能平衡

\[
\delta_C S_{\rm source}
\]

里的线性源。

也就是说，平方项适合作为“on-shell 不扰动”的诊断/正则结构，但未必能单独决定 \(C_{\mu\nu}\)。

因此需要额外的 constitutive potential：

\[
V(C,\chi,\ldots),
\]

使得

\[
\delta_C
\left[
S_{\rm source}+S_V+S_{\rm state}^{(2)}
\right]
=0
\]

在目标 \(C_{\mu\nu}\) 上成立。

这说明 action 化并没有被简单解决。

## 6. 物理分支传播条件

定义差变量：

\[
\delta\Phi=\Phi-S,
\]

\[
\delta\sigma=\sigma-\tilde\rho.
\]

若 shadow 方程和物质方程同型，并且没有来自 \(S_{\rm state}\) 的额外源，则：

\[
\delta\Phi|_\Sigma=0,
\qquad
\partial_n\delta\Phi|_\Sigma=0,
\qquad
\delta\sigma|_\Sigma=0
\]

会在适定唯一演化下保持为 0。

这个条件可以概括为：

\[
\boxed{
\left.
\mathcal E_\Phi^{\rm state}
\right|_{\rm phys}=0,
\qquad
\left.
\mathcal E_\sigma^{\rm state}
\right|_{\rm phys}=0.
}
\]

其中：

\[
\mathcal E_\Phi^{\rm state}
=
\frac{\delta S_{\rm state}}{\delta\Phi},
\qquad
\mathcal E_\sigma^{\rm state}
=
\frac{\delta S_{\rm state}}{\delta\sigma}.
\]

若这两个量不为 0，则差变量满足非齐次方程：

\[
\mathcal L_\Phi(\delta\Phi,\delta\sigma)
=
-\mathcal E_\Phi^{\rm state}|_{\rm phys},
\]

\[
\mathcal L_\sigma(\delta\Phi,\delta\sigma)
=
-\mathcal E_\sigma^{\rm state}|_{\rm phys}.
\]

即使初始相等，后续也会被源项推离：

\[
U\neq u,
\qquad
R\neq r.
\]

所以分支传播的必要条件非常明确：

\[
\boxed{
S_{\rm state}
\text{ 必须在物理分支上对 shadow fields “一阶静默”。}
}
\]

这里“一阶静默”定义为：不仅 \(S_{\rm state}\) 的约束值为零，而且它对 \(\Phi,\sigma\) 的一阶变分也为零。

## 7. 当前最小健康方案

综合上面分析，目前最健康的 action 化方向不是线性硬约束，而是：

1. 保留 \(S_{\rm source}\) 提供主要 stress；
2. 用平方型状态项或等价结构，使投影/无迹/守恒约束在目标分支上一阶静默；
3. 用一个尚待构造的 constitutive potential \(V(C,\chi,\ldots)\) 来平衡 \(C\) 的变分；
4. 暂时把 \(Q\to0\Rightarrow C\to0\) 当作分支/边界/正则性条件，而不是直接放入 \(\mu(Q_{\rm sh})C^2\) 变分权重；
5. 要求 shadow-field 的物理分支满足：
   \[
   \mathcal E_\Phi^{\rm state}|_{\rm phys}=0,
   \qquad
   \mathcal E_\sigma^{\rm state}|_{\rm phys}=0.
   \]

这可以压缩成：

\[
\boxed{
\text{shadow action 必须是 stealth-state action：}
\text{它能选择 }C\text{ 的状态，但在物理分支上不推动 }U,R\text{ 偏离 }u,r。
}
\]

这里的 stealth-state action 是本笔记定义的术语，意思是“状态项在目标约束面上对 shadow fields 的一阶变分为零”。

## 8. 本轮结论

shadow-field 分支不是已经完成的 action 理论。

它目前通过了一个重要检查：

\[
\boxed{
\text{它可以避免辅助 action 直接改写物质变分。}
}
\]

但它还必须通过两个新检查：

\[
\boxed{
\Theta_{\mu\nu}^{\rm state}
+\Theta_{\mu\nu}^{\rm sh}
\text{ 是否可接受？}
}
\]

\[
\boxed{
\mathcal E_\Phi^{\rm state}|_{\rm phys}
=
\mathcal E_\sigma^{\rm state}|_{\rm phys}
=0
\text{ 是否成立？}
}
\]

如果这两个条件不成立，物理分支会漂移，或者 metric 方程右边不再是我们想要的 trace=0 投影修正。

因此当前判断是：

\[
\boxed{
\text{shadow-field 分支仍可行，但需要 stealth-state 构造；}
\text{简单线性乘子或 }Q_{\rm sh}\text{ 权重会带来问题。}
}
\]

## 9. 下一步

下一步理论任务：

1. 构造一个具体的 stealth-state action 候选；
2. 检查它的 \(\delta_C\) 方程能否给出非零、无迹、守恒、位于 \(E_{\rm sh}\) 内的 \(C_{\mu\nu}\)；
3. 若第 2 步失败，说明 action 化要么需要额外自由度 \(\chi\)，要么应退回 equation-first。

下一步数值任务仍然是：

1. 在 `n=384` 上复查 trace=0 投影方程；
2. 若 trace=0 优势保持，再对 stealth-state 条件做 fixed-background residual 检验。
