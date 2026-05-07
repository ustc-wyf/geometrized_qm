# 144 放宽 pure `gtilde` 限制后的 `u,r` 引力作用量理论分析

## 1. 新目标

本轮路线调整是：

1. 仍取强约束，目标仍是让新表象的演化逐点接近、甚至等价于平直时空量子场论；
2. 暂时放开“引力作用量必须完全由 `gtilde` 衍生”的限制；
3. 允许 `gtilde` 表象下的引力作用量显含
   \[
   u_\mu=\partial_\mu S,\qquad r_\mu=\partial_\mu\sqrt{\rho}.
   \]

这一步把问题从“pure `gtilde` 几何 action 能否独立完成任务”改为“允许物质诱导方向进入 gravitational sector 后，能否构造强等价或强近似等价的 action”。

## 2. 必须修正的等价性口径

若当前 disformal 变换在某个局域 patch 上可逆，则存在反变换
\[
g_{\mu\nu}=G_{\mu\nu}[\tilde g,u,r,\rho,\ldots].
\]

形式上可以定义拉回的 EH 引力项
\[
S_{\rm grav}^{\rm pullback}[\tilde g,\rho,S]
=
\frac{M_P^2}{2}\int d^4x\,
\sqrt{-G[\tilde g,u,r,\rho,\ldots]}\,
R\!\left[G[\tilde g,u,r,\rho,\ldots]\right].
\]

但这里必须区分两件事。当前项目中的 `gtilde` 物质作用量不是把原 KG/Madelung 物质作用量逐项变量变换得到的；我们证明的是它变分后的物质方程与原物质方程等价。这是方程层面的等价，不是完整 action 的逐项拉回等价。

因此：

- 若把 A 支完整总作用量
  \[
  S_A[g,\rho,S]=S_{\rm EH}[g]+S_{\rm KG/Mad}[g,\rho,S]
  \]
  整体做可逆变量重写，则当然得到等价的总场方程；
- 但若只把 \(S_{\rm EH}[g]\) 拉回为 \(S_{\rm grav}^{\rm pullback}\)，再把物质项替换成另行构造的
  \[
  S_{\rm matter}^{(\tilde g)}[\tilde g,\tilde\rho,S],
  \]
  则这是一个 hybrid action，不是原总作用量的变量重写；
- 这个 hybrid action 不能自动推出与 A 支等价的 Einstein 方程，必须重新对 \(\tilde g,\rho,S,\tilde\rho\) 变分并检查。

密度映射例如
\[
\kappa\sqrt{|\tilde g|}\tilde\rho
=
X\sqrt{|g|}\rho
\]
或正测度分支版本，只能保证物质方程/守恒结构在相应口径下匹配；它本身不保证 metric variation 得到的源项与原 Einstein 方程一致。

修正后的结论是：

- “完整总作用量整体拉回”是严格等价但几乎同义反复的基准；
- “只拉回 EH + 使用新 `gtilde` 物质项”不是自动等价，而是需要检验的候选 action；
- 后续低阶 ansatz 应以 field-equation residual 和变分一致性为准，而不能只凭 \(S_{\rm EH}[g]\) 的拉回形式宣称等价。

## 3. 显含 `u,r` 后的物质变分

若
\[
S_{\rm grav}=S_{\rm grav}[\tilde g,u,r,\rho],
\]
则对 \(S\) 变分时有
\[
\delta u_\mu=\partial_\mu\delta S.
\]

因此连续性方程会出现额外项：
\[
\tilde\nabla_\mu J^\mu_{\rm matter}
+
\frac{1}{\sqrt{|\tilde g|}}
\tilde\nabla_\mu
\left(
\frac{\delta S_{\rm grav}}{\delta u_\mu}
\right)
=0
\]
的结构。

对 \(\rho\) 或 \(\sqrt{\rho}\) 变分时，也会从 \(r_\mu\) 和显含 \(\rho\) 的系数得到额外项。这些项不应自动看成错误；但也不能自动看成“保证等价”的补偿项。只有在完整总作用量整体拉回时，它们才必然与物质项共同重组为原 KG/Einstein 方程。对 hybrid action 来说，它们是必须显式检查的新贡献。

因此这条路线与“物质作用量单独对 \(u_\mu\) 变分给出纯 \(\tilde g\) 测地线”的目标不同。若仍要求物质 action 单独给测地线，则显含 \(u,r\) 的引力 action 泛型会破坏这一点；若目标是总作用量等价或近似等价于平直 QFT，则这是允许且必要的。

## 4. 首个低阶有效 ansatz

精确拉回作用量通常会很复杂。一个自然的低阶截断是各向异性 Ricci 耦合：
\[
S_{\rm grav}^{(u,r)}
=
\frac{M_P^2}{2}\int d^4x\sqrt{|\tilde g|}
\left[
A^{\mu\nu}(\mathcal I)\tilde R_{\mu\nu}
-2V(\mathcal I)
\right],
\]
其中
\[
A^{\mu\nu}
=
a_0(\mathcal I)\tilde g^{\mu\nu}
+a_1(\mathcal I)\hat u^\mu\hat u^\nu
+a_2(\mathcal I)\hat r^\mu\hat r^\nu
+a_3(\mathcal I)\hat u^{(\mu}\hat r^{\nu)}.
\]

\(\mathcal I\) 是无量纲标量集合，例如
\[
\tilde X=\tilde g^{\mu\nu}u_\mu u_\nu,\quad
\tilde Y=\tilde g^{\mu\nu}u_\mu r_\nu,\quad
\tilde Z=\tilde g^{\mu\nu}r_\mu r_\nu,\quad
\rho,
\]
以及必要时的 \(\tilde\nabla u,\tilde\nabla r\) 不变量。

这个 ansatz 的意义：

- \(a_0\tilde g^{\mu\nu}\tilde R_{\mu\nu}\) 是 EH-like 主项；
- \(u^\mu u^\nu\tilde R_{\mu\nu}\) 直接控制 Bohm 流方向上的曲率响应；
- \(r^\mu r^\nu\tilde R_{\mu\nu}\) 直接控制振幅梯度方向上的曲率响应；
- \(u^{(\mu}r^{\nu)}\tilde R_{\mu\nu}\) 控制相位流与振幅梯度混合方向；
- 这些正是 pure \(f(\tilde R)\) 和 pure \(f(\tilde R,I_2)\) 缺失的物质诱导张量方向。

## 5. metric 变分主结构

对
\[
\int\sqrt{|\tilde g|}\,A^{\mu\nu}\tilde R_{\mu\nu}
\]
做 metric 变分，若先把 \(A^{\mu\nu}\) 视作外给张量，则主结构为
\[
E_{\mu\nu}^{AR}
=
-\frac12\tilde g_{\mu\nu}A^{\alpha\beta}\tilde R_{\alpha\beta}
+\frac12
\left(
A_\mu{}^\alpha\tilde R_{\alpha\nu}
+A_\nu{}^\alpha\tilde R_{\alpha\mu}
\right)
\]
\[
+
\frac12
\left[
\tilde\Box A_{\mu\nu}
+\tilde g_{\mu\nu}\tilde\nabla_\alpha\tilde\nabla_\beta A^{\alpha\beta}
-\tilde\nabla_\alpha\tilde\nabla_\mu A^\alpha{}_\nu
-\tilde\nabla_\alpha\tilde\nabla_\nu A^\alpha{}_\mu
\right],
\]
再加上 \(A^{\mu\nu}\) 对 \(\tilde g_{\mu\nu}\) 的显式依赖变分项。

这点很重要：即使 action 只写成 \(A^{\mu\nu}\tilde R_{\mu\nu}\)，场方程中也自然出现 \(\tilde\nabla\tilde\nabla A\)。由于 \(A\) 由 \(u,r,\rho\) 组成，这些项会产生类似量子势的振幅/相位高阶导数结构。

## 6. 1+1d 与 3+1d 的差异

在 1+1d 或球对称径向约化中，若 \(u_\mu,r_\mu\) 线性无关，它们张成非平凡二维切空间。因此
\[
\tilde g_{\mu\nu},\quad
u_\mu u_\nu,\quad
r_\mu r_\nu,\quad
u_{(\mu}r_{\nu)}
\]
基本足以表达任意对称二张量。

但在一般 3+1d 中，\(u,r\) 只给出二维物质平面；横向二维子空间只通过 \(\tilde g_{\mu\nu}\) 被各向同性地看到。如果目标张量在横向平面有非平凡各向异性，仅靠代数 \(u,r\) basis 不够，需要加入：

- \(\tilde\nabla_\mu u_\nu\)、\(\tilde\nabla_\mu r_\nu\) 的对称/反对称部分；
- Hessian 型方向，例如 \(\tilde\nabla_\mu\tilde\nabla_\nu\rho\) 或 trace-free 部分；
- 更多曲率不变量或张量不变量。

## 7. 下一步

理论上建议的顺序：

1. 严格区分两条路线：完整总作用量整体拉回，以及“拉回/构造引力项 + 新 `gtilde` 物质项”的 hybrid action；
2. 先测试低阶
   \[
   A^{\mu\nu}(\mathcal I)\tilde R_{\mu\nu}-2V(\mathcal I)
   \]
   是否足够；
3. 在 1550nm 高斯干涉三切片上，用前面 pure \(f(R)\)、\(f(R,I_2)\) 相同的 fixed-background residual 口径拟合普适函数 \(a_i(\mathcal I),V(\mathcal I)\)；
4. 检验时必须包含 \(\tilde\nabla\tilde\nabla A^{\mu\nu}\) 导数项，不能只看代数项；
5. 若残差显著下降，再进入 Cauchy 主部、物质变分补偿项和 D 支投影初值求解；
6. 若仍失败，再加入 \(\nabla u,\nabla r\) 诱导的张量方向，而不是继续 pure \(f(\tilde R)\) 高阶多项式。

当前判断：

- 完整总作用量整体拉回可严格等价，但这不是当前 hybrid 物质 action 自动拥有的性质；
- 真正的非平凡问题是能否找到简单、普适、低阶、适定的 \(u,r\)-dependent effective gravitational action；
- 最自然的第一候选是 \(A^{\mu\nu}(u,r,\mathcal I)\tilde R_{\mu\nu}-2V(\mathcal I)\)。
