# 146. equation-first gBCD 闭合方程与守恒检验

日期：2026-05-07

## 1. 当前候选方程

当前从方程出发的最小可继续研究候选已经从 BCD 升级为 gBCD：

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

这里

- \(u_\mu=\partial_\mu S\)；
- \(r_\mu=\partial_\mu\ln \sqrt{\rho}\) 或等价的振幅梯度方向，具体数值脚本中来自 Bohm snapshot；
- \(A,B,C,D\) 是待闭合的标量系数字段；
- \(D u_{(\mu}r_{\nu)}\) 采用对称化约定，即 \(D(u_\mu r_\nu+r_\mu u_\nu)/2\)。

这个方程不是从某个已知局域作用量直接推出的，而是 equation-first 候选。作用量重建应后验检查。

## 2. 测地线条件与 full 守恒条件

Bianchi 恒等式给出

\[
\tilde\nabla^\mu
\left(
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu}
\right)=0.
\]

若只要求粒子轨迹仍是 \(\tilde g\) 测地线，则必要的方程级条件是横向力为零：

\[
h^\nu{}_\alpha \tilde\nabla^\mu \mathcal C_\mu{}^\alpha=0,
\]

其中 \(h^\nu{}_\alpha\) 是垂直于 \(u^\mu\) 的投影。

更强、更干净的条件是单独守恒：

\[
\tilde\nabla^\mu \mathcal C_{\mu\nu}=0.
\]

这会推出 \(\tilde\nabla^\mu\tilde T_{\mu\nu}=0\)，从而在 dust 型物质张量和连续性方程成立时给出测地线运动。数值上这比只约束横向力更强，但物理解释更清楚。

## 3. gBCD 守恒方程的显式形式

令所有升指标都由 \(\tilde g^{\mu\nu}\) 完成：

\[
u^\mu=\tilde g^{\mu\alpha}u_\alpha,\qquad
r^\mu=\tilde g^{\mu\alpha}r_\alpha.
\]

full 守恒方程展开为

\[
0=\tilde\nabla^\mu\mathcal C_{\mu\nu}.
\]

即

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

这个式子是 gBCD 系数闭合的核心 PDE。它给出 3+1d 中 4 个协变分量方程，或当前 2+1d 数值截面中的 3 个分量方程。

## 4. 为什么还没有完整动力学

gBCD 的关键点是：\(\mathcal C_{\mu\nu}\) 形式上像一个由 \(u,r\) 方向定义的各向异性有效应力张量。

但 \(A,B,C,D\) 四个标量不是自动由守恒方程唯一决定的：

- 在 2+1d 中，full 守恒给 3 个一阶 PDE；
- 在 3+1d 中，full 守恒给 4 个一阶 PDE；
- 方程本身还包括 Einstein 型的几何演化方程；
- 但若 \(A,B,C,D\) 被视为独立物理场，就仍需要类似“状态方程/本构关系”的额外规则。

因此当前结果不能说“完整 D 支动力学已经建立”。更准确的说法是：

\[
\boxed{
\text{gBCD 给出了一个能满足守恒与测地线要求的最小张量壳，}
\text{但还需要确定 } A,B,C,D \text{ 的本构闭合。}
}
\]

## 5. 数值检验结果

脚本：

`kg_examples/fit_equation_first_constrained_bcd.py`

新增能力：

- `--atoms g,uu,rr,ur`：使用 gBCD；
- `--hard-constraint`：用 KKT 系统强制约束；
- `--force-mode transverse`：硬约束 \(h\nabla C=0\)；
- `--force-mode full`：硬约束 \(\nabla C=0\)。

### 5.1 transverse hard constraint

输出：

`visualizations/equation_first_constrained_gbcd_hard_n96_3tau_summary/`

三时刻中心切片 weighted mean 代数残差：

- \(\tau=-3.5\)：`0.02858`
- \(\tau=0\)：`0.01177`
- \(\tau=+3.5\)：`0.04607`

横向力自然尺度 \(|J_\perp|/(|C|/h)\) 被压到 \(10^{-15}\sim10^{-11}\) 量级。

### 5.2 full divergence-free hard constraint

输出：

`visualizations/equation_first_constrained_gbcd_hard_full_n96_3tau_summary/`

三时刻中心切片 weighted mean 代数残差：

- \(\tau=-3.5\)：`0.02859`
- \(\tau=0\)：`0.01179`
- \(\tau=+3.5\)：`0.04641`

full 守恒约束也可被压到接近数值零。相邻 probe 时间层的代数残差升高到约 `0.14~0.23`，说明 full 守恒比 transverse 条件更强，会更明显地暴露系数时间闭合问题。

## 6. 直接推论

1. BCD 不是严格闭合的最小方程。

硬约束 \(h\nabla C=0\) 时，仅 \(Buu+Crr+D\,ur\) 会让中心代数残差暴涨，已经被当前 trusted 区域诊断排除为严格最小候选。

2. gBCD 是当前最小可继续候选。

加入 \(A\tilde g_{\mu\nu}\) 后，可以同时保持低代数残差和守恒/测地线约束。

3. full 守恒路线比 transverse 路线更适合后续理论化。

因为它自动使 \(\tilde T_{\mu\nu}\) 单独守恒，物质部分更接近标准最小耦合图景。但它需要更清楚的 \(A,B,C,D\) 本构闭合。

4. 逐切片投影不是动力学。

当前 KKT 解只说明：给定 A-reference 生成的 \(\tilde g,u,r\)，存在一组 \(A,B,C,D\) 可使方程近似成立并满足守恒。它尚未告诉我们真实 D 支演化时 \(A,B,C,D\) 如何从初值唯一推进。

## 7. 下一步

下一步应优先做两个检查：

1. 本构关系检查：

判断硬约束得到的 \(A,B,C,D\) 是否能写成局部标量函数，例如

\[
A=A(\rho,\tilde X,\tilde R,I_2,u\cdot r,r^2,\ldots),
\]

以及 \(B,C,D\) 的类似形式。

若可以，gBCD 可成为无新自由度的闭合方程。

2. 若不能，则把 \(A,B,C,D\) 视为辅助场：

需要给出额外状态方程或演化方程，并检查它是否可由作用量或 Helmholtz 条件支持。

在进入全动力学模拟前，必须先完成这一步，否则模拟器仍然会依赖人为投影规则。
