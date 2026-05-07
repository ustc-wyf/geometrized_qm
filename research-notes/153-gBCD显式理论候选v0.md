# 153. gBCD 显式理论候选 v0

日期：2026-05-07

## 1. 最终目标对齐

本项目最终目标不是在高斯干涉试验场上构造一个拟合器，而是得到一个可以写成理论的显式结构：

- 明确的场变量；
- 明确的场方程；
- 明确的闭合/选择原则；
- 在 \(Q\to0\) 或弱量子势极限下回到普通 Einstein 方程；
- 在高斯波包干涉试验场中不被数值残差否定。

本笔记把当前结果整理成第一个明确候选：

\[
\boxed{
\text{gBCD 显式理论候选 v0}
}
\]

它目前是 equation-first 理论，而不是已经完成的 action-first 理论。

## 2. 场变量

几何变量：

\[
\tilde g_{\mu\nu}.
\]

物质变量：

\[
\rho,\qquad S.
\]

定义：

\[
u_\mu=\partial_\mu S,
\qquad
r_\mu=\partial_\mu \ln\sqrt{\rho}
\]

或在当前数值实现中等价使用振幅梯度 covector。

辅助场变量：

\[
\lambda_I=(A,B,C,D).
\]

张量基底：

\[
E^I_{\mu\nu}
=
\left\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\right\}.
\]

修正张量：

\[
\mathcal C_{\mu\nu}
=
\lambda_I E^I_{\mu\nu}.
\]

## 3. 主场方程

候选主方程为

\[
\boxed{
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu}.
}
\]

即

\[
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C r_\mu r_\nu
+D u_{(\mu}r_{\nu)}.
\]

这里 \(\tilde T_{\mu\nu}\) 是变换后物质作用量对应的能动张量。物质部分仍按之前已完成的等价性分析来处理，使其给出质量壳方程、连续性方程和测地线解释。

## 4. 守恒条件

由 Bianchi identity 和物质守恒，应有

\[
\boxed{
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
}
\]

也就是

\[
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0.
\]

在 \(3+1d\) 中它有四个分量。需要强调：这四个分量是必要的，但由于当前 \(E^I_{\mu\nu}\) 的主符号结构通常落在 \(\mathrm{span}\{k,u,r\}\) 中，它不一定能独立、良态地决定全部 \(\lambda_I\)。因此还需要一个闭合/选择原则。

## 5. 闭合原则：约束解空间上的最小辅助泛函

v0 的核心闭合原则不是简单代数式，例如 trace=\(F(Q)\)，而是：

\[
\boxed{
\text{在满足主场方程和守恒条件的解空间中，使辅助场泛函最小。}
}
\]

当前最小版本为

\[
\boxed{
\mathcal A_\lambda
=
\frac12
\int d^4x\sqrt{|\tilde g|}
\left[
M^{IJ}\lambda_I\lambda_J
+K^{IJ}_{(t)}
\left(n^\alpha\tilde\nabla_\alpha\lambda_I\right)
\left(n^\beta\tilde\nabla_\beta\lambda_J\right)
\right].
}
\]

其中：

- \(M^{IJ}\) 是辅助场范数矩阵；
- \(K^{IJ}_{(t)}\) 是时间方向刚度矩阵；
- \(n^\mu\) 是当前 foliation 的单位时间方向。

当前数值默认对应：

\[
M^{IJ}\sim10^{-8}\delta^{IJ},
\qquad
K^{IJ}_{(t)}\sim10^{-5}\delta^{IJ},
\qquad
K_s=0.
\]

即：

\[
norm\_weight=10^{-8},
\qquad
time\_weight=10^{-5},
\qquad
space\_weight=0.
\]

## 6. 等价的受限变分写法

可以把 v0 写成带乘子的受限泛函：

\[
\mathfrak S_{\rm v0}
=
\mathcal A_\lambda
+
\int d^4x\sqrt{|\tilde g|}
\,
\Xi^{\mu\nu}
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
-\lambda_I E^I_{\mu\nu}
\right)
\]

\[
+
\int d^4x\sqrt{|\tilde g|}
\,
\xi^\nu
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu}).
\]

其中：

- \(\Xi^{\mu\nu}\) 强制主场方程；
- \(\xi^\nu\) 强制守恒条件。

对 \(\lambda_I\) 变分得到

\[
\boxed{
M^{IJ}\lambda_J
-\tilde\nabla_\alpha
\left(
K^{IJ\alpha\beta}\tilde\nabla_\beta\lambda_J
\right)
-\Xi^{\mu\nu}E^I_{\mu\nu}
-E^I_{\mu\nu}\tilde\nabla^\mu\xi^\nu
=0.
}
\]

这条方程就是缺失的闭合方程。它不是简单代数状态方程，而是由辅助场最小原则给出的 Euler-Lagrange 方程。

对 \(\Xi^{\mu\nu}\) 变分回到主场方程。

对 \(\xi^\nu\) 变分回到守恒条件。

这就是 v0 的显式方程系统。

## 7. 与当前数值试验的对应

当前数值器做的是 v0 的离散近似：

`kg_examples/fit_gbcd_auxiliary_gauge.py`

对应关系：

| 理论项 | 数值项 |
|---|---|
| 主场方程 | algebraic residual |
| 守恒条件 | nullspace hard constraint |
| \(M^{IJ}\lambda_I\lambda_J\) | `norm_weight` |
| \(K_t^{IJ}\nabla_t\lambda_I\nabla_t\lambda_J\) | `time_weight` |
| \(K_s^{IJ}\nabla_i\lambda_I\nabla_i\lambda_J\) | `space_weight` |

当前扫描给出默认候选：

\[
norm\_weight=10^{-8},
\quad
time\_weight=10^{-5},
\quad
space\_weight=0.
\]

三切片中心残差：

| 切片 | 中心残差 |
|---:|---:|
| \(\tau=-3.5\) | `0.03506` |
| \(\tau=0\) | `0.01241` |
| \(\tau=+3.5\) | `0.05017` |

因此 v0 在高斯干涉试验场中暂未被否定。

## 8. \(Q\to0\) 极限

本轮按用户要求暂时不把 \(Q\to0\) 条件放入默认候选。

但 v0 将来可自然加入该条件，只需把 \(M^{IJ}\) 替换为

\[
\mu(Q)M^{IJ},
\qquad
\mu(Q)\to\infty
\quad(Q\to0),
\]

即可让 \(\lambda_I\to0\)，从而

\[
\mathcal C_{\mu\nu}\to0,
\]

回到 Einstein 方程。

这说明 \(Q\to0\) 条件更适合作为辅助场范数权重或边界条件，而不是 trace=\(F(Q)\)。

## 9. 当前未完成点

v0 还不是最终理论，有三个未完成点：

1. 当前 \(n^\mu\) 引入了 foliation，理论上需要决定这是物理结构、规范选择，还是应改写成协变的 preferred direction，例如由 \(u^\mu\) 定义。

2. 当前空间刚度 \(K_s\) 的简单离散实现表现不好。若需要空间项，应改成协变 Laplacian 或弱形式。

3. 目前数值仍是三切片受限投影。下一步需要把 v0 的 Euler-Lagrange 方程改写成实际时间推进格式，而不是每个切片独立求最小。

## 10. 下一步

下一步目标：

\[
\boxed{
\text{从 v0 的受限变分方程推出可推进的 } \lambda_I \text{ 时间演化系统。}
}
\]

具体做法：

1. 固定默认参数 \(norm=10^{-8},time=10^{-5},space=0\)；
2. 从三层离散解中提取 \(\lambda_t,\lambda_{tt}\)；
3. 检查它是否满足 v0 的离散 Euler-Lagrange 关系；
4. 若成立，写出显式推进器；
5. 若不成立，说明仍缺少协变空间项、非齐次源项或更复杂的辅助场耦合。

## 11. 三层离散时间诊断

已对默认候选

\[
norm=10^{-8},\quad time=10^{-5},\quad space=0
\]

提取三层系数

\[
\lambda_-,\lambda_0,\lambda_+
\]

并计算归一化一阶差分与二阶差分：

\[
\Delta_1\lambda\sim \lambda_+-\lambda_0,\quad \lambda_0-\lambda_-,
\]

\[
\Delta_2\lambda=\lambda_+-2\lambda_0+\lambda_-.
\]

主要结果：

- \(\tau=0\) 干涉中心：所有分量时间差分都很小，例如 p95 二阶差分约 `0.002~0.007`；
- \(\tau=+3.5\)：二阶差分 p95 约 `0.05~0.12`；
- \(\tau=-3.5\)：`rr/ur` 方向最不平滑，二阶差分 p95 可到 `1.38` 与 `1.07`。

判断：

\[
\boxed{
\text{v0 若做显式时间推进，最容易出问题的是分离态的 }rr/ur\text{ 辅助自由度。}
}
\]

这提示下一步推进器不能只用简单显式更新。更合理的是：

- 对 \(\lambda_I\) 使用隐式或半隐式时间推进；
- 或者给 \(rr/ur\) 方向不同的 \(M^{IJ},K_t^{IJ}\) 权重；
- 或者加入协变空间弱形式，避免局部尖峰只靠时间项承担。
