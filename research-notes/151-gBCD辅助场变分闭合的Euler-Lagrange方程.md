# 151. gBCD 辅助场变分闭合的 Euler-Lagrange 方程

日期：2026-05-07

## 1. 目标

当前 equation-first 候选为

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2
+\lambda_I E^I_{\mu\nu},
\]

其中

\[
\lambda_I=(A,B,C,D),
\]

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

上一轮已经确认：

1. full conservation

\[
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0
\]

只有 3 条一阶方程，不能单独决定 4 个 \(\lambda_I\)；

2. 简单代数状态方程，例如 trace=\(F(Q)\) 或 \(D=0\)，会显著恶化残差；

3. \(Q\to0\) 回 Einstein 方程应作为极限一致性条件，而不是强行写成 trace=\(F(Q)\)。

因此下一步是把 \(\lambda_I\) 当作辅助场，用变分型闭合产生缺失的第 4 条方程。

## 2. 固定参考切片上的变分闭合泛函

在固定 \((\tilde g,u,r,\rho,S)\) 背景上，定义目标张量

\[
\mathcal R_{\mu\nu}
=
\tilde G_{\mu\nu}
-\tilde T_{\mu\nu}/M_P^2.
\]

希望

\[
\lambda_I E^I_{\mu\nu}\approx \mathcal R_{\mu\nu}.
\]

一个自然的投影闭合泛函是

\[
\mathcal J[\lambda,\xi]
=
\int d^dx\sqrt{|\tilde g|}
\left[
\frac12
\left(
\lambda_I E^I_{\mu\nu}-\mathcal R_{\mu\nu}
\right)
W^{\mu\nu\rho\sigma}
\left(
\lambda_J E^J_{\rho\sigma}-\mathcal R_{\rho\sigma}
\right)
\right.
\]

\[
\left.
+
\frac12 \mu(Q)M^{IJ}\lambda_I\lambda_J
+
\frac12 K^{IJ\alpha\beta}
\tilde\nabla_\alpha\lambda_I
\tilde\nabla_\beta\lambda_J
+
\xi^\nu\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})
\right].
\]

这里：

- \(W^{\mu\nu\rho\sigma}\)：张量残差的权重；
- \(M^{IJ}\)：辅助场质量矩阵；
- \(K^{IJ\alpha\beta}\)：辅助场梯度刚度；
- \(\xi^\nu\)：强制 full conservation 的 Lagrange multiplier；
- \(\mu(Q)\)：\(Q\to0\) 门控质量权重，例如

\[
\mu(Q)
=
\frac{1}{(\epsilon+|Q|/Q_0)^p}.
\]

当 \(Q\to0\) 时，\(\mu(Q)\) 变大，从而倾向 \(\lambda_I\to0\)，使方程回到 Einstein 形式。

## 3. 记号压缩

定义

\[
H_{IJ}
=
E^I_{\mu\nu}W^{\mu\nu\rho\sigma}E^J_{\rho\sigma},
\]

\[
b_I
=
E^I_{\mu\nu}W^{\mu\nu\rho\sigma}\mathcal R_{\rho\sigma}.
\]

则残差项对 \(\lambda_I\) 的变分给出

\[
H_{IJ}\lambda_J-b_I.
\]

## 4. Euler-Lagrange 方程

对 \(\lambda_I\) 变分，忽略边界项，得到

\[
\boxed{
H_{IJ}\lambda_J
-b_I
+\mu(Q)M_{IJ}\lambda_J
-\tilde\nabla_\alpha
\left(
K_{IJ}^{\alpha\beta}\tilde\nabla_\beta\lambda_J
\right)
-E^I_{\mu\nu}\tilde\nabla^\mu\xi^\nu
=0.
}
\]

对 \(\xi^\nu\) 变分，得到约束

\[
\boxed{
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0.
}
\]

所以固定背景上的辅助场闭合系统为

\[
\begin{cases}
H_{IJ}\lambda_J-b_I
+\mu M_{IJ}\lambda_J
-\tilde\nabla_\alpha(K_{IJ}^{\alpha\beta}\tilde\nabla_\beta\lambda_J)
-E^I_{\mu\nu}\tilde\nabla^\mu\xi^\nu
=0,\\
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0.
\end{cases}
\]

这是一组 saddle-point 方程：\(\lambda_I\) 是辅助应力系数，\(\xi^\nu\) 是 enforce conservation 的乘子。

## 5. 与当前数值实现的对应

当前脚本

`kg_examples/fit_gbcd_auxiliary_gauge.py`

已经实现了上述泛函的离散版本：

| 连续项 | 当前数值对应 |
|---|---|
| \(H_{IJ}\lambda_J-b_I\) | 代数残差最小化 |
| \(\tilde\nabla^\mu(\lambda_IE^I_{\mu\nu})=0\) | nullspace hard constraint |
| \(\mu(Q)M^{IJ}\lambda_I\lambda_J\) | `--q-gated-norm` |
| \(M^{IJ}\lambda_I\lambda_J\) | `--norm-weight` |
| 时间方向 \(K\nabla_t\lambda\nabla_t\lambda\) | `--time-weight` |
| 空间方向 \(K\nabla_i\lambda\nabla_i\lambda\) | `--space-weight` |

本轮数值结果说明：

- 纯代数 trace/Q 闭合不好；
- \(Q\)-门控范数项表现稳定，只小幅增加残差；
- 它同时降低系数粗糙度，符合辅助场 action 的直觉。

## 6. 作为真正动力学理论时的区别

上面的 \(\mathcal J[\lambda,\xi]\) 是“固定参考切片上的投影闭合”，用于从 A 参考生成的 \(\tilde g\) 上诊断可行性。

若要成为真正 D 支理论，应把 \(\tilde g\) 也作为自由场，写成某种总方程：

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2
+\lambda_I E^I_{\mu\nu},
\]

\[
\mathcal E_I[\lambda,\tilde g,u,r,Q]=0.
\]

其中 \(\mathcal E_I=0\) 应来自辅助场 action 对 \(\lambda_I\) 的变分。

如果完整 action 是协变的，且所有依赖的物质变量都按一致规则变分，那么 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\) 应该由 Bianchi identity 和物质方程自动协调，而不必人工每步投影。

但在当前 equation-first 阶段，保留 nullspace hard constraint 是安全做法，因为我们还没有完整 action。

## 7. 下一步

下一步应做两件事：

1. 数值层面：扫描 \(\mu(Q)\)、`time_weight`、`space_weight`，寻找残差与粗糙度的 Pareto 前沿。

2. 理论层面：尝试给出最小连续闭合方程

\[
\mathcal E_I
=
\mu(Q)M_{IJ}\lambda_J
-\tilde\nabla_\alpha(K_{IJ}^{\alpha\beta}\tilde\nabla_\beta\lambda_J)
+\text{可能的低阶耦合项}
=0,
\]

并检查它是否会把 \(\lambda_I\) 直接压成 0。若会，则必须保留场方程投影项或引入非齐次源项 \(b_I\)。

当前可行的最小模型是：

\[
\boxed{
\text{gBCD 张量壳}
+
\text{full conservation}
+
Q\text{-门控辅助场变分闭合}.
}
\]

补充：\(\mu(Q)\) 的首轮数值扫描已记录在
`research-notes/150-Q趋零条件作为辅助场变分权重.md` 第 7 节。
当前 tau=0 的温和折中点是 \(\epsilon=0.1,p=2\)，中心残差约 `1.28%`。
