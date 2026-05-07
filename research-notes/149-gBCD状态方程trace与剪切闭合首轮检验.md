# 149. gBCD 状态方程：trace 与剪切闭合首轮检验

日期：2026-05-07

## 1. 目的

上一轮确认：在当前 `2+1 txz` 约化中，

\[
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0,
\qquad
\lambda_I=(A,B,C,D)
\]

主符号 rank 为 3，而未知量有 4 个。因此 full conservation 不能独立产生全部 \(A,B,C,D\)，必须补一条状态方程/规范条件。

本轮先检验最小代数闭合：

1. trace closure：

\[
\mathcal C^\mu{}_\mu=F(Q),\qquad F(0)=0.
\]

2. shear-free closure：

\[
D=0,
\]

即不允许 \(u-r\) 剪切项。

## 2. 实现

新增脚本：

`kg_examples/fit_gbcd_trace_closure.py`

该脚本在原有 gBCD/full-conservation 的 nullspace hard constraint 基础上，再加入一组代数硬约束：

- `zero`：
  \[
  \mathcal C^\mu{}_\mu=0.
  \]
- `q1`：
  \[
  \mathcal C^\mu{}_\mu=\alpha\,Q/Q_0.
  \]
- `q2`：
  \[
  \mathcal C^\mu{}_\mu=\alpha\,Q/Q_0+\beta\,(Q/Q_0)^2.
  \]
- `shear0`：
  \[
  D=0.
  \]

其中

\[
Q=X_{\rm flat}-m^2,
\qquad
X_{\rm flat}=S_t^2-S_x^2-S_z^2,
\]

采用当前约定 \(X-m^2-Q=0\)。

所有 hard constraints 均通过 nullspace solve 保持。

## 3. trace closure 结果

正式输出：

- `visualizations/equation_first_gbcd_trace_closure_n96_tau0_trusted/`
- `visualizations/equation_first_gbcd_trace_closure_n96_taum3p5_trusted_zero/`
- `visualizations/equation_first_gbcd_trace_closure_n96_taup3p5_trusted_zero/`

tau=0：

| 闭合 | 中心代数 weighted mean | 中心 p95 | 结论 |
|---|---:|---:|---|
| `zero` | `0.02638` | `0.07935` | 比无 trace 的 `0.01187` 变差，但仍可研究 |
| `q1` | `2.01557` | `5.86496` | 失败 |
| `q2` | `0.78269` | `2.20743` | 失败 |

`trace=0` 跨时间切片：

| 中心切片 | 中心代数 weighted mean | 中心 p95 |
|---|---:|---:|
| \(\tau=-3.5\) | `0.11846` | `0.41823` |
| \(\tau=0\) | `0.02638` | `0.07935` |
| \(\tau=+3.5\) | `0.11031` | `0.38276` |

判断：

\[
\boxed{
\text{trace closure 能形式闭合，但跨时间残差明显恶化；简单 }F(Q)\text{ 型 trace 不成立。}
}
\]

特别是 `q1/q2` 的失败说明：用户提出的 \(Q\to0\) 回到 Einstein 方程是必要边界条件，但不能简单升级为 \(\mathcal C^\mu{}_\mu=F(Q)\) 这种只依赖 \(Q\) 的状态方程。

## 4. shear-free \(D=0\) 结果

正式输出：

- `visualizations/equation_first_gbcd_closure_shear0_n96_taum3p5_trusted/`
- `visualizations/equation_first_gbcd_closure_shear0_n96_tau0_trusted/`
- `visualizations/equation_first_gbcd_closure_shear0_n96_taup3p5_trusted/`

结果：

| 中心切片 | 中心代数 weighted mean | 中心 p95 | reduced 主符号 p95 条件数上界 |
|---|---:|---:|---:|
| \(\tau=-3.5\) | `0.10032` | `0.33380` | `1.43e9` |
| \(\tau=0\) | `0.03105` | `0.09533` | `4.62e9` |
| \(\tau=+3.5\) | `0.13447` | `0.57104` | `6.15e9` |

判断：

\[
\boxed{
D=0\text{ 也能形式闭合，但条件数很差，且分离态残差约 }10\%\text{ 以上。}
}
\]

因此简单“无剪切”不是好的状态方程。

## 5. 结论

本轮否定了两类最简单代数闭合：

- 单纯 trace closure；
- 单纯 shear-free closure。

这并不否定 gBCD 张量壳，因为无这些额外代数闭合时，gBCD/full-conservation/nullspace 的中心残差仍在 `1%~5%` 范围。

真正被否定的是：

\[
\boxed{
\text{用一条过于简单的代数关系直接固定缺失自由度。}
}
\]

## 6. 下一步

下一步应转向更自然的“变分型/辅助场型闭合”：

1. 把 \(\lambda_I=(A,B,C,D)\) 当作辅助场；
2. 用一个正定或半正定的辅助泛函选代表元，例如

\[
\int \sqrt{|\tilde g|}
\left[
\frac12 M^{IJ}\lambda_I\lambda_J
+\frac12 K^{IJ\mu\nu}\nabla_\mu\lambda_I\nabla_\nu\lambda_J
\right],
\]

再配合场方程投影和 full conservation；
3. 这样第 4 条闭合不是硬塞一个简单代数式，而是由“最小能量/最小范数/最小梯度”的 Euler-Lagrange 方程产生。

这也解释了为什么之前的 `norm/time` 代表元数值上表现更好：它们更接近这种辅助场变分闭合，而不是简单 trace 或剪切代数闭合。

