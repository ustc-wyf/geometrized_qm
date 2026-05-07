# 148. gBCD 辅助应力场代表元规范与 nullspace 硬约束

日期：2026-05-07

## 1. 背景

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

为了让物质测地线解释保持干净，本轮优先采用强条件

\[
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
\]

这时 \(\mathcal C_{\mu\nu}\) 更像一个辅助各向异性应力张量。真正需要补充的是：在满足场方程投影和 full conservation 的解空间中，如何选择唯一的 \(A,B,C,D\) 代表元。

## 2. KKT 实现问题与修正

旧实现用 KKT 方程的最小二乘解来实现 hard constraint。这个方法在条件数较差的切片会出现“名义硬约束、数值泄漏”：约束行被放进 KKT 系统，但由于病态线性代数，\(\tilde\nabla^\mu\mathcal C_{\mu\nu}\) 的原始范数可能异常放大。

本轮新增脚本：

`kg_examples/fit_gbcd_auxiliary_gauge.py`

并把 hard constraint 求解器升级为 nullspace 投影法：

1. 先对约束矩阵 \(F\) 做 SVD；
2. 写 \(y=N z\)，其中 \(FN=0\)；
3. 只在零空间内最小化代数残差与代表元正则；
4. 因此 \(Fy=0\) 被机器精度满足。

这个修正不改变物理方程，只改变数值求解 hard constraint 的实现方式。

## 3. 代表元规范测试

测试的代表元规范：

- `norm`：最小归一化系数范数；
- `space`：最小归一化空间邻点差；
- `time`：最小归一化相邻时间层差；
- `space_time`：空间与时间平滑同时加入。

第一轮 tau=0 的 KKT 扫描显示：

- 空间平滑权重 `1e-6` 会把中心代数残差从约 `1.2%` 恶化到约 `38%`；
- 空间平滑权重降到 `1e-7` 后仍会把中心残差推到约 `7%`；
- `1e-9~1e-10` 基本不改变结果，也没有真正稳定空间粗糙度；
- 因此空间平滑目前不能作为物理闭合规则，只能保留为待研究的数值 gauge。

时间平滑更温和：

- tau=0 时中心残差基本不变；
- 时间层间系数跳动可下降约一半；
- 但仍需用 nullspace 版本确认 hard conservation 没有泄漏。

## 4. nullspace 版三切片结果

正式输出：

- `visualizations/equation_first_gbcd_auxiliary_gauge_nullspace_n96_taum3p5_norm_time/`
- `visualizations/equation_first_gbcd_auxiliary_gauge_nullspace_n96_tau0_norm_time/`
- `visualizations/equation_first_gbcd_auxiliary_gauge_nullspace_n96_taup3p5_norm_time/`

三组切片均使用：

- `n=96`
- `fit_region=trusted`
- `force_region=trusted`
- `force_mode=full`
- `hard_solver=nullspace`
- `gauge_cases=norm,time`

核心数值：

| 切片 | 代表元 | 中心代数 weighted mean | 中心代数 p95 | 相邻层 weighted mean | hard constraint max |
|---|---:|---:|---:|---:|---:|
| \(\tau=-3.5\) | norm | `0.02946` | `0.14064` | `0.20495 / 0.16472` | `1.78e-15` |
| \(\tau=-3.5\) | time | `0.03252` | `0.14627` | `0.20683 / 0.16481` | `1.05e-15` |
| \(\tau=0\) | norm | `0.01187` | `0.03887` | `0.15589 / 0.13663` | `9.93e-17` |
| \(\tau=0\) | time | `0.01199` | `0.03887` | `0.15592 / 0.13666` | `8.79e-17` |
| \(\tau=+3.5\) | norm | `0.04647` | `0.22023` | `0.22788 / 0.16832` | `3.64e-16` |
| \(\tau=+3.5\) | time | `0.04692` | `0.22027` | `0.22801 / 0.16851` | `3.43e-16` |

这里 hard constraint max 是缩放后约束行 \(F_s y\) 的最大绝对值。它达到 \(10^{-16}\) 量级，说明 full conservation 在数值上确实是硬满足的。

## 5. 当前判断

本轮可以确定三点：

1. gBCD 张量壳在 full conservation 硬约束下仍然可行，中心代数残差维持在约 `1%~5%` 量级。
2. 之前 tau=±3.5 看到的 divergence 异常主要是 KKT 病态造成的数值泄漏，不应解释为物理失败。
3. `time` 代表元是温和候选，但不是最终物理本构；它只是在解空间中选更平滑的代表元，没有给出唯一的空间状态方程。

因此当前最稳妥的机制仍是：

\[
\boxed{
\mathcal C_{\mu\nu}
=\lambda_I E^I_{\mu\nu},
\quad
\lambda_I=(A,B,C,D),
\quad
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0,
}
\]

再补充一个真正的状态方程/代表元规范来唯一确定 \(\lambda_I\)。

## 6. 下一步

下一步不应再把“逐切片投影”误当成完整动力学。应当转向：

1. 写出 \(\lambda_I\) 的一阶守恒 PDE 主部，并判断哪些分量是约束、哪些分量可作演化。
2. 尝试一个最小状态方程，例如 trace condition、固定本征压力关系、或从辅助场作用量导出的代数约束。
3. 用该闭合系统从初值推进，而不是每个切片独立拟合。
4. 同时保留 nullspace hard constraint 作为每步投影/校正器，防止数值漂移破坏 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\)。

## 7. 守恒 PDE 主部秩检查

为判断 full conservation 是否能独立产生全部 \(A,B,C,D\)，新增脚本：

`kg_examples/diagnose_gbcd_conservation_principal_symbol.py`

对

\[
\tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})=0,
\qquad
\lambda_I=(A,B,C,D)
\]

只取 \(\lambda_I\) 导数的主部。对任意 covector \(k_\mu\)，主符号为

\[
P_{\nu I}(k)
=
\left[
k_\nu,\,
u_\nu(u^\mu k_\mu),\,
r_\nu(r^\mu k_\mu),\,
\frac12\{r_\nu(u^\mu k_\mu)+u_\nu(r^\mu k_\mu)\}
\right].
\]

在当前 `2+1 txz` 约化中，\(P\) 是 \(3\times 4\) 矩阵，因此最多秩为 3。数值检查输出：

`visualizations/equation_first_gbcd_conservation_principal_symbol_n96_trusted/`

三切片 `tau=-3.5,0,+3.5`，三个方向 `lab_t/lab_x/lab_z` 上：

- rank p50/p95/min/max 全部为 `3`；
- nullity p50 全部为 `1`；
- 条件数 p95 在 `6.4e5` 到 `1.86e7` 之间，说明主部不是低秩失败，但有明显病态方向。

结论：

\[
\boxed{
\text{在当前 2+1 约化中，full conservation 是必要条件，但不足以唯一演化四个 } \lambda_I。
}
\]

因此必须额外补充一个状态方程/规范条件。否则 \(A,B,C,D\) 的产生机制仍有一个自由函数，不能称为完整动力学理论。
