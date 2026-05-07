# trace0 硬守恒 KKT 与 ALM 原型

## 背景

上一轮已经确认，在 \(n=384\)、`core10`、三时间切片上，普通 trace=0 与完整守恒并不矛盾：

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu},
\qquad
C_{\mu\nu}\in
\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\},
\]

\[
\tilde g^{\mu\nu}C_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0.
\]

但 `force_weight` 扫描只把守恒作为 penalty 行加入，还不是等式硬约束。本轮目标是检查能否把完整守恒更接近地实现为硬约束，同时不破坏代数 gBCD 壳。

## 方法

更新脚本：

- `kg_examples/fit_gbcd_trace0_sparse_conservation.py`

新增三类求解模式：

- `--hard-project`：从逐点 trace0 初始代表元出发，只求 \(F\delta=-Fx_0\)。它能直接压守恒，但不最小化代数残差。
- `--hard-kkt`：求解约束最小二乘的 KKT 方程
  \[
  \min_\delta\|A(x_0+\delta)-b\|^2+\epsilon\|\delta\|^2,
  \qquad
  F(x_0+\delta)=0.
  \]
- `--hard-alm`：用增广拉格朗日/乘子法做同一个硬约束问题，把每个子步化成正定 least-squares：
  \[
  \min_\delta
  \|A\delta-(b-Ax_0)\|^2
  +\mu_k\|F\delta+Fx_0+\lambda_k/\mu_k\|^2,
  \]
  \[
  \lambda_{k+1}=\lambda_k+\mu_k F(x_0+\delta_k).
  \]

这里 \(A\) 是三切片代数场方程行，\(F\) 是完整守恒行。trace=0 已经在逐点系数 nullspace 中硬消元。

## 关键结果

小规模烟测 `n=96,tau=0`：

- `hard-kkt` 可以把 full divergence 压到约 \(3.6\times10^{-13}\)，中心代数残差保持在约 \(6.6\times10^{-4}\)。
- `hard-alm` 可以把 full divergence 压到约 \(4.8\times10^{-7}\)，中心代数残差同样保持在约 \(6.6\times10^{-4}\)。

高分辨率三切片 `n=384,core10,window=[-9,9]um`：

| tau | mode | center algebraic row residual | full divergence scale |
|---:|---|---:|---:|
| -3.5 | penalty:10 | \(5.478\times10^{-3}\) | \(7.740\times10^{-3}\) |
| -3.5 | hard-alm | \(6.009\times10^{-3}\) | \(1.987\times10^{-3}\) |
| 0 | penalty:10 | \(6.747\times10^{-4}\) | \(4.576\times10^{-4}\) |
| 0 | hard-alm | \(6.747\times10^{-4}\) | \(2.388\times10^{-5}\) |
| +3.5 | penalty:10 | \(7.698\times10^{-3}\) | \(6.502\times10^{-3}\) |
| +3.5 | hard-alm | \(8.172\times10^{-3}\) | \(1.940\times10^{-3}\) |

主要输出：

- `visualizations/equation_first_gbcd_trace0_hard_alm_n384_taum3p5_core10/`
- `visualizations/equation_first_gbcd_trace0_hard_alm_n384_tau0_core10/`
- `visualizations/equation_first_gbcd_trace0_hard_alm_n384_taup3p5_core10/`

## 解释

`hard-project` 不是可接受的物理/数值代表元选择：它只顾守恒，可能显著破坏代数场方程。

`hard-kkt` 在小规模和干涉中心表现很好，但在左侧分离态 `tau=-3.5,n=384` 中，当前无 SciPy/MINRES 的 LSQR saddle-point 实现仍然病态，可能为降低 KKT 残差牺牲代数项。因此不能把这一失败解释为理论矛盾。

`hard-alm` 更稳，因为它避免直接求解不定 KKT 系统。三切片结果显示：完整守恒可以比 penalty 扫描进一步压低约 \(3\sim20\) 倍，而代数残差只小幅增加。因此当前数据支持把 trace0 + full conservation 写成 equation-first 候选的硬约束结构，而不是仅作为 penalty 正则化。

## 当前理论判断

当前最小显式方程候选可继续写为：

\[
\boxed{
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+C_{\mu\nu},
\quad
C_{\mu\nu}\in E,
\quad
C^\mu{}_\mu=0,
\quad
\tilde\nabla^\mu C_{\mu\nu}=0,
}
\]

其中

\[
E=\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}.
\]

这一步仍不是 action 化成功，也不是完整演化理论完成；它说明在当前高斯干涉试验场上，`trace0 + full conservation` 作为显式 equation-first 场方程没有被数值反例否定。

## 下一步

1. 把本候选整理成一版可直接写入论文草稿的 equation-first proposal。
2. 明确 patch 条件：\(\Delta=u^2r^2-(u\cdot r)^2\neq0\)、\(r^\perp\) 不退化、以及 \(w^2\) 退化面上的分支处理。
3. 若继续求严格 saddle-point 数值解，需要 MINRES/Schur complement 或更强预条件器；当前 ALM 已足够作为方程可行性的三切片诊断。
4. 继续理论检查：该 equation-first 系统的约束传播、\(Q\to0\Rightarrow C\to0\) 分支条件、以及是否可能由更高层辅助场 action 产生。
