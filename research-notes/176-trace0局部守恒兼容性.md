# trace=0 局部守恒兼容性检查

日期：2026-05-08

## 目的

`175` 号笔记确认了在 `n=384` 高分辨率下，普通 trace=0 比 \(\mathsf q_E\)-trace 更适合作为当前最小标量闭合。但那只是逐点局部代数检查，没有施加完整守恒条件

\[
\tilde\nabla^\mu C_{\mu\nu}=0.
\]

本轮检查的问题是：

> 如果逐点求出满足 trace=0 的 \(C_{\mu\nu}\)，它是否会自动或近似满足完整守恒？

如果答案是“是”，说明局部闭合已经很接近完整方程；如果答案是“否”，说明完整守恒必须和 trace=0 联立求解，不能被当作事后诊断项。

## 方法

新增脚本：

`kg_examples/diagnose_gbcd_trace_local_conservation.py`

对每个中心切片 \(\tau=-3.5,0,+3.5\)，构造三张近邻时间切片

\[
\tau-\delta\tau,\qquad \tau,\qquad \tau+\delta\tau,
\]

其中 \(\delta\tau=2.5\times 10^{-7}\)。然后在每个点逐点求解

\[
C_{\mu\nu}
=
\lambda_g\tilde g_{\mu\nu}
+\lambda_u u_\mu u_\nu
+\lambda_r r_\mu r_\nu
+\lambda_{ur}u_{(\mu}r_{\nu)}
\]

并比较：

- `unconstrained`：不加 trace=0，只做局部最小二乘。
- `trace0`：逐点强制 \(\tilde g^{\mu\nu}C_{\mu\nu}=0\)。

然后用三切片有限差分直接计算中心切片的

\[
J_\nu=\tilde\nabla^\mu C_{\mu\nu}.
\]

主要无量纲诊断量为

\[
\frac{|J|}{|C|/h},
\]

其中 \(h=\min(\Delta t,\Delta x,\Delta z)\)。它表示散度相对于“一个网格导数尺度”的大小；若远小于 1，说明自然守恒性较好；若接近或大于 1，说明逐点代表元和守恒严重不兼容。

本轮使用 `n=384, core10, window=[-9,9] um`。

## 结果

加权平均结果如下：

| closure | tau | algebraic residual | full divergence scale | transverse divergence scale |
|---|---:|---:|---:|---:|
| trace0 | -3.5 | 0.004064 | 0.783756 | 1.385058 |
| trace0 | 0 | 0.000598 | 0.000937 | 0.000654 |
| trace0 | +3.5 | 0.006150 | 0.254636 | 1.080004 |
| unconstrained | -3.5 | 0.000648 | 0.023113 | 0.167660 |
| unconstrained | 0 | 0.000069 | 0.000236 | 0.000454 |
| unconstrained | +3.5 | 0.001647 | 0.007245 | 0.394867 |

输出文件：

- `visualizations/equation_first_gbcd_trace_local_conservation_n384_core10_trace0_3tau/summary.json`
- `visualizations/equation_first_gbcd_trace_local_conservation_n384_core10_trace0_3tau/trace_local_conservation_summary.png`
- `visualizations/equation_first_gbcd_trace_local_conservation_n384_core10_unconstrained_3tau/summary.json`
- `visualizations/equation_first_gbcd_trace_local_conservation_n384_core10_unconstrained_3tau/trace_local_conservation_summary.png`

## 解释

1. 干涉中点 \(\tau=0\) 很好：trace0 的代数残差约 \(6.0\times10^{-4}\)，自然 full divergence 约 \(9.4\times10^{-4}\)。这说明在干涉中心，局部 trace0 代表元已经几乎自动守恒。

2. 两个分离态不好：trace0 的 full divergence scale 在 \(\tau=-3.5\) 约 `0.78`，在 \(\tau=+3.5\) 约 `0.25`；transverse divergence scale 甚至接近或超过 1。这不是小误差。

3. `unconstrained` 的自然散度在分离态明显更小，尤其 full divergence scale：

\[
0.023\quad(\tau=-3.5),\qquad
0.0072\quad(\tau=+3.5).
\]

这说明分离态的散度问题不是“逐点局部拟合必然导致”，而是 trace=0 这个额外闭合把局部代表元推离了自然守恒方向。

4. 这不否定 trace=0。它说明 trace=0 不能以“先逐点求、再检查散度”的方式实现；正确方程必须把

\[
\tilde g^{\mu\nu}C_{\mu\nu}=0
\]

和

\[
\tilde\nabla^\mu C_{\mu\nu}=0
\]

联立起来共同选取 \(C_{\mu\nu}\)。

## 当前结论

trace=0 仍是当前主标量闭合，但完整理论必须写成联立方程组：

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2
+C_{\mu\nu},
\]

\[
C_{\mu\nu}\in
\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\},
\]

\[
\tilde g^{\mu\nu}C_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0.
\]

其中最后两个条件不能分步处理。特别是在分离态，完整守恒是确定 \(C_{\mu\nu}\) 的实质方程，而不是数值后处理。

## 下一步

旧的全局 hard-nullspace 求解器在 `n=384` 下太重。下一步应实现一个高分辨率可用的 sparse/global 求解器：

1. 对 trace=0 做逐点硬消元，即每个点把四个系数降到三个自由度。

2. 在消元后的自由度上联立代数残差和完整守恒残差。

3. 先用 sparse LSQR 做 force-weight 扫描，判断守恒能否压低而不显著破坏代数残差。

4. 如果扫描显示可行，再升级为真正 hard conservation 或 saddle-point/KKT 的矩阵自由求解器。

这一步是方程可行性检查，不是改回“构造 D 支模拟器”的目标。
