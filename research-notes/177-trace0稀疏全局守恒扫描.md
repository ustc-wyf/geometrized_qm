# trace=0 稀疏全局守恒扫描

日期：2026-05-08

## 目的

`176` 号笔记说明：逐点 trace=0 代表元在干涉中点几乎自然守恒，但在两个分离态不自然守恒。因此本轮要检查更关键的问题：

> trace=0 和完整守恒是否真的互相兼容？

如果二者兼容，那么完整方程组

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu},
\quad
C_{\mu\nu}\in\mathrm{span}\{\tilde g,uu,rr,ur\},
\]

\[
\tilde g^{\mu\nu}C_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0
\]

就仍是可行的 equation-first 主候选。

## 方法

新增脚本：

`kg_examples/fit_gbcd_trace0_sparse_conservation.py`

实现方式：

1. 对每个格点逐点构造 trace 行

\[
t_I=\tilde g^{\mu\nu}E^I_{\mu\nu},
\qquad
E^I_{\mu\nu}=\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}.
\]

2. 用 SVD 求 \(t_I\lambda^I=0\) 的三维 nullspace，从而把四个 \(\lambda_I\) 硬消元成三个自由度。

3. 在消元后的自由度上建立 sparse 线性系统，包含：

- 代数残差：\(C_{\mu\nu}\approx \tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2\)。
- 守恒残差：\(\tilde\nabla^\mu C_{\mu\nu}\approx0\)。

4. 扫描守恒行权重 `force_weight`。这一步仍是 penalty scan，不是最终 hard conservation saddle-point 求解。

本轮使用 `n=384, core10, window=[-9,9] um`，三切片：

\[
\tau=-3.5,\quad 0,\quad +3.5.
\]

主要报告两个量：

- `row_alg`：对称 6 分量 row-norm 的 floor-relative 代数残差，加权平均。
- `full_div`：\(|\tilde\nabla C|/(|C|/h)\) 的加权平均。

## 结果

左侧分离态 \(\tau=-3.5\)：

| force_weight | row_alg | full_div | transverse_div |
|---:|---:|---:|---:|
| 0 | 0.004797 | 0.783756 | 1.385058 |
| 0.3 | 0.005132 | 0.454681 | 1.918270 |
| 1 | 0.005380 | 0.271785 | 1.286534 |
| 3 | 0.005176 | 0.014175 | 0.215292 |
| 10 | 0.005478 | 0.007740 | 0.224111 |

干涉中点 \(\tau=0\)：

| force_weight | row_alg | full_div | transverse_div |
|---:|---:|---:|---:|
| 0 | 0.000675 | 0.000937 | 0.000654 |
| 1 | 0.000675 | 0.000925 | 0.000646 |
| 10 | 0.000675 | 0.000458 | 0.000341 |

右侧分离态 \(\tau=+3.5\)：

| force_weight | row_alg | full_div | transverse_div |
|---:|---:|---:|---:|
| 0 | 0.007307 | 0.254636 | 1.080004 |
| 0.3 | 0.007465 | 0.040685 | 0.578997 |
| 1 | 0.007429 | 0.040344 | 0.563219 |
| 3 | 0.007534 | 0.022098 | 0.424603 |
| 10 | 0.007698 | 0.006502 | 0.339096 |

输出文件：

- `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_taum3p5_core10_scan/summary.json`
- `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_taum3p5_core10_scan/trace0_sparse_conservation_scan.png`
- `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_tau0_core10_scan/summary.json`
- `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_tau0_core10_scan/trace0_sparse_conservation_scan.png`
- `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_taup3p5_core10_scan/summary.json`
- `visualizations/equation_first_gbcd_trace0_sparse_conservation_n384_taup3p5_core10_scan/trace0_sparse_conservation_scan.png`

## 解释

1. 在两个分离态，逐点 trace0 的自然守恒残差很大，但加入全局守恒 penalty 后可以被显著压低。

2. 左侧分离态最明显：

\[
{\rm full\_div}: 0.784\to0.00774,
\]

而中心切片代数残差只从

\[
0.00480\to0.00548.
\]

这说明 trace0 和完整守恒不是矛盾条件。

3. 右侧分离态也类似：

\[
{\rm full\_div}: 0.255\to0.00650,
\]

代数残差只从

\[
0.00731\to0.00770.
\]

4. 干涉中点本来就几乎自然守恒，加入守恒 penalty 后代数残差基本不变。

5. `transverse_div` 没有像 `full_div` 一样稳定压到很小，原因是这里的 transverse projection 使用 \(\tilde g\) 和 \(u\) 定义，而数值范数是普通欧氏分量范数；它不是严格正交投影范数，所以不应单独作为本轮否定标准。下一步 hard solve 时仍要直接监控完整 \(J_\nu=\tilde\nabla^\mu C_{\mu\nu}\) 的各分量。

## 当前结论

高分辨率三切片测试支持：

\[
\boxed{
\text{trace0 + full conservation 是兼容的。}
}
\]

也就是说，`176` 中看到的分离态大散度不是理论矛盾，而是因为“逐点 trace0 后处理”没有联立守恒。全局联立后，守恒残差可以压低，同时代数残差只轻微增加。

因此当前主候选仍是：

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu},
\]

\[
C_{\mu\nu}\in
\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\},
\]

\[
\tilde g^{\mu\nu}C_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0,
\]

再加 patch 条件和 \(Q\to0\Rightarrow C_{\mu\nu}\to0\) 的分支/边界条件。

## 下一步

1. 把当前 penalty scan 升级为真正的 hard conservation 或 saddle-point/KKT 矩阵自由求解器。

2. 记录 trace0 消元后的主部和 Cauchy 结构，确认 hard conservation 不会引入隐藏的病态自由度。

3. 从理论写作上，把当前方程组整理成显式 proposal：变量、方程、patch、分支条件、测地线兼容条件、与 action 化 no-go 的关系。

4. action 化仍不作为短期主线；若继续，必须走真实辅助场 \(\chi\)-sector 的有效应力路线。
