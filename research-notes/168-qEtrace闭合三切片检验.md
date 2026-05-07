# 168-qEtrace闭合三切片检验

日期：2026-05-07

## 一句话结论

\(\mathsf q_E\)-trace 是主符号上自然的标量闭合候选，但它的最简单代数版本在当前高斯干涉三切片数据上不如普通 trace=0。

具体说，下面三类闭合做了对照：

1. 普通 trace：
   \[
   \tilde g^{\mu\nu}\mathcal C_{\mu\nu}=F(Q).
   \]
2. \(u-r\) 二平面正定迹：
   \[
   \mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=F(Q).
   \]
3. \(F(Q)\) 取 \(0\)、线性 \(Q\)、二次 \(Q^2\) 三种最小形式。

结果是：

\[
\boxed{
\text{普通 trace}=0
\text{ 在三切片中给出最小代数残差；}
\mathsf q_E\text{-trace 的零源或简单 }Q\text{ 源版本被当前数据削弱。}
}
\]

因此，若继续保留 \(\mathsf q_E\) 思路，它更适合作为主符号/辅助场正则结构，而不是直接设成

\[
\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=0
\]

这样的硬代数状态方程。

## 1. 本轮检验的问题

上一份笔记发现：投影方程

\[
\Pi_E^\perp\mathcal R_{\mu\nu}=0
\]

加 harmonic gauge 后，仍有一个 \(E\)-方向主部模式没有被固定。

这个未定模式为

\[
\bar h^{(0)}_{\mu\nu}
=
-w_\mu w_\nu,
\]

\[
w_\mu=(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu.
\]

为了固定它，需要一个标量闭合条件。理论上自然的候选是

\[
\mathsf q_E^{\mu\nu}(\Pi_E\mathcal R)_{\mu\nu}
=
\chi(\mathcal Q)\Theta.
\]

本轮检验的问题是：

\[
\boxed{
\text{最简单的 }\mathsf q_E\text{-trace 闭合是否已经被高斯干涉数据支持？}
}
\]

## 2. 数值设置

修改脚本：

`kg_examples/fit_gbcd_trace_closure.py`

新增闭合：

- `qe0`：
  \[
  \mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=0.
  \]
- `qe_q1`：
  \[
  \mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=\alpha Q/Q_0.
  \]
- `qe_q2`：
  \[
  \mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=\alpha Q/Q_0+\beta(Q/Q_0)^2.
  \]

并和旧闭合比较：

- `zero`：
  \[
  \tilde g^{\mu\nu}\mathcal C_{\mu\nu}=0.
  \]
- `q1`：
  \[
  \tilde g^{\mu\nu}\mathcal C_{\mu\nu}=\alpha Q/Q_0.
  \]
- `q2`：
  \[
  \tilde g^{\mu\nu}\mathcal C_{\mu\nu}=\alpha Q/Q_0+\beta(Q/Q_0)^2.
  \]

共同参数：

- `full_resolution=96`
- `window_um=9`
- `fit_region=trusted`
- `force_region=trusted`
- `force_mask_erosion=1`
- `active_dilation=1`
- 三个中心切片：
  \[
  \tau=-3.5,\quad0,\quad+3.5.
  \]

输出目录：

- `visualizations/equation_first_gbcd_qe_trace_closure_n96_taum3p5_trusted/`
- `visualizations/equation_first_gbcd_qe_trace_closure_n96_tau0_trusted/`
- `visualizations/equation_first_gbcd_qe_trace_closure_n96_taup3p5_trusted/`

## 3. 指标定义

`central residual`：

当前中心切片上的代数场方程相对残差加权均值，越小越好。它衡量闭合条件加入后，\(\mathcal C_{\mu\nu}\) 还能不能解释

\[
\mathcal R_{\mu\nu}
=
\tilde G_{\mu\nu}
-\tilde T_{\mu\nu}/M_P^2.
\]

`probe residual`：

中心切片前后两个很近 probe 切片的 residual 平均值。它粗略衡量该闭合在短时间方向上是否稳定。

`rank`：

在采样主符号方向 `lab_t/lab_x/lab_z` 下，守恒主符号加上一个标量闭合后，对 \(\lambda_I=(A,B,C,D)\) 的秩。`rank=3` 表示标量闭合后，剩余三维 nullspace 被守恒方程看见。

## 4. 结果表

### \(\tau=-3.5\)

| closure | central residual | probe residual |
|---|---:|---:|
| trace=0 | 0.1185 | 0.2566 |
| trace=q1 | 0.1235 | 0.2604 |
| trace=q2 | 0.1229 | 0.2598 |
| qe=0 | 0.1876 | 0.2912 |
| qe=q1 | 0.1810 | 0.2848 |
| qe=q2 | 0.1810 | 0.2848 |

### \(\tau=0\)

| closure | central residual | probe residual |
|---|---:|---:|
| trace=0 | 0.0264 | 0.1859 |
| trace=q1 | 2.0156 | 1.6201 |
| trace=q2 | 0.7827 | 0.6677 |
| qe=0 | 0.2964 | 0.3860 |
| qe=q1 | 0.2908 | 0.3816 |
| qe=q2 | 0.2907 | 0.3814 |

### \(\tau=+3.5\)

| closure | central residual | probe residual |
|---|---:|---:|
| trace=0 | 0.1103 | 0.2598 |
| trace=q1 | 0.1170 | 0.2635 |
| trace=q2 | 0.1452 | 0.2906 |
| qe=0 | 0.2208 | 0.3498 |
| qe=q1 | 0.1942 | 0.3245 |
| qe=q2 | 0.1941 | 0.3245 |

## 5. 直接判断

### 5.1 普通 trace=0 目前最好

在三组中心切片中，普通 trace=0 都是当前最小 residual：

\[
0.1185,\quad0.0264,\quad0.1103.
\]

而最好的 \(\mathsf q_E\)-trace 版本为：

\[
0.1810,\quad0.2907,\quad0.1941.
\]

所以最小 \(\mathsf q_E\)-trace 硬闭合没有通过这一轮经验筛选。

### 5.2 \(F(Q)=Q,Q^2\) 不是好 RHS

普通 trace 的 `q1/q2` 在 \(\tau=0\) 明显恶化，尤其 `q1` 把 residual 拉到约 \(2.0\)。

\(\mathsf q_E\)-trace 的 `q1/q2` 只比 `qe0` 略有改善，且拟合出的 \(Q\) 多项式系数很小。这说明简单

\[
F(Q)=\alpha Q+\beta Q^2
\]

没有抓住合适的状态方程。

### 5.3 主符号 rank 没有否定这些闭合

所有测试闭合在采样的 `lab_t/lab_x/lab_z` 主方向上都给出 rank=3。

这意味着：

- 普通 trace=0 在这些方向上也能看见主符号缺口；
- \(\mathsf q_E\)-trace 的优势不是 rank，而是它对 \(-w_\mu w_\nu\) 的收缩符号更稳定；
- 但是这种主符号优势没有自动转化为更好的 algebraic closure。

## 6. 理论修正

上一轮说“普通 trace 未必可靠”是对的，但需要更精确：

普通 trace 对未定模式的收缩为

\[
\tilde g^{\mu\nu}(-w_\mu w_\nu)=-w^2.
\]

所以普通 trace 只在 \(w^2=0\) 时漏掉该模式。

\(\mathsf q_E\)-trace 对该模式为

\[
\mathsf q_E^{\mu\nu}(-w_\mu w_\nu)
=
-
\left[
(e_0\cdot w)^2+(e_1\cdot w)^2
\right],
\]

只要 \(w\neq0\) 就不会漏掉。

因此：

- 普通 trace 不是理论上必然坏；
- \(\mathsf q_E\)-trace 是更强的主符号安全条件；
- 但当前数据更支持普通 trace=0 作为最小代数闭合。

## 7. 当前结论

\[
\boxed{
\mathsf q_E\text{-trace 的最小硬闭合版本应暂时降级。}
}
\]

更准确地说：

1. \(\mathsf q_E\)-trace 仍是有价值的主符号诊断；
2. 它不应该直接取 \(0\) 或简单 \(Q,Q^2\) 多项式；
3. 如果继续使用它，应作为辅助场 action 中的正定范数或稳定性项，而不是唯一代数状态方程；
4. 当前最小 algebraic 状态方程反而是普通 trace=0 更强。

## 8. 普通 trace 的退化面检查

普通 trace=0 的理论风险是它在

\[
w^2=0
\]

时可能漏掉主符号未定模式，其中

\[
w_\mu=(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu.
\]

为此又检查了

\[
\frac{|w^2|}{\mathsf q_E(w,w)}
\]

在 `trusted` 区域的大小。这个比值越接近 0，说明 \(w\) 越接近 null，普通 trace 越容易失效。

输出：

`visualizations/equation_first_gbcd_trace_w2_degeneracy_n96_3tau/summary.json`

结果：

- \(\tau=-3.5\)：约 `12%~14%` 的 trusted 点在采样 `lab_t/lab_x/lab_z` 主方向上接近 \(w^2=0\)；
- \(\tau=0\)：约 `2.7%`；
- \(\tau=+3.5\)：约 `3.7%~5.5%`。

但这些点高度重合于 \(r^\perp\) 退化/近退化区域，也就是 \(\mathsf q_E\) 本身需要换 patch 的区域。

因此更准确的结论是：

\[
\boxed{
\text{普通 trace=0 是当前最小闭合，但必须把 }w^2=0\text{ 或 }r^\perp\text{ 退化区作为 patch 风险。}
}
\]

这不是立即否定 trace=0，而是要求后续理论把这些区域作为图册边界或分支面处理。

## 9. 下一步

下一步有两个合理方向：

1. 沿数据支持方向：把普通 trace=0 与 full conservation 组合成当前最小闭合，重新检查它的主符号退化面 \(w^2=0\) 在高斯干涉区域是否真实出现；
2. 沿作用量方向：把 \(\mathsf q_E\)-trace 放入辅助场正定范数，例如
   \[
   \int\sqrt{|\tilde g|}
   \left(\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}\right)^2,
   \]
   而不是把它硬设为零。
