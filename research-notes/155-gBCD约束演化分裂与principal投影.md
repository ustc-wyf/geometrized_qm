# gBCD 约束-演化分裂与 principal 投影

日期：2026-05-07

## 背景

上一轮 `lambda-only` 一步推进显示：

- 固定 A 支背景，只调整 \(\lambda_I=(A,B,C,D)\)，可以满足 full conservation；
- 但下一切片代数 gBCD 场方程残差仍约 `14%~17%`；
- 因此不能把 D 支动力学简化为固定背景上的 \(\lambda\) 推进。

重新检查二阶几何演化后发现：真正的 metric 方程不应把“下一切片代数场方程”作为显式一步目标。场方程在中心切片含有 \(\tilde g_+\) 的二阶时间差分：

\[
\tilde G_{\mu\nu}[\tilde g_-,\tilde g_0,\tilde g_+]
=
\tilde T_{\mu\nu}/M_P^2+\mathcal C_{\mu\nu}.
\]

因此正确分裂是：

- 当前切片的约束组合必须由当前数据和 \(\lambda\) 硬满足；
- 可演化组合由 \(\tilde g_+\) 的二阶时间导数决定。

## 主部矩阵

对中心切片线性化 \(\tilde G_{\mu\nu}\)，只保留 \(\delta\tilde g_+\) 通过二阶 lab-time 导数进入的主部：

\[
\delta G_{\mu\nu}^{(tt)}
=
\delta R_{\mu\nu}^{(tt)}
-\frac12\tilde g_{\mu\nu}\delta R^{(tt)}.
\]

其中

\[
\delta R_{\mu\nu}^{(tt)}
=\frac12\left[
\delta_{\mu0}\tilde g^{0\beta}\ddot h_{\beta\nu}
+\delta_{\nu0}\tilde g^{0\beta}\ddot h_{\beta\mu}
-\tilde g^{00}\ddot h_{\mu\nu}
-\delta_{\mu0}\delta_{\nu0}\tilde g^{\alpha\beta}\ddot h_{\alpha\beta}
\right],
\]

\[
\ddot h_{\mu\nu}\simeq \frac{\delta \tilde g_{+\,\mu\nu}}{\Delta t^2}.
\]

在当前 1550 nm 参数下：

\[
\Delta t/\Delta x \simeq 6\times10^{-7},
\]

所以该主部是显式时间推进中最强的局部项。

数值发现：该 \(6\times6\) 对称张量主部矩阵在三组切片的 trusted 区域逐点 rank=3。

解释：

- rank=3 的像空间是可由 \(\tilde g_+\) 加速度推进的演化组合；
- 剩下 3 个左零空间方向是约束组合；
- 约束组合不能靠未来一层 \(\tilde g_+\) 补救，必须在当前切片由初值投影/约束求解满足。

## 第一轮主部修正诊断

脚本：

- `kg_examples/diagnose_gbcd_metric_principal_correction.py`

输出：

- `visualizations/equation_first_gbcd_metric_principal_correction_n96_taum3p5/`
- `visualizations/equation_first_gbcd_metric_principal_correction_n96_tau0/`
- `visualizations/equation_first_gbcd_metric_principal_correction_n96_taup3p5/`

结果：

| tau | before weighted mean | after weighted mean | before p95 | after p95 |
|---:|---:|---:|---:|---:|
| -3.5 | 0.03506 | 0.02412 | 0.17693 | 0.10354 |
| 0 | 0.01241 | 0.01229 | 0.03928 | 0.03927 |
| +3.5 | 0.05017 | 0.04258 | 0.24125 | 0.15480 |

结论：若仍使用原来的 v0 \(\lambda\)，主部 \(\delta\tilde g_+\) 只能修正 rank=3 的演化组合，不能消掉约束残差。

## principal-constraint projection

于是改为先选取 \(\lambda\)，使中心残差本身落在可演化子空间：

\[
\mathcal R_{\mu\nu}
=
\mathcal C_{\mu\nu}
-
\left(\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2\right).
\]

令 \(N_A^{\mu\nu}\) 为主部矩阵的左零空间基，则硬约束为：

\[
N_A^{\mu\nu}\mathcal R_{\mu\nu}=0.
\]

同时保持 full conservation 硬约束：

\[
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
\]

然后在这些硬约束解空间中最小化原来的辅助场代表元：

\[
\lambda^2+(\nabla_t\lambda)^2.
\]

脚本：

- `kg_examples/fit_gbcd_principal_constraint_projection.py`

输出：

- `visualizations/equation_first_gbcd_principal_constraint_projection_n96_taum3p5/`
- `visualizations/equation_first_gbcd_principal_constraint_projection_n96_tau0/`
- `visualizations/equation_first_gbcd_principal_constraint_projection_n96_taup3p5/`

## principal-constraint projection 结果

| tau | hard equality max | before projection wmean | after projection wmean | after p95 | relative \(|\delta\tilde g_+|\) p95 |
|---:|---:|---:|---:|---:|---:|
| -3.5 | \(9.63\times10^{-14}\) | 1.16954 | \(3.57\times10^{-4}\) | 0.00140 | \(4.42\times10^{-11}\) |
| 0 | \(1.28\times10^{-14}\) | 1.31280 | \(4.92\times10^{-5}\) | 0.00024 | \(1.20\times10^{-9}\) |
| +3.5 | \(5.18\times10^{-14}\) | 0.90248 | \(2.79\times10^{-4}\) | 0.00162 | \(2.69\times10^{-10}\) |

这里：

- `before projection` 是 \(\lambda\) 投影后、但尚未用 \(\delta\tilde g_+\) 演化组合修正前的中心残差；
- `after projection` 是把可演化残差交给主部 \(\delta\tilde g_+\) 后剩余的约束残差；
- `relative |delta gtilde_+|` 是所需下一层度规修正相对当前度规范数的估计。

## 结论

1. 当前 gBCD 不是简单的每切片代数拟合问题，而是标准的“约束 + 演化”结构。

2. 之前 v0 \(\lambda\) 的中心残差小，并不等于适合作为初值；它还有不可由 \(\tilde g_+\) 推进消去的约束残差。

3. principal-constraint projection 给出了更合理的初值投影规则：
   - full conservation 硬满足；
   - metric 主部左零空间约束硬满足；
   - 剩余演化组合由 \(\delta\tilde g_+\) 负责。

4. 三组切片都跑通，硬约束残差约 \(10^{-13}\sim10^{-14}\)，投影后残差加权均值约 \(5\times10^{-5}\sim3.6\times10^{-4}\)。

5. 所需 \(\delta\tilde g_+\) 相对修正很小，p95 约 \(10^{-11}\sim10^{-9}\)。这说明当前路线数值上有希望。

## 完整几何代回复查

上表的 `after projection` 是只按 time-time 主部计算的残差。为了检查它是否已经等价于完整几何推进，我又把

\[
\tilde g_+^{\rm new}
=
\tilde g_+^{A}+\delta\tilde g_+
\]

直接代回完整的 `metric_jets_full -> geometry_data -> Einstein tensor` 计算。

结果：

| tau | exact after weighted mean | exact after p95 |
|---:|---:|---:|
| -3.5 | 0.34074 | 1.12551 |
| 0 | 0.15291 | 0.58780 |
| +3.5 | 0.03685 | 0.03253 |

解释：

- principal 投影本身是正确的约束/演化分裂诊断；
- 但逐点独立求出的 \(\delta\tilde g_+\) 并不是完整几何更新；
- 它会激发 mixed \(tx,tz\) 和低阶空间导数项；
- 简单邻点平滑会破坏主部方程，反而使 exact residual 变差。

因此不能把当前 pointwise principal correction 当作最终一步演化器。

下一步必须求解一个全局的 \(\delta\tilde g_+\) 方程，把：

- time-time 主部；
- mixed \(tx,tz\) 项；
- connection quadratic 的线性项；
- 必要的空间正则或边界条件；

放进同一个线性/半隐式系统中。

## 下一步

下一版真正一步演化器应使用以下流程：

1. 在当前切片求 \(\lambda_0\)，硬满足：
   \[
   \tilde\nabla^\mu\mathcal C_{\mu\nu}=0,
   \qquad
   N_A^{\mu\nu}\mathcal R_{\mu\nu}=0.
   \]

2. 用剩余可演化残差求全局光滑的 \(\delta\tilde g_+\)：
   \[
   P_{\rm dyn}\mathcal R_{\mu\nu}
   =
   \delta G_{\mu\nu}^{\rm full-linear}[\delta\tilde g_+].
   \]

3. 物质场 \(\rho,S\) 不能再固定为 A 支背景；应接入变换后连续性方程与质量壳方程，作为同一步系统中的物质推进。

4. 当前主部只保留 \(tt\) 项。下一步需要把 mixed \(tx,tz\) 与低阶 connection 项纳入全局 solve；它们不应改变 rank=3 的约束/演化基本分裂，但会影响真正可用的 \(\tilde g_+\)。
