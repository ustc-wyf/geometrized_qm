# 高分辨率口径与 full-linear metric update 首轮

日期：2026-05-07

## 问题

用户指出当前结果听不懂，并追问当前分辨率是否太低、难以处理 1550 nm 干涉条纹。

复查后确认：此前 `n=96` 诊断确实是低分辨率，只能做结构性算法测试，不能作为干涉条纹和二阶曲率残差的最终数值证据。

## 分辨率量化

旧设置：

- `full_resolution=96`
- `window=[-9,9] um`
- 实际局部窗口：`29 x 29`
- \(\Delta x\simeq0.617\,\mu m\)
- \(\lambda=1.55\,\mu m\)
- 每波长采样点数：约 `2.51`

这个采样对干涉条纹和二阶导数都明显不足。

新检查：

- 脚本：`kg_examples/render_a_branch_resolution_check.py`
- 输出：`visualizations/a_branch_resolution_check_n384_window9um/`
- 参数：`full_resolution=384`
- 实际局部窗口：`117 x 117`
- \(\Delta x\simeq0.154\,\mu m\)
- 每波长采样点数：约 `10.05`

图：

- `visualizations/a_branch_resolution_check_n384_window9um/a_branch_rho_resolution_check.png`

结论：

- `n=384` 是下一档可用的条纹分辨率；
- 若要更可靠地算曲率二阶导数，后续还应比较 `n=512/640`；
- 但当前 dense 约束/几何求解器不能直接上 `n=384`，否则矩阵规模和内存会爆。

## full-linear metric update

上一轮 pointwise principal correction 的问题是：

- 它逐点求 \(\delta\tilde g_+\)；
- 只精确处理 time-time 主部；
- 直接代回完整几何后会激发 mixed \(tx,tz\) 与空间导数项。

本轮实现全局 full-linear 低分辨率原型：

- 脚本：`kg_examples/solve_gbcd_full_linear_metric_update.py`
- 思路：对已有 Einstein tensor 离散 stencil 关于 \(\tilde g_+\) 做完整线性化；
- 包含：
  - local \(tt\) 项；
  - mixed \(tx,tz\) stencil 项；
  - connection quadratic 的线性项；
- 当前实现是 dense prototype，只能用于低分辨率诊断。

## tau=0 ridge 扫描

在 `n=96,tau=0` 上扫 ridge：

| ridge | linear residual wmean | exact residual wmean | exact residual p95 |
|---:|---:|---:|---:|
| \(10^{-10}\) | \(4.92\times10^{-5}\) | 0.22984 | 1.03152 |
| \(10^{-6}\) | \(4.94\times10^{-5}\) | 0.00661 | 0.04951 |
| \(10^{-4}\) | \(5.12\times10^{-5}\) | 0.00805 | 0.06984 |
| \(10^{-2}\) | 0.00218 | 0.02079 | 0.30999 |

结论：

- `ridge=1e-6` 当前最好；
- 太小的 ridge 会过拟合线性方程，导致 \(\delta\tilde g_+\) 太大，完整非线性代回变差；
- `ridge=1e-6` 能把 tau=0 的 exact residual 加权均值降到约 `0.66%`。

## 三切片 full-linear 结果

采用 `ridge=1e-6`：

| tau | before wmean | linear after wmean | exact after wmean | exact after p95 |
|---:|---:|---:|---:|---:|
| -3.5 | 1.16954 | \(4.10\times10^{-4}\) | 0.06792 | 0.86055 |
| 0 | 1.31280 | \(4.94\times10^{-5}\) | 0.00661 | 0.04951 |
| +3.5 | 0.90248 | \(2.87\times10^{-4}\) | 0.08244 | 0.70280 |

输出：

- `visualizations/equation_first_gbcd_full_linear_metric_update_n96_taum3p5_r1em6/`
- `visualizations/equation_first_gbcd_full_linear_metric_update_n96_tau0_r1em6/`
- `visualizations/equation_first_gbcd_full_linear_metric_update_n96_taup3p5_r1em6/`

解释：

- full-linear global update 明显优于 pointwise principal correction；
- tau=0 干涉中心表现最好；
- tau=\(\pm3.5\) 仍有较大 exact p95，可能来自低分辨率、边缘低密度区、以及当前 dense 原型的正则/边界条件不足。

## 当前判断

1. 当前 `n=96` 不能作为物理可信模拟，只能做算法结构测试。

2. `n=384` 可以解析条纹，是下一步可信诊断的最低档。

3. full-linear metric update 路线是对的，但不能继续用 dense 矩阵形式推进到高分辨率。

4. 下一步应把：
   - principal-constraint projection；
   - full-linear \(\delta\tilde g_+\) update；

   都改成 sparse/matrix-free 版本。

5. 完成后先在 `n=384` 的 `tau=0` 重跑，再扩展到 \(\tau=\pm3.5\)。
