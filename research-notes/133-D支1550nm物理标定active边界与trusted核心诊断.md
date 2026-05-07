# D支 1550nm 物理标定 active 边界与 trusted 核心诊断

日期：2026-05-02

## 背景

本轮继续研究 1550nm massive KG 近似下的 D 支三域全动力学原型：

- \(\lambda=1550\,{\rm nm}\)
- \(m=\omega_0/10\)
- KG 正频内积归一化
- \(\ell/l_P=10^{60}\)
- `n=64`
- `dt_old=2.5e-7`
- 几何闭合仍为 `inert_tridomain`

注意：`inert_tridomain` 仍只是第一版可执行闭合，不是最终完整张量几何求解器。

## Active Dilation 扫描

输出目录：

- `visualizations/d_physical_1550nm_kg_active_dilation_sweep_ell1e60_n64_dt2p5em7/`

新增决策图：

- `visualizations/d_physical_1550nm_kg_active_dilation_sweep_ell1e60_n64_dt2p5em7/active_dilation_sweep_decision_plot.png`

结果：

- `active_dilation=0`：`step=98` 判别式穿零，\(t=1.2096\times10^{-4}\,{\rm fs}\)
- `active_dilation=1`：完成 `180` 步
- `active_dilation=2`：完成 `180` 步，`measure_rel_l1_support=5.53e-3`，`n_cons_rel_l1_support=1.81e-4`
- `active_dilation=4`：与 `2` 基本相同
- `active_dilation=6`：step 14 非有限值，说明过度扩张会把病态低密度尾部纳入活跃演化

判断：

- 原先的早期质量壳穿零主要不是核心物理区失败，而是 active/support 截断边界造成的数值伪故障。
- 合理窗口是 `active_dilation=2` 或 `4`。
- 过度扩张不是好事，会引入低密度尾部病态。

## 原始 Support 停止条件的失败

用 `active_dilation=2` 延长到 `400` 步，仍以原始 support 作为停止条件：

- 输出目录：
  `visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_steps400_dt2p5em7/`
- 时间诊断：
  `visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_steps400_dt2p5em7/longrun_time_diagnostics.png`
- 位置诊断：
  `visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_steps400_dt2p5em7/failure_location_diagnostic.png`

运行停止在：

- `step=231`
- \(t=4.3318\times10^{-4}\,{\rm eV}^{-1}=2.8512\times10^{-4}\,{\rm fs}\)
- `disc_min_support=-1.535e-10`
- `measure_rel_l1_support=4.675`

但前一步：

- `step=230`
- `disc_min_support=1.742e-10`
- `measure_rel_l1_support=3.838e-2`
- `n_cons_rel_l1_support=2.553e-4`

最小判别式位置：

- grid index `(42,28)`
- \((x,z)=(9.251,-3.700)\,\mu{\rm m}\)
- 该点在原始 support 的边缘
- 它不在 active 边缘，说明 `active_dilation=2` 已经修复了活跃区截断，但原始 support 的阈值边缘仍会触发停止

在同一个失败时刻，若把 support 向内腐蚀一格得到 trusted core：

- `disc_min_trusted=3.865e-6`
- `measure_rel_l1_trusted=8.30e-5`

## 新增 Trusted Support 机制

修改脚本：

- `kg_examples/simulate_d_tridomain_full_dynamics.py`

新增参数：

- `--trusted-erosion N`：把原始 support 按 8 邻域向内腐蚀 `N` 格，得到 trusted core
- `--stop-mask support|trusted`：停止条件选择原始 support 或 trusted core

原始 support 指标不会删除；trusted 只用于把人为阈值边缘伪故障与主体演化分开。

## Trusted 核心区 400 步结果

运行设置：

- `active_dilation=2`
- `trusted_erosion=1`
- `stop_mask=trusted`
- `steps=400`

输出目录：

- `visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_trusted1_steps400_dt2p5em7/`

关键图：

- `visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_trusted1_steps400_dt2p5em7/trusted_vs_support_time_diagnostics.png`
- `visualizations/d_physical_1550nm_kg_ell1e60_n64_dilate2_trusted1_steps400_dt2p5em7/trusted_final_spatial_diagnostic.png`

最终：

- 完成 `400` 步
- \(t=7.50096\times10^{-4}\,{\rm eV}^{-1}=4.9372\times10^{-4}\,{\rm fs}\)
- `stopped_reason=completed`

原始 support 外圈：

- `measure_rel_l1_support=17.981`
- `n_cons_rel_l1_support=6.71e-4`
- `disc_min_support=-5.554e-8`

trusted 核心区：

- `measure_rel_l1_trusted=1.437e-4`
- `n_cons_rel_l1_trusted=9.622e-5`
- `disc_min_trusted=3.867e-6`
- `rho_tilde_max_trusted=1.16603e-3`，几乎保持初值尺度

空间定位：

- support 最小判别式仍在 `(42,28)`，\((9.251,-3.700)\,\mu{\rm m}\)，属于 support 边缘
- trusted 最小判别式在 `(41,27)`，\((8.326,-4.625)\,\mu{\rm m}\)，仍为正

## 判断

\[
\boxed{
\text{当前 D 支原型的早期失败主要来自 support/active 阈值边缘处理，而不是 trusted 核心区主体演化立即失效。}
}
\]

更细分地说：

- `active_dilation=0` 的失败是 active 截断边界伪故障
- `active_dilation=2` 修复 active 边界后，原始 support 的外圈阈值仍会出现局部判别式穿零和测度尖峰
- 把停止条件改成一格内缩 trusted core 后，主体区可以稳定完成 400 步，且 D-vs-A 误差仍在 \(10^{-4}\) 量级

这不是在物理上删去低密度区，而是在数值诊断上区分：

- 主体支撑区内可信演化
- 人为阈值边缘和低密度尾部的病态处理

## 下一步

建议下一步先做受控收敛：

- 固定 `active_dilation=2, trusted_erosion=1, stop_mask=trusted`
- 做 `n=64/96` 或 `64/96/128` 分辨率检查
- 做 `dt_old=2.5e-7/1.25e-7` 时间步检查
- 同时报告 support 与 trusted 两套指标
- 若 trusted 核心稳定而 support 外圈持续坏掉，则下一步应改成平滑权重/tapered support 或 body-fitted 低密度边界，而不是把 support 外圈当作 D 支物理失败
