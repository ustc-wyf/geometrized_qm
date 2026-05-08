# D 支 support 边界 active-set 规则

日期：2026-05-08

## 1. 问题

上一轮已经确认：

- `support` 级别可以短窗推进；
- `core10` 主支撑区没有负质量壳判别式和负 `rho_tilde`；
- 坏点集中在低密度 support 边缘带。

本轮继续判断：这些坏点是否来自人工 support 边界的密度通量，还是来自质量壳实根裕度本身耗尽。

## 2. zero-flux 边界对照

新增选项：

- `--matter-boundary-mode open`
- `--matter-boundary-mode zero_flux`

`zero_flux` 的含义：

- 只在数值通量层把 active 与非 active 之间的密度通量设为 0；
- 不裁剪 `rho_tilde`；
- 不修改质量壳方程；
- 不修改内部 active 区的连续性方程。

测试：

- `tau=0`
- `steps=10`
- `evolve-region=support`
- `global + trace0 + penalty + joint2`

结果：

- `zero_flux` 与 `open` 的最终坏点几乎完全相同；
- 最终仍有 support 负判别式 `2/4242`；
- 最终仍有 support 负 `rho_tilde` `1/4242`；
- 坏点坐标仍为 `(x,z)=(-19.5337,1.5627)`。

结论：

- 坏点不是由 active 外“真空通量”主导产生；
- 主要问题应来自低密度边缘带的质量壳实根裕度耗尽，或由此引发的 `j^t` 过零。

## 3. active-set mass-shell guard

新增选项：

- `--active-set-mode fixed`
- `--active-set-mode mass_shell_guard`
- `--active-set-disc-margin`
- `--active-set-rho-frac`

规则：

对满足下列条件的点，从后续 active 演化域中移除：

1. 该点属于当前 active；
2. 该点不属于 `core10`，即不在主支撑保护区；
3. `rho_tilde <= active_set_rho_frac * max(rho_tilde[active])`；
4. 质量壳判别式 `discriminant <= active_set_disc_margin`。

本轮测试取：

- `active_set_disc_margin = 1e-7`
- `active_set_rho_frac = 5e-3`

解释：

- 这是计算域 active-set 规则；
- 它不是裁剪 `rho_tilde`；
- 也不是把负判别式改成正；
- 它的物理口径是：低密度、非主支撑、质量壳实根裕度不足的点不再作为当前离散演化域的一部分。

## 4. active-set 测试结果

输出目录：

- `visualizations/d_harmonic_scan_tau0_steps10_support_joint2_penalty_active_guard/`

结果：

- 成功完成 10 步；
- active 总点数从 `4242` 降到 `4239`；
- 总共移除 `3` 个点；
- 被移除点全部不在 `core10`；
- `core10` 缺失点数为 `0`；
- active 区最终负判别式点数为 `0`；
- active 区最终负 `rho_tilde` 点数为 `0`。

关键数字：

- D residual core10 weighted mean `0.1002552`；
- D residual core10 p95 `0.3455522`；
- support `rho_pullback_weighted_l1 = 2.1369e-4`；
- core10 `rho_pullback_weighted_l1 = 2.0539e-4`；
- active 最小判别式 `8.8909e-8`；
- active 最小 `rho_tilde = 2.3082e-6`。

被移除的 3 个点：

- `(x,z)=(-19.5337,0.7813)`, `rho_A=4.069e-6`
- `(x,z)=(-19.5337,1.5627)`, `rho_A=4.755e-6`
- `(x,z)=(-2.3440,-22.6591)`, `rho_A=3.331e-6`

## 5. 判断

active-set mass-shell guard 是当前更合理的 support 边界/interface 规则：

- 它保持 `core10` 完整；
- 它不改变主物理区结果；
- 它清除了 active 区负判别式和负 `rho_tilde`；
- 它比单纯 `zero_flux` 更针对真正失败机制。

但它仍是数值边界规则，不是最终理论方程。

下一步：

- 把 `mass_shell_guard` 作为 support 长窗默认边界候选；
- 继续检查阈值 `active_set_disc_margin` 与 `active_set_rho_frac` 的敏感性；
- 同时优化全局 `C` 求解器，否则每轮 10 步仍约 8 到 10 分钟。
