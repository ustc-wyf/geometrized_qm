# D 支 support 边界短窗与判别式守卫

日期：2026-05-08

## 1. 测试目的

检验当前默认基线 `global + trace0 + penalty + joint2` 是否只能在 `core10` 上稳定，还是也能推进到完整 `support`。

## 2. 测试设置

- `tau = 0`
- `steps = 2` / `10`
- `evolve-region = support`
- 额外再跑一版 `--no-stop-on-negative-discriminant`

## 3. 结果

### support, 2 步

- D residual core10 weighted mean `0.1003`
- support pullback weighted L1 `4.34e-05`
- 没有负判别式

### support, 10 步, 默认硬停

- 在第 8 步因 `negative mass-shell discriminant: -1.037741e-08` 停止
- 负判别式比例只有 `2.36e-4`
- 说明不是大规模崩坏，而是一个很小的局部阈值触发

### support, 10 步, 关闭硬停

- 成功完成 10 步
- D residual core10 weighted mean `0.1002`
- support pullback weighted L1 `2.96e-04`
- 最终负判别式最小值 `-4.45e-08`
- 负判别式比例 `4.71e-4`

## 4. 解释

support 不是不能做，而是：

1. 当前硬停阈值太敏感，会把极小的负判别式当成失败；
2. 关闭硬停后，support 10 步仍保持同一量级的 residual；
3. 所以边界/interface 不是“立刻炸掉”的问题，而是“判别式守卫需要改成柔性门限”的问题。

## 5. 结论

- 当前默认骨架已经能从 `core10` 推到 `support`；
- 下一步 interface 规则不需要从零发明，但需要把负判别式硬停改成更温和的守卫；
- 这比继续加 metric corrector 更重要。

## 6. 代码更新：柔性守卫与时间层同步

已在 `kg_examples/simulate_d_harmonic_standalone.py` 中加入：

- `--disc-fraction-tolerance`：允许少量低密度/边界点出现负质量壳判别式；默认 `1e-3`。
- `--rho-tilde-fraction-tolerance`：允许少量边界点出现负 `rho_tilde`；默认 `1e-3`。
- `--stop-on-negative-rho-tilde/--no-stop-on-negative-rho-tilde`：控制是否让大面积负 `rho_tilde` 直接停止。
- `fields_final.npz` 现在保存 `n_cons` 与 `discriminant`，便于事后定位坏点。
- 修正 `global_auxiliary_solve` 的时间层同步：求下一层 `C` 时使用 `new_minus,new_center,metric_next_work,next_rho,next_current_on_center`，守恒的一步后向时间差分以前一层 `C` 为参照。

这里的“柔性守卫”只是可信性判定，不会裁剪 `rho_tilde`，也不会把负判别式点改成正值；所以它不是改变物理方程，而是避免极少数边界点过早中断整轮诊断。

## 7. 同步版验证

测试：

- `tau=0`
- `steps=10`
- `evolve-region=support`
- `global + trace0 + penalty + joint2`
- `metric-corrector=none`

输出目录：

- `visualizations/d_harmonic_scan_tau0_steps10_support_joint2_penalty_synced_softguard/`

结果：

- 成功完成 10 步，`stopped_reason = completed`。
- D residual core10 weighted mean `0.1002555`。
- D residual core10 p95 `0.3455607`。
- support 上 `rho_pullback_weighted_l1 = 2.9618e-4`。
- core10 上 `rho_pullback_weighted_l1 = 2.0538e-4`。
- support 上负判别式比例 `4.7148e-4`，低于默认容许比例 `1e-3`。
- core10 上负判别式点数为 `0`。
- support 上负 `rho_tilde` 点数为 `1/4242`，core10 上为 `0/1024`。

坏点定位：

- 唯一负 `rho_tilde` 点在 `(x,z)=(-19.5337, 1.5627)` 的 support 边缘；
- 该点 `rho_A=4.755e-6`，不属于 core10 主支撑区；
- 因此当前问题应归类为边界正性/interface 问题，而不是主物理区失控。

## 8. 下一步判断

- 不应降低分辨率来掩盖这个问题；
- 下一步应做正式的边界正性/interface 规则，例如 active-set 边界、低密度边界贴体处理或保守通量约束；
- 同时应优化全局 `C` 求解器，因为 `n=384,support,steps=10` 约需 8 分钟，主要瓶颈是单核全局最小二乘。
