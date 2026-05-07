# tau=-3.5 残差来源与 center-plus 初值投影检验

日期：2026-05-07

## 背景

上一轮 `n=384, core10` 三切片结果显示：

- `tau=0` 干涉中心：完整非线性残差加权均值约 `0.00644`。
- `tau=+3.5` 后分离态：约 `0.03864`。
- `tau=-3.5` 前分离态：约 `0.14577`，是当前主卡点。

本轮目标是判断 `tau=-3.5` 的高残差来自哪里，以及是否可以通过改变辅助张量基底或初值投影方式解决。

## 1. 残差来源分层

新增脚本：

`kg_examples/diagnose_gbcd_metric_update_residual_sources.py`

它把完整非线性残差按以下口径分层：

- `core50`：拟合 mask 内满足 `rho >= 0.5 rho_max` 的高密度主支撑区。
- `core10_only`：`0.1 rho_max <= rho < 0.5 rho_max`。
- `boundary_layer_0`：拟合 mask 内至少有一个 8-邻居在 mask 外的边缘格点。
- `det_abs` 分位区间。
- \(|\delta g|/|g|\) 分位区间。
- `tt, tx, tz, xx, xz, zz` 六个张量分量。

默认 `matter4`、`force-mask-erosion=1` 的结果：

- 全 mask exact residual weighted mean `0.14577`。
- `boundary_layer_0` weighted mean `0.3092`，明显最差。
- `core50` weighted mean `0.1124`，仍不够低。
- 低/中等 \(|det|\) 区更差，高 \(|det|\) 区更好。
- 各张量分量均有贡献，不是单一分量爆炸。
- 高残差不是由 \(|\delta g|/|g|\) 最大的点主导。

结论：边缘层确实是坏点集中区，但不是唯一根因；不能只靠裁 mask 或简单阻尼解释。

## 2. 边缘 force 约束测试

将 force mask erosion 从 `1` 改为 `0`，也就是把 `core10` 边缘点也纳入守恒约束。

结果：

- projection force rows 从 `6213` 增至 `7650`。
- projection central algebraic residual 几乎不变。
- full-linear update 后 exact residual weighted mean 从 `0.14577` 降到 `0.13917`。

结论：边缘守恒约束有小幅帮助，但不是根治。

## 3. time_weight=0 测试

将辅助场代表元里的时间平滑权重设为 `0`。

结果：

- projection system residual 与默认几乎相同。
- full-linear update 后 exact residual weighted mean 升到 `0.16319`。

结论：`tau=-3.5` 问题不是简单的时间方向代表元偏置。

## 4. 扩展辅助张量基底 matter6_normal

新增 `--atom-family`：

- `matter4`：\((g, uu, rr, ur)\)。
- `matter6_normal`：\((uu, rr, nn, ur, un, rn)\)，其中 \(n_\mu\) 是 2+1 维中由 \(u,r,\tilde g\) 构造的法向协向量。

相关脚本已支持保存和自动读取 `atom_family/atom_names`：

- `kg_examples/fit_gbcd_principal_constraint_projection_sparse.py`
- `kg_examples/solve_gbcd_full_linear_metric_update_sparse.py`

`matter6_normal` 结果：

- projection system residual 从 `5.39e-3` 降到 `3.60e-4`，说明旧 4 基底确实在投影层缺少张量方向。
- 但 full-linear update 后 exact residual weighted mean 为 `0.15269`，不如 `matter4` 默认。
- `matter6 + force edge` 为 `0.15187`，仍不如 `matter4 + force edge` 的 `0.13917`。

结论：补足 2+1 维局域张量基底能改善投影系统，但不能解决非线性几何闭合；主瓶颈不是单纯少了 `nn/un/rn`。

## 5. center-plus 初值投影检验

对 `kg_examples/solve_gbcd_full_linear_metric_update_sparse.py` 增加：

`--metric-variable-slices plus|center_plus|minus_center_plus`

含义：

- `plus`：旧模式，只修正未来切片 \(g_+\)。
- `center_plus`：同时修正当前切片 \(g_0\) 和未来切片 \(g_+\)。
- `minus_center_plus`：预留，允许三切片一起修正。

动机：

如果 A-derived 当前几何 \(g_0\) 本身不是 D 支一致初值，只修正 \(g_+\) 可能过度受限；`center_plus` 是初值投影求解器的第一步。

`matter4 + force edge + center_plus + ridge=1e-6`：

- linear residual weighted mean 从 `0.13917` 进一步压到 `0.08504`。
- 但 exact nonlinear residual weighted mean 暴涨到 `0.92042`。
- \(\delta g_0\) 的 p95 相对量级约 `1.25`，说明线性步过大。

线搜索：

- 新增 `kg_examples/scan_gbcd_metric_update_alpha.py`。
- 对 \(\alpha\delta g\) 扫描后，最佳仍为 \(\alpha=1\)，但残差仍 `0.92042`；小 \(\alpha\) 不能接近 plus-only 的 `0.13917`。

强正则化 `ridge=1e-3`：

- linear residual weighted mean 进一步降到 `0.05835`。
- exact nonlinear residual 仍为 `0.90481`。

结论：

- `center_plus` 说明“当前切片自由度”在线性层面确实能降低残差。
- 但当前实现只改几何，没有同步重算/约束 \(T_{\mu\nu}\)、\(\mathcal C_{\mu\nu}\) 与 \(u,r,\rho\) 对 \(g_0\) 的依赖，因此非线性回代不自洽。
- 这不是物理模型失败，而是初值投影器还不是完整 D 初值求解器。

## 当前判断

1. `tau=-3.5` 的主问题不是低分辨率、不是简单阻尼、不是单一张量分量爆炸。
2. 旧 `plus-only` 几何更新目前仍是最稳可用更新，最佳为 `matter4 + force edge`，exact residual weighted mean `0.13917`。
3. 真正的下一步应把初值投影从“只解几何”升级为“联合解几何 + 辅助应力 + 物质源一致性”。
4. 如果允许 \(g_0\) 变化，就必须同步处理：
   - \(\tilde T_{\mu\nu}[\tilde\rho,u,\tilde g]\)；
   - \(\mathcal C_{\mu\nu}[\lambda_I,u,r,\tilde g]\)；
   - 质量壳、连续性、守恒约束；
   - \(\rho^\tilde\leftrightarrow\rho\) 的拉回关系。

## 关键输出

- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_iter3000_coeff_cg6000/residual_sources/gbcd_residual_source_diagnostics.png`
- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_forceedge_coeff_cg6000/summary.json`
- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_matter6_coeff_cg6000/summary.json`
- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_forceedge_centerplus_coeff_cg6000/summary.json`
- `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_forceedge_centerplus_coeff_cg6000/alpha_scan/gbcd_metric_update_alpha_scan.png`
