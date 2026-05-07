# D支干涉时刻初态与 support 边缘表达说明

日期：2026-05-03

## 用户问题

用户指出两点：

1. `support` 边缘表达是否是新策略，它是否改变物理；
2. 参数确认不能只看 \(t=0\) 附近，还要看干涉发生时量子势变化，例如由 A 支在 \(t_{\rm old}=8\) 生成 D 支初态后，在该时刻附近的时间步长和分辨率是否仍平稳。

## support/trusted 的物理地位

当前脚本中：

- `trusted` 只影响诊断指标与停止条件；
- `support -> active` 定义实际演化窗口；
- `active` 外的 \(\tilde\rho,n_{\rm cons},u_t\) 被置零，metric tail 也被 sanitize。

因此：

- `trusted` 不是物理规则，也不改变物理方程；
- `support`/`active` 是数值计算区域与外边界处理；
- 它不能被解释为新的物理筛选机制；
- 如果未来用平滑权重/tapered support/body-fitted 边界，只要不改变方程源项，只是改变积分区域和低密度边界表示，就属于数值方法；
- 如果用它删除物质、修改源项或改变边界条件，则会改变物理，不应这样做。

本轮已把这个区分写入工作判断：后续必须报告 support/trusted 指标，但不能把 support 边缘表达当成新物理。

## 代码更新

修改：

- `kg_examples/mixed_tilde_initial_data.py`
- `kg_examples/simulate_d_tridomain_full_dynamics.py`

新增能力：

- `localized_direct_tilde_coordinate_initial(..., initial_time=...)`
- `simulate_d_tridomain_full_dynamics.py --initial-time`
- `simulate_d_tridomain_full_dynamics.py --initial-time-old-units`

含义：

- 先用 A 支自由 KG 精确谱演化到指定时刻 \(t_A\)；
- 再用该时刻的 \((\rho_A,S_A)\) 和混合变换构造 D 支初态；
- 随后只研究该初态附近的 D 支局部演化。

这解决了此前只能从 \(t=0\) 构造 D 初态的问题。

## \(t_{\rm old}=8\) 干涉窗口检查

物理换算：

- \(t_{\rm old}=8\) 对应 `initial_time_fs = 39.4977507 fs`；
- 旧比例下两个包相遇时间约 \(t_{\rm old}\sim 8.53\)，所以 \(t_{\rm old}=8\) 已经接近干涉/重叠区。

共同设置：

- 1550nm physical optical mode；
- KG norm；
- \(\ell/l_P=10^{60}\)；
- `dt_old=2.5e-7`；
- local evolution window `steps=200`，即 \(t_{\rm old}=5\times10^{-5}\)；
- `active_dilation=2`；
- `trusted_erosion=1`；
- `stop_mask=trusted`。

输出：

- `visualizations/d_physical_1550nm_kg_interference_time_check/interference_time_check_diagnostics.png`
- `visualizations/d_physical_1550nm_kg_interference_time_check/interference_time_check_summary.json`
- `visualizations/d_physical_1550nm_kg_interference_resolution_check/interference_resolution_diagnostics.png`
- `visualizations/d_physical_1550nm_kg_interference_resolution_check/interference_resolution_summary.json`

## 与 \(t_{\rm old}=0\) 的对比

`t_old=0, n=96` 半窗口：

- `measure_rel_l1_trusted = 1.448e-5`
- `n_cons_rel_l1_trusted = 2.012e-5`
- `disc_min_trusted = 3.190e-2`
- `measure_rel_l1_support = 1.648e-5`
- `disc_min_support = 8.080e-3`

`t_old=8, n=96` 半窗口：

- `measure_rel_l1_trusted = 1.099e-3`
- `n_cons_rel_l1_trusted = 2.360e-4`
- `disc_min_trusted = 4.198e-7`
- `measure_rel_l1_support = 1.126e-3`
- `disc_min_support = 1.685e-7`

判断：

- \(t_{\rm old}=8\) 明显更苛刻；
- 误差升高约两个数量级；
- 质量壳判别式裕度从 \(10^{-2}\) 降到 \(10^{-7}\)；
- support 和 trusted 仍同阶，说明这不是单纯 support 外圈问题。

## \(t_{\rm old}=8\) 的分辨率检查

`t_old=8, n=96`：

- `measure_rel_l1_trusted = 1.099e-3`
- `n_cons_rel_l1_trusted = 2.360e-4`
- `disc_min_trusted = 4.198e-7`
- `measure_rel_l1_support = 1.126e-3`
- `disc_min_support = 1.685e-7`

`t_old=8, n=128`：

- `measure_rel_l1_trusted = 1.562e-3`
- `n_cons_rel_l1_trusted = 3.186e-4`
- `disc_min_trusted = 6.064e-8`
- `measure_rel_l1_support = 1.571e-3`
- `disc_min_support = 6.064e-8`

判断：

- `n=128` 没有改善 \(t_{\rm old}=8\)；
- 判别式裕度反而更小；
- 因此干涉区困难不是单纯空间分辨率不足；
- 下一步应检查时间步或当前 `inert_tridomain` 几何闭合在干涉区的适用性。

## 当前结论

\[
\boxed{
\text{support 边缘表达是数值区域/边界表示问题，不是新物理；但干涉区稳定性必须单独校准，不能从 }t=0\text{ 外推。}
}
\]

更具体地：

- \(t=0\) 附近：`n=96, dt_old=2.5e-7` 已经相当稳定；
- \(t_{\rm old}=8\) 附近：同一参数可跑完短窗口，但判别式裕度很小，误差明显增大；
- `n=128` 没有修复干涉区问题；
- 因此下一步最小必要检查是 `t_old=8, n=96, dt_old=1.25e-7`，判断是否是时间步问题；
- 若时间步减半仍不改善，则问题更可能来自当前 `inert_tridomain` 闭合在干涉区的不足，或需要更真实的动态几何/边界匹配。

