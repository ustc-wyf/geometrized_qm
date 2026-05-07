# D支 \(t_{\rm old}=8\) 时间步复查与 support 路线调整

日期：2026-05-03

## 用户问题

用户要求：

1. 在确定干涉窗口的 `dt` 后，说明前面“优先研究 support 边缘表达：平滑权重、tapered support 或 body-fitted 低密度边界”的计划是否还要继续；
2. 明确 support 边缘表达是否改变物理，还是只是数值方法；
3. 参数确认不能只看 \(t=0\)，还要看干涉发生时由 A 支确定的 D 支初态，例如 \(t_{\rm old}=8\) 附近。

## 本轮复查设置

共同物理设置：

- 1550nm optical massive-KG 模式；
- \(m=\omega_0/10\)；
- KG 内积归一化；
- \(\ell/l_P=10^{60}\)；
- `active_dilation=2`；
- `trusted_erosion=1`；
- `stop_mask=trusted`；
- 初态由 A 支自由 KG 在 \(t_{\rm old}=8\) 的 \((\rho_A,S_A)\) 经 D 支混合变换生成。

物理时间：

- \(t_{\rm old}=8\) 对应 `initial_time_fs = 39.4977507 fs`；
- 旧比例相遇时间约 \(t_{\rm old}\simeq 8.53\)，即该点处在干涉/重叠窗口前沿。

对照组：

- `n96 dt`：`dt_old=2.5e-7, steps=200`；
- `n96 dt/2`：`dt_old=1.25e-7, steps=400`；
- `n128 dt`：`dt_old=2.5e-7, steps=200`。

三组都对应同一局部物理窗口：

- `t_completed_fs = 2.4686094e-4 fs`。

## 输出

- 图：`visualizations/d_physical_1550nm_kg_told8_dt_resolution_decision/told8_dt_resolution_decision.png`
- 数据：`visualizations/d_physical_1550nm_kg_told8_dt_resolution_decision/told8_dt_resolution_decision_summary.json`

## 数值结果

`n96 dt`：

- `measure_rel_l1_trusted = 1.098689981e-3`
- `n_cons_rel_l1_trusted = 2.360176039e-4`
- `disc_min_trusted = 4.198042102e-7`
- `measure_rel_l1_support = 1.126289735e-3`
- `disc_min_support = 1.684543303e-7`

`n96 dt/2`：

- `measure_rel_l1_trusted = 1.098689979e-3`
- `n_cons_rel_l1_trusted = 2.360176039e-4`
- `disc_min_trusted = 4.198042101e-7`
- `measure_rel_l1_support = 1.126299837e-3`
- `disc_min_support = 1.684543304e-7`

`n128 dt`：

- `measure_rel_l1_trusted = 1.561562759e-3`
- `n_cons_rel_l1_trusted = 3.186200611e-4`
- `disc_min_trusted = 6.063554930e-8`
- `measure_rel_l1_support = 1.570867391e-3`
- `disc_min_support = 6.063554930e-8`

## 判定

\[
\boxed{
\text{\(t_{\rm old}=8\) 干涉窗口的问题不是时间步太大，也不是普通加密到 \(n=128\) 就能解决的空间分辨率问题。}
}
\]

理由：

- 在 `n=96` 下把时间步减半，trusted/support 的测度误差、守恒误差和质量壳判别式几乎逐项不变；
- `n=128` 不但没有改善 trusted measure，判别式裕度还更小；
- support 与 trusted 在 \(t_{\rm old}=8\) 已经同阶，说明这不是早期 `n64` 那种纯 support 外圈爆点。

因此当前默认时间步仍可保留：

- `dt_old=2.5e-7`

但这个结论只针对当前 `inert_tridomain` 短窗口诊断，不等于完整 D 支几何动力学已经稳定。

## support 边缘表达计划的修正

原计划：

- 优先研究 support 边缘表达：平滑权重、tapered support 或 body-fitted 低密度边界。

修正后：

- 该方向仍需保留，但降级为数值边界卫生检查；
- 它的任务是确认低密度 support 边缘没有人为制造伪故障；
- 它不能作为主物理闭合，也不能通过删除物质、修改源项或改变场方程来稳定结果；
- 若只改变积分权重、边界表示或贴体网格采样，而连续极限中的方程和源项不变，则属于合法数值方法；
- 若改变 \(\tilde\rho\)、\(\tilde T_{\mu\nu}\)、质量壳方程、D 支作用量或边界条件的物理内容，则属于改模型，不能采用。

当前新的主线优先级：

1. 诊断 \(t_{\rm old}=8\) 干涉窗口中 `inert_tridomain` 闭合的失败来源；
2. 检查饱和区惯性 metric 代表、plateau algebraic 残差、质量壳判别式裕度和 D-vs-A 源项错配；
3. 把几何闭合从 `inert_tridomain` 升级到允许 metric/scalaron/interface 自洽调整的版本；
4. support 平滑或 body-fitted 低密度边界作为并行的数值收敛控制，而不是主物理机制。

## 当前结论

这轮结果支持如下判断：

\[
\boxed{
\text{早期 \(t=0\) 的 support 边缘问题确实存在；但干涉窗口的主要瓶颈已经转向几何闭合，而不是 support 边缘或 dt。}
}
\]

因此，后续不应继续盲目减小 `dt`，也不应单纯期待 support smoothing 修复干涉区；应把 `t_old=8` 作为新的强测试点，直接研究 D 支几何闭合如何在强干涉量子势结构下自洽演化。
