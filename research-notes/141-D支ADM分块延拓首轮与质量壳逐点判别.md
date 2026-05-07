# D 支 ADM 分块延拓首轮与质量壳逐点判别

日期：2026-05-03

## 问题

用户要求继续尝试更温和的饱和区 metric 代表选择，并确认“质量壳”诊断到底是逐点检查还是对支撑区做平均。

## 质量壳判别口径

当前代码中的质量壳判别是逐点的，不是平均的。

具体做法是，在每个网格点上由质量壳方程

\[
\tilde g^{\mu\nu}u_\mu u_\nu=m^2
\]

对 \(u_t\) 形成一个二次方程，并计算判别式

\[
D=b^2-a c .
\]

这在代码中作为 `current["discriminant"]` 存成逐点数组。

停止条件使用的是指定 mask 上的逐点最小值：

\[
\min_{\rm mask}D .
\]

因此只要支撑区或 trusted 区中有一个被纳入检查的网格点出现明显负判别式，就会触发停止；这不是支撑内平均后再判断。

## ADM 分块延拓

为避免直接对 \(g_{\mu\nu}\) 做 harmonic extension 破坏 Lorentz 性，已在

- `kg_examples/simulate_d_tridomain_full_dynamics.py`

中加入

- `--metric-extension-variable-mode covariant`
- `--metric-extension-variable-mode adm`

`adm` 模式尝试将协变度规分解为：

- log lapse；
- shift；
- 空间度规的 Cholesky 变量；

然后对这些变量做受限 harmonic extension，再重构 \(g_{\mu\nu}\)。

## 测试设置

- 物理模式：1550 nm，KG 归一化；
- \(m=\omega_0/10\)；
- \(\ell/l_P=10^{60}\)；
- `t_old=8`；
- `n=96`；
- `dt_old=2.5e-7`；
- `steps=20`；
- `geometry_closure=restricted_extension`；
- `metric_extension_variable_mode=adm`；
- flat boundary anchors 打开。

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_adm_flat_anchor_n96_steps20/summary.json`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_adm_flat_anchor_n96_steps20/d_tridomain_dynamics_restricted_extension_ell8.19075e+31_n96_t3.75048e-05.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_restricted_extension_adm_flat_anchor_n96_steps20/support_edge_tapered_diagnostics.png`

## 结果

关键诊断：

- `metric_extension_variable_mode = adm`
- `metric_extension_flat_anchor_count = 97`
- `metric_extension_solve_count = 5`
- `metric_extension_solve_fraction = 0.012165450121654502`
- `metric_extension_anchor_count = 97`
- `metric_extension_anchored_component_count = 1`
- `metric_extension_orphan_component_count = 0`
- `metric_extension_lorentz_fallback_count = 0`
- `metric_extension_adm_valid_fraction = 0.25790754257907544`
- `metric_extension_change_rel_p95 = 0.9851585826090669`
- `metric_extension_roughness_before_p95 = 368.9097786294534`
- `metric_extension_roughness_after_p95 = 1.6304060060012092e+53`
- `metric_extension_admissible_blend = 0`
- `metric_extension_admissible_attempts = 10`

最终物质诊断：

- `measure_rel_l1_trusted = 0.00011488058838534157`
- `disc_min_trusted = 4.0995731809267966e-07`
- `measure_rel_l1_support = 0.00011397749603909711`
- `n_cons_rel_l1_support = 2.339903038536604e-05`
- `saturated_fraction = 1.0`

## 解释

ADM 分块有两个正面结果：

- 没有触发 Lorentz fallback；
- flat boundary anchors 仍然使饱和连通块不再是 orphan。

但它仍未闭合饱和区 metric：

\[
\texttt{metric\_extension\_admissible\_blend}=0 .
\]

主要问题是当前 lab-time 切片下只有约 \(25.8\%\) active 点能被分解成实 ADM 变量；并且在可延拓区域上重构后的空间梯度粗糙度反而变得极大。

因此当前失败不是“没有 metric”，而是：

\[
\boxed{
\text{当前 full ADM 变量延拓在固定实验室时间切片上太强或太病态，不能作为饱和区代表选择。}
}
\]

## 下一步

更合理的下一步不是继续直接平滑全度规，而是收缩 ADM 更新自由度：

1. 只更新 lapse，固定 shift 与空间度规；
2. 或更新 lapse + shift，固定空间度规；
3. 或将 direct tensor interface / jump matching 作为额外物理锚点接入；
4. 同时检查 ADM invalid 点是否对应时间反向、signature 分支或当前实验室切片不适配。

