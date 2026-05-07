# D 支 ADM lapse 小步上限扫描

日期：2026-05-03

## 目的

在 `branch-preserving signature + ADM lapse-only` 闭合下，扫描饱和区代表 metric 小步松弛上限：

\[
\epsilon_g=\texttt{metric\_extension\_max\_rel\_change}.
\]

目标是判断这个闭合是否存在稳定窗口，而不是只依赖单点 `1e-3`。

## 固定设置

- 1550 nm 物理标定；
- KG 归一化；
- \(m=\omega_0/10\)；
- \(\ell/l_P=10^{60}\)；
- `t_old=8`；
- `n=96`；
- `dt_old=2.5e-7`；
- `steps=20`；
- `geometry_closure=restricted_extension`；
- `metric_extension_variable_mode=adm`；
- `metric_extension_adm_update_fields=lapse`。

## 输出

汇总输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_cap_scan/summary.json`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_physical_1550nm_kg_told8_adm_lapse_cap_scan/adm_lapse_cap_scan.png`

各组输出：

- `visualizations/d_physical_1550nm_kg_told8_adm_lapse_branch_preserve_cap1e4_n96_steps20/`
- `visualizations/d_physical_1550nm_kg_told8_adm_lapse_branch_preserve_cap3e4_n96_steps20/`
- `visualizations/d_physical_1550nm_kg_told8_adm_lapse_branch_preserve_cap1e3_n96_steps20/`
- `visualizations/d_physical_1550nm_kg_told8_adm_lapse_branch_preserve_cap3e3_n96_steps20/`

## 结果

| cap | trusted measure rel L1 | support measure rel L1 | trusted \(D_{\min}\) | accepted blend | actual change p95 | interface solved fraction |
|---:|---:|---:|---:|---:|---:|---:|
| inert | \(1.1488\times10^{-4}\) | \(1.1398\times10^{-4}\) | \(4.0996\times10^{-7}\) | 0 | 0 | 0.3762 |
| \(10^{-4}\) | \(1.3096\times10^{-4}\) | \(1.2891\times10^{-4}\) | \(4.0996\times10^{-7}\) | 1 | \(8.64\times10^{-5}\) | 0.4381 |
| \(3\times10^{-4}\) | \(1.6315\times10^{-4}\) | \(1.5864\times10^{-4}\) | \(4.0996\times10^{-7}\) | 1 | \(2.59\times10^{-4}\) | 0.4429 |
| \(10^{-3}\) | \(2.7536\times10^{-4}\) | \(2.6229\times10^{-4}\) | \(4.0996\times10^{-7}\) | 1 | \(8.64\times10^{-4}\) | 0.4498 |
| \(3\times10^{-3}\) | \(5.9205\times10^{-4}\) | \(5.5481\times10^{-4}\) | \(4.0996\times10^{-7}\) | 1 | \(2.59\times10^{-3}\) | 0.4498 |

## 判断

扫描显示：

- 所有 cap 组都被 admissible 检查接受；
- 质量壳逐点最小判别式几乎不变；
- `n_cons` 误差几乎不变；
- transformed measure 偏离随 cap 单调增大；
- interface solved fraction 在小 cap 下已经从 0.376 提升到约 0.44，之后收益逐渐饱和。

当前最好的折中是：

\[
\boxed{
\epsilon_g=3\times10^{-4}
}
\]

理由：

- 相比 inert，它引入了非平凡 metric 松弛；
- trusted measure 仍只有 \(1.63\times10^{-4}\)，接近 inert 的 \(1.15\times10^{-4}\)；
- interface solved fraction 从 0.376 提升到 0.443；
- 比 \(10^{-3}\) 更保守，避免过强边界锚点牵引。

## 下一步

建议把当前默认 working closure 更新为：

```text
--geometry-closure restricted_extension
--metric-extension-variable-mode adm
--metric-extension-adm-update-fields lapse
--metric-extension-max-rel-change 3e-4
```

下一步继续：

1. 更长步数测试；
2. 画 ADM invalid 区；
3. 与 direct tensor interface matching 联合；
4. 在 `t_old=0` 与其他干涉时刻复查同一 cap 是否仍稳定。

