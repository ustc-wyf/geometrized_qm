# D支 support 边缘表达后处理与固定参数判定

日期：2026-05-03

## 用户排队指令

用户要求在 \(dt\) 确定后，根据最近两轮分析判断：

1. 分辨率是否固定为 `n=96`；
2. 时间步是否固定；
3. 若没有提升，则不做动态参数；
4. 若固定参数足够，就继续推进 “support 边缘表达”。

## 参数判定

当前对照结果：

- `t_old=8,n96,dt_old=2.5e-7`
- `t_old=8,n96,dt_old=1.25e-7`
- `t_old=8,n128,dt_old=2.5e-7`

在 binary support、trusted core、support edge band、tapered support 四种口径下，`n96 dt` 与 `n96 dt/2` 几乎完全一致。

因此：

\[
\boxed{
n=96,\quad dt_{\rm old}=2.5\times 10^{-7}
}
\]

可作为当前 `inert_tridomain` 诊断阶段的固定工作参数。

不采用动态参数的理由：

- `dt` 减半无可见改善；
- `n=128` 无改善，甚至部分指标更差；
- 动态参数会增加变量，让后续 support/bulk closure 分析更难归因；
- 当前瓶颈不表现为 CFL 时间步误差，也不表现为普通空间分辨率误差。

## support 边缘表达实现

新增后处理脚本：

- `kg_examples/postprocess_d_support_edge.py`

它只读取已有完整场数据：

- `summary.json`
- `fields_final.npz`

并计算四类诊断口径：

1. binary support：原始二值 support；
2. trusted core：support 腐蚀一格后的核心；
3. support edge band：`support - trusted`；
4. tapered support：在 support 阈值以下按对数尺度平滑衰减的连续权重。

重要限定：

- 这是后处理；
- 不重新演化；
- 不改 \(\tilde\rho\)；
- 不改 \(\tilde T_{\mu\nu}\)；
- 不改质量壳方程；
- 不改作用量；
- 不改物理边界条件。

所以它是数值诊断方法，不是新物理。

## 输出

主输出目录：

- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/`

主图：

- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/support_edge_taper_overview.png`

逐 case 图：

- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/d_physical_1550nm_kg_told8_n96_dt2p5em7_steps200_support_edge_taper.png`
- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/d_physical_1550nm_kg_told8_n96_dt1p25em7_steps400_support_edge_taper.png`
- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/d_physical_1550nm_kg_told8_n128_dt2p5em7_steps200_support_edge_taper.png`

数据：

- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess/support_edge_taper_summary.json`

额外 taper 宽度检查：

- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess_taper0p5/`
- `visualizations/d_physical_1550nm_kg_told8_support_edge_postprocess_taper2p0/`

## 主要数值

### taper 宽度 0.5 decade

`n96 dt`：

- support `1.126290e-3`
- trusted `1.098690e-3`
- edge `1.346598e-3`
- taper `1.380657e-3`

`n96 dt/2`：

- support `1.126300e-3`
- trusted `1.098690e-3`
- edge `1.346689e-3`
- taper `1.380667e-3`

`n128 dt`：

- support `1.570867e-3`
- trusted `1.561563e-3`
- edge `1.749782e-3`
- taper `1.598716e-3`

### taper 宽度 1.0 decade

`n96 dt`：

- support `1.126290e-3`
- trusted `1.098690e-3`
- edge `1.346598e-3`
- taper `2.059287e-3`

`n96 dt/2`：

- support `1.126300e-3`
- trusted `1.098690e-3`
- edge `1.346689e-3`
- taper `2.059297e-3`

`n128 dt`：

- support `1.570867e-3`
- trusted `1.561563e-3`
- edge `1.749782e-3`
- taper `1.621821e-3`

### taper 宽度 2.0 decades

`n96 dt`：

- support `1.126290e-3`
- trusted `1.098690e-3`
- edge `1.346598e-3`
- taper `2.355253e-3`

`n96 dt/2`：

- support `1.126300e-3`
- trusted `1.098690e-3`
- edge `1.346689e-3`
- taper `2.355260e-3`

`n128 dt`：

- support `1.570867e-3`
- trusted `1.561563e-3`
- edge `1.749782e-3`
- taper `1.652621e-3`

## 解读

1. `dt` 固定：
   - 所有 support 表达口径下，`dt` 减半都几乎不改变结果。

2. `n=96` 固定：
   - `n128` 没有改善 measure relative L1；
   - 因此当前不应把 `n128` 作为默认动态升级。

3. support 边缘不是唯一主因：
   - support edge band 的误差确实略高于 trusted core；
   - 但 trusted core 本身也在 \(10^{-3}\) 量级；
   - 这说明边缘表达会影响诊断数值，但不能解释全部干涉窗口偏离。

4. tapered support 的作用：
   - taper 越宽，越多低密度外圈进入权重；
   - 误差会升高，但 `n96 dt` 与 `n96 dt/2` 仍保持一致；
   - 这说明 support 边缘表达不是时间步问题。

## 当前结论

\[
\boxed{
\text{当前固定 }n=96,\ dt_{\rm old}=2.5\times10^{-7}\text{；support 边缘表达继续作为后处理/诊断工具推进，不引入动态参数。}
}
\]

下一步应在固定参数下继续：

- 用 support/taper/body-fitted 口径清理低密度边界诊断；
- 但主物理瓶颈仍要回到几何闭合，尤其是 saturated bulk metric 代表与干涉区自洽调整。
