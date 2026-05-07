# 092 固定实验室坐标严格初值与冻结 `\tilde g` 物质原型

## 目标

在用户明确接受：

- `\tilde g` 可在实验室坐标下出现时间反向等现象；
- 不再强求 `\tilde g` 适合 `ADM/BSSN` 切片；
- 重点是变换 `(g,\rho,S)\mapsto \tilde g` 自洽；

之后，重新回到**固定实验室坐标** `(t,x,z)`，直接检查：

1. 严格由 `g_0=\eta,\rho_0,S_0` 生成的初始 `\tilde g_0` 是否自洽；
2. 在此基础上，物质守恒流写法是否可以直接数值推进。

## 1. 严格初值生成

新增文件：

- [mixed_tilde_initial_data.py](/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/mixed_tilde_initial_data.py)

其中新增了：

- `localized_direct_tilde_coordinate_snapshot(...)`
- `localized_direct_tilde_coordinate_initial(...)`

它们直接从局域双高斯波包的正频波函数导出：

\[
\rho_0,\quad S_0,\quad S_{t,0},\quad S_{x,0},\quad S_{z,0},
\]

以及

\[
\tilde g^{ab}_0,\quad \tilde g_{ab,0},\quad \partial_t \tilde g_{ab}|_{t=0}.
\]

这里不再经过 `ADM` 的 `lapse/shift` 提取。

## 2. 自洽性

在主支撑区

\[
\rho > 10^{-3}\rho_{\max}
\]

内，直接验证

\[
\tilde g^{ab}u_a u_b = m^2,
\qquad
u_a=(S_t,S_x,S_z)
\]

成立。

数值结果：

\[
\max_{\text{support}}
\left|
\tilde g^{ab}u_a u_b - m^2
\right|
\approx 2.08\times 10^{-9}.
\]

因此，这组局域双高斯波包初值在当前变换下是自洽的。

## 3. 型别与实验室时间

对 `256×256` 网格、默认局域双高斯波包参数：

- 主支撑区网格点数：
  \[
  3246
  \]
- 在其中，`\tilde g` 的 `t-x-z` 三维逆度规块保持
  \[
  (1+,2-)
  \]
  型别的点数：
  \[
  2414
  \]
- 其余点数：
  \[
  832
  \]

另外，即使在核心区

\[
\rho > 10^{-1}\rho_{\max},
\]

也有大量点满足

\[
\tilde g^{tt}\le 0.
\]

因此：

- “实验室时间 `t` 不是 `\tilde g` 的好时间函数”仍然成立；
- 但这不再被解释成理论否定，只被视为旧 `ADM/BSSN` 数值方法不适配。

## 4. 固定实验室坐标下的物质原型

新增文件：

- [coordinate_matter_evolution.py](/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/coordinate_matter_evolution.py)

其中新增了基于协变量 `u_i` 的版本：

- `solve_covector_time_component(...)`
- `current_components_from_rho_u(...)`
- `conservative_density_from_rho_u(...)`
- `recover_rho_from_conservative_density_and_u(...)`
- `coordinate_matter_rhs_covector(...)`

这里不再把 `S` 的 `unwrap(arg ψ)` 标量重构当成数值主变量，而是直接推进：

\[
N_{\mathrm{cons}},\quad u_x,\quad u_z.
\]

## 5. 冻结初始 `\tilde g_0` 的最小原型

新增脚本：

- [simulate_b_coordinate_frozen_metric.py](/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_b_coordinate_frozen_metric.py)

先把几何固定为严格生成的 `\tilde g_0`，只推进物质守恒流。

### 5.1 不做尾部处理时

直接在整个网格上反演 `\rho = N_{\mathrm{cons}} / (\sqrt{|\tilde g|} J^t/\rho)` 会在稀薄尾部立刻炸掉。

原因不是主支撑区的判别式，而是低密度尾部中 `J^t` 逼近零，导致密度反演病态。

### 5.2 把稀薄尾部按“零测度密度真空区”处理后

按用户先前物理解读，在主支撑区外不再强求速度定义；数值上对一个膨胀过的主支撑掩膜外设

\[
N_{\mathrm{cons}} = 0.
\]

重跑短时间窗：

\[
t=0.02,\qquad 80\ \text{步},\qquad dt=2.5\times 10^{-4}.
\]

输出文件：

- [summary.json](/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/b_coordinate_frozen_covector_masked_t002/summary.json)

结果：

- 主支撑区判别式始终保持正：
  \[
  \Delta_{\text{disc,min,support}}
  \approx 4.03
  \]
- `\rho` 有限，且未再出现前一版的 `NaN/Inf`
- 与初始 `\rho` 的相对 `L^1` 偏离：
  \[
  \approx 1.37\times 10^{-2}
  \]

因此，固定实验室坐标下的协变量守恒流写法，在“主支撑区 + 尾部真空化”处理下是可运行的。

## 6. 当前阶段性结论

1. 严格由 `(g_0,\rho_0,S_0)` 生成 `\tilde g_0` 是可行的。
2. 先前的主要数值假象来自：
   - `unwrap(arg ψ)` 相位梯度重构；
   - 稀薄尾部对 `J^t` 近零点的病态反演。
3. 在主支撑区内，当前局域双高斯波包初值与变换是自洽的。
4. 固定实验室坐标、以
   \[
   (N_{\mathrm{cons}},u_x,u_z)
   \]
   为基本物质量，再把尾部按真空区处理，是后续直接分量几何演化的合理起点。

## 7. 下一步

下一步不再回到 `ADM/BSSN`，而是：

1. 在固定实验室坐标下，给 `B` 支写直接分量几何演化原型；
2. 几何主变量直接取 `\tilde g_{ab}` 与其时间导数；
3. 物质继续使用当前已经稳定下来的
   \[
   (N_{\mathrm{cons}},u_x,u_z)
   \]
   写法；
4. 尾部继续按“零测度密度真空区”处理。
