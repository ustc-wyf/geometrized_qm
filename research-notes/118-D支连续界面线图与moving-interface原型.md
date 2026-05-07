## D支连续界面线图与 moving-interface 原型

日期：2026-05-02

### 背景

用户指出此前图中界面像随机小点，不利于判断物理成因，并要求继续推进完整模拟器。

本轮完成两件事：

1. 将界面可视化从“线段中点散点图”改为“连续线段图”。
2. 新增一个 moving-interface reduced prototype，用 jump-law 速度推进界面水平集，并每步重新抽取界面。

### 连续线段图

修改：

- `kg_examples/diagnose_d_trace_interface_weak_form.py`
  - `find_level_segments` 现在保留每条 marching-square 界面线段的端点：
    `x0,z0,x1,z1`。
- `kg_examples/diagnose_d_interface_speed_law.py`
  - 速度律图从散点改为 `LineCollection` 连续线段。
  - 未解点画为灰/黑色虚线，不再隐藏。

运行：

```bash
python3 kg_examples/diagnose_d_interface_speed_law.py \
  --output visualizations/d_interface_speed_law_lines_ell30_t16_128_fix \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 --interface-dt 0.0002
```

图：

`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_interface_speed_law_lines_ell30_t16_128_fix/d_interface_speed_law_ell30_n128_t16.png`

结果与上一轮一致：

- 总界面线段 `1361`
- 可解 `849`
- 未解 `512`
- 可解比例 `0.624`
- 修正速度 p95 `1.437`
- 修正法向范数 p95 `5.663e-3`

视觉结论：

- 界面不是随机点，而是连续曲线片段；
- 主要结构包括：上方竖向条纹、左右弧形叶片、中心缠结区、低密度边缘小弧；
- 这些结构大体贴着干涉后形成的量子曲率过渡带，而不是均匀随机撒点。

### moving-interface reduced prototype

新增脚本：

`kg_examples/prototype_d_moving_interface_reduced.py`

功能：

- 从参考切片的 `raw_y=ell^2 R_tilde` 抽取 `raw_y=±1` 界面；
- 用 leading trace/scalaron jump law 反求界面 \(F_t\)；
- 把可解线段的 \(F_t\) 用 cloud-in-cell 沉积到网格；
- 用邻近平均做 level-set velocity extension；
- 以

\[
raw_y^{n+1}=raw_y^n+\Delta t\,F_t
\]

做小步显式推进；
- 每一步重新从整张 `raw_y` 场抽取界面，因此形式上可以看到界面分裂、合并、新生、消失。

运行：

```bash
python3 kg_examples/prototype_d_moving_interface_reduced.py \
  --output visualizations/d_moving_interface_reduced_ell30_t16_96_fix \
  --ell 30 --resolution 96 --time 16 --levels 1 \
  --half-width-y 1 --samples 61 \
  --dt 0.0001 --steps 6 --extension-iterations 8 --render-every 1
```

输出：

`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_moving_interface_reduced_ell30_t16_96_fix/summary.json`

末帧：

`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_moving_interface_reduced_ell30_t16_96_fix/frames/moving_interface_step_006.png`

统计：

- step 0：`866` 段，可解 `555`，可解比例 `0.641`
- step 6：`958` 段，可解 `657`，可解比例 `0.686`
- 修正速度 p95 维持在 `~1.3-1.4`
- 修正法向范数 p95 维持在 `~5.66e-3`

### 重要限制

这个原型不是完整全动力学，也不是最终物理预测。它冻结了 `metric_inv`、`stress_trace`、`rho`，只推进 `raw_y` 水平集。

更重要的是，直接显式推进原始 `raw_y` 会迅速暴露 stiffness：因为 `raw_y` 在过渡层附近梯度极大，虽然界面速度 \(v_n\) 是 \(O(1)\)，但

\[
F_t=-v_n|\nabla raw_y|
\]

可以非常大。因此原始 `raw_y` 不是一个适合直接显式推进的 level-set 函数。

这不是物理上要削掉的发散，而是数值表示问题：下一版 moving-interface solver 应推进“界面位置”或“重初始化后的 signed-distance level set”，同时用 D 支 jump law 提供法向速度。

### 下一步

- 保留连续线段图作为默认界面展示方式；
- 完整模拟器下一步应改为：
  - 用 signed-distance / body-fitted interface 表示界面位置；
  - 每步由 D 支 jump/matching 条件给出速度；
  - 每步从全域水平集重新抽取界面，以允许生成/消失/合并/断裂；
  - 对未解段升级到 tensor jump/matching，而不是 clipping 或 damping；
  - 物质 \((\rho,S)\) 和 bulk 几何仍需接入自洽演化，当前原型尚未完成这一层。

### 线性 rho 图与断裂来源复查

新增脚本：

`kg_examples/plot_d_interface_linear_zooms.py`

运行：

```bash
python3 kg_examples/plot_d_interface_linear_zooms.py \
  --output visualizations/d_interface_linear_zooms_ell30_t16_192 \
  --ell 30 --resolution 192 --time 16 --levels 1 --rho-min-frac 1e-3
```

输出目录：

`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_interface_linear_zooms_ell30_t16_192/`

代表图：

- `d_interface_linear_global_ell30_n192_t16.png`
- `d_interface_linear_upper_transition_fan_ell30_n192_t16.png`
- `d_interface_linear_central_tangle_ell30_n192_t16.png`
- `d_interface_linear_left_leaf_ell30_n192_t16.png`
- `d_interface_linear_right_leaf_ell30_n192_t16.png`
- `d_interface_linear_low_density_edge_ell30_n192_t16.png`

图的三栏含义：

1. 线性 \(\rho\) + 全部 \(raw_y=\pm1\) 等值线，不加 \(\rho\) 支撑阈值；
2. 线性 \(\rho\) + 主支撑区保留界面，阈值为 \(\rho>10^{-3}\rho_{\max}\)；
3. \(raw_y\) 场 + 主支撑区保留界面。

关键统计：

- 全域不加支撑阈值的等值线段数：`47379`
- 主支撑区保留线段数：`2287`
- 被支撑阈值排除的线段数：`45092`
- 保留比例：`0.0483`
- 全部等值线段的中点 \(\rho\) 中位数：`2.27e-25`
- 主支撑区保留线段中点 \(\rho\) 中位数：`7.98e-3`

解释：

- 许多看似“断裂”的地方其实是低密度尾部被物理工作阈值排除了；
- \(raw_y=\pm1\) 是标量场水平集，不是物质边界，数学上可以是闭合曲线、局部弧段、分叉附近的小环，不必延伸到计算区域边界；
- 在中心和上方扇区，断裂还与多个水平集分支靠近、局部极值/鞍点和有限网格分辨率有关。

### signed-distance 预览

新增脚本：

`kg_examples/prototype_d_signed_distance_interface.py`

运行：

```bash
python3 kg_examples/prototype_d_signed_distance_interface.py \
  --output visualizations/d_signed_distance_interface_ell30_t16_128 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --rho-min-frac 1e-3 --near-distance 0.75
```

输出：

`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_interface_ell30_t16_128/d_signed_distance_interface_ell30_n128_t16.png`

定义：

- signed distance \(d(x,z)\) 在界面上为 `0`；
- \(|d|\) 是到最近界面线段的欧氏距离；
- 本轮约定 \(|raw_y|<1\) 一侧为负，\(|raw_y|\ge1\) 一侧为正。

目的：

- 它不是新物理场；
- 它是数值重初始化变量，用来在 moving-interface solver 中替代陡峭的原始 `raw_y`；
- 理想 signed distance 满足 \(|\nabla d|\simeq1\)。

首轮统计：

- 支撑区界面段数：`1361`
- 近界面支撑区网格点数：`1089`
- \(|\nabla d|\) 近界面 p95：`0.998`

结论：

signed-distance 表示可作为下一版 moving-interface solver 的数值变量；均值偏低来自多界面靠近、分支交汇和当前粗网格下最近距离函数不可光滑。
