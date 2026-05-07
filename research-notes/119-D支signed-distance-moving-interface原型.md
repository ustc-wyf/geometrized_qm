## D支 signed-distance moving-interface 原型

日期：2026-05-02

### 目的

上一版 `prototype_d_moving_interface_reduced.py` 直接推进

\[
raw_y=\ell^2\tilde R
\]

并暴露出强 stiffness：界面速度 \(v_n\) 本身是 \(O(1)\)，但由于
`raw_y` 在过渡层附近梯度很大，

\[
F_t=-v_n|\nabla raw_y|
\]

可变得极大。因此本轮把界面变量从原始 `raw_y` 改为 signed-distance level set。

### 新脚本

`kg_examples/prototype_d_signed_distance_moving_interface.py`

核心变量：

- `d_plus`：到 `raw_y=+1` 分支的 signed distance，初始化符号为 `sign(raw_y-1)`；
- `d_minus`：到 `raw_y=-1` 分支的 signed distance，初始化符号为 `sign(raw_y+1)`；
- 两个界面分别为 `d_plus=0` 与 `d_minus=0`。

每一步：

1. 从 `d_plus=0`、`d_minus=0` 抽取连续界面线段；
2. 在界面上使用 D 支 leading trace/scalaron jump law 反求法向速度；
3. 用 cloud-in-cell 把界面速度延拓到近界面网格；
4. 按 Hamilton-Jacobi 形式推进 signed distance：

   \[
   d_t+v_n|\nabla d|=0;
   \]

5. 从新的零水平集重新初始化 signed distance。

### 与物理方程的关系

本轮仍然只使用 D 支 trace/scalaron jump law 给界面速度，没有引入 clipping/damping。

但它还不是完整全动力学，因为：

- `rho` 冻结；
- `metric_inv` 冻结；
- `stress_trace` 冻结；
- `alpha=|\nabla raw_y_{\rm ref}|` 冻结，用来把 `q=raw_y` 的 jump law 转成 signed-distance 表示；
- 未解段仍只是标记为 dotted unresolved，没有用 tensor matching 修复。

因此本轮结果只能说明 signed-distance 表示比直接推进原始 `raw_y` 更适合作为 moving-interface 数值坐标，不能作为最终物理预测。

### 主要运行

较粗时间步：

```bash
python3 kg_examples/prototype_d_signed_distance_moving_interface.py \
  --output visualizations/d_signed_distance_moving_interface_ell30_t16_96 \
  --ell 30 --resolution 96 --time 16 --levels 1 \
  --half-width-y 1 --samples 61 \
  --dt 0.05 --steps 8 --extension-iterations 10 --render-every 2
```

小时间步对照：

```bash
python3 kg_examples/prototype_d_signed_distance_moving_interface.py \
  --output visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001 \
  --ell 30 --resolution 96 --time 16 --levels 1 \
  --half-width-y 1 --samples 61 \
  --dt 0.01 --steps 20 --extension-iterations 10 --render-every 5
```

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/summary.json`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/signed_distance_interface_overlay_ell30_n96_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/frames/signed_distance_interface_step_020.png`

### 小时间步结果

`ell=30,t=16,n=96,dt=0.01,steps=20`：

- step 0：
  - 界面段数 `864`
  - 可解段 `515`
  - 可解比例 `0.596`
  - 总界面长度 `307.053`
  - 修正速度 p95 `1.129`
- step 20：
  - 界面段数 `712`
  - 可解段 `445`
  - 可解比例 `0.625`
  - 总界面长度 `214.710`
  - 修正速度 p95 `3.105`

解释：

- signed-distance 重初始化没有导致界面整体崩溃；
- 界面总长度下降，说明部分小环/局部分支在 jump-law 推进下收缩或消失；
- 可解比例维持在 `0.60-0.63`，说明 trace-only jump law 的未解段不是时间步问题；
- 局部最大速度仍有尖峰，说明下一步必须补充 tensor jump/matching 或更完整的界面形状自由度，不能靠数值平滑掩盖。

### 当前结论

1. `raw_y` 是物理层变量，但不适合直接作为显式演化 level-set 数值变量。
2. signed-distance 是合规的数值坐标重置，因为它只改变界面表示，不改变 D 支 jump law。
3. 当前 moving-interface 路线已经从“画界面”推进到“界面可运动、可重初始化、可发生拓扑变化”的原型。
4. 真正全动力学的下一个瓶颈不是 signed-distance，而是完整 tensor matching 与 bulk 自洽演化。

### 下一步

- 在 `d_plus/d_minus` 界面上加入完整 tensor jump/matching 的最小闭合条件；
- 把 `rho,S` 的连续性/质量壳演化接入，而不是冻结参考物质；
- 把 `metric_inv/stress_trace/alpha` 从冻结参考量改成由当前 \(\tilde g,\rho,S\) 更新；
- 用局部贴体坐标或高阶界面 quadrature 改善 unresolved 段，而不是 clipping/damping。
