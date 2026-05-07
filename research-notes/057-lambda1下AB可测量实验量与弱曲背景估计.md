# 057 `lambda=1` 下 A/B 可测量实验量与弱曲背景估计

## 目标

把当前 `lambda_grav = 1` 的 A/B 差异，从单纯的 `L1` 误差，
进一步压成更接近实验读出的量：

```text
1. 条纹间距
2. 等效相位移
3. 可见度
4. 亮纹强度差
```

并在两类弱背景上做第一轮定量估计：

```text
1. 地球表面附近的 Schwarzschild 弱场
2. 已观测到的纳赫兹引力波背景
```

## 口径说明

这里仍然沿用当前最可信的受控 shared-baseline scaffold：

```text
A: 可信 PDE 骨架动态演化
B: 在同一组 A 的 rho,S 历史上推进 Einstein 几何并重构 rho_B
```

因此这里得到的“相位差”，准确说是：

```text
从最终条纹图样反推出来的等效 fringe-pattern shift
```

而不是一个“B 自己独立传播出来的额外相位场”。

## 脚本

主脚本：

```text
kg_examples/analyze_ab_measurable_observables.py
```

输出目录：

```text
visualizations/ab_measurable_observables_lambda1/
```

默认测量假设：

```text
示例物理条纹间距 d_phys = 1 micron
示例读出尺度 L_readout = 1 m
GW characteristic strain = 2.4e-15
```

后两项只是为了把极小效应换算成直观数量级；
平直背景下的核心结论并不依赖这些例子。

## 1. `lambda=1` 平直基线下的 A/B 实验量

取最终时刻 `t/T = 1`，
在中心横截线 `z = 0` 上抽取条纹。

### 1.1 条纹间距

得到：

```text
dominant fringe spacing (FFT)  ~ 0.44456
peak-to-peak spacing for A     ~ 0.44450
peak-to-peak spacing for B     ~ 0.44450
relative spacing shift         ~ 0
```

也就是说：

```text
在当前 `lambda=1` scaffold 下，
A 与 B 的条纹间距没有可分辨的变化。
```

### 1.2 等效相位移

从同一中心截线的主频模式提取：

```text
effective phase shift ~ -5.58e-12 rad
equivalent fringe shift / spacing ~ -8.88e-13
```

若拿一个示例物理条纹间距

```text
d_phys = 1 micron
```

来换算，则对应的等效条纹平移只有：

```text
Delta x_eff ~ -8.88e-19 m
```

这基本可以视为：

```text
没有可测的条纹整体平移。
```

### 1.3 可见度

中心亮纹可见度：

```text
V_A ~ 0.7452539705
V_B ~ 0.7452539969
relative visibility change ~ 3.54e-8
```

所以：

```text
可见度变化也极小。
```

### 1.4 亮纹强度差

真正显著的差别不在条纹位置，而在亮度重标定：

```text
central bright fringe relative intensity shift ~ 9.78e-5
mean bright-peak relative intensity shift      ~ 9.29e-5
support-region pointwise max relative diff     ~ 1.05e-4
```

因此当前 `lambda=1` 下最直接、最现实的 A/B 可测量差异是：

```text
亮纹/热点强度在 1e-4 量级上的系统性增亮
```

而不是：

```text
条纹间距改变
或
明显相位滑移
```

## 2. 地球表面 Schwarzschild 背景

### 2.1 理论层面

地球外部 Schwarzschild 真空区满足：

```text
G_mu nu = 0
```

因此按前面 `024` 的弱展开，
`S_EH[g~] - S_EH[g]` 的一阶 bulk 项

```text
Delta L^(1) ~ (Q/m^2) G_nn
```

在这里**精确为零**。

所以地球表面背景不会给当前 A/B 分裂再引入一个新的“一阶 Einstein 作用量 bulk 差异”。

### 2.2 只看局域读出尺度的潮汐调制

若把局域读出线段长度记为 `L_readout`，
则弱场潮汐量级可用

```text
kappa_E = GM_earth / (c^2 R_earth^3) ~ 1.72e-23 m^-2
```

估计。

在 `L_readout = 1 m` 时，
取几组代表取向，相对地球径向的角度 `theta`：

```text
radial (0 deg)                -> 3.43e-23
horizontal tangential (90 deg)-> 1.72e-23
tilted 45 deg                 -> 8.58e-24
magic angle 54.7 deg          -> 0
```

这里给出的就是潮汐读出调制的量级。

### 2.3 最大可观测差别的方向

最大方向是：

```text
让读出轴尽量沿地球径向
```

也就是实验平面里有一个主要读出方向“竖直”。

最小方向是：

```text
magic angle ~ 54.7 deg
```

在这个方向上，二极潮汐投影抵消。

### 2.4 实际影响

由于平直基线下的等效条纹位移已经只有

```text
~ 8.9e-19 m   (for d_phys = 1 micron)
```

再乘上地球潮汐调制 `~1e-23`，
差别完全可以忽略。

所以：

```text
地球表面 Schwarzschild 背景不会把当前 `lambda=1` 的 A/B 差异放大到实验可见水平。
```

## 3. 已观测纳赫兹引力波背景

### 3.1 采用的代表值

这里把“已观测背景”解释为 NANOGrav 15 年数据报告的纳赫兹背景信号。

其代表特征应变为：

```text
h_c = 2.4^{+0.7}_{-0.6} × 10^-15
at f_ref = 1 / year
```

### 3.2 读出调制

在横波 TT 规范里，最佳情况下的局域读出调制量级是：

```text
max |delta d / d| ~ h_c / 2 ~ 1.2e-15
```

代表取向：

```text
plus-aligned axis      -> max for + polarization
45 deg to plus axis    -> max for x polarization
orthogonal transverse  -> same magnitude, opposite sign for + mode
parallel to propagation-> 0
```

### 3.3 最大可观测差别的方向

若已知瞬时偏振主轴，最大方向是：

```text
让读出轴完全横向，并与 + 模偏振主轴对齐
```

若主要响应的是 `x` 偏振，则最佳方向改为：

```text
相对 + 轴旋转 45 deg
```

### 3.4 实际影响

对当前 `lambda=1` 的 A/B 基线差异，
即使取最佳方向，
也只是乘上

```text
1 + O(1e-15)
```

的调制。

因此不管看：

```text
条纹间距
等效相位移
还是等效位移
```

纳赫兹背景都不会把当前 A/B 差异显著放大。

## 4. 当前最核心结论

这一步最值得记住的是：

```text
1. `lambda=1` 下 A/B 的主要实验量差异，不是条纹平移，而是亮纹强度的 ~1e-4 重标定。
2. 在当前受控 scaffold 下，条纹间距变化 ~ 0，等效相位移只有 ~1e-12 rad。
3. 地球表面 Schwarzschild 真空背景的首阶 Einstein bulk 差异为零，残余潮汐调制 ~1e-23 (for 1 m readout)。
4. 已观测纳赫兹引力波背景的最佳取向调制也只有 ~1e-15。
5. 所以在这两类真实弱背景下，当前 `lambda=1` 的 A/B 差异并不会被放大成条纹间距或相位上的显著实验信号；最现实的可测入口仍是强度型 observables。
```

## 5. 后续最自然的下一步

若要继续追求“真正实验上最敏感的量”，最自然的下一步是：

```text
1. 把 A/B 的中心亮纹、旁瓣亮纹、热点积分计数做成系统时间序列
2. 再把这一套强度型 observables 搬到更强的 lambda，例如 3 或 10
3. 只有在 full backreacted B matter phase 被独立验证后，再认真追条纹相位差
```
