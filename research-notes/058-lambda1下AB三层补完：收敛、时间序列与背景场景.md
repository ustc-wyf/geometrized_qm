# 058 `lambda=1` 下 A/B 三层补完：收敛、时间序列与背景场景

## 目标

把前一轮还没有补完的三层内容真正落地：

```text
1. 直接对候选 A/B observable 做数值收敛性检查
2. 不只看最终时刻单条中心截线，而是补时间序列和积分型 readout
3. 把弱曲背景从“单点估计”升级成场景梯度 + 取向扫描
```

这一步仍然只讨论 `A/B`，并继续沿用当前最可信的 shared-backbone 口径：

```text
A: 可信 PDE 骨架动态演化
B: 在同一组 A 的 rho,S 历史上推进 Einstein 几何并重构 rho_B
```

## 脚本与输出

新增脚本：

```text
kg_examples/scan_ab_observable_convergence.py
kg_examples/analyze_ab_time_observables.py
kg_examples/compare_ab_background_scenarios.py
```

新增输出目录：

```text
visualizations/ab_observable_convergence/
visualizations/ab_time_observables_lambda1/
visualizations/ab_background_scenarios_lambda1/
```

## 1. 先把“不是数值误差”钉住

### 1.1 收敛扫描设置

时间步扫描：

```text
dt_target = 0.05, 0.025, 0.0125, 0.00625
nx = nz = 181
```

网格扫描：

```text
nx = nz = 141, 181, 221
dt_target = 0.0125
```

共同参数：

```text
lambda_grav = 1
alpha = 0.5
inner = 10
outer = 15
nkx = nkz = 81
```

### 1.2 关键结果

用最细时间步参考值：

```text
central_peak_rel_intensity_shift = 6.761437198341409e-05
mean_peak_rel_intensity_shift    = 6.561680454841884e-05
support_pointwise_rel_diff_max   = 7.820840064539984e-05
visibility_rel_diff              = 2.375198983714465e-08
effective_phase_shift_rad        = -1.9427825797490044e-12
equiv_shift_over_spacing         = -3.092034509198784e-13
```

数值地板取“最后两档 refinement 的差值上界”：

```text
central peak floor  ~ 3.57e-09
mean peak floor     ~ 6.70e-08
support max floor   ~ 2.21e-08
visibility floor    ~ 4.01e-09
phase floor         ~ 9.21e-11
```

于是信号与数值地板的比值约为：

```text
central peak  ~ 1.89e4
mean peak     ~ 9.8e2
support max   ~ 3.5e3
visibility    ~ 5.9
phase         ~ 2.1e-2
```

### 1.3 结论

当前最稳的 A/B observable 已经非常明确：

```text
最可信：中心亮纹强度差、中心热点积分强度差、主支撑区相对差
不适合作为主 observable：等效相位移、等效条纹位移
可见度：有信号，但没有强度型 observable 那么干净
```

这一步的重要含义是：

```text
“A/B 在 lambda=1 下确实存在非零强度型差异”
已经不再只是视觉印象，
而是一个通过 dt/grid refinement 直接支撑的数值结论。
```

## 2. 时间序列与积分型 readout

### 2.1 新 readout

除了最终时刻中心截线，还跟踪了：

```text
1. 全内区 B_vs_A_relL1
2. 中心点强度差
3. 中心圆盘积分强度差
4. 亮区总强度差
5. 干涉中心盒内积分强度差
6. 中心峰、平均峰、support max 的时间序列
```

默认积分窗：

```text
central_disk_radius = 1.5
bright_fraction = 0.5 * A_max
interference_box_halfwidth = 2.0
```

### 2.2 差异何时开始长出来

首次过阈时间：

```text
B_vs_A_relL1 >= 1e-6            at t/T ~ 0.3725
center_peak_rel_shift >= 1e-5   at t/T ~ 0.4066
central_disk_rel_shift >= 1e-5  at t/T ~ 0.4165
bright_region_rel_shift >= 1e-5 at t/T ~ 0.1297
```

这说明：

```text
1. A/B 的分裂不是初始时刻就有，而是在双束重叠过程中建立
2. 积分型亮区强度 readout 比“单个中心峰”更早开始响应
```

### 2.3 最终时刻的积分型读出

最终时刻得到：

```text
center_point_rel_shift      = 6.761080361107943e-05
central_disk_rel_shift      = 6.683999408006137e-05
bright_region_rel_shift     = 6.625028280014768e-05
interference_box_rel_shift  = 6.557413714432494e-05
central_peak_rel_shift      = 6.761080361123216e-05
mean_peak_rel_shift         = 6.561312672010778e-05
support_pointwise_rel_diff  = 7.818634974182434e-05
```

可以看到它们都稳定落在：

```text
6.6e-05 ~ 7.8e-05
```

这一窄范围内，说明“亮纹强度重标定”不是某个孤立像素的偶然产物，而是多种 readout 一致看到的结构。

## 3. 背景场景梯度与取向扫描

### 3.1 现在使用的平直基线

背景脚本已改为优先读取本轮 refinement 后的平直基线：

```text
central_peak_rel_intensity_shift = 6.761437198341409e-05
```

而不再沿用上一轮较粗口径的 `~9.8e-05`。

### 3.2 场景梯度

当前纳入的场景：

```text
1. Earth surface
2. LEO 400 km
3. GEO
4. Sun-Earth L1
5. 可调近源质量：1 tonne @ 10 cm, readout = 1 cm
6. nanohertz GW background
```

其中前五类都仍然是：

```text
beam region 内的 vacuum tidal 背景
```

所以它们仍然只是“平直 AB 分裂之上的读出调制”，不是新的 bulk Einstein 差异。

### 3.3 代表结果

最大取向一律是轴向读出：

```text
Earth surface axial modulation  ~ 3.43e-23
LEO axial modulation            ~ 2.86e-23
GEO axial modulation            ~ 1.18e-25
Sun-Earth L1 axial modulation   ~ 3.54e-30
1 tonne / 10 cm axial modulation~ 1.49e-25   (for 1 cm readout)
GW best-case modulation         ~ 1.2e-15
```

把这些调制乘到当前主 observable

```text
central_peak_rel_intensity_shift ~ 6.76e-05
```

上，得到的“背景诱导改变量”与当前数值地板相比：

```text
Earth surface   -> ~6.5e-19 of the numerical floor
LEO             -> ~5.4e-19
GEO             -> ~2.2e-21
Sun-Earth L1    -> ~6.7e-26
1 tonne / 10 cm -> ~2.8e-21
nanohertz GW    -> ~2.3e-11
```

### 3.4 结论

在当前真空弱场估计口径下：

```text
1. Earth/L1/GW 这些背景都不会把现有 A/B 强度差进一步放大成新的可观测入口
2. 新增的“可调近源质量”场景在实验设计上更灵活，但在当前 vacuum tidal 近似下仍然远低于数值地板
3. 若真要让背景本身成为主角，就不能再停留在 vacuum tidal 调制估计，而必须进入 full curved-background / non-vacuum overlap 的真正动力学
```

## 4. 当前阶段最扎实的结论

现在可以把 `lambda=1` 下 A/B 的现状压成一句更干净的话：

```text
A/B 在 lambda=1 下最稳的非零差异，不是条纹位置，也不是相位，
而是 6.7e-05 量级的亮纹/热点强度重标定；
它对 dt 和 grid refinement 稳定，
并且在多种积分型 readout 上彼此一致。
```

## 5. 还剩下什么没有做

这轮三层补完之后，仍然还没完成的是：

```text
1. full curved-background self-consistent PDE evolution
2. non-vacuum overlap 背景下的 bulk Einstein 差异
3. 真正实验噪声模型下的 signal-to-noise 估计
```

也就是说：

```text
“存在一个稳定非零 A/B observable”
这件事现在已经比较扎实；

但

“这个 observable 在哪些真实实验里最容易被看到”
还需要下一轮专门把噪声基准、器件系统误差和实验几何一起拉进来。
```
