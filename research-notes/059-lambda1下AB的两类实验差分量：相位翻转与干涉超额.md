# 059 `lambda=1` 下 A/B 的两类实验差分量：相位翻转与干涉超额

## 目标

把前一轮找到的“原始强度型 observable”

```text
central peak / central disk / bright-region intensity shift
```

进一步改写成更接近真实实验语言的“两种实验设置之间的差分量”：

```text
1. 相位翻转差分
2. 干涉超额计数
```

这类量的优点是：

```text
同一台装置内可切换
可做差分
对共模源强漂移与增益漂移更友好
```

## 定义

### 1. 相位翻转差分

取同样的双束干涉装置，只改相对相位：

```text
phi = 0
phi = pi
```

定义：

```text
D_phi = (N(phi=0) - N(phi=pi)) / (N(phi=0) + N(phi=pi))
```

### 2. 干涉超额计数

取同样的探测窗，比较：

```text
两束同时开
两束分别单开再相加
```

定义：

```text
E_int = (N_both - (N_beam1 + N_beam2)) / (N_beam1 + N_beam2)
```

### 3. 这里真正比较的量

对每一个实验差分量，
我们分别计算 `A` 和 `B` 的预测，
真正要看的理论分裂是：

```text
Delta_AB = observable_B - observable_A
```

## 脚本与输出

脚本：

```text
kg_examples/analyze_ab_differential_experiment_observables.py
```

输出目录：

```text
visualizations/ab_differential_experiment_observables_lambda1/
visualizations/ab_differential_experiment_observables_lambda1_dt025/
visualizations/ab_differential_experiment_observables_lambda1_dt00625/
```

使用的 refined 基线：

```text
lambda_grav = 1
dt_target = 0.0125
nx = nz = 181
nkx = nkz = 81
```

并额外补了：

```text
dt_target = 0.025
dt_target = 0.00625
```

做最小时间步收敛核对。

## 探测窗

为了不把结果绑死在一个任意孔径上，算了两个简单读出窗：

```text
1. central disk: 半径 1.5
2. interference box: |x|<=2, |z|<=2
```

## 1. 相位翻转差分

### 1.1 refined (`dt=0.0125`) 结果

中心圆盘：

```text
D_phi(A) final =  0.022489942607988086
D_phi(B) final =  0.022489940302998706
Delta_AB final = -2.304989379919853e-09
max |Delta_AB|  =  2.934160958992238e-08
t/T at max      =  0.8758241758241758
```

中心干涉盒：

```text
D_phi(A) final = -0.0023981217173558154
D_phi(B) final = -0.002398174913895193
Delta_AB final = -5.3196539377513286e-08
max |Delta_AB|  =  5.3196539377513286e-08
t/T at max      =  1.0
```

### 1.2 dt 稳定性

中心圆盘最终分裂：

```text
dt=0.025   -> -2.469119551179455e-09
dt=0.0125  -> -2.304989379919853e-09
dt=0.00625 -> -2.266545874596293e-09
```

中心干涉盒最终分裂：

```text
dt=0.025   -> -5.352088379905959e-08
dt=0.0125  -> -5.3196539377513286e-08
dt=0.00625 -> -5.311907348400438e-08
```

所以它是稳定的，但量级非常小。

## 2. 干涉超额计数

### 2.1 refined (`dt=0.0125`) 结果

中心圆盘：

```text
E_int(A) final = -0.48875502869600607
E_int(B) final = -0.4887550289351614
Delta_AB final = -2.39155306669403e-10
max |Delta_AB| =  3.3662228959840945e-07
t/T at max     =  0.5406593406593406
```

中心干涉盒：

```text
E_int(A) final = -0.5011990608586778
E_int(B) final = -0.501199178773787
Delta_AB final = -1.1791510912129155e-07
max |Delta_AB| =  1.2542084717459545e-06
t/T at max     =  0.5769230769230769
```

### 2.2 dt 稳定性

中心圆盘最大分裂：

```text
dt=0.025   -> 3.3435582863505786e-07
dt=0.0125  -> 3.3662228959840945e-07
dt=0.00625 -> 3.372462233364182e-07
```

中心干涉盒最大分裂：

```text
dt=0.025   -> 1.2465707220910005e-06
dt=0.0125  -> 1.2542084717459545e-06
dt=0.00625 -> 1.2563076258187422e-06
```

这说明干涉超额计数的 `AB` 分裂也稳定，而且明显强于相位翻转差分。

## 3. 哪一个更值得实验主打

从当前结果看：

```text
相位翻转差分：AB 分裂 ~ 1e-8
干涉超额计数：AB 分裂 ~ 1e-6
```

所以在这两类“实验双设置差分量”里，更值得优先主打的是：

```text
干涉超额计数
```

而且更好的读出窗是：

```text
interference box
```

因为当前最强的分裂出现在这里：

```text
max |Delta_AB(E_int)| ~ 1.26e-06
at t/T ~ 0.577
```

## 4. 为什么它们比原始强度差小这么多

这不是坏事，而是这些差分量的结构使然。

我们前一轮最强的原始信号是：

```text
~ 6.7e-05 的亮纹/热点强度重标定
```

而这里的两类实验差分量都带有：

```text
自归一化
两设置做差
```

所以大部分“共模强度重标定”会被抵消掉。

于是：

```text
实验上更干净
理论分裂也更小
```

这就是现在看到的 tradeoff。

## 5. 当前阶段的最准确结论

如果只问：

```text
有没有更接近真实实验流程的双设置 observable？
```

答案是：

```text
有，而且它们确实给出稳定非零的 A/B 分裂。
```

如果再问：

```text
哪一个最好？
```

答案是：

```text
干涉超额计数优于相位翻转差分；
中心干涉盒优于中心圆盘；
最佳量级大约是 1e-6。
```

但如果问：

```text
它们是否比原始强度型 observable 更强？
```

答案是否定的：

```text
不是。
它们更实验化、更抗共模噪声，
但会把理论信号本身也抵消掉一大截。
```

所以对下一轮实验讨论来说，最合理的口径应该是：

```text
1. 原始强度型 observable 给出最大的理论分裂量级：~6.7e-05
2. 双设置差分 observable 给出更实验友好的 clean channel：最好可到 ~1.3e-06
3. 后续真正要做的是把这两类 observable 一起拿去和实验噪声、漂移、归一化误差比较
```
