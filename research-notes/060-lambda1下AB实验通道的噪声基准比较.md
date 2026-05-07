# 060 `lambda=1` 下 A/B 实验通道的噪声基准比较

## 目标

把当前几条最重要的实验通道放进统一的噪声基准里比较：

```text
1. raw center-peak intensity shift
2. raw interference-box integrated intensity shift
3. phase-flip differential
4. interference-excess differential
```

并回答最实际的问题：

```text
哪一条通道在统计计数和系统误差上最划算？
```

## 用到的信号量级

来自前面 refined 结果：

```text
raw center peak              ~ 6.761080361107943e-05
raw interference box         ~ 6.557413714432494e-05
phase-flip differential box  ~ 5.3196539377513286e-08
interference excess box max  ~ 1.2542084717459545e-06
```

这里对差分通道都取了当前更有利的版本：

```text
phase flip: interference box
interference excess: interference box 且取时间峰值
```

## 1. 统计噪声模型

### 1.1 raw intensity

对相对强度差 `delta`，
若读出窗内总计数为 `N_window`，
Poisson shot noise 给出：

```text
sigma_rel ~ 1 / sqrt(N_window)
```

所以要达到 `k sigma`：

```text
N_window ~ (k / delta)^2
```

### 1.2 phase-flip differential

定义：

```text
D_phi = (N0 - Npi) / (N0 + Npi)
```

若 `N0, Npi` 都服从 Poisson 统计，则：

```text
Var(D_phi) ~ (1 - D_phi^2) / N_total
```

其中 `N_total = N0 + Npi`。

### 1.3 interference excess

定义：

```text
E_int = (N_both - (N1 + N2)) / (N1 + N2)
```

若三组计数都服从独立 Poisson 统计，则：

```text
Var(E_int) ~ (2 + 3E + E^2) / B
```

其中：

```text
E = E_int(A) at the operating point
B = N1 + N2
```

在当前最优点，

```text
E ~ -0.50035
```

因此系数约为：

```text
2 + 3E + E^2 ~ 0.7486
```

## 2. 统计计数要求

### 2.1 raw center peak

```text
1 sigma -> 2.19e8 counts
3 sigma -> 1.97e9 counts
5 sigma -> 5.47e9 counts
```

### 2.2 raw interference box

```text
1 sigma -> 2.33e8 counts
3 sigma -> 2.09e9 counts
5 sigma -> 5.81e9 counts
```

### 2.3 phase-flip differential (interference box)

```text
1 sigma -> 3.53e14 total counts across the two phase states
3 sigma -> 3.18e15
5 sigma -> 8.83e15
```

### 2.4 interference excess (interference box)

对分母计数 `B = N1 + N2`：

```text
1 sigma -> 4.76e11
3 sigma -> 4.29e12
5 sigma -> 1.19e13
```

若换成三组实验总计数：

```text
1 sigma -> 7.14e11 total counts
3 sigma -> 6.43e12
5 sigma -> 1.79e13
```

## 3. 系统误差容忍度

这里给两个预算：

```text
signal-size budget : 单个系统误差本身等于理论信号
10% budget         : 单个系统误差只吃掉理论信号的 10%
```

### 3.1 raw 通道的相对增益/源强残余漂移

raw center peak：

```text
signal-size budget -> 6.76e-05
10% budget         -> 6.76e-06
```

raw interference box：

```text
signal-size budget -> 6.56e-05
10% budget         -> 6.56e-06
```

### 3.2 phase-flip differential 的双相位态 pair imbalance

当前最优 `phase-flip box`：

```text
signal-size budget -> 5.32e-08
10% budget         -> 5.32e-09
```

也就是说：

```text
phase flip 若想拿来分辨 A/B，
两种相位设置之间的相对不平衡必须压到 1e-8 甚至 1e-9
```

这已经非常苛刻。

### 3.3 interference excess 的单束归一化误差

当前最优 `interference excess box`：

```text
signal-size budget -> 2.51e-06
10% budget         -> 2.51e-07
```

这意味着：

```text
若要把 interference excess 当主信号，
单束参考项 N1+N2 的相对归一化必须好到 ppm 甚至亚 ppm。
```

## 4. 结论：哪条通道最划算

### 4.1 phase flip 基本出局

原因很直接：

```text
统计上太贵：5 sigma 需要 ~8.8e15 counts
系统上也太苛刻：pair imbalance 要压到 ~1e-8
```

所以：

```text
phase-flip differential 更适合作为概念性 cross-check，
不适合作为当前项目的第一主通道。
```

### 4.2 interference excess 是“最干净的差分通道”，但不便宜

它的优点：

```text
确实是稳定非零
天然是实验双设置差分量
比 phase flip 强两个数量级
```

但代价是：

```text
5 sigma 仍需 ~1.8e13 total counts
单束归一化误差最好压到 ~2.5e-07
```

所以它更像：

```text
实验上的 clean companion channel
```

而不是最便宜的发现通道。

### 4.3 真正最划算的主通道仍是 raw intensity

在当前最简单的噪声模型下，
raw 通道虽然更依赖绝对归一化，
但它的信号大约大了 `50` 倍，
这个优势太明显了。

具体体现在：

```text
5 sigma counts only ~5e9
10% systematic budget only asks for ~6e-06 level residual drift
```

与 `interference excess` 比较：

```text
raw     : ~5.5e9 counts, ~6e-06 drift budget
E_int   : ~1.8e13 counts, ~2.5e-07 normalization budget
```

所以至少按当前这套 baseline：

```text
raw intensity channel is easier both statistically and systematically
```

## 5. 当前最合理的实验策略

如果现在就要给实验组一个分层建议，我会这样排：

### 第一优先级

```text
raw integrated intensity channel
```

更具体地说，
比起单像素中心峰，
我会更偏向：

```text
interference-box integrated intensity
```

因为它和中心峰有几乎一样大的信号量级，
但对像素配准和局域尖点更鲁棒。

### 第二优先级

```text
interference excess
```

把它作为：

```text
共模抑制更好的 companion channel
```

如果 raw 通道先给出 hint，
`E_int` 是很自然的后续验证线。

### 不建议当前主打

```text
phase-flip differential
```

## 6. 最后一句话

当前噪声基准比较之后，
我们已经可以把实验路线压成一句非常实用的话：

```text
若追求“最先看到信号”，优先做 raw integrated intensity；
若追求“更干净的差分交叉验证”，再做 interference excess；
phase flip 暂时不值得作为主通道投入。
```
