# 073 新混合变换下 shared-baseline 源项替换的第一步审计

## 目标

在真正重写 `A/B/C` shared-baseline 比较脚本之前，先单独检查：

```text
若把旧脚本里的 rho_A 源项替换成新混合变换要求的
N~ = sqrt(|g~|) rho~
或 rho~
，当前高斯波包 benchmark 上最关键的 X 场到底长什么样。
```

这一步的目的不是直接给出新物理结果，而是判断：

```text
新变换下能否立刻把旧脚本里的源项替掉，
还是必须先决定 X 过零时的 branch prescription。
```

## 1. 当前 shared-baseline 旧脚本的问题

旧脚本

```text
kg_examples/compare_abc_same_initial_evolution.py
```

在两个层面都把旧变换硬编码进去了：

1. 记录 / 作图时直接用

```text
rho_B = rho_A exp[-(beta_B + tau_B)]
rho_C = rho_A exp[-(beta_C + tau_C)]
```

2. 几何源项推进时，直接把

```text
rho_A, S_A
```

喂给 `rhs_branch_b` 与 `rhs_branch_c`

因此它不能被直接解释成“新变换下的数值结果”。

## 2. 新混合变换下最自然的第一替换对象

若沿当前已经统一的强匹配解释，
并暂取 `kappa` 为固定非零常数，则真正的正测度密度关系是：

```text
N~ := sqrt(|g~|) rho~ = (|X|/|kappa|) rho
```

这里 `X/kappa` 只能保留为一个 **branch diagnostic**：

```text
若要在单个正密度 patch 上做强匹配，必须有 kappa/X > 0 .
```

在当前 shared-baseline 中，原 `g` 表象的 backbone 就是平直 `A` 支，
所以第一步最自然的是先从 `A` benchmark 上审计：

```text
X
Q
N~_signed = (X/kappa) rho_A      （仅作 branch diagnostic）
N~_abs    = |X/kappa| rho_A      （正测度密度候选）
```

## 3. 审计脚本

新增：

```text
kg_examples/audit_mixed_transform_shared_backbone.py
```

它在五个代表时刻

```text
t/T = 0, 0.25, 0.5, 0.75, 1
```

上，从 exact `A` benchmark 计算：

```text
Q = Box sqrt(rho) / sqrt(rho)
X = m^2 + Q
```

并输出 `X/kappa`、`N~_signed`、`N~_abs` 在主支撑区的统计。
其中：

```text
N~_signed 不是物理正密度，只用来诊断 X 是否在固定 kappa 下改变符号；
N~_abs 才对应强匹配口径下的正测度密度大小。
```

## 4. 当前审计结果

使用默认参数：

```text
alpha = 0.5
kappa = 1
support mask: rho > 1e-4 * max(rho)
```

得到：

### t/T = 0

```text
X support min/max/mean = -2561.69 / 181.38 / -12.73
negative fraction      = 0.469
```

### t/T = 0.25

```text
X support min/max/mean = -3170.47 / 156.36 / -15.78
negative fraction      = 0.443
```

### t/T = 0.5

```text
X support min/max/mean = -6342.07 / 123.92 / -23.34
negative fraction      = 0.454
```

### t/T = 0.75

```text
X support min/max/mean = -6802.07 / 156.04 / -31.30
negative fraction      = 0.463
```

### t/T = 1

```text
X support min/max/mean = -4470.49 / 140.29 / -47.03
negative fraction      = 0.475
```

## 5. 当前结论

这一步说明：

```text
在当前高斯波包 benchmark 上，
X 在主支撑区内部就已经大量过零并改变符号。
```

因此：

1. 不能直接把旧脚本里的 `rho_A` 换成 `(X/kappa) rho_A`
   然后就宣布得到了“新变换下的 shared-baseline 数值结果”
2. 更强地说：

```text
当前 benchmark 的主支撑区内部就存在大量 X/kappa < 0 的区域，
因此不存在一个覆盖整个主支撑区的单一固定-kappa正密度强匹配 patch。
```

3. 在真正重跑几何源项之前，必须先决定：

```text
X < 0 的区域到底如何进源项：
- signed
- abs
- positive-part
- 或别的 branch prescription
```

## 6. 这一步的意义

这一审计并没有完成新的 A/B/C 数值结果，
但它完成了一个更基础的判断：

```text
旧 shared-baseline 管线不能被“机械替换源项”地升级到新混合变换；
在固定 kappa 下，必须先决定 branch / patch prescription。
```

先做 branch prescription，
然后才谈重写主脚本和正式重跑，
这是现在最稳的推进顺序。
