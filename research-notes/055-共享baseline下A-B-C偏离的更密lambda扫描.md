# 055 共享 baseline 下 A/B/C 偏离的更密 lambda 扫描

## 目标

在已经修正 `gravity -> 0` 极限并建立共享 baseline 之后，
对 `lambda_grav` 做一轮比之前更密的 continuation 扫描，
判断：

```text
1. B 从 A 开始显著偏离的大致量级
2. C 从 B 开始显著偏离的大致量级
3. 强刚度区需要的 C 子循环级别
```

## 1. 扫描设置

主脚本：

```text
kg_examples/scan_abc_comparison_backbone.py
```

底层 continuation solver：

```text
kg_examples/simulate_bc_small_gravity_backbone.py
```

当前参数：

```text
alpha = 0.5
mp = 300
ell = 0.02
dt_target = 0.025
```

更密的扫描网格：

```text
lambda_grav =
0,
1e-4, 3e-4,
1e-3, 3e-3,
1e-2, 3e-2,
1e-1, 3e-1,
1, 3, 10
```

分段子循环策略：

```text
lambda <= 1e-1    -> c_substeps = 1
lambda <= 3e-1    -> c_substeps = 2
lambda <= 1       -> c_substeps = 4
lambda > 1        -> c_substeps = 8
```

## 2. 主要结果

汇总文件：

```text
visualizations/abc_comparison_backbone_scan/summary.json
visualizations/abc_comparison_backbone_scan/thresholds.json
```

### 2.1 B vs A

`B_vs_A` 非常平滑，几乎线性随 `lambda_grav` 增长：

```text
lambda = 1e-4   -> B_vs_A ~ 3.13e-10
lambda = 1e-3   -> B_vs_A ~ 3.13e-9
lambda = 1e-2   -> B_vs_A ~ 3.13e-8
lambda = 1e-1   -> B_vs_A ~ 3.13e-7
lambda = 1      -> B_vs_A ~ 3.13e-6
lambda = 3      -> B_vs_A ~ 9.40e-6
lambda = 10     -> B_vs_A ~ 3.14e-5
```

当前最稳的阈值说法：

```text
B_vs_A >= 1e-6   首次出现在 lambda ~ 1
B_vs_A >= 1e-5   首次出现在 lambda ~ 10
```

### 2.2 C vs B

`C` 相对 `B` 的偏离在弱场区极小，
到中等 `lambda` 后进入 `10^-6 ~ 10^-5`，
到更强处进入几 `10^-5`：

```text
lambda = 3e-2   -> C_vs_B ~ 1.11e-6
lambda = 1e-1   -> C_vs_B ~ 9.90e-6
lambda = 3e-1   -> C_vs_B ~ 6.46e-6
lambda = 1      -> C_vs_B ~ 7.21e-6
lambda = 3      -> C_vs_B ~ 2.37e-5
lambda = 10     -> C_vs_B ~ 9.49e-5
```

当前最稳妥的阈值说法是：

```text
C_vs_B >= 1e-5 的稳健 crossing 出现在 lambda ~ 3
在 lambda <= 10 内，C_vs_B 仍未超过 1e-4
```

不直接把 `lambda = 0.1` 当成第一次 crossing，
是因为：

```text
0.1 -> 0.3 -> 1
```

这一段里 `C_vs_B` 并不严格单调，
而该区间又伴随着 `c_substeps` 的逐步提高，
所以不适合把某个单点放大直接读成纯物理结论。

### 2.3 C vs A

```text
lambda = 3e-2   -> C_vs_A ~ 1.21e-6
lambda = 1e-1   -> C_vs_A ~ 1.02e-5
lambda = 3e-1   -> C_vs_A ~ 7.40e-6
lambda = 1      -> C_vs_A ~ 1.03e-5
lambda = 3      -> C_vs_A ~ 3.31e-5
lambda = 10     -> C_vs_A ~ 1.26e-4
```

因此：

```text
C_vs_A >= 1e-4   首次出现在 lambda ~ 10
```

## 3. 解释

### 3.1 B 分支已经很干净

`B_vs_A` 的增长形状非常理想：

```text
平滑
单调
近似线性
```

这进一步支持：

```text
B branch numerics 现在已经足够可信，
可以作为后续物理比较的稳定支。
```

### 3.2 C 分支已经能做比较，但必须带着 solver-control 一起读

`C` 现在已经不是“完全坏掉”的支，
因为在受控子循环下它能稳定给出有限偏离。

但更准确的说法是：

```text
C 的物理偏离已经开始可比较，
不过 bridge 区间仍然和子循环如何开启有明显耦合。
```

也就是说：

```text
弱场区：C 和 B 先天地很接近
中场区：既有真实偏离，也仍有 solver-control 影响
强场区：只有在 subcycling 开足后，C 的偏离才值得读成物理
```

### 3.3 代表性比较点

当前最自然的代表点是：

```text
lambda ~ 3e-2   : C-B 刚进入 1e-6 量级
lambda ~ 3e-1   : bridge 区，开始需要轻度子循环
lambda ~ 3      : C-B 进入几 1e-5
lambda ~ 10     : C-A 进入 1e-4 量级
```

## 4. 计算代价

单独补跑：

```text
lambda = 3, c_substeps = 8
```

在当前机器上大约耗时：

```text
real time ~ 257 s
peak memory ~ 1.1e8 bytes
```

这说明：

```text
进入中强引力且需要较高子循环后，
dense continuation scan 已经不是轻量脚本。
```

因此 continuation 脚本做成可恢复/复用已有 case 是必要的。

## 5. 当前最核心结论

```text
B 分支的偏离是平滑、单调、近线性的；
C 分支在当前受控求解下也已经能比较，
但其“显著脱离 B”的稳健量级更接近 lambda ~ 3，
而不是更弱区间里那些还混着 solver 刚度效应的单点放大。
```

## 6. 下一步建议

最合理的下一步是：

```text
1. 从这轮 dense scan 里挑 3e-2 / 3e-1 / 3 / 10 四个代表点
2. 对比它们的 rho 图、差分图和中心线
3. 再决定“视觉上有意义的偏离”应把阈值落在哪一档
```
