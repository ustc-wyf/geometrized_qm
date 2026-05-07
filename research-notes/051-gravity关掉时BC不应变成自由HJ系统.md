# 051 gravity 关掉时 B/C 不应变成自由 HJ 系统

## 1. 当前真正的错误点

最近数值上最大的偏题不是某个时间步或边界系数，
而是我把

```text
gravity -> 0
```

错误地实现成了

```text
让 g~ 变成一个自由的平直背景，
然后单独去解
sqrt(-g~) rho~ ( g~^(μν) ∂_μS ∂_νS - m^2 )
对应的 HJ + continuity 系统。
```

这一步是错的。

## 2. 为什么错

在当前项目里，B/C 分支不是“另一个独立物质理论”，
而是对 A 分支物质理论的几何化重写。

因此：

```text
关掉引力
```

只意味着：

```text
去掉 / 压低 Einstein-Hilbert 或 f(R~) 那一部分的动力学回授，
不意味着把几何化约束一并丢掉。
```

也就是说，B/C 的正确无引力极限应当是：

```text
S_B^(grav off) = S_geom,matter + S_constraint
S_C^(grav off) = S_geom,matter + S_constraint
```

而不是：

```text
S_geom,matter
```

单独裸奔。

## 3. 作用量层面的结构

### A 分支

原始物质理论：

```text
S_A,matter[g, R, S]
= ∫ sqrt(-g) [ g^(μν) ∂_μR ∂_νR + ρ g^(μν) ∂_μS ∂_νS - m^2 ρ ]
```

或等价地给出：

```text
HJ: X + Q = m^2
Continuity: ∇_μ ( ρ ∂^μ S ) = 0
```

### B/C 分支

几何化后的物质作用量是：

```text
S_geom,matter[g~, rho~, S]
= ∫ sqrt(-g~) rho~ ( g~^(μν) ∂_μS ∂_νS - m^2 )
```

但它必须配合约束一起使用，例如：

```text
g~ = G[g, rho, S]
sqrt(-g~) rho~ (C + D X) = sqrt(-g) rho
```

或任何与项目主线等价的 reduced constitutive law。

只有：

```text
S_geom,matter + S_constraint
```

一起看，B/C 才是 A 的重写。

## 4. 因此 gravity -> 0 时真正应满足什么

正确的一致性要求是：

```text
去掉引力动力学后，
B/C 通过“物质项 + 几何化约束”
应当退回与 A 等价的同一物质理论。
```

所以此时不应该出现：

```text
B/C 仍然作为一个自由的、独立于 A 的平直 HJ+continuity 系统继续演化
```

那样做得到的不是项目要比较的理论。

## 5. 当前数值实现意味着什么

当前 `Bflat/Bfull/Cfull` 的 flat-limit 数值器，
本质上做的是：

```text
把约束去掉，
只保留几何化后的 HJ + continuity
```

因此它和 A 不一致，不应该被解释成：

```text
“物理上 B/C 本来就和 A 不同”
```

更准确的结论是：

```text
当前 flat-limit B/C solver 不是目标理论，
而只是一个偏掉的 surrogate system。
```

## 6. 后续数值路线应如何改

后面正确的求解顺序应当是：

```text
1. 先在无引力极限下，把 B/C 的约束重新带回去
2. 验证它们在这一极限下与 A 完全等价
3. 再打开 EH[g~] 或 f(R~) 的动力学回授
4. 这时三者的差异才可解释为引力部分导致的差异
```

## 7. 当前最核心结论

```text
我之前把 “gravity -> 0” 错实现成了 “free g~-HJ system”。
这正是当前 B/C 与 A 不一致的根源之一，
而且是一个概念层面的错误，不是单纯调参能修好的。
```
