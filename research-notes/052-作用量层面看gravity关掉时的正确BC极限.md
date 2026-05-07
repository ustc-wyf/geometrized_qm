# 052 作用量层面看 gravity 关掉时的正确 B/C 极限

## 1. 要回答的问题

用户指出的核心一致性要求是：

```text
当引力部分关掉时，
A / B / C 应该退回同一个物质作用量。
```

这是对的，不是额外要求。

因此，后续所有数值比较都必须以这个极限为准。

## 2. 当前最重要的澄清

此前的一个根本偏题是：

```text
把 gravity -> 0
错实现成了
“让 g~ 平直，然后单独解 HJ + continuity”
```

这不对应当前项目要比较的理论。

因为在本项目里：

```text
B/C 不是独立物质理论，
而是 A 的几何化重写。
```

所以关掉引力之后，B/C 仍然必须保留：

```text
S_geom,matter + S_constraint
```

而不是只留下

```text
S_geom,matter
```

## 3. 最简单几何化方案下的无引力极限

当前主线里的最简单 rank-1 几何化方案是：

```text
g~^(μν) = g^(μν) + (Q / X^2) u^μ u^ν
u_μ = ∂_μ S
X = g^(μν) ∂_μ S ∂_ν S
Q = - □ sqrt(rho) / sqrt(rho)
```

以及密度变换：

```text
sqrt(-g~) rho~ = sqrt(-g) rho X / m^2
```

在这个方案下有一个关键恒等式：

```text
g~^(μν) u_ν = (1 + Q/X) u^μ = (m^2 / X) u^μ
```

因此

```text
sqrt(-g~) rho~ g~^(μν) u_ν
= sqrt(-g) rho u^μ
```

这是作用量层面的真正无引力等价结构。

## 4. 这意味着什么

这说明：

```text
当 gravity -> 0 时，
正确的 B/C 极限不是“独立演化的平直 HJ+continuity 流”，
而是：
先有 A 的物质解，
再由约束代数地重构出 g~ 和 rho~。
```

换句话说：

```text
无引力极限下，
B/C 不该作为一个独立 PDE 系统单独往前跑。
```

它们在这个极限下应当是：

```text
A 的重写
```

## 5. 数值检验

新增脚本：

```text
kg_examples/check_constrained_gravity_off_equivalence.py
```

它做的是：

1. 用 `A` 的精确 free KG benchmark 生成 `rho, S`
2. 按约束重构 `g~` 与 `rho~`
3. 检查电流恒等式

输出在：

```text
kg_examples/outputs/gravity_off_equivalence_summary.json
```

最关键的结果不是 HJ 残差（那一步当前数值导数还比较粗糙），
而是电流恒等式：

### alpha = 0.5

```text
current_identity_residual_t ~ 1e-16
current_identity_residual_x ~ 1e-17
current_identity_residual_z ~ 1e-17
```

### alpha = 1.0

```text
current_identity_residual_t ~ 1e-16
current_identity_residual_x ~ 1e-17
current_identity_residual_z ~ 1e-17
```

这说明：

```text
按正确约束重构后，
B/C 的无引力极限电流与 A 完全一致，
精度上就是机器误差级别。
```

## 6. 当前应撤回的东西

因此，当前应正式撤回：

```text
把 Bflat / Cfull 这类“去掉引力后仍独立跑 HJ+continuity”的 flat-limit solver
当作目标理论的无引力极限。
```

它们只能算：

```text
一个偏掉的 surrogate system
```

不能再拿来判断

```text
B/C 和 A 的物理差异
```

## 7. 后续数值路线

后面正确的比较顺序应改成：

```text
1. 无引力极限：
   A 动力学求解
   + 由约束代数重构 B/C

2. 打开引力动力学后：
   再让 B/C 的几何场真正演化

3. 此时三者偏离，
   才能解释为引力动力学带来的偏离
```

## 8. 当前最核心结论

```text
作用量层面检查表明：
gravity -> 0 时，
正确的 B/C 极限应当与 A 完全一致；
我之前实现的 flat-limit B/C solver 目标错了。
```
