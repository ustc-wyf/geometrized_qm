# 018 最小混合逆度规下的连续性方程 no-go

## 目标

继续沿当前主线推进：

```text
在最简单消项解
g~^(mu nu) = g^(mu nu) + (Q/X^2) u^mu u^nu
```

已经成功消去拉氏量中的 `Q` 且成功重现 HJ 方程之后，
我们现在尝试加入“最少量的振幅结构”，
看能否把连续性方程也一起修回来。

本文件考虑的最小混合逆度规 ansatz 是：

```text
g~^(mu nu)
= C g^(mu nu)
 + D u^mu u^nu
 + E r^mu r^nu
 + F ( u^mu r^nu + r^mu u^nu )
```

其中

```text
u_mu = ∂_mu S
r_mu = ∂_mu sqrt(ρ)
X = g^(mu nu) u_mu u_nu
Y = g^(mu nu) u_mu r_nu
Z = g^(mu nu) r_mu r_nu
```

并继续采用密度变换约定：

```text
sqrt(-g) ρ = sqrt(-g~) ρ~
```

目标是同时要求：

1. 拉氏量中的 `Q` 精确消去
2. 连续性方程精确回到原形式

结果将表明：

```text
在这个最小混合 ansatz 下，
这两件事对泛型场构型不能同时成立，
除非 Q = 0。
```

## 1. 拉氏量消项条件

原拉氏量是

```text
L_orig = sqrt(-g) ρ ( X - m^2 + Q )
```

新拉氏量在当前约定下写成

```text
L_geom = sqrt(-g) ρ ( g~^(mu nu) u_mu u_nu - m^2 )
```

先计算

```text
g~^(mu nu) u_mu u_nu
```

逐项收缩得到：

```text
C g^(mu nu) u_mu u_nu = C X

D u^mu u^nu u_mu u_nu = D X^2

E r^mu r^nu u_mu u_nu = E Y^2

F ( u^mu r^nu + r^mu u^nu ) u_mu u_nu = 2 F X Y
```

所以

```text
g~^(mu nu) u_mu u_nu
= C X + D X^2 + E Y^2 + 2 F X Y
```

要使拉氏量精确消去 `Q`，必须满足

```text
C X + D X^2 + E Y^2 + 2 F X Y - m^2
= X - m^2 + Q
```

即

```text
(C - 1) X + D X^2 + E Y^2 + 2 F X Y = Q
```

记为

```text
(L)
```

## 2. 连续性方程回到原形式的条件

由新作用量对 `S` 变分得到的新连续性方程是

```text
∂_mu [ sqrt(-g) ρ g~^(mu nu) u_nu ] = 0
```

而原连续性方程是

```text
∂_mu [ sqrt(-g) ρ u^mu ] = 0
```

若要让新旧连续性方程对泛型场构型精确重合，
一个自然且足够强的要求是：

```text
g~^(mu nu) u_nu = u^mu
```

现在把它展开。

先逐项算

```text
C g^(mu nu) u_nu = C u^mu

D u^mu u^nu u_nu = D X u^mu

E r^mu r^nu u_nu = E Y r^mu

F ( u^mu r^nu + r^mu u^nu ) u_nu
= F ( Y u^mu + X r^mu )
```

因此

```text
g~^(mu nu) u_nu
= [ C + D X + F Y ] u^mu
 + [ E Y + F X ] r^mu
```

若要它等于 `u^mu`，则必须满足

```text
C + D X + F Y = 1
```

以及

```text
E Y + F X = 0
```

分别记为

```text
(C1)
(C2)
```

## 3. 联立三条条件

现在把 `(C1)` 与 `(C2)` 代回拉氏量消项条件 `(L)`。

由 `(C1)` 可得

```text
C - 1 = - D X - F Y
```

代入 `(L)`：

```text
(- D X - F Y) X + D X^2 + E Y^2 + 2 F X Y = Q
```

整理：

```text
- D X^2 - F X Y + D X^2 + E Y^2 + 2 F X Y = Q
```

得到

```text
E Y^2 + F X Y = Q
```

再利用 `(C2)`：

```text
F X = - E Y
```

所以

```text
E Y^2 + Y (F X) = E Y^2 - E Y^2 = 0
```

于是最终得到

```text
Q = 0
```

## 4. 结论

因此，在这个最小混合逆度规 ansatz 下：

```text
若要求
1. 拉氏量中的 Q 精确消去
2. 新连续性方程对泛型场构型精确回到原形式
```

那么代数上必然得到

```text
Q = 0
```

也就是说：

```text
除非量子势本来就不存在，
否则这两个要求不能同时满足。
```

这就是当前 ansatz 下的一个小型 no-go 定理。

## 5. 这个 no-go 的含义

它说明：

```text
把振幅结构以最线性的方式
E r^mu r^nu + F(u^mu r^nu + r^mu u^nu)
加进逆度规，
还不够同时完成两件事：
```

1. 让拉氏量中 `Q` 精确消失
2. 让连续性方程精确保持原形式

换句话说：

```text
最简单消项解失败的根源不是“还没加振幅项”，
而是“只加到这一级别还不够”。
```

## 6. 需要怎样突破这个 no-go

既然当前最小混合 ansatz 已经被代数上卡死，那么后续若要继续推进，只能放宽至少一项假设。

可能的方向有三类：

### 6.1 放宽“连续性方程必须逐字重合”

也许新几何下真正应守恒的对象不是原来的 `J^mu`，而是某个几何修正后的流。

这是最保守的放宽。

### 6.2 放宽逆度规 ansatz

例如让逆度规依赖更高阶结构，如：

```text
Q
∂_mu Q
Box sqrt(ρ)
```

或者更一般的非线性组合。

### 6.3 放宽作用量 ansatz

也许仅仅替换 kinetic metric 还不够，
还需要在几何侧作用量里加入额外几何标量或约束项。

## 7. 当前最稳的判断

到这一步，当前主线已经可以很清楚地收缩为：

```text
最简单消项解已经完成了 HJ / 测地线部分；
若要把守恒结构也一起修好，
必须超出这个最小混合逆度规 ansatz。
```

所以接下来的真正选择是：

```text
要么放宽“连续性方程必须原样重合”的要求，
要么升级 ansatz，
要么升级作用量。
```
