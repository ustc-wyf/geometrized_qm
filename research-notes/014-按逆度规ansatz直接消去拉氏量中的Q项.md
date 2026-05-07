# 014 按逆度规 ansatz 直接消去拉氏量中的 Q 项

## 目标

这一份文档完全按当前用户指定的路线做：

```text
不要先从 HJ 方程反推，
而是直接把逆度规 ansatz 代回拉氏量，
要求把量子势项 Q 消去。
```

采用的逆度规 ansatz 是：

```text
g~^(mu nu)
= C g^(mu nu)
 + D g^(mu alpha) g^(nu beta) ∂_alpha S ∂_beta S
```

定义

```text
u_mu = ∂_mu S
u^mu = g^(mu alpha) ∂_alpha S
X = g^(mu nu) ∂_mu S ∂_nu S = g^(mu nu) u_mu u_nu
Q = Box sqrt(ρ) / sqrt(ρ)
```

于是 ansatz 更紧凑地写成

```text
g~^(mu nu) = C g^(mu nu) + D u^mu u^nu
```

## 1. 原拉氏量

从已经得到的等价形式出发：

```text
L_orig = sqrt(-g) ρ ( X - m^2 - Q )
```

这正是当前要被几何化的对象。

## 2. 新拉氏量的最小写法

若采用新的逆度规，并继续采用当前项目已经接受的密度权重约定

```text
sqrt(-g) ρ = sqrt(-g~) ρ~
```

则最自然的新拉氏量写成

```text
L_geom = sqrt(-g~) ρ~ ( g~^(mu nu) ∂_mu S ∂_nu S - m^2 )
       = sqrt(-g) ρ ( g~^(mu nu) ∂_mu S ∂_nu S - m^2 )
```

因此关键只在于计算

```text
g~^(mu nu) ∂_mu S ∂_nu S
```

## 3. 把 ansatz 代回去

直接代入：

```text
g~^(mu nu) ∂_mu S ∂_nu S
= C g^(mu nu) ∂_mu S ∂_nu S
 + D u^mu u^nu ∂_mu S ∂_nu S
```

第一项就是

```text
C X
```

第二项因为

```text
u^mu ∂_mu S = g^(mu alpha) ∂_alpha S ∂_mu S = X
```

所以得到

```text
D X^2
```

于是

```text
g~^(mu nu) ∂_mu S ∂_nu S = C X + D X^2
```

从而新拉氏量变成

```text
L_geom = sqrt(-g) ρ ( C X + D X^2 - m^2 )
```

## 4. 消去 Q 项的条件

现在直接要求

```text
L_geom = L_orig
```

即

```text
C X + D X^2 - m^2 = X - m^2 - Q
```

所以唯一需要满足的条件是

```text
D X^2 + (C - 1) X + Q = 0
```

这就是按你要求“直接把逆度规代回拉氏量，消去 Q 项”后得到的核心方程。

## 5. 一个最简单解

最简单选择就是

```text
C = 1
```

此时

```text
D X^2 = -Q
```

即

```text
D = - Q / X^2
```

于是最简单逆度规 ansatz 变成

```text
g~^(mu nu)
= g^(mu nu)
 - (Q / X^2) u^mu u^nu
```

在这个选择下，

```text
g~^(mu nu) ∂_mu S ∂_nu S
= X - Q
```

所以新拉氏量立刻变成

```text
L_geom = sqrt(-g) ρ ( X - Q - m^2 )
       = sqrt(-g) ρ ( X - m^2 - Q )
       = L_orig
```

### 5.1 结论

因此：

```text
在 X ≠ 0 的区域上，
按你指定的逆度规 ansatz，
确实存在一个非常直接的最简单解，
它能把拉氏量里的 Q 项完全吸收到新逆度规里。
```

## 6. 一个更广泛的解族

更一般地，不必固定 `C = 1`。

核心条件

```text
D X^2 + (C - 1) X - Q = 0
```

说明：

```text
这不是唯一解，
而是一整族解。
```

例如：

### 6.1 任选 C，解 D

```text
D = [ -Q - (C - 1) X ] / X^2
```

### 6.2 任选 D，解 C

```text
C = 1 + ( -Q - D X^2 ) / X
```

所以结论可以直接写成：

```text
按这个逆度规 ansatz，
消去拉氏量中 Q 项的变换不是孤立的，
而是一个一函数自由度的解族。
```

## 7. 对应的 covariant 度规

若要进一步研究测地线、因果结构、可逆性，就需要把 covariant 度规也写出来。

从

```text
g~^(mu nu) = C g^(mu nu) + D u^mu u^nu
```

可得其逆矩阵为

```text
g~_mu nu
= C^(-1) g_mu nu
 - [ D / ( C (C + D X) ) ] u_mu u_nu
```

成立条件是

```text
C ≠ 0
C + D X ≠ 0
```

对于最简单解

```text
C = 1
D = - Q / X^2
```

得到

```text
g~_mu nu
= g_mu nu
 + [ Q / ( X (X - Q) ) ] u_mu u_nu
```

## 8. 这个结果说明什么

这一轮得到的是一个比我们前面绕 HJ 方程更直接的结果：

```text
在拉氏量层面，
Q 的吸收其实是非常直接的代数条件。
```

按你指定的类 II 风格逆度规 ansatz，
只要满足

```text
D X^2 + (C - 1) X + Q = 0
```

就能把

```text
X - m^2 - Q
```

改写成

```text
g~^(mu nu) ∂_mu S ∂_nu S - m^2
```

所以：

```text
至少在“拉氏量中消去 Q 项”这一步，
类 II 风格 ansatz 不但没有失败，
反而是成功且相当干净的。
```

## 9. 还没有解决的事

当然，这一步成功还不代表整个理论已经完成。

后面仍然必须继续检查：

1. 用这个新拉氏量变分出来的方程是否与原理论等价
2. 连续性方程在新几何下怎么解释
3. 轨迹是否真的能写成测地线
4. 经典极限、可逆性、因果结构是否良好

但是：

```text
“把 Q 从拉氏量里直接消去”这一步，
现在已经算出来了，
而且存在一个更广泛的解族。
```

## 10. 下一步建议

下一步最自然的是：

```text
固定这个最简单解
C = 1
D = - Q / X^2
```

然后在这个具体几何上继续检查：

1. 新拉氏量变分后的方程
2. 它与原 KG-Madelung 方程组的对应关系
3. 它是否真的支持测地线解释
