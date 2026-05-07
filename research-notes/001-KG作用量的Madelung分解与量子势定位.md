# 001 KG 作用量的 Madelung 分解与量子势定位

## 目标

在自然单位制 `\hbar=c=1`、闵氏度规 `\eta_{\mu\nu} = diag(1,-1,-1,-1)` 下，从自由复 Klein-Gordon 场的标准作用量出发，做

`ψ = \sqrt{\rho} e^{iS}`

的 Madelung 分解，明确：

- 作用量在 `\rho, S` 变量下的形式
- 连续性方程
- Hamilton-Jacobi 型方程
- 量子势 `Q` 的精确表达
- “量子势项”在作用量和方程中的地位

## 1. 起点：自由复 KG 场作用量

取复标量场 `ψ` 的自由作用量

```math
I[\psi,\psi^*]
= \int d^4x \left(
\partial_\mu \psi^* \partial^\mu \psi
- m^2 \psi^* \psi
\right).
```

对 `\psi^*` 变分可得 Klein-Gordon 方程

```math
(\Box + m^2)\psi = 0,
\qquad
\Box := \partial_\mu \partial^\mu .
```

这里选复标量场而不是实标量场，是因为 Bohm/Madelung 分解里需要非平凡相位 `S`。

## 2. Madelung 分解

写成

```math
\psi = \sqrt{\rho} e^{iS},
\qquad
\psi^* = \sqrt{\rho} e^{-iS}.
```

于是

```math
\partial_\mu \psi
= e^{iS}\left(
\partial_\mu \sqrt{\rho}
+ i \sqrt{\rho}\,\partial_\mu S
\right),
```

```math
\partial_\mu \psi^*
= e^{-iS}\left(
\partial_\mu \sqrt{\rho}
- i \sqrt{\rho}\,\partial_\mu S
\right).
```

两者相乘后交叉项抵消：

```math
\partial_\mu \psi^* \partial^\mu \psi
= \partial_\mu \sqrt{\rho}\,\partial^\mu \sqrt{\rho}
+ \rho \,\partial_\mu S \partial^\mu S.
```

因此作用量改写为

```math
I[\rho,S]
= \int d^4x
\left[
\partial_\mu \sqrt{\rho}\,\partial^\mu \sqrt{\rho}
+ \rho \,\partial_\mu S \partial^\mu S
- m^2 \rho
\right].
```

若进一步用

```math
\partial_\mu \sqrt{\rho}\,\partial^\mu \sqrt{\rho}
= \frac{1}{4\rho}\partial_\mu \rho \partial^\mu \rho,
```

则也可写成

```math
I[\rho,S]
= \int d^4x
\left[
\frac{1}{4\rho}\partial_\mu \rho \partial^\mu \rho
+ \rho \,\partial_\mu S \partial^\mu S
- m^2 \rho
\right].
```

## 3. 对相位 `S` 变分：连续性方程

作用量只通过 `\partial_\mu S` 依赖 `S`，所以

```math
\frac{\partial \mathcal L}{\partial(\partial_\mu S)}
= 2\rho \partial^\mu S.
```

Euler-Lagrange 方程给出

```math
\partial_\mu(\rho \partial^\mu S)=0.
```

这就是 KG 的守恒流连续性方程。若定义

```math
j^\mu = \rho \partial^\mu S,
```

则有

```math
\partial_\mu j^\mu = 0.
```

注意：这里的 `j^0 = \rho \partial^0 S` 不一定正定，这正是后续负粒子数密度问题的来源之一。

## 4. 对振幅 `\rho` 变分：Hamilton-Jacobi 型方程

更方便的做法是把 `R := \sqrt{\rho}` 当作基本变量。此时

```math
\mathcal L
= \partial_\mu R \partial^\mu R
+ R^2 \partial_\mu S \partial^\mu S
- m^2 R^2.
```

对 `R` 变分：

```math
\frac{\partial \mathcal L}{\partial R}
= 2R(\partial_\mu S \partial^\mu S - m^2),
```

```math
\frac{\partial \mathcal L}{\partial(\partial_\mu R)}
= 2\partial^\mu R.
```

因此

```math
2R(\partial_\mu S \partial^\mu S - m^2) - 2\Box R = 0,
```

即

```math
\partial_\mu S \partial^\mu S
= m^2 + \frac{\Box R}{R}.
```

写回 `R=\sqrt{\rho}`：

```math
\partial_\mu S \partial^\mu S
= m^2 + \frac{\Box \sqrt{\rho}}{\sqrt{\rho}}.
```

## 5. 量子势的定义

为了把上式写成“经典 Hamilton-Jacobi 方程 + 修正项”的形式，定义量子势

```math
Q := \frac{\Box \sqrt{\rho}}{\sqrt{\rho}}.
```

则 Hamilton-Jacobi 型方程变成

```math
\partial_\mu S \partial^\mu S
- m^2
- Q
= 0.
```

也可写为

```math
p_\mu p^\mu = m^2 + Q,
\qquad
p_\mu := \partial_\mu S.
```

这说明 Bohm 图景中的“有效质量平方”变成了

```math
m_{\mathrm{eff}}^2 = m^2 + Q.
```

当 `Q` 足够大时，`p_\mu` 的性质可能改变，这与超光速或类空四动量问题直接相关。

## 6. 量子势项在作用量中的地位

这一点需要特别澄清。

### 6.1 直接写在 `(\rho,S)` 作用量里时

原作用量中并没有一个已经显式写成 `\rho Q` 的独立项；它表现为

```math
\partial_\mu \sqrt{\rho}\,\partial^\mu \sqrt{\rho}
```

这一振幅梯度项。

因此，从“原始作用量密度”的角度看，量子势信息编码在振幅梯度能中，而不是单独写成名字叫 `Q` 的势能项。

### 6.2 经过分部积分后

因为

```math
\partial_\mu \sqrt{\rho}\,\partial^\mu \sqrt{\rho}
= \partial_\mu(\sqrt{\rho}\,\partial^\mu\sqrt{\rho})
- \sqrt{\rho}\,\Box\sqrt{\rho},
```

忽略边界项后，作用量可等价写成

```math
I[\rho,S]
\simeq
\int d^4x
\left[
\rho \,\partial_\mu S \partial^\mu S
- m^2 \rho
- \sqrt{\rho}\,\Box\sqrt{\rho}
\right].
```

又因为

```math
- \sqrt{\rho}\,\Box\sqrt{\rho}
= - \rho Q,
```

所以等价地可写成

```math
I[\rho,S]
\simeq
\int d^4x
\rho\left(
\partial_\mu S \partial^\mu S
- m^2
- Q
\right).
```

这里的 `\simeq` 表示差一个边界项。

所以如果后续要谈“把量子势项从作用量中吸收到新度规里”，最自然的出发点是这个分部积分后的等价形式，因为它把 `Q` 显式暴露了出来。

## 7. 当前结论

本步骤已经得到四个后续研究必需的对象：

1. 连续性方程

```math
\partial_\mu(\rho \partial^\mu S)=0.
```

2. Hamilton-Jacobi 型方程

```math
\partial_\mu S \partial^\mu S - m^2 - Q = 0.
```

3. 量子势

```math
Q = \frac{\Box\sqrt{\rho}}{\sqrt{\rho}}.
```

4. 显式量子势作用量形式（忽略边界项）

```math
I[\rho,S]
\simeq
\int d^4x\,
\rho\left(
\partial_\mu S \partial^\mu S - m^2 - Q
\right).
```

这正是下一步研究“能否通过 `g_{\mu\nu} \to \tilde g_{\mu\nu}` 吸收 `Q`”的最小出发点。

## 8. 未决问题

- 应把后续几何化瞄准原始梯度项 `\partial_\mu \sqrt{\rho}\partial^\mu\sqrt{\rho}`，还是瞄准显式写出的 `\rho Q` 项？
- 若采用 disformal transform，允许依赖的对象应是 `\partial_\mu S`、`\partial_\mu \rho`，还是二者都要进入？
- 若新度规吸收了 `Q`，连续性方程是否也应改写为 `\tilde g_{\mu\nu}` 下的协变守恒形式？

## 9. 对下一步的建议

下一步不宜立即写最一般的 disformal ansatz，而应先补上“成功判据”：

- 作用量层面成功是什么意思
- 轨迹层面成功是什么意思
- 观测层面不变具体指什么
- 哪些变换即使能消项也应判为不合格

这样做可以避免后面出现“找到了一个改写，但不知道它算不算物理上成功”的问题。
