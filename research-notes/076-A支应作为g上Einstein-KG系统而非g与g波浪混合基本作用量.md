# 076 A 支应作为 `g` 上 Einstein-KG 系统而非 `g` 与 `g~` 混合基本作用量

## 目标

明确修正当前项目对 `A` 支的定义：

- `A` 支的**基本总作用量**应全部用不带波浪号的 `g` 写；
- `g~` 不是 `A` 支基本作用量中的独立引力变量；
- `g~` 只是在求出 `(g,\rho,S)` 之后，由变换关系事后重构出来的有效几何，用来给出测地线解释。

这一定义与 `B/C` 支有本质区别：

- `B/C` 支的引力作用量直接写在 `g~` 上；
- `A` 支的引力作用量和物质作用量都写在 `g` 上。

## 1. A 支的正确基本作用量

从自由复 Klein-Gordon 场在弯曲时空 `g` 上的标准作用量出发：

```math
S_A[g,\psi]
=
\frac{M_P^2}{2}\int d^4x \sqrt{-g}\,R[g]
\int d^4x \sqrt{-g}
\left(
g^{\mu\nu}\partial_\mu\psi^*\partial_\nu\psi
-m^2\psi^*\psi
\right).
```

若写成 Madelung 变量

```math
\psi=\sqrt{\rho}\,e^{iS},
```

则等价为

```math
S_A[g,\rho,S]
=
\frac{M_P^2}{2}\int d^4x \sqrt{-g}\,R[g]
+
\int d^4x \sqrt{-g}
\left[
g^{\mu\nu}\partial_\mu\sqrt{\rho}\,\partial_\nu\sqrt{\rho}
+\rho\,g^{\mu\nu}\partial_\mu S\partial_\nu S
-m^2\rho
\right].
```

定义

```math
u_\mu:=\partial_\mu S,
\qquad
X:=g^{\mu\nu}u_\mu u_\nu,
\qquad
Q:=\frac{\Box_g\sqrt{\rho}}{\sqrt{\rho}},
```

则在分部积分并忽略边界项后，也可写成

```math
S_A[g,\rho,S]
\simeq
\frac{M_P^2}{2}\int d^4x \sqrt{-g}\,R[g]
+
\int d^4x \sqrt{-g}\,\rho\,(X-m^2-Q).
```

这里 `Q` 仍然是显式的量子势项；这正是 `A` 支与 `B/C` 支的根本区别。

## 2. A 支对应的场方程

对 `S` 变分得到连续性方程

```math
\nabla_\mu(\rho\,u^\mu)=0
\quad\Longleftrightarrow\quad
\partial_\mu(\sqrt{-g}\,\rho\,u^\mu)=0.
```

对 `\sqrt{\rho}` 变分得到 Hamilton-Jacobi 方程

```math
X=m^2+Q.
```

对 `g_{\mu\nu}` 变分得到 Einstein 方程

```math
M_P^2\,G_{\mu\nu}[g]=T_{\mu\nu}^{(A)},
```

其中物质能动张量是标准 Klein-Gordon 形式在 Madelung 变量下的改写：

```math
T_{\mu\nu}^{(A)}
=
\partial_\mu\sqrt{\rho}\,\partial_\nu\sqrt{\rho}
+\partial_\nu\sqrt{\rho}\,\partial_\mu\sqrt{\rho}
+2\rho\,\partial_\mu S\,\partial_\nu S
-g_{\mu\nu}
\left[
g^{\alpha\beta}\partial_\alpha\sqrt{\rho}\,\partial_\beta\sqrt{\rho}
+\rho\,X
-m^2\rho
\right].
```

因此，`A` 支本质上就是：

```text
Einstein 引力 + 在 g 上最小耦合的复 Klein-Gordon 标量场
```

只是我们后续用 `(\rho,S)` 变量而不是 `\psi` 变量来表达它。

## 3. `g~` 在 A 支中的地位

在 `A` 支中，`g~` 不应再作为基本作用量中的独立场。

更准确地说：

1. 先解 `A` 支基本方程，得到 `g,\rho,S`；
2. 再由选定的几何化变换

```math
g~^{\mu\nu}=g~^{\mu\nu}[g,\rho,S]
```

事后重构出 `g~`；
3. 然后检查在这个重构的 `g~` 中，

```math
g~^{\mu\nu}\partial_\mu S\partial_\nu S=\kappa
```

以及相应测地线解释是否成立。

所以：

```text
在 A 支中，g~ 是导出量，不是基本变分变量。
```

## 4. 为什么这不否定 A 支也可做 ADM 分解

`ADM` 分解不是“切换到 `g~`”。
它只是把某个 Lorentz 度规写成初值问题形式。

因此：

- `A` 支若做全动力学数值积分，应对 `g` 做 `ADM` 分解；
- `B/C` 支若做全动力学数值积分，应对 `g~` 做 `ADM` 分解。

两者的区别不在 `ADM` 本身，而在于：

- `A` 支基本引力度规是 `g`；
- `B/C` 支基本引力度规是 `g~`。

## 5. 对后续数值工作的直接影响

这一定义修正后，三支数值任务应明确分开：

### A 支

直接求解

```math
G_{\mu\nu}[g]=T_{\mu\nu}^{(A)}/M_P^2
```

以及 `(\rho,S)` 的连续性方程和 HJ 方程。

### B 支

求解以 `g~` 为基本引力度规的 Einstein-Bohm 系统。

### C 支

求解以 `g~` 为基本引力度规的 `f(R~)`-Bohm 系统。

因此：

```text
A/B/C 三支不是同一个基本作用量在不同记号下的写法，
而是三种不同的引力-物质耦合理论。
```

## 6. 当前结论

本项目后续若讨论 `A` 支，
应统一采用如下口径：

```text
A 支 = g 上的 Einstein-KG 系统；
g~ 只作为从 (g,ρ,S) 导出的有效几何。
```

不应再把

```text
S_EH[g] + L_matter[g~]
```

这种“引力在 `g` 上、物质在 `g~` 上”的混合作用量当作 `A` 支的基本定义。
