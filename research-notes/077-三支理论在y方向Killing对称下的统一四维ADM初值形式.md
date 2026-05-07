# 077 三支理论在 `y` 方向 Killing 对称下的统一四维 `ADM` 初值形式

## 目标

为避免把 `3+1 -> 2+1` 约化、共形重参数化和引力作用量选择混在一起，
本笔记把 `A/B/C` 三支统一写回同一个四维 `ADM` 初值问题形式。

这里只固定一条对称性假设：

```text
所有场对 y 无依赖，且 g_(iy)=0, N^y=0
```

其中

```text
i,j ∈ {x,z}.
```

于是四维度规统一写成

```math
ds^2
=
N^2 dt^2
- h_{ij}(dx^i+N^i dt)(dx^j+N^j dt)
- e^{2\beta}dy^2.
```

对应的三维空间度规是

```math
\gamma_{IJ}
=
\begin{pmatrix}
h_{xx} & h_{xz} & 0 \\
h_{xz} & h_{zz} & 0 \\
0      & 0      & e^{2\beta}
\end{pmatrix},
\qquad I,J\in\{x,z,y\}.
```

因此

```math
\sqrt{\gamma}=e^\beta\sqrt{h},
\qquad
h:=h_{xx}h_{zz}-h_{xz}^2.
```

## 1. 几何变量与外曲率

标准四维 `ADM` 变量为

```text
N,\quad N^x,\quad N^z,\quad h_{xx},h_{xz},h_{zz},\beta
```

三维空间外曲率定义为

```math
K_{IJ}
:=
-\frac{1}{2N}
\left(
\partial_t\gamma_{IJ}
-D_I N_J
-D_J N_I
\right).
```

在当前对称性下，可写成

```math
K_{IJ}
=
\begin{pmatrix}
K_{xx} & K_{xz} & 0 \\
K_{xz} & K_{zz} & 0 \\
0 & 0 & e^{2\beta}K_\beta
\end{pmatrix}.
```

于是

```math
K = h^{ij}K_{ij}+K_\beta,
```

```math
K_{IJ}K^{IJ}=K_{ij}K^{ij}+K_\beta^2.
```

几何演化方程统一是

```math
\partial_t\gamma_{IJ}
=
-2N K_{IJ}+\mathcal L_{\vec N}\gamma_{IJ}.
```

展开后：

```math
\partial_t h_{ij}=-2N K_{ij}+D_iN_j+D_jN_i,
```

```math
\partial_t\beta = -N K_\beta + N^i\partial_i\beta.
```

## 2. 统一的 Einstein 约束—演化结构

若某一分支的引力方程可写成

```math
M_P^2 G_{\mu\nu}=T_{\mu\nu}^{\rm eff},
```

则它的 `ADM` 形式统一为：

### Hamilton 约束

```math
{}^{(3)}R + K^2 - K_{IJ}K^{IJ}
=
\frac{2}{M_P^2}\,\mathcal E,
```

其中

```math
\mathcal E:=n^\mu n^\nu T_{\mu\nu}^{\rm eff}.
```

### Momentum 约束

```math
D_J\left(K^J{}_I-\delta^J{}_I K\right)
=
\frac{1}{M_P^2}\,\mathcal P_I,
```

其中

```math
\mathcal P_I:=-\gamma_I{}^\mu n^\nu T_{\mu\nu}^{\rm eff}.
```

### 外曲率演化

```math
(\partial_t-\mathcal L_{\vec N})K_{IJ}
=
-D_I D_J N
+N\Big[
{}^{(3)}R_{IJ}
+K K_{IJ}
-2K_{IK}K^K{}_J
-\frac{1}{M_P^2}
\Big(
\mathcal S_{IJ}
-\frac12\gamma_{IJ}(\mathcal S-\mathcal E)
\Big)
\Big],
```

其中

```math
\mathcal S_{IJ}:=\gamma_I{}^\mu\gamma_J{}^\nu T_{\mu\nu}^{\rm eff},
\qquad
\mathcal S:=\gamma^{IJ}\mathcal S_{IJ}.
```

因此三支的区别完全压缩到：

```text
各自的有效应力能量投影 (𝓔, 𝓟_I, 𝓢_IJ) 是什么。
```

## 3. A 支：g 上的 Einstein-Klein-Gordon 系统

### 3.1 基本作用量

```math
S_A[g,\psi]
=
\frac{M_P^2}{2}\int d^4x\sqrt{-g}R[g]
+
\int d^4x\sqrt{-g}
\left(
g^{\mu\nu}\partial_\mu\psi^*\partial_\nu\psi
-m^2\psi^*\psi
\right).
```

为数值稳定起见，建议把

```math
\psi=\phi_1+i\phi_2
```

作为基本演化变量，而不是直接推进 `(\rho,S)`。

### 3.2 标量场一阶化

定义

```math
\Pi_A
:=
\frac{1}{N}
\left(
\partial_t\phi_A - N^I\partial_I\phi_A
\right),
\qquad A=1,2.
```

则标量场演化方程是

```math
\partial_t\phi_A
=
N\Pi_A + N^I\partial_I\phi_A,
```

```math
\partial_t\Pi_A
=
N\left(
D^I D_I\phi_A + K\Pi_A - m^2\phi_A
\right)
+D^I N\,\partial_I\phi_A
+N^I D_I\Pi_A.
```

### 3.3 A 支应力能量投影

以当前复标量归一化，

```math
\mathcal E_A
=
\sum_{A=1}^2
\left(
\Pi_A^2 + D^I\phi_A D_I\phi_A + m^2\phi_A^2
\right),
```

```math
\mathcal P_I^{(A)}
=
-2\sum_{A=1}^2 \Pi_A D_I\phi_A,
```

```math
\mathcal S_{IJ}^{(A)}
=
2\sum_{A=1}^2 D_I\phi_A D_J\phi_A
+\gamma_{IJ}
\sum_{A=1}^2
\left(
\Pi_A^2 - D^K\phi_A D_K\phi_A - m^2\phi_A^2
\right).
```

因此 `A` 支现在已经是一个完全标准的：

```text
四维 Einstein + 复 Klein-Gordon
```

初值问题，只是再附加了 `y` 方向 Killing 对称。

### 3.4 从 `(ρ,S)` 初值到 `(\phi_1,\phi_2,\Pi_1,\Pi_2)` 初值

若给定初始 `(\rho,S)`，则

```math
\phi_1=\sqrt{\rho}\cos S,
\qquad
\phi_2=\sqrt{\rho}\sin S.
```

若同时给定初始 `\partial_t\rho,\partial_tS`，则

```math
\Pi_1
=
\frac{1}{N}
\left[
\frac{\partial_t\rho}{2\sqrt{\rho}}\cos S
-\sqrt{\rho}\sin S\,\partial_t S
-N^I\partial_I\phi_1
\right],
```

```math
\Pi_2
=
\frac{1}{N}
\left[
\frac{\partial_t\rho}{2\sqrt{\rho}}\sin S
+\sqrt{\rho}\cos S\,\partial_t S
-N^I\partial_I\phi_2
\right].
```

## 4. B 支：g~ 上的 Einstein-Bohm 系统

### 4.1 基本作用量

```math
S_B[g^\sim,\tilde\rho,S]
=
\frac{M_P^2}{2}\int d^4x\sqrt{-g^\sim}R[g^\sim]
+
\int d^4x\sqrt{-g^\sim}\,
\tilde\rho
\left(
g^{\sim\mu\nu}\partial_\mu S\partial_\nu S
-m^2
\right).
```

### 4.2 物质方程

定义

```math
p_I:=\partial_I S,\qquad
E:=\sqrt{m^2+\gamma^{IJ}p_Ip_J},
```

并选未来指向分支

```math
\partial_t S = N^I p_I - N E.
```

再定义守恒密度

```math
n:=\sqrt{\gamma}\,\tilde\rho\,E,
```

则连续性方程写成

```math
\partial_t n
=
-\partial_I\left[
nN^I + \frac{Nn}{E}\gamma^{IJ}p_J
\right].
```

最后由

```math
\tilde\rho = \frac{n}{\sqrt{\gamma}\,E}
```

回收物质密度。

### 4.3 B 支应力能量投影

在壳约束

```math
g^{\sim\mu\nu}\partial_\mu S\partial_\nu S = m^2
```

上，物质应力张量化为

```math
T_{\mu\nu}^{(B)} = 2\tilde\rho\,u_\mu u_\nu,
\qquad
u_\mu:=\partial_\mu S.
```

因此

```math
\mathcal E_B = 2\tilde\rho\,E^2,
```

```math
\mathcal P_I^{(B)} = 2\tilde\rho\,E\,p_I,
```

```math
\mathcal S_{IJ}^{(B)} = 2\tilde\rho\,p_I p_J.
```

所以 `B` 支就是

```text
四维 Einstein 方程
+ 一组 Bohm 型一阶物质方程
```

的完整耦合系统。

## 5. C 支：g~ 上的 f(R~) 标量-张量系统

### 5.1 辅助场形式

```math
S_C[g^\sim,\Phi,\tilde\rho,S]
=
\frac{M_P^2}{2}\int d^4x\sqrt{-g^\sim}
\left[
\Phi R[g^\sim] - U(\Phi)
\right]
+
\int d^4x\sqrt{-g^\sim}\,
\tilde\rho
\left(
g^{\sim\mu\nu}\partial_\mu S\partial_\nu S
-m^2
\right).
```

物质部分与 `B` 支相同，因此 `(\tilde\rho,S)` 的演化和 `B` 支完全一致。

### 5.2 协变场方程

```math
\Phi G_{\mu\nu}[g^\sim]
=
\frac{1}{M_P^2}T_{\mu\nu}^{(B)}
+\nabla_\mu\nabla_\nu\Phi
-g^\sim_{\mu\nu}\Box\Phi
-\frac12 U(\Phi)\,g^\sim_{\mu\nu},
```

并且

```math
3\Box\Phi + 2U(\Phi)-\Phi U_\Phi(\Phi)
=
\frac{1}{M_P^2}T^{(B)},
```

其中

```math
T^{(B)} = g^{\sim\mu\nu}T_{\mu\nu}^{(B)} = 2m^2\tilde\rho.
```

### 5.3 `ADM` 一阶化

定义

```math
\Pi_\Phi
:=
\frac{1}{N}
\left(
\partial_t\Phi - N^I\partial_I\Phi
\right).
```

则

```math
\partial_t\Phi = N\Pi_\Phi + N^I\partial_I\Phi.
```

而由标量波动方程得到

```math
\partial_t\Pi_\Phi
=
N\left(
D^I D_I\Phi + K\Pi_\Phi
+\frac{1}{3M_P^2}T^{(B)}
-\frac13\bigl(2U-\Phi U_\Phi\bigr)
\right)
+D^I N\,\partial_I\Phi
+N^I D_I\Pi_\Phi.
```

数值实现时，把

```math
\nabla_\mu\nabla_\nu\Phi
```

按 `ADM` 标准分解后并入有效源项，即可把 `C` 支也写成与上节相同的

```text
约束方程 + 外曲率演化 + 标量 Φ 演化
```

结构。

## 6. 当前阶段的直接结论

经过这一整理，三支现在已经处在同一个统一语言下：

### A 支

```text
四维 Einstein-Klein-Gordon，基本度规是 g
```

### B 支

```text
四维 Einstein-Bohm，基本度规是 g~
```

### C 支

```text
四维 f(R~)-Bohm 标量-张量系统，基本度规是 g~
```

并且三者都已经具备了：

1. 统一的 `ADM` 初值问题形式  
2. 明确的几何变量  
3. 明确的物质变量  
4. 可直接数值实现的约束—演化结构

下一步不再是继续改理论定义，而是：

```text
把这些方程逐项落成数值离散代码，并先跑 A/B 两支的最小全动力学验证，再把 C 支接上。
```
