# gauge-fixed 演化、约束传播与 patch 图册草案

日期：2026-05-08

## 目的

接续 `184-完整D支Cauchy闭合快查.md`。

本笔记处理三个直接影响“能否独立给出 D 支演化”的问题：

1. 采用什么 gauge，把方程写成演化系统；
2. 约束如何传播，哪些约束必须每步硬施加；
3. \(\Delta=0,r^\perp=0,w^2=0\) 等退化区域如何换图，而不是把数值 guard 当成物理方程。

## 1. 推荐的第一版 gauge：generalized harmonic

第一版理论分析建议采用 generalized harmonic gauge，而不是 ADM。

理由：

- 当前主部分析已经在 harmonic gauge 下完成；
- reduced Einstein 方程的约束传播结构最清楚；
- \(\tilde G_{\mu\nu}= \tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}\) 可直接写成 10 个 wave-like metric 方程；
- ADM 更适合最终数值器，但会先引入 lapse/shift、Hamiltonian/momentum 约束和更多离散选择。

定义 gauge covector：

\[
F_\nu
=
\tilde g_{\nu\alpha}\tilde\nabla_\mu\tilde\nabla^\mu x^\alpha
-H_\nu(\tilde g,S,\tilde\rho),
\]

最小 harmonic 取 \(H_\nu=0\)。

reduced Einstein-like 方程可写成示意形式：

\[
\tilde G^{(F)}_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2
+C_{\mu\nu},
\]

其中 \(\tilde G^{(F)}_{\mu\nu}\) 是把 \(\tilde G_{\mu\nu}\) 加上标准 harmonic gauge-fixing 主部后的 reduced Einstein tensor。其主部为

\[
-\frac12\tilde\square \bar h_{\mu\nu}.
\]

这一步不是改变物理方程，而是把坐标自由度固定后用于 Cauchy 演化。

## 2. 完整未知量

在 massive branch 的主 patch 中取变量：

\[
\tilde g_{\mu\nu},
\qquad
S,
\qquad
\tilde\rho,
\qquad
\lambda_I.
\]

其中

\[
C_{\mu\nu}=\lambda_I E^I_{\mu\nu},
\]

\[
E=\mathrm{span}
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

并定义：

\[
u_\mu=\partial_\mu S,
\qquad
r_\mu=\tilde\nabla_\mu\ln\sqrt{\tilde\rho}.
\]

为了做严格一阶化，也可以把 \(u_\mu,r_\mu\) 当作辅助变量，并加入：

\[
\partial_{[\mu}u_{\nu]}=0,
\qquad
r_\mu-\partial_\mu\ln\sqrt{\tilde\rho}=0.
\]

这能避免在主部讨论中混淆 \(\nabla C\) 对 \(r\) 产生的二阶 \(\tilde\rho\) 导数。

## 3. 方程组

### 3.1 reduced metric 方程

\[
\boxed{
\tilde G^{(F)}_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2
+C_{\mu\nu}.
}
\]

这是 10 个 gauge-fixed metric 演化方程。

### 3.2 物质方程

\[
\boxed{
\tilde g^{\mu\nu}u_\mu u_\nu=m^2.
}
\]

\[
\boxed{
\tilde\nabla_\mu(\tilde\rho\,u^\mu)=0.
}
\]

在给定 \(\tilde g\) 后，第一式选择 \(\partial_tS\) 的 Hamilton-Jacobi 分支，第二式以守恒形式推进 \(\tilde\rho\)。

### 3.3 辅助几何应力方程

\[
\boxed{
C_{\mu\nu}=\lambda_I E^I_{\mu\nu}.
}
\]

\[
\boxed{
\tilde g^{\mu\nu}C_{\mu\nu}=0.
}
\]

\[
\boxed{
\tilde\nabla^\mu C_{\mu\nu}=0.
}
\]

trace0 是状态方程；守恒是 Bianchi 相容性在 \(C\) sector 的传播方程。

## 4. 为什么这在主 patch 内闭合

在 \(3+1d\) 中，\(\lambda_I\) 有 4 个。

trace0 给出一个代数关系，消去 1 个系数。

守恒方程对 \(\lambda_I\) 的主符号为

\[
M_{\nu I}(\xi)=\xi^\mu E^I_{\mu\nu}.
\]

此前已算出泛型 rank 为 3。

因此，在 trace0 消元后，\(\tilde\nabla^\mu C_{\mu\nu}=0\) 正好给剩余 3 个 \(\lambda\) 自由度的一阶传播方程。

对 metric 部分，\(\Pi_E^\perp\mathcal R=0\) 加 harmonic gauge 后还剩一个 \(E\)-方向主部 null mode：

\[
\bar h^{(0)}_{\mu\nu}=-w_\mu w_\nu,
\]

\[
w_\mu=(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu.
\]

trace0 对该模式的收缩为

\[
-w^2.
\]

所以只要

\[
\boxed{w^2\neq0,}
\]

trace0 可以固定这个缺口。

因此，主 patch 条件是：

\[
\boxed{
d>2,\qquad
\Delta=u^2r^2-(u\cdot r)^2\neq0,
\qquad
w^2\neq0.
}
\]

再加物质非特征条件：

\[
\boxed{
\text{HJ 分支判别式非负，且 }u^t\neq0
\text{ 或等价的时间通量不退化。}
}
\]

在这些条件下，当前系统在 principal-symbol/counting 意义上是局部闭合候选。

## 5. 约束传播

### 5.1 harmonic gauge 约束

真实物理方程要求

\[
F_\nu=0.
\]

若 reduced 方程成立，且总源守恒：

\[
\tilde\nabla^\mu
\left(
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}
\right)=0,
\]

则由 Bianchi 恒等式得到 \(F_\nu\) 的齐次波动型传播方程：

\[
\tilde\square F_\nu+\tilde R_\nu{}^\mu F_\mu
=0
\]

加上低阶项。

所以若初始切片满足

\[
F_\nu=0,
\qquad
\partial_tF_\nu=0,
\]

则 harmonic gauge 约束传播。

这里总源守恒由两部分给出：

1. 物质方程推出 \(\tilde\nabla^\mu\tilde T_{\mu\nu}=0\)；
2. \(C\) sector 方程直接施加 \(\tilde\nabla^\mu C_{\mu\nu}=0\)。

### 5.2 Einstein constraint

在 full Einstein-like 写法中：

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}.
\]

四个 normal projection：

\[
n^\mu
\left(
\tilde G_{\mu\nu}
-\tilde T_{\mu\nu}/M_P^2
-C_{\mu\nu}
\right)=0
\]

是初值约束。

在 harmonic reduced 系统中，只要 gauge 约束和总源守恒传播，这四个 Einstein constraints 也按 GR 的标准机制传播。

因此它们不是每一步可任意重设的数值规则；它们必须在初值求解器中满足，并在数值器中用 constraint-preserving 方法监控。

### 5.3 trace0

trace0 不是由 Bianchi 推出的约束，而是状态方程：

\[
\tilde g^{\mu\nu}C_{\mu\nu}=0.
\]

它必须作为方程系统的一部分持续成立。

最干净的实现方式是：直接用 trace0 消元 \(\lambda_I\)，不要把 trace0 当作只在初始切片检查一次的约束。

### 5.4 integrability constraints

由于

\[
u_\mu=\partial_\mu S,
\]

需要

\[
\partial_{[\mu}u_{\nu]}=0.
\]

若直接推进 \(S\)，该条件自动成立。

若改成推进 \(u_\mu\)，必须把 curl-free 条件作为约束传播。

同理，若把 \(r_\mu\) 当独立变量，需要保证

\[
r_\mu=\partial_\mu\ln\sqrt{\tilde\rho}.
\]

## 6. patch 图册

### 6.1 主图 \(P_{ur}\)

条件：

\[
d>2,\qquad
\Delta\neq0,\qquad
w^2\neq0.
\]

使用四方向基底：

\[
E_{ur}=
\mathrm{span}\{
\tilde g,uu,rr,ur
\}.
\]

使用普通 trace0 闭合。

这是当前 proposal 的主工作区。

### 6.2 \(r\parallel u\) 或 \(r^\perp=0\) 图

若

\[
r_\mu=\alpha u_\mu,
\]

则

\[
rr=\alpha^2uu,
\qquad
ur=\alpha uu.
\]

四方向基底降为：

\[
E_u=\mathrm{span}\{\tilde g,uu\}.
\]

trace0 后只剩一个无迹方向：

\[
H^{(u)}_{\mu\nu}
=
u_\mu u_\nu-\frac{u^2}{d}\tilde g_{\mu\nu}.
\]

因此在该图内可写：

\[
C_{\mu\nu}=B\,H^{(u)}_{\mu\nu}.
\]

注意：守恒 \(\nabla^\mu(BH^{(u)}_{\mu\nu})=0\) 可能对一个标量 \(B\) 过定，因此该图常常会迫使：

- \(B=0\)，即回到 GR/branch 区；
- 或者要求两侧 bulk 几何满足额外边界匹配；
- 或者说明需要扩大基底，例如加入 \(\nabla_{(\mu}u_{\nu)}\)、\(\nabla_{(\mu}r_{\nu)}\)。

所以 \(E_u\) 是合法图册，但不是保证处处可行的“万能降维理论”。

### 6.3 \(w^2=0\) 主符号图

\(w^2=0\) 不是基底降秩，而是普通 trace0 看不见某个主部模式。

处理方式有三类：

1. 选择 Cauchy 面和 gauge，使物理解的相关主方向避开 \(w^2=0\)；
2. 在该图内切换到更强的 scalar closure，例如 \(\mathsf q_E\)-trace；
3. 把该区域作为特征/边界处理，要求从相邻非退化图延拓 \(C_{\mu\nu}\)。

第二种会改变状态方程，因此应标记为 proposal v1.2，而不是偷偷说仍是同一个 v1。

### 6.4 \(Q=0\) 或弱量子 branch 图

若选择 GR branch，则要求：

\[
Q=0\quad\Rightarrow\quad C=0
\]

或弱版本：

\[
\|C\|\ll
\|\tilde T/M_P^2\|+\|\tilde G\|.
\]

若不用硬写 \(C=\chi\hat C\)，则必须在 \(Q=0\) 界面施加无通量/正则性条件：

\[
n^\mu C_{\mu\nu}=0
\]

或有界延拓条件。

## 7. 图册切换规则

关键原则：

\[
\boxed{
C_{\mu\nu}\text{ 是主变量，}\lambda_I\text{ 只是局部坐标。}
}
\]

因此切换图册时要求连续的是 \(C_{\mu\nu}\)，不是各图中的系数。

从图 \(P_A\) 切到 \(P_B\) 时：

1. 用旧图重建张量 \(C_{\mu\nu}\)；
2. 在新图基底中重新投影/最小二乘得到新系数；
3. 若界面为真实薄层，需要满足
   \[
   n^\mu[C_{\mu\nu}]=0
   \]
   或对应的积分 jump law；
4. 若新图无法表达旧 \(C\)，则该点不是坐标奇点，而是当前最小基底失败，需要扩大 \(E\) 或进入边界/缺陷描述。

数值实现中的阈值，例如 \(\Delta_{\rm rel}<10^{-3}\)，只是 chart selection tolerance，不是物理常数。

## 8. 独立 D 支演化器骨架

第一版 standalone solver 应按以下流程写，而不是继续从 A 快照反推。

### 8.1 初值

在初始 Cauchy 面 \(\Sigma_0\) 给：

\[
\gamma_{ij},\quad K_{ij},
\quad S,\quad \tilde\rho,
\quad \lambda_I
\]

或 harmonic 坐标下等价的

\[
\tilde g_{\mu\nu},\quad \partial_t\tilde g_{\mu\nu}.
\]

初值必须满足：

1. HJ 壳条件；
2. Einstein constraints with source \(\tilde T/M_P^2+C\)；
3. trace0；
4. \(C\in E\)；
5. normal/tangential split of \(\nabla^\mu C_{\mu\nu}=0\)；
6. harmonic gauge constraints \(F_\nu=0,\partial_tF_\nu=0\)；
7. patch 条件或合法图册标签。

### 8.2 单步推进

每个时间步：

1. 由 \(S,\tilde\rho,\tilde g\) 计算 \(u,r,E\)；
2. 用 trace0 消元 \(\lambda_I\)；
3. 用 \(\nabla C=0\) 推进剩余 \(\lambda\)；
4. 用 HJ 方程推进 \(S\)，用连续性方程守恒推进 \(\tilde\rho\)；
5. 用 reduced Einstein 方程推进 \(\tilde g\)；
6. 监控 harmonic constraint、Einstein constraint、trace0、\(\nabla C\)、mass shell；
7. 若进入退化区域，按 patch 图册切换；
8. 必要时做 constraint projection，但投影必须保持物理方程，不可用软 penalty 假装成立。

## 9. 当前可写入文章的表述

可以写：

> In the non-degenerate massive patch, the matter Hamilton-Jacobi and continuity equations, together with the trace-free conserved \(E\)-sector correction, form a locally closed equation-first Cauchy candidate after generalized harmonic gauge fixing.

中文：

> 在非退化 massive patch 内，物质 HJ 方程、连续性方程与无迹守恒的 \(E\)-sector 修正张量结合后，在 generalized harmonic gauge 固定下形成一个局部闭合的 equation-first Cauchy 候选系统。

不能写：

> 我们已经证明该系统全局适定。

也不能写：

> 当前数值 hard-ALM 已经是独立 D 支演化器。

更准确的结论是：

\[
\boxed{
\text{当前方程已经具备独立 D 支演化理论的主部闭合框架；下一步是约束传播证明和 standalone solver。}
}
