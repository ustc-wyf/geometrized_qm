# D 支 harmonic 闭合定理草案与当前例子可演化性

日期：2026-05-08

## 0. 本轮要解决的问题

用户指出上一轮只是总结，没有继续处理以下问题：

1. 选定 harmonic gauge 或 ADM gauge；
2. 写出完整约束/演化分裂；
3. 证明约束传播；
4. 处理 \(\Delta=0,r^\perp=0,w^2=0\) 的 patch transition；
5. 判断当前例子能否独立求解 D 支演化。

本笔记给出一个明确版本。

结论先写：

\[
\boxed{
\text{第一版理论采用 generalized harmonic gauge，不采用 ADM 作为主表述。}
}
\]

\[
\boxed{
\text{在 massive、非退化主 patch 内，当前 equation-first D 支已形成局部 Cauchy 闭合候选。}
}
\]

\[
\boxed{
\text{当前高斯干涉例子的 core10 主支撑区几乎全落在主 patch 内，可作为独立 D 支演化器的局部初值测试场。}
}
\]

但：

\[
\boxed{
\text{这还不是全局适定性定理，也不是已经完成的生产级非线性 D 支演化器。}
}
\]

## 1. Gauge 选择：选 generalized harmonic，不选 ADM

第一版使用 generalized harmonic gauge：

\[
F_\nu
=\tilde g_{\nu\alpha}\tilde\nabla_\mu\tilde\nabla^\mu x^\alpha
-H_\nu(\tilde g,S,\tilde\rho).
\]

最小版本取 \(H_\nu=0\)。

理由：

- harmonic gauge 下 Einstein 方程的主部直接是波方程；
- Bianchi 恒等式给 harmonic 约束传播的标准证明；
- 前面投影主符号分析已经使用 harmonic gauge；
- ADM 适合最终数值实现，但会额外引入 lapse/shift、Hamiltonian/momentum split 和坐标选择，不适合作为当前理论闭合的第一版。

reduced 方程写成

\[
\boxed{
\tilde G^{(F)}_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}.
}
\]

其中 \(\tilde G^{(F)}\) 表示加上 harmonic gauge-fixing 后的 reduced Einstein tensor。其线性化主部为

\[
-\frac12\tilde\square\bar h_{\mu\nu}.
\]

这一步是坐标规范固定，不是修改物理方程。物理解还必须满足 \(F_\nu=0\)。

## 2. 完整未知量与方程

主 massive patch 中取未知量：

\[
\tilde g_{\mu\nu},\quad S,\quad \tilde\rho,\quad C_{\mu\nu}.
\]

在局部图册中把 \(C\) 展开为

\[
C_{\mu\nu}=\lambda_I E^I_{\mu\nu},
\]

\[
E=
\mathrm{span}
\left\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\right\},
\]

\[
u_\mu=\partial_\mu S,\qquad
r_\mu=\tilde\nabla_\mu\ln\sqrt{\tilde\rho}.
\]

方程组为

\[
\boxed{
\tilde G^{(F)}_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}.
}
\]

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

\[
\boxed{
C_{\mu\nu}\in E.
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

这里 \(\tilde T_{\mu\nu}\) 的归一化按当前物质作用量约定写；若代码里有整体因子 2，则整体吸收入 \(\tilde T\) 定义，不改变闭合结构。

## 3. 演化/约束分裂

令切片法向协向量为 \(\xi_\mu\)，实验室时间时 \(\xi_\mu=dt_\mu\)。

### 3.1 metric sector

演化方程：

\[
\tilde G^{(F)}_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}
\]

给出 10 个二阶波型方程，用于推进 \(\tilde g_{\mu\nu}\)。

约束：

\[
F_\nu=0,
\]

以及 full Einstein 方程相对 reduced 方程的 normal projection 初值约束。实际实现中等价于要求初始切片上

\[
F_\nu|_\Sigma=0,\qquad
\partial_tF_\nu|_\Sigma=0.
\]

### 3.2 matter sector

Hamilton-Jacobi 方程写成关于 \(S_t\) 的二次方程：

\[
\tilde g^{tt}S_t^2
+2\tilde g^{ti}S_tS_i
+\tilde g^{ij}S_iS_j
-m^2=0.
\]

若判别式非负，且选定分支后

\[
u^t=\tilde g^{t\mu}u_\mu\neq0,
\]

则可用该式推进 \(S\)。

连续性方程写成守恒形式：

\[
\partial_t(\sqrt{|\tilde g|}\tilde\rho u^t)
+\partial_i(\sqrt{|\tilde g|}\tilde\rho u^i)=0.
\]

这推进 \(\tilde\rho\)，并保持正密度需数值格式额外保证。

### 3.3 \(C\)-sector

在主图中

\[
\lambda_I=(\lambda_g,\lambda_{uu},\lambda_{rr},\lambda_{ur}).
\]

trace0 给出代数关系：

\[
\tau_I\lambda_I=0,
\]

\[
\tau_I=(d,u^2,r^2,u\cdot r).
\]

因此令

\[
\lambda_I=N_I{}^A a_A,\qquad A=1,2,3,
\]

其中 \(N_I{}^A\) 张成 \(\tau_I\lambda_I=0\) 的 nullspace。

守恒方程主部为

\[
\xi^\mu E^I_{\mu\nu}\partial_\xi\lambda_I.
\]

代入 trace0 后：

\[
M_{\nu A}(\xi)\partial_\xi a_A,\qquad
M_{\nu A}(\xi)=\xi^\mu E^I_{\mu\nu}N_I{}^A.
\]

主 patch 要求

\[
\boxed{
\mathrm{rank}\,M(\xi)=3.
}
\]

在 \(3+1d\) 中 \(\nu=0,1,2,3\)，所以这是 4 个方程对 3 个时间导数：其中 3 个独立组合是 \(a_A\) 的演化方程，剩下 1 个组合是 \(C\)-sector constraint，必须在初值中施加并在数值中监控/投影。

在当前 2+1 \(t,x,z\) 约化诊断中，\(\nu=0,1,2\)，rank=3 则直接给 3 个演化组合。

## 4. 约束传播证明

### 4.1 matter stress conservation

对 dust/HJ 型物质源

\[
\tilde T_{\mu\nu}=\tilde\rho u_\mu u_\nu
\]

有

\[
\tilde\nabla^\mu\tilde T_{\mu\nu}
=
u_\nu\tilde\nabla^\mu(\tilde\rho u_\mu)
+\tilde\rho u^\mu\tilde\nabla_\mu u_\nu.
\]

连续性方程给第一项为零。

又因为 \(u_\mu=\partial_\mu S\)，无挠 Levi-Civita 联络下

\[
u^\mu\tilde\nabla_\mu u_\nu
=
u^\mu\tilde\nabla_\nu u_\mu
=
\frac12\tilde\nabla_\nu(u^\mu u_\mu)
=
\frac12\tilde\nabla_\nu m^2
=0.
\]

因此

\[
\boxed{
\tilde\nabla^\mu\tilde T_{\mu\nu}=0.
}
\]

这也解释了为什么 HJ + continuity 能推出测地线运动。

### 4.2 total source conservation

我们把

\[
\tilde\nabla^\mu C_{\mu\nu}=0
\]

作为 \(C\)-sector 方程。因此

\[
\boxed{
\tilde\nabla^\mu
\left(
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}
\right)=0.
}
\]

### 4.3 harmonic constraint propagation

reduced 方程可写成

\[
\tilde G_{\mu\nu}
-\tilde\nabla_{(\mu}F_{\nu)}
+\frac12\tilde g_{\mu\nu}\tilde\nabla_\alpha F^\alpha
=
\tilde T_{\mu\nu}/M_P^2+C_{\mu\nu}.
\]

对两边取 \(\tilde\nabla^\mu\)。由 Bianchi 恒等式

\[
\tilde\nabla^\mu \tilde G_{\mu\nu}=0
\]

和总源守恒，得到

\[
\tilde\nabla^\mu
\left(
\tilde\nabla_{(\mu}F_{\nu)}
-\frac12\tilde g_{\mu\nu}\tilde\nabla_\alpha F^\alpha
\right)=0.
\]

标准化简给

\[
\boxed{
\tilde\square F_\nu+\tilde R_\nu{}^\mu F_\mu=0
}
\]

至多差一个整体因子，取决于 reduced tensor 的符号约定。

所以若初值满足

\[
F_\nu|_\Sigma=0,\qquad
\partial_tF_\nu|_\Sigma=0,
\]

则唯一性给

\[
F_\nu=0
\]

在依赖域内传播。

因此 reduced 解也是 full equation 解。

## 5. 主 patch 条件

定义

\[
\Delta=u^2r^2-(u\cdot r)^2.
\]

定义 metric 主符号缺口方向

\[
w_\mu=(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu.
\]

ordinary trace0 对该缺口的收缩为

\[
-w^2.
\]

因此主 patch 条件为

\[
\boxed{
d>2,\qquad
\Delta\neq0,\qquad
w^2\neq0,\qquad
u^t\neq0,\qquad
\mathrm{rank}\,M(\xi)=3.
}
\]

这些条件含义：

- \(d>2\)：严格 1+1d 需要单独理论，不能直接套四方向 \(E\) 图册；
- \(\Delta\neq0\)：\(u,r\) 张成的张量基底不降秩；
- \(w^2\neq0\)：trace0 确实关闭 metric principal gap；
- \(u^t\neq0\)：当前时间切片不是物质 HJ 的坏特征面；
- \(\mathrm{rank}\,M=3\)：trace0 后的 \(C\)-sector 守恒能推进 3 个独立 \(\lambda\) 自由度。

## 6. patch transition 规则

### 6.1 基本原则

主变量是

\[
\boxed{C_{\mu\nu}}
\]

而不是局部系数 \(\lambda_I\)。

图册切换时要求连续的是 \(C_{\mu\nu}\)。\(\lambda_I\) 可以发散或不连续，因为它只是局部坐标。

若出现真实界面，还要求法向通量匹配：

\[
\boxed{
n^\mu[C_{\mu\nu}]=0.
}
\]

### 6.2 \(\Delta=0\) 或 \(r^\perp=0\)

这是基底降秩。

若 \(r_\mu\parallel u_\mu\)，改用低维图

\[
E_u=\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu\}.
\]

trace0 后

\[
C_{\mu\nu}
=
B\left(
u_\mu u_\nu-\frac{u^2}{d}\tilde g_{\mu\nu}
\right).
\]

守恒

\[
\tilde\nabla^\mu C_{\mu\nu}=0
\]

可能过定 \(B\)。如果没有非零 \(B\) 解，则该 patch 的正则选择是

\[
B=0,\qquad C=0,
\]

或者必须扩大 basis。不能把 \(E_u\) 全局替代主图。

### 6.3 \(w^2=0\)

这不是 basis 降秩，而是 trace0 对 metric 主符号缺口失效。

有三种合法处理：

1. 换 Cauchy covector \(\xi\)，即换局部时间函数或 harmonic source，使 \(w^2\neq0\)；
2. 引入不同 scalar closure，例如 \(\mathsf q_E\)-trace，但这应标记为 proposal v1.2；
3. 把 \(w^2=0\) 面作为 characteristic/boundary surface 处理，施加界面匹配。

当前 proposal v1 不应假装 \(w^2=0\) 区域已由 ordinary trace0 自动解决。

### 6.4 \(Q\to0\) branch

\(Q\to0\Rightarrow C\to0\) 不是 trace0 + conservation 的自动推论。

若要保留 GR branch，应写成 branch 正则性：

\[
C_{\mu\nu}=\chi(\mathcal q)\hat C_{\mu\nu},
\]

\[
\chi(0)=\chi'(0)=0,\qquad
\hat C_{\mu\nu}\ \text{有界}.
\]

若 \(\chi\) 只有一阶零点，则必须补

\[
n^\mu\hat C_{\mu\nu}=0
\]

作为 \(\mathcal q=0\) 面上的无通量条件。

这属于 branch selection，不是 Cauchy 主部闭合本身。

## 7. 当前高斯干涉例子的 Cauchy patch 诊断

新增脚本：

```bash
python3 kg_examples/diagnose_d_cauchy_patch_atlas.py \
  --output visualizations/d_cauchy_patch_atlas_n384_core10 \
  --full-resolution 384 \
  --region core10 \
  --taus=-3.5,0,3.5
```

脚本检查：

- \(\Delta\) 是否非退化；
- \(w^2\) 是否非零；
- \(u^t\) 是否非零，即 HJ 是否可用当前时间推进；
- trace0 消元后的 \(\lambda\)-sector 时间主符号是否 rank=3；
- 并画出 patch label 图。

结果：

| slice | core10 点数 | 主 patch 点数 | 主 patch 占比 | 异常 |
|---|---:|---:|---:|---|
| \(\tau=-3.5\) | 2550 | 2550 | 100% | 无 |
| \(\tau=0\) | 1024 | 1023 | 99.902% | 1 个 \(\lambda\) rank bad 点 |
| \(\tau=+3.5\) | 3589 | 3589 | 100% | 无 |

图像：

- `visualizations/d_cauchy_patch_atlas_n384_core10/d_cauchy_patch_atlas_taum3p5.png`
- `visualizations/d_cauchy_patch_atlas_n384_core10/d_cauchy_patch_atlas_tau0.png`
- `visualizations/d_cauchy_patch_atlas_n384_core10/d_cauchy_patch_atlas_taup3p5.png`

解释：

在当前高斯干涉例子的主要物理支撑区，局部 Cauchy 前提基本成立。此前很多“界面/patch”焦虑主要属于低密度区、支撑边缘或数值代表元问题，不是 core10 主物理区的主障碍。

这说明：

\[
\boxed{
\text{当前例子可以作为独立 D 支演化器的初值测试场。}
}
\]

但这还不等于：

\[
\boxed{
\text{我们已经有完整独立 D 支演化结果。}
}
\]

因为还没有实现真正的非线性 generalized-harmonic D 支推进器。

## 8. 独立 D 支演化器的最小实现规格

若下一步写真正演化器，应按下面的最小系统实现，不再退回 A 支快照驱动：

### 8.1 初值

给定

\[
\tilde g_{\mu\nu}|_\Sigma,\quad
\partial_t\tilde g_{\mu\nu}|_\Sigma,\quad
S|_\Sigma,\quad
\tilde\rho|_\Sigma,\quad
a_A|_\Sigma.
\]

必须满足：

1. HJ mass shell；
2. \(\tilde\rho>0\)；
3. \(C=N^A a_A E_A\in E\)；
4. trace0；
5. \(C\)-sector constraint；
6. Einstein normal constraints；
7. harmonic constraints \(F_\nu=0,\partial_tF_\nu=0\)；
8. patch 条件。

### 8.2 每步推进

1. 用 HJ 分支求 \(S_t\)，推进 \(S\)；
2. 用守恒通量推进 \(\tilde\rho\)；
3. 用 \(\nabla C=0\) 的 3 个独立演化组合推进 \(a_A\)；
4. 用 reduced Einstein wave 方程推进 \(\tilde g_{\mu\nu}\)；
5. 对 trace0、harmonic、Einstein normal、\(C\)-constraint 做监控或约束投影；
6. 若进入 \(\Delta=0,w^2=0,u^t=0\) 等区域，按 patch 图册换图。

## 9. 本轮判断

五个问题的状态：

1. gauge：已选 generalized harmonic；
2. 完整分裂：已写出 metric/matter/\(C\)-sector 的演化与约束；
3. 约束传播：已由 matter conservation + \(C\) conservation + Bianchi 推出 harmonic 约束传播；
4. patch transition：已明确 \(\Delta=0,r^\perp=0,w^2=0,Q\to0\) 的处理原则；
5. 当前例子可演化性：core10 主物理区通过局部 Cauchy patch 诊断。

所以现在的真实下一步不是继续怀疑“方程是否闭合”，而是：

\[
\boxed{
\text{实现真正的 generalized-harmonic D 支独立演化器，并用当前三切片做初值/短时推进测试。}
}
\]

