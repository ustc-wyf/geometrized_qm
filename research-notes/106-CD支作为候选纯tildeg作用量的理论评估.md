# CD 支作为候选纯 `tilde g` 作用量的理论评估

时间：2026-05-01

## 问题

此前的 native 场方程残差检验只能说明：

```text
C/D 支没有被平直时空量子力学参考解立即排除。
```

它不能说明：

```text
C/D 支已经给出实验兼容的真实动力学预测。
```

真正要判断 C/D 是否被实验否定，下一步必须求解它们自己的全作用量演化，并把演化出的 `rho,S,tilde g` 可观测量与普通量子力学实验对比。

本笔记先只回答理论层面的问题：C/D 为什么有资格作为候选作用量，或者会因什么理由被理论上排除。

## 候选作用量

当前 C/D 都被理解为 native `tilde g` 表象里的 metric `f(tilde R)` 理论：

\[
S_{C/D}
=
\int d^4x\sqrt{-\tilde g}
\left[
\frac{M_P^2}{2} f(\tilde R)
+ \rho
\left(
\tilde g^{\mu\nu}\partial_\mu S\partial_\nu S-m^2
\right)
\right].
\]

其中

\[
f_C(\tilde R)
=
\frac{\tilde R}{\sqrt{1+(\ell^2\tilde R)^2}},
\qquad
f_D(\tilde R)
=
\frac{1}{\ell^2}\tanh(\ell^2\tilde R).
\]

令

\[
y=\ell^2\tilde R.
\]

则

\[
f_{C,R}=(1+y^2)^{-3/2},
\qquad
f_{D,R}=\operatorname{sech}^2 y.
\]

两者都满足

\[
f(0)=0,\qquad f_R(0)=1,\qquad f_R>0.
\]

## native 场方程

metric `f(R)` 的场方程是

\[
M_P^2
\left[
f_R\tilde R_{\mu\nu}
-\frac12 f\tilde g_{\mu\nu}
-
\left(
\tilde\nabla_\mu\tilde\nabla_\nu
-\tilde g_{\mu\nu}\tilde\Box
\right)f_R
\right]
=
\tilde T_{\mu\nu}.
\]

四维迹方程是

\[
M_P^2
\left(
f_R\tilde R-2f+3\tilde\Box f_R
\right)
=
\tilde T.
\]

物质项仍然最小耦合在 `tilde g` 上，因此对 `S,rho` 变分给出

\[
\tilde\nabla_\mu(\rho \tilde u^\mu)=0,
\qquad
\tilde g^{\mu\nu}u_\mu u_\nu=m^2.
\]

由于 \(u_\mu=\partial_\mu S\)，质量壳方程直接推出几何光学轨道沿 `tilde g` 测地线。因此 C/D 的一个重要优点是：引力作用量完全由 `tilde g` 衍生，不会在 `S` 方程里直接加入额外的 `u_\mu` 力项。

## 支持 C/D 继续作为候选的理由

### 1. 它们是自洽的协变作用量

C/D 不是 ad hoc 的场方程，而是明确来自 diffeomorphism invariant 的局域作用量。自由引力度规是 `tilde g`，物质也最小耦合在 `tilde g` 上，因此变分对象清楚。

这点比“把某个阻尼因子直接乘在 Einstein 方程左边”更干净，因为后者未必能从局域作用量推出，也未必自动满足 Bianchi 恒等式对应的守恒结构。

### 2. 弱曲率极限自动回到 Einstein-Hilbert

当

\[
|\ell^2\tilde R|\ll 1
\]

时，

\[
f_C(\tilde R)
=
\tilde R-\frac12\ell^4\tilde R^3+O(\tilde R^5),
\]

\[
f_D(\tilde R)
=
\tilde R-\frac13\ell^4\tilde R^3+O(\tilde R^5).
\]

因此二者都没有宇宙常数项，也没有改变线性化 Einstein-Hilbert 主项。若普通引力实验所在区域满足 `|ell^2 R_obs| << 1`，C/D 在这些区域可以自然接近 GR。

### 3. 强 `tilde R` 区域会被饱和

当

\[
|\ell^2\tilde R|\gg 1
\]

时，

\[
f_C\to \frac{\operatorname{sgn}\tilde R}{\ell^2},
\qquad
f_D\to \frac{\operatorname{sgn}\tilde R}{\ell^2}.
\]

同时

\[
f_{C,R}\sim |\ell^2\tilde R|^{-3},
\qquad
f_{D,R}\sim 4e^{-2|\ell^2\tilde R|}.
\]

所以在强 `tilde R` 体区，Einstein-like 响应 \(f_R\tilde R_{\mu\nu}\) 被压低。这个机制正好对应用户的物理动机：即使由量子势诱导的 `tilde R` 很大，也不希望它在普通平直量子力学实验里产生巨大的真实引力反作用。

### 4. 它们解释了为什么 B 支被快速否定

B 支是 \(f(\tilde R)=\tilde R\)，即 \(f_R=1\)。它没有强曲率饱和机制。

因此平直量子力学参考波函数生成的 `tilde g[rho,S]` 一旦具有很大的 `tilde R`，B 支 native Einstein 方程就必须直接响应这些曲率结构。此前残差检验中 B 支残差很大，理论上并不意外。

C/D 的残差显著下降，正是因为它们把大 `tilde R` 体区从 Einstein 响应中屏蔽掉了。

### 5. D 比 C 更贴近当前目标

C 的压制是幂律，D 的压制是指数型。

因此 D 在强 `tilde R` 体区更快进入饱和，更符合“尽量不偏离平直时空量子力学”的目标。此前数值残差中 D 的 bulk 代表性指标优于 C，与这个理论预期一致。

代价是：D 的过渡层通常更薄、更 stiff。

## 可能否定 C/D 的理论理由

### 1. 过渡层导数项不是可忽略修正

C/D 的危险项不是饱和体区的代数项，而是

\[
\left(
\tilde\nabla_\mu\tilde\nabla_\nu
-\tilde g_{\mu\nu}\tilde\Box
\right)f_R.
\]

当

\[
|\ell^2\tilde R|=O(1)
\]

时，链式法则会给出大的 \(f_{RR}\)、\(f_{RRR}\) 系数。即使 \(f\) 本身被 \(1/\ell^2\) 压低，过渡层里的导数项仍可能形成薄层尖峰。

在四维迹方程中，薄层法向积分近似给出类似边界条件：

\[
3M_P^2[\partial_n f_R]^+_-
\approx
\int_{\rm layer} \tilde T\,dn
-
M_P^2\int_{\rm layer}(f_R\tilde R-2f)\,dn.
\]

因此如果 \(f_R\) 的法向导数在过渡层两侧产生有限跳变，过渡层会表现得像有效 surface stress，而不是普通小扰动。

这件事若在实验可观测区域产生有限相位偏移、散射或轨迹偏折，就会否定 C/D。

### 2. \(f_R\to0\) 会带来退化或强耦合风险

在 metric `f(R)` 中，\(f_R\) 是 Einstein-like 张量项前的有效系数，也对应标量-张量重写里的

\[
\phi=f_R.
\]

C/D 在强饱和区都有

\[
\phi=f_R\to0.
\]

这可以解释为“引力响应被关小”，但也意味着场方程的主部可能退化。若在有限体积、有限密度的区域里 \(\phi\) 过小，理论可能出现：

- 度规方程对某些分量约束不足；
- 标量-张量形式强耦合；
- Cauchy 问题不适定；
- 小残差对应大解偏移。

因此“残差小”必须搭配线性响应/稳定性分析才有物理意义。

### 3. 标准 metric `f(R)` 稳定性条件对 C/D 不友好

标准 metric `f(R)` 作为普适引力理论时，常用健康条件包括

\[
f_R>0,
\qquad
f_{RR}\ge0
\]

至少在所关心的背景附近成立。前者避免有效 Planck 质量变号；后者与 scalaron 稳定性有关。

C/D 虽然满足 \(f_R>0\)，但

\[
f_{C,RR}
=
-3\ell^2 y(1+y^2)^{-5/2},
\]

\[
f_{D,RR}
=
-2\ell^2\operatorname{sech}^2y\,\tanh y.
\]

所以在 \(\tilde R>0\) 区域二者都有 \(f_{RR}<0\)。

四维背景上的 scalaron 质量量级为

\[
m_s^2
\sim
\frac{f_R-\tilde R f_{RR}}{3f_{RR}}.
\]

因此在正曲率背景上 C/D 会给出负的 \(m_s^2\) 倾向。这是一个很强的理论红旗。

更关键的是：只要我们要求 \(f_R\) 随正 \(\tilde R\) 从 1 单调下降到 0，某些正曲率区域里就不可避免有 \(f_{RR}<0\)。所以“正曲率饱和屏蔽”和“标准 metric f(R) 的 \(f_{RR}\ge0\) 健康条件”之间有结构性张力。

这不一定立刻排除 C/D，因为我们可能把它们看成受限有效理论，或未来改用 Palatini/更高结构理论。但若坚持把 C/D 当作普适 metric `f(R)` 引力理论，这一点是目前最强的否定理由。

### 4. `ell` 必须满足真实实验尺度分离

C/D 要同时做到两件事：

\[
|\ell^2 R_{\rm ordinary}|\ll1
\]

保证普通 GR/弱引力实验不变；

\[
|\ell^2 \tilde R_{\rm quantum}|\gg1
\]

保证量子势诱导的大 `tilde R` 被屏蔽。

也就是说必须存在尺度窗口：

\[
|R_{\rm ordinary}|\ll \ell^{-2}\ll |\tilde R_{\rm quantum}|.
\]

如果真实实验中不存在这样一个窗口，C/D 就会被否定。当前的 `ell=100` 只是无量纲数值实验参数，还不是物理标定。

### 5. 参考残差小不等于演化预测小

native 残差检验只说明平直量子力学参考解是候选场方程的近似解。

但从残差到演化误差还需要一个稳定性估计：

\[
\delta\Psi
\sim
\mathcal L^{-1}\mathcal E.
\]

若线性化算子 \(\mathcal L\) 在 \(f_R\to0\) 或过渡层附近条件数很大，则很小的残差也可能导致较大的真实演化偏离。

因此实验判断最终必须比较全作用量演化出的可观测量，而不能只停在残差图。

### 6. 单 KG 标量场还不能保证普适实验兼容

当前物质项只检验了一个 KG 标量场。若 `tilde g` 被宣称为所有物质共同耦合的物理度规，则还要说明标准模型场、电磁场、实验仪器、探测器读数如何耦合。

若只有这个 KG 标量场耦合到 `tilde g`，则理论会违反普适等效原理或至少引入强烈的物种依赖。这个问题不在当前数值残差检验范围内。

## 当前理论判断

### B 支

B 支可以基本排除为当前目标下的候选。

理由不是“数学不自洽”，而是它缺乏屏蔽机制，因此无法让平直量子力学参考构型成为 native `tilde g` Einstein 方程的近似解。

### C 支

C 支可以保留为“幂律饱和型对照候选”。

它的优点是过渡比 D 更缓，可能数值上稍温和；缺点是强曲率体区压制不如 D，残差和实验偏离可能更大。

### D 支

D 支是当前最值得优先推进的候选。

它最符合用户提出的物理图景：

```text
普通弱曲率区回到 EH；
量子势诱导的大 tilde R 区进入饱和；
bulk 中尽量贴近平直量子力学。
```

但 D 的理论风险也最清楚：

```text
更薄的过渡层；
更强的局部 stiffness；
f_R -> 0 的退化；
正曲率区 f_RR < 0 的 metric f(R) 稳定性红旗。
```

## 结论

当前最准确的表述是：

\[
\boxed{
C/D 支是有明确理论动机的有效候选作用量，
不是任意拼凑的方程。
}
\]

但也必须同时承认：

\[
\boxed{
它们还没有达到“理论上可接受且实验未排除”的层级。
}
\]

尤其是 \(f_{RR}<0\)、\(f_R\to0\)、过渡层边界条件这三点，是继续推进全作用量演化前必须牢记的理论风险。

因此下一步若继续研究 CD 支，最合理的顺序是：

1. 先把 C/D 的 scalar-tensor/`chi` 形式写成清楚的 Cauchy 问题，检查 principal part 与线性稳定性；
2. 对过渡层做积分形式的 jump/boundary condition 分析，判断薄层是否会产生有限实验效应；
3. 再用全作用量演化比较可观测量，而不是继续只看参考残差。

这也解释了为什么后续数值工作应优先 D 支，但不能把 D 支的残差小直接当成理论成功。
