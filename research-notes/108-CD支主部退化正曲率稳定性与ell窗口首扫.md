# CD 支主部退化、正曲率稳定性与 `ell` 窗口首扫

时间：2026-05-01

## 本轮任务安排

本轮按用户要求，明确执行以下任务：

1. 解释正 \(\tilde R\) 区域 \(f_{RR}<0\) 到底是灾难性、可能发散，还是可能安全；
2. 解释 `Cauchy 适定` 的含义；
3. 回答当 \(f_R\to0\) 时，是否可以把 \(f\to1/\ell^2\) 的代数项看成新的主部；
4. 实现并运行 `ell-window scan`；
5. 把理论分析和数值诊断写入项目笔记与长记忆。

## 对上一轮计划没有执行完的解释

上一轮提出了：

```text
scalar-tensor/chi Cauchy 主部分析
过渡层积分边界条件分析
线性响应稳定性估计
```

但实际完成的是概念解释与 `ell` 参数路线，尚未执行这三个技术分析，也没有展示对应计算或判据。

更准确的说法应该是：

```text
上一轮完成了这些技术分析的前置概念整理；
本轮开始把它们落到可检查的理论判据和 ell 扫描脚本上。
```

## 1. 正 \(\tilde R\) 区域 \(f_{RR}<0\) 是不是灾难

### 精确符号

对 C/D 支，令

\[
y=\ell^2\tilde R.
\]

C 支：

\[
f_{C,R}=(1+y^2)^{-3/2}>0,
\qquad
f_{C,RR}=-3\ell^2y(1+y^2)^{-5/2}.
\]

D 支：

\[
f_{D,R}=\operatorname{sech}^2y>0,
\qquad
f_{D,RR}=-2\ell^2\operatorname{sech}^2y\tanh y.
\]

因此：

\[
\tilde R>0 \Rightarrow f_{RR}<0,
\qquad
\tilde R<0 \Rightarrow f_{RR}>0.
\]

### scalaron 质量符号

metric `f(R)` 的 scalaron 质量平方近似为

\[
m_s^2
=
\frac{f_R-\tilde R f_{RR}}{3f_{RR}}.
\]

在 \(\tilde R>0\) 且 \(f_{RR}<0\) 时，

\[
f_R-\tilde R f_{RR}>0,
\qquad
3f_{RR}<0,
\]

所以

\[
m_s^2<0.
\]

这意味着正曲率区对应 tachyonic scalaron 倾向，也就是小扰动可能指数增长。

### 不是瞬间数学灾难

这个结论不是说：

```text
一出现 R>0，理论立刻发散。
```

它真正说明的是：

```text
若把 C/D 当作标准 metric f(R) 且 scalaron 是真实传播自由度，
正 R 区域存在一个线性不稳定方向。
```

是否真的灾难，取决于：

- 正 \(\tilde R\) 区域是否持续存在；
- 该区域大小是否能支撑不稳定模式；
- 增长时间 \(1/\sqrt{|m_s^2|}\) 是否短于实验/演化时间；
- 初始扰动是否会投影到 scalaron 模式上；
- \(f_R\to0\) 区域是否已经超出该 effective theory 的适用范围；
- 是否改用 Palatini 或受约束辅助场，使 metric scalaron 不再传播。

### 当前判断

因此正 \(\tilde R\) 区 \(f_{RR}<0\) 是一个强红旗，但不是自动定理式排除。

更精确地说：

\[
\boxed{
它表示“可能发散/可能不稳定”，不是“必然立即发散”。
}
\]

若后续线性响应分析显示增长率很大、并且该模式在实验时间内可被激发，则 C/D 作为 metric `f(R)` 会被理论上强烈否定。

若正曲率区只是极薄、瞬时、且 scalaron 模式被 cutoff 或约束掉，则它可能作为 effective theory 保留。

## 2. Cauchy 适定是什么意思

`Cauchy 适定` 通常包含三件事：

1. 存在性：给定合法初始数据，至少局域时间内有解；
2. 唯一性：同一组初始数据不会演化出多个不同解；
3. 连续依赖性：初始数据或数值误差的小扰动，只导致解的小扰动。

对数值相对论，第三点尤其关键。若连续依赖性失败，即使离散残差很小、步长很小，数值误差也可能被 PDE 本身无限放大。

在双曲 PDE 语言里，这通常要求主部给出强双曲或至少良好的能量估计。

因此 `Cauchy 主部分析` 的目的就是先问：

```text
这个方程有没有资格作为一个演化方程？
```

而不是直接问：

```text
某个参考场代进去残差是不是小？
```

## 3. \(f_R\to0\) 时能否改看 \(f\to1/\ell^2\) 项

用户提出的问题非常关键：

```text
当 f_R -> 0 时，它包含的高阶导数项已经不是最高阶的非零项了，
是不是观察 f -> 1/ell^2 的项作为最高阶导数项即可？
```

答案是：不能直接这样做。

原因是 \(f\to1/\ell^2\) 的项是

\[
-\frac12 f\,\tilde g_{\mu\nu},
\]

它是零阶代数项，不含导数。

在完全饱和且 \(f_R\) 空间时间常数的区域中，metric `f(R)` 方程会退化为近似

\[
-\frac12 f_{\rm sat}\tilde g_{\mu\nu}
\approx
\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

这不是一个演化方程，而是代数约束。

所以这里发生的不是：

```text
四阶/二阶主部消失后，零阶项成为健康的新主部。
```

而是：

```text
方程主部降秩，度规动力学被关掉或退化。
```

这可能是我们想要的“屏蔽/冻结”机制，但它必须按退化 PDE 或分区边界问题来处理，不能按普通 Cauchy 演化方程直接推进。

### 分区图景

更合理的结构是三块：

1. EH-like 区：\(|\ell^2\tilde R|\ll1\)，\(f_R\approx1\)，度规有正常 Einstein-like 主部；
2. saturated bulk：\(|\ell^2\tilde R|\gg1\)，\(f_R\approx0\)，度规主部退化，方程近似代数约束；
3. transition layer：\(|\ell^2\tilde R|=O(1)\)，\(\nabla\nabla f_R\) 可能形成薄层主导项。

所以后续应把 C/D 看成：

```text
EH-like 动力学区 + 饱和退化区 + 过渡层边界条件
```

而不是一个处处普通的 metric `f(R)` 演化系统。

## 4. 过渡层积分边界条件

四维迹方程为

\[
M_P^2(f_R\tilde R-2f+3\tilde\Box f_R)=\tilde T.
\]

令 \(\phi=f_R\)。跨过一个法向坐标为 \(n\) 的薄层积分，得到近似：

\[
3M_P^2[\partial_n\phi]^+_-
\approx
\int_{\rm layer}\tilde T\,dn
-M_P^2\int_{\rm layer}(f_R\tilde R-2f)\,dn
+\text{切向导数与几何修正}.
\]

这说明过渡层不是纯数值麻烦。若 \(\partial_n f_R\) 的跳变有限，它会表现为有效薄壳源。

因此后续必须判断：

- 这个跳变是否随网格收敛到有限值；
- 它是否给 \(S,\rho\) 演化带来可观测相位偏移；
- 它是否只局限在密度极低区；
- 它是否需要作为边界条件而不是普通体区 PDE 来处理。

## 5. 线性响应稳定性估计

参考残差检验实际检验的是：

\[
\mathcal E[\Psi_{\rm ref}]=r.
\]

真正的解写成

\[
\Psi=\Psi_{\rm ref}+\delta\Psi.
\]

线性化后：

\[
\mathcal L\,\delta\Psi=-r.
\]

因此偏离大小不只取决于残差 \(r\)，还取决于

\[
\mathcal L^{-1}.
\]

若 Cauchy 主部退化、scalaron tachyonic、或者过渡层形成近零模，则

\[
\|\mathcal L^{-1}\|
\]

可能很大。此时小残差也可能导致大偏离。

所以后续要检验的不是一句“残差小”，而是：

\[
\boxed{
\|\delta\Psi\|\lesssim \|\mathcal L^{-1}\|\,\|r\|
}
\]

是否仍小。

## 6. `ell-window scan` 实现

新脚本：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/scan_cd_ell_window.py
```

输出：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/cd_ell_window_scan_64/summary.json
```

设置：

- 网格：`64x64`
- 时刻：`t=0,8,16`
- 分支：`C,D`
- `ell = 1,3,10,30,100,300`
- 主支撑区：`\rho > 1e-3 rho_max`

输出指标包括：

- residual `p95/abs_max`
- derivative norm
- `f_R < 1e-4` 的支撑区比例
- `f_RR < 0` 的支撑区比例
- scalaron \(m_s^2\) 的正负/undefined 比例
- 严格过渡层比例
- 热点与粗过渡带重合率

## 7. 首扫结果

### 7.1 residual 会推动 `ell` 变大

取三个时刻中的最大支撑区 residual `p95`：

| branch | ell=1 | ell=10 | ell=30 | ell=100 | ell=300 |
|---|---:|---:|---:|---:|---:|
| C | 1.49e4 | 3.92e2 | 9.61 | 4.98e-2 | 5.55e-3 |
| D | 1.38e4 | 3.14e2 | 1.24 | 4.98e-2 | 5.55e-3 |

所以若只看 bulk residual，确实会倾向选大 `ell`。

这支持用户的直觉：为了不显著偏离平直量子力学，`ell` 不能太小。

### 7.2 但大 `ell` 会把支撑区推入饱和退化

同一扫描显示，`f_R<1e-4` 的支撑区比例随 `ell` 很快接近 1。

尤其对 `ell>=100`，C/D 的主支撑区基本全在

\[
f_R\ll1
\]

的饱和/退化区。

这意味着：

```text
大 ell 确实屏蔽了残差，
但它也把正常 Einstein-like 主部几乎关掉。
```

所以大 `ell` 不是免费午餐。

### 7.3 C 支的正曲率红旗持续存在

C 支中 \(f_{RR}<0\) 的支撑区比例在本扫描里基本跟正 \(\tilde R\) 区比例同步，最大约 `0.53`。

这说明 C 支的 metric scalaron 红旗不是少数孤立点，而是支撑区中相当大的正曲率部分。

### 7.4 D 支的“大 ell 安全感”有假象

D 支在大 `ell` 下显示 \(f_{RR}<0\) 比例很小，甚至 `ell=300` 时接近 0。

但这不是简单的健康信号。原因是 D 的指数饱和太快：

\[
f_R=\operatorname{sech}^2(\ell^2\tilde R)
\]

和

\[
f_{RR}=-2\ell^2 f_R\tanh(\ell^2\tilde R)
\]

会在双精度数值中直接下溢到接近 0。

此时 scalaron 质量公式

\[
m_s^2=\frac{f_R-\tilde R f_{RR}}{3f_{RR}}
\]

大量 undefined。

所以 D 大 `ell` 下的正确解读是：

```text
不是 scalaron 红旗消失，
而是系统进入更强的饱和/退化极限，metric scalaron 公式本身失去普通解释。
```

### 7.5 严格过渡层在 `64x64` 下继续欠采样

对于 `ell>=100`，严格过渡层

\[
0.5\le|\ell^2\tilde R|\le2
\]

在主支撑区中几乎没有采样点。

因此本扫描可以说明：

```text
大 ell 的 bulk residual 小；
大 ell 的主部退化强；
但大 ell 的真正过渡层仍未解析。
```

它还不能说明过渡层薄壳效应安全。

## 8. 当前结论

本轮比上一轮更清楚地说明：

\[
\boxed{
残差要求推动 ell 变大；
Cauchy/主部健康性要求不能盲目变大。
}
\]

因此 `ell` 的物理标定不是只找一个让 residual 小的值，而是找一个同时满足：

1. bulk residual 足够小；
2. \(f_R\to0\) 的退化区不导致不可适定；
3. 正曲率 \(f_{RR}<0\) 区不在实验时间内触发 scalaron 增长；
4. 过渡层 jump 条件不给出可观测薄壳效应；
5. 普通 GR 实验仍处于 \(|\ell^2R|\ll1\) 区域。

## 9. 下一步

下一步应该不是立刻全演化，而是把上面三块补成更硬的判据：

1. 对 `ell=30,100,300` 的 D 支做局部窗口高分辨率过渡层积分，估计 \([\partial_n f_R]^+_-\)；
2. 对正 \(\tilde R\) 区估计 scalaron 增长时间 \(1/\sqrt{|m_s^2|}\)，并区分 finite 与 undefined 饱和区；
3. 写出分区 PDE 观点：EH-like 区、saturated algebraic 区、transition layer 的匹配条件；
4. 在这些判据通过后，再选一个中间 `ell` 进入全作用量演化。
