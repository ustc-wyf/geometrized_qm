# 174-stealth-state-action候选的最小no-go

日期：2026-05-08

## 一句话结论

尝试构造具体 stealth-state action 后，出现一个非常明确的最小 no-go：

\[
\boxed{
\text{在最小局域 ansatz 下，}
C_{\mu\nu}
\text{ 不能同时作为可变分辅助场、无迹非零源、stealth-state 约束对象。}
}
\]

原因很简单但很致命：

1. 要让 \(C_{\mu\nu}\) 成为 metric 方程右边的源，需要线性项
   \[
   -\frac{M_P^2}{2}\int\sqrt{|\tilde g|}\,
   \tilde g^{\mu\nu}C_{\mu\nu}.
   \]
2. 但如果 \(C_{\mu\nu}\) 是自由变分场，这个线性项的 \(C\)-变分给出
   \[
   -\frac{M_P^2}{2}\tilde g^{\mu\nu},
   \]
   必须被其他项抵消。
3. 最自然的无迹乘子项正好抵消它；一旦抵消，产生 metric 源的线性项本身也在壳上消失。
4. 平方型 stealth-state 项在目标约束面上一阶变分为零，不能抵消这个线性源。

所以：

\[
\boxed{
\text{最小 stealth-state action 不能生成非零的 trace=0 投影源。}
}
\]

要继续 action 化，必须引入更非最小的补偿场/介质变量，或者接受当前理论先作为 equation-first 方程。

## 1. 要求回顾

我们希望 action 化的目标是：

\[
M_P^2\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}^{(m)}
+\Theta_{\mu\nu},
\]

并且在目标分支上：

\[
\Theta_{\mu\nu}/M_P^2=C_{\mu\nu},
\]

\[
C_{\mu\nu}\in E_{\rm sh},
\qquad
C^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0.
\]

其中：

\[
E_{\rm sh}
=
\mathrm{span}
\{
\tilde g_{\mu\nu},
U_\mu U_\nu,
R_\mu R_\nu,
U_{(\mu}R_{\nu)}
\}.
\]

还要求状态项在物理分支上 stealth：

\[
\left.\frac{\delta S_{\rm state}}{\delta\Phi}\right|_{\rm phys}=0,
\qquad
\left.\frac{\delta S_{\rm state}}{\delta\sigma}\right|_{\rm phys}=0.
\]

这里：

\[
U_\mu=\partial_\mu\Phi,
\qquad
R_\mu=\tilde\nabla_\mu\ln\sqrt{\sigma}.
\]

## 2. 最小候选 action

最小候选自然写为：

\[
S_{\rm aux}^{\rm trial}
=
S_{\rm source}
+S_{\rm tr}
+S_{\rm stealth}^{(2)}
+S_{\rm sh}.
\]

其中：

\[
S_{\rm source}
=
-\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\,
T_C,
\]

\[
T_C\equiv\tilde g^{\mu\nu}C_{\mu\nu}.
\]

无迹约束用乘子：

\[
S_{\rm tr}
=
\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\,
\alpha T_C.
\]

stealth-state 平方项取：

\[
S_{\rm stealth}^{(2)}
=
\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\left[
\kappa_\perp\|\Pi_{E_{\rm sh}}^\perp C\|^2
+\kappa_{\rm div}\|\tilde\nabla^\mu C_{\mu\nu}\|^2
\right].
\]

这个平方项的好处是：在

\[
\Pi_{E_{\rm sh}}^\perp C=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0
\]

上，它对 \(\Phi,\sigma\) 的一阶变分为零。因此它是 stealth 的。

## 3. \(C\)-变分的冲突

把线性部分合并：

\[
S_{\rm lin}
=
-\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
(1-\alpha)T_C.
\]

对 \(C_{\mu\nu}\) 变分得到：

\[
\frac{\delta S_{\rm lin}}{\delta C_{\mu\nu}}
=
-\frac{M_P^2}{2}
\sqrt{|\tilde g|}
(1-\alpha)\tilde g^{\mu\nu}.
\]

在目标约束面上：

\[
\Pi_{E_{\rm sh}}^\perp C=0,
\qquad
\tilde\nabla^\mu C_{\mu\nu}=0,
\]

平方型 stealth 项的一阶 \(C\)-变分也为零。

因此 \(C\)-方程要求：

\[
(1-\alpha)\tilde g^{\mu\nu}=0.
\]

也就是：

\[
\boxed{
\alpha=1.
}
\]

但一旦 \(\alpha=1\)，线性项整体变为：

\[
S_{\rm lin}=0.
\]

于是它对 metric 的变分也为零：

\[
\Theta_{\mu\nu}^{\rm lin}=0.
\]

这正是冲突：

\[
\boxed{
\text{为了让 }C\text{-变分成立，线性源必须被抵消；}
\text{但一旦抵消，metric 源也消失。}
}
\]

所以这个最小候选不能产生

\[
\Theta_{\mu\nu}=M_P^2 C_{\mu\nu}\neq0.
\]

## 4. 为什么平方型 stealth 项救不了

平方型状态项的设计目标是：

\[
S_{\rm stealth}^{(2)}=0,
\qquad
\delta S_{\rm stealth}^{(2)}=0
\quad
\text{在目标约束面上。}
\]

这正是它能保护 shadow 分支的原因。

但同一个性质也意味着：

\[
\left.
\frac{\delta S_{\rm stealth}^{(2)}}{\delta C_{\mu\nu}}
\right|_{\rm target}=0.
\]

所以它不能抵消 \(S_{\rm source}\) 的线性 \(C\)-变分。

如果把平方项改成非零一阶变分，它就可能决定 \(C\)，但也会通过

\[
\delta_{\Phi,\sigma}\Pi_{E_{\rm sh}}
\]

推动 shadow fields 偏离物理分支。

因此这里有一个清楚的二难：

\[
\boxed{
\text{stealth 的状态项不能决定 }C;
\qquad
\text{能决定 }C\text{ 的状态项通常不 stealth。}
}
\]

## 5. 一般局域势 \(V(C)\) 能否救场

可以尝试加入：

\[
S_V
=
\int\sqrt{|\tilde g|}\,V(C).
\]

如果 \(V\) 不依赖 \(U,R,\Phi,\sigma\)，它不会破坏 shadow 分支。

但为了抵消 \(S_{\rm source}\) 的 \(C\)-变分，需要：

\[
\frac{\partial V}{\partial C_{\mu\nu}}
\sim
\frac{M_P^2}{2}\tilde g^{\mu\nu}
\]

在目标 \(C\) 上成立。

若 \(C\) 被要求无迹，局域标量势通常分解为 trace 部分和 traceless 部分：

\[
C_{\mu\nu}
=
\widehat C_{\mu\nu}
+\frac1d \tilde g_{\mu\nu}T_C.
\]

在

\[
T_C=0
\]

处，能产生纯 \(\tilde g^{\mu\nu}\) 导数的正是 \(V\) 对 \(T_C\) 的线性导数。

但这等价于又加入一个线性 trace 项。它会重复第 3 节的问题：为了抵消 \(C\)-变分，线性 trace 系数会抵消 source，也会削掉 metric 源。

所以普通局域 \(V(C)\) 不能自然救场。

它可以做到两件事：

1. 决定 traceless 部分 \(\widehat C_{\mu\nu}\) 的大小；
2. 或者决定 trace 部分。

但它很难同时做到：

\[
\boxed{
C^\mu{}_\mu=0,
\qquad
\delta_C S=0,
\qquad
\Theta_{\mu\nu}=M_P^2 C_{\mu\nu}\neq0.
}
\]

## 6. 可能的逃逸路线

这个 no-go 不是完全否定 action 化，而是否定最小 ansatz。

可能逃逸路线有四类。

### 路线 A：\(C_{\mu\nu}\) 不作为自由变分场

把 \(C_{\mu\nu}\) 看成由其他自由场的有效应力自动生成：

\[
C_{\mu\nu}
=
\Theta_{\mu\nu}[\chi,\tilde g,U,R]/M_P^2.
\]

此时不再对 \(C\) 变分，也就没有上面的 \(C\)-方程冲突。

代价是必须找到具体 \(\chi\) sector，使 \(\Theta_{\mu\nu}\) 自动落在 \(E_{\rm sh}\)、无迹、守恒。

### 路线 B：接受 equation-first

保留方程：

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2
+C_{\mu\nu},
\]

\[
C\in E,\qquad C^\mu{}_\mu=0,\qquad \tilde\nabla^\mu C_{\mu\nu}=0.
\]

不再要求它来自普通局域 action。

这条路规格低一些，但目前最清楚，也最少引入新自由度。

### 路线 C：非局域或边界型分支选择

如果 \(Q\to0\Rightarrow C\to0\) 和 \(C\in E\) 不是由局域 action 强制，而是由边界条件、正则性或非局域投影选出，则可绕开局域 \(C\)-变分冲突。

代价是理论会变成非局域或带全局分支选择。

### 路线 D：带补偿的非最小乘子系统

引入额外补偿场，使乘子贡献的 shadow-source 和 metric stress 在物理分支上相互抵消。

这在形式上可能，但自由度会明显增加，且需要新的稳定性分析。

## 7. 当前结论

本轮真正得到的是一个负面但有用的结论：

\[
\boxed{
\text{最小 stealth-state action 不足以 action 化 trace=0 投影方程。}
}
\]

更具体地说：

1. 线性 source 项可以给 metric 方程提供 \(C_{\mu\nu}\)；
2. 但如果 \(C_{\mu\nu}\) 是自由场，它的变分必须被抵消；
3. 抵消项若是 trace 乘子或 trace 势，会同时取消 metric source；
4. stealth 平方项不会破坏分支，但也无法决定 \(C\)；
5. 非 stealth 项可以决定 \(C\)，但会推动 shadow fields 偏离物理分支。

因此当前 action 化路线已经收窄到：

\[
\boxed{
\text{寻找真实辅助场 }\chi\text{ 的有效应力 }\Theta_{\mu\nu},
\quad
\text{或承认 equation-first 是当前主理论形式。}
}
\]

## 8. 下一步

下一步有两个合理选择。

理论选择：

1. 尝试构造真实 \(\chi\)-sector，使其应力张量天然具有
   \[
   \Theta_{\mu\nu}/M_P^2\in E_{\rm sh},
   \qquad
   \Theta^\mu{}_\mu=0,
   \qquad
   \tilde\nabla^\mu\Theta_{\mu\nu}=0.
   \]
2. 如果做不到，就正式把 action 化降级，把 trace=0 投影方程作为 equation-first 理论继续推进。

数值选择：

1. 立刻做 `n=384` 复检，确认 trace=0 投影方程在可信条纹分辨率下仍然优于 \(\mathsf q_E\)-trace；
2. 若复检失败，则无需继续 action 化当前 trace=0 版本；
3. 若复检通过，再决定是否值得引入 \(\chi\)-sector。
