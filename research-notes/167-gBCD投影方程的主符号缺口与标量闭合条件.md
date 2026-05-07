# 167-gBCD投影方程的主符号缺口与标量闭合条件

日期：2026-05-07

## 一句话结论

投影型 gBCD 方程

\[
\Pi_E^\perp\mathcal R_{\mu\nu}=0
\]

加上 harmonic gauge 后，主符号仍然留下一条 \(E\)-方向的未定模式。因此它还不是完整的 Cauchy 动力学方程。

这不是数值问题，而是一个清楚的理论闭合要求：

\[
\boxed{
\text{还必须增加一个协变标量状态方程，来固定 }E\text{ 子空间中的剩余主部模式。}
}
\]

这也解释了为什么过去的 \(A,B,C,D\) 辅助场守恒方程总会出现一个 null direction：它不是算法偶然，而是方程结构本身需要一个额外标量闭合。

## 1. reduced principal symbol

记

\[
\mathcal R_{\mu\nu}
=
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

投影方程为

\[
\Pi_E^\perp\mathcal R_{\mu\nu}=0.
\]

取 trace-reversed 扰动

\[
\bar h_{\mu\nu}
=
h_{\mu\nu}-\frac12\tilde g_{\mu\nu}h.
\]

harmonic gauge 为

\[
F_\nu[h]
=
\tilde\nabla^\mu\bar h_{\mu\nu}=0.
\]

对 Fourier 主符号协向量 \(\xi_\mu\)，忽略常数因子后，reduced principal system 是：

\[
\boxed{
\xi^2\Pi_E^\perp\bar h_{\mu\nu}=0,
}
\]

\[
\boxed{
\xi^\mu\bar h_{\mu\nu}=0.
}
\]

这里

\[
\xi^2=\tilde g^{\mu\nu}\xi_\mu\xi_\nu.
\]

若 \(\xi^2\neq0\)，第一式给出

\[
\Pi_E^\perp\bar h=0,
\]

也就是

\[
\bar h_{\mu\nu}=\lambda_I E^I_{\mu\nu}.
\]

因此主符号可逆性的问题变成：

\[
\xi^\mu \lambda_I E^I_{\mu\nu}=0
\]

是否只允许 \(\lambda_I=0\)。

## 2. harmonic gauge 对 \(E\) 子空间的秩

定义矩阵

\[
M_{\nu I}(\xi)
=
\xi^\mu E^I_{\mu\nu}.
\]

四个列向量为

\[
M_{\nu 0}=\xi_\nu,
\]

\[
M_{\nu 1}=(\xi\cdot u)u_\nu,
\]

\[
M_{\nu 2}=(\xi\cdot r)r_\nu,
\]

\[
M_{\nu 3}
=
\frac12\left[
(\xi\cdot u)r_\nu
+(\xi\cdot r)u_\nu
\right].
\]

这些列向量都落在

\[
\mathrm{span}\{\xi_\nu,u_\nu,r_\nu\}
\]

中，所以

\[
\mathrm{rank}\,M\le3.
\]

在一般 \(3+1d\) 非退化情形下，若 \(\xi,u,r\) 线性无关且 \(\xi\cdot u,\xi\cdot r\) 不同时为零，则

\[
\mathrm{rank}\,M=3.
\]

于是 \(4\) 个 \(E\)-方向系数中仍有

\[
4-3=1
\]

个主部未定模式。

## 3. 剩余未定模式的显式形式

令

\[
p=\xi\cdot u,\qquad q=\xi\cdot r.
\]

一个核向量是

\[
N^I(\xi)
=
\left(
0,\,
-q^2,\,
-p^2,\,
2pq
\right).
\]

对应的张量扰动为

\[
\bar h^{(0)}_{\mu\nu}
=
-q^2 u_\mu u_\nu
-p^2 r_\mu r_\nu
+2pq\,u_{(\mu}r_{\nu)}.
\]

等价地，定义

\[
w_\mu=q\,u_\mu-p\,r_\mu,
\]

则

\[
\boxed{
\bar h^{(0)}_{\mu\nu}=-w_\mu w_\nu.
}
\]

并且

\[
\xi^\mu w_\mu
=
q(\xi\cdot u)-p(\xi\cdot r)
=0.
\]

所以这个模式正是 \(u-r\) 二平面中与主符号方向 \(\xi\) 横向的 rank-one 模式。

这说明：

\[
\boxed{
\Pi_E^\perp\mathcal R=0+\text{harmonic gauge}
\text{ 还不能固定 }w_\mu w_\nu\text{ 方向。}
}
\]

## 4. 为什么这和以前的 nullspace 结果一致

辅助应力守恒的主符号是

\[
\sigma_\xi(\tilde\nabla^\mu\mathcal C_{\mu\nu})
=
\xi^\mu \lambda_I E^I_{\mu\nu}
=
M_{\nu I}(\xi)\lambda_I.
\]

刚才已经看到 \(M\) 泛型 rank 为 \(3\)，所以它有一个 null direction。

因此过去数值中看到的“full conservation 主符号 rank=3、nullity=1”，不是 KKT 算法的偶然问题，而是这个方程族的真实结构。

结论是：

\[
\boxed{
\text{守恒方程能给出 3 个独立主部条件，还差 1 个标量状态方程。}
}
\]

## 5. 标量闭合条件的一般判据

假设补一个标量闭合

\[
\boxed{
\mathcal S(\lambda_I;\tilde g,u,r,\rho,S)=0.
}
\]

在线性主部上，它必须满足：

\[
\boxed{
\sigma_\xi(\mathcal S)[N(\xi)]\neq0
}
\]

对所有允许的非退化 Cauchy 主方向 \(\xi\) 成立。

如果闭合在最高阶上等价于

\[
s_I\lambda_I=\text{lower-order/source},
\]

则条件是

\[
\boxed{
s_I N^I(\xi)
\neq0.
}
\]

也就是

\[
-s_B(\xi\cdot r)^2
-s_C(\xi\cdot u)^2
+2s_D(\xi\cdot u)(\xi\cdot r)
\neq0.
\]

这给出一个非常直接的筛选标准：

- 任何标量状态方程如果对 \(N^I(\xi)\) 正交，就不能闭合主部；
- 单纯 trace 条件未必可靠，因为它可能在某些方向上与 \(N^I(\xi)\) 正交；
- 一个好的状态方程必须直接看见 \(u-r\) 平面中的这个 rank-one 横向模式。

## 6. 一个自然的协变候选：\(u-r\) 平面正定迹

在 massive branch 中

\[
u^2=m^2>0.
\]

令

\[
e^0_\mu=\frac{u_\mu}{\sqrt{u^2}}.
\]

取 \(r\) 垂直于 \(u\) 的部分：

\[
r^\perp_\mu
=
r_\mu-\frac{u\cdot r}{u^2}u_\mu.
\]

在 \(\Delta\neq0\) 的 massive branch 中通常有

\[
(r^\perp)^2<0.
\]

定义

\[
e^1_\mu
=
\frac{r^\perp_\mu}{\sqrt{-(r^\perp)^2}}.
\]

再定义 \(u-r\) 二平面上的正定收缩张量

\[
\mathsf q_E^{\mu\nu}
=
e_0^\mu e_0^\nu
+e_1^\mu e_1^\nu.
\]

这里 \(\mathsf q_E\) 不是 spacetime metric，而是只在 \(\{u,r\}\) 二平面上取正定范数的辅助张量。

一个自然标量闭合类是

\[
\boxed{
\mathsf q_E^{\mu\nu}
(\Pi_E\mathcal R)_{\mu\nu}
=
\chi(\mathcal Q)\,\Theta.
}
\]

其中

\[
\chi(0)=0,
\]

\(\Theta\) 是由低阶标量或辅助场给出的 regular 标量。

如果先取最小版本，可以测试

\[
\boxed{
\mathsf q_E^{\mu\nu}
(\Pi_E\mathcal R)_{\mu\nu}
=0.
}
\]

但这只是最小闭合候选，不应直接当作最终物理定律。

为什么这个闭合能看见缺口？

对未定模式

\[
N_{\mu\nu}=-w_\mu w_\nu,
\]

有

\[
\mathsf q_E^{\mu\nu}N_{\mu\nu}
=
-
\left[
(e_0^\mu w_\mu)^2
+(e_1^\mu w_\mu)^2
\right].
\]

只要 \(w_\mu\neq0\)，就有

\[
\mathsf q_E^{\mu\nu}N_{\mu\nu}<0.
\]

因此它不会漏掉这个主部 null mode。

## 7. 当前理论位置

目前可以明确说：

1. \(\Pi_E^\perp\mathcal R=0\) 是一个清楚的 6 方程投影型 Einstein-like 方程；
2. 它的退化 patch 条件已经明确为 \(d=2\) 或 \(\Delta=0\)；
3. 加 harmonic gauge 后，仍有一个 \(E\)-方向主部模式未定；
4. 因此完整理论必须再给出一个协变标量状态方程；
5. 这个标量方程最好同时承担 \(Q\to0\) 回 Einstein 的分支选择；
6. \(\mathsf q_E^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=\chi(\mathcal Q)\Theta\) 是当前最自然的首个候选闭合类。

这一步的意义是：我们没有把困难藏在数值器里，而是把投影理论缺的最后一个主部条件显式找出来了。

## 8. 下一步

下一步应做两件事：

1. 在高斯干涉三切片上检验 \(\mathsf q_E\)-trace 闭合是否比普通 trace 条件更符合已有数据；
2. 从 action/Helmholtz 角度检查这个正定平面迹闭合能否来自一个辅助场泛函，而不是人为添加。

