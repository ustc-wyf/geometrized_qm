# branch preservation 与 patch transition 首轮分析

日期：2026-05-08

## 一句话结论

当前 proposal v1 可以继续推进，但需要补两条“健康性规则”：

1. **branch preservation**：若要求
   \[
   Q\to0\Rightarrow C_{\mu\nu}\to0,
   \]
   则不能只在初始切片令 \(C=0\)。需要保证演化不会让
   \[
   C_{\mu\nu}=\chi(\mathcal q)\hat C_{\mu\nu}
   \]
   中的 \(\hat C_{\mu\nu}\) 在 \(\chi\to0\) 时发散。

2. **patch transition**：\(\Delta=0\) 是张量基底降秩；\(w^2=0\) 是 trace0 状态方程主符号退化。二者性质不同，不能混为一个“数值坏点”。

本轮结论是：branch preservation 可以通过“光滑二阶消失的 \(\chi\)”或“\(Q=0\) 面上的无通量边界条件”来控制；patch transition 则应把 \(C_{\mu\nu}\) 当作主变量、\(\lambda_I\) 只当作局部坐标，并在退化面切换到低维或替代闭合 patch。

## 1. branch preservation 的数学形式

设

\[
C_{\mu\nu}
=
\chi(\mathcal q)\hat C_{\mu\nu},
\qquad
\chi(0)=0.
\]

其中

\[
\mathcal q=\tilde Q_\rho/m^2
\]

是 massive branch 的无量纲量子强度。massless branch 需要另选尺度。

守恒条件为：

\[
\tilde\nabla^\mu C_{\mu\nu}=0.
\]

代入分支写法：

\[
\tilde\nabla^\mu
(\chi\hat C_{\mu\nu})
=
\chi\tilde\nabla^\mu\hat C_{\mu\nu}
+
(\tilde\nabla^\mu\chi)\hat C_{\mu\nu}
=0.
\]

如果直接除以 \(\chi\)，得到：

\[
\tilde\nabla^\mu\hat C_{\mu\nu}
+
(\tilde\nabla^\mu\ln\chi)\hat C_{\mu\nu}
=0.
\]

这里的问题很清楚：当 \(\chi\to0\) 时，\(\tilde\nabla\ln\chi\) 可能发散。因此 \(C\to0\) 是否保持，不只看 \(\chi(0)=0\)，还要看 \(\hat C\) 是否有界以及 \(\chi\) 的零点阶数。

## 2. \(\chi\) 的零点阶数

令局部

\[
\chi(\mathcal q)\sim |\mathcal q|^p.
\]

则

\[
\tilde\nabla_\mu\chi
\sim
p|\mathcal q|^{p-1}\tilde\nabla_\mu|\mathcal q|.
\]

若 \(\hat C_{\mu\nu}\) 有界，且 \(\tilde\nabla\mathcal q\) 有界，则：

- \(p>1\)：\(\tilde\nabla\chi\to0\)，守恒方程在 \(\mathcal q=0\) 面上不产生有限跳跃源；
- \(p=1\)：\(\tilde\nabla\chi\) 一般有限，守恒方程要求额外边界条件；
- \(0<p<1\)：\(\tilde\nabla\chi\) 发散，通常不适合作为正则分支。

因此当前最稳妥的选择是使用光滑的二阶消失函数，例如：

\[
\boxed{
\chi(\mathcal q)
=
\frac{\mathcal q^2}{\mathcal q^2+q_0^2}
}
\]

或任何满足

\[
\chi(0)=0,
\qquad
\chi'(0)=0,
\qquad
\chi>0\quad(\mathcal q\neq0)
\]

的光滑函数。

这比直接用 \(|\mathcal q|\) 更健康，因为 \(|\mathcal q|\) 在零点不可微。

## 3. 若使用一阶零点，需要什么边界条件

如果坚持使用

\[
\chi\sim|\mathcal q|,
\]

则在 \(\mathcal q=0\) 面上，守恒式

\[
\chi\tilde\nabla^\mu\hat C_{\mu\nu}
+
(\tilde\nabla^\mu\chi)\hat C_{\mu\nu}=0
\]

退化为：

\[
\boxed{
n^\mu\hat C_{\mu\nu}=0
\quad\text{on}\quad
\mathcal q=0,
}
\]

其中

\[
n_\mu\propto\tilde\nabla_\mu\mathcal q.
\]

这就是“无通量”边界条件：有效几何应力不能穿过 \(Q=0\) 分支边界向 GR 区域注入非零源。

所以一阶零点不是绝对不行，但要额外补边界匹配条件；二阶光滑零点则更自然。

## 4. 有界性条件

branch preservation 的核心不是 \(C\) 本身，而是

\[
\hat C_{\mu\nu}=C_{\mu\nu}/\chi.
\]

我们需要要求：

\[
\boxed{
\|\hat C\|_W<\infty
\quad
\text{near}\quad
\mathcal q=0.
}
\]

或者用辅助 action 的语言，要求有限加权能量：

\[
\int\sqrt{|\tilde g|}
\,
\mu(\mathcal q)\|C\|_W^2
<\infty,
\]

其中可取

\[
\mu(\mathcal q)\sim\chi(\mathcal q)^{-2}.
\]

这样有限能量就倾向于排除 \(C\) 不随 \(\chi\) 消失的分支。

但这仍只是分支选择/正则性条件；不是 trace0 和守恒自动推出的结论。

## 5. patch transition 的第一个原则：\(C\) 是主变量，\(\lambda_I\) 是坐标

在非退化 patch 内：

\[
C_{\mu\nu}=\lambda_I E^I_{\mu\nu},
\qquad
E=\mathrm{span}\{\tilde g,uu,rr,ur\}.
\]

当

\[
\Delta=u^2r^2-(u\cdot r)^2\to0
\]

时，Gram 矩阵

\[
\det H=\frac{d-2}{2}\Delta^3
\]

退化。此时发散的可能是系数 \(\lambda_I\)，不一定是张量 \(C_{\mu\nu}\) 本身。

因此 patch transition 的正确原则是：

\[
\boxed{
C_{\mu\nu}\text{ 是几何对象，}\lambda_I\text{ 只是局部坐标。}
}
\]

如果 \(\lambda_I\to\infty\) 但 \(C_{\mu\nu}\) 有有限极限，这更像坐标奇点；如果 \(C_{\mu\nu}\) 本身发散，才是物理/理论奇点。

## 6. \(\Delta=0\) 的低维 patch

当 \(r_\mu\) 与 \(u_\mu\) 线性相关，或 \(u,r\) 张成 null-degenerate 二平面时，\(E\) 降维。

最常见的光滑分支是：

\[
r_\mu=\alpha u_\mu.
\]

此时

\[
r_\mu r_\nu=\alpha^2 u_\mu u_\nu,
\qquad
u_{(\mu}r_{\nu)}=\alpha u_\mu u_\nu.
\]

于是

\[
E\to E_u
=
\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu\}.
\]

在这个 patch 中，应改写为：

\[
C_{\mu\nu}
=
A'\tilde g_{\mu\nu}
+
B'u_\mu u_\nu.
\]

trace0 给出：

\[
dA'+B'u^2=0.
\]

若 massive branch 上 \(u^2=m^2\)，则：

\[
A'=-\frac{m^2}{d}B'.
\]

所以低维 patch 中仍可有一个无迹方向：

\[
C_{\mu\nu}
=
B'
\left(
u_\mu u_\nu
-
\frac{m^2}{d}\tilde g_{\mu\nu}
\right).
\]

这说明 \(\Delta=0\) 不必自动判死；它可以是从四方向 patch 降到二维 patch。

但如果场方程要求的 \(\mathcal R_{\mu\nu}\) 在该点有无法落入 \(E_u\) 的分量，则当前最小 proposal 在该点失败，必须：

- 要么由 \(Q\to0\) 分支令 \(C\to0\)；
- 要么扩大张量基底，例如加入 \(\tilde\nabla_{(\mu}r_{\nu)}\)、\(\tilde\nabla_{(\mu}u_{\nu)}\) 等导数方向；
- 要么把该点作为真正的边界/缺陷面处理。

## 7. \(r^\perp=0\) 的意义

定义 \(u\)-正交的振幅梯度：

\[
r^\perp_\mu
=
r_\mu
-
\frac{u\cdot r}{u^2}u_\mu.
\]

当

\[
r^\perp_\mu=0
\]

时，\(r\) 与 \(u\) 平行，因此在 massive branch 中也有

\[
\Delta=0.
\]

所以 \(r^\perp=0\) 不是一个独立的新物理奇点，而是 \(\Delta=0\) 的一个常见可解释子情形。

在高斯波包中心，\(r_\mu\) 可能为零或接近平行于 \(u_\mu\)。这时如果 \(Q\) 同时很小，GR branch \(C\to0\) 最自然；如果 \(Q\) 不小，则需要低维 patch 或扩大基底。

## 8. \(w^2=0\) 不是同一种退化

此前主符号分析中定义：

\[
w_\mu
=
(\xi\cdot r)u_\mu
-
(\xi\cdot u)r_\mu.
\]

这里 \(\xi_\mu\) 是主符号/频率协向量，不是时空中的固定场。

trace0 对剩余未定主部模式的收缩正比于：

\[
\tilde g^{\mu\nu}(-w_\mu w_\nu)
=
-w^2.
\]

因此若

\[
w^2=0,
\]

普通 trace0 看不见这个主部模式。

这和 \(\Delta=0\) 不同：

- \(\Delta=0\)：张量基底 \(E\) 自身降秩，是 field-space/spacetime patch 问题；
- \(w^2=0\)：trace0 状态方程对某些主符号方向 \(\xi\) 失去控制，是 PDE 主部/特征方向问题。

所以 \(w^2=0\) 不能简单通过改写 \(\lambda_I\) 坐标解决。它需要检查状态方程本身。

## 9. \(w^2=0\) 的可能处理

候选处理有三类：

### 9.1 保持 trace0，但限制物理 patch

如果物理解的相关主符号方向从不进入 \(w^2=0\)，或者这些方向只出现在低密度/退化 patch 中，则 trace0 仍可作为主闭合。

这需要数值和解析共同检查，不能先验假定。

### 9.2 改用混合 trace 状态方程

把状态方程从普通 trace 改成：

\[
S^{\mu\nu}C_{\mu\nu}=0.
\]

其中 \(S^{\mu\nu}\) 需要满足：

\[
S^{\mu\nu}w_\mu w_\nu\neq0
\]

在目标 patch 内成立。

\(\mathsf q_E\)-trace 曾作为全局硬闭合表现较差，但它仍可能作为局部主符号补丁。更平滑的方式是定义：

\[
S^{\mu\nu}
=
(1-\eta)\tilde g^{\mu\nu}
+
\eta\,\mathsf q_E^{\mu\nu},
\]

其中 \(\mathsf q_E^{\mu\nu}\) 指此前构造的 \(u-r\) 二平面正定收缩张量，\(\eta\) 是只在 trace0 主符号退化附近打开的 bump function。

这会改变状态方程，因此不能悄悄当作“同一个理论”。若采用，proposal v1 就要升级为带 patchwise scalar closure 的 v2。

### 9.3 扩大方程组

如果既要保持普通 trace0，又要全局避免主符号退化，可能需要加入额外辅助场或高阶状态方程。这会改变最小理论规格，暂时不作为第一选择。

## 10. 当前阶段判断

### branch preservation

当前可接受的首版规则是：

\[
\boxed{
C_{\mu\nu}
=
\chi(\mathcal q)\hat C_{\mu\nu},
\qquad
\chi(0)=\chi'(0)=0,
\qquad
\hat C_{\mu\nu}\text{ 有界}.
}
\]

如果不用二阶消失的 \(\chi\)，则必须在 \(\mathcal q=0\) 面上补：

\[
\boxed{
n^\mu\hat C_{\mu\nu}=0.
}
\]

### patch transition

当前可接受的首版规则是：

\[
\boxed{
\text{以 }C_{\mu\nu}\text{ 为主变量；当 }E\text{ 降秩时，切到低维 }E_{\rm lim}\text{ patch。}
}
\]

对于

\[
r\parallel u,
\]

低维 patch 是：

\[
E_u=\mathrm{span}\{\tilde g,uu\}.
\]

若 \(\mathcal R\) 不能落入低维 patch，则当前最小 proposal 需要扩大基底或把该区域作为边界/缺陷处理。

对于 \(w^2=0\)，不能只靠基底换图，需要单独检查 trace0 状态方程的主符号。它是 proposal v1 的主要理论风险之一。

## 11. 下一步

下一步应做两个具体检查：

1. 在高斯干涉三切片上统计 \(\mathcal q\approx0\) 区域中 \(C/\chi\) 是否有界，测试二阶 \(\chi(\mathcal q)\) 的 branch preservation 口径。
2. 在 \(\Delta\approx0\) 与 \(r^\perp\approx0\) 区域，检查 \(\mathcal R\) 是否接近低维 patch \(E_u=\mathrm{span}\{\tilde g,uu\}\)；若不接近，说明最小基底必须扩展。
