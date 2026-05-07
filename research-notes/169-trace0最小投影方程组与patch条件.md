# 169-trace0最小投影方程组与patch条件

日期：2026-05-08

## 一句话结论

当前最小、且被三切片高斯干涉数据相对支持的 equation-first 候选不是
\(\mathsf q_E\)-trace 闭合，而是：

\[
\boxed{
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu},
\qquad
\mathcal C_{\mu\nu}\in E,
\qquad
\tilde g^{\mu\nu}\mathcal C_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
}
\]

其中

\[
E=\mathrm{span}
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

这可以等价写成投影形式：

\[
\boxed{
\Pi_E^\perp
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)=0,
}
\]

\[
\boxed{
\tilde g^{\mu\nu}
\Pi_E
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)=0,
}
\]

\[
\boxed{
\tilde\nabla^\mu
\Pi_E
\left(
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}
\right)=0.
}
\]

这组方程不是最终理论，只是目前最小显式候选。它还必须附带非退化 patch、\(Q\to0\) 回 Einstein 的分支条件，以及后续 Helmholtz/变分性检查。

## 1. 变量和口径

本笔记全程在 \(\tilde g\) 表象写方程。

基本变量取为：

\[
\tilde g_{\mu\nu},\qquad \tilde\rho,\qquad S.
\]

定义：

\[
u_\mu=\partial_\mu S.
\]

振幅方向写作：

\[
r_\mu=\tilde\nabla_\mu\ln\sqrt{\tilde\rho}.
\]

这里特别注意：

\[
\boxed{\rho\neq\tilde\rho.}
\]

若和 A 支平直量子力学快照比较，正确比较对象不是直接把 \(\tilde\rho\) 当成 \(\rho_A\)，而是先通过变换关系把 \(\tilde\rho\) 拉回 \(g\) 表象，再和 \(\rho_A\) 比较。

物质能动量记为

\[
\tilde T_{\mu\nu}.
\]

具体是否有 \(2\)、符号或 \(m^2\) 归一化因子，取决于物质作用量规范。本笔记只固定结构：

\[
\mathcal R_{\mu\nu}
\equiv
\tilde G_{\mu\nu}
-\frac{1}{M_P^2}\tilde T_{\mu\nu}.
\]

所有数值 residual 也按这个 Einstein 残差的口径理解。

## 2. 最小显式方程组

定义四个张量方向：

\[
E^0_{\mu\nu}=\tilde g_{\mu\nu},
\]

\[
E^1_{\mu\nu}=u_\mu u_\nu,
\]

\[
E^2_{\mu\nu}=r_\mu r_\nu,
\]

\[
E^3_{\mu\nu}=u_{(\mu}r_{\nu)}
=\frac12(u_\mu r_\nu+r_\mu u_\nu).
\]

修正张量写作：

\[
\mathcal C_{\mu\nu}
=
A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C r_\mu r_\nu
+D u_{(\mu}r_{\nu)}.
\]

候选场方程为：

\[
\boxed{
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T_{\mu\nu}
+\mathcal C_{\mu\nu}.
}
\]

第一条闭合条件是普通 trace=0：

\[
\boxed{
\tilde g^{\mu\nu}\mathcal C_{\mu\nu}=0.
}
\]

第二条是兼容/传播条件，即修正张量守恒：

\[
\boxed{
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
}
\]

它不是任意添加的数值规则。由于 Bianchi 恒等式给出

\[
\tilde\nabla^\mu\tilde G_{\mu\nu}=0,
\]

若物质部分仍满足

\[
\tilde\nabla^\mu\tilde T_{\mu\nu}=0,
\]

则完整方程自动要求 \(\mathcal C_{\mu\nu}\) 守恒。反过来，若先把 \(\mathcal C_{\mu\nu}\) 当作额外几何/介质张量来求解，守恒条件就是它和物质测地线解释兼容的传播条件。

## 3. 投影写法

在非退化 patch 内，定义投影：

\[
\Pi_E\mathcal R_{\mu\nu}
=
E^I_{\mu\nu}(H^{-1})_{IJ}
\langle E^J,\mathcal R\rangle,
\]

其中

\[
H_{IJ}=\langle E_I,E_J\rangle.
\]

于是

\[
\mathcal C_{\mu\nu}=\Pi_E\mathcal R_{\mu\nu}.
\]

场方程可以不显式引入 \(A,B,C,D\)，直接写成：

\[
\Pi_E^\perp\mathcal R_{\mu\nu}=0.
\]

普通 trace 闭合为：

\[
\tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0.
\]

在投影壳上 \(\mathcal R_{\mu\nu}=\Pi_E\mathcal R_{\mu\nu}\)，因此 trace 条件等价于

\[
\tilde g^{\mu\nu}\mathcal R_{\mu\nu}=0.
\]

在 \(3+1d\) 中：

\[
\tilde g^{\mu\nu}\tilde G_{\mu\nu}=-\tilde R,
\]

所以得到：

\[
\boxed{
-\tilde R-\frac{1}{M_P^2}\tilde T=0.
}
\]

这和普通 Einstein 方程取迹的形式一致；区别是本理论不要求整个 \(\mathcal R_{\mu\nu}=0\)，只要求它的非 \(E\) 分量为零，且 \(E\) 内剩余部分无迹并守恒。

## 4. 物质方程与测地线

物质部分仍应保留 Hamilton-Jacobi 壳方程和连续性方程：

\[
\boxed{
\tilde g^{\mu\nu}u_\mu u_\nu=m^2,
}
\]

\[
\boxed{
\tilde\nabla_\mu
\left(
\tilde\rho\,\tilde g^{\mu\nu}u_\nu
\right)=0.
}
\]

由于 \(u_\mu=\partial_\mu S\)，壳方程给出：

\[
u^\mu\tilde\nabla_\mu u_\nu
=
\frac12\tilde\nabla_\nu(u^\mu u_\mu)
=0.
\]

因此隐变量粒子仍沿 \(\tilde g\) 测地线运动。

这里要区分两件事：

1. equation-first 写法可以把 \(\mathcal C_{\mu\nu}\) 作为场方程侧的几何/介质修正；
2. 若以后要求它来自作用量，则引力作用量若显含 \(\tilde\rho,S\)，可能会反过来改变物质变分。

所以本候选还必须做 Helmholtz/变分性检查，不能直接宣称已经有 action。

## 5. \(Q\to0\) 回 Einstein 的分支条件

trace=0 不足以保证回到 Einstein。因为它只要求

\[
\mathcal C^\mu{}_\mu=0,
\]

并不强制

\[
\mathcal C_{\mu\nu}=0.
\]

因此还需要分支选择：

\[
\boxed{
Q\to0
\quad\Longrightarrow\quad
\mathcal C_{\mu\nu}\to0.
}
\]

这更像边界条件、正则性条件或辅助场 action 的真空选择，而不是再额外加一个逐点硬方程。否则容易把系统过约束。

一个可检验的局域版本是引入门控尺度：

\[
\mathcal C_{\mu\nu}
=
\chi(\mathcal Q)\,\widehat{\mathcal C}_{\mu\nu},
\qquad
\chi(0)=0,
\]

但前面对简单 \(F(Q)\) 的测试已经说明，不能把它粗暴简化为 trace=\(aQ+bQ^2\)。门控更适合作为分支选择或辅助泛函权重。

## 6. 非退化 patch 条件

投影方程只有在张量基底非退化时才是单图册表达。

Gram 行列式为：

\[
\det H=
\frac{d-2}{2}
\left(
u^2r^2-(u\cdot r)^2
\right)^3.
\]

因此至少需要：

\[
\boxed{
d>2,
\qquad
\Delta\equiv u^2r^2-(u\cdot r)^2\neq0.
}
\]

若 \(d=2\)，或者 \(r_\mu\) 与 \(u_\mu\) 线性相关，则

\[
\tilde g_{\mu\nu},\quad
u_\mu u_\nu,\quad
r_\mu r_\nu,\quad
u_{(\mu}r_{\nu)}
\]

不再是四个独立方向，必须换低维 basis 或切换 patch。

普通 trace=0 还有一个主符号风险。上一轮发现未定模式为

\[
\bar h^{(0)}_{\mu\nu}=-w_\mu w_\nu,
\]

\[
w_\mu=(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu.
\]

普通 trace 对它的收缩是：

\[
\tilde g^{\mu\nu}\bar h^{(0)}_{\mu\nu}=-w^2.
\]

所以当

\[
\boxed{
w^2=0
}
\]

时，普通 trace 会漏掉这个主部模式。

当前三切片检查显示，\(w^2\simeq0\) 点主要集中在 \(r^\perp\) 也退化或近退化的区域。因此这更像图册边界问题，而不是 trace=0 的全局否定。

## 7. patch 规则的当前版本

当前可采用的临时规则是：

1. 主 patch：
   \[
   d>2,\qquad \Delta\neq0,\qquad w^2\not\simeq0.
   \]
   使用 trace=0 最小闭合。
2. \(r^\perp\) 退化 patch：
   当 \(r_\mu\parallel u_\mu\) 或 \(r^\perp\to0\) 时，丢弃冗余 \(r\) 方向，改用较低维 tensor basis。
3. \(w^2\simeq0\) 主符号 patch：
   trace=0 可能漏掉一个模式，需要切换到更强的标量闭合或辅助场正定项，例如 \(\mathsf q_E\)-norm。
4. 严格 \(1+1d\)：
   不机械使用 3+1d 四方向投影；要单独写 \(1+1d\) 的独立 basis。

这些 patch 规则不是人为削峰，也不是数值阻尼。它们对应的是同一个协变方程在不同非退化图册中的表达方式。

## 8. 为什么分离态 residual 更大

上一轮三切片结果中，分离态 residual 比干涉中心更大。这并不表示“干涉时理论更好、分离时物理更坏”，而是说明当前最小闭合的误差来源主要不是干涉强度本身。

更具体地说：

1. 分离态 trusted 区域包含更多支撑边缘和尾部点；
2. 这些区域更容易出现 \(r^\perp\) 近退化、\(w^2\simeq0\) 或 basis condition number 变坏；
3. 干涉中心虽然振荡更强，但 \(u,r\) 张量方向更丰富，反而更容易让 \(\mathcal R_{\mu\nu}\) 落入 \(E\) 子空间；
4. 已有守恒 residual 很小，因此当前问题更像代数张量方向/patch 问题，而不是守恒方程本身失败。

所以 residual 变大提示的是：

\[
\boxed{
\text{trace=0 最小闭合在分离态的 chart/patch 条件更紧张。}
}
\]

## 9. 下一步

下一步不应回到“大量数值演化器调参”，而应继续完成方程本身的理论筛选：

1. Helmholtz/self-adjoint 检查：
   判断
   \[
   \Pi_E^\perp\mathcal R=0,\qquad
   \tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0
   \]
   是否可能来自某个局域作用量，或至少来自带辅助场的作用量。
2. patch 方程写法：
   把 \(\Delta=0\)、\(r^\perp=0\)、\(w^2=0\) 分别写成明确的 chart transition，而不是数值 `guard`。
3. 三切片复检：
   在 `n=384` 可信条纹分辨率上复查 trace=0 最小闭合，确认 `n=96` 结论不是低分辨率误导。
4. 分支条件：
   检查 \(Q\to0\Rightarrow\mathcal C\to0\) 能否由边界条件、辅助泛函或正则性唯一选出。

当前候选的地位是：

\[
\boxed{
\text{它是目前最清楚的显式 equation-first 方程组，但还不是最终理论。}
}
\]
