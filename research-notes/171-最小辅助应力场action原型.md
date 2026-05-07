# 171-最小辅助应力场action原型

日期：2026-05-08

## 一句话结论

如果要把 trace=0 投影方程 action 化，不能简单把

\[
\mathcal C_{\mu\nu}
=
A\tilde g_{\mu\nu}
+B u_\mu u_\nu
+C r_\mu r_\nu
+D u_{(\mu}r_{\nu)}
\]

直接塞进作用量，然后期待 metric 变分给出

\[
\tilde G_{\mu\nu}
=
\tilde T_{\mu\nu}/M_P^2+\mathcal C_{\mu\nu}.
\]

原因是：一旦 \(\mathcal C_{\mu\nu}\) 依赖 \(\tilde g,u,r,\tilde\rho,S\)，metric 变分和物质变分都会产生额外项。

因此 action 化时真正应该定义的是辅助 sector 的有效能动量：

\[
\boxed{
\Theta_{\mu\nu}
\equiv
-\frac{2}{\sqrt{|\tilde g|}}
\frac{\delta S_{\rm aux}}{\delta \tilde g^{\mu\nu}},
}
\]

然后要求

\[
\boxed{
\Theta_{\mu\nu}/M_P^2
\in
\mathrm{span}\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\},
\qquad
\Theta^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu\Theta_{\mu\nu}=0.
}
\]

也就是说：

\[
\boxed{
\lambda_I E^I_{\mu\nu}
\text{ 不能先验当作 action 里的裸源；}
\quad
\lambda_I
\text{ 更应当看成生成有效应力 }\Theta_{\mu\nu}
\text{ 的辅助变量。}
}
\]

## 1. 目标方程

当前 equation-first 候选为：

\[
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}\tilde T^{(m)}_{\mu\nu}
+\mathcal C_{\mu\nu},
\]

其中

\[
\mathcal C_{\mu\nu}
=
\lambda_I E^I_{\mu\nu},
\]

\[
E^I_{\mu\nu}
=
\{
\tilde g_{\mu\nu},
u_\mu u_\nu,
r_\mu r_\nu,
u_{(\mu}r_{\nu)}
\}.
\]

并附带：

\[
\mathcal C^\mu{}_\mu=0,
\]

\[
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0,
\]

\[
Q\to0\quad\Longrightarrow\quad \mathcal C_{\mu\nu}\to0.
\]

若来自普通作用量，最自然形式应为：

\[
S
=
S_{\rm EH}[\tilde g]
+S_m[\tilde g,\tilde\rho,S]
+S_{\rm aux}[\tilde g,\text{aux};\tilde\rho,S].
\]

metric 变分给出：

\[
M_P^2\tilde G_{\mu\nu}
=
\tilde T^{(m)}_{\mu\nu}
+\Theta_{\mu\nu},
\]

其中

\[
\Theta_{\mu\nu}
=
-\frac{2}{\sqrt{|\tilde g|}}
\frac{\delta S_{\rm aux}}{\delta\tilde g^{\mu\nu}}.
\]

因此 action 化的真实目标是：

\[
\boxed{
\mathcal C_{\mu\nu}
=
\Theta_{\mu\nu}/M_P^2,
}
\]

而不是把某个裸的 \(\lambda_I E^I_{\mu\nu}\) 直接叫作 \(\mathcal C_{\mu\nu}\)。

## 2. 一个看似可行但不够完整的源项

若把 \(C_{\mu\nu}\) 暂时当作独立协变张量，不让它依赖 \(\tilde g\)，考虑

\[
S_{\rm source}
=
-\frac{M_P^2}{2}
\int d^dx\sqrt{|\tilde g|}
\,
\tilde g^{\mu\nu}C_{\mu\nu}.
\]

对 \(\tilde g^{\mu\nu}\) 变分，得到

\[
\Theta_{\mu\nu}
=
M_P^2
\left(
C_{\mu\nu}
-\frac12\tilde g_{\mu\nu}C^\alpha{}_\alpha
\right).
\]

如果再有

\[
C^\alpha{}_\alpha=0,
\]

则

\[
\Theta_{\mu\nu}=M_P^2 C_{\mu\nu}.
\]

这说明一个重要事实：

\[
\boxed{
\text{独立、无迹的 }C_{\mu\nu}
\text{ 可以作为 metric 方程中的应力源。}
}
\]

但是这还不是完整理论，因为还没有解释：

1. 为什么 \(C_{\mu\nu}\) 必须落在 \(E\) 子空间；
2. 为什么 \(C^\mu{}_\mu=0\)；
3. 为什么 \(\tilde\nabla^\mu C_{\mu\nu}=0\)；
4. 为什么 \(Q\to0\) 时 \(C_{\mu\nu}\to0\)。

这些必须由辅助场方程或约束给出。

## 3. 为什么不能直接令 \(C_{\mu\nu}=\lambda_I E^I_{\mu\nu}\)

如果直接把

\[
C_{\mu\nu}=\lambda_I E^I_{\mu\nu}
\]

代入上一节的 \(S_{\rm source}\)，则

\[
S_{\rm source}
=
-\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\,
\lambda_I \tilde g^{\mu\nu}E^I_{\mu\nu}.
\]

此时 \(E^I_{\mu\nu}\) 显含：

\[
\tilde g_{\mu\nu},
\qquad
u_\mu=\partial_\mu S,
\qquad
r_\mu=\tilde\nabla_\mu\ln\sqrt{\tilde\rho}.
\]

于是 metric 变分不再只给

\[
M_P^2\lambda_I E^I_{\mu\nu}.
\]

它还会给：

1. \(\delta E^I_{\alpha\beta}/\delta\tilde g^{\mu\nu}\) 产生的额外 metric stress；
2. 投影算子、内积、升降指标变化产生的额外项；
3. 若 \(Q\)-门控权重 \(\mu(Q)\) 出现在 action 中，还会产生高阶 metric 变分项。

同样，物质变分也会改变。因为 \(u_\mu,r_\mu\) 依赖 \(S,\tilde\rho\)，所以：

\[
\frac{\delta S_{\rm aux}}{\delta S}\neq0,
\qquad
\frac{\delta S_{\rm aux}}{\delta\tilde\rho}\neq0
\]

一般成立。

这会把物质方程改成：

\[
\mathcal E_S^{(m)}
+\mathcal E_S^{({\rm aux})}
=0,
\]

\[
\mathcal E_{\tilde\rho}^{(m)}
+\mathcal E_{\tilde\rho}^{({\rm aux})}
=0.
\]

所以原本由物质 action 给出的 Hamilton-Jacobi 方程和连续性方程不再保持原样。

这正是用户此前担心的问题：

\[
\boxed{
\text{若引力/辅助 action 显含 }u,r,\tilde\rho,S,
\text{ 它通常会破坏“物质变分单独给测地线”的结构。}
}
\]

## 4. 最小辅助 sector 的正确写法

因此更稳妥的 action 原型不是先写

\[
\mathcal C_{\mu\nu}=\lambda_I E^I_{\mu\nu}
\]

进 metric 方程，而是写：

\[
S_{\rm aux}
=
S_{\rm source}[\tilde g,C]
+S_{\rm state}[C,\lambda,\zeta,\alpha,\ldots;\tilde g,u,r,Q].
\]

其中：

\[
S_{\rm source}
=
-\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\,
\tilde g^{\mu\nu}C_{\mu\nu}.
\]

\(S_{\rm state}\) 用来约束或选择 \(C_{\mu\nu}\) 的状态。

一个形式上的最小原型为：

\[
S_{\rm state}
=
\frac{M_P^2}{2}
\int\sqrt{|\tilde g|}
\left[
\alpha\,\tilde g^{\mu\nu}C_{\mu\nu}
+\zeta^\nu\tilde\nabla^\mu C_{\mu\nu}
+\Lambda_\perp^{\mu\nu}
(\Pi_E^\perp C)_{\mu\nu}
+\mu(Q) C_{\mu\nu}C^{\mu\nu}
\right].
\]

这里各符号含义如下：

\[
\alpha:
\text{ 无迹约束乘子；}
\]

\[
\zeta^\nu:
\text{ 守恒约束乘子；}
\]

\[
\Lambda_\perp^{\mu\nu}:
\text{ 强制 }C_{\mu\nu}\in E\text{ 的乘子；}
\]

\[
\mu(Q):
\text{ }Q\to0\text{ 时把 }C_{\mu\nu}\text{ 压向 0 的门控权重。}
\]

形式上，对乘子变分给出：

\[
C^\mu{}_\mu=0,
\]

\[
\tilde\nabla^\mu C_{\mu\nu}=0,
\]

\[
\Pi_E^\perp C_{\mu\nu}=0.
\]

当

\[
\mu(Q)\to\infty
\quad(Q\to0)
\]

时，有限 action 倾向要求

\[
C_{\mu\nu}\to0.
\]

这正是 trace=0 投影方程需要的结构。

## 5. 这个原型的关键 caveat

上一节的 \(S_{\rm state}\) 仍然不是最终作用量。

因为它依赖

\[
\Pi_E[\tilde g,u,r],
\qquad
Q[\tilde g,\tilde\rho],
\]

所以对 metric 和 matter 变分时，\(S_{\rm state}\) 自己也会贡献有效能动量和物质方程项。

因此真正进入 Einstein 方程右边的不是裸 \(C_{\mu\nu}\)，而是

\[
\Theta_{\mu\nu}
=
-\frac{2}{\sqrt{|\tilde g|}}
\frac{\delta (S_{\rm source}+S_{\rm state})}
{\delta\tilde g^{\mu\nu}}.
\]

于是完整 metric 方程为：

\[
\tilde G_{\mu\nu}
=
\frac{1}{M_P^2}
\left(
\tilde T^{(m)}_{\mu\nu}
+\Theta_{\mu\nu}
\right).
\]

若要和 equation-first 方程完全一致，必须要求的是：

\[
\boxed{
\Theta_{\mu\nu}/M_P^2
\in E,
\qquad
\Theta^\mu{}_\mu=0,
\qquad
\tilde\nabla^\mu\Theta_{\mu\nu}=0,
}
\]

而不只是要求裸 \(C_{\mu\nu}\) 满足这些式子。

这就是 action 化的核心修正。

## 6. 保持测地线方程的条件

物质部分要继续给出：

\[
\tilde g^{\mu\nu}u_\mu u_\nu=m^2,
\]

\[
\tilde\nabla_\mu(\tilde\rho u^\mu)=0,
\]

从而推出：

\[
u^\mu\tilde\nabla_\mu u_\nu=0.
\]

如果 \(S_{\rm aux}\) 显含 \(S,\tilde\rho\)，则必须满足：

\[
\boxed{
\frac{\delta S_{\rm aux}}{\delta S}=0,
\qquad
\frac{\delta S_{\rm aux}}{\delta\tilde\rho}=0
\quad
\text{在目标物质壳上成立。}
}
\]

否则物质方程会被辅助 sector 改写，隐变量测地线解释不再自动成立。

这给 action 化路线带来一个强约束：

\[
\boxed{
\text{要么 }S_{\rm aux}\text{ 不直接依赖 }S,\tilde\rho；
\quad
\text{要么必须设计补偿项，使其物质变分在壳上抵消。}
}
\]

但若 \(S_{\rm aux}\) 完全不依赖 \(S,\tilde\rho\)，又很难自然生成 \(u,r\) 张成的 \(E\) 子空间。

这就是当前 action 化路线最大的理论张力。

## 7. 两个可行分支

### 分支 A：equation-first 保持

保持当前方程组：

\[
\Pi_E^\perp\mathcal R=0,
\qquad
\tilde g^{\mu\nu}\Pi_E\mathcal R_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
\]

此时不强求它来自普通作用量。优点是直接保留测地线解释和显式方程；缺点是理论规格较低，需要单独证明适定性和实验一致性。

### 分支 B：引入独立 shadow fields

引入独立辅助方向：

\[
U_\mu,\qquad R_\mu,
\]

用它们构造

\[
E_{\rm aux}
=
\mathrm{span}
\{
\tilde g_{\mu\nu},
U_\mu U_\nu,
R_\mu R_\nu,
U_{(\mu}R_{\nu)}
\}.
\]

这样 \(S_{\rm aux}\) 可以不直接变分 \(S,\tilde\rho\)，从而不破坏物质方程。

然后需要额外机制让某个分支满足：

\[
U_\mu\approx u_\mu,
\qquad
R_\mu\approx r_\mu.
\]

这个机制不能简单用硬约束

\[
U_\mu-u_\mu=0,
\qquad
R_\mu-r_\mu=0
\]

直接加入 action，因为那又会通过约束乘子改写 \(S,\tilde\rho\) 的变分。

所以 shadow-field 分支的难点是：

\[
\boxed{
\text{如何让 }U,R\text{ 与 }u,r\text{ 选择同一物理分支，同时不破坏物质测地线方程。}
}
\]

## 8. 本轮结论

本轮不是证明 action 化失败，而是明确了正确目标。

错误目标是：

\[
\lambda_I E^I_{\mu\nu}
\text{ 直接作为 action 中的裸源，并且仍等于 metric 方程右边。}
\]

正确目标是：

\[
\boxed{
\Theta_{\mu\nu}
=
-\frac{2}{\sqrt{|\tilde g|}}
\frac{\delta S_{\rm aux}}{\delta\tilde g^{\mu\nu}},
\quad
\Theta_{\mu\nu}/M_P^2
\text{ 满足投影、无迹、守恒和 }Q\to0\text{ 分支条件。}
}
\]

同时必须检查：

\[
\boxed{
\delta S_{\rm aux}/\delta S
\text{ 和 }
\delta S_{\rm aux}/\delta\tilde\rho
\text{ 是否破坏物质壳方程。}
}
\]

因此当前路线图应调整为：

1. 保留 trace=0 投影方程作为 equation-first 候选；
2. action 化时改为寻找能产生目标有效应力 \(\Theta_{\mu\nu}\) 的辅助 sector；
3. 若辅助 sector 显含 \(u,r\)，必须同时处理物质变分补偿；
4. 若不想破坏物质变分，则考虑独立 \(U,R\) 的 shadow-field 分支。

## 9. 下一步

下一步有两个互补任务：

1. 理论任务：
   写出 shadow-field 分支的最小方程组，检查 \(U,R\) 与 \(u,r\) 的同分支选择是否可以作为边界/正则性条件，而不是硬变分约束。
2. 数值任务：
   在 `n=384` 可信条纹分辨率上复查 trace=0 投影方程，确认 `n=96` 中 trace=0 优于 \(\mathsf q_E\)-trace 的结论没有被低分辨率误导。

如果 `n=384` 仍支持 trace=0，那么 action 化应优先围绕辅助有效应力 \(\Theta_{\mu\nu}\) 展开，而不是回到纯 \(f(R)\) 或裸 \(\lambda_I E^I_{\mu\nu}\) 作用量。
