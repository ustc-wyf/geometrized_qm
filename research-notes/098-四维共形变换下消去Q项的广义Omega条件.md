# 四维共形变换下消去 `Q` 项的广义 `\Omega` 条件

## 目标

从 `A` 支完整拉氏量

\[
\mathcal L_A
=
\sqrt{-g}\left[
\frac{M_P^2}{2}R[g]
+
\rho\,(X-m^2-Q)
\right],
\qquad
X:=g^{\mu\nu}\partial_\mu S\,\partial_\nu S,
\qquad
Q:=\frac{\Box_g\sqrt{\rho}}{\sqrt{\rho}},
\]

出发，考虑四维纯共形变换

\[
\tilde g_{\mu\nu}=\Omega(\phi)^2 g_{\mu\nu},
\qquad
\phi:=\sqrt{\rho},
\]

并寻找什么样的 `\Omega(\phi)` 能在新表象里**消去** `Q` 型二阶导数项。

为避免记号堆叠，下文取

\[
\tilde\phi:=\phi,\qquad \tilde S:=S,
\]

只把度规换到 `\tilde g` 表象。

---

## 一、基本变换公式

在四维里：

\[
g_{\mu\nu}=\Omega^{-2}\tilde g_{\mu\nu},
\qquad
g^{\mu\nu}=\Omega^2\tilde g^{\mu\nu},
\qquad
\sqrt{-g}=\Omega^{-4}\sqrt{-\tilde g}.
\]

因此

\[
X=\Omega^2 \tilde X,
\qquad
\tilde X:=\tilde g^{\mu\nu}\partial_\mu \tilde S\,\partial_\nu \tilde S.
\]

对任意标量 `f`，四维 d'Alembert 算子满足

\[
\Box_g f
=
\Omega^2\left[
\tilde\Box f
-2\,\tilde g^{\mu\nu}(\partial_\mu\ln\Omega)\partial_\nu f
\right].
\]

令 `f=\phi`，并用

\[
\partial_\mu\ln\Omega
=
\frac{\Omega_\phi}{\Omega}\,\partial_\mu \phi,
\qquad
\Omega_\phi:=\frac{d\Omega}{d\phi},
\]

则

\[
Q
=
\frac{\Box_g\phi}{\phi}
=
\Omega^2\left[
\tilde Q_\phi
-2\frac{\Omega_\phi}{\Omega}\frac{(\tilde\nabla\phi)^2}{\phi}
\right],
\qquad
\tilde Q_\phi:=\frac{\tilde\Box\phi}{\phi}.
\]

---

## 二、物质部分在新表象中的一般形式

### 1. 动能项

\[
\sqrt{-g}\,\rho X
=
\Omega^{-4}\sqrt{-\tilde g}\,\phi^2\,(\Omega^2\tilde X)
=
\sqrt{-\tilde g}\,\phi^2\Omega^{-2}\tilde X.
\]

### 2. 质量项

\[
\sqrt{-g}\,\rho m^2
=
\sqrt{-\tilde g}\,\phi^2 m^2 \Omega^{-4}.
\]

### 3. `Q` 项

\[
-\sqrt{-g}\,\rho Q
=
-\Omega^{-4}\sqrt{-\tilde g}\,\phi^2\,
\Omega^2\left[
\tilde Q_\phi
-2\frac{\Omega_\phi}{\Omega}\frac{(\tilde\nabla\phi)^2}{\phi}
\right]
\]
\[
=
\sqrt{-\tilde g}
\left[
-\phi^2\Omega^{-2}\tilde Q_\phi
+2\phi\,\Omega^{-3}\Omega_\phi (\tilde\nabla\phi)^2
\right].
\]

因此物质部分整体为

\[
\boxed{
\mathcal L_{m}
=
\sqrt{-\tilde g}\left[
\phi^2\Omega^{-2}\tilde X
-\phi^2 m^2 \Omega^{-4}
-\phi^2\Omega^{-2}\tilde Q_\phi
+2\phi\,\Omega^{-3}\Omega_\phi (\tilde\nabla\phi)^2
\right].
}
\]

---

## 三、Einstein-Hilbert 项在新表象中的一般形式

四维共形变换下：

\[
\sqrt{-g}R[g]
=
\sqrt{-\tilde g}\left[
\Omega^{-2}\tilde R
+6\Omega^{-2}\tilde\Box\ln\Omega
-6\Omega^{-2}(\tilde\nabla\ln\Omega)^2
\right].
\]

将 `\ln\Omega` 对 `\phi` 展开，可得

\[
\sqrt{-g}R[g]
=
\sqrt{-\tilde g}\left[
\Omega^{-2}\tilde R
+6\Omega^{-3}\Omega_\phi\,\tilde\Box\phi
+6\left(
\Omega^{-3}\Omega_{\phi\phi}
-2\Omega^{-4}\Omega_\phi^2
\right)(\tilde\nabla\phi)^2
\right],
\]

其中

\[
\Omega_{\phi\phi}:=\frac{d^2\Omega}{d\phi^2}.
\]

利用 `\tilde\Box\phi=\phi \tilde Q_\phi`，可写成

\[
\boxed{
\mathcal L_{\rm EH}
=
\sqrt{-\tilde g}\left[
\frac{M_P^2}{2}\Omega^{-2}\tilde R
+3M_P^2\,\phi\,\Omega^{-3}\Omega_\phi\,\tilde Q_\phi
+3M_P^2\left(
\Omega^{-3}\Omega_{\phi\phi}
-2\Omega^{-4}\Omega_\phi^2
\right)(\tilde\nabla\phi)^2
\right].
}
\]

---

## 四、总拉氏量与消去 `Q` 的条件

把物质项与 EH 项相加，`Q` 型项 `\tilde Q_\phi` 的总系数为

\[
-\phi^2\Omega^{-2}
+3M_P^2\,\phi\,\Omega^{-3}\Omega_\phi.
\]

要它对任意场构型恒等消失，必须满足

\[
-\phi^2\Omega^{-2}
+3M_P^2\,\phi\,\Omega^{-3}\Omega_\phi
=
0.
\]

乘以 `\Omega^3/\phi` 后得到一阶常微分方程

\[
\boxed{
3M_P^2\,\Omega_\phi=\phi\,\Omega.
}
\]

解为

\[
\boxed{
\Omega(\phi)=\Omega_0\,\exp\!\left(\frac{\phi^2}{6M_P^2}\right)
=
\Omega_0\,\exp\!\left(\frac{\rho}{6M_P^2}\right).
}
\]

因此：

\[
\boxed{
\text{若要在四维纯共形变换下把 `Q` 项从全拉氏量中全局消去，}
\Omega\propto \sqrt{\rho}\text{ 不是正确答案；正确答案是指数型 }\Omega(\rho)\sim e^{\rho/(6M_P^2)}.
}
\]

---

## 五、代回后得到的简化全拉氏量

在满足

\[
3M_P^2\Omega_\phi=\phi\Omega
\]

的条件下，`Q` 项完全消失。此时梯度平方项的总系数可化简为

\[
\Omega^{-2}\left(1+\frac{\phi^2}{3M_P^2}\right).
\]

所以总拉氏量变成

\[
\boxed{
\mathcal L_A
=
\sqrt{-\tilde g}\left[
\frac{M_P^2}{2}\Omega^{-2}\tilde R
+\phi^2\Omega^{-2}\tilde X
-m^2\phi^2\Omega^{-4}
+\Omega^{-2}\left(1+\frac{\phi^2}{3M_P^2}\right)(\tilde\nabla\phi)^2
\right].
}
\]

其中

\[
\Omega(\phi)=\Omega_0\,\exp\!\left(\frac{\phi^2}{6M_P^2}\right).
\]

如果想把它完全改写成“全带 tilde”的变量，只需再定义

\[
\tilde\phi:=\phi,\qquad \tilde S:=S,
\]

并把 `\phi` 全部替换成 `\tilde\phi` 即可。

---

## 六、与 `\Omega\propto\sqrt{\rho}` 的比较

若取

\[
\Omega\propto\sqrt{\rho}=\phi,
\]

则虽然有某些项会变得漂亮（例如动能项可常数化），但 `Q` 项不会全局消失，最终仍留下

\[
(\text{常数}-3M_P^2 \times \text{场})\,\tilde Q
\]

这类系数。

所以：

\[
\boxed{
\Omega\propto\sqrt{\rho}
\]
\[
\text{适合做“部分几何化”，}
\qquad
\Omega\sim e^{\rho/(6M_P^2)}
\text{ 才能做“全局消去 }Q\text{”。}
}

---

## 七、解释

这个结果有两个重要含义：

1. **如果目标是把 `Q` 完全吸收到引力侧**，最自然的共形因子不是幂律型，而是指数型。
2. `Q` 虽然消失了，但代价不是零：
   - Einstein-Hilbert 项前面出现了场依赖因子 `\Omega^{-2}`;
   - 标量场的梯度平方项也留下了新的系数；
   - 所以它仍然不是“普通最小耦合标量 + 普通 Einstein-Hilbert”，而是一个新的非最小耦合标量-张量理论。
