# 指数型 `\Omega` 消去 `Q` 后的标量-张量形式

## 1. 起点

上一节已经得到：若取四维共形变换

\[
\tilde g_{\mu\nu}=\Omega(\phi)^2 g_{\mu\nu},
\qquad
\phi:=\sqrt{\rho},
\]

并要求原 `A` 支拉氏量中的

\[
Q=\frac{\Box_g\phi}{\phi}
\]

在新表象里全局消去，则唯一的一般解是

\[
\boxed{
\Omega(\phi)=\Omega_0\exp\!\left(\frac{\phi^2}{6M_P^2}\right).
}
\]

此时完整拉氏量可写成

\[
\boxed{
\mathcal L
=
\sqrt{-\tilde g}\left[
\frac{M_P^2}{2}\Omega^{-2}\tilde R
+\phi^2\Omega^{-2}\tilde X
-m^2\phi^2\Omega^{-4}
+\Omega^{-2}\left(1+\frac{\phi^2}{3M_P^2}\right)(\tilde\nabla\phi)^2
\right],
}
\]

其中

\[
\tilde X:=\tilde g^{\mu\nu}\partial_\mu S\,\partial_\nu S.
\]

目标是把它整理成更标准的“标量-张量”形式。

---

## 2. 定义非最小耦合标量

定义

\[
\boxed{
\Phi:=M_P^2\Omega^{-2}
=
\Phi_0\,\exp\!\left(-\frac{\phi^2}{3M_P^2}\right),
\qquad
\Phi_0:=M_P^2\Omega_0^{-2}.
}
\]

于是

\[
\frac{M_P^2}{2}\Omega^{-2}\tilde R=\frac12 \Phi\,\tilde R.
\]

反过来

\[
\boxed{
\phi^2 = -3M_P^2 \ln\!\frac{\Phi}{\Phi_0}.
}
\]

并且

\[
\partial_\mu \Phi
=
-\frac{2\phi}{3}\Omega^{-2}\partial_\mu\phi
=
-\frac{2\phi}{3M_P^2}\Phi\,\partial_\mu\phi.
}
\]

因此

\[
(\tilde\nabla\phi)^2
=
\frac{9M_P^4}{4\phi^2\Phi^2}(\tilde\nabla\Phi)^2.
\]

---

## 3. 梯度项写成 \(\Phi\) 的函数

原梯度项是

\[
\Omega^{-2}\left(1+\frac{\phi^2}{3M_P^2}\right)(\tilde\nabla\phi)^2
=
\frac{\Phi}{M_P^2}\left(1+\frac{\phi^2}{3M_P^2}\right)(\tilde\nabla\phi)^2.
\]

代入上式得

\[
\boxed{
\Omega^{-2}\left(1+\frac{\phi^2}{3M_P^2}\right)(\tilde\nabla\phi)^2
=
\frac{3M_P^2}{4\Phi}\,
\frac{1-\ln(\Phi/\Phi_0)}{-\ln(\Phi/\Phi_0)}
(\tilde\nabla\Phi)^2.
}
\]

更紧凑地，若记

\[
L:= -\ln(\Phi/\Phi_0)>0,
\]

则

\[
\boxed{
\Omega^{-2}\left(1+\frac{\phi^2}{3M_P^2}\right)(\tilde\nabla\phi)^2
=
\frac{3M_P^2}{4\Phi}\frac{1+L}{L}(\tilde\nabla\Phi)^2.
}
\]

---

## 4. 相位 \(S\) 的动能耦合

\[
\phi^2\Omega^{-2}\tilde X
=
\frac{\phi^2}{M_P^2}\Phi\,\tilde X
=
-3\Phi\ln\!\frac{\Phi}{\Phi_0}\,\tilde X.
}
\]

因此

\[
\boxed{
\phi^2\Omega^{-2}\tilde X
=
-3\Phi\ln\!\frac{\Phi}{\Phi_0}\,
\tilde g^{\mu\nu}\partial_\mu S\partial_\nu S.
}
\]

这说明 `S` 不再是最小耦合自由 KG 相位，而是带有一个场依赖前因子。

---

## 5. 质量势项

\[
-m^2\phi^2\Omega^{-4}
=
-m^2\phi^2\frac{\Phi^2}{M_P^4}
=
\frac{3m^2}{M_P^2}\Phi^2\ln\!\frac{\Phi}{\Phi_0}.
}
\]

所以可以定义有效势

\[
\boxed{
U(\Phi):=
-\frac{3m^2}{M_P^2}\Phi^2\ln\!\frac{\Phi}{\Phi_0},
}
\]

使得拉氏量里写成 `-U(\Phi)` 的标准符号。

由于 `\ln(\Phi/\Phi_0)\le 0`，这个 `U(\Phi)` 在允许区间内非负。

---

## 6. 标准标量-张量形式

把上面各项合并，得到

\[
\boxed{
\mathcal L
=
\sqrt{-\tilde g}\Bigg[
\frac12 \Phi\,\tilde R
+A(\Phi)\,\tilde g^{\mu\nu}\partial_\mu S\partial_\nu S
+B(\Phi)\,(\tilde\nabla\Phi)^2
-U(\Phi)
\Bigg],
}
\]

其中

\[
\boxed{
A(\Phi)=-3\Phi\ln\!\frac{\Phi}{\Phi_0},
}
\]

\[
\boxed{
B(\Phi)=\frac{3M_P^2}{4\Phi}\frac{1+L}{L},
\qquad
L:=-\ln(\Phi/\Phi_0),
}
\]

\[
\boxed{
U(\Phi)= -\frac{3m^2}{M_P^2}\Phi^2\ln\!\frac{\Phi}{\Phi_0}.
}
\]

这已经是一个非常标准的“非最小耦合标量-张量 + 另一标量 `S` 带场依赖动能前因子”的理论。

---

## 7. 它和通常 Brans-Dicke 形式的关系

通常 Jordan-frame 标量-张量理论写成

\[
\mathcal L_{\rm ST}
=
\sqrt{-\tilde g}
\left[
\frac12\Phi \tilde R
-\frac{\omega(\Phi)}{2\Phi}(\tilde\nabla\Phi)^2
-U(\Phi)
\right]
+\mathcal L_m.
\]

若硬要比照，则我们的梯度项对应

\[
-\frac{\omega(\Phi)}{2\Phi}=B(\Phi),
}
\]

因此

\[
\boxed{
\omega(\Phi)
=
-2\Phi B(\Phi)
=
-\frac{3M_P^2}{2}\frac{1+L}{L}.
}
\]

所以它并不是通常温和的 Brans-Dicke 参数，而是一个随 `\Phi` 变化、并且在这里取负的大函数。

这再次说明：这个理论不是普通 Brans-Dicke，而是一个较特殊的非最小耦合模型。

---

## 8. 与原来的 `B/C` 分支的关系

### 与 `B` 支

`B` 支的核心是把质量壳项几何化为

\[
\tilde g^{\mu\nu}u_\mu u_\nu-m^2,
}
\]

并把引力部分保持成纯 Einstein-Hilbert：

\[
\sqrt{-\tilde g}\,\tilde R.
}
\]

而这里的新理论：

- 的确把原来的 `Q` 消掉了；
- 但代价是引力项变成
  \[
  \frac12\Phi \tilde R
  \]
  的非最小耦合；
- 同时 `S` 的动能项也获得了 `A(\Phi)` 前因子。

所以：

\[
\boxed{
\text{它不像 B 支；它不是“纯 Einstein-Hilbert + 最小耦合物质”的 }\tilde g\text{ 理论。}
}
\]

### 与 `C` 支

`C` 支本身就是 `f(\tilde R)` 型修正重力；  
而这里得到的是

\[
\Phi \tilde R + B(\Phi)(\nabla\Phi)^2 - U(\Phi)
}
\]

这种标量-张量型修正重力。

因此：

\[
\boxed{
\text{它在结构上更接近标量-张量 / Jordan-frame 理论，
而不是我们当前的 B 支，也不是当前那个纯 }f(\tilde R)\text{ 型 C 支。}
}
\]

---

## 9. 总结

1. 放宽 `\Omega` 以后，确实可以找到一个精确消去 `Q` 的共形因子：

\[
\boxed{
\Omega(\rho)=\Omega_0 e^{\rho/(6M_P^2)}.
}
\]

2. 在这个选择下，原来的全拉氏量变成一个新的标量-张量理论：

\[
\boxed{
\mathcal L
=
\sqrt{-\tilde g}\Big[
\frac12\Phi\tilde R
+A(\Phi)(\partial S)^2
+B(\Phi)(\partial\Phi)^2
-U(\Phi)
\Big].
}
\]

3. 因此：
   - `Q` 被消掉了；
   - 但它不是无代价的；
   - 代价是引入了非最小耦合引力和新的梯度/势结构。

4. 从分类上说，这条路更像“把原 Bohm 量子势理论改写成一个特定 Jordan-frame 标量-张量理论”，而不是得到我们原来想要的 `B` 支那种“纯 \(\tilde g\)-Einstein + 最小耦合物质”形式。
