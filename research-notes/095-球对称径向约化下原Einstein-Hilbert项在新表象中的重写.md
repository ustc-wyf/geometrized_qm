# 球对称径向约化下原 Einstein-Hilbert 项在新表象中的重写

时间：2026-04-30

## 目标

在当前 `m\neq 0, C=1` 混合投影分支下，考虑 `3+1d` 球对称径向扇区，直接把联络之差
\[
\Delta\Gamma^\rho{}_{\mu\nu}
:=
\tilde\Gamma^\rho{}_{\mu\nu}-\Gamma^\rho{}_{\mu\nu}
\]
代回
\[
R[g]-\tilde R[\tilde g]
\]
的关系式，得到
\[
\sqrt{-g}\,R[g]
\]
在新表象中的显式重写公式。

这里的“新表象”是指：把原表象中的 `g` 全部重写成由 `(\tilde g,\rho,S)` 确定的对象，而不再把 `g` 当作独立基本变量。

## 1. 球对称径向扇区的度规与变换

取四维球对称度规
\[
ds^2
=
h_{ab}(x)\,dx^a dx^b
+
\mathcal R(x)^2 d\Omega_2^2,
\qquad
x^a=(t,r).
\]

这里：

- `h_{ab}` 是 `t-r` 平面上的二维 Lorentz 度规；
- `\mathcal R(t,r)` 是面积半径（areal radius）；
- `d\Omega_2^2` 是单位二球度规。

在球对称径向问题里，Bohm 平面正好就是 `t-r` 平面，因此当前变换在这个扇区上退化成：

\[
\tilde h_{ab}=\Omega^2 h_{ab},
\qquad
\Omega^2:=\frac{X}{m^2},
\]

而角向部分不变：
\[
\tilde g_{\theta\theta}=g_{\theta\theta}=\mathcal R^2,
\qquad
\tilde g_{\phi\phi}=g_{\phi\phi}=\mathcal R^2\sin^2\theta.
\]

因此
\[
\tilde g_{\mu\nu}
=
\begin{pmatrix}
\Omega^2 h_{ab} & 0 \\
0 & \mathcal R^2 \gamma_{ij}
\end{pmatrix},
\qquad
g_{\mu\nu}
=
\begin{pmatrix}
h_{ab} & 0 \\
0 & \mathcal R^2 \gamma_{ij}
\end{pmatrix}.
\]

体积元关系是
\[
\sqrt{-\tilde g}=\Omega^2\sqrt{-g},
\qquad
\sqrt{-g}=\Omega^{-2}\sqrt{-\tilde g}
=
\frac{m^2}{X}\sqrt{-\tilde g}.
\]

## 2. 联络之差的非零分量

记 `D_a` 为二维基底度规 `h_{ab}` 的 Levi-Civita 导数，`\tilde D_a` 为 `\tilde h_{ab}` 的 Levi-Civita 导数。

标准球对称 warped-product 联络的非零分量为：

### 原表象 `g`
\[
\Gamma^a{}_{bc}={}^{(h)}\Gamma^a{}_{bc},
\]
\[
\Gamma^a{}_{ij}
=
-\mathcal R\,D^a\mathcal R\,\gamma_{ij},
\]
\[
\Gamma^i{}_{aj}
=
\frac{D_a\mathcal R}{\mathcal R}\delta^i_j,
\]
\[
\Gamma^i{}_{jk}={}^{(\gamma)}\Gamma^i{}_{jk}.
\]

### 新表象 `\tilde g`
\[
\tilde\Gamma^a{}_{bc}={}^{(\tilde h)}\Gamma^a{}_{bc},
\]
\[
\tilde\Gamma^a{}_{ij}
=
-\mathcal R\,\tilde D^a\mathcal R\,\gamma_{ij}
=
-\Omega^{-2}\mathcal R\,D^a\mathcal R\,\gamma_{ij},
\]
\[
\tilde\Gamma^i{}_{aj}
=
\frac{D_a\mathcal R}{\mathcal R}\delta^i_j,
\]
\[
\tilde\Gamma^i{}_{jk}={}^{(\gamma)}\Gamma^i{}_{jk}.
\]

因此联络之差只有两类非零分量：

### 2.1 径向基底内部

令
\[
\varphi:=\ln\Omega=\frac12\ln\frac{X}{m^2}.
\]

则二维共形变换的标准公式给出
\[
\boxed{
\Delta\Gamma^a{}_{bc}
=
\delta^a_b\,D_c\varphi
+
\delta^a_c\,D_b\varphi
-
h_{bc}D^a\varphi.
}
\]

等价地，也可写成
\[
\boxed{
\Delta\Gamma^a{}_{bc}
=
\frac12
\Bigl[
\delta^a_b\,D_c\ln\Omega^2
+
\delta^a_c\,D_b\ln\Omega^2
-
h_{bc}D^a\ln\Omega^2
\Bigr].
}
\]

### 2.2 角向块的抬升项

\[
\boxed{
\Delta\Gamma^a{}_{ij}
=
\left(1-\Omega^{-2}\right)\mathcal R\,D^a\mathcal R\,\gamma_{ij}.
}
\]

利用 `g_{ij}=\tilde g_{ij}=\mathcal R^2\gamma_{ij}`，也可写成
\[
\boxed{
\Delta\Gamma^a{}_{ij}
=
\left(1-\Omega^{-2}\right)
\frac{D^a\mathcal R}{\mathcal R}\,g_{ij}.
}
\]

其余分量为零：
\[
\Delta\Gamma^i{}_{aj}=0,
\qquad
\Delta\Gamma^i{}_{jk}=0.
\]

所以球对称径向约化下，联络之差已经非常简洁：  
它由“二维共形块的联络差”与“角向纤维的抬升修正”两部分组成。

## 3. 由联络差得到 Ricci 标量关系

直接用 Ricci 张量差公式
\[
\tilde R_{\mu\nu}
=
R_{\mu\nu}
+
\nabla_\rho \Delta\Gamma^\rho{}_{\mu\nu}
-
\nabla_\nu \Delta\Gamma^\rho{}_{\mu\rho}
+
\Delta\Gamma^\rho{}_{\rho\lambda}\Delta\Gamma^\lambda{}_{\mu\nu}
-
\Delta\Gamma^\rho{}_{\nu\lambda}\Delta\Gamma^\lambda{}_{\mu\rho}
\]
逐项代入上面的非零分量，最后得到与 warped-product 标准公式一致的标量关系：

\[
\boxed{
R[g]
=
\Omega^2\,\tilde R[\tilde g]
+
\Box_h \ln\Omega^2
+
\frac{2}{\mathcal R^2}(1-\Omega^2).
}
\]

这里 `\Box_h := D^a D_a` 是原二维基底度规 `h_{ab}` 的 d'Alembertian。

进一步利用二维共形关系
\[
\tilde\Box_{\tilde h}\phi=\Omega^{-2}\Box_h\phi
\]
可把它改写成
\[
\boxed{
R[g]
=
\Omega^2
\left[
\tilde R[\tilde g]
+
\tilde\Box_{\tilde h}\ln\Omega^2
\right]
+
\frac{2}{\mathcal R^2}(1-\Omega^2).
}
\]

## 4. 原 Einstein-Hilbert 密度在新表象中的封闭公式

把
\[
\sqrt{-g}=\Omega^{-2}\sqrt{-\tilde g}
\]
代入上式，得到

\[
\boxed{
\sqrt{-g}\,R[g]
=
\sqrt{-\tilde g}
\left[
\tilde R[\tilde g]
+
\tilde\Box_{\tilde h}\ln\Omega^2
+
\frac{2}{\mathcal R^2}\left(\Omega^{-2}-1\right)
\right].
}
\]

再代回
\[
\Omega^2=\frac{X}{m^2},
\qquad
\Omega^{-2}=\frac{m^2}{X}
\]
得到最适合当前项目直接引用的形式：

\[
\boxed{
\sqrt{-g}\,R[g]
=
\sqrt{-\tilde g}
\left[
\tilde R[\tilde g]
+
\tilde\Box_{\tilde h}\ln\frac{X}{m^2}
+
\frac{2}{\mathcal R^2}\left(\frac{m^2}{X}-1\right)
\right].
}
\]

因此差值公式是

\[
\boxed{
\sqrt{-g}\,R[g]
-
\sqrt{-\tilde g}\,\tilde R[\tilde g]
=
\sqrt{-\tilde g}
\left[
\tilde\Box_{\tilde h}\ln\frac{X}{m^2}
+
\frac{2}{\mathcal R^2}\left(\frac{m^2}{X}-1\right)
\right].
}
\]

## 5. 解释

这条公式说明：在球对称径向约化下，原表象 Einstein-Hilbert 项重写到新表象后，只比 `\sqrt{-\tilde g}\tilde R` 多出两部分：

1. 一个纯二维基底上的共形修正项
\[
\tilde\Box_{\tilde h}\ln\frac{X}{m^2};
\]
2. 一个来自角向球面的几何修正项
\[
\frac{2}{\mathcal R^2}\left(\frac{m^2}{X}-1\right).
\]

这比一般 `3+1d` 情形简洁得多。

## 6. 两个极限

### 6.1 若 `X/m^2 \to 1`

则
\[
\ln\frac{X}{m^2}\to 0,
\qquad
\frac{m^2}{X}-1\to 0,
\]
所以
\[
\sqrt{-g}R[g]
\to
\sqrt{-\tilde g}\tilde R[\tilde g].
\]

也就是说，在弱量子势、缓变包络极限下，球对称径向扇区里 `A` 与 `B` 的引力项渐进一致。

### 6.2 若 `\mathcal R \to \infty`

则角向修正项
\[
\frac{2}{\mathcal R^2}\left(\frac{m^2}{X}-1\right)
\]
自动衰减。  
远离星体时，差异主要由二维共形项
\[
\tilde\Box_{\tilde h}\ln\frac{X}{m^2}
\]
控制。

## 7. 结论

球对称径向约化下，联络之差确实有非常简洁的结构，并且原 Einstein-Hilbert 密度在新表象中的重写公式也很漂亮：

\[
\boxed{
\sqrt{-g}\,R[g]
=
\sqrt{-\tilde g}
\left[
\tilde R[\tilde g]
+
\tilde\Box_{\tilde h}\ln\frac{X}{m^2}
+
\frac{2}{\mathcal R^2}\left(\frac{m^2}{X}-1\right)
\right].
}
\]

这说明对星体外部的径向问题，`A` 支原引力项与 `B` 支 Einstein-Hilbert 项之间的差异，已被压缩成：

- 一个二维基底上的共形 Laplacian 项；
- 一个角向曲率半径控制的球面修正项。

这正是后续研究球对称外部径向传播时最自然的起点。
