# 045 为什么 A 分支数值器无法贴近当前 benchmark

## 结论

当前用户认可的 `current_geometry_check` / `double_slit_beamlike_xz` 图样，
不是平直背景下自由 KG 全动力学初边值问题的精确解。

因此：

```text
任何试图用“全动力学时间演化器”去逐步逼近这张 benchmark 图的做法，
都不会稳定收敛到它。
```

问题不只是边界条件或网格，而是：

```text
benchmark 本身属于准单色束流/近似 beam 模型，
不是精确的 flat KG 解。
```

## 1. 检验方式

若

```text
ψ(x,z,t) = φ(x,z) e^{-iωt}
```

是平直 KG 的精确单频解，则 `φ` 必须满足：

```text
(∂_x^2 + ∂_z^2 + (ω^2 - m^2)) φ = 0
```

我把当前 `kg_double_slit_beamlike.py` 产生的 `ψ` 代回这条方程，
直接计算残差。

## 2. 数值结果

在当前主支撑区，
得到：

```text
abs residual mean main    ≈ 0.686
relative residual mean    ≈ 3.07e-2
relative residual max     ≈ 4.37
```

这说明：

```text
benchmark 与平直 KG 的单频 Helmholtz 方程只有近似相容，
而不是精确相容。
```

## 3. 物理与数值含义

这件事直接解释了前面的现象：

```text
1. 若把 benchmark 当作一个“实验上想看到的束流图样”，它是合理的。
2. 若把 benchmark 当作 A 分支全动力学 PDE 的精确目标解，它就不合理。
```

所以之前“让 A 分支数值器贴近 benchmark”的要求，
若坚持使用平直 KG 全动力学时间推进器，
在原则上就会卡住。

## 4. 后续正确路线

若要让三分支比较可信，有两条严格可行的路线：

### 路线 A

把 A 分支 benchmark 换成真正的精确 KG 解，
例如严格 Fourier 构造出来的解。

### 路线 B

把整个问题改写成与当前 benchmark 同类的稳态边界驱动问题，
也就是：

```text
不再拿“瞬态 Cauchy 演化”去逼近它，
而改求稳态/单频边界值问题。
```

对当前项目而言，
如果用户坚持使用当前这张几何上更符合直觉的 benchmark 图样，
那么更自然的是走路线 B。
