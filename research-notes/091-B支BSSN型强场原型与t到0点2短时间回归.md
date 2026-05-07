# B支 `BSSN` 型强场原型与 `t=0.2` 短时间回归

日期：2026-04-30

## 目标

在当前

- 局域双高斯波包初值
- `x,z\in[-20,20]`
- `B` 支强场几何演化

设定下，放弃先前的

- 原始 `ADM`
- 显式 `RK4`
- 零 `shift`
- 仅 `1+log` 切片

组合，改为一套 `BSSN` 型共形-无迹变量原型，并先检查它在短时间窗内是否能够重现旧 `B` 支的弱场早期行为。

## 新变量

新原型脚本：

`/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_b_bssn_driver.py`

以三维空间度规
\[
\gamma_{ij}=\mathrm{diag\ block}(h_{xx},h_{xz},h_{zz};e^{2\beta})
\]
为基础，引入
\[
\chi := (\det\gamma)^{-1/3},
\]
并定义二维共形块
\[
\tilde\gamma_{xx}=\chi h_{xx},\qquad
\tilde\gamma_{xz}=\chi h_{xz},\qquad
\tilde\gamma_{zz}=\chi h_{zz}.
\]

由于
\[
\det\tilde\gamma_{(2)}\cdot \tilde\gamma_{yy}=1,
\]
所以 `yy` 共形分量由二维共形块行列式代数决定，不作为独立变量推进。

外曲率改写为

- `K`：三维迹
- `\tilde A_{xx},\tilde A_{xz},\tilde A_{zz}`：共形无迹块

再加上

- `N`：lapse
- `\beta^x,\beta^z`：shift
- `n,S`：物质守恒密度与相位

共同组成当前原型状态变量。

## 规范条件

当前强场原型使用：

### 切片条件
\[
\partial_t N + \beta^i\partial_i N = -2 N K.
\]

### `shift` 驱动

不是完整的二阶 `Gamma-driver`，而是当前最小版本
\[
\partial_t \beta^i = \mu_\beta \tilde\Gamma^i - \eta_\beta \beta^i,
\]
其中
\[
\tilde\Gamma^x = -\partial_x \tilde\gamma^{xx}-\partial_z\tilde\gamma^{xz},
\qquad
\tilde\Gamma^z = -\partial_x \tilde\gamma^{xz}-\partial_z\tilde\gamma^{zz}.
\]

本轮测试用参数：
\[
\mu_\beta=0.5,\qquad \eta_\beta=0.5.
\]

## 物质与边界

物质部分继续使用原来的守恒系统
\[
\partial_t n+\partial_i(n v^i)=0,
\qquad
\partial_t S = \beta^i\partial_i S - N E.
\]

但边界条件不再沿用旧物理变量，而是在每一步先由 `BSSN` 状态重构物理几何，再按特征入流/出流条件更新 `n,S`。

几何边界条件则取弱场背景：
\[
\chi\to1,\quad
\tilde\gamma_{xx}\to1,\quad
\tilde\gamma_{xz}\to0,\quad
\tilde\gamma_{zz}\to1,
\]
\[
K\to0,\quad
\tilde A_{xx},\tilde A_{xz},\tilde A_{zz}\to0,\quad
N\to1,\quad
\beta^x,\beta^z\to0.
\]

## `t=0.2` 短时间回归

输出目录：

`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/b_bssn_driver_t02/`

运行参数：
\[
\text{steps}=800,\qquad dt=2.5\times10^{-4}.
\]

终态量级：
\[
\chi_{\min}\approx 0.9992891845,
\]
\[
N_{\min}\approx 0.9999883248,\qquad
N_{\max}\approx 1.0000012897,
\]
\[
\max|h_{xx}-1|\approx 1.0760\times10^{-3},
\]
\[
\max|h_{xz}|\approx 6.4893\times10^{-5},
\]
\[
\max|h_{zz}-1|\approx 1.0764\times10^{-3},
\]
\[
\max|\beta_{\text{geom}}|\approx 1.10\times10^{-7}.
\]

这些量级与旧 `B` 支在相同时间窗下的弱场结果保持一致，说明：

1. 新的 `BSSN` 型变量没有在早期引入明显的物理偏差；
2. 非零 `shift` 的引入目前仍保持在极小量级，
\[
\max|\beta^i| \sim 10^{-5},
\]
没有在早期破坏解；
3. 当前原型至少在短时间弱场区是自洽的。

## 当前结论

可以把这一步概括为：

\[
\boxed{
\text{B支的BSSN型强场原型已经建立，并已通过 }t=0.2\text{ 的短时间回归。}
}
\]

下一步不再停留在短时间窗，而是直接进行带检查点的长时间诊断，目标是检查它能否跨过旧方法在
\[
t_B\approx 6.84475
\]
附近的失稳时刻。
