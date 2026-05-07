# D 支 subcell 界面 toy 演化首轮

时间：2026-05-01

## 目的

用户要求在认识到过渡层可作为动态界面后，开始利用这种性质重新组织全动力学数值计算。

本轮实现第一版 `subcell/body-fitted interface toy evolution`。

这里的目标不是宣称已完成完整物理演化，而是先检验：

```text
动态界面 y=ell^2 R_tilde=±1 是否能作为子网格线源进入一个闭合数值校正器。
```

## 脚本

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_subcell_interface_evolution.py
```

## 变量定义

### y

\[
y=\ell^2\tilde R.
\]

对 D 支：

\[
f_R=\operatorname{sech}^2(y).
\]

### \(\chi\)

\[
\chi=\tilde R.
\]

本轮 toy evolution 的未知校正量是：

\[
\delta\chi.
\]

### trace residual

记作 `trace_residual`，是 D 支标量-张量形式在 `2+1d` 下的迹方程残差：

\[
\mathcal T_{\rm ref}
=
f_R\chi-\frac{3}{2}f+2\tilde\Box f_R-\frac{\tilde T}{M_P^2}.
\]

当前计算中固定参考度规与物质迹，只校正 \(\chi\)。因此它仍是 frozen-background toy，不是完整演化。

### subcell interface

界面为：

\[
y=+1,\qquad y=-1.
\]

脚本在每个穿越单元内用 marching-square 型线段表示界面，而不是把界面膨胀成网格厚带。

### line source

每条界面线段有：

- 中点 `(x,z)`;
- 线段长度 `length`;
- 法向短线积分诊断量；
- 线源强度。

本轮最有效的线源模式是 `source_mode=seed`：

```text
先用局部 IMEX 估计 delta_chi_seed，
再把 delta_chi_seed 采样到界面线段上作为线源方向。
```

线源用 cloud-in-cell 方式沉积到网格：

\[
S_\Gamma \sim \sum_\Gamma q_\Gamma\,\frac{\Delta s_\Gamma}{\Delta x\Delta z}.
\]

### toy equation

第一版用 screened-Poisson 型隐式校正：

\[
(1-L^2\Delta)\delta\chi=S_\Gamma.
\]

其中 `L=screen_length`。

当前实现为了快速闭环，使用周期谱解法；由于波包主支撑区远离边界，这一步只作为原型可以接受。后续若转正式演化，需要改成与边界条件一致的椭圆/IMEX solve。

## smoke test: 96x96

命令：

```text
python3 kg_examples/prototype_d_subcell_interface_evolution.py \
  --output visualizations/d_subcell_interface_toy_smoke_96_seed_local \
  --ell 30 --resolution 96 --time 16 --levels 1 \
  --samples 31 --screen-length 0.15 \
  --target-delta-p95 0.0002 --source-mode seed \
  --amplitudes=-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8
```

输出：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_subcell_interface_toy_smoke_96_seed_local/summary.json
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_subcell_interface_toy_smoke_96_seed_local/d_subcell_interface_toy_ell30_n96_t16.png
```

结果：

- support: `930`
- segments: `649`
- chosen amplitude: `-2`
- before abs_max: `3.86e5`
- after abs_max: `3.22e5`
- before p95: `0.492`
- after p95: `0.893`

解释：

界面线源能压低最极端尖峰，但会提高 p95。

## 128x128 复核

命令：

```text
python3 kg_examples/prototype_d_subcell_interface_evolution.py \
  --output visualizations/d_subcell_interface_toy_ell30_t16_128_seed_local \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --samples 41 --screen-length 0.15 \
  --target-delta-p95 0.0002 --source-mode seed \
  --amplitudes=-8,-4,-2,-1,-0.5,0,0.5,1,2,4,8
```

输出：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_subcell_interface_toy_ell30_t16_128_seed_local/summary.json
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_subcell_interface_toy_ell30_t16_128_seed_local/d_subcell_interface_toy_ell30_n128_t16.png
```

结果：

- support: `1669`
- segments: `888`
- chosen amplitude: `1`
- before abs_max: `3.62e5`
- after abs_max: `1.63e5`
- before abs_mean: `362.1`
- after abs_mean: `338.0`
- before p95: `1.014`
- after p95: `1.177`

解释：

`128x128` 下趋势更清楚：

```text
界面子网格校正确实能显著压低最极端尖峰；
但当前 frozen-background scalar chi toy 会把一部分误差扩散到更广区域，
使 p95 略升。
```

这说明 subcell interface 路线有数值价值，但第一版 toy 还不够成为完整演化器。

## 当前判断

### 已经成立

1. 动态界面 `y=±1` 可以被提取为子网格线段。
2. 线段源项可以沉积进网格并进入隐式校正方程。
3. 该校正能在 `96/128` 两个分辨率下压低最极端 trace residual。

### 尚未成立

1. 它还没有同时降低 bulk p95。
2. 它还没有推进 \((\rho,S)\)。
3. 它还没有更新完整 \(\tilde g\) 或约束系统。
4. 它仍使用 frozen metric/stress trace 和周期谱 screened solve。

## 结论

本轮说明：

\[
\boxed{
\text{subcell interface 方法可进入数值器，但不能只靠单个 }\delta\chi\text{ 校正完成全动力学。}
}
\]

下一步应把这个 toy 从“只改 \(\chi\)”升级为：

1. 界面源项只处理尖峰；
2. bulk 区域另用局部 IMEX/Newton 校正处理 p95；
3. saturated bulk 使用固定参考延拓或最小曲率延拓；
4. 最后把这三块耦合成一小步真正的 \((\rho,S,\tilde g)\) 演化。

