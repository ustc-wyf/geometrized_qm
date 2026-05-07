## D支 tensor jump 首轮诊断

日期：2026-05-02

### 目的

上一轮 moving-interface 原型只用了 trace/scalaron jump law：

\[
2\,\tilde g^{ab}F_aF_b\,\Delta(\partial_q f_R)
+\int_{\rm layer}
\left(f_R\tilde R-\frac32 f-\frac{\tilde T}{M_P^2}\right)dq=0.
\]

这只能决定界面速度的一个标量条件，不能保证完整张量场方程成立。本轮把 D 支 metric \(f(R)\) 方程的薄层主部提升到张量形式。

### 张量薄层主部

D 支场方程：

\[
f_R\tilde R_{ab}
-\frac12 f\tilde g_{ab}
-\tilde\nabla_a\tilde\nabla_b f_R
+\tilde g_{ab}\tilde\Box f_R
-\frac{\tilde T_{ab}}{M_P^2}=0.
\]

设界面坐标

\[
q=\ell^2\tilde R,\qquad F_a=\partial_a q.
\]

保留 \(f_R(q)\) 的主奇异二阶导数项，跨层积分得到首个张量 jump 诊断：

\[
\Delta(\partial_q f_R)
\left(
-F_aF_b+\tilde g_{ab}\tilde g^{cd}F_cF_d
\right)
+
\int_{\rm layer}
\left(
f_R\tilde R_{ab}
-\frac12 f\tilde g_{ab}
-\frac{\tilde T_{ab}}{M_P^2}
\right)dq
=0.
\]

取迹后回到前一轮 trace jump law。

### 新脚本

`kg_examples/diagnose_d_tensor_interface_jump.py`

关键处理：

- 界面仍取 `raw_y=ell^2 R_tilde=±1`；
- 在界面上插值 \(\tilde g_{ab}\)，再现场求逆，避免分别插值 \(\tilde g_{ab}\) 和 \(\tilde g^{ab}\) 导致迹条件被污染；
- Ricci 张量的无迹部分冻结，迹部强制按 \(R(q)=q/\ell^2\) 通过层变化；
- 先用 trace 条件解 \(F_t\)，再检查完整张量残差；
- 额外做 free-covector rank-one 兼容性检查。

### rank-one 兼容性含义

若目标张量

\[
H_{ab}:=-\frac{1}{\Delta(\partial_q f_R)}
\int_{\rm layer}(\cdots)_{ab}\,dq
\]

真的能由某个界面协向量 \(F_a\) 产生，则

\[
H_{ab}=-F_aF_b+\tilde g_{ab}F^2.
\]

在 \(2+1d\) 中有

\[
F^2=\frac12\operatorname{Tr}_{\tilde g}H,
\]

于是

\[
K_{ab}:=\tilde g_{ab}F^2-H_{ab}
\]

必须是 rank-one，并等于 \(F_aF_b\)。本轮用 `rank1_tail_relative` 度量偏离 rank-one 的程度。

这个检查不固定当前参考界面法向，因此它问的是“是否存在某个合适界面法向”，比 trace-only 速度律更接近真正的 tensor matching。

### 运行

`96x96`：

```bash
python3 kg_examples/diagnose_d_tensor_interface_jump.py \
  --output visualizations/d_tensor_interface_jump_ell30_t16_96_rankcheck_v2 \
  --ell 30 --resolution 96 --time 16 --levels 1 \
  --half-width-y 1 --samples 61
```

`128x128`：

```bash
python3 kg_examples/diagnose_d_tensor_interface_jump.py \
  --output visualizations/d_tensor_interface_jump_ell30_t16_128_direction \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81
```

输出图：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_interface_jump_ell30_t16_96_rankcheck_v2/d_tensor_interface_jump_ell30_n96_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_interface_jump_ell30_t16_128_direction/d_tensor_interface_jump_ell30_n128_t16.png`

### 结果：128x128

- 总界面段数：`1361`
- trace 速度通过数值验证的段数：`614`
- unresolved：`747`
- solved fraction：`0.451`
- trace 目标法向范数：
  - median `5.646e-3`
  - p95 `5.667e-3`
- trace residual after solve：
  - median `4.66e-10`
  - p95 `4.32e-6`
- 固定当前参考界面法向后，张量残差仍很大：
  - `tensor_residual_relative` median `4.37e5`
  - p95 `9.08e7`
  - 无迹部分同量级
- free-covector rank-one mismatch：
  - median `0.018`
  - p95 `0.431`
- free-covector 与当前参考界面空间法向对齐：
  - median `0.486`
  - p95 `0.995`

### 解释

1. trace 速度律本身是稳定的，目标法向范数仍是 \(O(10^{-3})\)。
2. 但把 trace 解代入完整张量 jump 后，固定参考界面方向的无迹残差极大。
3. free-covector rank-one mismatch 的中位数较小，说明目标张量并非普遍不可能由某个 \(F_a\) 表示。
4. 目标自由法向和当前 `raw_y=±1` 参考界面法向的对齐中位数只有 `0.486`，说明主要失败不是“速度大小没调好”，而是“界面形状/法向方向不对”。
5. p95 rank-one mismatch 仍达到 `0.431`，说明在尾部区域还需要更完整的张量 jump、切向项、两侧几何匹配或更高阶层内结构。

### 当前结论

仅推进 signed-distance 界面位置还不够。下一版 full-interface solver 至少需要：

- 让界面法向方向成为未知量，而不是完全继承平直参考 `raw_y` 的空间梯度；
- 在界面上同时求 trace 条件和无迹 tensor matching；
- 允许两侧 bulk 几何或外曲率类数据参与匹配；
- unresolved 段不能用 damping/clipping 处理。

### 下一步

最小下一步不是直接全演化，而是写一个 local tensor-matching interface solver：

1. 输入当前界面点的 \(\tilde g_{ab}\)、algebraic tensor integral；
2. 由 rank-one 近似反推出目标 covector \(F_a\)；
3. 把目标空间法向转成 signed-distance 界面速度和界面形状校正；
4. 在 signed-distance moving-interface 原型中用这个 tensor-projected normal 替代 trace-only normal；
5. 再检查界面运动是否从“固定参考法向巨大张量残差”下降到可控水平。
