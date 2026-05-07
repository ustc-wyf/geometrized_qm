## D支 direct tensor 目标法向重构预览

日期：2026-05-02

### 目的

用户指出 trace 估计不是充分条件，应直接使用张量条件。

本轮新增的重构预览不是新的物理规则，而是把 direct tensor 条件反推出的目标协向量

\[
F_a=(F_t,F_x,F_z)
\]

直接画出来，用来判断后续 moving-interface solver 应怎样更新界面法向/形状。

### 脚本

`kg_examples/plot_d_direct_tensor_normal_reconstruction.py`

它复用 `diagnose_d_direct_tensor_interface.py` 的 direct tensor 条件：

\[
H_{ab}=-F_aF_b+\tilde g_{ab}F^2,
\qquad
H_{ab}=-A_{ab}/\Delta(\partial_q f_R).
\]

验收条件仍是完整张量条件：

1. 可构造实 rank-one 协向量候选；
2. `rank1_tail_relative` 通过数值容差；
3. `rank1_negative_relative` 通过数值容差；
4. `direct_tensor_residual_relative` 通过数值容差。

没有先求 trace speed，也没有用 trace 作为 solved 条件。

### 图中变量

- `linear rho + current raw_y=+/-1`：线性坐标下的 \(\rho\)，叠加当前参考界面 `raw_y=ell^2 R_tilde=±1`。
- `direct tensor accepted`：黑线是通过 direct tensor 条件的界面段，灰色虚线是未通过段。
- `target linelets`：蓝色短线是与 direct tensor 目标空间协向量 \((F_x,F_z)\) 垂直的局部目标切线；红色短线是目标空间法向。
- `direct tensor residual / algebraic`：完整张量残差相对于代数层积分张量的大小。
- `target-current normal alignment`：direct tensor 目标空间法向与当前 `raw_y=±1` 空间法向的夹角对齐度，`1` 表示同向或反向，`0` 表示正交。
- `target coordinate speed`：\(-F_t/|F_{\rm spatial}|\)，只作为 direct tensor 目标协向量给出的坐标法向速度诊断。

### 运行

`10%` 数值验收：

```bash
python3 kg_examples/plot_d_direct_tensor_normal_reconstruction.py \
  --output visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.1 --tensor-rel-tol 0.1
```

`1%` 数值验收：

```bash
python3 kg_examples/plot_d_direct_tensor_normal_reconstruction.py \
  --output visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.01 --tensor-rel-tol 0.01
```

### 输出

`10%`：

- 图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/d_direct_tensor_normal_reconstruction_ell30_n128_t16.png`
- summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/summary.json`
- 每段数据：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/direct_tensor_interface_rows.jsonl`

`1%`：

- 图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/d_direct_tensor_normal_reconstruction_ell30_n128_t16.png`
- summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/summary.json`
- 每段数据：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/direct_tensor_interface_rows.jsonl`

### 结果

`ell=30,t=16,n=128`，总界面段 `1361`，实 covector candidate `1317`。

`10%` 验收：

- accepted：`456/1361`，fraction `0.335`
- accepted tensor residual median `0.0254`，p95 `0.0887`
- accepted target speed median `0.9957`，p95 `1.0124`
- accepted target/current normal alignment median `0.4186`

`1%` 验收：

- accepted：`62/1361`，fraction `0.0456`
- accepted tensor residual median `0.00775`，p95 `0.00972`
- accepted target speed median `0.9990`，p95 `1.0023`
- accepted target/current normal alignment median `0.4114`

### 结论

这次结果把问题锁定得更清楚：

1. trace-only 的速度律确实不能作为完整闭合；
2. direct tensor 条件给出的目标速度幅值很温和，接近 `1`；
3. 真正的障碍是当前 `raw_y=±1` 参考界面的空间法向通常不等于 tensor matching 要求的目标法向；
4. 下一步的全动力学界面模拟器应把界面形状/法向作为未知量，由 direct tensor 条件重构 signed-distance/body-fitted 界面；
5. 未通过 direct tensor 的段不能用 trace speed、damping 或 clipping 补掉，只能标为 tensor-unresolved，或进入更完整的两侧 bulk/tensor matching。
