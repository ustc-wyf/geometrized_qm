## D支 direct tensor interface 条件

日期：2026-05-02

### 背景

用户指出：trace jump law 只是必要条件，不是充分条件；既然目标是 tensor matching，就不应继续把 trace 估计当作闭合条件。

本轮已修正路线：直接从完整张量 jump 条件反推界面协向量 \(F_a\)，trace 只作为派生检查，不作为求解条件。

### 方程

张量薄层条件写成：

\[
\Delta(\partial_q f_R)
\left(-F_aF_b+\tilde g_{ab}F^2\right)+A_{ab}=0,
\]

其中

\[
A_{ab}:=
\int_{\rm layer}
\left(
f_R\tilde R_{ab}
-\frac12 f\tilde g_{ab}
-\frac{\tilde T_{ab}}{M_P^2}
\right)dq.
\]

定义

\[
H_{ab}:=-\frac{A_{ab}}{\Delta(\partial_q f_R)}.
\]

直接张量条件要求：

\[
H_{ab}=-F_aF_b+\tilde g_{ab}F^2.
\]

在 \(2+1d\) 中：

\[
F^2=\frac12\operatorname{Tr}_{\tilde g}H,
\]

于是

\[
F_aF_b=\tilde g_{ab}F^2-H_{ab}.
\]

所以直接判据是：右侧必须近似为一个 rank-one 的实协向量外积。

### 新脚本

`kg_examples/diagnose_d_direct_tensor_interface.py`

与上一轮 `diagnose_d_tensor_interface_jump.py` 的区别：

- 不先解 trace speed；
- 直接从 \(H_{ab}\) 构造 \(F_aF_b\)；
- 用 rank-one 分解得到目标 \(F_a\)；
- 用完整 tensor residual 检查是否接受；
- trace 只作为派生一致性检查。

`direct_tensor_solved` 的含义：

- 存在实的 rank-one 候选 covector；
- `rank1_tail_relative` 小于给定数值阈值；
- `rank1_negative_relative` 小于给定数值阈值；
- `direct_tensor_residual_relative` 小于给定数值阈值。

这些阈值只是数值验收精度，不是物理阻尼、裁剪或额外规则。

### 运行

`10%` 数值验收：

```bash
python3 kg_examples/diagnose_d_direct_tensor_interface.py \
  --output visualizations/d_direct_tensor_interface_ell30_t16_128_v2 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.1 --tensor-rel-tol 0.1
```

`1%` 数值验收：

```bash
python3 kg_examples/diagnose_d_direct_tensor_interface.py \
  --output visualizations/d_direct_tensor_interface_ell30_t16_128_strict001 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.01 --tensor-rel-tol 0.01
```

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_v2/d_direct_tensor_interface_ell30_n128_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_v2/summary.json`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_strict001/d_direct_tensor_interface_ell30_n128_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_interface_ell30_t16_128_strict001/summary.json`

### 结果

`ell=30,t=16,n=128`，总界面段数 `1361`。

候选层面：

- 可构造实 covector 候选：`1317`
- candidate fraction：`0.968`
- 全体候选 rank-one mismatch：
  - median `0.0354`
  - p95 `~1.0`
- 全体候选 direct tensor residual relative：
  - median `0.308`
  - p95 `4.318`

`10%` 验收：

- accepted：`456/1361`
- fraction：`0.335`
- accepted residual median：`0.0254`
- accepted residual p95：`0.0887`
- accepted speed median：`0.996`
- accepted speed p95：`1.012`
- accepted target/current spatial-normal alignment median：`0.419`

`1%` 验收：

- accepted：`62/1361`
- fraction：`0.0456`
- accepted residual median：`0.00775`
- accepted residual p95：`0.00972`
- accepted speed median：`0.999`
- accepted speed p95：`1.002`
- accepted target/current spatial-normal alignment median：`0.411`

### 解释

1. 用户的批评是正确的：trace jump law 不是充分条件。
2. 直接张量条件下，当前参考界面 `raw_y=±1` 只有一部分段能满足完整 tensor matching。
3. 通过直接张量条件得到的速度非常接近单位量级，说明速度幅值本身不是主要难点。
4. 即使 accepted 段，目标 \(F_a\) 的空间方向和当前参考界面法向的对齐中位数也只有 `~0.4`，说明真正需要改变的是界面形状/法向，而不是只调 \(F_t\)。
5. p95 尾部 rank-one mismatch 接近 `1`，说明部分区域不能由单一 rank-one 界面主部解释，可能需要：
   - 两侧 bulk 几何共同调整；
   - 更完整切向/外曲率匹配；
   - 层内结构而非单面近似；
   - 或在这些段上当前候选 D 支不可闭合。

### 结论

从现在开始，不能再把 trace-only moving-interface 当作主闭合。下一步应直接基于 tensor condition：

1. 构造由 \(H_{ab}\) 反推的目标 \(F_a\) 场；
2. 用目标 \(F_x,F_z\) 重构/校正界面形状，而不是沿当前 `raw_y` 法向推进；
3. 只在 direct tensor residual 合格的段上推进；
4. 对不合格段进入 body-fitted/tensor matching 或标记为局部闭合失败。
