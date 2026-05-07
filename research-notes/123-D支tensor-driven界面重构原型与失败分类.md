## D支 tensor-driven 界面重构原型与失败分类

日期：2026-05-02

### 目的

在用户指出 trace 估计不是充分条件后，当前路线已经改为直接使用完整 leading tensor condition。

本轮继续推进两件事：

1. 用 direct tensor 条件给出的目标空间协向量 \((F_x,F_z)\) 生成候选界面曲线；
2. 对 direct tensor 未通过段做失败类型分类，判断下一步是单纯重构界面法向，还是必须引入更完整的两侧 bulk/tensor matching。

### 物理口径

重构使用的条件仍是：

\[
H_{ab}=-F_aF_b+\tilde g_{ab}F^2.
\]

只有通过完整 direct tensor 验收的段才用于构造目标切线场。未通过段不会被 trace speed、damping、clipping 或插值强行补成可解段。

目标空间法向为：

\[
n_i \propto (F_x,F_z).
\]

目标界面切线取：

\[
\tau_i \propto (-F_z,F_x).
\]

这是把 tensor condition 给出的法向整理成连续界面表示，不是新增物理规则。

### 新脚本

#### 1. Tensor-driven 重构器

`kg_examples/prototype_d_tensor_reconstructed_interface.py`

功能：

- 重新计算 direct tensor rows；
- 分别处理 `raw_y=+1` 与 `raw_y=-1` 两个分支，避免混合两层界面；
- 只从 `direct_tensor_solved=1` 的段构造目标切线场；
- 沿目标切线场积分出候选界面曲线；
- 输出原始 rows、重构 curves 和图。

#### 2. Rejection mode 分类器

`kg_examples/diagnose_d_direct_tensor_rejection_modes.py`

分类：

- `accepted`：完整 direct tensor rank-one 与 residual 验收通过；
- `no_real_covector`：无法构造正的实 rank-one covector；
- `rank_tail`：\(g_{ab}F^2-H_{ab}\) 不够 rank-one；
- `negative`：rank-one 目标有过大负特征值权重；
- `tensor_residual`：rank-one 近似存在，但完整 tensor residual 超限。

组合失败用 `+` 连接。

### 运行

`10%` 验收：

```bash
python3 kg_examples/prototype_d_tensor_reconstructed_interface.py \
  --output visualizations/d_tensor_reconstructed_interface_ell30_t16_128 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.1 --tensor-rel-tol 0.1
```

```bash
python3 kg_examples/diagnose_d_direct_tensor_rejection_modes.py \
  --rows-jsonl visualizations/d_tensor_reconstructed_interface_ell30_t16_128/direct_tensor_interface_rows.jsonl \
  --output visualizations/d_direct_tensor_rejection_modes_ell30_t16_128 \
  --rank1-tol 0.1 --tensor-rel-tol 0.1 \
  --time 16 --resolution 128
```

`1%` 验收：

```bash
python3 kg_examples/prototype_d_tensor_reconstructed_interface.py \
  --output visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.01 --tensor-rel-tol 0.01
```

```bash
python3 kg_examples/diagnose_d_direct_tensor_rejection_modes.py \
  --rows-jsonl visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/direct_tensor_interface_rows.jsonl \
  --output visualizations/d_direct_tensor_rejection_modes_ell30_t16_128_strict001 \
  --rank1-tol 0.01 --tensor-rel-tol 0.01 \
  --time 16 --resolution 128
```

### 输出

`10%`：

- 重构图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128/d_tensor_reconstructed_interface_ell30_n128_t16.png`
- 重构 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128/summary.json`
- 重构曲线：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128/tensor_reconstructed_curves.json`
- 失败分类图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128/d_direct_tensor_rejection_modes_n128_t16.png`
- 失败分类 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128/summary.json`

`1%`：

- 重构图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/d_tensor_reconstructed_interface_ell30_n128_t16.png`
- 重构 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/summary.json`
- 重构曲线：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/tensor_reconstructed_curves.json`
- 失败分类图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128_strict001/d_direct_tensor_rejection_modes_n128_t16.png`
- 失败分类 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_rejection_modes_ell30_t16_128_strict001/summary.json`

### 重构结果

`ell=30,t=16,n=128`。

`10%` 验收：

- accepted：`456/1361`，fraction `0.335`
- 重构曲线：`29` 条
- 分支：`raw_y=-1` 有 `10` 条，`raw_y=+1` 有 `19` 条
- 总长度：`165.86`
- accepted 点覆盖率：`1.0`
- accepted 点到重构曲线最近距离 median `0.166`，p95 `0.376`

`1%` 验收：

- accepted：`62/1361`，fraction `0.0456`
- 重构曲线：`16` 条
- 分支：`raw_y=-1` 有 `4` 条，`raw_y=+1` 有 `12` 条
- 总长度：`44.15`
- accepted 点覆盖率：`1.0`
- accepted 点到重构曲线最近距离 median `0.0110`，p95 `0.331`

解释：accepted 区域本身可以被整理成连续目标曲线；限制不在重构算法，而在 direct tensor accepted 区域太少。

### 失败分类结果

`10%` 验收的 mode counts：

- `accepted`: `456`，`33.5%`
- `tensor_residual`: `378`，`27.8%`
- `rank_tail+negative+tensor_residual`: `451`，`33.1%`
- `no_real_covector`: `44`，`3.23%`
- 其他组合合计约 `2.35%`

`1%` 验收的 mode counts：

- `accepted`: `62`，`4.56%`
- `tensor_residual`: `277`，`20.4%`
- `rank_tail+tensor_residual`: `191`，`14.0%`
- `rank_tail+negative+tensor_residual`: `787`，`57.8%`
- `no_real_covector`: `44`，`3.23%`

### 解释

1. 当前 direct tensor accepted 区域主要位于较规则的外侧分支；中部/内侧复杂区大量 rejected。
2. 有一批段属于 `tensor_residual` 型：rank-one 条件近似可行，但完整 tensor residual 超限。这类段可能需要更完整的层内 profile、切向项、外曲率或两侧 bulk 几何共同调整。
3. 另一大批段属于 `rank_tail+negative+tensor_residual`：它们不是单一 rank-one 界面主部能直接表示的。若这些段不是数值误差，就意味着需要更高阶 matching 或当前 D 支闭合在这些局部失败。
4. `no_real_covector` 只占约 `3.2%`，不是主要问题。

### 下一步

下一步不应把 unresolved 段用插值补上，而应进入两条并行检查：

1. 对 `tensor_residual` 型段，加入更完整的 tensor matching 数据：两侧 bulk 几何、层内 Ricci profile、切向/外曲率项，检查是否能把完整 residual 降下来。
2. 对 `rank_tail+negative+tensor_residual` 型段，检查它们是否随分辨率、层宽采样、\(\ell\) 改变而消失；如果稳定存在，则可能是 D 支局部不可闭合的物理信号。
