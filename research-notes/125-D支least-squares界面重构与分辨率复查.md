# D支 least-squares 界面重构与分辨率复查

日期：2026-05-02

## 问题

上一轮 direct tensor 条件表明，原来的 rank-one 特征向量投影会把大量 `tensor_residual` 段判为失败。用户进一步确认当前计划，并追问“分辨率是否也需要提高”。本轮目标是：

- 把 `prototype_d_tensor_reconstructed_interface.py` 的界面重构器从 rank-one 投影升级为局部最小二乘 \(F_a\) 求解；
- 在同一物理设置 `ell=30,t=16` 下比较 `n=128` 与 `n=160`；
- 判断“提高分辨率”在当前问题里到底是主瓶颈还是辅助需求。

## 方法

物理条件保持不变，仍然是 direct tensor matching 的 leading 条件：

\[
\Delta(\partial_q f_R)\left(-F_aF_b+\tilde g_{ab}F^2\right)+A_{ab}=0 .
\]

其中 \(F_a\) 是时空界面协法向，\(A_{ab}\) 是沿过渡层短线积分得到的代数张量源。变化只在数值求解方式：

- 旧方法：先由 rank-one 投影给出 \(F_a\)，再检查完整张量残差；
- 新方法：以 rank-one 投影作为初值，用 Gauss-Newton 直接最小化完整张量残差；
- `accepted` 仍只由完整张量残差阈值决定，不使用 trace speed、damping、clipping 或人为补齐。

## 数值设置

共同设置：

- `branch = D`
- `ell = 30`
- `time = 16`
- `half_width_y = 1`
- `samples = 81`
- `tensor_rel_tol = 0.1`
- `solver = least_squares`
- `max_iter = 30`

运行命令：

```bash
python3 kg_examples/prototype_d_tensor_reconstructed_interface.py \
  --output visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_128 \
  --ell 30 --resolution 128 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.1 --tensor-rel-tol 0.1 \
  --solver least_squares --max-iter 30

python3 kg_examples/prototype_d_tensor_reconstructed_interface.py \
  --output visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_160 \
  --ell 30 --resolution 160 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --rank1-tol 0.1 --tensor-rel-tol 0.1 \
  --solver least_squares --max-iter 30
```

## 结果

`n=128`：

- 界面段 `row_count = 1361`
- candidate `1317`，candidate fraction `0.9677`
- accepted `867`，accepted fraction `0.6370`
- accepted residual median `0.01489`，p95 `0.07362`
- target speed median `0.99858`，p95 `1.02843`
- target/current normal alignment median `0.4449`
- 重构曲线 `30` 条，`-1:15,+1:15`
- accepted coverage `1.0`
- accepted nearest curve distance median `0.1435`，p95 `0.3385`

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_128/d_tensor_reconstructed_interface_ell30_n128_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_128/summary.json`

`n=160`：

- 界面段 `row_count = 1847`
- candidate `1769`，candidate fraction `0.9578`
- accepted `1214`，accepted fraction `0.6573`
- accepted residual median `0.01410`，p95 `0.06934`
- target speed median `0.99893`，p95 `1.02342`
- target/current normal alignment median `0.4600`
- 重构曲线 `40` 条，`-1:12,+1:28`
- accepted coverage `1.0`
- accepted nearest curve distance median `0.1210`，p95 `0.2628`

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_160/d_tensor_reconstructed_interface_ell30_n160_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_160/summary.json`

## 分辨率判断

分辨率确实需要提高，但其角色要严格限定：

- 提高分辨率能改善界面定位、插值和重构曲线的几何连续性；
- 在旧 rank-one 投影下，`n=96,128,160` 的 accepted fraction 为 `0.254 -> 0.335 -> 0.398`，说明解析度确实有影响；
- 但更大的改善来自求解器升级：在 `n=128` 下，least-squares 把 accepted fraction 从 `0.335` 提高到 `0.637`；
- 从 `n=128` 到 `n=160`，least-squares 只进一步从 `0.637` 提高到 `0.657`，属于必要但较小的数值收敛改善。

因此当前判断是：

\[
\text{主瓶颈是 tensor matching 的局部非线性求解；分辨率是后续收敛与界面几何解析的必要条件。}
\]

不能把加密网格当作替代张量条件、层内 profile 或两侧 bulk matching 的物理闭合。

## 下一步

- 后续 direct tensor 主求解器默认采用局部 least-squares \(F_a\)，rank-one 投影只保留为初值和诊断；
- 用 least-squares accepted 段作为 tensor-driven signed-distance/body-fitted interface 的输入；
- 对仍失败段，重点检查 `rank_tail+negative+tensor_residual` 与 `no_real_covector` 是否来自低密度、极陡层、局部 \(\tilde g\) 病态、缺失层内 profile 或真正局部不可闭合；
- 做 `n=192/240` 的局部窗口复查，而不是盲目全域加密。
