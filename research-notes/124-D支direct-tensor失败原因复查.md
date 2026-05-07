## D支 direct tensor 失败原因复查

日期：2026-05-02

### 背景

上一轮把 direct tensor 条件分成：

- `accepted`
- `tensor_residual`
- `rank_tail+negative+tensor_residual`
- `rank_tail+tensor_residual`
- `no_real_covector`

但这还不能说明失败是物理上的，还是数值/求解方式导致的。

本轮做三类复查：

1. 局部特征诊断：密度、`raw_y` 梯度、度规条件数等；
2. 稳定性检查：分辨率、\(\ell\)、层宽；
3. 局部最小二乘检查：直接最小化完整张量残差，判断 `tensor_residual` 是否只是初始投影算法不好。

### 新脚本

#### 特征诊断

`kg_examples/diagnose_d_direct_tensor_failure_features.py`

为每条界面段补充：

- `rho_relative = rho/rho_max`
- `raw_y_grad_norm = |\nabla raw_y|`
- `metric_condition_abs_eig`
- `metric_det_abs`
- `algebraic_tensor_norm`
- `target_norm_relative_error`
- 原始 direct tensor residual 与 rank-one 指标

#### 局部最小二乘检查

`kg_examples/diagnose_d_direct_tensor_least_squares.py`

直接把完整张量残差

\[
\Delta(\partial_q f_R)(-F_aF_b+\tilde g_{ab}F^2)+A_{ab}
\]

看成 \(F_t,F_x,F_z\) 的函数，用 Gauss-Newton 在三个未知量上最小化。

这不是添加新物理，而是更直接地求解用户要求的完整张量条件。它用于检查此前的 rank-one 特征向量投影是否过于粗糙。

### 局部特征结果

`ell=30,t=16,n=128,10%`。

按密度分箱：

- `rho/rho_max = 1e-3..1e-2`：accepted `0.217`
- `1e-2..5e-2`：accepted `0.225`
- `5e-2..1e-1`：accepted `0.310`
- `0.1..0.3`：accepted `0.422`
- `>0.3`：accepted `0.601`

按 `raw_y` 梯度分箱：

- `|grad raw_y|<1e3`：accepted `0.827`
- `1e3..1e5`：accepted `0.321`
- `1e5..1e8`：accepted `0.070`
- `1e8..1e12`：accepted `0.057`
- `>1e12`：accepted `0.182`

按度规条件数分箱：

- `cond<20`：accepted `0.485`
- `20..50`：accepted `0.523`
- `50..100`：accepted `0.353`
- `100..1e3`：accepted `0.030`
- `1e3..1e5`：accepted `0.030`
- `>1e5`：accepted `0.000`

解释：

失败明显集中在低密度、过渡层极陡、局部度规病态区域。

### 稳定性检查

固定 `ell=30,t=16,half-width=1`，改变分辨率：

- `n=96`：accepted `0.254`
- `n=128`：accepted `0.335`
- `n=160`：accepted `0.398`

提高分辨率会改善 accepted 比例，说明有一部分失败是解析度问题。

但主要失败类型没有消失：

- `n=96`：`rank_tail+negative+tensor_residual` 为 `0.400`
- `n=128`：`0.331`
- `n=160`：`0.287`

改变 \(\ell\)：

- `ell=10,n=128`：accepted `0.335`
- `ell=30,n=128`：accepted `0.335`

失败比例几乎不变，说明当前失败不是简单由 \(\ell\) 取值导致。

改变层宽：

- `half-width=0.5`：accepted `0.123`
- `half-width=1`：accepted `0.335`
- `half-width=2`：accepted `0.335`

解释：层宽太窄会漏掉完整过渡层；一旦覆盖到 `±1` 附近，继续加宽不改变结果。

### 局部最小二乘结果

`10%` 容差：

- 原始 direct tensor accepted：`456/1361 = 0.335`
- 最小二乘后 accepted：`867/1361 = 0.637`

分模式：

- 原本 `accepted`：`456/456` 仍 accepted；
- 原本 `tensor_residual`：`327/378` 被修复；
- 原本 `rank_tail+tensor_residual`：`6/28` 被修复；
- 原本 `rank_tail+negative+tensor_residual`：`78/453` 被修复；
- `no_real_covector`：`0/44` 被修复。

`1%` 容差：

- 原始 accepted：`62/1361 = 0.0456`
- 最小二乘后 accepted：`333/1361 = 0.2447`

分模式：

- 原本 `accepted`：`62/62` 仍 accepted；
- 原本 `tensor_residual`：`223/277` 被修复；
- 原本 `rank_tail+tensor_residual`：`4/191` 被修复；
- 原本 `rank_tail+negative+tensor_residual`：`44/787` 被修复；
- `no_real_covector`：`0/44` 被修复。

### 结论

失败原因至少分三层。

第一层：算法投影问题。

此前把 \(K_{ab}=\tilde g_{ab}F^2-H_{ab}\) 做普通 rank-one 特征向量投影，再代回完整张量方程。这对精确 rank-one 情况是对的，但对近似数据不一定是最优的完整张量残差解。

局部最小二乘说明：`tensor_residual` 型失败大部分可修复。因此这类不应再解释为物理失败，而应升级求解器：直接最小化完整 tensor residual。

第二层：数值解析度/局部病态问题。

失败强烈集中在：

- 低密度区；
- `raw_y` 梯度极大区；
- \(\tilde g\) 条件数大的区域。

提高分辨率能把 accepted 从 `0.254` 提到 `0.398`，说明至少一部分失败是界面未解析或局部插值病态。

第三层：可能的物理闭合问题。

`rank_tail+negative+tensor_residual` 中大多数即使局部最小二乘也不能通过：

- `10%` 下只有 `78/453` 被修复；
- `1%` 下只有 `44/787` 被修复。

这类更像真正困难：单一薄界面 rank-one 主部不够，需要补：

- 层内 Ricci profile；
- 切向/外曲率项；
- 两侧 bulk 几何共同调整；
- 或者它们确实是 D 支局部不可闭合的信号。

### 当前修正后的下一步

后续 direct tensor interface solver 应从“rank-one 特征向量投影”升级为：

\[
\min_{F_a}\left\|
\Delta(\partial_q f_R)(-F_aF_b+\tilde g_{ab}F^2)+A_{ab}
\right\|.
\]

然后：

1. 用最小二乘通过的 `63.7%` 段重构界面；
2. 对剩余失败段做分辨率/层内 profile/两侧 bulk matching 检查；
3. 不能再把 `tensor_residual` 型直接当作物理失败。

### 输出

- 特征图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_failure_features_ell30_t16_128/d_direct_tensor_failure_features_ell30_n128_t16.png`
- 特征 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_failure_features_ell30_t16_128/summary.json`
- 最小二乘图：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_least_squares_ell30_t16_128/d_direct_tensor_least_squares_ell30_n128_t16.png`
- 最小二乘 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_least_squares_ell30_t16_128/summary.json`
- 严格最小二乘 summary：`/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_direct_tensor_least_squares_ell30_t16_128_strict001/summary.json`
