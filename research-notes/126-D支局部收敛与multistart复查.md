# D支局部收敛与 multistart 复查

日期：2026-05-02

## 目的

在 least-squares direct tensor solver 已经显著优于 rank-one 投影后，需要判断剩余失败段的性质：

- 是否只是全域网格分辨率不够；
- 是否只是局部 Gauss-Newton 初值选择不佳；
- 还是已经进入需要补充层内 profile / 两侧 bulk matching 的物理闭合问题。

本轮只检查同一个 direct tensor 条件：

\[
\Delta(\partial_q f_R)\left(-F_aF_b+\tilde g_{ab}F^2\right)+A_{ab}=0 .
\]

未加入 trace speed、阻尼、裁剪、插值补齐或任何新的闭合规则。

## 新增脚本

新增局部收敛诊断：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_ls_local_resolution_convergence.py`

功能：

- 在 `n=160` 的 least-squares 未解段中自动选出最密集的物理窗口；
- 在同一物理窗口上比较 `n=128,160,192`；
- 输出全局和局部 accepted fraction、未解 residual p95、未解段局部特征。

同时增强 direct tensor solver：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_d_direct_tensor_least_squares.py`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_tensor_reconstructed_interface.py`

新增：

- `multistart_covector_seeds`
- `gauss_newton_refine_multistart`
- `--multistart`
- `--max-seeds`

multistart 只是对同一个张量残差换多个实协向量初值，不是新物理规则。

## 局部分辨率结果

运行：

```bash
python3 kg_examples/diagnose_d_ls_local_resolution_convergence.py \
  --output visualizations/d_ls_local_resolution_convergence_ell30_t16_128_160_192 \
  --ell 30 --time 16 \
  --resolutions 128,160,192 \
  --base-resolution 160 \
  --half-width-y 1 --samples 81 \
  --tensor-rel-tol 0.1 --max-iter 30 \
  --focus-rho-min-frac 1e-4 \
  --window-width-x 4 --window-width-z 4 \
  --top-windows 4
```

自动选出的四个窗口：

- `W1`: \(x\in[-4,0], z\in[4,8]\)
- `W2`: \(x\in[0,4], z\in[4,8]\)
- `W3`: \(x\in[-4,0], z\in[8,12]\)
- `W4`: \(x\in[0,4], z\in[8,12]\)

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_ls_local_resolution_convergence_ell30_t16_128_160_192/d_ls_local_resolution_convergence_ell30_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_ls_local_resolution_convergence_ell30_t16_128_160_192/summary.json`

全局 accepted fraction：

- `n=128`: `0.6370`
- `n=160`: `0.6573`
- `n=192`: `0.6458`

局部 accepted fraction：

- `W1`: `0.569 -> 0.601 -> 0.607`
- `W2`: `0.483 -> 0.596 -> 0.519`
- `W3`: `0.423 -> 0.453 -> 0.597`
- `W4`: `0.440 -> 0.486 -> 0.465`

未解段 residual p95：

- 全局：`1.008 -> 1.040 -> 1.010`
- `W1`: `1.001 -> 1.075 -> 1.011`
- `W2`: `1.008 -> 1.014 -> 1.006`
- `W3`: `1.031 -> 1.052 -> 1.014`
- `W4`: `1.083 -> 1.010 -> 1.038`

解释：

- 分辨率对某些局部窗口有改善，尤其 `W3`；
- 但改善不单调，且未解 residual p95 始终贴近 `~1`；
- 因此剩余失败不能主要解释为普通网格欠分辨。

## multistart 结果

运行：

```bash
python3 kg_examples/prototype_d_tensor_reconstructed_interface.py \
  --output visualizations/d_tensor_reconstructed_interface_ls_multistart_ell30_t16_160 \
  --ell 30 --resolution 160 --time 16 --levels 1 \
  --half-width-y 1 --samples 81 \
  --tensor-rel-tol 0.1 \
  --solver least_squares --max-iter 30 \
  --multistart --max-seeds 12
```

输出：

- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_multistart_ell30_t16_160/d_tensor_reconstructed_interface_ell30_n160_t16.png`
- `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_tensor_reconstructed_interface_ls_multistart_ell30_t16_160/summary.json`

对比普通 least-squares：

- 普通 least-squares：accepted `1214/1847 = 0.6573`
- multistart least-squares：accepted `1242/1847 = 0.6724`
- 只额外修复 `28` 段；
- 没有把大块 residual-fail 区域打通。

被修复的 `28` 段来源：

- `no_candidate`: `18`
- `rank_tail+negative+residual`: `6`
- `residual`: `4`

仍失败：

- `605` 段仍失败；
- 其中旧初始类型最多的是 `rank_tail+negative+residual`，约 `451` 段；
- 旧 residual 中位数约 `1.000005`，multistart 后仍约 `0.999989`。

## 判断

当前证据支持：

\[
\text{剩余主要失败不是单纯分辨率问题，也不是单初值 Gauss-Newton 局部极小问题。}
\]

更合理的解释是：

- 当前 leading thin-layer algebraic tensor integral 还缺层内 profile 信息；
- 或缺切向/外曲率项；
- 或缺两侧 EH-like / saturated bulk geometry matching；
- 或在低密度、极陡层、度规病态区域，D 支确实出现局部不可闭合。

## 下一步

不要继续盲目全域加密。下一步应在 `W1-W4` 四个窗口内做：

- 层内 `q` 方向 profile 采样，而不是只看积分后的 \(A_{ab}\)；
- 把 residual tensor 分解到 trace、normal-normal、normal-tangent、tangent-tangent 分量；
- 检查失败段是否集中在某个张量分量；
- 加入两侧 bulk matching 的候选项，判断 residual≈1 是否能被系统性消掉；
- 若仍不能消掉，再把这些段标记为 D 支局部闭合失败候选。
