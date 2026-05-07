# sparse/matrix-free 求解器首轮

日期：2026-05-07

## 目标

把上一轮 dense 的：

- principal-constraint projection；
- full-linear metric update；

改成不显式存储大矩阵的 sparse/matrix-free 版本，为 `n=384` 高分辨率做准备。

## full-linear metric update sparse 版本

新增脚本：

- `kg_examples/solve_gbcd_full_linear_metric_update_sparse.py`

做法：

- 不存 dense \(A\)；
- 对每个 \(\delta\tilde g_+\) 变量只存它影响的少数 stencil 行；
- 支持 `cg` normal equation 和 `lsqr` augmented least-squares。

### 算子一致性

在 `n=96,tau=0` 上，比较 dense \(A\) 与 sparse-column \(A\)：

- forward relative difference: \(1.75\times10^{-16}\)
- transpose relative difference: \(1.54\times10^{-16}\)

结论：稀疏装配本身正确。

### 迭代求解差距

和 dense `ridge=1e-6` 对照：

| solver | iterations | linear wmean | exact wmean | exact p95 |
|---|---:|---:|---:|---:|
| dense lstsq | direct | \(4.94\times10^{-5}\) | 0.00661 | 0.04951 |
| sparse CG | 2000 | \(4.71\times10^{-4}\) | 0.01156 | 0.17572 |
| sparse LSQR | 1000 | 0.00107 | 0.01413 | 0.16817 |

结论：

- sparse 算子正确；
- 当前迭代器还不能完全复现 dense direct solve；
- 但 exact residual 已经保持在百分之一级，说明方向仍可用。

## sparse principal-constraint projection

新增脚本：

- `kg_examples/fit_gbcd_principal_constraint_projection_sparse.py`

做法：

- algebraic rows；
- full conservation rows；
- principal left-null constraint rows；
- norm/time regularizer rows；

统一作为 sparse LSQR 最小二乘系统求解。

重要 caveat：

- 当前 sparse 版用 penalty rows 近似 hard constraints；
- 还不是 dense 版那种 nullspace hard solve。

### n=96 core10 对照

dense hard principal + dense full-linear：

- exact residual weighted mean: `0.02548`
- exact residual p95: `0.09238`

sparse soft-principal + sparse full-linear：

- exact residual weighted mean: `0.07885`
- exact residual p95: `0.53520`

结论：

- sparse pipeline 已跑通；
- 但精度明显差于 dense hard pipeline；
- 主要原因不是物理方程失败，而是 sparse principal projection 还只是 soft penalty，且 metric update 迭代解还未达到 dense direct solve 质量。

## 下一步

不能直接把当前 sparse soft pipeline 当作高分辨率物理结果。

下一步优先级：

1. 把 sparse principal projection 从 penalty 形式升级为约束保持形式：
   - nullspace/Schur complement；
   - 或 augmented Lagrangian 迭代逐步提高约束；
   - 并显式报告 conservation/principal constraint residual。

2. 改进 full-linear metric update 的迭代器：
   - 更好的预条件；
   - block/local preconditioner；
   - 或分区 direct + 全局迭代。

3. 只有当 `n=96` sparse pipeline 接近 dense pipeline 后，才进入 `n=384 core10`。

当前阶段结论：

- 稀疏算子装配正确；
- 高分辨率路线可行；
- 但高分辨率物理结论还不能产出，必须先完成约束保持和预条件改进。
