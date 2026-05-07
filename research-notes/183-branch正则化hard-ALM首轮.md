# branch 正则化 hard-ALM 首轮

日期：2026-05-08

## 目的

上一轮 `182-branch-patch数值检查.md` 发现：未显式加入 branch 正则性时，hard-ALM 得到的 \(C_{\mu\nu}\) 在 \(q\approx0\) 区域并不会自动满足 \(C\sim\chi(q)\)。

本轮把未知量改写为

\[
C_{\mu\nu}
=
\chi(\mathcal q)\hat C_{\mu\nu},
\]

并直接用同一套 trace0 + full conservation + ALM 求解 \(\hat C\)，测试带 branch 规则的 proposal v1.1 是否数值可行。

## 新增脚本

```text
kg_examples/fit_gbcd_trace0_branch_sparse_conservation.py
```

脚本逻辑：

1. 用三时间层构造 \(\mathcal q=(\tilde\nabla_\mu r^\mu+r_\mu r^\mu)/m^2\)；
2. 定义
   \[
   q_{\rm rel}=|\mathcal q|/q_{\rm scale},
   \qquad
   \chi=\frac{q_{\rm rel}^2}{q_{\rm rel}^2+q_0^2};
   \]
3. 以 \(\hat C\) 的系数为未知量，实际进入代数方程和守恒方程的是 \(C=\chi\hat C\)；
4. trace0 仍通过 nullspace 硬消元；
5. full conservation 仍通过 hard-ALM；
6. 可选加入 \(\hat C\) tensor norm 的 L2 boundedness rows，权重记为 `branch_bound_weight`。

注意：boundedness rows 是数值测试工具，不是新物理项本身。物理上对应的是“\(\hat C\) 有界”这一 branch 正则性条件。

## 输出

主要输出目录：

```text
visualizations/trace0_branch_hard_alm_n384_tau0_core10/
visualizations/trace0_branch_hard_alm_n384_tau0_core10_midscan/
visualizations/trace0_branch_hard_alm_n384_taum3p5_core10_bw3em4/
visualizations/trace0_branch_hard_alm_n384_taup3p5_core10_bw0/
visualizations/trace0_branch_hard_alm_n384_taup3p5_core10_bw1em4/
visualizations/trace0_branch_hard_alm_n384_taup3p5_core10_bw3em4/
```

每个目录包含：

```text
summary.json
trace0_branch_conservation_scan.png
trace0_branch_conservation_coefficients.npz
```

## \(\tau=0\) 扫描结果

`n=384, core10`，中心干涉切片。

| branch bound | algebraic wmean | algebraic p95 | full div wmean | \(\hat C\) wmean | \(\hat C\) p95 |
|---:|---:|---:|---:|---:|---:|
| 0 | 8.88e-3 | 3.22e-2 | 2.65e-5 | 1.76e6 | 8.77e4 |
| 1e-4 | 1.67e-2 | 4.26e-2 | 2.69e-5 | 1.57e4 | 6.86e4 |
| 2e-4 | 2.09e-2 | 4.96e-2 | 2.78e-5 | 1.10e4 | 6.00e4 |
| 3e-4 | 2.43e-2 | 6.28e-2 | 2.91e-5 | 8.69e3 | 5.71e4 |
| 5e-4 | 2.95e-2 | 9.22e-2 | 3.36e-5 | 6.16e3 | 5.00e4 |
| 1e-3 | 3.73e-2 | 1.46e-1 | 5.34e-5 | 3.76e3 | 3.24e4 |
| 1e-2 | 7.54e-2 | 9.24e-1 | 2.56e-3 | 8.24e2 | 4.94e3 |

结论：

- \(\tau=0\) 存在一段可用窗口，大约 `1e-4` 到 `5e-4`；
- 在此窗口中，守恒残差几乎不变，代数残差从约 0.009 上升到 0.017--0.029，但仍可控；
- \(\hat C\) 的 weighted mean 从 `1.76e6` 降到 `6e3--1.6e4`；
- 太强的有界性惩罚 `1e-2` 会把 algebraic p95 推到接近 1，不可接受。

这说明 branch-v1.1 不是立即失败；但 \(\hat C\) 有界性和拟合 A-reference 快照之间存在真实张力。

## \(\tau=-3.5\) 左侧分离态

测试权重：`branch_bound_weight=3e-4`。

结果：

- algebraic weighted mean：`2.21e-2`；
- algebraic p95：`7.59e-2`；
- full divergence weighted mean：`1.00e-3`；
- full divergence p95：`8.85e-4`；
- \(\hat C\) weighted mean：`5.48e5`；
- \(\hat C\) p95：`1.79e6`；
- \(\chi\) p50：`2.62e-3`。

判断：

- 代数残差仍可控；
- full divergence 与此前分离态 hard-ALM 的量级相近；
- 但因为大量点 \(\chi\) 极小，\(\hat C\) 仍非常大。

这说明分离态对 branch boundedness 更苛刻。若坚持同一个 \(\chi\) 归一化和全局 bound 权重，左侧并没有完全失败，但 \(\hat C\) 有界性还不够理想。

## \(\tau=+3.5\) 右侧分离态

测试了 `branch_bound_weight=0,1e-4,3e-4`。

| branch bound | algebraic wmean | algebraic p95 | full div wmean | full div p95 | \(\hat C\) wmean | \(\hat C\) p95 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1.33e-2 | 6.27e-2 | 3.27e-4 | 1.55e-3 | 2.35e8 | 3.54e6 |
| 1e-4 | 2.80e-2 | 1.07e-1 | 4.87e-3 | 3.18e-3 | 2.53e5 | 1.02e6 |
| 3e-4 | 3.51e-2 | 1.46e-1 | 3.25e-2 | 6.60e-3 | 1.23e5 | 6.90e5 |

判断：

- `bound=0` 时 full divergence 尚可，但 \(\hat C\) 极大；
- 一旦压 \(\hat C\)，full divergence 明显恶化；
- 右侧分离态比 \(\tau=0\) 和左侧分离态更紧张。

这说明当前 boundedness rows 的实现不能简单全局套用；右侧可能需要：

1. 更好的 ALM/KKT 线性代数；
2. 按 \(\chi\)、局部密度或 patch 类型自适应的 boundedness 权重；
3. 重新审视 \(\chi(\mathcal q)\) 的尺度选择，而不是每个切片只用 p95 归一化；
4. 或者承认强 \(C=\chi\hat C\) 正则分支与 A-reference 快照在分离态存在物理差异。

## 当前理论判断

本轮支持把 proposal 改写为 v1.1：

\[
\boxed{
\mathcal R_{\mu\nu}
=
\chi(\mathcal q)\hat C_{\mu\nu},
\qquad
\hat C_{\mu\nu}\in E,
\qquad
\tilde g^{\mu\nu}\hat C_{\mu\nu}=0,
\qquad
\tilde\nabla^\mu(\chi\hat C_{\mu\nu})=0.
}
\]

但需要明确：

- \(\hat C\) 有界不是自动结果，必须作为正则性条件；
- 数值实现中不能只用一个全局 constant penalty；
- 右侧分离态已经暴露出 branch boundedness 与 conservation 的张力；
- 因此 v1.1 仍是候选，不是完成理论。

## 下一步

1. 把 branch boundedness 从“普通 L2 惩罚”改成更接近理论的局部条件，例如：
   \[
   \|\hat C\|_W \le K
   \]
   的不等式/投影形式，而不是全局二次惩罚。
2. 检查右侧分离态的 full divergence outliers 位置，判断它们是否集中在 \(\chi\ll1\)、低密度边缘或 patch 退化区域。
3. 固定一个物理上更合理的 \(\chi\) 尺度，而不是每个切片用 p95 重新归一化。
4. 若右侧问题仍存在，应在理论上允许 A-reference 快照投影到 D-branch 初态后产生可观但受控的偏差，而不是要求 A 快照本身逐点满足强 branch 正则性。

