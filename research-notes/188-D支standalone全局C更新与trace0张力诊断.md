# D 支 standalone 全局 \(C\) 更新与 trace0 张力诊断

日期：2026-05-08

## 1. 做了什么

在 `kg_examples/simulate_d_harmonic_standalone.py` 中新增了 runtime 全局 \(C\)-update：

- `--auxiliary-update global`
- `--auxiliary-global-mode penalty|project|alm`
- `--auxiliary-trace0`
- `--global-conservation-weight`

目标不是冻结 \(C\)，而是每一步重新求 \(C_{\mu\nu}\)。

## 2. 诊断结果

固定算例：

- `tau=0`
- `n=384`
- `core10`
- 输入包：`visualizations/equation_first_gbcd_plus_initial_package_n384_tau0_core10_noactiveedge_guarded/gbcd_plus_initial_package.npz`

结果：

- `local + no trace0`：2 步后 D residual core10 weighted mean 约 `0.041`
- `global + trace0 + penalty`：D residual core10 weighted mean 约 `0.103`
- `global + trace0 + hard project`：守恒残差压到 `1.8e-5`，但代数残差炸到 `1e6` 量级，coeff 改变量约 `323`
- `global + trace0 + ALM`：守恒更好，但代数残差仍明显劣化
- `global + trace0 + linear-plus metric corrector`：把 D residual 和 \(\rho\) pullback 偏差都显著弄坏

## 3. 解释

这说明：

1. 冻结 \(C\) 只是诊断工具，不是物理规则。
2. 在当前“固定未来层度规 + 同步求 \(C\)”的离散结构里，`trace0 + full conservation` 不能被单独硬投影成唯一可行解。
3. 真正需要同步的是：
   - \(C\)-sector 更新，
   - metric 未来层更新，
   - 以及 trace0/守恒/代数方程之间的耦合。

换句话说，当前瓶颈已经从“有没有更新 \(C\)”变成“\(C\) 与 metric 是否要联合求解”。

## 4. 下一步

- 不再把冻结 \(C\) 当主路线；
- 不再单独硬压 trace0；
- 先做 \(C\)-metric 联合更新的最小原型，或在理论上明确 trace0 是否应改成 patch 条件而非全局硬约束。

