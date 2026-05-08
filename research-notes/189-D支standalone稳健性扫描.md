# D 支 standalone 稳健性扫描

日期：2026-05-08

## 1. 扫描设置

默认基线：

- `kg_examples/simulate_d_harmonic_standalone.py`
- `--auxiliary-update global`
- `--auxiliary-trace0`
- `--auxiliary-global-mode penalty`
- `--joint-outer-iterations 2`
- `--joint-metric-relaxation 1.0`
- `--joint-residual-gate-rel 0.0`
- `--metric-corrector none`

扫描切片：

- `tau = -3.5, 0, +3.5`
- `steps = 10`
- 额外补 `tau = 0, steps = 20`

## 2. 结果

### `tau = 0`

- 10 步：D residual core10 weighted mean `0.0987`
- 20 步：D residual core10 weighted mean `0.0977`
- 结论：长窗没有明显漂移，20 步与 10 步差别很小

### `tau = -3.5`

- 10 步：D residual core10 weighted mean `0.1176`
- \(\rho\) pullback core10 weighted L1 `1.94`
- 结论：残差仍在可控量级，但与平直/干涉中点相比误差更大

### `tau = +3.5`

- 10 步：D residual core10 weighted mean `0.0911`
- \(\rho\) pullback core10 weighted L1 `2.25`
- 结论：残差同样稳定，但 D-A 偏差在分离态更大

## 3. 解释

这个基线的关键信息是：

1. `trace0 + penalty + joint2` 是稳定的，不会在 10 到 20 步内自发爆掉。
2. `tau=0` 的长窗漂移很小，说明联立固定点已经可以当作短窗骨架。
3. 分离态的 \(\rho\) pullback 偏差更大，这更像“物理解偏差”而不是数值失控。

## 4. 结论

- 可以把这条基线当成当前默认稳定数值器骨架；
- 但它还不是最终生产版，因为：
  - `linear-plus` 仍不可靠；
  - 边界/interface 规则还未完全定型；
  - 长窗和更复杂初值还需继续验证。

