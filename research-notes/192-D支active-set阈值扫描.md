# D 支 active-set 阈值扫描

日期：2026-05-08

## 1. 目的

上一轮 `mass_shell_guard` 使用：

- `active_set_disc_margin = 1e-7`
- `active_set_rho_frac = 5e-3`

成功清除了 active 区负质量壳判别式和负 `rho_tilde`。本轮检查该结论是否依赖单个阈值。

## 2. 扫描设置

共同设置：

- `tau = 0`
- `steps = 10`
- `evolve-region = support`
- `active-set-mode = mass_shell_guard`
- `global + trace0 + penalty + joint2`
- `metric-corrector = none`

扫描组合：

| name | disc margin | rho frac |
|---|---:|---:|
| baseline | `1e-7` | `5e-3` |
| tighter disc | `5e-8` | `5e-3` |
| looser disc | `2e-7` | `5e-3` |
| stricter density | `1e-7` | `1e-3` |

## 3. 结果表

| run | dropped | core missing | active neg disc | active neg rho | D residual core10 wmean | core10 rho wL1 | support rho wL1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `1e-7,5e-3` | 3 | 0 | 0 | 0 | `0.1002552` | `2.0539e-4` | `2.1369e-4` |
| `5e-8,5e-3` | 2 | 0 | 0 | 0 | `0.1002546` | `2.0538e-4` | `2.1365e-4` |
| `2e-7,5e-3` | 12 | 0 | 0 | 0 | `0.1002570` | `2.0539e-4` | `2.1412e-4` |
| `1e-7,1e-3` | 2 | 0 | 1 | 1 | `0.1002781` | `2.0538e-4` | `2.3296e-4` |

## 4. 判断

1. `rho_frac=1e-3` 太严格，会留下 active 区负判别式和负 `rho_tilde`。
2. `rho_frac=5e-3` 是当前需要的低密度门限量级。
3. `disc_margin=5e-8` 已经足够清除 active 坏点，并且只移除 2 个点，比 `1e-7` 的 3 个点和 `2e-7` 的 12 个点更小。
4. 三个有效组合都保持 `core10` 完整，且 core10 残差和 `rho` pullback 几乎不变。

## 5. 当前默认建议

当前 support 边界/interface 候选参数更新为：

- `active_set_mode = mass_shell_guard`
- `active_set_disc_margin = 5e-8`
- `active_set_rho_frac = 5e-3`

这仍是数值计算域规则，不是理论方程本身。

## 6. 下一步

- 用该默认参数检查 `tau=-3.5,+3.5`；
- 如果三切片都稳定，再考虑更长时间窗；
- 同时开始优化全局 `C` 求解器，避免每 10 步都耗时数分钟。
