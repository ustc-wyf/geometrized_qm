# D 支分支保持 ADM 小步松弛闭合

日期：2026-05-03

## 背景

在前一轮 ADM 分块延拓中，候选 metric 没有被接受，表面诊断为：

\[
\texttt{metric\_extension\_admissible\_blend}=0 .
\]

进一步检查后发现，这个结论混合了两个问题：

1. ADM 无效点被 floor 后的伪 ADM 变量重构，导致粗糙度出现 \(10^{53}\) 量级伪影；
2. admissible 检查强制所有 trusted 点都是 \(+--\) 签名，而初态 \(\tilde g\) 在实验室时间切片下本来已有 7 个 trusted 点处在非 \(+--\) 惯性分支。

这与用户允许 \(\tilde g\) 出现时间反向/非标准分支的理论口径不一致。

## 代码修正

修改文件：

- `kg_examples/simulate_d_tridomain_full_dynamics.py`

新增或修正：

- `--metric-extension-adm-update-fields all|lapse|lapse_shift`
- `--metric-extension-max-rel-change`
- `inertia_signature`
- `branch_preserving_signature_mask`
- ADM 无效点保留原 metric，不再用 floor 后变量重构；
- admissible 签名检查从“强制 \(+--\)”改为“保持该点原有惯性分支不变”。

含义：

\[
\boxed{
\text{代表选择可以接受已有分支，但不能由数值更新制造新的签名/惯性类跳变。}
}
\]

## 诊断

初态 trusted 区：

- trusted 点数：128；
- 其中 121 个为 \(+--\)；
- 7 个点不是 \(+--\)，但更新候选在这些点 `delta_norm=0`；
- 因此原来的 `admissible_blend=0` 不是候选更新破坏这些点，而是验收逻辑把已有分支误判为候选失败。

回溯到 \(2^{-20}\) 时仍旧失败，也验证了旧逻辑的问题：

- 失败点数始终是 7；
- 失败原因都是 signature；
- 失败点的候选更新为 0。

## 小步松弛

直接接受平直边界锚点驱动的 ADM 延拓会使 D-A transformed measure 短时偏离过大：

- `steps=2`
- `lapse` uncapped：
  - `measure_rel_l1_trusted ≈ 3.396e-2`

因此加入候选 metric 相对变化上限：

\[
\frac{\|\delta \tilde g\|}{\|\tilde g\|} \le \epsilon_g .
\]

该上限只限制饱和区代表 metric 的选择步长，不改物质方程、质量壳方程、源项或 D 支作用量。

## 20 步结果

测试设置：

- 1550 nm；
- KG 归一化；
- \(m=\omega_0/10\)；
- \(\ell/l_P=10^{60}\)；
- `t_old=8`；
- `n=96`；
- `dt_old=2.5e-7`；
- `steps=20`；
- `geometry_closure=restricted_extension`；
- `metric_extension_variable_mode=adm`；
- `metric_extension_max_rel_change=1e-3`。

对照结果：

| run | trusted measure rel L1 | support measure rel L1 | trusted \(D_{\min}\) | accepted blend | interface solved fraction |
|---|---:|---:|---:|---:|---:|
| inert/orphan | \(1.1488\times10^{-4}\) | \(1.1398\times10^{-4}\) | \(4.0996\times10^{-7}\) | 0 | 0.376 |
| ADM lapse cap \(10^{-3}\) | \(2.7536\times10^{-4}\) | \(2.6229\times10^{-4}\) | \(4.0996\times10^{-7}\) | 1 | 0.450 |
| ADM lapse+shift cap \(10^{-3}\) | \(6.6356\times10^{-4}\) | \(6.2529\times10^{-4}\) | \(4.0996\times10^{-7}\) | 1 | 0.431 |

## 判断

当前最稳闭合候选是：

\[
\boxed{
\text{branch-preserving signature check}
+ \text{ADM lapse-only}
+ \text{max relative metric update }10^{-3}.
}
\]

它的优点：

- 更新被 admissible 检查接受；
- 不制造新的签名分支跳变；
- 质量壳逐点判别式裕度保持；
- D-A transformed measure 只从 inert 的 \(1.15\times10^{-4}\) 增至 \(2.75\times10^{-4}\)，没有爆炸；
- interface solved fraction 从 0.376 提高到约 0.450。

它的不足：

- 仍然是代表选择闭合，不是完整张量场方程的唯一解；
- `adm_valid_fraction≈0.258`，说明固定实验室时间 ADM 坐标只覆盖 active 区一部分；
- `lapse+shift` 比 `lapse` 偏离 A 更大，当前不应作为默认优先闭合。

## 下一步

下一步应以 `ADM lapse cap 1e-3` 为当前 working closure，继续做：

1. 更长时间步数测试；
2. `max_rel_change=1e-4,3e-4,1e-3,3e-3` 扫描；
3. 与 direct tensor interface matching 联合；
4. 对 ADM invalid 区画图，确认它们对应时间反向/非实验室切片区域还是数值尾部。

