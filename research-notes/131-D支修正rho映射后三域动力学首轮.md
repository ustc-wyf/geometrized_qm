# D支修正 rho 映射后三域动力学首轮

日期：2026-05-02

## 目的

在用户指出 `rho_A != rho_tilde` 后，重新启动 D 支同初态动力学计算：

- 初态不再使用 `rho_tilde=rho_A`；
- 使用正测度密度关系
  \[
  \sqrt{|\tilde g|}\tilde\rho={|X|\rho_A\over m^2}
  \]
  初始化 \(\tilde\rho\)；
- 比较对象改为变换后测度密度 \(\tilde N=\sqrt{|\tilde g|}\tilde\rho\) 和坐标守恒密度
  \[
  n_{\rm cons}=\sqrt{|\tilde g|}\tilde\rho\,\tilde g^{t\nu}u_\nu .
  \]

## 新脚本

- `kg_examples/simulate_d_tridomain_full_dynamics.py`

脚本内容：

- 用局域双高斯波包构造 A 参考初态；
- 用当前 disformal/mixed 变换构造初始 \(\tilde g_0\)；
- 用 \(\tilde N_0=|X_0|\rho_A/m^2\) 得到 \(\tilde\rho_0=\tilde N_0/\sqrt{|\tilde g_0|}\)；
- 以 \((\tilde\rho,u_i,\tilde g)\) 构造守恒量 \(n_{\rm cons}\)；
- 用 \(\tilde g\) 表象质量壳方程和守恒方程推进物质子系统；
- 计算 \(y=\ell^2\tilde R\)，按 weak / transition / saturated 三域统计；
- 用修正后的 \(\tilde\rho\) 计算 \(\tilde T_{\mu\nu}\) 与 plateau tensor diagnostic；
- 将 direct/trace jump 相关界面提取作为诊断输出。

## 本轮几何闭合

本轮使用 `geometry_closure=inert_tridomain`：

- weak/saturated bulk 在去掉高阶导数项后没有唯一局部几何演化；
- 因此本轮选择“惯性 bulk representative”：保持当前选定 \(\tilde g\)，只推进物质子系统；
- transition/interface 条件仍只作为诊断，还没有反过来唯一求解 bulk metric。

这不是最终完整张量几何求解器，而是修正 \(\tilde\rho\) 后的第一版 rule-closed 可执行系统。

## 运行结果

正式稳定输出：

- 输出目录：
  `visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/`
- 图：
  `visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/d_tridomain_dynamics_inert_tridomain_ell300_n64_t5.7e-05_fixed.png`
- 数据：
  `visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/fields_final.npz`
- 摘要：
  `visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/summary.json`

参数：

- `resolution=64`
- `dt=1e-6`
- `ell=300`
- `M_P=300`
- `steps=57`
- `t_completed=5.7e-05`

关键结果：

- `measure_rel_l1_support = 2.2905e-2`
- `n_cons_rel_l1_support = 5.9187e-5`
- `disc_min_support = 4.5188e-7`
- `rho_tilde_max_support = 0.476161`
- `mass_shell_defect_p95 = 3.95e-14`
- `weak_fraction = 0`
- `transition_fraction = 0`
- `saturated_fraction = 1`
- `interface_solved_fraction = 0.517`
- `plateau_lhs_p95 = 5.56e-3`
- `plateau_rhs_p95 = 3.48e-4`
- `plateau_relative_p95 = 0.999`

对照失败点：

- 同样参数尝试跑到 `steps=100` 时，在 `step=58` 出现
  `disc_min_support=-4.627e-8`；
- 因此本轮可信截止采用 `step=57`，不展示 `step=58` 后的病态尖峰。

## 解释

1. 这次 `t=0` 的 \(\tilde N\) 与 A 参考变换后测度密度吻合到机器精度，说明此前的 `rho_A=rho_tilde` 错误已被修正。

2. \(n_{\rm cons}\) 与 A 参考的相对差到截止时仍只有 \(O(10^{-5})\)，而 \(\tilde N\) 的相对差到 \(2.3\%\)。这说明坐标守恒流短时非常稳定，但测度密度对参考 \(\tilde g_A(t)\) 与当前惯性 \(\tilde g_D\) 的差异更敏感。

3. 对 `ell=300`，初始支撑区已经全部落在 saturated 区。快速检查显示即使 `ell=1`，支撑区仍有约 `81%` saturated；`ell>=10` 时基本全 saturated。这意味着当前局域双高斯的 \(\tilde R\) 在该变换下本来就非常大。

4. `plateau_relative_p95≈1` 不是数值爆炸，而是张量结构不匹配：\(-\frac12 f\tilde g_{\mu\nu}\) 是满秩 metric 型，而 \(\tilde T_{\mu\nu}\) 近似 rank-one。只有两边都足够小时，饱和区才可近似当作 \(0\simeq0\)。

5. 本轮 `ell=M_P=300` 给出 \(M_P^2/\ell^2=1\)，plateau 项并不被参数压低。如果要实现“饱和区近似无动力”的物理设想，需要后续扫描 \(\ell\gg M_P\) 或明确引入背景/常数项处理。

## 未决问题

- `inert_tridomain` 仍然不是唯一由作用量推出的完整几何演化；
- 饱和区 bulk metric 的内禀规则仍需理论指定；
- transition/interface 的 direct tensor jump 还没有作为强制边界条件反推 metric；
- 质量壳判别式在 `t≈5.8e-05` 穿零，需要换切片、body-fitted 坐标或界面匹配处理；
- 旧 residual/jump/tensor 诊断中使用 `rho_A` 的源项仍需系统性重算。
