# 129. D 支同初态 reduced 动态演化首轮

日期：2026-05-02

## 用户校正后的问题

用户指出，不能继续把 A 分支平直量子力学参考解静态代入 D 支 jump law，然后把不匹配解释成 D 支失败。真正要做的是：

- 以与 A 支相同的实验准备态作为初值；
- 按 D 支的物质方程演化 \(\rho_D,S_D\)；
- 让 \(\tilde g_D\) 与过渡界面按 D 支方程/jump law 共同演化；
- 再把 D 支结果与 A 支结果比较。

本轮实现的是第一步可执行动态原型，不是最终完整 D 支求解器。

## 新增与修改

- 新脚本：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/simulate_d_reduced_dynamic_same_initial.py`
- 修改：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/coordinate_matter_evolution.py`
- 修改内容：
  `coordinate_matter_rhs` 与 `coordinate_matter_rhs_covector` 现在返回已计算的 `flux_x/flux_z`，供守恒型密度通量使用。

脚本含义：

- 初态使用同一组局域双高斯波包 A 准备态；
- 物质变量用 D 支坐标度规下的 covector 形式演化；
- `frozen_initial` 模式固定初始 \(\tilde g\)，用于分离物质方程稳定性；
- `instant_transform` 模式每步用当前 \(\rho,u_\mu\) 代数重构 \(\tilde g[\rho,u]\)，作为反例/诊断；
- jump/interface law 当前作为动态诊断输出，还不是完整 tensor jump feedback。

## 关键实现修正

最初的守恒型 Rusanov 通量把 `n_cons=\sqrt{|\tilde g|}j^t` 误当作必须非负的密度做了裁剪，这是错误的。负频分支中 `n_cons` 可以带符号，真正应保持非负的是恢复出的 \(\rho\)，不是 `n_cons`。

修正后：

- Rusanov 速度使用 \(v^i=F^i/q\)，其中 \(q=n_\mathrm{cons}\)；
- 不再把 `n_cons` 裁剪成非负；
- 默认 `active_dilation=0`，避免把低密度、\(\tilde g\) 条件数病态的尾部强行纳入演化；
- 默认 `interface_every=5`，界面/jump 诊断降频，不改变物质演化；
- 默认在 support 上质量壳判别式变负时停止，避免继续用 `disc_floor` 伪造实数 \(u_t\)。

## 数值结果

### frozen_initial

命令：

```bash
python3 kg_examples/simulate_d_reduced_dynamic_same_initial.py \
  --output visualizations/d_reduced_dynamic_same_initial_frozen_n64_dt1e6_active0_valid74 \
  --metric-mode frozen_initial \
  --resolution 64 \
  --dt 1e-6 \
  --steps 74 \
  --ell 30 \
  --samples 21 \
  --active-dilation 0 \
  --interface-every 10
```

输出：

- 图：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_frozen_n64_dt1e6_active0_valid74/d_reduced_dynamic_vs_a_frozen_initial_ell30_n64_t7.4e-05.png`
- 摘要：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_frozen_n64_dt1e6_active0_valid74/summary.json`

最后可靠时间片：

- `t = 7.4e-05`
- `rho_rel_l1_support = 1.5836385e-3`
- `rho_max_abs_support = 1.4297401e-2`
- `phase_grad_rel_l1_support = 5.878061995e-5`
- `disc_min_support = 3.52259377e-8`
- `interface_solved_fraction = 0.50617`

若继续到 step 75，`disc_min_support` 变成 `-4.629266e-7`，即在当前实验室时间切片下质量壳方程无实数 \(u_t\)。因此 step 74 是本轮 frozen baseline 的最后可靠比较时间片。

### instant_transform

命令：

```bash
python3 kg_examples/simulate_d_reduced_dynamic_same_initial.py \
  --output visualizations/d_reduced_dynamic_same_initial_instant_n64_dt1e6_active0_stopdisc \
  --metric-mode instant_transform \
  --resolution 64 \
  --dt 1e-6 \
  --steps 74 \
  --ell 30 \
  --samples 21 \
  --active-dilation 0 \
  --interface-every 10
```

输出：

- 图：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_instant_n64_dt1e6_active0_stopdisc/d_reduced_dynamic_vs_a_instant_transform_ell30_n64_t2e-06.png`
- 摘要：
  `/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_reduced_dynamic_same_initial_instant_n64_dt1e6_active0_stopdisc/summary.json`

结果：

- step 1 仍接近 A，`rho_rel_l1_support = 2.50e-5`；
- step 2 立刻病态，`disc_min_support = -2.1399878545e11`，\(\rho\) 出现 \(10^{35}\) 量级异常。

解释：

每步代数重构 \(\tilde g[\rho,u]\) 需要 \(Q=\Box\sqrt{\rho}/\sqrt{\rho}\)，其中包含 \(\partial_t^2\sqrt{\rho}\)。用动态 \(\rho\) 的有限差分重构该二阶时间导数，会把微小时间离散误差按 \(1/\Delta t^2\) 放大。因此 `instant_transform` 不是可用的全动力学路线。

## 当前结论

1. 保符号守恒密度后，D 支物质方程短时可以从 A 初态稳定启动。
2. 在冻结初始 \(\tilde g\) 的基线下，D 与 A 的早期密度差先保持在 \(10^{-5}\sim10^{-4}\) 级，接近质量壳判别式穿零前增长到 \(1.6\times10^{-3}\)。
3. 实验室时间切片下的质量壳判别式穿零是当前 reduced 动态的硬边界，不能用 `disc_floor` 继续解释为物理演化。
4. 每步从 \(\rho\) 代数重构 \(\tilde g\) 的方法不可行；完整 D 求解器必须把 \(\tilde g\) 或界面/层变量作为独立演化对象，而不是用 \(Q\) 的二阶时间差分反复重构。
5. trace jump law 仍只作为诊断；根据此前决策，完整主线必须转向 direct tensor jump / body-fitted interface 与独立几何演化。

## 重要更正：本轮没有正确使用 \(\tilde\rho\) 初态

用户指出，本轮 reduced 动态中虽然从 A 支求出了 \(\rho_A,u_\mu\)，但代码并没有先把 \(\rho_A\) 变换为 \(\tilde\rho\)。复查确认该指出是正确的。

当前代码实际做的是：

\[
n_{\rm used}
=
\sqrt{|\tilde g|}\,\rho_A\,\tilde g^{t\nu}u_\nu .
\]

但在当前混合变换、非零质量分支 \(\kappa=m^2\) 的强匹配口径下，应优先使用测度密度关系：

\[
\kappa\sqrt{|{\tilde g}|}\,\tilde\rho
=
X\sqrt{|g|}\rho
\]

或者若坚持正测度密度：

\[
\sqrt{|{\tilde g}|}\,\tilde\rho
=
\frac{|X|}{m^2}\sqrt{|g|}\rho .
\]

在本轮平直 \(g=\eta\) 初态下，这意味着应先构造：

\[
\tilde\rho_0
=
\frac{|X_0|}{m^2}
\frac{\rho_A}{\sqrt{|\tilde g_0|}}
\]

再设

\[
n_{\rm correct}
=
\sqrt{|\tilde g_0|}\,\tilde\rho_0\,\tilde g_0^{t\nu}u_\nu
=
\frac{|X_0|}{m^2}\rho_A\,\tilde g_0^{t\nu}u_\nu .
\]

因此本轮 `frozen_initial` 的数值结果只能保留为“错误密度映射下的数值稳定性/方法学诊断”，不能作为 D 支同初态物理演化与 A 支的有效比较。下一轮必须修正 \(\rho_A\to\tilde\rho_0\) 和 observable 的回表象定义。

## 下一步

- 把当前 reduced matter evolution 保留为物质子系统；
- 不再使用 `instant_transform` 作为主算法；
- 修正初态密度映射：用 \(\tilde\rho_0\) 而不是 \(\rho_A\) 初始化 D 支物质方程；
- 比较时必须说明比较的是裸 \(\tilde\rho\)、测度密度 \(\sqrt{|\tilde g|}\tilde\rho\)、流 \(J^\mu\)，还是拉回到 \(g\) 表象后的有效密度；
- 进入 D 支一致初值/演化的下一版：
  - \(\tilde g\) 或 scalaron/interface 作为独立变量；
  - direct tensor jump law 给出界面法向/速度/生成消失条件；
  - 质量壳判别式穿零处需要更换时间切片、边界匹配或 body-fitted 坐标，不能简单 floor。
