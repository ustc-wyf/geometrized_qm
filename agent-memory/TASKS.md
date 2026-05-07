# 项目任务

## 总目标

- 构造并检验一个“量子势几何化”的最小理论框架：在不改变实验观测结果的前提下，把 Klein-Gordon 标量场中的量子势项通过合适的度规变换吸收到几何结构中

## 目标树

- 根目标：量子势几何化
- 分支 A：建立最小理论设定
- A1：固定记号、单位制、签名、Madelung 分解约定
- A2：写出自由 KG 场在背景时空中的作用量与场方程
- A3：把作用量改写为 `\rho` 与 `S` 的形式，明确量子势项的精确表达
- 分支 B：分析问题动机与目标判据
- B1：明确 Bohm/KG 图景下负粒子数密度与超光速问题出现在哪里
- B2：明确“几何化成功”的判据：作用量消去显含量子势项、轨迹变为测地线、观测不变
- B3：明确“物理可接受”的判据：协变性、可逆性、因果性、经典极限
- 分支 C：构造几何变换方案
- C1：确定 Bekenstein disformal transform 的一般形式与允许依赖的变量
- C2：判断变换应依赖 `\rho`、`S`、`\partial_\mu S`、`\partial_\mu \rho` 中哪些对象
- C3：求解使量子势项被吸收进新度规的条件
- 分支 D：检验变换后的动力学
- D1：推导新框架下的作用量与 Euler-Lagrange 方程
- D2：检查粒子轨迹能否写成新度规下的测地线方程
- D3：检查电流、密度、类时/类空性质如何在新几何中重解释
- 分支 E：一致性与极限
- E1：检查经典极限 `Q -> 0` 或缓变极限下 `\tilde g_{\mu\nu} -> g_{\mu\nu}`
- E2：检查实验观测等价性，不把纯数学重写误认为新物理预测
- E3：梳理未来与 Einstein-Hilbert 动力学衔接的条件，但暂不深入展开
- 分支 F：样例与反例
- F1：在闵氏背景上选最简单可算例子测试公式
- F2：寻找会导致变换失效、奇异或不可逆的反例
- F3：总结该方案适用范围

## 当前进展

- 已创建长期记忆目录与四个核心记忆文件
- 已将长记忆状态设为开启
- 已明确课题名称、研究动机、最小研究范围和候选变换类型
- 已形成第一版树状研究计划
- 已完成第一步正式推导：KG 作用量的 Madelung 分解与量子势定位
- 已完成第二步判据整理：成功、失败与物理可接受性标准
- 已完成第三步筛选：disformal ansatz 候选分类与第一轮优先级排序
- 已完成第四步比较：共形、相位单向、振幅单向、最小混合四类变换的初步并排分析
- 已完成第五步直接计算：四类变换代回作用量与 HJ 方程后的显式结构
- 已完成第六步：HJ 层面吸收 `Q` 的显式求解，并得到更广泛变换族的存在性结论
- 已完成第七步快速筛选：最简单 HJ 解的作用量兼容性与存在条件检查
- 已完成第八步方法论修正：确认“作用量逐项等价”不是必要条件，后续以物理方程等价为优先筛选标准
- 已完成第九步变分筛选：最简单共形方案破坏原连续性方程，类 II 的 `A=1` 最简单解导致零作用量退化
- 已完成第十步修正检验：在 `sqrt(-g)ρ = sqrt(-g~)ρ~` 约定下重算，类 II 退化结论保持不变，共形连续性问题仍存在
- 已完成第十一步方法修正：确认应“先变分再代回变换”，并撤回对类 II 最简单解的直接排除
- 已完成第十二步结构筛选：一般类 II 在当前作用量 ansatz 下泛型地不能同时满足原 HJ 方程与原连续性方程
- 已完成第十三步复核：确认在纯类 II 且只由 HJ 方程决定时，`A=1` 的 `B` 不被强制含 `Z`
- 已完成第十四步：按用户指定的逆度规 ansatz 直接代回拉氏量，得到消去 `Q` 项的最简单解与一般解族
- 已完成第十五步：检验最简单消项解，确认其成功重现 HJ 方程并支持测地线解释，但一般不能保留原连续性方程
- 已完成第十六步复核：把 `X=m^2-Q` 明确代入后，确认连续性方程差异等价于额外条件 `J^mu ∂_mu Q = 0`
- 已完成第十七步比较：确认 `J^mu ∂_mu Q = 0` 与新几何测地线方程不等价
- 已完成第十八步 no-go：最小混合逆度规 ansatz 下，若要求拉氏量精确消项且连续性方程原样重合，则必有 `Q=0`
- 已完成第十九步修正：去掉 `sqrt(-g)ρ` 不变假设后，纯 `u^mu u^nu` ansatz 可通过适当的 `ρ~` 变换同时重现 HJ 与连续性方程
- 已完成第二十步：最简单纯 ansatz 在 `X>0` 区域内同时满足可逆性、Lorentz 签名保持、实幅度与正确经典极限
- 已完成第二十一步分支分析：`X>0` 与 `X<0` 不对应“成功/失败”，而对应两条不同的纯几何分支
- 已完成第二十二步：对 Cursor 的对照分析做阶段性判断，确认关键差异在等价性定义与密度变换假设
- 已完成第二十三步：形成阶段性总结，并建立 `S_EH[g]` / `S_EH[g~]` 的 Einstein 作用量比较框架
- 已完成第二十四步：给出 `S_EH[g]` 与 `S_EH[g~]` 的一阶/二阶差异主项，并明确弱场与强量子势场景下的差异量级
- 已完成第二十五步：建立 1+1d 自由 KG 双高斯波包基线模型并跑出第一轮 `Q` 的量级结果
- 已完成第二十六步：用双高斯 KG 具体例子检验四种一致性机制，排除机制 1 和机制 2，保留机制 3 和机制 4
- 已完成第二十七步：提出 `Q` 依赖引力作用量的一般要求，并分析 `exp(-(κQ)^2)R~`、`1/(1+(κQ)^2)R~`、`ρ~R~` 的可行性
- 已完成第二十八步：在双高斯模型上具体测试 `F(Q)=1/(1+(κQ)^2)`，确认主支撑区较温和而节点导数项更危险
- 已完成第二十九步：在 1+1d 双高斯模型上直接比较 `EH[g]` 与 `EH[g~]`，确认其差异主要体现在局域几何量而非纯 EH 的 bulk 动力学预测
- 已完成第三十步：建立多版 3+1d 模型原型，并确认若干近似/错误几何不能作为正式主算例
- 已完成第三十一步：三种完整作用量已写成统一约束形式
- 已完成第三十二步：C 分支已完成辅助场重写
- 已完成第三十三步：3+1d→2+1d 一致约化框架已经建立
- 已完成第三十四步：在 `2+1d` conformal-Killing 约化下写出 `R[g~]` 与 `F(R~)R~` 的显式动力学方程，并实现第一版全动力学数值器
- 已建立独立于长期记忆的 `research-notes/` 留档目录
- 已确认本项目是从另一台电脑复制到当前机器；旧的 `agent-memory/` 作为权威历史保留，不做覆盖式重建
- 已在当前机器重新接管长记忆协议；后续每轮都应先按 `LOG.md` 中的最新状态执行 memory reload
- 已完成共享 baseline 下更密的 `lambda_grav` continuation 扫描，并得到新的代表区间：`B_vs_A` 约在 `lambda~1` 进入 `1e-6`，`C_vs_B` 约在 `lambda~3` 稳健进入 `1e-5`
- 已完成从 A 支 Einstein 方程变换得到 \(\mathcal R_{\mu\nu}\) 的分解路线，并修正反解度规记号为 \(\mathfrak g_{\mu\nu}\)，避免和 Einstein tensor \(G_{\mu\nu}\) 混淆。
- 已实现并运行 `kg_examples/diagnose_mathcal_r_pure_geometry.py`：
  - 高分辨率输出 `visualizations/mathcal_r_pure_geometry_n320_3tau/`；
  - 完整 \(\mathcal R_{\rm need}/M_P^2\) 的 pure-geometry 拟合成功是平凡 EH 主项；
  - 非平凡源 \(\tilde T_{\mu\nu}/M_P^2\) 的 pure-geometry 拟合 residual 约 `0.996`，说明不能据此恢复 pure-\(\tilde g\) action 路线。
- 已实现并运行 `kg_examples/diagnose_mathcal_r_gtilde_function.py`：
  - 高分辨率输出 `visualizations/mathcal_r_gtilde_function_n320_3tau/`；
  - 留一时间切片 kNN 检验显示：非平凡源 \(\tilde T_{\mu\nu}/M_P^2\) 不能由当前低阶局部 \(\tilde g\) 曲率/张量特征稳定预测；
  - 这不排除高阶或非局域 pure-\(\tilde g\) 泛函，但排除了“简单局部函数”的乐观读法。

## 推进策略

- 当前策略：深度优先
- 优先事项：先完成分支 A 与分支 B，确保研究对象、目标判据与数学对象都被严格写清，再进入变换构造
- 工作节奏：每次只推进一个小闭环，优先完成“定义清楚一个对象”再做下一个推导
- 暂缓事项：暂不处理多体情形、测量问题、完整引力动力学
- 计划修正：不再只盯住混合 ansatz，而是把共形、相位单向、振幅单向、混合四类放入统一比较框架

## 当前活跃分支

- 分支：D1-D2-E2-E3
- 状态：当前主线已从“继续加密参考残差图”转为“先评估 C/D pure-`\tilde g` 作用量的理论健康性，再决定是否进入全作用量演化”。阶段性判断是：B 支在当前目标下基本排除；C 支保留为幂律饱和对照；D 支作为优先候选，但需先处理 `f_R -> 0`、过渡层导数项和 `f_{RR}<0` 稳定性风险。

## 2026-05-07 当前 D/gBCD 数值主线更新

- 当前 equation-first gBCD 主线已推进到 `plus-only` 初始加速度求解器：
  - 保持 \(\tilde g_0,\rho^\tilde,u_i\) 不变；
  - 由 D/gBCD 场方程求 \(g_+\)，即初始加速度；
  - 不再通过大幅修改 \(\rho^\tilde\) 或中心切片 \(g_0\) 来硬配 D 方程。
- `tau=-3.5,n=384,core10` 左侧分离态当前最佳可推进包：
  - `metric_active_dilation=0` 的 plus-only 重新求解；
  - `auto_bad_zero` 局部可容许性守卫；
  - D 方程 exact residual weighted mean `0.13220`；
  - 一步 \(g_+\) 测度偏差 `1.09e-4`；
  - 负质量壳判别式比例 `0`。
- 关键输出：
  - `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_noactiveedge_coeff_cg6000/`
  - `visualizations/equation_first_gbcd_plus_local_guard_scan_n384_taum3p5_core10_noactiveedge_cg6000/`
  - `visualizations/equation_first_gbcd_plus_initial_package_n384_taum3p5_core10_noactiveedge_guarded/`
  - `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_taum3p5_core10_noactiveedge_guarded/`
- 新研究笔记：
  - `research-notes/162-plus-only初始加速度包与可容许性守卫.md`。

## 下一步

- 每一轮继续研究前，先按长记忆协议重读 `README.md`、`TASKS.md`、`LOG.md`、`DECISIONS.md` 与最新 handoff，再进入当前理论/数值主线
- 当前主线下一步：不要把完整 \(\mathcal R_{\rm need}\) 的平凡 pure-geometry 拟合当作成功；应继续分析剥离 EH 主项后的非平凡剩余，尤其是 \(T^A-T^{(\tilde m)}\)、\(u,r,\rho\)-dependent gravitational sector，以及其 action-level metric/S/rho 变分是否可闭合。
- 当前数值下一步：先把 `no-active-edge + auto_bad_zero` 流程扩展到 `tau=0,+3.5`，确认干涉中心和右侧分离态也能形成 guarded package；随后把每一步 metric solve 的局部可容许性守卫升级为正式强约束/patch 边界，而不是只做后处理。
- 已完成 `tau=0,+3.5` 扩展：
  - `tau=0` residual weighted mean `6.60e-5`，冻结点 `0`，一步 \(g_+\) 加权测度偏差 `3.00e-5`；
  - `tau=+3.5` residual weighted mean `0.01324`，冻结支撑点 `2` 且不在 core10，一步 \(g_+\) 加权测度偏差 `2.88e-5`；
  - 三切片均无负质量壳判别式。
- 当前数值下一步更新为：
  - 以三切片 guarded package 为标准验证集，实现多步推进原型；
  - 每一步执行 plus-only metric solve、matter step、可容许性守卫、D residual 与测度偏差记录；
  - 将 `auto_bad_zero` 从后处理升级为求解器内部 active-set/强约束；
  - 继续单独优化 `tau=-3.5` 的高 residual 来源。

## 目标层级纠偏

- 重要纠偏：当前目标不是“完成一个 D 支数值计算器”。
- 真正目标是：
  - 在当前高斯波包干涉例子上寻找一个能描述变换后几何的、类似 Einstein 方程的显式 equation-first 场方程；
  - 该方程应尽量一般化，而不是只作为逐点数值拟合；
  - 它还要与物质部分的等价性和隐变量粒子的 \(\tilde g\)-测地线解释兼容。
- 数值工作的位置：
  - 高斯干涉数值结果是“方程形式筛选器”和“反例检查器”；
  - 初值包、guarded package、多步推进原型都只是检验候选方程闭合性和稳定性的工具；
  - 不得把数值器本身当作最终研究路线。
- 当前理论下一步应优先写清：
  - gBCD 方程的最终候选形式；
  - \(A,B,C,D\) 是局部本构函数、辅助场，还是由某个变分/最小化原则产生；
  - \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\) 如何保证测地线；
  - 三切片数值结果对上述方程形式提供了什么支持、暴露了什么缺口。
- 若继续检查 pure-\(\tilde g\) 希望，应从 action-level 必要条件入手：\(\tilde\nabla\)-守恒、Helmholtz/self-adjoint integrability、以及可能的高阶曲率/非局域泛函，而不是只做低阶曲率多项式拟合。
- 若继续 pure-\(\tilde g\) 局域路线，优先级应改为：quadratic curvature basis -> 导数不变量 -> action-level integrability；不要再把单变量 \(f(\tilde R)\) 的高次多项式当主线。
- 当前对“需要到多少阶才可能较好拟合”的经验判断已经更新：标准 4 阶 quadratic local action 仍几乎失败，最轻量 6 阶 pure-\(f(\tilde R)\) 复查也几乎没有改善，且最小真 6 阶 finite-jet 局域基底 \(\{1,\tilde R,I_2,\tilde R^2,\tilde R I_2,\tilde R^3\}\) 仍给出总体加权残差 \(\sim0.9987\)。因此下一步若继续 pure-\(\tilde g\) 局域路线，应扩到 \(I_3\) 与显式导数不变量，而不是继续在这些最小基底内打转。
- 现在可转入 equation-first 路线：先构造满足协变性、对称性、\(\tilde\nabla\)-守恒和正确平直极限的 pure-\(\tilde g\)+matter 闭合方程，再在后验检查它是否满足 Helmholtz/self-adjoint integrability、能否重建作用量。
- 当前显式方程候选已提升为投影型 Einstein-like 方程：
  \[
  \Pi_E^\perp\left(\tilde G_{\mu\nu}-\frac{1}{M_P^2}\tilde T_{\mu\nu}\right)=0,
  \]
  其中 \(\Pi_E\) 投影到 \(\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\}\) 张量子空间。
- 该方程必须与守恒约束
  \[
  \tilde\nabla^\mu\Pi_E\left(\tilde G_{\mu\nu}-\frac{1}{M_P^2}\tilde T_{\mu\nu}\right)=0
  \]
  以及 \(Q\to0\) 时投影残差消失的分支条件共同使用；否则它只是“残差方向限制”，还不是完整理论。
- 数值任务现在只能服务于这个方程候选：检查投影残差、守恒、近退化 patch 和 \(Q\to0\) 极限，不能把多步数值器本身当成主目标。
- 已完成投影方程闭合条件首轮分析：
  \[
  \det H=\frac{d-2}{2}\left(u^2r^2-(u\cdot r)^2\right)^3.
  \]
  因此 \(d>2\) 且 \(u,r\) 张成的二平面非退化时，四方向投影良定义；严格 \(1+1d\) 必须单独处理。
- \(W\) 的定位已明确：它不改变在壳成员关系，只影响离壳残差、数值权重和辅助场代表元规范。
- 已完成 3+1d 方程计数与主部首轮分析：
  \[
  \dim\mathrm{Sym}^2=10,\qquad \dim E=4,\qquad \dim E^\perp=6.
  \]
  这与 metric theory 去掉 4 个坐标规范自由度后的计数一致；主部为
  \[
  \delta\mathcal E^{\rm prin}_{\mu\nu}
  =
  \Pi_E^\perp\delta\tilde G^{\rm prin}_{\mu\nu}.
  \]
  harmonic gauge 下近似为 \(-\frac12\Pi_E^\perp\tilde\square\bar h_{\mu\nu}\)。
- 已完成 harmonic gauge 后主符号缺口检查：
  \[
  M_{\nu I}(\xi)=\xi^\mu E^I_{\mu\nu}
  \]
  泛型 rank 为 3，因此还剩一个 \(E\)-方向主部模式
  \[
  N^I=(0,-q^2,-p^2,2pq),\qquad
  p=\xi\cdot u,\quad q=\xi\cdot r.
  \]
  对应张量是
  \[
  \bar h^{(0)}_{\mu\nu}=-[(\xi\cdot r)u_\mu-(\xi\cdot u)r_\mu]^2.
  \]
- 当前明确需求：完整理论必须补一个协变标量状态方程，且其主部必须满足 \(s_I N^I(\xi)\neq0\)。首个候选闭合类是
  \[
  \mathsf q_E^{\mu\nu}(\Pi_E\mathcal R)_{\mu\nu}
  =
  \chi(\mathcal Q)\Theta.
  \]
- 已完成 \(\mathsf q_E\)-trace 最小硬闭合三切片检验：
  - 普通 trace=0 的中心 residual 加权均值为 `tau=-3.5:0.1185, tau=0:0.0264, tau=+3.5:0.1103`；
  - 最好的 \(\mathsf q_E\)-trace 版本约为 `0.1810, 0.2907, 0.1941`；
  - 因此 \(\mathsf q_E^{\mu\nu}\mathcal C_{\mu\nu}=0\) 或简单 \(Q,Q^2\) RHS 暂不作为主状态方程。
- 当前标量闭合优先级更新：普通 trace=0 是当前最小 algebraic closure；\(\mathsf q_E\)-trace 保留为主符号诊断/辅助场正定范数候选。
- 已检查普通 trace=0 的主符号退化风险：
  - 接近 \(w^2=0\) 的 trusted 点比例约为 `tau=-3.5:12%~14%`、`tau=0:2.7%`、`tau=+3.5:3.7%~5.5%`；
  - 这些点主要对应 \(r^\perp\) 退化/近退化 patch；
  - 因此 trace=0 可作为当前最小闭合，但必须配合 patch/branch 规则，不能声称全局无退化。
- 历史上先测试过的最小三项候选是
  \[
  \mathcal C_{\mu\nu}=Buu+Crr+D\,u_{(\mu}r_{\nu)}
  \]
  并把 \(h^\nu{}_\alpha \tilde\nabla^\mu\mathcal C_\mu{}^\alpha=0\) 作为测地线保持约束；\(A\tilde g_{\mu\nu}\) 暂列为扩展项，只有三项 constrained fit 不够时再加入。
- 最新 constrained-fit 结论：三项 BCD 的硬横向力约束在 trusted 区域会导致中心代数残差灾难性升高，因此当前最小严格候选已升级为
  \[
  \mathcal C_{\mu\nu}
  =
  A\tilde g_{\mu\nu}
  +Buu+Crr+D\,u_{(\mu}r_{\nu)} ,
  \]
  且必须同时满足
  \(h^\nu{}_\alpha\tilde\nabla^\mu\mathcal C_\mu{}^\alpha=0\)。
- 下一步优先任务已经从“逐点系数闭合”升级为“投影场方程闭合”：明确 \(\Pi_E\)、守恒约束、\(Q\to0\) 分支和 action/integrability 条件如何共同定义一个理论。
- 数值任务：`kg_examples/fit_equation_first_constrained_bcd.py --atoms g,uu,rr,ur --hard-constraint` 等脚本只作为投影方程的诊断工具，不能把三切片系数拟合当作最终方程。
- 已完成 `n=96` 三时刻 sanity check：
  - transverse hard 与 full divergence-free hard 均可通过；
  - full hard 的中心代数残差约 `1%~5%`，但相邻 probe 时间层残差约 `14%~23%`；
  - 这确认 gBCD 张量壳有希望，但时间闭合问题是真问题。
- 本构关系下一步：
  - 当前简单局部标量函数 \(A,B,C,D=A,B,C,D(\log\rho,\tilde R,I_2,u^2,r^2,u\cdot r)\) 初检失败；
  - 下一步不要直接开全动力学，而应先处理 \(A,B,C,D\) 的本构闭合：
    1. 固定 KKT null-space/gauge，例如最小范数、平滑范数或物理正则项；
    2. 加入导数型特征，如 \(\nabla u,\nabla r,\nabla R\)；
    3. 或把 \(A,B,C,D\) 作为辅助场，提出状态方程/演化方程。
- 最新机制判断：
  - ridge/minimum-norm 正则可改善 KKT 代表元，但不能把 \(A,B,C,D\) 全部变成简单局部函数；
  - 低阶单势函数 \(L(\rho,u^2,r^2,u\cdot r)\) 通过 metric variation 产生 gBCD 的机制已经测试失败；
  - 当前最优先路线是辅助各向异性应力场机制：
    \(\mathcal C_{\mu\nu}=\lambda_I E^I_{\mu\nu}\)，\(\lambda_I=(A,B,C,D)\)。
- 下一步具体任务：
  - 为辅助应力场机制选择一个代表元规范，优先比较：
    1. 最小范数；
    2. 最小空间梯度；
    3. 最小时间梯度；
    4. trace 或某个能量条件；
  - 将 full 守恒 PDE 与该代表元规范组合成真正的 \(\lambda_I\) 演化/约束系统；
  - 再检查这个辅助场系统是否能由作用量或 Helmholtz 条件支持。
- 以 `research-notes/072-三种引力作用量在2+1d下的全动力学方程与当前数值状态.md` 为新的引力作用量统一索引，后续所有 A/B/C 讨论都先区分：
  - 形式全方程
  - 显式闭合 PDE
  - 当前可信数值口径
- 梳理负密度、超光速问题在 KG/Bohm 图景中具体出现在哪里
- 对各方案的几何作用量直接做变分，检查是否能导出等价的 HJ 方程与连续性方程
- 对一般共形类 `g~_mu nu = A g_mu nu` 直接变分，求连续性方程与 HJ 方程条件
- 对一般类 II `g~_mu nu = A g_mu nu + B u_mu u_nu` 直接变分，求连续性方程与 HJ 方程条件
- 在后续所有变分中统一使用 `sqrt(-g)ρ = sqrt(-g~)ρ~`
- 撤回由“先代回再变分”得到的错误排除结论，并以“先变分再代回”重做筛选
- 对类 III `g~_mu nu = A g_mu nu + C r_mu r_nu` 做同样的结构性分析
- 若类 III 也失败，则主线集中到混合类
- 检查这些方案下 Bohm 轨迹是否真的能写成测地线方程
- 如用户给出不同推导路径，继续对照符号和方程来源，避免把“额外条件引入的 Z”误当成纯 HJ 结论
- 对最简单拉氏量消项解 `C=1, D=Q/X^2` 直接做后续方程与几何检验
- 在最简单消项解的基础上加入最少量的振幅结构修正，检查能否同时修复连续性方程
- 若继续保留最简单消项解，则把 `J^mu ∂_mu Q = 0` 视为候选特殊解族条件进行检查
- 区分“HJ/测地线成功”与“连续性方程成功”，不要把两者混为一谈
- 在 no-go 基础上决定下一轮是放宽连续性要求，还是升级 ansatz / 作用量
- 对纯 ansatz 的最简单解继续检查：全局可逆性、签名保持、因果结构、经典极限
- 研究 `X = m^2 + Q > 0` 的物理含义，特别是它与原“超光速/类空”问题的关系；若翻阅旧笔记，始终注意其中常用的是 `Q_old = -Q`
- 研究 `X < 0` 的 spacelike branch 是否能把原来的“超光速/类空”问题重新解释为另一条良好的 Lorentzian 几何分支
- 若后续回到 `m=0`，不要直接把 `m≠0` 主线里的 `κ=m^2`、密度变换律和连续性恢复公式机械取极限；应先重写 null-shell 物质作用量
- 当前混合 `u-r` 路线已出现一个漂亮的投影型解：`D=Z, E=-Y, F=X`；下一步应优先讨论其退化面 `X=0`、`Q=0`、`Δ=XZ-Y^2=0`，并区分哪些是纯 patch 边界，哪些是更深层的几何/作用量退化
- 已完成这一步的第一轮文献分层：当前默认先把 `X=0`、`Q=0`、`Δ=0` 当作 atlas / 非可逆面候选；下一步应把这一分层落成项目内的局域判据，而不是只停在文献层面的分类
- 已补做纯共形分支的统一复核：当前应把“纯共形在方程级可行，但在 pulled-back 作用量级退化”作为后续与 rank-1 / 混合 ansatz 对照的基准结论
- 已补做另一组混合解 `D=Z, E=-aY, F=aX, C=aXZ` 的检查：该解在 `XΔ≠0` 的 patch 上也可同时恢复 HJ 与连续性方程，但其退化结构更复杂，后续应与投影型主解并排比较，而不宜直接取代
- 已抽出一般 `C,D,E,F` 混合 ansatz 的行列式 strict positivity 总公式；后续讨论任何新解时，优先先套这一总公式，再谈因果性或经典极限
- 已澄清 `ρ~ -> 0` 在 `X=0` 处能做什么、不能做什么；后续若要真正 regularize `X=0`，必须明确要 regularize 的是“电流”还是“壳约束/度规图册”
- 已把这一点加强成严格局域结论：在一般混合 ansatz 中，只要在正则开区域里同时恢复连续性方程和固定非零 HJ 壳长，控制 `u` 方向的组合就被锁成 `κ/X`；后续若要真正避开这一主奇异，必须明确放弃哪一条要求
- 整理“固定 (R,S) 的严格匹配路线”与“允许 ρ~ 变换的表象路线”的统一对照框架
- 比较 `S_EH[g]` 与 `S_EH[g~]` 的一阶/二阶差异项
- 估计两种引力作用量在“弱量子势 + 近似平直背景”中的差异量级
- 识别只在“强引力 + 强量子势”场景下出现的潜在判别量
- 进一步识别即使背景引力弱、但量子势大或变化快时也会放大的判别量
- 设计可用于区分 `S_EH[g]` 与 `S_EH[g~]` 的实验/观测思路
- 把双高斯 KG 基线模型放入弱 Schwarzschild 背景做第一轮推广
- 把双高斯 KG 基线模型放入引力波背景做第二轮推广
- 分析这个双高斯模型中哪些自然定义的实验信号会对节点附近的大 `Q` 尖峰不敏感，哪些会敏感
- 用双高斯基线模型具体测试 `F(Q)=1/(1+(κQ)^2) R~` 的主支撑区偏离与节点区偏离
- 检查 `∇F`、`∇∇F` 是否比 `F` 本身更危险
- 把 `F(Q)=1/(1+(κQ)^2) R~` 放入弱 Schwarzschild 背景做第一轮曲背景检验
- 在 1+1d 双高斯模型上继续沿 `F(Q)R~` 路线推进，获得真正的局域预测差异
- 为 A 分支选定一个明确的 `2+1d` reduced closure（不再含混地把 `g` 与 `g~` 的关系留空）
- 先把 `2+1d` 数值器的初值几何修回用户认可的 `z=±x` 两束相交束流，再重跑 B/C 与参考支
- 在同一初值下补做 A 分支的全动力学数值演化，使三支比较真正闭合
- 对 B/C 数值器做稳定性扫描：`M_P`、`ℓ`、网格、时间步长、边界海绵层
- 检查 `F(R~)R~` 分支出现强回授时，是物理结构还是数值 stiffness 主导
- 为 A 分支选定一个明确的 `2+1d` reduced closure，使 `sqrt(-g)R[g]` 这条支也能和 B/C 一样落到唯一的显式闭合 PDE
- 在 `A` 支的 `g` 上 Einstein-KG 全动力学写成一般 `2+1 ADM` 约束—演化系统后，再和 `B/C` 做真正的同初值全动力学比较
- 对新混合变换，shared-baseline 的下一步不再是直接重跑旧脚本，而是先确定 `X<0` 区域的 source prescription；在这一步明确前，不应机械地把旧 `rho_A` 替成 `(X/kappa) rho_A`，更不能把它误当成正测度密度；当前正密度口径应先写成 `sqrt(|g~|)ρ~ = (|X|/|κ|)ρ_A`
- 当前已经确认：对 `±45°` 高斯波包初值，当前混合变换会生成 `t-x`、`t-z`、`x-z` 非对角分量；因此真正的下一步不是继续修补旧对角 gauge 数值器，而是把三支全动力学统一搬到一般 `2+1 ADM` 度规下
- 新方法的第一步已经完成：共同物质部分已在一般 `ADM` 度规下重写成 `(n,S)` 的守恒系统；下一步是分别重推 B/C 的 `ADM` 引力演化方程，并为 A 支写出一般 `ADM` 下的 `g↔g~` 约束闭合
- 新方法的第二步已经完成：`B/C` 两支在一般 `2+1 ADM` 度规下的作用量和协变场方程已重写完成
- 当前最直接的下一步已修正为：把 `A` 支按“全部写在 `g` 上”的 Einstein-KG 系统重写成一般 `2+1 ADM` 约束—演化方程；`g~` 只保留为事后重构的有效几何
- 新方法已进一步统一到四维 `ADM` + `y` 方向 Killing 对称；当前不再把主要精力放在 `2+1` 共形重参数化上
- `A/B` 两支的最小同步规范原型已经启动并可短时间推进；当前最直接的下一步是：
  1. 以 `\Phi=1,\Pi_\Phi=0` 为校准点，校正 `C` 支的 `ADM` 有效源项与平直 `\Phi` 扣除
  2. 在此基础上为 `C` 支构造非平凡而一致的 `\Phi` 初值
  3. 为三支求解满足约束的初始几何，再做真正的同初值比较
  4. 显式计算共同初值下三支的动量约束残差，判断当前“共同 `g` 初值”近似是否已经够用
- 在当前认可的 `±45°` 初值几何下导出更稳定的单束/双束 sanity check 套件
- 围绕 `B -> A` 的弱耦合回归做最小测试，定位 B 分支离散实现中的 bug
- 在确认 B 分支能弱耦合退回 A 之后，再恢复 C 分支和更强耦合参数
- 重新定义“三支相同实验准备”的边界表述，转向共享 incoming characteristic data
- 物质边界的特征入流/出流条件已经完成；下一步不再继续修补物质边界，而是为 `C` 支构造真正的椭圆型 `Phi` 初值方程
- 在当前长时间稳定版本上，显式写出 `Phi` 初值的椭圆方程、离散形式与边界条件
- 原始完整椭圆型 `Phi` 初值已经探测过，结果过强并把 `C` 支立即打爆；下一步应改为“椭圆型初值的连续延拓/幅度控制”而不是整幅代入
- 以受控参数 `\lambda` 构造
  `\Phi_\lambda = 1 - \lambda(1-\Phi_{\text{elliptic}})`
  或等价的更温和辅助变量，扫描 `\lambda` 的可稳定区间与 `C-B` 分离量级
- 在引入椭圆型 `Phi` 初值后，再评估是否有必要把几何子系统从同步规范 `ADM` 切换到更严格的初边值问题形式
- 从 dense scan 中挑选 `lambda = 3e-2, 3e-1, 3, 10` 作为代表点，开始正式做 A/B/C 图样比较与物理解读
- 对 C/D 支在进入全作用量演化前新增理论检查：
  - 写清 scalar-tensor/`chi` Cauchy 问题和 principal part；
  - 对过渡层做积分形式 jump/boundary condition 分析；
  - 估计线性响应算子条件数，避免把“小残差”误读成“小演化偏离”；
  - 标定 `ell` 的物理尺度窗口 `|R_ordinary| << ell^-2 << |R_tilde_quantum|`。
- 新增近期任务：实现 `ell-window scan`，至少比较 `ell = 1, 3, 10, 30, 100, 300`，输出 `f_R` 退化区、`f_RR<0` 稳定性红旗区、scalaron 质量平方符号分布、过渡层宽度、热点重合率和残差统计。
- 已完成首轮 `ell-window scan`。下一步在进入全作用量演化前，应优先：
  - 对 `ell=30,100,300` 的 D 支做局部窗口高分辨率过渡层积分，估计 `[\partial_n f_R]^+_-`；
  - 对正 `R~` 区估计 scalaron 增长时间 `1/sqrt(|m_s^2|)`，并区分 finite 区与 `f_RR` 数值下溢/undefined 饱和区；
  - 把 C/D 写成 EH-like 区、saturated algebraic 区、transition layer 三分区匹配问题；
  - 然后再选择一个中间 `ell` 进入全作用量演化。
- 追问后修正的近期任务：
  - 对 `ell-window scan` 增补 `M_P^2 f_R`、`M_P^2 f`、`T_munu` 的同口径尺度比较；
  - 实现过渡层 jump/层强度诊断：`[\partial_n f_R]`、`\int \partial_n^2 f_R dn`、`\int |\partial_n^2 f_R|dn`；
  - 写出饱和区固定度规的候选选解原则，并判断它是否可由 `M_P -> infinity`、渐近平直边界条件或最小曲率原则推出。
- 已完成第一版 `M_P^2` 尺度比较和局部 jump proxy。下一步需要从 proxy 升级为真正法向短线积分：
  - 在左弧窗口上沿 `grad(ell^2 R_tilde)` 构造穿层短线；
  - 直接积分 `[\partial_n f_R]`、`\int \partial_n^2 f_R dn`、`\int |\partial_n^2 f_R|dn`；
  - 对 `ell=30,100` 做局部加密，确认严格层是否能被解析；
  - 同时写出 saturated bulk 固定度规的两个候选原则：最小曲率延拓与 EH-like 边界参考延拓。
- 路线升级：将过渡层作为 interface/jump 条件处理。
  - 下一步脚本目标：在 D 支 `ell=30,100` 的左弧窗口上提取 `y=ell^2 R_tilde` 的 `|y|=1` 水平集；
  - 沿水平集法向构造短线采样，计算 jump 与绝对层强度；
  - 做 `160/240/320` 分辨率对照，判断界面积分是否收敛；
  - 若收敛，开始把 full evolution 设计为 EH-like bulk、saturated bulk 与 interface matching 三块耦合。
- 已完成 D 支动态界面短线积分首轮。下一步：
  - 以 `ell=30` 为首个 interface toy evolution 参数，因为其界面积分已有初步收敛迹象；
  - 写出 bulk-interface toy system：EH-like bulk、saturated fixed/extended bulk、`y=±1` interface jump 源；
  - 对 `ell=100` 暂不硬推全演化，改做局部贴体坐标或更高分辨率界面采样。
- 已完成 `ell=30` bulk-interface toy system 首版。下一步：
  - 已修正首版脚本的左弧窗口限制，新增全域界面提取 `--bbox global`；
  - 已测试全域单界面 `|y|=1` 与全域多界面 `|y|=0.5,1,2`；
  - 结果显示全域 `|y|=1` 已覆盖 top 1% derivative hotspot 的 `100%`，多界面没有提高覆盖率，只增加层内源项采样；
  - 已写并跑通第一版 subcell/body-fitted interface toy evolution；
  - 第一版结论是：子网格界面线源能压低最极端尖峰，但单独 `δχ` 校正会抬高 p95；
  - 已按物理锁定口径检查 D 支迹方程空间法向弱形式，发现固定时间切片界面不够；
  - 下一步不再继续推进 screened-Poisson toy correction，而是推导时空动态界面 \(F(t,x,z)=ell^2 R_tilde\mp1=0\) 的 jump/matching 条件，再基于该弱形式设计数值求解器。
- 已完成 D 支 leading trace/scalaron 时空界面速度律首轮：
  - 已把 jump 条件改写成 `jump_law_residual_scaled=0`，可反求目标 `tilde g^{ab}F_aF_b`；
  - 已新增 `kg_examples/diagnose_d_interface_speed_law.py`，能在可解界面点上反求满足 jump law 的 `F_t` 与坐标法向速度；
  - `ell=30,t=16,128x128` 下可解比例约 `0.624`，修正速度 p95 约 `1.44`，修正后 scaled residual p95 约 `2.17e-19`；
  - 下一步应将该速度律升级成 moving-interface reduced simulator：界面点按 jump-law 速度推进，bulk 由 EH-like / saturated 区的合法弱形式闭合；
  - 对无实根点，不做数值削峰，应升级到完整 tensor jump/matching 或允许界面空间形状共同调整。
- 已完成连续线段界面图和第一版 moving-interface reduced prototype：
  - `diagnose_d_interface_speed_law.py` 已改为连续线段图，默认可用于分析界面物理结构；
  - 新脚本 `prototype_d_moving_interface_reduced.py` 每步重新抽取界面，形式上允许界面生成/消失/合并/断裂；
  - 首轮结果显示直接推进原始 `raw_y` 会暴露强 stiffness，因此下一步完整模拟器应改成 signed-distance / body-fitted interface 表示，而不是 raw_y 显式场推进；
  - 下一步：实现 signed-distance 界面重初始化版本，并开始接入 tensor jump/matching；未解段仍不得 clipping/damping。
- 已完成 signed-distance moving-interface 原型：
  - 新脚本 `prototype_d_signed_distance_moving_interface.py` 使用 `d_plus/d_minus` 表示 `raw_y=±1` 两个动态界面；
  - 每步用 D 支 leading trace/scalaron jump law 计算法向速度，并用 signed-distance 重初始化避免直接推进 `raw_y` 的 stiffness；
  - `ell=30,t=16,n=96,dt=0.01,steps=20` 下，界面段数 `864 -> 712`，总线长 `307.053 -> 214.710`，可解比例约 `0.60-0.63`；
  - 下一步应补完整 tensor jump/matching，使 unresolved 段有物理闭合；同时把冻结的 `rho/metric_inv/stress_trace/alpha` 改为由当前 \((\rho,S,\tilde g)\) 自洽更新。
- 已完成 D 支 tensor jump 首轮诊断：
  - 新脚本 `diagnose_d_tensor_interface_jump.py` 检查张量薄层主部；
  - `ell=30,t=16,n=128` 下 trace 条件稳定，但固定参考界面法向给出巨大无迹 tensor residual；
  - free-covector rank-one mismatch median `0.018`，说明很多段存在近似可行的目标 covector；
  - 目标自由法向与当前界面空间法向 alignment median `0.486`，说明下一步应调整界面形状/法向，而不是只调速度。
  - 下一步：实现 local tensor-matching interface solver，由 algebraic tensor integral 反推出目标 \(F_a\)，并把它接入 signed-distance moving-interface 原型。
- 已按用户纠正完成 direct tensor interface 条件首轮：
  - 新脚本 `diagnose_d_direct_tensor_interface.py` 不再先解 trace，而是直接由 \(H_{ab}=-F_aF_b+\tilde g_{ab}F^2\) 反推 \(F_a\)；
  - `ell=30,t=16,n=128` 下，`10%` 数值验收 accepted `456/1361`，`1%` 数值验收 accepted `62/1361`；
  - accepted 段速度 p95 接近 `1`，但目标/当前空间法向 alignment median 约 `0.41`；
  - 下一步应实现 tensor-driven interface reconstruction：用目标 \(F_x,F_z\) 重构 signed-distance/贴体界面形状，而不是沿当前 `raw_y` 法向做 trace-speed 推进。
- 已完成 direct tensor 目标法向/切线重构预览：
  - 新脚本 `plot_d_direct_tensor_normal_reconstruction.py` 输出 accepted/rejected 界面段、目标法向红线、目标切线蓝线、完整张量残差、alignment、target speed；
  - `accepted` 不使用 trace speed，只用 direct tensor rank-one + full tensor residual；
  - 输出目录：
    - `visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/`
    - `visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/`
  - 下一步：把目标 \(F_x,F_z\) 场接入真正的 level-set/body-fitted 重构器；未通过 direct tensor 的段必须保持 `tensor-unresolved`，不能退回 trace speed。
- 已完成 tensor-driven 界面重构原型：
  - 新脚本 `prototype_d_tensor_reconstructed_interface.py`；
  - 只使用 direct tensor accepted 段构造目标切线场，且 `raw_y=+1/-1` 分支分开；
  - `10%` 验收：accepted `456/1361`，重构曲线 `29` 条，总长度 `165.86`；
  - `1%` 验收：accepted `62/1361`，重构曲线 `16` 条，总长度 `44.15`；
  - 结论：accepted 区域可以整理成曲线，但 accepted 区域太少，不能构成完整全界面。
- 已完成 direct tensor rejected 模式分类：
  - 新脚本 `diagnose_d_direct_tensor_rejection_modes.py`；
  - `10%` 下：`tensor_residual` 约 `27.8%`，`rank_tail+negative+tensor_residual` 约 `33.1%`，`no_real_covector` 约 `3.23%`；
  - 下一步：
    - 对 `tensor_residual` 型段加入完整层内 profile、切向/外曲率项、两侧 bulk 几何匹配；
    - 对 `rank_tail+negative+tensor_residual` 型段做分辨率、层宽采样、ell 扫描，判断是数值假象还是 D 支局部不可闭合信号。
- 已完成 direct tensor 失败原因复查：
  - 新脚本 `diagnose_d_direct_tensor_failure_features.py` 显示失败强烈集中在低密度、高 `|grad raw_y|`、局部 \(\tilde g\) 条件数高的区域；
  - 分辨率提高会改善 accepted：`n=96,128,160` 下 `0.254 -> 0.335 -> 0.398`，但顽固失败类型不消失；
  - `ell=10` 与 `ell=30` 的失败比例基本相同；
  - `half-width=0.5` 太窄会坏，`half-width=1` 与 `2` 基本相同；
  - 新脚本 `diagnose_d_direct_tensor_least_squares.py` 显示原本 `tensor_residual` 型大多可由直接最小化完整张量残差修复；
  - `10%` 下最小二乘 accepted `867/1361=0.637`，`1%` 下 `333/1361=0.2447`。
- 已完成 least-squares 界面重构与分辨率复查：
  - `prototype_d_tensor_reconstructed_interface.py` 已支持 `--solver least_squares`；
  - `ell=30,t=16,n=128` 下 accepted `867/1361=0.6370`，重构曲线 `30` 条；
  - `ell=30,t=16,n=160` 下 accepted `1214/1847=0.6573`，重构曲线 `40` 条；
  - accepted residual p95 从 `0.0736` 降到 `0.0693`，nearest curve distance p95 从 `0.3385` 降到 `0.2628`；
  - 判断：分辨率需要继续提高用于界面几何收敛，但主瓶颈已从普通解析度转为完整 tensor matching 的局部求解和顽固失败段的物理闭合。
- 已完成局部分辨率收敛与 multistart 复查：
  - 新脚本 `diagnose_d_ls_local_resolution_convergence.py` 已自动选出 `W1-W4` 四个失败窗口；
  - `n=128,160,192` 全局 accepted fraction 为 `0.6370 -> 0.6573 -> 0.6458`，不是单调收敛；
  - `W1-W4` 中只有部分窗口改善，未解 residual p95 仍约 `1`；
  - `--multistart --max-seeds 12` 在 `n=160` 下只把 accepted 从 `1214` 提高到 `1242`，额外修复 `28` 段；
  - 判断：剩余失败不能主要归因于普通分辨率不足或单初值 Gauss-Newton 卡住。
- 已完成低密度/计算精度误差检查：
  - `rho_floor=1e-10,1e-12,1e-14` 完全不改变 `n=160` least-squares 结果；
  - `probe_dt=2.5e-4,5e-4,1e-3` 只带来小变化，未解 residual p95 仍约 `1.04`；
  - `samples=41,81,121,161` 基本不改变结果；
  - 高密度区仍有未解段，`rho/rho_max>=0.1` 时 unresolved fraction 约 `0.233`；
  - 判断：普通精度误差不是主因，但低密度区后续必须单独标记为不可靠支持区。
- 已按用户纠正修正 residual 诊断口径：
  - A 参考解上的 jump residual 不是 D 支最终失败率；
  - 它只衡量 A 分支平直量子力学参考解离 D 支约束/匹配流形有多远；
  - D 支真实演化中 \(\rho,S,\tilde g\) 和界面都应动态改变，并由 jump law / bulk 方程共同决定；
  - 若 A 参考解完全满足 D 支所有方程，才会说明两支在该问题上接近动力学等价；不完全满足是预期的 A/D 差异入口。
- 当前下一步：
  - 把 residual 诊断转化为 D 支一致初值构造问题：寻找最小 \(\delta\rho,\delta S,\delta\tilde g,\delta\)interface，使 D 支约束和 jump law 在初始切片上成立或尽量小；
  - 用 least-squares \(F_a\) 作为 direct tensor interface solver 的默认局部闭合；
  - 统计仍需分成 `support-trusted` 与 `low-density-untrusted`，不要用低密度失败段直接否定 D 支；
  - 在动态界面原型中让 \(\rho,S\) 也按 D 支物质方程改变，而不是继续冻结为 A 参考解；
  - 用 level-set/body-fitted 表示界面，并由 jump law 决定界面法向速度及生成/消失。
## 2026-04-26 Current Focus

- Keep A fixed as the trusted benchmark solver.
- Treat B matter transport as the shared core module for B/C and continue with the conservative Rusanov-style update.
- Debug C only through the `phi` dynamics now that `C(phi=1)` is back on the B scale.
- Latest isolation result: `phi` alone stays small, but `phi -> (tau,beta)` coupling still blows up even on fixed matter.
- Flat-limit consistency for `C(phi=1)` has now been repaired.
- Next concrete target: replace or reformulate the `phi` auxiliary variable itself, because current evidence says the remaining stiffness comes from the `phi` representation, not from the matter source.
- Single-beam B check is now essentially done:
  - `Bflat ≈ Bfull`
  - dt-halving changes little
  - both remain stably different from A
- Next comparative task: decide whether to accept this as a real physical split between A and B, or to go back one level and re-derive the flat-limit B matter equations from the action to confirm the implemented PDE is the intended one.
- C comparison result after the recent fixes:
  - shrinking `ell` still does not bring C closer to A
  - current priority is no longer line-width scanning but rethinking the C auxiliary-variable formulation itself
- Updated priority after user clarification:
  - treat “gravity switched off => A/B/C coincide” as the intended physical requirement unless re-derived otherwise
  - stop trusting the current flat-limit B/C matter PDE as the final target
  - next task is to re-derive the gravity-off / flat-limit equations from the intended action-level comparison before further solver tuning
- Immediate next task:
  - reintroduce the geometrization constraints into the gravity-off B/C limit
  - derive the reduced 2+1d solver target from `S_geom,matter + S_constraint`, not from a free flat HJ system
- Gravity-off target is now clarified:
  - A should be solved dynamically
  - B/C should first be reconstructed algebraically from A in the no-gravity limit
  - only after that should genuine B/C geometry dynamics be turned on
- Next implementation task:
  - refactor the comparison scripts so the gravity-off benchmark for B/C is constraint reconstruction from A, not `Bflat/Cflat` evolution
- New concrete next step:
  - use `kg_examples/simulate_bc_gravity_off_target.py` as the regression target
  - refactor the old B/C flat-limit solver so its `gravity=0` path converges to this target, not to the old free-HJ surrogate
- First working continuation/debug solver now exists:
  - `kg_examples/simulate_bc_small_gravity_backbone.py`
  - zero gravity and weak gravity limits are continuous and regression-clean
- Next task:
  - increase `lambda_grav` until visible but still controlled deviation appears
  - then decide how to replace the passive-geometry approximation with genuine B/C backreaction
- Current quick regression harness:
  - `kg_examples/simulate_abc_zero_gravity_regression.py` now produces the expected identical A/B/C `rho` maps at `gravity=0`
  - next real step is to reintroduce small gravity continuously on top of this shared backbone, instead of using the old incorrect flat-limit B/C solver
- Continuation solver status:
  - `B` looks numerically healthy on the shared backbone up to at least `lambda_grav ~ 10`
  - `C` still needs solver work; current harness shows dt-sensitivity around `lambda_grav ~ 1` and numerical re-pinning of `phi -> 1` by `lambda_grav ~ 10`
- Next concrete debug:
  - stop trying to rescue the current `phi=f_R` propagation variable
  - re-formulate C in a signed curvature-like auxiliary variable (`chi` or equivalent)
  - then repeat the continuation scan on that reformulated C solver
## 2026-04-26 C Solver Update

- `C` flat-limit source consistency has been repaired.
- Current status:
  - `lambda_grav=1e-3`: C still tracks A/B very closely
  - `lambda_grav=1.0`: C now separates from B, and dt-halving shows this regime is converging
  - `lambda_grav=10`: C improves strongly as dt shrinks and still appears to be converging downward
- Next immediate task:
  - use C-geometry subcycling as the default strong-gravity stabilization method
  - with the solver now stabilized, start the actual A/B/C comparative study on top of this unified baseline

## 2026-04-27 Current Focus

- Shared-baseline dense continuation scan is now available up to `lambda = 10`.
- B is now quantitatively confirmed to separate from A in a smooth near-linear way.
- C can be compared on the same baseline, but the bridge region must still be read together with the chosen subcycling schedule.
- Same-initial dynamic pattern comparison is now implemented in:
  - `kg_examples/compare_abc_same_initial_evolution.py`
- Completed representative cases:
  - `lambda = 1`, `c_substeps = 4`
  - `lambda = 3`, `c_substeps = 8`
- `lambda = 1` 的 A/B 实验量提取也已完成：
  - `kg_examples/analyze_ab_measurable_observables.py`
  - 当前结果表明 A/B 在 `lambda=1` 下几乎不出现条纹间距或相位滑移，最可测的入口是 `~1e-4` 的亮纹强度重标定
- `lambda = 1` 的 A/B observable 三层补完已完成：
  - `kg_examples/scan_ab_observable_convergence.py`
  - `kg_examples/analyze_ab_time_observables.py`
- 当前理论整理的最新优先级：
  - 先统一“类似 Bekenstein 的成功变换”应满足的三条主张：
    1. 保持因果序
    2. 原 `g` 表象下的 HJ 与连续性方程不变
    3. `g~` 表象下给出测地线解释
  - 然后再把 `X=0`、`Q=0`、`Δ=0` 作为边界/图册问题单独检查
- 当前最新交付：
  - 面向这条主线的 7 页短 deck 已重写并导出：
    `presentation-workspaces/x0-boundary-summary/output/output.pptx`
- 下一步：
  - 以这版 deck 的主线为准，继续检查 `X=0` 处的联络、曲率和弱延拓
  - `kg_examples/compare_ab_background_scenarios.py`
  - refined 后的主 observable 为 `central_peak_rel_intensity_shift ~ 6.76e-05`
  - 其数值地板约 `3.57e-09`，信号高出地板约 `1.9e4`
  - 中心圆盘、亮区总强度、干涉盒总强度等积分型 readout 在最终时刻都稳定落在 `6.6e-05 ~ 6.8e-05`
  - phase / fringe-shift 类型 observable 仍不够干净，不宜作为当前实验主入口
- Earth Schwarzschild 与纳赫兹 GW 背景的第一轮弱场取向估计也已完成：
  - 这两类真实弱背景都不会把当前 `lambda=1` 的 A/B 差异放大成显著的条纹位置或相位信号
- 背景层也已升级成显式场景梯度与取向扫描：
  - Earth surface / LEO / GEO / Sun-Earth L1 / 近源质量 / 纳赫兹 GW
  - 当前这些都仍是 vacuum tidal / readout modulation，而不是 full curved-background 自洽动力学
  - 以 refined 主 observable 为基线时，这些背景诱导改变量仍远低于当前数值地板
- 两类实验双设置差分量也已完成：
  - `kg_examples/analyze_ab_differential_experiment_observables.py`
  - `phase-flip differential` 的 A/B 分裂稳定在 `~1e-8`
  - `interference excess` 的 A/B 分裂可稳定到 `~1.3e-06`
  - 当前更优的差分通道是：
  - `interference excess`
  - `interference box`
- 实验通道噪声基准也已完成：
  - `kg_examples/analyze_ab_noise_baselines.py`
  - `raw intensity` 通道的 `5 sigma` 计数需求约 `5e9`
  - 最优差分通道 `interference excess` 的 `5 sigma` 总计数需求约 `1.8e13`
  - `phase flip` 的统计和系统误差代价都过高，当前不适合主打
  - 现阶段最合理的排序是：
  - `raw integrated intensity` first
  - `interference excess` as companion cross-check
- 用户已明确当前阶段策略：
  - 暂时挂起 `C` 与更一般作用量
  - 主线只讨论 `A/B`
  - 当前阶段的成功标准不是“全面证明理论安全”，而是找到至少一个可测的、非数值误差的 A/B observable，并把它和噪声基准比较
- Current controlled conclusion:
  - B departs from A later and more gently
  - C departs from B earlier and more strongly
  - the main separation builds up during beam overlap rather than being present at `t=0`
- Immediate next task:
  - stop prioritizing `lambda=3,10` as if they were direct real-world predictions; keep them only as sensitivity/debug references if needed
  - treat intensity-type observables as the canonical AB experimental channel:
  - `central_peak_rel_intensity_shift`
  - `central_disk_rel_shift`
  - `bright_region_rel_shift`
  - keep the new differential channels as the experiment-facing companion set:
  - `phase-flip differential`
  - `interference excess`
  - start the next phase of experimental discussion around explicit noise/systematics baselines rather than around numerical existence
  - this baseline comparison is now done; next phase should move from generic noise models to a more explicit experimental architecture if needed
  - compare:
  - Earth surface vs Sun-Earth L1
  - GW background vs flat
  - near-source-mass orientation/null-test scenarios
  - and also compare:
  - strongest raw signal scale `~6.7e-05`
  - cleanest differential-channel scale `~1.3e-06`
  - and now compare against:
  - shot-noise count budgets
  - gain/source drift budgets
  - single-beam normalization budgets
  - keep a clear boundary:
  - current background numbers are scenario ladders on top of the refined flat AB split
  - full curved-background self-consistent PDE evolution remains a later task, not a finished deliverable

## 2026-04-27 Presentation Deliverable

- 已完成一版完整技术汇报 deck：
  - `presentation-workspaces/quantum-potential-geometrization-ab/output/output.pptx`
- 当前 deck 任务状态：
  - `completed`：主线说明版 PPT（14 页，含数值模拟配图、A/B 结论、差分实验量、噪声基准、appendix）
  - `completed`：逐页 PNG 预览与 package QA
  - `completed`：图片资产稳定内嵌方案
- 当前建议的下一步分支：
  - `branch A`：把这份 deck 精简成论文答辩/组会版（更短、更强调 claim）
  - `branch B`：转写成论文初稿结构，准备 `.tex` / Overleaf 版本
  - `branch C`：继续扩实验架构细节，把 raw intensity 通道翻成更具体的 measurement protocol
- 当前 active branch 建议：
  - 优先 `branch B`
  - 原因：用户已经明确认为理论文章“在没有出错的前提下已经差不多可以发表”，而 deck 已把故事顺序整理好，最自然的后续就是论文初稿落地
- Ordering strategy：
  - `opportunistic`
  - 先用 deck 固化叙事，再根据用户选择切到 paper 或精简版 slides
- Next immediate actions：
  - 等用户指定是要：
    - 论文初稿结构 / `.tex`
    - deck 精简版
    - 继续实验 protocol
  - 在用户审阅当前 deck 后，优先修正方法口径与展示规范：
    - 明确区分“共享主干比较”与“全反作用自洽演化”
    - 全文统一改成中文主叙述，必要术语用“中文（English）”
    - 所有图和数值量第一次出现时都补上严格定义
    - 禁止继续使用未定义对象的 `(B-A)/A` 一类缩写

## 2026-04-28 Boundary-Interpretation Deck

- 已完成一版新的局部总结 deck：
  - `presentation-workspaces/x0-boundary-summary/output/output.pptx`
- 这版 deck 的目标：
  - 不讲全项目
  - 只统一 `X=0`、`ρ~`、`sqrt(|g~|)ρ~`、电流与速度发散之间的解释关系
- 当前状态：
  - `completed`：7 页小 deck
  - `completed`：逐页 PNG 预览
  - `completed`：package QA / output hygiene 检查
- 当前结论口径：
  - 应优先以 `sqrt(|g~|)ρ~` 作为 `g~` 表象里的密度解释对象
  - 在强连续性匹配且 `A_u=κ/X` 的分支里，统一有 `sqrt(|g~|)ρ~ -> 0`
  - 因而 `X=0` 可解释成 `g~` 表象中的零测度密度边界，速度样对象无需在该点有良好定义
- 下一步建议：
  - 用户先审阅这版小 deck，统一语言
  - 之后再继续检查 `X=0` 处联络、曲率与弱延拓问题

## 2026-04-29 Realistic Weak-Field Numerical Strategy

- 当前新的方法口径：
  - `A` 支先不作为主要全动力学数值对象
  - 先在平直背景上推进 `A` 支复 Klein-Gordon 场，生成初始 `(\rho,S)`
  - 主要全动力学数值工作放在 `B/C` 两支自身的完整耦合演化上
  - 最后用线性化 `A` 支弱场估计给出“把 A 当固定背景”的误差量级
- 已完成：
  - 新脚本 `kg_examples/simulate_bc_geometry_from_a_reference.py`
  - 参考 `A` 支已改成真正的平直背景复标量场演化，不再使用旧的 `(\rho,S)` 近似推进器
  - 第一轮结果已整理到 `research-notes/082-A支平直参考驱动BC几何的第一轮结果.md`
  - 更正后的完整耦合脚本 `kg_examples/simulate_bc_from_a_initial_data.py`
  - 对应记录 `research-notes/083-B与C两支从A支初始数据出发的完整耦合同步规范测试.md`
- 当前数值结论：
  - `A` 支弱场误差量级约为 `4e-3`
  - `B/C` 支的几何偏离已达 `1e-1` 量级
  - 把时间推进从二阶升到四阶后，时间步敏感性几乎不变
  - 把 `A` 的线性化弱场 `omega` 用作 `B/C` 初始几何后，结果也几乎不变
  - 即使改回“`A` 只给初始数据、`B/C` 完整耦合推进”的正确口径，这种增长仍然存在
  - 当前 `C` 仍取 `\Phi=1,\Pi_\Phi=0`，因此与 `B` 严格重合
- 下一步：
  - 优先检查同步规范本身是否导致了主要约束增长
  - 优先检查 `B` 支几何演化的约束传播与阻尼写法
  - 只有在 `B` 支达到可接受收敛后，再引入 `C` 支非平凡辅助场初值

## 2026-04-29 Boundary Fix After Weak-Field Strategy

- 已完成：
  - `B/C` 新原型中初始动量约束外曲率的补入
  - 动量约束残差监测
  - 几何边界条件修正：
    - 物质变量 `n,S` 做阻尼
    - 几何变量回到背景值而不是乘法压向 `0`
  - `C` 支 `r3_proxy` 初值在修正后原型上的长时测试
- 当前最新结论：
  - 先前的巨大 `B/C` 几何增长主要是边界条件伪影，不是时间推进阶数问题
  - 修正后 `B` 支在 `200` 步、`dt=2.5e-4` 下仍保持弱场稳定：
    - `max|h_xx-1| ≈ 7.74e-3`
    - `max|H| ≈ 2.12e-1`
    - `max|M| ≈ 4.82e-2`
  - `A` 支线性化弱场量级 `~3.9e-3` 与当前 `B` 支几何偏离同阶
  - 当前 `C` 支即使取 `r3_proxy` 初值，在 `ell=0.02` 下仍几乎和 `B` 重合：
    - `max|Phi-1| ≈ 2.6e-8`
- 当前 active branch：
  - `B/C` 弱场稳定演化与 `C` 支非平凡分裂
- Ordering strategy：
  - `depth-first`
  - 先固定当前稳定的 `B/C` 原型，再引入更强的 `C` 支初始 `Phi`
- Next immediate actions：
  - 延长 `B` 支的物理时间窗口，确认长时稳定性
  - 设计比 `r3_proxy` 更强、但仍自洽的 `C` 支初始 `Phi`
  - 用线性化 `A` 支解量化“把 `A` 当平直背景参考”的系统误差

## 2026-04-29 `C` 支 `ell` 扫描后的当前任务

- 已完成：
  - 在统一稳定版本上，为 `C` 支加入显式 `B-C` 终态差异摘要
  - 在 `ell = 0.02, 0.05, 0.1, 0.2, 0.5, 1.0` 上完成第一轮扫描
  - 形成结论：当前弱场双束问题里，默认 `ell=0.02` 时 `C` 数值上几乎退回 `B`
  - 形成结论：纯局域代数 `Phi` 初值不能满足当前迹方程的符号需求
- 当前活跃分支：
  - `B/C` 稳定弱场演化
  - `C` 支非平凡 `Phi` 椭圆型初值
- 下一步：
  - 以 `trace(T_B)` 和当前 `B` 支几何为输入，构造一个真正的椭圆型 `Phi` 初值方程
  - 比较该椭圆型 `Phi` 初值与当前 `r3_proxy` 的偏离量级和 `C-B` 分离量级
  - 只有在这一步后，才决定继续扩大 `ell`，还是固定 `ell` 转向更物理的 `Phi` 初值

## 2026-04-29 椭圆型 `Phi` 初值延拓扫描后的当前任务

- 已完成：
  - `elliptic_lambda = 1e-4, 3e-4, 1e-3, 3e-3`
    的小窗口扫描
- 当前结论：
  - 首选基准：`elliptic_lambda = 1e-4`
  - 激进对照：`elliptic_lambda = 3e-4`
  - `1e-3` 及以上当前约束代价已偏大
- 下一步：
  - 已完成：以 `elliptic_lambda = 1e-4` 和 `3e-4` 做 `steps=200` 的更长时间 `B/C` 比较
  - 当前新的默认顺序：
    - 默认基准：`elliptic_lambda = 1e-4`
    - 激进对照：`elliptic_lambda = 3e-4`
  - 复查 `ell = 0.2, 0.5, 1.0` 下该小窗口是否仍存在
  - 再评估是否有必要把几何子系统从同步规范 `ADM` 切换到更严格的初边值问题形式

## 2026-04-29 局域波包平直参考后的当前任务

- 已完成：
  - 新脚本：
    `kg_examples/simulate_flat_localized_crossing_packets.py`
  - 新笔记：
    `research-notes/089-固定平直时空下局域双高斯波包在±45度方向相交的演化.md`
  - 当前固定平直参考准备态已从束状初值改为真正局域的双高斯波包
  - 在区域 `[-20,20]^2`、时间窗 `t∈[0,16]` 上，边界相对密度最大值仅 `~1.3e-28`
- 当前活跃任务：
  - 以后续 `B/C` 全动力学计算为目标，把这组新的局域波包 `(\rho,S)` 初值接入几何演化
- 下一步：
  - 不再沿用旧的束状平直参考初值
  - 在新的局域波包初值上重新组织 `B/C` 的主数值实验
  - 已完成第一轮：
    - `C(flat)` 基线校验
    - `ell=0.2, elliptic_lambda=1e-4` 的 `40` 步与 `200` 步测试
  - 当前下一步：
    - 在这组局域波包初值上继续扩展 `C-B` 分离扫描
    - 视需要再比较 `ell=0.2` 与更大 `ell` 的局域波包结果

## 2026-04-29 A/B 同初值比较后的当前任务

- 已完成：
  - 新脚本：
    `kg_examples/analyze_ab_localized_reference_vs_b.py`
  - 已修正 `B` 支初值里守恒密度 `n` 的定义，使其严格由同一组 `(rho,S,g)` 直接构造
  - 已得到局域双高斯波包初值下 `A/B` 的第一轮定量差异和 `B` 支演化图示
- 当前结论：
  - 在 `t=0.05` 的早期时间窗内，`A/B` 密度差仍较小：
    `relative L1 ~ 1.57e-4`
    `max abs diff ~ 4.47e-5`
  - 同时 `B` 支几何偏离已达
    `max|h_xx-1| ~ 1.08e-3`
  - 因而当前口径下，`B` 的几何响应比其对 `A` 参考密度分布的早期反馈更显著
- 下一步：
  - 延长 `A/B` 比较的物理时间窗口，检查该 `1e-4` 级密度差是否继续积累
  - 继续在同一组局域波包初值上推进 `B/C` 的差异扫描
  - 已新增：在继续任何长跑前，先给 `B/C` 长时间演化加 checkpoint 进度输出，定位 `t>0.2` 后失稳发生的物理时间
  - 已完成一轮实现层面复查：`yy` 分量中严格可约掉的 `exp(±2 beta)` 已全部改写为 reduced 组合
  - 下一步：
    - 用新实现重新运行 `diagnose_bc_overflow.py`
    - 判断长时间溢出是否被显著推迟、以及剩余溢出是否主要来自物质守恒密度中的 `exp(beta)`
  - 进展：
    - 改写后新的溢出诊断已稳定推进到 `t=0.70`
    - 当前仍未出现溢出或非有限值
  - 下一步：
    - 已完成：继续运行到首次失稳点
  - 当前已知首次失稳在 `t=7.79925`
  - 下一步：
    - 定位是哪一个 `exp(...)` 在 `RK4` 中间试探态先溢出
    - 区分是 `n/rho` 反演中的 `exp(beta)`，还是同步规范几何子步本身导致的中间态爆炸

## 2026-04-29 规范切换后的当前任务

- 已完成：
  - 不再继续沿同步规范作为默认长时间演化规范
  - 已把 `B/C` 原型扩展为支持 `gauge_mode = synchronous | 1plog`
  - 已在零 shift 前提下接入 `1+log slicing`
    `lapse_t = -2 lapse K`
  - 已把几何演化中的 lapse Hessian 项接入
  - 已把物质特征边界条件改成使用当前 `lapse` 判定入流/出流
- 当前结论：
  - `1+log` 下 `C(flat)` 仍严格退回 `B`
  - 在 `ell=0.2, elliptic_lambda=1e-4`、局域双高斯波包初值下，
    `1+log` 到 `t=0.55` 的检查点仍稳定，且 `lapse` 已开始在 `C` 支自适应下降到 `~0.9969`
- 下一步：
  - 已完成：`1+log` 长时间稳定性测试
  - 当前结果：
    - `1+log` 在 `t=2.911` 首次失败，且失败早于同步规范的 `t=7.79925`
    - 失败主要来自非平凡 `C` 的 lapse 与空间度规先进入强非线性，而不是 `beta`
  - 当前下一步：
    - 已完成：`B` 支单独的 `1+log` 长时间诊断
    - 当前结果：
      - `C` 支在 `t=2.911` 首次失败
      - `B` 支在 `t=6.84475` 首次失败
      - 两者都无法在当前“零 shift + 1+log”下跑到 `t=16`
    - 当前下一步：
      - 判断是优先引入更稳的 shift 条件，还是直接改成更标准的约束保持几何演化形式

## 2026-05-02 D 支同初态动态演化后的当前任务

- 已完成：
  - 按用户校正，不再把 A 参考解静态 residual 当作 D 支失败率；
  - 新增 `kg_examples/simulate_d_reduced_dynamic_same_initial.py`，从同 A 初态出发演化 D 支 reduced 物质子系统；
  - 修正 `n_cons` 符号处理，确认保符号守恒更新后 frozen baseline 可短时接近 A；
  - 确认每步代数重构 \(\tilde g[\rho]\) 会因 \(Q\) 的二阶时间差分病态而在 step 2 失败。
- 已确认更正：
  - 上一轮 reduced 动态没有把 \(\rho_A\) 映射成 \(\tilde\rho\)，而是错误使用 `rho_tilde=rho_A` 初始化守恒流；
  - 因此上一轮 frozen baseline 不能作为 D 支同初态物理演化结论，只能作为方法学诊断。
- 当前结果：
  - frozen initial \(\tilde g\)：最后可靠时间片 `t=7.4e-05`，`rho_rel_l1_support=1.58e-3`，下一步质量壳判别式穿零；
  - instant transform：step 2 即 `disc_min_support=-2.14e11`，不可作为主算法。
- 当前活跃任务：
  - 把 reduced matter evolution 保留为物质子系统；
  - 构造真正的 D 支一致初值/演化器，其中 \(\tilde g\)、scalaron 或界面/层变量独立演化；
  - 用用户三域规则组织完整模拟：weak bulk、saturated bulk、transition interface；
  - direct tensor jump/body-fitted matching 必须进入主闭合，trace jump 只能作为诊断。
- 下一步：
  - 已完成首轮：修正 D 支物质初态，使用 `sqrt(|g~|)rho~ = |X|rho_A/m^2` 正测度关系初始化；
  - 已完成首轮：重新定义 A/D 比较 observable，输出测度密度 `sqrt(|g~|)rho~` 和守恒流密度 `n_cons`；
  - 第一优先：把所有 D 支 \(\tilde T_{\mu\nu}\) residual/jump/tensor 诊断中的 `rho_A` 源项替换为 \(\tilde\rho\) 或 \(\tilde N\)；
  - 已完成首轮：检查饱和区残留 \(-M_P^2 f\tilde g_{\mu\nu}/2\) 的尺度；`ell=M_P=300` 时 plateau 不被参数压低；
  - 设计独立几何/界面变量版本，避免用动态 \(\rho\) 的 \(\partial_t^2\sqrt{\rho}\) 重构 \(Q\)；
  - 在质量壳判别式接近零的区域，比较换时间切片、body-fitted 坐标和边界匹配三种处理；
  - 把同初态 D 演化结果与 A 支精确演化持续比较，优先输出密度差、相位梯度差、判别式和界面 solvability 的时间序列。

## 2026-05-02 修正 rho 映射后三域动力学首轮后的当前任务

- 已完成：
  - 新脚本 `kg_examples/simulate_d_tridomain_full_dynamics.py`；
  - 新笔记 `research-notes/131-D支修正rho映射后三域动力学首轮.md`；
  - 用 \(\tilde N=|X_A|\rho_A/m^2\) 初始化 \(\tilde\rho\)；
  - 在 `n=64,dt=1e-6,ell=M_P=300` 下跑到可信截止 `t=5.7e-05`；
  - 输出图和数据在 `visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/`。
- 当前结果：
  - `measure_rel_l1_support=2.29e-2`；
  - `n_cons_rel_l1_support=5.92e-5`；
  - `disc_min_support=4.52e-7`；
  - 下一步 `t=5.8e-05` 质量壳判别式变负；
  - 当前支撑区全为 saturated；
  - `M_P^2/ell^2=1`，plateau 项不小。
- 当前活跃任务：
  - 把旧 residual/jump/tensor 诊断的源项全部改成 \(\tilde\rho\) 或 \(\tilde N\) 后重算；
  - 扫描 \(\ell/M_P\)，特别是 \(\ell\gg M_P\)，检验饱和区 inert 近似是否改善；
  - 建立真正的 saturated bulk metric 规则：固定代表、最小曲率代表、边界匹配代表或 scalar-tensor 变量；
  - 把 direct tensor jump/body-fitted matching 从诊断升级为强边界条件；
  - 处理质量壳判别式穿零：换切片、body-fitted 坐标或界面匹配。
- 下一步建议顺序：
  - 先跑 `ell/M_P` 扫描：固定 `M_P=300`，比较 `ell=300,1000,3000` 的 plateau 尺度、判别式寿命和 D-vs-A observable；
  - 然后重算 D 支 tensor residual/jump 的 \(\tilde\rho\) 版本，检查旧 direct tensor accepted fraction 是否显著变化；
  - 最后进入真正 body-fitted/interface-strong solve。

## 2026-05-02 1550nm 物理标定后的当前任务

- 已完成：
  - 新增 `kg_examples/physical_units.py`；
  - 修改 `kg_examples/simulate_d_tridomain_full_dynamics.py`，支持 `--physical-optical` 模式；
  - 确定 1550nm、线宽约 \(k_0/10\)、旧波包比例缩放、自洽 \(m=\omega_0/10\) 的参数；
  - 记录到 `research-notes/132-1550nm物理标定参数准备.md`。
- 当前确定参数：
  - \(k_0=0.799898054\,{\rm eV}\)；
  - \(\omega_0=0.803927793\,{\rm eV}\)；
  - \(m=0.080392779\,{\rm eV}\)；
  - \(\sigma_\parallel=2.368\,\mu{\rm m}\)；
  - \(\sigma_\perp=1.776\,\mu{\rm m}\)；
  - 域半宽 \(29.603\,\mu{\rm m}\)；
  - 相遇时间 \(42.10\,{\rm fs}\)；
  - \(M_P=1.220890128\times10^{28}\,{\rm eV}\)。
- 阻塞确认：
  - 已确认：\(\ell/l_P\) 按推荐扫描；
  - 已确认：波函数按 KG 内积归一化。
- 已完成：
  - `--normalize-kg` 已加入；
  - 已完成 `ell/l_P=1e29,1e35,1e45,1e60` 的 `n=32` 短程扫描；
  - 已完成 `ell/l_P=1e60,n=64` 的物理标定首轮可信截止运行；
  - 已完成 active/support 边界诊断，并新增 `--trusted-erosion`、`--stop-mask support|trusted`。
- 当前结果：
  - 短程扫描图：`visualizations/d_physical_1550nm_kg_scan_n32_steps10/physical_1550nm_kg_ell_scan_summary.png`
  - `ell/l_P=1e60,n=64` 输出：`visualizations/d_physical_1550nm_kg_ell1e60_n64_stable17/`
  - 可信截止：`step=17`，`t=8.39e-5 fs`
  - `measure_rel_l1_support=1.542e-2`
  - `n_cons_rel_l1_support=1.809e-5`
  - `disc_min_support=2.60e-9`
  - active 边界修正后，`active_dilation=2,trusted_erosion=1,stop_mask=trusted` 可完成 400 步；
  - 400 步时原始 support 外圈坏掉：`measure_rel_l1_support=17.981`、`disc_min_support=-5.554e-8`；
  - 400 步时 trusted core 仍稳定：`measure_rel_l1_trusted=1.437e-4`、`n_cons_rel_l1_trusted=9.622e-5`、`disc_min_trusted=3.867e-6`。
- 下一步：
  - 当前默认下一轮数值设置：`ell/l_P=1e60`、`active_dilation=2`、`trusted_erosion=1`、`stop_mask=trusted`，同时继续报告原始 support 指标；
  - 已完成收敛检查：`dt_old=2.5e-7/1.25e-7` 几乎无差异，`n64->n96` 是决定性改善，`n96->n128` 半窗口只小幅变化；
  - 当前默认工作分辨率改为 `n=96`，默认时间步仍用 `dt_old=2.5e-7`；
  - 用户提醒后已确认：上述默认只适合 \(t=0\) 早期窗口，不能直接外推到干涉区；
  - 已完成 `t_old=8` 干涉窗口检查：`n=96` 可跑完但 trusted 判别式只剩 `4.198e-7`，误差约 `1e-3`；
  - `t_old=8,n=128` 没有改善，trusted 判别式降到 `6.064e-8`；
  - 已完成 `t_old=8,n=96,dt_old=1.25e-7` 时间步复查：与 `dt_old=2.5e-7` 几乎逐项相同，因此干涉区问题不是时间步太大；
  - `t_old=8,n=128` 也没有改善，因此不能把干涉区问题归为普通空间分辨率不足；
  - 已完成 support 边缘表达后处理：binary support、trusted core、support edge band、tapered support 四种口径均显示 `n96 dt` 与 `n96 dt/2` 几乎相同；
  - 当前固定 `n=96, dt_old=2.5e-7`，暂不引入动态参数；
  - support 边缘表达继续研究，但当前以 `kg_examples/postprocess_d_support_edge.py` 后处理诊断为主：平滑权重/tapered support、body-fitted 低密度边界或连续权重诊断均不能改变源项、方程或物理边界条件；
- 当前新的优先主线：诊断并升级 `inert_tridomain` 在干涉区的几何闭合，包括饱和区 metric 代表、plateau algebraic 残差、质量壳判别式裕度、interface/bulk matching 与 metric/scalaron 自洽调整；
- full tensor interface matching 已从诊断原型接入主闭合链，作为 saturated-bulk 代表延拓的额外物理锚点；
  - 已进一步明确：饱和区代表 metric 应作为接口约束下的受限延拓来选，而不是简单平均或固定平直度规；
  - 已把受限延拓准则落实为 `--geometry-closure restricted_extension`：
    - 无 weak/transition/interface 物理锚点的 saturated 连通块保持原代表；
    - 有锚点时做最小 metric-gradient / harmonic 延拓；
    - 延拓候选必须通过 Lorentz、质量壳判别式、正 \(\tilde\rho\) 检查；
    - 当前 `t_old=8,ell/l_P=1e60,n=96` 被识别为全 saturated orphan bulk，故 `restricted_extension` 与 `inert_tridomain` 完全一致；
    - 已加入用户提出的低物质边界平直锚点假设：`--flat-anchor-boundary-layers`、`--flat-anchor-rho-frac`、`--flat-anchor-measure-frac`；
    - 当前 `t_old=8,ell/l_P=1e60,n=96` 中 flat boundary anchors 有 `97` 个，使 saturated bulk 不再 orphan，但直接 \(g_{\mu\nu}\) harmonic update 被 admissible 检查拒绝，`metric_extension_admissible_blend=0`；
    - 已完成第一版 ADM 分块延拓：`--metric-extension-variable-mode adm`；
    - ADM 分块下 `lorentz_fallback_count=0`，但 `adm_valid_fraction≈0.258`，候选粗糙度爆到 `1.63e53`，仍被 admissible 检查拒绝，`metric_extension_admissible_blend=0`；
    - 已修正 ADM 无效点伪重构和强制 \(+--\) 验收问题；
    - 当前 working closure 初步改为 `branch-preserving signature + ADM lapse-only`；
    - 已完成 cap 扫描，当前最佳折中为 `--metric-extension-max-rel-change 3e-4`；
    - 该闭合在 `t_old=8,n=96,steps=20` 下被接受，trusted measure `1.6315e-4`、trusted disc `4.0996e-7`、interface solved fraction `0.4429`；
    - `ADM lapse+shift + cap 1e-3` 也被接受，但 trusted measure `6.6356e-4`，偏离更大，暂不作为默认优先闭合。
  - 已确认质量壳判别是逐点判别：
    - 每个网格点计算 `discriminant=b^2-ac`；
    - 停止条件取 support/trusted mask 上的最小值，不是支撑内平均；
    - 后续报告时必须明确这是“最坏点裕度”，不是平均裕度。
  - 重算 D 支 tensor residual/jump 的 \(\tilde\rho\)、KG 归一化版本；
  - 已完成 `metric_extension_max_rel_change=1e-4,3e-4,1e-3,3e-3` 扫描，优先默认 `3e-4`；
- 继续做更长步数测试，先用 `ADM lapse-only cap=3e-4`；
- 对 ADM invalid 区画图，确认它们是时间反向/非实验室切片区域还是数值尾部；
- 若要实验强度预测，再引入 photon number / pulse energy，而不是单粒子 KG norm。
- 下一步优先：
  - 物理时间窗已重定：旧 `[0,16]` 对应当前 `[0,T_max]`，其中 `T_max=120.01529382382436 eV^-1=78.99550140570793 fs`，旧 `t_old=8` 对应 `T_max/2=39.49775070285396 fs`；
  - 已完成 `T=0,T_max/2,T_max` 三切片局部短窗试运行；
  - 用户已修正当前优先级：先建立完整 D 支模型数值计算器，不以 D-A transformed-measure 偏差作为 full tensor anchor 的主要否决条件；
  - 下一步应先比较 `T=0,T_max/2,T_max` 三切片的 inert / ADM-lapse-only / tensor-anchor 三种闭合在 D 支内部自洽指标上的表现；
  - D 支内部自洽指标优先级：质量壳逐点判别式、\(\tilde\rho>0\)、守恒量稳定、D 支场方程/约束残差、full tensor interface residual、饱和区代表 metric 规则是否与界面匹配一致；
  - `T=0` 稀疏 tensor anchor 导致 D-A 偏差 `~1.95e-2` 不再自动判作失败，而是记录为 D 支可能与 A 支产生早期可观测差异；
  - 建议把 full tensor matching 暂时作为降频诊断或有条件锚点刷新，而不是每步无条件强锚点；
  - 全 `[0,T_max]` 生产演化前，需要重新处理长时间物质 RK4 稳定性；低分辨率粗步长 `dt_old=0.1/0.01` 全窗已在早期爆炸，不能作为可行全窗设置；
  - 若仍按 `dt_old=2.5e-7` 直接覆盖旧 `[0,16]`，理论上需约 `6.4e7` 步，当前脚本不可行，必须发展多尺度/重采样/更稳的全窗推进策略。

- 2026-05-04 新增近期任务：把 D 支全窗数值器从当前 `n_cons + u_i` 显式 RK4 升级为约束保持物质积分器。
  - 已完成底层性能优化：`--geometry-every`、`--metric-extension-every`、`--diagnostics-every`、`--no-render`；
  - 已确认优化版内存峰值约 `80~105MB`，16G 机器可承受；
  - 已完成第一版局部时间图册推进：
    - `--matter-variable-mode local_time`
    - `--local-time-patchwise`
    - `--local-time-stationary-candidates`
    - `--local-time-patch-bboxes`
  - 当前局部时间结果：
    - `t_old=0,n64,20` 步：trusted 覆盖 `100%`、`disc_min_trusted=3.86e-6`、用时约 `17.7s`；
    - `t_old=8,n32,20` 步：trusted 覆盖 `100%`、`disc_min_trusted=7.05e-3`、用时约 `13.6s`；
    - 内存仍远低于 16G，主瓶颈转为 patch 数、patch 边界通量和全窗步数。
  - 已完成 local_time 步长扫描：
    - `n32` 快扫对大步长过于乐观，不能作为生产依据；
    - `n64,t_old=8,max_tilt=5` 是当前关键约束；
    - `dt_old=1e-5` 作为安全候选，`dt_old=2.5e-5` 作为激进短窗探针；
    - `dt_old=5e-5` 在干涉区质量壳裕度约 `3e-11`，暂不作为生产候选；
    - `dt_old=2.5e-5,steps=100` 可跑完但 wall time 过高，说明算法结构还需优化。
  - 下一步不应继续盲目调分辨率或 support 阈值，而应优先实现/比较：
    - 把当前 mask/recombine patchwise 原型升级为 conservative patch-boundary flux / matching；
    - 对局部时间图册做 patch 合并或候选缓存，降低 `patch_count`；
    - 接入 `restricted_extension` 与 full tensor matching 的一致性门控，避免 metric 更新路径重新使用病态实验室 `n_cons`；
    - 扫描 `dt_old`，判断 local_time 模式是否允许比 `2.5e-7` 更大的可信步长；
    - 只在完成上述约束保持后，再测试完整旧 `[0,16]` / 物理 `[0,Tmax]` 是否能在 1 天内稳定跑完。
  - 当前更具体的执行顺序：
    - 优先减少 patch 数：按空间连通分量和相近 `(a,b)` 合并 patch，或缓存候选/metric_tau；
    - 实现真正 conservative patch-boundary flux，不再只做 mask/recombine；
    - 用 `dt_old=1e-5` 在 `t_old=0,8,16` 做 `n64` 更长短窗；
    - 只有在 local_time 物质推进稳定后，再接回几何 `restricted_extension` / full tensor matching。

## 2026-05-05 D 支 local_time 全窗优化下一步

- 已完成并保留的安全优化：
  - 解析 `tau=t+a x+b z` 的协变/逆变度规变换；
  - 每步复用 `metric_inv/det`；
  - `fast-profile` 尊重 `diagnostics/interface/metric-extension every=0`，可测纯演化；
  - `coordinate_matter_evolution.py` 的守恒量/密度恢复函数支持外传缓存的逆度规与行列式。
- 已否决的提速路线：
  - `local_time_tile_size>1`；
  - 粗 `local_time_step=1.0/1.25/2.0`；
  - `local_time_step=1.5` 即使细化 stationary quantization；
  - `local_time_atlas_every>1` 作为主优化；
  - 重写轻量 RHS 微优化。
- 当前可信性能：
  - `t_old=8,n64,dt_old=2.5e-5,steps=100` 纯演化约 `21.3s`、RSS 约 `151MB`；
  - 全旧窗 `[0,16]` 用 `dt_old=2.5e-5` 粗估 `30-38h`，尚未稳进 1 天；
  - `dt_old=1e-5` 更稳但仍需数天。
- 下一步优先级更新：
  - 不再优先粗化 chart 图册，因为覆盖诊断已否决；
  - 优先把 patch RK4 小数组循环迁移到编译/向量化后端，或改成真正 conservative multi-patch flux / characteristic update；
  - 继续保持 `local_time_uncovered_trusted_measure_fraction=0`、`disc_min_trusted>0`、`local_time_n_tau_rel_delta_sum` 小量作为进入长窗的硬门槛；
  - 若短期只做物理探索，可用 `dt_old=2.5e-5` 激进探针；若做保守生产，应仍以 `dt_old≈1e-5` 作为候选。

## 2026-05-05 并行调试与下一步提速任务

- 已新增可用的短窗调试工具：
  - `kg_examples/run_d_local_time_parallel_scan.py`
  - 默认 `--quick-diagnostics --skip-fields-npz`，只用于快速判断质量壳、正密度、局部时间 coverage 和 `n_tau` 守恒诊断，不用于最终物理图。
- 调试建议：
  - 默认并行 `--max-workers 2`；
  - 如果用户机器当前负载低，可临时开到 `4`；
  - 不建议开到 `8`，因为 4 进程已经显示明显资源竞争。
- 当前短窗最快安全口径：
  - `--quick-diagnostics --skip-fields-npz --no-render --fast-profile --diagnostics-every 0 --interface-every 0 --geometry-every 100000`
  - 对 `t_old=8,n64,dt_old=2.5e-5,steps=20` 可把耗时压到约 `3.8s`。
- 当前不可作为生产默认的调试开关：
  - `--local-time-atlas-every 20`：50 步短窗可提速，但 100 步 `n_tau` 守恒诊断恶化；
  - `--local-time-tile-size>1`、粗 `local-time-step` 方案仍保持否决。
- 下一步如果继续追求 1 天内完整 `[0,Tmax]`：
  - 第一优先：实现/试验编译化 patch RK4 内核，减少大量小 patch 上的 Python/NumPy `roll` 开销；
  - 第二优先：把 mask/recombine patchwise 原型升级为真正 conservative patch-boundary flux；
  - 第三优先：考虑特征推进/局部时间步策略，而不是继续微调诊断和图册复用。

## 2026-05-05 旧切片 5 分钟诊断后的近期任务

- 已完成旧 `t=0,4,8,12` 四切片约 5 分钟端点完整诊断，结果在：
  `visualizations/d_old_slices_5min_full_diagnostics/five_min_diagnostics_summary.json`
- 下一步优先级：
  - 先定位 `t_old=8` transformed measure 爆发：检查最终场中 `measure_density_D` 的最大值/位置、低密度 floor 区、干涉相消区、patch label 边界与 `u_t/j_tau` 的关系；
  - 再检查 `t_old=0` 的 `n_tau≈2.419` 和 coverage missed `0.4458%`，判断是否需要更完整局部时间候选或减小步长；
  - 对 `t_old=12` 检查 tapered/support-edge 偏差是否集中在低密度尾部，如果是，生产比较应更重 trusted 核心并改进支撑边缘处理；
  - 不应把“4 个任务都 completed”误读为可长窗生产，当前至少 `t_old=8` 已出现严重数值爆发。

## 2026-05-05 高分辨率 `[-9,9]um` D 局部窗口后的近期任务

- 当前必须采用的新初始化流程：
  - 在大域 A 正频 KG 参考上生成快照；
  - 裁剪到物理 `[-9,9]um` 窗口；
  - 再用裁剪后的 `rho,S,gtilde` 初始化 D 支；
  - 不得把小窗口本身当作周期 FFT 初始域。
- 当前脚本：
  `kg_examples/simulate_d_local_window_from_a_snapshot.py`
- 当前密度比较硬规则：
  - 主比较是 `rho_A` vs `rho_D_to_A=m^2 ntilde_D/|X_g[D]|`；
  - 不能直接把 `ntilde_D` 或 `rho_tilde_D` 与 `rho_A` 比。
- 当前可信分辨率基准：
  - `n_full=640` 大域裁剪后窗口实际 `195x195`；
  - `dx≈0.0925um`；
  - 每个干涉条纹约 `11.85` 点。
- 当前时间步判断：
  - `dt_old=2.5e-5` 在高分辨率干涉中心过大，一步后会破坏质量壳；
  - 局部短窗诊断应使用 `dt_old=2.5e-7` 或更小，除非新的积分器证明可放大步长。
- 当前诊断结果：
  - `tau=-3.5` 分离态可跑满 5 分钟，`rho_D_to_A` core/support 相对 L1 约 `5e-4`，质量壳判别式保持正；
  - `tau=0` 干涉中心只能跑约 `28-29` 步，trusted 区出现 `~1e-8` 量级负判别式；
  - 干涉中心 core 区 `rho_D_to_A` 相对 L1 仍约 `7e-6`，但 support/trusted L1 被低密度/近 `X_g[D]=0`/节点边缘主导。
- 已完成后续定位与修正：
  - 原 `tau=0` 失败/爆点主要来自低密度 support 边缘被当作 bulk 单元显式推进；
  - `simulate_d_local_window_from_a_snapshot.py` 已新增 `--evolve-mask active|support|trusted`；
  - 当前推荐短窗物理诊断使用 `--evolve-mask trusted --freeze-uncovered-chart --chart-boundary-halo 1`；
  - 新口径下 `tau=0` 和 `tau=-3.5` 均可跑满 5 分钟，质量壳判别式保持正，`rho_D_to_A` L1 约 `2e-4~4e-4`。
- 下一步优先：
  - 第一版 `trusted-buffer` conserved-current matching flux 已接入；下一步要验证它能否替代 `face` 作为默认局部窗口边界闭合；
  - 对 `tau=+3.5` 也跑同口径 5 分钟，确认分离后状态是否同样稳定；
  - 检查 `[-9,9]um` 在分离态 support 触边的影响，必要时用稍大窗口或只在 `tau=-2.5..+2.5` 做生产短窗；
  - 在边界 flux 明确后，再尝试增大 `dt_old` 或进入更长时间窗，不应回到 active 全域显式推进。

## 2026-05-05 算法状态核查后的当前任务

- 当前状态：
  - 局部短窗可信诊断算法已可用；
  - 完整 D 支全窗生产算法尚未闭合；
  - 关键缺口是正式 conservative boundary/interface flux，而不是再回到 active 全域显式推进。
- 已完成：
  - `rk4_patchwise_local_time_matter_step` / `rk4_local_time_matter_step` 已支持 `stencil_mask` 与 `write_mask` 分离；
  - `simulate_d_local_window_from_a_snapshot.py` 已支持 `--boundary-stencil-mask write|support|active`；
  - `tau=0` 的 20 秒、120 秒、165 步公平对照和 5 分钟 fixed-buffer 运行均完成。
- 当前第一优先：
  - fixed-buffer stencil 已升级出第一版 `matching` conserved-current interface flux；
  - 下一步不再是继续加数值耗散，而是比较 `matching_flux_weight=equal|abs_n`、`tau=±3.5` 和更长步数；
  - 继续跟踪 `disc_min_trusted`，因为 fixed-buffer 5 分钟虽未失败，但质量壳裕度已接近 `1e-11`。
- 已完成新增：
  - `--boundary-flux-mode central|face|rusanov`；
  - `face` 模式为显式有限体积面通量；
  - `rusanov` 模式为带数值耗散的保守通量候选；
  - `tau=0` 的 face 20 秒、Rusanov 20 秒、Rusanov 120 秒均已跑通。
  - `--boundary-flux-mode matching`；
  - `--matching-flux-weight equal|abs_n`；
  - `matching` 只在 `trusted-buffer` 面上用法向守恒流连续条件替换面通量，且用全局写回区识别物理界面，避免 patch/chart 边界误判；
  - `tau=0` 的 matching 20/60 步已跑通，60 步 `rho_pullback_rel_l1_trusted≈7.43e-05`、`disc_min_trusted≈8.50e-09`。
- 下一轮验证顺序：
  - Rusanov 耗散强度扫描已完成，结论是不作为默认物理算法；
  - 当前守恒通量基线应回到 `face` 或 `rusanov_strength=0`；
  - 第一版 matching flux 已完成，且 `tau=-3.5,0,+3.5` 初步对照已跑；
  - `tau=-3.5` 和 `tau=0` 在同 60 步下 matching 与 face 近乎同阶，`tau=+3.5` 两者均被 300s wall time 截断但未物理失败；
  - 下一步若继续数值闭合，应把 full tensor Einstein interface matching 作为耦合界面求解器接入，而不是只在物质连续性方程中继续微调标量通量；
  - 后续再把 full tensor Einstein interface matching 作为强边界条件接入，不能把当前物质通量 matching 误称为 full tensor matching；
  - 在上述闭合稳定后，再做 `dt_old` 放大扫描和全窗耗时估计。

## 2026-05-06 full tensor interface 接入后的当前任务

- 已完成：
  - `simulate_d_local_window_from_a_snapshot.py` 已能在高分辨率局部窗口上输出 full tensor interface rows、accepted/rejected 统计、accepted anchor mask 和诊断图；
  - `tau=0` 的 0 步与 20 步 full tensor 诊断已跑通，accepted fraction 约 `27.6%`；
  - 结果说明当前不是 full tensor 最小二乘“完全无解”，而是只有部分界面段满足当前局部张量 jump condition。
- 当前仍未完成：
  - full tensor accepted anchor 还没有作为强边界条件真正反馈到局部窗口的 metric/interface 演化中；
  - local-window 版本仍是固定初始 metric 的 matter-only short-window evolution 加 strong tensor diagnostic；
  - rejected tensor 段的处理规则仍需从“界面形状调整 / 两侧 bulk metric 调整 / 饱和区代表选择”三者中闭合。
- 下一步优先：
  - 把 `tensor_anchor_mask` 和 `tensor_interface_rows.jsonl` 作为几何闭合器输入，先实现 accepted 段参与 saturated-bulk metric representative selection；
  - 同时明确 rejected 段不强锚定，而是作为 unresolved interface diagnostic；
  - 若继续本地局部窗口路线，优先实现“accepted tensor anchors + flat low-matter boundary + weak/transition anchor”的 metric representative update，再检查质量壳和守恒量是否仍稳定。

## 2026-05-06 D 支初值投影任务

- 用户新目标：
  - 重新确定 D 支初值；
  - 目标函数：\(\tilde\rho\) 拉回 \(g\) 表象后的 \(\rho\) 与 A 支 \(\rho_A\) 的逐点加权偏离最小；
  - 约束：\(\tilde g\) 表象满足 D 支方程。
- 已完成首轮可行性检查：
  - 新脚本 `kg_examples/analyze_d_initial_projection_feasibility.py`；
  - 对固定 A-derived geometry 和 \(u_\mu\) 的 \(\tilde\rho\) 全局缩放做扫描；
  - 最优缩放仍为 `alpha=1`，而修复 tensor jump 所需物质量级约 `10^62` 倍，说明不能靠调 \(\tilde\rho\) 完成约束投影。
- 下一步：
  - 把初值投影变量从 \(\tilde\rho\) 扩展到 \(\tilde g\) 的 metric jets/界面几何；
  - 优先做低维参数化：保持边界与弱场区固定，只在饱和区/过渡层附近调整 `raw_y` 或等价 metric representative；
  - 目标函数应同时包含：
    `rho_pullback_weighted_error + lambda_tensor * full_tensor_residual + lambda_geom * metric_smoothness/branch-preservation penalty`；
  - 输出必须报告拉回 \(\rho\) 偏差、accepted tensor fraction、rejected residual、质量壳判别式和 Lorentz signature 保持情况。

## 2026-05-06 D 支方程认证后的任务修正

- 已新增认证脚本：
  `kg_examples/certify_d_branch_equations.py`
- 对当前 `tau=0`、20 步高分辨率局部窗口案例认证结果：
  - `certified=false`；
  - 通过：质量壳、正 \(\tilde\rho\)、拉回 \(\rho\) 加权目标；
  - 失败：full tensor interface 与 saturated-bulk algebraic equation。
- 当前不能声称：
  - 不能说当前初态满足 D 支完整场方程；
  - 不能说当前后续演化保证满足 D 支完整场方程；
  - 不能说当前数值方法是完整 D 支求解器。
- 下一步若用户要求“保证满足 D 支方程”，必须先实现：
  - 初值层：约束求解器，直接求解 full tensor interface 与 saturated-bulk algebraic/bulk equation；
  - 演化层：每步隐式或投影格式，把 D 支残差作为 hard acceptance criterion；
  - 输出层：每步写 `d_branch_equation_certificate`，若认证失败则拒绝该步。

## 2026-05-06 用户澄清后的保证性问题

- 用户澄清：
  - 不是问当前已有结果是否满足 D 支方程；
  - 是问此前提出的“大 \(\lambda\) 惩罚投影方案”能否保证初态和后续演化满足 D 支方程，若能则实施。
- 当前结论：
  - 大 \(\lambda\) 惩罚方案本身不能保证，只能近似压低残差；
  - 若要 guarantee，必须转成等式约束求解/隐式约束投影，并以 certificate 作为硬验收；
  - 当前可实施的是“certified constrained workflow”的骨架：候选初值/步进后必须通过 certificate，否则拒绝；但还缺少能把失败候选投影到约束面的 full D equality solver。
- 若继续推进，实现顺序应是：
  - 第一阶段：把 D 方程离散残差函数统一成 `E_D(state)`，不再散落在诊断脚本中；
  - 第二阶段：实现 coarse/local patch equality-constrained initial projection，变量至少包括 metric representative/jets；
  - 第三阶段：实现每步 implicit/projected update，step accept 条件是 `certificate=true`；
  - 第四阶段：再谈长时间演化和与 A 支比较。

## 2026-05-06 真正 D 初值求解器后的当前状态

- 已实施：
  - `kg_examples/solve_d_initial_data_constrained.py`
  - 它把 D 初值求解变成 hard-constrained gate，而不是大惩罚近似；
  - 若当前 D 方程必要条件不满足，则输出 `status=infeasible` 和图/JSON 报告，并拒绝进入演化。
- 当前结果：
  - `tau=0` 初始切片和旧候选 20 步切片都返回 `status=infeasible`；
  - 不可行不是局部少量网格点，而是支撑质量 `100%`；
  - 直接原因是 saturated bulk 代数张量方程要求满秩度规等于 on-shell 标量物质的秩一张量。
- 当前不能继续做的事：
  - 不能把已有 matter-only/local-window evolution 当成 D 支完整演化；
  - 不能在没有可行 D 初态时报告“D 初值后的偏差演化”。
- 下一步若继续 D 支路线，需要先做理论/模型层面的选择：
  - 重新审查 saturated bulk 近似是否过强，是否必须保留 \(f_R\) 导数项或全 f(R) PDE；
  - 或修改/减去 plateau cosmological-like 项，使饱和区不是 \(-\frac12 f\tilde g_{\mu\nu}=\tilde T_{\mu\nu}/M_P^2\) 的秩冲突形式；
  - 或改变物质应力张量/引力耦合，使 RHS 不是纯秩一 dust-like 张量；
  - 做出选择后再实现新的 equality-constrained 初值/演化求解器。

### 修正后的当前任务口径

- 上一段中 `status=infeasible` 只适用于 strict plateau algebraic gate，不再代表完整 D 支方程不可行。
- 当前最新状态：
  - `solve_d_initial_data_constrained.py` 已改为默认 `requires_full_fR_equation_solve`；
  - `diagnose_d_full_fr_equation_local_window.py` 已能直接计算完整 metric f(R) 方程项分解；
  - `tau=0` 诊断显示：当前 A-derived 初态在支撑区 \(\phi=f_R=0\) 且导数项为 `0`，所以该候选本身不满足完整 D 方程，但这不是完整 D 模型的数学不可行证明。
- 下一步更准确的任务：
  - 不能再用 pure plateau algebraic equation 做 D 初值硬约束；
  - 应把初值求解变量扩展为 \(\tilde g\) 及其时间 jets，使 \(\phi\) 及其导数项可以非零并参与补足张量结构；
  - 初值求解器应以完整 \(f(R)\) 残差作为 hard residual，而不是以 plateau rank condition 为 hard residual；
  - 需要考虑高精度/重参数化 raw_y/phi，因为当前 `ell/l_P=1e60` 让 \(\phi=\mathrm{sech}^2(raw_y)\) 在双精度下完全下溢。

## 2026-05-07 metric \(f(R)\) 方向匹配后的近期任务

- 已完成：
  - 新增 `kg_examples/diagnose_metric_fr_direction_matching.py`；
  - 对 `tau=-3.5,0,+3.5` 三个切片完成 direction matching 诊断；
  - 已确认必须同时看：
    \(\phi\tilde R_{\mu\nu}\)、
    \(\phi\tilde R_{\mu\nu}-\frac12 f\tilde g_{\mu\nu}\)、
    以及包含 \(\phi_R,\phi_{RR}\) 的完整 metric \(f(R)\) 局部结构。
- 当前判断：
  - \(-\frac12 f\tilde g_{\mu\nu}\) 项显著改善匹配；
  - 当前 D 函数不匹配；
  - 纯 \(f(\tilde R)\) 的最大疑点是反推出的 \(\phi,f\) 不像 \(\tilde R\) 的单值函数。
- 下一步优先：
  - 第一优先：对反推出的 \(\phi,f\) 做更干净的单值性分析，排除坐标 Frobenius 内积和低密度/节点区域造成的假离散；
  - 第二优先：把目标从 \(\tilde T_{\mu\nu}/M_P^2\) 扩展到非平直弱背景中的 \(R_{\mu\nu}\)，检查“逼近 A 支 Einstein 方程”是否比平直背景目标更容易；
  - 第三优先：如果单值性仍失败，转向更宽的几何作用量，例如依赖 \(\tilde R_{\mu\nu}\tilde R^{\mu\nu}\) 或其它张量不变量，而不是只依赖 \(\tilde R\)。

### 单值性复查后的更新

- 已完成第一优先：
  - 增加高密度核心区；
  - 增加 invariant trace/Ricci contraction 反推；
  - 输出 `metric_fr_single_value_tau_*` 三组目录。
- 当前结论：
  - 低密度/节点确实污染单值性；
  - 但在 `core_rho_10pct` 中，只有干涉中心 `tau=0` 显示出较接近单值的趋势，分离态仍明显不够；
  - invariant 反推通常残差更大，说明 Frobenius 局部拟合的乐观结果不应直接解释成存在纯 \(f(\tilde R)\)。
- 下一步优先级调整：
  - 若继续纯 \(f(\tilde R)\)：先做函数拟合并把拟合出来的 \(f(\tilde R)\) 重新回代全张量方程，而不是只看局部逐点最优；
  - 若允许扩展作用量：优先加入第二个标量不变量，例如 \(\tilde R_{\mu\nu}\tilde R^{\mu\nu}\)，检验二变量函数是否显著恢复单值性；
  - 仍需做非平直弱背景 \(R_{\mu\nu}\) 目标测试，以免平直背景 \(T/M_P^2\) 目标过于特殊。

### 保留二阶导数项后的普适 \(f(R)\) 拟合更新

- 已完成：
  - `kg_examples/fit_metric_fr_universal_function.py`；
  - 单一 Chebyshev/asinh \(f(\tilde R)\) 同时拟合 `tau=-3.5,0,+3.5` 的 `core10`；
  - 已包含 \(f_R,f_{RR},f_{RRR}\) 诱导的完整 metric \(f(R)\) 二阶导数项；
  - 已做列归一化，避免第一版矩阵尺度病态。
- 结果：
  - degree 6/10/14/18 均给出加权相对残差约 `1`；
  - core10 三切片 p50 residual 仍约 `1`；
  - 普适函数 LHS 量级只达到 RHS 的 `3e-5~5e-5`。
- 当前优先级：
  - 不应继续只提高纯 \(f(\tilde R)\) 多项式阶数；
  - 下一步若继续几何作用量路线，应加入第二不变量，例如 \(I_2=\tilde R_{\mu\nu}\tilde R^{\mu\nu}\)，做 \(f(\tilde R,I_2)\) 的同类普适函数检验；
  - 或者转向非平直弱背景 \(R_{\mu\nu}\) 目标，检查平直量子实验目标是否过强。
- 高阶复查：
  - 单变量 \(f(\tilde R)\) 已从 degree `18` 扫到 `22,26,30,34`；
  - degree `34` weighted residual 仍约 `0.99945`，core10 p50 residual 仍约 `1`；
  - 条件数已到 \(10^{14}\)，不建议继续单纯提高阶数。

### 单位与量纲待修正项

- 需要在后续脚本中显式区分：
  - 非约化 Planck mass \(1.2209e28 eV\)；
  - 约化 Planck mass \(2.4353e27 eV\)，适用于 \(M_P^2G=T\) 口径。
- 需要为 2+1d 截面到 3+1d 光束解释加入横向厚度/模式函数因子 \(L_y\)；
- 重新跑最终候选时，报告应同时给：
  - 采用的 \(M_P\) 约定；
  - \(L_y\) 选择；
  - 这些全局尺度重标定对绝对 \(f\) 系数与 RHS 量级的影响。
- 但当前判断：这些是绝对量纲修正，不太可能改变纯 \(f(\tilde R)\) 的相对残差/单值性结论。

### \(f(R,I_2)\) 候选更新

- 已完成：
  - `kg_examples/diagnose_metric_fr_ricci2_direction_matching.py`
  - `kg_examples/fit_metric_fr_ricci2_universal_algebraic.py`
  - `kg_examples/fit_metric_fr_ricci2_universal_full.py`
- 当前判断：
  - \(I_2=\tilde R_{\mu\nu}\tilde R^{\mu\nu}\) 明显改善局部代数张量方向匹配；
  - 但普适二维代数函数 \(f(R,I_2)\) 回代仍没有通过；
  - 完整 metric \(f(R,I_2)\) 导数项回代后，三切片共同拟合仍没有通过，`n=160` weighted residual 仍约 `0.994~0.998`。
- 下一步优先：
  - 暂时不建议继续单纯提高二维 Chebyshev 阶数；
  - 若继续几何作用量路线，应考虑更多不变量，例如 Riemann/Kretschmann、Weyl-like 结构、或与 \(u_\mu\) 诱导方向有关但仍能从 \(\tilde g\) 几何中抽取的标量/张量；
  - 若要继续 numerical scan，先向量化 `fit_metric_fr_ricci2_universal_full.py` 的张量二阶协变导数，否则高分辨率扫描会很慢。
- 最新高阶检查：
  - `degree=4:4` 三切片 `n=160` 已完成，无实质改善；
  - `degree=5:5` 原型运行过慢且 `4:4` 已无改善，已停止；
  - 若继续扫高阶，任务应先变成“优化完整导数项回代器”，再做 `5:5+` 扫描。

### 新路线：允许引力作用量显含 \(u,r\)

- 已修正理论基准：
  - 若 disformal 变换可逆，则“完整总作用量整体拉回”可与 A 支完全等价；
  - 但仅拉回 EH action \(S_{\rm EH}[G(\tilde g,u,r,\rho,\ldots)]\) 后配上本项目另行构造的 `gtilde` 物质 action，是 hybrid action，不自动等价；
  - 因此后续必须直接检验 hybrid action 的完整场方程，而不能只引用 EH 拉回。
- 下一步优先任务：
  - 推导并实现低阶
    \[
    A^{\mu\nu}(\mathcal I)\tilde R_{\mu\nu}-2V(\mathcal I)
    \]
    的 fixed-background residual fitter；
  - \(A^{\mu\nu}\) 首先取
    \(a_0\tilde g^{\mu\nu}+a_1\hat u^\mu\hat u^\nu+a_2\hat r^\mu\hat r^\nu+a_3\hat u^{(\mu}\hat r^{\nu)}\)；
  - 必须包含完整 metric variation 中的 \(\nabla\nabla A^{\mu\nu}\) 项；
  - 使用同一 1550nm Gaussian interference 三切片 `tau=-3.5,0,+3.5`、`core10/trusted` 口径，与 pure \(f(R)\)、\(f(R,I_2)\) 残差直接比较；
  - 若该低阶 ansatz 通过，再继续做 Cauchy 主部、物质变分补偿项和 D 支投影初值；
  - 若仍失败，再加入 \(\nabla u,\nabla r\) 诱导的张量方向。

### 新增优先任务：从变换后的 A 支 Einstein 方程提取剩余张量

- 新优先级高于直接猜 \(A^{\mu\nu}R_{\mu\nu}\)：
  - 先实现 \(C^\alpha{}_{\mu\nu}=\Gamma[G]-\tilde\Gamma\)；
  - 第一版直接用 A 参考中已知的 \(g\)，不要先尝试求解 \(g=G[\tilde g,u,r,\rho]\) 这个反问题；
  - 计算 \(G_{\mu\nu}[G]=\tilde G_{\mu\nu}+\mathcal H_{\mu\nu}\)；
  - 计算
    \[
    \mathcal R_{\mu\nu}
    =
    T^{A}_{\mu\nu}
    -T^{(\tilde m)}_{\mu\nu}
    -M_P^2\mathcal H_{\mu\nu};
    \]
  - 在 1550nm Gaussian interference 的 `tau=-3.5,0,+3.5` 三切片上分析 \(\mathcal R_{\mu\nu}\) 的张量方向、量级、单值性与主支撑区/节点区分布；
  - 再判断 \(\mathcal R_{\mu\nu}\) 是否可被 pure-\(\tilde g\) 几何项吸收，或是否必须引入 \(u,r\)-dependent gravitational sector。

### 2026-05-07 当前活跃任务：gBCD 辅助应力场闭合

- 已完成：
  - `kg_examples/fit_gbcd_auxiliary_gauge.py`；
  - nullspace hard constraint 求解器；
  - 三切片 `n=96, trusted, force_mode=full` 的 `norm/time` 代表元检验；
  - `kg_examples/diagnose_gbcd_conservation_principal_symbol.py`；
  - 确认 `2+1` 约化下 full conservation 主符号 rank=3、nullity=1；
  - `kg_examples/fit_gbcd_trace_closure.py`；
  - 否定简单 `trace=F(Q)` 与 `D=0` 代数闭合作为状态方程；
  - `kg_examples/fit_gbcd_auxiliary_gauge.py` 已升级 `--q-gated-norm`；
  - 三切片 \(Q\)-门控辅助范数测试已完成，确认该路线比简单 trace/Q 代数闭合更自然；
  - `research-notes/151-gBCD辅助场变分闭合的Euler-Lagrange方程.md`；
  - tau=0 \(Q\)-门控强度扫描已完成，输出 `visualizations/equation_first_gbcd_aux_qgated_scan_n96_tau0_summary/`；
  - 无 Q 门控最小辅助 action 扫描已完成，输出 `visualizations/equation_first_gbcd_aux_noq_time_scan_n96_3tau_summary/`；
  - 已写出 `research-notes/153-gBCD显式理论候选v0.md`；
  - 已写出 `research-notes/164-gBCD投影型Einstein-like方程候选.md`，把 gBCD 从 \(A,B,C,D\) 逐点代表元推进为显式投影型场方程。
  - 已做 v0 三层离散时间差分诊断，确认 `rr/ur` 是时间推进风险分量；
  - 研究笔记 `research-notes/148-gBCD辅助应力场代表元规范与nullspace硬约束.md`。
  - 研究笔记 `research-notes/149-gBCD状态方程trace与剪切闭合首轮检验.md`。
  - 研究笔记 `research-notes/150-Q趋零条件作为辅助场变分权重.md`。
  - 研究笔记 `research-notes/152-无Q门控的最小辅助场action扫描.md`。
- 当前结论：
  - gBCD/full-conservation 张量壳仍可行；
  - KKT 病态已经用 nullspace 投影解决；
  - 空间平滑代表元会污染代数场方程，暂不能作为物理状态方程；
  - 时间平滑只是温和 numerical gauge，不是最终本构。
  - full conservation 本身不能唯一产生四个 \(\lambda_I\)，还缺一个状态方程/规范条件。
  - 简单代数状态方程会显著恶化残差或引入病态条件数。
  - \(Q\to0\) 回 Einstein 应作为辅助场 action 的权重/边界条件，而不是 trace=\(F(Q)\)。
  - 当前闭合应理解为 saddle-point 辅助场变分问题；conservation 乘子对应 nullspace hard constraint。
  - 暂时默认不加 Q 门控、不加空间项，采用 `norm_weight=1e-8,time_weight=1e-5`。
  - gBCD v0 是 equation-first 显式理论候选，不是最终 action-first 理论。
  - lambda-only 一步推进已诊断：hard conservation 能满足，但固定 A 背景下下一切片代数场方程残差仍约 `14%~17%`，不能作为完整 D 支动力学。
  - metric time-principal 主部 rank=3，说明 gBCD 场方程自然分裂为 3 个约束组合 + 3 个演化组合。
  - principal-constraint projection 已跑通，可把约束残差压到 `1e-4` 以下量级，且所需相对 \(\delta\tilde g_+\) 很小。
  - 但 pointwise principal \(\delta\tilde g_+\) 完整代回后 exact residual 仍有 `3.7%~34%`，说明还需要全局 full-linear metric update。
  - `n=96` 分辨率太低，每波长只有 `2.51` 点；`n=384` 局部 `117x117`，每波长 `10.05` 点，是下一档最低可信条纹分辨率。
  - dense full-linear metric update 在 `n=96` 已跑通，`ridge=1e-6` 最稳，tau=0 exact residual wmean `0.00661`。
  - sparse full-linear 算子装配正确到 `1e-16`，但当前 CG/LSQR 迭代仍差于 dense direct solve。
  - sparse principal projection 当前是 soft penalty 版；`n=96 core10` sparse pipeline exact residual wmean `0.07885`，dense hard pipeline 为 `0.02548`。
  - `n=384,tau=0,core10` sparse pipeline 已跑通，exact residual weighted mean `0.00644`，p95 `0.04337`，corrected metric determinant 未退化。
  - `n=384,tau=-3.5,core10` 旧快速尝试未收敛，exact residual wmean `0.40698`；本轮向量化 solver + 3000 步 projection + 6000 步 metric CG 后降到 `0.14577`，但仍是三切片主卡点。
  - `n=384,tau=+3.5,core10` 已完成同口径 projection/update，exact residual wmean `0.03864`，p95 `0.10831`，corrected metric determinant 全正。
  - full-linear metric update 与 sparse projection 的 matvec 已改为 entry-array 向量化，当前高分辨率三切片能在数分钟内调试。
  - `tau=-3.5` 阻尼线搜索显示 `alpha=1` 最优，欠阻尼不能解决残差；左侧分离态问题不是简单 Newton 过冲。
- 下一步最高优先：
  - 暂停把 `plus-only` 包装成主目标，先完成投影型 Einstein-like 方程的理论闭合检查；
  - 已明确投影内积 \(W^{\mu\nu\rho\sigma}\) 对在壳解集不是物理自由参数，下一步只需记录不同 \(W\) 对离壳诊断/代表元的影响；
  - 已得到 Gram 矩阵 \(H_{IJ}\) 的退化条件，下一步要把 \(\Delta=0\) 写成 patch/branch 规则，避免把数值 `guard` 当成物理方程；
  - 检查 \(Q\to0\) 分支条件能否由 \(Q\)-门控辅助泛函、边界条件或正则化选择自然给出；
  - 已检验 \(\mathsf q_E\)-trace 最小硬闭合，结果不如普通 trace=0；也已确认 trace=0 的 \(w^2=0\) 风险主要集中在 \(r^\perp\) 退化 patch；
  - 已写出普通 trace=0 + full conservation + patch/branch 的显式最小方程组，见 `research-notes/169-trace0最小投影方程组与patch条件.md`；
  - 已完成 trace=0 投影方程的 Helmholtz/self-adjoint 必要条件检查，见 `research-notes/170-trace0投影方程的Helmholtz检查.md`；
  - 已写出最小辅助应力场 action 原型，见 `research-notes/171-最小辅助应力场action原型.md`；
  - 已分析 shadow-field 分支，见 `research-notes/172-shadow-field分支的最小方程组.md`；
  - 已检查 shadow stress 和物理分支传播条件，见 `research-notes/173-shadow-stress与物理分支传播条件.md`；
  - 已检查具体 stealth-state action 候选，见 `research-notes/174-stealth-state-action候选的最小no-go.md`，结论是最小局域 ansatz 不能生成非零 trace=0 投影源；
  - 已完成 `n=384` 可信条纹分辨率上的 trace=0 局部代数复检，见 `research-notes/175-n384-trace0局部复检.md`；
  - 复检结论是普通 trace=0 比 \(\mathsf q_E\)-trace 更稳定，继续作为主标量闭合；\(\mathsf q_E\) 暂降级为主符号/正定范数诊断工具；
  - 已完成 trace0 局部守恒兼容性检查，见 `research-notes/176-trace0局部守恒兼容性.md`；
  - 守恒检查结论是干涉中点局部 trace0 几乎自然守恒，但分离态自然散度很大，因此 trace0 必须和 full conservation 联立求解；
  - 已实现 trace0 硬消元后的 sparse/global 守恒 penalty 扫描，见 `research-notes/177-trace0稀疏全局守恒扫描.md`；
  - 三切片扫描结论是 trace0 + full conservation 兼容，`force_weight=10` 可把两个分离态 full divergence 压到 `~0.0065-0.0077`，代数残差只轻微增加；
  - 已把 sparse penalty 扫描升级为 hard-project、KKT 和 ALM 原型，见 `research-notes/178-trace0硬守恒KKT与ALM原型.md`；
  - 当前最稳口径是 `hard-alm`：三切片 full divergence 进一步压到 `~2e-5` 到 `~2e-3`，中心代数残差只小幅增加；
  - `hard-kkt` 在小规模有效，但当前 LSQR saddle-point 版本在左侧分离态病态，后续若追求严格 KKT 需 MINRES/Schur complement/预条件器；
  - 已整理 equation-first trace0 proposal v1，见 `research-notes/179-equation-first-trace0-proposal-v1.md`；
  - 已细化 \(Q\to0\Rightarrow C\to0\) 分支和约束传播逻辑，见 `research-notes/180-Qto0分支与约束传播.md`；
  - 下一步理论任务若继续 action 化，应转向真实辅助场 \(\chi\)-sector 的有效应力；否则应正式降级 action 化并继续 equation-first；
  - 下一步理论写作应检查 branch preservation：\(\lambda_I=\chi(\mathcal q)\hat\lambda_I\) 中 \(\hat\lambda_I\) 的有界性是否随演化保持；
  - 下一步还应写 patch transition：在 \(\Delta=0,w^2=0,r^\perp=0\) 退化面附近寻找替代张量基底，避免把投影坐标奇点误当物理奇点；
  - 已完成 branch preservation 与 patch transition 首轮分析，见 `research-notes/181-branch-preservation与patch-transition首轮.md`；
  - 已完成 branch / patch 数值检查，见 `research-notes/182-branch-patch数值检查.md` 和 `visualizations/trace0_branch_patch_check_n384_core10/`；
  - 数值检查结论是 trace0 + full conservation + hard-ALM 不自动保证 \(Q\to0\Rightarrow C\to0\)，必须把 \(C=\chi(\mathcal q)\hat C\)、\(\hat C\) 有界写成显式 branch 正则性；
  - 低维 \(E_u=\mathrm{span}\{\tilde g,uu\}\) 对分离态较好，但干涉中心仍需要完整四方向 patch，不能用 \(E_u\) 全局替代；
  - 已完成 branch 正则化 hard-ALM 首轮，见 `research-notes/183-branch正则化hard-ALM首轮.md`；
  - branch-v1.1 在 \(\tau=0\) 有可用窗口，但右侧分离态显示 \(\hat C\) 有界性与守恒存在张力；
  - 已完成完整 D 支 Cauchy 闭合快查，见 `research-notes/184-完整D支Cauchy闭合快查.md`；
  - 快查结论是物质 HJ + continuity 加上当前几何候选后，在主非退化 massive patch 内方程主部/计数已经接近闭合；下一步不是再加场方程，而是做 gauge-fixed 约束传播和 patch transition；
  - 下一步具体任务应把 \(\hat C\) 有界性从全局 L2 penalty 改成局部不等式/投影条件，并诊断右侧分离态 divergence outliers 的位置；
  - 另一个下一步理论任务是把 harmonic gauge 或 ADM gauge 下的完整演化/约束分裂写成文章级定理草案；
  - 下一步数值验证若继续，应使用 ALM 作为 hard conservation baseline，并只把 KKT 当作需要更好线性代数的严格求解方向；
  - 数值器只保留为验证工具：检验投影残差、守恒、近退化点和三切片高斯干涉反例。
- 数值注意：
  - 后续所有硬守恒/硬测地线投影默认用 nullspace 或等价约束保持算法；
  - 裸 KKT 只能做快速诊断；
  - 若要做全演化，每步投影校正也必须保持 \(\tilde\nabla^\mu\mathcal C_{\mu\nu}=0\)，不能用软 penalty 假装满足。
