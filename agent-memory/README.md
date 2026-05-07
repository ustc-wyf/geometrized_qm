# 项目记忆说明

## 项目概览

- 历史项目根目录（旧机）：`E:\量子势几何化_codex`
- 当前项目根目录（本机恢复）：`/Users/wangyunfei/Desktop/量子势几何化_codex`
- 路径说明：旧日志、旧 handoff 和研究笔记里若出现 Windows 路径，应视为与当前 macOS 路径指向同一份跨机器复制过来的项目
- 课题名称：量子势几何化
- 核心目标：尝试把 Klein-Gordon 标量场在 Bohm/Madelung 表述中的量子势项吸收到一个新的几何结构中，使粒子运动可在新度规下理解为测地线运动
- 当前研究范围：单个 Klein-Gordon 方程、自由标量场、背景时空与标量场的相互作用，不进入多体问题
- 当前约定：自然单位制 `\hbar=c=1`，闵氏度规 `diag(1,-1,-1,-1)`，Madelung 分解 `\psi=\sqrt{\rho}e^{iS}`
- 当前计算记号约定：`X=g^{\mu\nu}\partial_\mu S \partial_\nu S`，`Y=g^{\mu\nu}\partial_\mu S \partial_\nu \sqrt{\rho}`，`Z=g^{\mu\nu}\partial_\mu \sqrt{\rho} \partial_\nu \sqrt{\rho}`
- 长记忆用途：把项目结构、方法、注意事项等信息持续写入此目录，供后续轮回重载

## 长记忆协议

- 状态源：长记忆是否开启，以 `agent-memory/LOG.md` 中最近一次显式记录的 `ON/OFF` 为准
- 开关口令：支持“长记忆开/关”“进入长记忆”“agent memory on/off”“am on/off”等显式指令
- 开启行为：只要最新状态为 `ON`，每一轮正式工作前都必须先重读 `README.md`、`TASKS.md`、`LOG.md`、`DECISIONS.md`；必要时再补读 handoff、源码和研究笔记
- 重启判定：若当前 agent 不能可靠知道自己处于第几个轮回，则按 restart 处理；轮回号取 `LOG.md` 中上一次记录的轮回号加一
- 续接原则：跨机器、跨线程或重启式接管时，保留旧记忆原文，只追加新的轮回记录、路径映射和新的判断，不覆盖旧历史
- 写入分工：具体事件、时间线、失败尝试写入 `LOG.md`；抽象原则和可复用判断写入 `DECISIONS.md`

## 架构与分层

- 理论主线 1：从 KG 作用量出发做 Madelung 分解，明确量子势项和 Bohm 力学中的运动学对象
- 理论主线 2：构造 Bekenstein disformal transform 一类的度规变换，研究能否消去显含量子势能项
- 理论主线 3：检查变换后的作用量、场方程、轨迹方程与实验可观测量是否保持物理等价
- 理论主线 4：研究经典极限下新旧度规是否渐近一致，以及将来是否能与 Einstein-Hilbert 动力学衔接

## 重要文件

- `agent-memory/HANDOFF-2026-04-27.md`：给新 agent 的交接手册，概括当前理论结论、solver 状态、哪些结果可信、建议的下一步

- `agent-memory/README.md`：项目核心设定、主线、注意事项
- `agent-memory/TASKS.md`：树状目标、当前活跃分支、推进顺序
- `agent-memory/LOG.md`：轮回、日志、阶段输出
- `agent-memory/DECISIONS.md`：抽象经验、方法选择、风险判断
- `/Users/wangyunfei/.codex/skills/agent-memory/SKILL.md`：本机新增的自定义 `agent-memory` skill 定义，供未来会话复用这套长期记忆协议
- `research-notes/README.md`：人工查阅用的阶段成果目录
- `research-notes/001-KG作用量的Madelung分解与量子势定位.md`：第一份正式推导留档
- `research-notes/002-量子势几何化的成功失败与可接受判据.md`：后续筛选 ansatz 的统一判据
- `research-notes/003-disformal-ansatz候选与第一轮筛选.md`：候选度规族与当前主攻路线的筛选结果
- `research-notes/004-四类度规变换的初步比较分析.md`：共形、单向与混合路线的统一初步比较
- `research-notes/005-四类变换代回作用量与HJ方程的直接计算.md`：四类方案的直接代回计算结果
- `research-notes/006-HJ层面的吸收Q求解与更广泛存在性.md`：HJ 层面的显式吸收解与存在性结论
- `research-notes/007-快速筛选：作用量兼容性与存在条件.md`：最简单 HJ 解的作用量筛选与存在条件
- `research-notes/008-方法论修正：作用量等价不是必要条件.md`：把筛选标准从“作用量逐项等价”修正为“物理方程等价优先”
- `research-notes/009-变分筛选：最简单共形与类II方案.md`：直接检查最简单共形与类 II 方案的变分表现
- `research-notes/010-密度随度规变换后的重新检验.md`：在 `sqrt(-g)ρ = sqrt(-g~)ρ~` 约定下重算筛选结果
- `research-notes/011-方法修正后：先变分再代回变换.md`：按正确顺序重检共形与类 II，撤回类 II 的错误排除
- `research-notes/012-一般类II能否同时满足HJ与连续性方程.md`：一般类 II 在当前作用量 ansatz 下的结构性筛选结果
- `research-notes/013-类II中A=1时B是否必然含Z的复核.md`：复核 `A=1` 时 `B` 是否被方程强制含 `Z`
- `research-notes/014-按逆度规ansatz直接消去拉氏量中的Q项.md`：按用户指定方式直接在拉氏量中消去 `Q` 的计算
- `research-notes/015-最简单拉氏量消项解的后续检验.md`：检验最简单消项解对 HJ、连续性、测地线与可逆性的影响
- `research-notes/016-连续性方程中代入X=m2-Q后的精确比较.md`：把 `X=m^2-Q` 代入后精确比较新旧连续性方程
- `research-notes/017-J^mu∂_muQ与新几何测地线方程是否等价.md`：比较守恒条件与新几何测地线条件的关系
- `research-notes/018-最小混合逆度规下的连续性方程no-go.md`：最小混合振幅修正下无法同时完成拉氏量消项与原连续性方程重合
- `research-notes/019-去掉sqrt-g-rho不变假设后的纯u-ansatz翻盘.md`：放开 `ρ~` 变换后，纯 `u^mu u^nu` ansatz 同时重现 HJ 与连续性方程
- `research-notes/020-最简单解的可逆性签名因果结构与经典极限.md`：最简单完整候选在 `X>0` 区域内的几何与物理性质
- `research-notes/021-X正负对应的纯ansatz双分支结构.md`：`X>0` 与 `X<0` 对应两条不同的纯几何分支
- `research-notes/022-对Cursor对照分析的判断.md`：对 Cursor 对照意见的阶段性判断
- `research-notes/023-阶段性总结与Einstein作用量耦合框架.md`：当前阶段的正式总结，以及 `S_EH[g]` / `S_EH[g~]` 两个引力候选框架
- `research-notes/024-S_EH[g]与S_EH[g~]的弱展开与差异量级.md`：比较两种 Einstein 作用量在弱展开与强量子势场景下的差异
- `research-notes/025-1+1d自由KG双高斯波包基线模型.md`：1+1d 双高斯自由 KG 干涉模型的定义与模块化计算说明
- `research-notes/026-1+1d双高斯KG基线算例结果.md`：该模型在平直背景下的第一轮 `Q` 量级分析
- `research-notes/027-用双高斯KG例子检验四种一致性机制.md`：用具体双高斯模型检验四种“恢复实验预测”的思路
- `research-notes/028-Q依赖引力作用量的一般要求与例子分析.md`：分析 `F(Q)R~` / `ρ~R~` 类引力作用量的可行性与优先级
- `research-notes/029-双高斯模型上测试FQ=1-1pluskQ2.md`：在具体双高斯模型上测试 `F(Q)=1/(1+(κQ)^2)` 的大小与导数风险
- `research-notes/030-1+1d平直背景下EHg与EHgtilde的真正差异.md`：说明 1+1d 里纯 EH 比较的拓扑限制
- `research-notes/031-1+1d双高斯模型下两种总作用量的直接比较结果.md`：给出该模型下 `EH[g]` 与 `EH[g~]` 的直接数值比较
- `research-notes/032-3+1d自由KG双缝型双高斯波包模型.md`：第一版 3+1d 双缝/屏幕模型（后续已不作为主算例）
- `research-notes/033-3+1d对射双高斯波包模型.md`：对射高斯包几何设定（后续主要作为几何草稿参考）
- `research-notes/034-3+1d准单色双缝束流模型.md`：更贴近实验直觉的准单色双缝束流模型说明
- `research-notes/035-3+1d自由KG精确Fourier对射模型.md`：严格自由 KG 的 3+1d Fourier 表示与分块积分策略
- `research-notes/036-在当前3+1d模型上比较三种总作用量的差异入口.md`：说明三种总作用量在当前层级下真正可比较的“差异入口”
- `research-notes/037-三种总作用量在共同初始数据上的差异入口图.md`：共同初始数据上的入口量图示框架
- `research-notes/038-三种总作用量的最小闭合原型演化.md`：reduced closure 原型（已判定不可信，不作主结论）
- `research-notes/039-按F(Rtilde)Rtilde的三条reduced动力学方程.md`：`F(R~)R~` 的 reduced 原型（同样不作主结论）
- `research-notes/040-三种完整作用量的统一约束形式.md`：A/B/C 三种完整作用量的统一写法
- `research-notes/041-fR分支的辅助场重写.md`：C 分支辅助场重写
- `research-notes/042-三种作用量的3+1d到2+1d一致约化框架.md`：从 3+1d 到 2+1d 的一致约化框架
- `research-notes/044-2+1d下三分支动力学方程与数值闭合.md`：`2+1d` conformal-Killing 约化下的三分支动力学方程、B/C 的显式 PDE，以及 A 分支仍需额外 reduced closure 的说明
- `research-notes/072-三种引力作用量在2+1d下的全动力学方程与当前数值状态.md`：当前 A/B/C 引力作用量主线的统一索引，明确区分形式全方程、显式闭合 PDE 与当前数值可信层级
- `research-notes/055-共享baseline下A-B-C偏离的更密lambda扫描.md`：共享 baseline 上更密 `lambda_grav` continuation 的阶段结果与代表区间
- `research-notes/056-同初值下2+1d高斯波包干涉的ABC后续演化比较.md`：同一初值下 A/B/C 后续演化的受控动态图样比较，当前最直接回答“引力作用量差异如何显影”的阶段文档
- `research-notes/057-lambda1下AB可测量实验量与弱曲背景估计.md`：把 `lambda=1` 的 A/B 差异压成条纹间距、等效相位移、可见度和亮纹强度，并给出地球表面 Schwarzschild 与纳赫兹 GW 背景下的弱场取向估计
- `research-notes/058-lambda1下AB三层补完：收敛、时间序列与背景场景.md`：对 `lambda=1` 的 A/B 候选 observable 做 dt/grid 收敛性检查，补时间序列/积分型 readout，并把弱曲背景升级成场景梯度 + 取向扫描
- `research-notes/059-lambda1下AB的两类实验差分量：相位翻转与干涉超额.md`：把 A/B 信号改写成“相位翻转差分”和“干涉超额计数”两类双设置实验量，并比较它们各自的 A/B 分裂与最小 dt 收敛性
- `research-notes/060-lambda1下AB实验通道的噪声基准比较.md`：把 raw 强度通道、phase-flip differential 和 interference excess 放进统一的 shot-noise / 漂移 / 归一化误差基准里比较，并给出当前最实际的实验优先级
- `research-notes/062-混合投影型逆度规解的统一整理.md`：把用户提出的 `D=Z, E=-Y, F=X` 混合逆度规解整理成统一框架，并并排总结 `m≠0, C=1` 与 `m=0, C=1, κ≠0` 两条分支
- `research-notes/063-度规退化、不可逆变换与真奇点的文献分层.md`：整理 GR、disformal/Bekenstein、signature change 与 mimetic 语境里“坐标奇点 / 变量奇点 / 真曲率奇点 / 退化面”的分层，并给出对 `X=0`、`Q=0`、`Δ=0` 的项目内解释框架
- `research-notes/064-纯共形变换能否满足几何化要求.md`：重新用当前新约定和“先变分后代回”的口径检查纯共形变换，区分它在方程级的成功与在拉回旧变量作用量上的退化
- `research-notes/065-混合ansatz新解D等于Z-E等于负aY-F等于aX.md`：检查另一组混合逆度规解 `D=Z, E=-aY, F=aX, C=aXZ`，给出它恢复 HJ/连续性方程的条件、行列式比值和新增退化面
- `research-notes/066-一般混合逆度规ansatz的行列式恒正条件.md`：把一般 `C,D,E,F` 混合逆度规的行列式比值总公式和 strict positivity 条件单独整理出来，供后续所有具体解套用
- `research-notes/067-在X等于0处让rho波浪吸收奇异性的可能性.md`：系统检查“让 `ρ~` 在 `X=0` 趋于 0` 能否吸收奇异”，区分它对电流和对 HJ 主奇异的不同作用
- `research-notes/068-混合ansatz在X等于0附近的严格局域分析.md`：对一般混合 ansatz 做严格局域分析，证明在正则开区域里若同时恢复连续性方程与固定非零 HJ 壳长，则主组合必然是 `κ/X`
- `research-notes/074-新混合变换下的一般2加1ADM全动力学方法.md`：建立适配当前混合变换的一般 `2+1 ADM` 全动力学方法，放弃旧的对角 conformal gauge，并先把共同物质推进器重写成 `(n,S)` 的守恒系统
- `research-notes/075-B与C两支在一般2加1ADM度规下的引力动力学.md`：把 `B/C` 两支重新写到一般 `2+1 ADM` 语言下，明确它们分别对应 `2+1 Einstein + β + Bohm 物质` 与 `2+1` 标量-张量系统
- `research-notes/076-A支应作为g上Einstein-KG系统而非g与g波浪混合基本作用量.md`：正式纠正 `A` 支定义；`A` 支的基本总作用量应全部写在 `g` 上，本质是 `g` 上的 Einstein-Klein-Gordon 系统，`g~` 只作为由 `(g,\rho,S)` 事后重构出来的有效几何
- `research-notes/077-三支理论在y方向Killing对称下的统一四维ADM初值形式.md`：把 `A/B/C` 三支重新统一到四维 `ADM` 初值问题语言下；`A` 支用 `g` 上的 Einstein-KG，`B` 支用 `g~` 上 Einstein-Bohm，`C` 支用 `g~` 上辅助场 `f(R~)` 标量-张量系统
- `research-notes/078-A与B两支在统一四维ADM下的最小同步规范原型.md`：记录新的最小全动力学数值启动；在固定 `N=1, N^x=N^z=0` 的同步规范下，`A/B` 两支已经进入同一原型时间推进
- `research-notes/079-C支已接入统一原型但约束残差仍显著偏大.md`：记录 `C` 支已经接入同一同步规范原型，但其 Hamilton 约束残差当前明显高于 `A/B`，说明 `\nabla_\mu\nabla_\nu\Phi` 的 `ADM` 投影或有效源项重排仍需校正
- `kg_examples/simulate_2p1_full_dynamics.py`：`2+1d` 全动力学数值器，当前实际跑的是 `flat reference / EH[g~] / F(R~)R~` 三支
- `kg_examples/simulate_2p1_boundary_driven.py`：`2+1d` 边界持续入射版数值器，A 分支改为原始平直 KG，B/C 分别为 `EH[g~]` 与 `F(R~)R~`
- `kg_examples/check_benchmark_residual.py`：检查当前 `beamlike` benchmark 对平直 KG/Helmholtz 的残差
- `kg_examples/exact_kg_2p1_benchmark.py`：构造真正的 `2+1d` 正频 free KG 精确 benchmark，并输出 `alpha=0,0.5,1` 的 `x-z` 截面
- `kg_examples/calibrate_exact_kg_solver_2p1.py`：用精确 `2+1d` free KG benchmark 去校准 A 分支时间推进器的误差
- `kg_examples/solve_a_helmholtz_2p1.py`：A 分支的 Helmholtz 边界值尝试（已证实不适合当前非单频 benchmark）
- `kg_examples/calibrate_exact_kg_fft_2p1.py`：A 分支的 FFT 正频传播校准器，当前这是最可信的 A 数值基线
- `kg_examples/calibrate_exact_kg_fd_absorbing.py`：A 分支的“大域 + 吸收层 + PDE 时间推进”校准器
- `kg_examples/calibrate_exact_kg_ps_fd_absorbing.py`：A 分支的“谱 Laplacian + 吸收层 + 精确边界”校准器，当前最好的 PDE 基线
- `kg_examples/compare_abc_same_initial_evolution.py`：同一初值下的 A/B/C 动态比较脚本；A 用可信 PDE 骨架演化，B/C 在同一 A 的 `rho,S` 历史上推进几何并重构 `rho`
- `kg_examples/analyze_ab_measurable_observables.py`：只跑 A/B 的实验量提取脚本；当前默认已切到 refined 口径（`dt_target=0.0125, nx=nz=181, nkx=nkz=81`），输出 `lambda=1` 下的条纹间距、等效相位移、可见度、亮纹强度差，以及 Earth Schwarzschild / 纳赫兹 GW 背景的弱场取向估计
- `kg_examples/scan_ab_observable_convergence.py`：对 `lambda=1` 的 A/B 候选 observable 直接做时间步和网格收敛性扫描；当前已确认中心亮纹强度差是最稳的主 observable
- `kg_examples/analyze_ab_time_observables.py`：把 A/B 差异扩展到时间序列和积分型 readout，跟踪差异何时在束流重叠过程中长出来
- `kg_examples/compare_ab_background_scenarios.py`：把弱曲背景从单点估计升级成 Earth surface / LEO / GEO / Sun-Earth L1 / 近源质量 / 纳赫兹 GW 的场景梯度与取向扫描，并直接对比当前数值地板
- `kg_examples/analyze_ab_differential_experiment_observables.py`：把 A/B 差异改写成实验双设置差分量，当前已实现 `D_phi` 与 `E_int`，并对 `dt=0.025, 0.0125, 0.00625` 做最小时间步核对
- `kg_examples/analyze_ab_noise_baselines.py`：把当前候选实验通道放进统一的 shot-noise 与系统误差模型，输出所需计数、pair imbalance / normalization / drift 容忍度和 SNR 曲线
- `kg_examples/adm_matter_2p1.py`：新方法的第一块可运行代码；在一般 `2+1 ADM` 度规下实现共同物质作用量对应的 `(n,S)` 守恒演化与 `ρ` 回收
- `kg_examples/adm_ykilling_geometry.py`：`y` 方向 Killing 对称下的三维空间几何工具，包括二维度规逆、二维 Ricci 标量、`β` 的 Hessian/Laplacian 以及三维空间 Ricci 量
- `kg_examples/adm_sources_abc.py`：`A/B` 两支的应力能量投影工具；支持从 `(\phi_1,\phi_2,\Pi_1,\Pi_2)` 或 `(\rho,S)` 计算 `ADM` 源项
- `kg_examples/simulate_ab_adm_synchronous.py`：新的最小同步规范全动力学原型；已能在统一四维 `ADM` 变量下共同推进 `A/B` 两支
- `kg_examples/simulate_ab_adm_synchronous.py`：现已扩展到 `A/B/C` 三支共同原型；`C` 支可以短时间推进，但约束残差仍待压低
- `kg_examples/fit_equation_first_constrained_bcd.py`：当前 equation-first gBCD 主诊断脚本；支持 `--atoms g,uu,rr,ur`、`--hard-constraint`、`--force-mode transverse/full`，用于检查代数场方程残差与守恒/测地线约束
- `kg_examples/diagnose_gbcd_constitutive_coefficients.py`：读取 hard-constrained gBCD 的 \(A,B,C,D\) 系数场，检查它们能否由简单局部标量函数预测；当前初检结果是否定的
- `kg_examples/fit_gbcd_potential_constitutive.py`：检验单个局域势函数 \(L(\rho,u^2,r^2,u\cdot r)\) 经 metric variation 产生 gBCD 系数的机制；degree 2/3 当前失败
- `kg_examples/fit_gbcd_auxiliary_gauge.py`：gBCD 辅助各向异性应力场的代表元规范诊断；支持 nullspace hard constraint，当前用于比较 `norm/time/space` 等代表元选择；已加入 `--q-gated-norm` 用于把 \(Q\to0\) 回 Einstein 写成辅助场范数权重
- `kg_examples/diagnose_gbcd_conservation_principal_symbol.py`：检查 gBCD full-conservation 方程对 \(\lambda_I=(A,B,C,D)\) 的一阶主符号；当前确认 `2+1` 约化下 rank=3、nullity=1，需要额外状态方程
- `kg_examples/fit_gbcd_trace_closure.py`：在 gBCD/full-conservation/nullspace 上测试额外代数状态方程；当前 `trace=F(Q)` 与 `D=0` 均不理想
- `research-notes/146-equation-first-gBCD闭合方程与守恒检验.md`：gBCD 方程、守恒展开式、三时刻 hard constraint 检验和本构闭合问题的最新理论索引
- `research-notes/147-gBCD系数ABCD的可能产生机制.md`：整理 gBCD 系数的三类产生机制：简单局部本构函数、单势函数 metric variation、辅助各向异性应力场；当前支持第三类
- `research-notes/148-gBCD辅助应力场代表元规范与nullspace硬约束.md`：记录 KKT hard constraint 的数值泄漏问题、nullspace 修正、三切片 `norm/time` 代表元结果与下一步闭合任务
- `research-notes/149-gBCD状态方程trace与剪切闭合首轮检验.md`：记录 `trace=0`、`trace=aQ+bQ^2` 与 `D=0` 的首轮状态方程检验；结论是简单代数闭合会恶化残差
- `research-notes/150-Q趋零条件作为辅助场变分权重.md`：记录 \(Q\to0\) 回 Einstein 条件的正确实现方式：作为辅助场 action 的权重/边界条件，而非 trace=\(F(Q)\)
- `research-notes/151-gBCD辅助场变分闭合的Euler-Lagrange方程.md`：推导固定参考切片上的 saddle-point 辅助场闭合泛函和 Euler-Lagrange 方程，连接数值正则项与连续作用量
- `research-notes/152-无Q门控的最小辅助场action扫描.md`：按用户要求暂不考虑 \(Q\)，扫描 norm/time/space 最小辅助 action；当前默认候选为 `norm_weight=1e-8,time_weight=1e-5,space_weight=0`
- `research-notes/153-gBCD显式理论候选v0.md`：把当前结果整理成 equation-first 显式理论候选 v0：主场方程、守恒、辅助场最小泛函和受限变分方程
- `research-notes/164-gBCD投影型Einstein-like方程候选.md`：当前 equation-first 主线的显式方程候选；把 Einstein 残差投影到 \(\{\tilde g,uu,rr,ur\}\) 张量子空间，并要求投影修正张量守恒和 \(Q\to0\) 回 Einstein 分支
- `research-notes/165-gBCD投影方程闭合条件.md`：投影方程的首轮理论闭合检查；说明 \(W\) 不改变在壳成员关系，给出 \(\det H=\frac{d-2}{2}(u^2r^2-(u\cdot r)^2)^3\)，并明确 \(1+1d\)、\(\Delta=0\)、守恒/测地线和 Helmholtz 问题
- `research-notes/166-gBCD投影方程的3加1计数与主部.md`：投影方程的 3+1d 方程计数与主部检查；说明 6 个投影方程和 4 个坐标规范自由度的关系，以及 harmonic gauge 下 \(\Pi_E^\perp\delta G\) 主部
- `research-notes/167-gBCD投影方程的主符号缺口与标量闭合条件.md`：投影方程加 harmonic gauge 后仍有一个 \(E\)-方向主部 null mode；给出 \(N^I=(0,-q^2,-p^2,2pq)\)，并提出 \(u-r\) 平面正定迹作为首个标量闭合候选
- `research-notes/168-qEtrace闭合三切片检验.md`：\(\mathsf q_E\)-trace 标量闭合的三切片检验；当前最小硬闭合版本不如普通 trace=0，下一步应检查 trace=0 的 \(w^2=0\) 退化面或把 \(\mathsf q_E\) 改作辅助场范数
- `research-notes/159-n384三切片gBCD求解器升级与物理主线判断.md`：记录 `n=384,core10,tau=-3.5,0,+3.5` 三切片高分辨率 gBCD 求解器升级、残差结果和下一步物理主线判断
- `kg_examples/plot_gbcd_metric_update_diagnostics.py`：读取 full-linear metric update 的 `.npz`，生成包含 `rho`、残差、\(\delta g_+\)、`det(corrected g+)` 的诊断图，并在图内说明 white contour 与 `core10` 定义
- `visualizations/full_dynamics_2p1/`：`2+1d` 全动力学数值结果图和 `summary.json`
- `visualizations/boundary_driven_2p1/`：边界驱动版 A/B/C 三支结果
- `visualizations/abc_same_initial_evolution_lambda1/`：`lambda=1` 的同初值 A/B/C 动态比较图与汇总
- `visualizations/abc_same_initial_evolution_lambda3/`：`lambda=3` 的同初值 A/B/C 动态比较图与汇总
- `visualizations/ab_measurable_observables_lambda1/`：`lambda=1` 下 A/B 的中心条纹实验量提取结果与弱曲背景估计
- `visualizations/ab_observable_convergence/`：`lambda=1` 下 A/B observable 的 dt/grid 收敛性图与汇总
- `visualizations/ab_time_observables_lambda1/`：A/B 差异的时间序列、积分型 readout 曲线和相对差 montage
- `visualizations/ab_background_scenarios_lambda1/`：弱曲背景场景梯度、取向扫描与“背景诱导改变量 vs 数值地板”的汇总
- `visualizations/ab_differential_experiment_observables_lambda1/`：refined 基线下两类实验差分量的时间序列与汇总
- `visualizations/ab_differential_experiment_observables_lambda1_dt025/`：较粗时间步核对
- `visualizations/ab_differential_experiment_observables_lambda1_dt00625/`：较细时间步核对
- `visualizations/ab_noise_baselines_lambda1/`：各实验通道的 required counts 与 shot-noise SNR 曲线
- `presentation-workspaces/quantum-potential-geometrization-ab/`：当前完整技术汇报 deck 的 workspace
  - `src/deck.mjs`：PPT 源码
  - `output/output.pptx`：最终交付版
  - `scratch/previews/`：逐页 PNG 预览
  - `scratch/quality-report.json`：PPT package QA 结果

## 常用方法

- 先固定最小模型，再逐步增加结构，避免一开始同时处理多体、量子场论、引力动力学
- 优先从作用量层面而非仅从运动方程层面构造理论，这样更容易检查变换前后物理等价性
- 每一步都要区分三件事：数学上可改写、物理上可解释、实验上可等价
- 筛选时优先检查变分得到的物理方程是否等价；作用量逐项等价只算强充分条件，不算必要条件
- 每得到一个新公式，都要检查经典极限、因果结构、可逆性与符号问题
- 每个相对完整的小步骤都单独留档到 `research-notes/`，供人工查阅
- 当前主攻 ansatz：同时包含 `u_\mu=\partial_\mu S` 与 `r_\mu=\partial_\mu \sqrt{\rho}` 的混合 disformal 子类
- 当前展示习惯：优先使用纯文本公式，减少窗口中的 LaTeX 渲染负担
- 当前新增主线：优先采用“把逆度规 ansatz 直接代回拉氏量消去 `Q` 项”的检验方式
- 当前新增主线：允许 `ρ~` 真正作为变换对象后，重新评估纯 `u^mu u^nu` ansatz 的完整可行性
- 当前新增主线：以最简单纯 ansatz 为基底，比较 `S_EH[g]` 与 `S_EH[g~]` 的弱差异
- 当前新增主线：在混合 `u-r` ansatz 中，`D=Z, E=-Y, F=X` 这组系数可被重写成作用在 `span{u,r}` 上的投影型混合块；当前最自然的下一步是单独分析其退化面 `X=0`、`Q=0`、`Δ=XZ-Y^2=0`
- 当前新增主线：在 `2+1d` conformal-Killing 约化下，把 `R[g~]` 和 `F(R~)R~` 写成显式波动系统再数值积分；`R[g]` 分支若要与之同级比较，还需要额外指定 `g↔g~` 的 reduced closure
- 当前新增主线：若目标是重现用户认可的 `±45°` 相干束干涉图样，数值问题应优先写成边界持续入射问题，而不是自由 Cauchy 衰减问题
- 当前新增主线：当前 `beamlike` 图样不是精确 flat KG 解；若坚持把它当目标实验图样，则三分支比较应转向稳态边界值问题，而不是要求 A 分支时间推进器贴近它
- 当前新增主线：已经建立 `2+1d` 精确 free KG benchmark；后续若要让 A 数值器“贴基准”，应贴这个精确 benchmark，而不是贴 beamlike 近似图样
- 当前新增主线：已把精确 benchmark 的观察窗扩大到 `[-10,10]` 并把传播轴动量展宽压到 `sigma_parallel/k0 = 0.01`；当前 A 分支时间推进器的快速误差探针显示仍需进一步校准
- 当前新增主线：A 分支的 `leapfrog + 全边界 Dirichlet` 校准测试已经表明该方法不收敛；后续若要让 A 真正贴近精确 benchmark，需要改边界处理或直接转向频域/Helmholtz 型求解
- 当前新增主线：A 分支的 Helmholtz 尝试已判定不适合当前非单频精确 benchmark；FFT 正频传播校准则在中心干涉区贴合良好，当前最适合作为 A 的数值基线
- 当前新增主线：A 分支的 PDE 校准已推进到“谱 Laplacian + 吸收层 + 精确边界”版本；这版在内核区把相对 `L1` 压到约 `0.22%`，相对 `L∞` 压到约 `4.7%`
- 当前新增主线：同初值 A/B/C 图样比较已经有了一个受控“共享主干”数值支架；当前做法是让 A 可信演化、把 A 的 `rho,S` 历史作为共享物质演化主干，再让 B/C 在这同一组历史上推进几何并重构 `rho`，用来隔离“引力作用量改变后几何响应会怎样变”
- 术语约定：`共享主干` 指“多分支共享同一条物质演化历史”，`数值支架` 指“为检验某个理论问题而搭建的受控计算框架”；二者都不等于“全反作用自洽演化”
- 当前新增主线：`lambda=1` 的 A/B observable 三层补完已经完成；当前最可信的主 observable 是亮纹/热点强度约 `6.7e-05` 的系统性重标定，而不是条纹间距或条纹相位的明显平移
- 当前新增主线：中心亮纹强度差 `central_peak_rel_intensity_shift` 的 refinement 后数值为 `~6.76e-05`，其数值地板约 `3.6e-09`，信号高出地板约四个数量级
- 当前新增主线：积分型 readout（中心圆盘、亮区总强度、干涉盒总强度）在最终时刻都落在同一 `6.6e-05 ~ 6.8e-05` 量级，并且亮区积分信号比单个中心峰更早开始响应
- 当前新增主线：当前的 Earth / L1 / GW / 近源质量比较已经升级成显式场景梯度与取向扫描；但它们仍只是真空弱场下的 tidal/readout modulation，不是 full curved-background 自洽动力学
- 当前新增主线：更实验化的“双设置差分 observable”也已经算完；`相位翻转差分` 的 A/B 分裂只有 `~1e-8`，`干涉超额计数` 的 A/B 分裂可到 `~1.3e-06`
- 当前新增主线：在这两类实验差分量里，`interference excess` 明显优于 `phase flip differential`，而 `interference box` 读出窗优于 `central disk`
- 当前新增主线：噪声基准比较已经完成；按当前最简单的 shot-noise / drift / normalization 模型，最实际的主通道不是差分量，而仍是 raw integrated intensity
- 当前新增主线：`phase flip` 在当前口径下统计和系统误差代价都过高，适合当概念性 cross-check，不适合当前主打
- 当前新增主线：已经整理出一份 14 页的完整技术汇报 deck，叙事顺序与当前论文主线一致，可直接作为后续写 paper 或精简组会版的骨架
- 当前主线纠偏：数值器不是目标；目标是写出可一般化的变换后 Einstein-like 方程。当前候选是投影型方程
  \[
  \Pi_E^\perp(\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2)=0,
  \]
  配合投影修正张量守恒、物质壳方程/连续性方程和 \(Q\to0\) 回 Einstein 分支条件。

## 注意事项

- 不要一开始就把目标扩大成完整量子引力理论；当前目标是最小模型中的“量子势几何化”
- 负粒子数密度与超光速问题是动机，不是自动已经被解释的结论；必须靠具体推导来检验
- `disformal transform` 的选取必须满足协变性、可逆性、物理解释性，不能只为消项而消项
- 若未来引入引力作用量，必须清楚区分“固定背景上的有效几何改写”和“真正动态度规理论”
- 当前 `2+1d` 数值器里的 `flat reference` 不是严格的 A 分支完整闭合，只是当前在未额外指定 `g↔g~` reduced closure 前，最干净的参考支
- 当前关于三种引力作用量的讨论必须始终区分三层：
  - 形式作用量与形式场方程
  - `2+1d` 约化后是否已有唯一显式闭合 PDE
  - 当前数值结果是否只是 shared-baseline，而非三支完全自洽全反作用
- 当前最可信的“同初值 A/B/C 后续演化比较”仍是受控 shared-baseline 版本，而不是 full backreacted B/C matter solver；若要把结论升级成真正三支完全自洽闭合比较，必须先重新验证 B/C 的自反馈物质求解器
- 当前新混合变换下，正测度密度必须写成 `sqrt(|g~|)ρ~ = (|X|/|κ|) sqrt(|g|)ρ`；带符号的 `X/κ` 只能作为 branch diagnostic 使用，不能直接当作物理正密度
- 因而在固定 `κ` 的强匹配口径下，只能在 `κ/X > 0` 的单个 patch 上工作；若 `X` 在主支撑区内部变号，就不能再机械地把旧 shared-baseline 源项写成 `(X/κ)ρ_A`
- 当前更关键的新判断是：对 `±45°` 高斯波包初值，当前混合变换会 generically 生成 `t-x`、`t-z`、`x-z` 非对角分量，所以旧的对角 `τ,β` 全动力学数值器已经不对应当前模型；后续必须改用一般 `2+1 ADM` 度规
- 当前 `B/C` 两支已经重新写成一般 `ADM` 初值问题；但 `A` 支的定义也已进一步澄清：它不是“引力在 `g` 上、物质在 `g~` 上”的混合作用量，而是 `g` 上的 Einstein-Klein-Gordon 系统
- 作图规范新增硬规则：任何图中出现的轮廓线、掩膜、窗口名、变量符号，都必须在图注、图内文字或配套摘要中明确给出数学定义；不能默认用户知道白线/青线/热点/`chi`/`\delta chi`/`\mathcal T` 等符号含义。
- 因此当前 A 分支的真正未完成项，不再是“把混合作用量闭合”，而是：在一般 `2+1 ADM` 变量下，把 `g` 上的 Einstein-KG 系统写成可数值推进的约束—演化方程，并在每一步演化后由当前混合变换事后重构 `g~`
- 当前方法已再次收束：不再依赖 `2+1` 共形重参数化来组织三支，而是统一改用四维 `ADM` 加 `y` 方向 Killing 对称；这样 `A/B/C` 的比较语言更直接，也避免约化共形因子带来的歧义
- 当前数值阶段已真正开始：`A/B/C` 三支都已进入同一最小同步规范原型；其中 `A/B` 的短时间行为较一致，而 `C` 支虽然能推进，但 Hamilton 约束残差显著偏大，下一步应优先校正 `C` 支的 `ADM` 投影
- 当前 `C` 支校准已进一步完成：在统一四维 `ADM` 原型中，只要取 `\Phi=1,\Pi_\Phi=0`，`C` 支的有效能量、动量、空间应力、`r3` 历史和 Hamilton 约束残差都会逐项、逐步严格退回 `B` 支；这说明 `C` 支的平直辅助场归一化和有效 Einstein 源的平直扣除已经对齐
- 当前因此不再把 `C` 支的主要问题理解成“源项投影整体错误”，而是收缩成两个真正剩余问题：
  - 如何为 `C` 支构造非平凡而一致的 `\Phi` 初值
  - 如何为三支共同求解满足约束的初始几何，而不是从“平直几何 + 非零物质”的非约束态硬启动
- 当前又新增一层更精确的初值判断：在同一组 `\rho,S` 上，`A` 支与 `B` 支的约束源差异实际上很大；因此“三支完全同一组 `g` 初值”一般不能同时严格满足三支约束，只能在弱场下作为近似共同初值
- 当前第一版共同初始几何已建立：
  - 共同空间度规：时间对称、共形平坦、由 `(\mathcal E_A+\mathcal E_B)/2` 驱动
  - 共同外曲率：由 `(\mathcal P^{(A)}+\mathcal P^{(B)})/2` 的最小 York 型向量势近似给出
  - 这使 `t=0` 的三支 Hamilton 残差从此前 `O(1)` 量级降到 `~8e-2`
- 当前 `C` 支的第一条非平凡初值入口也已接入：
  \[
  \Phi_0 = (1+(\ell^2 R_{\text{proxy}})^2)^{-3/2},
  \]
  其中 `R_proxy` 取共同初始几何的空间曲率；但在当前 `M_P=300, \ell=0.02` 下，这个偏离仍极小，说明当前弱场参数下 `C` 支若要显著分开，还需要更强曲率或更有针对性的 `\Phi` 初值选择
- 当前方法论已再次调整到更贴近现实弱场的顺序：
  - `A` 支先作为平直背景上的参考支，用来提供初始 `(\rho,S)` 与弱场量级估计
  - `B/C` 支作为主要全动力学对象，应该从这组初始 `(\rho,S)` 出发各自完整耦合推进
  - `A` 支的引力回授先以线性化弱场估计给出误差量级，而不再与 `B/C` 同步重火力推进
- 当前这一新版方法已落成：
  `kg_examples/simulate_bc_geometry_from_a_reference.py`
  其中参考 `A` 支不再用旧的 `(\rho,S)` 近似推进器，而是直接在平直背景上推进复 Klein-Gordon 场
- 当前边界条件也已进一步分层：
  - 物质子系统 `n,S` 已从经验海绵层升级到特征入流/出流边界条件，入流数据来自同一高斯束准备下的解析参考解
  - 几何子系统在当前同步规范 `ADM` 下仍只使用弱场背景边界值回归，不能称为严格的约束保持边界条件
- 这一步的长时间数值检验已经完成：物质边界升级后，`B/C` 的弱场稳定性保留，且 `ell=0.2` 时 `C-B` 分离量级几乎不变；因此当前真正限制 `C` 支可分辨性的，已不再是物质边界，而是 `Phi` 初值过弱
- 当前第一轮结果表明：
  - `A` 支的线性化弱场几何偏离量级约为 `4e-3`
  - `B/C` 支在同一参考物质驱动下的几何偏离可达 `10^{-1}` 量级
  - 但 `B/C` 的几何更新仍有明显时间步敏感性，因此还不能直接做最终物理解读
- 当前又新增了更正后的完整耦合版本：
  `kg_examples/simulate_bc_from_a_initial_data.py`
  - 这版不再把整条 `A` 历史当作 `B/C` 的外部源
  - 而是只用 `A` 支提供初始 `(\rho,S)` 和初始弱场 `omega`
  - `B/C` 自己完整推进 `(n,S,g~)`；这里 `n = e^\beta\sqrt{h}\rho E` 作为真正的守恒密度变量
- 当前这一更正后的 `B/C` 完整耦合同步规范测试说明：
  - 初始 Hamilton 残差可压到 `~7.5e-2`
  - 但一进入演化，残差增长与几何偏离量级仍与上一轮几乎相同
  - 因而当前主要问题已经收缩到同步规范和 `B` 支约束传播本身，而不再主要归咎于参考 `A` 的使用方式
- 当前又新增一个更关键的数值判断：
  - 先前 `B/C` 原型里看似灾难性的几何增长，主因并不是时间推进阶数，也不是初始外曲率，而是**几何边界条件实现错误**
  - 错误做法是把 `h_xx,h_xz,h_zz,beta,k_xx,k_xz,k_zz,k_beta` 与物质变量一起做乘法阻尼；这会把几何边界从背景值硬压向 `0`
  - 修正后应采用：
    - `n,S` 做边界阻尼
    - `h_xx -> 1, h_xz -> 0, h_zz -> 1, beta -> 0`
    - `k_xx,k_xz,k_zz,k_beta -> 0`
    - `Phi -> 1, Pi_Phi -> 0`
  - 修正后 `B/C` 弱场演化在 `10/20/200` 步测试里都保持稳定，`B` 支几何偏离维持在 `~7.7e-3` 而不是此前伪造出来的 `10^{-1}` 量级
  - 因此当前 `B/C` 原型已经从“明显伪不稳定”进入“弱场稳定演化”阶段
- 当前 `C` 支非平凡初值的最新状态：
  - `r3_proxy` 入口已经重新接入修正后的稳定原型
  - 但在当前 `ell=0.02`、弱场高斯束参数下，`max|Phi-1|` 到 `200` 步时也只有 `~2.6e-8`
  - 因而 `C` 若要与 `B` 明显分开，下一步必须设计更强、但仍自洽的 `Phi` 初值，而不是继续沿用过弱的 `r3_proxy`
- 当前 `C` 支参数扫描的经验结论：
  - 在固定 `steps=200, dt=2.5e-4`、`r3_proxy` 初值下，`ell=0.02` 到 `0.2` 时 `C-B` 几何差异仍只在 `10^{-8}` 到 `10^{-6}` 量级
  - `ell` 提高到 `0.5` 左右后，`C-B` 差异才开始进入 `10^{-6}` 到 `10^{-5}` 量级
  - `ell=1.0` 时，`C-B` 差异已到 `10^{-4}` 量级，但约束残差差距也同步增大
  - 因此默认 `ell=0.02` 下 `C` 几乎退回 `B` 应被视作模型在弱曲率下的真实性质，而不是单纯数值器未分开
- 关于 `Phi` 初值的当前方法边界：
  - 纯局域代数条件不能生成更强的 `Phi` 初值，因为当前模型中 `0<Phi<=1` 时 `2u(Phi)-Phi u_Phi(Phi) <= 0`，而物质迹 `trace(T_B)` 为正
  - 因而下一步应改为椭圆型 `Phi` 初值方程，而不是继续寻找只依赖局域曲率标量的闭式公式
- 当前若对外展示数值结果，必须把“共享主干比较”和“全反作用自洽演化”明确分开；前者不能写成后者的完成版
- 当前对 Earth Schwarzschild 和纳赫兹 GW 背景的数字，只是建立在弱场真空近似上的“实验读出调制估计”；它们不是 full curved-background B solver 的最终答案
- 当前 refined 基线下，phase / fringe-shift 类型 observable 仍然不够干净：其信号量级没有显著高于空间离散导致的数值地板；后续实验讨论应优先围绕强度型 observable，而不是相位型 observable
- 当前 refined 基线下，双设置差分量虽然实验上更抗共模噪声，但会把理论分裂本身也部分抵消掉；后续实验设计应同时保留“最大理论分裂量级”和“最干净实验差分通道”这两个口径
- 当前最简单噪声基准下，`raw intensity` 约需 `5e9` 级计数即可到 `5 sigma`，而最优差分通道 `interference excess` 约需 `1.8e13` 总计数；后续若要翻盘差分路线，必须引入更具体的实验架构来证明它能极大降低现实系统误差
- 当前对图、表、PPT 的书写规范应统一为：正文以中文为主，必要术语写成“中文（English）”；所有画出的量都必须在图题、图注或正文第一次出现处给出明确定义，不能只写 `(B-A)/A` 这类未指明对象的缩写
- 当前写作规范再加一条：任何自创的工作名词、窗口名、区域名、掩膜名都必须在第一次出现处立即给出明确定义，不能先用后补；例如若使用“主支撑区”“核心区”“过渡层”这类词，必须同时写出其数学定义或数值阈值
- 当前 `Q` 记号已经切换到新约定：`Q = + \Box sqrt(ρ) / sqrt(ρ)`；若翻阅较早笔记，需注意很多地方使用的是旧约定 `Q_old = -Q`。统一迁移说明见 `research-notes/061-Q约定改为正BoxR除以R.md`
- 当前混合投影型解的现状是：在 `Δ ≠ 0` 的 patch 上，`m≠0, C=1` 可得到自动恢复经典极限的主分支；`m=0, C=1, κ≠0` 可得到行列式严格正的局域 patch 分支，但无经典恢复。真正尚未处理的是退化面 `X=0`、`Q=0`、`Δ=0`
- 当前对 `X=0`、`Q=0`、`Δ=0` 的默认方法论判断是：先把它们当作变换图册边界或非可逆面候选，而不是直接当作原 `g` 系的真物理奇点；只有在原不变量或 geodesic completeness 也失控时，才升级为“真奇点”判断
- 当前对纯共形变换的最新判断是：若允许 `ρ~` 一起变换并按“先变分后代回”检查，则 pure conformal 在方程级可恢复 HJ、连续性和测地线解释；但若把匹配条件直接拉回旧变量作用量，则 integrand 恒等塌成零，因此它不能同时充当一个非退化的 pulled-back 有效作用量
- 当前还需特别区分 `m≠0` 与 `m=0` 两类壳约束：前者在当前 rank-1 物质作用量里是正则分支；后者对应 `null-shell` 退化分支，不能简单看成把 massive 公式里 `m -> 0`。若后续认真处理 `m=0`，应优先改用带辅助场或拉格朗日乘子的作用量，而不是继续沿用 `∫ sqrt(-g~)ρ~(g~^(μν)u_μu_ν-κ)` 这一 massive-like 写法
- 当前机器为了补跑汇总图，新增了项目本地依赖目录 `.agent-deps/`（含 `matplotlib` 等）；它只是本地运行辅助，不是理论结果的一部分
- 若继续使用 `artifact-tool` 生成 PPT，图片资产不要只靠 `path/uri` 引用；应先拷到 deck workspace 的本地资产目录，再转成 `data URL` 以避免 `zero-byte media`
- 当前 `D` 支 bulk-interface 路线已修正为全域动态界面口径：
  - `kg_examples/prototype_d_bulk_interface_toy.py` 支持 `--bbox global` 与 `--levels 0.5,1,2`
  - 水平集连通分支属于拓扑问题，但选择多少个 `|y|` level 主要是过渡层源项的数值采样问题
  - `ell=30,t=16,n=240` 下，全域单界面 `|y|=1` 已覆盖 top 1% derivative hotspot 的 `100%`
  - 后续瓶颈是 subcell/body-fitted interface 表示，而不是继续增加 rasterized interface band 的 level 数
- 当前必须严格区分两种数值任务：
  - reference residual 诊断：把平直 QM 参考历史 \((\rho_{\rm ref},S_{\rm ref})\) 代入候选 native 方程，检查 residual；此时物质分布与平直 QM 差异按定义为零
  - full evolution：只用平直 QM 给初态/入流边界，随后由候选作用量自洽推进 \((\rho,S,\tilde g)\)，再和实验或平直 QM 预测比较
- 当前第一版 subcell interface toy evolution 已新增：
  - 脚本 `kg_examples/prototype_d_subcell_interface_evolution.py`
  - 研究笔记 `research-notes/115-D支subcell界面toy演化首轮.md`
  - 它用 `y=±1` 子网格线段作为界面源，不再把界面当厚带；当前能压低最极端 trace residual，但会抬高 p95，因此下一步必须与 bulk IMEX/Newton 校正和 saturated bulk 延拓耦合
- 当前已修正界面变量命名：界面变量必须写成 `raw_y = ell^2 R_tilde`；旧 `y_ref` 在部分脚本中实际是 `tanh(raw_y)`，不得再用它定义界面
- 当前物理锁定诊断显示：固定时间切片上的空间界面弱形式不够，D 支过渡层必须作为时空动态界面 \(F(t,x,z)=ell^2 R_tilde\mp1=0\) 处理；相关脚本为 `kg_examples/diagnose_d_trace_interface_weak_form.py` 与 `kg_examples/diagnose_d_interface_spacetime_kinematics.py`，笔记为 `research-notes/116-D支迹方程界面弱形式与时空界面必要性.md`
- 当前 D 支时空动态界面已推进到速度律核心模块：
  - `kg_examples/diagnose_d_spacetime_interface_jump.py` 现在输出 `jump_law_residual_scaled`，用于从 trace/scalaron 薄层弱形式反求界面速度
  - `kg_examples/diagnose_d_interface_speed_law.py` 从 leading jump law 反求 `F_t` 与坐标法向速度，是 moving-interface reduced simulator 的核心模块
  - `research-notes/117-D支时空界面速度律首轮.md` 记录首轮结果：`ell=30,t=16,128x128` 可解比例约 `0.624`，可解点修正后 scaled residual p95 约 `2.17e-19`
  - 该模块仍不是完整 tensor matching 全动力学；无实根点必须升级到完整张量 jump 或允许界面空间形状共同调整，不能用 clipping/damping 补掉
- 当前界面图展示口径已升级：
  - 以后默认用连续线段展示过渡层界面，而不是散点；
  - `kg_examples/diagnose_d_trace_interface_weak_form.py` 的 `find_level_segments` 已保留线段端点；
  - `kg_examples/diagnose_d_interface_speed_law.py` 已改为连续线段速度图；
  - `kg_examples/prototype_d_moving_interface_reduced.py` 是第一版 moving-interface reduced prototype，每步重抽取 `raw_y=±1` 界面，但它冻结 `metric_inv/stress_trace/rho`，不是完整自洽动力学；
  - 该原型显示不能直接显式推进原始 `raw_y`，下一版应采用 signed-distance / body-fitted interface 表示并接入 tensor jump/matching；
  - 研究笔记：`research-notes/118-D支连续界面线图与moving-interface原型.md`
- 当前 signed-distance moving-interface 原型已完成：
  - 脚本 `kg_examples/prototype_d_signed_distance_moving_interface.py`
  - `d_plus=0` 表示 `raw_y=+1`，`d_minus=0` 表示 `raw_y=-1`
  - 每步用 D 支 leading trace/scalaron jump law 计算界面速度，并按 \(d_t+v_n|\nabla d|=0\) 推进后重初始化
  - 主输出：`visualizations/d_signed_distance_moving_interface_ell30_t16_96_dt001/summary.json`
  - 结果显示 signed-distance 表示比直接推进 `raw_y` 更稳，但 trace-only jump law 仍有约三分之一 unresolved 段
  - 研究笔记：`research-notes/119-D支signed-distance-moving-interface原型.md`
- 当前 D 支 tensor jump 首轮诊断已完成：
  - 脚本 `kg_examples/diagnose_d_tensor_interface_jump.py`
  - 研究笔记 `research-notes/120-D支tensor-jump首轮诊断.md`
  - 结论：trace-only speed law 只能保证取迹；固定参考 `raw_y=±1` 界面法向不能满足无迹 tensor matching
  - `ell=30,t=16,n=128` 下 rank-one mismatch median `0.018`，但目标自由法向与当前空间法向 alignment median `0.486`
  - 下一步应实现 local tensor-matching interface solver，让界面法向/形状成为未知量，而不是继续只调界面速度
- 当前已按用户要求转为 direct tensor interface 条件：
  - 脚本 `kg_examples/diagnose_d_direct_tensor_interface.py`
  - 研究笔记 `research-notes/121-D支direct-tensor-interface条件.md`
  - 该脚本直接检查 \(H_{ab}=-F_aF_b+\tilde g_{ab}F^2\)，trace 只是派生检查
  - `ell=30,t=16,n=128` 下，`10%` 验收 accepted `456/1361`，`1%` 验收 accepted `62/1361`
  - 后续不能再把 trace-only signed-distance moving interface 当作主闭合；必须做 tensor-driven interface reconstruction
- 当前 direct tensor 目标法向重构预览已完成：
  - 脚本 `kg_examples/plot_d_direct_tensor_normal_reconstruction.py`
  - 研究笔记 `research-notes/122-D支direct-tensor目标法向重构预览.md`
  - 输出 `visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128/` 与 `visualizations/d_direct_tensor_normal_reconstruction_ell30_t16_128_strict001/`
  - 图中蓝色短线是 direct tensor 目标切线，红色短线是目标法向；它们只是重构诊断，不是新物理规则
  - 结论：target speed 约为 `1`，但 target/current normal alignment median 约 `0.41`，说明下一步必须重构界面法向/形状
- 当前 tensor-driven 界面重构原型已完成：
  - 脚本 `kg_examples/prototype_d_tensor_reconstructed_interface.py`
  - 研究笔记 `research-notes/123-D支tensor-driven界面重构原型与失败分类.md`
  - 输出 `visualizations/d_tensor_reconstructed_interface_ell30_t16_128/` 与 `visualizations/d_tensor_reconstructed_interface_ell30_t16_128_strict001/`
  - `10%` 验收下 direct tensor accepted 区域可重构为 `29` 条曲线；`1%` 下为 `16` 条短曲线
  - 该结果只说明 accepted 区域可积分成曲线，不说明全界面已闭合
- 当前 direct tensor rejected 模式分类已完成：
  - 脚本 `kg_examples/diagnose_d_direct_tensor_rejection_modes.py`
  - 输出 `visualizations/d_direct_tensor_rejection_modes_ell30_t16_128/` 与 `visualizations/d_direct_tensor_rejection_modes_ell30_t16_128_strict001/`
  - `10%` 下 rejected 主类为 `tensor_residual` 和 `rank_tail+negative+tensor_residual`
  - 后续应分型处理 unresolved：前者补完整 tensor matching，后者做分辨率/ell/层宽稳定性检查
- 当前 direct tensor 失败原因复查已完成：
  - 脚本 `kg_examples/diagnose_d_direct_tensor_failure_features.py`
  - 脚本 `kg_examples/diagnose_d_direct_tensor_least_squares.py`
  - 研究笔记 `research-notes/124-D支direct-tensor失败原因复查.md`
  - 失败明显集中在低密度、高 `|grad raw_y|`、局部 \(\tilde g\) 条件数高的区域；
  - 提高分辨率会改善但不消除顽固失败；
  - `tensor_residual` 型大多可由局部最小二乘 \(F_a\) 修复，因此不应再当作物理失败；
  - 后续 direct tensor solver 应改用完整张量残差最小二乘，真正困难集中到 `rank_tail+negative+tensor_residual` 与 `no_real_covector`
- 当前 least-squares 界面重构与分辨率复查已完成：
  - 脚本 `kg_examples/prototype_d_tensor_reconstructed_interface.py` 现支持 `--solver least_squares`
  - 研究笔记 `research-notes/125-D支least-squares界面重构与分辨率复查.md`
  - 输出 `visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_128/` 与 `visualizations/d_tensor_reconstructed_interface_ls_ell30_t16_160/`
  - `ell=30,t=16,n=128` 下 accepted `867/1361=0.6370`
  - `ell=30,t=16,n=160` 下 accepted `1214/1847=0.6573`
  - 结论：分辨率确实需要提高，用于薄层解析和界面几何收敛；但主突破来自把 \(F_a\) 求解从 rank-one 投影改为完整张量残差最小二乘，分辨率不能替代 tensor matching / 层内 profile / 两侧 bulk matching
- 当前局部分辨率收敛与 multistart 复查已完成：
  - 脚本 `kg_examples/diagnose_d_ls_local_resolution_convergence.py`
  - 研究笔记 `research-notes/126-D支局部收敛与multistart复查.md`
  - 输出 `visualizations/d_ls_local_resolution_convergence_ell30_t16_128_160_192/`
  - multistart 输出 `visualizations/d_tensor_reconstructed_interface_ls_multistart_ell30_t16_160/`
  - `W1-W4` 四个默认失败窗口为 `[-4,0]x[4,8]`、`[0,4]x[4,8]`、`[-4,0]x[8,12]`、`[0,4]x[8,12]`
  - `n=128/160/192` 全局 accepted 不单调，未解 residual p95 仍约 `1`
  - multistart 只额外修复 `28` 段，因此剩余失败更应转入层内 profile / 两侧 bulk matching 检查
- 当前低密度/计算精度误差检查已完成：
  - 研究笔记 `research-notes/127-D支低密度区是否为计算精度误差.md`
  - `rho_floor=1e-10,1e-12,1e-14` 不改变 `n=160` direct tensor least-squares 统计
  - `probe_dt` 与层内 `samples` 扫描也没有系统性压低 residual≈1
  - 高密度支撑区仍有未解段，因此普通精度误差不是主因
  - 但低密度区仍应标记为 `low-density-untrusted`，后续不能用低密度失败直接否定 D 支
- 当前 residual 诊断口径已按用户纠正修正：
  - 研究笔记 `research-notes/128-A参考残差与D支动态匹配的关系.md`
  - A 分支参考解上的 jump residual 不是 D 支静态失败率；
  - 它只衡量 A 参考解距离 D 支约束/匹配流形有多远；
  - D 支真实演化中 \(\rho,S,\tilde g\) 和界面都应改变；
  - 后续主线应转为 D 支一致初值构造 + dynamic level-set/body-fitted interface evolution，而不是继续静态检验 A 参考解是否逐点满足 D 支方程
- 当前 D 支同初态 reduced 动态原型已完成首轮：
  - 脚本 `kg_examples/simulate_d_reduced_dynamic_same_initial.py`
  - 研究笔记 `research-notes/129-D支同初态reduced动态演化首轮.md`
  - `coordinate_matter_evolution.py` 现在返回 `flux_x/flux_z`，用于守恒型通量
  - `n_cons=sqrt(|g~|)j^t` 必须保符号，不能按普通密度裁剪为非负
  - frozen initial \(\tilde g\) 下，`n=64,dt=1e-6` 可到最后可靠时间片 `t=7.4e-05`，此时 `rho_rel_l1_support=1.58e-3`，下一步质量壳判别式穿零
  - 每步代数重构 \(\tilde g[\rho]\) 的 `instant_transform` 模式在 step 2 即失败，原因是 \(Q\) 中 \(\partial_t^2\sqrt{\rho}\) 的有限差分 \(1/dt^2\) 放大
  - 后续不要把 `instant_transform` 当作主路线；应把 \(\tilde g\)、scalaron 或界面/层变量作为独立演化对象，并用 direct tensor jump/body-fitted matching 闭合
- 当前 D 支修正 \(\rho_A\to\tilde\rho\) 后的三域动力学原型已完成首轮：
  - 脚本 `kg_examples/simulate_d_tridomain_full_dynamics.py`
  - 研究笔记 `research-notes/131-D支修正rho映射后三域动力学首轮.md`
  - 初态使用正测度关系 `sqrt(|g~|)rho~=|X|rho_A/m^2`，不再使用 `rho_tilde=rho_A`
  - 默认输出比较对象：
    - `ntilde_measure=sqrt(|g~|)rho~`
    - `n_cons=sqrt(|g~|)rho~ g~^{tν}u_ν`
  - 首轮稳定输出：
    - `visualizations/d_tridomain_full_n64_dt1e6_ell300_mp300_stable57/`
    - `t=5.7e-05`
    - `measure_rel_l1_support=2.29e-2`
    - `n_cons_rel_l1_support=5.92e-5`
    - `disc_min_support=4.52e-7`
  - 同参数下一步 `t=5.8e-05` 质量壳判别式变负，不能继续把 floor 后的 \(u_t\) 当物理解
  - 当前 `geometry_closure=inert_tridomain` 只是可执行闭合：bulk metric 惯性保持，interface/jump 作为诊断；不能声称是最终完整张量几何演化
  - 当前局域双高斯在该 D 变换下支撑区天然深处 saturated 区：`ell=300` 全 saturated，快速检查 `ell=1` 仍约 `81%` saturated
  - `ell=M_P=300` 时 \(M_P^2/\ell^2=1\)，plateau 项未被参数压低；若要饱和区近似 inert，需扫 \(\ell\gg M_P\) 或明确处理 plateau 常数项
- 当前 1550nm 物理标定支持已加入：
  - 模块 `kg_examples/physical_units.py`
  - `simulate_d_tridomain_full_dynamics.py` 支持 `--physical-optical`
  - 内部单位为自然单位 eV / eV\(^{-1}\)，输出记录 μm/fs
  - 1550nm 且自洽 `m=omega0/10` 时：
    - \(k_0=0.799898054\,{\rm eV}\)
    - \(\omega_0=0.803927793\,{\rm eV}\)
    - \(m=0.080392779\,{\rm eV}\)
    - \(\sigma_\parallel=2.368\,\mu{\rm m}\)
    - \(\sigma_\perp=1.776\,\mu{\rm m}\)
    - 域半宽 \(29.603\,\mu{\rm m}\)
    - 相遇时间 \(42.10\,{\rm fs}\)
  - 物理跑数值前必须确认 \(\ell/l_P\) 和波函数/光强归一化；否则 \(M_P\) 源项比较没有真实实验强度含义
  - 用户已确认首轮采用 KG 内积归一化：
    \(N_{\rm KG}=\int i(\psi^*\dot\psi-\psi\dot\psi^*)dx dz=1\)
  - 用户已确认 \(\ell/l_P\) 按推荐扫描；首轮扫描 `1e29,1e35,1e45,1e60`
  - 当前首轮正式物理标定运行采用 `ell/l_P=1e60,n=64`，输出在 `visualizations/d_physical_1550nm_kg_ell1e60_n64_stable17/`
  - 可信截止 `step=17,t=8.39e-5 fs`；下一步质量壳判别式穿零，因此后续重点是切片/body-fitted 或判别式边界处理
  - 后续已完成 active/support 边界诊断：
    - `simulate_d_tridomain_full_dynamics.py` 支持 `--trusted-erosion` 与 `--stop-mask support|trusted`
    - `trusted` 定义为原始 support 按 8 邻域向内腐蚀，用于区分 support 阈值边缘伪故障和主体演化
    - 当前默认下一轮设置为 `active_dilation=2, trusted_erosion=1, stop_mask=trusted`，但仍必须报告原始 support 指标
    - `ell/l_P=1e60,n=64,dt_old=2.5e-7` 下 trusted core 可完成 400 步，`measure_rel_l1_trusted=1.437e-4`、`disc_min_trusted=3.867e-6`
    - 同一时刻原始 support 外圈已坏掉，`measure_rel_l1_support=17.981`、`disc_min_support=-5.554e-8`，说明当前主要病灶在 support 阈值边缘/低密度尾部
    - 后续时空收敛检查显示：`dt_old=2.5e-7` 减半几乎不改变 `n64` 结果，而 `n64->n96` 是数量级改善；当前默认工作分辨率应提升到 `n=96`
    - 用户提醒后新增 `--initial-time` / `--initial-time-old-units`，可由 A 支任意时刻生成 D 初态
    - `t_old=8` 干涉窗口明显更苛刻：`n96` 半窗口 trusted measure 约 `1.10e-3`、trusted disc 约 `4.20e-7`；`n128` 未改善
    - 已完成 `t_old=8,n96,dt_old=1.25e-7` 时间步复查，结果与 `dt_old=2.5e-7` 几乎逐项相同
    - 因此 `t_old=8` 干涉窗口问题不是时间步太大，也不是普通 `n128` 加密能修复
    - 已新增 `kg_examples/postprocess_d_support_edge.py`，以后处理方式比较 binary support、trusted core、support edge band、tapered support 四种口径
    - support/trusted/taper 是数值区域/诊断表示，不是新物理；support 边缘表达继续作为数值边界卫生检查，但主线瓶颈已转向 `inert_tridomain` 几何闭合
    - 当前固定工作参数：`n=96, dt_old=2.5e-7`，暂不引入动态参数
    - `simulate_d_tridomain_full_dynamics.py` 现支持 `--geometry-closure restricted_extension` 与 `--metric-extension-variable-mode covariant|adm`
    - 质量壳判别式是逐点数组；停止条件取 support/trusted mask 上的最小值，不是支撑内平均
    - 低物质边界平直锚点已实现为 `--flat-anchor-boundary-layers`、`--flat-anchor-rho-frac`、`--flat-anchor-measure-frac`
    - 直接 \(g_{\mu\nu}\) 延拓和平直边界锚点首轮均被 admissible 检查拒绝，`metric_extension_admissible_blend=0`
    - ADM 分块延拓首轮也被拒绝：`adm_valid_fraction≈0.258`、`roughness_after_p95≈1.63e53`、`admissible_blend=0`
    - 后续已修正 ADM 无效点伪重构和强制 \(+--\) 签名验收，新增 branch-preserving signature 检查
    - `--metric-extension-adm-update-fields all|lapse|lapse_shift` 与 `--metric-extension-max-rel-change` 已加入
    - 当前最稳 working closure：`restricted_extension + adm + lapse + max_rel_change=3e-4`
    - 该闭合在 `t_old=8,n=96,steps=20,ell/l_P=1e60` 下完成，trusted measure `1.6315e-4`、trusted disc `4.0996e-7`、accepted blend `1.0`
    - cap 扫描汇总在 `visualizations/d_physical_1550nm_kg_told8_adm_lapse_cap_scan/summary.json` 与 `adm_lapse_cap_scan.png`
    - 旧时间窗 `[0,16]` 在当前 1550nm 标定下对应 `[0,T_max]`，`T_max=120.01529382382436 eV^-1=78.99550140570793 fs`；旧 `t_old=8` 对应 `T_max/2=39.49775070285396 fs`
    - 已跑 `T=0,T_max/2,T_max` 三切片局部短窗：`T_max/2` 和 `T_max` 的 D-A trusted transformed-measure 偏差分别约 `1.15e-4` 与 `3.58e-5`，但 `T=0` 的稀疏 full tensor anchor 把偏差抬到约 `1.95e-2`
    - 因此 full tensor matching 虽已可作为物理锚点来源，但需要门控；不能默认每步、每段 accepted row 都强锚定 metric 代表
    - 低分辨率全窗粗步长 `dt_old=0.1/0.01` 已在早期 D matter RK4 爆炸；而若用可信局部步长 `dt_old=2.5e-7` 跑完整 `[0,T_max]` 需要约 `6.4e7` 步，当前脚本不可直接生产全窗
- 新发现的重要结构、方法、约束后，及时更新本文件

## 2026-04-29 椭圆型 `Phi` 初值延拓窗口

- 当前受控延拓
  `Phi_lambda = 1 - lambda (1 - Phi_elliptic)`
  在固定
  `steps=20, dt=2.5e-4, ell=0.2`
  的短时测试下，已经显示出一个小而可用的窗口：
  `elliptic_lambda ~ 1e-4` 到 `3e-4`
- 其中：
  - `1e-4`：`C-B` 分离已可见，约束代价最温和
  - `3e-4`：分离增强，但约束代价明显上升
  - `1e-3` 及以上：当前同步规范弱场原型下已偏强
- 因而后续若要先得到可信的 `C-B` 分离，再控制住约束代价，应优先把 `elliptic_lambda = 1e-4` 当作新的默认基准。
- 这一判断已经过长时间复查：
  在 `ell=0.2, dt=2.5e-4, steps=200` 下，
  `elliptic_lambda = 1e-4`
  时 `C-B` 几何差异会继续积累到 `10^{-4}` 量级，而 Hamilton 约束差异仍保持在 `~1e-1` 量级；因此它不是只适合短时探针的参数，而是当前可用于正式长时间比较的默认基准。

## 2026-04-29 局域波包平直参考

- 当前固定平直参考初值已从“沿传播方向延伸到边界的准单色高斯束”切换为“真正空间局域的双高斯波包”
- 当前默认局域波包参数：
  - `m=1`
  - `k0=6`
  - `sigma_parallel=1.6`
  - `sigma_perp=1.2`
  - `alpha=0.5`
  - `phi0=0`
  - 初始中心 `(-6,-6)` 与 `(6,-6)`
  - 传播方向分别为 `+45^\circ` 与 `-45^\circ`
- 当前固定参考区域：
  `x,z \in [-20,20]`
- 当前数值确认：
  - 理论相遇时间 `t_meet ≈ 8.60`
  - 采样时间窗 `t ∈ [0,16]`
  - 边界相对密度最大值仅 `~1.3e-28`
- 当前结论：
  这套局域波包初值已经真正满足“边界近似真空平直背景”的要求，应作为后续 `B/C` 全动力学研究的默认平直参考准备态
- 当前在这组局域波包初值上，`B/C` 原型已完成第一轮实际接入：
  - `C(flat)` 仍逐项严格退回 `B`
  - 在 `ell=0.2, elliptic_lambda=1e-4` 下，`40` 步时已有 `10^{-6}` 量级的 `C-B` 几何分离
  - 延长到 `200` 步后，`C-B` 分离可积累到 `10^{-5}~10^{-4}` 量级，而边界相对密度仍仅 `~10^{-9}`

## 2026-04-29 `B/C` 几何演化的当前规范

- 当前 `B/C` 全动力学原型不再只支持同步规范；现已支持两种 `lapse` 规范：
  - `gauge_mode = synchronous`
  - `gauge_mode = 1plog`
- 当前 `1+log` 的实现保持零 shift，只把 `lapse` 作为额外演化变量加入：
  `lapse_t = -2 lapse K`
- 几何演化中的 `K_{ij}` 与 `K_beta` 方程已补上 lapse 的 Hessian 项与 `-grad(beta)·grad(lapse)` 项，因此这不是简单把 `N` 从 `1` 改成变量，而是完整接入了零 shift 下的 `1+log` 切片
- 当前物质特征边界条件也已经改成使用**当前 `lapse`** 来判定法向入流/出流，不再默认 `N=1`
- 当前判断：
  - 同步规范已知会在 `t=7.79925` 首次失败
  - `1+log` 目前作为新的默认长时间候选规范，长时间诊断仍在进行中

## 2026-05-04 D 支全窗优化状态

- `kg_examples/simulate_d_tridomain_full_dynamics.py` 已新增生产调度开关：
  - `--geometry-every N`
  - `--metric-extension-every N`
  - `--diagnostics-every N`
  - `--no-render`
- 这些开关只改变昂贵几何/诊断刷新频率，不改变 D 支物质方程本身。
- 当前优化版在 `n=96` 下短跑峰值 RSS 约 `100MB`，16G MacBook 的内存不是主要瓶颈。
- 关键 caveat：当前 `n_cons + u_i` 显式 RK4 变量在 `j^t≈0` 或时间切片翻转区域会病态；它会在 `t_old≈0.00137` 附近触碰质量壳/正密度边界，不能直接用于完整 `[0,T_max]` 生产演化。
- 实验性开关 `--matter-projection`、`--matter-positivity-limiter` 仅用于诊断约束保持需求，不能把其结果当成正式物理预测，除非后续证明 limiter 对守恒量和场方程残差的影响可控。
- `kg_examples/simulate_d_tridomain_full_dynamics.py` 现有第一版局部时间物质推进：
  - `--matter-variable-mode local_time`
  - `--local-time-patchwise`
  - `--local-time-stationary-candidates`
  - `--local-time-patch-bboxes`
- 该模式用 `measure_density=sqrt(|gtilde|)rho_tilde` 作主变量，在局部 `tau=t+a x+b z` chart 中推进，再把结果重构回实验室 `t` 切片。
- `--local-time-stationary-candidates` 通过解析 stationary point 补足粗 `(a,b)` 网格漏掉的窄时间锥；它是坐标图册改进，不是新物理参数。
- `--local-time-patch-bboxes` 只把每个 patch 的计算裁到一层 halo 包围盒，保持中心差分 stencil；这是性能优化，不是削峰/阻尼。
- 当前局部时间模式已在 `t_old=0,n64,20` 步和 `t_old=8,n32,20` 步达到 trusted 覆盖 `100%`、正密度和质量壳裕度正常。
- 仍需牢记 caveat：当前 patchwise 实现还是 mask/recombine 原型，不是最终 conservative patch-boundary flux；接入 `restricted_extension`/full tensor matching 前还要检查 metric 更新路径不要重新引入实验室 `n_cons/j^t` 病态。
- 2026-05-04 步长扫描后更新：
  - `--local-time-max-tilt` 默认从 `3` 改为 `5`，因为 `t_old=8,n64` 中 `max_tilt=3` 会漏掉约 `2.67%` trusted transformed measure，而 `max_tilt=5` 可恢复 `100%` 覆盖；
  - `local_time_measure_rel_delta` 不是局部时间守恒误差；真正局部坐标守恒量诊断是 `local_time_n_tau_rel_delta_sum`，对应 `sqrt(|gtilde|)rho_tilde j^tau`；
  - 当前安全候选步长暂定 `dt_old≈1e-5`；`dt_old=2.5e-5` 只作激进短窗探针；
  - `t_old=8,n64,max_tilt=5,dt_old=2.5e-5,steps=100` 已完成，质量壳仍正，但 wall time 过高，说明全窗前必须继续优化 patch 数/候选缓存/保守 patch flux。

## 2026-05-05 local_time 性能与生产 caveat 更新

- `simulate_d_tridomain_full_dynamics.py` 的 local_time 路径现在复用预计算 `metric_inv/det`，并用解析公式做 `tau=t+a x+b z` 的协变/逆变度规变换。
- 这些优化不改变 D 支方程；它们只减少重复矩阵逆和局部坐标变换开销。
- `--fast-profile` 现在尊重用户显式设置的 `--diagnostics-every 0`、`--interface-every 0`、`--metric-extension-every 0`，可用于测纯物质推进成本。
- 新增实验开关：
  - `--local-time-tile-size`：每个 tile 共用一个 chart；当前扫描显示不合格，不作为生产默认；
  - `--local-time-atlas-every`：复用 chart atlas 若干步；当前速度收益很小，长复用会恶化 `n_tau` 守恒诊断，不作为主线。
- 当前可信性能基准：
  - `t_old=8,n64,dt_old=2.5e-5,steps=100` 纯演化约 `21.3s`、RSS `~151MB`、coverage `1.0`、`disc_min_trusted≈1.60e-8`；
  - 完整旧 `[0,16]` 在 `dt_old=2.5e-5` 下粗估约 `30-38h`，尚未稳进 1 天；在 `dt_old=1e-5` 下仍约数天。
- 后续重启时不要再优先尝试粗化 chart 图册来提速；已经被覆盖和质量壳诊断否决。若继续追求 1 天全窗，应转向编译化 patch RK4、conservative multi-patch 通量、或更严谨的局部时间步/特征推进。

## 2026-05-05 并行调试工具

- 新增 `kg_examples/run_d_local_time_parallel_scan.py`：
  - 用独立子进程并行跑多个 D 支 local_time 短窗；
  - 汇总 `summary.json`、`run.log`、质量壳、正密度、local-time coverage、`local_time_n_tau_rel_delta_sum`；
  - 默认启用 `--quick-diagnostics --skip-fields-npz`，因此适合快速调试，不适合最终物理图输出。
- `simulate_d_tridomain_full_dynamics.py` 新增：
  - `--quick-diagnostics`：短窗调试时不强制端点曲率、interface、plateau 和 A-reference 重型诊断；
  - `--skip-fields-npz`：跳过压缩保存完整场；
  - `matplotlib` 延迟导入，`--no-render` 时不再付出绘图库启动成本。
- 推荐短窗调试命令骨架：
  `python3 kg_examples/run_d_local_time_parallel_scan.py --output <out> --max-workers 2 --initial-times-old 0,8 --dt-olds 2.5e-5 --resolutions 64 --steps 50`
- 当前建议：
  - 默认 `--max-workers 2`，机器空闲时可试 `4`；
  - 不建议开满 `8`；
  - 最终/生产口径仍应关闭 quick-only shortcuts，至少做一次完整诊断对照。

## 2026-05-05 D 支高分辨率局部窗口当前入口

- 当前主要脚本：
  `kg_examples/simulate_d_local_window_from_a_snapshot.py`
- 它的正确流程是：
  - 在大域 `[-20,20]` 旧坐标、`n_full=640` 上构造 A 支正频 KG 快照；
  - 裁剪到物理 `[-9,9]um`；
  - 用裁剪得到的 `rho,S,gtilde` 初始化 D 支局部窗口演化；
  - 禁止把小窗口直接当作周期 FFT 初始域。
- 当前推荐短窗诊断命令参数：
  - `--dt-old 2.5e-7`
  - `--full-resolution 640`
  - `--window-um 9`
  - `--local-time-max-tilt 5`
  - `--local-time-step 0.5`
  - `--evolve-mask trusted`
  - `--boundary-stencil-mask active`
  - `--boundary-flux-mode matching --matching-flux-weight abs_n` 作为当前最新 trusted-buffer 守恒流匹配候选
  - `--boundary-flux-mode face` 仍保留为无耗散守恒通量基线；`central` 保留为旧对照，`rusanov` 只作耗散诊断
  - `--freeze-uncovered-chart`
  - `--chart-boundary-halo 1`
- 当前新增掩码含义：
  - `support`：初态密度和 transformed measure 双阈值支撑；
  - `trusted`：`support` 内缩一层，用作 bulk 演化和主要物理诊断；
  - `active`：`support` 外扩两层，用作局部窗口 halo；
  - `evolve_base`：实际显式推进区，当前推荐为 `trusted`；
  - `boundary_stencil_mask`：固定缓冲/ghost 区；当前推荐为 `active`，即差分通量可读取 active halo，但只把 `trusted` 写回；
  - `boundary_flux_mode`：`central` 为旧中心差分基线；`face` 为显式有限体积面通量；`rusanov` 为带数值黏性的诊断通量，不是新物理源项；
  - `boundary_flux_mode=matching`：只在 `trusted-buffer` 网格面上把两侧 one-sided flux 替换为单一法向 conserved-current flux；当前是物质连续性方程的界面通量匹配，不是 full tensor Einstein matching；
  - `chart_boundary`：`evolve_base` 内当前局部时间图册覆盖不到的临时界面单元；
  - `evolved_trusted_no_chart_halo`：排除 chart boundary 邻域后的 bulk 统计口径。
- 当前可信图/数据输出：
  - `visualizations/d_local_window_snapshot_tau0_trusted_evolve_5min/`
  - `visualizations/d_local_window_snapshot_tau_m3p5_trusted_evolve_5min/`
  - `visualizations/d_local_window_snapshot_tau0_fixed_buffer_steps165/`
  - `visualizations/d_local_window_snapshot_tau0_fixed_buffer_5min/`
  - `visualizations/d_local_window_snapshot_tau0_face_flux_smoke20/`
  - `visualizations/d_local_window_snapshot_tau0_rusanov_flux_smoke120/`
  - `visualizations/d_local_window_tau0_matching_flux_steps20/`
  - `visualizations/d_local_window_tau0_matching_flux_steps60/`
- 当前解释 caveat：
  - 正确密度对比仍是 `rho_A` 与 `rho_D_to_A=m^2*ntilde_D/|X_g[D]|`；
  - 低密度 support 边缘不能再被直接当作 bulk 演化失败证据；
  - fixed-buffer stencil 不是阻尼、削峰或改源项；第一版 matching flux 已把 trusted-buffer 物质通量由匹配条件决定；
  - Rusanov 强度扫描后，当前不采用 Rusanov 作为默认物理闭合；
  - 当前 matching flux 仍不是 full tensor Einstein interface matching，后续不能把这两层混同；
  - `tau=-3.5` 在 `[-9,9]um` 中 support 会触边，因此该运行是稳定性诊断，不是最终边界条件证明。

## 2026-05-07 \(f(R)\) 与 \(f(R,I_2)\) 作用量诊断入口

- 纯 metric \(f(R)\) 局部方向/单值性脚本：
  `kg_examples/diagnose_metric_fr_direction_matching.py`
- 纯 metric \(f(R)\) 普适函数完整二阶导数项回代：
  `kg_examples/fit_metric_fr_universal_function.py`
- \(f(R,I_2)\) 局部代数方向诊断：
  `kg_examples/diagnose_metric_fr_ricci2_direction_matching.py`
- \(f(R,I_2)\) 普适二维函数代数项回代：
  `kg_examples/fit_metric_fr_ricci2_universal_algebraic.py`
- \(f(R,I_2)\) 普适二维函数完整 metric 方程回代：
  `kg_examples/fit_metric_fr_ricci2_universal_full.py`
- 口径：
  - \(I_2=\tilde R_{\mu\nu}\tilde R^{\mu\nu}\)；
  - \(Q_{\mu\nu}=\tilde R_{\mu\alpha}\tilde R^\alpha{}_\nu\)；
  - \(f(R,I_2)\) 的代数主项是
    \(f_R R_{\mu\nu}+2f_{I_2}Q_{\mu\nu}-\frac12 f g_{\mu\nu}\)。
- Caveat：
  - `fit_metric_fr_ricci2_universal_algebraic.py` 尚未包含完整 metric variation 的导数项；
  - 因此它只能否定“单靠代数普适二维函数已足够”，不能否定完整 \(f(R,I_2)\)。
  - `fit_metric_fr_ricci2_universal_full.py` 已包含完整导数项，但当前是 fixed-background 回代/拟合诊断，不是完整 D 支动力学求解器。

## 2026-05-07 gBCD 高分辨率当前入口

- 当前最低可信条纹分辨率仍是 `n=384`、窗口 `[-9,9]um`、`core10`，每 1550nm 波长约 `10.05` 点。
- 主要脚本：
  - `kg_examples/fit_gbcd_principal_constraint_projection_sparse.py`
  - `kg_examples/solve_gbcd_full_linear_metric_update_sparse.py`
  - `kg_examples/diagnose_gbcd_metric_update_residual_sources.py`
  - `kg_examples/scan_gbcd_metric_update_alpha.py`
  - `kg_examples/export_gbcd_plus_initial_package.py`
  - `kg_examples/test_gbcd_plus_initial_package_one_step.py`
  - `kg_examples/scan_gbcd_plus_alpha_admissibility.py`
  - `kg_examples/scan_gbcd_plus_local_guard.py`
- 新增选项：
  - `--atom-family matter4|matter4_plus_normal|matter6_normal`：用于辅助张量 \(\mathcal C_{\mu\nu}\) 的局域张量基底诊断；
  - `--metric-variable-slices plus|center_plus|minus_center_plus`：用于 metric 初值投影诊断，默认仍为 `plus`。
- 当前最佳 `tau=-3.5` 单独 metric update 结果：
  - `matter4 + force-mask-erosion=0 + plus-only`；
  - exact nonlinear residual weighted mean `0.13917`；
  - 输出在 `visualizations/equation_first_gbcd_sparse_pipeline_n384_taum3p5_core10_pen1e4_forceedge_coeff_cg6000/`。
- 当前最佳可推进 `tau=-3.5` 初始加速度包：
  - `matter4 + plus-only + metric_active_dilation=0 + auto_bad_zero`；
  - exact nonlinear residual weighted mean `0.13220`；
  - 初始 \(\rho^\tilde\) 拉回偏差约 `1.4e-17`；
  - 一步 \(g_+\) 测度偏差约 `1.09e-4`；
  - 负质量壳判别式比例 `0`；
  - 输出在 `visualizations/equation_first_gbcd_plus_initial_package_n384_taum3p5_core10_noactiveedge_guarded/` 和 `visualizations/equation_first_gbcd_plus_initial_package_one_step_n384_taum3p5_core10_noactiveedge_guarded/`。
- 三切片 guarded package 验证：
  - 汇总输出在 `visualizations/equation_first_gbcd_guarded_three_tau_summary/`；
  - `tau=-3.5` residual `0.13220`，一步加权测度偏差 `1.04e-4`；
  - `tau=0` residual `6.60e-5`，一步加权测度偏差 `3.00e-5`；
  - `tau=+3.5` residual `0.01324`，一步加权测度偏差 `2.88e-5`；
  - 三者一步 \(g_+\) 重构负质量壳判别式比例均为 `0`。
- 重要结论：
  - `matter6_normal` 改善 projection system residual，但不改善完整非线性回代；
  - `center_plus` 只改几何会出现“线性残差低、非线性残差高”的不自洽；
  - `center_plus + source` 已被联合投影诊断否定为不健康路线；
  - 当前主线应视为 `plus-only` 初始加速度求解器，而不是大幅修改中心切片或 \(\rho^\tilde\)；
  - 旧 active-edge 解的可容许性失败主要来自边界/退化点，下一步要把局部可容许性守卫写成正式强约束。

## 2026-05-08 trace=0 最小投影方程组当前入口

- 当前主线不是构造 D 支数值器，而是寻找显式、可推广的 equation-first Einstein-like 方程；数值器只是检验工具。
- 最新理论笔记：
  `research-notes/169-trace0最小投影方程组与patch条件.md`
- Helmholtz 检查笔记：
  `research-notes/170-trace0投影方程的Helmholtz检查.md`
- 最小辅助 action 原型：
  `research-notes/171-最小辅助应力场action原型.md`
- shadow-field 分支：
  `research-notes/172-shadow-field分支的最小方程组.md`
- shadow stress 与分支传播：
  `research-notes/173-shadow-stress与物理分支传播条件.md`
- stealth-state 最小 no-go：
  `research-notes/174-stealth-state-action候选的最小no-go.md`
- n=384 trace=0 局部复检：
  `research-notes/175-n384-trace0局部复检.md`
- trace0 局部守恒兼容性：
  `research-notes/176-trace0局部守恒兼容性.md`
- trace0 稀疏全局守恒扫描：
  `research-notes/177-trace0稀疏全局守恒扫描.md`
- trace0 硬守恒 KKT 与 ALM 原型：
  `research-notes/178-trace0硬守恒KKT与ALM原型.md`
- 当前最小候选：
  \[
  \tilde G_{\mu\nu}
  =
  \tilde T_{\mu\nu}/M_P^2
  +\mathcal C_{\mu\nu},
  \quad
  \mathcal C_{\mu\nu}\in
  \mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu,r_\mu r_\nu,u_{(\mu}r_{\nu)}\},
  \]
  \[
  \tilde g^{\mu\nu}\mathcal C_{\mu\nu}=0,
  \qquad
  \tilde\nabla^\mu\mathcal C_{\mu\nu}=0.
  \]
- 同一候选的投影写法：
  \[
  \Pi_E^\perp
  \left(\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2\right)=0,
  \quad
  \tilde g^{\mu\nu}\Pi_E
  \left(\tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2\right)=0.
  \]
- 重要 caveat：
  - 不可混用 \(\rho\) 和 \(\tilde\rho\)；
  - \(Q\to0\Rightarrow\mathcal C\to0\) 仍是分支/边界/正则性条件；
  - 非退化 patch 需要 \(d>2\)、\(\Delta=u^2r^2-(u\cdot r)^2\neq0\)，并处理 \(w^2=0\)、\(r^\perp=0\) 和严格 \(1+1d\)；
  - Helmholtz 检查显示裸投影方程不宜声称为普通纯 metric 作用量方程，因为 \((\delta\Pi_E)\mathcal R\) 一般破坏自伴性，且 7 条 metric 条件存在 Noether 计数压力；
  - action 化时真正重要的是辅助 sector 的有效能动量 \(\Theta_{\mu\nu}\)，不是裸 \(\lambda_I E^I_{\mu\nu}\)；
  - 若辅助 action 显含 \(u,r,\tilde\rho,S\)，通常会改写物质变分；保持测地线需要补偿或独立 \(U,R\) shadow fields；
  - shadow-field 分支引入独立 \((\Phi,\sigma)\)，用 \(U=d\Phi,R=d\ln\sqrt\sigma\) 构造 \(E_{\rm sh}\)，并通过同型动力学和初始/边界条件选择 \(U=u,R=r\) 的物理分支；
  - shadow stress 检查显示线性乘子和 \(Q_{\rm sh}\)-门控 variational weight 会产生额外 source；物理分支传播要求状态项对 \(\Phi,\sigma\) 一阶静默；
  - 具体 stealth-state action 候选检查显示最小局域 ansatz 有 no-go：抵消 \(C\)-变分会同时取消 metric 源；
  - `n=384` 局部复检支持普通 trace=0 继续作为主标量闭合，\(\mathsf q_E\) 暂降级为主符号/正定范数诊断工具；
  - trace0 局部守恒检查显示干涉中点几乎自然守恒，但分离态自然散度很大；
  - trace0 稀疏全局守恒扫描显示 `trace0 + full conservation` 兼容；
  - trace0 硬守恒原型显示 ALM 在三切片上比 penalty 进一步压低 full divergence，且只小幅增加代数残差；
  - 下一步应整理 equation-first 显式 proposal；action 化若继续则转向真实辅助场 \(\chi\)-sector；严格 KKT 数值求解若继续则需要 MINRES/Schur complement/更强预条件器。
