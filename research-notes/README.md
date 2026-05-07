# 研究留档

这个目录用于保存可人工查阅的阶段性成果，与 `agent-memory/` 分离。

使用规则：

- 每个相对完整的小步骤单独存成一个 Markdown 文件
- 文件尽量按编号递增，便于追踪逻辑发展
- 文件内容优先记录：问题、设定、推导、结论、未决问题
- 即使推导未完成，也应保留中间版本，避免丢失思路

当前文件：

- `001-KG作用量的Madelung分解与量子势定位.md`
- `002-量子势几何化的成功失败与可接受判据.md`
- `003-disformal-ansatz候选与第一轮筛选.md`
- `004-四类度规变换的初步比较分析.md`
- `005-四类变换代回作用量与HJ方程的直接计算.md`
- `006-HJ层面的吸收Q求解与更广泛存在性.md`
- `007-快速筛选：作用量兼容性与存在条件.md`
- `008-方法论修正：作用量等价不是必要条件.md`
- `009-变分筛选：最简单共形与类II方案.md`
- `010-密度随度规变换后的重新检验.md`
- `011-方法修正后：先变分再代回变换.md`
- `012-一般类II能否同时满足HJ与连续性方程.md`
- `013-类II中A=1时B是否必然含Z的复核.md`
- `014-按逆度规ansatz直接消去拉氏量中的Q项.md`
- `015-最简单拉氏量消项解的后续检验.md`
- `016-连续性方程中代入X=m2-Q后的精确比较.md`
- `017-J^mu∂_muQ与新几何测地线方程是否等价.md`
- `018-最小混合逆度规下的连续性方程no-go.md`
- `019-去掉sqrt-g-rho不变假设后的纯u-ansatz翻盘.md`
- `020-最简单解的可逆性签名因果结构与经典极限.md`
- `021-X正负对应的纯ansatz双分支结构.md`
- `022-对Cursor对照分析的判断.md`
- `023-阶段性总结与Einstein作用量耦合框架.md`
- `024-S_EH[g]与S_EH[g~]的弱展开与差异量级.md`
- `025-1+1d自由KG双高斯波包基线模型.md`
- `026-1+1d双高斯KG基线算例结果.md`
- `027-用双高斯KG例子检验四种一致性机制.md`
- `028-Q依赖引力作用量的一般要求与例子分析.md`
- `029-双高斯模型上测试FQ=1-1pluskQ2.md`
- `030-1+1d平直背景下EHg与EHgtilde的真正差异.md`
- `031-1+1d双高斯模型下两种总作用量的直接比较结果.md`
- `061-Q约定改为正BoxR除以R.md`
- `062-混合投影型逆度规解的统一整理.md`
- `063-度规退化、不可逆变换与真奇点的文献分层.md`
- `064-纯共形变换能否满足几何化要求.md`
- `065-混合ansatz新解D等于Z-E等于负aY-F等于aX.md`
- `066-一般混合逆度规ansatz的行列式恒正条件.md`
- `067-在X等于0处让rho波浪吸收奇异性的可能性.md`
- `068-混合ansatz在X等于0附近的严格局域分析.md`
- `069-混合ansatz下散度等价的主方程与绕开kappa除以X的条件.md`
- `070-关于rho波浪为零时速度无定义的物理解读与混合投影反例.md`
- `071-测度密度趋零与X等于0处速度无定义的统一解释.md`
- `117-D支时空界面速度律首轮.md`：把 D 支 trace/scalaron 的时空薄层弱形式改写为界面速度律，并给出 `ell=30` 的首轮数值结果。
- `118-D支连续界面线图与moving-interface原型.md`：把界面图从散点改为连续线段，并记录第一版 moving-interface reduced prototype 及其限制。
- `119-D支signed-distance-moving-interface原型.md`：把 moving-interface 原型从直接推进 `raw_y` 升级为 `d_plus/d_minus` signed-distance level set，并记录 `ell=30,t=16,n=96` 的短程结果。
- `120-D支tensor-jump首轮诊断.md`：把 D 支 trace jump law 升级到张量薄层主部，发现固定参考界面法向不能满足无迹 tensor matching，但 free-covector rank-one 兼容性中位数较好。
- `121-D支direct-tensor-interface条件.md`：按用户要求不再以 trace 为闭合条件，直接由完整 tensor jump 条件反推 \(F_a\)，并给出 `ell=30,t=16,n=128` 的 `10%/1%` 数值验收结果。
- `122-D支direct-tensor目标法向重构预览.md`：把 direct tensor 条件反推出的目标 \(F_a\) 画成目标法向/切线，确认后续模拟器必须重构界面形状而不是继续沿 trace speed 推进。
- `123-D支tensor-driven界面重构原型与失败分类.md`：用 direct tensor accepted 段积分出目标界面曲线，并把 rejected 段分成 tensor residual、rank-one/negative 等失败类型。
- `124-D支direct-tensor失败原因复查.md`：通过局部特征、分辨率/ell/层宽检查和最小二乘 \(F_a\) 求解，确认 `tensor_residual` 大多是投影算法问题，而顽固失败集中在低密度/极陡界面/可能局部不可闭合段。
- `125-D支least-squares界面重构与分辨率复查.md`：把 tensor-driven 界面重构器升级为局部 least-squares \(F_a\) 求解，并确认分辨率需要提高但不是替代 tensor matching 的物理闭合。
- `126-D支局部收敛与multistart复查.md`：在 `W1-W4` 失败窗口上做 `128/160/192` 局部分辨率对照，并用 multistart \(F_a\) 复查剩余失败是否只是初值问题。
- `127-D支低密度区是否为计算精度误差.md`：通过密度分箱、`rho_floor`、`probe_dt` 和层内 `samples` 敏感性检查，判断普通计算精度误差不是 residual≈1 的主要来源，但低密度区需单独标记。
- `128-A参考残差与D支动态匹配的关系.md`：修正 residual 诊断口径，明确 A 参考解不严格满足 D 支 jump law 是预期内的，真正任务是构造 D 支一致初值和动态界面演化。
- `129-D支同初态reduced动态演化首轮.md`：实现同 A 初态出发的 D 支 reduced 动态原型，确认保符号守恒密度后 frozen baseline 可短时演化，而每步代数重构 \(\tilde g[\rho]\) 会因 \(Q\) 的二阶时间差分立刻病态。
- `130-D支三域规则可模拟性与rho映射错误审计.md`：确认弱曲率区、饱和区、过渡界面三域规则可作为 bulk-interface 模拟路线，但必须处理饱和区残留 \(f\) 项，并系统修正 \(\rho_A\to\tilde\rho\) 后重算 D 支源项、support 和 observable。
- `131-D支修正rho映射后三域动力学首轮.md`：实现修正 \(\tilde\rho\) 初态后的 D 支三域动力学原型，跑通 `n=64,dt=1e-6,ell=M_P=300` 到可信截止 `t=5.7e-05`，并确认饱和区几何闭合仍是主要理论缺口。
- `132-1550nm物理标定参数准备.md`：按用户要求加入 1550nm 光学波包物理标定，确定 \(k_0,\omega_0,m,\sigma\)、真实 \(M_P\) 和自然单位换算，并标记仍需确认 \(\ell/l_P\) 与波函数归一化。
- `144-放宽纯gtilde限制后的ur引力作用量理论分析.md`：放宽 pure `gtilde` 限制后，建立精确拉回 EH action 的存在性基准，并提出首个低阶 \(A^{\mu\nu}(u,r)\tilde R_{\mu\nu}-2V\) 检验路线。
- `145-从g表象Einstein方程变换得到gtilde方程的分解路线.md`：从 A 支 Einstein 方程直接变换，定义 \(\mathcal R_{\mu\nu}=T^A-T^{(\tilde m)}-M_P^2\mathcal H\)，把 action 猜测问题转化为剩余张量的结构分解问题；已补入 1550nm 高斯干涉上 \(\mathcal R_{\rm need}\) 的 pure-geometry 诊断和留一切片 \(\tilde g\)-函数性诊断，结论是完整 \(\mathcal R_{\rm need}\) 的拟合成功主要是平凡 EH 主项，非平凡 \(\tilde T/M_P^2\) 仍不能由简单局部 pure geometry 特征稳定表示。
- `146-equation-first-gBCD闭合方程与守恒检验.md`：从方程出发提出 gBCD 张量壳，比较 transverse/full 守恒硬约束，并把问题重写为辅助各向异性应力张量闭合。
- `147-gBCD系数ABCD的可能产生机制.md`：比较简单局部本构函数、单个局域势函数和辅助应力场三类机制；当前支持辅助应力场机制。
- `148-gBCD辅助应力场代表元规范与nullspace硬约束.md`：修正 KKT hard constraint 的数值泄漏，采用 nullspace 投影，并记录三切片 `norm/time` 代表元结果。
- `149-gBCD状态方程trace与剪切闭合首轮检验.md`：检验 `trace=0`、`trace=aQ+bQ^2` 和 `D=0` 三类最小代数状态方程，结论是它们能形式闭合但会明显恶化残差或引入病态。
- `150-Q趋零条件作为辅助场变分权重.md`：把 \(Q\to0\) 回 Einstein 条件改写为辅助场 action 中的 \(Q\)-门控范数权重，并完成三切片数值检验。
- `151-gBCD辅助场变分闭合的Euler-Lagrange方程.md`：推导 gBCD 辅助场 saddle-point 泛函、Euler-Lagrange 方程和其与 `norm/time/space/q-gated` 数值正则的对应。
- `152-无Q门控的最小辅助场action扫描.md`：修正 3+1d 守恒方程数量理解，并在无 \(Q\) 门控下扫描 norm/time/space 辅助场 action，得到当前默认 `time_weight=1e-5` 候选。
- `153-gBCD显式理论候选v0.md`：给出当前最明确的 equation-first 显式理论候选，包括场变量、主场方程、守恒条件、辅助场最小泛函和下一步推进器风险。
- `154-gBCD-v0一步lambda推进诊断.md`：实现并诊断固定 A 背景上的 v0 一步 \(\lambda\) 推进器，确认 hard conservation 可机器精度满足，但 lambda-only 不能让下一切片代数 gBCD 方程保持低残差。
- `155-gBCD约束演化分裂与principal投影.md`：发现 metric time-principal 主部 rank=3，建立 gBCD 的约束/演化分裂，并实现 principal-constraint projection，三组切片投影后残差降到 \(10^{-4}\) 以下量级。
- `156-高分辨率口径与full-linear-metric-update首轮.md`：确认 `n=96` 每波长仅约 2.5 点、只能作算法诊断；生成 `n=384` A 支条纹图，并实现 dense full-linear metric update 原型，确定 `ridge=1e-6` 为首轮较稳参数。
- `157-sparse-matrixfree求解器首轮.md`：实现 sparse full-linear metric update 和 sparse principal projection 首版；确认 full-linear 稀疏算子与 dense 一致，但当前迭代器和 soft constraint projection 尚未达到 dense hard pipeline 精度。
- `158-n384-core10首个高分辨率D支几何更新结果.md`：跑通 `n=384,tau=0,core10` 高分辨率 sparse pipeline，完整非线性残差加权均值约 `0.00644`，并确认 corrected metric 在核心区未退化；分离态快速尝试仍受求解器收敛限制。
- `159-n384三切片gBCD求解器升级与物理主线判断.md`：把 full-linear metric update 和 principal projection 的稀疏 matvec 改为 entry-array 向量化，完成 `n=384,core10,tau=-3.5,0,+3.5` 三切片检查；干涉中心残差约 `0.00644`，右侧分离态约 `0.0386`，左侧分离态约 `0.146`，当前主卡点定位为左侧分离态的约束/代表元闭合而非低分辨率或简单阻尼问题。
- `160-tau负3p5残差来源与center-plus初值投影检验.md`：对左侧分离态做残差分层，确认边缘层最差但非唯一根因；测试 edge-force、`time_weight=0`、`matter6_normal` 和 `center_plus` 初值投影，结论是只改几何的 center-plus 在线性层面有用但非线性不自洽，下一步必须联合处理几何、辅助应力、物质源和 \(\rho^\tilde\) 拉回。
- `161-D初值联合投影与plus-only主自由度.md`：实现联合投影原型与 trust-region 扫描，确认 `center_plus+source` 会要求荒唐的 \(\eta=\delta\log(\sqrt{|\tilde g|}\tilde\rho)\) 且破坏 \(\rho\) pullback；当前健康自由度是 `plus-only`，即保持初始物质/中心度规不变，由 D 方程确定下一切片或初始加速度。
- `162-plus-only初始加速度包与可容许性守卫.md`：把 `plus-only` 升级成可导出的 D 初始加速度包；确认旧 active-edge 失败来自局部边界/退化点可容许性问题，`metric_active_dilation=0 + auto_bad_zero` 在 `tau=-3.5,n=384,core10` 上给出 residual 约 `0.132`、一步测度偏差约 `1.1e-4`、负判别式比例 `0` 的当前最佳可推进初态。
- `163-三切片guarded-plus-only初始包验证.md`：把同一 guarded plus-only 流程扩展到 `tau=-3.5,0,+3.5`，确认三切片均能生成可推进初始包；干涉中心 residual 约 `6.6e-5`，右侧分离态约 `0.013`，左侧分离态约 `0.132`，三者一步物质重构均无负判别式。
- `164-gBCD投影型Einstein-like方程候选.md`：把 gBCD 从逐点 \(A,B,C,D\) 拟合提升为显式投影型 Einstein-like 方程；主方程是 Einstein 残差落在 \(\{\tilde g,uu,rr,ur\}\) 张量子空间，配合投影修正张量守恒、物质壳方程/连续性方程和 \(Q\to0\) 回 Einstein 分支条件。
- `165-gBCD投影方程闭合条件.md`：进一步分析投影方程的闭合条件，证明 \(W\) 不改变在壳解集，给出 Gram 行列式 \(\det H=\frac{d-2}{2}(u^2r^2-(u\cdot r)^2)^3\)，并区分守恒约束、测地线来源、\(Q\to0\) 分支和 action/Helmholtz 检查。
- `166-gBCD投影方程的3加1计数与主部.md`：分析 3+1d 下投影方程的方程计数和主部；确认它给出 6 个独立张量条件，主部是 harmonic gauge 下 Einstein 波算子的 \(\Pi_E^\perp\) 投影，但 Cauchy 适定性仍需 gauge-fixed principal symbol 检查。
- `167-gBCD投影方程的主符号缺口与标量闭合条件.md`：检查 harmonic gauge 后的 reduced principal symbol，发现仍有一个 \(E\)-方向主部模式 \(-[(\xi\cdot r)u-(\xi\cdot u)r]^2\) 未定；因此完整理论还需一个协变标量状态方程，并提出 \(u-r\) 平面正定迹闭合类作为首个候选。
- `168-qEtrace闭合三切片检验.md`：在高斯干涉三切片上检验 \(\mathsf q_E\)-trace 标量闭合；结果显示它的零源或简单 \(Q,Q^2\) 硬闭合版本均不如普通 trace=0，因此 \(\mathsf q_E\) 暂降级为主符号/辅助场正定范数候选。
- `169-trace0最小投影方程组与patch条件.md`：把当前数据支持的最小 equation-first 候选写成显式方程组：Einstein 残差等于一个落在 \(\{\tilde g,uu,rr,ur\}\) 子空间内的守恒、无迹修正张量，并明确 \(\rho\neq\tilde\rho\)、\(Q\to0\) 分支和 \(\Delta,w^2,r^\perp\) patch 条件。

补充说明：

- 自 `061-Q约定改为正BoxR除以R.md` 起，项目默认采用的新约定是 `Q = + Box sqrt(ρ) / sqrt(ρ)`
- 若查阅更早的历史笔记，需注意其中许多地方使用的是旧约定 `Q_old = - Box sqrt(ρ) / sqrt(ρ)`，应先按 `Q_old = -Q` 换回当前约定再读
