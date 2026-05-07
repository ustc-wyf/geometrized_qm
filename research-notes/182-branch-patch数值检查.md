# branch / patch 数值检查：\(C/\chi\) 有界性与低维 \(E_u\) patch

日期：2026-05-08

## 目的

接续 `181-branch-preservation与patch-transition首轮.md` 的两个具体检查：

1. 在高斯波包干涉三切片上检查 \(\mathcal q\approx0\) 区域中 \(C/\chi\) 是否有界；
2. 在 \(\Delta\approx0\)、\(r^\perp\approx0\) 附近检查
   \[
   \mathcal R_{\mu\nu}
   =
   \tilde G_{\mu\nu}-\tilde T_{\mu\nu}/M_P^2
   \]
   是否接近低维 patch
   \[
   E_u=\mathrm{span}\{\tilde g_{\mu\nu},u_\mu u_\nu\}.
   \]

这里使用已经得到的 `hard-alm` 解作为 trace0 + full conservation 的基准代表元。

## 新增脚本与输出

新增脚本：

```text
kg_examples/diagnose_trace0_branch_patch.py
```

运行输出：

```text
visualizations/trace0_branch_patch_check_n384_core10/summary.json
visualizations/trace0_branch_patch_check_n384_core10/branch_patch_summary.png
visualizations/trace0_branch_patch_check_n384_core10/branch_patch_taum3p5.png
visualizations/trace0_branch_patch_check_n384_core10/branch_patch_tau0.png
visualizations/trace0_branch_patch_check_n384_core10/branch_patch_taup3p5.png
```

参数：

- `full_resolution=384`
- `window_um=9`
- `fit_region=core10`
- `force_region=core10`
- `diagnostic_region=core10`
- 三切片：\(\tau=-3.5,0,+3.5\)

## 定义

脚本中使用

\[
\tilde Q_\rho
=
\tilde\nabla_\mu r^\mu+r_\mu r^\mu,
\qquad
\mathcal q=\tilde Q_\rho/m^2.
\]

由于物理单位下 \(\mathcal q\) 的绝对量级随切片变化明显，本轮同时使用相对强度

\[
q_{\rm rel}
=
\frac{|\mathcal q|}
{\mathrm{p95}(|\mathcal q|)\text{ in diagnostic region}}.
\]

分支函数取诊断型

\[
\chi(q_{\rm rel})
=
\frac{q_{\rm rel}^2}{q_{\rm rel}^2+q_0^2},
\qquad q_0=0.1.
\]

这不是最终物理归一化，只用于测试“当 \(q\) 接近零时，现有 hard-ALM 代表元是否自然满足 \(C\sim\chi\)”。

低维 patch 检查包含两种残差：

1. `Eu2_residual`：逐点拟合到 \(\mathrm{span}\{\tilde g,uu\}\)；
2. `Eu_trace0_residual`：逐点拟合到单个 trace0 方向
   \[
   u_\mu u_\nu-\frac{u^2}{3}\tilde g_{\mu\nu},
   \]
   这里是 \(2+1\) 数值切片，所以维数取 \(d=3\)。

## 主要数值结果

### \(\tau=-3.5\)

- `core10` 点数：2550；
- \(q_{\rm rel}\le0.1\) 点数：2171；
- \(|C|\) weighted p95：`3.50e3`；
- \(|C|/\chi\) weighted p95：`1.32e6`；
- `Eu2_residual` weighted p95：`3.73e-2`；
- `Eu_trace0_residual` weighted p95：`2.22e-1`；
- \(\Delta_{\rm rel}\le0.1\) 点数：81，低维 patch `Eu2` p95：`2.63e-2`。

### \(\tau=0\)

- `core10` 点数：1024；
- \(q_{\rm rel}\le0.1\) 点数：533；
- \(|C|\) weighted p95：`8.37e2`；
- \(|C|/\chi\) weighted p95：`6.94e4`；
- `Eu2_residual` weighted p95：`2.44e-1`；
- `Eu_trace0_residual` weighted p95：`4.42e-1`；
- \(\Delta_{\rm rel}\le0.1\) 点数只有 1。

### \(\tau=+3.5\)

- `core10` 点数：3589；
- \(q_{\rm rel}\le0.1\) 点数：2915；
- \(|C|\) weighted p95：`4.64e3`；
- \(|C|/\chi\) weighted p95：`2.36e6`；
- `Eu2_residual` weighted p95：`9.28e-2`；
- `Eu_trace0_residual` weighted p95：`3.47e-1`；
- \(\Delta_{\rm rel}\le0.1\) 点数：55，低维 patch `Eu2` p95：`3.69e-1`。

## 结论 1：\(Q\to0\Rightarrow C\to0\) 不是 hard-ALM 自动结果

三个切片都有大量 \(q_{\rm rel}\le0.1\) 的点，但 \(|C|/\chi\) 的 p95 仍很大。

因此当前最小条件：

\[
C\in E,\qquad
\tilde g^{\mu\nu}C_{\mu\nu}=0,\qquad
\tilde\nabla^\mu C_{\mu\nu}=0
\]

不会自动推出

\[
Q\to0\Rightarrow C\to0.
\]

这支持前一轮理论判断：GR branch 必须作为正则性/边界条件显式加入，例如

\[
\boxed{
C_{\mu\nu}
=
\chi(\mathcal q)\hat C_{\mu\nu},
\qquad
\chi(0)=\chi'(0)=0,
\qquad
\hat C_{\mu\nu}\text{ 有界}.
}
\]

也就是说，\(\chi\)-branch 不是装饰项，而是 proposal v1 必须补上的健康性规则。

## 结论 2：低维 \(E_u\) patch 对分离态有效，对干涉中心不足

低维 patch \(\mathrm{span}\{\tilde g,uu\}\) 对三切片的表现：

- \(\tau=-3.5\)：weighted p95 `3.73e-2`，很好；
- \(\tau=0\)：weighted p95 `2.44e-1`，明显变差；
- \(\tau=+3.5\)：weighted p95 `9.28e-2`，仍可接受但比左侧分离态差。

这说明：

1. 分离态中 \(\mathcal R\) 很大程度上可被低维 \(E_u\) patch 表示；
2. 干涉中心需要完整 \(E=\mathrm{span}\{\tilde g,uu,rr,ur\}\)，不能退化成只有 \(g,uu\)；
3. 单个 trace0 方向 \(uu-(u^2/d)g\) 通常不够，特别是 \(\tau=0\) 和 \(\tau=+3.5\)。

## 结论 3：core10 内真正 \(\Delta\approx0\) 点不多

在 `core10` 区域，\(\Delta_{\rm rel}\) 的 p50/p95 基本是 1。

这意味着这组三切片的主要问题并不是“处处接近 \(r\parallel u\) 的低维退化”，而是：

- 大部分 core10 点处于非退化四方向 patch；
- 只有少量点进入 \(\Delta_{\rm rel}\le0.1\)；
- 因此 patch transition 规则仍需要写入理论，但它不是当前三切片残差/branch 问题的主因。

## 对 proposal 的修正

当前 equation-first proposal 不应只写：

\[
\mathcal R_{\mu\nu}=C_{\mu\nu},
\quad
C\in E,
\quad
\mathrm{tr}_{\tilde g}C=0,
\quad
\tilde\nabla C=0.
\]

更健康的首版应写成：

\[
\boxed{
\mathcal R_{\mu\nu}
=
\chi(\mathcal q)\hat C_{\mu\nu}
}
\]

其中

\[
\hat C_{\mu\nu}\in E,\qquad
\tilde g^{\mu\nu}\hat C_{\mu\nu}=0\quad(\mathcal q\neq0),
\]

并且完整守恒条件是

\[
\tilde\nabla^\mu(\chi\hat C_{\mu\nu})=0.
\]

这里 trace0 可以等价施加在 \(C\) 或 \(\hat C\) 上，只要 \(\chi\neq0\)；在 \(\chi=0\) 面上应改用正则延拓/无通量条件。

## 下一步

1. 把 proposal v1 升级成带 \(\chi\)-branch 的 v1.1：
   \[
   C=\chi(\mathcal q)\hat C,\quad \hat C\text{ 有界},\quad \tilde\nabla(\chi\hat C)=0.
   \]
2. 在数值上重新求解带 branch 权重/有界性惩罚的 hard-ALM：
   - 未加 branch 时 \(C/\chi\) 明显不有界；
   - 加入 branch 后检查代数 residual 和守恒 residual 是否仍可同时保持小。
3. patch 方面保留 \(C\) 为主变量、\(\lambda_I\) 为局部坐标的规则；
   低维 \(E_u\) patch 只作为 \(\Delta\to0\) 局部图，而不是替代完整四方向 patch。

