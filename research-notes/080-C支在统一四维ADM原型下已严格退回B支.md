# 080 `C` 支在统一四维 `ADM` 原型下已严格退回 `B` 支

## 目标

确认统一四维 `ADM` 同步规范原型中，`C` 支在
\[
\Phi=1,\qquad \Pi_\Phi=0
\]
时，是否严格退回 `B` 支。

这一步是 `C` 支数值器进入正式比较前的最低一致性校准条件。

## 1. 修正内容

在

```text
kg_examples/adm_sources_abc.py
kg_examples/simulate_ab_adm_synchronous.py
```

中做了两类修正。

### 1.1 势函数与导数

将 `C` 支的势函数改成与旧稳定实现一致的平直极限口径：

\[
U(\Phi)=R(\Phi)\bigl(\Phi-\Phi^{1/3}\bigr),
\qquad
R(\Phi)=\frac{1}{\ell^2}\sqrt{\Phi^{-2/3}-1},
\]

并取

\[
U_\Phi(\Phi)=R(\Phi).
\]

于是自动满足

\[
U(1)=0,\qquad U_\Phi(1)=0.
\]

### 1.2 平直 `\Phi=1` 伪源扣除

在 `C` 支的有效 Einstein 源中，原始

\[
\Box\Phi
=
\frac{T^{(B)}}{3M_P^2}
\;-\;
\frac{2U-\Phi U_\Phi}{3}
\]

在 `\Phi=1` 时仍保留了

\[
\frac{T^{(B)}}{3M_P^2}
\]

这一项，会使“冻结的辅助场”无法数值退回 `B` 支。

因此定义

\[
\Box\Phi_{\mathrm{flat}}

=
\frac{T^{(B)}}{3M_P^2},
\qquad
\Box\Phi_{\mathrm{corr}}
=
\Box\Phi-\Box\Phi_{\mathrm{flat}}.
\]

后续 `C` 支的空间应力投影与 `\Pi_\Phi` 演化都只使用
\[
\Box\Phi_{\mathrm{corr}}.
\]

## 2. 逐项退回检查

在标准 `\pm45^\circ` 双束高斯波包初值上，直接比较

\[
\Phi=1,\qquad \Pi_\Phi=0
\]

时 `B` 与 `C` 的有效源项，得到：

- `energy` 差值最大值：`0`
- `mom_x` 差值最大值：`0`
- `mom_z` 差值最大值：`0`
- `s_xx` 差值最大值：`0`
- `s_xz` 差值最大值：`0`
- `s_zz` 差值最大值：`0`
- `s_yy` 差值最大值：`0`

同时

\[
U(1)=0,\qquad U_\Phi(1)=0,\qquad \Box\Phi_{\mathrm{corr}}=0.
\]

因此此时 `C` 支的几何有效源与 `B` 支完全一致。

## 3. 时间推进检查

运行：

```text
python3 kg_examples/simulate_ab_adm_synchronous.py --steps 3 --dt 0.0005
python3 kg_examples/simulate_ab_adm_synchronous.py --steps 10 --dt 0.0005
```

结果表明：

- `C` 支的 `r3` 历史与 `B` 支逐步完全重合
- `C` 支的 Hamilton 约束残差历史与 `B` 支逐步完全重合
- `C` 支的 `rho` 最大值与 `B` 支完全重合
- `C` 支的 `phi` 始终保持 `1`

因此在当前统一四维 `ADM` 原型中，

\[
\boxed{
\Phi=1,\ \Pi_\Phi=0
\quad\Longrightarrow\quad
C\text{ 支严格退回 }B\text{ 支。}
}
\]

## 4. 这一步的意义

这一步只说明：

1. `C` 支数值器的平直辅助场归一化已经正确；
2. `C` 支的有效 Einstein 源投影不再带有错误的平直伪源；
3. 当前 `C` 支原型已具备做“非平凡 `\Phi` 初值”演化的最低一致性条件。

它还**不**意味着：

1. `C` 支已经进入与 `B` 明显分开的正式物理解；
2. 当前初值已经满足三支各自的 Hamilton / momentum 约束；
3. 三支真正同一物理初值的 full dynamics 比较已经完成。

## 5. 下一步

校准完成后，后续真正的问题已经收缩成两件：

1. 如何为 `C` 支构造非平凡而一致的 `\Phi` 初值；
2. 如何为三支求解各自的初值约束，使当前双束高斯波包不再从“平直几何 + 非零物质”的非约束态硬启动。
