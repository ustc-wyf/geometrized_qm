# \(M_P^2\) 有效刚度与过渡层 jump 首个数值诊断

时间：2026-05-01

## 目的

承接用户追问，本轮不再只看 \(f_R\)，而是检查：

1. \(M_P^2 f_R\) 是否仍足够大，能否解释为“几何被钉住”；
2. 过渡层是否能在当前网格上采样；
3. 局部 jump 代理量是否显示强薄层；
4. C/D 的退化区能否类比普通平直量子力学中的 \(g=\eta\) 选解原则。

## 脚本与输出

新脚本：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/diagnose_cd_transition_jump_mp2.py
```

输出：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_transition_jump_mp2_t16_160/summary.json
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/c_transition_jump_mp2_t16_160/summary.json
```

设置：

- 时间：`t=16`
- 网格：`160x160`
- 窗口：左弧窗口
  \[
  -7.4\le x\le -2.2,\qquad 2.6\le z\le6.2
  \]
- 分支：`C,D`
- `ell = 10,30,100,300`
- \(M_P=300\)

## 定义

### 支撑区

\[
\rho > 10^{-3}\rho_{\max}.
\]

### 严格过渡层

\[
0.5\le |\ell^2\tilde R|\le2.
\]

### 粗过渡带

\[
1\le |\ell^2\tilde R|\le100.
\]

粗过渡带只用于定位，不等同于真正薄层。

### 有效刚度

度规主部前的有效系数是

\[
M_{\rm eff}^2=M_P^2f_R.
\]

因此 \(f_R\ll1\) 不自动意味着主部无效；真正要看 \(M_P^2f_R\) 与物质源、代数项和导数项的相对大小。

### jump 代理量

令

\[
\phi=f_R,\qquad y=\ell^2\tilde R.
\]

用空间切片上的

\[
n_i=\frac{\partial_i y}{|\nabla y|}
\]

作为过渡层法向，计算：

\[
\partial_n\phi,
\qquad
\partial_n^2\phi.
\]

局部宽度代理量为

\[
\Delta n_{\rm proxy}=\frac{1}{|\nabla y|}.
\]

局部 signed jump 代理量：

\[
J_{\rm proxy}=\partial_n^2\phi\,\Delta n_{\rm proxy}.
\]

局部绝对层强度：

\[
A_{\rm layer}=|\partial_n^2\phi|\,\Delta n_{\rm proxy}.
\]

注意：这不是完整线积分，只是第一版局部诊断。真正 jump 分析还需要沿法向穿层积分。

## 结果摘要

### D 支

左弧窗口，`t=16`：

| ell | strict count | broad count | p95 \(M_P^2f_R\) | p95 \(T\) | broad p95 derivative/T | strict p95 \(A_{\rm layer}\) | signed/abs ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 31 | 188 | 5.77e4 | 3.23e1 | 4.31e7 | 1.31 | 0.128 |
| 30 | 5 | 139 | 1.56 | 3.23e1 | 9.11e6 | 5.07e-2 | 0.204 |
| 100 | 0 | 20 | 0 | 3.23e1 | 8.51e7 | 0 | 0 |
| 300 | 0 | 3 | 0 | 3.23e1 | 5.33e-5 | 0 | 0 |

### C 支

左弧窗口，`t=16`：

| ell | strict count | broad count | p95 \(M_P^2f_R\) | p95 \(T\) | broad p95 derivative/T | strict p95 \(A_{\rm layer}\) | signed/abs ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 31 | 188 | 5.01e4 | 3.23e1 | 3.36e7 | 1.17 | 0.076 |
| 30 | 5 | 139 | 3.62e2 | 3.23e1 | 3.16e7 | 4.64e-2 | 0.188 |
| 100 | 0 | 20 | 2.74e-1 | 3.23e1 | 4.82e8 | 0 | 0 |
| 300 | 0 | 3 | 3.76e-4 | 3.23e1 | 2.29e6 | 0 | 0 |

## 解读 1：\(M_P^2 f_R\) 改变了退化判断

如果只看 \(f_R\)，`ell=10` 已经很小。

但乘上 \(M_P^2=9\times10^4\) 后，`ell=10` 的左弧窗口中

\[
M_P^2f_R
\]

仍远大于 \(T\) 的 p95 尺度。

这意味着 `ell=10` 更像“几何仍被有效钉住”，不是主部已经完全失效。

`ell=30` 开始出现差异：

- C 支 p95 \(M_P^2f_R\approx3.62\times10^2\)，仍大于 \(T\) p95；
- D 支 p95 \(M_P^2f_R\approx1.56\)，已经低于 \(T\) p95。

所以 D 支因指数压制更快进入真正退化区。

`ell=100,300` 下，C/D 都明显进入强退化或数值下溢区。

## 解读 2：过渡层导数项相对 \(T\) 极大

在粗过渡带中，`ell=10,30,100` 的

\[
\frac{M_P^2\|\text{derivative term}\|}{\|\tilde T\|}
\]

p95 可达 \(10^6\sim10^8\) 量级。

这说明过渡带里的导数项不是小修正；它正是旧全动力学推进困难的主要 stiffness 来源。

`ell=300` 的 D 支这个比例很小，不应立即解读为安全，因为严格层和粗带采样点都很少，且 \(f_R,f_{RR}\) 已经强烈下溢。它更像“当前网格没有解析到层”。

## 解读 3：当前 jump 诊断还不是完整 jump

对 `ell=10`，严格层在窗口中有 31 个采样点；对 `ell=30` 只有 5 个；对 `ell=100,300` 为 0。

因此：

\[
\boxed{
ell=100,300 的真正严格过渡层在 160x160 下仍未解析。
}
\]

当前只能说：

- `ell=10` 的层强度可初步估计；
- `ell=30` 的层强度只能粗略参考；
- `ell=100,300` 的 jump 安全性不能由当前网格判断。

`signed/abs ratio` 在 `ell=10,30` 大约 `0.08` 到 `0.20`，说明局部 signed jump 代理量存在明显正负抵消。

但绝对层强度并不为零，因此不能说过渡层没有物理/数值效应。

## 解读 4：用户的“固定度规”思路是可行但需要提升为选解原则

标准平直量子力学里选 \(g=\eta\) 的理由是：

\[
M_P^2G_{\mu\nu}=T_{\mu\nu},
\qquad
G_{\mu\nu}=O(M_P^{-2}),
\]

再加上渐近平直边界条件和无自由引力波初值。

C/D 饱和区如果 \(f_R\to0\)，并不自动给普通 Cauchy 演化。但可以尝试定义一个类似的选解原则：

```text
在 saturated bulk 中不独立演化 tilde g；
tilde g 由参考平直度规、最小曲率原则或边界匹配条件钉住；
EH-like 区负责正常几何响应；
transition layer 负责匹配两个区域。
```

这会把 C/D 从普通 metric `f(R)` 演化理论，变成一种 singular-limit / matched-boundary effective theory。

这不是坏事，但必须明确写进理论定义里。

## 当前结论

1. 用户关于 \(M_P^2\) 的提醒是正确的：退化判断不能只看 \(f_R\)，必须看 \(M_P^2 f_R\)。
2. 在当前参数 \(M_P=300\) 下，`ell=10` 仍像有效钉住区，`ell=30` 是过渡，`ell=100,300` 才是明显强退化。
3. 过渡带导数项相对物质源极大，旧全显式动力学失稳并不奇怪。
4. 当前网格不能判断 `ell=100,300` 的 jump 安全性，因为严格层没有采样。
5. 下一步若继续严谨推进，应做真正的局部法向线积分，而不是只看局部 proxy。

## 下一步

建议下一步做两个小闭环：

1. 在左弧窗口上构造沿 \(\nabla(\ell^2\tilde R)\) 的局部法向短线，直接积分
   \[
   [\partial_n f_R],\quad
   \int \partial_n^2 f_R\,dn,\quad
   \int|\partial_n^2 f_R|dn.
   \]
2. 写出 saturated bulk 的候选选解原则：
   \[
   \tilde g_{\rm sat}=\operatorname*{argmin}_{\tilde g}
   \int_{\rm sat}|\tilde R|^2
   \]
   或
   \[
   \tilde g_{\rm sat}=\text{reference extension from EH-like boundary}.
   \]

然后检查这些选解原则是否与物质守恒和过渡层 jump 条件相容。
