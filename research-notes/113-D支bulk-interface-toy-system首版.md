# D 支 bulk-interface toy system 首版

时间：2026-05-01

## 目的

用户要求继续执行：

```text
以 ell=30 写 bulk-interface toy system：
EH-like bulk、saturated fixed/extended bulk、y=±1 interface jump 源项。
```

本笔记记录第一版可运行 toy system。

## 脚本与输出

脚本：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/kg_examples/prototype_d_bulk_interface_toy.py
```

输出：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240/summary.json
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240/d_bulk_interface_toy_ell30_n240_t16.png
```

设置：

- 分支：`D`
- `ell=30`
- `t=16`
- 分辨率：`240x240`
- \(M_P=300\)

## 动态界面定义

令

\[
y(t,x,z)=\ell^2\tilde R(t,x,z).
\]

动态界面为

\[
y=+1,\qquad y=-1.
\]

当前脚本在固定参考切片 \(t=16\) 上生成界面；全演化时应在每个时间步从当前 \(\tilde R(t,x,z)\) 重新生成。

## 区域划分

### EH-like bulk

\[
\rho>10^{-3}\rho_{\max},\qquad |y|<1.
\]

第一版策略：

```text
用 D/scalar-tensor bulk 方程或其弱曲率近似推进。
```

### saturated bulk

\[
\rho>10^{-3}\rho_{\max},\qquad |y|>1.
\]

第一版策略：

```text
不作为普通 Cauchy PDE 自主演化；
先用固定参考延拓或最小曲率延拓。
```

### interface

\[
y=\pm1.
\]

第一版策略：

```text
由法向短线积分给出 jump/matching 源项。
```

## 结果摘要

### 区域计数

- 主支撑区格点数：`5837`
- EH-like bulk：`9`
- saturated bulk：`5828`
- interface lines：`420`
- interface pixels：`171`
- interface near support：`712`
- derivative hotspots：`59`

### 支撑区比例

- EH-like bulk：`0.00154`
- saturated bulk：`0.99846`
- interface near：`0.12198`
- derivative hotspot near interface：`0.40678`

这里 `interface near` 指界面像素向外膨胀 `3` 个网格格点后的邻域。

### 界面源强度

界面短线积分统计：

- signed integral p95：
  \[
  8.07\times10^3
  \]
- absolute integral p95：
  \[
  4.23\times10^4
  \]
- global signed/abs ratio：
  \[
  5.31\times10^{-2}.
  \]

这说明界面层内部有强正负抵消；净 jump 比绝对层强度小得多。

## 图示

图文件：

```text
/Users/wangyunfei/Desktop/量子势几何化_codex/visualizations/d_bulk_interface_toy_ell30_t16_240/d_bulk_interface_toy_ell30_n240_t16.png
```

四个面板分别为：

1. \(\log_{10}\rho\) 与动态界面 \(|y|=1\)；
2. bulk-interface 区域分类：
   - `0`: EH-like bulk
   - `1`: saturated bulk
   - `2`: interface near
3. 原始导数项范数，lime 为 top 1% hotspot，cyan 为 interface near；
4. rasterized interface absolute source proxy。

## 解读

### 1. 当前参考切片几乎全部处于 saturated bulk

`ell=30` 下，主支撑区中 `99.846%` 已经满足 \(|y|>1\)。

这说明如果用 interface 路线，真正需要普通 bulk PDE 推进的 EH-like 区很少；主要问题变成：

```text
saturated bulk 如何固定/延拓；
interface 如何匹配。
```

### 2. 界面邻域抓住了相当一部分 stiff 区，但不是全部

interface near 覆盖主支撑区约 `12.2%`，覆盖 derivative top 1% hotspot 的约 `40.7%`。

这说明：

```text
界面源确实抓住了一部分 stiff 结构；
但还有一部分导数项热点不在当前 y=±1 的简单界面邻域里。
```

可能原因：

- hotspot 来自更宽的粗过渡带；
- 需要同时追踪多个 level set，例如 \(|y|=0.5,1,2\)；
- 界面邻域半径还不够；
- 部分热点来自低密度区或几何导数噪声，而不是主界面。

### 3. 第一版 toy system 可运行，但还不是完整演化器

当前已经形成三个可计算对象：

1. EH-like bulk mask；
2. saturated bulk mask；
3. interface jump source proxy。

但还没有真正推进 \(\tilde g(t)\)。下一步要把这些对象接入一个简化演化规则。

## 下一步

建议下一步做两个分支：

1. 改进界面定义：
   - 同时测试 \(|y|=0.5,1,2\) 多界面；
   - 比较 hotspot 覆盖率是否显著提高。
2. 写 toy evolution：
   - EH-like bulk：用弱曲率 bulk 方程；
   - saturated bulk：固定为参考延拓；
   - interface：加入 jump 源项；
   - 检查一小步预测是否比原始 reference residual 更稳定。

当前结论：

\[
\boxed{
bulk-interface toy system 已经落成第一版；
ell=30 可作为首个工作参数。
}
\]

但还需要改进界面定义，才能覆盖更多 stiff hotspot。
