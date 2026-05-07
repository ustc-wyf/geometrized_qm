# 101 D 支标量-张量参考校正 IMEX 原型

## 目标

把当前更合理的新路线先落成一个**最小可运行原型**：

```text
D 支
+ 带符号曲率代理 chi = R~_ref
+ 标量-张量重写
+ 平直量子力学参考解预测
+ 局部 IMEX 校正种子
+ 局部加密指标
```

这里的重点不是一步到位实现完整全几何积分器，而是先确认：

1. `D` 支能否摆脱 `phi=f_R` 的 flat-limit 病态；
2. 当前基准上哪些区域真正决定 stiff 时间步；
3. 参考解校正的第一层对象是否数值上稳定、可解释。

## 变量选择

对

```text
f_D(R~) = tanh(ell^2 R~) / ell^2
```

取带符号曲率代理

```text
chi := R~_ref
```

再定义

```text
y = tanh(ell^2 chi)
phi = f_R = 1 - y^2 = sech^2(ell^2 chi)
U(chi) = chi * phi - f_D(chi)
```

这样做的目的，是避免直接把 `phi=f_R` 当主变量时在 `phi -> 1` 附近出现的病态压缩。

## 当前实现

脚本：

```text
kg_examples/prototype_d_reference_imex.py
```

它做了四件事：

1. 用当前局域双高斯平直量子力学基准，在 `t = 0, 8, 16` 生成 native `tilde g_ref[rho,S]`；
2. 计算 `D` 支 metric-only 场方程残差与 `2+1d` 迹残差；
3. 生成带符号 `chi_ref, phi_ref` 与一个局部 IMEX 校正种子 `delta_chi_seed`；
4. 生成局部加密掩膜。

## 区域定义

### 主支撑区

```text
rho > 1e-3 * rho_max
```

### 高密度区

```text
rho > 1e-1 * rho_max
```

### 近过渡带

```text
| |ell^2 R~| - 1 | <= 1
```

这里“近过渡带”只是第一版局部加密指标，不是严格理论上的过渡层定义。

### 局部加密掩膜

第一版取

```text
support
&
(
  near_transition
  or residual_norm > p99
  or |grad phi| > p99
  or derivative_norm > p99
)
```

也就是说：只加密真正的尾部尖点和过渡带附近，而不再像旧的粗糙阈值那样把整个主支撑区都卷进去。

## 第一轮结果（96x96）

输出目录：

```text
visualizations/d_scalar_tensor_reference_imex_prototype_96/
```

摘要：

```text
summary.json
```

### 观察 1

原型本身能稳定跑完，不再像之前的纯 metric 显式推进那样一上来就炸。

### 观察 2

主支撑区里的 `D` 支残差代表性量仍然维持在我们前面看到的“小 bulk + 薄层尖峰”结构：

- `t=8` 时 `p95` 仍然小；
- `t=16` 时 `p95` 仍然小，但 `abs_max` 由薄层尖峰控制。

### 观察 3

新的局部加密掩膜不再覆盖整个主支撑区，第一轮大约只覆盖：

- `t=0`：约 `1.1%`
- `t=8`：约 `1.0%`
- `t=16`：约 `1.7%`

这说明：

```text
最大局部刚性确实来自极少数局部区域，
而不是整个 bulk。
```

### 观察 4

`delta_chi_seed` 的代表性尺度在主支撑区里是 `O(1)`，没有再出现平直极限下的病态爆大。

这说明：

```text
用带符号 chi 而不是直接用 phi，
至少在“第一层校正种子”这一步已经明显更健康。
```

## 当前结论

这个原型还不是完整全几何积分器，但它已经确认了新路线里的三个关键判断：

1. `D` 支比 `C/B` 更适合作为下一阶段基底；
2. `phi=f_R` 不是好主变量，带符号 `chi` 更自然；
3. 旧求解器真正被拖死的是极小比例的局部 stiff 区，而不是整个参考 bulk。

## 下一步

下一步不应回头修旧的 pure metric 显式推进器，而应继续沿这条路：

```text
D 支
+ chi/y 变量
+ 参考解预测
+ 局部 IMEX 校正
+ 局部网格加密
```

更具体地说，下一步应把当前的 `delta_chi_seed` 与局部加密掩膜接进真正的时间推进器，先做一个：

```text
reference-predictor + local implicit chi corrector
```

的时间步原型。
