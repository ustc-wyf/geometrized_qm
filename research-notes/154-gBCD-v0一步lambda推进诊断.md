# gBCD v0 一步 \(\lambda\) 推进诊断

日期：2026-05-07

## 目的

检验当前 equation-first gBCD v0 是否可以从已知两层
\[
\lambda_-,\lambda_0
\]
推出下一层 \(\lambda_+\)，从而替代“每个时间切片都用 A 支快照重新做三切片投影”的做法。

这里的 \(\lambda_I=(A,B,C,D)\) 是辅助各向异性应力张量
\[
\mathcal C_{\mu\nu}
=A\tilde g_{\mu\nu}+B u_\mu u_\nu+C r_\mu r_\nu+D u_{(\mu}r_{\nu)}
\]
的四个系数。

## 方程

一步原型固定 A 支背景上的 \(\rho,S,\tilde g,u,r\)，只求解 \(\lambda_+\)。

约束/目标如下：

1. 下一层代数场方程尽量满足：
   \[
   \lambda_I^+E^I_{\mu\nu}(t+\Delta t)
   \simeq
   \tilde G_{\mu\nu}(t+\Delta t)-\tilde T_{\mu\nu}(t+\Delta t)/M_P^2 .
   \]

2. 中心切片 full conservation 作为硬约束：
   \[
   \tilde\nabla^\mu(\lambda_I E^I_{\mu\nu})\big|_{t}=0 .
   \]

3. 在硬约束解空间中加入 v0 的弱辅助 action 代表元：
   \[
   w_0|\lambda_+|^2+w_t|\lambda_+-\lambda_{\rm target}|^2 .
   \]

测试了两种目标：

- `hold`：\(\lambda_{\rm target}=\lambda_0\)；
- `inertial`：\(\lambda_{\rm target}=2\lambda_0-\lambda_-\)。

## 数值结果

脚本：

- `kg_examples/test_gbcd_v0_one_step_lambda.py`

输出：

- `visualizations/equation_first_gbcd_v0_one_step_n96_taum3p5/`
- `visualizations/equation_first_gbcd_v0_one_step_n96_tau0/`
- `visualizations/equation_first_gbcd_v0_one_step_n96_taup3p5/`

参数：

- `full_resolution=96`
- `fit_region=trusted`
- `force_region=trusted`
- `norm_weight=1e-8`
- `time_weight=1e-5`
- 不使用 Q 门控

三组结果摘要：

| center tau | mode | hard conservation max | plus algebraic weighted mean | plus algebraic p50 | plus algebraic p95 | \(\lambda_+\) vs full projection p95 avg |
|---:|---|---:|---:|---:|---:|---:|
| -3.5 | hold | \(2.22\times10^{-16}\) | 0.16768 | 0.03041 | 0.60743 | 0.01872 |
| -3.5 | inertial | \(2.22\times10^{-16}\) | 0.16926 | 0.03275 | 0.60743 | 0.06529 |
| 0 | hold | \(4.86\times10^{-17}\) | 0.13704 | 0.01703 | 0.56949 | 0.00072 |
| 0 | inertial | \(5.09\times10^{-17}\) | 0.13683 | 0.01661 | 0.56949 | 0.00102 |
| +3.5 | hold | \(9.58\times10^{-16}\) | 0.16966 | 0.05651 | 0.68655 | 0.00078 |
| +3.5 | inertial | \(1.07\times10^{-15}\) | 0.17029 | 0.05642 | 0.68655 | 0.00356 |

其中：

- `hard conservation max` 是中心切片守恒硬约束残差的最大绝对值；
- `plus algebraic weighted mean` 是下一切片代数场方程的密度加权相对残差；
- `\(\lambda_+\) vs full projection p95 avg` 是一步预测的 \(\lambda_+\) 与完整三切片投影中 \(\lambda_+\) 的归一化 p95 差异。

## 结论

1. 硬守恒约束可以满足到机器精度。

2. 一步推进器能够基本复现完整三切片投影中的 \(\lambda_+\)：
   - tau=0 和 tau=+3.5 的 p95 平均差异小于 \(0.4\%\)；
   - tau=-3.5 的 `hold` 约 \(1.9\%\)，`inertial` 较差。

3. 但这不是完整动力学成功。固定 A 支背景、只更新 \(\lambda\) 时，下一切片代数场方程残差仍为 \(14\%\sim17\%\)，p95 可到 \(57\%\sim69\%\)。

4. 该结果说明：
   - 当前三切片投影中的大部分结构确实由 hard conservation 决定；
   - 但若 \(\rho,S,\tilde g\) 固定为 A 支背景，只靠 \(\lambda\) 不能让 D 支方程在下一切片继续成立；
   - 因此真正的 D 支推进必须允许物质场和/或几何场一起改变。

## 下一步

不要继续只优化 \(\lambda\)-only 推进器。

下一版最小一步系统应至少同时更新：

\[
(\delta\rho,\delta S,\delta\lambda_I)
\]

或者等价地更新守恒密度和相位梯度，使：

1. 中心或下一切片的 full conservation 硬满足；
2. 下一切片代数 gBCD 场方程低残差；
3. 物质部分仍满足变换后连续性方程和质量壳/测地线条件；
4. 与 A 支的差异只作为软目标，而不是硬约束。

这也是从“投影 A 快照”转向“真正 D 支时间演化”的必要步骤。
