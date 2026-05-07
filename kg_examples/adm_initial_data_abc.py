from __future__ import annotations

import numpy as np

from adm_matter_2p1 import matter_rhs
from adm_sources_abc import (
    a_branch_sources_from_complex,
    b_branch_sources_from_ns,
    complex_initial_data_from_rho_s,
)
from adm_ykilling_geometry import beta_hessian_and_spatial_ricci


def ddx(f: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def ddz(f: np.ndarray, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out


def a_complex_state_from_rho_s_time_symmetric(
    rho: np.ndarray,
    s: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    """
    从同一组 (rho, S, g) 构造 A 支的复标量初值。

    这里采用时间对称几何假设：
    - N = 1
    - N^x = N^z = 0
    - K_ij = 0, K_beta = 0

    在这一假设下，可由共同物质守恒系统计算
    S_t, n_t, E, rho_t，
    再重建
    phi_1, phi_2, Pi_1, Pi_2。
    """
    lapse = np.ones_like(rho)
    zero = np.zeros_like(rho)
    mrhs = matter_rhs(rho, s, lapse, zero, zero, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    sx = mrhs["sx"]
    sz = mrhs["sz"]
    e = mrhs["energy"]
    s_t = mrhs["s_t"]
    n_t = mrhs["n_t"]
    sqrt_h = mrhs["sqrt_h"]

    s_t_x = ddx(s_t, dx)
    s_t_z = ddz(s_t, dz)
    e_t = (sx * s_t_x + sz * s_t_z) / np.maximum(e, 1.0e-12)
    rho_t = (n_t - np.exp(beta) * sqrt_h * rho * e_t) / np.maximum(np.exp(beta) * sqrt_h * e, 1.0e-12)

    amp = np.sqrt(np.clip(rho, 1.0e-12, None))
    phi1 = amp * np.cos(s)
    phi2 = amp * np.sin(s)
    amp_t = 0.5 * rho_t / amp
    pi1 = amp_t * np.cos(s) - amp * np.sin(s) * s_t
    pi2 = amp_t * np.sin(s) + amp * np.cos(s) * s_t

    return {
        "phi1": phi1,
        "phi2": phi2,
        "pi1": pi1,
        "pi2": pi2,
        "rho_t": rho_t,
        "s_t": s_t,
        "energy": e,
    }


def solve_time_symmetric_conformal_geometry(
    energy: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
    max_iters: int = 2000,
    omega_tol: float = 1.0e-10,
    relax: float = 0.6,
) -> dict[str, np.ndarray | float | int]:
    """
    在 ansatz
      h_ij = e^{2 omega} delta_ij,  h_xz = 0,  beta = 0,
      K_ij = 0, K_beta = 0
    下，解时间对称 Hamilton 约束
      R^(3) = 2 E / M_P^2.

    对 beta=0 的当前 ansatz，有
      R^(3) = -2 e^{-2 omega} Delta omega,
    从而得到非线性 Poisson 型方程
      Delta omega = - e^{2 omega} E / M_P^2.

    边界条件取 omega = 0。
    """
    if not np.allclose(dx, dz, rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("当前实现假设 dx = dz。")

    nx, nz = energy.shape
    h2 = dx * dx
    omega = np.zeros_like(energy)
    rhs_scale = -energy / (mp * mp)

    for it in range(max_iters):
        source = rhs_scale * np.exp(np.clip(2.0 * omega, -20.0, 20.0))
        omega_new = omega.copy()
        omega_new[1:-1, 1:-1] = 0.25 * (
            omega[2:, 1:-1]
            + omega[:-2, 1:-1]
            + omega[1:-1, 2:]
            + omega[1:-1, :-2]
            - h2 * source[1:-1, 1:-1]
        )
        omega_new = relax * omega_new + (1.0 - relax) * omega
        omega_new[0, :] = 0.0
        omega_new[-1, :] = 0.0
        omega_new[:, 0] = 0.0
        omega_new[:, -1] = 0.0
        delta = float(np.max(np.abs(omega_new - omega)))
        omega = omega_new
        if delta < omega_tol:
            break

    conformal = np.exp(np.clip(2.0 * omega, -20.0, 20.0))
    h_xx = conformal.copy()
    h_xz = np.zeros_like(conformal)
    h_zz = conformal.copy()
    beta = np.zeros_like(conformal)
    geo = beta_hessian_and_spatial_ricci(h_xx, h_xz, h_zz, beta, dx, dz)
    ham_res = geo["r3_scalar"] - 2.0 * energy / (mp * mp)
    return {
        "omega": omega,
        "h_xx": h_xx,
        "h_xz": h_xz,
        "h_zz": h_zz,
        "beta": beta,
        "r3_scalar": geo["r3_scalar"],
        "hamiltonian_residual": ham_res,
        "iterations": it + 1,
        "final_delta": delta,
    }


def solve_average_momentum_extrinsic_curvature(
    mom_x: np.ndarray,
    mom_z: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
    max_iters: int = 2000,
    tol: float = 1.0e-10,
    relax: float = 0.7,
) -> dict[str, np.ndarray | float | int]:
    """
    在平坦二维背景近似下，用最简单的 York 型向量势 ansatz 解动量约束：

      K = 0,
      K_beta = 0,
      K_ij = d_i W_j + d_j W_i - delta_ij d_k W_k,
      Delta W_i = mom_i / M_P^2.

    这是一个近似共同初值构造，用于把非零动量密度纳入同一组初始几何。
    """
    if not np.allclose(dx, dz, rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("当前实现假设 dx = dz。")

    h2 = dx * dx
    wx = np.zeros_like(mom_x)
    wz = np.zeros_like(mom_z)

    rhs_x = mom_x / (mp * mp)
    rhs_z = mom_z / (mp * mp)

    delta = np.inf
    for it in range(max_iters):
        wx_new = wx.copy()
        wz_new = wz.copy()
        wx_new[1:-1, 1:-1] = 0.25 * (
            wx[2:, 1:-1] + wx[:-2, 1:-1] + wx[1:-1, 2:] + wx[1:-1, :-2] - h2 * rhs_x[1:-1, 1:-1]
        )
        wz_new[1:-1, 1:-1] = 0.25 * (
            wz[2:, 1:-1] + wz[:-2, 1:-1] + wz[1:-1, 2:] + wz[1:-1, :-2] - h2 * rhs_z[1:-1, 1:-1]
        )
        wx_new = relax * wx_new + (1.0 - relax) * wx
        wz_new = relax * wz_new + (1.0 - relax) * wz
        for arr in (wx_new, wz_new):
            arr[0, :] = 0.0
            arr[-1, :] = 0.0
            arr[:, 0] = 0.0
            arr[:, -1] = 0.0
        delta = max(float(np.max(np.abs(wx_new - wx))), float(np.max(np.abs(wz_new - wz))))
        wx = wx_new
        wz = wz_new
        if delta < tol:
            break

    div_w = ddx(wx, dx) + ddz(wz, dz)
    k_xx = 2.0 * ddx(wx, dx) - div_w
    k_zz = 2.0 * ddz(wz, dz) - div_w
    k_xz = ddx(wz, dx) + ddz(wx, dz)

    return {
        "wx": wx,
        "wz": wz,
        "k_xx": k_xx,
        "k_xz": k_xz,
        "k_zz": k_zz,
        "iterations": it + 1,
        "final_delta": delta,
    }


def branch_a_sources_from_rho_s_time_symmetric(
    rho: np.ndarray,
    s: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    a_init = a_complex_state_from_rho_s_time_symmetric(rho, s, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    src = a_branch_sources_from_complex(
        a_init["phi1"],
        a_init["phi2"],
        a_init["pi1"],
        a_init["pi2"],
        h_xx,
        h_xz,
        h_zz,
        beta,
        mass,
        dx,
        dz,
    )
    return {**a_init, **src}


def branch_b_sources_from_rho_s(
    rho: np.ndarray,
    s: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    one = np.ones_like(rho)
    zero = np.zeros_like(rho)
    return b_branch_sources_from_ns(rho, s, one, zero, zero, h_xx, h_xz, h_zz, beta, mass, dx, dz)


def compare_a_b_sources_common_rho_s(
    rho: np.ndarray,
    s: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> dict[str, float]:
    a_src = branch_a_sources_from_rho_s_time_symmetric(rho, s, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    b_src = branch_b_sources_from_rho_s(rho, s, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    out: dict[str, float] = {}
    for key in ("energy", "mom_x", "mom_z", "s_xx", "s_xz", "s_zz", "s_yy_reduced"):
        a = a_src[key]
        b = b_src[key]
        out[f"{key}_max_abs_diff"] = float(np.max(np.abs(a - b)))
        out[f"{key}_rel_l1_diff"] = float(np.sum(np.abs(a - b)) / max(np.sum(np.abs(a)), 1.0e-12))
    return out
