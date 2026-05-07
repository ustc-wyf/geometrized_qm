from __future__ import annotations

import numpy as np

from adm_matter_2p1 import matter_observables
from adm_ykilling_geometry import inverse_2metric


def u_of_phi(phi: np.ndarray, ell: float) -> np.ndarray:
    phi_safe = np.clip(phi, 1.0e-6, 1.0)
    r_phi = np.sqrt(np.maximum(np.power(phi_safe, -2.0 / 3.0) - 1.0, 0.0)) / (ell * ell)
    return r_phi * (phi_safe - np.power(phi_safe, 1.0 / 3.0))


def u_phi_of_phi(phi: np.ndarray, ell: float) -> np.ndarray:
    phi_safe = np.clip(phi, 1.0e-6, 1.0)
    return np.sqrt(np.maximum(np.power(phi_safe, -2.0 / 3.0) - 1.0, 0.0)) / (ell * ell)


def a_branch_sources_from_complex(
    phi1: np.ndarray,
    phi2: np.ndarray,
    pi1: np.ndarray,
    pi2: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    from adm_ykilling_geometry import ddx, ddz

    hxx_inv, hxz_inv, hzz_inv, det_h = inverse_2metric(h_xx, h_xz, h_zz)

    phi1_x = ddx(phi1, dx)
    phi1_z = ddz(phi1, dz)
    phi2_x = ddx(phi2, dx)
    phi2_z = ddz(phi2, dz)

    grad1_sq = hxx_inv * phi1_x * phi1_x + 2.0 * hxz_inv * phi1_x * phi1_z + hzz_inv * phi1_z * phi1_z
    grad2_sq = hxx_inv * phi2_x * phi2_x + 2.0 * hxz_inv * phi2_x * phi2_z + hzz_inv * phi2_z * phi2_z

    energy = pi1 * pi1 + pi2 * pi2 + grad1_sq + grad2_sq + mass * mass * (phi1 * phi1 + phi2 * phi2)
    mom_x = -2.0 * (pi1 * phi1_x + pi2 * phi2_x)
    mom_z = -2.0 * (pi1 * phi1_z + pi2 * phi2_z)

    common = pi1 * pi1 + pi2 * pi2 - grad1_sq - grad2_sq - mass * mass * (phi1 * phi1 + phi2 * phi2)

    s_xx = 2.0 * (phi1_x * phi1_x + phi2_x * phi2_x) + h_xx * common
    s_xz = 2.0 * (phi1_x * phi1_z + phi2_x * phi2_z) + h_xz * common
    s_zz = 2.0 * (phi1_z * phi1_z + phi2_z * phi2_z) + h_zz * common
    return {
        "energy": energy,
        "mom_x": mom_x,
        "mom_z": mom_z,
        "s_xx": s_xx,
        "s_xz": s_xz,
        "s_zz": s_zz,
        "s_yy_reduced": common,
        "det_h": det_h,
    }


def complex_initial_data_from_rho_s(
    rho: np.ndarray,
    s_phase: np.ndarray,
) -> dict[str, np.ndarray]:
    amp = np.sqrt(np.clip(rho, 0.0, None))
    return {
        "phi1": amp * np.cos(s_phase),
        "phi2": amp * np.sin(s_phase),
    }


def b_branch_sources_from_ns(
    rho: np.ndarray,
    s: np.ndarray,
    lapse: np.ndarray,
    shift_x: np.ndarray,
    shift_z: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    obs = matter_observables(
        rho=rho,
        s=s,
        lapse=lapse,
        shift_x=shift_x,
        shift_z=shift_z,
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        mass=mass,
        dx=dx,
        dz=dz,
    )
    rho_eff = rho
    energy = obs["energy"]
    sx = obs["sx"]
    sz = obs["sz"]

    return {
        "energy": 2.0 * rho_eff * energy * energy,
        "mom_x": 2.0 * rho_eff * energy * sx,
        "mom_z": 2.0 * rho_eff * energy * sz,
        "s_xx": 2.0 * rho_eff * sx * sx,
        "s_xz": 2.0 * rho_eff * sx * sz,
        "s_zz": 2.0 * rho_eff * sz * sz,
        "s_yy_reduced": np.zeros_like(rho_eff),
        "rho": rho_eff,
        "phase_t": obs["s_t"],
    }


def c_branch_effective_sources(
    rho: np.ndarray,
    s: np.ndarray,
    phi: np.ndarray,
    pi_phi: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    k_xx: np.ndarray,
    k_xz: np.ndarray,
    k_zz: np.ndarray,
    k_beta: np.ndarray,
    lapse: np.ndarray,
    shift_x: np.ndarray,
    shift_z: np.ndarray,
    mass: float,
    mp: float,
    ell: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    from adm_ykilling_geometry import inverse_2metric, scalar_hessian_ykilling, scalar_laplacian_ykilling

    bsrc = b_branch_sources_from_ns(rho, s, lapse, shift_x, shift_z, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    hxx_inv, hxz_inv, hzz_inv, _ = inverse_2metric(h_xx, h_xz, h_zz)
    k_trace = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz + k_beta

    hphi = scalar_hessian_ykilling(h_xx, h_xz, h_zz, phi, dx, dz)
    lap3_phi = scalar_laplacian_ykilling(h_xx, h_xz, h_zz, beta, phi, dx, dz)
    beta_x = ddx_noimport(beta, dx)
    beta_z = ddz_noimport(beta, dz)
    grad_beta_dot_grad_phi = (
        hxx_inv * beta_x * hphi["fx"]
        + hxz_inv * (beta_x * hphi["fz"] + beta_z * hphi["fx"])
        + hzz_inv * beta_z * hphi["fz"]
    )

    u = u_of_phi(phi, ell)
    u_phi = u_phi_of_phi(phi, ell)
    trace_b = 2.0 * mass * mass * bsrc["rho"]
    box_phi = trace_b / (3.0 * mp * mp) - (2.0 * u - phi * u_phi) / 3.0
    box_phi_flat = trace_b / (3.0 * mp * mp)
    box_phi_corr = box_phi - box_phi_flat

    a_e = lap3_phi + k_trace * pi_phi - 0.5 * u

    kx_mixed_x = hxx_inv * k_xx + hxz_inv * k_xz
    kx_mixed_z = hxz_inv * k_xx + hzz_inv * k_xz
    kz_mixed_x = hxx_inv * k_xz + hxz_inv * k_zz
    kz_mixed_z = hxz_inv * k_xz + hzz_inv * k_zz

    p_phi_x = -(ddx_noimport(pi_phi, dx) + kx_mixed_x * hphi["fx"] + kx_mixed_z * hphi["fz"])
    p_phi_z = -(ddz_noimport(pi_phi, dz) + kz_mixed_x * hphi["fx"] + kz_mixed_z * hphi["fz"])

    a_s_xx = hphi["hess_xx"] + k_xx * pi_phi + h_xx * box_phi_corr + 0.5 * u * h_xx
    a_s_xz = hphi["hess_xz"] + k_xz * pi_phi + h_xz * box_phi_corr + 0.5 * u * h_xz
    a_s_zz = hphi["hess_zz"] + k_zz * pi_phi + h_zz * box_phi_corr + 0.5 * u * h_zz
    a_s_yy_reduced = grad_beta_dot_grad_phi + k_beta * pi_phi + box_phi_corr + 0.5 * u

    phi_safe = np.clip(phi, 1.0e-6, None)
    energy_eff = bsrc["energy"] / phi_safe + (mp * mp / phi_safe) * a_e
    mom_x_eff = bsrc["mom_x"] / phi_safe + (mp * mp / phi_safe) * p_phi_x
    mom_z_eff = bsrc["mom_z"] / phi_safe + (mp * mp / phi_safe) * p_phi_z
    s_xx_eff = bsrc["s_xx"] / phi_safe + (mp * mp / phi_safe) * a_s_xx
    s_xz_eff = bsrc["s_xz"] / phi_safe + (mp * mp / phi_safe) * a_s_xz
    s_zz_eff = bsrc["s_zz"] / phi_safe + (mp * mp / phi_safe) * a_s_zz
    s_yy_eff_reduced = bsrc["s_yy_reduced"] / phi_safe + (mp * mp / phi_safe) * a_s_yy_reduced

    return {
        "energy": energy_eff,
        "mom_x": mom_x_eff,
        "mom_z": mom_z_eff,
        "s_xx": s_xx_eff,
        "s_xz": s_xz_eff,
        "s_zz": s_zz_eff,
        "s_yy_reduced": s_yy_eff_reduced,
        "rho": bsrc["rho"],
        "phase_t": bsrc["phase_t"],
        "trace_b": trace_b,
        "u": u,
        "u_phi": u_phi,
        "box_phi_corr": box_phi_corr,
        "box_phi_flat": box_phi_flat,
        "lap3_phi": lap3_phi,
        "k_trace": k_trace,
    }


def ddx_noimport(f: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def ddz_noimport(f: np.ndarray, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out
