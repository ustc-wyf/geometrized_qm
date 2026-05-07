from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from calibrate_exact_kg_ps_fd_absorbing import (
    apply_boundary,
    boundary_exact,
    make_sponge,
    spectral_lap,
)
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz


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


def d2x(f: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dx * dx)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def d2z(f: np.ndarray, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / (dz * dz)
    out[:, 0] = out[:, 1]
    out[:, -1] = out[:, -2]
    return out


def reconstruct_from_psi(psi_prev: np.ndarray, psi_now: np.ndarray, psi_next: np.ndarray, dx: float, dz: float, dt: float, m: float):
    rho = np.abs(psi_now) ** 2
    rho_safe = np.maximum(rho, 1e-12)
    sqrt_rho = np.sqrt(rho_safe)
    sqrt_rho_prev = np.sqrt(np.maximum(np.abs(psi_prev) ** 2, 1e-12))
    sqrt_rho_next = np.sqrt(np.maximum(np.abs(psi_next) ** 2, 1e-12))

    dpsi_t = (psi_next - psi_prev) / (2.0 * dt)
    dpsi_x = ddx(psi_now, dx)
    dpsi_z = ddz(psi_now, dz)

    s = np.unwrap(np.unwrap(np.angle(psi_now), axis=0), axis=1)
    s_t = np.imag(np.conj(psi_now) * dpsi_t) / rho_safe
    s_x = np.imag(np.conj(psi_now) * dpsi_x) / rho_safe
    s_z = np.imag(np.conj(psi_now) * dpsi_z) / rho_safe

    X = s_t * s_t - s_x * s_x - s_z * s_z
    d2_sqrt_dt2 = (sqrt_rho_next - 2.0 * sqrt_rho + sqrt_rho_prev) / (dt * dt)
    d2_sqrt_dx2 = d2x(sqrt_rho, dx)
    d2_sqrt_dz2 = d2z(sqrt_rho, dz)
    Q = (d2_sqrt_dt2 - d2_sqrt_dx2 - d2_sqrt_dz2) / np.maximum(sqrt_rho, 1e-12)

    Xt = np.maximum(X, 1e-12)
    sqrt_det_ratio = np.sqrt(np.maximum(Xt / (m * m), 1e-12))
    rho_tilde = rho * sqrt_det_ratio

    return {
        "rho": rho,
        "S": s,
        "X": X,
        "Q": Q,
        "rho_tilde": rho_tilde,
        "sqrt_det_ratio": sqrt_det_ratio,
    }


def run(outdir: str | Path, alpha: float = 0.5):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    inner = 10.0
    outer = 20.0
    params = Exact2p1Params(alpha=alpha, x_half_range=outer, z_half_range=outer, nx=181, nz=181, nkx=61, nkz=61)
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    dt = 0.025
    t_final = params.overlap_time
    steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, params)
    _, _, _, rho_exact = integrate_xz(alpha, t_final, params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)

    for n in range(steps):
        t_next = (n + 1) * dt
        psi_b = boundary_exact(alpha, t_next, x, z, params)
        rhs = spectral_lap(psi_now, dx, dz) - params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

    recon = reconstruct_from_psi(psi_prev, psi_now, psi_next, dx, dz, dt, params.m)
    rho_num = recon["rho"]

    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    A = rho_exact[np.ix_(mask_x, mask_z)]
    A_num = rho_num[np.ix_(mask_x, mask_z)]
    summary = {
        "alpha": alpha,
        "relative_l1_A_solver": float(np.mean(np.abs(A_num - A)) / np.max(A)),
        "X_support_min": float(np.min(recon["X"][rho_num > 1e-3 * np.max(rho_num)])),
        "X_support_max": float(np.max(recon["X"][rho_num > 1e-3 * np.max(rho_num)])),
        "Q_support_min": float(np.min(recon["Q"][rho_num > 1e-3 * np.max(rho_num)])),
        "Q_support_max": float(np.max(recon["Q"][rho_num > 1e-3 * np.max(rho_num)])),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    np.savez_compressed(out / f"gravity_off_target_alpha_{alpha:.1f}.npz", x=x, z=z, rho=rho_num, rho_tilde=recon["rho_tilde"], X=recon["X"], Q=recon["Q"], S=recon["S"])
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    run(Path(__file__).resolve().parent.parent / "visualizations" / "bc_gravity_off_target", alpha=0.5)
