from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz


def d_dx(f: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def d_dz(f: np.ndarray, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out


def d2_dx2(f: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dx * dx)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def d2_dz2(f: np.ndarray, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / (dz * dz)
    out[:, 0] = out[:, 1]
    out[:, -1] = out[:, -2]
    return out


def summarize(field: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.abs(field[mask])
    return {
        "linf": float(np.max(vals)),
        "l1_mean": float(np.mean(vals)),
        "l2_rms": float(np.sqrt(np.mean(vals * vals))),
    }


def main():
    out = Path(__file__).resolve().parent.parent / "kg_examples" / "outputs" / "gravity_off_equivalence_summary.json"
    params = Exact2p1Params(x_half_range=20.0, z_half_range=20.0, nx=121, nz=121, nkx=81, nkz=81)
    t0 = params.overlap_time
    dt = 0.01

    summary: dict[str, dict[str, object]] = {"params": {"t0": t0, "dt_fd": dt}}

    for alpha in (0.5, 1.0):
        x, z, psi0, rho = integrate_xz(alpha, t0, params)
        dx = float(x[1] - x[0])
        dz = float(z[1] - z[0])
        Xg, Zg = np.meshgrid(x, z, indexing="ij")

        psi_p = evaluate_points(alpha, Xg, Zg, t0 + dt, params)
        psi_m = evaluate_points(alpha, Xg, Zg, t0 - dt, params)

        dpsi_t = (psi_p - psi_m) / (2.0 * dt)
        dpsi_x = d_dx(psi0, dx)
        dpsi_z = d_dz(psi0, dz)

        rho_safe = np.maximum(rho, 1e-12)
        sqrt_rho = np.sqrt(rho_safe)
        sqrt_rho_p = np.sqrt(np.maximum(np.abs(psi_p) ** 2, 1e-12))
        sqrt_rho_m = np.sqrt(np.maximum(np.abs(psi_m) ** 2, 1e-12))

        s_t = np.imag(np.conj(psi0) * dpsi_t) / rho_safe
        s_x = np.imag(np.conj(psi0) * dpsi_x) / rho_safe
        s_z = np.imag(np.conj(psi0) * dpsi_z) / rho_safe

        X = s_t * s_t - s_x * s_x - s_z * s_z

        d2_sqrt_dt2 = (sqrt_rho_p - 2.0 * sqrt_rho + sqrt_rho_m) / (dt * dt)
        d2_sqrt_dx2 = d2_dx2(sqrt_rho, dx)
        d2_sqrt_dz2 = d2_dz2(sqrt_rho, dz)
        Q = (d2_sqrt_dt2 - d2_sqrt_dx2 - d2_sqrt_dz2) / np.maximum(sqrt_rho, 1e-12)

        hj_residual_A = X - Q - params.m * params.m

        # Gravity-off constrained reconstruction from the main rank-1 solution.
        D = -Q / np.maximum(X * X, 1e-12)
        Xt = np.maximum(X, 1e-12)
        sqrt_det_ratio = np.sqrt(np.maximum(Xt / (params.m * params.m), 1e-12))
        rho_tilde = rho * sqrt_det_ratio

        # Current identity implied by the constraint:
        # sqrt(-g~) rho~ g~^{mu nu} \partial_nu S = rho \partial^mu S
        current_A_0 = rho * s_t
        current_A_x = -rho * s_x
        current_A_z = -rho * s_z

        weight = sqrt_det_ratio * rho_tilde * (params.m * params.m / Xt)
        current_tilde_0 = weight * s_t
        current_tilde_x = weight * (-s_x)
        current_tilde_z = weight * (-s_z)

        current_identity_0 = current_tilde_0 - current_A_0
        current_identity_x = current_tilde_x - current_A_x
        current_identity_z = current_tilde_z - current_A_z

        support_mask = rho > (1e-3 * np.max(rho))

        summary[str(alpha)] = {
            "support_fraction": float(np.mean(support_mask)),
            "A_HJ_residual": summarize(hj_residual_A, support_mask),
            "constrained_tilde_HJ_residual": summarize((X - Q - params.m * params.m), support_mask),
            "current_identity_residual_t": summarize(current_identity_0, support_mask),
            "current_identity_residual_x": summarize(current_identity_x, support_mask),
            "current_identity_residual_z": summarize(current_identity_z, support_mask),
            "X_stats": {
                "min_support": float(np.min(X[support_mask])),
                "max_support": float(np.max(X[support_mask])),
            },
            "Q_stats": {
                "min_support": float(np.min(Q[support_mask])),
                "max_support": float(np.max(Q[support_mask])),
            },
        }

    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
