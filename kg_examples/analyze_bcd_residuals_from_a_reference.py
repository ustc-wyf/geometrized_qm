from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_c_terms_from_a_reference import (
    metric_jets_full,
    scalar_jets_full,
)
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def geometry_data(
    metric_cov: np.ndarray,
    dg: np.ndarray,
    d2g: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nx, nz = metric_cov.shape[:2]
    ginv = np.zeros_like(metric_cov)
    gamma2 = np.zeros(metric_cov.shape[:2] + (3, 3, 3), dtype=float)
    ricci = np.zeros(metric_cov.shape[:2] + (3, 3), dtype=float)
    r_scalar = np.zeros(metric_cov.shape[:2], dtype=float)

    for i in range(nx):
        for j in range(nz):
            g = metric_cov[i, j]
            ginv_ij = np.linalg.pinv(g, rcond=1.0e-12, hermitian=True)
            ginv[i, j] = ginv_ij

            gamma1 = np.zeros((3, 3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        gamma1[a, b, c] = 0.5 * (
                            dg[i, j, b, a, c] + dg[i, j, c, a, b] - dg[i, j, a, b, c]
                        )
            gamma2_ij = np.einsum("ad,dbc->abc", ginv_ij, gamma1, optimize=True)
            gamma2[i, j] = gamma2_ij

            dginv = np.zeros((3, 3, 3), dtype=float)
            for e in range(3):
                dginv[e] = -ginv_ij @ dg[i, j, e] @ ginv_ij

            dgamma1 = np.zeros((3, 3, 3, 3), dtype=float)
            for e in range(3):
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            dgamma1[e, a, b, c] = 0.5 * (
                                d2g[i, j, e, b, a, c]
                                + d2g[i, j, e, c, a, b]
                                - d2g[i, j, e, a, b, c]
                            )

            dgamma2 = np.zeros((3, 3, 3, 3), dtype=float)
            for e in range(3):
                dgamma2[e] = np.einsum("ad,dbc->abc", dginv[e], gamma1, optimize=True) + np.einsum(
                    "ad,dbc->abc", ginv_ij, dgamma1[e], optimize=True
                )

            ric = np.zeros((3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    term1 = sum(dgamma2[c, c, a, b] for c in range(3))
                    term2 = sum(dgamma2[b, c, a, c] for c in range(3))
                    quad1 = 0.0
                    quad2 = 0.0
                    for c in range(3):
                        for d in range(3):
                            quad1 += gamma2_ij[c, a, b] * gamma2_ij[d, c, d]
                            quad2 += gamma2_ij[c, a, d] * gamma2_ij[d, b, c]
                    ric[a, b] = term1 - term2 + quad1 - quad2
            ricci[i, j] = ric
            r_scalar[i, j] = float(np.einsum("ab,ab->", ginv_ij, ric))
    return ginv, gamma2, ricci, r_scalar


def covariant_hessian_from_jets(
    scalar_grad: np.ndarray,
    scalar_hess: np.ndarray,
    gamma2: np.ndarray,
) -> np.ndarray:
    nx, nz = scalar_grad.shape[:2]
    cov_hess = np.zeros((nx, nz, 3, 3), dtype=float)
    for i in range(nx):
        for j in range(nz):
            for a in range(3):
                for b in range(3):
                    val = scalar_hess[i, j, a, b]
                    for c in range(3):
                        val -= gamma2[i, j, c, a, b] * scalar_grad[i, j, c]
                    cov_hess[i, j, a, b] = val
    return cov_hess


def tensor_frobenius(tensor: np.ndarray) -> np.ndarray:
    return np.sqrt(np.sum(np.abs(tensor) ** 2, axis=(-2, -1)))


def scalar_stats(field: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.abs(field[mask])
    if vals.size == 0:
        return {"abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0}
    return {
        "abs_max": float(np.max(vals)),
        "abs_mean": float(np.mean(vals)),
        "abs_median": float(np.median(vals)),
        "p95": float(np.percentile(vals, 95.0)),
    }


def stress_tensor_tilde(
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    s_t: np.ndarray,
    s_x: np.ndarray,
    s_z: np.ndarray,
    m: float,
) -> tuple[np.ndarray, np.ndarray]:
    u_cov = np.stack([s_t, s_x, s_z], axis=-1)
    tilde_x = np.einsum("...ab,...a,...b->...", metric_inv, u_cov, u_cov, optimize=True)
    outer = np.einsum("...a,...b->...ab", u_cov, u_cov, optimize=True)
    t_tensor = 2.0 * rho[..., None, None] * outer - metric_cov * (
        rho * (tilde_x - m * m)
    )[..., None, None]
    return np.real_if_close(t_tensor), np.real_if_close(tilde_x)


def branch_f_and_fr(branch: str, r_scalar: np.ndarray, ell: float) -> tuple[np.ndarray, np.ndarray]:
    if branch == "B":
        return np.real_if_close(r_scalar), np.ones_like(r_scalar)
    l2 = ell * ell
    if branch == "C":
        f = r_scalar / np.sqrt(1.0 + (l2 * r_scalar) ** 2)
        fr = 1.0 / np.power(1.0 + (l2 * r_scalar) ** 2, 1.5)
        return np.real_if_close(f), np.real_if_close(fr)
    if branch == "D":
        tanh_arg = np.tanh(l2 * r_scalar)
        f = tanh_arg / l2
        fr = 1.0 - tanh_arg * tanh_arg
        return np.real_if_close(f), np.real_if_close(fr)
    raise ValueError(f"Unsupported branch: {branch}")


def branch_residual_snapshot(
    branch: str,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    probe_dt: float,
    ell: float,
    mp: float,
    mass: float,
    rho_floor: float,
) -> dict[str, object]:
    snap_m = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - probe_dt, rho_floor=rho_floor)
    snap_0 = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t, rho_floor=rho_floor)
    snap_p = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + probe_dt, rho_floor=rho_floor)

    metric_m = np.real_if_close(snap_m["cov_txz"])
    metric_0 = np.real_if_close(snap_0["cov_txz"])
    metric_p = np.real_if_close(snap_p["cov_txz"])

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, probe_dt, dx, dz)
    ginv_0, gamma2_0, ricci_0, r_tilde_0 = geometry_data(metric_0, dg, d2g)

    rho = np.real_if_close(snap_0["bohm"]["rho"])
    s_t = np.real_if_close(snap_0["bohm"]["s_t"])
    s_x = np.real_if_close(snap_0["bohm"]["s_x"])
    s_z = np.real_if_close(snap_0["bohm"]["s_z"])

    t_tensor, tilde_x = stress_tensor_tilde(metric_0, ginv_0, rho, s_t, s_x, s_z, m=mass)
    rhs_tensor = t_tensor / (mp * mp)

    f0, fr0 = branch_f_and_fr(branch, r_tilde_0, ell)
    if branch == "B":
        cov_hess_fr = np.zeros(metric_0.shape[:2] + (3, 3), dtype=float)
        box_fr = np.zeros(metric_0.shape[:2], dtype=float)
    else:
        dg_m, d2g_m = metric_jets_full(
            np.real_if_close(
                localized_direct_tilde_coordinate_snapshot(
                    psi0, psi0_hat, omega, x, z, t - 2.0 * probe_dt, rho_floor=rho_floor
                )["cov_txz"]
            ),
            metric_m,
            metric_0,
            probe_dt,
            dx,
            dz,
        )
        _, _, _, r_tilde_m = geometry_data(metric_m, dg_m, d2g_m)
        dg_p, d2g_p = metric_jets_full(
            metric_0,
            metric_p,
            np.real_if_close(
                localized_direct_tilde_coordinate_snapshot(
                    psi0, psi0_hat, omega, x, z, t + 2.0 * probe_dt, rho_floor=rho_floor
                )["cov_txz"]
            ),
            probe_dt,
            dx,
            dz,
        )
        _, _, _, r_tilde_p = geometry_data(metric_p, dg_p, d2g_p)
        _, fr_m = branch_f_and_fr(branch, r_tilde_m, ell)
        _, fr_p = branch_f_and_fr(branch, r_tilde_p, ell)
        dfr, d2fr = scalar_jets_full(fr_m, fr0, fr_p, probe_dt, dx, dz)
        cov_hess_fr = covariant_hessian_from_jets(dfr, d2fr, gamma2_0)
        box_fr = np.einsum("...ab,...ab->...", ginv_0, cov_hess_fr, optimize=True)

    algebraic_tensor = fr0[..., None, None] * ricci_0 - 0.5 * f0[..., None, None] * metric_0
    derivative_tensor = -(cov_hess_fr - metric_0 * box_fr[..., None, None])
    lhs_tensor = algebraic_tensor + derivative_tensor
    residual_tensor = lhs_tensor - rhs_tensor

    support_mask = rho > 1.0e-3 * float(np.max(rho))
    dense_core_mask = rho > 1.0e-1 * float(np.max(rho))

    lhs_norm = tensor_frobenius(lhs_tensor)
    rhs_norm = tensor_frobenius(rhs_tensor)
    residual_norm = tensor_frobenius(residual_tensor)
    algebraic_norm = tensor_frobenius(algebraic_tensor)
    derivative_norm = tensor_frobenius(derivative_tensor)

    relative_to_rhs = residual_norm / np.maximum(rhs_norm, 1.0e-30)
    relative_to_lhs = residual_norm / np.maximum(lhs_norm, 1.0e-30)
    mass_shell_defect = np.abs(tilde_x - mass * mass)

    return {
        "branch": branch,
        "time": float(t),
        "R_tilde_support": scalar_stats(r_tilde_0, support_mask),
        "R_tilde_dense_core": scalar_stats(r_tilde_0, dense_core_mask),
        "tilde_X_minus_m2_support": scalar_stats(tilde_x - mass * mass, support_mask),
        "tilde_X_minus_m2_dense_core": scalar_stats(tilde_x - mass * mass, dense_core_mask),
        "lhs_norm_support": scalar_stats(lhs_norm, support_mask),
        "lhs_norm_dense_core": scalar_stats(lhs_norm, dense_core_mask),
        "rhs_norm_support": scalar_stats(rhs_norm, support_mask),
        "rhs_norm_dense_core": scalar_stats(rhs_norm, dense_core_mask),
        "residual_norm_support": scalar_stats(residual_norm, support_mask),
        "residual_norm_dense_core": scalar_stats(residual_norm, dense_core_mask),
        "algebraic_norm_support": scalar_stats(algebraic_norm, support_mask),
        "algebraic_norm_dense_core": scalar_stats(algebraic_norm, dense_core_mask),
        "derivative_norm_support": scalar_stats(derivative_norm, support_mask),
        "derivative_norm_dense_core": scalar_stats(derivative_norm, dense_core_mask),
        "relative_residual_to_rhs_support": scalar_stats(relative_to_rhs, support_mask),
        "relative_residual_to_rhs_dense_core": scalar_stats(relative_to_rhs, dense_core_mask),
        "relative_residual_to_lhs_support": scalar_stats(relative_to_lhs, support_mask),
        "relative_residual_to_lhs_dense_core": scalar_stats(relative_to_lhs, dense_core_mask),
        "mass_shell_defect_support": scalar_stats(mass_shell_defect, support_mask),
        "mass_shell_defect_dense_core": scalar_stats(mass_shell_defect, dense_core_mask),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--nx", type=int, default=96)
    parser.add_argument("--nz", type=int, default=96)
    parser.add_argument("--times", type=str, default="0,8,16")
    args = parser.parse_args()

    params = FlatLocalizedCrossingParams(nx=args.nx, nz=args.nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    times = [float(tok) for tok in args.times.split(",") if tok.strip()]
    snapshots: list[dict[str, object]] = []
    for branch in ("B", "C", "D"):
        for t in times:
            snapshots.append(
                branch_residual_snapshot(
                    branch=branch,
                    psi0=psi0,
                    psi0_hat=psi0_hat,
                    omega=omega,
                    x=x,
                    z=z,
                    t=t,
                    probe_dt=args.probe_dt,
                    ell=args.ell,
                    mp=args.mp,
                    mass=params.m,
                    rho_floor=args.rho_floor,
                )
            )

    by_branch: dict[str, dict[str, float]] = {}
    for branch in ("B", "C", "D"):
        entries = [s for s in snapshots if s["branch"] == branch]
        by_branch[branch] = {
            "peak_support_residual_abs_max": max(s["residual_norm_support"]["abs_max"] for s in entries),
            "peak_dense_core_residual_abs_max": max(s["residual_norm_dense_core"]["abs_max"] for s in entries),
            "peak_support_relative_to_rhs": max(s["relative_residual_to_rhs_support"]["abs_max"] for s in entries),
            "peak_dense_core_relative_to_rhs": max(s["relative_residual_to_rhs_dense_core"]["abs_max"] for s in entries),
            "peak_support_derivative_abs_max": max(s["derivative_norm_support"]["abs_max"] for s in entries),
            "peak_dense_core_derivative_abs_max": max(s["derivative_norm_dense_core"]["abs_max"] for s in entries),
        }

    summary = {
        "params": {
            "ell": args.ell,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "nx": args.nx,
            "nz": args.nz,
            "times": times,
            "mass": params.m,
        },
        "definitions": {
            "support_mask": "rho > 1e-3 * rho_max",
            "dense_core_mask": "rho > 1e-1 * rho_max",
            "branch_B": "f_B(R_tilde) = R_tilde",
            "branch_C": "f_C(R_tilde) = R_tilde / sqrt(1 + (ell^2 R_tilde)^2)",
            "branch_D": "f_D(R_tilde) = tanh(ell^2 R_tilde) / ell^2",
            "residual_tensor": "E_munu = [f_R R_munu - 1/2 f g_munu - (nabla_munu - g_munu Box) f_R] - T_munu / M_P^2",
        },
        "summary": by_branch,
        "snapshots": snapshots,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
