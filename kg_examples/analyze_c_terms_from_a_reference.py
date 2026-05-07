from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from mixed_tilde_initial_data import (
    exact_localized_wave_derivatives,
    localized_direct_tilde_coordinate_snapshot,
)
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def d1x(f: np.ndarray, dx: float) -> np.ndarray:
    return np.gradient(f, dx, axis=0, edge_order=2)


def d1z(f: np.ndarray, dz: float) -> np.ndarray:
    return np.gradient(f, dz, axis=1, edge_order=2)


def d2xx(f: np.ndarray, dx: float) -> np.ndarray:
    return np.gradient(np.gradient(f, dx, axis=0, edge_order=2), dx, axis=0, edge_order=2)


def d2zz(f: np.ndarray, dz: float) -> np.ndarray:
    return np.gradient(np.gradient(f, dz, axis=1, edge_order=2), dz, axis=1, edge_order=2)


def d2xz(f: np.ndarray, dx: float, dz: float) -> np.ndarray:
    return np.gradient(np.gradient(f, dx, axis=0, edge_order=2), dz, axis=1, edge_order=2)


def metric_jets_full(
    metric_m: np.ndarray,
    metric_0: np.ndarray,
    metric_p: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    nx, nz = metric_0.shape[:2]
    dg = np.zeros((nx, nz, 3, 3, 3), dtype=float)
    d2g = np.zeros((nx, nz, 3, 3, 3, 3), dtype=float)

    metric_t = (metric_p - metric_m) / (2.0 * dt)
    metric_tt = (metric_p - 2.0 * metric_0 + metric_m) / (dt * dt)

    dg[..., 0, :, :] = metric_t
    dg[..., 1, :, :] = d1x(metric_0, dx)
    dg[..., 2, :, :] = d1z(metric_0, dz)

    d2g[..., 0, 0, :, :] = metric_tt
    d2g[..., 0, 1, :, :] = d1x(metric_t, dx)
    d2g[..., 1, 0, :, :] = d2g[..., 0, 1, :, :]
    d2g[..., 0, 2, :, :] = d1z(metric_t, dz)
    d2g[..., 2, 0, :, :] = d2g[..., 0, 2, :, :]
    d2g[..., 1, 1, :, :] = d2xx(metric_0, dx)
    d2g[..., 1, 2, :, :] = d2xz(metric_0, dx, dz)
    d2g[..., 2, 1, :, :] = d2g[..., 1, 2, :, :]
    d2g[..., 2, 2, :, :] = d2zz(metric_0, dz)
    return dg, d2g


def scalar_jets_full(
    f_m: np.ndarray,
    f_0: np.ndarray,
    f_p: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    df = np.zeros(f_0.shape + (3,), dtype=float)
    d2f = np.zeros(f_0.shape + (3, 3), dtype=float)

    f_t = (f_p - f_m) / (2.0 * dt)
    f_tt = (f_p - 2.0 * f_0 + f_m) / (dt * dt)

    df[..., 0] = f_t
    df[..., 1] = d1x(f_0, dx)
    df[..., 2] = d1z(f_0, dz)

    d2f[..., 0, 0] = f_tt
    d2f[..., 0, 1] = d1x(f_t, dx)
    d2f[..., 1, 0] = d2f[..., 0, 1]
    d2f[..., 0, 2] = d1z(f_t, dz)
    d2f[..., 2, 0] = d2f[..., 0, 2]
    d2f[..., 1, 1] = d2xx(f_0, dx)
    d2f[..., 1, 2] = d2xz(f_0, dx, dz)
    d2f[..., 2, 1] = d2f[..., 1, 2]
    d2f[..., 2, 2] = d2zz(f_0, dz)
    return df, d2f


def raw_ricci_tensor_scalar(metric_cov: np.ndarray, dg: np.ndarray, d2g: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    nx, nz = metric_cov.shape[:2]
    ricci = np.zeros((nx, nz, 3, 3), dtype=float)
    r_scalar = np.zeros((nx, nz), dtype=float)
    for i in range(nx):
        for j in range(nz):
            g = metric_cov[i, j]
            ginv = np.linalg.pinv(g, rcond=1.0e-12, hermitian=True)

            gamma1 = np.zeros((3, 3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        gamma1[a, b, c] = 0.5 * (
                            dg[i, j, b, a, c] + dg[i, j, c, a, b] - dg[i, j, a, b, c]
                        )
            gamma2 = np.einsum("ad,dbc->abc", ginv, gamma1, optimize=True)

            dginv = np.zeros((3, 3, 3), dtype=float)
            for e in range(3):
                dginv[e] = -ginv @ dg[i, j, e] @ ginv

            dgamma1 = np.zeros((3, 3, 3, 3), dtype=float)
            for e in range(3):
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            dgamma1[e, a, b, c] = 0.5 * (
                                d2g[i, j, e, b, a, c] + d2g[i, j, e, c, a, b] - d2g[i, j, e, a, b, c]
                            )

            dgamma2 = np.zeros((3, 3, 3, 3), dtype=float)
            for e in range(3):
                dgamma2[e] = np.einsum("ad,dbc->abc", dginv[e], gamma1, optimize=True) + np.einsum(
                    "ad,dbc->abc", ginv, dgamma1[e], optimize=True
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
                            quad1 += gamma2[c, a, b] * gamma2[d, c, d]
                            quad2 += gamma2[c, a, d] * gamma2[d, b, c]
                    ric[a, b] = term1 - term2 + quad1 - quad2

            ricci[i, j] = ric
            r_scalar[i, j] = float(np.einsum("ab,ab->", ginv, ric))
    return ricci, r_scalar


def box_scalar(metric_cov: np.ndarray, dg: np.ndarray, scalar_grad: np.ndarray, scalar_hess: np.ndarray) -> np.ndarray:
    nx, nz = metric_cov.shape[:2]
    out = np.zeros((nx, nz), dtype=float)
    for i in range(nx):
        for j in range(nz):
            g = metric_cov[i, j]
            ginv = np.linalg.pinv(g, rcond=1.0e-12, hermitian=True)
            gamma1 = np.zeros((3, 3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        gamma1[a, b, c] = 0.5 * (
                            dg[i, j, b, a, c] + dg[i, j, c, a, b] - dg[i, j, a, b, c]
                        )
            gamma2 = np.einsum("ad,dbc->abc", ginv, gamma1, optimize=True)
            val = 0.0
            for a in range(3):
                for b in range(3):
                    nabla_ab = scalar_hess[i, j, a, b]
                    for c in range(3):
                        nabla_ab -= gamma2[c, a, b] * scalar_grad[i, j, c]
                    val += ginv[a, b] * nabla_ab
            out[i, j] = val
    return out


def support_stats(field: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.abs(field[mask])
    if vals.size == 0:
        return {"abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0}
    return {
        "abs_max": float(np.max(vals)),
        "abs_mean": float(np.mean(vals)),
        "abs_median": float(np.median(vals)),
        "p95": float(np.percentile(vals, 95.0)),
    }


def analyze_snapshot(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    probe_dt: float,
    ell: float,
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
    _, r_tilde = raw_ricci_tensor_scalar(metric_0, dg, d2g)

    l2 = ell * ell
    f0 = r_tilde / np.sqrt(1.0 + (l2 * r_tilde) ** 2)
    fr0 = 1.0 / np.power(1.0 + (l2 * r_tilde) ** 2, 1.5)

    # Recompute f_R at neighboring times for time derivatives.
    dg_m, d2g_m = metric_jets_full(
        np.real_if_close(localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - 2.0 * probe_dt, rho_floor=rho_floor)["cov_txz"]),
        metric_m,
        metric_0,
        probe_dt,
        dx,
        dz,
    )
    _, r_m = raw_ricci_tensor_scalar(metric_m, dg_m, d2g_m)
    dg_p, d2g_p = metric_jets_full(
        metric_0,
        metric_p,
        np.real_if_close(localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + 2.0 * probe_dt, rho_floor=rho_floor)["cov_txz"]),
        probe_dt,
        dx,
        dz,
    )
    _, r_p = raw_ricci_tensor_scalar(metric_p, dg_p, d2g_p)

    fr_m = 1.0 / np.power(1.0 + (l2 * r_m) ** 2, 1.5)
    fr_p = 1.0 / np.power(1.0 + (l2 * r_p) ** 2, 1.5)

    dfr, d2fr = scalar_jets_full(fr_m, fr0, fr_p, probe_dt, dx, dz)
    box_fr = box_scalar(metric_0, dg, dfr, d2fr)

    sqrt_abs_g = np.sqrt(np.abs(np.linalg.det(metric_0.reshape(-1, 3, 3))).reshape(metric_0.shape[:2]))
    rho = snap_0["bohm"]["rho"]
    support = rho > 1.0e-3 * float(np.max(rho))
    core = rho > 1.0e-1 * float(np.max(rho))

    term_frR = fr0 * r_tilde
    term_alg = -2.0 * f0
    term_box = 3.0 * box_fr

    ratio_support = float(support_stats(term_box, support)["abs_max"] / max(support_stats(term_alg, support)["abs_max"], 1.0e-30))
    ratio_core = float(support_stats(term_box, core)["abs_max"] / max(support_stats(term_alg, core)["abs_max"], 1.0e-30))

    return {
        "time": float(t),
        "support_fraction": float(np.mean(support)),
        "core_fraction": float(np.mean(core)),
        "R_tilde_support": support_stats(r_tilde, support),
        "R_tilde_core": support_stats(r_tilde, core),
        "sqrtg_R_tilde_support": support_stats(sqrt_abs_g * r_tilde, support),
        "sqrtg_R_tilde_core": support_stats(sqrt_abs_g * r_tilde, core),
        "f_support": support_stats(f0, support),
        "f_core": support_stats(f0, core),
        "fR_support": support_stats(fr0, support),
        "fR_core": support_stats(fr0, core),
        "term_fR_R_support": support_stats(term_frR, support),
        "term_fR_R_core": support_stats(term_frR, core),
        "term_minus2f_support": support_stats(term_alg, support),
        "term_minus2f_core": support_stats(term_alg, core),
        "term_3boxfR_support": support_stats(term_box, support),
        "term_3boxfR_core": support_stats(term_box, core),
        "boxfR_support": support_stats(box_fr, support),
        "boxfR_core": support_stats(box_fr, core),
        "derivative_to_algebraic_ratio_support": ratio_support,
        "derivative_to_algebraic_ratio_core": ratio_core,
        "ell2_R_support_max": float((l2) * support_stats(r_tilde, support)["abs_max"]),
        "ell2_R_core_max": float((l2) * support_stats(r_tilde, core)["abs_max"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=100.0)
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
    snapshots = [
        analyze_snapshot(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t=float(t),
            probe_dt=args.probe_dt,
            ell=args.ell,
            rho_floor=args.rho_floor,
        )
        for t in times
    ]

    peak_support_ratio = max(s["derivative_to_algebraic_ratio_support"] for s in snapshots)
    peak_core_ratio = max(s["derivative_to_algebraic_ratio_core"] for s in snapshots)
    peak_support_absR = max(s["R_tilde_support"]["abs_max"] for s in snapshots)
    peak_core_absR = max(s["R_tilde_core"]["abs_max"] for s in snapshots)

    summary = {
        "params": {
            "ell": args.ell,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "nx": args.nx,
            "nz": args.nz,
            "times": [float(t) for t in times],
        },
        "summary": {
            "peak_support_abs_R_tilde": peak_support_absR,
            "peak_core_abs_R_tilde": peak_core_absR,
            "peak_support_derivative_to_algebraic_ratio": peak_support_ratio,
            "peak_core_derivative_to_algebraic_ratio": peak_core_ratio,
            "derivative_term_dominates_support": bool(peak_support_ratio > 1.0),
            "derivative_term_dominates_core": bool(peak_core_ratio > 1.0),
        },
        "snapshots": snapshots,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
