from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map
from simulate_bc_psfd_imex import spectral_k2, spectral_lap


def compute_q_and_x(
    alpha: float,
    t: float,
    dt: float,
    params: Exact2p1Params,
    rho_floor: float,
):
    x, z, psi_prev, _ = integrate_xz(alpha, t - dt, params)
    _, _, psi_now, _ = integrate_xz(alpha, t, params)
    _, _, psi_next, _ = integrate_xz(alpha, t + dt, params)

    rho_prev = np.abs(psi_prev) ** 2
    rho_now = np.abs(psi_now) ** 2
    rho_next = np.abs(psi_next) ** 2

    sqrt_prev = np.sqrt(np.maximum(rho_prev, rho_floor))
    sqrt_now = np.sqrt(np.maximum(rho_now, rho_floor))
    sqrt_next = np.sqrt(np.maximum(rho_next, rho_floor))

    dx = 2.0 * params.x_half_range / (params.nx - 1)
    dz = 2.0 * params.z_half_range / (params.nz - 1)
    k2 = spectral_k2(params.nx, params.nz, dx, dz)

    sqrt_tt = (sqrt_next - 2.0 * sqrt_now + sqrt_prev) / (dt * dt)
    sqrt_lap = spectral_lap(sqrt_now, k2)
    q_field = (sqrt_tt - sqrt_lap) / np.maximum(sqrt_now, rho_floor)
    x_field = params.m * params.m + q_field
    return x, z, rho_now, q_field, x_field


def support_stats(field: np.ndarray, rho: np.ndarray, support_frac: float) -> dict[str, float]:
    mask = rho > support_frac * float(np.max(rho))
    vals = field[mask]
    neg_frac = float(np.mean(vals < 0.0)) if vals.size else 0.0
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)) if vals.size else 0.0,
        "max": float(np.max(vals)) if vals.size else 0.0,
        "mean": float(np.mean(vals)) if vals.size else 0.0,
        "negative_fraction": neg_frac,
    }


def run(
    outdir: str | Path,
    alpha: float = 0.5,
    kappa: float = 1.0,
    dt_target: float = 0.025,
    x_half_range: float = 15.0,
    z_half_range: float = 15.0,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 61,
    nkz: int = 61,
    support_frac: float = 1.0e-4,
    rho_floor: float = 1.0e-14,
):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = Exact2p1Params(
        alpha=alpha,
        x_half_range=x_half_range,
        z_half_range=z_half_range,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )
    total_time = params.overlap_time
    steps = int(round(total_time / dt_target))
    dt = total_time / steps
    times = [0.0, 0.25 * total_time, 0.5 * total_time, 0.75 * total_time, total_time]

    summary: dict[str, object] = {
        "alpha": alpha,
        "kappa": kappa,
        "dt": dt,
        "support_fraction": support_frac,
        "times_over_T": [],
        "snapshots": [],
    }

    for idx, t in enumerate(times):
        x, z, rho, q_field, x_field = compute_q_and_x(alpha, t, dt, params, rho_floor)
        ratio_signed = x_field / kappa
        n_signed = ratio_signed * rho
        n_abs = np.abs(ratio_signed) * rho

        snap = {
            "time_over_T": float(t / total_time),
            "rho_support_stats": support_stats(rho, rho, support_frac),
            "Q_support_stats": support_stats(q_field, rho, support_frac),
            "X_support_stats": support_stats(x_field, rho, support_frac),
            "X_over_kappa_support_stats": support_stats(ratio_signed, rho, support_frac),
            "N_signed_support_stats": support_stats(n_signed, rho, support_frac),
            "N_abs_support_stats": support_stats(n_abs, rho, support_frac),
            "global": {
                "Q_min": float(np.min(q_field)),
                "Q_max": float(np.max(q_field)),
                "X_min": float(np.min(x_field)),
                "X_max": float(np.max(x_field)),
            },
        }
        summary["times_over_T"].append(float(t / total_time))
        summary["snapshots"].append(snap)

        if idx in (0, 2, 4):
            tag = f"t{idx}"
            render_map(out / f"rho_{tag}.png", x, z, rho, f"rho, t/T={t / total_time:.2f}")
            render_map(out / f"Q_{tag}.png", x, z, q_field, f"Q, t/T={t / total_time:.2f}")
            render_map(out / f"X_{tag}.png", x, z, x_field, f"X, t/T={t / total_time:.2f}")
            render_map(out / f"N_signed_{tag}.png", x, z, n_signed, f"N_signed, t/T={t / total_time:.2f}")
            render_map(out / f"N_abs_{tag}.png", x, z, n_abs, f"N_abs, t/T={t / total_time:.2f}")

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--kappa", type=float, default=1.0)
    parser.add_argument("--dt-target", type=float, default=0.025)
    parser.add_argument("--x-half-range", type=float, default=15.0)
    parser.add_argument("--z-half-range", type=float, default=15.0)
    parser.add_argument("--nx", type=int, default=181)
    parser.add_argument("--nz", type=int, default=181)
    parser.add_argument("--nkx", type=int, default=61)
    parser.add_argument("--nkz", type=int, default=61)
    parser.add_argument("--support-frac", type=float, default=1.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-14)
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(
            Path(__file__).resolve().parent.parent
            / "visualizations"
            / "mixed_transform_shared_backbone_audit"
        ),
    )
    args = parser.parse_args()
    run(
        args.outdir,
        alpha=args.alpha,
        kappa=args.kappa,
        dt_target=args.dt_target,
        x_half_range=args.x_half_range,
        z_half_range=args.z_half_range,
        nx=args.nx,
        nz=args.nz,
        nkx=args.nkx,
        nkz=args.nkz,
        support_frac=args.support_frac,
        rho_floor=args.rho_floor,
    )
