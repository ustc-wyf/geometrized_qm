from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_d_trace_interface_weak_form import find_level_segments, parse_positive_levels, real_array
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from prototype_d_reference_imex import reference_snapshot_data
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def d_phi_d_raw_y(raw_y: np.ndarray) -> np.ndarray:
    tanh_y = np.tanh(raw_y)
    phi = 1.0 - tanh_y * tanh_y
    return -2.0 * phi * tanh_y


def d_branch_f(raw_y: np.ndarray, ell: float) -> np.ndarray:
    return np.tanh(raw_y) / (ell * ell)


def d_branch_phi(raw_y: np.ndarray) -> np.ndarray:
    tanh_y = np.tanh(raw_y)
    return 1.0 - tanh_y * tanh_y


def sign_fractions(values: list[float], scale_floor: float = 1.0e-12) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"positive": 0.0, "negative": 0.0, "near_zero": 0.0}
    scale = max(float(np.percentile(np.abs(vals), 95.0)), scale_floor)
    eps = 1.0e-8 * scale
    return {
        "positive": float(np.mean(vals > eps)),
        "negative": float(np.mean(vals < -eps)),
        "near_zero": float(np.mean(np.abs(vals) <= eps)),
    }


def interface_jump_rows(
    raw_y: np.ndarray,
    raw_y_m: np.ndarray,
    raw_y_p: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    stress_trace: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    ell: float,
    probe_dt: float,
    half_width_y: float,
    samples: int,
) -> list[dict[str, float]]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_t = (raw_y_p - raw_y_m) / (2.0 * probe_dt)
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)

    rows: list[dict[str, float]] = []
    q = np.linspace(level - half_width_y, level + half_width_y, samples)
    phi_prime = d_phi_d_raw_y(q)
    delta_phi_prime = float(phi_prime[-1] - phi_prime[0])
    phi_q = d_branch_phi(q)
    f_q = d_branch_f(q, ell)
    chi_q = q / (ell * ell)
    # Stress is slowly varying across the thin layer in this leading estimate.
    for seg in find_level_segments(raw_y, rho, x, z, level):
        xp = np.array([seg["x"]])
        zp = np.array([seg["z"]])
        ft = float(bilinear_sample(raw_y_t, x, z, xp, zp)[0])
        fx = float(bilinear_sample(raw_y_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(raw_y_z, x, z, xp, zp)[0])
        ginv_here = np.zeros((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                ginv_here[a, b] = float(bilinear_sample(metric_inv[..., a, b], x, z, xp, zp)[0])
        covector = np.array([ft, fx, fz], dtype=float)
        normal_norm = float(covector @ ginv_here @ covector)
        abs_norm_sqrt = float(np.sqrt(abs(normal_norm)))
        if abs_norm_sqrt <= 1.0e-300:
            continue
        sigma = float(np.sign(normal_norm))
        stress_here = float(bilinear_sample(stress_trace, x, z, xp, zp)[0])

        # Thin-layer spacetime normal estimate:
        #   2 Box(phi) -> 2 sigma [d_phi/d_eta]
        # with d_eta = d(raw_y)/sqrt(|grad raw_y|_g^2).
        normal_flux_jump = float(2.0 * sigma * abs_norm_sqrt * delta_phi_prime)
        algebraic_q = phi_q * chi_q - 1.5 * f_q - stress_here
        algebraic_integral_raw_y = float(np.trapezoid(algebraic_q, q))
        algebraic_abs_integral_raw_y = float(np.trapezoid(np.abs(algebraic_q), q))
        bulk_integral = algebraic_integral_raw_y / abs_norm_sqrt
        bulk_abs_integral = algebraic_abs_integral_raw_y / abs_norm_sqrt
        weak_residual = normal_flux_jump + bulk_integral
        jump_law_residual_scaled = float(2.0 * sigma * abs_norm_sqrt * abs_norm_sqrt * delta_phi_prime + algebraic_integral_raw_y)
        target_same_signature = np.nan
        denom = 2.0 * sigma * delta_phi_prime
        if abs(denom) > 1.0e-300:
            target_abs_norm_sq = -algebraic_integral_raw_y / denom
            if target_abs_norm_sq >= 0.0:
                target_same_signature = float(sigma * target_abs_norm_sq)
        rows.append(
            {
                "x": float(seg["x"]),
                "z": float(seg["z"]),
                "level": float(level),
                "normal_norm_tilde_inverse": normal_norm,
                "normal_norm_sign": sigma,
                "sqrt_abs_normal_norm": abs_norm_sqrt,
                "delta_phi_prime": delta_phi_prime,
                "normal_flux_jump": normal_flux_jump,
                "algebraic_integral_raw_y": algebraic_integral_raw_y,
                "algebraic_abs_integral_raw_y": algebraic_abs_integral_raw_y,
                "bulk_integral": bulk_integral,
                "bulk_abs_integral": bulk_abs_integral,
                "weak_residual": weak_residual,
                "jump_law_residual_scaled": jump_law_residual_scaled,
                "target_normal_norm_same_signature": target_same_signature,
                "weak_to_abs_ratio": float(abs(weak_residual) / max(abs(normal_flux_jump) + bulk_abs_integral, 1.0e-300)),
                "stress_trace": stress_here,
                "raw_y_t": ft,
                "raw_y_x": fx,
                "raw_y_z": fz,
            }
        )
    return rows


def summarize(rows: list[dict[str, float]]) -> dict[str, object]:
    return {
        "count": len(rows),
        "normal_flux_jump": scalar_stats([row["normal_flux_jump"] for row in rows]),
        "bulk_integral": scalar_stats([row["bulk_integral"] for row in rows]),
        "algebraic_integral_raw_y": scalar_stats([row["algebraic_integral_raw_y"] for row in rows]),
        "bulk_abs_integral": scalar_stats([row["bulk_abs_integral"] for row in rows]),
        "weak_residual": scalar_stats([row["weak_residual"] for row in rows]),
        "jump_law_residual_scaled": scalar_stats([row["jump_law_residual_scaled"] for row in rows]),
        "weak_to_abs_ratio": scalar_stats([row["weak_to_abs_ratio"] for row in rows]),
        "sqrt_abs_normal_norm": scalar_stats([row["sqrt_abs_normal_norm"] for row in rows]),
        "target_normal_norm_same_signature": scalar_stats([row["target_normal_norm_same_signature"] for row in rows]),
        "normal_norm_sign_fractions": sign_fractions([row["normal_norm_tilde_inverse"] for row in rows]),
        "global_weak_to_abs_ratio": float(
            abs(sum(row["weak_residual"] for row in rows))
            / max(sum(abs(row["normal_flux_jump"]) + row["bulk_abs_integral"] for row in rows), 1.0e-300)
        )
        if rows
        else 0.0,
    }


def render_plot(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    rows: list[dict[str, float]],
    levels: list[float],
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)
    im0 = axes[0].pcolormesh(xg, zg, np.log10(np.maximum(rho, 1.0e-16)), shading="auto", cmap="viridis")
    axes[0].contour(xg, zg, np.abs(raw_y), levels=levels, colors="white", linewidths=0.8)
    axes[0].set_title("log10 rho; white: |raw_y| levels")
    fig.colorbar(im0, ax=axes[0])

    if rows:
        xs = [row["x"] for row in rows]
        zs = [row["z"] for row in rows]
        flux = [row["normal_flux_jump"] for row in rows]
        vmax1 = max(float(np.percentile(np.abs(flux), 95.0)), 1.0e-12)
        sc1 = axes[1].scatter(xs, zs, c=flux, s=5.0, cmap="coolwarm", vmin=-vmax1, vmax=vmax1)
        fig.colorbar(sc1, ax=axes[1])
        weak = [row["weak_residual"] for row in rows]
        vmax2 = max(float(np.percentile(np.abs(weak), 95.0)), 1.0e-12)
        sc2 = axes[2].scatter(xs, zs, c=weak, s=5.0, cmap="coolwarm", vmin=-vmax2, vmax=vmax2)
        fig.colorbar(sc2, ax=axes[2])
    axes[1].set_title("spacetime normal flux jump")
    axes[2].set_title("spacetime weak residual")
    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    params = FlatLocalizedCrossingParams(nx=args.resolution, nz=args.resolution)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    ref = reference_snapshot_data(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=args.time,
        probe_dt=args.probe_dt,
        ell=args.ell,
        mp=args.mp,
        mass=params.m,
        rho_floor=args.rho_floor,
    )
    raw_y = args.ell * args.ell * real_array(ref["chi_ref"])
    raw_y_m = args.ell * args.ell * real_array(ref["chi_m"])
    raw_y_p = args.ell * args.ell * real_array(ref["chi_p"])
    rho = real_array(ref["rho"])
    metric_inv = real_array(ref["metric_inv"])
    stress_trace = real_array(ref["stress_trace"])
    levels = parse_positive_levels(args.levels)
    rows: list[dict[str, float]] = []
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        rows.extend(
            interface_jump_rows(
                raw_y=raw_y,
                raw_y_m=raw_y_m,
                raw_y_p=raw_y_p,
                metric_inv=metric_inv,
                rho=rho,
                stress_trace=stress_trace,
                x=x,
                z=z,
                level=level,
                ell=args.ell,
                probe_dt=args.probe_dt,
                half_width_y=args.half_width_y,
                samples=args.samples,
            )
        )

    fig_path = out / f"d_spacetime_interface_jump_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_plot(fig_path, x=x, z=z, rho=rho, raw_y=raw_y, rows=rows, levels=levels)
    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "levels_abs_raw_y": levels,
            "half_width_y": args.half_width_y,
            "samples": args.samples,
        },
        "definitions": {
            "raw_y": "raw_y=ell^2 R_tilde; interface is raw_y=+/-level",
            "normal_flux_jump": "Leading thin-layer contribution of 2*Box(f_R): 2*sigma*sqrt(|g^ab F_a F_b|)*(d_phi/d_raw_y|right - d_phi/d_raw_y|left)",
            "bulk_integral": "Integral over raw_y of (phi*R_tilde - 3/2*f - T/Mp^2)/sqrt(|g^ab F_a F_b|)",
            "weak_residual": "normal_flux_jump + bulk_integral. This is the leading spacetime-interface trace weak form, not yet the full tensor jump condition.",
            "jump_law_residual_scaled": "Equivalent jump-law residual after multiplying by sqrt(|g^ab F_a F_b|): 2*sigma*|g^ab F_a F_b|*Delta(d_phi/d_raw_y)+Integral(algebraic d raw_y). This is the form used to solve for interface speed.",
            "target_normal_norm_same_signature": "The g^ab F_a F_b value that would satisfy the leading trace jump law if the interface kept the same normal signature as the reference slice; NaN means no same-signature real solution.",
            "half_width_y": "Integrates raw_y from level-half_width_y to level+half_width_y.",
        },
        "summary": summarize(rows),
        "files": {"figure": str(fig_path.resolve())},
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolution", type=int, default=128)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=101)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
