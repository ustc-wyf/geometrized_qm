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


def evaluate_interface_rows(
    raw_y: np.ndarray,
    raw_y_m: np.ndarray,
    raw_y_p: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    probe_dt: float,
) -> list[dict[str, float]]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_t = (raw_y_p - raw_y_m) / (2.0 * probe_dt)
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)

    rows: list[dict[str, float]] = []
    for seg in find_level_segments(raw_y, rho, x, z, level):
        xp = np.array([seg["x"]])
        zp = np.array([seg["z"]])
        ft = float(bilinear_sample(raw_y_t, x, z, xp, zp)[0])
        fx = float(bilinear_sample(raw_y_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(raw_y_z, x, z, xp, zp)[0])
        grad_spatial = float(np.hypot(fx, fz))
        v_normal = float(-ft / max(grad_spatial, 1.0e-300))

        ginv_here = np.zeros((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                ginv_here[a, b] = float(bilinear_sample(metric_inv[..., a, b], x, z, xp, zp)[0])
        covector = np.array([ft, fx, fz], dtype=float)
        normal_norm = float(covector @ ginv_here @ covector)
        rows.append(
            {
                "x": float(seg["x"]),
                "z": float(seg["z"]),
                "level": float(level),
                "segment_length": float(seg["length"]),
                "raw_y_t": ft,
                "raw_y_x": fx,
                "raw_y_z": fz,
                "spatial_grad_norm": grad_spatial,
                "coordinate_normal_speed": v_normal,
                "normal_norm_tilde_inverse": normal_norm,
            }
        )
    return rows


def sign_fractions(values: list[float]) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"positive": 0.0, "negative": 0.0, "near_zero": 0.0}
    scale = max(float(np.percentile(np.abs(vals), 95.0)), 1.0e-12)
    eps = 1.0e-8 * scale
    return {
        "positive": float(np.mean(vals > eps)),
        "negative": float(np.mean(vals < -eps)),
        "near_zero": float(np.mean(np.abs(vals) <= eps)),
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

    valid = rows
    if valid:
        xs = [row["x"] for row in valid]
        zs = [row["z"] for row in valid]
        speed = [row["coordinate_normal_speed"] for row in valid]
        vmax = max(float(np.percentile(np.abs(speed), 95.0)), 1.0e-12)
        sc1 = axes[1].scatter(xs, zs, c=speed, s=5.0, cmap="coolwarm", vmin=-vmax, vmax=vmax)
        fig.colorbar(sc1, ax=axes[1])
        norm = [row["normal_norm_tilde_inverse"] for row in valid]
        vmax2 = max(float(np.percentile(np.abs(norm), 95.0)), 1.0e-12)
        sc2 = axes[2].scatter(xs, zs, c=norm, s=5.0, cmap="coolwarm", vmin=-vmax2, vmax=vmax2)
        fig.colorbar(sc2, ax=axes[2])
    axes[1].set_title("coordinate normal speed -F_t/|grad F|")
    axes[2].set_title("tilde inverse norm g^ab F_a F_b")

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
    levels = parse_positive_levels(args.levels)
    rows: list[dict[str, float]] = []
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        rows.extend(
            evaluate_interface_rows(
                raw_y=raw_y,
                raw_y_m=raw_y_m,
                raw_y_p=raw_y_p,
                metric_inv=metric_inv,
                rho=rho,
                x=x,
                z=z,
                level=level,
                probe_dt=args.probe_dt,
            )
        )

    fig_path = out / f"d_interface_spacetime_kinematics_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_plot(fig_path, x=x, z=z, rho=rho, raw_y=raw_y, rows=rows, levels=levels)
    speeds = [row["coordinate_normal_speed"] for row in rows]
    norms = [row["normal_norm_tilde_inverse"] for row in rows]
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
        },
        "definitions": {
            "raw_y": "raw_y=ell^2 R_tilde; interface is raw_y=+/-level",
            "coordinate_normal_speed": "v_n=-F_t/sqrt(F_x^2+F_z^2) for F=raw_y-level on the t,x,z coordinate grid",
            "normal_norm_tilde_inverse": "g_tilde^{ab} F_a F_b. With +-- signature, positive means timelike normal/spacelike interface; negative means spacelike normal/timelike interface.",
        },
        "counts": {
            "rows": len(rows),
        },
        "stats": {
            "coordinate_normal_speed": scalar_stats(speeds),
            "normal_norm_tilde_inverse": scalar_stats(norms),
            "normal_norm_sign_fractions": sign_fractions(norms),
            "spatial_grad_norm": scalar_stats([row["spatial_grad_norm"] for row in rows]),
            "raw_y_t": scalar_stats([row["raw_y_t"] for row in rows]),
        },
        "files": {
            "figure": str(fig_path.resolve()),
        },
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
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
