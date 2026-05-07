from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from plot_d_interface_linear_zooms import extract_segments
from diagnose_d_trace_interface_weak_form import parse_positive_levels, real_array
from prototype_d_reference_imex import reference_snapshot_data
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def point_segment_distance_squared(
    px: np.ndarray,
    pz: np.ndarray,
    x0: float,
    z0: float,
    x1: float,
    z1: float,
) -> np.ndarray:
    vx = x1 - x0
    vz = z1 - z0
    denom = max(vx * vx + vz * vz, 1.0e-300)
    tau = ((px - x0) * vx + (pz - z0) * vz) / denom
    tau = np.clip(tau, 0.0, 1.0)
    cx = x0 + tau * vx
    cz = z0 + tau * vz
    return (px - cx) * (px - cx) + (pz - cz) * (pz - cz)


def signed_distance_to_segments(
    raw_y: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    segments: list[dict[str, float]],
    level: float,
) -> np.ndarray:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    dist2 = np.full_like(raw_y, np.inf, dtype=float)
    for seg in segments:
        d2 = point_segment_distance_squared(xg, zg, seg["x0"], seg["z0"], seg["x1"], seg["z1"])
        dist2 = np.minimum(dist2, d2)
    distance = np.sqrt(dist2)
    # Negative: EH-like side |raw_y|<level. Positive: saturated side |raw_y|>level.
    sign = np.where(np.abs(raw_y) >= level, 1.0, -1.0)
    return sign * distance


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_segments(ax: plt.Axes, rows: list[dict[str, float]], colors: str | list[str], linewidth: float = 0.85) -> None:
    if not rows:
        return
    ax.add_collection(LineCollection(row_segments(rows), colors=colors, linewidths=linewidth))


def render(
    path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    signed_distance: np.ndarray,
    segments: list[dict[str, float]],
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8), constrained_layout=True)
    colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in segments]

    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))
    im0 = axes[0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_segments(axes[0], segments, colors)
    axes[0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0].set_title("linear rho + support interface")
    fig.colorbar(im0, ax=axes[0])

    raw_clip = max(float(np.percentile(np.abs(raw_y[support]), 95.0)), 1.0) if np.any(support) else 1.0
    im1 = axes[1].pcolormesh(xg, zg, raw_y, shading="auto", cmap="coolwarm", vmin=-raw_clip, vmax=raw_clip)
    add_segments(axes[1], segments, colors)
    axes[1].set_title("raw_y = ell^2 R_tilde")
    fig.colorbar(im1, ax=axes[1])

    sd_clip = max(float(np.percentile(np.abs(signed_distance[support]), 95.0)), 1.0) if np.any(support) else 1.0
    im2 = axes[2].pcolormesh(xg, zg, signed_distance, shading="auto", cmap="coolwarm", vmin=-sd_clip, vmax=sd_clip)
    add_segments(axes[2], segments, colors)
    axes[2].set_title("signed distance to |raw_y|=1")
    fig.colorbar(im2, ax=axes[2])

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
    fig.savefig(path, dpi=180)
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
    rho = real_array(ref["rho"])
    raw_y = args.ell * args.ell * real_array(ref["chi_ref"])
    levels = parse_positive_levels(args.levels)
    if len(levels) != 1:
        raise ValueError("This signed-distance preview currently expects a single absolute level, usually 1.")
    level = levels[0]
    segments: list[dict[str, float]] = []
    for signed_level in (-level, level):
        segments.extend(extract_segments(raw_y, rho, x, z, level=signed_level, rho_min_frac=args.rho_min_frac))
    signed_distance = signed_distance_to_segments(raw_y, x, z, segments, level=level)

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    sd_x, sd_z = np.gradient(signed_distance, dx, dz, edge_order=2)
    grad_norm = np.sqrt(sd_x * sd_x + sd_z * sd_z)
    near = np.abs(signed_distance) < args.near_distance
    support = rho > args.rho_min_frac * float(np.max(rho))
    near_support = near & support

    fig_path = out / f"d_signed_distance_interface_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render(fig_path, x=x, z=z, rho=rho, raw_y=raw_y, signed_distance=signed_distance, segments=segments, rho_min_frac=args.rho_min_frac)

    vals = grad_norm[near_support]
    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolution": args.resolution,
            "rho_min_frac": args.rho_min_frac,
            "level_abs_raw_y": level,
            "near_distance": args.near_distance,
        },
        "definitions": {
            "signed_distance": "Negative where |raw_y|<level, positive where |raw_y|>=level; magnitude is the approximate Euclidean distance to the support-kept |raw_y|=level interface segments.",
            "purpose": "Numerical reinitialization preview for the moving-interface solver. This is not a new physical field.",
            "body_fitted_relation": "A body-fitted solver would place grid/finite-volume faces on these zero-distance curves instead of representing them on a background Cartesian grid.",
        },
        "stats": {
            "segment_count": len(segments),
            "near_support_count": int(np.sum(near_support)),
            "grad_signed_distance_near_support_abs_mean": float(np.mean(np.abs(vals))) if vals.size else 0.0,
            "grad_signed_distance_near_support_abs_median": float(np.median(np.abs(vals))) if vals.size else 0.0,
            "grad_signed_distance_near_support_abs_p95": float(np.percentile(np.abs(vals), 95.0)) if vals.size else 0.0,
        },
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
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    parser.add_argument("--near-distance", type=float, default=0.75)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
