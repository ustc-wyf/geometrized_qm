from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_bcd_transition_layers import make_fields, mask_centroid
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def bounding_box_from_mask(
    mask: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    margin_x: float,
    margin_z: float,
) -> tuple[float, float, float, float] | None:
    if not np.any(mask):
        return None
    ii, jj = np.where(mask)
    xmin = float(x[np.min(ii)] - margin_x)
    xmax = float(x[np.max(ii)] + margin_x)
    zmin = float(z[np.min(jj)] - margin_z)
    zmax = float(z[np.max(jj)] + margin_z)
    return xmin, xmax, zmin, zmax


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


def render_zoom_figure(
    branch: str,
    t: float,
    x: np.ndarray,
    z: np.ndarray,
    fields: dict[str, np.ndarray],
    strict_band: tuple[float, float],
    broad_band: tuple[float, float],
    out_path: Path,
) -> dict[str, object]:
    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    rho = fields["rho"]
    r_tilde = fields["R_tilde"]
    ell2r = fields["ell2R_abs"]
    derivative = fields["derivative_norm"]
    residual = fields["residual_norm"]

    support_mask = rho > 1.0e-3 * float(np.max(rho))
    dense_core_mask = rho > 1.0e-1 * float(np.max(rho))
    strict_mask = (ell2r >= strict_band[0]) & (ell2r <= strict_band[1])
    broad_mask = (ell2r >= broad_band[0]) & (ell2r <= broad_band[1])

    hotspot_threshold = float(np.percentile(residual[support_mask], 99.0)) if np.any(support_mask) else float(np.max(residual))
    hotspot_mask = residual >= hotspot_threshold

    focus_mask = (broad_mask & support_mask) | (hotspot_mask & support_mask) | dense_core_mask
    bbox = bounding_box_from_mask(focus_mask, x, z, margin_x=2.0, margin_z=2.0)
    if bbox is None:
        bbox = (-12.0, 12.0, -12.0, 12.0)
    xmin, xmax, zmin, zmax = bbox

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.2), constrained_layout=True)

    rho_vals = rho[focus_mask] if np.any(focus_mask) else rho.reshape(-1)
    rho_vals = rho_vals[np.isfinite(rho_vals)]
    rho_vmin = float(np.percentile(rho_vals, 1.0)) if rho_vals.size else 0.0
    rho_vmax = float(np.percentile(rho_vals, 99.0)) if rho_vals.size else 1.0
    if rho_vmax <= rho_vmin:
        rho_vmax = rho_vmin + 1.0

    r_vals = r_tilde[focus_mask] if np.any(focus_mask) else r_tilde.reshape(-1)
    r_vals = r_vals[np.isfinite(r_vals)]
    r_vmax = float(np.percentile(np.abs(r_vals), 99.0)) if r_vals.size else 1.0
    if r_vmax <= 0.0:
        r_vmax = 1.0

    d_vals = derivative[focus_mask] if np.any(focus_mask) else derivative.reshape(-1)
    d_vals = d_vals[np.isfinite(d_vals)]
    d_vmin = float(np.percentile(d_vals, 1.0)) if d_vals.size else 0.0
    d_vmax = float(np.percentile(d_vals, 99.0)) if d_vals.size else 1.0
    if d_vmax <= d_vmin:
        d_vmax = d_vmin + 1.0

    e_vals = residual[focus_mask] if np.any(focus_mask) else residual.reshape(-1)
    e_vals = e_vals[np.isfinite(e_vals)]
    e_vmin = float(np.percentile(e_vals, 1.0)) if e_vals.size else 0.0
    e_vmax = float(np.percentile(e_vals, 99.0)) if e_vals.size else 1.0
    if e_vmax <= e_vmin:
        e_vmax = e_vmin + 1.0

    im0 = axes[0, 0].pcolormesh(Xg, Zg, rho, shading="auto", cmap="viridis", vmin=rho_vmin, vmax=rho_vmax)
    axes[0, 0].contour(Xg, Zg, broad_mask.astype(float), levels=[0.5], colors="w", linewidths=1.0)
    axes[0, 0].contour(Xg, Zg, strict_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
    axes[0, 0].set_title(fr"{branch} branch, $t={t:.1f}$: rho")
    fig.colorbar(im0, ax=axes[0, 0], label="rho")

    im1 = axes[0, 1].pcolormesh(Xg, Zg, r_tilde, shading="auto", cmap="coolwarm", vmin=-r_vmax, vmax=r_vmax)
    axes[0, 1].contour(Xg, Zg, broad_mask.astype(float), levels=[0.5], colors="w", linewidths=1.0)
    axes[0, 1].contour(Xg, Zg, strict_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
    axes[0, 1].set_title(r"$\tilde R$")
    fig.colorbar(im1, ax=axes[0, 1], label=r"$\tilde R$")

    im2 = axes[1, 0].pcolormesh(Xg, Zg, derivative, shading="auto", cmap="inferno", vmin=d_vmin, vmax=d_vmax)
    axes[1, 0].contour(Xg, Zg, broad_mask.astype(float), levels=[0.5], colors="w", linewidths=1.0)
    axes[1, 0].contour(Xg, Zg, strict_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
    axes[1, 0].contour(Xg, Zg, hotspot_mask.astype(float), levels=[0.5], colors="lime", linewidths=1.0)
    axes[1, 0].set_title("Derivative-term norm with hotspots")
    fig.colorbar(im2, ax=axes[1, 0], label=r"$\|D(f_R)\|_F$")

    im3 = axes[1, 1].pcolormesh(Xg, Zg, residual, shading="auto", cmap="plasma", vmin=e_vmin, vmax=e_vmax)
    axes[1, 1].contour(Xg, Zg, broad_mask.astype(float), levels=[0.5], colors="w", linewidths=1.0)
    axes[1, 1].contour(Xg, Zg, strict_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
    axes[1, 1].contour(Xg, Zg, dense_core_mask.astype(float), levels=[0.5], colors="lime", linewidths=1.0, linestyles="--")
    axes[1, 1].set_title("Residual norm")
    fig.colorbar(im3, ax=axes[1, 1], label=r"$\|E\|_F$")

    for ax in axes.ravel():
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(zmin, zmax)
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")

    fig.savefig(out_path, dpi=200)
    plt.close(fig)

    strict_support = strict_mask & support_mask
    broad_support = broad_mask & support_mask
    hotspot_support = hotspot_mask & support_mask
    return {
        "time": float(t),
        "strict_band": [float(strict_band[0]), float(strict_band[1])],
        "broad_band": [float(broad_band[0]), float(broad_band[1])],
        "strict_count_global": int(np.sum(strict_mask)),
        "strict_count_support": int(np.sum(strict_support)),
        "broad_count_global": int(np.sum(broad_mask)),
        "broad_count_support": int(np.sum(broad_support)),
        "hotspot_count_support": int(np.sum(hotspot_support)),
        "strict_fraction_global": float(np.mean(strict_mask)),
        "strict_fraction_support": float(np.mean(strict_support[support_mask])) if np.any(support_mask) else 0.0,
        "broad_fraction_global": float(np.mean(broad_mask)),
        "broad_fraction_support": float(np.mean(broad_support[support_mask])) if np.any(support_mask) else 0.0,
        "hotspot_inside_strict_fraction": float(np.mean((hotspot_mask & strict_mask)[hotspot_mask])) if np.any(hotspot_mask) else 0.0,
        "hotspot_inside_broad_fraction": float(np.mean((hotspot_mask & broad_mask)[hotspot_mask])) if np.any(hotspot_mask) else 0.0,
        "strict_centroid_global": mask_centroid(strict_mask, x, z),
        "strict_centroid_support": mask_centroid(strict_support, x, z),
        "broad_centroid_global": mask_centroid(broad_mask, x, z),
        "broad_centroid_support": mask_centroid(broad_support, x, z),
        "hotspot_centroid_support": mask_centroid(hotspot_support, x, z),
        "hotspot_threshold": hotspot_threshold,
        "focus_bbox": [xmin, xmax, zmin, zmax],
        "residual_norm_strict_support": scalar_stats(residual, strict_support),
        "residual_norm_broad_support": scalar_stats(residual, broad_support),
        "figure": str(out_path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--nx", type=int, default=192)
    parser.add_argument("--nz", type=int, default=192)
    parser.add_argument("--times", type=str, default="8,16")
    parser.add_argument("--strict-band", type=str, default="0.5,2.0")
    parser.add_argument("--broad-band", type=str, default="1,100")
    args = parser.parse_args()

    strict_band = tuple(float(tok) for tok in args.strict_band.split(","))
    broad_band = tuple(float(tok) for tok in args.broad_band.split(","))
    params = FlatLocalizedCrossingParams(nx=args.nx, nz=args.nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    times = [float(tok) for tok in args.times.split(",") if tok.strip()]

    args.output.mkdir(parents=True, exist_ok=True)
    records: dict[str, list[dict[str, object]]] = {"C": [], "D": []}
    for branch in ("C", "D"):
        for t in times:
            fields = make_fields(
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
            rec = render_zoom_figure(
                branch=branch,
                t=t,
                x=x,
                z=z,
                fields=fields,
                strict_band=strict_band,
                broad_band=broad_band,
                out_path=args.output / f"{branch.lower()}_zoom_t{str(t).replace('.', 'p')}.png",
            )
            records[branch].append(rec)

    summary = {
        "params": {
            "ell": args.ell,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "nx": args.nx,
            "nz": args.nz,
            "times": times,
            "strict_band": list(strict_band),
            "broad_band": list(broad_band),
        },
        "definitions": {
            "support_mask": "rho > 1e-3 * rho_max",
            "dense_core_mask": "rho > 1e-1 * rho_max",
            "strict_transition_mask": f"{strict_band[0]} <= |ell^2 R_tilde| <= {strict_band[1]}",
            "broad_transition_band": f"{broad_band[0]} <= |ell^2 R_tilde| <= {broad_band[1]}",
            "hotspot_mask": "residual_norm in top 1% within support_mask",
        },
        "records": records,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
