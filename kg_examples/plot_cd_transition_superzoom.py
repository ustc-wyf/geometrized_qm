from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_bcd_transition_layers import make_fields
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def _color_limits_positive(field: np.ndarray, mask: np.ndarray | None = None, qlo: float = 1.0, qhi: float = 99.0) -> tuple[float, float]:
    vals = field[mask] if mask is not None else field.reshape(-1)
    vals = vals[np.isfinite(vals)]
    vals = vals[vals >= 0.0]
    if vals.size == 0:
        return 0.0, 1.0
    vmin = float(np.percentile(vals, qlo))
    vmax = float(np.percentile(vals, qhi))
    if vmax <= vmin:
        vmax = float(np.max(vals))
        vmin = float(np.min(vals))
        if vmax <= vmin:
            vmax = vmin + 1.0
    return vmin, vmax


def _color_limits_signed(field: np.ndarray, mask: np.ndarray | None = None, qhi: float = 99.0) -> tuple[float, float]:
    vals = field[mask] if mask is not None else field.reshape(-1)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return -1.0, 1.0
    vmax = float(np.percentile(np.abs(vals), qhi))
    if vmax <= 0.0:
        vmax = 1.0
    return -vmax, vmax


def render_superzoom_figure(
    branch: str,
    t: float,
    x: np.ndarray,
    z: np.ndarray,
    fields: dict[str, np.ndarray],
    strict_band: tuple[float, float],
    broad_band: tuple[float, float],
    windows: list[tuple[str, tuple[float, float, float, float]]],
    out_path: Path,
) -> dict[str, object]:
    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    rho = fields["rho"]
    r_tilde = fields["R_tilde"]
    ell2r = fields["ell2R_abs"]
    derivative = fields["derivative_norm"]
    residual = fields["residual_norm"]

    support = rho > 1.0e-3 * float(np.max(rho))
    strict = (ell2r >= strict_band[0]) & (ell2r <= strict_band[1])
    broad = (ell2r >= broad_band[0]) & (ell2r <= broad_band[1])
    hotspot_thr = float(np.percentile(residual[support], 99.0)) if np.any(support) else float(np.max(residual))
    hotspot = residual >= hotspot_thr

    summary = {"hotspot_threshold": hotspot_thr, "windows": {}}
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 12.8), constrained_layout=True)

    for row, (name, (xmin, xmax, zmin, zmax)) in enumerate(windows):
        win = (Xg >= xmin) & (Xg <= xmax) & (Zg >= zmin) & (Zg <= zmax)
        strict_support = strict & support & win
        broad_support = broad & support & win
        hotspot_support = hotspot & support & win
        hotspot_inside_broad = hotspot & broad & support & win
        hotspot_inside_strict = hotspot & strict & support & win

        summary["windows"][name] = {
            "bbox": [xmin, xmax, zmin, zmax],
            "strict_support_count": int(np.sum(strict_support)),
            "broad_support_count": int(np.sum(broad_support)),
            "hotspot_support_count": int(np.sum(hotspot_support)),
            "hotspot_inside_broad_count": int(np.sum(hotspot_inside_broad)),
            "hotspot_inside_strict_count": int(np.sum(hotspot_inside_strict)),
            "hotspot_inside_broad_fraction": float(np.sum(hotspot_inside_broad) / np.sum(hotspot_support)) if np.any(hotspot_support) else 0.0,
            "hotspot_inside_strict_fraction": float(np.sum(hotspot_inside_strict) / np.sum(hotspot_support)) if np.any(hotspot_support) else 0.0,
            "residual_p95_support": float(np.percentile(residual[support & win], 95.0)) if np.any(support & win) else 0.0,
            "residual_max_support": float(np.max(residual[support & win])) if np.any(support & win) else 0.0,
        }

        rho_vmin, rho_vmax = _color_limits_positive(rho, win)
        r_vmin, r_vmax = _color_limits_signed(r_tilde, win)
        d_vmin, d_vmax = _color_limits_positive(derivative, win)
        e_vmin, e_vmax = _color_limits_positive(residual, win)

        ax = axes[row, 0]
        im = ax.pcolormesh(Xg, Zg, rho, shading="auto", cmap="viridis", vmin=rho_vmin, vmax=rho_vmax)
        ax.contour(Xg, Zg, broad.astype(float), levels=[0.5], colors="w", linewidths=1.0)
        ax.contour(Xg, Zg, strict.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
        ax.set_title(f"{branch} {name}: rho")
        fig.colorbar(im, ax=ax, label="rho")

        ax = axes[row, 1]
        im = ax.pcolormesh(Xg, Zg, r_tilde, shading="auto", cmap="coolwarm", vmin=r_vmin, vmax=r_vmax)
        ax.contour(Xg, Zg, broad.astype(float), levels=[0.5], colors="w", linewidths=1.0)
        ax.contour(Xg, Zg, strict.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
        ax.set_title(f"{branch} {name}: R_tilde")
        fig.colorbar(im, ax=ax, label="R_tilde")

        ax = axes[row, 2]
        im = ax.pcolormesh(Xg, Zg, residual, shading="auto", cmap="plasma", vmin=e_vmin, vmax=e_vmax)
        ax.contour(Xg, Zg, broad.astype(float), levels=[0.5], colors="w", linewidths=1.0)
        ax.contour(Xg, Zg, strict.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
        ax.contour(Xg, Zg, hotspot.astype(float), levels=[0.5], colors="lime", linewidths=1.0)
        ax.set_title(f"{branch} {name}: residual")
        fig.colorbar(im, ax=ax, label="||E||_F")

        for col in range(3):
            ax = axes[row, col]
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(zmin, zmax)
            ax.set_xlabel("x")
            ax.set_ylabel("z")
            ax.set_aspect("equal")

    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    summary["figure"] = str(out_path.resolve())
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--nx", type=int, default=256)
    parser.add_argument("--nz", type=int, default=256)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--strict-band", type=str, default="0.5,2.0")
    parser.add_argument("--broad-band", type=str, default="1,100")
    args = parser.parse_args()

    strict_band = tuple(float(tok) for tok in args.strict_band.split(","))
    broad_band = tuple(float(tok) for tok in args.broad_band.split(","))
    windows = [
        ("left_arc", (-7.4, -2.2, 2.6, 6.2)),
        ("center_column", (-2.4, 2.4, 0.0, 10.2)),
        ("right_arc", (2.2, 7.4, 2.6, 6.2)),
    ]

    params = FlatLocalizedCrossingParams(nx=args.nx, nz=args.nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    args.output.mkdir(parents=True, exist_ok=True)

    summary = {
        "params": {
            "ell": args.ell,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "nx": args.nx,
            "nz": args.nz,
            "time": args.time,
            "strict_band": list(strict_band),
            "broad_band": list(broad_band),
        },
        "definitions": {
            "support_mask": "rho > 1e-3 * rho_max",
            "strict_transition_mask": f"{strict_band[0]} <= |ell^2 R_tilde| <= {strict_band[1]}",
            "broad_transition_band": f"{broad_band[0]} <= |ell^2 R_tilde| <= {broad_band[1]}",
            "hotspot_mask": "residual_norm in top 1% within support_mask",
        },
        "records": {},
    }

    for branch in ("C", "D"):
        fields = make_fields(
            branch=branch,
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
        summary["records"][branch] = render_superzoom_figure(
            branch,
            args.time,
            x,
            z,
            fields,
            strict_band,
            broad_band,
            windows,
            args.output / f"{branch.lower()}_superzoom_t{str(args.time).replace('.', 'p')}.png",
        )

    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
