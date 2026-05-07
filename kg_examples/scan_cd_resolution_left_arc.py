from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np

from plot_bcd_transition_layers import make_fields
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


LEFT_ARC_BBOX = (-7.4, -2.2, 2.6, 6.2)


def mask_centroid(mask: np.ndarray, xg: np.ndarray, zg: np.ndarray) -> tuple[float, float] | None:
    if not np.any(mask):
        return None
    w = mask.astype(float)
    norm = float(np.sum(w))
    return float(np.sum(w * xg) / norm), float(np.sum(w * zg) / norm)


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


def signed_limits(field: np.ndarray, mask: np.ndarray | None = None, qhi: float = 99.0) -> tuple[float, float]:
    vals = field[mask] if mask is not None else field.reshape(-1)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return -1.0, 1.0
    vmax = float(np.percentile(np.abs(vals), qhi))
    vmax = max(vmax, 1.0e-12)
    return -vmax, vmax


def positive_limits(field: np.ndarray, mask: np.ndarray | None = None, qlo: float = 1.0, qhi: float = 99.0) -> tuple[float, float]:
    vals = field[mask] if mask is not None else field.reshape(-1)
    vals = vals[np.isfinite(vals)]
    vals = vals[vals >= 0.0]
    if vals.size == 0:
        return 0.0, 1.0
    vmin = float(np.percentile(vals, qlo))
    vmax = float(np.percentile(vals, qhi))
    if vmax <= vmin:
        vmax = max(float(np.max(vals)), vmin + 1.0)
    return vmin, vmax


def compute_branch_resolution_record(
    branch: str,
    nx: int,
    nz: int,
    t: float,
    ell: float,
    mp: float,
    probe_dt: float,
    rho_floor: float,
    strict_band: tuple[float, float],
    broad_band: tuple[float, float],
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    params = FlatLocalizedCrossingParams(nx=nx, nz=nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    fields = make_fields(
        branch=branch,
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=t,
        probe_dt=probe_dt,
        ell=ell,
        mp=mp,
        mass=params.m,
        rho_floor=rho_floor,
    )

    xg, zg = np.meshgrid(x, z, indexing="ij")
    xmin, xmax, zmin, zmax = LEFT_ARC_BBOX
    win = (xg >= xmin) & (xg <= xmax) & (zg >= zmin) & (zg <= zmax)

    rho = fields["rho"]
    residual = fields["residual_norm"]
    r_tilde = fields["R_tilde"]
    ell2r = fields["ell2R_abs"]

    support = rho > 1.0e-3 * float(np.max(rho))
    dense = rho > 1.0e-1 * float(np.max(rho))
    strict = (ell2r >= strict_band[0]) & (ell2r <= strict_band[1])
    broad = (ell2r >= broad_band[0]) & (ell2r <= broad_band[1])
    hotspot_thr = float(np.percentile(residual[support], 99.0)) if np.any(support) else float(np.max(residual))
    hotspot = residual >= hotspot_thr

    win_support = win & support
    win_dense = win & dense
    win_strict = win & support & strict
    win_broad = win & support & broad
    win_hotspot = win & support & hotspot

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    area = dx * dz
    l1_abs = float(np.sum(np.abs(residual[win_support])) * area) if np.any(win_support) else 0.0

    max_idx = np.unravel_index(np.argmax(np.where(win_support, np.abs(residual), -1.0)), residual.shape) if np.any(win_support) else None
    max_point = (float(xg[max_idx]), float(zg[max_idx])) if max_idx is not None else None

    record = {
        "branch": branch,
        "nx": int(nx),
        "nz": int(nz),
        "time": float(t),
        "bbox": list(LEFT_ARC_BBOX),
        "definitions": {
            "support": "rho > 1e-3 * rho_max",
            "dense": "rho > 1e-1 * rho_max",
            "strict_transition": f"{strict_band[0]} <= |ell^2 R_tilde| <= {strict_band[1]}",
            "broad_transition": f"{broad_band[0]} <= |ell^2 R_tilde| <= {broad_band[1]}",
            "hotspot": "residual_norm in top 1% within support",
        },
        "window_counts": {
            "support": int(np.sum(win_support)),
            "dense": int(np.sum(win_dense)),
            "strict_support": int(np.sum(win_strict)),
            "broad_support": int(np.sum(win_broad)),
            "hotspot_support": int(np.sum(win_hotspot)),
        },
        "fractions": {
            "hotspot_inside_broad": float(np.sum(win_hotspot & broad) / np.sum(win_hotspot)) if np.any(win_hotspot) else 0.0,
            "hotspot_inside_strict": float(np.sum(win_hotspot & strict) / np.sum(win_hotspot)) if np.any(win_hotspot) else 0.0,
            "broad_in_support": float(np.sum(win_broad) / np.sum(win_support)) if np.any(win_support) else 0.0,
            "strict_in_support": float(np.sum(win_strict) / np.sum(win_support)) if np.any(win_support) else 0.0,
        },
        "residual_window_support": support_stats(residual, win_support),
        "residual_window_l1_abs": l1_abs,
        "max_residual_point_in_window_support": max_point,
        "centroids": {
            "broad_support": mask_centroid(win_broad, xg, zg),
            "hotspot_support": mask_centroid(win_hotspot, xg, zg),
        },
        "grid_spacing": {"dx": dx, "dz": dz},
    }
    return record, {
        "xg": xg,
        "zg": zg,
        "win": win,
        "support": support,
        "strict": strict,
        "broad": broad,
        "hotspot": hotspot,
        **fields,
    }


def render_branch_resolution_figure(
    branch: str,
    rows: list[tuple[int, dict[str, np.ndarray]]],
    out_path: Path,
) -> None:
    nrows = len(rows)
    fig, axes = plt.subplots(nrows, 4, figsize=(16.0, 4.0 * nrows), constrained_layout=True)
    if nrows == 1:
        axes = np.array([axes])

    note = (
        "Window = left arc: -7.4 <= x <= -2.2, 2.6 <= z <= 6.2. "
        "White contour = broad transition band 1 <= |ell^2 R_tilde| <= 100. "
        "Cyan contour = strict transition layer 0.5 <= |ell^2 R_tilde| <= 2. "
        "Lime contour = residual hotspot = top 1% within support. "
        "Support means rho > 1e-3 * rho_max."
    )
    fig.text(0.02, 0.99, note, va="top", fontsize=9)

    for row_idx, (nx, data) in enumerate(rows):
        xg = data["xg"]
        zg = data["zg"]
        win = data["win"]
        rho = data["rho"]
        r_tilde = data["R_tilde"]
        residual = data["residual_norm"]
        support = data["support"]
        strict = data["strict"]
        broad = data["broad"]
        hotspot = data["hotspot"]

        win_support = win & support

        rho_vmin, rho_vmax = positive_limits(rho, win_support)
        r_vmin, r_vmax = signed_limits(r_tilde, win_support)
        e_vmin, e_vmax = positive_limits(residual, win_support)
        max_abs = max(float(np.nanmax(np.abs(residual[win_support]))) if np.any(win_support) else 1.0, 1.0e-12)
        linthresh = max(1.0e-8, 1.0e-3 * max_abs)

        panels = [
            ("rho (linear)", rho, {"cmap": "viridis", "vmin": rho_vmin, "vmax": rho_vmax}),
            (r"$\tilde R$ (linear)", r_tilde, {"cmap": "coolwarm", "vmin": r_vmin, "vmax": r_vmax}),
            ("residual (linear)", residual, {"cmap": "plasma", "vmin": e_vmin, "vmax": e_vmax}),
            (
                "residual (symlog)",
                residual,
                {"cmap": "plasma", "norm": SymLogNorm(linthresh=linthresh, vmin=-max_abs, vmax=max_abs)},
            ),
        ]

        for col_idx, (title, field, kwargs) in enumerate(panels):
            ax = axes[row_idx, col_idx]
            im = ax.pcolormesh(xg, zg, field, shading="auto", **kwargs)
            ax.contour(xg, zg, broad.astype(float), levels=[0.5], colors="white", linewidths=0.9)
            ax.contour(xg, zg, strict.astype(float), levels=[0.5], colors="cyan", linewidths=1.0)
            if "residual" in title:
                ax.contour(xg, zg, hotspot.astype(float), levels=[0.5], colors="lime", linewidths=0.9)
            xmin, xmax, zmin, zmax = LEFT_ARC_BBOX
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(zmin, zmax)
            ax.set_aspect("equal")
            ax.set_xlabel("x")
            ax.set_ylabel("z")
            ax.set_title(f"{branch}, {nx}x{nx}: {title}")
            fig.colorbar(im, ax=ax, shrink=0.85)

    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--resolutions", type=str, default="64,96,128,192")
    parser.add_argument("--branches", type=str, default="C,D")
    parser.add_argument("--strict-band", type=str, default="0.5,2.0")
    parser.add_argument("--broad-band", type=str, default="1,100")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    resolutions = [int(tok) for tok in args.resolutions.split(",") if tok.strip()]
    branches = [tok.strip().upper() for tok in args.branches.split(",") if tok.strip()]
    strict_band = tuple(float(tok) for tok in args.strict_band.split(","))
    broad_band = tuple(float(tok) for tok in args.broad_band.split(","))

    summary: dict[str, object] = {
        "params": {
            "time": args.time,
            "ell": args.ell,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "resolutions": resolutions,
            "bbox": list(LEFT_ARC_BBOX),
            "strict_band": list(strict_band),
            "broad_band": list(broad_band),
        },
        "records": {branch: [] for branch in branches},
        "files": {},
    }

    for branch in branches:
        plot_rows: list[tuple[int, dict[str, np.ndarray]]] = []
        for n in resolutions:
            record, plot_data = compute_branch_resolution_record(
                branch=branch,
                nx=n,
                nz=n,
                t=args.time,
                ell=args.ell,
                mp=args.mp,
                probe_dt=args.probe_dt,
                rho_floor=args.rho_floor,
                strict_band=strict_band,
                broad_band=broad_band,
            )
            summary["records"][branch].append(record)
            plot_rows.append((n, plot_data))
        fig_path = args.output / f"{branch.lower()}_left_arc_resolution_t{str(args.time).replace('.', 'p')}.png"
        render_branch_resolution_figure(branch, plot_rows, fig_path)
        summary["files"][branch] = str(fig_path.resolve())

    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
