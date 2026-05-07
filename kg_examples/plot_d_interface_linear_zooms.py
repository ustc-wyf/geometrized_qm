from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from diagnose_d_trace_interface_weak_form import parse_positive_levels, real_array
from prototype_d_reference_imex import reference_snapshot_data
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


DEFAULT_WINDOWS = {
    "global": (-20.0, 20.0, -20.0, 20.0),
    "upper_transition_fan": (-8.0, 8.0, 1.0, 11.0),
    "left_leaf": (-8.5, -0.5, 1.0, 7.0),
    "right_leaf": (0.5, 8.5, 1.0, 7.0),
    "central_tangle": (-3.5, 3.5, 3.0, 10.5),
    "low_density_edge": (-10.0, 10.0, -2.0, 2.0),
}


def bilinear_sample(field: np.ndarray, x: np.ndarray, z: np.ndarray, xp: float, zp: float) -> float:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    u = (xp - float(x[0])) / dx
    v = (zp - float(z[0])) / dz
    i0 = int(np.clip(np.floor(u), 0, len(x) - 2))
    j0 = int(np.clip(np.floor(v), 0, len(z) - 2))
    tx = float(np.clip(u - i0, 0.0, 1.0))
    tz = float(np.clip(v - j0, 0.0, 1.0))
    return float(
        (1.0 - tx) * (1.0 - tz) * field[i0, j0]
        + tx * (1.0 - tz) * field[i0 + 1, j0]
        + (1.0 - tx) * tz * field[i0, j0 + 1]
        + tx * tz * field[i0 + 1, j0 + 1]
    )


def extract_segments(
    raw_y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    rho_min_frac: float,
) -> list[dict[str, float]]:
    rho_threshold = rho_min_frac * float(np.max(rho))
    segments: list[dict[str, float]] = []
    for i in range(len(x) - 1):
        for j in range(len(z) - 1):
            cell_y = raw_y[i : i + 2, j : j + 2]
            if not (np.nanmin(cell_y) <= level <= np.nanmax(cell_y)):
                continue
            corners = [
                (float(x[i]), float(z[j]), float(raw_y[i, j])),
                (float(x[i + 1]), float(z[j]), float(raw_y[i + 1, j])),
                (float(x[i + 1]), float(z[j + 1]), float(raw_y[i + 1, j + 1])),
                (float(x[i]), float(z[j + 1]), float(raw_y[i, j + 1])),
            ]
            crossings: list[tuple[float, float]] = []
            for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
                x0, z0, y0 = corners[a]
                x1, z1, y1 = corners[b]
                dy = y1 - y0
                if abs(dy) <= 1.0e-300:
                    continue
                tau = (level - y0) / dy
                if 0.0 <= tau <= 1.0:
                    crossings.append((x0 + tau * (x1 - x0), z0 + tau * (z1 - z0)))
            if len(crossings) < 2:
                continue
            p0, p1 = crossings[0], crossings[1]
            if len(crossings) > 2:
                best_len = -1.0
                for a in range(len(crossings)):
                    for b in range(a + 1, len(crossings)):
                        dist = float(np.hypot(crossings[a][0] - crossings[b][0], crossings[a][1] - crossings[b][1]))
                        if dist > best_len:
                            best_len = dist
                            p0, p1 = crossings[a], crossings[b]
            xp = 0.5 * (p0[0] + p1[0])
            zp = 0.5 * (p0[1] + p1[1])
            rho_mid = bilinear_sample(rho, x, z, xp, zp)
            if rho_mid < rho_threshold:
                continue
            length = float(np.hypot(p1[0] - p0[0], p1[1] - p0[1]))
            if length <= 0.0:
                continue
            segments.append(
                {
                    "x0": float(p0[0]),
                    "z0": float(p0[1]),
                    "x1": float(p1[0]),
                    "z1": float(p1[1]),
                    "x": float(xp),
                    "z": float(zp),
                    "rho_mid": float(rho_mid),
                    "level": float(level),
                    "length": length,
                }
            )
    return segments


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_segments(ax: plt.Axes, rows: list[dict[str, float]], colors: str | list[str], linewidth: float, alpha: float, linestyle: str = "solid") -> None:
    if not rows:
        return
    lc = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth, alpha=alpha, linestyles=linestyle)
    ax.add_collection(lc)


def crop_rows(rows: list[dict[str, float]], bbox: tuple[float, float, float, float]) -> list[dict[str, float]]:
    xmin, xmax, zmin, zmax = bbox
    return [row for row in rows if xmin <= row["x"] <= xmax and zmin <= row["z"] <= zmax]


def render_window(
    path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    all_rows: list[dict[str, float]],
    support_rows: list[dict[str, float]],
    bbox: tuple[float, float, float, float],
    title: str,
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    xmin, xmax, zmin, zmax = bbox
    support = rho > rho_min_frac * float(np.max(rho))

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[(xg >= xmin) & (xg <= xmax) & (zg >= zmin) & (zg <= zmax)], 99.5)), 1.0e-16)
    im0 = axes[0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_segments(axes[0], crop_rows(all_rows, bbox), colors="0.15", linewidth=0.55, alpha=0.45)
    axes[0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0].set_title("linear rho + all raw_y=+/-1 contours")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    support_crop = crop_rows(support_rows, bbox)
    colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in support_crop]
    add_segments(axes[1], support_crop, colors=colors, linewidth=0.95, alpha=0.95)
    axes[1].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[1].set_title("linear rho + support-kept interface")
    fig.colorbar(im1, ax=axes[1])

    raw_clip = max(float(np.percentile(np.abs(raw_y[(xg >= xmin) & (xg <= xmax) & (zg >= zmin) & (zg <= zmax)]), 95.0)), 1.0)
    im2 = axes[2].pcolormesh(xg, zg, raw_y, shading="auto", cmap="coolwarm", vmin=-raw_clip, vmax=raw_clip)
    add_segments(axes[2], support_crop, colors=colors, linewidth=0.9, alpha=0.95)
    axes[2].set_title("raw_y field + support-kept interface")
    fig.colorbar(im2, ax=axes[2])

    for ax in axes:
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(zmin, zmax)
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
    fig.suptitle(f"{title}; dashed white: rho>{rho_min_frac:g} rho_max", fontsize=12)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def stats_for_rows(rows: list[dict[str, float]]) -> dict[str, float]:
    if not rows:
        return {"segment_count": 0, "total_length": 0.0, "rho_mid_median": 0.0, "rho_mid_p05": 0.0}
    rhos = np.array([row["rho_mid"] for row in rows], dtype=float)
    return {
        "segment_count": len(rows),
        "total_length": float(sum(row["length"] for row in rows)),
        "rho_mid_median": float(np.median(rhos)),
        "rho_mid_p05": float(np.percentile(rhos, 5.0)),
    }


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
    signed_levels = [sign * lev for lev in levels for sign in (-1.0, 1.0)]
    all_rows: list[dict[str, float]] = []
    support_rows: list[dict[str, float]] = []
    for level in signed_levels:
        all_rows.extend(extract_segments(raw_y, rho, x, z, level=level, rho_min_frac=0.0))
        support_rows.extend(extract_segments(raw_y, rho, x, z, level=level, rho_min_frac=args.rho_min_frac))

    files: dict[str, str] = {}
    for name, bbox in DEFAULT_WINDOWS.items():
        path = out / f"d_interface_linear_{name}_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
        render_window(
            path,
            x=x,
            z=z,
            rho=rho,
            raw_y=raw_y,
            all_rows=all_rows,
            support_rows=support_rows,
            bbox=bbox,
            title=name,
            rho_min_frac=args.rho_min_frac,
        )
        files[name] = str(path.resolve())

    support_fraction = len(support_rows) / max(len(all_rows), 1)
    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "rho_min_frac": args.rho_min_frac,
            "levels_abs_raw_y": levels,
        },
        "definitions": {
            "raw_y": "raw_y=ell^2 R_tilde. Orange lines are raw_y=+level and green lines are raw_y=-level in support-kept panels.",
            "all_contours": "All marching-square raw_y=+/-level segments, no rho support threshold.",
            "support_kept_interface": "Segments whose midpoint satisfies rho >= rho_min_frac*rho_max.",
            "dashed_white": "Boundary of support rho>rho_min_frac*rho_max.",
            "linear_rho": "rho plotted on a linear color scale, not log scale.",
        },
        "segment_stats": {
            "all": stats_for_rows(all_rows),
            "support_kept": stats_for_rows(support_rows),
            "support_kept_fraction_by_segment_count": float(support_fraction),
            "removed_by_support_threshold": int(len(all_rows) - len(support_rows)),
        },
        "window_stats": {
            name: {
                "all": stats_for_rows(crop_rows(all_rows, bbox)),
                "support_kept": stats_for_rows(crop_rows(support_rows, bbox)),
            }
            for name, bbox in DEFAULT_WINDOWS.items()
        },
        "files": files,
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolution", type=int, default=192)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
