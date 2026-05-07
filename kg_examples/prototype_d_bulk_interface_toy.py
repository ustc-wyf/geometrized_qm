from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from integrate_d_interface_jumps import (
    bilinear_sample,
    integrate_lines,
    scalar_stats,
    spatial_derivatives,
)
from plot_bcd_transition_layers import make_fields
from scan_cd_ell_window import branch_f_derivatives
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def parse_positive_levels(raw: str) -> list[float]:
    levels = sorted({abs(float(part.strip())) for part in raw.split(",") if part.strip()})
    if not levels or levels[0] <= 0.0:
        raise ValueError("--levels must contain positive numbers, e.g. 0.5,1,2")
    return levels


def levels_slug(levels: list[float]) -> str:
    return "_".join(f"{level:g}".replace(".", "p") for level in levels)


def parse_bbox(raw: str) -> tuple[float, float, float, float] | None:
    if raw.strip().lower() in {"global", "all", "none"}:
        return None
    parts = [float(part.strip()) for part in raw.split(",") if part.strip()]
    if len(parts) != 4:
        raise ValueError("--bbox must be 'global' or four numbers: xmin,xmax,zmin,zmax")
    xmin, xmax, zmin, zmax = parts
    if xmin >= xmax or zmin >= zmax:
        raise ValueError("--bbox must satisfy xmin < xmax and zmin < zmax")
    return xmin, xmax, zmin, zmax


def dilate_mask(mask: np.ndarray, radius: int) -> np.ndarray:
    out = mask.copy()
    for _ in range(radius):
        padded = np.pad(out, 1, mode="constant", constant_values=False)
        nxt = np.zeros_like(out)
        for di in range(3):
            for dj in range(3):
                nxt |= padded[di : di + out.shape[0], dj : dj + out.shape[1]]
        out = nxt
    return out


def add_to_grid(accum: np.ndarray, x: np.ndarray, z: np.ndarray, xp: float, zp: float, value: float) -> None:
    i = int(np.argmin(np.abs(x - xp)))
    j = int(np.argmin(np.abs(z - zp)))
    if 0 <= i < accum.shape[0] and 0 <= j < accum.shape[1]:
        accum[i, j] += value


def find_interface_points_for_bbox(
    y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    bbox: tuple[float, float, float, float] | None,
) -> list[tuple[float, float]]:
    if bbox is None:
        xmin, xmax, zmin, zmax = float(x[0]), float(x[-1]), float(z[0]), float(z[-1])
    else:
        xmin, xmax, zmin, zmax = bbox

    rho_threshold = 1.0e-3 * float(np.max(rho))
    y_der = spatial_derivatives(y, x, z)
    points: list[tuple[float, float]] = []

    for i in range(len(x) - 1):
        for j in range(len(z) - 1):
            cell_y = y[i : i + 2, j : j + 2]
            if not (np.nanmin(cell_y) <= level <= np.nanmax(cell_y)):
                continue
            corners = [
                (float(x[i]), float(z[j]), float(y[i, j])),
                (float(x[i + 1]), float(z[j]), float(y[i + 1, j])),
                (float(x[i + 1]), float(z[j + 1]), float(y[i + 1, j + 1])),
                (float(x[i]), float(z[j + 1]), float(y[i, j + 1])),
            ]
            edge_points: list[tuple[float, float]] = []
            for a, b in ((0, 1), (1, 2), (2, 3), (3, 0)):
                x0, z0, y0 = corners[a]
                x1, z1, y1 = corners[b]
                dy = y1 - y0
                if abs(dy) <= 1.0e-300:
                    continue
                tau = (level - y0) / dy
                if 0.0 <= tau <= 1.0:
                    xp = x0 + tau * (x1 - x0)
                    zp = z0 + tau * (z1 - z0)
                    if xmin <= xp <= xmax and zmin <= zp <= zmax:
                        edge_points.append((float(xp), float(zp)))
            if not edge_points:
                continue

            # A cell-crossing level segment is represented by its midpoint.
            xp = float(np.mean([point[0] for point in edge_points]))
            zp = float(np.mean([point[1] for point in edge_points]))
            rho_p = float(bilinear_sample(rho, x, z, np.array([xp]), np.array([zp]))[0])
            if rho_p <= rho_threshold:
                continue
            gx = float(bilinear_sample(y_der["fx"], x, z, np.array([xp]), np.array([zp]))[0])
            gz = float(bilinear_sample(y_der["fz"], x, z, np.array([xp]), np.array([zp]))[0])
            if gx * gx + gz * gz <= 1.0e-300:
                continue
            points.append((xp, zp))
    return points


def build_interface_rows(
    y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    positive_levels: list[float],
    half_width_y: float,
    samples: int,
    bbox: tuple[float, float, float, float] | None,
) -> list[dict[str, float]]:
    all_rows: list[dict[str, float]] = []
    signed_levels = [sign * level for level in positive_levels for sign in (-1.0, 1.0)]
    for level in signed_levels:
        points = find_interface_points_for_bbox(y, rho, x, z, level=level, bbox=bbox)
        rows = integrate_lines(y, x, z, points, half_width_y=half_width_y, samples=samples)
        for row in rows:
            row["level"] = level
        all_rows.extend(rows)
    return all_rows


def render_summary_figure(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    y: np.ndarray,
    derivative_norm: np.ndarray,
    mode_map: np.ndarray,
    interface_abs_source: np.ndarray,
    hotspot_mask: np.ndarray,
    interface_near_mask: np.ndarray,
    ell: float,
    t: float,
    positive_levels: list[float],
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 10.5), constrained_layout=True)

    im0 = axes[0, 0].pcolormesh(xg, zg, np.log10(np.maximum(rho, 1.0e-16)), shading="auto", cmap="viridis")
    axes[0, 0].contour(xg, zg, np.abs(y), levels=positive_levels, colors="white", linewidths=0.9)
    levels_label = ",".join(f"{level:g}" for level in positive_levels)
    axes[0, 0].set_title(r"$\log_{10}\rho$ with dynamic interfaces $|\ell^2\tilde R|=$" + levels_label)
    fig.colorbar(im0, ax=axes[0, 0])

    cmap = matplotlib.colors.ListedColormap(["#2f6f9f", "#d9a441", "#c7362f"])
    im1 = axes[0, 1].pcolormesh(xg, zg, mode_map, shading="auto", cmap=cmap, vmin=0, vmax=2)
    axes[0, 1].set_title("bulk-interface mode map: 0=EH-like, 1=saturated, 2=interface")
    fig.colorbar(im1, ax=axes[0, 1], ticks=[0, 1, 2])

    im2 = axes[1, 0].pcolormesh(
        xg,
        zg,
        np.log10(np.maximum(derivative_norm, 1.0e-30)),
        shading="auto",
        cmap="inferno",
    )
    axes[1, 0].contour(xg, zg, hotspot_mask.astype(float), levels=[0.5], colors="lime", linewidths=1.0)
    axes[1, 0].contour(xg, zg, interface_near_mask.astype(float), levels=[0.5], colors="cyan", linewidths=1.0)
    axes[1, 0].set_title("original derivative-term norm; lime=top 1%, cyan=near interface")
    fig.colorbar(im2, ax=axes[1, 0])

    im3 = axes[1, 1].pcolormesh(
        xg,
        zg,
        np.log10(np.maximum(interface_abs_source, 1.0e-30)),
        shading="auto",
        cmap="magma",
    )
    axes[1, 1].contour(xg, zg, hotspot_mask.astype(float), levels=[0.5], colors="lime", linewidths=1.0)
    axes[1, 1].set_title("interface absolute source proxy")
    fig.colorbar(im3, ax=axes[1, 1])

    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")

    fig.suptitle(f"D bulk-interface toy, ell={ell:g}, t={t:g}, |y| levels={levels_label}", fontsize=13)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolution", type=int, default=240)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--half-width-y", type=float, default=2.0)
    parser.add_argument("--samples", type=int, default=101)
    parser.add_argument("--near-radius", type=int, default=2)
    parser.add_argument(
        "--levels",
        type=str,
        default="1",
        help="Comma-separated positive |y| levels. Each level adds both y=+level and y=-level interfaces.",
    )
    parser.add_argument(
        "--bbox",
        type=str,
        default="global",
        help="Interface extraction box: 'global' or xmin,xmax,zmin,zmax. Earlier left-arc diagnostics used -7.4,-2.2,2.6,6.2.",
    )
    args = parser.parse_args()
    positive_levels = parse_positive_levels(args.levels)
    bbox = parse_bbox(args.bbox)

    params = FlatLocalizedCrossingParams(nx=args.resolution, nz=args.resolution)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    fields = make_fields(
        branch="D",
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

    rho = np.real_if_close(fields["rho"])
    r_tilde = np.real_if_close(fields["R_tilde"])
    derivative_norm = np.real_if_close(fields["derivative_norm"])
    residual_norm = np.real_if_close(fields["residual_norm"])
    y = np.real_if_close(branch_f_derivatives("D", r_tilde, args.ell)["y"])

    support = rho > 1.0e-3 * float(np.max(rho))
    eh_like = support & (np.abs(y) < 1.0)
    saturated = support & (np.abs(y) > 1.0)

    rows = build_interface_rows(
        y=y,
        rho=rho,
        x=x,
        z=z,
        positive_levels=positive_levels,
        half_width_y=args.half_width_y,
        samples=args.samples,
        bbox=bbox,
    )

    interface_signed_source = np.zeros_like(rho, dtype=float)
    interface_abs_source = np.zeros_like(rho, dtype=float)
    interface_mask = np.zeros_like(rho, dtype=bool)
    for row in rows:
        add_to_grid(interface_signed_source, x, z, row["x"], row["z"], row["signed_integral_d2phi"])
        add_to_grid(interface_abs_source, x, z, row["x"], row["z"], row["abs_integral_d2phi"])
        i = int(np.argmin(np.abs(x - row["x"])))
        j = int(np.argmin(np.abs(z - row["z"])))
        interface_mask[i, j] = True

    interface_near = dilate_mask(interface_mask, args.near_radius)
    mode_map = np.zeros_like(rho, dtype=float)
    mode_map[saturated] = 1.0
    mode_map[interface_near & support] = 2.0

    hotspot_thr = float(np.percentile(derivative_norm[support], 99.0)) if np.any(support) else float(np.max(derivative_norm))
    hotspot = support & (derivative_norm >= hotspot_thr)

    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    fig_path = (
        output
        / f"d_bulk_interface_toy_ell{args.ell:g}_n{args.resolution}_t{args.time:g}_levels{levels_slug(positive_levels)}.png"
    )
    render_summary_figure(
        out_path=fig_path,
        x=x,
        z=z,
        rho=rho,
        y=y,
        derivative_norm=derivative_norm,
        mode_map=mode_map,
        interface_abs_source=interface_abs_source,
        hotspot_mask=hotspot,
        interface_near_mask=interface_near & support,
        ell=args.ell,
        t=args.time,
        positive_levels=positive_levels,
    )

    signed_vals = np.array([row["signed_integral_d2phi"] for row in rows], dtype=float)
    abs_vals = np.array([row["abs_integral_d2phi"] for row in rows], dtype=float)
    rows_by_level = {}
    for level in sorted({float(row["level"]) for row in rows}):
        level_rows = [row for row in rows if float(row["level"]) == level]
        level_signed = np.array([row["signed_integral_d2phi"] for row in level_rows], dtype=float)
        level_abs = np.array([row["abs_integral_d2phi"] for row in level_rows], dtype=float)
        rows_by_level[f"{level:g}"] = {
            "line_count": len(level_rows),
            "signed_integral": scalar_stats(level_signed),
            "absolute_integral": scalar_stats(level_abs),
        }
    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "half_width_y": args.half_width_y,
            "samples": args.samples,
            "near_radius": args.near_radius,
            "positive_levels": positive_levels,
            "bbox": "global" if bbox is None else list(bbox),
        },
        "definitions": {
            "y": "y=ell^2 R_tilde",
            "dynamic_interfaces": "level sets y=+/-level for each positive level in params.positive_levels, inside params.bbox; regenerated from current fields each time slice",
            "EH_like_bulk": "support & |y| < 1",
            "saturated_bulk": "support & |y| > 1",
            "interface_near": "interface pixels dilated by near_radius cells",
            "hotspot": "top 1% of original derivative_norm within support",
            "interface_abs_source": "rasterized sum of line integral int |d_n^2 f_R| dn; proxy, not final continuum source",
        },
        "counts": {
            "support": int(np.sum(support)),
            "EH_like_bulk": int(np.sum(eh_like)),
            "saturated_bulk": int(np.sum(saturated)),
            "interface_lines": len(rows),
            "interface_pixels": int(np.sum(interface_mask)),
            "interface_near_support": int(np.sum(interface_near & support)),
            "derivative_hotspots": int(np.sum(hotspot)),
        },
        "fractions_support": {
            "EH_like_bulk": float(np.sum(eh_like) / max(np.sum(support), 1)),
            "saturated_bulk": float(np.sum(saturated) / max(np.sum(support), 1)),
            "interface_near": float(np.sum(interface_near & support) / max(np.sum(support), 1)),
            "hotspot_near_interface": float(np.sum(hotspot & interface_near) / max(np.sum(hotspot), 1)),
        },
        "line_integral_stats": {
            "signed_integral": scalar_stats(signed_vals),
            "absolute_integral": scalar_stats(abs_vals),
            "global_signed_to_abs_ratio": float(abs(np.sum(signed_vals)) / max(np.sum(abs_vals), 1.0e-300)) if rows else 0.0,
            "by_signed_level": rows_by_level,
        },
        "field_stats_support": {
            "derivative_norm": scalar_stats(derivative_norm[support]),
            "residual_norm": scalar_stats(residual_norm[support]),
            "interface_abs_source_nonzero": scalar_stats(interface_abs_source[interface_abs_source > 0.0]),
        },
        "toy_policy": {
            "EH_like_bulk": "evolve with ordinary D/scalar-tensor bulk equations or their weak-curvature approximation",
            "saturated_bulk": "do not evolve as ordinary Cauchy PDE in the first toy model; use fixed/reference extension or minimum-curvature extension",
            "interface": "apply jump/matching source calibrated by line integrals across y=+/-1",
        },
        "figure": str(fig_path.resolve()),
    }

    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
