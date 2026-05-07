from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from integrate_d_interface_jumps import bilinear_sample, scalar_stats, spatial_derivatives
from prototype_d_reference_imex import reference_snapshot_data
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def real_array(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values)
    if np.iscomplexobj(arr):
        return np.asarray(np.real(arr), dtype=float)
    return np.asarray(arr, dtype=float)


def parse_positive_levels(raw: str) -> list[float]:
    levels = sorted({abs(float(part.strip())) for part in raw.split(",") if part.strip()})
    if not levels or levels[0] <= 0.0:
        raise ValueError("--levels must contain positive numbers, e.g. 1 or 0.5,1,2")
    return levels


def find_level_segments(
    raw_y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
) -> list[dict[str, float]]:
    rho_threshold = 1.0e-3 * float(np.max(rho))
    y_der = spatial_derivatives(raw_y, x, z)
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
            length = float(np.hypot(p1[0] - p0[0], p1[1] - p0[1]))
            if length <= 0.0:
                continue
            xp = 0.5 * (p0[0] + p1[0])
            zp = 0.5 * (p0[1] + p1[1])
            rho_p = float(bilinear_sample(rho, x, z, np.array([xp]), np.array([zp]))[0])
            if rho_p <= rho_threshold:
                continue
            gx = float(bilinear_sample(y_der["fx"], x, z, np.array([xp]), np.array([zp]))[0])
            gz = float(bilinear_sample(y_der["fz"], x, z, np.array([xp]), np.array([zp]))[0])
            grad = float(np.hypot(gx, gz))
            if grad <= 1.0e-300:
                continue
            segments.append(
                {
                    "x0": float(p0[0]),
                    "z0": float(p0[1]),
                    "x1": float(p1[0]),
                    "z1": float(p1[1]),
                    "x": xp,
                    "z": zp,
                    "length": length,
                    "level": level,
                    "grad_raw_y": grad,
                    "nx": gx / grad,
                    "nz": gz / grad,
                }
            )
    return segments


def line_integrals_for_segment(
    seg: dict[str, float],
    raw_y: np.ndarray,
    chi: np.ndarray,
    phi: np.ndarray,
    f: np.ndarray,
    stress_trace: np.ndarray,
    trace_residual: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    half_width_y: float,
    samples: int,
) -> dict[str, float]:
    grad = max(float(seg["grad_raw_y"]), 1.0e-300)
    half_width = half_width_y / grad
    s = np.linspace(-half_width, half_width, samples)
    xs = float(seg["x"]) + s * float(seg["nx"])
    zs = float(seg["z"]) + s * float(seg["nz"])
    valid = (xs >= float(x[0])) & (xs <= float(x[-1])) & (zs >= float(z[0])) & (zs <= float(z[-1]))
    if int(np.sum(valid)) < max(7, samples // 2):
        return {"valid": 0.0}
    xs = xs[valid]
    zs = zs[valid]
    s = s[valid]

    raw_y_vals = bilinear_sample(raw_y, x, z, xs, zs)
    chi_vals = bilinear_sample(chi, x, z, xs, zs)
    phi_vals = bilinear_sample(phi, x, z, xs, zs)
    f_vals = bilinear_sample(f, x, z, xs, zs)
    stress_vals = bilinear_sample(stress_trace, x, z, xs, zs)
    residual_vals = bilinear_sample(trace_residual, x, z, xs, zs)

    algebraic_vals = phi_vals * chi_vals - 1.5 * f_vals - stress_vals
    dphi_dn = np.gradient(phi_vals, s, edge_order=2)
    d2phi_dnn = np.gradient(dphi_dn, s, edge_order=2)

    normal_jump = float(2.0 * (dphi_dn[-1] - dphi_dn[0]))
    normal_abs = float(np.trapezoid(np.abs(2.0 * d2phi_dnn), s))
    bulk_integral = float(np.trapezoid(algebraic_vals, s))
    bulk_abs = float(np.trapezoid(np.abs(algebraic_vals), s))
    trace_integral = float(np.trapezoid(residual_vals, s))
    trace_abs = float(np.trapezoid(np.abs(residual_vals), s))
    weak_residual = normal_jump + bulk_integral

    return {
        "valid": 1.0,
        "x": float(seg["x"]),
        "z": float(seg["z"]),
        "level": float(seg["level"]),
        "segment_length": float(seg["length"]),
        "grad_raw_y": float(seg["grad_raw_y"]),
        "half_width": float(half_width),
        "raw_y_min": float(np.min(raw_y_vals)),
        "raw_y_max": float(np.max(raw_y_vals)),
        "normal_jump": normal_jump,
        "normal_abs_integral": normal_abs,
        "bulk_integral": bulk_integral,
        "bulk_abs_integral": bulk_abs,
        "weak_residual": weak_residual,
        "trace_integral_direct": trace_integral,
        "trace_abs_integral_direct": trace_abs,
        "direct_minus_weak": trace_integral - weak_residual,
        "normal_to_bulk_abs_ratio": float(abs(normal_jump) / max(abs(bulk_integral), 1.0e-300)),
        "weak_to_abs_ratio": float(abs(weak_residual) / max(normal_abs + bulk_abs, 1.0e-300)),
    }


def summarize_rows(rows: list[dict[str, float]]) -> dict[str, object]:
    valid = [row for row in rows if row.get("valid", 0.0) > 0.5]
    return {
        "count": len(valid),
        "normal_jump": scalar_stats([row["normal_jump"] for row in valid]),
        "normal_abs_integral": scalar_stats([row["normal_abs_integral"] for row in valid]),
        "bulk_integral": scalar_stats([row["bulk_integral"] for row in valid]),
        "bulk_abs_integral": scalar_stats([row["bulk_abs_integral"] for row in valid]),
        "weak_residual": scalar_stats([row["weak_residual"] for row in valid]),
        "trace_integral_direct": scalar_stats([row["trace_integral_direct"] for row in valid]),
        "direct_minus_weak": scalar_stats([row["direct_minus_weak"] for row in valid]),
        "normal_to_bulk_abs_ratio": scalar_stats([row["normal_to_bulk_abs_ratio"] for row in valid]),
        "weak_to_abs_ratio": scalar_stats([row["weak_to_abs_ratio"] for row in valid]),
        "global_weak_to_abs_ratio": float(
            abs(sum(row["weak_residual"] for row in valid))
            / max(sum(row["normal_abs_integral"] + row["bulk_abs_integral"] for row in valid), 1.0e-300)
        )
        if valid
        else 0.0,
    }


def render_plot(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    trace_residual: np.ndarray,
    rows: list[dict[str, float]],
    levels: list[float],
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > 1.0e-3 * float(np.max(rho))
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6), constrained_layout=True)

    im0 = axes[0].pcolormesh(xg, zg, np.log10(np.maximum(rho, 1.0e-16)), shading="auto", cmap="viridis")
    axes[0].contour(xg, zg, np.abs(raw_y), levels=levels, colors="white", linewidths=0.8)
    axes[0].set_title("log10 rho; white: |raw_y| levels")
    fig.colorbar(im0, ax=axes[0])

    vmax = float(np.percentile(np.abs(trace_residual[support]), 99.0)) if np.any(support) else 1.0
    vmax = max(vmax, 1.0e-12)
    im1 = axes[1].pcolormesh(xg, zg, trace_residual, shading="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    axes[1].contour(xg, zg, np.abs(raw_y), levels=levels, colors="black", linewidths=0.6)
    axes[1].set_title("D trace residual on reference slice")
    fig.colorbar(im1, ax=axes[1])

    valid = [row for row in rows if row.get("valid", 0.0) > 0.5]
    if valid:
        xs = [row["x"] for row in valid]
        zs = [row["z"] for row in valid]
        vals = [row["weak_residual"] for row in valid]
        vmax2 = max(float(np.percentile(np.abs(vals), 95.0)), 1.0e-12)
        sc = axes[2].scatter(xs, zs, c=vals, s=5.0, cmap="coolwarm", vmin=-vmax2, vmax=vmax2)
        fig.colorbar(sc, ax=axes[2])
    axes[2].contour(xg, zg, support.astype(float), levels=[0.5], colors="gray", linewidths=0.5)
    axes[2].set_title("weak residual on subcell interface samples")

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

    chi = real_array(ref["chi_ref"])
    raw_y = args.ell * args.ell * chi
    tanh_y = np.tanh(raw_y)
    phi = 1.0 - tanh_y * tanh_y
    f = tanh_y / (args.ell * args.ell)
    rho = real_array(ref["rho"])
    stress_trace = real_array(ref["stress_trace"])
    trace_residual = real_array(ref["trace_residual_2p1"])

    levels = parse_positive_levels(args.levels)
    signed_levels = [sign * level for level in levels for sign in (-1.0, 1.0)]
    rows: list[dict[str, float]] = []
    for level in signed_levels:
        segments = find_level_segments(raw_y, rho, x, z, level)
        for seg in segments:
            rows.append(
                line_integrals_for_segment(
                    seg,
                    raw_y=raw_y,
                    chi=chi,
                    phi=phi,
                    f=f,
                    stress_trace=stress_trace,
                    trace_residual=trace_residual,
                    x=x,
                    z=z,
                    half_width_y=args.half_width_y,
                    samples=args.samples,
                )
            )

    fig_path = out / f"d_trace_interface_weak_form_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_plot(fig_path, x=x, z=z, rho=rho, raw_y=raw_y, trace_residual=trace_residual, rows=rows, levels=levels)

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
            "raw_y": "raw_y = ell^2 * R_tilde. This is the interface variable.",
            "phi": "phi = f_R = sech^2(raw_y)",
            "trace_equation": "phi*R_tilde - 3/2*f + 2*Box(phi) - T/Mp^2 = 0 in 2+1d",
            "normal_jump": "2 * [partial_n phi] across a short spatial normal line; proxy for the singular normal part of 2*Box(phi)",
            "bulk_integral": "integral of phi*R_tilde - 3/2*f - T/Mp^2 along the same normal line",
            "weak_residual": "normal_jump + bulk_integral. In the ideal thin-layer weak trace equation this should vanish up to tangential/time/connection contributions.",
            "direct_minus_weak": "direct integral of the full trace residual minus weak_residual; estimates omitted time/tangential/connection pieces and interpolation error",
        },
        "counts": {
            "rows": len(rows),
            "valid_rows": int(sum(1 for row in rows if row.get("valid", 0.0) > 0.5)),
        },
        "summary": summarize_rows(rows),
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
    parser.add_argument("--half-width-y", type=float, default=2.0)
    parser.add_argument("--samples", type=int, default=81)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
