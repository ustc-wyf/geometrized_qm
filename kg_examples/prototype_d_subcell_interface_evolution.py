from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from integrate_d_interface_jumps import bilinear_sample, integrate_lines, scalar_stats, spatial_derivatives
from prototype_d_reference_imex import (
    array_stats,
    imex_local_seed,
    reference_snapshot_data,
    trace_residual_2p1_from_chi,
)
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


def cloud_in_cell_add(
    accum: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    xp: float,
    zp: float,
    value: float,
) -> None:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    u = (xp - float(x[0])) / dx
    v = (zp - float(z[0])) / dz
    i0 = int(np.floor(u))
    j0 = int(np.floor(v))
    if i0 < 0 or i0 >= len(x) - 1 or j0 < 0 or j0 >= len(z) - 1:
        return
    tx = float(np.clip(u - i0, 0.0, 1.0))
    tz = float(np.clip(v - j0, 0.0, 1.0))
    weights = (
        ((1.0 - tx) * (1.0 - tz), i0, j0),
        (tx * (1.0 - tz), i0 + 1, j0),
        ((1.0 - tx) * tz, i0, j0 + 1),
        (tx * tz, i0 + 1, j0 + 1),
    )
    for weight, i, j in weights:
        accum[i, j] += weight * value


def extract_interface_segments(
    y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
) -> list[dict[str, float]]:
    rho_threshold = 1.0e-3 * float(np.max(rho))
    y_der = spatial_derivatives(y, x, z)
    segments: list[dict[str, float]] = []
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

            # For ambiguous marching-square cells this picks the longest chord,
            # which is adequate for this first subcell source prototype.
            p0 = crossings[0]
            p1 = crossings[1]
            if len(crossings) > 2:
                best_len = -1.0
                for a in range(len(crossings)):
                    for b in range(a + 1, len(crossings)):
                        dist = float(np.hypot(crossings[a][0] - crossings[b][0], crossings[a][1] - crossings[b][1]))
                        if dist > best_len:
                            best_len = dist
                            p0 = crossings[a]
                            p1 = crossings[b]
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
            segments.append({"x": xp, "z": zp, "length": length, "level": level, "grad_y": grad})
    return segments


def solve_screened_periodic(source: np.ndarray, x: np.ndarray, z: np.ndarray, screen_length: float) -> np.ndarray:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    kx = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(len(z), d=dz)
    kx2, kz2 = np.meshgrid(kx * kx, kz * kz, indexing="ij")
    denom = 1.0 + screen_length * screen_length * (kx2 + kz2)
    return real_array(np.fft.ifft2(np.fft.fft2(source) / denom))


def build_subcell_line_source(
    ref: dict[str, np.ndarray],
    x: np.ndarray,
    z: np.ndarray,
    positive_levels: list[float],
    half_width_y: float,
    samples: int,
    source_mode: str,
    ell: float,
) -> dict[str, object]:
    y = real_array(ref.get("raw_y_ref", ref["y_ref"]))
    rho = real_array(ref["rho"])
    trace_residual = real_array(ref["trace_residual_2p1"])
    seed = imex_local_seed(
        chi_ref=real_array(ref["chi_ref"]),
        r_tilde=real_array(ref["R_tilde"]),
        trace_residual=trace_residual,
        derivative_norm=real_array(ref["derivative_norm"]),
        ell=ell,
        dx=float(ref["dx"]),
        dz=float(ref["dz"]),
    )
    delta_seed = real_array(seed["delta_chi_seed"])
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    cell_area = dx * dz

    signed_levels = [sign * level for level in positive_levels for sign in (-1.0, 1.0)]
    segments: list[dict[str, float]] = []
    for level in signed_levels:
        segments.extend(extract_interface_segments(y, rho, x, z, level))

    points = [(seg["x"], seg["z"]) for seg in segments]
    line_rows = integrate_lines(y, x, z, points, half_width_y=half_width_y, samples=samples)
    for seg, row in zip(segments, line_rows):
        seg.update(row)

    source = np.zeros_like(y, dtype=float)
    line_strengths: list[float] = []
    for seg in segments:
        residual_here = float(bilinear_sample(trace_residual, x, z, np.array([seg["x"]]), np.array([seg["z"]]))[0])
        if source_mode == "trace":
            strength = -residual_here
        elif source_mode == "seed":
            strength = float(bilinear_sample(delta_seed, x, z, np.array([seg["x"]]), np.array([seg["z"]]))[0])
        elif source_mode == "signed_jump":
            strength = -float(seg.get("signed_integral_d2phi", 0.0))
        elif source_mode == "hybrid":
            jump_sign = np.sign(float(seg.get("signed_integral_d2phi", 0.0)))
            strength = -abs(residual_here) * (jump_sign if jump_sign != 0.0 else np.sign(residual_here))
        else:
            raise ValueError(f"unknown source_mode: {source_mode}")
        line_strengths.append(strength)
        cloud_in_cell_add(source, x, z, seg["x"], seg["z"], strength * seg["length"] / cell_area)

    return {
        "segments": segments,
        "source": source,
        "line_strengths": np.array(line_strengths, dtype=float),
    }


def evaluate_trace_residual(
    ref: dict[str, np.ndarray],
    chi_0: np.ndarray,
    ell: float,
    probe_dt: float,
) -> np.ndarray:
    return real_array(
        trace_residual_2p1_from_chi(
            chi_m=ref["chi_m"],
            chi_0=chi_0,
            chi_p=ref["chi_p"],
            metric_m=ref["metric_cov_m"],
            metric_0=ref["metric_cov"],
            metric_p=ref["metric_cov_p"],
            stress_trace=ref["stress_trace"],
            ell=ell,
            probe_dt=probe_dt,
            dx=float(ref["dx"]),
            dz=float(ref["dz"]),
        )["trace_residual_2p1"]
    )


def choose_line_search_update(
    ref: dict[str, np.ndarray],
    delta_unit: np.ndarray,
    support: np.ndarray,
    ell: float,
    probe_dt: float,
    amplitudes: list[float],
) -> dict[str, object]:
    before = real_array(ref["trace_residual_2p1"])
    results: list[dict[str, object]] = []
    best: dict[str, object] | None = None
    for amp in amplitudes:
        chi_trial = real_array(ref["chi_ref"]) + amp * delta_unit
        residual = evaluate_trace_residual(ref, chi_trial, ell=ell, probe_dt=probe_dt)
        stats = array_stats(residual, support)
        score = stats["p95"] + 1.0e-5 * stats["abs_max"]
        row = {"amplitude": float(amp), "score": float(score), "support_trace": stats}
        results.append(row)
        if best is None or score < float(best["score"]):
            best = {"amplitude": float(amp), "score": float(score), "residual": residual, "chi": chi_trial, "stats": stats}
    assert best is not None
    best["line_search"] = results
    best["before_stats"] = array_stats(before, support)
    return best


def render_figure(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    ref: dict[str, np.ndarray],
    source: np.ndarray,
    delta: np.ndarray,
    residual_after: np.ndarray,
    support: np.ndarray,
    ell: float,
    t: float,
    levels: list[float],
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    source = real_array(source)
    delta = real_array(delta)
    residual_after = real_array(residual_after)
    y_abs = np.abs(real_array(ref.get("raw_y_ref", ref["y_ref"])))
    residual_before = real_array(ref["trace_residual_2p1"])
    change = np.abs(residual_after) - np.abs(residual_before)

    panels = [
        (np.log10(np.maximum(real_array(ref["rho"]), 1.0e-16)), r"$\log_{10}\rho$", "viridis", False),
        (source, r"subcell line source $S_\Gamma$", "coolwarm", True),
        (delta, r"chosen $\delta\chi$ from screened interface solve", "coolwarm", True),
        (residual_before, "trace residual before T_ref", "coolwarm", True),
        (residual_after, "trace residual after T_toy", "coolwarm", True),
        (change, "|T_toy|-|T_ref|", "coolwarm", True),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.8), constrained_layout=True)
    for ax, (field, title, cmap, symmetric) in zip(axes.ravel(), panels):
        if symmetric:
            vmax = float(np.percentile(np.abs(field[np.isfinite(field)]), 99.0)) if np.any(np.isfinite(field)) else 1.0
            vmax = max(vmax, 1.0e-12)
            im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
        else:
            im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap)
        ax.contour(xg, zg, y_abs, levels=levels, colors="white", linewidths=0.7)
        ax.contour(xg, zg, support.astype(float), levels=[0.5], colors="black", linewidths=0.5, linestyles="--")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, shrink=0.86)
    fig.suptitle(
        "D subcell interface toy correction: white contours are |y| levels; dashed black is support "
        r"($\rho>10^{-3}\rho_{\max}$).",
        fontsize=12,
    )
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
    rho = real_array(ref["rho"])
    support = rho > 1.0e-3 * float(np.max(rho))
    levels = parse_positive_levels(args.levels)
    line_source = build_subcell_line_source(
        ref=ref,
        x=x,
        z=z,
        positive_levels=levels,
        half_width_y=args.half_width_y,
        samples=args.samples,
        source_mode=args.source_mode,
        ell=args.ell,
    )
    source = real_array(line_source["source"])
    delta_raw = solve_screened_periodic(source, x=x, z=z, screen_length=args.screen_length)

    raw_p95 = float(np.percentile(np.abs(delta_raw[support]), 95.0)) if np.any(support) else 0.0
    if raw_p95 > 0.0:
        delta_unit = delta_raw * (args.target_delta_p95 / raw_p95)
    else:
        delta_unit = delta_raw

    amplitudes = [float(v) for v in args.amplitudes.split(",") if v.strip()]
    chosen = choose_line_search_update(
        ref=ref,
        delta_unit=delta_unit,
        support=support,
        ell=args.ell,
        probe_dt=args.probe_dt,
        amplitudes=amplitudes,
    )
    amplitude = float(chosen["amplitude"])
    delta = amplitude * delta_unit
    residual_after = real_array(chosen["residual"])

    fig_path = out / f"d_subcell_interface_toy_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_figure(
        out_path=fig_path,
        x=x,
        z=z,
        ref=ref,
        source=source,
        delta=delta,
        residual_after=residual_after,
        support=support,
        ell=args.ell,
        t=args.time,
        levels=levels,
    )

    source_nonzero = source[np.abs(source) > 0.0]
    segments = line_source["segments"]
    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "levels_abs_y": levels,
            "source_mode": args.source_mode,
            "screen_length": args.screen_length,
            "target_delta_p95": args.target_delta_p95,
            "chosen_amplitude": amplitude,
        },
        "definitions": {
            "y": "y=ell^2 R_tilde. For D branch, f_R=sech^2(y).",
            "subcell_interface": "line segments from the dynamic level sets y=+/-level; not a raster thick band",
            "line_source": "cloud-in-cell deposition of interface segment strengths times segment length divided by cell area",
            "toy_equation": "(1 - screen_length^2 Delta) delta_chi = line_source, solved spectrally with periodic boundary as a first prototype",
            "support": "rho > 1e-3 * rho_max",
            "trace_residual": "2+1d scalar-tensor D-branch trace residual evaluated at fixed reference metric and stress trace",
        },
        "counts": {
            "support": int(np.sum(support)),
            "segments": len(segments),
            "source_nonzero_pixels": int(source_nonzero.size),
        },
        "line_source_stats": {
            "source_nonzero": scalar_stats(source_nonzero),
            "segment_strength": scalar_stats(line_source["line_strengths"]),
            "segment_length": scalar_stats([seg["length"] for seg in segments]),
            "signed_integral_d2phi": scalar_stats([seg.get("signed_integral_d2phi", 0.0) for seg in segments]),
            "abs_integral_d2phi": scalar_stats([seg.get("abs_integral_d2phi", 0.0) for seg in segments]),
        },
        "correction_stats": {
            "delta_raw_support_p95": raw_p95,
            "delta_unit_support": array_stats(delta_unit, support),
            "delta_chosen_support": array_stats(delta, support),
        },
        "residual_stats": {
            "before_support": array_stats(real_array(ref["trace_residual_2p1"]), support),
            "after_support": array_stats(residual_after, support),
            "abs_change_support": array_stats(np.abs(residual_after) - np.abs(real_array(ref["trace_residual_2p1"])), support),
            "line_search": chosen["line_search"],
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
    parser.add_argument("--source-mode", choices=["trace", "seed", "signed_jump", "hybrid"], default="seed")
    parser.add_argument("--half-width-y", type=float, default=2.0)
    parser.add_argument("--samples", type=int, default=51)
    parser.add_argument("--screen-length", type=float, default=0.75)
    parser.add_argument("--target-delta-p95", type=float, default=0.02)
    parser.add_argument("--amplitudes", type=str, default="-2,-1,-0.5,-0.25,0,0.25,0.5,1,2")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
