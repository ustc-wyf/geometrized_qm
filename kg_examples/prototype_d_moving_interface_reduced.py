from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from diagnose_d_interface_speed_law import (
    solve_ft_for_target_norm,
    target_normal_norm_candidates,
)
from diagnose_d_spacetime_interface_jump import d_branch_f, d_branch_phi, d_phi_d_raw_y
from diagnose_d_trace_interface_weak_form import find_level_segments, parse_positive_levels, real_array
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from prototype_d_reference_imex import reference_snapshot_data
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def cloud_in_cell_add(
    accum: np.ndarray,
    weights: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    xp: float,
    zp: float,
    value: float,
    weight_scale: float,
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
    for weight, i, j in (
        ((1.0 - tx) * (1.0 - tz), i0, j0),
        (tx * (1.0 - tz), i0 + 1, j0),
        ((1.0 - tx) * tz, i0, j0 + 1),
        (tx * tz, i0 + 1, j0 + 1),
    ):
        scaled = weight * weight_scale
        accum[i, j] += scaled * value
        weights[i, j] += scaled


def neighbor_average_extend(values: np.ndarray, known: np.ndarray, iterations: int) -> np.ndarray:
    out = values.copy()
    mask = known.copy()
    for _ in range(iterations):
        neighbor_sum = np.zeros_like(out)
        neighbor_count = np.zeros_like(out)
        for shifted_values, shifted_mask in (
            (np.roll(out, 1, axis=0), np.roll(mask, 1, axis=0)),
            (np.roll(out, -1, axis=0), np.roll(mask, -1, axis=0)),
            (np.roll(out, 1, axis=1), np.roll(mask, 1, axis=1)),
            (np.roll(out, -1, axis=1), np.roll(mask, -1, axis=1)),
        ):
            neighbor_sum += np.where(shifted_mask, shifted_values, 0.0)
            neighbor_count += shifted_mask.astype(float)
        fill = (~mask) & (neighbor_count > 0.0)
        if not np.any(fill):
            break
        out[fill] = neighbor_sum[fill] / neighbor_count[fill]
        mask[fill] = True
    return out


def compute_interface_rows(
    raw_y: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    stress_trace: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    ell: float,
    half_width_y: float,
    samples: int,
    ft_guess_grid: np.ndarray,
) -> list[dict[str, float]]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)
    q = np.linspace(level - half_width_y, level + half_width_y, samples)
    phi_prime = d_phi_d_raw_y(q)
    delta_phi_prime = float(phi_prime[-1] - phi_prime[0])
    phi_q = d_branch_phi(q)
    f_q = d_branch_f(q, ell)
    chi_q = q / (ell * ell)

    rows: list[dict[str, float]] = []
    for seg in find_level_segments(raw_y, rho, x, z, level):
        xp = np.array([seg["x"]])
        zp = np.array([seg["z"]])
        fx = float(bilinear_sample(raw_y_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(raw_y_z, x, z, xp, zp)[0])
        grad_spatial = float(np.hypot(fx, fz))
        if grad_spatial <= 1.0e-300:
            continue

        ft_guess = float(bilinear_sample(ft_guess_grid, x, z, xp, zp)[0])
        ginv_here = np.zeros((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                ginv_here[a, b] = float(bilinear_sample(metric_inv[..., a, b], x, z, xp, zp)[0])
        covector_guess = np.array([ft_guess, fx, fz], dtype=float)
        normal_norm_guess = float(covector_guess @ ginv_here @ covector_guess)

        stress_here = float(bilinear_sample(stress_trace, x, z, xp, zp)[0])
        algebraic_q = phi_q * chi_q - 1.5 * f_q - stress_here
        algebraic_integral = float(np.trapezoid(algebraic_q, q))
        candidates = target_normal_norm_candidates(delta_phi_prime, algebraic_integral)

        best: dict[str, float] | None = None
        for cand in candidates:
            ft_corr, disc = solve_ft_for_target_norm(
                ginv=ginv_here,
                fx=fx,
                fz=fz,
                target_norm=cand["normal_norm"],
                ft_reference=ft_guess,
            )
            if not np.isfinite(ft_corr):
                continue
            speed_corr = float(-ft_corr / grad_spatial)
            row = {
                **cand,
                "ft_corrected": ft_corr,
                "quadratic_discriminant": disc,
                "coordinate_normal_speed_corrected": speed_corr,
                "abs_ft_change": abs(ft_corr - ft_guess),
            }
            if best is None or row["abs_ft_change"] < best["abs_ft_change"]:
                best = row

        base = {
            "x0": float(seg["x0"]),
            "z0": float(seg["z0"]),
            "x1": float(seg["x1"]),
            "z1": float(seg["z1"]),
            "x": float(seg["x"]),
            "z": float(seg["z"]),
            "level": float(level),
            "segment_length": float(seg["length"]),
            "spatial_grad_norm": grad_spatial,
            "ft_guess": ft_guess,
            "normal_norm_guess": normal_norm_guess,
            "algebraic_integral_raw_y": algebraic_integral,
        }
        if best is None:
            rows.append({**base, "solved": 0.0})
            continue
        rows.append(
            {
                **base,
                "solved": 1.0,
                "ft_corrected": best["ft_corrected"],
                "coordinate_normal_speed_corrected": best["coordinate_normal_speed_corrected"],
                "target_signature": best["sigma"],
                "normal_norm_corrected": best["normal_norm"],
            }
        )
    return rows


def rows_to_ft_grid(
    rows: list[dict[str, float]],
    x: np.ndarray,
    z: np.ndarray,
    extension_iterations: int,
) -> tuple[np.ndarray, np.ndarray]:
    accum = np.zeros((len(x), len(z)), dtype=float)
    weights = np.zeros_like(accum)
    for row in rows:
        if row.get("solved", 0.0) <= 0.5:
            continue
        cloud_in_cell_add(
            accum,
            weights,
            x,
            z,
            row["x"],
            row["z"],
            row["ft_corrected"],
            max(row["segment_length"], 1.0e-12),
        )
    known = weights > 0.0
    ft_grid = np.zeros_like(accum)
    ft_grid[known] = accum[known] / weights[known]
    if extension_iterations > 0:
        ft_grid = neighbor_average_extend(ft_grid, known, iterations=extension_iterations)
    return ft_grid, known


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_line_collection(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    values: list[float] | np.ndarray | None = None,
    colors: str | list[str] = "white",
    cmap: str = "coolwarm",
    vlim: float | None = None,
    linewidth: float = 0.9,
    linestyle: str = "solid",
) -> LineCollection | None:
    if not rows:
        return None
    if values is None:
        lc = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth, linestyles=linestyle)
    else:
        vals = np.asarray(values, dtype=float)
        lc = LineCollection(row_segments(rows), cmap=cmap, linewidths=linewidth, linestyles=linestyle)
        lc.set_array(vals)
        if vlim is not None:
            lc.set_clim(-vlim, vlim)
    ax.add_collection(lc)
    return lc


def render_frame(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    rows: list[dict[str, float]],
    step: int,
    dt: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    unsolved = [row for row in rows if row.get("solved", 0.0) <= 0.5]
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.7), constrained_layout=True)

    im0 = axes[0].pcolormesh(xg, zg, np.log10(np.maximum(rho, 1.0e-16)), shading="auto", cmap="viridis")
    colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in rows]
    add_line_collection(axes[0], rows, colors=colors, linewidth=0.85)
    axes[0].set_title("rho background; orange +1 / green -1")
    fig.colorbar(im0, ax=axes[0])

    vmax_y = max(float(np.percentile(np.abs(raw_y[np.isfinite(raw_y)]), 95.0)), 1.0)
    im1 = axes[1].pcolormesh(xg, zg, raw_y, shading="auto", cmap="coolwarm", vmin=-vmax_y, vmax=vmax_y)
    add_line_collection(axes[1], solved, colors="#111111", linewidth=0.85)
    add_line_collection(axes[1], unsolved, colors="#777777", linewidth=0.65, linestyle="dotted")
    axes[1].set_title("raw_y; solid solved, dotted unresolved")
    fig.colorbar(im1, ax=axes[1])

    speeds = [row["coordinate_normal_speed_corrected"] for row in solved]
    vlim = max(float(np.percentile(np.abs(speeds), 95.0)), 1.0e-12) if speeds else 1.0
    lc = add_line_collection(axes[2], solved, values=speeds, vlim=vlim, linewidth=1.05)
    add_line_collection(axes[2], unsolved, colors="0.3", linewidth=0.65, linestyle="dotted")
    if lc is not None:
        fig.colorbar(lc, ax=axes[2])
    axes[2].set_title("jump-law interface speed")

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
        ax.set_aspect("equal")
    fig.suptitle(f"D moving-interface reduced prototype, step={step}, pseudo time={step * dt:.3e}", fontsize=12)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def summarize_step(step: int, rows: list[dict[str, float]]) -> dict[str, object]:
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    unsolved = len(rows) - len(solved)
    return {
        "step": step,
        "segment_count": len(rows),
        "solved_count": len(solved),
        "unresolved_count": unsolved,
        "solved_fraction": float(len(solved) / max(len(rows), 1)),
        "total_line_length": float(sum(row["segment_length"] for row in rows)),
        "solved_line_length": float(sum(row["segment_length"] for row in solved)),
        "speed_corrected": scalar_stats([row["coordinate_normal_speed_corrected"] for row in solved]),
        "normal_norm_corrected": scalar_stats([row["normal_norm_corrected"] for row in solved]),
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    frame_dir = out / "frames"
    frame_dir.mkdir(parents=True, exist_ok=True)

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
    rho = real_array(ref["rho"])
    metric_inv = real_array(ref["metric_inv"])
    stress_trace = real_array(ref["stress_trace"])
    levels = parse_positive_levels(args.levels)
    ft_guess_grid = np.zeros_like(raw_y)

    step_summaries: list[dict[str, object]] = []
    frame_paths: list[str] = []
    for step in range(args.steps + 1):
        rows: list[dict[str, float]] = []
        for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
            rows.extend(
                compute_interface_rows(
                    raw_y=raw_y,
                    metric_inv=metric_inv,
                    rho=rho,
                    stress_trace=stress_trace,
                    x=x,
                    z=z,
                    level=level,
                    ell=args.ell,
                    half_width_y=args.half_width_y,
                    samples=args.samples,
                    ft_guess_grid=ft_guess_grid,
                )
            )
        step_summaries.append(summarize_step(step, rows))
        if step % args.render_every == 0 or step == args.steps:
            frame_path = frame_dir / f"moving_interface_step_{step:03d}.png"
            render_frame(frame_path, x=x, z=z, rho=rho, raw_y=raw_y, rows=rows, step=step, dt=args.dt)
            frame_paths.append(str(frame_path.resolve()))
        if step == args.steps:
            break

        ft_grid, known = rows_to_ft_grid(rows, x=x, z=z, extension_iterations=args.extension_iterations)
        raw_y = raw_y + args.dt * ft_grid
        ft_guess_grid = ft_grid

    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "reference_time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "dt": args.dt,
            "steps": args.steps,
            "levels_abs_raw_y": levels,
            "half_width_y": args.half_width_y,
            "samples": args.samples,
            "extension_iterations": args.extension_iterations,
        },
        "definitions": {
            "scope": "Reduced moving-interface prototype. It evolves raw_y by the leading trace/scalaron jump-law interface F_t with frozen metric/stress/rho coefficients; it is not the final full tensor-matching dynamics.",
            "topology_handling": "At every step interfaces are re-extracted from the full raw_y field, so existing curves may split, merge, appear, or disappear if the evolved level set crosses +/-level.",
            "unresolved_segments": "Segments where the trace-only jump law has no real F_t root. These are drawn as dotted lines and are not clipped or damped into fake solutions.",
            "velocity_extension": "Numerical level-set extension of solved interface F_t from line segments to nearby grid points by cloud-in-cell deposition and neighbor averaging; this is an algorithmic extension, not an extra physical source.",
        },
        "steps": step_summaries,
        "files": {"frames": frame_paths},
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=61)
    parser.add_argument("--dt", type=float, default=1.0e-4)
    parser.add_argument("--steps", type=int, default=6)
    parser.add_argument("--extension-iterations", type=int, default=8)
    parser.add_argument("--render-every", type=int, default=1)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
