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
from prototype_d_moving_interface_reduced import cloud_in_cell_add, neighbor_average_extend
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
    sign_source: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    segments: list[dict[str, float]],
) -> np.ndarray:
    """Return sign(sign_source) times distance to the nearest segment."""
    if not segments:
        return sign_source.copy()
    xg, zg = np.meshgrid(x, z, indexing="ij")
    dist2 = np.full_like(sign_source, np.inf, dtype=float)
    for seg in segments:
        dist2 = np.minimum(
            dist2,
            point_segment_distance_squared(xg, zg, seg["x0"], seg["z0"], seg["x1"], seg["z1"]),
        )
    sign = np.where(sign_source >= 0.0, 1.0, -1.0)
    return sign * np.sqrt(dist2)


def initialize_branch_distance(
    raw_y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
) -> tuple[np.ndarray, list[dict[str, float]]]:
    segments = find_level_segments(raw_y, rho, x, z, level)
    return signed_distance_to_segments(raw_y - level, x, z, segments), segments


def reinitialize_distance(
    distance_field: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, list[dict[str, float]], bool]:
    segments = find_level_segments(distance_field, rho, x, z, 0.0)
    if not segments:
        return distance_field.copy(), segments, False
    return signed_distance_to_segments(distance_field, x, z, segments), segments, True


def compute_distance_interface_rows(
    distance_field: np.ndarray,
    raw_y_ref_grad_norm: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    stress_trace: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    branch_level: float,
    ell: float,
    half_width_y: float,
    samples: int,
    dt_reference_grid: np.ndarray,
) -> list[dict[str, float]]:
    """Compute jump-law speeds for a signed-distance level set.

    The D-branch thin-layer law is naturally written in the physical layer
    coordinate q=ell^2 R_tilde.  This prototype uses the local frozen layer
    compression alpha=|grad q_ref| to convert the target norm for dq into a
    target norm for the signed-distance covector dd.  This is a numerical
    representation test, not a replacement for full tensor matching.
    """
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dist_x, dist_z = np.gradient(distance_field, dx, dz, edge_order=2)

    q = np.linspace(branch_level - half_width_y, branch_level + half_width_y, samples)
    phi_prime = d_phi_d_raw_y(q)
    delta_phi_prime = float(phi_prime[-1] - phi_prime[0])
    phi_q = d_branch_phi(q)
    f_q = d_branch_f(q, ell)
    chi_q = q / (ell * ell)

    rows: list[dict[str, float]] = []
    for seg in find_level_segments(distance_field, rho, x, z, 0.0):
        xp = np.array([seg["x"]])
        zp = np.array([seg["z"]])
        fx = float(bilinear_sample(dist_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(dist_z, x, z, xp, zp)[0])
        grad_dist = float(np.hypot(fx, fz))
        if grad_dist <= 1.0e-300:
            continue

        alpha = float(bilinear_sample(raw_y_ref_grad_norm, x, z, xp, zp)[0])
        if alpha <= 1.0e-300:
            rows.append({**seg, "level": branch_level, "solved": 0.0, "failure_reason": "zero_layer_compression"})
            continue

        ginv_here = np.zeros((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                ginv_here[a, b] = float(bilinear_sample(metric_inv[..., a, b], x, z, xp, zp)[0])

        stress_here = float(bilinear_sample(stress_trace, x, z, xp, zp)[0])
        algebraic_q = phi_q * chi_q - 1.5 * f_q - stress_here
        algebraic_integral = float(np.trapezoid(algebraic_q, q))
        candidates = target_normal_norm_candidates(delta_phi_prime, algebraic_integral)
        dt_ref = float(bilinear_sample(dt_reference_grid, x, z, xp, zp)[0])

        best: dict[str, float] | None = None
        for cand in candidates:
            target_norm_distance = cand["normal_norm"] / (alpha * alpha)
            dt_corr, disc = solve_ft_for_target_norm(
                ginv=ginv_here,
                fx=fx,
                fz=fz,
                target_norm=target_norm_distance,
                ft_reference=dt_ref,
            )
            if not np.isfinite(dt_corr):
                continue
            speed_corr = float(-dt_corr / grad_dist)
            row = {
                **cand,
                "target_norm_distance": target_norm_distance,
                "dt_corrected": dt_corr,
                "quadratic_discriminant": disc,
                "coordinate_normal_speed_corrected": speed_corr,
                "abs_dt_change": abs(dt_corr - dt_ref),
            }
            if best is None or row["abs_dt_change"] < best["abs_dt_change"]:
                best = row

        base = {
            "x0": float(seg["x0"]),
            "z0": float(seg["z0"]),
            "x1": float(seg["x1"]),
            "z1": float(seg["z1"]),
            "x": float(seg["x"]),
            "z": float(seg["z"]),
            "level": float(branch_level),
            "segment_length": float(seg["length"]),
            "grad_distance_norm": grad_dist,
            "alpha_raw_y_per_distance": alpha,
            "algebraic_integral_raw_y": algebraic_integral,
            "dt_reference": dt_ref,
        }
        if best is None:
            rows.append({**base, "solved": 0.0, "failure_reason": "no_real_speed_root"})
            continue
        rows.append(
            {
                **base,
                "solved": 1.0,
                "dt_corrected": best["dt_corrected"],
                "coordinate_normal_speed_corrected": best["coordinate_normal_speed_corrected"],
                "target_signature": best["sigma"],
                "normal_norm_raw_y_target": best["normal_norm"],
                "normal_norm_distance_target": best["target_norm_distance"],
            }
        )
    return rows


def rows_to_speed_grid(
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
            row["coordinate_normal_speed_corrected"],
            max(row["segment_length"], 1.0e-12),
        )
    known = weights > 0.0
    speed = np.zeros_like(accum)
    speed[known] = accum[known] / weights[known]
    if extension_iterations > 0:
        speed = neighbor_average_extend(speed, known, iterations=extension_iterations)
    return speed, known


def evolve_distance(distance_field: np.ndarray, speed_grid: np.ndarray, x: np.ndarray, z: np.ndarray, dt: float) -> np.ndarray:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dist_x, dist_z = np.gradient(distance_field, dx, dz, edge_order=2)
    grad_norm = np.sqrt(dist_x * dist_x + dist_z * dist_z)
    return distance_field - dt * speed_grid * grad_norm


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_lines(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    colors: str | list[str] = "white",
    values: np.ndarray | list[float] | None = None,
    linewidth: float = 0.9,
    linestyle: str = "solid",
    cmap: str = "coolwarm",
    vlim: float | None = None,
) -> LineCollection | None:
    if not rows:
        return None
    if values is None:
        lc = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth, linestyles=linestyle)
    else:
        lc = LineCollection(row_segments(rows), cmap=cmap, linewidths=linewidth, linestyles=linestyle)
        lc.set_array(np.asarray(values, dtype=float))
        if vlim is not None:
            lc.set_clim(-vlim, vlim)
    ax.add_collection(lc)
    return lc


def render_frame(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    d_minus: np.ndarray,
    d_plus: np.ndarray,
    rows: list[dict[str, float]],
    step: int,
    dt: float,
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    unsolved = [row for row in rows if row.get("solved", 0.0) <= 0.5]
    line_colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in rows]
    nearest_abs_distance = np.minimum(np.abs(d_minus), np.abs(d_plus))

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 10.0), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))
    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_lines(axes[0, 0], rows, colors=line_colors, linewidth=0.85)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0, 0].set_title("linear rho + signed-distance interfaces")
    fig.colorbar(im0, ax=axes[0, 0])

    d_clip = max(float(np.percentile(nearest_abs_distance[support], 95.0)), 1.0) if np.any(support) else 1.0
    im1 = axes[0, 1].pcolormesh(xg, zg, nearest_abs_distance, shading="auto", cmap="magma", vmin=0.0, vmax=d_clip)
    add_lines(axes[0, 1], rows, colors=line_colors, linewidth=0.85)
    axes[0, 1].set_title("nearest |signed distance|")
    fig.colorbar(im1, ax=axes[0, 1])

    speeds = [row["coordinate_normal_speed_corrected"] for row in solved]
    vlim = max(float(np.percentile(np.abs(speeds), 95.0)), 1.0e-12) if speeds else 1.0
    lc = add_lines(axes[1, 0], solved, values=speeds, linewidth=1.05, vlim=vlim)
    add_lines(axes[1, 0], unsolved, colors="0.25", linewidth=0.65, linestyle="dotted")
    if lc is not None:
        fig.colorbar(lc, ax=axes[1, 0])
    axes[1, 0].set_title("jump-law normal speed; dotted unresolved")

    im3 = axes[1, 1].pcolormesh(xg, zg, d_plus - d_minus, shading="auto", cmap="coolwarm")
    add_lines(axes[1, 1], solved, colors="black", linewidth=0.75)
    add_lines(axes[1, 1], unsolved, colors="0.35", linewidth=0.65, linestyle="dotted")
    axes[1, 1].set_title("d_plus - d_minus branch separation")
    fig.colorbar(im3, ax=axes[1, 1])

    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
        ax.set_aspect("equal")
    fig.suptitle(f"D signed-distance moving-interface prototype, step={step}, pseudo time={step * dt:.3e}", fontsize=12)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_overlay(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    initial_rows: list[dict[str, float]],
    final_rows: list[dict[str, float]],
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    fig, ax = plt.subplots(figsize=(7.6, 7.0), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))
    im = ax.pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    fig.colorbar(im, ax=ax)
    add_lines(ax, initial_rows, colors="white", linewidth=0.7, linestyle="dashed")
    final_colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in final_rows]
    add_lines(ax, final_rows, colors=final_colors, linewidth=1.0)
    ax.contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dotted")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_aspect("equal")
    ax.set_title("initial dashed white; final orange(+1)/green(-1)")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def summarize_step(step: int, rows: list[dict[str, float]], reinit_ok: dict[str, bool]) -> dict[str, object]:
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    return {
        "step": step,
        "segment_count": len(rows),
        "solved_count": len(solved),
        "unresolved_count": len(rows) - len(solved),
        "solved_fraction": float(len(solved) / max(len(rows), 1)),
        "total_line_length": float(sum(row["segment_length"] for row in rows)),
        "solved_line_length": float(sum(row["segment_length"] for row in solved)),
        "speed_corrected": scalar_stats([row["coordinate_normal_speed_corrected"] for row in solved]),
        "normal_norm_raw_y_target": scalar_stats([row["normal_norm_raw_y_target"] for row in solved]),
        "normal_norm_distance_target": scalar_stats([row["normal_norm_distance_target"] for row in solved]),
        "alpha_raw_y_per_distance": scalar_stats([row["alpha_raw_y_per_distance"] for row in rows if "alpha_raw_y_per_distance" in row]),
        "reinitialized": reinit_ok,
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
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)
    raw_y_ref_grad_norm = np.sqrt(raw_y_x * raw_y_x + raw_y_z * raw_y_z)

    levels = parse_positive_levels(args.levels)
    if levels != [1.0]:
        raise ValueError("This prototype currently expects --levels 1 so d_minus/d_plus represent raw_y=-1/+1.")

    d_minus, initial_minus = initialize_branch_distance(raw_y, rho, x, z, -1.0)
    d_plus, initial_plus = initialize_branch_distance(raw_y, rho, x, z, 1.0)
    dt_guess_minus = np.zeros_like(raw_y)
    dt_guess_plus = np.zeros_like(raw_y)

    step_summaries: list[dict[str, object]] = []
    frame_paths: list[str] = []
    initial_rows_for_overlay: list[dict[str, float]] = []
    final_rows_for_overlay: list[dict[str, float]] = []

    for step in range(args.steps + 1):
        rows_minus = compute_distance_interface_rows(
            distance_field=d_minus,
            raw_y_ref_grad_norm=raw_y_ref_grad_norm,
            metric_inv=metric_inv,
            rho=rho,
            stress_trace=stress_trace,
            x=x,
            z=z,
            branch_level=-1.0,
            ell=args.ell,
            half_width_y=args.half_width_y,
            samples=args.samples,
            dt_reference_grid=dt_guess_minus,
        )
        rows_plus = compute_distance_interface_rows(
            distance_field=d_plus,
            raw_y_ref_grad_norm=raw_y_ref_grad_norm,
            metric_inv=metric_inv,
            rho=rho,
            stress_trace=stress_trace,
            x=x,
            z=z,
            branch_level=1.0,
            ell=args.ell,
            half_width_y=args.half_width_y,
            samples=args.samples,
            dt_reference_grid=dt_guess_plus,
        )
        rows = rows_minus + rows_plus
        if step == 0:
            initial_rows_for_overlay = rows
        final_rows_for_overlay = rows
        step_summaries.append(summarize_step(step, rows, {"minus": True, "plus": True}))

        if step % args.render_every == 0 or step == args.steps:
            frame_path = frame_dir / f"signed_distance_interface_step_{step:03d}.png"
            render_frame(
                frame_path,
                x=x,
                z=z,
                rho=rho,
                d_minus=d_minus,
                d_plus=d_plus,
                rows=rows,
                step=step,
                dt=args.dt,
                rho_min_frac=args.rho_min_frac,
            )
            frame_paths.append(str(frame_path.resolve()))
        if step == args.steps:
            break

        speed_minus, _ = rows_to_speed_grid(rows_minus, x=x, z=z, extension_iterations=args.extension_iterations)
        speed_plus, _ = rows_to_speed_grid(rows_plus, x=x, z=z, extension_iterations=args.extension_iterations)
        d_minus_tentative = evolve_distance(d_minus, speed_minus, x=x, z=z, dt=args.dt)
        d_plus_tentative = evolve_distance(d_plus, speed_plus, x=x, z=z, dt=args.dt)
        d_minus, _, ok_minus = reinitialize_distance(d_minus_tentative, rho, x, z)
        d_plus, _, ok_plus = reinitialize_distance(d_plus_tentative, rho, x, z)
        dt_guess_minus = -speed_minus
        dt_guess_plus = -speed_plus
        if not (ok_minus and ok_plus):
            break

    overlay_path = out / f"signed_distance_interface_overlay_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_overlay(
        overlay_path,
        x=x,
        z=z,
        rho=rho,
        initial_rows=initial_rows_for_overlay,
        final_rows=final_rows_for_overlay,
        rho_min_frac=args.rho_min_frac,
    )

    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "reference_time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "dt": args.dt,
            "steps_requested": args.steps,
            "steps_completed": len(step_summaries) - 1,
            "levels_abs_raw_y": levels,
            "half_width_y": args.half_width_y,
            "samples": args.samples,
            "rho_min_frac": args.rho_min_frac,
            "extension_iterations": args.extension_iterations,
        },
        "definitions": {
            "scope": "Reduced signed-distance moving-interface prototype with frozen rho, metric_inv, stress_trace and frozen raw_y layer-compression alpha. It tests the interface representation; it is not the final self-consistent tensor-matching dynamics.",
            "d_minus": "Signed distance to the raw_y=-1 branch, with sign convention sign(raw_y+1) at initialization and preserved by reinitialization.",
            "d_plus": "Signed distance to the raw_y=+1 branch, with sign convention sign(raw_y-1) at initialization and preserved by reinitialization.",
            "speed_law": "D-branch leading trace/scalaron jump law is converted from q=raw_y to signed-distance by alpha=|grad raw_y_ref|, then solved for d_t and normal speed.",
            "reinitialization": "After each Hamilton-Jacobi interface update, d is reset to distance from its zero contour. This is a numerical coordinate reset, not a new physical damping/source.",
            "unresolved_segments": "Dotted segments have no real root in the trace-only jump law and are not clipped or damped into artificial solutions.",
        },
        "initial_raw_y_segment_counts": {
            "minus": len(initial_minus),
            "plus": len(initial_plus),
            "total": len(initial_minus) + len(initial_plus),
        },
        "steps": step_summaries,
        "files": {
            "frames": frame_paths,
            "overlay": str(overlay_path.resolve()),
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
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=61)
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    parser.add_argument("--dt", type=float, default=5.0e-2)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--extension-iterations", type=int, default=10)
    parser.add_argument("--render-every", type=int, default=2)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
