from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

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


def solve_ft_for_target_norm(
    ginv: np.ndarray,
    fx: float,
    fz: float,
    target_norm: float,
    ft_reference: float,
) -> tuple[float, float]:
    """Solve g^{ab} F_a F_b = target_norm for F_t.

    The returned root is the one closest to the measured/reference F_t. This
    keeps the jump correction minimal when both quadratic roots are real.
    """
    a = float(ginv[0, 0])
    b = 2.0 * float(ginv[0, 1] * fx + ginv[0, 2] * fz)
    c = float(ginv[1, 1] * fx * fx + 2.0 * ginv[1, 2] * fx * fz + ginv[2, 2] * fz * fz - target_norm)
    if abs(a) <= 1.0e-300:
        if abs(b) <= 1.0e-300:
            return np.nan, np.nan
        root = -c / b
        return float(root), 0.0
    discriminant = b * b - 4.0 * a * c
    if discriminant < 0.0:
        return np.nan, float(discriminant)
    sqrt_disc = float(np.sqrt(discriminant))
    roots = np.array([(-b - sqrt_disc) / (2.0 * a), (-b + sqrt_disc) / (2.0 * a)], dtype=float)
    root = float(roots[np.argmin(np.abs(roots - ft_reference))])
    return root, float(discriminant)


def target_normal_norm_candidates(delta_phi_prime: float, algebraic_integral: float) -> list[dict[str, float]]:
    candidates: list[dict[str, float]] = []
    for sigma in (-1.0, 1.0):
        denom = 2.0 * sigma * delta_phi_prime
        if abs(denom) <= 1.0e-300:
            continue
        abs_norm_sq = -algebraic_integral / denom
        if abs_norm_sq < 0.0 or not np.isfinite(abs_norm_sq):
            continue
        candidates.append(
            {
                "sigma": float(sigma),
                "abs_normal_norm": float(abs_norm_sq),
                "normal_norm": float(sigma * abs_norm_sq),
            }
        )
    return candidates


def interface_speed_rows(
    raw_y: np.ndarray,
    raw_y_m: np.ndarray,
    raw_y_p: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    stress_trace: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    ell: float,
    probe_dt: float,
    half_width_y: float,
    samples: int,
) -> list[dict[str, float]]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_t = (raw_y_p - raw_y_m) / (2.0 * probe_dt)
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
        ft_ref = float(bilinear_sample(raw_y_t, x, z, xp, zp)[0])
        fx = float(bilinear_sample(raw_y_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(raw_y_z, x, z, xp, zp)[0])
        grad_spatial = float(np.hypot(fx, fz))
        if grad_spatial <= 1.0e-300:
            continue

        ginv_here = np.zeros((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                ginv_here[a, b] = float(bilinear_sample(metric_inv[..., a, b], x, z, xp, zp)[0])
        covector_ref = np.array([ft_ref, fx, fz], dtype=float)
        normal_norm_ref = float(covector_ref @ ginv_here @ covector_ref)
        speed_ref = float(-ft_ref / grad_spatial)

        stress_here = float(bilinear_sample(stress_trace, x, z, xp, zp)[0])
        algebraic_q = phi_q * chi_q - 1.5 * f_q - stress_here
        algebraic_integral = float(np.trapezoid(algebraic_q, q))
        algebraic_abs_integral = float(np.trapezoid(np.abs(algebraic_q), q))
        candidates = target_normal_norm_candidates(delta_phi_prime, algebraic_integral)

        best: dict[str, float] | None = None
        for cand in candidates:
            ft_corr, disc = solve_ft_for_target_norm(
                ginv=ginv_here,
                fx=fx,
                fz=fz,
                target_norm=cand["normal_norm"],
                ft_reference=ft_ref,
            )
            if not np.isfinite(ft_corr):
                continue
            speed_corr = float(-ft_corr / grad_spatial)
            row = {
                **cand,
                "ft_corrected": ft_corr,
                "quadratic_discriminant": disc,
                "coordinate_normal_speed_corrected": speed_corr,
                "abs_speed_change": abs(speed_corr - speed_ref),
                "abs_ft_change": abs(ft_corr - ft_ref),
            }
            if best is None or row["abs_ft_change"] < best["abs_ft_change"]:
                best = row

        if best is None:
            rows.append(
                {
                    "x0": float(seg["x0"]),
                    "z0": float(seg["z0"]),
                    "x1": float(seg["x1"]),
                    "z1": float(seg["z1"]),
                    "x": float(seg["x"]),
                    "z": float(seg["z"]),
                    "level": float(level),
                    "solved": 0.0,
                    "ft_reference": ft_ref,
                    "coordinate_normal_speed_reference": speed_ref,
                    "normal_norm_reference": normal_norm_ref,
                    "spatial_grad_norm": grad_spatial,
                    "delta_phi_prime": delta_phi_prime,
                    "algebraic_integral_raw_y": algebraic_integral,
                    "algebraic_abs_integral_raw_y": algebraic_abs_integral,
                }
            )
            continue

        normal_norm_corr = best["normal_norm"]
        scaled_residual_corr = float(2.0 * best["sigma"] * abs(normal_norm_corr) * delta_phi_prime + algebraic_integral)
        rows.append(
            {
                "x0": float(seg["x0"]),
                "z0": float(seg["z0"]),
                "x1": float(seg["x1"]),
                "z1": float(seg["z1"]),
                "x": float(seg["x"]),
                "z": float(seg["z"]),
                "level": float(level),
                "solved": 1.0,
                "ft_reference": ft_ref,
                "ft_corrected": best["ft_corrected"],
                "coordinate_normal_speed_reference": speed_ref,
                "coordinate_normal_speed_corrected": best["coordinate_normal_speed_corrected"],
                "abs_speed_change": best["abs_speed_change"],
                "abs_ft_change": best["abs_ft_change"],
                "normal_norm_reference": normal_norm_ref,
                "normal_norm_corrected": normal_norm_corr,
                "target_signature": best["sigma"],
                "reference_to_target_norm_ratio": float(
                    abs(normal_norm_ref) / max(abs(normal_norm_corr), 1.0e-300)
                ),
                "spatial_grad_norm": grad_spatial,
                "delta_phi_prime": delta_phi_prime,
                "algebraic_integral_raw_y": algebraic_integral,
                "algebraic_abs_integral_raw_y": algebraic_abs_integral,
                "scaled_residual_corrected": scaled_residual_corr,
                "nx_spatial": fx / grad_spatial,
                "nz_spatial": fz / grad_spatial,
                "x_euler": float(seg["x"] + best["coordinate_normal_speed_corrected"] * fx / grad_spatial),
                "z_euler": float(seg["z"] + best["coordinate_normal_speed_corrected"] * fz / grad_spatial),
            }
        )
    return rows


def summarize(rows: list[dict[str, float]]) -> dict[str, object]:
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    unsolved = len(rows) - len(solved)
    return {
        "count": len(rows),
        "solved_count": len(solved),
        "unsolved_count": unsolved,
        "solved_fraction": float(len(solved) / max(len(rows), 1)),
        "coordinate_normal_speed_reference": scalar_stats([row["coordinate_normal_speed_reference"] for row in solved]),
        "coordinate_normal_speed_corrected": scalar_stats([row["coordinate_normal_speed_corrected"] for row in solved]),
        "abs_speed_change": scalar_stats([row["abs_speed_change"] for row in solved]),
        "normal_norm_reference": scalar_stats([row["normal_norm_reference"] for row in solved]),
        "normal_norm_corrected": scalar_stats([row["normal_norm_corrected"] for row in solved]),
        "reference_to_target_norm_ratio": scalar_stats([row["reference_to_target_norm_ratio"] for row in solved]),
        "algebraic_integral_raw_y": scalar_stats([row["algebraic_integral_raw_y"] for row in solved]),
        "scaled_residual_corrected": scalar_stats([row["scaled_residual_corrected"] for row in solved]),
        "target_signature_positive_fraction": float(
            np.mean([row["target_signature"] > 0.0 for row in solved]) if solved else 0.0
        ),
    }


def render_plot(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    rows: list[dict[str, float]],
    levels: list[float],
    interface_dt: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    unsolved = [row for row in rows if row.get("solved", 0.0) <= 0.5]
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 9.0), constrained_layout=True)

    im0 = axes[0, 0].pcolormesh(xg, zg, np.log10(np.maximum(rho, 1.0e-16)), shading="auto", cmap="viridis")
    add_interface_lines(axes[0, 0], rows, color_by="level")
    axes[0, 0].set_title("log10 rho; lines: raw_y=+/-levels")
    fig.colorbar(im0, ax=axes[0, 0])

    if solved:
        speed_ref = np.array([row["coordinate_normal_speed_reference"] for row in solved])
        speed_corr = np.array([row["coordinate_normal_speed_corrected"] for row in solved])
        speed_delta = speed_corr - speed_ref
        vmax_ref = max(float(np.percentile(np.abs(speed_ref), 95.0)), 1.0e-12)
        vmax_corr = max(float(np.percentile(np.abs(speed_corr), 95.0)), 1.0e-12)
        vmax_delta = max(float(np.percentile(np.abs(speed_delta), 95.0)), 1.0e-12)
        lc1 = add_interface_lines(
            axes[0, 1],
            solved,
            values=speed_ref,
            cmap="coolwarm",
            vmin=-vmax_ref,
            vmax=vmax_ref,
            linewidth=1.15,
        )
        add_unsolved_lines(axes[0, 1], unsolved)
        fig.colorbar(lc1, ax=axes[0, 1])
        lc2 = add_interface_lines(
            axes[1, 0],
            solved,
            values=speed_corr,
            cmap="coolwarm",
            vmin=-vmax_corr,
            vmax=vmax_corr,
            linewidth=1.15,
        )
        add_unsolved_lines(axes[1, 0], unsolved)
        fig.colorbar(lc2, ax=axes[1, 0])
        lc3 = add_interface_lines(
            axes[1, 1],
            solved,
            values=speed_delta,
            cmap="coolwarm",
            vmin=-vmax_delta,
            vmax=vmax_delta,
            linewidth=1.15,
        )
        add_unsolved_lines(axes[1, 1], unsolved)
        fig.colorbar(lc3, ax=axes[1, 1])
        if interface_dt > 0.0:
            xs = np.array([row["x"] for row in solved])
            zs = np.array([row["z"] for row in solved])
            nx = np.array([row["nx_spatial"] for row in solved])
            nz = np.array([row["nz_spatial"] for row in solved])
            axes[1, 0].quiver(xs, zs, interface_dt * speed_corr * nx, interface_dt * speed_corr * nz, color="black", width=0.002, scale_units="xy", angles="xy", scale=1.0, alpha=0.55)

    axes[0, 1].set_title("reference coordinate normal speed")
    axes[1, 0].set_title("jump-law corrected speed; arrows: one Euler step")
    axes[1, 1].set_title("corrected speed - reference speed")
    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
        ax.set_aspect("equal")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_unsolved_lines(ax: plt.Axes, rows: list[dict[str, float]]) -> None:
    if not rows:
        return
    lc = LineCollection(row_segments(rows), colors="0.25", linewidths=0.75, linestyles="dotted", alpha=0.8)
    ax.add_collection(lc)


def add_interface_lines(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    values: np.ndarray | None = None,
    cmap: str = "coolwarm",
    vmin: float | None = None,
    vmax: float | None = None,
    linewidth: float = 0.9,
    color_by: str | None = None,
) -> LineCollection | None:
    if not rows:
        return None
    segments = row_segments(rows)
    if values is not None:
        lc = LineCollection(segments, cmap=cmap, linewidths=linewidth)
        lc.set_array(np.asarray(values, dtype=float))
        if vmin is not None and vmax is not None:
            lc.set_clim(vmin, vmax)
    elif color_by == "level":
        colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in rows]
        lc = LineCollection(segments, colors=colors, linewidths=linewidth, alpha=0.95)
    else:
        lc = LineCollection(segments, colors="white", linewidths=linewidth, alpha=0.95)
    ax.add_collection(lc)
    return lc


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
    raw_y = args.ell * args.ell * real_array(ref["chi_ref"])
    raw_y_m = args.ell * args.ell * real_array(ref["chi_m"])
    raw_y_p = args.ell * args.ell * real_array(ref["chi_p"])
    rho = real_array(ref["rho"])
    metric_inv = real_array(ref["metric_inv"])
    stress_trace = real_array(ref["stress_trace"])
    levels = parse_positive_levels(args.levels)

    rows: list[dict[str, float]] = []
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        rows.extend(
            interface_speed_rows(
                raw_y=raw_y,
                raw_y_m=raw_y_m,
                raw_y_p=raw_y_p,
                metric_inv=metric_inv,
                rho=rho,
                stress_trace=stress_trace,
                x=x,
                z=z,
                level=level,
                ell=args.ell,
                probe_dt=args.probe_dt,
                half_width_y=args.half_width_y,
                samples=args.samples,
            )
        )

    fig_path = out / f"d_interface_speed_law_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_plot(fig_path, x=x, z=z, rho=rho, raw_y=raw_y, rows=rows, levels=levels, interface_dt=args.interface_dt)
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
            "interface_dt_for_arrows": args.interface_dt,
        },
        "definitions": {
            "raw_y": "raw_y=ell^2 R_tilde; interfaces are raw_y=+/-level.",
            "speed_law": "Solve 2*sigma*|g^ab F_a F_b|*Delta(d_phi/d raw_y)+Integral(algebraic d raw_y)=0 for g^ab F_a F_b, then solve the quadratic relation for F_t.",
            "coordinate_normal_speed": "v_n=-F_t/sqrt(F_x^2+F_z^2) on the current t,x,z coordinate grid.",
            "corrected_speed": "The minimal-change root, chosen closest to the reference F_t among real roots.",
            "scope": "Leading trace/scalaron spacetime-interface condition only; this is the interface-dynamics core for the reduced solver, not the final full tensor matching system.",
        },
        "summary": summarize(rows),
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
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument("--interface-dt", type=float, default=2.0e-4)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
