from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from diagnose_d_spacetime_interface_jump import d_phi_d_raw_y
from diagnose_d_trace_interface_weak_form import find_level_segments, parse_positive_levels, real_array
from diagnose_d_tensor_interface_jump import (
    algebraic_tensor_integral,
    reference_tensor_data,
    sample_matrix,
    tensor_norm,
    trace_with_inverse,
)
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def direct_tensor_covector(
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    algebraic_tensor: np.ndarray,
    delta_phi_prime: float,
    current_spatial_normal: np.ndarray,
) -> dict[str, float | np.ndarray]:
    """Solve the leading tensor jump condition directly for a covector.

    The tensor equation is
        Delta(phi_q) * (-F_a F_b + g_ab F^2) + A_ab = 0.
    Define H_ab=-A_ab/Delta(phi_q). Then in 2+1 dimensions
        F^2 = Tr_g(H)/2,
        F_a F_b = g_ab F^2 - H_ab.
    A real covector exists only when the right-hand side is rank-one positive
    semi-definite, up to numerical error.
    """
    if abs(delta_phi_prime) <= 1.0e-300:
        return {"candidate": 0.0, "failure_reason": "zero_delta_phi_prime"}

    target_h = -algebraic_tensor / delta_phi_prime
    target_norm = 0.5 * trace_with_inverse(metric_inv, target_h)
    k_cov = metric_cov * target_norm - target_h
    k_sym = 0.5 * (k_cov + k_cov.T)
    evals, evecs = np.linalg.eigh(k_sym)
    positive = np.where(evals > 0.0)[0]
    if positive.size == 0:
        svals = np.linalg.svd(k_sym, compute_uv=False)
        return {
            "candidate": 0.0,
            "failure_reason": "no_positive_rank1_eigenvalue",
            "direct_target_norm": float(target_norm),
            "rank1_tail_relative": float(np.sqrt(np.sum(svals * svals)) / max(float(np.max(svals)), 1.0e-300)),
        }

    lead_idx = int(positive[np.argmax(evals[positive])])
    lead_eval = float(evals[lead_idx])
    covector = np.sqrt(max(lead_eval, 0.0)) * np.asarray(evecs[:, lead_idx], dtype=float)
    spatial = covector[1:3]
    if np.linalg.norm(spatial) > 1.0e-300 and np.dot(spatial, current_spatial_normal) < 0.0:
        covector = -covector
        spatial = -spatial

    k_rank1 = np.outer(covector, covector)
    rank_residual = k_sym - k_rank1
    svals = np.linalg.svd(k_sym, compute_uv=False)
    rank_tail = float(np.sqrt(np.sum(np.delete(evals, lead_idx) ** 2)))
    rank_scale = max(float(np.sqrt(np.sum(evals * evals))), 1.0e-300)
    negative_weight = float(np.sum(np.abs(evals[evals < 0.0])))
    negative_relative = negative_weight / max(abs(lead_eval), 1.0e-300)

    norm_from_covector = float(covector @ metric_inv @ covector)
    singular_tensor = delta_phi_prime * (-np.outer(covector, covector) + metric_cov * norm_from_covector)
    tensor_residual = singular_tensor + algebraic_tensor
    algebraic_norm = max(tensor_norm(algebraic_tensor), 1.0e-300)
    residual_relative = tensor_norm(tensor_residual) / algebraic_norm

    spatial_norm = float(np.linalg.norm(spatial))
    speed = float(-covector[0] / spatial_norm) if spatial_norm > 1.0e-300 else np.nan
    alignment = np.nan
    cur_norm = float(np.linalg.norm(current_spatial_normal))
    if spatial_norm > 1.0e-300 and cur_norm > 1.0e-300:
        alignment = float(abs(np.dot(spatial, current_spatial_normal)) / (spatial_norm * cur_norm))

    return {
        "candidate": 1.0,
        "covector": covector,
        "direct_target_norm": float(target_norm),
        "direct_covector_norm": norm_from_covector,
        "direct_norm_error": norm_from_covector - target_norm,
        "direct_speed": speed,
        "direct_spatial_norm": spatial_norm,
        "direct_tensor_residual_norm": tensor_norm(tensor_residual),
        "direct_tensor_residual_relative": residual_relative,
        "rank1_tail_relative": rank_tail / rank_scale,
        "rank1_negative_relative": negative_relative,
        "rank1_residual_relative": tensor_norm(rank_residual) / max(tensor_norm(k_sym), 1.0e-300),
        "rank1_leading_eigenvalue": lead_eval,
        "rank1_s0": float(svals[0]) if svals.size else 0.0,
        "rank1_s1": float(svals[1]) if svals.size > 1 else 0.0,
        "rank1_s2": float(svals[2]) if svals.size > 2 else 0.0,
        "spatial_alignment_with_current": alignment,
        "target_covector_t": float(covector[0]),
        "target_covector_x": float(covector[1]),
        "target_covector_z": float(covector[2]),
    }


def direct_tensor_rows(
    raw_y: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    ricci: np.ndarray,
    rho: np.ndarray,
    rhs_tensor: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
    ell: float,
    half_width_y: float,
    samples: int,
    rank1_tol: float,
    tensor_rel_tol: float,
) -> list[dict[str, float]]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)
    q = np.linspace(level - half_width_y, level + half_width_y, samples)
    delta_phi_prime = float(d_phi_d_raw_y(q)[-1] - d_phi_d_raw_y(q)[0])

    rows: list[dict[str, float]] = []
    for seg in find_level_segments(raw_y, rho, x, z, level):
        xp = np.array([seg["x"]])
        zp = np.array([seg["z"]])
        fx = float(bilinear_sample(raw_y_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(raw_y_z, x, z, xp, zp)[0])
        current_spatial = np.array([fx, fz], dtype=float)
        g_cov = sample_matrix(metric_cov, x, z, xp, zp)
        g_inv = np.linalg.pinv(g_cov, rcond=1.0e-12, hermitian=True)
        algebraic_tensor, trace_from_tensor, trace_direct = algebraic_tensor_integral(
            q=q,
            level=level,
            ell=ell,
            metric_cov=g_cov,
            metric_inv=g_inv,
            ricci=sample_matrix(ricci, x, z, xp, zp),
            rhs_tensor=sample_matrix(rhs_tensor, x, z, xp, zp),
        )
        result = direct_tensor_covector(
            metric_cov=g_cov,
            metric_inv=g_inv,
            algebraic_tensor=algebraic_tensor,
            delta_phi_prime=delta_phi_prime,
            current_spatial_normal=current_spatial,
        )
        base: dict[str, float] = {
            "x0": float(seg["x0"]),
            "z0": float(seg["z0"]),
            "x1": float(seg["x1"]),
            "z1": float(seg["z1"]),
            "x": float(seg["x"]),
            "z": float(seg["z"]),
            "level": float(level),
            "segment_length": float(seg["length"]),
            "delta_phi_prime": delta_phi_prime,
            "trace_algebraic_from_tensor": trace_from_tensor,
            "trace_algebraic_direct": trace_direct,
            "trace_integral_mismatch": trace_from_tensor - trace_direct,
            "algebraic_tensor_norm": tensor_norm(algebraic_tensor),
        }
        for key, value in result.items():
            if key == "covector" or isinstance(value, str):
                continue
            base[key] = float(value)
        base["direct_tensor_solved"] = float(
            result.get("candidate", 0.0) > 0.5
            and base.get("rank1_tail_relative", np.inf) <= rank1_tol
            and base.get("rank1_negative_relative", np.inf) <= rank1_tol
            and base.get("direct_tensor_residual_relative", np.inf) <= tensor_rel_tol
        )
        rows.append(base)
    return rows


def summarize(rows: list[dict[str, float]]) -> dict[str, object]:
    candidates = [row for row in rows if row.get("candidate", 0.0) > 0.5]
    solved = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    return {
        "count": len(rows),
        "candidate_count": len(candidates),
        "candidate_fraction": float(len(candidates) / max(len(rows), 1)),
        "direct_tensor_solved_count": len(solved),
        "direct_tensor_solved_fraction": float(len(solved) / max(len(rows), 1)),
        "rank1_tail_relative": scalar_stats([row["rank1_tail_relative"] for row in candidates]),
        "rank1_negative_relative": scalar_stats([row["rank1_negative_relative"] for row in candidates]),
        "rank1_residual_relative": scalar_stats([row["rank1_residual_relative"] for row in candidates]),
        "direct_tensor_residual_relative": scalar_stats([row["direct_tensor_residual_relative"] for row in candidates]),
        "direct_norm_error": scalar_stats([row["direct_norm_error"] for row in candidates]),
        "direct_speed": scalar_stats([row["direct_speed"] for row in candidates]),
        "spatial_alignment_with_current": scalar_stats([row["spatial_alignment_with_current"] for row in candidates]),
        "direct_target_norm": scalar_stats([row["direct_target_norm"] for row in candidates]),
        "trace_integral_mismatch": scalar_stats([row["trace_integral_mismatch"] for row in rows]),
        "accepted": {
            "rank1_tail_relative": scalar_stats([row["rank1_tail_relative"] for row in solved]),
            "rank1_negative_relative": scalar_stats([row["rank1_negative_relative"] for row in solved]),
            "direct_tensor_residual_relative": scalar_stats([row["direct_tensor_residual_relative"] for row in solved]),
            "direct_norm_error": scalar_stats([row["direct_norm_error"] for row in solved]),
            "direct_speed": scalar_stats([row["direct_speed"] for row in solved]),
            "spatial_alignment_with_current": scalar_stats([row["spatial_alignment_with_current"] for row in solved]),
            "direct_target_norm": scalar_stats([row["direct_target_norm"] for row in solved]),
        },
    }


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_lines(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    colors: str | list[str] = "white",
    values: list[float] | np.ndarray | None = None,
    linewidth: float = 0.9,
    linestyle: str = "solid",
    cmap: str = "magma",
    vmin: float | None = None,
    vmax: float | None = None,
) -> LineCollection | None:
    if not rows:
        return None
    if values is None:
        lc = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth, linestyles=linestyle)
    else:
        lc = LineCollection(row_segments(rows), cmap=cmap, linewidths=linewidth, linestyles=linestyle)
        lc.set_array(np.asarray(values, dtype=float))
        if vmin is not None and vmax is not None:
            lc.set_clim(vmin, vmax)
    ax.add_collection(lc)
    return lc


def render_plot(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    rows: list[dict[str, float]],
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    candidates = [row for row in rows if row.get("candidate", 0.0) > 0.5]
    solved = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    rejected = [row for row in rows if row.get("direct_tensor_solved", 0.0) <= 0.5]
    colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in rows]

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.8), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))
    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_lines(axes[0, 0], rows, colors=colors, linewidth=0.85)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0, 0].set_title("linear rho + raw_y=+/-1")
    fig.colorbar(im0, ax=axes[0, 0])

    add_lines(axes[0, 1], solved, colors="black", linewidth=0.9)
    add_lines(axes[0, 1], rejected, colors="0.35", linewidth=0.65, linestyle="dotted")
    axes[0, 1].set_title("direct tensor accepted; dotted rejected")

    if candidates:
        vals = np.array([row["rank1_tail_relative"] for row in candidates], dtype=float)
        vmax = max(float(np.percentile(vals[np.isfinite(vals)], 95.0)), 1.0e-12)
        lc = add_lines(axes[0, 2], candidates, values=vals, vmin=0.0, vmax=vmax, linewidth=1.0)
        fig.colorbar(lc, ax=axes[0, 2])
    axes[0, 2].set_title("rank-one mismatch")

    if candidates:
        vals = np.array([row["direct_tensor_residual_relative"] for row in candidates], dtype=float)
        vmax = max(float(np.percentile(vals[np.isfinite(vals)], 95.0)), 1.0e-12)
        lc = add_lines(axes[1, 0], candidates, values=vals, vmin=0.0, vmax=vmax, linewidth=1.0)
        fig.colorbar(lc, ax=axes[1, 0])

        vals = np.array([row["spatial_alignment_with_current"] for row in candidates], dtype=float)
        lc = add_lines(axes[1, 1], candidates, values=vals, vmin=0.0, vmax=1.0, linewidth=1.0, cmap="viridis")
        fig.colorbar(lc, ax=axes[1, 1])

        vals = np.array([row["direct_speed"] for row in candidates], dtype=float)
        vmax = max(float(np.percentile(np.abs(vals[np.isfinite(vals)]), 95.0)), 1.0e-12)
        lc = add_lines(axes[1, 2], candidates, values=vals, vmin=-vmax, vmax=vmax, linewidth=1.0, cmap="coolwarm")
        fig.colorbar(lc, ax=axes[1, 2])
    axes[1, 0].set_title("direct tensor residual / algebraic")
    axes[1, 1].set_title("target spatial normal alignment")
    axes[1, 2].set_title("direct tensor speed")

    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
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
    data = reference_tensor_data(
        x=x,
        z=z,
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        t=args.time,
        probe_dt=args.probe_dt,
        ell=args.ell,
        mp=args.mp,
        mass=params.m,
        rho_floor=args.rho_floor,
    )
    rows: list[dict[str, float]] = []
    levels = parse_positive_levels(args.levels)
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        rows.extend(
            direct_tensor_rows(
                raw_y=real_array(data["raw_y"]),
                metric_cov=real_array(data["metric_cov"]),
                metric_inv=real_array(data["metric_inv"]),
                ricci=real_array(data["ricci"]),
                rho=real_array(data["rho"]),
                rhs_tensor=real_array(data["rhs_tensor"]),
                x=x,
                z=z,
                level=level,
                ell=args.ell,
                half_width_y=args.half_width_y,
                samples=args.samples,
                rank1_tol=args.rank1_tol,
                tensor_rel_tol=args.tensor_rel_tol,
            )
        )

    fig_path = out / f"d_direct_tensor_interface_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render_plot(fig_path, x=x, z=z, rho=real_array(data["rho"]), rows=rows, rho_min_frac=args.rho_min_frac)
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
            "rho_min_frac": args.rho_min_frac,
            "rank1_tol": args.rank1_tol,
            "tensor_rel_tol": args.tensor_rel_tol,
        },
        "definitions": {
            "direct_tensor_condition": "Solves the leading tensor jump condition directly: H_ab=-F_a F_b+g_ab F^2 with H_ab=-A_ab/Delta(d f_R/dq). Trace is only a derived check, not the solve condition.",
            "direct_tensor_solved": "A numerical acceptance flag requiring a real rank-one covector and residuals below the configured tolerances. The tolerance is numerical, not a new physical damping rule.",
            "rank1_tail_relative": "How far g_ab F^2-H_ab is from rank one after removing the best rank-one piece.",
            "spatial_alignment_with_current": "Alignment between the direct tensor target spatial covector and the current raw_y interface spatial normal.",
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
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    parser.add_argument("--rank1-tol", type=float, default=0.1)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
