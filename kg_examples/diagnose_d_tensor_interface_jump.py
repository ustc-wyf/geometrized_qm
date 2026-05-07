from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from diagnose_d_interface_speed_law import solve_ft_for_target_norm
from diagnose_d_spacetime_interface_jump import d_branch_f, d_branch_phi, d_phi_d_raw_y
from diagnose_d_trace_interface_weak_form import find_level_segments, parse_positive_levels, real_array
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def tensor_norm(tensor: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.asarray(tensor, dtype=float) ** 2)))


def sample_matrix(field: np.ndarray, x: np.ndarray, z: np.ndarray, xp: np.ndarray, zp: np.ndarray) -> np.ndarray:
    out = np.zeros(field.shape[-2:], dtype=float)
    for a in range(field.shape[-2]):
        for b in range(field.shape[-1]):
            out[a, b] = float(bilinear_sample(field[..., a, b], x, z, xp, zp)[0])
    return out


def trace_with_inverse(metric_inv: np.ndarray, tensor: np.ndarray) -> float:
    return float(np.einsum("ab,ab->", metric_inv, tensor, optimize=True))


def traceless_part(metric_cov: np.ndarray, metric_inv: np.ndarray, tensor: np.ndarray) -> np.ndarray:
    dim = float(metric_cov.shape[0])
    return tensor - (trace_with_inverse(metric_inv, tensor) / dim) * metric_cov


def rank_one_compatibility(metric_cov: np.ndarray, metric_inv: np.ndarray, target_h: np.ndarray) -> dict[str, float]:
    """Check whether target_h can be written as -F_a F_b + g_ab F^2.

    If H_ab=-F_aF_b+g_ab N and n=3, then Tr_g H=2N and
    K_ab=g_ab N-H_ab must equal F_aF_b, hence be rank one.
    """
    target_n = 0.5 * trace_with_inverse(metric_inv, target_h)
    k_cov = metric_cov * target_n - target_h
    svals = np.linalg.svd(k_cov, compute_uv=False)
    evals, evecs = np.linalg.eigh(0.5 * (k_cov + k_cov.T))
    lead_idx = int(np.argmax(np.abs(evals)))
    lead_eval = float(evals[lead_idx])
    lead_vec = np.asarray(evecs[:, lead_idx], dtype=float)
    negative_weight = float(np.sum(np.abs(evals[evals < 0.0])))
    s0 = max(float(svals[0]), 1.0e-300)
    rank_tail = float(np.sqrt(np.sum(svals[1:] * svals[1:])))
    metric_norm_error = trace_with_inverse(metric_inv, k_cov) - target_n
    return {
        "free_covector_target_norm": target_n,
        "rank1_s0": float(svals[0]),
        "rank1_s1": float(svals[1]) if svals.size > 1 else 0.0,
        "rank1_s2": float(svals[2]) if svals.size > 2 else 0.0,
        "rank1_tail_relative": rank_tail / s0,
        "rank1_leading_eigenvalue": lead_eval,
        "rank1_negative_weight_relative": negative_weight / s0,
        "free_covector_t_unit": float(lead_vec[0]),
        "free_covector_x_unit": float(lead_vec[1]),
        "free_covector_z_unit": float(lead_vec[2]),
        "rank1_metric_norm_error": metric_norm_error,
    }


def reference_tensor_data(
    x: np.ndarray,
    z: np.ndarray,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    t: float,
    probe_dt: float,
    ell: float,
    mp: float,
    mass: float,
    rho_floor: float,
) -> dict[str, np.ndarray]:
    snap_mm = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - 2.0 * probe_dt, rho_floor=rho_floor)
    snap_m = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - probe_dt, rho_floor=rho_floor)
    snap_0 = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t, rho_floor=rho_floor)
    snap_p = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + probe_dt, rho_floor=rho_floor)
    snap_pp = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + 2.0 * probe_dt, rho_floor=rho_floor)

    metric_mm = real_array(snap_mm["cov_txz"])
    metric_m = real_array(snap_m["cov_txz"])
    metric_0 = real_array(snap_0["cov_txz"])
    metric_p = real_array(snap_p["cov_txz"])
    metric_pp = real_array(snap_pp["cov_txz"])

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, probe_dt, dx, dz)
    metric_inv, _, ricci, r_scalar = geometry_data(metric_0, dg, d2g)

    dg_m, d2g_m = metric_jets_full(metric_mm, metric_m, metric_0, probe_dt, dx, dz)
    _, _, _, r_scalar_m = geometry_data(metric_m, dg_m, d2g_m)
    dg_p, d2g_p = metric_jets_full(metric_0, metric_p, metric_pp, probe_dt, dx, dz)
    _, _, _, r_scalar_p = geometry_data(metric_p, dg_p, d2g_p)

    bohm = snap_0["bohm"]
    rho = real_array(bohm["rho"])
    s_t = real_array(bohm["s_t"])
    s_x = real_array(bohm["s_x"])
    s_z = real_array(bohm["s_z"])
    stress_tensor, _ = stress_tensor_tilde(metric_0, metric_inv, rho, s_t, s_x, s_z, m=mass)

    l2 = ell * ell
    return {
        "metric_cov": metric_0,
        "metric_inv": metric_inv,
        "ricci": real_array(ricci),
        "raw_y": l2 * real_array(r_scalar),
        "raw_y_m": l2 * real_array(r_scalar_m),
        "raw_y_p": l2 * real_array(r_scalar_p),
        "rho": rho,
        "rhs_tensor": real_array(stress_tensor) / (mp * mp),
    }


def algebraic_tensor_integral(
    q: np.ndarray,
    level: float,
    ell: float,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    ricci: np.ndarray,
    rhs_tensor: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    """Thin-layer algebraic tensor integral with trace-consistent Ricci scaling.

    The sampled Ricci tensor is first projected to its traceless part using the
    sampled metric.  The pure-trace Ricci part is then rebuilt so that the scalar
    curvature follows R=q/ell^2 across the layer.  This keeps the tensor
    diagnostic consistent with the earlier trace jump law while avoiding an
    arbitrary full Ricci profile through the layer.
    """
    l2 = ell * ell
    phi_q = d_branch_phi(q)
    f_q = d_branch_f(q, ell)
    width = float(q[-1] - q[0])
    ricci_traceless = traceless_part(metric_cov, metric_inv, ricci)
    int_phi = float(np.trapezoid(phi_q, q))
    int_phi_trace = float(np.trapezoid(phi_q * q / (3.0 * l2), q))
    int_f = float(np.trapezoid(f_q, q))
    tensor = int_phi * ricci_traceless + int_phi_trace * metric_cov - 0.5 * int_f * metric_cov - width * rhs_tensor

    trace_integral_from_tensor = trace_with_inverse(metric_inv, tensor)
    stress_trace = trace_with_inverse(metric_inv, rhs_tensor)
    trace_integral_direct = float(np.trapezoid(phi_q * (q / l2) - 1.5 * f_q - stress_trace, q))
    return tensor, trace_integral_from_tensor, trace_integral_direct


def tensor_jump_rows(
    raw_y: np.ndarray,
    raw_y_m: np.ndarray,
    raw_y_p: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    ricci: np.ndarray,
    rho: np.ndarray,
    rhs_tensor: np.ndarray,
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

    rows: list[dict[str, float]] = []
    for seg in find_level_segments(raw_y, rho, x, z, level):
        xp = np.array([seg["x"]])
        zp = np.array([seg["z"]])
        ft_ref = float(bilinear_sample(raw_y_t, x, z, xp, zp)[0])
        fx = float(bilinear_sample(raw_y_x, x, z, xp, zp)[0])
        fz = float(bilinear_sample(raw_y_z, x, z, xp, zp)[0])
        g_cov = sample_matrix(metric_cov, x, z, xp, zp)
        # In a body-fitted/subcell calculation the inverse should be computed
        # from the interpolated metric itself.  Interpolating g_ab and g^ab
        # separately breaks g^ab g_bc=delta^a_c and contaminates the tensor
        # trace condition.
        g_inv = np.linalg.pinv(g_cov, rcond=1.0e-12, hermitian=True)
        ricci_here = sample_matrix(ricci, x, z, xp, zp)
        rhs_here = sample_matrix(rhs_tensor, x, z, xp, zp)
        algebraic_tensor, trace_from_tensor, trace_direct = algebraic_tensor_integral(
            q=q,
            level=level,
            ell=ell,
            metric_cov=g_cov,
            metric_inv=g_inv,
            ricci=ricci_here,
            rhs_tensor=rhs_here,
        )

        target_norm = np.nan
        if abs(delta_phi_prime) > 1.0e-300:
            target_norm = float(-trace_from_tensor / (2.0 * delta_phi_prime))
        ft_corr, disc = solve_ft_for_target_norm(
            ginv=g_inv,
            fx=fx,
            fz=fz,
            target_norm=target_norm,
            ft_reference=ft_ref,
        )

        base = {
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
            "target_normal_norm": target_norm,
            "quadratic_discriminant": float(disc),
            "algebraic_tensor_norm": tensor_norm(algebraic_tensor),
            "algebraic_traceless_norm": tensor_norm(traceless_part(g_cov, g_inv, algebraic_tensor)),
            "ft_reference": ft_ref,
        }
        if abs(delta_phi_prime) > 1.0e-300:
            target_h = -algebraic_tensor / delta_phi_prime
            base.update(rank_one_compatibility(g_cov, g_inv, target_h))
            free_spatial = np.array([base["free_covector_x_unit"], base["free_covector_z_unit"]], dtype=float)
            current_spatial = np.array([fx, fz], dtype=float)
            denom = float(np.linalg.norm(free_spatial) * np.linalg.norm(current_spatial))
            base["free_to_current_spatial_alignment"] = (
                float(abs(np.dot(free_spatial, current_spatial)) / denom) if denom > 1.0e-300 else np.nan
            )
        if not np.isfinite(ft_corr):
            rows.append({**base, "solved": 0.0})
            continue

        covector = np.array([ft_corr, fx, fz], dtype=float)
        normal_norm = float(covector @ g_inv @ covector)
        normal_norm_error = normal_norm - target_norm
        root_tolerance = 1.0e-4 * max(1.0, abs(target_norm))
        if not np.isfinite(normal_norm_error) or abs(normal_norm_error) > root_tolerance:
            rows.append(
                {
                    **base,
                    "solved": 0.0,
                    "failure_reason": "root_verification_failed",
                    "ft_candidate": float(ft_corr),
                    "normal_norm_candidate": normal_norm,
                    "normal_norm_error": normal_norm_error,
                }
            )
            continue
        singular_tensor = delta_phi_prime * (-np.outer(covector, covector) + g_cov * normal_norm)
        residual_tensor = singular_tensor + algebraic_tensor
        residual_traceless = traceless_part(g_cov, g_inv, residual_tensor)
        singular_traceless = traceless_part(g_cov, g_inv, singular_tensor)
        trace_residual = trace_with_inverse(g_inv, residual_tensor)
        alg_norm = max(tensor_norm(algebraic_tensor), 1.0e-300)
        alg_tl_norm = max(tensor_norm(traceless_part(g_cov, g_inv, algebraic_tensor)), 1.0e-300)
        rows.append(
            {
                **base,
                "solved": 1.0,
                "ft_corrected": float(ft_corr),
                "coordinate_normal_speed_corrected": float(-ft_corr / max(np.hypot(fx, fz), 1.0e-300)),
                "normal_norm_corrected": normal_norm,
                "normal_norm_error": normal_norm_error,
                "singular_tensor_norm": tensor_norm(singular_tensor),
                "singular_traceless_norm": tensor_norm(singular_traceless),
                "tensor_residual_norm": tensor_norm(residual_tensor),
                "traceless_residual_norm": tensor_norm(residual_traceless),
                "tensor_residual_relative": tensor_norm(residual_tensor) / alg_norm,
                "traceless_residual_relative": tensor_norm(residual_traceless) / alg_tl_norm,
                "trace_residual_after_solve": trace_residual,
            }
        )
    return rows


def summarize(rows: list[dict[str, float]]) -> dict[str, object]:
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    return {
        "count": len(rows),
        "solved_count": len(solved),
        "unresolved_count": len(rows) - len(solved),
        "solved_fraction": float(len(solved) / max(len(rows), 1)),
        "target_normal_norm": scalar_stats([row["target_normal_norm"] for row in rows]),
        "quadratic_discriminant": scalar_stats([row["quadratic_discriminant"] for row in rows]),
        "trace_integral_mismatch": scalar_stats([row["trace_integral_mismatch"] for row in rows]),
        "algebraic_tensor_norm": scalar_stats([row["algebraic_tensor_norm"] for row in rows]),
        "algebraic_traceless_norm": scalar_stats([row["algebraic_traceless_norm"] for row in rows]),
        "rank1_tail_relative": scalar_stats([row["rank1_tail_relative"] for row in rows if "rank1_tail_relative" in row]),
        "rank1_metric_norm_error": scalar_stats(
            [row["rank1_metric_norm_error"] for row in rows if "rank1_metric_norm_error" in row]
        ),
        "rank1_negative_weight_relative": scalar_stats(
            [row["rank1_negative_weight_relative"] for row in rows if "rank1_negative_weight_relative" in row]
        ),
        "free_to_current_spatial_alignment": scalar_stats(
            [row["free_to_current_spatial_alignment"] for row in rows if "free_to_current_spatial_alignment" in row]
        ),
        "free_covector_target_norm": scalar_stats(
            [row["free_covector_target_norm"] for row in rows if "free_covector_target_norm" in row]
        ),
        "coordinate_normal_speed_corrected": scalar_stats(
            [row["coordinate_normal_speed_corrected"] for row in solved]
        ),
        "normal_norm_corrected": scalar_stats([row["normal_norm_corrected"] for row in solved]),
        "normal_norm_error": scalar_stats([row["normal_norm_error"] for row in solved]),
        "tensor_residual_norm": scalar_stats([row["tensor_residual_norm"] for row in solved]),
        "traceless_residual_norm": scalar_stats([row["traceless_residual_norm"] for row in solved]),
        "tensor_residual_relative": scalar_stats([row["tensor_residual_relative"] for row in solved]),
        "traceless_residual_relative": scalar_stats([row["traceless_residual_relative"] for row in solved]),
        "trace_residual_after_solve": scalar_stats([row["trace_residual_after_solve"] for row in solved]),
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
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    unsolved = [row for row in rows if row.get("solved", 0.0) <= 0.5]
    colors = ["#d95f02" if row["level"] > 0.0 else "#1b9e77" for row in rows]

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.8), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))
    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_lines(axes[0, 0], rows, colors=colors, linewidth=0.85)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0, 0].set_title("linear rho + raw_y=+/-1")
    fig.colorbar(im0, ax=axes[0, 0])

    add_lines(axes[0, 1], solved, colors="black", linewidth=0.85)
    add_lines(axes[0, 1], unsolved, colors="0.35", linewidth=0.7, linestyle="dotted")
    axes[0, 1].set_title("trace solvability; dotted unresolved")

    rank_rows = [row for row in rows if "rank1_tail_relative" in row]
    if rank_rows:
        rank_vals = np.array([row["rank1_tail_relative"] for row in rank_rows], dtype=float)
        rank_vmax = max(float(np.percentile(rank_vals[np.isfinite(rank_vals)], 95.0)), 1.0e-12)
        lc_rank = add_lines(axes[0, 2], rank_rows, values=rank_vals, vmin=0.0, vmax=rank_vmax, linewidth=1.0)
        fig.colorbar(lc_rank, ax=axes[0, 2])
    axes[0, 2].set_title("free-covector rank-one mismatch")

    if solved:
        rel = np.array([row["tensor_residual_relative"] for row in solved], dtype=float)
        rel_vmax = max(float(np.percentile(rel[np.isfinite(rel)], 95.0)), 1.0e-12)
        lc2 = add_lines(axes[1, 0], solved, values=rel, vmin=0.0, vmax=rel_vmax, linewidth=1.05)
        add_lines(axes[1, 0], unsolved, colors="0.35", linewidth=0.7, linestyle="dotted")
        fig.colorbar(lc2, ax=axes[1, 0])

        tl = np.array([row["traceless_residual_relative"] for row in solved], dtype=float)
        tl_vmax = max(float(np.percentile(tl[np.isfinite(tl)], 95.0)), 1.0e-12)
        lc3 = add_lines(axes[1, 1], solved, values=tl, vmin=0.0, vmax=tl_vmax, linewidth=1.05)
        add_lines(axes[1, 1], unsolved, colors="0.35", linewidth=0.7, linestyle="dotted")
        fig.colorbar(lc3, ax=axes[1, 1])

        speeds = np.array([row["coordinate_normal_speed_corrected"] for row in solved], dtype=float)
        speed_vmax = max(float(np.percentile(np.abs(speeds[np.isfinite(speeds)]), 95.0)), 1.0e-12)
        lc4 = add_lines(axes[1, 2], solved, values=speeds, vmin=-speed_vmax, vmax=speed_vmax, linewidth=1.05, cmap="coolwarm")
        add_lines(axes[1, 2], unsolved, colors="0.35", linewidth=0.7, linestyle="dotted")
        fig.colorbar(lc4, ax=axes[1, 2])
    axes[1, 0].set_title("full tensor residual / algebraic tensor")
    axes[1, 1].set_title("traceless residual / traceless algebraic")
    axes[1, 2].set_title("trace-law normal speed")

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
    levels = parse_positive_levels(args.levels)
    rows: list[dict[str, float]] = []
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        rows.extend(
            tensor_jump_rows(
                raw_y=real_array(data["raw_y"]),
                raw_y_m=real_array(data["raw_y_m"]),
                raw_y_p=real_array(data["raw_y_p"]),
                metric_cov=real_array(data["metric_cov"]),
                metric_inv=real_array(data["metric_inv"]),
                ricci=real_array(data["ricci"]),
                rho=real_array(data["rho"]),
                rhs_tensor=real_array(data["rhs_tensor"]),
                x=x,
                z=z,
                level=level,
                ell=args.ell,
                probe_dt=args.probe_dt,
                half_width_y=args.half_width_y,
                samples=args.samples,
            )
        )

    fig_path = out / f"d_tensor_interface_jump_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
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
        },
        "definitions": {
            "tensor_jump_law": "Leading thin-layer tensor condition: Delta(phi_q)*(-F_a F_b + g_ab g^cd F_c F_d) + Integral(algebraic_ab dq)=0.",
            "algebraic_tensor_ansatz": "Ricci traceless part is frozen at the interface; the pure trace part is shifted so that R(q)=q/ell^2 through the layer. This is a first tensor diagnostic, not the final tensor matching closure.",
            "trace_solvability": "No-real-F_t segments cannot satisfy even the trace of the tensor jump law for the current fixed spatial interface normal.",
            "rank1_tail_relative": "Free-covector algebraic compatibility check. Zero means the target tensor can be written as -F_a F_b + g_ab F^2 for some covector F_a, before imposing the current interface normal.",
            "free_to_current_spatial_alignment": "Absolute Euclidean alignment between the best free covector spatial direction and the current raw_y interface spatial normal. One means same direction up to sign; zero means orthogonal.",
            "tensor_residual_relative": "Frobenius norm of the tensor jump residual after enforcing the trace speed, divided by the algebraic tensor norm.",
            "traceless_residual_relative": "Same diagnostic after removing the metric trace; this measures the genuinely tensorial mismatch left by trace-only matching.",
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
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=61)
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
