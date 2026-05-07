from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from diagnose_d_direct_tensor_interface import direct_tensor_covector, direct_tensor_rows
from diagnose_d_direct_tensor_least_squares import gauss_newton_refine, gauss_newton_refine_multistart
from diagnose_d_spacetime_interface_jump import d_phi_d_raw_y
from diagnose_d_trace_interface_weak_form import find_level_segments, parse_positive_levels, real_array
from diagnose_d_tensor_interface_jump import algebraic_tensor_integral, reference_tensor_data, sample_matrix, tensor_norm
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


@dataclass
class TargetField:
    coords: np.ndarray
    tangents: np.ndarray
    residuals: np.ndarray
    speeds: np.ndarray
    levels: np.ndarray


def least_squares_tensor_rows(
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
    tensor_rel_tol: float,
    max_iter: int,
    multistart: bool = False,
    max_seeds: int = 24,
) -> list[dict[str, float]]:
    """Solve the direct tensor condition by minimizing the full tensor residual.

    The older rank-one projection is still used only to get an initial real
    covector. Acceptance here is based on the full tensor residual after
    local least-squares refinement, not on the trace equation.
    """
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
        algebraic_norm = max(tensor_norm(algebraic_tensor), 1.0e-300)
        initial = direct_tensor_covector(
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
            "candidate": float(initial.get("candidate", 0.0)),
            "delta_phi_prime": delta_phi_prime,
            "trace_algebraic_from_tensor": trace_from_tensor,
            "trace_algebraic_direct": trace_direct,
            "trace_integral_mismatch": trace_from_tensor - trace_direct,
            "algebraic_tensor_norm": tensor_norm(algebraic_tensor),
            "initial_direct_tensor_residual_relative": float(initial.get("direct_tensor_residual_relative", np.nan)),
            "rank1_tail_relative": float(initial.get("rank1_tail_relative", np.nan)),
            "rank1_negative_relative": float(initial.get("rank1_negative_relative", np.nan)),
            "solver": 1.0,
        }
        covector = initial.get("covector")
        if not isinstance(covector, np.ndarray) and not multistart:
            base["direct_tensor_solved"] = 0.0
            rows.append(base)
            continue

        if multistart:
            refined, refined_norm, iterations, seed_count = gauss_newton_refine_multistart(
                covector if isinstance(covector, np.ndarray) else None,
                g_cov,
                g_inv,
                algebraic_tensor,
                delta_phi_prime,
                max_iter=max_iter,
                current_spatial_normal=current_spatial,
                max_seeds=max_seeds,
            )
            base["ls_seed_count"] = float(seed_count)
            if refined is None:
                base["direct_tensor_solved"] = 0.0
                rows.append(base)
                continue
        else:
            refined, refined_norm, iterations = gauss_newton_refine(
                covector,
                g_cov,
                g_inv,
                algebraic_tensor,
                delta_phi_prime,
                max_iter=max_iter,
            )
            base["ls_seed_count"] = 1.0
        residual_relative = float(refined_norm / algebraic_norm)
        spatial_norm = float(np.linalg.norm(refined[1:3]))
        alignment = np.nan
        cur_norm = float(np.linalg.norm(current_spatial))
        if spatial_norm > 1.0e-300 and cur_norm > 1.0e-300:
            alignment = float(abs(np.dot(refined[1:3], current_spatial)) / (spatial_norm * cur_norm))
        base.update(
            {
                "direct_tensor_solved": float(residual_relative <= tensor_rel_tol),
                "direct_tensor_residual_relative": residual_relative,
                "ls_iterations": float(iterations),
                "ls_improvement_factor": float(
                    base["initial_direct_tensor_residual_relative"] / max(residual_relative, 1.0e-300)
                )
                if np.isfinite(base["initial_direct_tensor_residual_relative"])
                else np.nan,
                "direct_covector_norm": float(refined @ g_inv @ refined),
                "direct_speed": float(-refined[0] / spatial_norm) if spatial_norm > 1.0e-300 else np.nan,
                "direct_spatial_norm": spatial_norm,
                "spatial_alignment_with_current": alignment,
                "target_covector_t": float(refined[0]),
                "target_covector_x": float(refined[1]),
                "target_covector_z": float(refined[2]),
            }
        )
        rows.append(base)
    return rows


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def curve_segments(curves: list[dict[str, object]]) -> list[list[tuple[float, float]]]:
    segments: list[list[tuple[float, float]]] = []
    for curve in curves:
        points = np.asarray(curve["points"], dtype=float)
        if points.shape[0] < 2:
            continue
        for i in range(points.shape[0] - 1):
            segments.append([(float(points[i, 0]), float(points[i, 1])), (float(points[i + 1, 0]), float(points[i + 1, 1]))])
    return segments


def add_collection(
    ax: plt.Axes,
    segments: list[list[tuple[float, float]]],
    *,
    colors: str | list[str] = "black",
    values: np.ndarray | None = None,
    cmap: str = "viridis",
    linewidth: float = 0.8,
    linestyle: str = "solid",
    vmin: float | None = None,
    vmax: float | None = None,
) -> LineCollection | None:
    if not segments:
        return None
    if values is None:
        collection = LineCollection(segments, colors=colors, linewidths=linewidth, linestyles=linestyle)
    else:
        collection = LineCollection(segments, cmap=cmap, linewidths=linewidth, linestyles=linestyle)
        collection.set_array(values)
        if vmin is not None and vmax is not None:
            collection.set_clim(vmin, vmax)
    ax.add_collection(collection)
    return collection


def target_field_from_rows(rows: list[dict[str, float]]) -> TargetField:
    accepted = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    coords: list[list[float]] = []
    tangents: list[list[float]] = []
    residuals: list[float] = []
    speeds: list[float] = []
    levels: list[float] = []
    for row in accepted:
        normal = np.array([row.get("target_covector_x", np.nan), row.get("target_covector_z", np.nan)], dtype=float)
        norm = float(np.linalg.norm(normal))
        if not np.isfinite(norm) or norm <= 1.0e-300:
            continue
        tangent = np.array([-normal[1], normal[0]], dtype=float) / norm
        coords.append([row["x"], row["z"]])
        tangents.append([float(tangent[0]), float(tangent[1])])
        residuals.append(row["direct_tensor_residual_relative"])
        speeds.append(row["direct_speed"])
        levels.append(row["level"])
    return TargetField(
        coords=np.asarray(coords, dtype=float),
        tangents=np.asarray(tangents, dtype=float),
        residuals=np.asarray(residuals, dtype=float),
        speeds=np.asarray(speeds, dtype=float),
        levels=np.asarray(levels, dtype=float),
    )


def query_tangent(
    field: TargetField,
    point: np.ndarray,
    prev_dir: np.ndarray | None,
    *,
    level: float,
    radius: float,
    k_nearest: int,
) -> tuple[np.ndarray, float, int] | None:
    if field.coords.size == 0:
        return None
    same_level = np.abs(field.levels - level) <= 1.0e-12
    if not np.any(same_level):
        return None
    coords = field.coords[same_level]
    tangents = field.tangents[same_level]
    residuals = field.residuals[same_level]
    d2 = np.sum((coords - point[None, :]) ** 2, axis=1)
    order = np.argsort(d2)
    order = order[: max(1, min(k_nearest, order.size))]
    within = order[d2[order] <= radius * radius]
    if within.size == 0:
        return None

    base = tangents[within[0]]
    reference = base if prev_dir is None or np.linalg.norm(prev_dir) <= 1.0e-300 else prev_dir
    sigma2 = max((0.55 * radius) ** 2, 1.0e-300)
    weighted = np.zeros(2, dtype=float)
    for idx in within:
        tangent = np.array(tangents[idx], dtype=float)
        if float(np.dot(tangent, reference)) < 0.0:
            tangent = -tangent
        # Lower residual tensor matches get more influence, but no failed
        # segment is ever promoted into the field.
        quality = 1.0 / max(float(residuals[idx]), 1.0e-4)
        weight = float(np.exp(-d2[idx] / sigma2) * quality)
        weighted += weight * tangent
    norm = float(np.linalg.norm(weighted))
    if not np.isfinite(norm) or norm <= 1.0e-300:
        return None
    return weighted / norm, float(np.sqrt(d2[within[0]])), int(within.size)


def integrate_branch(
    field: TargetField,
    seed: np.ndarray,
    init_tangent: np.ndarray,
    *,
    level: float,
    ds: float,
    max_steps: int,
    radius: float,
    k_nearest: int,
    bounds: tuple[float, float, float, float],
) -> tuple[np.ndarray, dict[str, float]]:
    xmin, xmax, zmin, zmax = bounds

    def one_side(sign: float) -> list[np.ndarray]:
        points: list[np.ndarray] = []
        point = np.array(seed, dtype=float)
        prev = sign * np.array(init_tangent, dtype=float)
        for _ in range(max_steps):
            first = query_tangent(field, point, prev, level=level, radius=radius, k_nearest=k_nearest)
            if first is None:
                break
            direction, _, _ = first
            if float(np.dot(direction, prev)) < 0.0:
                direction = -direction
            midpoint = point + 0.5 * ds * direction
            second = query_tangent(field, midpoint, direction, level=level, radius=radius, k_nearest=k_nearest)
            if second is not None:
                direction = second[0]
                if float(np.dot(direction, prev)) < 0.0:
                    direction = -direction
            new_point = point + ds * direction
            if not (xmin <= new_point[0] <= xmax and zmin <= new_point[1] <= zmax):
                break
            if points and np.linalg.norm(new_point - points[-1]) <= 1.0e-12:
                break
            points.append(new_point)
            point = new_point
            prev = direction
        return points

    backward = one_side(-1.0)
    forward = one_side(1.0)
    points = list(reversed(backward)) + [np.array(seed, dtype=float)] + forward
    curve = np.asarray(points, dtype=float)
    if curve.shape[0] >= 2:
        length = float(np.sum(np.linalg.norm(np.diff(curve, axis=0), axis=1)))
    else:
        length = 0.0
    return curve, {"point_count": float(curve.shape[0]), "length": length}


def reconstruct_curves(
    field: TargetField,
    *,
    ds: float,
    max_steps: int,
    radius: float,
    k_nearest: int,
    seed_separation: float,
    cover_radius: float,
    max_curves_per_level: int,
    bounds: tuple[float, float, float, float],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    curves: list[dict[str, object]] = []
    if field.coords.size == 0:
        return curves, {"accepted_coverage_fraction": 0.0}

    covered = np.zeros(field.coords.shape[0], dtype=bool)
    # Start with the best local tensor matches, then let each streamline cover
    # nearby accepted samples. This avoids seeding many duplicate curves.
    order = np.argsort(field.residuals)
    level_counts: dict[float, int] = {}
    for idx in order:
        level = float(field.levels[idx])
        if level_counts.get(level, 0) >= max_curves_per_level:
            continue
        if covered[idx]:
            continue
        seed = field.coords[idx]
        curve, stats = integrate_branch(
            field,
            seed,
            field.tangents[idx],
            level=level,
            ds=ds,
            max_steps=max_steps,
            radius=radius,
            k_nearest=k_nearest,
            bounds=bounds,
        )
        if curve.shape[0] < 3 or stats["length"] < seed_separation:
            continue
        curves.append(
            {
                "level": level,
                "seed_x": float(seed[0]),
                "seed_z": float(seed[1]),
                "seed_residual": float(field.residuals[idx]),
                "seed_speed": float(field.speeds[idx]),
                "point_count": int(curve.shape[0]),
                "length": float(stats["length"]),
                "points": curve.tolist(),
            }
        )
        level_counts[level] = level_counts.get(level, 0) + 1
        d2 = np.sum((field.coords[:, None, :] - curve[None, :, :]) ** 2, axis=2)
        covered |= np.min(d2, axis=1) <= cover_radius * cover_radius

    nearest_distances: list[float] = []
    if curves:
        all_curve_points = np.vstack([np.asarray(curve["points"], dtype=float) for curve in curves])
        d2 = np.sum((field.coords[:, None, :] - all_curve_points[None, :, :]) ** 2, axis=2)
        nearest_distances = np.sqrt(np.min(d2, axis=1)).tolist()

    total_length = float(sum(float(curve["length"]) for curve in curves))
    summary = {
        "curve_count": len(curves),
        "curve_count_by_level": {str(k): v for k, v in sorted(level_counts.items())},
        "total_length": total_length,
        "accepted_coverage_fraction": float(np.mean(covered)) if covered.size else 0.0,
        "accepted_nearest_curve_distance": scalar_stats(nearest_distances),
        "curve_length": scalar_stats([float(curve["length"]) for curve in curves]),
        "curve_point_count": scalar_stats([float(curve["point_count"]) for curve in curves]),
    }
    return curves, summary


def summarize_rows(rows: list[dict[str, float]], field: TargetField, curves_summary: dict[str, object]) -> dict[str, object]:
    candidates = [row for row in rows if row.get("candidate", 0.0) > 0.5]
    accepted = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    return {
        "row_count": len(rows),
        "candidate_count": len(candidates),
        "candidate_fraction": float(len(candidates) / max(len(rows), 1)),
        "accepted_count": len(accepted),
        "accepted_fraction": float(len(accepted) / max(len(rows), 1)),
        "target_field_count": int(field.coords.shape[0]),
        "accepted_direct_tensor_residual_relative": scalar_stats(
            [row["direct_tensor_residual_relative"] for row in accepted]
        ),
        "accepted_target_speed": scalar_stats([row["direct_speed"] for row in accepted]),
        "accepted_target_current_alignment": scalar_stats(
            [row["spatial_alignment_with_current"] for row in accepted]
        ),
        "reconstruction": curves_summary,
    }


def render(
    out_path: Path,
    *,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    rows: list[dict[str, float]],
    curves: list[dict[str, object]],
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    accepted = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    rejected = [row for row in rows if row.get("direct_tensor_solved", 0.0) <= 0.5]
    current_colors = ["#1b9e77" if row["level"] < 0.0 else "#d95f02" for row in rows]
    curve_colors = ["#0077b6" if float(curve["level"]) < 0.0 else "#c1121f" for curve in curves]

    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.8), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))

    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_collection(axes[0, 0], row_segments(rows), colors=current_colors, linewidth=0.7)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0, 0].set_title("linear rho + current raw_y=+/-1")
    fig.colorbar(im0, ax=axes[0, 0])

    add_collection(axes[0, 1], row_segments(rejected), colors="0.7", linewidth=0.5, linestyle="dotted")
    add_collection(axes[0, 1], row_segments(accepted), colors="black", linewidth=0.9)
    axes[0, 1].set_title("direct tensor accepted/rejected")

    add_collection(axes[0, 2], row_segments(rows), colors="0.82", linewidth=0.45)
    add_collection(axes[0, 2], curve_segments(curves), colors=curve_colors, linewidth=1.25)
    axes[0, 2].set_title("reconstructed target curves")

    add_collection(axes[1, 0], row_segments(accepted), colors="0.72", linewidth=0.65)
    add_collection(axes[1, 0], curve_segments(curves), colors=curve_colors, linewidth=1.25)
    axes[1, 0].set_title("accepted current segments vs target curves")

    if accepted:
        vals = np.asarray([row["spatial_alignment_with_current"] for row in accepted], dtype=float)
        lc = add_collection(
            axes[1, 1],
            row_segments(accepted),
            values=vals,
            cmap="viridis",
            linewidth=1.0,
            vmin=0.0,
            vmax=1.0,
        )
        fig.colorbar(lc, ax=axes[1, 1])
    axes[1, 1].set_title("target-current normal alignment")

    if accepted:
        vals = np.asarray([row["direct_speed"] for row in accepted], dtype=float)
        finite = vals[np.isfinite(vals)]
        vmax = max(float(np.percentile(np.abs(finite), 95.0)), 1.0e-12) if finite.size else 1.0
        lc = add_collection(
            axes[1, 2],
            row_segments(accepted),
            values=vals,
            cmap="coolwarm",
            linewidth=1.0,
            vmin=-vmax,
            vmax=vmax,
        )
        fig.colorbar(lc, ax=axes[1, 2])
    axes[1, 2].set_title("target speed from full tensor condition")

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
        if args.solver == "least_squares":
            rows.extend(
                least_squares_tensor_rows(
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
                    tensor_rel_tol=args.tensor_rel_tol,
                    max_iter=args.max_iter,
                    multistart=args.multistart,
                    max_seeds=args.max_seeds,
                )
            )
        else:
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

    dx = float(abs(x[1] - x[0]))
    dz = float(abs(z[1] - z[0]))
    grid_step = float(np.hypot(dx, dz))
    ds = args.ds if args.ds > 0.0 else 0.45 * grid_step
    radius = args.field_radius if args.field_radius > 0.0 else 2.5 * grid_step
    seed_separation = args.seed_separation if args.seed_separation > 0.0 else 1.25 * grid_step
    cover_radius = args.cover_radius if args.cover_radius > 0.0 else 1.0 * grid_step

    field = target_field_from_rows(rows)
    curves, curves_summary = reconstruct_curves(
        field,
        ds=ds,
        max_steps=args.max_steps,
        radius=radius,
        k_nearest=args.k_nearest,
        seed_separation=seed_separation,
        cover_radius=cover_radius,
        max_curves_per_level=args.max_curves_per_level,
        bounds=(float(x[0]), float(x[-1]), float(z[0]), float(z[-1])),
    )

    fig_path = out / f"d_tensor_reconstructed_interface_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render(
        fig_path,
        x=x,
        z=z,
        rho=real_array(data["rho"]),
        rows=rows,
        curves=curves,
        rho_min_frac=args.rho_min_frac,
    )

    rows_path = out / "direct_tensor_interface_rows.jsonl"
    with rows_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
    curves_path = out / "tensor_reconstructed_curves.json"
    curves_path.write_text(json.dumps(curves, ensure_ascii=False, indent=2), encoding="utf-8")

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
            "rank1_tol": args.rank1_tol,
            "tensor_rel_tol": args.tensor_rel_tol,
            "solver": args.solver,
            "max_iter": args.max_iter,
            "multistart": args.multistart,
            "max_seeds": args.max_seeds,
            "ds": ds,
            "field_radius": radius,
            "seed_separation": seed_separation,
            "cover_radius": cover_radius,
            "k_nearest": args.k_nearest,
            "max_steps": args.max_steps,
            "max_curves_per_level": args.max_curves_per_level,
        },
        "definitions": {
            "scope": "Prototype reconstruction of interface curves from the full direct-tensor target covector field. It does not solve trace speed and does not promote unresolved segments.",
            "target_tangent": "The local curve tangent is perpendicular to the direct-tensor target spatial covector (F_x,F_z).",
            "least_squares_solver": "When --solver least_squares is used, F_a is locally refined by minimizing the full tensor residual; this is a numerical solve of the same tensor condition, not a new physical rule.",
            "multistart_solver": "When --multistart is used, several real covector seeds are tried for the same least-squares residual. It is a solver diagnostic, not an extra closure condition.",
            "accepted_coverage_fraction": "Fraction of accepted direct-tensor samples within cover_radius of a reconstructed curve.",
            "tensor_unresolved": "Rows not accepted by the direct tensor rank-one/residual test remain unresolved and are not filled by trace speed, damping, or clipping.",
        },
        "summary": summarize_rows(rows, field, curves_summary),
        "files": {
            "figure": str(fig_path.resolve()),
            "rows_jsonl": str(rows_path.resolve()),
            "curves_json": str(curves_path.resolve()),
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
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    parser.add_argument("--rank1-tol", type=float, default=0.1)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument("--solver", choices=["rank_projection", "least_squares"], default="rank_projection")
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--multistart", action="store_true")
    parser.add_argument("--max-seeds", type=int, default=24)
    parser.add_argument("--ds", type=float, default=0.0)
    parser.add_argument("--field-radius", type=float, default=0.0)
    parser.add_argument("--seed-separation", type=float, default=0.0)
    parser.add_argument("--cover-radius", type=float, default=0.0)
    parser.add_argument("--k-nearest", type=int, default=12)
    parser.add_argument("--max-steps", type=int, default=80)
    parser.add_argument("--max-curves-per-level", type=int, default=36)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
