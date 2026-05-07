from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from diagnose_d_direct_tensor_interface import direct_tensor_covector
from diagnose_d_spacetime_interface_jump import d_phi_d_raw_y
from diagnose_d_trace_interface_weak_form import find_level_segments, parse_positive_levels, real_array
from diagnose_d_tensor_interface_jump import (
    algebraic_tensor_integral,
    reference_tensor_data,
    sample_matrix,
    tensor_norm,
)
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


SYMMETRIC_COMPONENTS = [(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0), (0, 1, np.sqrt(2.0)), (0, 2, np.sqrt(2.0)), (1, 2, np.sqrt(2.0))]


def classify(row: dict[str, float], *, rank1_tol: float, tensor_rel_tol: float) -> str:
    if row.get("direct_tensor_solved", 0.0) > 0.5:
        return "accepted"
    if row.get("candidate", 0.0) <= 0.5:
        return "no_real_covector"
    failures: list[str] = []
    if row.get("rank1_tail_relative", np.inf) > rank1_tol:
        failures.append("rank_tail")
    if row.get("rank1_negative_relative", np.inf) > rank1_tol:
        failures.append("negative")
    if row.get("direct_tensor_residual_relative", np.inf) > tensor_rel_tol:
        failures.append("tensor_residual")
    return "+".join(failures) if failures else "unclassified_numeric"


def residual_tensor(covector: np.ndarray, metric_cov: np.ndarray, metric_inv: np.ndarray, algebraic_tensor: np.ndarray, delta_phi_prime: float) -> np.ndarray:
    normal_norm = float(covector @ metric_inv @ covector)
    return delta_phi_prime * (-np.outer(covector, covector) + metric_cov * normal_norm) + algebraic_tensor


def residual_vector_and_jacobian(
    covector: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    algebraic_tensor: np.ndarray,
    delta_phi_prime: float,
) -> tuple[np.ndarray, np.ndarray]:
    residual = residual_tensor(covector, metric_cov, metric_inv, algebraic_tensor, delta_phi_prime)
    metric_f = metric_inv @ covector
    vec = np.zeros(6, dtype=float)
    jac = np.zeros((6, 3), dtype=float)
    for row_index, (a, b, factor) in enumerate(SYMMETRIC_COMPONENTS):
        vec[row_index] = factor * residual[a, b]
        for c in range(3):
            deriv = delta_phi_prime * (
                -(1.0 if a == c else 0.0) * covector[b]
                - covector[a] * (1.0 if b == c else 0.0)
                + 2.0 * metric_cov[a, b] * metric_f[c]
            )
            jac[row_index, c] = factor * deriv
    return vec, jac


def gauss_newton_refine(
    covector: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    algebraic_tensor: np.ndarray,
    delta_phi_prime: float,
    *,
    max_iter: int,
) -> tuple[np.ndarray, float, int]:
    current = np.asarray(covector, dtype=float).copy()
    best_vec, _ = residual_vector_and_jacobian(current, metric_cov, metric_inv, algebraic_tensor, delta_phi_prime)
    best_norm = float(np.linalg.norm(best_vec))
    for iteration in range(max_iter):
        vec, jac = residual_vector_and_jacobian(current, metric_cov, metric_inv, algebraic_tensor, delta_phi_prime)
        current_norm = float(np.linalg.norm(vec))
        if current_norm < best_norm:
            best_norm = current_norm
        scale = max(float(np.linalg.norm(jac)), 1.0)
        lhs = jac.T @ jac + (1.0e-12 * scale * scale) * np.eye(3)
        rhs = -jac.T @ vec
        try:
            step = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(lhs, rhs, rcond=1.0e-12)[0]
        if not np.all(np.isfinite(step)) or float(np.linalg.norm(step)) <= 1.0e-14 * max(float(np.linalg.norm(current)), 1.0):
            return current, current_norm, iteration
        accepted = False
        alpha = 1.0
        for _ in range(12):
            trial = current + alpha * step
            trial_vec, _ = residual_vector_and_jacobian(trial, metric_cov, metric_inv, algebraic_tensor, delta_phi_prime)
            trial_norm = float(np.linalg.norm(trial_vec))
            if np.isfinite(trial_norm) and trial_norm < current_norm:
                current = trial
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return current, current_norm, iteration
    final_vec, _ = residual_vector_and_jacobian(current, metric_cov, metric_inv, algebraic_tensor, delta_phi_prime)
    return current, float(np.linalg.norm(final_vec)), max_iter


def _align_with_current_spatial(covector: np.ndarray, current_spatial_normal: np.ndarray | None) -> np.ndarray:
    out = np.asarray(covector, dtype=float).copy()
    if current_spatial_normal is None:
        return out
    spatial = out[1:3]
    if np.linalg.norm(spatial) > 1.0e-300 and np.linalg.norm(current_spatial_normal) > 1.0e-300:
        if float(np.dot(spatial, current_spatial_normal)) < 0.0:
            out = -out
    return out


def multistart_covector_seeds(
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    algebraic_tensor: np.ndarray,
    delta_phi_prime: float,
    *,
    primary: np.ndarray | None = None,
    current_spatial_normal: np.ndarray | None = None,
    max_seeds: int = 24,
) -> list[np.ndarray]:
    """Generate real covector seeds for the same full tensor residual.

    This does not add a new physical rule. It only probes whether the local
    nonlinear least-squares solve was trapped by the rank-one projection seed.
    """
    seeds: list[np.ndarray] = []

    def add(seed: np.ndarray) -> None:
        seed = _align_with_current_spatial(seed, current_spatial_normal)
        if not np.all(np.isfinite(seed)):
            return
        norm = float(np.linalg.norm(seed))
        if norm <= 1.0e-300:
            return
        for old in seeds:
            if float(np.linalg.norm(seed - old)) <= 1.0e-8 * max(norm, float(np.linalg.norm(old)), 1.0):
                return
        seeds.append(seed)

    if primary is not None:
        primary = np.asarray(primary, dtype=float)
        for factor in (0.5, 1.0, 2.0, -1.0):
            add(factor * primary)

    if abs(delta_phi_prime) > 1.0e-300:
        target_h = -algebraic_tensor / delta_phi_prime
        target_norm = 0.5 * float(np.sum(metric_inv * target_h))
        k_cov = metric_cov * target_norm - target_h
        k_sym = 0.5 * (k_cov + k_cov.T)
        evals, evecs = np.linalg.eigh(k_sym)
        order = np.argsort(np.abs(evals))[::-1]
        for idx in order:
            scale = float(np.sqrt(max(abs(evals[idx]), 1.0e-300)))
            direction = np.asarray(evecs[:, idx], dtype=float)
            for factor in (0.5, 1.0, 2.0):
                add(factor * scale * direction)

    seed_scale = 1.0
    if seeds:
        seed_scale = float(np.median([np.linalg.norm(seed) for seed in seeds]))
    else:
        seed_scale = max(float(np.sqrt(np.linalg.norm(algebraic_tensor))), 1.0)

    add(np.array([seed_scale, 0.0, 0.0], dtype=float))
    if current_spatial_normal is not None and np.linalg.norm(current_spatial_normal) > 1.0e-300:
        spatial = np.asarray(current_spatial_normal, dtype=float) / float(np.linalg.norm(current_spatial_normal))
        add(np.array([0.0, seed_scale * spatial[0], seed_scale * spatial[1]], dtype=float))
        add(np.array([seed_scale, seed_scale * spatial[0], seed_scale * spatial[1]], dtype=float))
        add(np.array([-seed_scale, seed_scale * spatial[0], seed_scale * spatial[1]], dtype=float))

    return seeds[:max_seeds]


def gauss_newton_refine_multistart(
    primary_covector: np.ndarray | None,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    algebraic_tensor: np.ndarray,
    delta_phi_prime: float,
    *,
    max_iter: int,
    current_spatial_normal: np.ndarray | None = None,
    max_seeds: int = 24,
) -> tuple[np.ndarray | None, float, int, int]:
    seeds = multistart_covector_seeds(
        metric_cov,
        metric_inv,
        algebraic_tensor,
        delta_phi_prime,
        primary=primary_covector,
        current_spatial_normal=current_spatial_normal,
        max_seeds=max_seeds,
    )
    best_covector: np.ndarray | None = None
    best_norm = np.inf
    best_iterations = 0
    for seed in seeds:
        refined, residual_norm, iterations = gauss_newton_refine(
            seed,
            metric_cov,
            metric_inv,
            algebraic_tensor,
            delta_phi_prime,
            max_iter=max_iter,
        )
        if np.isfinite(residual_norm) and residual_norm < best_norm:
            best_covector = refined
            best_norm = float(residual_norm)
            best_iterations = int(iterations)
    if best_covector is None:
        return None, np.inf, 0, len(seeds)
    best_covector = _align_with_current_spatial(best_covector, current_spatial_normal)
    return best_covector, best_norm, best_iterations, len(seeds)


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def render(out_path: Path, rows: list[dict[str, float]], x: np.ndarray, z: np.ndarray, rho: np.ndarray) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.8), constrained_layout=True)
    im0 = axes[0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=max(float(np.percentile(rho, 99.5)), 1.0e-16))
    axes[0].set_title("linear rho")
    fig.colorbar(im0, ax=axes[0])

    for ax, key, title, vmax_percentile in [
        (axes[1], "initial_direct_tensor_residual_relative", "initial direct tensor residual", 95.0),
        (axes[2], "ls_direct_tensor_residual_relative", "after local least-squares", 95.0),
    ]:
        candidate_rows = [row for row in rows if np.isfinite(row.get(key, np.nan))]
        if candidate_rows:
            vals = np.asarray([row[key] for row in candidate_rows], dtype=float)
            vmax = max(float(np.percentile(vals[np.isfinite(vals)], vmax_percentile)), 1.0e-12)
            collection = LineCollection(row_segments(candidate_rows), cmap="magma", linewidths=0.95)
            collection.set_array(vals)
            collection.set_clim(0.0, vmax)
            ax.add_collection(collection)
            fig.colorbar(collection, ax=ax)
        ax.set_title(title)

    for ax in axes:
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
        ax.set_aspect("equal")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def summarize(rows: list[dict[str, float]], *, tensor_rel_tol: float) -> dict[str, object]:
    by_mode: dict[str, list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        by_mode[row["initial_mode"]].append(row)
    mode_summary: dict[str, object] = {}
    for mode, group in sorted(by_mode.items(), key=lambda item: item[0]):
        mode_summary[mode] = {
            "count": len(group),
            "ls_accepted_count": int(sum(row.get("ls_direct_tensor_residual_relative", np.inf) <= tensor_rel_tol for row in group)),
            "initial_direct_tensor_residual_relative": scalar_stats(
                [row.get("initial_direct_tensor_residual_relative", np.nan) for row in group]
            ),
            "ls_direct_tensor_residual_relative": scalar_stats(
                [row.get("ls_direct_tensor_residual_relative", np.nan) for row in group]
            ),
            "ls_improvement_factor": scalar_stats([row.get("ls_improvement_factor", np.nan) for row in group]),
            "ls_iterations": scalar_stats([row.get("ls_iterations", np.nan) for row in group]),
        }
    return {
        "count": len(rows),
        "initial_mode_counts": dict(Counter(row["initial_mode"] for row in rows)),
        "ls_accepted_count": int(sum(row.get("ls_direct_tensor_residual_relative", np.inf) <= tensor_rel_tol for row in rows)),
        "ls_accepted_fraction": float(
            np.mean([row.get("ls_direct_tensor_residual_relative", np.inf) <= tensor_rel_tol for row in rows])
        )
        if rows
        else 0.0,
        "mode_summary": mode_summary,
    }


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
    raw_y = real_array(data["raw_y"])
    rho = real_array(data["rho"])
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)

    rows: list[dict[str, float]] = []
    levels = parse_positive_levels(args.levels)
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        q = np.linspace(level - args.half_width_y, level + args.half_width_y, args.samples)
        delta_phi_prime = float(d_phi_d_raw_y(q)[-1] - d_phi_d_raw_y(q)[0])
        for seg in find_level_segments(raw_y, rho, x, z, level):
            xp = np.array([seg["x"]])
            zp = np.array([seg["z"]])
            current_spatial = np.array(
                [
                    float(bilinear_sample(raw_y_x, x, z, xp, zp)[0]),
                    float(bilinear_sample(raw_y_z, x, z, xp, zp)[0]),
                ],
                dtype=float,
            )
            g_cov = sample_matrix(real_array(data["metric_cov"]), x, z, xp, zp)
            g_inv = np.linalg.pinv(g_cov, rcond=1.0e-12, hermitian=True)
            algebraic_tensor, _, _ = algebraic_tensor_integral(
                q=q,
                level=level,
                ell=args.ell,
                metric_cov=g_cov,
                metric_inv=g_inv,
                ricci=sample_matrix(real_array(data["ricci"]), x, z, xp, zp),
                rhs_tensor=sample_matrix(real_array(data["rhs_tensor"]), x, z, xp, zp),
            )
            algebraic_norm = max(tensor_norm(algebraic_tensor), 1.0e-300)
            result = direct_tensor_covector(
                metric_cov=g_cov,
                metric_inv=g_inv,
                algebraic_tensor=algebraic_tensor,
                delta_phi_prime=delta_phi_prime,
                current_spatial_normal=current_spatial,
            )
            row: dict[str, float] = {
                "x0": float(seg["x0"]),
                "z0": float(seg["z0"]),
                "x1": float(seg["x1"]),
                "z1": float(seg["z1"]),
                "x": float(seg["x"]),
                "z": float(seg["z"]),
                "level": float(level),
                "candidate": float(result.get("candidate", 0.0)),
                "initial_direct_tensor_residual_relative": float(result.get("direct_tensor_residual_relative", np.nan)),
                "rank1_tail_relative": float(result.get("rank1_tail_relative", np.nan)),
                "rank1_negative_relative": float(result.get("rank1_negative_relative", np.nan)),
            }
            row["direct_tensor_solved"] = float(
                row["candidate"] > 0.5
                and row.get("rank1_tail_relative", np.inf) <= args.rank1_tol
                and row.get("rank1_negative_relative", np.inf) <= args.rank1_tol
                and row.get("initial_direct_tensor_residual_relative", np.inf) <= args.tensor_rel_tol
            )
            row["initial_mode"] = classify(row, rank1_tol=args.rank1_tol, tensor_rel_tol=args.tensor_rel_tol)
            covector = result.get("covector")
            if isinstance(covector, np.ndarray):
                refined, refined_norm, iterations = gauss_newton_refine(
                    covector,
                    g_cov,
                    g_inv,
                    algebraic_tensor,
                    delta_phi_prime,
                    max_iter=args.max_iter,
                )
                row["ls_direct_tensor_residual_relative"] = float(refined_norm / algebraic_norm)
                initial_rel = row["initial_direct_tensor_residual_relative"]
                row["ls_improvement_factor"] = float(initial_rel / max(row["ls_direct_tensor_residual_relative"], 1.0e-300))
                row["ls_iterations"] = float(iterations)
                row["ls_speed"] = float(-refined[0] / np.linalg.norm(refined[1:3])) if np.linalg.norm(refined[1:3]) > 1.0e-300 else np.nan
            rows.append(row)

    fig_path = out / f"d_direct_tensor_least_squares_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render(fig_path, rows, x, z, rho)
    rows_path = out / "direct_tensor_least_squares_rows.jsonl"
    with rows_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
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
            "max_iter": args.max_iter,
        },
        "definitions": {
            "least_squares_check": "Local Gauss-Newton minimization of the full tensor residual over the three covector components F_a, initialized from the rank-one direct-tensor covector.",
            "purpose": "Tests whether tensor_residual failures are an artifact of the initial rank-one projection rather than a genuine full tensor mismatch.",
        },
        "summary": summarize(rows, tensor_rel_tol=args.tensor_rel_tol),
        "files": {"figure": str(fig_path.resolve()), "rows_jsonl": str(rows_path.resolve())},
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
    parser.add_argument("--rank1-tol", type=float, default=0.1)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument("--max-iter", type=int, default=30)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
