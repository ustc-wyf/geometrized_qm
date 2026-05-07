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

from diagnose_d_direct_tensor_interface import direct_tensor_rows
from diagnose_d_trace_interface_weak_form import parse_positive_levels, real_array
from diagnose_d_tensor_interface_jump import reference_tensor_data, sample_matrix
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


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


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_lines(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    *,
    values: np.ndarray | None = None,
    colors: str | list[str] = "black",
    cmap: str = "viridis",
    linewidth: float = 0.85,
    vmin: float | None = None,
    vmax: float | None = None,
) -> LineCollection | None:
    if not rows:
        return None
    if values is None:
        collection = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth)
    else:
        collection = LineCollection(row_segments(rows), cmap=cmap, linewidths=linewidth)
        collection.set_array(values)
        if vmin is not None and vmax is not None:
            collection.set_clim(vmin, vmax)
    ax.add_collection(collection)
    return collection


def sample_scalar(field: np.ndarray, x: np.ndarray, z: np.ndarray, row: dict[str, float]) -> float:
    return float(bilinear_sample(field, x, z, np.array([row["x"]]), np.array([row["z"]]))[0])


def add_features(
    rows: list[dict[str, float]],
    *,
    data: dict[str, np.ndarray],
    x: np.ndarray,
    z: np.ndarray,
) -> list[dict[str, float]]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_y = real_array(data["raw_y"])
    rho = real_array(data["rho"])
    raw_y_x, raw_y_z = np.gradient(raw_y, dx, dz, edge_order=2)
    raw_y_xx = np.gradient(raw_y_x, dx, axis=0, edge_order=2)
    raw_y_zz = np.gradient(raw_y_z, dz, axis=1, edge_order=2)
    rho_max = max(float(np.max(rho)), 1.0e-300)

    enriched: list[dict[str, float]] = []
    for row in rows:
        xp = np.array([row["x"]])
        zp = np.array([row["z"]])
        g_cov = sample_matrix(real_array(data["metric_cov"]), x, z, xp, zp)
        g_inv = np.linalg.pinv(g_cov, rcond=1.0e-12, hermitian=True)
        eig_abs = np.sort(np.abs(np.linalg.eigvalsh(0.5 * (g_cov + g_cov.T))))
        cond = float(eig_abs[-1] / max(eig_abs[0], 1.0e-300))
        det_abs = float(abs(np.linalg.det(g_cov)))
        grad_x = sample_scalar(raw_y_x, x, z, row)
        grad_z = sample_scalar(raw_y_z, x, z, row)
        grad_norm = float(np.hypot(grad_x, grad_z))
        lap_abs = abs(sample_scalar(raw_y_xx + raw_y_zz, x, z, row))
        rho_here = sample_scalar(rho, x, z, row)
        direct_target_norm = row.get("direct_target_norm", np.nan)
        direct_norm_error = row.get("direct_norm_error", np.nan)
        target_norm_relative_error = (
            float(abs(direct_norm_error) / max(abs(direct_target_norm), 1.0e-300))
            if np.isfinite(direct_norm_error) and np.isfinite(direct_target_norm)
            else np.nan
        )
        enriched.append(
            {
                **row,
                "rho_here": rho_here,
                "rho_relative": float(rho_here / rho_max),
                "log10_rho_relative": float(np.log10(max(rho_here / rho_max, 1.0e-300))),
                "raw_y_grad_norm": grad_norm,
                "raw_y_laplacian_abs": float(lap_abs),
                "metric_condition_abs_eig": cond,
                "metric_det_abs": det_abs,
                "target_norm_relative_error": target_norm_relative_error,
                "ginv_fro_norm": float(np.sqrt(np.sum(g_inv * g_inv))),
            }
        )
    return enriched


def grouped_stats(rows: list[dict[str, float]], modes: list[str], keys: list[str]) -> dict[str, object]:
    grouped: dict[str, list[dict[str, float]]] = defaultdict(list)
    for row, mode in zip(rows, modes):
        grouped[mode].append(row)
    out: dict[str, object] = {}
    for mode in sorted(grouped, key=lambda m: (m != "accepted", m)):
        group = grouped[mode]
        out[mode] = {"count": len(group), "fraction": float(len(group) / max(len(rows), 1))}
        for key in keys:
            out[mode][key] = scalar_stats([row.get(key, np.nan) for row in group])
    return out


def render(
    out_path: Path,
    *,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    rows: list[dict[str, float]],
    modes: list[str],
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.8), constrained_layout=True)

    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=max(float(np.percentile(rho, 99.5)), 1.0e-16))
    add_lines(axes[0, 0], rows, colors=["#1b9e77" if row["level"] < 0.0 else "#d95f02" for row in rows], linewidth=0.65)
    axes[0, 0].set_title("linear rho + raw_y=+/-1")
    fig.colorbar(im0, ax=axes[0, 0])

    feature_specs = [
        ("log10_rho_relative", "log10(rho/rho_max)", "viridis", None, 0.0),
        ("raw_y_grad_norm", "|grad raw_y|", "magma", 0.0, None),
        ("metric_condition_abs_eig", "metric abs-eigen condition", "magma", 1.0, None),
        ("target_norm_relative_error", "|F^2-target|/|target|", "magma", 0.0, None),
        ("direct_tensor_residual_relative", "direct tensor residual", "magma", 0.0, None),
    ]
    for ax, (key, title, cmap, vmin, vmax) in zip(axes.ravel()[1:], feature_specs):
        vals = np.asarray([row.get(key, np.nan) for row in rows], dtype=float)
        finite = vals[np.isfinite(vals)]
        if finite.size == 0:
            ax.set_title(title)
            continue
        if vmax is None:
            vmax_plot = max(float(np.percentile(finite, 95.0)), 1.0e-12)
        else:
            vmax_plot = vmax
        vmin_plot = float(np.percentile(finite, 5.0)) if vmin is None else vmin
        lc = add_lines(ax, rows, values=vals, cmap=cmap, linewidth=0.9, vmin=vmin_plot, vmax=vmax_plot)
        fig.colorbar(lc, ax=ax)
        ax.set_title(title)

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

    rows = add_features(rows, data=data, x=x, z=z)
    modes = [classify(row, rank1_tol=args.rank1_tol, tensor_rel_tol=args.tensor_rel_tol) for row in rows]
    keys = [
        "rho_relative",
        "log10_rho_relative",
        "raw_y_grad_norm",
        "raw_y_laplacian_abs",
        "metric_condition_abs_eig",
        "metric_det_abs",
        "ginv_fro_norm",
        "algebraic_tensor_norm",
        "direct_target_norm",
        "direct_norm_error",
        "target_norm_relative_error",
        "rank1_tail_relative",
        "rank1_negative_relative",
        "direct_tensor_residual_relative",
        "spatial_alignment_with_current",
        "direct_speed",
    ]
    counts = Counter(modes)
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
        },
        "definitions": {
            "metric_condition_abs_eig": "Condition number from absolute eigenvalues of the local interpolated covariant metric.",
            "target_norm_relative_error": "Relative mismatch between the extracted covector norm F^2 and the tensor-implied target norm.",
            "raw_y_grad_norm": "Euclidean spatial gradient magnitude of raw_y=ell^2 R_tilde on the reference slice.",
        },
        "summary": {
            "count": len(rows),
            "mode_counts": dict(counts),
            "mode_fractions": {key: float(value / max(len(rows), 1)) for key, value in counts.items()},
            "mode_feature_stats": grouped_stats(rows, modes, keys),
        },
    }

    fig_path = out / f"d_direct_tensor_failure_features_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render(fig_path, x=x, z=z, rho=real_array(data["rho"]), rows=rows, modes=modes)
    rows_path = out / "direct_tensor_failure_features_rows.jsonl"
    with rows_path.open("w", encoding="utf-8") as handle:
        for row, mode in zip(rows, modes):
            handle.write(json.dumps({"mode": mode, **row}, ensure_ascii=False, sort_keys=True) + "\n")
    summary["files"] = {"figure": str(fig_path.resolve()), "rows_jsonl": str(rows_path.resolve())}
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
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
