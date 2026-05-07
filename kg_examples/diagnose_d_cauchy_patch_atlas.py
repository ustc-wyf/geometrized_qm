from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import build_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


def tau_tag(tau: float) -> str:
    if abs(tau) < 5.0e-13:
        return "tau0"
    sign = "p" if tau > 0.0 else "m"
    return "tau" + sign + f"{abs(tau):g}".replace(".", "p")


def scalar_contract(metric_inv: np.ndarray, a_cov: np.ndarray, b_cov: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...a,...b->...", metric_inv, a_cov, b_cov, optimize=True)


def raise_covector(metric_inv: np.ndarray, cov: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...b->...a", metric_inv, cov, optimize=True)


def weighted_stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {
            "count": 0,
            "mean": 0.0,
            "p50": 0.0,
            "p95": 0.0,
            "p99": 0.0,
            "min": 0.0,
            "max": 0.0,
        }
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "p99": float(np.percentile(vals, 99.0)),
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
    }
    if weights is not None:
        all_w = np.asarray(weights, dtype=float)
        w = np.maximum(all_w[mask], 0.0)
        wsum = float(np.sum(w))
        if wsum > 0.0:
            order = np.argsort(vals)
            sorted_vals = vals[order]
            sorted_w = w[order]
            cdf = np.cumsum(sorted_w) / wsum
            out["weighted_mean"] = float(np.sum(vals * w) / wsum)
            out["weighted_p50"] = float(sorted_vals[np.searchsorted(cdf, 0.50, side="left")])
            out["weighted_p95"] = float(sorted_vals[np.searchsorted(cdf, 0.95, side="left")])
    return out


def nullspace_of_row(row: np.ndarray, rcond: float = 1.0e-12) -> np.ndarray:
    mat = np.asarray(row, dtype=float)[None, :]
    _, s, vh = np.linalg.svd(mat, full_matrices=True)
    tol = rcond * (float(s[0]) if s.size else 0.0)
    rank = int(np.count_nonzero(s > tol))
    return vh[rank:].T


def lambda_time_symbol_after_trace(case) -> tuple[np.ndarray, np.ndarray]:
    """Time-principal matrix for lambda after trace0 elimination.

    In the 2+1 txz reduction the four basis tensors are
    g, uu, rr, ur.  Trace0 removes one algebraic direction.  The time
    part of div C=0 should have rank three on the remaining subspace.
    """

    metric_inv = case.metric_inv
    u_cov = case.u_cov
    r_cov = case.r_cov
    u_up = raise_covector(metric_inv, u_cov)
    r_up = raise_covector(metric_inv, r_cov)
    xi_cov = np.array([1.0, 0.0, 0.0])
    u_dot_xi = u_up[..., 0]
    r_dot_xi = r_up[..., 0]
    rows = np.zeros(u_cov.shape[:-1] + (3, 4), dtype=float)
    rows[..., :, 0] = xi_cov
    rows[..., :, 1] = u_cov * u_dot_xi[..., None]
    rows[..., :, 2] = r_cov * r_dot_xi[..., None]
    rows[..., :, 3] = 0.5 * (r_cov * u_dot_xi[..., None] + u_cov * r_dot_xi[..., None])

    u2 = scalar_contract(metric_inv, u_cov, u_cov)
    r2 = scalar_contract(metric_inv, r_cov, r_cov)
    ur = scalar_contract(metric_inv, u_cov, r_cov)
    trace_row = np.stack([np.full_like(u2, 3.0), u2, r2, ur], axis=-1)

    smin = np.full(u2.shape, np.nan, dtype=float)
    rank = np.zeros(u2.shape, dtype=float)
    for i_raw, j_raw in np.ndindex(u2.shape):
        i = int(i_raw)
        j = int(j_raw)
        ns = nullspace_of_row(trace_row[i, j])
        reduced = rows[i, j] @ ns
        if reduced.size == 0 or not np.all(np.isfinite(reduced)):
            continue
        s = np.linalg.svd(reduced, compute_uv=False)
        if s.size == 0:
            continue
        tol = max(reduced.shape) * np.finfo(float).eps * float(s[0])
        rank[i, j] = float(np.count_nonzero(s > tol))
        smin[i, j] = float(s[-1])
    return rank, smin


def analyze_case(args: argparse.Namespace, ref: dict[str, object], tau: float) -> dict[str, object]:
    case = build_case(args, ref, tau)
    metric_inv = case.metric_inv
    u_cov = case.u_cov
    r_cov = case.r_cov
    u_up = raise_covector(metric_inv, u_cov)
    r_up = raise_covector(metric_inv, r_cov)

    u2 = scalar_contract(metric_inv, u_cov, u_cov)
    r2 = scalar_contract(metric_inv, r_cov, r_cov)
    ur = scalar_contract(metric_inv, u_cov, r_cov)
    delta = u2 * r2 - ur * ur
    delta_scale = np.maximum(np.abs(u2 * r2) + np.abs(ur * ur), 1.0e-300)
    delta_rel = np.abs(delta) / delta_scale

    r_parallel = (ur / np.maximum(u2, 1.0e-300))[..., None] * u_cov
    r_perp = r_cov - r_parallel
    rperp2 = scalar_contract(metric_inv, r_perp, r_perp)
    rperp_rel = np.sqrt(np.abs(rperp2)) / np.maximum(np.sqrt(np.abs(r2)), 1.0e-300)

    xi_dot_u = u_up[..., 0]
    xi_dot_r = r_up[..., 0]
    w_cov = xi_dot_r[..., None] * u_cov - xi_dot_u[..., None] * r_cov
    w2 = scalar_contract(metric_inv, w_cov, w_cov)
    w2_scale = np.maximum(np.abs(xi_dot_r * xi_dot_r * u2) + np.abs(xi_dot_u * xi_dot_u * r2), 1.0e-300)
    w2_rel = np.abs(w2) / w2_scale

    hj_nonchar = np.abs(u_up[..., 0]) / np.maximum(np.sqrt(np.abs(u2)), 1.0e-300)
    lambda_rank, lambda_smin = lambda_time_symbol_after_trace(case)

    region = getattr(case, args.region)
    weight = np.sqrt(np.maximum(case.rho_a, 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    main = (
        region
        & (delta_rel > float(args.delta_rel_min))
        & (w2_rel > float(args.w2_rel_min))
        & (hj_nonchar > float(args.hj_nonchar_min))
        & (lambda_rank >= 3.0)
    )
    low_dim = region & ((delta_rel <= float(args.delta_rel_min)) | (rperp_rel <= float(args.rperp_rel_min)))
    w_characteristic = region & (delta_rel > float(args.delta_rel_min)) & (w2_rel <= float(args.w2_rel_min))
    hj_degenerate = region & (hj_nonchar <= float(args.hj_nonchar_min))
    lambda_bad = region & (lambda_rank < 3.0)
    unresolved = region & ~(main | low_dim | w_characteristic | hj_degenerate | lambda_bad)

    labels = np.full(case.rho_a.shape, -1, dtype=int)
    labels[region] = 0
    labels[main] = 1
    labels[low_dim] = 2
    labels[w_characteristic] = 3
    labels[hj_degenerate] = 4
    labels[lambda_bad] = 5
    labels[unresolved] = 6

    counts = {
        "region": int(np.count_nonzero(region)),
        "main_patch": int(np.count_nonzero(main)),
        "low_dimensional_basis_patch": int(np.count_nonzero(low_dim)),
        "w_characteristic_patch": int(np.count_nonzero(w_characteristic)),
        "hj_degenerate_patch": int(np.count_nonzero(hj_degenerate)),
        "lambda_time_rank_bad": int(np.count_nonzero(lambda_bad)),
        "unresolved": int(np.count_nonzero(unresolved)),
    }
    denom = max(counts["region"], 1)
    fractions = {name: float(value / denom) for name, value in counts.items() if name != "region"}

    plot_path = args.output / f"d_cauchy_patch_atlas_{tau_tag(tau)}.png"
    render_case_plot(
        plot_path,
        case,
        region,
        labels,
        delta_rel,
        w2_rel,
        hj_nonchar,
        lambda_smin,
    )
    return {
        "tau": float(tau),
        "plot_png": str(plot_path.resolve()),
        "shape": [int(case.rho_a.shape[0]), int(case.rho_a.shape[1])],
        "counts": counts,
        "fractions": fractions,
        "stats": {
            "delta_rel": weighted_stats(delta_rel[region], weight[region]),
            "rperp_rel": weighted_stats(rperp_rel[region], weight[region]),
            "w2_rel": weighted_stats(w2_rel[region], weight[region]),
            "hj_nonchar_abs_ut_over_sqrt_u2": weighted_stats(hj_nonchar[region], weight[region]),
            "lambda_trace_reduced_time_rank": weighted_stats(lambda_rank[region], weight[region]),
            "lambda_trace_reduced_time_smin": weighted_stats(lambda_smin[region], weight[region]),
        },
    }


def render_case_plot(
    out_path: Path,
    case,
    region: np.ndarray,
    labels: np.ndarray,
    delta_rel: np.ndarray,
    w2_rel: np.ndarray,
    hj_nonchar: np.ndarray,
    lambda_smin: np.ndarray,
) -> None:
    x_um = case.x * HBAR_C_EV_M * 1.0e6
    z_um = case.z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.0), constrained_layout=True)
    fields = [
        ("rho_A", case.rho_a, "viridis", False),
        ("patch labels", labels.astype(float), "tab10", False),
        ("log10 Delta_rel", np.log10(np.maximum(delta_rel, 1.0e-16)), "magma", True),
        ("log10 |w^2|_rel", np.log10(np.maximum(w2_rel, 1.0e-16)), "magma", True),
        ("log10 HJ nonchar", np.log10(np.maximum(hj_nonchar, 1.0e-16)), "cividis", True),
        ("log10 lambda smin", np.log10(np.maximum(lambda_smin, 1.0e-300)), "plasma", True),
    ]
    for ax, (title, data, cmap, masked_percentile) in zip(axes.ravel(), fields):
        arr = np.asarray(data, dtype=float)
        if title == "patch labels":
            im = ax.pcolormesh(xg, zg, arr, shading="auto", cmap=cmap, vmin=0, vmax=6)
        else:
            if masked_percentile and np.any(region):
                vals = arr[region & np.isfinite(arr)]
            else:
                vals = arr[np.isfinite(arr)]
            if vals.size:
                vmin = float(np.percentile(vals, 1.0))
                vmax = float(np.percentile(vals, 99.0))
                if vmax <= vmin:
                    vmax = vmin + 1.0
            else:
                vmin, vmax = 0.0, 1.0
            im = ax.pcolormesh(xg, zg, arr, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour(xg, zg, region.astype(float), levels=[0.5], colors="white", linewidths=0.5)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(
        "labels: 1 main, 2 low-dimensional E_u, 3 w^2=0 characteristic, "
        "4 HJ degenerate, 5 lambda-rank bad, 6 unresolved",
        fontsize=10,
    )
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    taus = [float(x.strip()) for x in args.taus.split(",") if x.strip()]
    cases = [analyze_case(args, ref, tau) for tau in taus]
    report = {
        "parameters": {
            "taus": taus,
            "region": str(args.region),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "delta_rel_min": float(args.delta_rel_min),
            "rperp_rel_min": float(args.rperp_rel_min),
            "w2_rel_min": float(args.w2_rel_min),
            "hj_nonchar_min": float(args.hj_nonchar_min),
        },
        "definitions": {
            "main_patch": "Region where Delta, w^2, HJ time characteristic, and trace-reduced lambda time symbol are all non-degenerate.",
            "Delta": "u^2 r^2-(u.r)^2. Delta=0 means the {u,r} tensor chart loses rank.",
            "w2": "w_mu=(xi.r)u_mu-(xi.u)r_mu with xi=dt. w^2=0 means ordinary trace0 does not close the metric principal gap on that Cauchy surface.",
            "HJ_nonchar": "|u^t|/sqrt(|u^2|). Small values mean lab-time slices are bad for solving the HJ branch.",
            "lambda_smin": "Smallest singular value of the time-principal div(C)=0 symbol after trace0 elimination.",
            "scope": "This is a Cauchy-patch admissibility diagnostic, not a full nonlinear D evolution.",
        },
        "cases": cases,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plots": [case["plot_png"] for case in cases],
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument("--delta-rel-min", type=float, default=1.0e-8)
    parser.add_argument("--rperp-rel-min", type=float, default=1.0e-6)
    parser.add_argument("--w2-rel-min", type=float, default=1.0e-8)
    parser.add_argument("--hj-nonchar-min", type=float, default=1.0e-8)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
