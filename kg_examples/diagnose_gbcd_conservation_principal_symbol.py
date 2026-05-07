from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import build_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


def stats(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "min": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "min": float(np.min(vals)),
        "max": float(np.max(vals)),
    }


def principal_matrix(case, k_cov: np.ndarray) -> np.ndarray:
    """Principal symbol for derivatives of lambda=(A,B,C,D).

    Conservation of C_mn=lambda_I E^I_mn has principal part

        k_n A + u_n (u.k) B + r_n (r.k) C
        + 1/2 (r_n (u.k) + u_n (r.k)) D.

    In the current 2+1 reduction this is a 3 x 4 matrix at each grid point.
    """

    ginv = case.metric_inv
    u_cov = case.u_cov
    r_cov = case.r_cov
    u_up = np.einsum("...ab,...b->...a", ginv, u_cov, optimize=True)
    r_up = np.einsum("...ab,...b->...a", ginv, r_cov, optimize=True)
    uk = np.einsum("...a,a->...", u_up, k_cov, optimize=True)
    rk = np.einsum("...a,a->...", r_up, k_cov, optimize=True)
    p = np.zeros(u_cov.shape[:-1] + (3, 4), dtype=float)
    p[..., :, 0] = k_cov
    p[..., :, 1] = u_cov * uk[..., None]
    p[..., :, 2] = r_cov * rk[..., None]
    p[..., :, 3] = 0.5 * (r_cov * uk[..., None] + u_cov * rk[..., None])
    return p


def analyze_symbol(case, region: str, k_name: str, k_cov: np.ndarray) -> dict[str, object]:
    mask = getattr(case, region)
    p = principal_matrix(case, np.asarray(k_cov, dtype=float))
    mats = p[mask]
    ranks = []
    smins = []
    smaxs = []
    conds = []
    nullity = []
    for mat in mats:
        s = np.linalg.svd(mat, compute_uv=False)
        tol = max(mat.shape) * np.finfo(float).eps * (float(s[0]) if s.size else 0.0)
        rank = int(np.count_nonzero(s > tol))
        ranks.append(rank)
        smins.append(float(s[-1]) if s.size else 0.0)
        smaxs.append(float(s[0]) if s.size else 0.0)
        conds.append(float(s[0] / max(s[-1], 1.0e-300)) if s.size else 0.0)
        nullity.append(4 - rank)
    return {
        "k_name": k_name,
        "k_cov": [float(x) for x in k_cov],
        "rank": stats(np.asarray(ranks, dtype=float)),
        "nullity": stats(np.asarray(nullity, dtype=float)),
        "singular_min": stats(np.asarray(smins, dtype=float)),
        "singular_max": stats(np.asarray(smaxs, dtype=float)),
        "condition": stats(np.asarray(conds, dtype=float)),
    }


def render_plot(path: Path, records: list[dict[str, object]]) -> None:
    labels = [f"{r['tau_label']} {r['k_name']}" for r in records]
    rank_p50 = [r["rank"]["p50"] for r in records]
    null_p50 = [r["nullity"]["p50"] for r in records]
    cond_p95 = [r["condition"]["p95"] for r in records]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    axes[0].bar(x, rank_p50)
    axes[0].set_title("principal rank p50")
    axes[1].bar(x, null_p50)
    axes[1].set_title("lambda nullity p50")
    axes[2].bar(x, cond_p95)
    axes[2].set_yscale("log")
    axes[2].set_title("condition p95")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    taus = [float(x.strip()) for x in args.taus.split(",") if x.strip()]
    k_list = {
        "lab_t": np.array([1.0, 0.0, 0.0]),
        "lab_x": np.array([0.0, 1.0, 0.0]),
        "lab_z": np.array([0.0, 0.0, 1.0]),
    }
    records = []
    for tau in taus:
        case = build_case(args, ref, tau)
        for k_name, k_cov in k_list.items():
            rec = analyze_symbol(case, args.region, k_name, k_cov)
            rec["tau"] = tau
            rec["tau_label"] = f"tau={tau:g}"
            rec["region"] = args.region
            rec["count"] = int(np.count_nonzero(getattr(case, args.region)))
            records.append(rec)
    plot_path = args.output / "gbcd_conservation_principal_symbol.png"
    render_plot(plot_path, records)
    report = {
        "parameters": {
            "taus": taus,
            "region": args.region,
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
        },
        "definition": {
            "equation": "Principal part of nabla^mu(lambda_I E^I_mn)=0 for lambda=(A,B,C,D).",
            "columns": "A, B, C, D; rows are covector index nu in 2+1 txz.",
            "principal_symbol": "P_nI(k)=[k_n, u_n(u.k), r_n(r.k), 1/2(r_n(u.k)+u_n(r.k))].",
            "interpretation": "In 2+1, P is 3x4, so conservation alone cannot determine all four lambda derivatives; nullity >= 1 means an additional state equation/gauge is structurally required.",
        },
        "records": records,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
