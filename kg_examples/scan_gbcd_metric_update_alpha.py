from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from solve_gbcd_full_linear_metric_update_sparse import load_lambda_cov
from test_gbcd_v0_one_step_lambda import weighted_stats


def parse_alphas(text: str) -> list[float]:
    out: list[float] = []
    for part in text.split(","):
        item = part.strip()
        if item:
            out.append(float(item))
    if not out:
        raise ValueError("empty alpha list")
    return out


def optional_array(data: np.lib.npyio.NpzFile, key: str, shape_like: np.ndarray) -> np.ndarray:
    if key in data.files:
        return np.asarray(data[key], dtype=float)
    return np.zeros_like(shape_like)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    case = build_case(args, ref, float(args.tau))
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(case.x[1] - case.x[0])
    dz = float(case.z[1] - case.z[0])

    coeff_path = args.coefficients
    cov0, coeff_mask = load_lambda_cov({0: case}, coeff_path, str(args.atom_family))
    data = np.load(args.data)
    mask = np.asarray(data["mask"], dtype=bool)
    if not np.array_equal(mask, coeff_mask):
        raise ValueError("data mask and coefficient mask differ")

    delta_m = optional_array(data, "delta_metric_minus", full_case.geom_m.metric_cov)
    delta_0 = optional_array(data, "delta_metric_center", full_case.geom_0.metric_cov)
    delta_p = optional_array(data, "delta_metric_plus", full_case.geom_p.metric_cov)
    weights = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))

    rows = []
    for alpha in parse_alphas(args.alphas):
        metric_m = full_case.geom_m.metric_cov + alpha * delta_m
        metric_0 = full_case.geom_0.metric_cov + alpha * delta_0
        metric_p = full_case.geom_p.metric_cov + alpha * delta_p
        geom, _ = make_time_geometry(metric_0, metric_m, metric_p, dt, dx, dz, with_dgamma=False)
        ein = geom.ricci - 0.5 * metric_0 * geom.r_scalar[..., None, None]
        exact_r_need = ein - case.source_tilde
        rel = tensor_norm(cov0 - exact_r_need) / np.maximum(tensor_norm(exact_r_need), 1.0e-300)
        stats = weighted_stats(rel[mask], weights)
        rows.append(
            {
                "alpha": float(alpha),
                "weighted_mean": float(stats["weighted_mean"]),
                "p50": float(stats["p50"]),
                "p95": float(stats["p95"]),
                "max": float(stats["max"]),
            }
        )

    best = min(rows, key=lambda row: row["weighted_mean"])
    alpha_values = np.asarray([row["alpha"] for row in rows], dtype=float)
    weighted_mean = np.asarray([row["weighted_mean"] for row in rows], dtype=float)
    p95 = np.asarray([row["p95"] for row in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    positive = alpha_values > 0.0
    if np.any(positive):
        ax.set_xscale("symlog", linthresh=max(float(np.min(alpha_values[positive])) * 0.5, 1.0e-12))
    ax.plot(alpha_values, weighted_mean, "o-", label="weighted mean")
    ax.plot(alpha_values, p95, "s--", label="p95")
    ax.axvline(best["alpha"], color="black", linewidth=0.8, linestyle=":", label=f"best alpha={best['alpha']:.3g}")
    ax.set_xlabel("alpha multiplying the solved metric update")
    ax.set_ylabel("exact nonlinear relative residual")
    ax.set_title("Damped metric-update line search")
    ax.grid(True, alpha=0.25)
    ax.legend()
    plot_path = args.output / "gbcd_metric_update_alpha_scan.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "atom_family": str(args.atom_family),
            "data": str(args.data.resolve()),
            "coefficients": str(coeff_path.resolve()),
            "alphas": [float(x) for x in alpha_values],
        },
        "definition": {
            "alpha": "multiplies the already solved metric update before exact nonlinear field-equation back-substitution",
            "weighted_mean": "rho-weighted mean of ||C - (G - T)|| / ||G - T|| over the fitting mask",
        },
        "best": best,
        "rows": rows,
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
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
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
    parser.add_argument(
        "--atom-family",
        choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"],
        default="auto",
    )
    parser.add_argument(
        "--alphas",
        default="0,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.2,0.3,0.5,0.7,1.0",
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
