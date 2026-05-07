from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_bcd_residuals_from_a_reference import stress_tensor_tilde
from coordinate_matter_evolution import solve_covector_time_component
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det
from solve_gbcd_full_linear_metric_update_sparse import load_lambda_cov
from solve_gbcd_joint_initial_projection_sparse import weighted_pullback_l1_from_eta
from test_gbcd_v0_one_step_lambda import weighted_stats


def evaluate_alpha(
    *,
    alpha: float,
    data: np.lib.npyio.NpzFile,
    cov0: np.ndarray,
    case0,
    full_case,
    dt: float,
    dx: float,
    dz: float,
    mass: float,
    mp: float,
    fit_mask: np.ndarray,
    rho_floor: float,
) -> dict[str, float]:
    metric_m = full_case.geom_m.metric_cov + float(alpha) * np.asarray(data["delta_metric_minus"], dtype=float)
    metric_0 = full_case.geom_0.metric_cov + float(alpha) * np.asarray(data["delta_metric_center"], dtype=float)
    metric_p = full_case.geom_p.metric_cov + float(alpha) * np.asarray(data["delta_metric_plus"], dtype=float)
    eta = float(alpha) * np.asarray(data["eta_log_measure"], dtype=float)

    geom, _ = make_time_geometry(metric_0, metric_m, metric_p, dt, dx, dz, with_dgamma=False)
    ein = geom.ricci - 0.5 * metric_0 * geom.r_scalar[..., None, None]
    measure = case0.measure_tilde * np.exp(np.clip(eta, -50.0, 50.0))
    sqrt_abs_g = safe_sqrt_abs_det(metric_0)
    rho_tilde = np.maximum(measure / np.maximum(sqrt_abs_g, 1.0e-300), rho_floor)
    phase = solve_covector_time_component(
        geom.metric_inv,
        case0.u_cov[..., 1],
        case0.u_cov[..., 2],
        mass,
        branch="negative_frequency",
        u_t_reference=case0.u_cov[..., 0],
    )
    stress, tilde_x = stress_tensor_tilde(
        metric_0,
        geom.metric_inv,
        rho_tilde,
        phase["u_t"],
        case0.u_cov[..., 1],
        case0.u_cov[..., 2],
        m=mass,
    )
    source = stress / (mp * mp)
    r_need = ein - source
    rel = tensor_norm(cov0 - r_need) / np.maximum(tensor_norm(r_need), 1.0e-300)
    weights = np.sqrt(np.maximum(case0.rho_a[fit_mask], 0.0) / max(float(np.max(case0.rho_a)), 1.0e-300))
    det = np.linalg.det(metric_0.reshape(-1, 3, 3)).reshape(metric_0.shape[:2])
    return {
        "alpha": float(alpha),
        "exact_weighted_mean": float(weighted_stats(rel[fit_mask], weights)["weighted_mean"]),
        "exact_p50": float(weighted_stats(rel[fit_mask], weights)["p50"]),
        "exact_p95": float(weighted_stats(rel[fit_mask], weights)["p95"]),
        "rho_pullback_weighted_l1": weighted_pullback_l1_from_eta(case0.rho_a, eta, fit_mask),
        "eta_abs_max": float(np.max(np.abs(eta[fit_mask]))) if np.any(fit_mask) else 0.0,
        "negative_discriminant_fraction": float(np.count_nonzero(phase["discriminant"][fit_mask] < 0.0) / max(np.count_nonzero(fit_mask), 1)),
        "mass_shell_abs_p95": float(np.percentile(np.abs((tilde_x - mass * mass)[fit_mask]), 95.0)) if np.any(fit_mask) else 0.0,
        "det_min": float(np.min(det[fit_mask])) if np.any(fit_mask) else 0.0,
    }


def render(out_path: Path, rows: list[dict[str, float]]) -> None:
    alpha = np.asarray([r["alpha"] for r in rows], dtype=float)
    exact = np.asarray([r["exact_weighted_mean"] for r in rows], dtype=float)
    exact_p95 = np.asarray([r["exact_p95"] for r in rows], dtype=float)
    rho = np.asarray([r["rho_pullback_weighted_l1"] for r in rows], dtype=float)
    neg_disc = np.asarray([r["negative_discriminant_fraction"] for r in rows], dtype=float)
    det_min = np.asarray([r["det_min"] for r in rows], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.8), constrained_layout=True)
    axes[0, 0].plot(alpha, exact, marker="o", label="weighted mean")
    axes[0, 0].plot(alpha, exact_p95, marker=".", label="p95")
    axes[0, 0].set_xlabel("alpha")
    axes[0, 0].set_ylabel("exact residual")
    axes[0, 0].set_title("D tensor residual after scaled joint update")
    axes[0, 0].legend()

    axes[0, 1].plot(alpha, rho, marker="o", color="tab:orange")
    axes[0, 1].set_xlabel("alpha")
    axes[0, 1].set_ylabel("weighted relative L1")
    axes[0, 1].set_title("rho pullback deviation")

    axes[1, 0].plot(alpha, neg_disc, marker="o", color="tab:red")
    axes[1, 0].set_xlabel("alpha")
    axes[1, 0].set_ylabel("fraction")
    axes[1, 0].set_title("negative mass-shell discriminant fraction")

    axes[1, 1].plot(alpha, det_min, marker="o", color="tab:green")
    axes[1, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[1, 1].set_xlabel("alpha")
    axes[1, 1].set_ylabel("min det(gtilde)")
    axes[1, 1].set_title("metric nondegeneracy check")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    case0 = build_case(args, ref, float(args.tau))
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(case0.x[1] - case0.x[0])
    dz = float(case0.z[1] - case0.z[0])
    cov0, fit_mask = load_lambda_cov({0: case0}, args.coefficients, str(args.atom_family))
    data = np.load(args.projection_data)
    alphas = np.unique(
        np.concatenate(
            [
                np.linspace(0.0, 1.0, int(args.linear_count)),
                np.asarray([float(x) for x in args.extra_alphas.split(",") if x.strip()], dtype=float),
            ]
        )
    )
    rows = [
        evaluate_alpha(
            alpha=float(alpha),
            data=data,
            cov0=cov0,
            case0=case0,
            full_case=full_case,
            dt=dt,
            dx=dx,
            dz=dz,
            mass=float(ref["params"].m),
            mp=float(args.mp),
            fit_mask=fit_mask,
            rho_floor=float(args.rho_floor),
        )
        for alpha in alphas
    ]
    best_exact = min(rows, key=lambda row: row["exact_weighted_mean"])
    feasible = [
        row
        for row in rows
        if row["rho_pullback_weighted_l1"] <= float(args.max_rho_deviation)
        and row["negative_discriminant_fraction"] <= float(args.max_negative_discriminant_fraction)
        and row["det_min"] > 0.0
    ]
    best_feasible = min(feasible, key=lambda row: row["exact_weighted_mean"]) if feasible else None
    plot_path = args.output / "gbcd_joint_projection_alpha_scan.png"
    render(plot_path, rows)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "projection_data": str(args.projection_data.resolve()),
            "coefficients": str(args.coefficients.resolve()),
            "max_rho_deviation": float(args.max_rho_deviation),
            "max_negative_discriminant_fraction": float(args.max_negative_discriminant_fraction),
        },
        "definition": {
            "alpha": "Trust-region scale applied to the joint correction from the linearized solve: g -> g + alpha*delta_g, eta -> alpha*eta.",
            "purpose": "Check whether a smaller physically admissible step improves the exact nonlinear D residual without large rho pullback drift.",
        },
        "best_exact": best_exact,
        "best_feasible_under_user_like_bounds": best_feasible,
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
    parser.add_argument("--projection-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--atom-family", choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"], default="auto")
    parser.add_argument("--linear-count", type=int, default=41)
    parser.add_argument("--extra-alphas", type=str, default="0.001,0.002,0.005,0.01,0.02,0.05,0.1")
    parser.add_argument("--max-rho-deviation", type=float, default=0.05)
    parser.add_argument("--max-negative-discriminant-fraction", type=float, default=0.0)
    parser.add_argument("--full-resolution", type=int, default=96)
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
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
