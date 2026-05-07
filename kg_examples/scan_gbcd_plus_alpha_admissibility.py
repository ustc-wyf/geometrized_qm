from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from coordinate_matter_evolution import conservative_density_from_rho_u, coordinate_matter_rhs_covector
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_reduced_dynamic_same_initial import rk4_matter_step
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det
from solve_gbcd_full_linear_metric_update_sparse import load_lambda_cov
from test_gbcd_v0_one_step_lambda import weighted_stats


def parse_alphas(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise ValueError("empty alpha list")
    return values


def relative_l1(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(np.abs(candidate[valid] - reference[valid])) / denom)


def optional_array(data: np.lib.npyio.NpzFile, key: str, shape_like: np.ndarray) -> np.ndarray:
    if key in data.files:
        return np.asarray(data[key], dtype=float)
    return np.zeros_like(shape_like)


def render(out_path: Path, rows: list[dict[str, float]]) -> None:
    alpha = np.asarray([row["alpha"] for row in rows], dtype=float)
    residual = np.asarray([row["field_residual_weighted_mean"] for row in rows], dtype=float)
    rho_l1 = np.asarray([row["one_step_plus_vs_A_relative_l1_support"] for row in rows], dtype=float)
    neg_disc = np.asarray([row["negative_discriminant_fraction_plus"] for row in rows], dtype=float)
    det_min = np.asarray([row["det_plus_min_support"] for row in rows], dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.8), constrained_layout=True)
    axes[0, 0].plot(alpha, residual, marker="o")
    axes[0, 0].set_title("D field-equation residual")
    axes[0, 0].set_xlabel("alpha")
    axes[0, 0].set_ylabel("weighted mean residual")
    axes[0, 0].grid(alpha=0.25)

    axes[0, 1].plot(alpha, rho_l1, marker="o", color="tab:orange")
    axes[0, 1].set_title("one-step rho/measure deviation from A")
    axes[0, 1].set_xlabel("alpha")
    axes[0, 1].set_ylabel("relative L1 on support")
    axes[0, 1].grid(alpha=0.25)

    axes[1, 0].plot(alpha, neg_disc, marker="o", color="tab:red")
    axes[1, 0].set_title("negative mass-shell discriminant fraction")
    axes[1, 0].set_xlabel("alpha")
    axes[1, 0].set_ylabel("fraction on support")
    axes[1, 0].grid(alpha=0.25)

    axes[1, 1].plot(alpha, det_min, marker="o", color="tab:green")
    axes[1, 1].axhline(0.0, color="black", linewidth=0.8)
    axes[1, 1].set_title("min det(g_plus) on support")
    axes[1, 1].set_xlabel("alpha")
    axes[1, 1].set_ylabel("min det")
    axes[1, 1].grid(alpha=0.25)

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    case = build_case(args, ref, float(args.tau))
    case_next = build_case(args, ref, float(args.tau) + float(args.force_probe_dt_old))
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(case.x[1] - case.x[0])
    dz = float(case.z[1] - case.z[0])
    cov0, coeff_mask = load_lambda_cov({0: case}, args.coefficients, str(args.atom_family))
    data = np.load(args.data)
    mask = np.asarray(data["mask"], dtype=bool)
    if not np.array_equal(mask, coeff_mask):
        raise ValueError("data mask and coefficient mask differ")

    delta_m = optional_array(data, "delta_metric_minus", full_case.geom_m.metric_cov)
    delta_0 = optional_array(data, "delta_metric_center", full_case.geom_0.metric_cov)
    delta_p = optional_array(data, "delta_metric_plus", full_case.geom_p.metric_cov)
    rho_tilde0 = np.maximum(case.measure_tilde / np.maximum(safe_sqrt_abs_det(full_case.geom_0.metric_cov), 1.0e-300), float(args.rho_floor))
    cons0 = conservative_density_from_rho_u(
        rho=np.where(case.support, rho_tilde0, 0.0),
        u_x=case.u_cov[..., 1],
        u_z=case.u_cov[..., 2],
        metric_cov_txz=full_case.geom_0.metric_cov,
        mass=float(ref["params"].m),
        branch="negative_frequency",
        u_t_reference=case.u_cov[..., 0],
        rho_floor=float(args.rho_floor),
    )
    n_next, ux_next, uz_next, current_center = rk4_matter_step(
        n_cons=cons0["n_cons"],
        u_x=case.u_cov[..., 1],
        u_z=case.u_cov[..., 2],
        metric_cov=full_case.geom_0.metric_cov,
        mass=float(ref["params"].m),
        dx=dx,
        dz=dz,
        dt=dt,
        u_t_reference=case.u_cov[..., 0],
        active_mask=case.support,
    )
    weights = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    rows = []
    for alpha in parse_alphas(args.alphas):
        metric_m = full_case.geom_m.metric_cov + float(alpha) * delta_m
        metric_0 = full_case.geom_0.metric_cov + float(alpha) * delta_0
        metric_p = full_case.geom_p.metric_cov + float(alpha) * delta_p
        geom, _ = make_time_geometry(metric_0, metric_m, metric_p, dt, dx, dz, with_dgamma=False)
        ein = geom.ricci - 0.5 * metric_0 * geom.r_scalar[..., None, None]
        rel = tensor_norm(cov0 - (ein - case.source_tilde)) / np.maximum(tensor_norm(ein - case.source_tilde), 1.0e-300)
        stats = weighted_stats(rel[mask], weights)
        current_plus = coordinate_matter_rhs_covector(
            n_cons=n_next,
            u_x=ux_next,
            u_z=uz_next,
            metric_cov_txz=metric_p,
            mass=float(ref["params"].m),
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=current_center["u_t"],
        )
        measure_plus = safe_sqrt_abs_det(metric_p) * current_plus["rho"]
        det_plus = np.linalg.det(metric_p.reshape(-1, 3, 3)).reshape(metric_p.shape[:2])
        rows.append(
            {
                "alpha": float(alpha),
                "field_residual_weighted_mean": float(stats["weighted_mean"]),
                "field_residual_p95": float(stats["p95"]),
                "one_step_plus_vs_A_relative_l1_support": relative_l1(case_next.measure_tilde, measure_plus, case.support),
                "negative_discriminant_fraction_plus": float(np.count_nonzero(current_plus["discriminant"][case.support] < 0.0) / max(np.count_nonzero(case.support), 1)),
                "disc_plus_min_support": float(np.min(current_plus["discriminant"][case.support])),
                "disc_plus_p01_support": float(np.percentile(current_plus["discriminant"][case.support], 1.0)),
                "det_plus_min_support": float(np.min(det_plus[case.support])),
            }
        )

    feasible = [
        row
        for row in rows
        if row["negative_discriminant_fraction_plus"] <= float(args.max_negative_discriminant_fraction)
        and row["one_step_plus_vs_A_relative_l1_support"] <= float(args.max_measure_deviation)
    ]
    best_residual = min(rows, key=lambda row: row["field_residual_weighted_mean"])
    best_feasible = min(feasible, key=lambda row: row["field_residual_weighted_mean"]) if feasible else None
    plot_path = args.output / "gbcd_plus_alpha_admissibility_scan.png"
    render(plot_path, rows)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "dt_old": float(args.force_probe_dt_old),
            "data": str(args.data.resolve()),
            "coefficients": str(args.coefficients.resolve()),
            "max_measure_deviation": float(args.max_measure_deviation),
            "max_negative_discriminant_fraction": float(args.max_negative_discriminant_fraction),
        },
        "meaning": {
            "alpha": "Scale applied to the solved future-slice metric correction delta g_+. alpha=1 is the full field-equation update; alpha=0 is the A-reference future slice.",
            "field_residual": "How well the D/gBCD field equation is satisfied at the initial slice.",
            "admissibility": "Whether the one-step matter state can be reconstructed on the candidate future metric without negative mass-shell discriminant or huge measure drift.",
        },
        "best_residual": best_residual,
        "best_feasible": best_feasible,
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
    parser.add_argument("--tau", type=float, default=-3.5)
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
    parser.add_argument("--atom-family", choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"], default="auto")
    parser.add_argument("--alphas", default="0,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.05,0.075,0.1,0.15,0.2,0.3,0.5,0.7,1.0")
    parser.add_argument("--max-measure-deviation", type=float, default=0.05)
    parser.add_argument("--max-negative-discriminant-fraction", type=float, default=0.0)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
