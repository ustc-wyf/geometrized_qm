from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from analyze_bcd_residuals_from_a_reference import stress_tensor_tilde
from coordinate_matter_evolution import conservative_density_from_rho_u, coordinate_matter_rhs_covector
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_reduced_dynamic_same_initial import rk4_matter_step
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det
from solve_gbcd_full_linear_metric_update_sparse import load_lambda_cov
from test_gbcd_v0_one_step_lambda import weighted_stats


def relative_l1(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(np.abs(candidate[valid] - reference[valid])) / denom)


def apply_plus_guard(
    *,
    args: argparse.Namespace,
    ref: dict[str, object],
    case0,
    full_case,
    metric_center: np.ndarray,
    metric_plus: np.ndarray,
    delta_plus: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, dict[str, object], np.ndarray]:
    if args.plus_guard == "none":
        return metric_plus, {"mode": "none"}, np.ones(case0.rho_a.shape, dtype=float)

    rho_tilde0 = np.maximum(case0.measure_tilde / np.maximum(safe_sqrt_abs_det(metric_center), 1.0e-300), float(args.rho_floor))
    cons0 = conservative_density_from_rho_u(
        rho=np.where(case0.support, rho_tilde0, 0.0),
        u_x=case0.u_cov[..., 1],
        u_z=case0.u_cov[..., 2],
        metric_cov_txz=metric_center,
        mass=float(ref["params"].m),
        branch="negative_frequency",
        u_t_reference=case0.u_cov[..., 0],
        rho_floor=float(args.rho_floor),
    )
    n_next, ux_next, uz_next, current_center = rk4_matter_step(
        n_cons=cons0["n_cons"],
        u_x=case0.u_cov[..., 1],
        u_z=case0.u_cov[..., 2],
        metric_cov=metric_center,
        mass=float(ref["params"].m),
        dx=dx,
        dz=dz,
        dt=dt,
        u_t_reference=case0.u_cov[..., 0],
        active_mask=case0.support,
    )
    current_plus = coordinate_matter_rhs_covector(
        n_cons=n_next,
        u_x=ux_next,
        u_z=uz_next,
        metric_cov_txz=metric_plus,
        mass=float(ref["params"].m),
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=current_center["u_t"],
    )
    det_plus = np.linalg.det(metric_plus.reshape(-1, 3, 3)).reshape(metric_plus.shape[:2])
    bad = (current_plus["discriminant"] < 0.0) | (det_plus <= float(args.plus_guard_det_floor))
    bad &= case0.support
    if int(args.plus_guard_dilation) > 0:
        bad = base.dilate_mask_8(bad, int(args.plus_guard_dilation)) & case0.support
    scale = np.where(bad, 0.0, 1.0)
    guarded_plus = full_case.geom_p.metric_cov + scale[..., None, None] * delta_plus
    return (
        guarded_plus,
        {
            "mode": str(args.plus_guard),
            "dilation": int(args.plus_guard_dilation),
            "det_floor": float(args.plus_guard_det_floor),
            "frozen_support_points": int(np.count_nonzero(bad & case0.support)),
            "frozen_core10_points": int(np.count_nonzero(bad & case0.core10)),
            "bad_before_negative_discriminant_support_points": int(
                np.count_nonzero((current_plus["discriminant"] < 0.0) & case0.support)
            ),
            "bad_before_nonpositive_det_support_points": int(np.count_nonzero((det_plus <= float(args.plus_guard_det_floor)) & case0.support)),
        },
        scale,
    )


def render_summary(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho_a: np.ndarray,
    measure_a: np.ndarray,
    exact_rel: np.ndarray,
    support: np.ndarray,
    fit_mask: np.ndarray,
    raw_y: np.ndarray,
) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 10.0), constrained_layout=True)
    panels = [
        (rho_a, "A rho on initial slice", "viridis", False),
        (measure_a, "initial transformed measure", "viridis", False),
        (exact_rel, "D field-equation residual", "magma", True),
        (np.sign(raw_y) * np.log10(1.0 + np.abs(raw_y)), "sign(raw_y) log10(1+|raw_y|)", "coolwarm", False),
    ]
    for ax, (field, title, cmap, log_scale) in zip(axes.flat, panels):
        data = np.asarray(field, dtype=float)
        if log_scale:
            data = np.log10(1.0 + np.maximum(data, 0.0))
            title = "log10(1+" + title + ")"
        im = ax.pcolormesh(xg, zg, data, shading="auto", cmap=cmap)
        ax.contour(xg, zg, support.astype(float), levels=[0.5], colors=["white"], linewidths=0.6)
        ax.contour(xg, zg, fit_mask.astype(float), levels=[0.5], colors=["cyan"], linewidths=0.6)
        ax.set_title(title + "; white=support, cyan=fit mask")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax)
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

    projection = np.load(args.plus_projection_data)
    metric_minus = np.asarray(
        projection["corrected_metric_minus"] if "corrected_metric_minus" in projection.files else full_case.geom_m.metric_cov,
        dtype=float,
    )
    metric_center = np.asarray(
        projection["corrected_metric_center"] if "corrected_metric_center" in projection.files else full_case.geom_0.metric_cov,
        dtype=float,
    )
    metric_plus = np.asarray(projection["corrected_metric_plus"], dtype=float)
    delta_plus = np.asarray(
        projection["delta_metric_plus"] if "delta_metric_plus" in projection.files else metric_plus - full_case.geom_p.metric_cov,
        dtype=float,
    )
    cov0, fit_mask = load_lambda_cov({0: case0}, args.coefficients, str(args.atom_family))

    metric_plus, guard_info, guard_scale = apply_plus_guard(
        args=args,
        ref=ref,
        case0=case0,
        full_case=full_case,
        metric_center=metric_center,
        metric_plus=metric_plus,
        delta_plus=delta_plus,
        dt=dt,
        dx=dx,
        dz=dz,
    )

    geom, _ = make_time_geometry(metric_center, metric_minus, metric_plus, dt, dx, dz, with_dgamma=False)
    ein = geom.ricci - 0.5 * metric_center * geom.r_scalar[..., None, None]
    sqrt_abs_g = safe_sqrt_abs_det(metric_center)
    rho_tilde = np.maximum(case0.measure_tilde / np.maximum(sqrt_abs_g, 1.0e-300), float(args.rho_floor))
    stress, tilde_x = stress_tensor_tilde(
        metric_center,
        geom.metric_inv,
        rho_tilde,
        case0.u_cov[..., 0],
        case0.u_cov[..., 1],
        case0.u_cov[..., 2],
        m=float(ref["params"].m),
    )
    source = stress / (float(args.mp) * float(args.mp))
    r_need = ein - source
    exact_rel = tensor_norm(cov0 - r_need) / np.maximum(tensor_norm(r_need), 1.0e-300)
    weights = np.sqrt(np.maximum(case0.rho_a[fit_mask], 0.0) / max(float(np.max(case0.rho_a)), 1.0e-300))

    cons = conservative_density_from_rho_u(
        rho=np.where(case0.support, rho_tilde, 0.0),
        u_x=case0.u_cov[..., 1],
        u_z=case0.u_cov[..., 2],
        metric_cov_txz=metric_center,
        mass=float(ref["params"].m),
        branch="negative_frequency",
        u_t_reference=case0.u_cov[..., 0],
        rho_floor=float(args.rho_floor),
    )
    ell = float(args.ell if args.ell is not None else args.ell_over_planck / float(args.mp))
    raw_y = ell * ell * geom.r_scalar
    plot_path = args.output / "gbcd_plus_initial_package_summary.png"
    render_summary(plot_path, case0.x, case0.z, case0.rho_a, case0.measure_tilde, exact_rel, case0.support, fit_mask, raw_y)

    package_path = args.output / "gbcd_plus_initial_package.npz"
    np.savez_compressed(
        package_path,
        x=case0.x,
        z=case0.z,
        support=case0.support,
        trusted=case0.trusted,
        core10=case0.core10,
        fit_mask=fit_mask,
        metric_minus=metric_minus,
        metric_center=metric_center,
        metric_plus=metric_plus,
        metric_inv_center=geom.metric_inv,
        ricci_center=geom.ricci,
        r_scalar_center=geom.r_scalar,
        raw_y=raw_y,
        rho_A=case0.rho_a,
        measure_tilde=case0.measure_tilde,
        rho_tilde=rho_tilde,
        u_t=case0.u_cov[..., 0],
        u_x=case0.u_cov[..., 1],
        u_z=case0.u_cov[..., 2],
        n_cons=cons["n_cons"],
        flux_x=cons["flux_x"],
        flux_z=cons["flux_z"],
        discriminant=cons["discriminant"],
        auxiliary_cov=cov0,
        source_tilde=source,
        r_need=r_need,
        exact_relative_residual=exact_rel,
        plus_guard_scale=guard_scale,
    )

    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "dt_ev_inv": float(dt),
            "dt_old": float(args.force_probe_dt_old),
            "dx_ev_inv": float(dx),
            "dz_ev_inv": float(dz),
            "mass_ev": float(ref["params"].m),
            "mp_ev": float(args.mp),
            "ell_ev_inv": float(ell),
            "plus_projection_data": str(args.plus_projection_data.resolve()),
            "coefficients": str(args.coefficients.resolve()),
            "plus_guard": str(args.plus_guard),
        },
        "plus_guard": guard_info,
        "meaning": {
            "metric_minus_center_plus": "Three tilde-metric time slices (g_-, g_0, g_+) centered on the initial slice.  In plus-only mode g_0 and matter are unchanged; g_+ is chosen by the D/gBCD field equation.",
            "measure_tilde": "sqrt(|gtilde|) * rho_tilde, initialized from the A/KG state by |X| rho_A / m^2.",
            "auxiliary_cov": "The fitted auxiliary tensor C_mn used in the equation-first gBCD field equation.",
            "exact_relative_residual": "Pointwise relative error after substituting this package into Gtilde - Ttilde/Mp^2 = C. Smaller is better.",
        },
        "diagnostics": {
            "fit_points": int(np.count_nonzero(fit_mask)),
            "support_points": int(np.count_nonzero(case0.support)),
            "exact_relative_residual": weighted_stats(exact_rel[fit_mask], weights),
            "rho_pullback_relative_l1_support": relative_l1(case0.measure_tilde, sqrt_abs_g * rho_tilde, case0.support),
            "mass_shell_defect_abs_p95_fit": float(np.percentile(np.abs((tilde_x - float(ref["params"].m) ** 2)[fit_mask]), 95.0)),
            "negative_discriminant_fraction_fit": float(np.count_nonzero(cons["discriminant"][fit_mask] < 0.0) / max(np.count_nonzero(fit_mask), 1)),
            "det_center_min_fit": float(np.min(np.linalg.det(metric_center.reshape(-1, 3, 3)).reshape(metric_center.shape[:2])[fit_mask])),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "package_npz": str(package_path.resolve()),
            "plot_png": str(plot_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plus-projection-data", type=Path, required=True)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=-3.5)
    parser.add_argument("--atom-family", choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"], default="auto")
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
    parser.add_argument("--plus-guard", choices=["none", "auto_bad_zero"], default="none")
    parser.add_argument("--plus-guard-dilation", type=int, default=0)
    parser.add_argument("--plus-guard-det-floor", type=float, default=0.0)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
