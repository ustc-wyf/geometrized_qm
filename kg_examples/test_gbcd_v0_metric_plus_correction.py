from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from test_gbcd_v0_one_step_lambda import ATOM_NAMES, load_full_coeff, weighted_stats


def metric_plus_variable_index(mask: np.ndarray) -> tuple[list[tuple[int, int, int, int]], np.ndarray]:
    variables: list[tuple[int, int, int, int]] = []
    col = -np.ones(mask.shape + (3, 3), dtype=int)
    for i, j in np.argwhere(mask):
        for a, b in base.SYMMETRIC_COMPONENTS:
            n = len(variables)
            variables.append((int(i), int(j), int(a), int(b)))
            col[int(i), int(j), int(a), int(b)] = n
            col[int(i), int(j), int(b), int(a)] = n
    return variables, col


def center_einstein_from_metric_plus(full_case, metric_plus: np.ndarray, dt: float, dx: float, dz: float) -> np.ndarray:
    geom, _ = make_time_geometry(
        full_case.geom_0.metric_cov,
        full_case.geom_m.metric_cov,
        metric_plus,
        dt,
        dx,
        dz,
        with_dgamma=False,
    )
    return geom.ricci - 0.5 * full_case.geom_0.metric_cov * geom.r_scalar[..., None, None]


def solve_scaled_min_norm(a: np.ndarray, b: np.ndarray, ridge: float) -> tuple[np.ndarray, dict[str, object]]:
    col_scale = np.linalg.norm(a, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_s = a / col_scale[None, :]
    if ridge > 0.0:
        a_solve = np.vstack([a_s, float(ridge) * np.eye(a_s.shape[1])])
        b_solve = np.concatenate([b, np.zeros(a_s.shape[1], dtype=float)])
    else:
        a_solve = a_s
        b_solve = b
    y, residuals, rank, singular = np.linalg.lstsq(a_solve, b_solve, rcond=1.0e-12)
    x = y / col_scale
    residual = a @ x - b
    return x, {
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "condition_scaled": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_relative_residual": float(np.linalg.norm(residual) / max(np.linalg.norm(b), 1.0e-300)),
        "weighted_residual_norm": float(np.linalg.norm(residual)),
        "weighted_rhs_norm": float(np.linalg.norm(b)),
        "ridge": float(ridge),
        "lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
    }


def render_maps(
    path: Path,
    x: np.ndarray,
    z: np.ndarray,
    mask: np.ndarray,
    before_rel: np.ndarray,
    after_rel: np.ndarray,
    delta_g_norm: np.ndarray,
) -> None:
    x_um = x * 1.973269804593025e-7
    z_um = z * 1.973269804593025e-7
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), constrained_layout=True)
    fields = [
        ("central algebraic residual before", before_rel, "magma"),
        ("central algebraic residual after metric-plus correction", after_rel, "magma"),
        (r"$|\delta \\tilde g_+|$", delta_g_norm, "viridis"),
    ]
    for ax, (title, field, cmap) in zip(axes, fields):
        finite = field[mask & np.isfinite(field)]
        vmax = max(float(np.percentile(finite, 95.0)) if finite.size else 1.0, 1.0e-16)
        im = ax.pcolormesh(xg, zg, np.clip(field, 0.0, vmax), shading="auto", cmap=cmap, vmin=0.0, vmax=vmax)
        ax.contour(xg, zg, mask.astype(float), levels=[0.5], colors="white", linewidths=0.7)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ATOM_NAMES
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.dt_old) * old_scale

    cases = {
        key: build_case(args, ref, float(args.tau) + key * float(args.force_probe_dt_old))
        for key in base.TIME_KEYS
    }
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])

    a_sys, b_sys, meta = base.assemble_system(
        cases,
        full_case,
        region=args.fit_region,
        force_region=args.force_region,
        target_floor_frac=float(args.target_floor_frac),
        force_weight=1.0,
        force_mode="full",
        force_erosion=int(args.force_mask_erosion),
        active_dilation=int(args.active_dilation),
        dt=dt,
        dx=dx,
        dz=dz,
    )
    del a_sys, b_sys
    coeff_path = args.coefficients
    if coeff_path is None:
        if abs(float(args.tau)) < 1.0e-12:
            coeff_path = Path("visualizations/equation_first_gbcd_aux_noq_time_scan_n96_tau0_tw_1em5/gbcd_auxiliary_gauge_coefficients.npz")
        elif float(args.tau) < 0:
            coeff_path = Path("visualizations/equation_first_gbcd_aux_noq_time1em5_n96_taum3p5/gbcd_auxiliary_gauge_coefficients.npz")
        else:
            coeff_path = Path("visualizations/equation_first_gbcd_aux_noq_time1em5_n96_taup3p5/gbcd_auxiliary_gauge_coefficients.npz")
    coeff = load_full_coeff(coeff_path, meta, "time", cases[0].rho_a.shape)
    fields = base.coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
    cov = base.coeff_fields_to_cov(cases, fields)

    mask = getattr(cases[0], args.fit_region)
    idx = np.where(mask.reshape(-1))[0]
    target_rows = symmetric_rows(cases[0].r_need, idx)
    residual_tensor = cov[0] - cases[0].r_need
    residual_rows = symmetric_rows(residual_tensor, idx)
    target_norm = np.sqrt(np.sum(target_rows**2, axis=1))
    finite_target = target_norm[np.isfinite(target_norm)]
    floor = float(args.target_floor_frac) * max(
        float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0,
        1.0e-300,
    )
    weights_point = np.sqrt(np.maximum(cases[0].rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    weights_point = weights_point / np.maximum(target_norm, floor)
    b = (residual_rows * weights_point[:, None]).reshape(-1)

    metric_mask = base.dilate_mask_8(mask, int(args.metric_active_dilation))
    metric_mask[0, :] = False
    metric_mask[-1, :] = False
    metric_mask[:, 0] = False
    metric_mask[:, -1] = False
    variables, _ = metric_plus_variable_index(metric_mask)
    baseline_g = center_einstein_from_metric_plus(full_case, full_case.geom_p.metric_cov, dt, dx, dz)
    jac = np.zeros((b.size, len(variables)), dtype=float)
    eps_base = float(args.perturbation)
    metric_plus_base = np.array(full_case.geom_p.metric_cov, dtype=float, copy=True)
    for col, (i, j, a, c) in enumerate(variables):
        perturbed = metric_plus_base.copy()
        eps = eps_base * max(1.0, abs(float(metric_plus_base[i, j, a, c])))
        perturbed[i, j, a, c] += eps
        if a != c:
            perturbed[i, j, c, a] += eps
        dg = center_einstein_from_metric_plus(full_case, perturbed, dt, dx, dz) - baseline_g
        rows = symmetric_rows(dg, idx)
        jac[:, col] = (rows * weights_point[:, None]).reshape(-1) / eps

    delta, solve = solve_scaled_min_norm(jac, b, float(args.ridge))
    delta_metric_plus = np.zeros_like(metric_plus_base)
    for value, (i, j, a, c) in zip(delta, variables):
        delta_metric_plus[i, j, a, c] += value
        if a != c:
            delta_metric_plus[i, j, c, a] += value
    corrected_g = center_einstein_from_metric_plus(full_case, metric_plus_base + delta_metric_plus, dt, dx, dz)
    corrected_r_need = cases[0].r_need + (corrected_g - baseline_g)

    before_rel = tensor_norm(cov[0] - cases[0].r_need) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    after_rel = tensor_norm(cov[0] - corrected_r_need) / np.maximum(tensor_norm(corrected_r_need), 1.0e-300)
    delta_g_norm = tensor_norm(delta_metric_plus)
    metric_norm = tensor_norm(metric_plus_base)
    weights_eval = np.sqrt(np.maximum(cases[0].rho_a[mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))

    plot_path = args.output / "gbcd_v0_metric_plus_correction.png"
    render_maps(plot_path, cases[0].x, cases[0].z, mask, before_rel, after_rel, delta_g_norm)
    coeff_path_out = args.output / "metric_plus_correction.npz"
    np.savez_compressed(
        coeff_path_out,
        delta_metric_plus=delta_metric_plus,
        metric_plus_base=metric_plus_base,
        before_rel=before_rel,
        after_rel=after_rel,
        mask=mask,
    )
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "metric_active_dilation": int(args.metric_active_dilation),
            "perturbation": float(args.perturbation),
            "ridge": float(args.ridge),
            "coefficients": str(coeff_path.resolve()),
            "dt_ev_inv": float(dt),
            "dx_ev_inv": float(dx),
            "dz_ev_inv": float(dz),
        },
        "definition": {
            "goal": "Linearize the center-slice Einstein tensor with respect to the next metric slice gtilde_+ and solve the smallest delta gtilde_+ that makes the center algebraic gBCD field equation closer to exact.",
            "before_equation": "C_0 = G_0[g_- , g_0, g_+]-T_0/Mp^2.",
            "correction_equation": "delta G_0[delta g_+] ≈ C_0 - (G_0-T_0/Mp^2).",
        },
        "counts": {
            "fit_points": int(np.count_nonzero(mask)),
            "metric_variables": int(len(variables)),
            "rows": int(jac.shape[0]),
        },
        "solve": solve,
        "evaluation": {
            "before_relative_residual": weighted_stats(before_rel[mask], weights_eval),
            "after_relative_residual": weighted_stats(after_rel[mask], weights_eval),
            "delta_metric_plus_norm": weighted_stats(delta_g_norm[metric_mask]),
            "relative_delta_metric_plus_norm": weighted_stats(
                delta_g_norm[metric_mask] / np.maximum(metric_norm[metric_mask], 1.0e-300)
            ),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "correction_npz": str(coeff_path_out.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--coefficients", type=Path, default=None)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
    parser.add_argument("--active-dilation", type=int, default=1)
    parser.add_argument("--metric-active-dilation", type=int, default=1)
    parser.add_argument("--perturbation", type=float, default=1.0e-8)
    parser.add_argument("--ridge", type=float, default=1.0e-10)
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
