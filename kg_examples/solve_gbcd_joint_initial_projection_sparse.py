from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from analyze_bcd_residuals_from_a_reference import stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from coordinate_matter_evolution import solve_covector_time_component
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det
from solve_gbcd_full_linear_metric_update_sparse import (
    COMPONENTS,
    apply_a,
    apply_b_entries,
    apply_bt_entries,
    build_entry_arrays,
    build_sparse_columns,
    load_lambda_cov,
    lsqr_solve,
    metric_variable_index_for_slices,
    parse_metric_variable_slices,
    sym_to_vec,
)
from test_gbcd_v0_one_step_lambda import weighted_stats


JointVariable = Union[tuple[str, int, int, int, int], tuple[str, int, int]]


def stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }
    if weights is not None:
        w = np.asarray(weights, dtype=float)[mask]
        w = np.maximum(w, 0.0)
        if float(np.sum(w)) > 0.0:
            out["weighted_mean"] = float(np.sum(w * vals) / np.sum(w))
    return out


def column_norms_from_column_blocks(col_vals: list[np.ndarray]) -> np.ndarray:
    norms = np.asarray([float(np.linalg.norm(vals)) for vals in col_vals], dtype=float)
    finite = norms[np.isfinite(norms) & (norms > 0.0)]
    floor = 1.0e-30 if finite.size == 0 else max(float(np.percentile(finite, 5.0)) * 1.0e-12, 1.0e-30)
    return np.where(norms > floor, norms, floor)


def add_source_log_measure_columns(
    *,
    col_rows: list[np.ndarray],
    col_vals: list[np.ndarray],
    variables: list[JointVariable],
    source_tilde: np.ndarray,
    row_points: list[tuple[int, int]],
    row_weights: np.ndarray,
) -> tuple[int, int]:
    """Append one local source-amplitude variable per fitted point.

    The linearized equation is
        delta G_mn - delta(T_mn/Mp^2) = residual_mn.
    For a log-measure perturbation eta, with metric and u fixed at the
    linearization point, delta(T/Mp^2)=eta*(T/Mp^2).  Therefore the source
    column enters with a minus sign.
    """
    start = len(variables)
    for n, (i, j) in enumerate(row_points):
        variables.append(("source_log_measure", int(i), int(j)))
        col_rows.append(np.asarray([n], dtype=np.int32))
        col_vals.append((-row_weights[n] * sym_to_vec(source_tilde[i, j]))[None, :])
    return start, len(variables)


def render_summary(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    fit_mask: np.ndarray,
    rho_a: np.ndarray,
    before_rel: np.ndarray,
    linear_after_rel: np.ndarray,
    exact_after_rel: np.ndarray,
    eta: np.ndarray,
    discriminant: np.ndarray,
) -> None:
    extent = [
        float(x[0] * HBAR_C_EV_M * 1.0e6),
        float(x[-1] * HBAR_C_EV_M * 1.0e6),
        float(z[0] * HBAR_C_EV_M * 1.0e6),
        float(z[-1] * HBAR_C_EV_M * 1.0e6),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.8), constrained_layout=True)
    panels = [
        (rho_a, "A-branch rho", "viridis", False),
        (before_rel, "relative residual before", "magma", True),
        (linear_after_rel, "linear joint residual after", "magma", True),
        (exact_after_rel, "exact back-substitution residual", "magma", True),
        (eta, "eta = delta log transformed measure", "coolwarm", False),
        (discriminant, "mass-shell discriminant after", "coolwarm", False),
    ]
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    for ax, (field, title, cmap, log_scale) in zip(axes.flat, panels):
        data = np.asarray(field, dtype=float)
        if log_scale:
            data = np.log10(1.0 + np.maximum(data, 0.0))
            title = "log10(1+ " + title + ")"
        im = ax.imshow(data.T, origin="lower", extent=extent, aspect="equal", cmap=cmap)
        try:
            ax.contour(x_um, z_um, fit_mask.T.astype(float), levels=[0.5], colors=["white"], linewidths=0.7)
        except ValueError:
            pass
        ax.set_title(title)
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def weighted_pullback_l1_from_eta(rho_a: np.ndarray, eta: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(rho_a) & np.isfinite(eta)
    if not np.any(valid):
        return 0.0
    weight = np.maximum(rho_a[valid], 0.0)
    denom = max(float(np.sum(weight * np.abs(rho_a[valid]))), 1.0e-300)
    candidate = rho_a[valid] * np.exp(np.clip(eta[valid], -50.0, 50.0))
    return float(np.sum(weight * np.abs(candidate - rho_a[valid])) / denom)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV

    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    cases = {key: build_case(args, ref, float(args.tau) + key * float(args.force_probe_dt_old)) for key in base.TIME_KEYS}
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])

    coeff_path = args.coefficients
    if coeff_path is None:
        raise ValueError("--coefficients is required for the joint projection prototype")
    cov0, fit_mask = load_lambda_cov(cases, coeff_path, str(args.atom_family))
    residual = cov0 - cases[0].r_need

    active_metric = base.dilate_mask_8(fit_mask, int(args.metric_active_dilation))
    active_metric[0, :] = False
    active_metric[-1, :] = False
    active_metric[:, 0] = False
    active_metric[:, -1] = False
    metric_slices = parse_metric_variable_slices(str(args.metric_variable_slices))
    metric_variables = metric_variable_index_for_slices(active_metric, metric_slices)
    row_points = [(int(i), int(j)) for i, j in np.argwhere(fit_mask)]

    target_norm = np.asarray([max(float(tensor_norm(cases[0].r_need)[i, j]), 1.0e-300) for i, j in row_points])
    finite_target = target_norm[np.isfinite(target_norm)]
    floor = float(args.target_floor_frac) * max(float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0, 1.0e-300)
    row_weights = np.asarray(
        [
            np.sqrt(max(cases[0].rho_a[i, j], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
            / max(target_norm[n], floor)
            for n, (i, j) in enumerate(row_points)
        ],
        dtype=float,
    )
    rhs = np.concatenate([row_weights[n] * sym_to_vec(residual[i, j]) for n, (i, j) in enumerate(row_points)])

    dg_base, d2g_base = metric_jets_full(
        full_case.geom_m.metric_cov,
        full_case.geom_0.metric_cov,
        full_case.geom_p.metric_cov,
        dt,
        dx,
        dz,
    )
    col_rows, col_vals, _old_col_norms = build_sparse_columns(
        metric_cov=full_case.geom_0.metric_cov,
        dg_base=dg_base,
        d2g_base=d2g_base,
        variables=metric_variables,
        row_points=row_points,
        row_weights=row_weights,
        dt=dt,
        dx=dx,
        dz=dz,
        eps=float(args.linear_eps),
    )
    variables: list[JointVariable] = list(metric_variables)
    source_start = source_end = len(variables)
    if args.source_variable_mode == "log_measure":
        source_start, source_end = add_source_log_measure_columns(
            col_rows=col_rows,
            col_vals=col_vals,
            variables=variables,
            source_tilde=cases[0].source_tilde,
            row_points=row_points,
            row_weights=row_weights,
        )
    elif args.source_variable_mode != "none":
        raise ValueError(f"unknown source variable mode: {args.source_variable_mode}")

    col_norms = column_norms_from_column_blocks(col_vals)
    entry_row, entry_col, entry_val, entry_n_rows = build_entry_arrays(col_rows, col_vals, col_norms, len(row_points))

    def b_apply(y: np.ndarray) -> np.ndarray:
        return apply_b_entries(entry_row, entry_col, entry_val, y, entry_n_rows)

    def bt_apply(v: np.ndarray) -> np.ndarray:
        return apply_bt_entries(entry_row, entry_col, entry_val, v, len(variables))

    nvars = len(variables)
    metric_reg = np.zeros(nvars, dtype=float)
    metric_reg[: len(metric_variables)] = float(args.metric_ridge)
    source_reg = np.zeros(nvars, dtype=float)
    source_close = np.zeros(nvars, dtype=float)
    if source_end > source_start:
        source_reg[source_start:source_end] = float(args.source_ridge)
        for col, var in enumerate(variables[source_start:source_end], start=source_start):
            _, i, j = var
            rho_w = np.sqrt(max(float(cases[0].rho_a[i, j]), 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
            source_close[col] = np.sqrt(max(float(args.source_closeness_weight), 0.0)) * rho_w / col_norms[col]

    def aprod_aug(y: np.ndarray) -> np.ndarray:
        return np.concatenate(
            [
                b_apply(y),
                metric_reg * y,
                source_reg * y,
                source_close * y,
            ]
        )

    def atprod_aug(u: np.ndarray) -> np.ndarray:
        m0 = entry_n_rows
        m1 = m0 + nvars
        m2 = m1 + nvars
        m3 = m2 + nvars
        return (
            bt_apply(u[:m0])
            + metric_reg * u[m0:m1]
            + source_reg * u[m1:m2]
            + source_close * u[m2:m3]
        )

    y, solve_info = lsqr_solve(
        aprod_aug,
        atprod_aug,
        np.concatenate([rhs, np.zeros(nvars * 3, dtype=float)]),
        nvars,
        tol=float(args.lsqr_tol),
        maxiter=int(args.lsqr_maxiter),
    )
    delta_vec = y / col_norms
    delta_vec_eval = delta_vec.copy()

    delta_metrics = {
        "minus": np.zeros_like(full_case.geom_m.metric_cov),
        "center": np.zeros_like(full_case.geom_0.metric_cov),
        "plus": np.zeros_like(full_case.geom_p.metric_cov),
    }
    eta = np.zeros_like(cases[0].rho_a, dtype=float)
    eta_raw = np.zeros_like(cases[0].rho_a, dtype=float)
    eta_clip_count = 0
    for col, (value, var) in enumerate(zip(delta_vec, variables)):
        if var[0] == "source_log_measure":
            _, i, j = var
            raw_value = float(value)
            clipped_value = float(np.clip(raw_value, -float(args.eta_clip), float(args.eta_clip)))
            eta_raw[int(i), int(j)] = raw_value
            eta[int(i), int(j)] = clipped_value
            delta_vec_eval[col] = clipped_value
            if clipped_value != raw_value:
                eta_clip_count += 1
            continue
        slice_name, i, j, a, b = var
        delta_metrics[str(slice_name)][int(i), int(j), int(a), int(b)] += value
        if int(a) != int(b):
            delta_metrics[str(slice_name)][int(i), int(j), int(b), int(a)] += value
    weighted_pred = apply_a(col_rows, col_vals, delta_vec_eval, len(row_points))

    corrected_metric_minus = full_case.geom_m.metric_cov + delta_metrics["minus"]
    corrected_metric_center = full_case.geom_0.metric_cov + delta_metrics["center"]
    corrected_metric_plus = full_case.geom_p.metric_cov + delta_metrics["plus"]
    geom_corr, _ = make_time_geometry(
        corrected_metric_center,
        corrected_metric_minus,
        corrected_metric_plus,
        dt,
        dx,
        dz,
        with_dgamma=False,
    )
    ein_corr = geom_corr.ricci - 0.5 * corrected_metric_center * geom_corr.r_scalar[..., None, None]

    measure_corr = cases[0].measure_tilde * np.exp(np.clip(eta, -float(args.eta_clip), float(args.eta_clip)))
    sqrt_abs_g_corr = safe_sqrt_abs_det(corrected_metric_center)
    rho_tilde_corr = np.maximum(measure_corr / np.maximum(sqrt_abs_g_corr, 1.0e-300), float(args.rho_floor))
    phase_corr = solve_covector_time_component(
        geom_corr.metric_inv,
        cases[0].u_cov[..., 1],
        cases[0].u_cov[..., 2],
        float(ref["params"].m),
        branch="negative_frequency",
        u_t_reference=cases[0].u_cov[..., 0],
    )
    u_t_corr = phase_corr["u_t"]
    stress_corr, tilde_x_corr = stress_tensor_tilde(
        corrected_metric_center,
        geom_corr.metric_inv,
        rho_tilde_corr,
        u_t_corr,
        cases[0].u_cov[..., 1],
        cases[0].u_cov[..., 2],
        m=float(ref["params"].m),
    )
    source_corr = stress_corr / (float(args.mp) * float(args.mp))
    exact_r_need = ein_corr - source_corr

    predicted = np.zeros_like(residual)
    width = len(COMPONENTS)
    for n, (i, j) in enumerate(row_points):
        vals = weighted_pred[n * width : (n + 1) * width] / max(row_weights[n], 1.0e-300)
        for value, (a, b) in zip(vals, COMPONENTS):
            predicted[i, j, a, b] = value
            predicted[i, j, b, a] = value

    before_rel = tensor_norm(residual) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    linear_after_rel = tensor_norm(residual - predicted) / np.maximum(tensor_norm(cases[0].r_need + predicted), 1.0e-300)
    exact_after_rel = tensor_norm(cov0 - exact_r_need) / np.maximum(tensor_norm(exact_r_need), 1.0e-300)
    metric_norm_center = tensor_norm(full_case.geom_0.metric_cov)
    metric_norm_plus = tensor_norm(full_case.geom_p.metric_cov)
    delta_norm_center = tensor_norm(delta_metrics["center"])
    delta_norm_plus = tensor_norm(delta_metrics["plus"])

    weights_eval = np.sqrt(np.maximum(cases[0].rho_a[fit_mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    plot_path = args.output / "gbcd_joint_initial_projection_sparse.png"
    render_summary(
        plot_path,
        cases[0].x,
        cases[0].z,
        fit_mask,
        cases[0].rho_a,
        before_rel,
        linear_after_rel,
        exact_after_rel,
        eta,
        phase_corr["discriminant"],
    )
    data_path = args.output / "joint_initial_projection_sparse.npz"
    np.savez_compressed(
        data_path,
        delta_metric_minus=delta_metrics["minus"],
        delta_metric_center=delta_metrics["center"],
        delta_metric_plus=delta_metrics["plus"],
        corrected_metric_minus=corrected_metric_minus,
        corrected_metric_center=corrected_metric_center,
        corrected_metric_plus=corrected_metric_plus,
        eta_log_measure=eta,
        corrected_measure=measure_corr,
        corrected_rho_tilde=rho_tilde_corr,
        corrected_u_t=u_t_corr,
        before_rel=before_rel,
        linear_after_rel=linear_after_rel,
        exact_after_rel=exact_after_rel,
        eta_raw_log_measure=eta_raw,
        mass_shell_defect=tilde_x_corr - float(ref["params"].m) ** 2,
        discriminant=phase_corr["discriminant"],
        mask=fit_mask,
    )

    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "metric_variable_slices": str(args.metric_variable_slices),
            "source_variable_mode": str(args.source_variable_mode),
            "metric_ridge": float(args.metric_ridge),
            "source_ridge": float(args.source_ridge),
            "source_closeness_weight": float(args.source_closeness_weight),
            "dt_ev_inv": float(dt),
            "dx_ev_inv": float(dx),
            "dz_ev_inv": float(dz),
            "coefficients": str(coeff_path.resolve()),
        },
        "definition": {
            "goal": "Joint initial projection for equation-first gBCD: solve one sparse linearized system for metric jets and transformed-measure source amplitude.",
            "linearized_equation": "delta G_mn - delta(T_mn/Mp^2) = C_mn - (G_mn - T_mn/Mp^2).",
            "source_variable": "eta = delta log(sqrt(|gtilde|) rho_tilde). After solving, the script recomputes rho_tilde from the corrected metric determinant and re-solves the mass shell for u_t.",
            "caveat": "This is a first joint-projection prototype. It enforces the D tensor equation by linearized least squares, then checks exact nonlinear back-substitution; it does not yet impose continuity/integrability as hard rows.",
        },
        "counts": {
            "fit_points": int(np.count_nonzero(fit_mask)),
            "metric_active_points": int(np.count_nonzero(active_metric)),
            "metric_variables": int(len(metric_variables)),
            "source_variables": int(source_end - source_start),
            "total_variables": int(nvars),
            "field_rows": int(entry_n_rows),
            "scalar_nonzeros": int(entry_val.size),
        },
        "solve": solve_info,
        "evaluation": {
            "before_relative_residual": weighted_stats(before_rel[fit_mask], weights_eval),
            "linear_after_relative_residual": weighted_stats(linear_after_rel[fit_mask], weights_eval),
            "exact_after_relative_residual": weighted_stats(exact_after_rel[fit_mask], weights_eval),
            "eta_raw_abs": stats(np.abs(eta_raw[fit_mask]), weights_eval),
            "eta_abs": stats(np.abs(eta[fit_mask]), weights_eval),
            "eta_clip_fraction": float(eta_clip_count / max(source_end - source_start, 1)),
            "rho_pullback_weighted_relative_l1_linearized": weighted_pullback_l1_from_eta(cases[0].rho_a, eta, fit_mask),
            "relative_delta_metric_center_norm": weighted_stats(
                delta_norm_center[active_metric] / np.maximum(metric_norm_center[active_metric], 1.0e-300)
            ),
            "relative_delta_metric_plus_norm": weighted_stats(
                delta_norm_plus[active_metric] / np.maximum(metric_norm_plus[active_metric], 1.0e-300)
            ),
            "mass_shell_defect_abs": stats(np.abs((tilde_x_corr - float(ref["params"].m) ** 2)[fit_mask]), weights_eval),
            "mass_shell_discriminant": stats(phase_corr["discriminant"][fit_mask], weights_eval),
            "negative_discriminant_fraction": float(np.count_nonzero(phase_corr["discriminant"][fit_mask] < 0.0) / max(np.count_nonzero(fit_mask), 1)),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "data_npz": str(data_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--metric-active-dilation", type=int, default=1)
    parser.add_argument("--metric-variable-slices", choices=["plus", "center_plus", "minus_center_plus"], default="center_plus")
    parser.add_argument("--source-variable-mode", choices=["none", "log_measure"], default="log_measure")
    parser.add_argument("--linear-eps", type=float, default=1.0e-8)
    parser.add_argument("--metric-ridge", type=float, default=1.0e-6)
    parser.add_argument("--source-ridge", type=float, default=1.0e-8)
    parser.add_argument("--source-closeness-weight", type=float, default=1.0)
    parser.add_argument("--eta-clip", type=float, default=2.0)
    parser.add_argument("--lsqr-tol", type=float, default=1.0e-6)
    parser.add_argument("--lsqr-maxiter", type=int, default=1000)
    parser.add_argument("--atom-family", choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"], default="auto")
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
