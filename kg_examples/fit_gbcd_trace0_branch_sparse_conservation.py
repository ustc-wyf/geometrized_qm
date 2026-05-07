from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_equation_first_transverse_force import (
    covector_norm_euclidean,
    divergence_cov2,
    derivative_tensor,
    transverse_covector,
)
from diagnose_gbcd_trace_local_closure import ATOM_NAMES, trace_coefficients, weighted_stats
from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_gbcd_principal_constraint_projection_sparse import append_row
from fit_gbcd_trace0_sparse_conservation import (
    TIME_KEYS,
    add_projected,
    augmented_lagrangian_update,
    build_basis,
    build_column_index,
    coeff_to_lambda_fields,
    make_active_mask,
)
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def scalar_contract(metric_inv: np.ndarray, a_cov: np.ndarray, b_cov: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...a,...b->...", metric_inv, a_cov, b_cov, optimize=True)


def raise_covector(metric_inv: np.ndarray, cov: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...b->...a", metric_inv, cov, optimize=True)


def qtilde_from_triplet(case_m, case_0, case_p, gamma: np.ndarray, dt: float, dx: float, dz: float) -> np.ndarray:
    r_up_m = raise_covector(case_m.metric_inv, case_m.r_cov)
    r_up_0 = raise_covector(case_0.metric_inv, case_0.r_cov)
    r_up_p = raise_covector(case_p.metric_inv, case_p.r_cov)
    div = (r_up_p[..., 0] - r_up_m[..., 0]) / (2.0 * dt)
    div = div + np.gradient(r_up_0[..., 1], dx, axis=0, edge_order=2)
    div = div + np.gradient(r_up_0[..., 2], dz, axis=1, edge_order=2)
    conn = np.zeros_like(div)
    for mu in range(3):
        for lam in range(3):
            conn += gamma[..., mu, mu, lam] * r_up_0[..., lam]
    r2 = scalar_contract(case_0.metric_inv, case_0.r_cov, case_0.r_cov)
    return div + conn + r2


def build_chi_fields(
    args: argparse.Namespace,
    ref: dict[str, object],
    cases: dict[int, object],
    tau: float,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[dict[int, np.ndarray], dict[str, object]]:
    probe = float(args.force_probe_dt_old)
    ext = dict(cases)
    ext[-2] = build_case(args, ref, tau - 2.0 * probe)
    ext[2] = build_case(args, ref, tau + 2.0 * probe)
    full = {key: build_full_case(args, ref, tau + key * probe) for key in TIME_KEYS}
    m2 = max(float(ref["params"].m) ** 2, 1.0e-300)
    q = {
        key: qtilde_from_triplet(ext[key - 1], ext[key], ext[key + 1], full[key].geom_0.gamma, dt, dx, dz) / m2
        for key in TIME_KEYS
    }
    q_abs_values = []
    for key in TIME_KEYS:
        mask = getattr(cases[key], args.fit_region)
        vals = np.abs(q[key][mask & np.isfinite(q[key])])
        if vals.size:
            q_abs_values.append(vals)
    if q_abs_values:
        q_all = np.concatenate(q_abs_values)
        q_scale = float(np.percentile(q_all, float(args.branch_q_scale_percentile)))
    else:
        q_scale = 1.0
    q_scale = max(q_scale, 1.0e-300)
    chi = {}
    for key in TIME_KEYS:
        q_rel = np.abs(q[key]) / q_scale
        chi[key] = q_rel * q_rel / (q_rel * q_rel + float(args.branch_q0) ** 2)
        chi[key] = np.maximum(chi[key], float(args.branch_chi_floor))
    return chi, {
        "q_scale": float(q_scale),
        "q_scale_percentile": float(args.branch_q_scale_percentile),
        "branch_q0": float(args.branch_q0),
        "branch_chi_floor": float(args.branch_chi_floor),
        "q_abs_stats": {
            str(key): weighted_stats(
                np.abs(q[key][getattr(cases[key], args.fit_region)]),
                np.sqrt(
                    np.maximum(cases[key].rho_a[getattr(cases[key], args.fit_region)], 0.0)
                    / max(float(np.max(cases[key].rho_a)), 1.0e-300)
                ),
            )
            for key in TIME_KEYS
        },
        "chi_stats": {
            str(key): weighted_stats(
                chi[key][getattr(cases[key], args.fit_region)],
                np.sqrt(
                    np.maximum(cases[key].rho_a[getattr(cases[key], args.fit_region)], 0.0)
                    / max(float(np.max(cases[key].rho_a)), 1.0e-300)
                ),
            )
            for key in TIME_KEYS
        },
    }


def scaled_atoms(case, chi: np.ndarray) -> dict[str, np.ndarray]:
    atoms = tensor_atoms(case, "matter4")
    return {name: chi[..., None, None] * value for name, value in atoms.items()}


def add_branch_algebraic_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    chi_fields: dict[int, np.ndarray],
    col_index: np.ndarray,
    basis: dict[int, np.ndarray],
    region: str,
    target_floor_frac: float,
) -> int:
    count = 0
    atom_cache = {key: scaled_atoms(cases[key], chi_fields[key]) for key in TIME_KEYS}
    for tpos, key in enumerate(TIME_KEYS):
        case = cases[key]
        mask = getattr(case, region)
        idx = np.where(mask.reshape(-1))[0]
        target_rows = symmetric_rows(case.r_need, idx)
        weights = base.relative_row_weight(case, idx, target_rows, target_floor_frac)
        for p, (i_raw, j_raw) in enumerate(np.argwhere(mask)):
            i = int(i_raw)
            j = int(j_raw)
            for comp_pos, (mu, nu) in enumerate(SYMMETRIC_COMPONENTS):
                entries: list[tuple[int, float]] = []
                for apos, atom_name in enumerate(ATOM_NAMES):
                    value = weights[p] * atom_cache[key][atom_name][i, j, mu, nu]
                    add_projected(entries, col_index, basis, tpos, key, apos, i, j, value)
                append_row(row_cols, row_vals, rhs, entries, float(weights[p] * target_rows[p, comp_pos]))
                count += 1
    return count


def add_branch_bound_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    col_index: np.ndarray,
    basis: dict[int, np.ndarray],
    region: str,
    target_floor_frac: float,
    bound_weight: float,
) -> int:
    if bound_weight <= 0.0:
        return 0
    count = 0
    atom_cache = {key: tensor_atoms(cases[key], "matter4") for key in TIME_KEYS}
    for tpos, key in enumerate(TIME_KEYS):
        case = cases[key]
        mask = getattr(case, region)
        idx = np.where(mask.reshape(-1))[0]
        target_rows = symmetric_rows(case.r_need, idx)
        target_norm = np.sqrt(np.sum(target_rows**2, axis=1))
        finite = target_norm[np.isfinite(target_norm)]
        target_scale = max(
            float(target_floor_frac) * (float(np.percentile(finite, 95.0)) if finite.size else 0.0),
            1.0e-300,
        )
        rho_w = np.sqrt(np.maximum(case.rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
        for p, (i_raw, j_raw) in enumerate(np.argwhere(mask)):
            i = int(i_raw)
            j = int(j_raw)
            w = float(bound_weight) * float(rho_w[p]) / target_scale
            for mu, nu in SYMMETRIC_COMPONENTS:
                entries: list[tuple[int, float]] = []
                for apos, atom_name in enumerate(ATOM_NAMES):
                    value = w * atom_cache[key][atom_name][i, j, mu, nu]
                    add_projected(entries, col_index, basis, tpos, key, apos, i, j, value)
                append_row(row_cols, row_vals, rhs, entries, 0.0)
                count += 1
    return count


def add_branch_force_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    chi_fields: dict[int, np.ndarray],
    full_gamma_case,
    col_index: np.ndarray,
    basis: dict[int, np.ndarray],
    force_region: str,
    target_floor_frac: float,
    force_erosion: int,
    dt: float,
    dx: float,
    dz: float,
) -> int:
    atom_cache = {key: scaled_atoms(cases[key], chi_fields[key]) for key in TIME_KEYS}
    force_mask = base.make_force_mask(cases[0], force_region, force_erosion)
    force_weights = base.force_row_weight(cases[0], force_mask, min(abs(dt), abs(dx), abs(dz)), target_floor_frac, 1.0)
    metric_inv = full_gamma_case.geom_0.metric_inv
    gamma = full_gamma_case.geom_0.gamma
    count = 0
    for p, (i_raw, j_raw) in enumerate(np.argwhere(force_mask)):
        i = int(i_raw)
        j = int(j_raw)
        for sigma in range(3):
            entries: list[tuple[int, float]] = []
            for mu in range(3):
                for alpha in range(3):
                    pref = metric_inv[i, j, mu, alpha]
                    if not np.isfinite(pref):
                        continue
                    if alpha == 0:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            add_projected(
                                entries,
                                col_index,
                                basis,
                                2,
                                1,
                                apos,
                                i,
                                j,
                                pref * atom_cache[1][atom_name][i, j, mu, sigma] / (2.0 * dt),
                            )
                            add_projected(
                                entries,
                                col_index,
                                basis,
                                0,
                                -1,
                                apos,
                                i,
                                j,
                                -pref * atom_cache[-1][atom_name][i, j, mu, sigma] / (2.0 * dt),
                            )
                    elif alpha == 1:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            add_projected(
                                entries,
                                col_index,
                                basis,
                                1,
                                0,
                                apos,
                                i + 1,
                                j,
                                pref * atom_cache[0][atom_name][i + 1, j, mu, sigma] / (2.0 * dx),
                            )
                            add_projected(
                                entries,
                                col_index,
                                basis,
                                1,
                                0,
                                apos,
                                i - 1,
                                j,
                                -pref * atom_cache[0][atom_name][i - 1, j, mu, sigma] / (2.0 * dx),
                            )
                    else:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            add_projected(
                                entries,
                                col_index,
                                basis,
                                1,
                                0,
                                apos,
                                i,
                                j + 1,
                                pref * atom_cache[0][atom_name][i, j + 1, mu, sigma] / (2.0 * dz),
                            )
                            add_projected(
                                entries,
                                col_index,
                                basis,
                                1,
                                0,
                                apos,
                                i,
                                j - 1,
                                -pref * atom_cache[0][atom_name][i, j - 1, mu, sigma] / (2.0 * dz),
                            )
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, mu]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom_name in enumerate(ATOM_NAMES):
                                add_projected(
                                    entries,
                                    col_index,
                                    basis,
                                    1,
                                    0,
                                    apos,
                                    i,
                                    j,
                                    conn * atom_cache[0][atom_name][i, j, lam, sigma],
                                )
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, sigma]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom_name in enumerate(ATOM_NAMES):
                                add_projected(
                                    entries,
                                    col_index,
                                    basis,
                                    1,
                                    0,
                                    apos,
                                    i,
                                    j,
                                    conn * atom_cache[0][atom_name][i, j, mu, lam],
                                )
            entries = [(c, float(force_weights[p]) * v) for c, v in entries]
            append_row(row_cols, row_vals, rhs, entries, 0.0)
            count += 1
    return count


def branch_tensor_fields(cases: dict[int, object], lambda_fields: dict[int, np.ndarray], chi_fields: dict[int, np.ndarray]) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray]]:
    c_out: dict[int, np.ndarray] = {}
    hat_out: dict[int, np.ndarray] = {}
    for key in TIME_KEYS:
        atoms = tensor_atoms(cases[key], "matter4")
        hat = np.zeros_like(cases[key].r_need)
        for apos, atom_name in enumerate(ATOM_NAMES):
            hat += lambda_fields[key][..., apos, None, None] * atoms[atom_name]
        hat_out[key] = hat
        c_out[key] = chi_fields[key][..., None, None] * hat
    return c_out, hat_out


def evaluate_branch_solution(
    cases: dict[int, object],
    chi_fields: dict[int, np.ndarray],
    full_gamma_case,
    coeff: np.ndarray,
    col_index: np.ndarray,
    basis: dict[int, np.ndarray],
    region: str,
    force_region: str,
    force_erosion: int,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, object]:
    lambda_fields = coeff_to_lambda_fields(coeff, col_index, basis, cases[0].rho_a.shape)
    cov, hat = branch_tensor_fields(cases, lambda_fields, chi_fields)
    by_time = {}
    for key in TIME_KEYS:
        mask = getattr(cases[key], region)
        err_norm = tensor_norm(cov[key] - cases[key].r_need)
        target_norm = tensor_norm(cases[key].r_need)
        idx = np.where(mask.reshape(-1))[0]
        err_rows = symmetric_rows(cov[key] - cases[key].r_need, idx)
        target_rows = symmetric_rows(cases[key].r_need, idx)
        err_row_norm = np.sqrt(np.sum(err_rows**2, axis=1))
        target_row_norm = np.sqrt(np.sum(target_rows**2, axis=1))
        finite = target_row_norm[np.isfinite(target_row_norm)]
        row_floor = 0.05 * max(float(np.percentile(finite, 95.0)) if finite.size else 0.0, 1.0e-300)
        rel_row_floor = err_row_norm / np.maximum(target_row_norm, row_floor)
        weights = np.sqrt(np.maximum(cases[key].rho_a[mask], 0.0) / max(float(np.max(cases[key].rho_a)), 1.0e-300))
        trace_value = np.einsum("...ab,...ab->...", cases[key].metric_inv, cov[key], optimize=True)
        by_time[str(key)] = {
            "tau": float(cases[key].tau_old),
            "algebraic_row_floor_relative_residual": weighted_stats(rel_row_floor, weights),
            "trace_abs": weighted_stats(trace_value[mask], weights),
            "C_norm": weighted_stats(tensor_norm(cov[key])[mask], weights),
            "Chat_norm": weighted_stats(tensor_norm(hat[key])[mask], weights),
            "chi": weighted_stats(chi_fields[key][mask], weights),
        }
    dcov = derivative_tensor(cov[-1], cov[0], cov[1], dt, dx, dz)
    j_cov = divergence_cov2(cov[0], dcov, full_gamma_case.geom_0.metric_inv, full_gamma_case.geom_0.gamma)
    j_perp = transverse_covector(j_cov, cases[0].u_cov, full_gamma_case.geom_0.metric_inv)
    total = covector_norm_euclidean(j_cov)
    perp = covector_norm_euclidean(j_perp)
    c_norm = tensor_norm(cov[0])
    h = min(abs(dt), abs(dx), abs(dz))
    force_mask = base.make_force_mask(cases[0], force_region, force_erosion)
    weights = np.sqrt(np.maximum(cases[0].rho_a[force_mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    return {
        "by_time": by_time,
        "central_divergence": {
            "force_count": int(np.count_nonzero(force_mask)),
            "full_over_derivative_scale": weighted_stats(
                total[force_mask] / np.maximum(c_norm[force_mask] / max(h, 1.0e-300), 1.0e-300),
                weights,
            ),
            "transverse_over_derivative_scale": weighted_stats(
                perp[force_mask] / np.maximum(c_norm[force_mask] / max(h, 1.0e-300), 1.0e-300),
                weights,
            ),
        },
    }


def render(path: Path, records: list[dict[str, object]]) -> None:
    labels = [str(r["label"]) for r in records]
    alg = [float(r["evaluation"]["by_time"]["0"]["algebraic_row_floor_relative_residual"]["weighted_mean"]) for r in records]
    div = [float(r["evaluation"]["central_divergence"]["full_over_derivative_scale"]["weighted_mean"]) for r in records]
    chat = [float(r["evaluation"]["by_time"]["0"]["Chat_norm"]["p95"]) for r in records]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), constrained_layout=True)
    axes[0].plot(x, alg, "o-")
    axes[0].set_yscale("log")
    axes[0].set_title("algebraic residual")
    axes[1].plot(x, div, "o-")
    axes[1].set_yscale("log")
    axes[1].set_title("full divergence scale")
    axes[2].plot(x, chat, "o-")
    axes[2].set_yscale("log")
    axes[2].set_title("p95 |Chat|")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ATOM_NAMES
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    tau = float(args.tau)
    probe = float(args.force_probe_dt_old)
    cases = {
        -1: build_case(args, ref, tau - probe),
        0: build_case(args, ref, tau),
        1: build_case(args, ref, tau + probe),
    }
    full = build_full_case(args, ref, tau)
    scale = ref["scale"]
    dt = probe * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    chi_fields, chi_report = build_chi_fields(args, ref, cases, tau, dt, dx, dz)
    active = make_active_mask(cases, args.fit_region, args.force_region, int(args.active_dilation))
    col_index, ncols = build_column_index(active)
    basis = build_basis(cases, active)
    alg_cols: list[np.ndarray] = []
    alg_vals: list[np.ndarray] = []
    alg_rhs: list[float] = []
    force_cols: list[np.ndarray] = []
    force_vals: list[np.ndarray] = []
    force_rhs: list[float] = []
    n_alg = add_branch_algebraic_rows(
        alg_cols,
        alg_vals,
        alg_rhs,
        cases,
        chi_fields,
        col_index,
        basis,
        args.fit_region,
        float(args.target_floor_frac),
    )
    n_force = add_branch_force_rows(
        force_cols,
        force_vals,
        force_rhs,
        cases,
        chi_fields,
        full,
        col_index,
        basis,
        args.force_region,
        float(args.target_floor_frac),
        int(args.force_mask_erosion),
        dt,
        dx,
        dz,
    )
    initial = np.zeros(ncols, dtype=float)
    records = []
    coeff_arrays = {}
    for bound_weight in [float(x.strip()) for x in args.branch_bound_weights.split(",") if x.strip()]:
        row_cols = list(alg_cols)
        row_vals = list(alg_vals)
        rhs = list(alg_rhs)
        n_bound = add_branch_bound_rows(
            row_cols,
            row_vals,
            rhs,
            cases,
            col_index,
            basis,
            args.fit_region,
            float(args.target_floor_frac),
            float(bound_weight),
        )
        coeff, solve = augmented_lagrangian_update(
            initial,
            row_cols,
            row_vals,
            rhs,
            force_cols,
            force_vals,
            ncols,
            mu0=float(args.hard_alm_mu0),
            growth=float(args.hard_alm_growth),
            outer=int(args.hard_alm_outer),
            tol=float(args.hard_alm_tol),
            maxiter=int(args.hard_alm_maxiter),
        )
        evaluation = evaluate_branch_solution(
            cases,
            chi_fields,
            full,
            coeff,
            col_index,
            basis,
            args.fit_region,
            args.force_region,
            int(args.force_mask_erosion),
            dt,
            dx,
            dz,
        )
        label = f"bound:{bound_weight:g}"
        records.append(
            {
                "label": label,
                "branch_bound_weight": float(bound_weight),
                "n_bound_rows": int(n_bound),
                "solve": solve,
                "evaluation": evaluation,
            }
        )
        coeff_arrays[f"coeff_branch_bound_{bound_weight:g}"] = coeff
    plot_path = args.output / "trace0_branch_conservation_scan.png"
    coeff_path = args.output / "trace0_branch_conservation_coefficients.npz"
    render(plot_path, records)
    np.savez_compressed(coeff_path, **coeff_arrays)
    report = {
        "parameters": {
            "tau": tau,
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "force_probe_dt_old": float(args.force_probe_dt_old),
            "ncols": int(ncols),
            "n_algebraic_rows": int(n_alg),
            "n_force_rows": int(n_force),
            "active_count": int(np.count_nonzero(active)),
            "branch_bound_weights": [float(x.strip()) for x in args.branch_bound_weights.split(",") if x.strip()],
            "branch_definition": "Solve for Chat in C=chi(q) Chat, with trace0 imposed on Chat/C through the same nullspace and full conservation imposed on C.",
            "bound_rows": "Optional L2 penalty on Chat tensor rows. This is a numerical boundedness test, not a new physical term by itself.",
        },
        "chi": chi_report,
        "records": records,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "coefficients_npz": str(coeff_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
    parser.add_argument("--active-dilation", type=int, default=1)
    parser.add_argument("--branch-q0", type=float, default=0.1)
    parser.add_argument("--branch-q-scale-percentile", type=float, default=95.0)
    parser.add_argument("--branch-chi-floor", type=float, default=1.0e-8)
    parser.add_argument("--branch-bound-weights", type=str, default="0,1e-4,1e-3,1e-2")
    parser.add_argument("--hard-alm-mu0", type=float, default=100.0)
    parser.add_argument("--hard-alm-growth", type=float, default=10.0)
    parser.add_argument("--hard-alm-outer", type=int, default=3)
    parser.add_argument("--hard-alm-tol", type=float, default=1.0e-6)
    parser.add_argument("--hard-alm-maxiter", type=int, default=300)
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
    run(parser.parse_args())


if __name__ == "__main__":
    main()
