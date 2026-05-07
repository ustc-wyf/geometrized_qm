from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_gbcd_metric_principal_correction import (
    COMPONENTS,
    local_principal_matrix,
    sym_to_vec,
    vec_to_sym,
)
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from test_gbcd_v0_one_step_lambda import ATOM_NAMES, weighted_stats


def add_coeff(row: np.ndarray, col_index: np.ndarray, tpos: int, apos: int, i: int, j: int, value: float) -> None:
    c = int(col_index[tpos, apos, i, j])
    if c >= 0 and np.isfinite(value):
        row[c] += float(value)


def atom_scales_from_coeff_guess(coeff: np.ndarray, meta: dict[str, object]) -> dict[str, float]:
    col = np.asarray(meta["column_index"], dtype=int)
    scales: dict[str, float] = {}
    for apos, atom in enumerate(ATOM_NAMES):
        vals = []
        for tpos in range(3):
            ids = col[tpos, apos]
            mask = ids >= 0
            vals.append(coeff[ids[mask]])
        arr = np.concatenate(vals) if vals else np.asarray([1.0])
        finite = arr[np.isfinite(arr)]
        scales[atom] = max(float(np.percentile(np.abs(finite), 95.0)) if finite.size else 1.0, 1.0e-30)
    return scales


def add_regularizer_rows(
    meta: dict[str, object],
    scales: dict[str, float],
    *,
    norm_weight: float,
    time_weight: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    col = np.asarray(meta["column_index"], dtype=int)
    active = np.asarray(meta["active_mask"], dtype=bool)
    ncols = int(meta["ncols"])
    rows: list[np.ndarray] = []
    rhs: list[float] = []
    counts = {"norm_rows": 0, "time_rows": 0}
    sqrt_norm = float(np.sqrt(max(norm_weight, 0.0)))
    sqrt_time = float(np.sqrt(max(time_weight, 0.0)))
    for apos, atom in enumerate(ATOM_NAMES):
        scale = scales[atom]
        for i, j in np.argwhere(active):
            i = int(i)
            j = int(j)
            for tpos in range(3):
                c = int(col[tpos, apos, i, j])
                if c >= 0 and sqrt_norm > 0.0:
                    row = np.zeros(ncols, dtype=float)
                    row[c] = sqrt_norm / scale
                    rows.append(row)
                    rhs.append(0.0)
                    counts["norm_rows"] += 1
            if sqrt_time > 0.0:
                for t0, t1 in [(0, 1), (1, 2)]:
                    c0 = int(col[t0, apos, i, j])
                    c1 = int(col[t1, apos, i, j])
                    if c0 >= 0 and c1 >= 0:
                        row = np.zeros(ncols, dtype=float)
                        row[c1] = sqrt_time / scale
                        row[c0] = -sqrt_time / scale
                        rows.append(row)
                        rhs.append(0.0)
                        counts["time_rows"] += 1
    if not rows:
        return np.zeros((0, ncols), dtype=float), np.zeros(0, dtype=float), counts
    return np.vstack(rows), np.asarray(rhs, dtype=float), counts


def split_alg_force(a: np.ndarray, b: np.ndarray, meta: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_force = 3 * int(meta["force_count"])
    if n_force <= 0:
        return a, b, np.zeros((0, a.shape[1]), dtype=float), np.zeros(0, dtype=float)
    return a[:-n_force], b[:-n_force], a[-n_force:], b[-n_force:]


def principal_constraint_rows(cases: dict[int, object], meta: dict[str, object], region: str, dt: float) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    col_index = np.asarray(meta["column_index"], dtype=int)
    mask = getattr(cases[0], region)
    atoms = tensor_atoms(cases[0], "matter4")
    rows: list[np.ndarray] = []
    rhs: list[float] = []
    ranks = []
    for i_raw, j_raw in np.argwhere(mask):
        i = int(i_raw)
        j = int(j_raw)
        mat = local_principal_matrix(cases[0].metric_cov[i, j], cases[0].metric_inv[i, j], dt)
        u, singular, _ = np.linalg.svd(mat, full_matrices=True)
        tol = max(mat.shape) * np.finfo(float).eps * (float(singular[0]) if singular.size else 0.0)
        rank = int(np.count_nonzero(singular > tol)) if singular.size else 0
        ranks.append(rank)
        left_null = u[:, rank:]
        target_vec = sym_to_vec(cases[0].r_need[i, j])
        target_norm = max(float(np.linalg.norm(target_vec)), 1.0e-300)
        rho_w = float(np.sqrt(max(cases[0].rho_a[i, j], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300)))
        weight = rho_w / target_norm
        for null_vec in left_null.T:
            row = np.zeros(int(meta["ncols"]), dtype=float)
            for apos, atom in enumerate(ATOM_NAMES):
                atom_vec = sym_to_vec(atoms[atom][i, j])
                add_coeff(row, col_index, 1, apos, i, j, weight * float(np.dot(null_vec, atom_vec)))
            rows.append(row)
            rhs.append(weight * float(np.dot(null_vec, target_vec)))
    info = {
        "constraint_rows": int(len(rows)),
        "rank_counts": {str(k): int(ranks.count(k)) for k in sorted(set(ranks))},
    }
    if not rows:
        return np.zeros((0, int(meta["ncols"])), dtype=float), np.zeros(0, dtype=float), info
    return np.vstack(rows), np.asarray(rhs, dtype=float), info


def solve_nonhomogeneous_hard(a: np.ndarray, b: np.ndarray, f: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    col_scale = np.linalg.norm(a, axis=0) + np.linalg.norm(f, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_s = a / col_scale[None, :]
    f_s = f / col_scale[None, :]
    u, singular, vh = np.linalg.svd(f_s, full_matrices=True)
    if singular.size:
        tol = max(f_s.shape) * np.finfo(float).eps * float(singular[0])
        rank = int(np.count_nonzero(singular > tol))
    else:
        tol = 0.0
        rank = 0
    if rank:
        y0 = vh[:rank].T @ ((u[:, :rank].T @ d) / singular[:rank])
    else:
        y0 = np.zeros(a.shape[1], dtype=float)
    null = vh[rank:].T
    if null.shape[1] > 0:
        lhs = a_s @ null
        rhs = b - a_s @ y0
        z = np.linalg.lstsq(lhs, rhs, rcond=1.0e-12)[0]
        y = y0 + null @ z
    else:
        y = y0
    x = y / col_scale
    eq = f @ x - d
    obj = a @ x - b
    return x, {
        "hard_rank": rank,
        "hard_nullity": int(null.shape[1]),
        "hard_singular_min": float(singular[-1]) if singular.size else 0.0,
        "hard_singular_max": float(singular[0]) if singular.size else 0.0,
        "hard_singular_tol": float(tol),
        "hard_equality_norm": float(np.linalg.norm(eq)),
        "hard_equality_max_abs": float(np.max(np.abs(eq))) if eq.size else 0.0,
        "objective_norm": float(np.linalg.norm(obj)),
    }


def principal_project_residual(cases: dict[int, object], residual: np.ndarray, mask: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    predicted = np.zeros_like(residual)
    delta_metric = np.zeros_like(residual)
    for i_raw, j_raw in np.argwhere(mask):
        i = int(i_raw)
        j = int(j_raw)
        mat = local_principal_matrix(cases[0].metric_cov[i, j], cases[0].metric_inv[i, j], dt)
        rhs = sym_to_vec(residual[i, j])
        sol, *_ = np.linalg.lstsq(mat, rhs, rcond=1.0e-12)
        predicted[i, j] = vec_to_sym(mat @ sol)
        delta_metric[i, j] = vec_to_sym(sol)
    return predicted, delta_metric


def determinant_stats(metric: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    det = np.linalg.det(metric)
    vals = det[mask & np.isfinite(det)]
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95_abs": 0.0, "max": 0.0, "sign_changes": 0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95_abs": float(np.percentile(np.abs(vals), 95.0)),
        "max": float(np.max(vals)),
        "negative_count": int(np.count_nonzero(vals < 0.0)),
        "positive_count": int(np.count_nonzero(vals > 0.0)),
    }


def render(path: Path, x: np.ndarray, z: np.ndarray, mask: np.ndarray, before: np.ndarray, after: np.ndarray, hnorm: np.ndarray) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), constrained_layout=True)
    fields = [
        ("before projection residual", before, "magma"),
        ("after principal projection residual", after, "magma"),
        (r"$|\delta \tilde g_+|$ needed", hnorm, "viridis"),
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
    dt = float(args.force_probe_dt_old) * old_scale
    cases = {
        key: build_case(args, ref, float(args.tau) + key * float(args.force_probe_dt_old))
        for key in base.TIME_KEYS
    }
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    a, b, meta = base.assemble_system(
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
    a_alg, b_alg, f_force, d_force = split_alg_force(a, b, meta)

    prelim = np.linalg.lstsq(a_alg, b_alg, rcond=1.0e-12)[0]
    scales = atom_scales_from_coeff_guess(prelim, meta)
    a_reg, b_reg, reg_counts = add_regularizer_rows(
        meta,
        scales,
        norm_weight=float(args.norm_weight),
        time_weight=float(args.time_weight),
    )
    f_constraint, d_constraint, constraint_info = principal_constraint_rows(cases, meta, args.fit_region, dt)
    f_hard = np.vstack([f_force, f_constraint])
    d_hard = np.concatenate([d_force, d_constraint])
    a_obj = np.vstack([a_alg, a_reg])
    b_obj = np.concatenate([b_alg, b_reg])
    coeff, solve = solve_nonhomogeneous_hard(a_obj, b_obj, f_hard, d_hard)

    fields = base.coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
    cov = base.coeff_fields_to_cov(cases, fields)
    mask = getattr(cases[0], args.fit_region)
    residual = cov[0] - cases[0].r_need
    projected, delta_metric = principal_project_residual(cases, residual, mask, dt)
    corrected_metric_plus = cases[1].metric_cov + delta_metric
    after = residual - projected
    before_rel = tensor_norm(residual) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    after_rel = tensor_norm(after) / np.maximum(tensor_norm(cases[0].r_need + projected), 1.0e-300)
    hnorm = tensor_norm(delta_metric)
    metric_norm = tensor_norm(cases[0].metric_cov)
    weights = np.sqrt(np.maximum(cases[0].rho_a[mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    plot_path = args.output / "gbcd_principal_constraint_projection.png"
    render(plot_path, cases[0].x, cases[0].z, mask, before_rel, after_rel, hnorm)
    coeff_path = args.output / "principal_constraint_coefficients.npz"
    save = {
        "active_mask": np.asarray(meta["active_mask"], dtype=bool),
        "mask": mask,
        "before_rel": before_rel,
        "after_rel": after_rel,
        "projected_delta_einstein": projected,
        "delta_metric_plus": delta_metric,
        "corrected_metric_plus": corrected_metric_plus,
    }
    for key in base.TIME_KEYS:
        suffix = {-1: "m", 0: "0", 1: "p"}[key]
        for apos, atom in enumerate(ATOM_NAMES):
            save[f"{atom}_{suffix}"] = fields[key][..., apos]
    np.savez_compressed(coeff_path, **save)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "norm_weight": float(args.norm_weight),
            "time_weight": float(args.time_weight),
            "dt_ev_inv": float(dt),
        },
        "definition": {
            "goal": "Choose lambda so that the center algebraic residual lies in the image of the metric time-principal Einstein operator, while full conservation remains hard.",
            "hard_constraints": "full conservation plus left-null/principal constraints N^T(C_0-r_need_0)=0.",
        },
        "counts": {
            "ncols": int(meta["ncols"]),
            "algebraic_rows": int(a_alg.shape[0]),
            "force_rows": int(f_force.shape[0]),
            "principal_constraint_rows": int(f_constraint.shape[0]),
            "regularizer_counts": reg_counts,
            **constraint_info,
        },
        "solve": solve,
        "evaluation": {
            "before_principal_projection_relative_residual": weighted_stats(before_rel[mask], weights),
            "after_principal_projection_relative_residual": weighted_stats(after_rel[mask], weights),
            "delta_metric_plus_norm": weighted_stats(hnorm[mask], weights),
            "relative_delta_metric_plus_norm": weighted_stats(
                hnorm[mask] / np.maximum(metric_norm[mask], 1.0e-300),
                weights,
            ),
            "metric_plus_det_base": determinant_stats(cases[1].metric_cov, mask),
            "metric_plus_det_corrected": determinant_stats(corrected_metric_plus, mask),
        },
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
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--norm-weight", type=float, default=1.0e-8)
    parser.add_argument("--time-weight", type=float, default=1.0e-5)
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
    parser.add_argument("--active-dilation", type=int, default=1)
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
