from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_gbcd_metric_principal_correction import local_principal_matrix, sym_to_vec
from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from solve_gbcd_full_linear_metric_update_sparse import lsqr_solve
from test_gbcd_v0_one_step_lambda import ATOM_NAMES as DEFAULT_ATOM_NAMES
from test_gbcd_v0_one_step_lambda import weighted_stats


ATOM_NAMES = tuple(DEFAULT_ATOM_NAMES)
ATOM_FAMILY = "matter4"


def add_entry(entries: list[tuple[int, float]], col_index: np.ndarray, tpos: int, apos: int, i: int, j: int, value: float) -> None:
    c = int(col_index[tpos, apos, i, j])
    if c >= 0 and np.isfinite(value) and value != 0.0:
        entries.append((c, float(value)))


def append_row(row_cols: list[np.ndarray], row_vals: list[np.ndarray], rhs: list[float], entries: list[tuple[int, float]], b: float) -> None:
    if not entries:
        return
    acc: dict[int, float] = {}
    for c, v in entries:
        acc[c] = acc.get(c, 0.0) + float(v)
    cols = np.fromiter(acc.keys(), dtype=np.int32)
    vals = np.fromiter(acc.values(), dtype=float)
    good = np.isfinite(vals) & (vals != 0.0)
    if not np.any(good):
        return
    row_cols.append(cols[good])
    row_vals.append(vals[good])
    rhs.append(float(b))


def sparse_apply(row_cols: list[np.ndarray], row_vals: list[np.ndarray], x: np.ndarray) -> np.ndarray:
    y = np.zeros(len(row_cols), dtype=float)
    for r, (cols, vals) in enumerate(zip(row_cols, row_vals)):
        y[r] = float(np.dot(vals, x[cols]))
    return y


def sparse_apply_t(row_cols: list[np.ndarray], row_vals: list[np.ndarray], y: np.ndarray, ncols: int) -> np.ndarray:
    out = np.zeros(ncols, dtype=float)
    for yr, cols, vals in zip(y, row_cols, row_vals):
        if yr != 0.0:
            np.add.at(out, cols, yr * vals)
    return out


def column_norms_from_rows(row_cols: list[np.ndarray], row_vals: list[np.ndarray], ncols: int) -> np.ndarray:
    norms2 = np.zeros(ncols, dtype=float)
    for cols, vals in zip(row_cols, row_vals):
        np.add.at(norms2, cols, vals * vals)
    norms = np.sqrt(norms2)
    finite = norms[np.isfinite(norms) & (norms > 0.0)]
    floor = 1.0e-30 if finite.size == 0 else max(float(np.percentile(finite, 5.0)) * 1.0e-12, 1.0e-30)
    return np.where(norms > floor, norms, floor)


def build_entry_arrays_from_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    col_norms: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    val_parts: list[np.ndarray] = []
    for row, (cols, vals) in enumerate(zip(row_cols, row_vals)):
        if cols.size == 0:
            continue
        scaled = vals / col_norms[cols]
        good = np.isfinite(scaled) & (scaled != 0.0)
        if np.any(good):
            row_parts.append(np.full(int(np.count_nonzero(good)), row, dtype=np.int64))
            col_parts.append(cols[good].astype(np.int64, copy=False))
            val_parts.append(scaled[good].astype(float, copy=False))
    if not row_parts:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0, dtype=float)
    return np.concatenate(row_parts), np.concatenate(col_parts), np.concatenate(val_parts)


def entry_apply(row_idx: np.ndarray, col_idx: np.ndarray, vals: np.ndarray, x: np.ndarray, nrows: int) -> np.ndarray:
    if vals.size == 0:
        return np.zeros(nrows, dtype=float)
    return np.bincount(row_idx, weights=vals * x[col_idx], minlength=nrows).astype(float, copy=False)


def entry_apply_t(row_idx: np.ndarray, col_idx: np.ndarray, vals: np.ndarray, y: np.ndarray, ncols: int) -> np.ndarray:
    if vals.size == 0:
        return np.zeros(ncols, dtype=float)
    return np.bincount(col_idx, weights=vals * y[row_idx], minlength=ncols).astype(float, copy=False)


def add_algebraic_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    col_index: np.ndarray,
    region: str,
    target_floor_frac: float,
) -> int:
    count = 0
    atom_cache = {key: tensor_atoms(cases[key], ATOM_FAMILY) for key in base.TIME_KEYS}
    for tpos, key in enumerate(base.TIME_KEYS):
        case = cases[key]
        mask = getattr(case, region)
        idx = np.where(mask.reshape(-1))[0]
        target_rows = symmetric_rows(case.r_need, idx)
        weights = base.relative_row_weight(case, idx, target_rows, target_floor_frac)
        for p, (i_raw, j_raw) in enumerate(np.argwhere(mask)):
            i = int(i_raw)
            j = int(j_raw)
            for comp_pos, (mu, nu) in enumerate(base.SYMMETRIC_COMPONENTS):
                entries = []
                for apos, atom in enumerate(ATOM_NAMES):
                    add_entry(entries, col_index, tpos, apos, i, j, weights[p] * atom_cache[key][atom][i, j, mu, nu])
                append_row(row_cols, row_vals, rhs, entries, weights[p] * target_rows[p, comp_pos])
                count += 1
    return count


def add_force_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    full_case,
    col_index: np.ndarray,
    force_region: str,
    target_floor_frac: float,
    force_weight: float,
    force_erosion: int,
    dt: float,
    dx: float,
    dz: float,
) -> int:
    force_mask = base.make_force_mask(cases[0], force_region, force_erosion)
    force_weights = base.force_row_weight(cases[0], force_mask, min(abs(dt), abs(dx), abs(dz)), target_floor_frac, force_weight)
    metric_inv = full_case.geom_0.metric_inv
    gamma = full_case.geom_0.gamma
    atom_cache = {key: tensor_atoms(cases[key], ATOM_FAMILY) for key in base.TIME_KEYS}
    count = 0
    for p, (i_raw, j_raw) in enumerate(np.argwhere(force_mask)):
        i = int(i_raw)
        j = int(j_raw)
        row_sigma: list[list[tuple[int, float]]] = [[], [], []]
        for sigma in range(3):
            for mu in range(3):
                for alpha in range(3):
                    pref = metric_inv[i, j, mu, alpha]
                    if not np.isfinite(pref):
                        continue
                    if alpha == 0:
                        for apos, atom in enumerate(ATOM_NAMES):
                            add_entry(row_sigma[sigma], col_index, 2, apos, i, j, pref * atom_cache[1][atom][i, j, mu, sigma] / (2.0 * dt))
                            add_entry(row_sigma[sigma], col_index, 0, apos, i, j, -pref * atom_cache[-1][atom][i, j, mu, sigma] / (2.0 * dt))
                    elif alpha == 1:
                        for apos, atom in enumerate(ATOM_NAMES):
                            add_entry(row_sigma[sigma], col_index, 1, apos, i + 1, j, pref * atom_cache[0][atom][i + 1, j, mu, sigma] / (2.0 * dx))
                            add_entry(row_sigma[sigma], col_index, 1, apos, i - 1, j, -pref * atom_cache[0][atom][i - 1, j, mu, sigma] / (2.0 * dx))
                    else:
                        for apos, atom in enumerate(ATOM_NAMES):
                            add_entry(row_sigma[sigma], col_index, 1, apos, i, j + 1, pref * atom_cache[0][atom][i, j + 1, mu, sigma] / (2.0 * dz))
                            add_entry(row_sigma[sigma], col_index, 1, apos, i, j - 1, -pref * atom_cache[0][atom][i, j - 1, mu, sigma] / (2.0 * dz))
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, mu]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom in enumerate(ATOM_NAMES):
                                add_entry(row_sigma[sigma], col_index, 1, apos, i, j, conn * atom_cache[0][atom][i, j, lam, sigma])
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, sigma]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom in enumerate(ATOM_NAMES):
                                add_entry(row_sigma[sigma], col_index, 1, apos, i, j, conn * atom_cache[0][atom][i, j, mu, lam])
        for sigma in range(3):
            entries = [(c, force_weights[p] * v) for c, v in row_sigma[sigma]]
            append_row(row_cols, row_vals, rhs, entries, 0.0)
            count += 1
    return count


def add_principal_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    col_index: np.ndarray,
    region: str,
    dt: float,
    weight: float,
) -> tuple[int, dict[str, int]]:
    mask = getattr(cases[0], region)
    atoms = tensor_atoms(cases[0], ATOM_FAMILY)
    ranks = []
    count = 0
    for i_raw, j_raw in np.argwhere(mask):
        i = int(i_raw)
        j = int(j_raw)
        mat = local_principal_matrix(cases[0].metric_cov[i, j], cases[0].metric_inv[i, j], dt)
        u, singular, _ = np.linalg.svd(mat, full_matrices=True)
        tol = max(mat.shape) * np.finfo(float).eps * (float(singular[0]) if singular.size else 0.0)
        rank = int(np.count_nonzero(singular > tol)) if singular.size else 0
        ranks.append(rank)
        target_vec = sym_to_vec(cases[0].r_need[i, j])
        target_norm = max(float(np.linalg.norm(target_vec)), 1.0e-300)
        rho_w = float(np.sqrt(max(cases[0].rho_a[i, j], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300)))
        row_weight = float(weight) * rho_w / target_norm
        for null_vec in u[:, rank:].T:
            entries = []
            for apos, atom in enumerate(ATOM_NAMES):
                add_entry(entries, col_index, 1, apos, i, j, row_weight * float(np.dot(null_vec, sym_to_vec(atoms[atom][i, j]))))
            append_row(row_cols, row_vals, rhs, entries, row_weight * float(np.dot(null_vec, target_vec)))
            count += 1
    return count, {str(k): int(ranks.count(k)) for k in sorted(set(ranks))}


def add_regularizer_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    meta: dict[str, object],
    atom_scales: dict[str, float],
    norm_weight: float,
    time_weight: float,
) -> dict[str, int]:
    col = np.asarray(meta["column_index"], dtype=int)
    active = np.asarray(meta["active_mask"], dtype=bool)
    counts = {"norm_rows": 0, "time_rows": 0}
    sqrt_norm = float(np.sqrt(max(norm_weight, 0.0)))
    sqrt_time = float(np.sqrt(max(time_weight, 0.0)))
    for apos, atom in enumerate(ATOM_NAMES):
        scale = max(float(atom_scales[atom]), 1.0e-30)
        for i, j in np.argwhere(active):
            i = int(i)
            j = int(j)
            if sqrt_norm > 0.0:
                for tpos in range(3):
                    c = int(col[tpos, apos, i, j])
                    if c >= 0:
                        append_row(row_cols, row_vals, rhs, [(c, sqrt_norm / scale)], 0.0)
                        counts["norm_rows"] += 1
            if sqrt_time > 0.0:
                for t0, t1 in [(0, 1), (1, 2)]:
                    c0 = int(col[t0, apos, i, j])
                    c1 = int(col[t1, apos, i, j])
                    if c0 >= 0 and c1 >= 0:
                        append_row(row_cols, row_vals, rhs, [(c1, sqrt_time / scale), (c0, -sqrt_time / scale)], 0.0)
                        counts["time_rows"] += 1
    return counts


def estimate_atom_scales(cases: dict[int, object], active_mask: np.ndarray) -> dict[str, float]:
    vals = {atom: [] for atom in ATOM_NAMES}
    for key in base.TIME_KEYS:
        atoms = tensor_atoms(cases[key], ATOM_FAMILY)
        for i, j in np.argwhere(active_mask):
            a = np.stack([sym_to_vec(atoms[atom][int(i), int(j)]) for atom in ATOM_NAMES], axis=1)
            b = sym_to_vec(cases[key].r_need[int(i), int(j)])
            try:
                sol = np.linalg.lstsq(a, b, rcond=1.0e-12)[0]
            except np.linalg.LinAlgError:
                sol = np.zeros(len(ATOM_NAMES), dtype=float)
            for apos, atom in enumerate(ATOM_NAMES):
                if np.isfinite(sol[apos]):
                    vals[atom].append(float(sol[apos]))
    scales = {}
    for atom in ATOM_NAMES:
        arr = np.asarray(vals[atom], dtype=float)
        finite = arr[np.isfinite(arr)]
        scales[atom] = max(float(np.percentile(np.abs(finite), 95.0)) if finite.size else 1.0, 1.0e-30)
    return scales


def solve_sparse_lsqr(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: np.ndarray,
    ncols: int,
    tol: float,
    maxiter: int,
) -> tuple[np.ndarray, dict[str, object]]:
    col_norms = column_norms_from_rows(row_cols, row_vals, ncols)
    row_idx, col_idx, vals = build_entry_arrays_from_rows(row_cols, row_vals, col_norms)
    nrows = len(row_cols)

    def aprod(y: np.ndarray) -> np.ndarray:
        return entry_apply(row_idx, col_idx, vals, y, nrows)

    def atprod(u: np.ndarray) -> np.ndarray:
        return entry_apply_t(row_idx, col_idx, vals, u, ncols)

    y, info = lsqr_solve(aprod, atprod, rhs, ncols, tol=tol, maxiter=maxiter)
    x = y / col_norms
    residual = sparse_apply(row_cols, row_vals, x) - rhs
    info["system_relative_residual"] = float(np.linalg.norm(residual) / max(float(np.linalg.norm(rhs)), 1.0e-300))
    info["operator_backend"] = "entry-arrays"
    info["scalar_nonzeros"] = int(vals.size)
    return x, info


def run(args: argparse.Namespace) -> dict[str, object]:
    global ATOM_NAMES, ATOM_FAMILY
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    cases = {key: build_case(args, ref, float(args.tau) + key * float(args.force_probe_dt_old)) for key in base.TIME_KEYS}
    full_case = build_full_case(args, ref, float(args.tau))
    ATOM_FAMILY = str(args.atom_family)
    ATOM_NAMES = tuple(tensor_atoms(cases[0], ATOM_FAMILY).keys())
    base.ATOM_NAMES = ATOM_NAMES
    base.ATOM_FAMILY = ATOM_FAMILY
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    active_mask = base.make_active_mask(cases, args.fit_region, args.force_region, int(args.active_dilation))
    col_index, ncols = base.build_column_index(active_mask)
    meta = {
        "active_mask": active_mask,
        "column_index": col_index,
        "ncols": int(ncols),
        "force_count": int(np.count_nonzero(base.make_force_mask(cases[0], args.force_region, int(args.force_mask_erosion)))),
    }
    row_cols: list[np.ndarray] = []
    row_vals: list[np.ndarray] = []
    rhs_list: list[float] = []
    algebraic_rows = add_algebraic_rows(row_cols, row_vals, rhs_list, cases, col_index, args.fit_region, float(args.target_floor_frac))
    force_rows = add_force_rows(
        row_cols,
        row_vals,
        rhs_list,
        cases,
        full_case,
        col_index,
        args.force_region,
        float(args.target_floor_frac),
        float(args.force_penalty),
        int(args.force_mask_erosion),
        dt,
        dx,
        dz,
    )
    principal_rows, rank_counts = add_principal_rows(
        row_cols,
        row_vals,
        rhs_list,
        cases,
        col_index,
        args.fit_region,
        dt,
        float(args.principal_penalty),
    )
    atom_scales = estimate_atom_scales(cases, active_mask)
    reg_counts = add_regularizer_rows(row_cols, row_vals, rhs_list, meta, atom_scales, float(args.norm_weight), float(args.time_weight))
    rhs = np.asarray(rhs_list, dtype=float)
    coeff, solve = solve_sparse_lsqr(row_cols, row_vals, rhs, ncols, float(args.lsqr_tol), int(args.lsqr_maxiter))
    fields = base.coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
    cov = base.coeff_fields_to_cov(cases, fields)
    mask = getattr(cases[0], args.fit_region)
    weights = np.sqrt(np.maximum(cases[0].rho_a[mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    rel0 = tensor_norm(cov[0] - cases[0].r_need) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    out_path = args.output / "principal_constraint_coefficients.npz"
    save = {
        "active_mask": active_mask,
        "mask": mask,
        "central_relative_residual": rel0,
        "atom_family": np.asarray(ATOM_FAMILY),
        "atom_names": np.asarray(ATOM_NAMES),
    }
    for key in base.TIME_KEYS:
        suffix = {-1: "m", 0: "0", 1: "p"}[key]
        for apos, atom in enumerate(ATOM_NAMES):
            save[f"{atom}_{suffix}"] = fields[key][..., apos]
    np.savez_compressed(out_path, **save)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "force_penalty": float(args.force_penalty),
            "principal_penalty": float(args.principal_penalty),
            "norm_weight": float(args.norm_weight),
            "time_weight": float(args.time_weight),
            "atom_family": ATOM_FAMILY,
            "atom_names": list(ATOM_NAMES),
            "dt_ev_inv": float(dt),
        },
        "definition": {
            "goal": "Sparse LSQR prototype for principal-constraint lambda projection.",
            "caveat": "Force and principal constraints are enforced by large penalty rows in this prototype, not exact nullspace hard constraints.",
        },
        "counts": {
            "ncols": int(ncols),
            "rows": int(len(row_cols)),
            "algebraic_rows": int(algebraic_rows),
            "force_rows": int(force_rows),
            "principal_rows": int(principal_rows),
            "regularizer_counts": reg_counts,
            "atom_scales": atom_scales,
            "rank_counts": rank_counts,
            "nnz": int(sum(len(c) for c in row_cols)),
        },
        "solve": solve,
        "evaluation": {
            "central_algebraic_relative_residual": weighted_stats(rel0[mask], weights),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "coefficients_npz": str(out_path.resolve()),
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
    parser.add_argument("--force-penalty", type=float, default=100.0)
    parser.add_argument("--principal-penalty", type=float, default=100.0)
    parser.add_argument("--norm-weight", type=float, default=1.0e-8)
    parser.add_argument("--time-weight", type=float, default=1.0e-5)
    parser.add_argument(
        "--atom-family",
        choices=["matter4", "matter4_plus_normal", "matter6_normal"],
        default="matter4",
    )
    parser.add_argument("--lsqr-tol", type=float, default=1.0e-6)
    parser.add_argument("--lsqr-maxiter", type=int, default=1000)
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
