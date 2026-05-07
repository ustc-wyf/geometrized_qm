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
from diagnose_gbcd_trace_local_closure import ATOM_NAMES, solve_local, trace_coefficients, weighted_stats
from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_gbcd_principal_constraint_projection_sparse import append_row, solve_sparse_lsqr, sparse_apply
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


TIME_KEYS = (-1, 0, 1)
SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def trace_nullspace(row: np.ndarray) -> np.ndarray | None:
    row = np.asarray(row, dtype=float)
    if not np.all(np.isfinite(row)):
        return None
    norm = float(np.linalg.norm(row))
    if norm <= 1.0e-300:
        return None
    _, s, vh = np.linalg.svd(row.reshape(1, -1), full_matrices=True)
    tol = row.size * np.finfo(float).eps * (float(s[0]) if s.size else 0.0)
    rank = int(np.count_nonzero(s > tol))
    n = vh[rank:].T
    return n if n.shape == (len(ATOM_NAMES), len(ATOM_NAMES) - 1) else None


def make_active_mask(cases: dict[int, object], region: str, force_region: str, active_dilation: int) -> np.ndarray:
    base.ATOM_NAMES = ATOM_NAMES
    return base.make_active_mask(cases, region, force_region, active_dilation)


def build_column_index(active_mask: np.ndarray) -> tuple[np.ndarray, int]:
    free_dim = len(ATOM_NAMES) - 1
    col = -np.ones((len(TIME_KEYS), free_dim) + active_mask.shape, dtype=int)
    ncols = 0
    for tpos, _ in enumerate(TIME_KEYS):
        for k in range(free_dim):
            for i_raw, j_raw in np.argwhere(active_mask):
                col[tpos, k, int(i_raw), int(j_raw)] = ncols
                ncols += 1
    return col, ncols


def build_basis(cases: dict[int, object], active_mask: np.ndarray) -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    for key in TIME_KEYS:
        rows = trace_coefficients(cases[key])
        basis = np.zeros(cases[key].rho_a.shape + (len(ATOM_NAMES), len(ATOM_NAMES) - 1), dtype=float)
        for i_raw, j_raw in np.argwhere(active_mask):
            i = int(i_raw)
            j = int(j_raw)
            n = trace_nullspace(rows[i, j])
            if n is not None:
                basis[i, j, :, :] = n
        out[key] = basis
    return out


def add_projected(
    entries: list[tuple[int, float]],
    col_index: np.ndarray,
    basis: dict[int, np.ndarray],
    tpos: int,
    key: int,
    apos: int,
    i: int,
    j: int,
    value: float,
) -> None:
    if not np.isfinite(value) or value == 0.0:
        return
    for k in range(len(ATOM_NAMES) - 1):
        c = int(col_index[tpos, k, i, j])
        if c < 0:
            continue
        v = float(value) * float(basis[key][i, j, apos, k])
        if np.isfinite(v) and v != 0.0:
            entries.append((c, v))


def add_algebraic_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
    col_index: np.ndarray,
    basis: dict[int, np.ndarray],
    region: str,
    target_floor_frac: float,
) -> int:
    count = 0
    atom_cache = {key: tensor_atoms(cases[key], "matter4") for key in TIME_KEYS}
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


def add_force_rows(
    row_cols: list[np.ndarray],
    row_vals: list[np.ndarray],
    rhs: list[float],
    cases: dict[int, object],
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
    atom_cache = {key: tensor_atoms(cases[key], "matter4") for key in TIME_KEYS}
    force_mask = base.make_force_mask(cases[0], force_region, force_erosion)
    force_weights = base.force_row_weight(
        cases[0],
        force_mask,
        min(abs(dt), abs(dx), abs(dz)),
        target_floor_frac,
        1.0,
    )
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
                            val_p = pref * atom_cache[1][atom_name][i, j, mu, sigma] / (2.0 * dt)
                            val_m = -pref * atom_cache[-1][atom_name][i, j, mu, sigma] / (2.0 * dt)
                            add_projected(entries, col_index, basis, 2, 1, apos, i, j, val_p)
                            add_projected(entries, col_index, basis, 0, -1, apos, i, j, val_m)
                    elif alpha == 1:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            val_p = pref * atom_cache[0][atom_name][i + 1, j, mu, sigma] / (2.0 * dx)
                            val_m = -pref * atom_cache[0][atom_name][i - 1, j, mu, sigma] / (2.0 * dx)
                            add_projected(entries, col_index, basis, 1, 0, apos, i + 1, j, val_p)
                            add_projected(entries, col_index, basis, 1, 0, apos, i - 1, j, val_m)
                    else:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            val_p = pref * atom_cache[0][atom_name][i, j + 1, mu, sigma] / (2.0 * dz)
                            val_m = -pref * atom_cache[0][atom_name][i, j - 1, mu, sigma] / (2.0 * dz)
                            add_projected(entries, col_index, basis, 1, 0, apos, i, j + 1, val_p)
                            add_projected(entries, col_index, basis, 1, 0, apos, i, j - 1, val_m)
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, mu]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom_name in enumerate(ATOM_NAMES):
                                value = conn * atom_cache[0][atom_name][i, j, lam, sigma]
                                add_projected(entries, col_index, basis, 1, 0, apos, i, j, value)
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, sigma]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom_name in enumerate(ATOM_NAMES):
                                value = conn * atom_cache[0][atom_name][i, j, mu, lam]
                                add_projected(entries, col_index, basis, 1, 0, apos, i, j, value)
            entries = [(c, float(force_weights[p]) * v) for c, v in entries]
            append_row(row_cols, row_vals, rhs, entries, 0.0)
            count += 1
    return count


def scaled_rows(
    alg_cols: list[np.ndarray],
    alg_vals: list[np.ndarray],
    alg_rhs: list[float],
    force_cols: list[np.ndarray],
    force_vals: list[np.ndarray],
    force_rhs: list[float],
    force_weight: float,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray]:
    cols = list(alg_cols)
    vals = list(alg_vals)
    rhs = list(alg_rhs)
    w = float(force_weight)
    for c, v, b in zip(force_cols, force_vals, force_rhs):
        cols.append(c)
        vals.append(w * v)
        rhs.append(w * float(b))
    return cols, vals, np.asarray(rhs, dtype=float)


def coeff_to_lambda_fields(coeff: np.ndarray, col_index: np.ndarray, basis: dict[int, np.ndarray], shape: tuple[int, int]) -> dict[int, np.ndarray]:
    fields: dict[int, np.ndarray] = {}
    for tpos, key in enumerate(TIME_KEYS):
        arr = np.zeros(shape + (len(ATOM_NAMES),), dtype=float)
        for i_raw, j_raw in np.argwhere(col_index[tpos, 0] >= 0):
            i = int(i_raw)
            j = int(j_raw)
            y = np.asarray([coeff[int(col_index[tpos, k, i, j])] for k in range(len(ATOM_NAMES) - 1)], dtype=float)
            arr[i, j, :] = basis[key][i, j] @ y
        fields[key] = arr
    return fields


def local_initial_coeff(cases: dict[int, object], col_index: np.ndarray, basis: dict[int, np.ndarray]) -> np.ndarray:
    ncols = int(np.max(col_index)) + 1 if np.any(col_index >= 0) else 0
    coeff = np.zeros(ncols, dtype=float)
    atom_cache = {key: tensor_atoms(cases[key], "matter4") for key in TIME_KEYS}
    trace_cache = {key: trace_coefficients(cases[key]) for key in TIME_KEYS}
    for tpos, key in enumerate(TIME_KEYS):
        for i_raw, j_raw in np.argwhere(col_index[tpos, 0] >= 0):
            i = int(i_raw)
            j = int(j_raw)
            a = np.zeros((len(SYMMETRIC_COMPONENTS), len(ATOM_NAMES)), dtype=float)
            for apos, atom_name in enumerate(ATOM_NAMES):
                for cpos, (mu, nu) in enumerate(SYMMETRIC_COMPONENTS):
                    a[cpos, apos] = atom_cache[key][atom_name][i, j, mu, nu]
            b = np.asarray([cases[key].r_need[i, j, mu, nu] for mu, nu in SYMMETRIC_COMPONENTS], dtype=float)
            lam, _ = solve_local(a, b, trace_cache[key][i, j])
            if not np.all(np.isfinite(lam)):
                continue
            y = basis[key][i, j].T @ lam
            for k in range(len(ATOM_NAMES) - 1):
                c = int(col_index[tpos, k, i, j])
                if c >= 0 and np.isfinite(y[k]):
                    coeff[c] = float(y[k])
    return coeff


def tensor_fields(cases: dict[int, object], lambda_fields: dict[int, np.ndarray]) -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    for key in TIME_KEYS:
        atoms = tensor_atoms(cases[key], "matter4")
        cov = np.zeros_like(cases[key].r_need)
        for apos, atom_name in enumerate(ATOM_NAMES):
            cov += lambda_fields[key][..., apos, None, None] * atoms[atom_name]
        out[key] = cov
    return out


def evaluate_solution(
    cases: dict[int, object],
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
    cov = tensor_fields(cases, lambda_fields)
    by_time = {}
    for key in TIME_KEYS:
        mask = getattr(cases[key], region)
        active = col_index[TIME_KEYS.index(key), 0] >= 0
        err_norm = tensor_norm(cov[key] - cases[key].r_need)
        target_norm = tensor_norm(cases[key].r_need)
        idx = np.where(mask.reshape(-1))[0]
        err_rows = symmetric_rows(cov[key] - cases[key].r_need, idx)
        target_rows = symmetric_rows(cases[key].r_need, idx)
        err_row_norm = np.sqrt(np.sum(err_rows**2, axis=1))
        target_row_norm = np.sqrt(np.sum(target_rows**2, axis=1))
        active_idx = np.where(active.reshape(-1))[0]
        active_target_rows = symmetric_rows(cases[key].r_need, active_idx)
        active_target_row_norm = np.sqrt(np.sum(active_target_rows**2, axis=1))
        finite_target = target_norm[active & np.isfinite(target_norm)]
        finite_target_rows = active_target_row_norm[np.isfinite(active_target_row_norm)]
        floor = 0.05 * max(float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0, 1.0e-300)
        row_floor = 0.05 * max(
            float(np.percentile(finite_target_rows, 95.0)) if finite_target_rows.size else 0.0,
            1.0e-300,
        )
        rel = err_norm / np.maximum(target_norm, 1.0e-300)
        rel_floor = err_norm / np.maximum(target_norm, floor)
        rel_row_floor = err_row_norm / np.maximum(target_row_norm, row_floor)
        weights = np.sqrt(np.maximum(cases[key].rho_a[mask], 0.0) / max(float(np.max(cases[key].rho_a)), 1.0e-300))
        trace_value = np.einsum("...ab,...ab->...", cases[key].metric_inv, cov[key], optimize=True)
        by_time[str(key)] = {
            "tau": float(cases[key].tau_old),
            "algebraic_relative_residual": weighted_stats(rel[mask], weights),
            "algebraic_floor_relative_residual": weighted_stats(rel_floor[mask], weights),
            "algebraic_row_floor_relative_residual": weighted_stats(rel_row_floor, weights),
            "trace_abs": weighted_stats(trace_value[mask], weights),
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
    weights = [float(r["force_weight"]) for r in records]
    alg = [float(r["evaluation"]["by_time"]["0"]["algebraic_row_floor_relative_residual"]["weighted_mean"]) for r in records]
    div = [float(r["evaluation"]["central_divergence"]["full_over_derivative_scale"]["weighted_mean"]) for r in records]
    x = np.arange(len(weights))
    fig, ax = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
    ax.plot(x, alg, "o-", label="algebraic residual")
    ax.plot(x, div, "o-", label="full divergence scale")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{w:g}" for w in weights], rotation=25, ha="right")
    ax.set_xlabel("force weight")
    ax.set_title("trace0 hard-eliminated sparse conservation scan")
    ax.legend()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ATOM_NAMES
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    probe = float(args.force_probe_dt_old)
    tau = float(args.tau)
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
    active = make_active_mask(cases, args.fit_region, args.force_region, int(args.active_dilation))
    col_index, ncols = build_column_index(active)
    basis = build_basis(cases, active)

    alg_cols: list[np.ndarray] = []
    alg_vals: list[np.ndarray] = []
    alg_rhs: list[float] = []
    force_cols: list[np.ndarray] = []
    force_vals: list[np.ndarray] = []
    force_rhs: list[float] = []
    n_alg = add_algebraic_rows(alg_cols, alg_vals, alg_rhs, cases, col_index, basis, args.fit_region, float(args.target_floor_frac))
    n_force = add_force_rows(
        force_cols,
        force_vals,
        force_rhs,
        cases,
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
    initial = local_initial_coeff(cases, col_index, basis)
    records = []
    coeff_arrays = {}
    for force_weight in [float(x.strip()) for x in args.force_weights.split(",") if x.strip()]:
        print(f"[trace0-sparse] tau={tau:g} force_weight={force_weight:g}", flush=True)
        row_cols, row_vals, rhs = scaled_rows(alg_cols, alg_vals, alg_rhs, force_cols, force_vals, force_rhs, force_weight)
        rhs_delta = rhs - sparse_apply(row_cols, row_vals, initial)
        delta, solve = solve_sparse_lsqr(row_cols, row_vals, rhs_delta, ncols, float(args.lsqr_tol), int(args.lsqr_maxiter))
        coeff = initial + delta
        evaluation = evaluate_solution(
            cases,
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
        record = {"force_weight": force_weight, "solve": solve, "evaluation": evaluation}
        records.append(record)
        coeff_arrays[f"coeff_force_{force_weight:g}"] = coeff
    plot_path = args.output / "trace0_sparse_conservation_scan.png"
    coeff_path = args.output / "trace0_sparse_conservation_coefficients.npz"
    render(plot_path, records)
    np.savez_compressed(coeff_path, **coeff_arrays)
    report = {
        "parameters": {
            "tau": tau,
            "force_weights": [float(r["force_weight"]) for r in records],
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "force_probe_dt_old": float(args.force_probe_dt_old),
            "ncols": int(ncols),
            "n_algebraic_rows": int(n_alg),
            "n_force_rows": int(n_force),
            "active_count": int(np.count_nonzero(active)),
        },
        "definition": {
            "purpose": "Sparse high-resolution feasibility scan for trace0 + full conservation.",
            "trace0": "Imposed exactly by pointwise nullspace elimination of g^{mn} C_mn=0.",
            "force_weight": "Penalty weight on full conservation rows; this scan is not yet a hard-conservation saddle-point solve.",
            "next_if_successful": "Upgrade to matrix-free hard conservation/KKT after finding a stable force-weight regime.",
        },
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
    parser.add_argument("--force-weights", type=str, default="0,0.03,0.1,0.3,1,3,10")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
    parser.add_argument("--active-dilation", type=int, default=1)
    parser.add_argument("--lsqr-tol", type=float, default=1.0e-6)
    parser.add_argument("--lsqr-maxiter", type=int, default=1000)
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
