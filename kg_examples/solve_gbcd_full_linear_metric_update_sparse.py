from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import fit_equation_first_constrained_bcd as base
from analyze_c_terms_from_a_reference import metric_jets_full
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from solve_gbcd_full_linear_metric_update import (
    COMPONENTS,
    add_symmetric_component,
    affected_points,
    einstein_point,
    metric_variable_index,
    perturb_jets_for_metric_plus_variable,
    render,
    sym_to_vec,
)
from test_gbcd_v0_one_step_lambda import ATOM_NAMES as DEFAULT_ATOM_NAMES
from test_gbcd_v0_one_step_lambda import weighted_stats


ATOM_NAMES = tuple(DEFAULT_ATOM_NAMES)
MetricVariable = tuple[str, int, int, int, int]


def parse_metric_variable_slices(mode: str) -> tuple[str, ...]:
    if mode == "plus":
        return ("plus",)
    if mode == "center_plus":
        return ("center", "plus")
    if mode == "minus_center_plus":
        return ("minus", "center", "plus")
    raise ValueError(f"unknown metric variable slice mode: {mode}")


def metric_variable_index_for_slices(mask: np.ndarray, slices: tuple[str, ...]) -> list[MetricVariable]:
    variables: list[MetricVariable] = []
    active_points = np.argwhere(mask)
    for slice_name in slices:
        for i, j in active_points:
            for a, b in COMPONENTS:
                variables.append((slice_name, int(i), int(j), int(a), int(b)))
    return variables


def affected_points_for_metric_variable(var: MetricVariable, shape: tuple[int, int]) -> list[tuple[int, int]]:
    slice_name, i, j, _, _ = var
    if slice_name in {"minus", "plus"}:
        return affected_points(i, j, shape)
    pts = {(i, j)}
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)):
        ii = i + di
        jj = j + dj
        if 0 <= ii < shape[0] and 0 <= jj < shape[1]:
            pts.add((ii, jj))
    return sorted(pts)


def perturb_jets_for_metric_variable(
    point: tuple[int, int],
    var: MetricVariable,
    *,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    slice_name, qi, qj, a, b = var
    pi, pj = point
    dmetric = np.zeros((3, 3), dtype=float)
    ddg = np.zeros((3, 3, 3), dtype=float)
    dd2g = np.zeros((3, 3, 3, 3), dtype=float)

    if slice_name in {"minus", "plus"}:
        time_sign = 1.0 if slice_name == "plus" else -1.0
        if pi == qi and pj == qj:
            add_symmetric_component(ddg[0], a, b, time_sign / (2.0 * dt))
            add_symmetric_component(dd2g[0, 0], a, b, 1.0 / (dt * dt))
        if pj == qj and abs(pi - qi) == 1:
            space_sign = 1.0 if qi == pi + 1 else -1.0
            coeff = time_sign * space_sign / (4.0 * dx * dt)
            add_symmetric_component(dd2g[0, 1], a, b, coeff)
            add_symmetric_component(dd2g[1, 0], a, b, coeff)
        if pi == qi and abs(pj - qj) == 1:
            space_sign = 1.0 if qj == pj + 1 else -1.0
            coeff = time_sign * space_sign / (4.0 * dz * dt)
            add_symmetric_component(dd2g[0, 2], a, b, coeff)
            add_symmetric_component(dd2g[2, 0], a, b, coeff)
        return dmetric, ddg, dd2g

    if slice_name != "center":
        raise ValueError(f"unknown metric variable slice: {slice_name}")

    if pi == qi and pj == qj:
        add_symmetric_component(dmetric, a, b, 1.0)
        add_symmetric_component(dd2g[0, 0], a, b, -2.0 / (dt * dt))
        add_symmetric_component(dd2g[1, 1], a, b, -2.0 / (dx * dx))
        add_symmetric_component(dd2g[2, 2], a, b, -2.0 / (dz * dz))
    if pj == qj and abs(pi - qi) == 1:
        space_sign = 1.0 if qi == pi + 1 else -1.0
        add_symmetric_component(ddg[1], a, b, space_sign / (2.0 * dx))
        add_symmetric_component(dd2g[1, 1], a, b, 1.0 / (dx * dx))
    if pi == qi and abs(pj - qj) == 1:
        space_sign = 1.0 if qj == pj + 1 else -1.0
        add_symmetric_component(ddg[2], a, b, space_sign / (2.0 * dz))
        add_symmetric_component(dd2g[2, 2], a, b, 1.0 / (dz * dz))
    if abs(pi - qi) == 1 and abs(pj - qj) == 1:
        sign_x = 1.0 if qi == pi + 1 else -1.0
        sign_z = 1.0 if qj == pj + 1 else -1.0
        coeff = sign_x * sign_z / (4.0 * dx * dz)
        add_symmetric_component(dd2g[1, 2], a, b, coeff)
        add_symmetric_component(dd2g[2, 1], a, b, coeff)
    return dmetric, ddg, dd2g


def build_sparse_columns(
    *,
    metric_cov: np.ndarray,
    dg_base: np.ndarray,
    d2g_base: np.ndarray,
    variables: list[MetricVariable],
    row_points: list[tuple[int, int]],
    row_weights: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    eps: float,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray]:
    point_to_row = {pt: n for n, pt in enumerate(row_points)}
    col_rows: list[np.ndarray] = []
    col_vals: list[np.ndarray] = []
    col_norms = np.zeros(len(variables), dtype=float)
    for col, var in enumerate(variables):
        rows = []
        vals = []
        for pt in affected_points_for_metric_variable(var, metric_cov.shape[:2]):
            row_index = point_to_row.get(pt)
            if row_index is None:
                continue
            dmetric, ddg, dd2g = perturb_jets_for_metric_variable(pt, var, dt=dt, dx=dx, dz=dz)
            i, j = pt
            e_plus = einstein_point(
                metric_cov[i, j] + eps * dmetric,
                dg_base[i, j] + eps * ddg,
                d2g_base[i, j] + eps * dd2g,
            )
            e_minus = einstein_point(
                metric_cov[i, j] - eps * dmetric,
                dg_base[i, j] - eps * ddg,
                d2g_base[i, j] - eps * dd2g,
            )
            deriv = row_weights[row_index] * sym_to_vec((e_plus - e_minus) / (2.0 * eps))
            if np.any(np.isfinite(deriv)) and float(np.linalg.norm(deriv)) > 0.0:
                rows.append(row_index)
                vals.append(deriv)
        if rows:
            row_arr = np.asarray(rows, dtype=np.int32)
            val_arr = np.vstack(vals).astype(float)
            col_norms[col] = float(np.linalg.norm(val_arr))
        else:
            row_arr = np.zeros(0, dtype=np.int32)
            val_arr = np.zeros((0, len(COMPONENTS)), dtype=float)
            col_norms[col] = 0.0
        col_rows.append(row_arr)
        col_vals.append(val_arr)
    finite = col_norms[np.isfinite(col_norms) & (col_norms > 0.0)]
    floor = 1.0e-30 if finite.size == 0 else max(float(np.percentile(finite, 5.0)) * 1.0e-12, 1.0e-30)
    col_norms = np.where(col_norms > floor, col_norms, floor)
    return col_rows, col_vals, col_norms


def apply_a(col_rows: list[np.ndarray], col_vals: list[np.ndarray], x: np.ndarray, n_points: int) -> np.ndarray:
    y = np.zeros(n_points * len(COMPONENTS), dtype=float)
    width = len(COMPONENTS)
    for j, xj in enumerate(x):
        if xj == 0.0:
            continue
        rows = col_rows[j]
        vals = col_vals[j]
        for local, row in enumerate(rows):
            start = int(row) * width
            y[start : start + width] += xj * vals[local]
    return y


def apply_at(col_rows: list[np.ndarray], col_vals: list[np.ndarray], y: np.ndarray) -> np.ndarray:
    out = np.zeros(len(col_rows), dtype=float)
    width = len(COMPONENTS)
    for j, rows in enumerate(col_rows):
        vals = col_vals[j]
        acc = 0.0
        for local, row in enumerate(rows):
            start = int(row) * width
            acc += float(np.dot(vals[local], y[start : start + width]))
        out[j] = acc
    return out


def build_entry_arrays(
    col_rows: list[np.ndarray],
    col_vals: list[np.ndarray],
    col_norms: np.ndarray,
    n_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Build scalar COO-style arrays for the scaled operator B=A/column_norm.

    The original sparse-column representation is compact and easy to verify,
    but its matvecs spend most time in Python loops.  These arrays let NumPy do
    the repeated scatter/gather operations with bincount while preserving the
    same linear operator.
    """
    width = len(COMPONENTS)
    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    val_parts: list[np.ndarray] = []
    comp = np.arange(width, dtype=np.int64)
    for col, (rows, vals) in enumerate(zip(col_rows, col_vals)):
        if rows.size == 0:
            continue
        scaled = vals / col_norms[col]
        scalar_rows = (rows.astype(np.int64)[:, None] * width + comp[None, :]).reshape(-1)
        scalar_cols = np.full(scalar_rows.size, col, dtype=np.int64)
        scalar_vals = scaled.reshape(-1)
        good = np.isfinite(scalar_vals) & (scalar_vals != 0.0)
        if np.any(good):
            row_parts.append(scalar_rows[good])
            col_parts.append(scalar_cols[good])
            val_parts.append(scalar_vals[good].astype(float, copy=False))
    if not row_parts:
        return (
            np.zeros(0, dtype=np.int64),
            np.zeros(0, dtype=np.int64),
            np.zeros(0, dtype=float),
            int(n_points * width),
        )
    return (
        np.concatenate(row_parts),
        np.concatenate(col_parts),
        np.concatenate(val_parts),
        int(n_points * width),
    )


def apply_b_entries(row_idx: np.ndarray, col_idx: np.ndarray, vals: np.ndarray, x: np.ndarray, n_rows: int) -> np.ndarray:
    if vals.size == 0:
        return np.zeros(n_rows, dtype=float)
    return np.bincount(row_idx, weights=vals * x[col_idx], minlength=n_rows).astype(float, copy=False)


def apply_bt_entries(row_idx: np.ndarray, col_idx: np.ndarray, vals: np.ndarray, y: np.ndarray, n_cols: int) -> np.ndarray:
    if vals.size == 0:
        return np.zeros(n_cols, dtype=float)
    return np.bincount(col_idx, weights=vals * y[row_idx], minlength=n_cols).astype(float, copy=False)


def cg_solve(matvec, rhs: np.ndarray, *, tol: float, maxiter: int, preconditioner=None) -> tuple[np.ndarray, dict[str, object]]:
    x = np.zeros_like(rhs)
    r = rhs - matvec(x)
    z = preconditioner(r) if preconditioner is not None else r.copy()
    p = z.copy()
    rz_old = float(np.dot(r, z))
    rhs_norm = max(float(np.linalg.norm(rhs)), 1.0e-300)
    history = []
    r_norm = float(np.linalg.norm(r))
    if r_norm / rhs_norm < tol:
        return x, {"iterations": 0, "relative_residual": float(r_norm / rhs_norm), "history": history}
    for it in range(1, int(maxiter) + 1):
        ap = matvec(p)
        denom = float(np.dot(p, ap))
        if abs(denom) < 1.0e-300:
            break
        alpha = rz_old / denom
        x += alpha * p
        r -= alpha * ap
        r_norm = float(np.linalg.norm(r))
        rel = r_norm / rhs_norm
        if it <= 10 or it % 10 == 0:
            history.append([it, rel])
        if rel < tol:
            return x, {"iterations": it, "relative_residual": rel, "history": history}
        z = preconditioner(r) if preconditioner is not None else r.copy()
        rz_new = float(np.dot(r, z))
        beta = rz_new / max(rz_old, 1.0e-300)
        p = z + beta * p
        rz_old = rz_new
    return x, {"iterations": it if "it" in locals() else 0, "relative_residual": rel if "rel" in locals() else 0.0, "history": history}


def make_block_jacobi_preconditioner(
    col_rows: list[np.ndarray],
    col_vals: list[np.ndarray],
    col_norms: np.ndarray,
    variables: list[MetricVariable],
    ridge: float,
):
    blocks: list[tuple[np.ndarray, np.ndarray]] = []
    n = len(variables)
    start = 0
    while start < n:
        s0, i0, j0 = variables[start][0], variables[start][1], variables[start][2]
        end = start
        while end < n and variables[end][0] == s0 and variables[end][1] == i0 and variables[end][2] == j0:
            end += 1
        ids = np.arange(start, end, dtype=int)
        k = len(ids)
        block = (float(ridge) ** 2) * np.eye(k)
        row_maps = []
        for col in ids:
            mapping = {int(r): col_vals[col][m] / col_norms[col] for m, r in enumerate(col_rows[col])}
            row_maps.append(mapping)
        for a in range(k):
            for b in range(a, k):
                common = set(row_maps[a]).intersection(row_maps[b])
                val = 0.0
                for row in common:
                    val += float(np.dot(row_maps[a][row], row_maps[b][row]))
                block[a, b] += val
                if a != b:
                    block[b, a] += val
        inv = np.linalg.pinv(block, rcond=1.0e-12, hermitian=True)
        blocks.append((ids, inv))
        start = end

    def apply(vec: np.ndarray) -> np.ndarray:
        out = np.zeros_like(vec)
        for ids, inv in blocks:
            out[ids] = inv @ vec[ids]
        return out

    return apply, {"block_count": len(blocks), "block_size_min": int(min(len(ids) for ids, _ in blocks)), "block_size_max": int(max(len(ids) for ids, _ in blocks))}


def lsqr_solve(aprod, atprod, b: np.ndarray, n: int, *, tol: float, maxiter: int) -> tuple[np.ndarray, dict[str, object]]:
    x = np.zeros(n, dtype=float)
    u = b.astype(float).copy()
    beta = float(np.linalg.norm(u))
    if beta > 0.0:
        u /= beta
    v = atprod(u)
    alpha = float(np.linalg.norm(v))
    if alpha > 0.0:
        v /= alpha
    w = v.copy()
    phibar = beta
    rhobar = alpha
    bnorm = max(beta, 1.0e-300)
    history = []
    rel = 1.0
    for it in range(1, int(maxiter) + 1):
        u = aprod(v) - alpha * u
        beta = float(np.linalg.norm(u))
        if beta > 0.0:
            u /= beta
        v = atprod(u) - beta * v
        alpha = float(np.linalg.norm(v))
        if alpha > 0.0:
            v /= alpha
        rho = float(np.hypot(rhobar, beta))
        if rho <= 1.0e-300:
            break
        c = rhobar / rho
        s = beta / rho
        theta = s * alpha
        rhobar = -c * alpha
        phi = c * phibar
        phibar = s * phibar
        x += (phi / rho) * w
        w = v - (theta / rho) * w
        rel = abs(phibar) / bnorm
        if it <= 10 or it % 10 == 0:
            history.append([it, float(rel)])
        if rel < tol:
            break
    return x, {"iterations": it if "it" in locals() else 0, "relative_residual_estimate": float(rel), "history": history}


def _npz_scalar_string(data: np.lib.npyio.NpzFile, key: str, default: str) -> str:
    if key not in data.files:
        return default
    value = np.asarray(data[key])
    if value.shape == ():
        return str(value.item())
    if value.size == 1:
        return str(value.reshape(-1)[0])
    return default


def _npz_atom_names(data: np.lib.npyio.NpzFile, case, atom_family: str) -> tuple[str, ...]:
    if "atom_names" in data.files:
        return tuple(str(item) for item in np.asarray(data["atom_names"]).tolist())
    return tuple(tensor_atoms(case, atom_family).keys())


def load_lambda_cov(cases: dict[int, object], coeff_path: Path, atom_family: str = "auto") -> tuple[np.ndarray, np.ndarray]:
    coeff_data = np.load(coeff_path)
    if atom_family == "auto":
        atom_family = _npz_scalar_string(coeff_data, "atom_family", "matter4")
    atom_names = _npz_atom_names(coeff_data, cases[0], atom_family)
    fields = {}
    for key, suffix in [(-1, "m"), (0, "0"), (1, "p")]:
        arr = np.zeros(cases[0].rho_a.shape + (len(atom_names),), dtype=float)
        for apos, atom in enumerate(atom_names):
            arr[..., apos] = coeff_data[f"{atom}_{suffix}"]
        fields[key] = arr
    atoms = tensor_atoms(cases[0], atom_family)
    cov0 = np.zeros_like(cases[0].r_need)
    for apos, atom in enumerate(atom_names):
        cov0 += fields[0][..., apos, None, None] * atoms[atom]
    return cov0, np.asarray(coeff_data["mask"], dtype=bool)


def default_coeff_path(tau: float, n: int) -> Path:
    if n != 96:
        raise ValueError("default coefficient paths are currently available only for n=96")
    if abs(tau) < 1.0e-12:
        tag = "tau0"
    elif tau < 0.0:
        tag = "taum3p5"
    else:
        tag = "taup3p5"
    return Path(f"visualizations/equation_first_gbcd_principal_constraint_projection_n96_{tag}/principal_constraint_coefficients.npz")


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

    coeff_path = args.coefficients or default_coeff_path(float(args.tau), int(args.full_resolution))
    cov0, fit_mask = load_lambda_cov(cases, coeff_path, str(args.atom_family))
    residual = cov0 - cases[0].r_need

    active_metric = base.dilate_mask_8(fit_mask, int(args.metric_active_dilation))
    active_metric[0, :] = False
    active_metric[-1, :] = False
    active_metric[:, 0] = False
    active_metric[:, -1] = False
    metric_slices = parse_metric_variable_slices(str(args.metric_variable_slices))
    variables = metric_variable_index_for_slices(active_metric, metric_slices)
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
    dg_base, d2g_base = metric_jets_full(full_case.geom_m.metric_cov, full_case.geom_0.metric_cov, full_case.geom_p.metric_cov, dt, dx, dz)
    col_rows, col_vals, col_norms = build_sparse_columns(
        metric_cov=full_case.geom_0.metric_cov,
        dg_base=dg_base,
        d2g_base=d2g_base,
        variables=variables,
        row_points=row_points,
        row_weights=row_weights,
        dt=dt,
        dx=dx,
        dz=dz,
        eps=float(args.linear_eps),
    )

    entry_arrays_info: dict[str, object] | None = None
    if args.operator_backend == "entry-arrays":
        entry_row, entry_col, entry_val, entry_n_rows = build_entry_arrays(col_rows, col_vals, col_norms, len(row_points))
        entry_arrays_info = {
            "scalar_nonzeros": int(entry_val.size),
            "rows": int(entry_n_rows),
            "cols": int(len(variables)),
        }

        def b_apply(y: np.ndarray) -> np.ndarray:
            return apply_b_entries(entry_row, entry_col, entry_val, y, entry_n_rows)

        def bt_apply(v: np.ndarray) -> np.ndarray:
            return apply_bt_entries(entry_row, entry_col, entry_val, v, len(variables))

    else:

        def b_apply(y: np.ndarray) -> np.ndarray:
            return apply_a(col_rows, col_vals, y / col_norms, len(row_points))

        def bt_apply(v: np.ndarray) -> np.ndarray:
            return apply_at(col_rows, col_vals, v) / col_norms

    bt_rhs = bt_apply(rhs)

    def normal_matvec(y: np.ndarray) -> np.ndarray:
        return bt_apply(b_apply(y)) + float(args.ridge) ** 2 * y

    preconditioner = None
    preconditioner_info = {"type": "none"}
    if args.preconditioner == "diagonal":
        if args.operator_backend == "entry-arrays":
            diag = np.bincount(entry_col, weights=entry_val * entry_val, minlength=len(variables)).astype(float, copy=False)
        else:
            diag = np.zeros(len(variables), dtype=float)
            for col, vals in enumerate(col_vals):
                if vals.size:
                    diag[col] = float(np.sum((vals / col_norms[col]) ** 2))
        diag += float(args.ridge) ** 2
        floor_diag = max(float(np.percentile(diag[np.isfinite(diag) & (diag > 0.0)], 5.0)) * 1.0e-12, 1.0e-300)
        inv_diag = 1.0 / np.maximum(diag, floor_diag)

        def preconditioner(vec: np.ndarray) -> np.ndarray:
            return inv_diag * vec

        preconditioner_info = {
            "type": "diagonal",
            "diag_min": float(np.min(diag)),
            "diag_p50": float(np.percentile(diag, 50.0)),
            "diag_p95": float(np.percentile(diag, 95.0)),
        }
    elif args.preconditioner == "block-jacobi":
        preconditioner, info = make_block_jacobi_preconditioner(col_rows, col_vals, col_norms, variables, float(args.ridge))
        preconditioner_info = {"type": "block-jacobi", **info}
    if args.solver == "cg":
        y, solve_info = cg_solve(
            normal_matvec,
            bt_rhs,
            tol=float(args.cg_tol),
            maxiter=int(args.cg_maxiter),
            preconditioner=preconditioner,
        )
        solve_info["solver"] = "cg_normal_equations"
    elif args.solver == "lsqr":
        m = rhs.size
        n = bt_rhs.size

        def aprod_aug(y_in: np.ndarray) -> np.ndarray:
            return np.concatenate([b_apply(y_in), float(args.ridge) * y_in])

        def atprod_aug(u_in: np.ndarray) -> np.ndarray:
            return bt_apply(u_in[:m]) + float(args.ridge) * u_in[m : m + n]

        y, solve_info = lsqr_solve(
            aprod_aug,
            atprod_aug,
            np.concatenate([rhs, np.zeros(n, dtype=float)]),
            n,
            tol=float(args.cg_tol),
            maxiter=int(args.cg_maxiter),
        )
        solve_info["solver"] = "lsqr_augmented"
    else:
        raise ValueError(f"unknown solver: {args.solver}")
    delta_vec = y / col_norms
    weighted_pred = apply_a(col_rows, col_vals, delta_vec, len(row_points))

    delta_metrics = {
        "minus": np.zeros_like(full_case.geom_m.metric_cov),
        "center": np.zeros_like(full_case.geom_0.metric_cov),
        "plus": np.zeros_like(full_case.geom_p.metric_cov),
    }
    for value, (slice_name, i, j, a, b) in zip(delta_vec, variables):
        delta_metrics[slice_name][i, j, a, b] += value
        if a != b:
            delta_metrics[slice_name][i, j, b, a] += value
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
    exact_r_need = ein_corr - cases[0].source_tilde

    predicted = np.zeros_like(residual)
    width = len(COMPONENTS)
    for n, (i, j) in enumerate(row_points):
        predicted[i, j] = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]
        )
        vals = weighted_pred[n * width : (n + 1) * width] / max(row_weights[n], 1.0e-300)
        for value, (a, b) in zip(vals, COMPONENTS):
            predicted[i, j, a, b] = value
            predicted[i, j, b, a] = value

    before_rel = tensor_norm(residual) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    linear_after_rel = tensor_norm(residual - predicted) / np.maximum(tensor_norm(cases[0].r_need + predicted), 1.0e-300)
    exact_after_rel = tensor_norm(cov0 - exact_r_need) / np.maximum(tensor_norm(exact_r_need), 1.0e-300)
    metric_norm_plus = tensor_norm(full_case.geom_p.metric_cov)
    delta_norm_plus = tensor_norm(delta_metrics["plus"])
    metric_norm_center = tensor_norm(full_case.geom_0.metric_cov)
    delta_norm_center = tensor_norm(delta_metrics["center"])
    weights_eval = np.sqrt(np.maximum(cases[0].rho_a[fit_mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    plot_path = args.output / "gbcd_full_linear_metric_update_sparse.png"
    render(plot_path, cases[0].x, cases[0].z, fit_mask, before_rel, linear_after_rel, exact_after_rel)
    data_path = args.output / "full_linear_metric_update_sparse.npz"
    np.savez_compressed(
        data_path,
        delta_metric_minus=delta_metrics["minus"],
        delta_metric_center=delta_metrics["center"],
        delta_metric_plus=delta_metrics["plus"],
        corrected_metric_minus=corrected_metric_minus,
        corrected_metric_center=corrected_metric_center,
        corrected_metric_plus=corrected_metric_plus,
        before_rel=before_rel,
        linear_after_rel=linear_after_rel,
        exact_after_rel=exact_after_rel,
        mask=fit_mask,
    )
    normal_residual = normal_matvec(y) - bt_rhs
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "metric_active_dilation": int(args.metric_active_dilation),
            "dt_ev_inv": float(dt),
            "dx_ev_inv": float(dx),
            "dz_ev_inv": float(dz),
            "linear_eps": float(args.linear_eps),
            "ridge": float(args.ridge),
            "solver": args.solver,
            "preconditioner": args.preconditioner,
            "operator_backend": args.operator_backend,
            "cg_tol": float(args.cg_tol),
            "cg_maxiter": int(args.cg_maxiter),
            "atom_family": str(args.atom_family),
            "metric_variable_slices": str(args.metric_variable_slices),
            "coefficients": str(coeff_path.resolve()),
        },
        "definition": {
            "goal": "Sparse-column / matrix-free version of full-linear metric update. It avoids storing the dense rows x variables matrix.",
            "normal_equation": "(B^T B + ridge^2 I)y = B^T b with B=A/column_norm and delta_metric=y/column_norm.",
        },
        "counts": {
            "fit_points": int(np.count_nonzero(fit_mask)),
            "metric_active_points": int(np.count_nonzero(active_metric)),
            "metric_slice_count": int(len(metric_slices)),
            "variables": int(len(variables)),
            "rows": int(len(row_points) * len(COMPONENTS)),
            "nonzero_column_blocks": int(sum(len(r) for r in col_rows)),
        },
        "solve": {
            **solve_info,
            "preconditioner": preconditioner_info,
            "operator_backend_info": entry_arrays_info or {"type": "python-column-loops"},
            "normal_relative_residual": float(np.linalg.norm(normal_residual) / max(float(np.linalg.norm(bt_rhs)), 1.0e-300)),
        },
        "evaluation": {
            "before_relative_residual": weighted_stats(before_rel[fit_mask], weights_eval),
            "linear_after_relative_residual": weighted_stats(linear_after_rel[fit_mask], weights_eval),
            "exact_after_relative_residual": weighted_stats(exact_after_rel[fit_mask], weights_eval),
            "relative_delta_metric_center_norm": weighted_stats(
                delta_norm_center[active_metric] / np.maximum(metric_norm_center[active_metric], 1.0e-300)
            ),
            "relative_delta_metric_plus_norm": weighted_stats(
                delta_norm_plus[active_metric] / np.maximum(metric_norm_plus[active_metric], 1.0e-300)
            ),
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
    parser.add_argument("--coefficients", type=Path, default=None)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--metric-active-dilation", type=int, default=1)
    parser.add_argument(
        "--metric-variable-slices",
        choices=["plus", "center_plus", "minus_center_plus"],
        default="plus",
    )
    parser.add_argument("--linear-eps", type=float, default=1.0e-8)
    parser.add_argument("--ridge", type=float, default=1.0e-6)
    parser.add_argument("--cg-tol", type=float, default=1.0e-8)
    parser.add_argument("--cg-maxiter", type=int, default=300)
    parser.add_argument("--solver", choices=["cg", "lsqr"], default="lsqr")
    parser.add_argument("--operator-backend", choices=["entry-arrays", "python-column-loops"], default="entry-arrays")
    parser.add_argument(
        "--atom-family",
        choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"],
        default="auto",
    )
    parser.add_argument("--preconditioner", choices=["none", "diagonal", "block-jacobi"], default="none")
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
