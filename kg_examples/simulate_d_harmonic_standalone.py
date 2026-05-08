from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_c_terms_from_a_reference import metric_jets_full
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference, real_array
from simulate_d_reduced_dynamic_same_initial import rk4_matter_step
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det, stress_tensor_tilde
from solve_gbcd_full_linear_metric_update import COMPONENTS, sym_to_vec
from solve_gbcd_full_linear_metric_update_sparse import (
    apply_b_entries,
    apply_bt_entries,
    build_entry_arrays,
    build_sparse_columns,
    cg_solve,
    metric_variable_index_for_slices,
)
from fit_equation_first_tensor_couplings import tensor_atoms
from diagnose_gbcd_trace_local_closure import nullspace_basis, solve_local, trace_coefficients
from fit_gbcd_principal_constraint_projection_sparse import append_row, solve_sparse_lsqr, sparse_apply


def relative_l1(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(np.abs(candidate[valid] - reference[valid])) / denom)


def weighted_relative_l1(reference: np.ndarray, candidate: np.ndarray, weight: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(weight[valid] * np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(weight[valid] * np.abs(candidate[valid] - reference[valid])) / denom)


def discriminant_guard(discriminant: np.ndarray, mask: np.ndarray, *, tolerance: float, fraction_tolerance: float) -> dict[str, float | bool]:
    if not np.any(mask):
        return {
            "stop": False,
            "min": 0.0,
            "negative_fraction": 0.0,
            "large_negative_fraction": 0.0,
        }
    values = np.asarray(discriminant[mask], dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {
            "stop": False,
            "min": 0.0,
            "negative_fraction": 0.0,
            "large_negative_fraction": 0.0,
        }
    tol = abs(float(tolerance))
    neg = values < 0.0
    large_neg = values < -tol
    negative_fraction = float(np.count_nonzero(neg) / values.size)
    large_negative_fraction = float(np.count_nonzero(large_neg) / values.size)
    return {
        "stop": bool(large_negative_fraction > max(float(fraction_tolerance), 0.0)),
        "min": float(np.min(values)),
        "negative_fraction": negative_fraction,
        "large_negative_fraction": large_negative_fraction,
    }


def lower_bound_guard(values: np.ndarray, mask: np.ndarray, *, tolerance: float, fraction_tolerance: float) -> dict[str, float | bool]:
    if not np.any(mask):
        return {
            "stop": False,
            "min": 0.0,
            "negative_fraction": 0.0,
            "large_negative_fraction": 0.0,
        }
    selected = np.asarray(values[mask], dtype=float)
    selected = selected[np.isfinite(selected)]
    if selected.size == 0:
        return {
            "stop": False,
            "min": 0.0,
            "negative_fraction": 0.0,
            "large_negative_fraction": 0.0,
        }
    tol = abs(float(tolerance))
    neg = selected < 0.0
    large_neg = selected < -tol
    negative_fraction = float(np.count_nonzero(neg) / selected.size)
    large_negative_fraction = float(np.count_nonzero(large_neg) / selected.size)
    return {
        "stop": bool(large_negative_fraction > max(float(fraction_tolerance), 0.0)),
        "min": float(np.min(selected)),
        "negative_fraction": negative_fraction,
        "large_negative_fraction": large_negative_fraction,
    }


def stats(values: np.ndarray, *, abs_value: bool = True) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if abs_value:
        vals = np.abs(vals)
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "p99": float(np.percentile(vals, 99.0)),
        "max": float(np.max(vals)),
    }


def infer_tau_from_package(path: Path) -> float:
    text = str(path)
    if "taum3p5" in text:
        return -3.5
    if "taup3p5" in text:
        return 3.5
    if "tau0" in text:
        return 0.0
    return 0.0


def pullback_rho_to_a(measure_tilde: np.ndarray, u_t: np.ndarray, u_x: np.ndarray, u_z: np.ndarray, mass: float, x_floor: float, mask: np.ndarray) -> np.ndarray:
    x_g = u_t * u_t - u_x * u_x - u_z * u_z
    return np.where(mask & (np.abs(x_g) > x_floor), mass * mass * measure_tilde / np.maximum(np.abs(x_g), x_floor), 0.0)


def current_measure(metric_cov: np.ndarray, rho_tilde: np.ndarray) -> np.ndarray:
    return safe_sqrt_abs_det(metric_cov) * rho_tilde


def recompute_r_need(
    metric_minus: np.ndarray,
    metric_center: np.ndarray,
    metric_plus: np.ndarray,
    rho_tilde: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    *,
    mass: float,
    mp: float,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    dg, d2g = metric_jets_full(metric_minus, metric_center, metric_plus, dt, dx, dz)
    geom, _ = make_time_geometry(metric_center, metric_minus, metric_plus, dt, dx, dz, with_dgamma=False)
    # make_time_geometry and metric_jets_full use the same finite-difference
    # convention; dg/d2g is kept here only to make this function's inputs clear.
    del dg, d2g
    ein = geom.ricci - 0.5 * metric_center * geom.r_scalar[..., None, None]
    stress, tilde_x = stress_tensor_tilde(metric_center, geom.metric_inv, rho_tilde, u_t, u_x, u_z, m=mass)
    r_need = ein - stress / (mp * mp)
    return r_need, geom.r_scalar, geom.metric_inv, real_array(tilde_x)


def trace_project_auxiliary(aux: np.ndarray, metric_inv: np.ndarray, metric_cov: np.ndarray) -> np.ndarray:
    trace = np.einsum("...ab,...ab->...", metric_inv, aux, optimize=True)
    dim = aux.shape[-1]
    return aux - (trace / float(dim))[..., None, None] * metric_cov


class RuntimeCase:
    def __init__(
        self,
        *,
        metric_cov: np.ndarray,
        metric_inv: np.ndarray,
        r_need: np.ndarray,
        u_cov: np.ndarray,
        r_cov: np.ndarray,
        rho_a: np.ndarray,
    ) -> None:
        self.metric_cov = metric_cov
        self.metric_inv = metric_inv
        self.r_need = r_need
        self.u_cov = u_cov
        self.r_cov = r_cov
        self.rho_a = rho_a


def ddx(field: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * dx)


def ddz(field: np.ndarray, dz: float) -> np.ndarray:
    return (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dz)


def local_auxiliary_solve(
    *,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    r_need: np.ndarray,
    rho_tilde: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    dx: float,
    dz: float,
    mask: np.ndarray,
    trace0: bool,
) -> np.ndarray:
    root = np.sqrt(np.maximum(rho_tilde, 1.0e-300))
    r_cov = np.stack([np.zeros_like(root), ddx(root, dx) / np.maximum(root, 1.0e-300), ddz(root, dz) / np.maximum(root, 1.0e-300)], axis=-1)
    case = RuntimeCase(
        metric_cov=metric_cov,
        metric_inv=metric_inv,
        r_need=r_need,
        u_cov=np.stack([u_t, u_x, u_z], axis=-1),
        r_cov=r_cov,
        rho_a=rho_tilde,
    )
    atoms = tensor_atoms(case, "matter4")
    trace_rows = trace_coefficients(case) if trace0 else None
    out = np.zeros_like(r_need)
    components = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    for i_raw, j_raw in np.argwhere(mask):
        i = int(i_raw)
        j = int(j_raw)
        a = np.zeros((len(components), 4), dtype=float)
        b = np.zeros(len(components), dtype=float)
        for cpos, (mu, nu) in enumerate(components):
            b[cpos] = r_need[i, j, mu, nu]
            for apos, name in enumerate(("g", "uu", "rr", "ur")):
                a[cpos, apos] = atoms[name][i, j, mu, nu]
        constraint = trace_rows[i, j] if trace0 and trace_rows is not None else None
        coeff, res = solve_local(a, b, constraint)
        if not np.all(np.isfinite(coeff)):
            continue
        for apos, name in enumerate(("g", "uu", "rr", "ur")):
            out[i, j] += coeff[apos] * atoms[name][i, j]
    return out


def _auxiliary_basis(metric_inv: np.ndarray, atoms: dict[str, np.ndarray], mask: np.ndarray, *, trace0: bool) -> tuple[np.ndarray, np.ndarray]:
    atom_names = ("g", "uu", "rr", "ur")
    free_dim = len(atom_names) - 1 if trace0 else len(atom_names)
    basis = np.zeros(mask.shape + (len(atom_names), free_dim), dtype=float)
    valid = mask.copy()
    if not trace0:
        eye = np.eye(len(atom_names), dtype=float)
        basis[mask] = eye
        return basis, valid

    trace_rows = np.zeros(mask.shape + (len(atom_names),), dtype=float)
    for apos, name in enumerate(atom_names):
        trace_rows[..., apos] = np.einsum("...ab,...ab->...", metric_inv, atoms[name], optimize=True)
    valid[:, :] = False
    for i_raw, j_raw in np.argwhere(mask):
        i = int(i_raw)
        j = int(j_raw)
        n = nullspace_basis(trace_rows[i, j])
        if n is None or n.shape != (len(atom_names), free_dim):
            continue
        basis[i, j] = n
        valid[i, j] = True
    return basis, valid


def _auxiliary_column_index(valid: np.ndarray, free_dim: int) -> tuple[np.ndarray, int]:
    col = -np.ones(valid.shape + (free_dim,), dtype=int)
    ncols = 0
    for i_raw, j_raw in np.argwhere(valid):
        i = int(i_raw)
        j = int(j_raw)
        for k in range(free_dim):
            col[i, j, k] = ncols
            ncols += 1
    return col, ncols


def _add_auxiliary_component_entries(
    entries: list[tuple[int, float]],
    *,
    col: np.ndarray,
    basis: np.ndarray,
    atoms: dict[str, np.ndarray],
    i: int,
    j: int,
    mu: int,
    nu: int,
    prefactor: float,
) -> None:
    if not np.isfinite(prefactor) or prefactor == 0.0:
        return
    atom_names = ("g", "uu", "rr", "ur")
    for k in range(col.shape[-1]):
        c = int(col[i, j, k])
        if c < 0:
            continue
        value = 0.0
        for apos, name in enumerate(atom_names):
            value += float(basis[i, j, apos, k]) * float(atoms[name][i, j, mu, nu])
        value *= float(prefactor)
        if np.isfinite(value) and value != 0.0:
            entries.append((c, value))


def _coefficients_to_auxiliary(
    coeff: np.ndarray,
    *,
    col: np.ndarray,
    basis: np.ndarray,
    atoms: dict[str, np.ndarray],
    shape: tuple[int, int],
) -> np.ndarray:
    atom_names = ("g", "uu", "rr", "ur")
    out = np.zeros(shape + (3, 3), dtype=float)
    for i_raw, j_raw in np.argwhere(col[..., 0] >= 0):
        i = int(i_raw)
        j = int(j_raw)
        lam = np.zeros(len(atom_names), dtype=float)
        for k in range(col.shape[-1]):
            c = int(col[i, j, k])
            if c >= 0:
                lam += float(coeff[c]) * basis[i, j, :, k]
        for apos, name in enumerate(atom_names):
            out[i, j] += lam[apos] * atoms[name][i, j]
    return out


def _auxiliary_project_conservation(
    initial: np.ndarray,
    force_cols: list[np.ndarray],
    force_vals: list[np.ndarray],
    force_rhs: list[float],
    ncols: int,
    *,
    tol: float,
    maxiter: int,
) -> tuple[np.ndarray, dict[str, object]]:
    target = np.asarray(force_rhs, dtype=float)
    initial_residual = sparse_apply(force_cols, force_vals, initial) - target
    delta_rhs = -initial_residual
    delta, solve_info = solve_sparse_lsqr(force_cols, force_vals, delta_rhs, ncols, tol, maxiter)
    coeff = initial + delta
    final_residual = sparse_apply(force_cols, force_vals, coeff) - target
    return coeff, {
        **solve_info,
        "mode": "hard_conservation_projection",
        "initial_force_norm": float(np.linalg.norm(initial_residual)),
        "projected_force_norm": float(np.linalg.norm(final_residual)),
        "projected_force_relative_norm": float(
            np.linalg.norm(final_residual) / max(float(np.linalg.norm(initial_residual)), 1.0e-300)
        ),
        "delta_norm": float(np.linalg.norm(delta)),
        "initial_norm": float(np.linalg.norm(initial)),
        "delta_over_initial_norm": float(np.linalg.norm(delta) / max(float(np.linalg.norm(initial)), 1.0e-300)),
    }


def _auxiliary_alm_conservation(
    initial: np.ndarray,
    alg_cols: list[np.ndarray],
    alg_vals: list[np.ndarray],
    alg_rhs: list[float],
    force_cols: list[np.ndarray],
    force_vals: list[np.ndarray],
    force_rhs: list[float],
    ncols: int,
    *,
    mu0: float,
    growth: float,
    outer: int,
    tol: float,
    maxiter: int,
) -> tuple[np.ndarray, dict[str, object]]:
    alg_target = np.asarray(alg_rhs, dtype=float)
    force_target = np.asarray(force_rhs, dtype=float)
    alg_delta_rhs = alg_target - sparse_apply(alg_cols, alg_vals, initial)
    force_delta_rhs = force_target - sparse_apply(force_cols, force_vals, initial)
    lam = np.zeros(len(force_cols), dtype=float)
    coeff = initial.copy()
    history: list[dict[str, float]] = []
    mu = float(mu0)
    for outer_it in range(1, int(outer) + 1):
        sqrt_mu = float(np.sqrt(max(mu, 0.0)))
        row_cols = list(alg_cols)
        row_vals = list(alg_vals)
        rhs = list(alg_delta_rhs)
        shifted_force_rhs = force_delta_rhs - lam / max(mu, 1.0e-300)
        for cols, vals, b in zip(force_cols, force_vals, shifted_force_rhs):
            row_cols.append(cols)
            row_vals.append(sqrt_mu * vals)
            rhs.append(sqrt_mu * float(b))
        delta, solve_info = solve_sparse_lsqr(row_cols, row_vals, np.asarray(rhs, dtype=float), ncols, tol, maxiter)
        coeff = initial + delta
        force_residual = sparse_apply(force_cols, force_vals, coeff) - force_target
        alg_residual = sparse_apply(alg_cols, alg_vals, coeff) - alg_target
        lam = lam + mu * force_residual
        history.append(
            {
                "outer": float(outer_it),
                "mu": float(mu),
                "lsqr_iterations": float(solve_info.get("iterations", 0)),
                "force_norm": float(np.linalg.norm(force_residual)),
                "algebraic_relative_norm": float(
                    np.linalg.norm(alg_residual) / max(float(np.linalg.norm(alg_target)), 1.0e-300)
                ),
                "delta_over_initial_norm": float(
                    np.linalg.norm(delta) / max(float(np.linalg.norm(initial)), 1.0e-300)
                ),
            }
        )
        mu *= float(growth)
    final_force = sparse_apply(force_cols, force_vals, coeff) - force_target
    final_alg = sparse_apply(alg_cols, alg_vals, coeff) - alg_target
    initial_force = sparse_apply(force_cols, force_vals, initial) - force_target
    return coeff, {
        "mode": "augmented_lagrangian_conservation",
        "outer_iterations": int(outer),
        "mu0": float(mu0),
        "growth": float(growth),
        "projected_force_norm": float(np.linalg.norm(final_force)),
        "projected_force_relative_norm": float(
            np.linalg.norm(final_force) / max(float(np.linalg.norm(initial_force)), 1.0e-300)
        ),
        "algebraic_relative_residual_norm": float(
            np.linalg.norm(final_alg) / max(float(np.linalg.norm(alg_target)), 1.0e-300)
        ),
        "delta_over_initial_norm": float(
            np.linalg.norm(coeff - initial) / max(float(np.linalg.norm(initial)), 1.0e-300)
        ),
        "history": history,
    }


def _initial_auxiliary_coefficients(
    *,
    col: np.ndarray,
    basis: np.ndarray,
    atoms: dict[str, np.ndarray],
    r_need: np.ndarray,
) -> np.ndarray:
    ncols = int(np.max(col)) + 1 if np.any(col >= 0) else 0
    coeff = np.zeros(ncols, dtype=float)
    atom_names = ("g", "uu", "rr", "ur")
    components = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    for i_raw, j_raw in np.argwhere(col[..., 0] >= 0):
        i = int(i_raw)
        j = int(j_raw)
        a = np.zeros((len(components), col.shape[-1]), dtype=float)
        b = np.zeros(len(components), dtype=float)
        for cpos, (mu, nu) in enumerate(components):
            b[cpos] = r_need[i, j, mu, nu]
            for k in range(col.shape[-1]):
                for apos, name in enumerate(atom_names):
                    a[cpos, k] += basis[i, j, apos, k] * atoms[name][i, j, mu, nu]
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        sol = np.linalg.lstsq(a, b, rcond=1.0e-12)[0]
        if not np.all(np.isfinite(sol)):
            continue
        for k, value in enumerate(sol):
            c = int(col[i, j, k])
            if c >= 0:
                coeff[c] = float(value)
    return coeff


def global_auxiliary_solve(
    *,
    metric_minus: np.ndarray,
    metric_center: np.ndarray,
    metric_plus: np.ndarray,
    metric_inv: np.ndarray,
    r_need: np.ndarray,
    rho_tilde: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    previous_auxiliary: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    mask: np.ndarray,
    trace0: bool,
    solve_mode: str,
    target_floor_frac: float,
    conservation_weight: float,
    lsqr_tol: float,
    lsqr_maxiter: int,
    alm_mu0: float,
    alm_growth: float,
    alm_outer: int,
) -> tuple[np.ndarray, dict[str, object]]:
    """Runtime C-sector solve.

    This is the first standalone update that does not freeze C.  It enforces the
    trace-free condition by eliminating one local coefficient and enforces
    conservation as a weighted backward-time/centered-space covariant-divergence
    equation.  The time derivative is intentionally one-sided, so the update can
    be used causally during evolution.
    """
    root = np.sqrt(np.maximum(rho_tilde, 1.0e-300))
    r_cov = np.stack(
        [
            np.zeros_like(root),
            ddx(root, dx) / np.maximum(root, 1.0e-300),
            ddz(root, dz) / np.maximum(root, 1.0e-300),
        ],
        axis=-1,
    )
    case = RuntimeCase(
        metric_cov=metric_center,
        metric_inv=metric_inv,
        r_need=r_need,
        u_cov=np.stack([u_t, u_x, u_z], axis=-1),
        r_cov=r_cov,
        rho_a=rho_tilde,
    )
    atoms = tensor_atoms(case, "matter4")
    basis, valid = _auxiliary_basis(metric_inv, atoms, mask, trace0=trace0)
    free_dim = basis.shape[-1]
    col, ncols = _auxiliary_column_index(valid, free_dim)
    if ncols == 0:
        return np.zeros_like(r_need), {"enabled": True, "reason": "empty global C solve"}

    row_cols: list[np.ndarray] = []
    row_vals: list[np.ndarray] = []
    rhs: list[float] = []
    alg_cols: list[np.ndarray] = []
    alg_vals: list[np.ndarray] = []
    alg_rhs: list[float] = []
    force_cols: list[np.ndarray] = []
    force_vals: list[np.ndarray] = []
    force_rhs: list[float] = []
    components = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
    row_target = np.sqrt(np.sum(np.asarray([r_need[..., mu, nu] for mu, nu in components]) ** 2, axis=0))
    finite_target = row_target[valid & np.isfinite(row_target)]
    target_floor = float(target_floor_frac) * max(
        float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0,
        1.0e-300,
    )
    max_rho = max(float(np.nanmax(np.maximum(rho_tilde[mask], 0.0))) if np.any(mask) else 0.0, 1.0e-300)
    for i_raw, j_raw in np.argwhere(valid):
        i = int(i_raw)
        j = int(j_raw)
        weight = np.sqrt(max(float(rho_tilde[i, j]), 0.0) / max_rho) / max(float(row_target[i, j]), target_floor)
        for mu, nu in components:
            entries: list[tuple[int, float]] = []
            _add_auxiliary_component_entries(
                entries,
                col=col,
                basis=basis,
                atoms=atoms,
                i=i,
                j=j,
                mu=mu,
                nu=nu,
                prefactor=weight,
            )
            append_row(row_cols, row_vals, rhs, entries, weight * float(r_need[i, j, mu, nu]))
            append_row(alg_cols, alg_vals, alg_rhs, entries, weight * float(r_need[i, j, mu, nu]))

    geom, _ = make_time_geometry(metric_center, metric_minus, metric_plus, dt, dx, dz, with_dgamma=False)
    h = max(min(abs(float(dt)), abs(float(dx)), abs(float(dz))), 1.0e-300)
    force_scale = max(target_floor / h, 1.0e-300)
    force_mask = valid.copy()
    force_mask[0, :] = False
    force_mask[-1, :] = False
    force_mask[:, 0] = False
    force_mask[:, -1] = False
    force_mask &= np.roll(valid, 1, 0) & np.roll(valid, -1, 0) & np.roll(valid, 1, 1) & np.roll(valid, -1, 1)
    for i_raw, j_raw in np.argwhere(force_mask):
        i = int(i_raw)
        j = int(j_raw)
        weight = float(conservation_weight) * np.sqrt(max(float(rho_tilde[i, j]), 0.0) / max_rho) / force_scale
        if not np.isfinite(weight) or weight == 0.0:
            continue
        for nu in range(3):
            entries: list[tuple[int, float]] = []
            rhs_value = 0.0
            for mu in range(3):
                for alpha in range(3):
                    pref = float(geom.metric_inv[i, j, mu, alpha])
                    if not np.isfinite(pref) or pref == 0.0:
                        continue
                    if alpha == 0:
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i,
                            j=j,
                            mu=mu,
                            nu=nu,
                            prefactor=weight * pref / dt,
                        )
                        rhs_value += weight * pref * float(previous_auxiliary[i, j, mu, nu]) / dt
                    elif alpha == 1:
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i + 1,
                            j=j,
                            mu=mu,
                            nu=nu,
                            prefactor=weight * pref / (2.0 * dx),
                        )
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i - 1,
                            j=j,
                            mu=mu,
                            nu=nu,
                            prefactor=-weight * pref / (2.0 * dx),
                        )
                    else:
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i,
                            j=j + 1,
                            mu=mu,
                            nu=nu,
                            prefactor=weight * pref / (2.0 * dz),
                        )
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i,
                            j=j - 1,
                            mu=mu,
                            nu=nu,
                            prefactor=-weight * pref / (2.0 * dz),
                        )
                    for lam in range(3):
                        conn = -pref * float(geom.gamma[i, j, lam, alpha, mu])
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i,
                            j=j,
                            mu=lam,
                            nu=nu,
                            prefactor=weight * conn,
                        )
                    for lam in range(3):
                        conn = -pref * float(geom.gamma[i, j, lam, alpha, nu])
                        _add_auxiliary_component_entries(
                            entries,
                            col=col,
                            basis=basis,
                            atoms=atoms,
                            i=i,
                            j=j,
                            mu=mu,
                            nu=lam,
                            prefactor=weight * conn,
                        )
            append_row(row_cols, row_vals, rhs, entries, rhs_value)
            append_row(force_cols, force_vals, force_rhs, entries, rhs_value)

    initial = _initial_auxiliary_coefficients(col=col, basis=basis, atoms=atoms, r_need=r_need)
    if solve_mode == "penalty":
        coeff, solve_info = solve_sparse_lsqr(row_cols, row_vals, np.asarray(rhs, dtype=float), ncols, lsqr_tol, lsqr_maxiter)
    elif solve_mode == "project":
        coeff, solve_info = _auxiliary_project_conservation(
            initial,
            force_cols,
            force_vals,
            force_rhs,
            ncols,
            tol=lsqr_tol,
            maxiter=lsqr_maxiter,
        )
    elif solve_mode == "alm":
        coeff, solve_info = _auxiliary_alm_conservation(
            initial,
            alg_cols,
            alg_vals,
            alg_rhs,
            force_cols,
            force_vals,
            force_rhs,
            ncols,
            mu0=alm_mu0,
            growth=alm_growth,
            outer=alm_outer,
            tol=lsqr_tol,
            maxiter=lsqr_maxiter,
        )
    else:
        raise ValueError(f"unknown global auxiliary solve mode: {solve_mode}")
    auxiliary = _coefficients_to_auxiliary(coeff, col=col, basis=basis, atoms=atoms, shape=mask.shape)
    alg_res = sparse_apply(alg_cols, alg_vals, coeff) - np.asarray(alg_rhs, dtype=float)
    force_res = sparse_apply(force_cols, force_vals, coeff) - np.asarray(force_rhs, dtype=float) if force_cols else np.zeros(0)
    initial_force = sparse_apply(force_cols, force_vals, initial) - np.asarray(force_rhs, dtype=float) if force_cols else np.zeros(0)
    info = {
        **solve_info,
        "enabled": True,
        "trace0": bool(trace0),
        "valid_points": int(np.count_nonzero(valid)),
        "conservation_points": int(np.count_nonzero(force_mask)),
        "unknowns": int(ncols),
        "rows_total": int(len(row_cols)),
        "rows_algebraic": int(len(alg_cols)),
        "rows_conservation": int(len(force_cols)),
        "solve_mode": str(solve_mode),
        "conservation_weight": float(conservation_weight),
        "alm_mu0": float(alm_mu0),
        "alm_growth": float(alm_growth),
        "alm_outer": int(alm_outer),
        "target_floor": float(target_floor),
        "force_scale": float(force_scale),
        "algebraic_relative_system_residual": float(
            np.linalg.norm(alg_res) / max(float(np.linalg.norm(np.asarray(alg_rhs, dtype=float))), 1.0e-300)
        ),
        "conservation_weighted_norm": float(np.linalg.norm(force_res)),
        "conservation_weighted_relative_to_initial": float(
            np.linalg.norm(force_res) / max(float(np.linalg.norm(initial_force)), 1.0e-300)
        )
        if force_res.size
        else 0.0,
        "initial_coeff_norm": float(np.linalg.norm(initial)),
        "coeff_norm": float(np.linalg.norm(coeff)),
        "coeff_change_over_initial": float(np.linalg.norm(coeff - initial) / max(float(np.linalg.norm(initial)), 1.0e-300)),
        "conservation_discretization": "backward-time, centered-space covariant divergence; causal runtime approximation",
    }
    return auxiliary, info


def linear_correct_metric_plus(
    *,
    metric_minus: np.ndarray,
    metric_center: np.ndarray,
    metric_plus_guess: np.ndarray,
    auxiliary: np.ndarray,
    r_need: np.ndarray,
    active: np.ndarray,
    rho_weight: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    ridge: float,
    eps: float,
    cg_tol: float,
    cg_maxiter: int,
    target_floor_frac: float,
) -> tuple[np.ndarray, dict[str, object]]:
    row_points = [(int(i), int(j)) for i, j in np.argwhere(active)]
    if not row_points:
        return metric_plus_guess, {"enabled": True, "fit_points": 0, "reason": "empty active mask"}
    target_norm_grid = tensor_norm(r_need)
    target_norm = np.asarray([max(float(target_norm_grid[i, j]), 1.0e-300) for i, j in row_points], dtype=float)
    finite_target = target_norm[np.isfinite(target_norm)]
    floor = float(target_floor_frac) * max(float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0, 1.0e-300)
    row_weights = np.asarray(
        [
            np.sqrt(max(float(rho_weight[i, j]), 0.0) / max(float(np.max(rho_weight)), 1.0e-300))
            / max(target_norm[n], floor)
            for n, (i, j) in enumerate(row_points)
        ],
        dtype=float,
    )
    residual = auxiliary - r_need
    rhs = np.concatenate([row_weights[n] * sym_to_vec(residual[i, j]) for n, (i, j) in enumerate(row_points)])
    active_metric = active.copy()
    active_metric[0, :] = False
    active_metric[-1, :] = False
    active_metric[:, 0] = False
    active_metric[:, -1] = False
    variables = metric_variable_index_for_slices(active_metric, ("plus",))
    dg_base, d2g_base = metric_jets_full(metric_minus, metric_center, metric_plus_guess, dt, dx, dz)
    col_rows, col_vals, col_norms = build_sparse_columns(
        metric_cov=metric_center,
        dg_base=dg_base,
        d2g_base=d2g_base,
        variables=variables,
        row_points=row_points,
        row_weights=row_weights,
        dt=dt,
        dx=dx,
        dz=dz,
        eps=eps,
    )
    entry_row, entry_col, entry_val, entry_n_rows = build_entry_arrays(col_rows, col_vals, col_norms, len(row_points))

    def b_apply(y: np.ndarray) -> np.ndarray:
        return apply_b_entries(entry_row, entry_col, entry_val, y, entry_n_rows)

    def bt_apply(v: np.ndarray) -> np.ndarray:
        return apply_bt_entries(entry_row, entry_col, entry_val, v, len(variables))

    bt_rhs = bt_apply(rhs)

    def normal_matvec(y: np.ndarray) -> np.ndarray:
        return bt_apply(b_apply(y)) + float(ridge) ** 2 * y

    y, solve_info = cg_solve(normal_matvec, bt_rhs, tol=float(cg_tol), maxiter=int(cg_maxiter))
    delta_vec = y / col_norms
    delta_metric = np.zeros_like(metric_plus_guess)
    for value, (_, i, j, a, b) in zip(delta_vec, variables):
        delta_metric[i, j, a, b] += value
        if a != b:
            delta_metric[i, j, b, a] += value
    corrected = metric_plus_guess + delta_metric
    corrected = 0.5 * (corrected + np.swapaxes(corrected, -1, -2))
    return corrected, {
        "enabled": True,
        "fit_points": int(np.count_nonzero(active)),
        "variables": int(len(variables)),
        "rows": int(rhs.size),
        "scalar_nonzeros": int(entry_val.size),
        "ridge": float(ridge),
        "linear_eps": float(eps),
        "solve": solve_info,
        "relative_delta_metric_p95": float(
            np.percentile(
                tensor_norm(delta_metric[active_metric]) / np.maximum(tensor_norm(metric_plus_guess[active_metric]), 1.0e-300),
                95.0,
            )
        )
        if np.any(active_metric)
        else 0.0,
    }


def render_summary(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    support: np.ndarray,
    rho_a: np.ndarray,
    rho_d_a: np.ndarray,
    residual: np.ndarray,
    records: list[dict[str, float]],
) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    vmax = max(
        float(np.percentile(rho_a[support], 99.5)) if np.any(support) else 0.0,
        float(np.percentile(rho_d_a[support], 99.5)) if np.any(support) else 0.0,
        1.0e-300,
    )
    diff = rho_d_a - rho_a
    diff_v = max(float(np.percentile(np.abs(diff[support]), 99.0)) if np.any(support) else 0.0, 1.0e-300)
    res_v = max(float(np.percentile(np.log10(1.0 + residual[support]), 99.0)) if np.any(support) else 1.0, 1.0)

    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.2), constrained_layout=True)
    panels = [
        (rho_a, "A rho at final time", "viridis", 0.0, vmax),
        (rho_d_a, "D pulled-back rho at final time", "viridis", 0.0, vmax),
        (diff, "D - A pulled-back rho", "coolwarm", -diff_v, diff_v),
        (np.log10(1.0 + np.maximum(residual, 0.0)), "log10(1+D residual)", "magma", 0.0, res_v),
    ]
    for ax, (field, title, cmap, vmin, vmax_panel) in zip(axes.flat[:4], panels):
        im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax_panel)
        ax.contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.55)
        ax.set_title(title + "; white=support")
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    times = np.asarray([r["tau_old"] for r in records], dtype=float)
    axes[1, 1].plot(times, [r["rho_pullback_rel_l1_support"] for r in records], marker="o", label="support")
    axes[1, 1].plot(times, [r["rho_pullback_rel_l1_core10"] for r in records], marker="o", label="core10")
    axes[1, 1].set_title("D-A rho relative L1")
    axes[1, 1].set_xlabel("tau_old")
    axes[1, 1].set_yscale("log")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend()
    axes[1, 2].plot(times, [r["d_equation_residual_wmean_core10"] for r in records], marker="o", label="D residual")
    axes[1, 2].plot(times, [max(r["disc_min_support"], 1.0e-300) for r in records], marker="o", label="disc min")
    axes[1, 2].set_title("equation / mass-shell diagnostics")
    axes[1, 2].set_xlabel("tau_old")
    axes[1, 2].set_yscale("log")
    axes[1, 2].grid(alpha=0.25)
    axes[1, 2].legend()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    package_path = args.package
    package = np.load(package_path)
    tau0 = float(args.tau if args.tau is not None else infer_tau_from_package(package_path))

    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    mass = float(ref["params"].m)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.dt_old) * old_scale

    x = np.asarray(package["x"], dtype=float)
    z = np.asarray(package["z"], dtype=float)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    support = np.asarray(package["support"], dtype=bool)
    core10 = np.asarray(package["core10"], dtype=bool)
    fit_mask = np.asarray(package["fit_mask"], dtype=bool)
    active = support.copy()
    if args.evolve_region == "core10":
        active = core10.copy()
    elif args.evolve_region == "fit_mask":
        active = fit_mask.copy()
    for _ in range(int(args.active_dilation)):
        active = active | np.roll(active, 1, 0) | np.roll(active, -1, 0) | np.roll(active, 1, 1) | np.roll(active, -1, 1)
    active &= support

    metric_minus = np.asarray(package["metric_minus"], dtype=float)
    metric_center = np.asarray(package["metric_center"], dtype=float)
    metric_plus = np.asarray(package["metric_plus"], dtype=float)
    auxiliary = np.asarray(package["auxiliary_cov"], dtype=float)
    n_cons = np.asarray(package["n_cons"], dtype=float)
    u_x = np.asarray(package["u_x"], dtype=float)
    u_z = np.asarray(package["u_z"], dtype=float)
    u_t_ref = np.asarray(package["u_t"], dtype=float)
    measure0 = np.asarray(package["measure_tilde"], dtype=float)
    sqrt_g0 = safe_sqrt_abs_det(metric_center)
    rho_tilde = np.where(sqrt_g0 > 0.0, measure0 / np.maximum(sqrt_g0, 1.0e-300), 0.0)
    rho_tilde = np.where(active, rho_tilde, 0.0)
    auxiliary_prev = auxiliary.copy()

    records: list[dict[str, float]] = []
    current = None
    stopped_reason = "completed"
    steps_completed = 0
    for step in range(int(args.steps) + 1):
        steps_completed = step
        tau = tau0 + step * float(args.dt_old)
        metric_inv_center = np.linalg.pinv(metric_center.reshape(-1, 3, 3), rcond=1.0e-12, hermitian=True).reshape(metric_center.shape)
        det_center = np.linalg.det(metric_center.reshape(-1, 3, 3)).reshape(metric_center.shape[:2])
        current = None
        from coordinate_matter_evolution import coordinate_matter_rhs_covector

        current = coordinate_matter_rhs_covector(
            n_cons=np.where(active, n_cons, 0.0),
            u_x=u_x,
            u_z=u_z,
            metric_cov_txz=metric_center,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=u_t_ref,
            metric_inv_txz=metric_inv_center,
            det_cov_txz=det_center,
        )
        rho_tilde = np.where(active, current["rho"], 0.0)
        measure_tilde = np.where(active, current_measure(metric_center, rho_tilde), 0.0)
        r_need, r_scalar, geom_inv, tilde_x = recompute_r_need(
            metric_minus,
            metric_center,
            metric_plus,
            rho_tilde,
            current["u_t"],
            u_x,
            u_z,
            mass=mass,
            mp=float(args.mp),
            dt=dt,
            dx=dx,
            dz=dz,
        )
        if bool(args.trace_project_auxiliary):
            auxiliary = trace_project_auxiliary(auxiliary, geom_inv, metric_center)
        aux_info: dict[str, object] = {"enabled": False}
        if args.auxiliary_update == "local":
            auxiliary = local_auxiliary_solve(
                metric_cov=metric_center,
                metric_inv=geom_inv,
                r_need=r_need,
                rho_tilde=rho_tilde,
                u_t=current["u_t"],
                u_x=u_x,
                u_z=u_z,
                dx=dx,
                dz=dz,
                mask=active,
                trace0=bool(args.auxiliary_trace0),
            )
            aux_info = {"enabled": True, "mode": "local"}
        elif args.auxiliary_update == "global":
            auxiliary, aux_info = global_auxiliary_solve(
                metric_minus=metric_minus,
                metric_center=metric_center,
                metric_plus=metric_plus,
                metric_inv=geom_inv,
                r_need=r_need,
                rho_tilde=rho_tilde,
                u_t=current["u_t"],
                u_x=u_x,
                u_z=u_z,
                previous_auxiliary=auxiliary_prev,
                dt=dt,
                dx=dx,
                dz=dz,
                mask=active,
                trace0=bool(args.auxiliary_trace0),
                solve_mode=str(args.auxiliary_global_mode),
                target_floor_frac=float(args.target_floor_frac),
                conservation_weight=float(args.global_conservation_weight),
                lsqr_tol=float(args.aux_lsqr_tol),
                lsqr_maxiter=int(args.aux_lsqr_maxiter),
                alm_mu0=float(args.aux_alm_mu0),
                alm_growth=float(args.aux_alm_growth),
                alm_outer=int(args.aux_alm_outer),
            )
        d_residual = tensor_norm(auxiliary - r_need) / np.maximum(tensor_norm(r_need), 1.0e-300)
        case_a = build_case(args, ref, tau)
        rho_pull = pullback_rho_to_a(measure_tilde, current["u_t"], u_x, u_z, mass, float(args.x_floor), active)
        weight = np.maximum(case_a.rho_a, 0.0)
        wsum_core = np.maximum(weight[core10], 0.0)
        stop_mask = active
        guard = discriminant_guard(
            current["discriminant"],
            stop_mask,
            tolerance=float(args.disc_tolerance),
            fraction_tolerance=float(args.disc_fraction_tolerance),
        )
        support_guard = discriminant_guard(
            current["discriminant"],
            support,
            tolerance=float(args.disc_tolerance),
            fraction_tolerance=float(args.disc_fraction_tolerance),
        )
        rho_guard = lower_bound_guard(
            rho_tilde,
            stop_mask,
            tolerance=float(args.rho_tilde_tolerance),
            fraction_tolerance=float(args.rho_tilde_fraction_tolerance),
        )
        rho_support_guard = lower_bound_guard(
            rho_tilde,
            support,
            tolerance=float(args.rho_tilde_tolerance),
            fraction_tolerance=float(args.rho_tilde_fraction_tolerance),
        )
        disc_min = float(guard["min"])
        rec = {
            "step": int(step),
            "tau_old": float(tau),
            "rho_pullback_rel_l1_support": relative_l1(case_a.rho_a, rho_pull, support),
            "rho_pullback_weighted_l1_support": weighted_relative_l1(case_a.rho_a, rho_pull, weight, support),
            "rho_pullback_rel_l1_core10": relative_l1(case_a.rho_a, rho_pull, core10),
            "rho_pullback_weighted_l1_core10": weighted_relative_l1(case_a.rho_a, rho_pull, weight, core10),
            "d_equation_residual_wmean_core10": float(
                np.sum(d_residual[core10] * wsum_core) / max(float(np.sum(wsum_core)), 1.0e-300)
            )
            if np.any(core10)
            else 0.0,
            "d_equation_residual_p95_core10": float(np.percentile(d_residual[core10], 95.0)) if np.any(core10) else 0.0,
            "disc_min_active": disc_min,
            "disc_min_support": float(support_guard["min"]),
            "negative_discriminant_fraction_active": float(guard["negative_fraction"]),
            "large_negative_discriminant_fraction_active": float(guard["large_negative_fraction"]),
            "negative_discriminant_fraction_support": float(support_guard["negative_fraction"]),
            "large_negative_discriminant_fraction_support": float(support_guard["large_negative_fraction"]),
            "rho_tilde_min_active": float(np.nanmin(rho_tilde[active])) if np.any(active) else 0.0,
            "rho_tilde_max_active": float(np.nanmax(rho_tilde[active])) if np.any(active) else 0.0,
            "rho_tilde_negative_fraction_active": float(rho_guard["negative_fraction"]),
            "rho_tilde_large_negative_fraction_active": float(rho_guard["large_negative_fraction"]),
            "rho_tilde_negative_fraction_support": float(rho_support_guard["negative_fraction"]),
            "rho_tilde_large_negative_fraction_support": float(rho_support_guard["large_negative_fraction"]),
            "raw_y_abs_p95_core10": float(np.percentile(np.abs(r_scalar[core10]), 95.0)) if np.any(core10) else 0.0,
        }
        if aux_info.get("enabled", False):
            for key in (
                "valid_points",
                "conservation_points",
                "unknowns",
                "rows_total",
                "rows_algebraic",
                "rows_conservation",
                "iterations",
                "system_relative_residual",
                "algebraic_relative_system_residual",
                "conservation_weighted_norm",
                "conservation_weighted_relative_to_initial",
                "coeff_change_over_initial",
                "projected_force_relative_norm",
                "algebraic_relative_residual_norm",
                "delta_over_initial_norm",
            ):
                if key in aux_info:
                    rec[f"auxiliary_{key}"] = float(aux_info[key])  # type: ignore[arg-type]
        records.append(rec)
        if step == int(args.steps):
            break
        if args.stop_on_negative_discriminant and bool(guard["stop"]):
            stopped_reason = (
                f"negative mass-shell discriminant: min={disc_min:.6e}, "
                f"large_fraction={float(guard['large_negative_fraction']):.6e}, "
                f"tolerance={float(args.disc_tolerance):.6e}"
            )
            break
        if args.stop_on_negative_rho_tilde and bool(rho_guard["stop"]):
            stopped_reason = (
                f"negative rho_tilde: min={float(rho_guard['min']):.6e}, "
                f"large_fraction={float(rho_guard['large_negative_fraction']):.6e}, "
                f"tolerance={float(args.rho_tilde_tolerance):.6e}"
            )
            break

        auxiliary_current = auxiliary.copy()
        try:
            n_next, ux_next, uz_next, current_next = rk4_matter_step(
                n_cons=np.where(active, n_cons, 0.0),
                u_x=u_x,
                u_z=u_z,
                metric_cov=metric_center,
                mass=mass,
                dx=dx,
                dz=dz,
                dt=dt,
                u_t_reference=current["u_t"],
                active_mask=active,
            )
        except FloatingPointError as exc:
            stopped_reason = f"matter step failed: {exc}"
            break

        # Advance the three-level metric history correctly:
        #   old: (g_{n-1}, g_n, g_{n+1})
        #   new center is the already available g_{n+1};
        #   now solve/extrapolate a fresh g_{n+2}.  The first version of this
        #   prototype accidentally reused g_{n+1} as both center and future,
        #   which artificially injected a huge second time derivative.
        new_minus = metric_center
        new_center = metric_plus
        metric_next = 2.0 * new_center - new_minus
        metric_next = np.where(active[..., None, None], metric_next, new_center)
        metric_next = 0.5 * (metric_next + np.swapaxes(metric_next, -1, -2))

        next_metric_inv = np.linalg.pinv(new_center.reshape(-1, 3, 3), rcond=1.0e-12, hermitian=True).reshape(new_center.shape)
        next_det = np.linalg.det(new_center.reshape(-1, 3, 3)).reshape(new_center.shape[:2])
        from coordinate_matter_evolution import coordinate_matter_rhs_covector

        next_current_on_center = coordinate_matter_rhs_covector(
            n_cons=np.where(active, n_next, 0.0),
            u_x=ux_next,
            u_z=uz_next,
            metric_cov_txz=new_center,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=current_next["u_t"],
            metric_inv_txz=next_metric_inv,
            det_cov_txz=next_det,
        )
        next_rho = np.where(active, next_current_on_center["rho"], 0.0)
        metric_next_work = metric_next
        corrector_info = {"enabled": False}
        aux_info = {"enabled": False}
        joint_outer = max(1, int(args.joint_outer_iterations))
        metric_relax = float(args.joint_metric_relaxation)
        residual_gate_rel = float(args.joint_residual_gate_rel)
        last_accepted_metric = metric_next_work.copy()
        last_d_residual_mean: float | None = None
        for _joint_it in range(joint_outer):
            r_need, r_scalar, geom_inv, tilde_x = recompute_r_need(
                new_minus,
                new_center,
                metric_next_work,
                next_rho,
                next_current_on_center["u_t"],
                ux_next,
                uz_next,
                mass=mass,
                mp=float(args.mp),
                dt=dt,
                dx=dx,
                dz=dz,
            )
            if bool(args.trace_project_auxiliary):
                auxiliary = trace_project_auxiliary(auxiliary, geom_inv, new_center)
            if args.auxiliary_update == "global":
                auxiliary, aux_info = global_auxiliary_solve(
                    metric_minus=new_minus,
                    metric_center=new_center,
                    metric_plus=metric_next_work,
                    metric_inv=geom_inv,
                    r_need=r_need,
                    rho_tilde=next_rho,
                    u_t=next_current_on_center["u_t"],
                    u_x=ux_next,
                    u_z=uz_next,
                    previous_auxiliary=auxiliary_current,
                    dt=dt,
                    dx=dx,
                    dz=dz,
                    mask=active,
                    trace0=bool(args.auxiliary_trace0),
                    solve_mode=str(args.auxiliary_global_mode),
                    target_floor_frac=float(args.target_floor_frac),
                    conservation_weight=float(args.global_conservation_weight),
                    lsqr_tol=float(args.aux_lsqr_tol),
                    lsqr_maxiter=int(args.aux_lsqr_maxiter),
                    alm_mu0=float(args.aux_alm_mu0),
                    alm_growth=float(args.aux_alm_growth),
                    alm_outer=int(args.aux_alm_outer),
                )
            elif args.auxiliary_update == "local":
                auxiliary = local_auxiliary_solve(
                    metric_cov=new_center,
                    metric_inv=geom_inv,
                    r_need=r_need,
                    rho_tilde=next_rho,
                    u_t=next_current_on_center["u_t"],
                    u_x=ux_next,
                    u_z=uz_next,
                    dx=dx,
                    dz=dz,
                    mask=active,
                    trace0=bool(args.auxiliary_trace0),
                )
                aux_info = {"enabled": True, "mode": "local"}
            d_residual = tensor_norm(auxiliary - r_need) / np.maximum(tensor_norm(r_need), 1.0e-300)
            d_residual_mean = float(
                np.sum(d_residual[core10] * wsum_core) / max(float(np.sum(wsum_core)), 1.0e-300)
            )
            if last_d_residual_mean is not None and d_residual_mean > last_d_residual_mean * (1.0 + residual_gate_rel):
                metric_next_work = last_accepted_metric
                break
            next_r_need, _, _, _ = recompute_r_need(
                new_minus,
                new_center,
                metric_next_work,
                next_rho,
                next_current_on_center["u_t"],
                ux_next,
                uz_next,
                mass=mass,
                mp=float(args.mp),
                dt=dt,
                dx=dx,
                dz=dz,
            )
            if args.metric_corrector == "linear-plus":
                metric_next_candidate, corrector_info = linear_correct_metric_plus(
                    metric_minus=new_minus,
                    metric_center=new_center,
                    metric_plus_guess=metric_next_work,
                    auxiliary=auxiliary,
                    r_need=next_r_need,
                    active=active,
                    rho_weight=case_a.rho_a,
                    dt=dt,
                    dx=dx,
                    dz=dz,
                    ridge=float(args.corrector_ridge),
                    eps=float(args.corrector_eps),
                    cg_tol=float(args.corrector_cg_tol),
                    cg_maxiter=int(args.corrector_cg_maxiter),
                    target_floor_frac=float(args.target_floor_frac),
                )
                metric_next_candidate = (1.0 - metric_relax) * metric_next_work + metric_relax * metric_next_candidate
                metric_next_candidate = 0.5 * (metric_next_candidate + np.swapaxes(metric_next_candidate, -1, -2))
            else:
                metric_next_candidate = metric_next_work
            delta_metric_rel = float(
                np.percentile(
                    tensor_norm(metric_next_candidate - metric_next_work)
                    / np.maximum(tensor_norm(metric_next_work), 1.0e-300),
                    95.0,
                )
            ) if np.any(active) else 0.0
            last_d_residual_mean = d_residual_mean
            last_accepted_metric = metric_next_work.copy()
            metric_next_work = metric_next_candidate
            if delta_metric_rel < 1.0e-8:
                break

        metric_next = metric_next_work
        if records:
            records[-1]["metric_corrector_enabled"] = float(bool(corrector_info.get("enabled", False)))
            records[-1]["metric_corrector_delta_p95"] = float(corrector_info.get("relative_delta_metric_p95", 0.0))
            records[-1]["auxiliary_solve_mode"] = str(aux_info.get("mode", aux_info.get("solve_mode", "")))
            records[-1]["auxiliary_joint_outer"] = float(joint_outer)
            records[-1]["auxiliary_joint_metric_relaxation"] = float(metric_relax)
            records[-1]["d_equation_residual_wmean_core10"] = float(
                np.sum(d_residual[core10] * wsum_core) / max(float(np.sum(wsum_core)), 1.0e-300)
            )
            records[-1]["d_equation_residual_p95_core10"] = float(np.percentile(d_residual[core10], 95.0)) if np.any(core10) else 0.0
            records[-1]["raw_y_abs_p95_core10"] = float(np.percentile(np.abs(r_scalar[core10]), 95.0)) if np.any(core10) else 0.0
            for key in (
                "valid_points",
                "conservation_points",
                "unknowns",
                "rows_total",
                "rows_algebraic",
                "rows_conservation",
                "iterations",
                "system_relative_residual",
                "algebraic_relative_system_residual",
                "conservation_weighted_norm",
                "conservation_weighted_relative_to_initial",
                "coeff_change_over_initial",
                "projected_force_relative_norm",
                "algebraic_relative_residual_norm",
                "delta_over_initial_norm",
            ):
                if key in aux_info:
                    records[-1][f"auxiliary_{key}"] = float(aux_info[key])  # type: ignore[arg-type]

        auxiliary_prev = auxiliary_current.copy()
        metric_minus, metric_center, metric_plus = new_minus, new_center, metric_next
        n_cons = np.where(active, n_next, 0.0)
        u_x = np.where(active, ux_next, 0.0)
        u_z = np.where(active, uz_next, 0.0)
        u_t_ref = current_next["u_t"]

    final_tau = records[-1]["tau_old"]
    final_a = build_case(args, ref, final_tau)
    final_rho_pull = pullback_rho_to_a(
        current_measure(metric_center, rho_tilde),
        current["u_t"],
        u_x,
        u_z,
        mass,
        float(args.x_floor),
        active,
    )
    final_residual = d_residual
    plot_path = args.output / "d_harmonic_standalone_vs_a.png"
    render_summary(plot_path, x, z, support, final_a.rho_a, final_rho_pull, final_residual, records)
    fields_path = args.output / "fields_final.npz"
    np.savez_compressed(
        fields_path,
        x=x,
        z=z,
        support=support,
        core10=core10,
        active=active,
        rho_A=final_a.rho_a,
        rho_D_to_A=final_rho_pull,
        rho_tilde=rho_tilde,
        measure_tilde=current_measure(metric_center, rho_tilde),
        n_cons=n_cons,
        discriminant=current["discriminant"],
        u_t=current["u_t"],
        u_x=u_x,
        u_z=u_z,
        metric_minus=metric_minus,
        metric_center=metric_center,
        metric_plus=metric_plus,
        auxiliary_cov=auxiliary,
        d_equation_residual=final_residual,
    )
    report = {
        "parameters": {
            "package": str(package_path.resolve()),
            "tau_initial": float(tau0),
            "steps_requested": int(args.steps),
            "steps_completed": int(steps_completed),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(dt),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "evolve_region": str(args.evolve_region),
            "active_dilation": int(args.active_dilation),
            "mp_ev": float(args.mp),
            "mass_ev": float(mass),
            "stopped_reason": stopped_reason,
            "recompute_first_plus": bool(args.recompute_first_plus),
            "trace_project_auxiliary": bool(args.trace_project_auxiliary),
            "metric_corrector": str(args.metric_corrector),
            "auxiliary_update": str(args.auxiliary_update),
            "auxiliary_trace0": bool(args.auxiliary_trace0),
            "auxiliary_global_mode": str(args.auxiliary_global_mode),
            "global_conservation_weight": float(args.global_conservation_weight),
            "aux_lsqr_tol": float(args.aux_lsqr_tol),
            "aux_lsqr_maxiter": int(args.aux_lsqr_maxiter),
            "aux_alm_mu0": float(args.aux_alm_mu0),
            "aux_alm_growth": float(args.aux_alm_growth),
            "aux_alm_outer": int(args.aux_alm_outer),
            "disc_tolerance": float(args.disc_tolerance),
            "disc_fraction_tolerance": float(args.disc_fraction_tolerance),
            "rho_tilde_tolerance": float(args.rho_tilde_tolerance),
            "rho_tilde_fraction_tolerance": float(args.rho_tilde_fraction_tolerance),
            "stop_on_negative_rho_tilde": bool(args.stop_on_negative_rho_tilde),
        },
        "scope": {
            "what_this_is": "Minimal standalone D-branch generalized-harmonic principal prototype.  It advances D matter variables and metric history from a D initial package; A/KG is used only as a posterior comparison.",
            "what_is_not_finished": "The metric update is an explicit principal-part wave update with residual feedback, not yet the full implicit nonlinear generalized-harmonic Einstein solve with simultaneous C-sector conservation.",
        },
        "final": records[-1],
        "records": records,
        "diagnostics": {
            "final_residual_core10": stats(final_residual[core10]),
            "final_rho_pullback_point_rel_core10": stats(
                np.abs(final_rho_pull[core10] - final_a.rho_a[core10]) / np.maximum(np.abs(final_a.rho_a[core10]), 1.0e-300),
                abs_value=False,
            ),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "fields_npz": str(fields_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=None)
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--evolve-region", choices=["support", "core10", "fit_mask"], default="core10")
    parser.add_argument("--active-dilation", type=int, default=0)
    parser.add_argument("--stop-on-negative-discriminant", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--recompute-first-plus", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--trace-project-auxiliary", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--auxiliary-update", choices=["frozen", "local", "global"], default="frozen")
    parser.add_argument("--auxiliary-trace0", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--auxiliary-global-mode", choices=["penalty", "project", "alm"], default="penalty")
    parser.add_argument("--global-conservation-weight", type=float, default=1.0)
    parser.add_argument("--aux-lsqr-tol", type=float, default=1.0e-6)
    parser.add_argument("--aux-lsqr-maxiter", type=int, default=800)
    parser.add_argument("--aux-alm-mu0", type=float, default=1.0)
    parser.add_argument("--aux-alm-growth", type=float, default=10.0)
    parser.add_argument("--aux-alm-outer", type=int, default=3)
    parser.add_argument("--metric-corrector", choices=["none", "linear-plus"], default="none")
    parser.add_argument("--corrector-ridge", type=float, default=1.0e-6)
    parser.add_argument("--corrector-eps", type=float, default=1.0e-8)
    parser.add_argument("--corrector-cg-tol", type=float, default=1.0e-7)
    parser.add_argument("--corrector-cg-maxiter", type=int, default=600)
    parser.add_argument("--joint-outer-iterations", type=int, default=2)
    parser.add_argument("--joint-metric-relaxation", type=float, default=1.0)
    parser.add_argument("--joint-residual-gate-rel", type=float, default=0.0)
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--disc-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--disc-fraction-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--stop-on-negative-rho-tilde", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-tilde-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--rho-tilde-fraction-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
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
