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
from diagnose_gbcd_trace_local_closure import solve_local, trace_coefficients


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
        d_residual = tensor_norm(auxiliary - r_need) / np.maximum(tensor_norm(r_need), 1.0e-300)
        case_a = build_case(args, ref, tau)
        rho_pull = pullback_rho_to_a(measure_tilde, current["u_t"], u_x, u_z, mass, float(args.x_floor), active)
        weight = np.maximum(case_a.rho_a, 0.0)
        wsum_core = np.maximum(weight[core10], 0.0)
        stop_mask = active
        disc_min = float(np.nanmin(current["discriminant"][stop_mask])) if np.any(stop_mask) else 0.0
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
            "disc_min_support": float(np.nanmin(current["discriminant"][support])) if np.any(support) else 0.0,
            "negative_discriminant_fraction_active": float(
                np.count_nonzero(current["discriminant"][stop_mask] < 0.0) / max(int(np.count_nonzero(stop_mask)), 1)
            ),
            "negative_discriminant_fraction_support": float(
                np.count_nonzero(current["discriminant"][support] < 0.0) / max(int(np.count_nonzero(support)), 1)
            ),
            "rho_tilde_min_active": float(np.nanmin(rho_tilde[active])) if np.any(active) else 0.0,
            "rho_tilde_max_active": float(np.nanmax(rho_tilde[active])) if np.any(active) else 0.0,
            "raw_y_abs_p95_core10": float(np.percentile(np.abs(r_scalar[core10]), 95.0)) if np.any(core10) else 0.0,
        }
        records.append(rec)
        if step == int(args.steps):
            break
        if args.stop_on_negative_discriminant and disc_min < -abs(float(args.disc_tolerance)):
            stopped_reason = f"negative mass-shell discriminant: {disc_min:.6e}"
            break

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
        next_r_need, _, _, _ = recompute_r_need(
            new_minus,
            new_center,
            metric_next,
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
        corrector_info = {"enabled": False}
        if args.metric_corrector == "linear-plus":
            metric_next, corrector_info = linear_correct_metric_plus(
                metric_minus=new_minus,
                metric_center=new_center,
                metric_plus_guess=metric_next,
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
        if records:
            records[-1]["metric_corrector_enabled"] = float(bool(corrector_info.get("enabled", False)))
            records[-1]["metric_corrector_delta_p95"] = float(corrector_info.get("relative_delta_metric_p95", 0.0))

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
    parser.add_argument("--auxiliary-update", choices=["frozen", "local"], default="frozen")
    parser.add_argument("--auxiliary-trace0", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--metric-corrector", choices=["none", "linear-plus"], default="none")
    parser.add_argument("--corrector-ridge", type=float, default=1.0e-6)
    parser.add_argument("--corrector-eps", type=float, default=1.0e-8)
    parser.add_argument("--corrector-cg-tol", type=float, default=1.0e-7)
    parser.add_argument("--corrector-cg-maxiter", type=int, default=600)
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--disc-tolerance", type=float, default=1.0e-12)
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
