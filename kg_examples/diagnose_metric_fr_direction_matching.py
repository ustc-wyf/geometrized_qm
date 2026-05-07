from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_bcd_residuals_from_a_reference import (
    covariant_hessian_from_jets,
    geometry_data,
    stress_tensor_tilde,
)
from analyze_c_terms_from_a_reference import metric_jets_full, scalar_jets_full
from diagnose_d_full_fr_equation_local_window import branch_d_f_phi
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import (
    build_physical_reference,
    make_cropped_snapshot,
    real_array,
)
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det


def stats_abs(values: np.ndarray | list[float]) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = np.abs(vals[np.isfinite(vals)])
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }


def stats_signed(values: np.ndarray | list[float]) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {
            "count": 0,
            "min": 0.0,
            "p05": 0.0,
            "p50": 0.0,
            "p95": 0.0,
            "max": 0.0,
            "abs_p50": 0.0,
            "abs_p95": 0.0,
        }
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p05": float(np.percentile(vals, 5.0)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
        "abs_p50": float(np.percentile(np.abs(vals), 50.0)),
        "abs_p95": float(np.percentile(np.abs(vals), 95.0)),
    }


def tensor_norm(tensor: np.ndarray) -> np.ndarray:
    return np.sqrt(np.einsum("...ab,...ab->...", tensor, tensor, optimize=True))


def tensor_dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...ab->...", a, b, optimize=True)


def cosine_to_target(tensor: np.ndarray, target: np.ndarray) -> np.ndarray:
    denom = np.maximum(tensor_norm(tensor) * tensor_norm(target), 1.0e-300)
    return tensor_dot(tensor, target) / denom


def make_snapshot_pack(args: argparse.Namespace, time_ev_inv: float, ix: np.ndarray, iz: np.ndarray, ref: dict[str, object]) -> dict[str, object]:
    return make_cropped_snapshot(
        psi0=np.asarray(ref["psi0"]),
        psi0_hat=np.asarray(ref["psi0_hat"]),
        omega=np.asarray(ref["omega"]),
        x_full=np.asarray(ref["x_full"]),
        z_full=np.asarray(ref["z_full"]),
        ix=ix,
        iz=iz,
        t=time_ev_inv,
        rho_floor=float(args.rho_floor),
        x_floor=float(args.x_floor),
        pinv_rcond=float(args.pinv_rcond),
        with_metric=True,
    )


def fit_local_basis(
    basis: np.ndarray,
    target: np.ndarray,
    mask: np.ndarray,
    rcond: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Least-squares target ~= sum_k coeff_k basis_k in fixed txz components."""
    nx, nz = target.shape[:2]
    k = basis.shape[2]
    coeff = np.full((nx, nz, k), np.nan, dtype=float)
    residual_rel = np.full((nx, nz), np.nan, dtype=float)
    fit_norm = np.full((nx, nz), np.nan, dtype=float)
    target_norm = tensor_norm(target)
    flat_mask = mask.reshape(-1)
    basis_flat = basis.reshape(-1, k, 9)
    target_flat = target.reshape(-1, 9)
    coeff_flat = coeff.reshape(-1, k)
    residual_flat = residual_rel.reshape(-1)
    fit_norm_flat = fit_norm.reshape(-1)
    target_norm_flat = target_norm.reshape(-1)

    for idx in np.where(flat_mask)[0]:
        mat = basis_flat[idx].T
        rhs = target_flat[idx]
        if not np.all(np.isfinite(mat)) or not np.all(np.isfinite(rhs)):
            continue
        if np.linalg.norm(mat) <= 1.0e-300:
            continue
        sol, *_ = np.linalg.lstsq(mat, rhs, rcond=rcond)
        fitted = mat @ sol
        coeff_flat[idx] = sol
        fit_norm_flat[idx] = float(np.linalg.norm(fitted))
        residual_flat[idx] = float(np.linalg.norm(fitted - rhs) / max(target_norm_flat[idx], 1.0e-300))

    return coeff, residual_rel, fit_norm


def single_value_diagnostic(x: np.ndarray, y: np.ndarray, mask: np.ndarray, bins: int = 50, min_count: int = 20) -> dict[str, float]:
    xv = np.asarray(x[mask], dtype=float)
    yv = np.asarray(y[mask], dtype=float)
    good = np.isfinite(xv) & np.isfinite(yv)
    xv = xv[good]
    yv = yv[good]
    if xv.size < min_count:
        return {"usable_bins": 0, "median_relative_mad": np.nan, "p95_relative_mad": np.nan, "global_abs_p95": np.nan}
    sx = np.sign(xv) * np.log10(1.0 + np.abs(xv))
    edges = np.linspace(float(np.min(sx)), float(np.max(sx)), int(bins) + 1)
    global_scale = max(float(np.percentile(np.abs(yv), 95.0)), 1.0e-300)
    rel_mads: list[float] = []
    global_mads: list[float] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        bin_mask = (sx >= lo) & (sx < hi)
        vals = yv[bin_mask]
        if vals.size < min_count:
            continue
        med = float(np.median(vals))
        mad = float(np.median(np.abs(vals - med)))
        rel_mads.append(mad / max(abs(med), 1.0e-12 * global_scale, 1.0e-300))
        global_mads.append(mad / max(global_scale, 1.0e-300))
    if not rel_mads:
        return {
            "usable_bins": 0,
            "median_relative_mad": np.nan,
            "p95_relative_mad": np.nan,
            "median_global_mad": np.nan,
            "p95_global_mad": np.nan,
            "global_abs_p95": global_scale,
        }
    return {
        "usable_bins": int(len(rel_mads)),
        "median_relative_mad": float(np.median(rel_mads)),
        "p95_relative_mad": float(np.percentile(rel_mads, 95.0)),
        "median_global_mad": float(np.median(global_mads)),
        "p95_global_mad": float(np.percentile(global_mads, 95.0)),
        "global_abs_p95": global_scale,
    }


def invariant_contract2(metric_inv: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.einsum("...ac,...bd,...ab,...cd->...", metric_inv, metric_inv, a, b, optimize=True)


def invariant_trace(metric_inv: np.ndarray, tensor: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...ab->...", metric_inv, tensor, optimize=True)


def invariant_algebraic_coefficients(
    metric_inv: np.ndarray,
    ricci: np.ndarray,
    metric_cov: np.ndarray,
    target: np.ndarray,
    mask: np.ndarray,
    rcond_floor: float = 1.0e-14,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve scalar contractions of phi*R_ab - 0.5*f*g_ab = target_ab.

    The two scalar equations are the trace and Ricci contraction.  This gives a
    coordinate-independent check of the local coefficients before the full
    tensor residual is inspected.
    """
    n_dim = ricci.shape[-1]
    r_scalar = invariant_trace(metric_inv, ricci)
    t_trace = invariant_trace(metric_inv, target)
    ricci2 = invariant_contract2(metric_inv, ricci, ricci)
    ricci_target = invariant_contract2(metric_inv, ricci, target)
    det = -0.5 * r_scalar * r_scalar + 0.5 * n_dim * ricci2

    coeff = np.full(ricci.shape[:2] + (2,), np.nan, dtype=float)
    residual = np.full(ricci.shape[:2], np.nan, dtype=float)
    rhs_norm = tensor_norm(target)

    good = mask & np.isfinite(det) & (np.abs(det) > rcond_floor * np.maximum(1.0, np.abs(r_scalar * r_scalar) + np.abs(ricci2)))
    phi = np.full_like(r_scalar, np.nan, dtype=float)
    f = np.full_like(r_scalar, np.nan, dtype=float)
    phi[good] = ((-0.5 * r_scalar[good]) * t_trace[good] + (0.5 * n_dim) * ricci_target[good]) / det[good]
    f[good] = ((-ricci2[good]) * t_trace[good] + r_scalar[good] * ricci_target[good]) / det[good]
    coeff[..., 0] = phi
    coeff[..., 1] = f
    fitted = phi[..., None, None] * ricci - 0.5 * f[..., None, None] * metric_cov
    residual[good] = tensor_norm(fitted - target)[good] / np.maximum(rhs_norm[good], 1.0e-300)
    return coeff, residual


def render_maps(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    support: np.ndarray,
    ricci_resid: np.ndarray,
    alg_resid: np.ndarray,
    full_resid: np.ndarray,
    d_full_resid: np.ndarray,
) -> None:
    extent = [
        float(x[0] * HBAR_C_EV_M * 1.0e6),
        float(x[-1] * HBAR_C_EV_M * 1.0e6),
        float(z[0] * HBAR_C_EV_M * 1.0e6),
        float(z[-1] * HBAR_C_EV_M * 1.0e6),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.0), constrained_layout=True)
    im0 = axes[0, 0].imshow(rho.T, origin="lower", extent=extent, aspect="equal", cmap="magma")
    axes[0, 0].contour(x * HBAR_C_EV_M * 1.0e6, z * HBAR_C_EV_M * 1.0e6, support.T.astype(float), levels=[0.5], colors=["white"], linewidths=0.6)
    axes[0, 0].set_title("rho; white=support")
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046)

    panels = [
        (ricci_resid, "best scalar * Ricci residual"),
        (alg_resid, "best local [phi Ricci - 1/2 f g] residual"),
        (full_resid, "best local full f(R) jet residual"),
        (d_full_resid, "current D full lhs residual"),
    ]
    for ax, (field, title) in zip(axes.ravel()[1:], panels):
        im = ax.imshow(np.log10(1.0e-16 + np.clip(field, 0.0, 1.0e16)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        ax.set_title("log10 " + title)
        fig.colorbar(im, ax=ax, fraction=0.046)
    axes[1, 2].axis("off")
    for ax in axes.ravel():
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_scatter(
    out_path: Path,
    r_scalar: np.ndarray,
    coeff_alg: np.ndarray,
    coeff_full: np.ndarray,
    mask: np.ndarray,
    sample_count: int,
) -> None:
    r = np.asarray(r_scalar[mask], dtype=float)
    ca = np.asarray(coeff_alg[mask], dtype=float)
    cf = np.asarray(coeff_full[mask], dtype=float)
    good = np.isfinite(r) & np.all(np.isfinite(ca), axis=1) & np.all(np.isfinite(cf), axis=1)
    r = r[good]
    ca = ca[good]
    cf = cf[good]
    if r.size > sample_count:
        rng = np.random.default_rng(12345)
        take = rng.choice(r.size, size=sample_count, replace=False)
        r = r[take]
        ca = ca[take]
        cf = cf[take]

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.2), constrained_layout=True)
    panels = [
        (ca[:, 0] if ca.size else np.array([]), "alg fit phi=f_R"),
        (ca[:, 1] if ca.size else np.array([]), "alg fit f"),
        (cf[:, 0] if cf.size else np.array([]), "full fit phi"),
        (cf[:, 1] if cf.size else np.array([]), "full fit f"),
        (cf[:, 2] if cf.size else np.array([]), "full fit phi_R"),
        (cf[:, 3] if cf.size else np.array([]), "full fit phi_RR"),
    ]
    for ax, (y, title) in zip(axes.ravel(), panels):
        if r.size and y.size:
            ax.scatter(r, y, s=2, alpha=0.18, linewidths=0)
        ax.set_xscale("symlog", linthresh=1.0e-30)
        ax.set_yscale("symlog", linthresh=1.0e-80)
        ax.set_xlabel("R_tilde")
        ax.set_ylabel(title)
        ax.set_title(title + " vs R_tilde")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def summarize_region(
    name: str,
    mask: np.ndarray,
    r_scalar: np.ndarray,
    rhs_norm: np.ndarray,
    ricci_cos: np.ndarray,
    d_alg_cos: np.ndarray,
    d_full_cos: np.ndarray,
    ricci_resid: np.ndarray,
    alg_resid: np.ndarray,
    full_resid: np.ndarray,
    d_alg_resid: np.ndarray,
    d_full_resid: np.ndarray,
    coeff_ricci: np.ndarray,
    coeff_alg: np.ndarray,
    coeff_full: np.ndarray,
    coeff_alg_inv: np.ndarray,
    resid_alg_inv: np.ndarray,
) -> dict[str, object]:
    return {
        "name": name,
        "count": int(np.count_nonzero(mask)),
        "R_tilde": stats_signed(r_scalar[mask]),
        "rhs_norm": stats_abs(rhs_norm[mask]),
        "cos_Ricci_to_T": stats_signed(ricci_cos[mask]),
        "cos_current_D_algebraic_to_T": stats_signed(d_alg_cos[mask]),
        "cos_current_D_full_lhs_to_T": stats_signed(d_full_cos[mask]),
        "residual_best_scalar_Ricci_to_T": stats_abs(ricci_resid[mask]),
        "residual_best_algebraic_local_phi_f_to_T": stats_abs(alg_resid[mask]),
        "residual_invariant_trace_ricci_phi_f_to_T": stats_abs(resid_alg_inv[mask]),
        "residual_best_full_local_f_jets_to_T": stats_abs(full_resid[mask]),
        "residual_current_D_algebraic_to_T": stats_abs(d_alg_resid[mask]),
        "residual_current_D_full_lhs_to_T": stats_abs(d_full_resid[mask]),
        "coeff_best_scalar_phi_Ricci": stats_signed(coeff_ricci[..., 0][mask]),
        "coeff_algebraic_phi": stats_signed(coeff_alg[..., 0][mask]),
        "coeff_algebraic_f": stats_signed(coeff_alg[..., 1][mask]),
        "coeff_invariant_algebraic_phi": stats_signed(coeff_alg_inv[..., 0][mask]),
        "coeff_invariant_algebraic_f": stats_signed(coeff_alg_inv[..., 1][mask]),
        "coeff_full_phi": stats_signed(coeff_full[..., 0][mask]),
        "coeff_full_f": stats_signed(coeff_full[..., 1][mask]),
        "coeff_full_phi_R": stats_signed(coeff_full[..., 2][mask]),
        "coeff_full_phi_RR": stats_signed(coeff_full[..., 3][mask]),
        "single_value_alg_phi_vs_R": single_value_diagnostic(r_scalar, coeff_alg[..., 0], mask),
        "single_value_alg_f_vs_R": single_value_diagnostic(r_scalar, coeff_alg[..., 1], mask),
        "single_value_invariant_alg_phi_vs_R": single_value_diagnostic(r_scalar, coeff_alg_inv[..., 0], mask),
        "single_value_invariant_alg_f_vs_R": single_value_diagnostic(r_scalar, coeff_alg_inv[..., 1], mask),
        "single_value_full_phi_vs_R": single_value_diagnostic(r_scalar, coeff_full[..., 0], mask),
        "single_value_full_f_vs_R": single_value_diagnostic(r_scalar, coeff_full[..., 1], mask),
        "single_value_full_phi_R_vs_R": single_value_diagnostic(r_scalar, coeff_full[..., 2], mask),
        "single_value_full_phi_RR_vs_R": single_value_diagnostic(r_scalar, coeff_full[..., 3], mask),
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    if args.ell is None:
        args.ell = float(args.ell_over_planck) * PLANCK_LENGTH_EV_INV

    ref = build_physical_reference(args)
    params = ref["params"]
    scale = ref["scale"]
    x_full = np.asarray(ref["x_full"], dtype=float)
    z_full = np.asarray(ref["z_full"], dtype=float)
    old_scale = float(scale.old_dimensionless_scale_ev_inv)
    t_meet_old = float(params.t_meet / old_scale)
    old_time = t_meet_old + float(args.tau_old)
    time_ev_inv = old_time * old_scale
    dt = float(args.dt_old) * old_scale

    window_ev_inv = float(args.window_um) * 1.0e-6 / HBAR_C_EV_M
    ix = np.where(np.abs(x_full) <= window_ev_inv)[0]
    iz = np.where(np.abs(z_full) <= window_ev_inv)[0]
    x = x_full[ix]
    z = z_full[iz]
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    snaps = {
        -2: make_snapshot_pack(args, time_ev_inv - 2.0 * dt, ix, iz, ref),
        -1: make_snapshot_pack(args, time_ev_inv - dt, ix, iz, ref),
        0: make_snapshot_pack(args, time_ev_inv, ix, iz, ref),
        1: make_snapshot_pack(args, time_ev_inv + dt, ix, iz, ref),
        2: make_snapshot_pack(args, time_ev_inv + 2.0 * dt, ix, iz, ref),
    }
    metric_m2 = real_array(snaps[-2]["cov_txz"])
    metric_m = real_array(snaps[-1]["cov_txz"])
    metric_0 = real_array(snaps[0]["cov_txz"])
    metric_p = real_array(snaps[1]["cov_txz"])
    metric_p2 = real_array(snaps[2]["cov_txz"])

    dg_0, d2g_0 = metric_jets_full(metric_m, metric_0, metric_p, dt, dx, dz)
    metric_inv, gamma2, ricci, r_scalar = geometry_data(metric_0, dg_0, d2g_0)

    dg_m, d2g_m = metric_jets_full(metric_m2, metric_m, metric_0, dt, dx, dz)
    _, _, _, r_scalar_m = geometry_data(metric_m, dg_m, d2g_m)
    dg_p, d2g_p = metric_jets_full(metric_0, metric_p, metric_p2, dt, dx, dz)
    _, _, _, r_scalar_p = geometry_data(metric_p, dg_p, d2g_p)

    dR, d2R = scalar_jets_full(r_scalar_m, r_scalar, r_scalar_p, dt, dx, dz)
    cov_hess_R = covariant_hessian_from_jets(dR, d2R, gamma2)
    box_R = np.einsum("...ab,...ab->...", metric_inv, cov_hess_R, optimize=True)
    grad_R_sq = np.einsum("...ab,...a,...b->...", metric_inv, dR, dR, optimize=True)

    bohm = snaps[0]["bohm"]
    rho_a = real_array(bohm["rho"])
    x_field = real_array(bohm["X"])
    measure = np.abs(x_field) * rho_a / (float(params.m) * float(params.m))
    sqrt_abs_g = safe_sqrt_abs_det(metric_0)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)
    stress, tilde_x = stress_tensor_tilde(
        metric_0,
        metric_inv,
        rho_tilde,
        real_array(bohm["s_t"]),
        real_array(bohm["s_x"]),
        real_array(bohm["s_z"]),
        m=float(params.m),
    )
    target = real_array(stress) / (float(args.mp) * float(args.mp))
    rhs_norm = tensor_norm(target)

    support = (rho_a > float(args.support_rho_frac) * float(np.max(rho_a))) & (
        measure > float(args.support_measure_frac) * float(np.max(measure))
    )
    trusted = erode_mask_8(support, int(args.trusted_erosion))
    if not np.any(trusted):
        trusted = support.copy()

    b_ricci = ricci[:, :, None, :, :]
    b_alg = np.stack([ricci, -0.5 * metric_0], axis=2)
    b_phi_r = metric_0 * box_R[..., None, None] - cov_hess_R
    b_phi_rr = metric_0 * grad_R_sq[..., None, None] - np.einsum("...a,...b->...ab", dR, dR, optimize=True)
    b_full = np.stack([ricci, -0.5 * metric_0, b_phi_r, b_phi_rr], axis=2)

    coeff_ricci, resid_ricci, _ = fit_local_basis(b_ricci, target, support, float(args.lstsq_rcond))
    coeff_alg, resid_alg, _ = fit_local_basis(b_alg, target, support, float(args.lstsq_rcond))
    coeff_full, resid_full, _ = fit_local_basis(b_full, target, support, float(args.lstsq_rcond))
    coeff_alg_inv, resid_alg_inv = invariant_algebraic_coefficients(metric_inv, ricci, metric_0, target, support)

    f_d, phi_d, raw_y = branch_d_f_phi(r_scalar, float(args.ell))
    _, phi_m, _ = branch_d_f_phi(r_scalar_m, float(args.ell))
    _, phi_p, _ = branch_d_f_phi(r_scalar_p, float(args.ell))
    dphi_d, d2phi_d = scalar_jets_full(phi_m, phi_d, phi_p, dt, dx, dz)
    cov_hess_phi_d = covariant_hessian_from_jets(dphi_d, d2phi_d, gamma2)
    box_phi_d = np.einsum("...ab,...ab->...", metric_inv, cov_hess_phi_d, optimize=True)
    d_algebraic = phi_d[..., None, None] * ricci - 0.5 * f_d[..., None, None] * metric_0
    d_derivative = metric_0 * box_phi_d[..., None, None] - cov_hess_phi_d
    d_lhs = d_algebraic + d_derivative
    d_alg_resid = tensor_norm(d_algebraic - target) / np.maximum(rhs_norm, 1.0e-300)
    d_full_resid = tensor_norm(d_lhs - target) / np.maximum(rhs_norm, 1.0e-300)

    ricci_cos = cosine_to_target(ricci, target)
    d_alg_cos = cosine_to_target(d_algebraic, target)
    d_full_cos = cosine_to_target(d_lhs, target)

    maps_path = out / "metric_fr_direction_residual_maps.png"
    scatter_path = out / "metric_fr_candidate_coefficients_vs_R.png"
    render_maps(maps_path, x, z, rho_a, support, resid_ricci, resid_alg, resid_full, d_full_resid)
    render_scatter(scatter_path, r_scalar, coeff_alg, coeff_full, trusted, int(args.scatter_sample_count))

    np.savez_compressed(
        out / "metric_fr_direction_matching_fields.npz",
        x=x,
        z=z,
        rho_A=rho_a,
        support=support,
        trusted=trusted,
        R_tilde=r_scalar,
        raw_y=raw_y,
        rhs_norm=rhs_norm,
        coeff_ricci=coeff_ricci,
        coeff_alg=coeff_alg,
        coeff_full=coeff_full,
        coeff_alg_inv=coeff_alg_inv,
        resid_alg_inv=resid_alg_inv,
        resid_ricci=resid_ricci,
        resid_alg=resid_alg,
        resid_full=resid_full,
        current_D_algebraic_resid=d_alg_resid,
        current_D_full_resid=d_full_resid,
        cos_Ricci_to_T=ricci_cos,
        cos_current_D_algebraic_to_T=d_alg_cos,
        cos_current_D_full_lhs_to_T=d_full_cos,
        mass_shell_defect=real_array(tilde_x) - float(params.m) * float(params.m),
    )

    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "tau_old": float(args.tau_old),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(dt),
            "mass": float(params.m),
            "mp": float(args.mp),
            "ell": float(args.ell),
            "ell_over_planck_length": float(args.ell) / PLANCK_LENGTH_EV_INV,
        },
        "definitions": {
            "target": "lower-index Ttilde_mn/M_P^2 in the txz coordinate chart; with the current pullback convention this is the matter-side tensor to match.",
            "residual_best_scalar_Ricci_to_T": "min_phi ||phi*Rtilde_mn - target||/||target|| pointwise.",
            "residual_best_algebraic_local_phi_f_to_T": "min_phi,f ||phi*Rtilde_mn - 0.5*f*gtilde_mn - target||/||target|| pointwise.",
            "residual_invariant_trace_ricci_phi_f_to_T": "phi,f solved from invariant trace and Ricci-contraction equations, then checked against the full tensor target.",
            "residual_best_full_local_f_jets_to_T": "min_phi,f,phi_R,phi_RR ||phi*Rmn - 0.5*f*gmn + phi_R*(g BoxR-HessR)+phi_RR*(g gradR^2-dR dR) - target||/||target|| pointwise.",
            "single_value_diagnostic": "Bins by sign(R)*log10(1+|R|). Small within-bin MAD, especially median_global_mad, is a necessary condition for coefficients to be functions of R only.",
            "inner_product": "Euclidean Frobenius product of lower-index txz components; used as a stable coordinate diagnostic, not as a final invariant statement.",
        },
        "summary": {
            "support": summarize_region(
                "support",
                support,
                r_scalar,
                rhs_norm,
                ricci_cos,
                d_alg_cos,
                d_full_cos,
                resid_ricci,
                resid_alg,
                resid_full,
                d_alg_resid,
                d_full_resid,
                coeff_ricci,
                coeff_alg,
                coeff_full,
                coeff_alg_inv,
                resid_alg_inv,
            ),
            "trusted": summarize_region(
                "trusted",
                trusted,
                r_scalar,
                rhs_norm,
                ricci_cos,
                d_alg_cos,
                d_full_cos,
                resid_ricci,
                resid_alg,
                resid_full,
                d_alg_resid,
                d_full_resid,
                coeff_ricci,
                coeff_alg,
                coeff_full,
                coeff_alg_inv,
                resid_alg_inv,
            ),
            "core_rho_1pct": summarize_region(
                "core_rho_1pct",
                support & (rho_a > 1.0e-2 * float(np.max(rho_a))),
                r_scalar,
                rhs_norm,
                ricci_cos,
                d_alg_cos,
                d_full_cos,
                resid_ricci,
                resid_alg,
                resid_full,
                d_alg_resid,
                d_full_resid,
                coeff_ricci,
                coeff_alg,
                coeff_full,
                coeff_alg_inv,
                resid_alg_inv,
            ),
            "core_rho_10pct": summarize_region(
                "core_rho_10pct",
                support & (rho_a > 1.0e-1 * float(np.max(rho_a))),
                r_scalar,
                rhs_norm,
                ricci_cos,
                d_alg_cos,
                d_full_cos,
                resid_ricci,
                resid_alg,
                resid_full,
                d_alg_resid,
                d_full_resid,
                coeff_ricci,
                coeff_alg,
                coeff_full,
                coeff_alg_inv,
                resid_alg_inv,
            ),
        },
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "fields_npz": str((out / "metric_fr_direction_matching_fields.npz").resolve()),
            "residual_maps_png": str(maps_path.resolve()),
            "candidate_coefficients_scatter_png": str(scatter_path.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau-old", type=float, default=0.0)
    parser.add_argument("--full-resolution", type=int, default=640)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument("--lstsq-rcond", type=float, default=1.0e-12)
    parser.add_argument("--scatter-sample-count", type=int, default=20000)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
