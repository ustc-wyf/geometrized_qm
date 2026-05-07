from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from diagnose_metric_fr_direction_matching import (
    fit_local_basis,
    single_value_diagnostic,
    stats_abs,
    stats_signed,
    tensor_norm,
)
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference, make_cropped_snapshot, real_array
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det


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


def invariant_contract2(metric_inv: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.einsum("...ac,...bd,...ab,...cd->...", metric_inv, metric_inv, a, b, optimize=True)


def symlog_values(values: np.ndarray) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    return np.sign(vals) * np.log10(1.0 + np.abs(vals))


def single_value_2d_diagnostic(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    mask: np.ndarray,
    bins: int = 32,
    min_count: int = 20,
) -> dict[str, float]:
    xv = np.asarray(x[mask], dtype=float)
    yv = np.asarray(y[mask], dtype=float)
    zv = np.asarray(z[mask], dtype=float)
    good = np.isfinite(xv) & np.isfinite(yv) & np.isfinite(zv)
    xv = xv[good]
    yv = yv[good]
    zv = zv[good]
    if xv.size < min_count:
        return {
            "usable_bins": 0,
            "median_relative_mad": np.nan,
            "p95_relative_mad": np.nan,
            "median_global_mad": np.nan,
            "p95_global_mad": np.nan,
            "global_abs_p95": np.nan,
        }

    sx = symlog_values(xv)
    sy = symlog_values(yv)
    x_edges = np.linspace(float(np.min(sx)), float(np.max(sx)), int(bins) + 1)
    y_edges = np.linspace(float(np.min(sy)), float(np.max(sy)), int(bins) + 1)
    global_scale = max(float(np.percentile(np.abs(zv), 95.0)), 1.0e-300)

    rel_mads: list[float] = []
    global_mads: list[float] = []
    for i in range(int(bins)):
        x_sel = (sx >= x_edges[i]) & (sx < x_edges[i + 1])
        if not np.any(x_sel):
            continue
        for j in range(int(bins)):
            sel = x_sel & (sy >= y_edges[j]) & (sy < y_edges[j + 1])
            vals = zv[sel]
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


def invariant_algebraic_coefficients_rp(
    metric_inv: np.ndarray,
    ricci: np.ndarray,
    metric_cov: np.ndarray,
    target: np.ndarray,
    mask: np.ndarray,
    rcond_floor: float = 1.0e-14,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve pointwise from invariant contractions.

    We solve the 3x3 system built from contractions with g^{ab}, R^{ab}, and
    Q^{ab}=R^a{}_c R^{cb}.  This is a stricter invariant diagnostic than a
    coordinate Frobenius fit.
    """
    nx, nz = ricci.shape[:2]
    coeff = np.full((nx, nz, 3), np.nan, dtype=float)
    residual = np.full((nx, nz), np.nan, dtype=float)
    cond = np.full((nx, nz), np.nan, dtype=float)
    rhs_norm = tensor_norm(target)

    ricci_up = np.einsum("...ac,...cb->...ab", metric_inv, ricci, optimize=True)
    q_tensor = np.einsum("...ac,...cb->...ab", ricci, ricci_up, optimize=True)

    r_scalar = np.einsum("...ab,...ab->...", metric_inv, ricci, optimize=True)
    p_scalar = invariant_contract2(metric_inv, ricci, ricci)
    q_scalar = invariant_contract2(metric_inv, ricci, q_tensor)
    s_scalar = invariant_contract2(metric_inv, q_tensor, q_tensor)

    t0 = np.einsum("...ab,...ab->...", metric_inv, target, optimize=True)
    t1 = invariant_contract2(metric_inv, ricci, target)
    t2 = invariant_contract2(metric_inv, q_tensor, target)

    flat_mask = mask.reshape(-1)
    coeff_flat = coeff.reshape(-1, 3)
    residual_flat = residual.reshape(-1)
    cond_flat = cond.reshape(-1)

    r_flat = r_scalar.reshape(-1)
    p_flat = p_scalar.reshape(-1)
    q_flat = q_scalar.reshape(-1)
    s_flat = s_scalar.reshape(-1)
    t0_flat = t0.reshape(-1)
    t1_flat = t1.reshape(-1)
    t2_flat = t2.reshape(-1)
    target_flat = target.reshape(-1, 9)
    rhs_flat = rhs_norm.reshape(-1)
    ricci_flat = ricci.reshape(-1, 9)
    q_flat_tensor = q_tensor.reshape(-1, 9)
    metric_flat = metric_cov.reshape(-1, 9)

    for idx in np.where(flat_mask)[0]:
        a = np.array(
            [
                [r_flat[idx], 2.0 * p_flat[idx], -0.5 * 3.0],
                [p_flat[idx], 2.0 * q_flat[idx], -0.5 * r_flat[idx]],
                [q_flat[idx], 2.0 * s_flat[idx], -0.5 * p_flat[idx]],
            ],
            dtype=float,
        )
        b = np.array([t0_flat[idx], t1_flat[idx], t2_flat[idx]], dtype=float)
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        if np.linalg.norm(a) <= 1.0e-300:
            continue
        sol, *_ = np.linalg.lstsq(a, b, rcond=rcond_floor)
        coeff_flat[idx] = sol
        fitted = (
            sol[0] * ricci_flat[idx].reshape(3, 3)
            + 2.0 * sol[1] * q_flat_tensor[idx].reshape(3, 3)
            - 0.5 * sol[2] * metric_flat[idx].reshape(3, 3)
        )
        residual_flat[idx] = float(np.linalg.norm(fitted - target_flat[idx].reshape(3, 3)) / max(rhs_flat[idx], 1.0e-300))
        svals = np.linalg.svd(a, compute_uv=False)
        cond_flat[idx] = float(np.max(svals) / max(np.min(svals), 1.0e-300))

    return coeff, residual, cond


def render_maps(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    support: np.ndarray,
    r_scalar: np.ndarray,
    i2_scalar: np.ndarray,
    resid_r_only: np.ndarray,
    resid_rp: np.ndarray,
) -> None:
    extent = [
        float(x[0] * HBAR_C_EV_M * 1.0e6),
        float(x[-1] * HBAR_C_EV_M * 1.0e6),
        float(z[0] * HBAR_C_EV_M * 1.0e6),
        float(z[-1] * HBAR_C_EV_M * 1.0e6),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15.8, 9.2), constrained_layout=True)

    im0 = axes[0, 0].imshow(rho.T, origin="lower", extent=extent, aspect="equal", cmap="magma")
    axes[0, 0].contour(x * HBAR_C_EV_M * 1.0e6, z * HBAR_C_EV_M * 1.0e6, support.T.astype(float), levels=[0.5], colors=["white"], linewidths=0.6)
    axes[0, 0].set_title("rho; white=support")
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046)

    im1 = axes[0, 1].imshow(symlog_values(r_scalar).T, origin="lower", extent=extent, aspect="equal", cmap="coolwarm")
    axes[0, 1].set_title("sign(R) log10(1+|R|)")
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046)

    im2 = axes[0, 2].imshow(symlog_values(i2_scalar).T, origin="lower", extent=extent, aspect="equal", cmap="coolwarm")
    axes[0, 2].set_title("sign(I2) log10(1+|I2|)")
    fig.colorbar(im2, ax=axes[0, 2], fraction=0.046)

    im3 = axes[1, 0].imshow(np.log10(1.0e-300 + np.clip(resid_r_only, 0.0, 1.0e16)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
    axes[1, 0].set_title("log10 residual pure f(R)")
    fig.colorbar(im3, ax=axes[1, 0], fraction=0.046)

    im4 = axes[1, 1].imshow(np.log10(1.0e-300 + np.clip(resid_rp, 0.0, 1.0e16)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
    axes[1, 1].set_title("log10 residual f(R,I2)")
    fig.colorbar(im4, ax=axes[1, 1], fraction=0.046)

    ratio = np.maximum(resid_r_only, 1.0e-300) / np.maximum(resid_rp, 1.0e-300)
    im5 = axes[1, 2].imshow(np.log10(1.0e-300 + ratio).T, origin="lower", extent=extent, aspect="equal", cmap="plasma")
    axes[1, 2].set_title("log10 (pure / RP)")
    fig.colorbar(im5, ax=axes[1, 2], fraction=0.046)

    for ax in axes.ravel():
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_scatter(
    out_path: Path,
    r_scalar: np.ndarray,
    i2_scalar: np.ndarray,
    coeff_r_only: np.ndarray,
    coeff_rp_local: np.ndarray,
    coeff_rp_inv: np.ndarray,
    mask: np.ndarray,
    sample_count: int,
) -> None:
    r = np.asarray(r_scalar[mask], dtype=float)
    i2 = np.asarray(i2_scalar[mask], dtype=float)
    c_r = np.asarray(coeff_r_only[mask], dtype=float)
    c_lp = np.asarray(coeff_rp_local[mask], dtype=float)
    c_iv = np.asarray(coeff_rp_inv[mask], dtype=float)
    good = (
        np.isfinite(r)
        & np.isfinite(i2)
        & np.all(np.isfinite(c_r), axis=1)
        & np.all(np.isfinite(c_lp), axis=1)
        & np.all(np.isfinite(c_iv), axis=1)
    )
    r = r[good]
    i2 = i2[good]
    c_r = c_r[good]
    c_lp = c_lp[good]
    c_iv = c_iv[good]
    if r.size > sample_count:
        rng = np.random.default_rng(12345)
        take = rng.choice(r.size, size=sample_count, replace=False)
        r = r[take]
        i2 = i2[take]
        c_r = c_r[take]
        c_lp = c_lp[take]
        c_iv = c_iv[take]

    fig, axes = plt.subplots(2, 3, figsize=(16.4, 9.6), constrained_layout=True)
    panels = [
        (c_r[:, 0] if c_r.size else np.array([]), "pure phi", axes[0, 0]),
        (c_r[:, 1] if c_r.size else np.array([]), "pure f", axes[0, 1]),
        (c_lp[:, 0] if c_lp.size else np.array([]), "RP local phi", axes[0, 2]),
        (c_lp[:, 1] if c_lp.size else np.array([]), "RP local psi", axes[1, 0]),
        (c_lp[:, 2] if c_lp.size else np.array([]), "RP local f", axes[1, 1]),
        (c_iv[:, 1] if c_iv.size else np.array([]), "RP invariant psi", axes[1, 2]),
    ]
    for y, title, ax in panels:
        if r.size and y.size:
            sc = ax.scatter(symlog_values(r), symlog_values(i2), c=y, s=2, alpha=0.18, linewidths=0, cmap="coolwarm")
            fig.colorbar(sc, ax=ax, fraction=0.046)
        ax.set_xlabel("sign(R) log10(1+|R|)")
        ax.set_ylabel("sign(I2) log10(1+|I2|)")
        ax.set_title(title)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def summarize_region(
    name: str,
    mask: np.ndarray,
    r_scalar: np.ndarray,
    i2_scalar: np.ndarray,
    rhs_norm: np.ndarray,
    coeff_r_only: np.ndarray,
    resid_r_only: np.ndarray,
    coeff_rp_local: np.ndarray,
    resid_rp_local: np.ndarray,
    coeff_rp_inv: np.ndarray,
    resid_rp_inv: np.ndarray,
) -> dict[str, object]:
    return {
        "name": name,
        "count": int(np.count_nonzero(mask)),
        "R_tilde": stats_signed(r_scalar[mask]),
        "I2_tilde": stats_signed(i2_scalar[mask]),
        "rhs_norm": stats_abs(rhs_norm[mask]),
        "residual_pure_R_local": stats_abs(resid_r_only[mask]),
        "residual_RP_local": stats_abs(resid_rp_local[mask]),
        "residual_RP_invariant": stats_abs(resid_rp_inv[mask]),
        "coeff_pure_R_phi": stats_signed(coeff_r_only[..., 0][mask]),
        "coeff_pure_R_f": stats_signed(coeff_r_only[..., 1][mask]),
        "coeff_RP_local_phi": stats_signed(coeff_rp_local[..., 0][mask]),
        "coeff_RP_local_psi": stats_signed(coeff_rp_local[..., 1][mask]),
        "coeff_RP_local_f": stats_signed(coeff_rp_local[..., 2][mask]),
        "coeff_RP_inv_phi": stats_signed(coeff_rp_inv[..., 0][mask]),
        "coeff_RP_inv_psi": stats_signed(coeff_rp_inv[..., 1][mask]),
        "coeff_RP_inv_f": stats_signed(coeff_rp_inv[..., 2][mask]),
        "single_value_pure_R_phi_vs_R": single_value_diagnostic(r_scalar, coeff_r_only[..., 0], mask),
        "single_value_pure_R_f_vs_R": single_value_diagnostic(r_scalar, coeff_r_only[..., 1], mask),
        "single_value_RP_local_phi_vs_RI2": single_value_2d_diagnostic(r_scalar, i2_scalar, coeff_rp_local[..., 0], mask),
        "single_value_RP_local_psi_vs_RI2": single_value_2d_diagnostic(r_scalar, i2_scalar, coeff_rp_local[..., 1], mask),
        "single_value_RP_local_f_vs_RI2": single_value_2d_diagnostic(r_scalar, i2_scalar, coeff_rp_local[..., 2], mask),
        "single_value_RP_inv_phi_vs_RI2": single_value_2d_diagnostic(r_scalar, i2_scalar, coeff_rp_inv[..., 0], mask),
        "single_value_RP_inv_psi_vs_RI2": single_value_2d_diagnostic(r_scalar, i2_scalar, coeff_rp_inv[..., 1], mask),
        "single_value_RP_inv_f_vs_RI2": single_value_2d_diagnostic(r_scalar, i2_scalar, coeff_rp_inv[..., 2], mask),
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
        -1: make_snapshot_pack(args, time_ev_inv - dt, ix, iz, ref),
        0: make_snapshot_pack(args, time_ev_inv, ix, iz, ref),
        1: make_snapshot_pack(args, time_ev_inv + dt, ix, iz, ref),
    }

    metric_m = real_array(snaps[-1]["cov_txz"])
    metric_0 = real_array(snaps[0]["cov_txz"])
    metric_p = real_array(snaps[1]["cov_txz"])

    dg_0, d2g_0 = metric_jets_full(metric_m, metric_0, metric_p, dt, dx, dz)
    metric_inv, _, ricci, r_scalar = geometry_data(metric_0, dg_0, d2g_0)

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

    ricci_sq = invariant_contract2(metric_inv, ricci, ricci)
    ricci_up = np.einsum("...ac,...cb->...ab", metric_inv, ricci, optimize=True)
    q_tensor = np.einsum("...ac,...cb->...ab", ricci, ricci_up, optimize=True)
    q_scalar = invariant_contract2(metric_inv, ricci, q_tensor)
    q2_scalar = invariant_contract2(metric_inv, q_tensor, q_tensor)

    b_pure_r = np.stack([ricci, -0.5 * metric_0], axis=2)
    b_rp = np.stack([ricci, 2.0 * q_tensor, -0.5 * metric_0], axis=2)

    coeff_r_only, resid_r_only, _ = fit_local_basis(b_pure_r, target, support, float(args.lstsq_rcond))
    coeff_rp_local, resid_rp_local, _ = fit_local_basis(b_rp, target, support, float(args.lstsq_rcond))
    coeff_rp_inv, resid_rp_inv, inv_cond = invariant_algebraic_coefficients_rp(metric_inv, ricci, metric_0, target, support)

    resid_r_only_abs = resid_r_only
    resid_rp_local_abs = resid_rp_local

    maps_path = out / "metric_fr_rp_residual_maps.png"
    scatter_path = out / "metric_fr_rp_candidate_coefficients_vs_invariants.png"
    render_maps(maps_path, x, z, rho_a, support, r_scalar, ricci_sq, resid_r_only_abs, resid_rp_local_abs)
    render_scatter(scatter_path, r_scalar, ricci_sq, coeff_r_only, coeff_rp_local, coeff_rp_inv, trusted, int(args.scatter_sample_count))

    np.savez_compressed(
        out / "metric_fr_rp_direction_matching_fields.npz",
        x=x,
        z=z,
        rho_A=rho_a,
        support=support,
        trusted=trusted,
        R_tilde=r_scalar,
        I2_tilde=ricci_sq,
        Q_scalar=q_scalar,
        Q2_scalar=q2_scalar,
        rhs_norm=rhs_norm,
        coeff_pure_R=coeff_r_only,
        coeff_RP_local=coeff_rp_local,
        coeff_RP_inv=coeff_rp_inv,
        resid_pure_R=resid_r_only_abs,
        resid_RP_local=resid_rp_local_abs,
        resid_RP_inv=resid_rp_inv,
        inv_condition=inv_cond,
        mass_shell_defect=real_array(tilde_x) - float(params.m) * float(params.m),
    )

    summary = {
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
            "target": "lower-index Ttilde_mn/M_P^2 in the txz chart.",
            "pure_fR_local": "min_phi,f ||phi*Rmn - 1/2 f gmn - target||/||target|| pointwise.",
            "f_R_I2_local": "min_phi,psi,f ||phi*Rmn + 2 psi*Qmn - 1/2 f gmn - target||/||target|| pointwise.",
            "f_R_I2_invariant": "phi,psi,f solved from trace/Ricci/Q contractions, then checked against the full tensor target.",
            "single_value_2d": "Bins by sign(R)log10(1+|R|) and sign(I2)log10(1+|I2|); small within-bin MAD is necessary for a 2-variable function.",
        },
        "summary": {
            "support": summarize_region(
                "support",
                support,
                r_scalar,
                ricci_sq,
                rhs_norm,
                coeff_r_only,
                resid_r_only_abs,
                coeff_rp_local,
                resid_rp_local_abs,
                coeff_rp_inv,
                resid_rp_inv,
            ),
            "trusted": summarize_region(
                "trusted",
                trusted,
                r_scalar,
                ricci_sq,
                rhs_norm,
                coeff_r_only,
                resid_r_only_abs,
                coeff_rp_local,
                resid_rp_local_abs,
                coeff_rp_inv,
                resid_rp_inv,
            ),
            "core_rho_1pct": summarize_region(
                "core_rho_1pct",
                support & (rho_a > 1.0e-2 * float(np.max(rho_a))),
                r_scalar,
                ricci_sq,
                rhs_norm,
                coeff_r_only,
                resid_r_only_abs,
                coeff_rp_local,
                resid_rp_local_abs,
                coeff_rp_inv,
                resid_rp_inv,
            ),
            "core_rho_10pct": summarize_region(
                "core_rho_10pct",
                support & (rho_a > 1.0e-1 * float(np.max(rho_a))),
                r_scalar,
                ricci_sq,
                rhs_norm,
                coeff_r_only,
                resid_r_only_abs,
                coeff_rp_local,
                resid_rp_local_abs,
                coeff_rp_inv,
                resid_rp_inv,
            ),
        },
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "fields_npz": str((out / "metric_fr_rp_direction_matching_fields.npz").resolve()),
            "residual_maps_png": str(maps_path.resolve()),
            "candidate_coefficients_scatter_png": str(scatter_path.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return summary


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
