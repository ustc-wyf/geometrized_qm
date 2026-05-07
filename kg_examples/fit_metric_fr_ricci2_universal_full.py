from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev as cheb

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import d1x, d1z, d2xx, d2xz, d2zz, metric_jets_full
from diagnose_metric_fr_direction_matching import stats_abs, tensor_norm
from diagnose_metric_fr_ricci2_direction_matching import invariant_contract2
from fit_metric_fr_ricci2_universal_algebraic import basis_values_2d, make_fit_scales
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference, make_cropped_snapshot, real_array
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det


@dataclass
class TimeGeometry:
    metric_cov: np.ndarray
    metric_inv: np.ndarray
    gamma: np.ndarray
    ricci: np.ndarray
    r_scalar: np.ndarray
    p_scalar: np.ndarray
    ricci_square_tensor: np.ndarray
    ricci_mixed: np.ndarray


@dataclass
class CaseFullData:
    tau_old: float
    x: np.ndarray
    z: np.ndarray
    rho: np.ndarray
    target: np.ndarray
    rhs_norm: np.ndarray
    masks: dict[str, np.ndarray]
    geom_m: TimeGeometry
    geom_0: TimeGeometry
    geom_p: TimeGeometry
    dgamma_0: np.ndarray


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


def geometry_data_with_dgamma(
    metric_cov: np.ndarray,
    dg: np.ndarray,
    d2g: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nx, nz = metric_cov.shape[:2]
    ginv = np.zeros_like(metric_cov)
    gamma2 = np.zeros(metric_cov.shape[:2] + (3, 3, 3), dtype=float)
    dgamma2 = np.zeros(metric_cov.shape[:2] + (3, 3, 3, 3), dtype=float)
    ricci = np.zeros(metric_cov.shape[:2] + (3, 3), dtype=float)
    r_scalar = np.zeros(metric_cov.shape[:2], dtype=float)

    for i in range(nx):
        for j in range(nz):
            g = metric_cov[i, j]
            ginv_ij = np.linalg.pinv(g, rcond=1.0e-12, hermitian=True)
            ginv[i, j] = ginv_ij

            gamma1 = np.zeros((3, 3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        gamma1[a, b, c] = 0.5 * (
                            dg[i, j, b, a, c] + dg[i, j, c, a, b] - dg[i, j, a, b, c]
                        )
            gamma2_ij = np.einsum("ad,dbc->abc", ginv_ij, gamma1, optimize=True)
            gamma2[i, j] = gamma2_ij

            dginv = np.zeros((3, 3, 3), dtype=float)
            for e in range(3):
                dginv[e] = -ginv_ij @ dg[i, j, e] @ ginv_ij

            dgamma1 = np.zeros((3, 3, 3, 3), dtype=float)
            for e in range(3):
                for a in range(3):
                    for b in range(3):
                        for c in range(3):
                            dgamma1[e, a, b, c] = 0.5 * (
                                d2g[i, j, e, b, a, c]
                                + d2g[i, j, e, c, a, b]
                                - d2g[i, j, e, a, b, c]
                            )

            dgamma2_ij = np.zeros((3, 3, 3, 3), dtype=float)
            for e in range(3):
                dgamma2_ij[e] = np.einsum("ad,dbc->abc", dginv[e], gamma1, optimize=True) + np.einsum(
                    "ad,dbc->abc", ginv_ij, dgamma1[e], optimize=True
                )
            dgamma2[i, j] = dgamma2_ij

            ric = np.zeros((3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    term1 = sum(dgamma2_ij[c, c, a, b] for c in range(3))
                    term2 = sum(dgamma2_ij[b, c, a, c] for c in range(3))
                    quad1 = 0.0
                    quad2 = 0.0
                    for c in range(3):
                        for d in range(3):
                            quad1 += gamma2_ij[c, a, b] * gamma2_ij[d, c, d]
                            quad2 += gamma2_ij[c, a, d] * gamma2_ij[d, b, c]
                    ric[a, b] = term1 - term2 + quad1 - quad2
            ricci[i, j] = ric
            r_scalar[i, j] = float(np.einsum("ab,ab->", ginv_ij, ric))
    return ginv, gamma2, dgamma2, ricci, r_scalar


def make_time_geometry(metric_cov: np.ndarray, metric_m: np.ndarray, metric_p: np.ndarray, dt: float, dx: float, dz: float, with_dgamma: bool) -> tuple[TimeGeometry, np.ndarray | None]:
    dg, d2g = metric_jets_full(metric_m, metric_cov, metric_p, dt, dx, dz)
    if with_dgamma:
        metric_inv, gamma, dgamma, ricci, r_scalar = geometry_data_with_dgamma(metric_cov, dg, d2g)
    else:
        metric_inv, gamma, ricci, r_scalar = geometry_data(metric_cov, dg, d2g)
        dgamma = None
    ricci_mixed = np.einsum("...ac,...cb->...ab", metric_inv, ricci, optimize=True)
    ricci_square_tensor = np.einsum("...ac,...cb->...ab", ricci, ricci_mixed, optimize=True)
    p_scalar = invariant_contract2(metric_inv, ricci, ricci)
    return (
        TimeGeometry(
            metric_cov=metric_cov,
            metric_inv=metric_inv,
            gamma=gamma,
            ricci=ricci,
            r_scalar=r_scalar,
            p_scalar=p_scalar,
            ricci_square_tensor=ricci_square_tensor,
            ricci_mixed=ricci_mixed,
        ),
        dgamma,
    )


def field_jets_full(f_m: np.ndarray, f_0: np.ndarray, f_p: np.ndarray, dt: float, dx: float, dz: float) -> tuple[np.ndarray, np.ndarray]:
    tail = f_0.shape[2:]
    df = np.zeros(f_0.shape[:2] + (3,) + tail, dtype=float)
    d2f = np.zeros(f_0.shape[:2] + (3, 3) + tail, dtype=float)
    tail_slice = (slice(None),) * len(tail)

    f_t = (f_p - f_m) / (2.0 * dt)
    f_tt = (f_p - 2.0 * f_0 + f_m) / (dt * dt)

    df[(slice(None), slice(None), 0) + tail_slice] = f_t
    df[(slice(None), slice(None), 1) + tail_slice] = d1x(f_0, dx)
    df[(slice(None), slice(None), 2) + tail_slice] = d1z(f_0, dz)

    d2f[(slice(None), slice(None), 0, 0) + tail_slice] = f_tt
    d2f[(slice(None), slice(None), 0, 1) + tail_slice] = d1x(f_t, dx)
    d2f[(slice(None), slice(None), 1, 0) + tail_slice] = d2f[(slice(None), slice(None), 0, 1) + tail_slice]
    d2f[(slice(None), slice(None), 0, 2) + tail_slice] = d1z(f_t, dz)
    d2f[(slice(None), slice(None), 2, 0) + tail_slice] = d2f[(slice(None), slice(None), 0, 2) + tail_slice]
    d2f[(slice(None), slice(None), 1, 1) + tail_slice] = d2xx(f_0, dx)
    d2f[(slice(None), slice(None), 1, 2) + tail_slice] = d2xz(f_0, dx, dz)
    d2f[(slice(None), slice(None), 2, 1) + tail_slice] = d2f[(slice(None), slice(None), 1, 2) + tail_slice]
    d2f[(slice(None), slice(None), 2, 2) + tail_slice] = d2zz(f_0, dz)
    return df, d2f


def covariant_hessian_scalar(df: np.ndarray, d2f: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    return d2f - np.einsum("...cab,...c->...ab", gamma, df, optimize=True)


def covariant_hessian_cov2(
    tensor: np.ndarray,
    dtensor: np.ndarray,
    d2tensor: np.ndarray,
    gamma: np.ndarray,
    dgamma: np.ndarray,
) -> np.ndarray:
    nx, nz = tensor.shape[:2]
    hess = np.zeros(tensor.shape[:2] + (3, 3, 3, 3), dtype=float)
    for i in range(nx):
        for j in range(nz):
            a0 = tensor[i, j]
            da = dtensor[i, j]
            d2a = d2tensor[i, j]
            g = gamma[i, j]
            dg = dgamma[i, j]
            first = np.zeros((3, 3, 3), dtype=float)
            for beta in range(3):
                for mu in range(3):
                    for nu in range(3):
                        val = da[beta, mu, nu]
                        for lam in range(3):
                            val -= g[lam, beta, mu] * a0[lam, nu]
                            val -= g[lam, beta, nu] * a0[mu, lam]
                        first[beta, mu, nu] = val
            for alpha in range(3):
                for beta in range(3):
                    for mu in range(3):
                        for nu in range(3):
                            val = d2a[alpha, beta, mu, nu]
                            for lam in range(3):
                                val -= dg[alpha, lam, beta, mu] * a0[lam, nu]
                                val -= g[lam, beta, mu] * da[alpha, lam, nu]
                                val -= dg[alpha, lam, beta, nu] * a0[mu, lam]
                                val -= g[lam, beta, nu] * da[alpha, mu, lam]
                                val -= g[lam, alpha, beta] * first[lam, mu, nu]
                                val -= g[lam, alpha, mu] * first[beta, lam, nu]
                                val -= g[lam, alpha, nu] * first[beta, mu, lam]
                            hess[i, j, alpha, beta, mu, nu] = val
    return hess


def covariant_hessian_mixed11(
    tensor: np.ndarray,
    dtensor: np.ndarray,
    d2tensor: np.ndarray,
    gamma: np.ndarray,
    dgamma: np.ndarray,
) -> np.ndarray:
    nx, nz = tensor.shape[:2]
    hess = np.zeros(tensor.shape[:2] + (3, 3, 3, 3), dtype=float)
    for i in range(nx):
        for j in range(nz):
            c0 = tensor[i, j]
            dc = dtensor[i, j]
            d2c = d2tensor[i, j]
            g = gamma[i, j]
            dg = dgamma[i, j]
            first = np.zeros((3, 3, 3), dtype=float)
            for beta in range(3):
                for up in range(3):
                    for down in range(3):
                        val = dc[beta, up, down]
                        for lam in range(3):
                            val += g[up, beta, lam] * c0[lam, down]
                            val -= g[lam, beta, down] * c0[up, lam]
                        first[beta, up, down] = val
            for alpha in range(3):
                for beta in range(3):
                    for up in range(3):
                        for down in range(3):
                            val = d2c[alpha, beta, up, down]
                            for lam in range(3):
                                val += dg[alpha, up, beta, lam] * c0[lam, down]
                                val += g[up, beta, lam] * dc[alpha, lam, down]
                                val -= dg[alpha, lam, beta, down] * c0[up, lam]
                                val -= g[lam, beta, down] * dc[alpha, up, lam]
                                val -= g[lam, alpha, beta] * first[lam, up, down]
                                val += g[up, alpha, lam] * first[beta, lam, down]
                                val -= g[lam, alpha, down] * first[beta, up, lam]
                            hess[i, j, alpha, beta, up, down] = val
    return hess


def full_lhs_from_fields(
    geom_m: TimeGeometry,
    geom_0: TimeGeometry,
    geom_p: TimeGeometry,
    dgamma_0: np.ndarray,
    f_m: np.ndarray,
    f_0: np.ndarray,
    f_p: np.ndarray,
    fr_m: np.ndarray,
    fr_0: np.ndarray,
    fr_p: np.ndarray,
    fy_m: np.ndarray,
    fy_0: np.ndarray,
    fy_p: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    dfr, d2fr = field_jets_full(fr_m, fr_0, fr_p, dt, dx, dz)
    hess_fr = covariant_hessian_scalar(dfr, d2fr, geom_0.gamma)
    box_fr = np.einsum("...ab,...ab->...", geom_0.metric_inv, hess_fr, optimize=True)

    a_m = fy_m[..., None, None] * geom_m.ricci
    a_0 = fy_0[..., None, None] * geom_0.ricci
    a_p = fy_p[..., None, None] * geom_p.ricci
    da, d2a = field_jets_full(a_m, a_0, a_p, dt, dx, dz)
    hess_a = covariant_hessian_cov2(a_0, da, d2a, geom_0.gamma, dgamma_0)
    box_a = np.einsum("...rs,...rsmn->...mn", geom_0.metric_inv, hess_a, optimize=True)
    divdiv_a = np.einsum("...ar,...bs,...rsab->...", geom_0.metric_inv, geom_0.metric_inv, hess_a, optimize=True)

    c_m = fy_m[..., None, None] * geom_m.ricci_mixed
    c_0 = fy_0[..., None, None] * geom_0.ricci_mixed
    c_p = fy_p[..., None, None] * geom_p.ricci_mixed
    dc, d2c = field_jets_full(c_m, c_0, c_p, dt, dx, dz)
    hess_c = covariant_hessian_mixed11(c_0, dc, d2c, geom_0.gamma, dgamma_0)
    last = np.zeros_like(geom_0.metric_cov)
    for mu in range(3):
        for nu in range(3):
            val = np.zeros(geom_0.r_scalar.shape, dtype=float)
            for alpha in range(3):
                val += hess_c[..., alpha, nu, alpha, mu]
                val += hess_c[..., alpha, mu, alpha, nu]
            last[..., mu, nu] = -val

    algebraic = (
        fr_0[..., None, None] * geom_0.ricci
        + 2.0 * fy_0[..., None, None] * geom_0.ricci_square_tensor
        - 0.5 * f_0[..., None, None] * geom_0.metric_cov
    )
    scalar_derivative = geom_0.metric_cov * box_fr[..., None, None] - hess_fr
    ricci_square_derivative = box_a + geom_0.metric_cov * divdiv_a[..., None, None] + last
    lhs = algebraic + scalar_derivative + ricci_square_derivative
    return {
        "lhs": lhs,
        "algebraic": algebraic,
        "scalar_derivative": scalar_derivative,
        "ricci_square_derivative": ricci_square_derivative,
    }


def build_case(args: argparse.Namespace, ref: dict[str, object], tau_old: float) -> CaseFullData:
    params = ref["params"]
    scale = ref["scale"]
    x_full = np.asarray(ref["x_full"], dtype=float)
    z_full = np.asarray(ref["z_full"], dtype=float)
    old_scale = float(scale.old_dimensionless_scale_ev_inv)
    t_meet_old = float(params.t_meet / old_scale)
    old_time = t_meet_old + float(tau_old)
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
    metrics = {key: real_array(snap["cov_txz"]) for key, snap in snaps.items()}

    geom_m, _ = make_time_geometry(metrics[-1], metrics[-2], metrics[0], dt, dx, dz, with_dgamma=False)
    geom_0, dgamma_0 = make_time_geometry(metrics[0], metrics[-1], metrics[1], dt, dx, dz, with_dgamma=True)
    geom_p, _ = make_time_geometry(metrics[1], metrics[0], metrics[2], dt, dx, dz, with_dgamma=False)
    if dgamma_0 is None:
        raise RuntimeError("central dgamma was not computed")

    bohm = snaps[0]["bohm"]
    rho = real_array(bohm["rho"])
    x_field = real_array(bohm["X"])
    measure = np.abs(x_field) * rho / (float(params.m) * float(params.m))
    sqrt_abs_g = safe_sqrt_abs_det(geom_0.metric_cov)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)
    stress, _ = stress_tensor_tilde(
        geom_0.metric_cov,
        geom_0.metric_inv,
        rho_tilde,
        real_array(bohm["s_t"]),
        real_array(bohm["s_x"]),
        real_array(bohm["s_z"]),
        m=float(params.m),
    )
    target = real_array(stress) / (float(args.mp) * float(args.mp))
    rhs_norm = tensor_norm(target)

    support = (rho > float(args.support_rho_frac) * float(np.max(rho))) & (
        measure > float(args.support_measure_frac) * float(np.max(measure))
    )
    trusted = erode_mask_8(support, int(args.trusted_erosion))
    if not np.any(trusted):
        trusted = support.copy()
    masks = {
        "support": support,
        "trusted": trusted,
        "core1": support & (rho > 1.0e-2 * float(np.max(rho))),
        "core10": support & (rho > 1.0e-1 * float(np.max(rho))),
    }
    return CaseFullData(
        tau_old=float(tau_old),
        x=x,
        z=z,
        rho=rho,
        target=target,
        rhs_norm=rhs_norm,
        masks=masks,
        geom_m=geom_m,
        geom_0=geom_0,
        geom_p=geom_p,
        dgamma_0=dgamma_0,
    )


def case_basis_values(case: CaseFullData, degree_r: int, degree_p: int, scales: dict[str, float]) -> dict[str, np.ndarray]:
    values: dict[str, np.ndarray] = {}
    for key, geom in [("m", case.geom_m), ("0", case.geom_0), ("p", case.geom_p)]:
        f, fr, fy = basis_values_2d(
            geom.r_scalar.reshape(-1),
            geom.p_scalar.reshape(-1),
            degree_r,
            degree_p,
            scales["r0"],
            scales["p0"],
            scales["r_s_min"],
            scales["r_s_max"],
            scales["p_s_min"],
            scales["p_s_max"],
        )
        shape = geom.r_scalar.shape + (f.shape[1],)
        values[f"f_{key}"] = f.reshape(shape)
        values[f"fr_{key}"] = fr.reshape(shape)
        values[f"fy_{key}"] = fy.reshape(shape)
    return values


def design_rows_for_case_full(
    case: CaseFullData,
    mask: np.ndarray,
    degree_r: int,
    degree_p: int,
    scales: dict[str, float],
    target_floor_frac: float,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    basis = case_basis_values(case, degree_r, degree_p, scales)
    k = basis["f_0"].shape[-1]
    idx = np.where(mask.reshape(-1))[0]
    rows = np.zeros((idx.size * 9, k), dtype=float)
    target = case.target.reshape(-1, 9)[idx]
    target_rows = target.reshape(idx.size * 9)

    for col in range(k):
        parts = full_lhs_from_fields(
            case.geom_m,
            case.geom_0,
            case.geom_p,
            case.dgamma_0,
            basis["f_m"][..., col],
            basis["f_0"][..., col],
            basis["f_p"][..., col],
            basis["fr_m"][..., col],
            basis["fr_0"][..., col],
            basis["fr_p"][..., col],
            basis["fy_m"][..., col],
            basis["fy_0"][..., col],
            basis["fy_p"][..., col],
            dt,
            dx,
            dz,
        )
        rows[:, col] = parts["lhs"].reshape(-1, 9)[idx].reshape(idx.size * 9)

    rho = case.rho.reshape(-1)[idx]
    rhs_norm = case.rhs_norm.reshape(-1)[idx]
    rho_w = np.sqrt(np.maximum(rho, 0.0) / max(float(np.max(case.rho)), 1.0e-300))
    floor = float(target_floor_frac) * max(float(np.percentile(rhs_norm[np.isfinite(rhs_norm)], 95.0)), 1.0e-300)
    rel_w = 1.0 / np.maximum(rhs_norm, floor)
    w = rho_w * rel_w
    w_rows = np.repeat(w, 9)
    return rows * w_rows[:, None], target_rows * w_rows


def fit_universal_full(
    cases: list[CaseFullData],
    fit_region: str,
    degree_r: int,
    degree_p: int,
    target_floor_frac: float,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, object]:
    pseudo_cases = []
    for case in cases:
        pseudo = type("PseudoCase", (), {})()
        pseudo.r_scalar = case.geom_0.r_scalar
        pseudo.p_scalar = case.geom_0.p_scalar
        pseudo.masks = case.masks
        pseudo_cases.append(pseudo)
    scales = make_fit_scales(pseudo_cases, fit_region)
    matrices = [
        design_rows_for_case_full(case, case.masks[fit_region], degree_r, degree_p, scales, target_floor_frac, dt, dx, dz)
        for case in cases
    ]
    a = np.vstack([m[0] for m in matrices])
    b = np.concatenate([m[1] for m in matrices])
    column_scale = np.linalg.norm(a, axis=0)
    column_scale = np.where(column_scale > 0.0, column_scale, 1.0)
    a_scaled = a / column_scale[None, :]
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_scaled, b, rcond=1.0e-12)
    coeff = coeff_scaled / column_scale
    residual_norm = float(np.linalg.norm(a @ coeff - b))
    target_norm = float(np.linalg.norm(b))
    return {
        "degree_r": int(degree_r),
        "degree_p": int(degree_p),
        "fit_region": fit_region,
        "coeff": coeff,
        "column_scale": column_scale,
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "scaled_condition": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
        "weighted_residual_norm": residual_norm,
        "weighted_target_norm": target_norm,
        "weighted_relative_residual": residual_norm / max(target_norm, 1.0e-300),
        **scales,
    }


def evaluate_case_full(case: CaseFullData, fit: dict[str, object], region: str, dt: float, dx: float, dz: float) -> dict[str, object]:
    degree_r = int(fit["degree_r"])
    degree_p = int(fit["degree_p"])
    scales = {key: float(fit[key]) for key in ["r0", "p0", "r_s_min", "r_s_max", "p_s_min", "p_s_max"]}
    basis = case_basis_values(case, degree_r, degree_p, scales)
    coeff = np.asarray(fit["coeff"], dtype=float)
    fields = {}
    for name, value in basis.items():
        fields[name] = value @ coeff
    parts = full_lhs_from_fields(
        case.geom_m,
        case.geom_0,
        case.geom_p,
        case.dgamma_0,
        fields["f_m"],
        fields["f_0"],
        fields["f_p"],
        fields["fr_m"],
        fields["fr_0"],
        fields["fr_p"],
        fields["fy_m"],
        fields["fy_0"],
        fields["fy_p"],
        dt,
        dx,
        dz,
    )
    lhs = parts["lhs"]
    residual = tensor_norm(lhs - case.target) / np.maximum(case.rhs_norm, 1.0e-300)
    lhs_norm = tensor_norm(lhs)
    rel_lhs_rhs = tensor_norm(lhs - case.target) / np.maximum(lhs_norm + case.rhs_norm, 1.0e-300)
    mask = case.masks[region]
    return {
        "count": int(np.count_nonzero(mask)),
        "residual_to_rhs": stats_abs(residual[mask]),
        "residual_to_lhs_plus_rhs": stats_abs(rel_lhs_rhs[mask]),
        "lhs_norm": stats_abs(lhs_norm[mask]),
        "rhs_norm": stats_abs(case.rhs_norm[mask]),
        "algebraic_norm": stats_abs(tensor_norm(parts["algebraic"])[mask]),
        "scalar_derivative_norm": stats_abs(tensor_norm(parts["scalar_derivative"])[mask]),
        "ricci_square_derivative_norm": stats_abs(tensor_norm(parts["ricci_square_derivative"])[mask]),
        "f": stats_abs(fields["f_0"][mask]),
        "f_R": stats_abs(fields["fr_0"][mask]),
        "f_I2": stats_abs(fields["fy_0"][mask]),
    }


def render_residual_plot(out_path: Path, cases: list[CaseFullData], fit: dict[str, object], region: str, dt: float, dx: float, dz: float) -> None:
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.8), constrained_layout=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        degree_r = int(fit["degree_r"])
        degree_p = int(fit["degree_p"])
        scales = {key: float(fit[key]) for key in ["r0", "p0", "r_s_min", "r_s_max", "p_s_min", "p_s_max"]}
        basis = case_basis_values(case, degree_r, degree_p, scales)
        coeff = np.asarray(fit["coeff"], dtype=float)
        fields = {name: value @ coeff for name, value in basis.items()}
        parts = full_lhs_from_fields(
            case.geom_m,
            case.geom_0,
            case.geom_p,
            case.dgamma_0,
            fields["f_m"],
            fields["f_0"],
            fields["f_p"],
            fields["fr_m"],
            fields["fr_0"],
            fields["fr_p"],
            fields["fy_m"],
            fields["fy_0"],
            fields["fy_p"],
            dt,
            dx,
            dz,
        )
        residual = tensor_norm(parts["lhs"] - case.target) / np.maximum(case.rhs_norm, 1.0e-300)
        extent = [
            float(case.x[0] * HBAR_C_EV_M * 1.0e6),
            float(case.x[-1] * HBAR_C_EV_M * 1.0e6),
            float(case.z[0] * HBAR_C_EV_M * 1.0e6),
            float(case.z[-1] * HBAR_C_EV_M * 1.0e6),
        ]
        im = ax.imshow(np.log10(1.0e-6 + np.clip(residual, 0.0, 1.0e6)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        ax.contour(case.x * HBAR_C_EV_M * 1.0e6, case.z * HBAR_C_EV_M * 1.0e6, case.masks[region].T.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        ax.set_title(f"tau={case.tau_old:g}; log10 full residual")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    if args.ell is None:
        args.ell = float(args.ell_over_planck) * PLANCK_LENGTH_EV_INV

    ref = build_physical_reference(args)
    cases = [build_case(args, ref, tau) for tau in [float(x) for x in args.taus.split(",") if x.strip()]]
    scale = ref["scale"]
    dt = float(args.dt_old) * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    degrees = [tuple(int(part) for part in item.split(":")) for item in args.degrees.split(",") if item.strip()]
    fits = [
        fit_universal_full(cases, args.fit_region, degree_r, degree_p, float(args.target_floor_frac), dt, dx, dz)
        for degree_r, degree_p in degrees
    ]

    report_fits: list[dict[str, object]] = []
    for fit in fits:
        entry: dict[str, object] = {k: v for k, v in fit.items() if k not in {"coeff", "column_scale"}}
        entry["coeff_abs"] = stats_abs(np.asarray(fit["coeff"]))
        entry["column_scale_abs"] = stats_abs(np.asarray(fit["column_scale"]))
        entry["by_case"] = {}
        for case in cases:
            entry["by_case"][f"tau={case.tau_old:g}"] = {
                region: evaluate_case_full(case, fit, region, dt, dx, dz)
                for region in ["trusted", "core1", "core10"]
            }
        report_fits.append(entry)

    best_fit = fits[-1]
    residual_plot = out / "universal_f_of_R_I2_full_residual_maps.png"
    render_residual_plot(residual_plot, cases, best_fit, args.fit_region, dt, dx, dz)

    np.savez_compressed(
        out / "universal_f_of_R_I2_full_fit_coefficients.npz",
        **{f"degree_{fit['degree_r']}_{fit['degree_p']}_coeff": np.asarray(fit["coeff"]) for fit in fits},
        **{f"degree_{fit['degree_r']}_{fit['degree_p']}_column_scale": np.asarray(fit["column_scale"]) for fit in fits},
        degrees=np.asarray(degrees),
    )
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(dt),
            "taus": [case.tau_old for case in cases],
            "fit_region": args.fit_region,
            "degrees": [list(d) for d in degrees],
            "mp": float(args.mp),
            "ell": float(args.ell),
            "ell_over_planck_length": float(args.ell) / PLANCK_LENGTH_EV_INV,
            "target_floor_frac": float(args.target_floor_frac),
        },
        "definition": {
            "f_model": "f(R,I2)=sum c_nm T_n(t(asinh(R/R0))) T_m(t(asinh(I2/I20))). f_R and f_I2 are analytic derivatives of the same fitted function.",
            "I2": "I2=R_mn R^mn, with contractions using gtilde inverse.",
            "full_metric_lhs": "f_R R_mn - 1/2 f g_mn + (g_mn Box - nabla_m nabla_n)f_R + 2 f_I2 R_m^a R_an + Box(f_I2 R_mn) + g_mn nabla^a nabla^b(f_I2 R_ab) - 2 nabla_a nabla_b(f_I2 R^a_(m delta^b_n)).",
            "diagnostic_scope": "Fixed-background回代/拟合诊断；尚不是完整D支动力学求解器。",
        },
        "fits": report_fits,
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "coefficients_npz": str((out / "universal_f_of_R_I2_full_fit_coefficients.npz").resolve()),
            "residual_plot": str(residual_plot.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--fit-region", choices=["trusted", "core1", "core10"], default="core10")
    parser.add_argument("--degrees", type=str, default="2:2,3:3,4:4")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
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
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
