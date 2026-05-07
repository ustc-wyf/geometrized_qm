from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_bcd_residuals_from_a_reference import stress_tensor_tilde
from coordinate_matter_evolution import (
    coordinate_matter_rhs,
    conservative_density_from_rho_u,
    coordinate_matter_rhs_covector,
    inverse_metric_block,
)
from mixed_tilde_initial_data import (
    bohm_fields_from_wave_derivatives,
    exact_localized_wave_derivatives,
    localized_direct_tilde_coordinate_initial,
    localized_direct_tilde_coordinate_snapshot,
)
from prototype_d_tensor_reconstructed_interface import least_squares_tensor_rows
from simulate_d_reduced_dynamic_same_initial import (
    add_lines,
    interface_jump_diagnostics,
    metric_curvature_from_history,
    relative_l1,
    rk4_matter_step,
    sanitize_metric_tail,
)
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    phase_field,
    spectral_omega,
)
from physical_units import (
    HBAR_C_EV_M,
    PLANCK_LENGTH_EV_INV,
    PLANCK_MASS_EV,
    ev_inv_to_fs,
    optical_crossing_params_from_physical_scale,
)


def load_pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def real_array(value: np.ndarray) -> np.ndarray:
    return np.asarray(np.real_if_close(value), dtype=float)


def tensor_norm(tensor: np.ndarray) -> np.ndarray:
    return np.sqrt(np.einsum("...ab,...ab->...", tensor, tensor, optimize=True))


def scalar_stats(field: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.asarray(field, dtype=float)[mask]
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    abs_vals = np.abs(vals)
    return {
        "mean": float(np.mean(abs_vals)),
        "p50": float(np.percentile(abs_vals, 50.0)),
        "p95": float(np.percentile(abs_vals, 95.0)),
        "max": float(np.max(abs_vals)),
    }


def list_scalar_stats(values: list[float] | np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0}
    abs_vals = np.abs(vals)
    return {
        "count": int(vals.size),
        "abs_max": float(np.max(abs_vals)),
        "abs_mean": float(np.mean(abs_vals)),
        "abs_median": float(np.median(abs_vals)),
        "p95": float(np.percentile(abs_vals, 95.0)),
    }


def smoothstep01(value: np.ndarray) -> np.ndarray:
    x = np.clip(np.asarray(value, dtype=float), 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def log_taper_weight(field: np.ndarray, high_threshold: float, decades: float, floor: float) -> np.ndarray:
    """Smooth diagnostic weight that reaches one at the binary support threshold.

    This is deliberately diagnostic-only.  It does not change the evolved
    equations or sources; it just replaces a sharp support integral by a
    continuous low-density edge integral for convergence checks.
    """
    high = max(float(high_threshold), floor)
    if decades <= 0.0:
        return (np.asarray(field, dtype=float) >= high).astype(float)
    low = max(high * (10.0 ** (-float(decades))), floor)
    denom = max(np.log(high) - np.log(low), 1.0e-300)
    coord = (np.log(np.maximum(np.asarray(field, dtype=float), floor)) - np.log(low)) / denom
    return smoothstep01(coord)


def support_taper_weights(
    rho: np.ndarray,
    measure: np.ndarray,
    rho_frac: float,
    measure_frac: float,
    decades: float,
    floor: float,
) -> np.ndarray:
    rho_threshold = float(rho_frac) * float(np.max(rho))
    measure_threshold = float(measure_frac) * float(np.max(measure))
    rho_weight = log_taper_weight(rho, rho_threshold, decades, floor)
    measure_weight = log_taper_weight(measure, measure_threshold, decades, floor)
    return np.minimum(rho_weight, measure_weight)


def weighted_relative_l1(a: np.ndarray, b: np.ndarray, weights: np.ndarray) -> float:
    w = np.asarray(weights, dtype=float)
    valid = np.isfinite(a) & np.isfinite(b) & np.isfinite(w) & (w > 0.0)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(w[valid] * np.abs(a[valid]))), 1.0e-300)
    return float(np.sum(w[valid] * np.abs(a[valid] - b[valid])) / denom)


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    vals = np.asarray(values, dtype=float).ravel()
    w = np.asarray(weights, dtype=float).ravel()
    valid = np.isfinite(vals) & np.isfinite(w) & (w > 0.0)
    if not np.any(valid):
        return 0.0
    vals = vals[valid]
    w = w[valid]
    order = np.argsort(vals)
    vals = vals[order]
    w = w[order]
    cdf = np.cumsum(w)
    target = float(np.clip(quantile, 0.0, 1.0)) * float(cdf[-1])
    return float(vals[min(int(np.searchsorted(cdf, target, side="left")), vals.size - 1)])


def point_segment_distance_squared(px: np.ndarray, pz: np.ndarray, seg: dict[str, float]) -> np.ndarray:
    vx = seg["x1"] - seg["x0"]
    vz = seg["z1"] - seg["z0"]
    denom = max(vx * vx + vz * vz, 1.0e-300)
    tau = ((px - seg["x0"]) * vx + (pz - seg["z0"]) * vz) / denom
    tau = np.clip(tau, 0.0, 1.0)
    cx = seg["x0"] + tau * vx
    cz = seg["z0"] + tau * vz
    return (px - cx) ** 2 + (pz - cz) ** 2


def tensor_interface_anchor_mask(rows: list[dict[str, float]], x: np.ndarray, z: np.ndarray, band_radius: float) -> np.ndarray:
    if band_radius <= 0.0 or not rows:
        return np.zeros((len(x), len(z)), dtype=bool)
    xg, zg = np.meshgrid(x, z, indexing="ij")
    mask = np.zeros_like(xg, dtype=bool)
    r2 = float(band_radius) * float(band_radius)
    for row in rows:
        if row.get("direct_tensor_solved", 0.0) <= 0.5:
            continue
        mask |= point_segment_distance_squared(xg, zg, row) <= r2
    return mask


def full_tensor_interface_diagnostics(
    *,
    raw_y: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    ricci: np.ndarray,
    rho: np.ndarray,
    rhs_tensor: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    ell: float,
    half_width_y: float,
    samples: int,
    tensor_rel_tol: float,
    max_iter: int,
    multistart: bool,
    max_seeds: int,
    band_factor: float,
) -> dict[str, object]:
    rows: list[dict[str, float]] = []
    for level in (-1.0, 1.0):
        rows.extend(
            least_squares_tensor_rows(
                raw_y=raw_y,
                metric_cov=metric_cov,
                metric_inv=metric_inv,
                ricci=ricci,
                rho=rho,
                rhs_tensor=rhs_tensor,
                x=x,
                z=z,
                level=level,
                ell=ell,
                half_width_y=half_width_y,
                samples=samples,
                tensor_rel_tol=tensor_rel_tol,
                max_iter=max_iter,
                multistart=multistart,
                max_seeds=max_seeds,
            )
        )

    solved = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    anchor_mask = tensor_interface_anchor_mask(rows=solved, x=x, z=z, band_radius=band_factor * float(np.hypot(dx, dz)))
    return {
        "rows": rows,
        "accepted_rows": solved,
        "segment_count": len(rows),
        "candidate_count": int(sum(row.get("candidate", 0.0) > 0.5 for row in rows)),
        "solved_count": len(solved),
        "solved_fraction": float(len(solved) / max(len(rows), 1)),
        "anchor_count": int(np.count_nonzero(anchor_mask)),
        "anchor_fraction": float(np.count_nonzero(anchor_mask) / max(anchor_mask.size, 1)),
        "direct_tensor_residual_relative": list_scalar_stats([row.get("direct_tensor_residual_relative", np.nan) for row in solved]),
        "tensor_anchor_mask": anchor_mask,
    }


def build_crossing_params(args: argparse.Namespace) -> tuple[FlatLocalizedCrossingParams, dict[str, float] | None]:
    if args.physical_optical:
        params, scale = optical_crossing_params_from_physical_scale(
            wavelength_nm=args.wavelength_nm,
            mass_over_omega=args.mass_over_omega,
            resolution=args.resolution,
            alpha=args.alpha,
            phi0=args.phi0,
            rho_floor=args.rho_floor,
        )
        return params, scale.to_json()
    return (
        FlatLocalizedCrossingParams(
            nx=args.resolution,
            nz=args.resolution,
            alpha=args.alpha,
            phi0=args.phi0,
            rho_floor=args.rho_floor,
        ),
        None,
    )


def scale_interface_rows(rows: list[dict[str, float]], factor: float) -> list[dict[str, float]]:
    scaled: list[dict[str, float]] = []
    for row in rows:
        item = dict(row)
        for key in ("x0", "x1", "z0", "z1"):
            if key in item:
                item[key] = float(item[key]) * factor
        scaled.append(item)
    return scaled


def evolve_spectral_state(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    initial_time: float,
) -> tuple[np.ndarray, np.ndarray]:
    if initial_time == 0.0:
        return psi0, psi0_hat
    psi_hat_initial = psi0_hat * np.exp(-1j * omega * float(initial_time))
    psi_initial = np.fft.ifft2(psi_hat_initial)
    return psi_initial, psi_hat_initial


def shift_mask_no_wrap(mask: np.ndarray, di: int, dj: int) -> np.ndarray:
    shifted = np.zeros_like(mask, dtype=bool)
    nx, nz = mask.shape
    src_i0 = max(0, -di)
    src_i1 = min(nx, nx - di)
    dst_i0 = max(0, di)
    dst_i1 = min(nx, nx + di)
    src_j0 = max(0, -dj)
    src_j1 = min(nz, nz - dj)
    dst_j0 = max(0, dj)
    dst_j1 = min(nz, nz + dj)
    if src_i0 < src_i1 and src_j0 < src_j1:
        shifted[dst_i0:dst_i1, dst_j0:dst_j1] = mask[src_i0:src_i1, src_j0:src_j1]
    return shifted


def erode_mask_8(mask: np.ndarray, iterations: int) -> np.ndarray:
    eroded = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(iterations, 0)):
        next_mask = eroded.copy()
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                next_mask &= shift_mask_no_wrap(eroded, di, dj)
        eroded = next_mask
    return eroded


def dilate_mask_4(mask: np.ndarray, iterations: int) -> np.ndarray:
    dilated = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(iterations, 0)):
        dilated = (
            dilated
            | shift_mask_no_wrap(dilated, 1, 0)
            | shift_mask_no_wrap(dilated, -1, 0)
            | shift_mask_no_wrap(dilated, 0, 1)
            | shift_mask_no_wrap(dilated, 0, -1)
        )
    return dilated


def kg_positive_frequency_norm(psi: np.ndarray, psi_t: np.ndarray, dx: float, dz: float) -> float:
    """Positive-frequency KG norm for the (+--) convention used here.

    For a mode exp(-i omega t), this is integral 2 omega |psi|^2 dx dz.
    """
    density = np.real(1j * (np.conj(psi) * psi_t - psi * np.conj(psi_t)))
    return float(np.sum(density) * dx * dz)


def safe_sqrt_abs_det(metric_cov: np.ndarray, floor: float = 1.0e-300) -> np.ndarray:
    det = np.linalg.det(metric_cov.reshape(-1, 3, 3)).reshape(metric_cov.shape[:2])
    return np.sqrt(np.maximum(np.abs(det), floor))


def tilde_measure_from_bohm(bohm: dict[str, np.ndarray], mass: float) -> np.ndarray:
    """Positive transformed matter measure density sqrt(|gtilde|) rho_tilde.

    We use the positive branch discussed in the notes:
        sqrt(|gtilde|) rho_tilde = |X| rho / m^2
    for flat original g.  The sign of X is still saved in diagnostics.
    """
    return np.abs(real_array(bohm["X"])) * real_array(bohm["rho"]) / (mass * mass)


def tilde_density_from_measure(
    ntilde_measure: np.ndarray,
    metric_cov: np.ndarray,
    rho_floor: float,
) -> np.ndarray:
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    return np.maximum(ntilde_measure / sqrt_abs_g, rho_floor)


def current_from_measure_and_u(
    measure_density: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov: np.ndarray,
    mass: float,
    u_t_reference: np.ndarray | None,
    active_mask: np.ndarray,
    rho_floor: float,
    metric_inv: np.ndarray | None = None,
    det_cov: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Build lab-slice diagnostics from positive transformed measure density.

    Unlike ``coordinate_matter_rhs_covector`` this does not recover rho by
    dividing by ``j^t``.  It is therefore suitable for output diagnostics when
    the actual update is performed in a better local time chart.
    """
    if metric_inv is None or det_cov is None:
        metric_inv, det_cov = inverse_metric_block(metric_cov)
    sqrt_abs_g = np.sqrt(np.maximum(np.abs(det_cov), 1.0e-300))
    rho = np.where(active_mask, np.maximum(measure_density / sqrt_abs_g, rho_floor), 0.0)
    current = conservative_density_from_rho_u(
        rho=rho,
        u_x=np.where(active_mask, u_x, 0.0),
        u_z=np.where(active_mask, u_z, 0.0),
        metric_cov_txz=metric_cov,
        mass=mass,
        branch="negative_frequency",
        u_t_reference=u_t_reference,
        rho_floor=rho_floor,
        metric_inv_txz=metric_inv,
        det_cov_txz=det_cov,
    )
    for key in ("u_t", "u_x", "u_z", "discriminant", "j_t", "j_x", "j_z", "n_cons", "flux_x", "flux_z", "sqrt_abs_g"):
        if key in current:
            current[key] = np.where(active_mask, current[key], 0.0)
    current["rho"] = rho
    current["u_x"] = np.where(active_mask, u_x, 0.0)
    current["u_z"] = np.where(active_mask, u_z, 0.0)
    return np.where(active_mask, current["n_cons"], 0.0), current


def flat_metric_cov_like(mask: np.ndarray) -> np.ndarray:
    eta = np.zeros(mask.shape + (3, 3), dtype=float)
    eta[..., 0, 0] = 1.0
    eta[..., 1, 1] = -1.0
    eta[..., 2, 2] = -1.0
    return eta


def unwrap_phase_like(candidate: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Unwrap a phase field and choose the 2*pi branch closest to a reference."""
    unwrapped = np.unwrap(np.unwrap(np.asarray(candidate, dtype=float), axis=0), axis=1)
    ref = np.asarray(reference, dtype=float)
    offset = np.round(np.median((ref - unwrapped) / (2.0 * np.pi)))
    return unwrapped + 2.0 * np.pi * offset


def gradient_from_phase(s_phase: np.ndarray, active_mask: np.ndarray, dx: float, dz: float) -> tuple[np.ndarray, np.ndarray]:
    sx = (np.roll(s_phase, -1, axis=0) - np.roll(s_phase, 1, axis=0)) / (2.0 * dx)
    sz = (np.roll(s_phase, -1, axis=1) - np.roll(s_phase, 1, axis=1)) / (2.0 * dz)
    return np.where(active_mask, sx, 0.0), np.where(active_mask, sz, 0.0)


def phase_current_from_state(
    n_cons: np.ndarray,
    s_phase: np.ndarray,
    metric_cov: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    active_mask: np.ndarray,
    s_t_reference: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    current = coordinate_matter_rhs(
        n_cons=np.where(active_mask, n_cons, 0.0),
        s=s_phase,
        metric_cov_txz=metric_cov,
        mass=mass,
        dx=dx,
        dz=dz,
        phase_branch="negative_frequency",
        s_t_reference=s_t_reference,
    )
    ux, uz = gradient_from_phase(s_phase, active_mask, dx, dz)
    for key in ("rho", "s_t", "discriminant", "j_t", "j_x", "j_z", "flux_x", "flux_z"):
        if key in current:
            current[key] = np.where(active_mask, current[key], 0.0)
    current["u_t"] = current["s_t"]
    current["u_x"] = ux
    current["u_z"] = uz
    return np.where(active_mask, n_cons, 0.0), ux, uz, current


def local_time_covariant_metric(metric_cov: np.ndarray, tau_a: float, tau_b: float) -> np.ndarray:
    """Covariant metric in the local time chart tau=t+a x+b z.

    This is the exact tensor transform for the linear chart change, but it is
    used here as a numerical preconditioner for the matter update.
    """
    a = float(tau_a)
    b = float(tau_b)
    g = np.asarray(metric_cov, dtype=float)
    out = np.empty_like(g)
    gtt = g[..., 0, 0]
    gtx = g[..., 0, 1]
    gtz = g[..., 0, 2]
    gxx = g[..., 1, 1]
    gxz = g[..., 1, 2]
    gzz = g[..., 2, 2]
    out[..., 0, 0] = gtt
    out[..., 0, 1] = gtx - a * gtt
    out[..., 1, 0] = out[..., 0, 1]
    out[..., 0, 2] = gtz - b * gtt
    out[..., 2, 0] = out[..., 0, 2]
    out[..., 1, 1] = gxx - 2.0 * a * gtx + a * a * gtt
    out[..., 1, 2] = gxz - a * gtz - b * gtx + a * b * gtt
    out[..., 2, 1] = out[..., 1, 2]
    out[..., 2, 2] = gzz - 2.0 * b * gtz + b * b * gtt
    return out


def local_time_contravariant_metric(metric_inv: np.ndarray, tau_a: float, tau_b: float) -> np.ndarray:
    """Contravariant metric in the local chart tau=t+a x+b z.

    The chart has unit Jacobian determinant, so the covariant determinant is
    unchanged.  Using this direct inverse-metric transform avoids inverting the
    same patch metric repeatedly during local-time RK stages.
    """
    a = float(tau_a)
    b = float(tau_b)
    g = np.asarray(metric_inv, dtype=float)
    out = np.empty_like(g)
    gtt = g[..., 0, 0]
    gtx = g[..., 0, 1]
    gtz = g[..., 0, 2]
    gxx = g[..., 1, 1]
    gxz = g[..., 1, 2]
    gzz = g[..., 2, 2]
    out[..., 0, 0] = gtt + 2.0 * a * gtx + 2.0 * b * gtz + a * a * gxx + 2.0 * a * b * gxz + b * b * gzz
    out[..., 0, 1] = gtx + a * gxx + b * gxz
    out[..., 1, 0] = out[..., 0, 1]
    out[..., 0, 2] = gtz + a * gxz + b * gzz
    out[..., 2, 0] = out[..., 0, 2]
    out[..., 1, 1] = gxx
    out[..., 1, 2] = gxz
    out[..., 2, 1] = gxz
    out[..., 2, 2] = gzz
    return out


def local_time_covector_from_lab(u_t: np.ndarray, u_x: np.ndarray, u_z: np.ndarray, tau_a: float, tau_b: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_tau = np.asarray(u_t, dtype=float)
    return u_tau, np.asarray(u_x, dtype=float) - float(tau_a) * u_tau, np.asarray(u_z, dtype=float) - float(tau_b) * u_tau


def lab_covector_from_local(u_tau: np.ndarray, u_x_tau: np.ndarray, u_z_tau: np.ndarray, tau_a: float, tau_b: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u_tau = np.asarray(u_tau, dtype=float)
    return u_tau, np.asarray(u_x_tau, dtype=float) + float(tau_a) * u_tau, np.asarray(u_z_tau, dtype=float) + float(tau_b) * u_tau


def local_time_candidate_arrays(max_tilt: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    if max_tilt < 0.0:
        raise ValueError("--local-time-max-tilt must be nonnegative")
    if step <= 0.0:
        raise ValueError("--local-time-step must be positive")
    values = np.arange(-float(max_tilt), float(max_tilt) + 0.5 * float(step), float(step))
    if values.size == 0:
        values = np.array([0.0], dtype=float)
    return values, values.copy()


def local_time_stationary_candidates(
    *,
    metric_inv: np.ndarray,
    active: np.ndarray,
    norm_floor: float,
    max_abs_tilt: float,
    quantization: float,
) -> tuple[list[tuple[float, float]], dict[str, float]]:
    """Analytic tau=t+a x+b z candidates where tau is most timelike.

    For a fixed cell the inverse-metric norm of d tau is a quadratic form in
    (a,b).  If its spatial Hessian is negative definite, the stationary point
    is the local maximum.  Adding those points to the coarse grid prevents a
    legitimate but narrow time cone from being skipped by the coarse scan.
    """
    gtt = metric_inv[..., 0, 0]
    gtx = metric_inv[..., 0, 1]
    gtz = metric_inv[..., 0, 2]
    gxx = metric_inv[..., 1, 1]
    gxz = metric_inv[..., 1, 2]
    gzz = metric_inv[..., 2, 2]

    tr = gxx + gzz
    det = gxx * gzz - gxz * gxz
    eig_disc = np.maximum((gxx - gzz) ** 2 + 4.0 * gxz * gxz, 0.0)
    max_eig = 0.5 * (tr + np.sqrt(eig_disc))
    valid = np.asarray(active, dtype=bool) & (max_eig < 0.0) & (np.abs(det) > 1.0e-300)

    h0 = (gzz * gtx - gxz * gtz) / det
    h1 = (-gxz * gtx + gxx * gtz) / det
    a_star = -h0
    b_star = -h1
    q_max = gtt - (gtx * h0 + gtz * h1)
    valid &= (q_max > norm_floor) & np.isfinite(a_star) & np.isfinite(b_star)
    valid &= (np.abs(a_star) <= max_abs_tilt) & (np.abs(b_star) <= max_abs_tilt)

    candidates: list[tuple[float, float]] = []
    if np.any(valid):
        if quantization > 0.0:
            a_vals = np.round(a_star[valid] / quantization) * quantization
            b_vals = np.round(b_star[valid] / quantization) * quantization
        else:
            a_vals = a_star[valid]
            b_vals = b_star[valid]
        for a, b in zip(a_vals, b_vals):
            candidates.append((round(float(a), 8), round(float(b), 8)))
    unique = list(dict.fromkeys(candidates))
    diag = {
        "stationary_candidate_raw_count": float(np.count_nonzero(valid)),
        "stationary_candidate_unique_count": float(len(unique)),
    }
    return unique, diag


def select_local_time_chart(
    *,
    current: dict[str, np.ndarray],
    metric_inv: np.ndarray,
    measure: np.ndarray,
    trusted: np.ndarray,
    a_values: np.ndarray,
    b_values: np.ndarray,
    norm_floor: float,
    admissible_fraction_floor: float,
) -> tuple[float, float, dict[str, float]]:
    jt = np.asarray(current["j_t"], dtype=float)
    jx = np.asarray(current["j_x"], dtype=float)
    jz = np.asarray(current["j_z"], dtype=float)
    measure = np.asarray(measure, dtype=float)
    trusted = np.asarray(trusted, dtype=bool)

    gtt = metric_inv[..., 0, 0]
    gtx = metric_inv[..., 0, 1]
    gtz = metric_inv[..., 0, 2]
    gxx = metric_inv[..., 1, 1]
    gxz = metric_inv[..., 1, 2]
    gzz = metric_inv[..., 2, 2]
    trusted_measure = max(float(np.sum(measure[trusted])), 1.0e-300)

    def evaluate(a: float, b: float) -> dict[str, float] | None:
        tau_norm = gtt + 2.0 * a * gtx + 2.0 * b * gtz + (a * a) * gxx + 2.0 * a * b * gxz + (b * b) * gzz
        good = trusted & (tau_norm > norm_floor)
        good_measure = float(np.sum(measure[good]))
        if good_measure <= 0.0:
            return None
        admissible_fraction = good_measure / trusted_measure
        if admissible_fraction <= 0.0:
            return None
        j_tau = jt + a * jx + b * jz
        abs_j = np.abs(j_tau[good])
        score = float(np.sum(measure[good] * abs_j) / max(good_measure, 1.0e-300))
        score *= float(admissible_fraction)
        return {
            "a": float(a),
            "b": float(b),
            "score": score,
            "admissible_fraction": float(admissible_fraction),
            "abs_jtau_mean": float(np.mean(abs_j)),
            "abs_jtau_p05": float(np.percentile(abs_j, 5.0)),
            "abs_jtau_p50": float(np.percentile(abs_j, 50.0)),
            "abs_jtau_p95": float(np.percentile(abs_j, 95.0)),
            "tau_norm_min": float(np.nanmin(tau_norm[trusted])),
            "tau_norm_p05": float(np.percentile(tau_norm[trusted], 5.0)),
        }

    best: dict[str, float] | None = None
    for a in a_values:
        for b in b_values:
            item = evaluate(float(a), float(b))
            if item is None:
                continue
            if best is None or item["score"] > best["score"]:
                best = item

    if best is None:
        best = {
            "a": 0.0,
            "b": 0.0,
            "score": 0.0,
            "admissible_fraction": 0.0,
            "abs_jtau_mean": 0.0,
            "abs_jtau_p05": 0.0,
            "abs_jtau_p50": 0.0,
            "abs_jtau_p95": 0.0,
            "tau_norm_min": 0.0,
            "tau_norm_p05": 0.0,
        }

    if best["admissible_fraction"] < admissible_fraction_floor:
        fallback = evaluate(0.0, 0.0)
        if fallback is not None and fallback["score"] >= best["score"]:
            best = fallback
        best["fallback_to_lab"] = 1.0
        return 0.0, 0.0, best

    refine_step = max(float(a_values[1] - a_values[0]) if a_values.size > 1 else 0.25, 0.25) / 2.0
    a_ref = np.arange(max(-float(np.max(np.abs(a_values))), best["a"] - refine_step), min(float(np.max(np.abs(a_values))), best["a"] + refine_step) + 0.5 * refine_step, refine_step)
    b_ref = np.arange(max(-float(np.max(np.abs(b_values))), best["b"] - refine_step), min(float(np.max(np.abs(b_values))), best["b"] + refine_step) + 0.5 * refine_step, refine_step)
    for a in a_ref:
        for b in b_ref:
            item = evaluate(float(a), float(b))
            if item is None:
                continue
            if item["score"] > best["score"]:
                best = item

    best["fallback_to_lab"] = 0.0
    return float(best["a"]), float(best["b"]), best


def select_local_time_patch_labels(
    *,
    current: dict[str, np.ndarray],
    metric_inv: np.ndarray,
    measure: np.ndarray,
    active: np.ndarray,
    trusted: np.ndarray,
    a_values: np.ndarray,
    b_values: np.ndarray,
    norm_floor: float,
    include_stationary_candidates: bool = False,
    stationary_candidate_quantization: float = 0.05,
    tile_size: int = 1,
) -> tuple[np.ndarray, list[tuple[float, float]], dict[str, float]]:
    """Choose a local time chart per cell.

    This is a numerical chart atlas for the matter update.  It does not add a
    force or damping term; it only avoids recovering rho from a bad lab-time
    projection when another local time covector is non-characteristic.
    """
    candidates = [(float(a), float(b)) for a in a_values for b in b_values]
    stationary_diag = {
        "stationary_candidate_raw_count": 0.0,
        "stationary_candidate_unique_count": 0.0,
    }
    if include_stationary_candidates:
        max_abs_tilt = max(
            float(np.max(np.abs(a_values))) if a_values.size else 0.0,
            float(np.max(np.abs(b_values))) if b_values.size else 0.0,
        )
        extra, stationary_diag = local_time_stationary_candidates(
            metric_inv=metric_inv,
            active=active,
            norm_floor=norm_floor,
            max_abs_tilt=max_abs_tilt,
            quantization=stationary_candidate_quantization,
        )
        candidates.extend(extra)
        candidates = list(dict.fromkeys(candidates))
    if (0.0, 0.0) not in candidates:
        candidates.insert(0, (0.0, 0.0))
    jt = np.asarray(current["j_t"], dtype=float)
    jx = np.asarray(current["j_x"], dtype=float)
    jz = np.asarray(current["j_z"], dtype=float)
    gtt = metric_inv[..., 0, 0]
    gtx = metric_inv[..., 0, 1]
    gtz = metric_inv[..., 0, 2]
    gxx = metric_inv[..., 1, 1]
    gxz = metric_inv[..., 1, 2]
    gzz = metric_inv[..., 2, 2]

    score = np.full(active.shape + (len(candidates),), -np.inf, dtype=float)
    tau_norm_stack = np.full_like(score, np.nan)
    jtau_stack = np.zeros_like(score)
    for idx, (a, b) in enumerate(candidates):
        tau_norm = gtt + 2.0 * a * gtx + 2.0 * b * gtz + (a * a) * gxx + 2.0 * a * b * gxz + (b * b) * gzz
        j_tau = jt + a * jx + b * jz
        admissible = active & (tau_norm > norm_floor)
        tau_norm_stack[..., idx] = tau_norm
        jtau_stack[..., idx] = j_tau
        score[..., idx] = np.where(admissible, np.abs(j_tau), -np.inf)

    lab_idx = candidates.index((0.0, 0.0))
    tile = max(int(tile_size), 1)
    if tile > 1:
        label = np.zeros(active.shape, dtype=int)
        pool_label = np.full(active.shape, lab_idx, dtype=int)
        patch_candidates: list[tuple[float, float]] = []
        measure_arr = np.asarray(measure, dtype=float)
        nx, nz = active.shape
        for i0 in range(0, nx, tile):
            for j0 in range(0, nz, tile):
                sl = (slice(i0, min(i0 + tile, nx)), slice(j0, min(j0 + tile, nz)))
                tile_active = active[sl]
                if not np.any(tile_active):
                    continue
                scope = trusted[sl] & tile_active
                if not np.any(scope):
                    scope = tile_active
                weights = np.where(scope, measure_arr[sl], 0.0)
                total_weight = max(float(np.sum(weights)), 1.0e-300)
                tile_score = score[sl]
                finite = np.isfinite(tile_score) & scope[..., None]
                good_weight = np.sum(np.where(finite, weights[..., None], 0.0), axis=(0, 1))
                coverage = good_weight / total_weight
                safe_tile_score = np.where(finite, tile_score, 0.0)
                weighted_abs = np.sum(weights[..., None] * safe_tile_score, axis=(0, 1))
                mean_abs = weighted_abs / np.maximum(good_weight, 1.0e-300)
                objective = coverage + 1.0e-6 * mean_abs
                if not np.any(np.isfinite(objective)) or float(np.max(objective)) <= 0.0:
                    chosen = lab_idx
                else:
                    chosen = int(np.argmax(objective))
                patch_label = len(patch_candidates)
                patch_candidates.append(candidates[chosen])
                label_view = label[sl]
                label_view[tile_active] = patch_label
                pool_view = pool_label[sl]
                pool_view[tile_active] = chosen
        if not patch_candidates:
            patch_candidates = [(0.0, 0.0)]
        candidates_out = patch_candidates
        selected_pool = pool_label
    else:
        label = np.argmax(score, axis=-1).astype(int)
        selected_pool = label.copy()
        best = np.take_along_axis(score, selected_pool[..., None], axis=-1)[..., 0]
        no_candidate = active & (~np.isfinite(best))
        label[no_candidate] = lab_idx
        selected_pool[no_candidate] = lab_idx
        candidates_out = candidates
    best = np.take_along_axis(score, selected_pool[..., None], axis=-1)[..., 0]
    best_jtau = np.take_along_axis(np.abs(jtau_stack), selected_pool[..., None], axis=-1)[..., 0]
    best_norm = np.take_along_axis(tau_norm_stack, selected_pool[..., None], axis=-1)[..., 0]

    total_measure = max(float(np.sum(np.asarray(measure, dtype=float)[trusted])), 1.0e-300)
    uncovered = trusted & (~np.isfinite(best))
    used_labels = sorted(int(x) for x in np.unique(label[active]))
    label_is_lab = np.asarray([candidate == (0.0, 0.0) for candidate in candidates_out], dtype=bool)
    label_edges = 0
    if np.any(active):
        x_edges = active[:-1, :] & active[1:, :] & (label[:-1, :] != label[1:, :])
        z_edges = active[:, :-1] & active[:, 1:] & (label[:, :-1] != label[:, 1:])
        label_edges = int(np.count_nonzero(x_edges) + np.count_nonzero(z_edges))
    diag = {
        "patch_count": float(len(used_labels)),
        "candidate_count": float(len(candidates)),
        "tile_size": float(tile),
        **stationary_diag,
        "patch_label_edges": float(label_edges),
        "uncovered_trusted_measure_fraction": float(np.sum(np.asarray(measure, dtype=float)[uncovered]) / total_measure),
        "abs_jtau_p05": float(np.percentile(best_jtau[trusted], 5.0)) if np.any(trusted) else 0.0,
        "abs_jtau_p50": float(np.percentile(best_jtau[trusted], 50.0)) if np.any(trusted) else 0.0,
        "abs_jtau_p95": float(np.percentile(best_jtau[trusted], 95.0)) if np.any(trusted) else 0.0,
        "tau_norm_min": float(np.nanmin(best_norm[trusted])) if np.any(trusted) else 0.0,
        "tau_norm_p05": float(np.percentile(best_norm[trusted], 5.0)) if np.any(trusted) else 0.0,
        "lab_label_fraction": float(np.count_nonzero(active & label_is_lab[label]) / max(int(np.count_nonzero(active)), 1)),
    }
    return label, candidates_out, diag


def local_time_label_diagnostics(
    *,
    current: dict[str, np.ndarray],
    metric_inv: np.ndarray,
    measure: np.ndarray,
    active: np.ndarray,
    trusted: np.ndarray,
    labels: np.ndarray,
    candidates: list[tuple[float, float]],
    norm_floor: float,
) -> dict[str, float]:
    """Re-evaluate a cached local-time atlas on the current fields."""
    jt = np.asarray(current["j_t"], dtype=float)
    jx = np.asarray(current["j_x"], dtype=float)
    jz = np.asarray(current["j_z"], dtype=float)
    gtt = metric_inv[..., 0, 0]
    gtx = metric_inv[..., 0, 1]
    gtz = metric_inv[..., 0, 2]
    gxx = metric_inv[..., 1, 1]
    gxz = metric_inv[..., 1, 2]
    gzz = metric_inv[..., 2, 2]
    selected_tau_norm = np.full(active.shape, np.nan, dtype=float)
    selected_abs_jtau = np.full(active.shape, np.nan, dtype=float)
    label_arr = np.asarray(labels, dtype=int)
    used_labels = sorted(int(x) for x in np.unique(label_arr[active]))
    for label_value in used_labels:
        a, b = candidates[label_value]
        mask = active & (label_arr == label_value)
        if not np.any(mask):
            continue
        tau_norm = gtt + 2.0 * a * gtx + 2.0 * b * gtz + (a * a) * gxx + 2.0 * a * b * gxz + (b * b) * gzz
        j_tau = jt + a * jx + b * jz
        selected_tau_norm[mask] = tau_norm[mask]
        selected_abs_jtau[mask] = np.abs(j_tau[mask])
    total_measure = max(float(np.sum(np.asarray(measure, dtype=float)[trusted])), 1.0e-300)
    covered = trusted & np.isfinite(selected_tau_norm) & (selected_tau_norm > norm_floor)
    uncovered = trusted & (~covered)
    label_edges = 0
    if np.any(active):
        x_edges = active[:-1, :] & active[1:, :] & (label_arr[:-1, :] != label_arr[1:, :])
        z_edges = active[:, :-1] & active[:, 1:] & (label_arr[:, :-1] != label_arr[:, 1:])
        label_edges = int(np.count_nonzero(x_edges) + np.count_nonzero(z_edges))
    label_is_lab = np.asarray([candidate == (0.0, 0.0) for candidate in candidates], dtype=bool)
    return {
        "patch_count": float(len(used_labels)),
        "candidate_count": float(len(candidates)),
        "patch_label_edges": float(label_edges),
        "uncovered_trusted_measure_fraction": float(np.sum(np.asarray(measure, dtype=float)[uncovered]) / total_measure),
        "abs_jtau_p05": float(np.nanpercentile(selected_abs_jtau[trusted], 5.0)) if np.any(trusted) else 0.0,
        "abs_jtau_p50": float(np.nanpercentile(selected_abs_jtau[trusted], 50.0)) if np.any(trusted) else 0.0,
        "abs_jtau_p95": float(np.nanpercentile(selected_abs_jtau[trusted], 95.0)) if np.any(trusted) else 0.0,
        "tau_norm_min": float(np.nanmin(selected_tau_norm[trusted])) if np.any(trusted) else 0.0,
        "tau_norm_p05": float(np.nanpercentile(selected_tau_norm[trusted], 5.0)) if np.any(trusted) else 0.0,
        "lab_label_fraction": float(np.count_nonzero(active & label_is_lab[label_arr]) / max(int(np.count_nonzero(active)), 1)),
    }


def rk4_patchwise_local_time_matter_step(
    *,
    measure_density: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    metric_det: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    dt: float,
    u_t_reference: np.ndarray,
    active_mask: np.ndarray,
    labels: np.ndarray,
    candidates: list[tuple[float, float]],
    stencil_mask: np.ndarray | None = None,
    flux_divergence_mode: str = "central",
    rusanov_strength: float = 1.0,
    matching_flux_weight: str = "abs_n",
    use_bounding_boxes: bool = True,
    split_components: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    write_mask = np.asarray(active_mask, dtype=bool)
    if stencil_mask is None:
        stencil_mask_arr = write_mask
    else:
        stencil_mask_arr = np.asarray(stencil_mask, dtype=bool) | write_mask
    measure_acc = np.zeros_like(measure_density, dtype=float)
    n_acc = np.zeros_like(measure_density, dtype=float)
    ux_acc = np.zeros_like(u_x, dtype=float)
    uz_acc = np.zeros_like(u_z, dtype=float)
    ut_ref = np.asarray(u_t_reference, dtype=float)
    total_input = max(float(np.sum(measure_density[write_mask])), 1.0e-300)

    def component_masks(mask: np.ndarray) -> list[np.ndarray]:
        """Split a patch label into 8-connected components for tighter stencils."""
        remaining = np.asarray(mask, dtype=bool).copy()
        comps: list[np.ndarray] = []
        nx, nz = remaining.shape
        while np.any(remaining):
            start = tuple(int(v) for v in np.argwhere(remaining)[0])
            stack = [start]
            remaining[start] = False
            coords: list[tuple[int, int]] = []
            while stack:
                i, j = stack.pop()
                coords.append((i, j))
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        if di == 0 and dj == 0:
                            continue
                        ii = i + di
                        jj = j + dj
                        if 0 <= ii < nx and 0 <= jj < nz and remaining[ii, jj]:
                            remaining[ii, jj] = False
                            stack.append((ii, jj))
            comp = np.zeros_like(mask, dtype=bool)
            ii, jj = np.asarray(coords, dtype=int).T
            comp[ii, jj] = True
            comps.append(comp)
        return comps

    used_count = 0
    component_count = 0
    bbox_area_sum = 0
    full_area_sum = 0
    bbox_fallback_count = 0
    n_tau_input_abs_sum = 0.0
    n_tau_delta_abs_sum = 0.0
    fv_diag_sum_keys = (
        "face_flux_fv_face_x_count",
        "face_flux_fv_face_z_count",
        "face_flux_fv_divergence_sum_write",
        "face_flux_fv_matching_interface_x_count",
        "face_flux_fv_matching_interface_z_count",
        "face_flux_fv_matching_flux_jump_l1",
        "face_flux_fv_matching_flux_abs_sum",
    )
    fv_diag_max_keys = (
        "face_flux_fv_rusanov_max_speed_x",
        "face_flux_fv_rusanov_max_speed_z",
        "face_flux_fv_matching_flux_jump_p95",
    )
    fv_diag_acc: dict[str, float] = {key: 0.0 for key in fv_diag_sum_keys}
    fv_diag_acc.update({key: 0.0 for key in fv_diag_max_keys})

    def accumulate_patch_diag(patch_diag: dict[str, float]) -> None:
        for key in fv_diag_sum_keys:
            fv_diag_acc[key] += float(patch_diag.get(key, 0.0))
        for key in fv_diag_max_keys:
            fv_diag_acc[key] = max(fv_diag_acc[key], float(patch_diag.get(key, 0.0)))

    for label_value in sorted(int(x) for x in np.unique(labels[write_mask])):
        patch = write_mask & (labels == label_value)
        if not np.any(patch):
            continue
        a, b = candidates[label_value]

        pieces = component_masks(patch) if (use_bounding_boxes and split_components) else [patch]
        component_count += len(pieces)
        for patch_piece in pieces:
            use_full_grid = not use_bounding_boxes
            if use_bounding_boxes:
                idx = np.argwhere(patch_piece)
                i0, j0 = idx.min(axis=0)
                i1, j1 = idx.max(axis=0) + 1
                if i0 <= 0 or j0 <= 0 or i1 >= patch.shape[0] or j1 >= patch.shape[1]:
                    use_full_grid = True
                    bbox_fallback_count += 1
                else:
                    sl = (slice(i0 - 1, i1 + 1), slice(j0 - 1, j1 + 1))
                    patch_local = patch_piece[sl]
                    stencil_local = stencil_mask_arr[sl]
                    full_area_sum += int(measure_density.size)
                    bbox_area_sum += int(stencil_local.size)
                    measure_p, n_p, ux_p, uz_p, _current_p, patch_diag = rk4_local_time_matter_step(
                        measure_density=measure_density[sl],
                        u_x=u_x[sl],
                        u_z=u_z[sl],
                        metric_cov=metric_cov[sl],
                        metric_inv=metric_inv[sl],
                        metric_det=metric_det[sl],
                        mass=mass,
                        dx=dx,
                        dz=dz,
                        dt=dt,
                        u_t_reference=ut_ref[sl],
                        active_mask=stencil_local,
                        write_mask=patch_local,
                        tau_a=a,
                        tau_b=b,
                        flux_divergence_mode=flux_divergence_mode,
                        rusanov_strength=rusanov_strength,
                        matching_flux_weight=matching_flux_weight,
                        matching_write_mask=write_mask[sl],
                    )
                    accumulate_patch_diag(patch_diag)
                    n_tau_input_abs_sum += abs(float(patch_diag.get("n_tau_input_sum", 0.0)))
                    n_tau_delta_abs_sum += abs(float(patch_diag.get("n_tau_output_sum", 0.0)) - float(patch_diag.get("n_tau_input_sum", 0.0)))
                    view = measure_acc[sl]
                    view[patch_local] = measure_p[patch_local]
                    view = n_acc[sl]
                    view[patch_local] = n_p[patch_local]
                    view = ux_acc[sl]
                    view[patch_local] = ux_p[patch_local]
                    view = uz_acc[sl]
                    view[patch_local] = uz_p[patch_local]

            if use_full_grid:
                full_area_sum += int(measure_density.size)
                bbox_area_sum += int(measure_density.size)
                measure_p, n_p, ux_p, uz_p, _current_p, patch_diag = rk4_local_time_matter_step(
                    measure_density=measure_density,
                    u_x=u_x,
                    u_z=u_z,
                    metric_cov=metric_cov,
                    metric_inv=metric_inv,
                    metric_det=metric_det,
                    mass=mass,
                    dx=dx,
                    dz=dz,
                    dt=dt,
                    u_t_reference=ut_ref,
                    active_mask=stencil_mask_arr,
                    write_mask=patch_piece,
                    tau_a=a,
                    tau_b=b,
                    flux_divergence_mode=flux_divergence_mode,
                    rusanov_strength=rusanov_strength,
                    matching_flux_weight=matching_flux_weight,
                    matching_write_mask=write_mask,
                )
                accumulate_patch_diag(patch_diag)
                n_tau_input_abs_sum += abs(float(patch_diag.get("n_tau_input_sum", 0.0)))
                n_tau_delta_abs_sum += abs(float(patch_diag.get("n_tau_output_sum", 0.0)) - float(patch_diag.get("n_tau_input_sum", 0.0)))
                measure_acc[patch_piece] = measure_p[patch_piece]
                n_acc[patch_piece] = n_p[patch_piece]
                ux_acc[patch_piece] = ux_p[patch_piece]
                uz_acc[patch_piece] = uz_p[patch_piece]
        used_count += 1

    measure_next = np.where(write_mask, measure_acc, 0.0)
    ux_next = np.where(write_mask, ux_acc, 0.0)
    uz_next = np.where(write_mask, uz_acc, 0.0)
    n_next, current_next = current_from_measure_and_u(
        measure_density=measure_next,
        u_x=ux_next,
        u_z=uz_next,
        metric_cov=metric_cov,
        mass=mass,
        u_t_reference=ut_ref,
        active_mask=write_mask,
        rho_floor=1.0e-14,
        metric_inv=metric_inv,
        det_cov=metric_det,
    )
    output_total = float(np.sum(measure_next[write_mask]))
    diag = {
        "patches_used": float(used_count),
        "patch_components_used": float(component_count),
        "patch_bbox_area_sum": float(bbox_area_sum),
        "patch_full_area_sum": float(full_area_sum),
        "patch_bbox_area_fraction": float(bbox_area_sum / max(full_area_sum, 1)),
        "patch_bbox_fallback_count": float(bbox_fallback_count),
        "write_cell_count": float(np.count_nonzero(write_mask)),
        "stencil_cell_count": float(np.count_nonzero(stencil_mask_arr)),
        "fixed_buffer_cell_count": float(np.count_nonzero(stencil_mask_arr & (~write_mask))),
        "flux_divergence_mode_face": 1.0 if flux_divergence_mode == "face" else 0.0,
        "flux_divergence_mode_rusanov": 1.0 if flux_divergence_mode == "rusanov" else 0.0,
        "flux_divergence_mode_matching": 1.0 if flux_divergence_mode == "matching" else 0.0,
        "rusanov_strength": float(rusanov_strength),
        "matching_flux_weight_abs_n": 1.0 if matching_flux_weight == "abs_n" else 0.0,
        "measure_rel_delta": float(abs(output_total - total_input) / total_input),
        "n_tau_rel_delta_sum": float(n_tau_delta_abs_sum / max(n_tau_input_abs_sum, 1.0e-300)),
    }
    diag.update(fv_diag_acc)
    return measure_next, n_next, ux_next, uz_next, current_next, diag


def finite_volume_flux_divergence(
    n_cons: np.ndarray,
    flux_x: np.ndarray,
    flux_z: np.ndarray,
    dx: float,
    dz: float,
    stencil_mask: np.ndarray,
    write_mask: np.ndarray,
    mode: str = "face",
    rusanov_strength: float = 1.0,
    matching_write_mask: np.ndarray | None = None,
    matching_flux_weight: str = "abs_n",
) -> tuple[np.ndarray, dict[str, float]]:
    """Compute non-periodic finite-volume divergence from cell-face fluxes."""
    if mode not in {"face", "rusanov", "matching"}:
        raise ValueError(f"unknown finite-volume flux mode: {mode}")
    if matching_flux_weight not in {"equal", "abs_n"}:
        raise ValueError(f"unknown matching_flux_weight: {matching_flux_weight}")
    n = np.asarray(n_cons, dtype=float)
    fx = np.asarray(flux_x, dtype=float)
    fz = np.asarray(flux_z, dtype=float)
    stencil = np.asarray(stencil_mask, dtype=bool)
    write = np.asarray(write_mask, dtype=bool)
    match_write = write if matching_write_mask is None else (np.asarray(matching_write_mask, dtype=bool) | write)
    match_write = match_write & stencil
    nx, nz = fx.shape
    face_x = np.zeros((nx + 1, nz), dtype=float)
    face_z = np.zeros((nx, nz + 1), dtype=float)
    valid_x = stencil[:-1, :] & stencil[1:, :]
    valid_z = stencil[:, :-1] & stencil[:, 1:]
    fx_face = 0.5 * (fx[:-1, :] + fx[1:, :])
    fz_face = 0.5 * (fz[:, :-1] + fz[:, 1:])
    ax = np.zeros_like(valid_x, dtype=float)
    az = np.zeros_like(valid_z, dtype=float)
    if mode == "rusanov":
        n_lx = n[:-1, :]
        n_rx = n[1:, :]
        n_lz = n[:, :-1]
        n_rz = n[:, 1:]
        vx_l = np.abs(fx[:-1, :] / np.where(np.abs(n_lx) > 1.0e-300, n_lx, np.sign(n_lx) * 1.0e-300 + (n_lx == 0.0) * 1.0e-300))
        vx_r = np.abs(fx[1:, :] / np.where(np.abs(n_rx) > 1.0e-300, n_rx, np.sign(n_rx) * 1.0e-300 + (n_rx == 0.0) * 1.0e-300))
        vz_l = np.abs(fz[:, :-1] / np.where(np.abs(n_lz) > 1.0e-300, n_lz, np.sign(n_lz) * 1.0e-300 + (n_lz == 0.0) * 1.0e-300))
        vz_r = np.abs(fz[:, 1:] / np.where(np.abs(n_rz) > 1.0e-300, n_rz, np.sign(n_rz) * 1.0e-300 + (n_rz == 0.0) * 1.0e-300))
        ax = np.maximum(vx_l, vx_r)
        az = np.maximum(vz_l, vz_r)
        strength = max(float(rusanov_strength), 0.0)
        fx_face = 0.5 * (fx[:-1, :] + fx[1:, :]) - 0.5 * strength * ax * (n_rx - n_lx)
        fz_face = 0.5 * (fz[:, :-1] + fz[:, 1:]) - 0.5 * strength * az * (n_rz - n_lz)
    affected_x = write[:-1, :] | write[1:, :]
    affected_z = write[:, :-1] | write[:, 1:]
    matching_interface_x = valid_x & affected_x & (match_write[:-1, :] ^ match_write[1:, :])
    matching_interface_z = valid_z & affected_z & (match_write[:, :-1] ^ match_write[:, 1:])
    matching_flux_jump_l1 = 0.0
    matching_flux_jump_p95 = 0.0
    matching_flux_abs_sum = 0.0
    if mode == "matching":
        if matching_flux_weight == "equal":
            wx_l = np.ones_like(fx_face)
            wx_r = np.ones_like(fx_face)
            wz_l = np.ones_like(fz_face)
            wz_r = np.ones_like(fz_face)
        else:
            n_scale = max(float(np.nanmax(np.abs(n[stencil]))) if np.any(stencil) else 0.0, 1.0e-300)
            weight_floor = max(1.0e-12 * n_scale, 1.0e-300)
            wx_l = np.abs(n[:-1, :]) + weight_floor
            wx_r = np.abs(n[1:, :]) + weight_floor
            wz_l = np.abs(n[:, :-1]) + weight_floor
            wz_r = np.abs(n[:, 1:]) + weight_floor
        fx_match = (wx_l * fx[:-1, :] + wx_r * fx[1:, :]) / np.maximum(wx_l + wx_r, 1.0e-300)
        fz_match = (wz_l * fz[:, :-1] + wz_r * fz[:, 1:]) / np.maximum(wz_l + wz_r, 1.0e-300)
        fx_face = np.where(matching_interface_x, fx_match, fx_face)
        fz_face = np.where(matching_interface_z, fz_match, fz_face)
        jumps = []
        if np.any(matching_interface_x):
            x_jump = np.abs(fx[1:, :] - fx[:-1, :])[matching_interface_x]
            jumps.append(x_jump)
            matching_flux_abs_sum += float(np.sum(np.abs(fx_face[matching_interface_x])))
        if np.any(matching_interface_z):
            z_jump = np.abs(fz[:, 1:] - fz[:, :-1])[matching_interface_z]
            jumps.append(z_jump)
            matching_flux_abs_sum += float(np.sum(np.abs(fz_face[matching_interface_z])))
        if jumps:
            jump_values = np.concatenate(jumps)
            matching_flux_jump_l1 = float(np.sum(jump_values))
            matching_flux_jump_p95 = float(np.percentile(jump_values, 95.0))
    face_x[1:nx, :] = np.where(valid_x, fx_face, 0.0)
    face_z[:, 1:nz] = np.where(valid_z, fz_face, 0.0)
    div = (face_x[1:, :] - face_x[:-1, :]) / dx + (face_z[:, 1:] - face_z[:, :-1]) / dz
    div = np.where(write, div, 0.0)
    diag = {
        "fv_face_x_count": float(np.count_nonzero(valid_x)),
        "fv_face_z_count": float(np.count_nonzero(valid_z)),
        "fv_divergence_sum_write": float(np.sum(div[write])) if np.any(write) else 0.0,
        "fv_rusanov_max_speed_x": float(np.nanmax(ax[valid_x])) if mode == "rusanov" and np.any(valid_x) else 0.0,
        "fv_rusanov_max_speed_z": float(np.nanmax(az[valid_z])) if mode == "rusanov" and np.any(valid_z) else 0.0,
        "fv_rusanov_strength": float(rusanov_strength) if mode == "rusanov" else 0.0,
        "fv_matching_interface_x_count": float(np.count_nonzero(matching_interface_x)) if mode == "matching" else 0.0,
        "fv_matching_interface_z_count": float(np.count_nonzero(matching_interface_z)) if mode == "matching" else 0.0,
        "fv_matching_flux_jump_l1": matching_flux_jump_l1,
        "fv_matching_flux_jump_p95": matching_flux_jump_p95,
        "fv_matching_flux_abs_sum": matching_flux_abs_sum,
        "fv_matching_weight_abs_n": 1.0 if mode == "matching" and matching_flux_weight == "abs_n" else 0.0,
    }
    return div, diag


def rk4_local_time_matter_step(
    measure_density: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    metric_det: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    dt: float,
    u_t_reference: np.ndarray,
    active_mask: np.ndarray,
    tau_a: float,
    tau_b: float,
    write_mask: np.ndarray | None = None,
    flux_divergence_mode: str = "central",
    rusanov_strength: float = 1.0,
    matching_flux_weight: str = "abs_n",
    matching_write_mask: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    if flux_divergence_mode not in {"central", "face", "rusanov", "matching"}:
        raise ValueError(f"unknown flux_divergence_mode: {flux_divergence_mode}")
    if matching_flux_weight not in {"equal", "abs_n"}:
        raise ValueError(f"unknown matching_flux_weight: {matching_flux_weight}")
    stencil_mask = np.asarray(active_mask, dtype=bool)
    if write_mask is None:
        write_mask_arr = stencil_mask
    else:
        write_mask_arr = np.asarray(write_mask, dtype=bool)
        stencil_mask = stencil_mask | write_mask_arr
    if matching_write_mask is None:
        matching_write_mask_arr = write_mask_arr
    else:
        matching_write_mask_arr = (np.asarray(matching_write_mask, dtype=bool) | write_mask_arr) & stencil_mask
    metric_tau = local_time_covariant_metric(metric_cov, tau_a, tau_b)
    metric_tau_inv = local_time_contravariant_metric(metric_inv, tau_a, tau_b)
    metric_tau_det = metric_det
    u_t_lab = np.asarray(u_t_reference, dtype=float)
    u_tau = u_t_lab
    u_x_tau = np.asarray(u_x, dtype=float) - float(tau_a) * u_t_lab
    u_z_tau = np.asarray(u_z, dtype=float) - float(tau_b) * u_t_lab
    rho_lab = np.maximum(np.asarray(measure_density, dtype=float) / np.sqrt(np.maximum(np.abs(metric_det), 1.0e-300)), 1.0e-14)
    cons_tau0 = conservative_density_from_rho_u(
        rho=rho_lab,
        u_x=u_x_tau,
        u_z=u_z_tau,
        metric_cov_txz=metric_tau,
        mass=mass,
        branch="negative_frequency",
        u_t_reference=u_tau,
        rho_floor=1.0e-14,
        metric_inv_txz=metric_tau_inv,
        det_cov_txz=metric_tau_det,
    )
    n_tau = np.where(stencil_mask, cons_tau0["n_cons"], 0.0)
    input_measure_sum = max(float(np.sum(np.asarray(measure_density, dtype=float)[write_mask_arr])), 1.0e-300)
    input_n_tau_sum = float(np.sum(n_tau[write_mask_arr]))
    input_n_tau_stencil_sum = float(np.sum(n_tau[stencil_mask]))
    fv_diag_last: dict[str, float] = {}

    def rhs(n0: np.ndarray, ux0: np.ndarray, uz0: np.ndarray, ut_ref: np.ndarray):
        nonlocal fv_diag_last
        out = coordinate_matter_rhs_covector(
            n_cons=np.where(stencil_mask, n0, 0.0),
            u_x=np.where(stencil_mask, ux0, 0.0),
            u_z=np.where(stencil_mask, uz0, 0.0),
            metric_cov_txz=metric_tau,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=ut_ref,
            metric_inv_txz=metric_tau_inv,
            det_cov_txz=metric_tau_det,
        )
        if flux_divergence_mode in {"face", "rusanov", "matching"}:
            div, fv_diag_last = finite_volume_flux_divergence(
                np.where(stencil_mask, n0, 0.0),
                out["flux_x"],
                out["flux_z"],
                dx,
                dz,
                stencil_mask=stencil_mask,
                write_mask=write_mask_arr,
                mode=flux_divergence_mode,
                rusanov_strength=rusanov_strength,
                matching_write_mask=matching_write_mask_arr,
                matching_flux_weight=matching_flux_weight,
            )
            n_t = -div
        else:
            fv_diag_last = {}
            n_t = np.where(write_mask_arr, out["n_t"], 0.0)
        ux_t = np.where(write_mask_arr, out["u_x_t"], 0.0)
        uz_t = np.where(write_mask_arr, out["u_z_t"], 0.0)
        return n_t, ux_t, uz_t, out

    k1_n, k1_x, k1_z, o1 = rhs(n_tau, u_x_tau, u_z_tau, u_tau)
    k2_n, k2_x, k2_z, _ = rhs(n_tau + 0.5 * dt * k1_n, u_x_tau + 0.5 * dt * k1_x, u_z_tau + 0.5 * dt * k1_z, o1["u_t"])
    k3_n, k3_x, k3_z, _ = rhs(n_tau + 0.5 * dt * k2_n, u_x_tau + 0.5 * dt * k2_x, u_z_tau + 0.5 * dt * k2_z, o1["u_t"])
    k4_n, k4_x, k4_z, _ = rhs(n_tau + dt * k3_n, u_x_tau + dt * k3_x, u_z_tau + dt * k3_z, o1["u_t"])

    n_tau_next = n_tau + (dt / 6.0) * (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n)
    u_x_tau_next = u_x_tau + (dt / 6.0) * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x)
    u_z_tau_next = u_z_tau + (dt / 6.0) * (k1_z + 2.0 * k2_z + 2.0 * k3_z + k4_z)
    n_tau_next = np.where(write_mask_arr, n_tau_next, n_tau)
    u_x_tau_next = np.where(write_mask_arr, u_x_tau_next, u_x_tau)
    u_z_tau_next = np.where(write_mask_arr, u_z_tau_next, u_z_tau)

    final_local = coordinate_matter_rhs_covector(
        n_cons=np.where(stencil_mask, n_tau_next, 0.0),
        u_x=np.where(stencil_mask, u_x_tau_next, 0.0),
        u_z=np.where(stencil_mask, u_z_tau_next, 0.0),
        metric_cov_txz=metric_tau,
        mass=mass,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=o1["u_t"],
        metric_inv_txz=metric_tau_inv,
        det_cov_txz=metric_tau_det,
    )
    for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z"):
        if key in final_local:
            final_local[key] = np.where(write_mask_arr, final_local[key], 0.0)
    final_local["u_x"] = np.where(write_mask_arr, u_x_tau_next, 0.0)
    final_local["u_z"] = np.where(write_mask_arr, u_z_tau_next, 0.0)
    final_local["u_tau"] = final_local["u_t"]

    u_t_lab_next, u_x_lab_next, u_z_lab_next = lab_covector_from_local(
        final_local["u_tau"],
        final_local["u_x"],
        final_local["u_z"],
        tau_a,
        tau_b,
    )
    lab_state = conservative_density_from_rho_u(
        rho=np.where(write_mask_arr, final_local["rho"], 0.0),
        u_x=u_x_lab_next,
        u_z=u_z_lab_next,
        metric_cov_txz=metric_cov,
        mass=mass,
        branch="negative_frequency",
        u_t_reference=u_t_lab_next,
        rho_floor=1.0e-14,
        metric_inv_txz=metric_inv,
        det_cov_txz=metric_det,
    )
    for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z", "flux_x", "flux_z"):
        if key in lab_state:
            lab_state[key] = np.where(write_mask_arr, lab_state[key], 0.0)
    lab_state["u_x"] = np.where(write_mask_arr, u_x_lab_next, 0.0)
    lab_state["u_z"] = np.where(write_mask_arr, u_z_lab_next, 0.0)
    lab_state["n_cons"] = np.where(write_mask_arr, lab_state["n_cons"], 0.0)
    measure_next = np.where(write_mask_arr, final_local["sqrt_abs_g"] * np.maximum(final_local["rho"], 0.0), 0.0)
    output_measure_sum = float(np.sum(measure_next[write_mask_arr]))
    output_n_tau_sum = float(np.sum(n_tau_next[write_mask_arr]))
    output_n_tau_stencil_sum = float(np.sum(n_tau_next[stencil_mask]))
    diag = {
        "tau_a": float(tau_a),
        "tau_b": float(tau_b),
        "tau_norm_min_local": float(np.nanmin(final_local["discriminant"][write_mask_arr])) if np.any(write_mask_arr) else 0.0,
        "write_cell_count": float(np.count_nonzero(write_mask_arr)),
        "stencil_cell_count": float(np.count_nonzero(stencil_mask)),
        "fixed_buffer_cell_count": float(np.count_nonzero(stencil_mask & (~write_mask_arr))),
        "flux_divergence_mode_face": 1.0 if flux_divergence_mode == "face" else 0.0,
        "flux_divergence_mode_rusanov": 1.0 if flux_divergence_mode == "rusanov" else 0.0,
        "flux_divergence_mode_matching": 1.0 if flux_divergence_mode == "matching" else 0.0,
        "rusanov_strength": float(rusanov_strength),
        "matching_flux_weight_abs_n": 1.0 if matching_flux_weight == "abs_n" else 0.0,
        "measure_rel_delta": float(abs(output_measure_sum - input_measure_sum) / input_measure_sum),
        "n_tau_input_sum": input_n_tau_sum,
        "n_tau_output_sum": output_n_tau_sum,
        "n_tau_stencil_input_sum": input_n_tau_stencil_sum,
        "n_tau_stencil_output_sum": output_n_tau_stencil_sum,
        "n_tau_rel_delta": float(abs(output_n_tau_sum - input_n_tau_sum) / max(abs(input_n_tau_sum), 1.0e-300)),
    }
    diag.update({f"face_flux_{key}": value for key, value in fv_diag_last.items()})
    return measure_next, lab_state["n_cons"], lab_state["u_x"], lab_state["u_z"], lab_state, diag


def rk4_phase_matter_step(
    n_cons: np.ndarray,
    s_phase: np.ndarray,
    metric_cov: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    dt: float,
    s_t_reference: np.ndarray,
    active_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    def rhs(n0: np.ndarray, s0: np.ndarray, st_ref: np.ndarray):
        out = coordinate_matter_rhs(
            n_cons=np.where(active_mask, n0, 0.0),
            s=s0,
            metric_cov_txz=metric_cov,
            mass=mass,
            dx=dx,
            dz=dz,
            phase_branch="negative_frequency",
            s_t_reference=st_ref,
        )
        n_t = np.where(active_mask, out["n_t"], 0.0)
        s_t = np.where(active_mask, out["s_t"], 0.0)
        return n_t, s_t, out

    k1_n, k1_s, o1 = rhs(n_cons, s_phase, s_t_reference)
    k2_n, k2_s, _ = rhs(n_cons + 0.5 * dt * k1_n, s_phase + 0.5 * dt * k1_s, o1["s_t"])
    k3_n, k3_s, _ = rhs(n_cons + 0.5 * dt * k2_n, s_phase + 0.5 * dt * k2_s, o1["s_t"])
    k4_n, k4_s, _ = rhs(n_cons + dt * k3_n, s_phase + dt * k3_s, o1["s_t"])

    n_next = n_cons + (dt / 6.0) * (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n)
    s_next = s_phase + (dt / 6.0) * (k1_s + 2.0 * k2_s + 2.0 * k3_s + k4_s)
    n_next = np.where(active_mask, n_next, 0.0)
    s_next = unwrap_phase_like(s_next, s_phase)
    ux_next, uz_next = gradient_from_phase(s_next, active_mask, dx, dz)
    final = coordinate_matter_rhs(
        n_cons=n_next,
        s=s_next,
        metric_cov_txz=metric_cov,
        mass=mass,
        dx=dx,
        dz=dz,
        phase_branch="negative_frequency",
        s_t_reference=o1["s_t"],
    )
    for key in ("rho", "s_t", "discriminant", "j_t", "j_x", "j_z", "flux_x", "flux_z"):
        if key in final:
            final[key] = np.where(active_mask, final[key], 0.0)
    final["u_t"] = final["s_t"]
    final["u_x"] = ux_next
    final["u_z"] = uz_next
    if not (
        np.all(np.isfinite(n_next[active_mask]))
        and np.all(np.isfinite(s_next[active_mask]))
        and np.all(np.isfinite(ux_next[active_mask]))
        and np.all(np.isfinite(uz_next[active_mask]))
        and np.all(np.isfinite(final["rho"][active_mask]))
        and np.all(np.isfinite(final["u_t"][active_mask]))
    ):
        raise FloatingPointError("non-finite value produced by D phase-mode RK4 step")
    return n_next, s_next, ux_next, uz_next, final


def lorentzian_cov_mask(metric_cov: np.ndarray, det_floor: float = 1.0e-14) -> np.ndarray:
    sym = 0.5 * (metric_cov + np.swapaxes(metric_cov, -1, -2))
    mats = sym.reshape(-1, 3, 3)
    eigvals = np.linalg.eigvalsh(mats).reshape(sym.shape[:-2] + (3,))
    det = np.linalg.det(mats).reshape(sym.shape[:-2])
    pos_count = np.count_nonzero(eigvals > 0.0, axis=-1)
    neg_count = np.count_nonzero(eigvals < 0.0, axis=-1)
    finite = np.all(np.isfinite(sym), axis=(-1, -2))
    return finite & (np.abs(det) > det_floor) & (pos_count == 1) & (neg_count == 2)


def inertia_signature(metric_cov: np.ndarray, det_floor: float = 1.0e-14) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sym = 0.5 * (metric_cov + np.swapaxes(metric_cov, -1, -2))
    mats = sym.reshape(-1, 3, 3)
    eigvals = np.linalg.eigvalsh(mats).reshape(sym.shape[:-2] + (3,))
    det = np.linalg.det(mats).reshape(sym.shape[:-2])
    pos_count = np.count_nonzero(eigvals > 0.0, axis=-1)
    neg_count = np.count_nonzero(eigvals < 0.0, axis=-1)
    finite_nondegenerate = np.all(np.isfinite(sym), axis=(-1, -2)) & (np.abs(det) > det_floor)
    return pos_count, neg_count, finite_nondegenerate


def branch_preserving_signature_mask(
    base_metric_cov: np.ndarray,
    candidate_metric_cov: np.ndarray,
    det_floor: float = 1.0e-14,
) -> np.ndarray:
    """Accept nondegenerate updates that keep each point on its signature branch."""
    base_pos, base_neg, base_ok = inertia_signature(base_metric_cov, det_floor=det_floor)
    cand_pos, cand_neg, cand_ok = inertia_signature(candidate_metric_cov, det_floor=det_floor)
    return base_ok & cand_ok & (base_pos == cand_pos) & (base_neg == cand_neg)


def metric_neighbor_mean_no_wrap(metric_cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    total = np.zeros_like(metric_cov)
    count = np.zeros(metric_cov.shape[:2], dtype=float)

    total[1:, :, :, :] += metric_cov[:-1, :, :, :]
    count[1:, :] += 1.0
    total[:-1, :, :, :] += metric_cov[1:, :, :, :]
    count[:-1, :] += 1.0
    total[:, 1:, :, :] += metric_cov[:, :-1, :, :]
    count[:, 1:] += 1.0
    total[:, :-1, :, :] += metric_cov[:, 1:, :, :]
    count[:, :-1] += 1.0
    return total, count


def restricted_metric_neighbor_mean_no_wrap(
    metric_cov: np.ndarray,
    valid_neighbor_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    total = np.zeros_like(metric_cov)
    count = np.zeros(metric_cov.shape[:2], dtype=float)
    valid = np.asarray(valid_neighbor_mask, dtype=bool)

    src = valid[:-1, :]
    total[1:, :, :, :] += metric_cov[:-1, :, :, :] * src[..., None, None]
    count[1:, :] += src
    src = valid[1:, :]
    total[:-1, :, :, :] += metric_cov[1:, :, :, :] * src[..., None, None]
    count[:-1, :] += src
    src = valid[:, :-1]
    total[:, 1:, :, :] += metric_cov[:, :-1, :, :] * src[..., None, None]
    count[:, 1:] += src
    src = valid[:, 1:]
    total[:, :-1, :, :] += metric_cov[:, 1:, :, :] * src[..., None, None]
    count[:, :-1] += src
    return total, count


def restricted_field_neighbor_mean_no_wrap(
    field: np.ndarray,
    valid_neighbor_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    total = np.zeros_like(field)
    count = np.zeros(field.shape[:2], dtype=float)
    valid = np.asarray(valid_neighbor_mask, dtype=bool)
    tail = (None,) * max(field.ndim - 2, 0)

    src = valid[:-1, :]
    total[1:, ...] += field[:-1, ...] * src[(slice(None), slice(None)) + tail]
    count[1:, :] += src
    src = valid[1:, :]
    total[:-1, ...] += field[1:, ...] * src[(slice(None), slice(None)) + tail]
    count[:-1, :] += src
    src = valid[:, :-1]
    total[:, 1:, ...] += field[:, :-1, ...] * src[(slice(None), slice(None)) + tail]
    count[:, 1:] += src
    src = valid[:, 1:]
    total[:, :-1, ...] += field[:, 1:, ...] * src[(slice(None), slice(None)) + tail]
    count[:, :-1] += src
    return total, count


def active_boundary_mask(active_mask: np.ndarray) -> np.ndarray:
    active = np.asarray(active_mask, dtype=bool)
    interior = active.copy()
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        interior &= shift_mask_no_wrap(active, di, dj)
    return active & (~interior)


def low_matter_flat_anchor_mask(
    rho: np.ndarray,
    measure: np.ndarray,
    active_mask: np.ndarray,
    rho_frac: float,
    measure_frac: float,
    boundary_layers: int,
) -> np.ndarray:
    """Asymptotically flat Dirichlet anchor on the low-matter outer boundary.

    The mask deliberately grows inward only from the numerical active boundary.
    It should not select interior destructive-interference nodes, because those
    can be physical features rather than asymptotic vacuum.
    """
    if boundary_layers <= 0:
        return np.zeros_like(active_mask, dtype=bool)
    active = np.asarray(active_mask, dtype=bool)
    boundary = active_boundary_mask(active)
    shell = boundary.copy()
    grow = boundary.copy()
    for _ in range(max(int(boundary_layers) - 1, 0)):
        grow = dilate_mask_4(grow, 1) & active
        shell |= grow
    rho_threshold = float(rho_frac) * float(np.max(rho))
    measure_threshold = float(measure_frac) * float(np.max(measure))
    return shell & (rho <= rho_threshold) & (measure <= measure_threshold)


def adm_variables_from_metric_cov(
    metric_cov: np.ndarray,
    floor: float = 1.0e-14,
) -> tuple[np.ndarray, np.ndarray]:
    g = 0.5 * (np.asarray(metric_cov, dtype=float) + np.swapaxes(metric_cov, -1, -2))
    h_xx = -g[..., 1, 1]
    h_xz = -g[..., 1, 2]
    h_zz = -g[..., 2, 2]
    det_h = h_xx * h_zz - h_xz * h_xz
    valid_h = np.isfinite(det_h) & (h_xx > floor) & (det_h > floor)

    beta_x = -(h_zz * g[..., 0, 1] - h_xz * g[..., 0, 2]) / np.maximum(det_h, floor)
    beta_z = -((-h_xz) * g[..., 0, 1] + h_xx * g[..., 0, 2]) / np.maximum(det_h, floor)
    beta_h_beta = h_xx * beta_x * beta_x + 2.0 * h_xz * beta_x * beta_z + h_zz * beta_z * beta_z
    lapse_sq = g[..., 0, 0] + beta_h_beta
    valid = valid_h & np.isfinite(lapse_sq) & (lapse_sq > floor)

    h_xx_safe = np.maximum(h_xx, floor)
    det_h_safe = np.maximum(det_h, floor)
    lapse = np.sqrt(np.maximum(lapse_sq, floor))
    ell_x = np.sqrt(h_xx_safe)
    chol_b = h_xz / np.maximum(ell_x, np.sqrt(floor))
    ell_z_sq = det_h_safe / h_xx_safe

    variables = np.stack(
        [
            np.log(np.maximum(lapse, floor)),
            beta_x,
            beta_z,
            np.log(np.maximum(ell_x, floor)),
            chol_b,
            0.5 * np.log(np.maximum(ell_z_sq, floor)),
        ],
        axis=-1,
    )
    finite = np.all(np.isfinite(variables), axis=-1)
    return np.real_if_close(variables), valid & finite


def metric_cov_from_adm_variables(variables: np.ndarray) -> np.ndarray:
    vals = np.asarray(variables, dtype=float)
    lapse = np.exp(vals[..., 0])
    beta_x = vals[..., 1]
    beta_z = vals[..., 2]
    ell_x = np.exp(vals[..., 3])
    chol_b = vals[..., 4]
    ell_z = np.exp(vals[..., 5])

    h_xx = ell_x * ell_x
    h_xz = ell_x * chol_b
    h_zz = chol_b * chol_b + ell_z * ell_z
    beta_h_beta = h_xx * beta_x * beta_x + 2.0 * h_xz * beta_x * beta_z + h_zz * beta_z * beta_z

    g = np.zeros(vals.shape[:-1] + (3, 3), dtype=float)
    g[..., 0, 0] = lapse * lapse - beta_h_beta
    g[..., 0, 1] = -(h_xx * beta_x + h_xz * beta_z)
    g[..., 1, 0] = g[..., 0, 1]
    g[..., 0, 2] = -(h_xz * beta_x + h_zz * beta_z)
    g[..., 2, 0] = g[..., 0, 2]
    g[..., 1, 1] = -h_xx
    g[..., 1, 2] = -h_xz
    g[..., 2, 1] = -h_xz
    g[..., 2, 2] = -h_zz
    return np.real_if_close(g)


def anchored_solve_components(solve_mask: np.ndarray, anchor_mask: np.ndarray) -> tuple[np.ndarray, int, int]:
    solve = np.asarray(solve_mask, dtype=bool)
    anchor = np.asarray(anchor_mask, dtype=bool)
    seen = np.zeros_like(solve, dtype=bool)
    eligible = np.zeros_like(solve, dtype=bool)
    anchored_count = 0
    orphan_count = 0
    nx, nz = solve.shape

    for i0, j0 in np.argwhere(solve):
        if seen[i0, j0]:
            continue
        stack = [(int(i0), int(j0))]
        seen[i0, j0] = True
        cells: list[tuple[int, int]] = []
        touches_anchor = False
        while stack:
            i, j = stack.pop()
            cells.append((i, j))
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ni = i + di
                nj = j + dj
                if ni < 0 or ni >= nx or nj < 0 or nj >= nz:
                    continue
                if anchor[ni, nj]:
                    touches_anchor = True
                elif solve[ni, nj] and not seen[ni, nj]:
                    seen[ni, nj] = True
                    stack.append((ni, nj))
        if touches_anchor:
            anchored_count += 1
            for i, j in cells:
                eligible[i, j] = True
        else:
            orphan_count += 1
    return eligible, anchored_count, orphan_count


def metric_roughness(metric_cov: np.ndarray, dx: float, dz: float) -> np.ndarray:
    grad_x, grad_z = np.gradient(metric_cov, dx, dz, axis=(0, 1), edge_order=1)
    return np.sqrt(
        np.sum(grad_x * grad_x + grad_z * grad_z, axis=(-2, -1))
    )


def empty_metric_extension_diagnostic() -> dict[str, float]:
    return {
        "metric_extension_solve_count": 0.0,
        "metric_extension_solve_fraction": 0.0,
        "metric_extension_anchor_count": 0.0,
        "metric_extension_flat_anchor_count": 0.0,
        "metric_extension_anchored_component_count": 0.0,
        "metric_extension_orphan_component_count": 0.0,
        "metric_extension_lorentz_fallback_count": 0.0,
        "metric_extension_change_rel_p95": 0.0,
        "metric_extension_roughness_before_p95": 0.0,
        "metric_extension_roughness_after_p95": 0.0,
        "metric_extension_adm_valid_fraction": 0.0,
        "metric_extension_admissible_blend": 0.0,
        "metric_extension_admissible_attempts": 0.0,
        "metric_extension_admissible_disc_min": 0.0,
        "metric_extension_admissible_rho_min": 0.0,
        "metric_extension_last_rejected_disc_min": 0.0,
        "metric_extension_last_rejected_rho_min": 0.0,
        "metric_extension_admissible_lorentz_failed_attempts": 0.0,
        "metric_extension_admissible_finite_failed_attempts": 0.0,
        "metric_extension_admissible_disc_failed_attempts": 0.0,
        "metric_extension_admissible_rho_failed_attempts": 0.0,
        "metric_extension_last_rejected_blend": 0.0,
        "metric_extension_last_rejected_reason_code": 0.0,
    }


def empty_matter_projection_diagnostic() -> dict[str, float]:
    return {
        "matter_projection_enabled": 0.0,
        "matter_projection_initial_bad_count": 0.0,
        "matter_projection_final_bad_count": 0.0,
        "matter_projection_blend_iterations": 0.0,
        "matter_projection_disc_bad_count": 0.0,
        "matter_projection_rho_bad_count": 0.0,
        "matter_projection_finite_bad_count": 0.0,
        "matter_projection_min_disc": 0.0,
        "matter_projection_min_rho": 0.0,
        "matter_limiter_applied_count": 0.0,
        "matter_limiter_max_abs_n_delta": 0.0,
    }


def empty_local_time_diagnostic() -> dict[str, float]:
    return {
        "local_time_tau_a": 0.0,
        "local_time_tau_b": 0.0,
        "local_time_score": 0.0,
        "local_time_admissible_fraction": 0.0,
        "local_time_atlas_reused": 0.0,
        "local_time_tile_size": 1.0,
        "local_time_patch_count": 0.0,
        "local_time_candidate_count": 0.0,
        "local_time_stationary_candidate_raw_count": 0.0,
        "local_time_stationary_candidate_unique_count": 0.0,
        "local_time_patch_label_edges": 0.0,
        "local_time_patch_lab_fraction": 0.0,
        "local_time_uncovered_trusted_measure_fraction": 0.0,
        "local_time_abs_jtau_p05": 0.0,
        "local_time_abs_jtau_p50": 0.0,
        "local_time_abs_jtau_p95": 0.0,
        "local_time_tau_norm_min": 0.0,
        "local_time_tau_norm_p05": 0.0,
        "local_time_fallback_to_lab": 0.0,
        "local_time_patch_bbox_area_fraction": 0.0,
        "local_time_patch_bbox_fallback_count": 0.0,
        "local_time_patch_component_count": 0.0,
        "local_time_measure_rel_delta": 0.0,
        "local_time_n_tau_rel_delta_sum": 0.0,
    }


def project_covector_matter_update(
    *,
    base_n_cons: np.ndarray,
    base_u_x: np.ndarray,
    base_u_z: np.ndarray,
    n_cons: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    u_t_reference: np.ndarray,
    projection_mask: np.ndarray,
    disc_floor: float,
    rho_floor: float,
    max_blends: int,
    positivity_limiter: bool,
    limiter_rho_floor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    """Project an explicit matter update back into the algebraic constraint set.

    This is a local constrained-integrator step: it does not change the
    continuum equations, but prevents an explicit RK proposal from crossing the
    mass-shell/positive-density boundary in one grid cell.
    """
    mask = np.asarray(projection_mask, dtype=bool)
    n_out = np.asarray(n_cons, dtype=float).copy()
    ux_out = np.asarray(u_x, dtype=float).copy()
    uz_out = np.asarray(u_z, dtype=float).copy()

    def recover() -> dict[str, np.ndarray]:
        rec = coordinate_matter_rhs_covector(
            n_cons=n_out,
            u_x=ux_out,
            u_z=uz_out,
            metric_cov_txz=metric_cov,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=u_t_reference,
        )
        for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z"):
            if key in rec:
                rec[key] = np.where(mask, rec[key], 0.0)
        return rec

    def bad_mask(rec: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        finite_bad = mask & (
            (~np.isfinite(n_out))
            | (~np.isfinite(ux_out))
            | (~np.isfinite(uz_out))
            | (~np.isfinite(rec["rho"]))
            | (~np.isfinite(rec["u_t"]))
            | (~np.isfinite(rec["discriminant"]))
        )
        disc_bad = mask & (rec["discriminant"] < float(disc_floor))
        rho_bad = mask & (rec["rho"] < float(rho_floor))
        bad = finite_bad | disc_bad | rho_bad
        return bad, finite_bad, disc_bad, rho_bad

    rec = recover()
    bad, finite_bad, disc_bad, rho_bad = bad_mask(rec)
    initial_bad = int(np.count_nonzero(bad))
    iterations = 0
    while np.any(bad) and iterations < max(int(max_blends), 0):
        n_out[bad] = 0.5 * (n_out[bad] + base_n_cons[bad])
        ux_out[bad] = 0.5 * (ux_out[bad] + base_u_x[bad])
        uz_out[bad] = 0.5 * (uz_out[bad] + base_u_z[bad])
        rec = recover()
        bad, finite_bad, disc_bad, rho_bad = bad_mask(rec)
        iterations += 1

    limiter_count = 0
    limiter_max_abs_n_delta = 0.0
    if positivity_limiter and np.any(bad):
        limiter_mask = bad & rho_bad & (~finite_bad) & (~disc_bad)
        if np.any(limiter_mask):
            old_n = n_out[limiter_mask].copy()
            rho_target = max(float(limiter_rho_floor), 0.0)
            rho_safe = np.where(np.abs(rec["rho"]) > 1.0e-300, rec["rho"], np.sign(rec["rho"]) * 1.0e-300 + (rec["rho"] == 0.0) * 1.0e-300)
            jt_unit = n_out / (np.maximum(rec["sqrt_abs_g"], 1.0e-300) * rho_safe)
            finite_jt = limiter_mask & np.isfinite(jt_unit)
            n_out[finite_jt] = rho_target * rec["sqrt_abs_g"][finite_jt] * jt_unit[finite_jt]
            limiter_count = int(np.count_nonzero(finite_jt))
            if limiter_count:
                limiter_max_abs_n_delta = float(np.nanmax(np.abs(n_out[limiter_mask] - old_n)))
            rec = recover()
            bad, finite_bad, disc_bad, rho_bad = bad_mask(rec)

    diag = empty_matter_projection_diagnostic()
    diag.update(
        {
            "matter_projection_enabled": 1.0,
            "matter_projection_initial_bad_count": float(initial_bad),
            "matter_projection_final_bad_count": float(np.count_nonzero(bad)),
            "matter_projection_blend_iterations": float(iterations),
            "matter_projection_disc_bad_count": float(np.count_nonzero(disc_bad)),
            "matter_projection_rho_bad_count": float(np.count_nonzero(rho_bad)),
            "matter_projection_finite_bad_count": float(np.count_nonzero(finite_bad)),
            "matter_projection_min_disc": float(np.nanmin(rec["discriminant"][mask])) if np.any(mask) else 0.0,
            "matter_projection_min_rho": float(np.nanmin(rec["rho"][mask])) if np.any(mask) else 0.0,
            "matter_limiter_applied_count": float(limiter_count),
            "matter_limiter_max_abs_n_delta": float(limiter_max_abs_n_delta),
        }
    )
    return n_out, ux_out, uz_out, rec, diag


def restricted_harmonic_metric_extension(
    metric_cov: np.ndarray,
    solve_mask: np.ndarray,
    anchor_mask: np.ndarray,
    flat_anchor_mask: np.ndarray | None,
    active_mask: np.ndarray,
    dx: float,
    dz: float,
    iterations: int,
    anchor_weight: float,
    enforce_lorentz: bool,
    variable_mode: str,
    adm_update_fields: str,
    max_rel_change: float,
) -> tuple[np.ndarray, dict[str, float]]:
    """Choose a saturated-bulk representative by constrained harmonic extension.

    In the saturated region f_R is effectively zero, so the local highest-order
    metric equation no longer selects a unique representative.  This closure
    fixes the representative by minimizing a discrete metric-gradient energy
    with Dirichlet data inherited from the non-saturated/interface/outer cells.
    It does not alter the matter equations, source tensor, or D-branch action.
    """
    base = np.asarray(metric_cov, dtype=float)
    active = np.asarray(active_mask, dtype=bool)
    eta = flat_metric_cov_like(active_mask)
    anchor = np.asarray(anchor_mask, dtype=bool) & active
    flat_anchor = np.zeros_like(active, dtype=bool) if flat_anchor_mask is None else (np.asarray(flat_anchor_mask, dtype=bool) & anchor)
    base_with_flat = base.copy()
    base_with_flat[flat_anchor] = eta[flat_anchor]
    if variable_mode == "adm":
        _, valid_adm = adm_variables_from_metric_cov(base_with_flat)
        anchor &= valid_adm
    elif variable_mode != "covariant":
        raise ValueError(f"unknown metric extension variable mode: {variable_mode}")
    solve = np.asarray(solve_mask, dtype=bool) & active & (~anchor)
    solve &= np.all(np.isfinite(base), axis=(-1, -2))
    anchor &= np.all(np.isfinite(base), axis=(-1, -2))
    flat_anchor &= anchor
    solve, anchored_components, orphan_components = anchored_solve_components(solve, anchor)
    if not np.any(solve):
        diag = empty_metric_extension_diagnostic()
        diag.update(
            {
                "metric_extension_anchor_count": float(np.count_nonzero(anchor)),
                "metric_extension_flat_anchor_count": float(np.count_nonzero(flat_anchor)),
                "metric_extension_anchored_component_count": float(anchored_components),
                "metric_extension_orphan_component_count": float(orphan_components),
            }
        )
        return base.copy(), diag

    out = base.copy()
    out[~active_mask] = eta[~active_mask]
    out[flat_anchor] = eta[flat_anchor]
    rough_before = metric_roughness(out, dx, dz)
    valid_neighbor_mask = solve | anchor
    adm_valid_fraction = 0.0

    omega = max(float(anchor_weight), 0.0)
    if variable_mode == "adm":
        variables, valid_adm = adm_variables_from_metric_cov(out)
        adm_valid_fraction = float(np.count_nonzero(valid_adm & active) / max(int(np.count_nonzero(active)), 1))
        flat_vars, _ = adm_variables_from_metric_cov(eta)
        variables[flat_anchor] = flat_vars[flat_anchor]
        anchor_values = variables.copy()
        if adm_update_fields == "all":
            update_components = np.ones(variables.shape[-1], dtype=bool)
        elif adm_update_fields == "lapse":
            update_components = np.array([True, False, False, False, False, False], dtype=bool)
        elif adm_update_fields == "lapse_shift":
            update_components = np.array([True, True, True, False, False, False], dtype=bool)
        else:
            raise ValueError(f"unknown ADM update field set: {adm_update_fields}")
        solve &= valid_adm
        valid_neighbor_mask = solve | anchor
        for _ in range(max(int(iterations), 0)):
            neighbor_sum, neighbor_count = restricted_field_neighbor_mean_no_wrap(variables, valid_neighbor_mask)
            denom = np.maximum(neighbor_count[..., None] + omega, 1.0)
            proposal = (neighbor_sum + omega * anchor_values) / denom
            variables[solve, :] = np.where(update_components, proposal[solve, :], anchor_values[solve, :])
            variables[flat_anchor] = flat_vars[flat_anchor]
        out = metric_cov_from_adm_variables(variables)
        out[active & (~valid_adm)] = base_with_flat[active & (~valid_adm)]
        out[~active_mask] = eta[~active_mask]
        out[flat_anchor] = eta[flat_anchor]
    else:
        anchor_values = out.copy()
        for _ in range(max(int(iterations), 0)):
            neighbor_sum, neighbor_count = restricted_metric_neighbor_mean_no_wrap(out, valid_neighbor_mask)
            denom = np.maximum(neighbor_count[..., None, None] + omega, 1.0)
            proposal = (neighbor_sum + omega * anchor_values) / denom
            out[solve] = proposal[solve]
            out[~active_mask] = eta[~active_mask]
            out[flat_anchor] = eta[flat_anchor]
            out = 0.5 * (out + np.swapaxes(out, -1, -2))

    fallback_count = 0
    if enforce_lorentz:
        good = branch_preserving_signature_mask(base, out)
        bad = solve & (~good)
        fallback_count = int(np.count_nonzero(bad))
        out[bad] = base[bad]

    diff_norm = tensor_norm(out - base)
    base_norm = np.maximum(tensor_norm(base), 1.0e-300)
    rel_change = diff_norm / base_norm
    rel_max = float(max_rel_change)
    cap_blend = 1.0
    if rel_max > 0.0 and np.any(solve):
        solve_rel_max = float(np.nanmax(rel_change[solve]))
        if np.isfinite(solve_rel_max) and solve_rel_max > rel_max:
            cap_blend = rel_max / solve_rel_max
            out = base + cap_blend * (out - base)
            out[~active_mask] = eta[~active_mask]
            out[flat_anchor] = eta[flat_anchor]
            out = 0.5 * (out + np.swapaxes(out, -1, -2))
            diff_norm = tensor_norm(out - base)
            rel_change = diff_norm / base_norm
    rough_after = metric_roughness(out, dx, dz)
    solve_count = int(np.count_nonzero(solve))
    active_count = max(int(np.count_nonzero(active_mask)), 1)
    diag = {
        "metric_extension_solve_count": float(solve_count),
        "metric_extension_solve_fraction": float(solve_count / active_count),
        "metric_extension_anchor_count": float(np.count_nonzero(valid_neighbor_mask & (~solve))),
        "metric_extension_flat_anchor_count": float(np.count_nonzero(flat_anchor)),
        "metric_extension_anchored_component_count": float(anchored_components),
        "metric_extension_orphan_component_count": float(orphan_components),
        "metric_extension_lorentz_fallback_count": float(fallback_count),
        "metric_extension_change_rel_p95": scalar_stats(rel_change, solve)["p95"],
        "metric_extension_roughness_before_p95": scalar_stats(rough_before, solve)["p95"],
        "metric_extension_roughness_after_p95": scalar_stats(rough_after, solve)["p95"],
        "metric_extension_adm_valid_fraction": adm_valid_fraction,
        "metric_extension_adm_update_component_count": float(np.count_nonzero(update_components))
        if variable_mode == "adm"
        else 0.0,
        "metric_extension_adm_preserved_invalid_count": float(np.count_nonzero(active & (~valid_adm)))
        if variable_mode == "adm"
        else 0.0,
        "metric_extension_rel_change_cap_blend": float(cap_blend),
    }
    return out, diag


def select_admissible_metric_update(
    base_metric: np.ndarray,
    target_metric: np.ndarray,
    n_cons: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    u_t_reference: np.ndarray,
    active_mask: np.ndarray,
    check_mask: np.ndarray,
    disc_tolerance: float,
    rho_floor: float,
    max_halvings: int,
) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    check = np.asarray(check_mask, dtype=bool) & np.asarray(active_mask, dtype=bool)
    if not np.any(check):
        check = np.asarray(active_mask, dtype=bool)

    base = np.asarray(base_metric, dtype=float)
    target = np.asarray(target_metric, dtype=float)
    if not np.any(np.abs(target - base) > 0.0):
        rec = coordinate_matter_rhs_covector(
            n_cons=n_cons,
            u_x=u_x,
            u_z=u_z,
            metric_cov_txz=base,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=u_t_reference,
        )
        for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z"):
            if key in rec:
                rec[key] = np.where(active_mask, rec[key], 0.0)
        diag = {
            "metric_extension_admissible_blend": 0.0,
            "metric_extension_admissible_attempts": 0.0,
            "metric_extension_admissible_disc_min": float(np.nanmin(rec["discriminant"][check])),
            "metric_extension_admissible_rho_min": float(np.nanmin(rec["rho"][check])),
        }
        return base.copy(), rec, diag

    blends = [1.0 / (2.0 ** k) for k in range(max(int(max_halvings), 0) + 1)] + [0.0]
    attempts = 0
    last_rec: dict[str, np.ndarray] | None = None
    last_disc_min = 0.0
    last_rho_min = 0.0
    last_rejected_blend = 0.0
    last_rejected_reason_code = 0.0
    lorentz_failed_attempts = 0
    finite_failed_attempts = 0
    disc_failed_attempts = 0
    rho_failed_attempts = 0

    eta = flat_metric_cov_like(active_mask)
    for blend in blends:
        attempts += 1
        candidate = base + float(blend) * (target - base)
        candidate[~active_mask] = eta[~active_mask]
        candidate = 0.5 * (candidate + np.swapaxes(candidate, -1, -2))
        if not np.all(branch_preserving_signature_mask(base, candidate)[check]):
            if blend > 0.0:
                lorentz_failed_attempts += 1
                last_rejected_blend = float(blend)
                last_rejected_reason_code = 1.0
            continue
        rec = coordinate_matter_rhs_covector(
            n_cons=n_cons,
            u_x=u_x,
            u_z=u_z,
            metric_cov_txz=candidate,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=u_t_reference,
        )
        for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z"):
            if key in rec:
                rec[key] = np.where(active_mask, rec[key], 0.0)
        finite = (
            np.all(np.isfinite(candidate[check]))
            and np.all(np.isfinite(rec["rho"][check]))
            and np.all(np.isfinite(rec["u_t"][check]))
            and np.all(np.isfinite(rec["discriminant"][check]))
        )
        last_rec = rec
        last_disc_min = float(np.nanmin(rec["discriminant"][check]))
        last_rho_min = float(np.nanmin(rec["rho"][check]))
        if finite and last_disc_min >= -abs(float(disc_tolerance)) and last_rho_min >= -abs(float(rho_floor)):
            diag = {
                "metric_extension_admissible_blend": float(blend),
                "metric_extension_admissible_attempts": float(attempts),
                "metric_extension_admissible_disc_min": last_disc_min,
                "metric_extension_admissible_rho_min": last_rho_min,
                "metric_extension_admissible_lorentz_failed_attempts": float(lorentz_failed_attempts),
                "metric_extension_admissible_finite_failed_attempts": float(finite_failed_attempts),
                "metric_extension_admissible_disc_failed_attempts": float(disc_failed_attempts),
                "metric_extension_admissible_rho_failed_attempts": float(rho_failed_attempts),
                "metric_extension_last_rejected_blend": float(last_rejected_blend),
                "metric_extension_last_rejected_reason_code": float(last_rejected_reason_code),
            }
            return candidate, rec, diag
        if blend > 0.0:
            last_rejected_blend = float(blend)
            if not finite:
                finite_failed_attempts += 1
                last_rejected_reason_code = 2.0
            elif last_disc_min < -abs(float(disc_tolerance)):
                disc_failed_attempts += 1
                last_rejected_reason_code = 3.0
            elif last_rho_min < -abs(float(rho_floor)):
                rho_failed_attempts += 1
                last_rejected_reason_code = 4.0

    rec = coordinate_matter_rhs_covector(
        n_cons=n_cons,
        u_x=u_x,
        u_z=u_z,
        metric_cov_txz=base,
        mass=mass,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=u_t_reference,
    )
    for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z"):
        if key in rec:
            rec[key] = np.where(active_mask, rec[key], 0.0)
    diag = {
        "metric_extension_admissible_blend": 0.0,
        "metric_extension_admissible_attempts": float(attempts),
        "metric_extension_admissible_disc_min": float(np.nanmin(rec["discriminant"][check]))
        if np.any(check)
        else last_disc_min,
        "metric_extension_admissible_rho_min": float(np.nanmin(rec["rho"][check])) if np.any(check) else last_rho_min,
        "metric_extension_admissible_lorentz_failed_attempts": float(lorentz_failed_attempts),
        "metric_extension_admissible_finite_failed_attempts": float(finite_failed_attempts),
        "metric_extension_admissible_disc_failed_attempts": float(disc_failed_attempts),
        "metric_extension_admissible_rho_failed_attempts": float(rho_failed_attempts),
        "metric_extension_last_rejected_blend": float(last_rejected_blend),
        "metric_extension_last_rejected_reason_code": float(last_rejected_reason_code),
    }
    if last_rec is not None:
        diag["metric_extension_last_rejected_disc_min"] = last_disc_min
        diag["metric_extension_last_rejected_rho_min"] = last_rho_min
    return base.copy(), rec, diag


def reference_a_snapshot(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    mass: float,
    rho_floor: float,
    active: np.ndarray,
) -> dict[str, np.ndarray]:
    snap = localized_direct_tilde_coordinate_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=t,
        rho_floor=rho_floor,
    )
    bohm = snap["bohm"]
    metric_cov = sanitize_metric_tail(real_array(snap["cov_txz"]), active)
    measure = tilde_measure_from_bohm(bohm, mass)
    rho_tilde = tilde_density_from_measure(measure, metric_cov, rho_floor)
    return {
        "rho": real_array(bohm["rho"]),
        "rho_tilde": rho_tilde,
        "ntilde_measure": measure,
        "s_t": real_array(bohm["s_t"]),
        "s_x": real_array(bohm["s_x"]),
        "s_z": real_array(bohm["s_z"]),
        "X": real_array(bohm["X"]),
        "metric_cov": metric_cov,
        "sqrt_abs_g": safe_sqrt_abs_det(metric_cov),
    }


def domain_masks(raw_y: np.ndarray, support: np.ndarray, weak_y: float, saturation_phi: float) -> dict[str, np.ndarray]:
    abs_y = np.abs(raw_y)
    # D branch phi=f_R=sech^2(y).  Avoid overflow by using cosh only in the moderate range.
    clipped = np.clip(abs_y, 0.0, 350.0)
    phi = 1.0 / (np.cosh(clipped) ** 2)
    weak = support & (abs_y <= weak_y)
    saturated = support & (phi <= saturation_phi)
    transition = support & (~weak) & (~saturated)
    return {
        "weak": weak,
        "saturated": saturated,
        "transition": transition,
        "phi": phi,
    }


def plateau_tensor_diagnostic(
    raw_y: np.ndarray,
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    rho_tilde: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    mass: float,
    ell: float,
    mp: float,
    mask: np.ndarray,
) -> dict[str, float]:
    stress, tilde_x = stress_tensor_tilde(metric_cov, metric_inv, rho_tilde, u_t, u_x, u_z, m=mass)
    f = np.tanh(raw_y) / (ell * ell)
    lhs = -0.5 * f[..., None, None] * metric_cov
    rhs = real_array(stress) / (mp * mp)
    residual = lhs - rhs
    lhs_norm = tensor_norm(lhs)
    rhs_norm = tensor_norm(rhs)
    residual_norm = tensor_norm(residual)
    denom = np.maximum(lhs_norm + rhs_norm, 1.0e-300)
    rel = residual_norm / denom
    return {
        "plateau_lhs_p95": scalar_stats(lhs_norm, mask)["p95"],
        "plateau_rhs_p95": scalar_stats(rhs_norm, mask)["p95"],
        "plateau_residual_p95": scalar_stats(residual_norm, mask)["p95"],
        "plateau_relative_p95": scalar_stats(rel, mask)["p95"],
        "mass_shell_defect_p95": scalar_stats(tilde_x - mass * mass, mask)["p95"],
    }


def render_summary(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    ref: dict[str, np.ndarray],
    d_measure: np.ndarray,
    rho_tilde_d: np.ndarray,
    raw_y: np.ndarray | None,
    masks: dict[str, np.ndarray] | None,
    interface_rows: list[dict[str, float]],
    records: list[dict[str, float]],
    support: np.ndarray,
    axis_unit: str = "natural length unit",
) -> None:
    plt = load_pyplot()
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.4), constrained_layout=True)

    vmax_measure = max(
        float(np.percentile(ref["ntilde_measure"][support], 99.5)),
        float(np.percentile(d_measure[support], 99.5)),
        1.0e-16,
    )
    im0 = axes[0, 0].pcolormesh(
        xg,
        zg,
        ref["ntilde_measure"],
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax_measure,
    )
    axes[0, 0].set_title("A reference: sqrt(|g~|) rho~")
    fig.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].pcolormesh(
        xg,
        zg,
        d_measure,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax_measure,
    )
    add_lines(axes[0, 1], interface_rows, colors="white", linewidth=0.7)
    axes[0, 1].set_title("D evolved: sqrt(|g~|) rho~")
    fig.colorbar(im1, ax=axes[0, 1])

    diff = d_measure - ref["ntilde_measure"]
    dmax = max(float(np.percentile(np.abs(diff[support]), 99.0)), 1.0e-16)
    im2 = axes[0, 2].pcolormesh(xg, zg, diff, shading="auto", cmap="coolwarm", vmin=-dmax, vmax=dmax)
    axes[0, 2].set_title("D - A transformed measure")
    fig.colorbar(im2, ax=axes[0, 2])

    if raw_y is not None:
        clip = max(float(np.percentile(np.abs(raw_y[support]), 95.0)), 1.0)
        im3 = axes[1, 0].pcolormesh(
            xg,
            zg,
            np.clip(raw_y, -clip, clip),
            shading="auto",
            cmap="coolwarm",
            vmin=-clip,
            vmax=clip,
        )
        add_lines(axes[1, 0], interface_rows, colors="black", linewidth=0.6)
        axes[1, 0].set_title("raw_y = ell^2 R~ with y=+-1 lines")
        fig.colorbar(im3, ax=axes[1, 0])
    else:
        axes[1, 0].set_axis_off()

    times = np.asarray([r["t"] for r in records], dtype=float)
    axes[1, 1].plot(times, [r["measure_rel_l1_support"] for r in records], marker="o", label="measure rel L1")
    axes[1, 1].plot(times, [r["n_cons_rel_l1_support"] for r in records], marker="o", label="conserved n rel L1")
    axes[1, 1].set_title("D vs A transformed observables")
    axes[1, 1].set_xlabel("t")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.25)

    axes[1, 2].plot(times, [r["weak_fraction"] for r in records], marker="o", label="weak")
    axes[1, 2].plot(times, [r["transition_fraction"] for r in records], marker="o", label="transition")
    axes[1, 2].plot(times, [r["saturated_fraction"] for r in records], marker="o", label="saturated")
    axes[1, 2].set_title("three-domain fractions on support")
    axes[1, 2].set_xlabel("t")
    axes[1, 2].legend()
    axes[1, 2].grid(alpha=0.25)

    for ax in axes[0, :].tolist() + [axes[1, 0]]:
        if ax.has_data():
            ax.set_xlabel(f"x [{axis_unit}]")
            ax.set_ylabel(f"z [{axis_unit}]")
            ax.set_aspect("equal")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def contour_mask(ax, xg: np.ndarray, zg: np.ndarray, mask: np.ndarray, *, color: str, label: str) -> None:
    if np.any(mask) and np.any(~mask):
        ax.contour(xg, zg, mask.astype(float), levels=[0.5], colors=color, linewidths=0.8)
        ax.plot([], [], color=color, linewidth=0.8, label=label)


def render_support_edge_summary(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    ref: dict[str, np.ndarray],
    d_measure: np.ndarray,
    records: list[dict[str, float]],
    support: np.ndarray,
    trusted: np.ndarray,
    active: np.ndarray,
    taper_weight: np.ndarray,
    axis_unit: str = "natural length unit",
) -> None:
    plt = load_pyplot()
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.4), constrained_layout=True)

    vmax = max(float(np.percentile(ref["ntilde_measure"][support], 99.5)), 1.0e-16)
    im0 = axes[0, 0].pcolormesh(xg, zg, ref["ntilde_measure"], shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
    contour_mask(axes[0, 0], xg, zg, support, color="white", label="binary support")
    contour_mask(axes[0, 0], xg, zg, trusted, color="cyan", label="trusted core")
    axes[0, 0].set_title("A transformed measure with support contours")
    axes[0, 0].legend(loc="upper right", fontsize=8)
    fig.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].pcolormesh(xg, zg, taper_weight, shading="auto", cmap="magma", vmin=0.0, vmax=1.0)
    contour_mask(axes[0, 1], xg, zg, active, color="white", label="active evolution window")
    axes[0, 1].set_title("diagnostic tapered-support weight")
    axes[0, 1].legend(loc="upper right", fontsize=8)
    fig.colorbar(im1, ax=axes[0, 1])

    diff = d_measure - ref["ntilde_measure"]
    dmax = max(float(np.percentile(np.abs(diff[support]), 99.0)), 1.0e-16)
    im2 = axes[0, 2].pcolormesh(xg, zg, diff, shading="auto", cmap="coolwarm", vmin=-dmax, vmax=dmax)
    contour_mask(axes[0, 2], xg, zg, support, color="black", label="binary support")
    axes[0, 2].set_title("D - A measure with support contour")
    axes[0, 2].legend(loc="upper right", fontsize=8)
    fig.colorbar(im2, ax=axes[0, 2])

    category = np.zeros_like(taper_weight)
    category[active] = 1.0
    category[support] = 2.0
    category[trusted] = 3.0
    im3 = axes[1, 0].pcolormesh(xg, zg, category, shading="auto", cmap="cividis", vmin=0.0, vmax=3.0)
    axes[1, 0].set_title("mask categories: 1 active, 2 support, 3 trusted")
    fig.colorbar(im3, ax=axes[1, 0])

    times = np.asarray([r["t"] for r in records], dtype=float)
    axes[1, 1].plot(times, [r["measure_rel_l1_support"] for r in records], label="binary support")
    axes[1, 1].plot(times, [r["measure_rel_l1_trusted"] for r in records], label="trusted core")
    axes[1, 1].plot(times, [r["measure_rel_l1_tapered"] for r in records], label="tapered weight")
    axes[1, 1].plot(times, [r["measure_rel_l1_support_edge"] for r in records], label="support edge band")
    axes[1, 1].set_title("measure relative L1 by support expression")
    axes[1, 1].set_xlabel("t")
    axes[1, 1].set_yscale("log")
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].grid(alpha=0.25)

    axes[1, 2].plot(times, [r["disc_min_support"] for r in records], label="disc min support")
    axes[1, 2].plot(times, [r["disc_min_trusted"] for r in records], label="disc min trusted")
    axes[1, 2].plot(times, [r["disc_p01_tapered"] for r in records], label="disc p01 tapered")
    axes[1, 2].set_title("mass-shell discriminant diagnostics")
    axes[1, 2].set_xlabel("t")
    axes[1, 2].legend(fontsize=8)
    axes[1, 2].grid(alpha=0.25)

    for ax in axes.ravel():
        if ax.has_data():
            ax.set_xlabel(f"x [{axis_unit}]")
            ax.set_ylabel(f"z [{axis_unit}]")
            ax.set_aspect("equal" if ax in axes[0, :].tolist() + [axes[1, 0]] else "auto")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    params, physical_scale = build_crossing_params(args)
    if args.physical_optical:
        if args.mp is None:
            args.mp = PLANCK_MASS_EV
        if args.ell_over_planck is not None:
            args.ell = args.ell_over_planck * PLANCK_LENGTH_EV_INV
        elif args.ell is None:
            raise ValueError("physical optical mode requires --ell-over-planck or explicit --ell in eV^-1")
        if args.dt_old_units:
            if physical_scale is None:
                raise ValueError("--dt-old-units requires physical optical mode")
            args.dt = args.dt * physical_scale["old_dimensionless_scale_ev_inv"]
    else:
        if args.mp is None:
            args.mp = 300.0
        if args.ell is None:
            args.ell = 300.0

    x, z, X_grid, Z_grid = make_grid(params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(X_grid, Z_grid, params)
    psi0_norm_before = float(np.sum(np.abs(psi0) ** 2) * dx * dz)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_before = kg_positive_frequency_norm(psi0, psi0_t, dx, dz)
    if args.normalize_probability and args.normalize_kg:
        raise ValueError("Choose at most one of --normalize-probability and --normalize-kg")
    if args.normalize_probability:
        psi0 = psi0 / np.sqrt(max(psi0_norm_before, 1.0e-300))
        psi0_hat = np.fft.fft2(psi0)
        psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    if args.normalize_kg:
        if kg_norm_before <= 0.0:
            raise ValueError(f"KG norm must be positive for positive-frequency normalization, got {kg_norm_before}")
        psi0 = psi0 / np.sqrt(kg_norm_before)
        psi0_hat = np.fft.fft2(psi0)
        psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    psi0_norm_after = float(np.sum(np.abs(psi0) ** 2) * dx * dz)
    kg_norm_after = kg_positive_frequency_norm(psi0, psi0_t, dx, dz)

    init = localized_direct_tilde_coordinate_initial(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        initial_time=args.initial_time,
        probe_dt=args.dt,
        rho_floor=args.rho_floor,
    )
    psi_evolve0, psi_evolve0_hat = evolve_spectral_state(psi0, psi0_hat, omega, args.initial_time)
    diag = init["diagnostics"]
    rho_a0 = real_array(init["rho0"])
    s_phase = unwrap_phase_like(real_array(init["s0"]), real_array(init["s0"]))
    if args.matter_variable_mode == "phase":
        ux, uz = gradient_from_phase(s_phase, active_mask=np.ones_like(rho_a0, dtype=bool), dx=dx, dz=dz)
    else:
        ux = real_array(diag["s_x0"])
        uz = real_array(diag["s_z0"])
    ut = real_array(diag["s_t0"])
    x_field_0 = real_array(diag["X0"])
    measure0 = np.abs(x_field_0) * rho_a0 / (params.m * params.m)

    support = (rho_a0 > args.support_rho_frac * float(np.max(rho_a0))) & (
        measure0 > args.support_measure_frac * float(np.max(measure0))
    )
    trusted = erode_mask_8(support, args.trusted_erosion)
    if not np.any(trusted):
        trusted = support.copy()
    stop_mask = trusted if args.stop_mask == "trusted" else support
    active = dilate_mask_4(support, args.active_dilation)
    projection_mask = {
        "stop": stop_mask,
        "support": support,
        "active": active,
    }[args.matter_projection_mask]
    taper_weight = support_taper_weights(
        rho=rho_a0,
        measure=measure0,
        rho_frac=args.support_rho_frac,
        measure_frac=args.support_measure_frac,
        decades=args.support_taper_decades,
        floor=args.rho_floor,
    )
    taper_weight_active = np.where(active, taper_weight, 0.0)
    support_edge = support & (~trusted)
    if not np.any(support_edge):
        support_edge = support.copy()
    taper_weight_total = max(float(np.sum(taper_weight)), 1.0e-300)
    taper_weight_outside_active_fraction = float(np.sum(taper_weight[~active]) / taper_weight_total)

    metric_cov = sanitize_metric_tail(real_array(init["metric_cov_txz_0"]), active)
    metric_t0 = real_array(init["metric_cov_txz_t0"])
    metric_pprev = sanitize_metric_tail(metric_cov - args.dt * metric_t0, active)
    metric_prev = metric_cov.copy()
    metric_next_linear = sanitize_metric_tail(metric_cov + args.dt * metric_t0, active)

    rho_tilde0 = tilde_density_from_measure(measure0, metric_cov, args.rho_floor)
    rho_tilde0 = np.where(active, rho_tilde0, 0.0)
    cons = conservative_density_from_rho_u(
        rho=rho_tilde0,
        u_x=ux,
        u_z=uz,
        metric_cov_txz=metric_cov,
        mass=params.m,
        branch="negative_frequency",
        u_t_reference=ut,
        rho_floor=args.rho_floor,
    )
    n_cons = np.where(active, cons["n_cons"], 0.0)
    current = coordinate_matter_rhs_covector(
        n_cons=n_cons,
        u_x=ux,
        u_z=uz,
        metric_cov_txz=metric_cov,
        mass=params.m,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=ut,
    )
    for key in ("rho", "u_t", "discriminant"):
        current[key] = np.where(active, current[key], 0.0)
    if args.matter_variable_mode == "phase":
        n_cons, ux, uz, current = phase_current_from_state(
            n_cons=n_cons,
            s_phase=s_phase,
            metric_cov=metric_cov,
            mass=params.m,
            dx=dx,
            dz=dz,
            active_mask=active,
            s_t_reference=ut,
        )
    measure_density = np.where(active, measure0, 0.0)
    if args.matter_variable_mode == "local_time":
        n_cons, current = current_from_measure_and_u(
            measure_density=measure_density,
            u_x=ux,
            u_z=uz,
            metric_cov=metric_cov,
            mass=params.m,
            u_t_reference=ut,
            active_mask=active,
            rho_floor=args.rho_floor,
        )

    records: list[dict[str, float]] = []
    last_raw_y: np.ndarray | None = None
    last_rows: list[dict[str, float]] = []
    last_tensor_rows: list[dict[str, float]] = []
    last_tensor_interface_diag: dict[str, object] = {
        "rows": [],
        "accepted_rows": [],
        "segment_count": 0,
        "candidate_count": 0,
        "solved_count": 0,
        "solved_fraction": 0.0,
        "anchor_count": 0,
        "anchor_fraction": 0.0,
        "direct_tensor_residual_relative": {"count": 0, "abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0},
        "tensor_anchor_mask": np.zeros_like(rho_a0, dtype=bool),
    }
    last_masks: dict[str, np.ndarray] | None = None
    last_extension_diag = empty_metric_extension_diagnostic()
    last_matter_projection_diag = empty_matter_projection_diagnostic()
    last_local_time_diag = empty_local_time_diagnostic()
    stopped_reason = "completed"
    completed_step = 0
    geometry_every = max(int(args.geometry_every), 1)
    metric_extension_every = int(args.metric_extension_every)
    diagnostics_every = int(args.diagnostics_every)
    geom: dict[str, np.ndarray] | None = None
    raw_y: np.ndarray | None = None
    masks: dict[str, np.ndarray] | None = None
    geom_step = -1
    metric_inv_current: np.ndarray | None = None
    metric_det_current: np.ndarray | None = None
    cached_local_time_labels: np.ndarray | None = None
    cached_local_time_candidates: list[tuple[float, float]] | None = None
    cached_local_time_age = 0

    for step in range(args.steps + 1):
        completed_step = step
        t = step * args.dt

        metric_inv_current = None
        metric_det_current = None
        force_endpoint_diagnostics = not bool(args.quick_diagnostics)
        diagnostics_due = (
            (force_endpoint_diagnostics and (step == 0 or step == args.steps))
            or (diagnostics_every > 0 and step % diagnostics_every == 0)
        )
        interface_due = args.interface_every > 0 and (
            step % args.interface_every == 0 or (force_endpoint_diagnostics and step == args.steps)
        )
        metric_extension_due = (
            args.geometry_closure == "restricted_extension"
            and step < args.steps
            and metric_extension_every > 0
            and step % metric_extension_every == 0
        )
        geometry_due = (
            diagnostics_due
            or interface_due
            or metric_extension_due
            or (args.render and step == args.steps)
            or (
                force_endpoint_diagnostics
                and (geom is None or step == 0 or step == args.steps or step % geometry_every == 0)
            )
        )

        if args.matter_variable_mode == "local_time":
            metric_inv_current, metric_det_current = inverse_metric_block(metric_cov)
            n_cons, current = current_from_measure_and_u(
                measure_density=measure_density,
                u_x=ux,
                u_z=uz,
                metric_cov=metric_cov,
                mass=params.m,
                u_t_reference=current.get("u_t", ut),
                active_mask=active,
                rho_floor=args.rho_floor,
                metric_inv=metric_inv_current,
                det_cov=metric_det_current,
            )
            rho_tilde_d = np.maximum(measure_density / np.sqrt(np.maximum(np.abs(metric_det_current), 1.0e-300)), args.rho_floor)
            rho_tilde_d = np.where(active, rho_tilde_d, 0.0)
            current["rho"] = rho_tilde_d
            measure_d = np.where(active, measure_density, 0.0)
        else:
            rho_tilde_d = real_array(current["rho"])
            sqrt_abs_g_d = safe_sqrt_abs_det(metric_cov)
            measure_d = sqrt_abs_g_d * rho_tilde_d

        if geometry_due:
            if step == 0:
                geom = metric_curvature_from_history(metric_pprev, metric_prev, metric_next_linear, args.dt, dx, dz)
            else:
                geom = metric_curvature_from_history(metric_pprev, metric_prev, metric_cov, args.dt, dx, dz)
            raw_y = args.ell * args.ell * geom["R_tilde"]
            masks = domain_masks(raw_y, support, args.weak_y, args.saturation_phi)
            last_raw_y = raw_y
            last_masks = masks
            geom_step = step
        if (diagnostics_due or interface_due or metric_extension_due or args.render) and (
            geom is None or raw_y is None or masks is None
        ):
            raise RuntimeError("internal error: geometry cache was not initialized")

        interface = {"rows": [], "segment_count": 0, "solved_count": 0, "solved_fraction": 0.0}
        tensor_interface = {
            "rows": [],
            "accepted_rows": [],
            "segment_count": 0,
            "candidate_count": 0,
            "solved_count": 0,
            "solved_fraction": 0.0,
            "anchor_count": 0,
            "anchor_fraction": 0.0,
            "direct_tensor_residual_relative": {"count": 0, "abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0},
            "tensor_anchor_mask": np.zeros_like(rho_a0, dtype=bool),
        }
        if interface_due:
            stress, _ = stress_tensor_tilde(
                geom["metric_cov"],
                geom["metric_inv"],
                rho_tilde_d,
                current["u_t"],
                ux,
                uz,
                m=params.m,
            )
            rhs = real_array(stress) / (args.mp * args.mp)
            stress_trace = np.einsum("...ab,...ab->...", geom["metric_inv"], rhs, optimize=True)
            interface = interface_jump_diagnostics(
                raw_y=raw_y,
                metric_inv=geom["metric_inv"],
                rho=measure_d,
                stress_trace=stress_trace,
                x=x,
                z=z,
                ell=args.ell,
                half_width_y=args.half_width_y,
                samples=args.samples,
            )
            last_rows = interface["rows"]
            tensor_interface = full_tensor_interface_diagnostics(
                raw_y=raw_y,
                metric_cov=geom["metric_cov"],
                metric_inv=geom["metric_inv"],
                ricci=geom["ricci"],
                rho=measure_d,
                rhs_tensor=rhs,
                x=x,
                z=z,
                ell=args.ell,
                half_width_y=args.half_width_y,
                samples=args.samples,
                tensor_rel_tol=args.tensor_rel_tol,
                max_iter=args.tensor_max_iter,
                multistart=args.tensor_multistart,
                max_seeds=args.tensor_max_seeds,
                band_factor=args.tensor_anchor_band_factor,
            )
            last_tensor_rows = tensor_interface["rows"]
            last_tensor_interface_diag = tensor_interface

        if diagnostics_due:
            ref = reference_a_snapshot(psi_evolve0, psi_evolve0_hat, omega, x, z, t, params.m, args.rho_floor, active)
            ref_cons = conservative_density_from_rho_u(
                rho=np.where(active, ref["rho_tilde"], 0.0),
                u_x=ref["s_x"],
                u_z=ref["s_z"],
                metric_cov_txz=ref["metric_cov"],
                mass=params.m,
                branch="negative_frequency",
                u_t_reference=ref["s_t"],
                rho_floor=args.rho_floor,
            )
            sat_diag = plateau_tensor_diagnostic(
                raw_y=raw_y,
                metric_cov=geom["metric_cov"],
                metric_inv=geom["metric_inv"],
                rho_tilde=rho_tilde_d,
                u_t=current["u_t"],
                u_x=ux,
                u_z=uz,
                mass=params.m,
                ell=args.ell,
                mp=args.mp,
                mask=masks["saturated"],
            )
            support_count = max(int(np.count_nonzero(support)), 1)
            records.append(
                {
                    "step": int(step),
                    "t": float(t),
                    "geometry_cache_step": int(geom_step),
                    "measure_rel_l1_support": relative_l1(ref["ntilde_measure"], measure_d, support),
                    "measure_max_abs_support": float(np.max(np.abs((measure_d - ref["ntilde_measure"])[support]))),
                    "n_cons_rel_l1_support": relative_l1(ref_cons["n_cons"], n_cons, support),
                    "measure_rel_l1_support_edge": relative_l1(ref["ntilde_measure"], measure_d, support_edge),
                    "n_cons_rel_l1_support_edge": relative_l1(ref_cons["n_cons"], n_cons, support_edge),
                    "disc_min_support_edge": float(np.nanmin(current["discriminant"][support_edge])),
                    "measure_rel_l1_tapered": weighted_relative_l1(ref["ntilde_measure"], measure_d, taper_weight_active),
                    "n_cons_rel_l1_tapered": weighted_relative_l1(ref_cons["n_cons"], n_cons, taper_weight_active),
                    "disc_min_tapered": float(np.nanmin(current["discriminant"][taper_weight_active > 0.0]))
                    if np.any(taper_weight_active > 0.0)
                    else 0.0,
                    "disc_p01_tapered": weighted_quantile(current["discriminant"], taper_weight_active, 0.01),
                    "rho_tilde_min_support": float(np.nanmin(rho_tilde_d[support])),
                    "rho_tilde_max_support": float(np.nanmax(rho_tilde_d[support])),
                    "disc_min_support": float(np.nanmin(current["discriminant"][support])),
                    "measure_rel_l1_trusted": relative_l1(ref["ntilde_measure"], measure_d, trusted),
                    "measure_max_abs_trusted": float(np.max(np.abs((measure_d - ref["ntilde_measure"])[trusted]))),
                    "n_cons_rel_l1_trusted": relative_l1(ref_cons["n_cons"], n_cons, trusted),
                    "rho_tilde_min_trusted": float(np.nanmin(rho_tilde_d[trusted])),
                    "rho_tilde_max_trusted": float(np.nanmax(rho_tilde_d[trusted])),
                    "disc_min_trusted": float(np.nanmin(current["discriminant"][trusted])),
                    "weak_fraction": float(np.count_nonzero(masks["weak"]) / support_count),
                    "transition_fraction": float(np.count_nonzero(masks["transition"]) / support_count),
                    "saturated_fraction": float(np.count_nonzero(masks["saturated"]) / support_count),
                    "trusted_count": int(np.count_nonzero(trusted)),
                    "support_count": int(np.count_nonzero(support)),
                    "support_edge_count": int(np.count_nonzero(support_edge)),
                    "taper_weight_sum": float(np.sum(taper_weight_active)),
                    "taper_edge_count": int(np.count_nonzero((taper_weight > 0.0) & (taper_weight < 1.0))),
                    "taper_weight_outside_active_fraction": taper_weight_outside_active_fraction,
                    "interface_segment_count": int(interface["segment_count"]),
                    "interface_solved_count": int(interface["solved_count"]),
                    "interface_solved_fraction": float(interface["solved_fraction"]),
                    "tensor_interface_segment_count": int(tensor_interface["segment_count"]),
                    "tensor_interface_candidate_count": int(tensor_interface["candidate_count"]),
                    "tensor_interface_solved_count": int(tensor_interface["solved_count"]),
                    "tensor_interface_solved_fraction": float(tensor_interface["solved_fraction"]),
                    "tensor_interface_anchor_count": int(tensor_interface["anchor_count"]),
                    "tensor_interface_anchor_fraction": float(tensor_interface["anchor_fraction"]),
                    "tensor_interface_direct_tensor_residual_relative_p95": float(tensor_interface["direct_tensor_residual_relative"]["p95"]),
                    **sat_diag,
                    **last_extension_diag,
                    **last_matter_projection_diag,
                    **last_local_time_diag,
                }
            )

        disc_min_stop = float(np.nanmin(current["discriminant"][stop_mask]))
        if args.stop_on_negative_discriminant and disc_min_stop < -abs(args.disc_tolerance):
            stopped_reason = f"negative mass-shell discriminant on {args.stop_mask} mask: {disc_min_stop:.6e}"
            break
        if step == args.steps:
            break

        try:
            base_n_cons = n_cons.copy()
            base_ux = ux.copy()
            base_uz = uz.copy()
            base_ut = current["u_t"].copy()
            if args.matter_variable_mode == "phase":
                n_cons, s_phase, ux, uz, current = rk4_phase_matter_step(
                    n_cons=n_cons,
                    s_phase=s_phase,
                    metric_cov=metric_cov,
                    mass=params.m,
                    dx=dx,
                    dz=dz,
                    dt=args.dt,
                    s_t_reference=current["u_t"],
                    active_mask=active,
                )
                last_matter_projection_diag = empty_matter_projection_diagnostic()
                last_local_time_diag = empty_local_time_diagnostic()
            elif args.matter_variable_mode == "local_time":
                if metric_inv_current is None or metric_det_current is None:
                    metric_inv_current, metric_det_current = inverse_metric_block(metric_cov)
                metric_inv_for_clock = metric_inv_current
                metric_det_for_clock = metric_det_current
                a_values, b_values = local_time_candidate_arrays(args.local_time_max_tilt, args.local_time_step)
                if args.local_time_patchwise:
                    reuse_atlas = (
                        args.local_time_atlas_every > 1
                        and cached_local_time_labels is not None
                        and cached_local_time_candidates is not None
                        and cached_local_time_age < args.local_time_atlas_every
                    )
                    if reuse_atlas:
                        labels = cached_local_time_labels
                        candidates = cached_local_time_candidates
                        selector_diag = {
                            **local_time_label_diagnostics(
                                current=current,
                                metric_inv=metric_inv_for_clock,
                                measure=measure_density,
                                active=active,
                                trusted=stop_mask,
                                labels=labels,
                                candidates=candidates,
                                norm_floor=args.local_time_norm_floor,
                            ),
                            "stationary_candidate_raw_count": 0.0,
                            "stationary_candidate_unique_count": 0.0,
                            "tile_size": float(args.local_time_tile_size),
                        }
                        cached_local_time_age += 1
                    else:
                        labels, candidates, selector_diag = select_local_time_patch_labels(
                            current=current,
                            metric_inv=metric_inv_for_clock,
                            measure=measure_density,
                            active=active,
                            trusted=stop_mask,
                            a_values=a_values,
                            b_values=b_values,
                            norm_floor=args.local_time_norm_floor,
                            include_stationary_candidates=args.local_time_stationary_candidates,
                            stationary_candidate_quantization=args.local_time_stationary_quantization,
                            tile_size=args.local_time_tile_size,
                        )
                        cached_local_time_labels = labels
                        cached_local_time_candidates = candidates
                        cached_local_time_age = 1
                    measure_density, n_cons, ux, uz, current, step_diag = rk4_patchwise_local_time_matter_step(
                        measure_density=measure_density,
                        u_x=ux,
                        u_z=uz,
                        metric_cov=metric_cov,
                        metric_inv=metric_inv_for_clock,
                        metric_det=metric_det_for_clock,
                        mass=params.m,
                        dx=dx,
                        dz=dz,
                        dt=args.dt,
                        u_t_reference=current["u_t"],
                        active_mask=active,
                        labels=labels,
                        candidates=candidates,
                        use_bounding_boxes=args.local_time_patch_bboxes,
                        split_components=args.local_time_split_components,
                    )
                    tau_a = 0.0
                    tau_b = 0.0
                    local_diag = {
                        "local_time_tau_a": 0.0,
                        "local_time_tau_b": 0.0,
                        "local_time_score": 0.0,
                        "local_time_admissible_fraction": 1.0 - float(selector_diag["uncovered_trusted_measure_fraction"]),
                        "local_time_atlas_reused": 1.0 if reuse_atlas else 0.0,
                        "local_time_tile_size": float(selector_diag.get("tile_size", args.local_time_tile_size)),
                        "local_time_patch_count": float(selector_diag["patch_count"]),
                        "local_time_candidate_count": float(selector_diag["candidate_count"]),
                        "local_time_stationary_candidate_raw_count": float(selector_diag["stationary_candidate_raw_count"]),
                        "local_time_stationary_candidate_unique_count": float(selector_diag["stationary_candidate_unique_count"]),
                        "local_time_patch_label_edges": float(selector_diag["patch_label_edges"]),
                        "local_time_patch_lab_fraction": float(selector_diag["lab_label_fraction"]),
                        "local_time_uncovered_trusted_measure_fraction": float(selector_diag["uncovered_trusted_measure_fraction"]),
                        "local_time_abs_jtau_p05": float(selector_diag["abs_jtau_p05"]),
                        "local_time_abs_jtau_p50": float(selector_diag["abs_jtau_p50"]),
                        "local_time_abs_jtau_p95": float(selector_diag["abs_jtau_p95"]),
                        "local_time_tau_norm_min": float(selector_diag["tau_norm_min"]),
                        "local_time_tau_norm_p05": float(selector_diag["tau_norm_p05"]),
                        "local_time_fallback_to_lab": 0.0,
                    }
                    local_diag["local_time_patch_bbox_area_fraction"] = float(step_diag.get("patch_bbox_area_fraction", 0.0))
                    local_diag["local_time_patch_bbox_fallback_count"] = float(step_diag.get("patch_bbox_fallback_count", 0.0))
                    local_diag["local_time_patch_component_count"] = float(step_diag.get("patch_components_used", 0.0))
                    local_diag["local_time_n_tau_rel_delta_sum"] = float(step_diag.get("n_tau_rel_delta_sum", 0.0))
                else:
                    tau_a, tau_b, selector_diag = select_local_time_chart(
                        current=current,
                        metric_inv=metric_inv_for_clock,
                        measure=measure_density,
                        trusted=stop_mask,
                        a_values=a_values,
                        b_values=b_values,
                        norm_floor=args.local_time_norm_floor,
                        admissible_fraction_floor=args.local_time_admissible_fraction_floor,
                    )
                    measure_density, n_cons, ux, uz, current, step_diag = rk4_local_time_matter_step(
                        measure_density=measure_density,
                        u_x=ux,
                        u_z=uz,
                        metric_cov=metric_cov,
                        metric_inv=metric_inv_for_clock,
                        metric_det=metric_det_for_clock,
                        mass=params.m,
                        dx=dx,
                        dz=dz,
                        dt=args.dt,
                        u_t_reference=current["u_t"],
                        active_mask=active,
                        tau_a=tau_a,
                        tau_b=tau_b,
                    )
                    local_diag = {
                        "local_time_tau_a": float(tau_a),
                        "local_time_tau_b": float(tau_b),
                        "local_time_score": float(selector_diag["score"]),
                        "local_time_admissible_fraction": float(selector_diag["admissible_fraction"]),
                        "local_time_patch_count": 1.0,
                        "local_time_candidate_count": float(a_values.size * b_values.size),
                        "local_time_stationary_candidate_raw_count": 0.0,
                        "local_time_stationary_candidate_unique_count": 0.0,
                        "local_time_patch_label_edges": 0.0,
                        "local_time_patch_lab_fraction": 1.0 if (tau_a == 0.0 and tau_b == 0.0) else 0.0,
                        "local_time_uncovered_trusted_measure_fraction": 1.0 - float(selector_diag["admissible_fraction"]),
                        "local_time_abs_jtau_p05": float(selector_diag["abs_jtau_p05"]),
                        "local_time_abs_jtau_p50": float(selector_diag["abs_jtau_p50"]),
                        "local_time_abs_jtau_p95": float(selector_diag["abs_jtau_p95"]),
                        "local_time_tau_norm_min": float(selector_diag["tau_norm_min"]),
                        "local_time_tau_norm_p05": float(selector_diag["tau_norm_p05"]),
                        "local_time_fallback_to_lab": float(selector_diag.get("fallback_to_lab", 0.0)),
                        "local_time_patch_bbox_area_fraction": 1.0,
                        "local_time_patch_bbox_fallback_count": 0.0,
                        "local_time_patch_component_count": 1.0,
                        "local_time_n_tau_rel_delta_sum": float(step_diag.get("n_tau_rel_delta", 0.0)),
                    }
                last_matter_projection_diag = empty_matter_projection_diagnostic()
                last_local_time_diag = {
                    **local_diag,
                    "local_time_measure_rel_delta": float(step_diag.get("measure_rel_delta", 0.0)),
                }
            else:
                n_cons, ux, uz, current = rk4_matter_step(
                    n_cons=n_cons,
                    u_x=ux,
                    u_z=uz,
                    metric_cov=metric_cov,
                    mass=params.m,
                    dx=dx,
                    dz=dz,
                    dt=args.dt,
                    u_t_reference=current["u_t"],
                    active_mask=active,
                )
                if args.matter_projection:
                    n_cons, ux, uz, current, last_matter_projection_diag = project_covector_matter_update(
                        base_n_cons=base_n_cons,
                        base_u_x=base_ux,
                        base_u_z=base_uz,
                        n_cons=n_cons,
                        u_x=ux,
                        u_z=uz,
                        metric_cov=metric_cov,
                        mass=params.m,
                        dx=dx,
                        dz=dz,
                        u_t_reference=base_ut,
                        projection_mask=projection_mask,
                        disc_floor=args.matter_projection_disc_floor,
                        rho_floor=args.matter_projection_rho_floor,
                        max_blends=args.matter_projection_max_blends,
                        positivity_limiter=args.matter_positivity_limiter,
                        limiter_rho_floor=args.matter_limiter_rho_floor,
                    )
                    if last_matter_projection_diag["matter_projection_final_bad_count"] > 0.0:
                        completed_step = step + 1
                        stopped_reason = (
                            "matter projection failed on "
                            f"{args.matter_projection_mask} mask: "
                            f"{int(last_matter_projection_diag['matter_projection_final_bad_count'])} bad cells"
                        )
                        break
                else:
                    last_matter_projection_diag = empty_matter_projection_diagnostic()
        except FloatingPointError as exc:
            stopped_reason = str(exc)
            break

        if args.geometry_closure == "inert_tridomain":
            # This is not a damping term.  It is the explicit three-domain closure:
            # weak/saturated bulks have no unique local metric evolution after
            # high-derivative terms are dropped, so the first executable closure
            # keeps the current bulk geometry and lets the matter equations move.
            metric_pprev = metric_prev
            metric_prev = metric_cov
            metric_cov = metric_cov
        elif args.geometry_closure == "restricted_extension" and metric_extension_due:
            # Saturated-bulk representative rule:
            # keep weak/transition/interface data fixed, and extend the metric
            # through saturated interiors by a minimal discrete-curvature
            # representative.  This selects a representative only where f_R has
            # already removed the local highest-order metric equation.
            physical_anchor_mask = support & (~masks["saturated"])
            tensor_anchor_mask = np.asarray(last_tensor_interface_diag["tensor_anchor_mask"], dtype=bool)
            flat_anchor_mask = low_matter_flat_anchor_mask(
                rho=rho_a0,
                measure=measure0,
                active_mask=active,
                rho_frac=args.flat_anchor_rho_frac,
                measure_frac=args.flat_anchor_measure_frac,
                boundary_layers=args.flat_anchor_boundary_layers,
            )
            anchor_mask = physical_anchor_mask | flat_anchor_mask | tensor_anchor_mask
            extended_metric, extension_diag = restricted_harmonic_metric_extension(
                metric_cov=metric_cov,
                solve_mask=masks["saturated"],
                anchor_mask=anchor_mask,
                flat_anchor_mask=flat_anchor_mask,
                active_mask=active,
                dx=dx,
                dz=dz,
                iterations=args.metric_extension_iterations,
                anchor_weight=args.metric_extension_anchor_weight,
                enforce_lorentz=args.metric_extension_enforce_lorentz,
                variable_mode=args.metric_extension_variable_mode,
                adm_update_fields=args.metric_extension_adm_update_fields,
                max_rel_change=args.metric_extension_max_rel_change,
            )
            extended_metric, updated_current, admissible_diag = select_admissible_metric_update(
                base_metric=metric_cov,
                target_metric=extended_metric,
                n_cons=n_cons,
                u_x=ux,
                u_z=uz,
                mass=params.m,
                dx=dx,
                dz=dz,
                u_t_reference=current["u_t"],
                active_mask=active,
                check_mask=stop_mask,
                disc_tolerance=args.disc_tolerance,
                rho_floor=args.rho_floor,
                max_halvings=args.metric_extension_max_halvings,
            )
            current = updated_current
            extension_diag.update(admissible_diag)
            last_extension_diag = extension_diag
            metric_pprev = metric_prev
            metric_prev = metric_cov
            metric_cov = extended_metric
        elif args.geometry_closure == "restricted_extension":
            metric_pprev = metric_prev
            metric_prev = metric_cov
            metric_cov = metric_cov
        else:
            raise ValueError(args.geometry_closure)

    t_completed = completed_step * args.dt
    need_final_reference = (not args.quick_diagnostics) or args.render or (not args.skip_fields_npz)
    final_ref = (
        reference_a_snapshot(psi_evolve0, psi_evolve0_hat, omega, x, z, t_completed, params.m, args.rho_floor, active)
        if need_final_reference
        else None
    )
    final_rho_tilde_d = real_array(current["rho"])
    final_measure_d = safe_sqrt_abs_det(metric_cov) * final_rho_tilde_d
    if args.matter_variable_mode == "local_time":
        final_measure_d = np.where(active, measure_density, 0.0)
        final_rho_tilde_d = np.where(active, tilde_density_from_measure(final_measure_d, metric_cov, args.rho_floor), 0.0)
        current["rho"] = final_rho_tilde_d
    if final_ref is None:
        final_reference_metrics = {
            "measure_rel_l1_support": None,
            "measure_max_abs_support": None,
            "measure_rel_l1_support_edge": None,
            "measure_rel_l1_tapered": None,
            "measure_rel_l1_trusted": None,
            "measure_max_abs_trusted": None,
        }
    else:
        final_reference_metrics = {
            "measure_rel_l1_support": relative_l1(final_ref["ntilde_measure"], final_measure_d, support),
            "measure_max_abs_support": float(np.max(np.abs((final_measure_d - final_ref["ntilde_measure"])[support]))),
            "measure_rel_l1_support_edge": relative_l1(final_ref["ntilde_measure"], final_measure_d, support_edge),
            "measure_rel_l1_tapered": weighted_relative_l1(final_ref["ntilde_measure"], final_measure_d, taper_weight_active),
            "measure_rel_l1_trusted": relative_l1(final_ref["ntilde_measure"], final_measure_d, trusted),
            "measure_max_abs_trusted": float(np.max(np.abs((final_measure_d - final_ref["ntilde_measure"])[trusted]))),
        }
    terminal_record = {
        "step": int(completed_step),
        "t": float(t_completed),
        "geometry_cache_step": int(geom_step),
        "last_diagnostic_step": int(records[-1]["step"]) if records else -1,
        **final_reference_metrics,
        "rho_tilde_min_support": float(np.nanmin(final_rho_tilde_d[support])),
        "rho_tilde_max_support": float(np.nanmax(final_rho_tilde_d[support])),
        "disc_min_support": float(np.nanmin(current["discriminant"][support])),
        "rho_tilde_min_trusted": float(np.nanmin(final_rho_tilde_d[trusted])),
        "rho_tilde_max_trusted": float(np.nanmax(final_rho_tilde_d[trusted])),
        "disc_min_trusted": float(np.nanmin(current["discriminant"][trusted])),
        "disc_min_support_edge": float(np.nanmin(current["discriminant"][support_edge])),
        "metric_extension_admissible_blend": float(last_extension_diag.get("metric_extension_admissible_blend", 0.0)),
        "metric_extension_solve_fraction": float(last_extension_diag.get("metric_extension_solve_fraction", 0.0)),
        "trusted_count": int(np.count_nonzero(trusted)),
        "support_count": int(np.count_nonzero(support)),
        "active_count": int(np.count_nonzero(active)),
        **last_matter_projection_diag,
        **last_local_time_diag,
    }
    if last_masks is not None:
        support_count = max(int(np.count_nonzero(support)), 1)
        terminal_record.update(
            {
                "weak_fraction": float(np.count_nonzero(last_masks["weak"]) / support_count),
                "transition_fraction": float(np.count_nonzero(last_masks["transition"]) / support_count),
                "saturated_fraction": float(np.count_nonzero(last_masks["saturated"]) / support_count),
            }
        )
    else:
        terminal_record.update(
            {
                "weak_fraction": None,
                "transition_fraction": None,
                "saturated_fraction": None,
            }
        )

    fig_path = out / f"d_tridomain_dynamics_{args.geometry_closure}_ell{args.ell:g}_n{args.resolution}_t{t_completed:g}.png"
    if args.physical_optical:
        plot_factor = HBAR_C_EV_M * 1.0e6
        x_plot = x * plot_factor
        z_plot = z * plot_factor
        rows_plot = scale_interface_rows(last_rows, plot_factor)
        axis_unit = "um"
    else:
        x_plot = x
        z_plot = z
        rows_plot = last_rows
        axis_unit = "code"
    edge_fig_path = out / "support_edge_tapered_diagnostics.png"
    if args.render:
        if final_ref is None:
            final_ref = reference_a_snapshot(
                psi_evolve0, psi_evolve0_hat, omega, x, z, t_completed, params.m, args.rho_floor, active
            )
        if last_raw_y is None or last_masks is None:
            geom_final = metric_curvature_from_history(metric_pprev, metric_prev, metric_cov, args.dt, dx, dz)
            last_raw_y = args.ell * args.ell * geom_final["R_tilde"]
            last_masks = domain_masks(last_raw_y, support, args.weak_y, args.saturation_phi)
        render_summary(
            fig_path,
            x=x_plot,
            z=z_plot,
            ref=final_ref,
            d_measure=final_measure_d,
            rho_tilde_d=final_rho_tilde_d,
            raw_y=last_raw_y,
            masks=last_masks,
            interface_rows=rows_plot,
            records=records,
            support=support,
            axis_unit=axis_unit,
        )
        render_support_edge_summary(
            edge_fig_path,
            x=x_plot,
            z=z_plot,
            ref=final_ref,
            d_measure=final_measure_d,
            records=records,
            support=support,
            trusted=trusted,
            active=active,
            taper_weight=taper_weight,
            axis_unit=axis_unit,
        )

    npz_path = out / "fields_final.npz"
    if not args.skip_fields_npz:
        if last_raw_y is None:
            geom_final = metric_curvature_from_history(metric_pprev, metric_prev, metric_cov, args.dt, dx, dz)
            last_raw_y = args.ell * args.ell * geom_final["R_tilde"]
        np.savez_compressed(
            npz_path,
            x=x,
            z=z,
            support=support,
            active=active,
            trusted=trusted,
            taper_weight=taper_weight,
            rho_A=final_ref["rho"],
            X_A=final_ref["X"],
            ntilde_measure_A=final_ref["ntilde_measure"],
            rho_tilde_A=final_ref["rho_tilde"],
            ntilde_measure_D=final_measure_d,
            rho_tilde_D=final_rho_tilde_d,
            measure_density_D=measure_density,
            n_cons_D=n_cons,
            u_t_D=current["u_t"],
            u_x_D=ux,
            u_z_D=uz,
            s_phase_D=s_phase,
            metric_cov_D=metric_cov,
            raw_y=np.zeros_like(rho_a0) if last_raw_y is None else last_raw_y,
        )

    summary = {
        "params": {
            "branch": "D",
            "scope": "Corrected three-domain D dynamics prototype for the 2+1d crossing Gaussian scene. Matter evolves with the tilde-metric mass shell and conservation law, initialized with sqrt(|gtilde|)rho_tilde=|X|rho_A/m^2. Geometry closure is explicit: inert_tridomain keeps the selected bulk representative fixed; restricted_extension selects a saturated-bulk representative by constrained minimal-curvature extension while preserving weak/transition/interface data. This is still a rule-closed representative solve, not yet a unique full tensor metric solve in the underdetermined saturated bulk.",
            "geometry_closure": args.geometry_closure,
            "matter_variable_mode": args.matter_variable_mode,
            "local_time_patchwise": bool(args.local_time_patchwise),
            "local_time_atlas_every": args.local_time_atlas_every,
            "local_time_tile_size": args.local_time_tile_size,
            "local_time_max_tilt": args.local_time_max_tilt,
            "local_time_step": args.local_time_step,
            "local_time_stationary_candidates": bool(args.local_time_stationary_candidates),
            "local_time_stationary_quantization": args.local_time_stationary_quantization,
            "local_time_patch_bboxes": bool(args.local_time_patch_bboxes),
            "local_time_split_components": bool(args.local_time_split_components),
            "local_time_norm_floor": args.local_time_norm_floor,
            "local_time_admissible_fraction_floor": args.local_time_admissible_fraction_floor,
            "matter_projection": bool(args.matter_projection),
            "matter_projection_mask": args.matter_projection_mask,
            "matter_projection_disc_floor": args.matter_projection_disc_floor,
            "matter_projection_rho_floor": args.matter_projection_rho_floor,
            "matter_projection_max_blends": args.matter_projection_max_blends,
            "matter_positivity_limiter": bool(args.matter_positivity_limiter),
            "matter_limiter_rho_floor": args.matter_limiter_rho_floor,
            "metric_extension_iterations": args.metric_extension_iterations,
            "metric_extension_anchor_weight": args.metric_extension_anchor_weight,
            "metric_extension_enforce_lorentz": bool(args.metric_extension_enforce_lorentz),
            "metric_extension_variable_mode": args.metric_extension_variable_mode,
            "metric_extension_adm_update_fields": args.metric_extension_adm_update_fields,
            "metric_extension_max_rel_change": args.metric_extension_max_rel_change,
            "metric_extension_max_halvings": args.metric_extension_max_halvings,
            "flat_anchor_boundary_layers": args.flat_anchor_boundary_layers,
            "flat_anchor_rho_frac": args.flat_anchor_rho_frac,
            "flat_anchor_measure_frac": args.flat_anchor_measure_frac,
            "geometry_every": geometry_every,
            "metric_extension_every": metric_extension_every,
            "diagnostics_every": diagnostics_every,
            "quick_diagnostics": bool(args.quick_diagnostics),
            "skip_fields_npz": bool(args.skip_fields_npz),
            "render": bool(args.render),
            "units": "hbar=c=1, energies in eV and lengths/times in eV^-1" if args.physical_optical else "dimensionless code units",
            "physical_optical": bool(args.physical_optical),
            "physical_scale": physical_scale,
            "psi0_norm_before": psi0_norm_before,
            "psi0_norm_after": psi0_norm_after,
            "kg_norm_before": kg_norm_before,
            "kg_norm_after": kg_norm_after,
            "normalize_probability": bool(args.normalize_probability),
            "normalize_kg": bool(args.normalize_kg),
            "ell": args.ell,
            "ell_over_planck_length": None if args.ell is None else float(args.ell / PLANCK_LENGTH_EV_INV),
            "mp": args.mp,
            "plateau_scale_mp2_over_ell2": float((args.mp * args.mp) / (args.ell * args.ell)),
            "resolution": args.resolution,
            "dt": args.dt,
            "dt_fs": ev_inv_to_fs(args.dt) if args.physical_optical else None,
            "initial_time": args.initial_time,
            "initial_time_fs": ev_inv_to_fs(args.initial_time) if args.physical_optical else None,
            "initial_time_old_units": (
                None
                if (not args.physical_optical or physical_scale is None)
                else float(args.initial_time / physical_scale["old_dimensionless_scale_ev_inv"])
            ),
            "steps_requested": args.steps,
            "steps_completed": completed_step,
            "t_requested": args.steps * args.dt,
            "t_completed": t_completed,
            "t_completed_fs": ev_inv_to_fs(t_completed) if args.physical_optical else None,
            "rho_floor": args.rho_floor,
            "support_rho_frac": args.support_rho_frac,
            "support_measure_frac": args.support_measure_frac,
            "active_dilation": args.active_dilation,
            "trusted_erosion": args.trusted_erosion,
            "stop_mask": args.stop_mask,
            "support_taper_decades": args.support_taper_decades,
            "taper_weight_outside_active_fraction": taper_weight_outside_active_fraction,
            "weak_y": args.weak_y,
            "saturation_phi": args.saturation_phi,
            "interface_every": args.interface_every,
            "stopped_reason": stopped_reason,
        },
        "definitions": {
            "ntilde_measure": "sqrt(|tilde g|) tilde rho.  A-reference value is |X_A| rho_A/m^2 on flat original g.",
            "n_cons": "coordinate conserved density sqrt(|tilde g|) tilde rho tilde g^{t nu} u_nu.",
            "matter_variable_mode": "covector evolves u_x,u_z directly; phase evolves S and reconstructs u_i=partial_i S each step, which preserves phase-gradient compatibility better.",
            "local_time_mode": "Experimental matter integrator: evolves the positive transformed measure density in selected linear time charts tau=t+a x+b z, then reports fields back on the lab-t slice.  Patchwise mode selects a chart per cell/patch.  This avoids using lab n_cons/j^t as the primary rho recovery variable.",
            "local_time_measure_rel_delta": "Relative change of sqrt(|gtilde|)rho_tilde across one local-time patchwise step. This is a diagnostic measure, not the conserved density in a tilted tau chart.",
            "local_time_n_tau_rel_delta_sum": "Patch-summed relative change of the actual coordinate conserved density sqrt(|gtilde|)rho_tilde j^tau during the local tau update.",
            "matter_projection": "Optional local constrained-integrator projection after an explicit covector RK step: cells in projection_mask that violate disc>=disc_floor, rho>=rho_floor, or finiteness are blended back toward the previous step until algebraically admissible.",
            "matter_positivity_limiter": "Optional last-stage positivity limiter for projection cells whose only remaining violation is tiny negative reconstructed rho; it resets rho to matter_limiter_rho_floor by rebuilding n_cons from the local j_t unit factor.",
            "raw_y": "ell^2 times the scalar curvature R_tilde computed from the current selected tilde metric history.",
            "weak_fraction": "fraction of initial support where |raw_y| <= weak_y.",
            "saturated_fraction": "fraction of initial support where f_R=sech^2(raw_y) <= saturation_phi.",
            "transition_fraction": "support points that are neither weak nor saturated.",
            "plateau_relative_p95": "p95 of ||-0.5 f g - T/M_P^2||/(||-0.5 f g||+||T/M_P^2||) in the saturated mask.",
            "trusted": "support eroded by trusted_erosion grid cells with 8-neighbor connectivity.  It separates arbitrary support-threshold edge artifacts from the interior support diagnostics.",
            "tapered_support_weight": "diagnostic-only smooth weight: each rho/measure threshold is softened over support_taper_decades below the binary support threshold; weights do not alter evolution equations or sources.",
            "support_edge": "binary support minus trusted core.  If the erosion is zero or empty this falls back to support.",
            "restricted_extension": "Representative metric rule for saturated bulk only: hold weak/transition/interface cells fixed, then minimize a discrete metric-gradient energy inside saturated cells, with optional anchoring to the prior disformal/inert representative and Lorentzian fallback.",
            "metric_extension_solve_fraction": "active-mask fraction where the restricted_extension rule updated the saturated-bulk representative.",
            "metric_extension_anchor_count": "number of active cells that can serve as physical Dirichlet anchors for saturated-bulk extension; by default these are support cells outside the saturated mask, i.e. weak or transition cells.",
            "metric_extension_flat_anchor_count": "number of low-matter active-boundary cells fixed to flat tilde metric as an asymptotically-flat/vacuum-equivalence boundary condition.",
            "metric_extension_orphan_component_count": "number of saturated connected components that did not touch a physical anchor and were therefore left on the prior representative rather than being driven by the artificial support boundary.",
            "metric_extension_change_rel_p95": "p95 of ||g_extended-g_old||/||g_old|| on the extension solve mask.",
            "metric_extension_roughness_before_after_p95": "p95 of the metric spatial-gradient norm on the extension solve mask before/after restricted extension.",
            "metric_extension_variable_mode": "covariant directly extends g_munu; adm extends log lapse, shift, and Cholesky spatial metric variables before reconstructing g_munu.",
            "metric_extension_adm_update_fields": "ADM-only field subset updated by restricted_extension: all updates lapse, shift, and spatial metric Cholesky variables; lapse updates only log lapse; lapse_shift updates log lapse and shift while holding the spatial metric fixed.",
            "metric_extension_adm_valid_fraction": "active-mask fraction where the current metric could be decomposed into real ADM variables with positive lapse and positive spatial metric.",
            "metric_extension_adm_update_component_count": "number of ADM components updated by the selected metric_extension_adm_update_fields closure.",
            "metric_extension_rel_change_cap_blend": "extra cap applied before admissibility checks so the proposed representative update has bounded relative metric change on the solve mask; one means no cap was needed.",
            "metric_extension_admissible_blend": "accepted fraction of the proposed restricted extension after Lorentz/signature, mass-shell discriminant, and positive-rho checks; zero means the prior representative was retained.",
            "geometry_every": "Production scheduling knob: recompute the expensive Ricci/scalar-curvature geometry cache every N matter microsteps, and always at step 0/final or when diagnostics/interface/metric-extension need fresh geometry.",
            "metric_extension_every": "Production scheduling knob: apply the saturated-bulk representative extension every N matter microsteps.  Zero disables extension refresh and keeps the previous representative between geometry checks.",
            "diagnostics_every": "Production scheduling knob: compute A-reference comparison, plateau algebraic residuals, and record rows every N matter microsteps.  Zero records only step 0 and final.",
            "quick_diagnostics": "Debug-speed mode: do not force expensive endpoint curvature/interface/plateau diagnostics.  It leaves the evolution equations unchanged and still reports algebraic final matter diagnostics.",
            "skip_fields_npz": "Debug-speed mode: skip writing the compressed final field archive.  This changes only output, not the evolved fields.",
            "important_caveat": "The bare rho_A=rho_tilde error is fixed here.  The remaining caveat is not a coding shortcut but a theoretical closure choice: saturated bulk f_R->0 does not uniquely evolve the full metric unless an internal/boundary representative is specified.",
        },
        "final": terminal_record,
        "last_diagnostic_record": records[-1] if records else None,
        "records": records,
        "files": {
            "figure": str(fig_path.resolve()) if args.render else None,
            "support_edge_figure": str(edge_fig_path.resolve()) if args.render else None,
            "fields_npz": None if args.skip_fields_npz else str(npz_path.resolve()),
            "summary": str((out / "summary.json").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--geometry-closure", choices=["inert_tridomain", "restricted_extension"], default="inert_tridomain")
    parser.add_argument(
        "--matter-variable-mode",
        choices=["covector", "phase", "local_time"],
        default="covector",
        help="Matter variables: covector evolves lab n_cons,u_x,u_z; phase evolves lab n_cons,S; local_time evolves positive measure density in a selected tau=t+a x+b z chart and reports back on lab-t slices.",
    )
    parser.add_argument(
        "--local-time-max-tilt",
        type=float,
        default=5.0,
        help="Maximum |a|,|b| searched for local_time chart tau=t+a x+b z.",
    )
    parser.add_argument(
        "--local-time-patchwise",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="In local_time mode, select a local time chart per cell/patch instead of using one global tau for the whole active mask.",
    )
    parser.add_argument(
        "--local-time-step",
        type=float,
        default=0.5,
        help="Coarse grid step in a,b for selecting the local_time chart.",
    )
    parser.add_argument(
        "--local-time-stationary-candidates",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Patchwise local_time only: add analytic per-cell stationary candidates for tau=t+a x+b z, so narrow time cones are not skipped by the coarse a,b grid.",
    )
    parser.add_argument(
        "--local-time-stationary-quantization",
        type=float,
        default=0.1,
        help="Quantization used to merge nearby analytic stationary local-time candidates. This is a numerical atlas compression, not a physical parameter.",
    )
    parser.add_argument(
        "--local-time-patch-bboxes",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Patchwise local_time only: evaluate each patch on its one-cell-halo bounding box instead of the whole grid. This is a performance optimization preserving the same centered-difference stencil away from global boundaries.",
    )
    parser.add_argument(
        "--local-time-atlas-every",
        type=int,
        default=1,
        help="Patchwise local_time only: recompute the local chart atlas every N matter steps and reuse it in between. Values >1 are experimental and must keep coverage diagnostics acceptable.",
    )
    parser.add_argument(
        "--local-time-tile-size",
        type=int,
        default=1,
        help="Patchwise local_time only: select one local time chart per tile_size x tile_size block. A value >1 is an experimental atlas coarsening that must be validated by coverage and mass-shell diagnostics.",
    )
    parser.add_argument(
        "--local-time-split-components",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Patchwise local_time only: split each chart label into connected components before building bounding boxes. This may reduce stencil area but can be slower in pure Python when many tiny components are present.",
    )
    parser.add_argument(
        "--local-time-norm-floor",
        type=float,
        default=-1.0e300,
        help="Minimum gtilde^{mu nu} tau_mu tau_nu allowed when selecting a local_time chart. The default does not impose a timelike-only filter.",
    )
    parser.add_argument(
        "--local-time-admissible-fraction-floor",
        type=float,
        default=0.0,
        help="If the selected local chart covers less than this transformed-measure fraction of the stop mask, fall back to lab time.",
    )
    parser.add_argument(
        "--matter-projection",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="After covector RK steps, locally blend cells back toward the previous state until mass-shell discriminant and positive density constraints are satisfied.",
    )
    parser.add_argument(
        "--matter-projection-mask",
        choices=["stop", "support", "active"],
        default="stop",
        help="Mask on which matter projection enforces algebraic admissibility.",
    )
    parser.add_argument(
        "--matter-projection-disc-floor",
        type=float,
        default=0.0,
        help="Minimum allowed mass-shell discriminant on the projection mask.",
    )
    parser.add_argument(
        "--matter-projection-rho-floor",
        type=float,
        default=0.0,
        help="Minimum allowed reconstructed tilde rho on the projection mask.",
    )
    parser.add_argument(
        "--matter-projection-max-blends",
        type=int,
        default=16,
        help="Maximum local halvings toward the previous matter state during matter projection.",
    )
    parser.add_argument(
        "--matter-positivity-limiter",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="After projection backtracking, reset remaining rho-only negative cells to a small nonnegative rho floor by rebuilding n_cons.",
    )
    parser.add_argument(
        "--matter-limiter-rho-floor",
        type=float,
        default=0.0,
        help="rho floor used by --matter-positivity-limiter.",
    )
    parser.add_argument(
        "--metric-extension-iterations",
        type=int,
        default=40,
        help="Jacobi-style iterations for the saturated-bulk restricted harmonic metric extension.",
    )
    parser.add_argument(
        "--metric-extension-anchor-weight",
        type=float,
        default=0.05,
        help="Small positive weight tying the saturated-bulk extension to the prior representative.",
    )
    parser.add_argument(
        "--metric-extension-enforce-lorentz",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fallback to the prior local representative if the extension violates Lorentz signature.",
    )
    parser.add_argument(
        "--metric-extension-variable-mode",
        choices=["covariant", "adm"],
        default="covariant",
        help="Variables used for restricted metric extension. adm is usually gentler than extending g_munu directly.",
    )
    parser.add_argument(
        "--metric-extension-adm-update-fields",
        choices=["all", "lapse", "lapse_shift"],
        default="all",
        help="When --metric-extension-variable-mode=adm, restrict which ADM fields are updated by the saturated-bulk extension.",
    )
    parser.add_argument(
        "--metric-extension-max-rel-change",
        type=float,
        default=0.0,
        help="Cap the proposed restricted_extension relative metric update on the solve mask before admissibility checks. Set <=0 to disable.",
    )
    parser.add_argument(
        "--metric-extension-max-halvings",
        type=int,
        default=8,
        help="Backtracking halvings used to keep restricted_extension inside the admissible Lorentz/mass-shell/positive-rho set.",
    )
    parser.add_argument(
        "--flat-anchor-boundary-layers",
        type=int,
        default=2,
        help="Low-matter active-boundary layers fixed to flat metric as an asymptotically-flat anchor for restricted_extension. Set 0 to disable.",
    )
    parser.add_argument(
        "--flat-anchor-rho-frac",
        type=float,
        default=1.0e-3,
        help="rho/rho_max threshold for low-matter flat boundary anchors.",
    )
    parser.add_argument(
        "--flat-anchor-measure-frac",
        type=float,
        default=1.0e-3,
        help="transformed-measure/max threshold for low-matter flat boundary anchors.",
    )
    parser.add_argument(
        "--geometry-every",
        type=int,
        default=1,
        help="Recompute expensive Ricci/scalar-curvature geometry every N matter steps. Forced at step 0/final and when diagnostics/interface/extension run.",
    )
    parser.add_argument(
        "--metric-extension-every",
        type=int,
        default=1,
        help="Apply restricted saturated-bulk metric extension every N matter steps. Set 0 to keep the representative fixed between geometry checks.",
    )
    parser.add_argument(
        "--diagnostics-every",
        type=int,
        default=1,
        help="Compute A-reference comparison and plateau diagnostics every N matter steps. Set 0 to record only step 0 and final.",
    )
    parser.add_argument(
        "--render",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Render summary figures. Disable for timing scans to reduce memory and runtime.",
    )
    parser.add_argument(
        "--fast-profile",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Use a production timing profile with sparse diagnostics, sparse interface checks, and no rendering unless explicitly re-enabled.",
    )
    parser.add_argument(
        "--quick-diagnostics",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Debug-speed mode: do not force endpoint curvature/interface/plateau diagnostics. Evolution equations are unchanged.",
    )
    parser.add_argument(
        "--skip-fields-npz",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Debug-speed mode: skip writing the compressed final fields archive.",
    )
    parser.add_argument("--ell", type=float, default=None, help="D-branch length scale. In physical mode this is in eV^-1.")
    parser.add_argument("--ell-over-planck", type=float, default=None, help="Set ell = this factor times 1/M_P in physical mode.")
    parser.add_argument("--mp", type=float, default=None, help="Planck mass parameter. In physical mode default is the physical M_P in eV.")
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument("--dt", type=float, default=1.0e-6)
    parser.add_argument(
        "--dt-old-units",
        action="store_true",
        help="In physical optical mode, interpret --dt as old dimensionless dt and scale it by old_k0/k0_physical.",
    )
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument(
        "--initial-time",
        type=float,
        default=0.0,
        help="A-reference time used to build the D initial data. In physical mode this is eV^-1.",
    )
    parser.add_argument(
        "--initial-time-old-units",
        action="store_true",
        help="In physical optical mode, interpret --initial-time as old dimensionless time and scale it by old_k0/k0_physical.",
    )
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument(
        "--support-taper-decades",
        type=float,
        default=1.0,
        help="Diagnostic-only smooth support edge width in decades below the binary rho/measure thresholds.",
    )
    parser.add_argument("--active-dilation", type=int, default=0)
    parser.add_argument("--trusted-erosion", type=int, default=0)
    parser.add_argument("--stop-mask", choices=["support", "trusted"], default="support")
    parser.add_argument("--weak-y", type=float, default=0.1)
    parser.add_argument("--saturation-phi", type=float, default=1.0e-3)
    parser.add_argument("--interface-every", type=int, default=10)
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=61)
    parser.add_argument(
        "--tensor-rel-tol",
        type=float,
        default=0.1,
        help="Full tensor interface residual tolerance used to accept a tensor-matched interface segment.",
    )
    parser.add_argument(
        "--tensor-max-iter",
        type=int,
        default=30,
        help="Maximum Gauss-Newton refinement iterations for the full tensor interface solve.",
    )
    parser.add_argument(
        "--tensor-max-seeds",
        type=int,
        default=24,
        help="Maximum multistart seeds for the full tensor interface solve.",
    )
    parser.add_argument(
        "--tensor-multistart",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use multistart least-squares refinement for the full tensor interface solve.",
    )
    parser.add_argument(
        "--tensor-anchor-band-factor",
        type=float,
        default=1.0,
        help="Band width in grid steps used to rasterize accepted tensor-matched interface rows as saturated-bulk anchors.",
    )
    parser.add_argument("--disc-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--stop-on-negative-discriminant", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--physical-optical", action="store_true", help="Use 1550nm-style physical optical packet in eV natural units.")
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1, help="Use m = ratio * omega0 with omega0^2=k0^2+m^2.")
    parser.add_argument("--normalize-probability", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    if args.fast_profile:
        if args.render:
            args.render = False
        args.geometry_every = max(int(args.geometry_every), 8)
        if args.geometry_closure == "restricted_extension":
            if args.metric_extension_every > 0:
                args.metric_extension_every = max(int(args.metric_extension_every), 8)
        else:
            args.metric_extension_every = 0
        if args.diagnostics_every > 0:
            args.diagnostics_every = max(int(args.diagnostics_every), 8)
        if args.interface_every > 0:
            args.interface_every = max(int(args.interface_every), 8)
        args.tensor_multistart = False
        args.tensor_max_seeds = min(int(args.tensor_max_seeds), 8)
        args.tensor_max_iter = min(int(args.tensor_max_iter), 12)
    if args.physical_optical and args.initial_time_old_units:
        _, scale = build_crossing_params(args)
        if scale is None:
            raise ValueError("--initial-time-old-units requires physical optical mode")
        args.initial_time = args.initial_time * scale["old_dimensionless_scale_ev_inv"]
    run(args)


if __name__ == "__main__":
    main()
