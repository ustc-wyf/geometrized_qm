from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from coordinate_matter_evolution import inverse_metric_block
from mixed_tilde_initial_data import (
    bohm_fields_from_wave_derivatives,
    exact_localized_wave_derivatives,
    safe_covariant_from_inverse,
    tilde_inverse_from_flat_mixed_transform,
)
from physical_units import (
    HBAR_C_EV_M,
    PLANCK_LENGTH_EV_INV,
    PLANCK_MASS_EV,
    ev_inv_to_fs,
    optical_crossing_params_from_physical_scale,
)
from simulate_d_reduced_dynamic_same_initial import relative_l1, sanitize_metric_tail
from simulate_d_tridomain_full_dynamics import (
    current_from_measure_and_u,
    dilate_mask_4,
    erode_mask_8,
    full_tensor_interface_diagnostics,
    kg_positive_frequency_norm,
    local_time_candidate_arrays,
    rk4_patchwise_local_time_matter_step,
    safe_sqrt_abs_det,
    select_local_time_patch_labels,
    support_taper_weights,
    weighted_relative_l1,
)
from simulate_flat_localized_crossing_packets import initial_wavefunction, make_grid, spectral_omega


def real_array(value: np.ndarray) -> np.ndarray:
    return np.asarray(np.real_if_close(value), dtype=float)


def crop_array(field: np.ndarray, ix: np.ndarray, iz: np.ndarray) -> np.ndarray:
    return np.asarray(field)[np.ix_(ix, iz)]


def crop_wave(wave: dict[str, np.ndarray], ix: np.ndarray, iz: np.ndarray) -> dict[str, np.ndarray]:
    return {key: crop_array(value, ix, iz) for key, value in wave.items()}


def scalar_stats(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p05": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p05": float(np.percentile(vals, 5.0)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }


def touches_boundary(mask: np.ndarray) -> bool:
    if mask.size == 0:
        return False
    return bool(np.any(mask[0, :]) or np.any(mask[-1, :]) or np.any(mask[:, 0]) or np.any(mask[:, -1]))


def edge_max_fraction(field: np.ndarray) -> float:
    edge = np.concatenate([field[0, :], field[-1, :], field[:, 0], field[:, -1]])
    return float(np.max(edge) / max(float(np.max(field)), 1.0e-300))


def point_relative_p95(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return 0.0
    rel = np.abs(candidate[valid] - reference[valid]) / np.maximum(np.abs(reference[valid]), 1.0e-300)
    return float(np.percentile(rel, 95.0))


def relative_l1_or_none(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float | None:
    if not np.any(mask):
        return None
    return float(relative_l1(reference, candidate, mask))


def point_relative_stats(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return {"count": 0, "p50": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0}
    rel = np.abs(candidate[valid] - reference[valid]) / np.maximum(np.abs(reference[valid]), 1.0e-300)
    return {
        "count": int(rel.size),
        "p50": float(np.percentile(rel, 50.0)),
        "p95": float(np.percentile(rel, 95.0)),
        "p99": float(np.percentile(rel, 99.0)),
        "max": float(np.max(rel)),
    }


def outlier_mass_fractions(
    reference: np.ndarray,
    candidate: np.ndarray,
    mask: np.ndarray,
    weight: np.ndarray,
    thresholds: tuple[float, ...] = (1.0e-2, 1.0e-1, 1.0),
) -> dict[str, dict[str, float]]:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate) & np.isfinite(weight)
    total_weight = max(float(np.sum(weight[valid])), 1.0e-300)
    result: dict[str, dict[str, float]] = {}
    if not np.any(valid):
        for threshold in thresholds:
            result[f"rel_gt_{threshold:g}"] = {"count": 0, "weight_fraction": 0.0}
        return result
    rel = np.abs(candidate - reference) / np.maximum(np.abs(reference), 1.0e-300)
    for threshold in thresholds:
        bad = valid & (rel > threshold)
        result[f"rel_gt_{threshold:g}"] = {
            "count": int(np.count_nonzero(bad)),
            "weight_fraction": float(np.sum(weight[bad]) / total_weight),
        }
    return result


def format_optional(value: float | None, precision: str = ".3e") -> str:
    if value is None:
        return "n/a"
    return format(value, precision)


def compact_stats(values: np.ndarray | list[float]) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    abs_vals = np.abs(vals)
    return {
        "count": int(vals.size),
        "min": float(np.min(abs_vals)),
        "p50": float(np.percentile(abs_vals, 50.0)),
        "p95": float(np.percentile(abs_vals, 95.0)),
        "max": float(np.max(abs_vals)),
    }


def write_tensor_rows_jsonl(rows: list[dict[str, float]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            serial = {key: float(value) for key, value in row.items() if np.isscalar(value) and np.isfinite(value)}
            handle.write(json.dumps(serial, ensure_ascii=False, sort_keys=True) + "\n")


def tensor_interface_summary(diag: dict[str, object]) -> dict[str, object]:
    rows = list(diag.get("rows", []))
    accepted = [row for row in rows if float(row.get("direct_tensor_solved", 0.0)) > 0.5]
    rejected = [row for row in rows if float(row.get("direct_tensor_solved", 0.0)) <= 0.5]
    residual_all = [float(row.get("direct_tensor_residual_relative", np.nan)) for row in rows]
    residual_accepted = [float(row.get("direct_tensor_residual_relative", np.nan)) for row in accepted]
    initial_residual_all = [float(row.get("initial_direct_tensor_residual_relative", np.nan)) for row in rows]
    return {
        "segment_count": int(diag.get("segment_count", 0)),
        "candidate_count": int(diag.get("candidate_count", 0)),
        "solved_count": int(diag.get("solved_count", 0)),
        "rejected_count": int(len(rejected)),
        "solved_fraction": float(diag.get("solved_fraction", 0.0)),
        "anchor_count": int(diag.get("anchor_count", 0)),
        "anchor_fraction": float(diag.get("anchor_fraction", 0.0)),
        "direct_tensor_residual_relative_all": compact_stats(residual_all),
        "direct_tensor_residual_relative_accepted": compact_stats(residual_accepted),
        "initial_direct_tensor_residual_relative_all": compact_stats(initial_residual_all),
        "definition": (
            "full tensor interface diagnostic: each raw_y=ell^2 Rtilde level segment at raw_y=+/-1 is tested "
            "against the leading thin-layer tensor jump condition using local least-squares covector reconstruction."
        ),
        "accepted_meaning": (
            "Accepted means the local segment admits a real covector satisfying the full tensor residual tolerance. "
            "It is a candidate strong interface anchor, not yet a complete coupled geometry solve."
        ),
    }


def compute_full_tensor_interface_for_state(
    *,
    metric_m: np.ndarray,
    metric_0: np.ndarray,
    metric_p: np.ndarray,
    metric_cov: np.ndarray,
    measure_density: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    dt_probe: float,
    mass: float,
    ell: float,
    mp: float,
    half_width_y: float,
    samples: int,
    tensor_rel_tol: float,
    tensor_max_iter: int,
    tensor_multistart: bool,
    tensor_max_seeds: int,
    tensor_anchor_band_factor: float,
) -> dict[str, object]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, dt_probe, dx, dz)
    metric_inv_geom, _, ricci, r_scalar = geometry_data(metric_cov, dg, d2g)
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure_density / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)
    stress, tilde_x = stress_tensor_tilde(
        metric_cov,
        metric_inv_geom,
        rho_tilde,
        u_t,
        u_x,
        u_z,
        m=mass,
    )
    rhs = real_array(stress) / (float(mp) * float(mp))
    raw_y = float(ell) * float(ell) * real_array(r_scalar)
    diag = full_tensor_interface_diagnostics(
        raw_y=raw_y,
        metric_cov=metric_cov,
        metric_inv=metric_inv_geom,
        ricci=real_array(ricci),
        rho=measure_density,
        rhs_tensor=rhs,
        x=x,
        z=z,
        ell=float(ell),
        half_width_y=float(half_width_y),
        samples=int(samples),
        tensor_rel_tol=float(tensor_rel_tol),
        max_iter=int(tensor_max_iter),
        multistart=bool(tensor_multistart),
        max_seeds=int(tensor_max_seeds),
        band_factor=float(tensor_anchor_band_factor),
    )
    diag["raw_y"] = raw_y
    diag["rho_tilde"] = rho_tilde
    diag["tilde_x"] = real_array(tilde_x)
    diag["r_scalar"] = real_array(r_scalar)
    diag["rhs_tensor_norm"] = np.sqrt(np.einsum("...ab,...ab->...", rhs, rhs, optimize=True))
    return diag


def render_tensor_interface_case(
    *,
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho_a: np.ndarray,
    raw_y: np.ndarray,
    tensor_anchor_mask: np.ndarray,
    rows: list[dict[str, float]],
    summary_text: str,
) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    accepted = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    rejected = [row for row in rows if row.get("direct_tensor_solved", 0.0) <= 0.5]

    def row_segments(selected: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
        return [
            [
                (row["x0"] * HBAR_C_EV_M * 1.0e6, row["z0"] * HBAR_C_EV_M * 1.0e6),
                (row["x1"] * HBAR_C_EV_M * 1.0e6, row["z1"] * HBAR_C_EV_M * 1.0e6),
            ]
            for row in selected
        ]

    raw_vmax = max(float(np.percentile(np.abs(raw_y[np.isfinite(raw_y)]), 99.0)), 1.0)
    rho_vmax = max(float(np.max(rho_a)), 1.0e-300)
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.6), constrained_layout=True)

    im0 = axes[0].pcolormesh(xg, zg, rho_a, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    axes[0].set_title("A rho background")
    fig.colorbar(im0, ax=axes[0])

    im1 = axes[1].pcolormesh(xg, zg, np.clip(raw_y, -raw_vmax, raw_vmax), shading="auto", cmap="coolwarm")
    axes[1].contour(xg, zg, raw_y, levels=[-1.0, 1.0], colors=["black", "white"], linewidths=0.8)
    axes[1].set_title("raw_y=ell^2 Rtilde; contours at -1,+1")
    fig.colorbar(im1, ax=axes[1])

    im2 = axes[2].pcolormesh(xg, zg, tensor_anchor_mask.astype(float), shading="auto", cmap="Greys", vmin=0.0, vmax=1.0)
    axes[2].add_collection(LineCollection(row_segments(rejected), colors="tomato", linewidths=0.75, alpha=0.75))
    axes[2].add_collection(LineCollection(row_segments(accepted), colors="lime", linewidths=1.15, alpha=0.95))
    axes[2].set_title("tensor interface: green accepted, red rejected")
    fig.colorbar(im2, ax=axes[2])

    for ax in axes:
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
    axes[2].text(0.0, -0.14, summary_text, transform=axes[2].transAxes, va="top", ha="left", fontsize=8)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def selected_local_time_admissible_mask(
    *,
    metric_inv: np.ndarray,
    labels: np.ndarray,
    candidates: list[tuple[float, float]],
    active: np.ndarray,
    norm_floor: float,
) -> np.ndarray:
    gtt = metric_inv[..., 0, 0]
    gtx = metric_inv[..., 0, 1]
    gtz = metric_inv[..., 0, 2]
    gxx = metric_inv[..., 1, 1]
    gxz = metric_inv[..., 1, 2]
    gzz = metric_inv[..., 2, 2]
    out = np.zeros_like(active, dtype=bool)
    label_arr = np.asarray(labels, dtype=int)
    for label_value in sorted(int(x) for x in np.unique(label_arr[active])):
        if label_value < 0 or label_value >= len(candidates):
            continue
        a, b = candidates[label_value]
        tau_norm = gtt + 2.0 * a * gtx + 2.0 * b * gtz + (a * a) * gxx + 2.0 * a * b * gxz + (b * b) * gzz
        mask = active & (label_arr == label_value) & np.isfinite(tau_norm) & (tau_norm > norm_floor)
        out |= mask
    return out


def make_cropped_snapshot(
    *,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x_full: np.ndarray,
    z_full: np.ndarray,
    ix: np.ndarray,
    iz: np.ndarray,
    t: float,
    rho_floor: float,
    x_floor: float,
    pinv_rcond: float,
    with_metric: bool,
) -> dict[str, object]:
    wave_full = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x_full, z_full, t)
    rho_full = np.abs(wave_full["psi"]) ** 2
    wave = crop_wave(wave_full, ix, iz)
    bohm = bohm_fields_from_wave_derivatives(wave, rho_floor=rho_floor)

    dx_full = float(x_full[1] - x_full[0])
    dz_full = float(z_full[1] - z_full[0])
    rho_window_fraction = float(
        np.sum(real_array(bohm["rho"])) * dx_full * dz_full / max(float(np.sum(rho_full) * dx_full * dz_full), 1.0e-300)
    )
    result: dict[str, object] = {
        "bohm": bohm,
        "rho_window_fraction": rho_window_fraction,
    }
    if with_metric:
        inv = tilde_inverse_from_flat_mixed_transform(bohm, x_floor=x_floor, pinv_rcond=pinv_rcond)
        cov = safe_covariant_from_inverse(
            inv["tilde_ginv_txz"],
            bohm["rho"],
            det_floor=1.0e-12,
            tail_rel_cut=1.0e-8,
            pinv_rcond=pinv_rcond,
        )
        result.update(
            {
                "inverse": inv,
                "cov_txz": real_array(cov["cov_txz"]),
                "det_txz": real_array(cov["det_txz"]),
                "singular_mask": cov["singular_mask"],
                "tail_mask": cov["tail_mask"],
            }
        )
    return result


def build_physical_reference(args: argparse.Namespace) -> dict[str, object]:
    params, scale = optical_crossing_params_from_physical_scale(
        wavelength_nm=args.wavelength_nm,
        mass_over_omega=args.mass_over_omega,
        resolution=args.full_resolution,
        alpha=args.alpha,
        phi0=args.phi0,
        rho_floor=args.rho_floor,
    )
    x_full, z_full, x_grid, z_grid = make_grid(params)
    dx_full = float(x_full[1] - x_full[0])
    dz_full = float(z_full[1] - z_full[0])
    psi0 = initial_wavefunction(x_grid, z_grid, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x_full, z_full, params.m)
    psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_before = kg_positive_frequency_norm(psi0, psi0_t, dx_full, dz_full)
    if args.normalize_kg:
        if kg_norm_before <= 0.0:
            raise ValueError(f"KG norm must be positive, got {kg_norm_before}")
        psi0 = psi0 / np.sqrt(kg_norm_before)
        psi0_hat = np.fft.fft2(psi0)
        psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_after = kg_positive_frequency_norm(psi0, psi0_t, dx_full, dz_full)
    return {
        "params": params,
        "scale": scale,
        "x_full": x_full,
        "z_full": z_full,
        "psi0": psi0,
        "psi0_hat": psi0_hat,
        "omega": omega,
        "kg_norm_before": float(kg_norm_before),
        "kg_norm_after": float(kg_norm_after),
    }


def render_case(
    *,
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho_a: np.ndarray,
    rho_d_to_a: np.ndarray,
    support: np.ndarray,
    trusted: np.ndarray,
    active: np.ndarray,
    chart_boundary: np.ndarray | None,
    summary_text: str,
) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    diff = rho_d_to_a - rho_a
    valid_support = support & np.isfinite(rho_a) & np.isfinite(rho_d_to_a)
    color_vmax = max(float(np.max(rho_a)), 1.0e-300)
    diff_vmax = max(float(np.percentile(np.abs(diff[valid_support]), 99.0)) if np.any(valid_support) else 0.0, 1.0e-300)
    point_rel = np.abs(diff) / np.maximum(np.abs(rho_a), 1.0e-300)
    rel_vmax = max(float(np.percentile(point_rel[valid_support], 99.0)) if np.any(valid_support) else 0.0, 1.0e-12)

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.0), constrained_layout=True)

    im0 = axes[0, 0].pcolormesh(xg, zg, rho_a, shading="auto", cmap="viridis", vmin=0.0, vmax=color_vmax)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7)
    axes[0, 0].set_title("A rho, white=initial support")
    fig.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].pcolormesh(
        xg,
        zg,
        np.clip(rho_d_to_a, 0.0, color_vmax),
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=color_vmax,
    )
    axes[0, 1].contour(xg, zg, trusted.astype(float), levels=[0.5], colors="cyan", linewidths=0.7)
    if chart_boundary is not None and np.any(chart_boundary):
        axes[0, 1].contour(xg, zg, chart_boundary.astype(float), levels=[0.5], colors="red", linewidths=0.7)
    axes[0, 1].set_title("D pulled back to A rho, same A color range")
    fig.colorbar(im1, ax=axes[0, 1])

    im2 = axes[0, 2].pcolormesh(xg, zg, diff, shading="auto", cmap="coolwarm", vmin=-diff_vmax, vmax=diff_vmax)
    axes[0, 2].contour(xg, zg, support.astype(float), levels=[0.5], colors="black", linewidths=0.6)
    axes[0, 2].set_title("D_to_A - A rho")
    fig.colorbar(im2, ax=axes[0, 2])

    mask_field = np.zeros_like(rho_a)
    mask_field[active] = 1.0
    mask_field[support] = 2.0
    mask_field[trusted] = 3.0
    if chart_boundary is not None:
        mask_field[chart_boundary] = 4.0
    im3 = axes[1, 0].pcolormesh(xg, zg, mask_field, shading="auto", cmap="cividis", vmin=0.0, vmax=4.0)
    axes[1, 0].set_title("mask: 1 active, 2 support, 3 trusted, 4 chart boundary")
    fig.colorbar(im3, ax=axes[1, 0])

    im4 = axes[1, 1].pcolormesh(
        xg,
        zg,
        np.clip(point_rel, 0.0, rel_vmax),
        shading="auto",
        cmap="inferno",
        vmin=0.0,
        vmax=rel_vmax,
    )
    axes[1, 1].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.6)
    if chart_boundary is not None and np.any(chart_boundary):
        axes[1, 1].contour(xg, zg, chart_boundary.astype(float), levels=[0.5], colors="cyan", linewidths=0.6)
    axes[1, 1].set_title("|D_to_A-A|/|A|, clipped at support p99")
    fig.colorbar(im4, ax=axes[1, 1])

    axes[1, 2].axis("off")
    axes[1, 2].text(0.0, 1.0, summary_text, va="top", ha="left", family="monospace", fontsize=9)

    for ax in axes.ravel()[:5]:
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run_case(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    if args.ell_over_planck is not None:
        args.ell = args.ell_over_planck * PLANCK_LENGTH_EV_INV
    elif args.ell is None:
        args.ell = 1.0e60 * PLANCK_LENGTH_EV_INV

    ref = build_physical_reference(args)
    params = ref["params"]
    scale = ref["scale"]
    x_full = ref["x_full"]
    z_full = ref["z_full"]
    psi0 = ref["psi0"]
    psi0_hat = ref["psi0_hat"]
    omega = ref["omega"]

    old_scale = float(scale.old_dimensionless_scale_ev_inv)
    t_meet_old = float(params.t_meet / old_scale)
    old_time = t_meet_old + float(args.tau_old)
    initial_time = old_time * old_scale
    dt = float(args.dt_old) * old_scale
    probe_dt = dt

    window_ev_inv = float(args.window_um) * 1.0e-6 / HBAR_C_EV_M
    ix = np.where(np.abs(x_full) <= window_ev_inv)[0]
    iz = np.where(np.abs(z_full) <= window_ev_inv)[0]
    if ix.size < 8 or iz.size < 8:
        raise ValueError("window contains too few grid points")
    x = x_full[ix]
    z = z_full[iz]
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    snap0 = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=initial_time,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=True,
    )
    snapp = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=initial_time + probe_dt,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=True,
    )
    snapm = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=initial_time - probe_dt,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=True,
    )

    bohm0 = snap0["bohm"]
    rho_a0 = real_array(bohm0["rho"])
    ux = real_array(bohm0["s_x"])
    uz = real_array(bohm0["s_z"])
    ut = real_array(bohm0["s_t"])
    x_field_0 = real_array(bohm0["X"])
    measure0 = np.abs(x_field_0) * rho_a0 / (float(params.m) * float(params.m))
    support = (rho_a0 > args.support_rho_frac * float(np.max(rho_a0))) & (
        measure0 > args.support_measure_frac * float(np.max(measure0))
    )
    trusted = erode_mask_8(support, args.trusted_erosion)
    if not np.any(trusted):
        trusted = support.copy()
    active = dilate_mask_4(support, args.active_dilation)
    stop_mask = trusted if args.stop_mask == "trusted" else support
    if args.evolve_mask == "active":
        evolve_base = active.copy()
    elif args.evolve_mask == "support":
        evolve_base = support.copy()
    elif args.evolve_mask == "trusted":
        evolve_base = trusted.copy()
    else:
        raise ValueError(f"unknown evolve_mask: {args.evolve_mask}")
    if args.boundary_stencil_mask == "active":
        boundary_stencil_base = active.copy()
    elif args.boundary_stencil_mask == "support":
        boundary_stencil_base = support.copy()
    elif args.boundary_stencil_mask == "write":
        boundary_stencil_base = evolve_base.copy()
    else:
        raise ValueError(f"unknown boundary_stencil_mask: {args.boundary_stencil_mask}")
    boundary_stencil_base = boundary_stencil_base | evolve_base
    taper = support_taper_weights(
        rho=rho_a0,
        measure=measure0,
        rho_frac=args.support_rho_frac,
        measure_frac=args.support_measure_frac,
        decades=args.support_taper_decades,
        floor=args.rho_floor,
    )
    taper_active = np.where(active, taper, 0.0)

    metric_m_ref = sanitize_metric_tail(real_array(snapm["cov_txz"]), active)
    metric_0_ref = sanitize_metric_tail(real_array(snap0["cov_txz"]), active)
    metric_p_ref = sanitize_metric_tail(real_array(snapp["cov_txz"]), active)
    metric_cov = metric_0_ref.copy()
    metric_t0 = (metric_p_ref - metric_m_ref) / (2.0 * probe_dt)
    metric_pprev = sanitize_metric_tail(metric_cov - dt * metric_t0, active)
    metric_prev = metric_cov.copy()
    metric_next_linear = sanitize_metric_tail(metric_cov + dt * metric_t0, active)
    del metric_pprev, metric_prev, metric_next_linear

    metric_inv, metric_det = inverse_metric_block(metric_cov)
    measure_density = np.where(active, measure0, 0.0)
    n_cons, current = current_from_measure_and_u(
        measure_density=measure_density,
        u_x=ux,
        u_z=uz,
        metric_cov=metric_cov,
        mass=float(params.m),
        u_t_reference=ut,
        active_mask=active,
        rho_floor=args.rho_floor,
        metric_inv=metric_inv,
        det_cov=metric_det,
    )
    initial_x_g = ut * ut - ux * ux - uz * uz
    initial_rho_pullback = np.where(np.abs(initial_x_g) > args.x_floor, float(params.m) ** 2 * measure0 / np.abs(initial_x_g), 0.0)
    initial_pullback_rel_l1_support = relative_l1(rho_a0, initial_rho_pullback, support)

    a_values, b_values = local_time_candidate_arrays(args.local_time_max_tilt, args.local_time_step)
    cached_labels: np.ndarray | None = None
    cached_candidates: list[tuple[float, float]] | None = None
    cached_age = 0
    last_selector_diag: dict[str, float] = {}
    last_step_diag: dict[str, float] = {}
    stopped_reason = "wall_time_budget_reached"
    step = 0
    advanced_time = 0.0
    min_accepted_dt = float("inf")
    max_halvings_used = 0
    rejected_attempts = 0
    chart_boundary_max_count = 0
    chart_boundary_max_measure_fraction = 0.0
    chart_boundary_max_trusted_measure_fraction = 0.0
    last_chart_boundary_count = 0
    last_step_active_count = int(np.count_nonzero(evolve_base))
    last_chart_boundary_mask = np.zeros_like(active, dtype=bool)
    last_step_active_mask = evolve_base.copy()
    start = time.perf_counter()
    while step < args.max_steps and (time.perf_counter() - start) < args.wall_seconds:
        metric_inv, metric_det = inverse_metric_block(metric_cov)
        n_cons, current = current_from_measure_and_u(
            measure_density=measure_density,
            u_x=ux,
            u_z=uz,
            metric_cov=metric_cov,
            mass=float(params.m),
            u_t_reference=current.get("u_t", ut),
            active_mask=active,
            rho_floor=args.rho_floor,
            metric_inv=metric_inv,
            det_cov=metric_det,
        )
        disc_min = float(np.nanmin(current["discriminant"][stop_mask])) if np.any(stop_mask) else 0.0
        if args.stop_on_negative_discriminant and disc_min < -abs(args.disc_tolerance):
            stopped_reason = f"negative mass-shell discriminant on {args.stop_mask}: {disc_min:.6e}"
            break

        reuse_atlas = (
            args.local_time_atlas_every > 1
            and cached_labels is not None
            and cached_candidates is not None
            and cached_age < args.local_time_atlas_every
        )
        if reuse_atlas:
            labels = cached_labels
            candidates = cached_candidates
            cached_age += 1
            selector_diag = dict(last_selector_diag)
            selector_diag["atlas_reused"] = 1.0
        else:
            labels, candidates, selector_diag = select_local_time_patch_labels(
                current=current,
                metric_inv=metric_inv,
                measure=measure_density,
                active=evolve_base,
                trusted=stop_mask,
                a_values=a_values,
                b_values=b_values,
                norm_floor=args.local_time_norm_floor,
                include_stationary_candidates=args.local_time_stationary_candidates,
                stationary_candidate_quantization=args.local_time_stationary_quantization,
                tile_size=args.local_time_tile_size,
            )
            selector_diag["atlas_reused"] = 0.0
            cached_labels = labels
            cached_candidates = candidates
            cached_age = 1

        step_active = evolve_base
        chart_boundary = np.zeros_like(active, dtype=bool)
        if args.freeze_uncovered_chart:
            chart_admissible = selected_local_time_admissible_mask(
                metric_inv=metric_inv,
                labels=labels,
                candidates=candidates,
                active=evolve_base,
                norm_floor=args.local_time_norm_floor,
            )
            chart_boundary = evolve_base & (~chart_admissible)
            if args.chart_boundary_halo > 0:
                chart_boundary = dilate_mask_4(chart_boundary, args.chart_boundary_halo) & evolve_base
            step_active = evolve_base & chart_admissible
            if args.chart_boundary_halo > 0:
                step_active = evolve_base & (~chart_boundary)
            if not np.any(step_active):
                stopped_reason = "no admissible local-time chart cells remain"
                break
        last_chart_boundary_count = int(np.count_nonzero(chart_boundary))
        last_step_active_count = int(np.count_nonzero(step_active))
        last_chart_boundary_mask = chart_boundary.copy()
        last_step_active_mask = step_active.copy()
        total_measure_active = max(float(np.sum(measure_density[active])), 1.0e-300)
        total_measure_trusted = max(float(np.sum(measure_density[trusted])), 1.0e-300)
        chart_measure_fraction = float(np.sum(measure_density[chart_boundary]) / total_measure_active)
        chart_trusted_measure_fraction = float(np.sum(measure_density[chart_boundary & trusted]) / total_measure_trusted)
        chart_boundary_max_count = max(chart_boundary_max_count, last_chart_boundary_count)
        chart_boundary_max_measure_fraction = max(chart_boundary_max_measure_fraction, chart_measure_fraction)
        chart_boundary_max_trusted_measure_fraction = max(chart_boundary_max_trusted_measure_fraction, chart_trusted_measure_fraction)

        base_measure = measure_density
        base_ux = ux
        base_uz = uz
        base_current = current
        accepted = False
        best_reject_reason = "not attempted"
        max_halvings = int(args.max_step_halvings) if args.adaptive_step else 0
        for halving in range(max_halvings + 1):
            dt_try = dt / (2.0**halving)
            try:
                trial_measure, trial_n, trial_ux, trial_uz, trial_current, step_diag = rk4_patchwise_local_time_matter_step(
                    measure_density=base_measure,
                    u_x=base_ux,
                    u_z=base_uz,
                    metric_cov=metric_cov,
                    metric_inv=metric_inv,
                    metric_det=metric_det,
                    mass=float(params.m),
                    dx=dx,
                    dz=dz,
                    dt=dt_try,
                    u_t_reference=base_current["u_t"],
                    active_mask=step_active,
                    stencil_mask=boundary_stencil_base,
                    labels=labels,
                    candidates=candidates,
                    flux_divergence_mode=args.boundary_flux_mode,
                    rusanov_strength=args.rusanov_strength,
                    matching_flux_weight=args.matching_flux_weight,
                    use_bounding_boxes=args.local_time_patch_bboxes,
                    split_components=args.local_time_split_components,
                )
            except FloatingPointError as exc:
                best_reject_reason = f"non-finite trial: {exc}"
                rejected_attempts += 1
                continue

            trial_measure = np.where(step_active, trial_measure, base_measure)
            trial_ux = np.where(step_active, trial_ux, base_ux)
            trial_uz = np.where(step_active, trial_uz, base_uz)
            trial_n, trial_current = current_from_measure_and_u(
                measure_density=trial_measure,
                u_x=trial_ux,
                u_z=trial_uz,
                metric_cov=metric_cov,
                mass=float(params.m),
                u_t_reference=base_current["u_t"],
                active_mask=active,
                rho_floor=args.rho_floor,
                metric_inv=metric_inv,
                det_cov=metric_det,
            )
            trial_disc_min = float(np.nanmin(trial_current["discriminant"][stop_mask])) if np.any(stop_mask) else 0.0
            trial_measure_delta = float(step_diag.get("measure_rel_delta", 0.0))
            finite_ok = (
                np.all(np.isfinite(trial_measure[active]))
                and np.all(np.isfinite(trial_ux[active]))
                and np.all(np.isfinite(trial_uz[active]))
                and np.all(np.isfinite(trial_current["u_t"][active]))
            )
            disc_ok = (not args.stop_on_negative_discriminant) or trial_disc_min >= -abs(args.disc_tolerance)
            measure_ok = trial_measure_delta <= float(args.max_measure_rel_delta)
            if finite_ok and disc_ok and measure_ok:
                measure_density, n_cons, ux, uz, current = trial_measure, trial_n, trial_ux, trial_uz, trial_current
                step_diag = dict(step_diag)
                step_diag["accepted_dt"] = float(dt_try)
                step_diag["accepted_dt_old"] = float(dt_try / old_scale)
                step_diag["adaptive_halvings"] = float(halving)
                step_diag["trial_disc_min"] = trial_disc_min
                step_diag["chart_boundary_count"] = float(last_chart_boundary_count)
                step_diag["step_active_count"] = float(last_step_active_count)
                step_diag["chart_boundary_measure_fraction"] = chart_measure_fraction
                step_diag["chart_boundary_trusted_measure_fraction"] = chart_trusted_measure_fraction
                accepted = True
                max_halvings_used = max(max_halvings_used, halving)
                min_accepted_dt = min(min_accepted_dt, dt_try)
                break
            rejected_attempts += 1
            best_reject_reason = (
                f"disc_min={trial_disc_min:.6e}, measure_rel_delta={trial_measure_delta:.6e}, finite_ok={finite_ok}"
            )
        if not accepted:
            stopped_reason = f"adaptive step failed after {max_halvings + 1} attempts: {best_reject_reason}"
            break
        last_selector_diag = selector_diag
        last_step_diag = step_diag
        advanced_time += float(step_diag.get("accepted_dt", dt))
        step += 1
    else:
        if step >= args.max_steps:
            stopped_reason = "max_steps_reached"

    evolution_wall_seconds = time.perf_counter() - start
    metric_inv, metric_det = inverse_metric_block(metric_cov)
    n_cons, current = current_from_measure_and_u(
        measure_density=measure_density,
        u_x=ux,
        u_z=uz,
        metric_cov=metric_cov,
        mass=float(params.m),
        u_t_reference=current.get("u_t", ut),
        active_mask=active,
        rho_floor=args.rho_floor,
        metric_inv=metric_inv,
        det_cov=metric_det,
    )

    final_time = initial_time + advanced_time
    final_snap = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=final_time,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=False,
    )
    bohm_final = final_snap["bohm"]
    rho_a_final = real_array(bohm_final["rho"])
    x_a_final = real_array(bohm_final["X"])
    measure_a_final = np.abs(x_a_final) * rho_a_final / (float(params.m) * float(params.m))

    x_g_d = current["u_t"] * current["u_t"] - ux * ux - uz * uz
    rho_d_to_a = np.where(
        active & (np.abs(x_g_d) > args.x_floor),
        float(params.m) ** 2 * measure_density / np.maximum(np.abs(x_g_d), args.x_floor),
        0.0,
    )
    measure_d = np.where(active, measure_density, 0.0)
    rho_core = rho_a_final > args.core_rho_frac * float(np.max(rho_a_final))
    if not np.any(rho_core):
        rho_core = trusted.copy()
    chart_boundary_halo_mask = dilate_mask_4(last_chart_boundary_mask, 1) & active
    evolved_active = active & (~last_chart_boundary_mask)
    evolved_support = support & evolved_active
    evolved_trusted = trusted & evolved_active
    evolved_trusted_no_halo = trusted & (~chart_boundary_halo_mask)
    evolved_core = rho_core & evolved_active
    chart_boundary_support = support & last_chart_boundary_mask
    chart_boundary_trusted = trusted & last_chart_boundary_mask
    support_rel_l1 = relative_l1(rho_a_final, rho_d_to_a, support)
    trusted_rel_l1 = relative_l1(rho_a_final, rho_d_to_a, trusted)
    core_rel_l1 = relative_l1(rho_a_final, rho_d_to_a, rho_core)
    evolved_trusted_rel_l1 = relative_l1_or_none(rho_a_final, rho_d_to_a, evolved_trusted)
    evolved_trusted_no_halo_rel_l1 = relative_l1_or_none(rho_a_final, rho_d_to_a, evolved_trusted_no_halo)

    tensor_anchor_mask = np.zeros_like(active, dtype=bool)
    tensor_summary: dict[str, object] = {
        "enabled": False,
        "segment_count": 0,
        "candidate_count": 0,
        "solved_count": 0,
        "solved_fraction": 0.0,
        "anchor_count": 0,
        "anchor_fraction": 0.0,
    }
    tensor_rows_path: Path | None = None
    tensor_figure_path: Path | None = None
    tensor_raw_y: np.ndarray | None = None
    if args.tensor_interface_diagnostics:
        tensor_start = time.perf_counter()
        tensor_diag = compute_full_tensor_interface_for_state(
            metric_m=metric_m_ref,
            metric_0=metric_0_ref,
            metric_p=metric_p_ref,
            metric_cov=metric_cov,
            measure_density=measure_density,
            u_t=current["u_t"],
            u_x=ux,
            u_z=uz,
            x=x,
            z=z,
            dt_probe=probe_dt,
            mass=float(params.m),
            ell=float(args.ell),
            mp=float(args.mp),
            half_width_y=float(args.half_width_y),
            samples=int(args.samples),
            tensor_rel_tol=float(args.tensor_rel_tol),
            tensor_max_iter=int(args.tensor_max_iter),
            tensor_multistart=bool(args.tensor_multistart),
            tensor_max_seeds=int(args.tensor_max_seeds),
            tensor_anchor_band_factor=float(args.tensor_anchor_band_factor),
        )
        tensor_anchor_mask = np.asarray(tensor_diag["tensor_anchor_mask"], dtype=bool)
        tensor_raw_y = np.asarray(tensor_diag["raw_y"], dtype=float)
        rows = list(tensor_diag.get("rows", []))
        tensor_rows_path = out / "tensor_interface_rows.jsonl"
        write_tensor_rows_jsonl(rows, tensor_rows_path)
        tensor_summary = tensor_interface_summary(tensor_diag)
        tensor_summary.update(
            {
                "enabled": True,
                "wall_seconds": float(time.perf_counter() - tensor_start),
                "raw_y_stats_active": scalar_stats(tensor_raw_y[active]),
                "raw_y_stats_support": scalar_stats(tensor_raw_y[support]),
                "rho_tilde_stats_support": scalar_stats(np.asarray(tensor_diag["rho_tilde"], dtype=float)[support]),
                "tilde_x_minus_m2_stats_support": scalar_stats(
                    (np.asarray(tensor_diag["tilde_x"], dtype=float) - float(params.m) ** 2)[support]
                ),
                "rhs_tensor_norm_stats_support": scalar_stats(
                    np.asarray(tensor_diag["rhs_tensor_norm"], dtype=float)[support]
                ),
                "geometry_history_scope": (
                    "The current local-window metric is fixed to the initial A-derived metric. "
                    "The tensor diagnostic therefore uses the initial A-derived metric jets "
                    "(t0-dt,t0,t0+dt) and the evolved D matter state. It is a strong interface "
                    "diagnostic/anchor candidate, not yet a full coupled metric evolution."
                ),
            }
        )
        tensor_figure_path = out / "tensor_interface_diagnostic.png"
        render_tensor_interface_case(
            out_path=tensor_figure_path,
            x=x,
            z=z,
            rho_a=rho_a_final,
            raw_y=tensor_raw_y,
            tensor_anchor_mask=tensor_anchor_mask,
            rows=rows,
            summary_text=(
                f"segments={tensor_summary['segment_count']}, "
                f"accepted={tensor_summary['solved_count']}, "
                f"fraction={tensor_summary['solved_fraction']:.3f}\n"
                f"accepted residual p95="
                f"{tensor_summary['direct_tensor_residual_relative_accepted']['p95']:.3e}"
            ),
        )

    figure_path = out / "rho_A_vs_D_pulled_back.png"
    summary_text = "\n".join(
        [
            f"tau_old={args.tau_old:g}, old_t0={old_time:.6g}",
            f"steps={step}, wall={evolution_wall_seconds:.1f}s",
            f"advanced_old={advanced_time / old_scale:.6g}",
            f"rho relL1 support={support_rel_l1:.3e}",
            f"rho relL1 trusted={trusted_rel_l1:.3e}",
            f"rho relL1 evolved_trusted={format_optional(evolved_trusted_rel_l1)}",
            f"rho relL1 no_boundary_halo={format_optional(evolved_trusted_no_halo_rel_l1)}",
            f"rho relL1 core={core_rel_l1:.3e}",
            f"chart boundary cells={last_chart_boundary_count}",
            f"disc_min trusted={float(np.nanmin(current['discriminant'][trusted])):.3e}",
            f"active touches boundary={touches_boundary(active)}",
        ]
    )
    render_case(
        out_path=figure_path,
        x=x,
        z=z,
        rho_a=rho_a_final,
        rho_d_to_a=rho_d_to_a,
        support=support,
        trusted=trusted,
        active=active,
        chart_boundary=last_chart_boundary_mask,
        summary_text=summary_text,
    )

    fields_path = out / "fields_final.npz"
    np.savez_compressed(
        fields_path,
        x=x,
        z=z,
        support=support,
        trusted=trusted,
        active=active,
        boundary_stencil_base=boundary_stencil_base,
        chart_boundary=last_chart_boundary_mask,
        step_active=last_step_active_mask,
        chart_boundary_halo=chart_boundary_halo_mask,
        evolve_base=evolve_base,
        evolved_active=evolved_active,
        evolved_support=evolved_support,
        evolved_trusted=evolved_trusted,
        evolved_trusted_no_halo=evolved_trusted_no_halo,
        rho_A=rho_a_final,
        rho_D_to_A=rho_d_to_a,
        ntilde_measure_A=measure_a_final,
        ntilde_measure_D=measure_d,
        u_t_D=current["u_t"],
        u_x_D=ux,
        u_z_D=uz,
        X_g_D=x_g_d,
        metric_cov_D=metric_cov,
        taper=taper,
        tensor_anchor_mask=tensor_anchor_mask,
        tensor_raw_y=np.zeros_like(rho_a_final) if tensor_raw_y is None else tensor_raw_y,
    )

    dx_old = dx / old_scale
    fringe_period_old = 2.0 * np.pi / (np.sqrt(2.0) * 6.0)
    summary = {
        "params": {
            "scope": "D local-window diagnostic initialized from a high-resolution full-domain A/KG snapshot cropped to the physical observation window. It does not reinitialize a small periodic wave packet.",
            "units": "hbar=c=1, eV/eV^-1 internally; plots use micrometers.",
            "wavelength_nm": args.wavelength_nm,
            "mass_over_omega": args.mass_over_omega,
            "normalize_kg": bool(args.normalize_kg),
            "kg_norm_before": ref["kg_norm_before"],
            "kg_norm_after": ref["kg_norm_after"],
            "full_resolution": args.full_resolution,
            "window_um": args.window_um,
            "window_grid_shape": [int(ix.size), int(iz.size)],
            "dx_um": float(dx * HBAR_C_EV_M * 1.0e6),
            "dx_old": float(dx_old),
            "points_per_interference_fringe": float(fringe_period_old / dx_old),
            "old_dimensionless_scale_ev_inv": old_scale,
            "t_meet_old": t_meet_old,
            "tau_old": float(args.tau_old),
            "initial_time_old": old_time,
            "initial_time_ev_inv": float(initial_time),
            "initial_time_fs": ev_inv_to_fs(initial_time),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(dt),
            "dt_fs": ev_inv_to_fs(dt),
            "adaptive_step": bool(args.adaptive_step),
            "max_step_halvings": int(args.max_step_halvings),
            "max_measure_rel_delta": float(args.max_measure_rel_delta),
            "min_accepted_dt_old": None if not np.isfinite(min_accepted_dt) else float(min_accepted_dt / old_scale),
            "max_halvings_used": int(max_halvings_used),
            "adaptive_rejected_attempts": int(rejected_attempts),
            "freeze_uncovered_chart": bool(args.freeze_uncovered_chart),
            "chart_boundary_halo": int(args.chart_boundary_halo),
            "chart_boundary_max_count": int(chart_boundary_max_count),
            "chart_boundary_max_measure_fraction": float(chart_boundary_max_measure_fraction),
            "chart_boundary_max_trusted_measure_fraction": float(chart_boundary_max_trusted_measure_fraction),
            "last_chart_boundary_count": int(last_chart_boundary_count),
            "last_step_active_count": int(last_step_active_count),
            "wall_seconds_requested": float(args.wall_seconds),
            "max_steps": int(args.max_steps),
            "steps_completed": int(step),
            "old_time_advanced": float(advanced_time / old_scale),
            "t_completed_ev_inv": float(advanced_time),
            "t_completed_fs": ev_inv_to_fs(advanced_time),
            "final_time_old": float(old_time + advanced_time / old_scale),
            "evolution_wall_seconds": float(evolution_wall_seconds),
            "steps_per_second": float(step / max(evolution_wall_seconds, 1.0e-300)),
            "stopped_reason": stopped_reason,
            "support_rho_frac": args.support_rho_frac,
            "support_measure_frac": args.support_measure_frac,
            "trusted_erosion": args.trusted_erosion,
            "active_dilation": args.active_dilation,
            "stop_mask": args.stop_mask,
            "evolve_mask": args.evolve_mask,
            "boundary_stencil_mask": args.boundary_stencil_mask,
            "boundary_flux_mode": args.boundary_flux_mode,
            "rusanov_strength": float(args.rusanov_strength),
            "matching_flux_weight": args.matching_flux_weight,
            "ell": float(args.ell),
            "ell_over_planck_length": float(args.ell / PLANCK_LENGTH_EV_INV),
            "mp": float(args.mp),
            "boundary_stencil_count": int(np.count_nonzero(boundary_stencil_base)),
            "fixed_buffer_count": int(np.count_nonzero(boundary_stencil_base & (~evolve_base))),
            "local_time_max_tilt": args.local_time_max_tilt,
            "local_time_step": args.local_time_step,
            "local_time_stationary_candidates": bool(args.local_time_stationary_candidates),
            "local_time_stationary_quantization": args.local_time_stationary_quantization,
            "local_time_tile_size": args.local_time_tile_size,
            "local_time_atlas_every": args.local_time_atlas_every,
            "local_time_patch_bboxes": bool(args.local_time_patch_bboxes),
            "local_time_split_components": bool(args.local_time_split_components),
            "tensor_interface_diagnostics": bool(args.tensor_interface_diagnostics),
            "half_width_y": float(args.half_width_y),
            "samples": int(args.samples),
            "tensor_rel_tol": float(args.tensor_rel_tol),
            "tensor_max_iter": int(args.tensor_max_iter),
            "tensor_max_seeds": int(args.tensor_max_seeds),
            "tensor_multistart": bool(args.tensor_multistart),
            "tensor_anchor_band_factor": float(args.tensor_anchor_band_factor),
        },
        "definitions": {
            "rho_D_to_A": "D branch transformed matter measure pulled back to the flat/A representation: rho_D_to_A = m^2 * (sqrt(|gtilde|) rho_tilde)_D / |X_g[D]|.",
            "X_g_D": "Flat g-frame mass-shell scalar of the evolved D covector, X_g[D]=u_t^2-u_x^2-u_z^2.",
            "support": "Initial binary support: rho_A0 > support_rho_frac*max(rho_A0) and ntilde_A0 > support_measure_frac*max(ntilde_A0).",
            "trusted": "support eroded by trusted_erosion grid cells; used to avoid support-edge stencil artifacts.",
            "active": "support dilated by active_dilation grid cells; outside active the D matter variables are not evolved.",
            "evolve_mask": "Which region is explicitly advanced by the local-time matter solver. Cells inside active but outside evolve_mask are boundary/buffer cells for this reduced local-window solve.",
            "boundary_stencil_mask": "Fixed buffer/ghost region used by the local-time finite-difference stencil. Only evolve_mask cells are written back; buffer cells outside evolve_mask provide boundary flux data and are not bulk-evolved.",
            "boundary_flux_mode": "central uses the legacy centered derivative; face uses explicit non-periodic centered finite-volume face fluxes; rusanov adds a local Lax-Friedrichs/Rusanov penalty; matching replaces trusted-buffer face fluxes by a single normal conserved-current flux from the one-sided matching condition. This matching mode is not yet full tensor matching.",
            "rusanov_strength": "Penalty multiplier for boundary_flux_mode=rusanov. 0 is centered face flux, 1 is standard Rusanov.",
            "matching_flux_weight": "For boundary_flux_mode=matching, equal uses an unweighted least-squares normal-flux continuity condition, while abs_n weights the two one-sided fluxes by |n_tau| so low-density halo cells do not dominate the matched face flux.",
            "full_tensor_interface": "Direct test of the leading thin-layer tensor jump condition on raw_y=ell^2 Rtilde = +/-1 segments. It is stricter than trace-only and separate from the material conserved-current matching flux.",
            "tensor_anchor_mask": "Grid cells lying within tensor_anchor_band_factor*sqrt(dx^2+dz^2) of accepted full-tensor interface segments. This is saved for later metric/interface closure; it is not yet used by this local-window matter-only evolution.",
            "chart_boundary": "Cells inside evolve_mask that are not covered by the currently selected admissible local-time chart, optionally dilated by chart_boundary_halo. They are treated as temporary chart/interface boundary cells for that step, not forced through a bad explicit update.",
            "evolved_trusted": "trusted cells excluding the final chart_boundary mask.",
            "evolved_trusted_no_halo": "trusted cells excluding a one-cell halo around the final chart_boundary mask; used to separate bulk evolution from interface-neighbor stencil artifacts.",
            "white_contour": "White contour in the figure is support on A rho panel and support on point-relative panel.",
            "cyan_contour": "Cyan contour in the figure is trusted on D pulled-back rho panel, and chart-boundary on the point-relative panel when present.",
            "red_contour": "Red contour on the D pulled-back rho panel is the final chart_boundary/interface mask.",
        },
        "initial": {
            "rho_window_fraction": float(snap0["rho_window_fraction"]),
            "rho_edge_max_fraction": edge_max_fraction(rho_a0),
            "initial_pullback_rel_l1_support": float(initial_pullback_rel_l1_support),
            "support_count": int(np.count_nonzero(support)),
            "trusted_count": int(np.count_nonzero(trusted)),
            "active_count": int(np.count_nonzero(active)),
            "evolve_base_count": int(np.count_nonzero(evolve_base)),
            "support_touches_boundary": touches_boundary(support),
            "trusted_touches_boundary": touches_boundary(trusted),
            "active_touches_boundary": touches_boundary(active),
            "singular_count_active": int(np.count_nonzero(np.asarray(snap0["singular_mask"], dtype=bool) & active)),
            "tail_count_active": int(np.count_nonzero(np.asarray(snap0["tail_mask"], dtype=bool) & active)),
        },
        "final": {
            "rho_window_fraction": float(final_snap["rho_window_fraction"]),
            "rho_edge_max_fraction": edge_max_fraction(rho_a_final),
            "rho_pullback_rel_l1_support": relative_l1(rho_a_final, rho_d_to_a, support),
            "rho_pullback_rel_l1_trusted": relative_l1(rho_a_final, rho_d_to_a, trusted),
            "rho_pullback_rel_l1_core": relative_l1(rho_a_final, rho_d_to_a, rho_core),
            "rho_pullback_rel_l1_evolved_support": relative_l1_or_none(rho_a_final, rho_d_to_a, evolved_support),
            "rho_pullback_rel_l1_evolved_trusted": evolved_trusted_rel_l1,
            "rho_pullback_rel_l1_evolved_trusted_no_chart_halo": evolved_trusted_no_halo_rel_l1,
            "rho_pullback_rel_l1_evolved_core": relative_l1_or_none(rho_a_final, rho_d_to_a, evolved_core),
            "rho_pullback_rel_l1_chart_boundary_support": relative_l1_or_none(rho_a_final, rho_d_to_a, chart_boundary_support),
            "rho_pullback_rel_l1_chart_boundary_trusted": relative_l1_or_none(rho_a_final, rho_d_to_a, chart_boundary_trusted),
            "rho_pullback_rel_l1_tapered": weighted_relative_l1(rho_a_final, rho_d_to_a, taper_active),
            "rho_pullback_point_rel_p95_support": point_relative_p95(rho_a_final, rho_d_to_a, support),
            "rho_pullback_point_rel_p95_trusted": point_relative_p95(rho_a_final, rho_d_to_a, trusted),
            "rho_pullback_point_rel_p95_core": point_relative_p95(rho_a_final, rho_d_to_a, rho_core),
            "rho_pullback_point_rel_stats_support": point_relative_stats(rho_a_final, rho_d_to_a, support),
            "rho_pullback_point_rel_stats_trusted": point_relative_stats(rho_a_final, rho_d_to_a, trusted),
            "rho_pullback_point_rel_stats_evolved_trusted": point_relative_stats(rho_a_final, rho_d_to_a, evolved_trusted),
            "rho_pullback_outlier_mass_fractions_support": outlier_mass_fractions(
                rho_a_final, rho_d_to_a, support, rho_a_final
            ),
            "rho_pullback_outlier_mass_fractions_trusted": outlier_mass_fractions(
                rho_a_final, rho_d_to_a, trusted, rho_a_final
            ),
            "rho_pullback_outlier_mass_fractions_evolved_trusted": outlier_mass_fractions(
                rho_a_final, rho_d_to_a, evolved_trusted, rho_a_final
            ),
            "chart_boundary_support_count": int(np.count_nonzero(chart_boundary_support)),
            "chart_boundary_trusted_count": int(np.count_nonzero(chart_boundary_trusted)),
            "chart_boundary_support_rho_mass_fraction": float(
                np.sum(rho_a_final[chart_boundary_support]) / max(float(np.sum(rho_a_final[support])), 1.0e-300)
            ),
            "chart_boundary_trusted_rho_mass_fraction": float(
                np.sum(rho_a_final[chart_boundary_trusted]) / max(float(np.sum(rho_a_final[trusted])), 1.0e-300)
            ),
            "measure_rel_l1_support": relative_l1(measure_a_final, measure_d, support),
            "measure_rel_l1_trusted": relative_l1(measure_a_final, measure_d, trusted),
            "measure_rel_l1_evolved_trusted": relative_l1_or_none(measure_a_final, measure_d, evolved_trusted),
            "measure_rel_l1_evolved_trusted_no_chart_halo": relative_l1_or_none(
                measure_a_final, measure_d, evolved_trusted_no_halo
            ),
            "measure_rel_l1_tapered": weighted_relative_l1(measure_a_final, measure_d, taper_active),
            "disc_min_support": float(np.nanmin(current["discriminant"][support])) if np.any(support) else 0.0,
            "disc_min_trusted": float(np.nanmin(current["discriminant"][trusted])) if np.any(trusted) else 0.0,
            "rho_D_to_A_stats_support": scalar_stats(rho_d_to_a[support]),
            "rho_A_stats_support": scalar_stats(rho_a_final[support]),
            "X_g_D_stats_support": scalar_stats(x_g_d[support]),
            "full_tensor_interface": tensor_summary,
            "local_time": {
                **{key: float(value) for key, value in last_selector_diag.items() if np.isscalar(value)},
                **{f"step_{key}": float(value) for key, value in last_step_diag.items() if np.isscalar(value)},
            },
        },
        "artifacts": {
            "figure": str(figure_path.resolve()),
            "fields_npz": str(fields_path.resolve()),
            "summary": str((out / "summary.json").resolve()),
            "tensor_interface_figure": None if tensor_figure_path is None else str(tensor_figure_path.resolve()),
            "tensor_interface_rows_jsonl": None if tensor_rows_path is None else str(tensor_rows_path.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau-old", type=float, required=True, help="Shifted old time tau = old_t - t_meet_old.")
    parser.add_argument("--wall-seconds", type=float, default=300.0)
    parser.add_argument("--max-steps", type=int, default=10_000_000)
    parser.add_argument("--full-resolution", type=int, default=640)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-5)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--ell", type=float, default=None, help="D-branch length scale in eV^-1; default is ell_over_planck/M_P.")
    parser.add_argument(
        "--ell-over-planck",
        type=float,
        default=1.0e60,
        help="Set ell = this factor times the Planck length 1/M_P. Use with the physical Planck mass default.",
    )
    parser.add_argument("--mp", type=float, default=None, help="Planck mass parameter in eV; default is physical M_P.")
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-taper-decades", type=float, default=1.0)
    parser.add_argument("--core-rho-frac", type=float, default=0.1)
    parser.add_argument("--active-dilation", type=int, default=2)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument("--stop-mask", choices=["support", "trusted"], default="trusted")
    parser.add_argument(
        "--evolve-mask",
        choices=["active", "support", "trusted"],
        default="active",
        help="Region explicitly advanced by the local-time solver. Use trusted to treat the low-density support edge as a fixed boundary/buffer.",
    )
    parser.add_argument(
        "--boundary-stencil-mask",
        choices=["write", "support", "active"],
        default="active",
        help=(
            "Region used as fixed buffer/ghost data for local-time finite-difference stencils. "
            "Only evolve_mask cells are written back; support/active cells outside evolve_mask provide boundary flux data."
        ),
    )
    parser.add_argument(
        "--boundary-flux-mode",
        choices=["central", "face", "rusanov", "matching"],
        default="central",
        help=(
            "How the local-time conserved-density flux divergence is computed. "
            "central keeps the legacy centered derivative; face uses explicit non-periodic finite-volume centered face fluxes; "
            "rusanov adds a local Lax-Friedrichs/Rusanov penalty to the face flux; matching uses a single normal conserved-current flux on trusted-buffer faces."
        ),
    )
    parser.add_argument(
        "--rusanov-strength",
        type=float,
        default=1.0,
        help="Penalty strength for --boundary-flux-mode rusanov. 0 equals centered face flux; 1 is standard Rusanov.",
    )
    parser.add_argument(
        "--matching-flux-weight",
        choices=["equal", "abs_n"],
        default="abs_n",
        help=(
            "Weight used by --boundary-flux-mode matching. equal enforces unweighted one-sided normal-flux matching; "
            "abs_n weights the two one-sided fluxes by |n_tau|."
        ),
    )
    parser.add_argument("--disc-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--stop-on-negative-discriminant", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--adaptive-step",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Halve a local-time RK step until algebraic admissibility checks pass. This changes step size only, not the equations.",
    )
    parser.add_argument("--max-step-halvings", type=int, default=12)
    parser.add_argument(
        "--max-measure-rel-delta",
        type=float,
        default=1.0e-2,
        help="Reject a single explicit step if sqrt(|gtilde|)rho_tilde changes by more than this relative amount on active cells.",
    )
    parser.add_argument(
        "--freeze-uncovered-chart",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Treat cells not covered by the selected local-time chart as temporary chart-boundary cells for that step instead of forcing a bad lab-time update.",
    )
    parser.add_argument(
        "--chart-boundary-halo",
        type=int,
        default=0,
        help="Dilate temporary chart-boundary cells by this many 4-neighbor layers before freezing them. This is a chart/interface treatment, not damping.",
    )
    parser.add_argument("--local-time-max-tilt", type=float, default=5.0)
    parser.add_argument("--local-time-step", type=float, default=0.5)
    parser.add_argument("--local-time-norm-floor", type=float, default=1.0e-12)
    parser.add_argument("--local-time-stationary-candidates", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--local-time-stationary-quantization", type=float, default=0.1)
    parser.add_argument("--local-time-tile-size", type=int, default=1)
    parser.add_argument("--local-time-atlas-every", type=int, default=1)
    parser.add_argument("--local-time-patch-bboxes", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--local-time-split-components", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument(
        "--tensor-interface-diagnostics",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Compute direct full-tensor interface jump diagnostics on raw_y=ell^2 Rtilde=+/-1 after the local-window run. "
            "This does not yet enforce the tensor condition as an evolution boundary."
        ),
    )
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=61)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument("--tensor-max-iter", type=int, default=30)
    parser.add_argument("--tensor-max-seeds", type=int, default=24)
    parser.add_argument("--tensor-multistart", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--tensor-anchor-band-factor", type=float, default=1.0)
    args = parser.parse_args()
    run_case(args)


if __name__ == "__main__":
    main()
