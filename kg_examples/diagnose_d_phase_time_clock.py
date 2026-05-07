from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_d_jt_degeneracy import parse_float_list, real_array, scalar_stats, sign_change_edge_count
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
from phase_clock_projection import project_to_gradient, project_to_weighted_gradient
from physical_units import HBAR_C_EV_M, ev_inv_to_fs
from simulate_d_tridomain_full_dynamics import (
    build_crossing_params,
    erode_mask_8,
    kg_positive_frequency_norm,
    safe_sqrt_abs_det,
    tilde_density_from_measure,
    tilde_measure_from_bohm,
)
from simulate_flat_localized_crossing_packets import initial_wavefunction, make_grid, spectral_omega


def measure_fraction(mask: np.ndarray, measure: np.ndarray, base: np.ndarray) -> float:
    denom = max(float(np.sum(measure[base])), 1.0e-300)
    return float(np.sum(measure[base & mask]) / denom)


def summarize_mask(
    *,
    mask_name: str,
    mask: np.ndarray,
    measure: np.ndarray,
    jt_unit: np.ndarray,
    js_unit: np.ndarray,
    mass_sq: float,
    thresholds: list[float],
) -> dict[str, object]:
    rel_js = (js_unit - mass_sq) / max(abs(mass_sq), 1.0e-300)
    threshold_rows = {}
    for threshold in thresholds:
        threshold_rows[f"abs_jt_lt_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask & (np.abs(jt_unit) < threshold))),
            "measure_fraction": measure_fraction(np.abs(jt_unit) < threshold, measure, mask),
        }
        threshold_rows[f"abs_js_minus_m2_over_m2_gt_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask & (np.abs(rel_js) > threshold))),
            "measure_fraction": measure_fraction(np.abs(rel_js) > threshold, measure, mask),
        }
    return {
        "name": mask_name,
        "count": int(np.count_nonzero(mask)),
        "measure_sum": float(np.sum(measure[mask])),
        "jt_unit": scalar_stats(jt_unit[mask]),
        "abs_jt_unit": scalar_stats(np.abs(jt_unit[mask])),
        "jt_positive_count": int(np.count_nonzero(mask & (jt_unit > 0.0))),
        "jt_negative_count": int(np.count_nonzero(mask & (jt_unit < 0.0))),
        "jt_sign_change_edge_count": sign_change_edge_count(jt_unit, mask),
        "js_unit": scalar_stats(js_unit[mask]),
        "relative_js_minus_m2": scalar_stats(rel_js[mask]),
        "abs_relative_js_minus_m2": scalar_stats(np.abs(rel_js[mask])),
        "n_phase_measure": scalar_stats((mass_sq * measure)[mask]),
        "thresholds": threshold_rows,
    }


def analyze_slice(
    *,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t_old: float,
    t: float,
    mass: float,
    rho_floor: float,
    support_rho_frac: float,
    support_measure_frac: float,
    trusted_erosion: int,
    thresholds: list[float],
    projection_mode: str,
) -> dict[str, object]:
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
    rho = real_array(bohm["rho"])
    measure = tilde_measure_from_bohm(bohm, mass)
    metric_inv = real_array(snap["inverse"]["tilde_ginv_txz"])
    metric_cov = real_array(snap["cov_txz"])
    rho_tilde = tilde_density_from_measure(measure, metric_cov, rho_floor)
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    u_cov = np.stack(
        [
            real_array(bohm["s_t"]),
            real_array(bohm["s_x"]),
            real_array(bohm["s_z"]),
        ],
        axis=-1,
    )
    j_unit = np.einsum("...mn,...n->...m", metric_inv, u_cov)
    jt_unit = j_unit[..., 0]
    js_unit = np.einsum("...m,...m->...", u_cov, j_unit)
    support = (rho > support_rho_frac * float(np.max(rho))) & (
        measure > support_measure_frac * float(np.max(measure))
    )
    trusted = erode_mask_8(support, trusted_erosion)
    if not np.any(trusted):
        trusted = support.copy()
    mass_sq = float(mass * mass)
    if projection_mode == "weighted":
        phase_proj = project_to_weighted_gradient(
            real_array(bohm["s_x"]),
            real_array(bohm["s_z"]),
            float(x[1] - x[0]),
            float(z[1] - z[0]),
            weight=np.where(trusted, measure, 0.0),
            mean_mode="weighted",
        )
    elif projection_mode == "periodic":
        phase_proj = project_to_gradient(
            real_array(bohm["s_x"]),
            real_array(bohm["s_z"]),
            float(x[1] - x[0]),
            float(z[1] - z[0]),
        )
    else:
        raise ValueError("projection_mode must be 'weighted' or 'periodic'")
    return {
        "t_old": float(t_old),
        "t": float(t),
        "t_fs": ev_inv_to_fs(t),
        "rho": rho,
        "measure": measure,
        "rho_tilde": rho_tilde,
        "sqrt_abs_g": sqrt_abs_g,
        "jt_unit": jt_unit,
        "js_unit": js_unit,
        "phase_projection_rel_resid": phase_proj["rel_resid"],
        "phase_projection_weighted_rel_resid": float(phase_proj.get("weighted_rel_resid", np.nan)),
        "support": support,
        "trusted": trusted,
        "summary": {
            "t_old": float(t_old),
            "t": float(t),
            "t_fs": ev_inv_to_fs(t),
            "support": summarize_mask(
                mask_name="support",
                mask=support,
                measure=measure,
                jt_unit=jt_unit,
                js_unit=js_unit,
                mass_sq=mass_sq,
                thresholds=thresholds,
            ),
            "phase_projection_rel_resid_support": scalar_stats(phase_proj["rel_resid"][support]),
            "phase_projection_weighted_rel_resid": float(phase_proj.get("weighted_rel_resid", np.nan)),
            "trusted": summarize_mask(
                mask_name="trusted",
                mask=trusted,
                measure=measure,
                jt_unit=jt_unit,
                js_unit=js_unit,
                mass_sq=mass_sq,
                thresholds=thresholds,
            ),
            "phase_projection_rel_resid_trusted": scalar_stats(phase_proj["rel_resid"][trusted]),
        },
    }


def render_plot(out_path: Path, x: np.ndarray, z: np.ndarray, slices: list[dict[str, object]], mass_sq: float) -> None:
    plot_factor = HBAR_C_EV_M * 1.0e6
    x_plot = x * plot_factor
    z_plot = z * plot_factor
    xg, zg = np.meshgrid(x_plot, z_plot, indexing="ij")
    n = len(slices)
    fig, axes = plt.subplots(n, 4, figsize=(20.0, 4.6 * n), constrained_layout=True)
    if n == 1:
        axes = axes[None, :]
    for row, item in enumerate(slices):
        rho = np.asarray(item["rho"], dtype=float)
        support = np.asarray(item["support"], dtype=bool)
        trusted = np.asarray(item["trusted"], dtype=bool)
        jt = np.asarray(item["jt_unit"], dtype=float)
        js = np.asarray(item["js_unit"], dtype=float)
        phase_resid = np.asarray(item["phase_projection_rel_resid"], dtype=float)
        title = f"t_old={float(item['t_old']):g}, t={float(item['t_fs']):.3g} fs"

        vmax_rho = max(float(np.percentile(rho[support], 99.5)) if np.any(support) else float(np.max(rho)), 1.0e-16)
        im0 = axes[row, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax_rho)
        axes[row, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7)
        axes[row, 0].contour(xg, zg, trusted.astype(float), levels=[0.5], colors="cyan", linewidths=0.7)
        axes[row, 0].set_title(f"{title}: rho_A")
        fig.colorbar(im0, ax=axes[row, 0])

        clip_jt = max(float(np.percentile(np.abs(jt[support]), 95.0)) if np.any(support) else float(np.max(np.abs(jt))), 1.0e-12)
        im1 = axes[row, 1].pcolormesh(
            xg,
            zg,
            np.clip(jt, -clip_jt, clip_jt),
            shading="auto",
            cmap="coolwarm",
            vmin=-clip_jt,
            vmax=clip_jt,
        )
        axes[row, 1].contour(xg, zg, jt, levels=[0.0], colors="black", linewidths=0.8)
        axes[row, 1].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        axes[row, 1].set_title("lab projection j^t/rho_tilde")
        fig.colorbar(im1, ax=axes[row, 1])

        log_abs_jt = np.log10(np.maximum(np.abs(jt), 1.0e-16))
        im2 = axes[row, 2].pcolormesh(xg, zg, log_abs_jt, shading="auto", cmap="magma")
        axes[row, 2].contour(xg, zg, jt, levels=[0.0], colors="white", linewidths=0.8)
        axes[row, 2].contour(xg, zg, support.astype(float), levels=[0.5], colors="cyan", linewidths=0.6)
        axes[row, 2].set_title("log10 |j^t/rho_tilde|")
        fig.colorbar(im2, ax=axes[row, 2])

        rel_js = (js - mass_sq) / max(abs(mass_sq), 1.0e-300)
        vmax_rel = max(float(np.percentile(np.abs(rel_js[support]), 99.0)) if np.any(support) else float(np.max(np.abs(rel_js))), 1.0e-14)
        im3 = axes[row, 3].pcolormesh(
            xg,
            zg,
            np.log10(np.maximum(phase_resid, 1.0e-16)),
            shading="auto",
            cmap="magma",
        )
        axes[row, 3].contour(xg, zg, support.astype(float), levels=[0.5], colors="black", linewidths=0.6)
        axes[row, 3].set_title("log10 phase projection residual")
        fig.colorbar(im3, ax=axes[row, 3])

        for ax in axes[row, :]:
            ax.set_xlabel("x [um]")
            ax.set_ylabel("z [um]")
            ax.set_aspect("equal")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    params, physical_scale = build_crossing_params(args)
    if not args.physical_optical or physical_scale is None:
        raise ValueError("This diagnostic expects --physical-optical for old-time labels.")

    x, z, x_grid, z_grid = make_grid(params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(x_grid, z_grid, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_before = kg_positive_frequency_norm(psi0, psi0_t, dx, dz)
    if args.normalize_kg:
        psi0 = psi0 / np.sqrt(max(kg_norm_before, 1.0e-300))
        psi0_hat = np.fft.fft2(psi0)
        psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_after = kg_positive_frequency_norm(psi0, psi0_t, dx, dz)

    old_scale = physical_scale["old_dimensionless_scale_ev_inv"]
    times_old = parse_float_list(args.times_old)
    plot_times_old = parse_float_list(args.plot_times_old)
    thresholds = parse_float_list(args.thresholds)
    plot_set = set(float(t) for t in plot_times_old)
    slices = [
        analyze_slice(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t_old=t_old,
            t=t_old * old_scale,
            mass=params.m,
            rho_floor=args.rho_floor,
            support_rho_frac=args.support_rho_frac,
            support_measure_frac=args.support_measure_frac,
            trusted_erosion=args.trusted_erosion,
            thresholds=thresholds,
            projection_mode=args.projection_mode,
        )
        for t_old in times_old
    ]
    plot_slices = [item for item in slices if float(item["t_old"]) in plot_set]
    if not plot_slices:
        plot_slices = slices[: min(len(slices), 3)]
    fig_path = out / "d_phase_time_clock_summary.png"
    render_plot(fig_path, x, z, plot_slices, float(params.m * params.m))
    summary = {
        "params": {
            "resolution": args.resolution,
            "physical_optical": bool(args.physical_optical),
            "physical_scale": physical_scale,
            "normalize_kg": bool(args.normalize_kg),
            "kg_norm_before": kg_norm_before,
            "kg_norm_after": kg_norm_after,
            "support_rho_frac": args.support_rho_frac,
            "support_measure_frac": args.support_measure_frac,
            "trusted_erosion": args.trusted_erosion,
            "times_old": times_old,
            "plot_times_old": plot_times_old,
            "thresholds": thresholds,
            "mass_sq": float(params.m * params.m),
        },
        "definitions": {
            "jt_unit": "j^t/rho_tilde = gtilde^{t nu} u_nu.  Lab-time conserved density sqrt(|gtilde|)rho_tilde jt_unit is ill-conditioned where this crosses zero.",
            "js_unit": "j^S/rho_tilde = u_mu gtilde^{mu nu}u_nu for tau=S.  On the D mass shell this equals m^2, so the phase-time conserved density is m^2 sqrt(|gtilde|)rho_tilde.",
            "phase_projection_rel_resid": "Helmholtz projection residual for reconstructing a single-valued phase potential from (u_x,u_z). Small values mean the phase field is a good local clock variable on the slice.",
            "phase_time_clock_interpretation": "This is a coordinate/variable diagnostic, not a new force or damping term.  It suggests replacing a global lab-time Cauchy density by a phase-time or characteristic matter integrator.",
        },
        "slices": [item["summary"] for item in slices],
        "files": {
            "figure": str(fig_path.resolve()),
            "summary": str((out / "summary.json").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--physical-optical", action="store_true", default=True)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--times-old", default="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16")
    parser.add_argument("--plot-times-old", default="0,8,16")
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument(
        "--projection-mode",
        choices=["weighted", "periodic"],
        default="weighted",
        help="Phase reconstruction mode. weighted is the recommended local-support default.",
    )
    parser.add_argument("--thresholds", default="1e-12,1e-9,1e-6,1e-3,1e-2")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
