from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
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


def real_array(value: np.ndarray) -> np.ndarray:
    return np.asarray(np.real_if_close(value), dtype=float)


def parse_float_list(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def scalar_stats(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p01": 0.0, "p05": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p01": float(np.percentile(vals, 1.0)),
        "p05": float(np.percentile(vals, 5.0)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }


def sign_change_edge_count(field: np.ndarray, mask: np.ndarray) -> int:
    f = np.asarray(field, dtype=float)
    m = np.asarray(mask, dtype=bool)
    x_edges = m[:-1, :] & m[1:, :] & (f[:-1, :] * f[1:, :] <= 0.0)
    z_edges = m[:, :-1] & m[:, 1:] & (f[:, :-1] * f[:, 1:] <= 0.0)
    return int(np.count_nonzero(x_edges) + np.count_nonzero(z_edges))


def mass_shell_discriminant(metric_inv: np.ndarray, ux: np.ndarray, uz: np.ndarray, mass: float) -> np.ndarray:
    a = metric_inv[..., 0, 0]
    b = metric_inv[..., 0, 1] * ux + metric_inv[..., 0, 2] * uz
    c = (
        metric_inv[..., 1, 1] * ux * ux
        + 2.0 * metric_inv[..., 1, 2] * ux * uz
        + metric_inv[..., 2, 2] * uz * uz
        - mass * mass
    )
    return b * b - a * c


def measure_fraction(values: np.ndarray, measure: np.ndarray, mask: np.ndarray) -> float:
    denom = max(float(np.sum(measure[mask])), 1.0e-300)
    return float(np.sum(measure[mask & values]) / denom)


def analyze_slice(
    *,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    time_label: str,
    mass: float,
    rho_floor: float,
    support_rho_frac: float,
    support_measure_frac: float,
    trusted_erosion: int,
    jt_thresholds: list[float],
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
    metric_inv = real_array(snap["inverse"]["tilde_ginv_txz"])
    metric_cov = real_array(snap["cov_txz"])
    measure = tilde_measure_from_bohm(bohm, mass)
    rho_tilde = tilde_density_from_measure(measure, metric_cov, rho_floor)
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    ux = real_array(bohm["s_x"])
    uz = real_array(bohm["s_z"])
    ut = real_array(bohm["s_t"])
    jt_unit = metric_inv[..., 0, 0] * ut + metric_inv[..., 0, 1] * ux + metric_inv[..., 0, 2] * uz
    n_cons = sqrt_abs_g * rho_tilde * jt_unit
    disc = mass_shell_discriminant(metric_inv, ux, uz, mass)

    support = (rho > support_rho_frac * float(np.max(rho))) & (
        measure > support_measure_frac * float(np.max(measure))
    )
    trusted = erode_mask_8(support, trusted_erosion)
    if not np.any(trusted):
        trusted = support.copy()

    def mask_summary(name: str, mask: np.ndarray) -> dict[str, object]:
        total_measure = max(float(np.sum(measure[mask])), 1.0e-300)
        total_rho = max(float(np.sum(rho[mask])), 1.0e-300)
        threshold_rows = {}
        for threshold in jt_thresholds:
            near = np.abs(jt_unit) < threshold
            threshold_rows[f"abs_jt_lt_{threshold:g}"] = {
                "count": int(np.count_nonzero(mask & near)),
                "count_fraction": float(np.count_nonzero(mask & near) / max(int(np.count_nonzero(mask)), 1)),
                "measure_fraction": float(np.sum(measure[mask & near]) / total_measure),
                "rho_fraction": float(np.sum(rho[mask & near]) / total_rho),
            }
        return {
            "name": name,
            "count": int(np.count_nonzero(mask)),
            "measure_sum": float(np.sum(measure[mask])),
            "rho_sum": float(np.sum(rho[mask])),
            "jt_unit": scalar_stats(jt_unit[mask]),
            "abs_jt_unit": scalar_stats(np.abs(jt_unit[mask])),
            "n_cons": scalar_stats(n_cons[mask]),
            "rho_tilde": scalar_stats(rho_tilde[mask]),
            "disc": scalar_stats(disc[mask]),
            "jt_positive_count": int(np.count_nonzero(mask & (jt_unit > 0.0))),
            "jt_negative_count": int(np.count_nonzero(mask & (jt_unit < 0.0))),
            "jt_positive_measure_fraction": measure_fraction(jt_unit > 0.0, measure, mask),
            "jt_negative_measure_fraction": measure_fraction(jt_unit < 0.0, measure, mask),
            "jt_sign_change_edge_count": sign_change_edge_count(jt_unit, mask),
            "thresholds": threshold_rows,
        }

    return {
        "time_label": time_label,
        "t": float(t),
        "t_fs": ev_inv_to_fs(t),
        "rho": rho,
        "measure": measure,
        "rho_tilde": rho_tilde,
        "jt_unit": jt_unit,
        "n_cons": n_cons,
        "disc": disc,
        "support": support,
        "trusted": trusted,
        "summary": {
            "time_label": time_label,
            "t": float(t),
            "t_fs": ev_inv_to_fs(t),
            "global": mask_summary("global", np.ones_like(rho, dtype=bool)),
            "support": mask_summary("support", support),
            "trusted": mask_summary("trusted", trusted),
        },
    }


def render_summary_plot(out_path: Path, x: np.ndarray, z: np.ndarray, slices: list[dict[str, object]]) -> None:
    plot_factor = HBAR_C_EV_M * 1.0e6
    x_plot = x * plot_factor
    z_plot = z * plot_factor
    xg, zg = np.meshgrid(x_plot, z_plot, indexing="ij")
    n = len(slices)
    fig, axes = plt.subplots(n, 3, figsize=(15.5, 4.4 * n), constrained_layout=True)
    if n == 1:
        axes = axes[None, :]

    for row, item in enumerate(slices):
        rho = item["rho"]
        measure = item["measure"]
        jt = item["jt_unit"]
        support = item["support"]
        trusted = item["trusted"]
        title = str(item["time_label"])

        vmax_rho = max(float(np.percentile(rho[support], 99.5)) if np.any(support) else float(np.max(rho)), 1.0e-16)
        im0 = axes[row, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax_rho)
        axes[row, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7)
        axes[row, 0].contour(xg, zg, trusted.astype(float), levels=[0.5], colors="cyan", linewidths=0.7)
        axes[row, 0].set_title(f"{title}: rho_A")
        fig.colorbar(im0, ax=axes[row, 0])

        clip = max(float(np.percentile(np.abs(jt[support]), 95.0)) if np.any(support) else float(np.max(np.abs(jt))), 1.0e-12)
        im1 = axes[row, 1].pcolormesh(
            xg,
            zg,
            np.clip(jt, -clip, clip),
            shading="auto",
            cmap="coolwarm",
            vmin=-clip,
            vmax=clip,
        )
        axes[row, 1].contour(xg, zg, jt, levels=[0.0], colors="black", linewidths=0.9)
        axes[row, 1].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        axes[row, 1].set_title(f"{title}: jt_unit = gtilde^{{t nu}} u_nu")
        fig.colorbar(im1, ax=axes[row, 1])

        log_abs = np.log10(np.maximum(np.abs(jt), 1.0e-16))
        im2 = axes[row, 2].pcolormesh(xg, zg, log_abs, shading="auto", cmap="magma")
        axes[row, 2].contour(xg, zg, jt, levels=[0.0], colors="white", linewidths=0.8)
        axes[row, 2].contour(xg, zg, support.astype(float), levels=[0.5], colors="cyan", linewidths=0.6)
        axes[row, 2].set_title(f"{title}: log10 |jt_unit|")
        fig.colorbar(im2, ax=axes[row, 2])

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
        raise ValueError("This diagnostic currently expects --physical-optical so old-time labels are meaningful.")

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
    time_old_values = parse_float_list(args.times_old)
    thresholds = parse_float_list(args.jt_thresholds)
    slices = [
        analyze_slice(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t=t_old * old_scale,
            time_label=f"t_old={t_old:g}",
            mass=params.m,
            rho_floor=args.rho_floor,
            support_rho_frac=args.support_rho_frac,
            support_measure_frac=args.support_measure_frac,
            trusted_erosion=args.trusted_erosion,
            jt_thresholds=thresholds,
        )
        for t_old in time_old_values
    ]

    fig_path = out / "jt_degeneracy_summary.png"
    render_summary_plot(fig_path, x, z, slices)
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
            "jt_thresholds": thresholds,
            "times_old": time_old_values,
        },
        "definitions": {
            "jt_unit": "gtilde^{t nu} u_nu.  The coordinate conserved density is sqrt(|gtilde|) rho_tilde jt_unit; jt_unit near zero makes rho_tilde recovery from n_cons ill-conditioned.",
            "support": "rho and transformed-measure binary support at the requested fractions of their maxima.",
            "trusted": "support eroded by trusted_erosion cells.",
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
    parser.add_argument("--times-old", default="0,8,16")
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument("--jt-thresholds", default="1e-6,1e-4,1e-3,1e-2")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
