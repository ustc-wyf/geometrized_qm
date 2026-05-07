from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_bcd_residuals_from_a_reference import (
    branch_f_and_fr,
    covariant_hessian_from_jets,
    geometry_data,
    metric_jets_full,
    scalar_jets_full,
    stress_tensor_tilde,
    tensor_frobenius,
)
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def make_fields(
    branch: str,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    probe_dt: float,
    ell: float,
    mp: float,
    mass: float,
    rho_floor: float,
) -> dict[str, np.ndarray]:
    snap_m = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - probe_dt, rho_floor=rho_floor)
    snap_0 = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t, rho_floor=rho_floor)
    snap_p = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + probe_dt, rho_floor=rho_floor)

    metric_m = np.real_if_close(snap_m["cov_txz"])
    metric_0 = np.real_if_close(snap_0["cov_txz"])
    metric_p = np.real_if_close(snap_p["cov_txz"])
    rho = np.real_if_close(snap_0["bohm"]["rho"])
    s_t = np.real_if_close(snap_0["bohm"]["s_t"])
    s_x = np.real_if_close(snap_0["bohm"]["s_x"])
    s_z = np.real_if_close(snap_0["bohm"]["s_z"])

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, probe_dt, dx, dz)
    ginv_0, gamma2_0, ricci_0, r_tilde_0 = geometry_data(metric_0, dg, d2g)

    t_tensor, tilde_x = stress_tensor_tilde(metric_0, ginv_0, rho, s_t, s_x, s_z, m=mass)
    rhs_tensor = t_tensor / (mp * mp)

    f0, fr0 = branch_f_and_fr(branch, r_tilde_0, ell)

    if branch == "B":
        cov_hess_fr = np.zeros(metric_0.shape[:2] + (3, 3), dtype=float)
        box_fr = np.zeros(metric_0.shape[:2], dtype=float)
    else:
        dg_m, d2g_m = metric_jets_full(
            np.real_if_close(
                localized_direct_tilde_coordinate_snapshot(
                    psi0, psi0_hat, omega, x, z, t - 2.0 * probe_dt, rho_floor=rho_floor
                )["cov_txz"]
            ),
            metric_m,
            metric_0,
            probe_dt,
            dx,
            dz,
        )
        _, _, _, r_tilde_m = geometry_data(metric_m, dg_m, d2g_m)
        dg_p, d2g_p = metric_jets_full(
            metric_0,
            metric_p,
            np.real_if_close(
                localized_direct_tilde_coordinate_snapshot(
                    psi0, psi0_hat, omega, x, z, t + 2.0 * probe_dt, rho_floor=rho_floor
                )["cov_txz"]
            ),
            probe_dt,
            dx,
            dz,
        )
        _, _, _, r_tilde_p = geometry_data(metric_p, dg_p, d2g_p)
        _, fr_m = branch_f_and_fr(branch, r_tilde_m, ell)
        _, fr_p = branch_f_and_fr(branch, r_tilde_p, ell)
        dfr, d2fr = scalar_jets_full(fr_m, fr0, fr_p, probe_dt, dx, dz)
        cov_hess_fr = covariant_hessian_from_jets(dfr, d2fr, gamma2_0)
        box_fr = np.einsum("...ab,...ab->...", ginv_0, cov_hess_fr, optimize=True)

    algebraic_tensor = fr0[..., None, None] * ricci_0 - 0.5 * f0[..., None, None] * metric_0
    derivative_tensor = -(cov_hess_fr - metric_0 * box_fr[..., None, None])
    lhs_tensor = algebraic_tensor + derivative_tensor
    residual_tensor = lhs_tensor - rhs_tensor

    return {
        "rho": rho,
        "R_tilde": r_tilde_0,
        "ell2R_abs": np.abs((ell * ell) * r_tilde_0),
        "fR": fr0,
        "rhs_norm": tensor_frobenius(rhs_tensor),
        "lhs_norm": tensor_frobenius(lhs_tensor),
        "residual_norm": tensor_frobenius(residual_tensor),
        "algebraic_norm": tensor_frobenius(algebraic_tensor),
        "derivative_norm": tensor_frobenius(derivative_tensor),
        "tilde_X_minus_m2": tilde_x - mass * mass,
    }


def mask_centroid(mask: np.ndarray, x: np.ndarray, z: np.ndarray) -> tuple[float, float] | None:
    if not np.any(mask):
        return None
    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    weights = mask.astype(float)
    norm = float(np.sum(weights))
    return float(np.sum(weights * Xg) / norm), float(np.sum(weights * Zg) / norm)


def render_branch_time_plot(
    branch: str,
    t: float,
    x: np.ndarray,
    z: np.ndarray,
    fields: dict[str, np.ndarray],
    ell2_band: tuple[float, float],
    out_path: Path,
) -> dict[str, object]:
    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    rho = fields["rho"]
    log_rho = np.log10(np.maximum(rho, 1.0e-16))
    residual = fields["residual_norm"]
    derivative = fields["derivative_norm"]
    ell2r = fields["ell2R_abs"]

    support_mask = rho > 1.0e-3 * float(np.max(rho))
    dense_core_mask = rho > 1.0e-1 * float(np.max(rho))
    transition_mask = (ell2r >= ell2_band[0]) & (ell2r <= ell2_band[1])

    hotspot_threshold = float(np.percentile(residual[support_mask], 99.0)) if np.any(support_mask) else float(np.max(residual))
    hotspot_mask = residual >= hotspot_threshold
    support_transition = transition_mask & support_mask
    hotspot_transition = hotspot_mask & transition_mask

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.2), constrained_layout=True)

    im0 = axes[0, 0].pcolormesh(Xg, Zg, log_rho, shading="auto", cmap="viridis")
    axes[0, 0].contour(Xg, Zg, transition_mask.astype(float), levels=[0.5], colors="w", linewidths=1.2)
    axes[0, 0].set_title(fr"{branch} branch, $t={t:.1f}$: log10 rho with transition layer")
    fig.colorbar(im0, ax=axes[0, 0], label=r"$\log_{10}\rho$")

    im1 = axes[0, 1].pcolormesh(Xg, Zg, np.log10(np.maximum(ell2r, 1.0e-16)), shading="auto", cmap="magma")
    axes[0, 1].contour(Xg, Zg, transition_mask.astype(float), levels=[0.5], colors="c", linewidths=1.2)
    axes[0, 1].set_title(fr"$\log_{{10}}|\ell^2\tilde R|$")
    fig.colorbar(im1, ax=axes[0, 1], label=r"$\log_{10}|\ell^2\tilde R|$")

    im2 = axes[1, 0].pcolormesh(Xg, Zg, np.log10(np.maximum(derivative, 1.0e-30)), shading="auto", cmap="inferno")
    axes[1, 0].contour(Xg, Zg, transition_mask.astype(float), levels=[0.5], colors="w", linewidths=1.0)
    axes[1, 0].contour(Xg, Zg, hotspot_mask.astype(float), levels=[0.5], colors="lime", linewidths=1.0)
    axes[1, 0].set_title("Derivative-term norm with residual hotspots")
    fig.colorbar(im2, ax=axes[1, 0], label=r"$\log_{10}\|D(f_R)\|_F$")

    im3 = axes[1, 1].pcolormesh(Xg, Zg, np.log10(np.maximum(residual, 1.0e-30)), shading="auto", cmap="plasma")
    axes[1, 1].contour(Xg, Zg, transition_mask.astype(float), levels=[0.5], colors="w", linewidths=1.0)
    axes[1, 1].contour(Xg, Zg, dense_core_mask.astype(float), levels=[0.5], colors="lime", linewidths=1.0, linestyles="--")
    axes[1, 1].set_title("Residual norm")
    fig.colorbar(im3, ax=axes[1, 1], label=r"$\log_{10}\|E\|_F$")

    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)

    return {
        "time": float(t),
        "transition_fraction_global": float(np.mean(transition_mask)),
        "transition_fraction_in_support": float(np.mean(support_transition[support_mask])) if np.any(support_mask) else 0.0,
        "hotspot_fraction_in_support": float(np.mean(hotspot_mask[support_mask])) if np.any(support_mask) else 0.0,
        "hotspot_inside_transition_fraction": float(np.mean(hotspot_transition[hotspot_mask])) if np.any(hotspot_mask) else 0.0,
        "transition_centroid_global": mask_centroid(transition_mask, x, z),
        "transition_centroid_in_support": mask_centroid(support_transition, x, z),
        "hotspot_centroid": mask_centroid(hotspot_mask & support_mask, x, z),
        "hotspot_threshold": hotspot_threshold,
        "ell2_band": [float(ell2_band[0]), float(ell2_band[1])],
        "figure": str(out_path.resolve()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--nx", type=int, default=96)
    parser.add_argument("--nz", type=int, default=96)
    parser.add_argument("--times", type=str, default="8,16")
    parser.add_argument("--ell2-band", type=str, default="0.5,2.0")
    args = parser.parse_args()

    ell2_band = tuple(float(tok) for tok in args.ell2_band.split(","))
    params = FlatLocalizedCrossingParams(nx=args.nx, nz=args.nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    times = [float(tok) for tok in args.times.split(",") if tok.strip()]

    args.output.mkdir(parents=True, exist_ok=True)
    records: dict[str, list[dict[str, object]]] = {"C": [], "D": []}
    for branch in ("C", "D"):
        for t in times:
            fields = make_fields(
                branch=branch,
                psi0=psi0,
                psi0_hat=psi0_hat,
                omega=omega,
                x=x,
                z=z,
                t=t,
                probe_dt=args.probe_dt,
                ell=args.ell,
                mp=args.mp,
                mass=params.m,
                rho_floor=args.rho_floor,
            )
            rec = render_branch_time_plot(
                branch=branch,
                t=t,
                x=x,
                z=z,
                fields=fields,
                ell2_band=ell2_band,
                out_path=args.output / f"{branch.lower()}_transition_t{str(t).replace('.', 'p')}.png",
            )
            records[branch].append(rec)

    summary = {
        "params": {
            "ell": args.ell,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "nx": args.nx,
            "nz": args.nz,
            "times": times,
            "ell2_band": list(ell2_band),
        },
        "definitions": {
            "support_mask": "rho > 1e-3 * rho_max",
            "dense_core_mask": "rho > 1e-1 * rho_max",
            "transition_mask": f"{ell2_band[0]} <= |ell^2 R_tilde| <= {ell2_band[1]}",
            "hotspot_mask": "residual_norm in top 1% within support_mask",
        },
        "records": records,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
