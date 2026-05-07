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
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import (
    build_physical_reference,
    make_cropped_snapshot,
    real_array,
)
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det


def stats(values: np.ndarray | list[float]) -> dict[str, float]:
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


def tensor_norm(tensor: np.ndarray) -> np.ndarray:
    return np.sqrt(np.einsum("...ab,...ab->...", tensor, tensor, optimize=True))


def branch_d_f_phi(r_scalar: np.ndarray, ell: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw_y = (ell * ell) * r_scalar
    tanh_y = np.tanh(raw_y)
    f = tanh_y / (ell * ell)
    phi = 1.0 - tanh_y * tanh_y
    return np.real_if_close(f), np.real_if_close(phi), np.real_if_close(raw_y)


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


def render_terms(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    raw_y: np.ndarray,
    algebraic_norm: np.ndarray,
    derivative_norm: np.ndarray,
    residual_relative: np.ndarray,
    support: np.ndarray,
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
    axes[0, 0].set_title("rho, white=support")
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046)

    im1 = axes[0, 1].imshow((np.sign(raw_y) * np.log10(1.0 + np.abs(raw_y))).T, origin="lower", extent=extent, aspect="equal", cmap="coolwarm")
    axes[0, 1].set_title("sign(raw_y) log10(1+|raw_y|)")
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046)

    im2 = axes[0, 2].imshow(np.log10(1.0e-300 + algebraic_norm).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
    axes[0, 2].set_title("log10 ||phi R - 1/2 f g||")
    fig.colorbar(im2, ax=axes[0, 2], fraction=0.046)

    im3 = axes[1, 0].imshow(np.log10(1.0e-300 + derivative_norm).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
    axes[1, 0].set_title("log10 ||(g Box - nabla nabla)phi||")
    fig.colorbar(im3, ax=axes[1, 0], fraction=0.046)

    ratio = derivative_norm / np.maximum(algebraic_norm, 1.0e-300)
    im4 = axes[1, 1].imshow(np.log10(1.0e-300 + ratio).T, origin="lower", extent=extent, aspect="equal", cmap="cividis")
    axes[1, 1].set_title("log10 derivative/algebraic")
    fig.colorbar(im4, ax=axes[1, 1], fraction=0.046)

    im5 = axes[1, 2].imshow(np.log10(1.0e-300 + residual_relative).T, origin="lower", extent=extent, aspect="equal", cmap="inferno")
    axes[1, 2].set_title("log10 residual/(lhs+rhs)")
    fig.colorbar(im5, ax=axes[1, 2], fraction=0.046)

    for ax in axes.ravel():
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
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
    f0, phi0, raw_y = branch_d_f_phi(r_scalar, float(args.ell))

    dg_m, d2g_m = metric_jets_full(metric_m2, metric_m, metric_0, dt, dx, dz)
    _, _, _, r_scalar_m = geometry_data(metric_m, dg_m, d2g_m)
    _, phi_m, _ = branch_d_f_phi(r_scalar_m, float(args.ell))
    dg_p, d2g_p = metric_jets_full(metric_0, metric_p, metric_p2, dt, dx, dz)
    _, _, _, r_scalar_p = geometry_data(metric_p, dg_p, d2g_p)
    _, phi_p, _ = branch_d_f_phi(r_scalar_p, float(args.ell))

    dphi, d2phi = scalar_jets_full(phi_m, phi0, phi_p, dt, dx, dz)
    cov_hess_phi = covariant_hessian_from_jets(dphi, d2phi, gamma2)
    box_phi = np.einsum("...ab,...ab->...", metric_inv, cov_hess_phi, optimize=True)

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

    algebraic_tensor = phi0[..., None, None] * ricci - 0.5 * f0[..., None, None] * metric_0
    derivative_tensor = metric_0 * box_phi[..., None, None] - cov_hess_phi
    lhs_tensor = algebraic_tensor + derivative_tensor
    rhs_tensor = real_array(stress) / (float(args.mp) * float(args.mp))
    residual_tensor = lhs_tensor - rhs_tensor

    algebraic_norm = tensor_norm(algebraic_tensor)
    derivative_norm = tensor_norm(derivative_tensor)
    lhs_norm = tensor_norm(lhs_tensor)
    rhs_norm = tensor_norm(rhs_tensor)
    residual_norm = tensor_norm(residual_tensor)
    relative_residual = residual_norm / np.maximum(lhs_norm + rhs_norm, 1.0e-300)
    derivative_to_algebraic = derivative_norm / np.maximum(algebraic_norm, 1.0e-300)

    support = (rho_a > float(args.support_rho_frac) * float(np.max(rho_a))) & (
        measure > float(args.support_measure_frac) * float(np.max(measure))
    )
    trusted = erode_mask_8(support, int(args.trusted_erosion))
    if not np.any(trusted):
        trusted = support.copy()
    clipped = np.clip(np.abs(raw_y), 0.0, 350.0)
    f_r = 1.0 / (np.cosh(clipped) ** 2)
    saturated = support & (f_r <= float(args.saturation_phi))

    plot_path = out / "d_full_fr_equation_terms.png"
    render_terms(plot_path, x, z, rho_a, raw_y, algebraic_norm, derivative_norm, relative_residual, support)
    np.savez_compressed(
        out / "d_full_fr_equation_fields.npz",
        x=x,
        z=z,
        rho_A=rho_a,
        support=support,
        trusted=trusted,
        saturated=saturated,
        raw_y=raw_y,
        R_tilde=r_scalar,
        phi=phi0,
        f=f0,
        algebraic_norm=algebraic_norm,
        derivative_norm=derivative_norm,
        lhs_norm=lhs_norm,
        rhs_norm=rhs_norm,
        residual_norm=residual_norm,
        relative_residual=relative_residual,
        derivative_to_algebraic=derivative_to_algebraic,
        mass_shell_defect=real_array(tilde_x) - float(params.m) * float(params.m),
    )

    masks = {
        "support": support,
        "trusted": trusted,
        "saturated": saturated,
    }
    term_summary: dict[str, object] = {}
    for name, mask in masks.items():
        term_summary[name] = {
            "count": int(np.count_nonzero(mask)),
            "raw_y_abs": stats(raw_y[mask]),
            "phi": stats(phi0[mask]),
            "R_tilde": stats(r_scalar[mask]),
            "algebraic_norm": stats(algebraic_norm[mask]),
            "derivative_norm": stats(derivative_norm[mask]),
            "derivative_to_algebraic": stats(derivative_to_algebraic[mask]),
            "lhs_norm": stats(lhs_norm[mask]),
            "rhs_norm": stats(rhs_norm[mask]),
            "residual_norm": stats(residual_norm[mask]),
            "relative_residual_lhs_plus_rhs": stats(relative_residual[mask]),
            "mass_shell_defect": stats((real_array(tilde_x) - float(params.m) * float(params.m))[mask]),
        }

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
            "saturation_phi": float(args.saturation_phi),
        },
        "equation": {
            "lhs": "phi*Ricci_mn - 1/2*f*g_mn + (g_mn*Box - nabla_m nabla_n)phi",
            "rhs": "T_mn/M_P^2",
            "phi": "f_R = sech^2(ell^2 R_tilde)",
            "f": "tanh(ell^2 R_tilde)/ell^2",
            "residual_relative": "||lhs-rhs||/(||lhs||+||rhs||)",
        },
        "summary": term_summary,
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "fields_npz": str((out / "d_full_fr_equation_fields.npz").resolve()),
            "terms_png": str(plot_path.resolve()),
        },
        "interpretation": (
            "This checks the full metric f(R) equation on the A-derived transformed metric history. "
            "It is a residual diagnostic, not yet an implicit equality-constrained solve."
        ),
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
    parser.add_argument("--saturation-phi", type=float, default=1.0e-3)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
