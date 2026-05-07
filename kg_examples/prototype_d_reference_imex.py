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
    stress_tensor_tilde,
    tensor_frobenius,
)
from analyze_c_terms_from_a_reference import metric_jets_full, scalar_jets_full
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def d_branch_signed_fields(chi: np.ndarray, ell: float) -> dict[str, np.ndarray]:
    l2 = ell * ell
    raw_y = l2 * chi
    tanh_raw_y = np.tanh(raw_y)
    phi = 1.0 - tanh_raw_y * tanh_raw_y
    f = tanh_raw_y / l2
    # Auxiliary-potential density for the scalar-tensor rewrite:
    #   U(chi) = chi * phi(chi) - f(chi)
    u = chi * phi - f
    phi_chi = -2.0 * l2 * tanh_raw_y * phi
    return {
        # Keep "y" for backward compatibility with older diagnostic scripts,
        # but the actual interface variable is raw_y = ell^2 * R_tilde.
        "y": np.real_if_close(tanh_raw_y),
        "raw_y": np.real_if_close(raw_y),
        "tanh_raw_y": np.real_if_close(tanh_raw_y),
        "phi": np.real_if_close(phi),
        "f": np.real_if_close(f),
        "U": np.real_if_close(u),
        "phi_chi": np.real_if_close(phi_chi),
    }


def scalar_covariant_box_from_phi(
    phi_m: np.ndarray,
    phi_0: np.ndarray,
    phi_p: np.ndarray,
    metric_m: np.ndarray,
    metric_0: np.ndarray,
    metric_p: np.ndarray,
    probe_dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, probe_dt, dx, dz)
    ginv_0, gamma2_0, _, _ = geometry_data(metric_0, dg, d2g)
    dphi, d2phi = scalar_jets_full(phi_m, phi_0, phi_p, probe_dt, dx, dz)
    cov_hess_phi = covariant_hessian_from_jets(dphi, d2phi, gamma2_0)
    box_phi = np.einsum("...ab,...ab->...", ginv_0, cov_hess_phi, optimize=True)
    return np.real_if_close(box_phi), np.real_if_close(cov_hess_phi)


def trace_residual_2p1_from_chi(
    chi_m: np.ndarray,
    chi_0: np.ndarray,
    chi_p: np.ndarray,
    metric_m: np.ndarray,
    metric_0: np.ndarray,
    metric_p: np.ndarray,
    stress_trace: np.ndarray,
    ell: float,
    probe_dt: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    d_m = d_branch_signed_fields(chi_m, ell)
    d_0 = d_branch_signed_fields(chi_0, ell)
    d_p = d_branch_signed_fields(chi_p, ell)
    box_phi, cov_hess_phi = scalar_covariant_box_from_phi(
        d_m["phi"],
        d_0["phi"],
        d_p["phi"],
        metric_m,
        metric_0,
        metric_p,
        probe_dt,
        dx,
        dz,
    )
    trace_residual = d_0["phi"] * chi_0 - 1.5 * d_0["f"] + 2.0 * box_phi - stress_trace
    return {
        "trace_residual_2p1": np.real_if_close(trace_residual),
        "box_phi": np.real_if_close(box_phi),
        "cov_hess_phi": np.real_if_close(cov_hess_phi),
        "phi": d_0["phi"],
        "f": d_0["f"],
        "U": d_0["U"],
        "phi_chi": d_0["phi_chi"],
    }


def reference_snapshot_data(
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
    snap_mm = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - 2.0 * probe_dt, rho_floor=rho_floor)
    snap_m = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t - probe_dt, rho_floor=rho_floor)
    snap_0 = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t, rho_floor=rho_floor)
    snap_p = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + probe_dt, rho_floor=rho_floor)
    snap_pp = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, t + 2.0 * probe_dt, rho_floor=rho_floor)

    metric_mm = np.real_if_close(snap_mm["cov_txz"])
    metric_m = np.real_if_close(snap_m["cov_txz"])
    metric_0 = np.real_if_close(snap_0["cov_txz"])
    metric_p = np.real_if_close(snap_p["cov_txz"])
    metric_pp = np.real_if_close(snap_pp["cov_txz"])
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, probe_dt, dx, dz)
    ginv_0, gamma2_0, ricci_0, r_tilde_0 = geometry_data(metric_0, dg, d2g)

    bohm = snap_0["bohm"]
    rho = np.real_if_close(bohm["rho"])
    s_t = np.real_if_close(bohm["s_t"])
    s_x = np.real_if_close(bohm["s_x"])
    s_z = np.real_if_close(bohm["s_z"])

    t_tensor, tilde_x = stress_tensor_tilde(metric_0, ginv_0, rho, s_t, s_x, s_z, m=mass)
    rhs_tensor = t_tensor / (mp * mp)

    f_0, fr_0 = branch_f_and_fr("D", r_tilde_0, ell)
    _, _, _, r_tilde_m = geometry_data(metric_m, *metric_jets_full(metric_mm, metric_m, metric_0, probe_dt, dx, dz))
    _, _, _, r_tilde_p = geometry_data(metric_p, *metric_jets_full(metric_0, metric_p, metric_pp, probe_dt, dx, dz))
    _, fr_m = branch_f_and_fr("D", r_tilde_m, ell)
    _, fr_p = branch_f_and_fr("D", r_tilde_p, ell)
    dfr, d2fr = scalar_jets_full(fr_m, fr_0, fr_p, probe_dt, dx, dz)
    cov_hess_fr = covariant_hessian_from_jets(dfr, d2fr, gamma2_0)
    box_fr = np.einsum("...ab,...ab->...", ginv_0, cov_hess_fr, optimize=True)

    algebraic_tensor = fr_0[..., None, None] * ricci_0 - 0.5 * f_0[..., None, None] * metric_0
    derivative_tensor = -(cov_hess_fr - metric_0 * box_fr[..., None, None])
    lhs_tensor = algebraic_tensor + derivative_tensor
    residual_tensor = lhs_tensor - rhs_tensor
    residual_norm = tensor_frobenius(residual_tensor)
    derivative_norm = tensor_frobenius(derivative_tensor)

    trace_residual = np.einsum("...ab,...ab->...", ginv_0, residual_tensor, optimize=True)

    d_signed = d_branch_signed_fields(r_tilde_0, ell)
    d_signed_m = d_branch_signed_fields(r_tilde_m, ell)
    d_signed_p = d_branch_signed_fields(r_tilde_p, ell)
    box_phi, _ = scalar_covariant_box_from_phi(
        d_signed_m["phi"],
        d_signed["phi"],
        d_signed_p["phi"],
        metric_m,
        metric_0,
        metric_p,
        probe_dt,
        dx,
        dz,
    )

    # In 2+1d the trace equation is
    #   f_R R - 3/2 f + 2 □ f_R = T / M_P^2.
    trace_rhs = np.einsum("...ab,...ab->...", ginv_0, rhs_tensor, optimize=True)
    trace_residual_2p1 = d_signed["phi"] * r_tilde_0 - 1.5 * d_signed["f"] + 2.0 * box_phi - trace_rhs

    return {
        "metric_cov": metric_0,
        "metric_cov_m": metric_m,
        "metric_cov_p": metric_p,
        "metric_inv": ginv_0,
        "rho": rho,
        "R_tilde": np.real_if_close(r_tilde_0),
        "tilde_X": np.real_if_close(tilde_x),
        "stress_trace": np.real_if_close(trace_rhs),
        "residual_tensor": np.real_if_close(residual_tensor),
        "residual_norm": np.real_if_close(residual_norm),
        "derivative_norm": np.real_if_close(derivative_norm),
        "trace_residual_tensor": np.real_if_close(trace_residual),
        "trace_residual_2p1": np.real_if_close(trace_residual_2p1),
        "chi_m": np.real_if_close(r_tilde_m),
        "chi_ref": np.real_if_close(r_tilde_0),
        "chi_p": np.real_if_close(r_tilde_p),
        "raw_y_ref": d_signed["raw_y"],
        "y_ref": d_signed["y"],
        "phi_ref": d_signed["phi"],
        "f_ref": d_signed["f"],
        "U_ref": d_signed["U"],
        "phi_chi_ref": d_signed["phi_chi"],
        "box_phi_ref": np.real_if_close(box_phi),
        "dx": dx,
        "dz": dz,
    }


def build_refine_mask(
    rho: np.ndarray,
    r_tilde: np.ndarray,
    phi_ref: np.ndarray,
    residual_norm: np.ndarray,
    derivative_norm: np.ndarray,
    dx: float,
    dz: float,
    ell: float,
) -> dict[str, np.ndarray | float]:
    rho_max = float(np.max(rho))
    support_mask = rho > 1.0e-3 * rho_max
    high_density_mask = rho > 1.0e-1 * rho_max
    l2r = np.abs((ell * ell) * r_tilde)
    broad_transition_mask = (l2r >= 1.0) & (l2r <= 100.0)
    near_transition_mask = np.abs(l2r - 1.0) <= 1.0

    grad_phi_x = np.gradient(phi_ref, dx, axis=0, edge_order=2)
    grad_phi_z = np.gradient(phi_ref, dz, axis=1, edge_order=2)
    grad_phi_norm = np.sqrt(grad_phi_x * grad_phi_x + grad_phi_z * grad_phi_z)

    support_vals = residual_norm[support_mask]
    residual_p95 = float(np.percentile(support_vals, 95.0)) if support_vals.size else 0.0
    residual_p99 = float(np.percentile(support_vals, 99.0)) if support_vals.size else 0.0
    grad_support_vals = grad_phi_norm[support_mask]
    grad_p95 = float(np.percentile(grad_support_vals, 95.0)) if grad_support_vals.size else 0.0
    grad_p99 = float(np.percentile(grad_support_vals, 99.0)) if grad_support_vals.size else 0.0
    deriv_support_vals = derivative_norm[support_mask]
    deriv_p95 = float(np.percentile(deriv_support_vals, 95.0)) if deriv_support_vals.size else 0.0
    deriv_p99 = float(np.percentile(deriv_support_vals, 99.0)) if deriv_support_vals.size else 0.0

    hotspot_mask = support_mask & (residual_norm > residual_p99)
    high_grad_mask = support_mask & (grad_phi_norm > grad_p99)
    high_derivative_mask = support_mask & (derivative_norm > deriv_p99)
    refine_mask = hotspot_mask | high_grad_mask | high_derivative_mask | (support_mask & near_transition_mask)

    return {
        "support_mask": support_mask,
        "high_density_mask": high_density_mask,
        "broad_transition_mask": broad_transition_mask,
        "near_transition_mask": near_transition_mask,
        "grad_phi_norm": grad_phi_norm,
        "hotspot_mask": hotspot_mask,
        "high_grad_mask": high_grad_mask,
        "high_derivative_mask": high_derivative_mask,
        "refine_mask": refine_mask,
        "residual_p95_support": residual_p95,
        "residual_p99_support": residual_p99,
        "grad_phi_p95_support": grad_p95,
        "grad_phi_p99_support": grad_p99,
        "derivative_p95_support": deriv_p95,
        "derivative_p99_support": deriv_p99,
    }


def imex_local_seed(
    chi_ref: np.ndarray,
    r_tilde: np.ndarray,
    trace_residual: np.ndarray,
    derivative_norm: np.ndarray,
    ell: float,
    dx: float,
    dz: float,
    safety: float = 0.2,
) -> dict[str, np.ndarray | float]:
    l2 = ell * ell
    y = np.tanh(l2 * chi_ref)
    phi = 1.0 - y * y
    phi_chi = -2.0 * l2 * y * phi
    # Local algebraic stiffness for the 2+1d trace equation at fixed reference metric:
    #   A(chi) = phi(chi) * R_ref - 3/2 f(chi)
    # so A'(chi) = phi_chi * R_ref - 3/2 phi.
    lambda_alg = np.abs(phi_chi * r_tilde - 1.5 * phi)
    # Explicit derivative piece is approximated by the already measured derivative norm
    # plus a coarse grid penalty that grows as the feature scale approaches dx,dz.
    grid_penalty = np.abs(phi_chi) * (1.0 / max(dx * dx, 1.0e-12) + 1.0 / max(dz * dz, 1.0e-12))
    lambda_exp = np.abs(derivative_norm) + grid_penalty
    # Avoid the flat-limit denominator collapse: the D branch still has a finite
    # curvature transition scale set by 1/ell^2 even when the local algebraic
    # derivative happens to vanish.
    lambda_floor = np.full_like(lambda_alg, 1.0 / max(ell * ell, 1.0e-12))
    lambda_total = np.maximum(lambda_alg + lambda_exp, lambda_floor)
    dt_local = safety / lambda_total
    delta_chi_seed = -trace_residual / lambda_total
    return {
        "lambda_alg": np.real_if_close(lambda_alg),
        "lambda_exp": np.real_if_close(lambda_exp),
        "lambda_floor": np.real_if_close(lambda_floor),
        "lambda_total": np.real_if_close(lambda_total),
        "dt_local": np.real_if_close(dt_local),
        "dt_global_suggested": float(np.min(dt_local)),
        "delta_chi_seed": np.real_if_close(delta_chi_seed),
    }


def array_stats(field: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.abs(field[mask])
    if vals.size == 0:
        return {"abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0}
    return {
        "abs_max": float(np.max(vals)),
        "abs_mean": float(np.mean(vals)),
        "abs_median": float(np.median(vals)),
        "p95": float(np.percentile(vals, 95.0)),
    }


def render_maps(
    outdir: Path,
    branch_data: dict[str, np.ndarray],
    masks: dict[str, np.ndarray | float],
    imex_seed: dict[str, np.ndarray | float],
    x: np.ndarray,
    z: np.ndarray,
    t: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.0), constrained_layout=True)
    panels = [
        (branch_data["rho"], r"$\rho$", "viridis"),
        (branch_data["R_tilde"], r"$\tilde R_{\rm ref}$", "coolwarm"),
        (branch_data["residual_norm"], r"$||\mathcal{E}_{ab}^{(D)}||$", "magma"),
        (branch_data["phi_ref"], r"$\phi_{\rm ref}=f_R(\chi_{\rm ref})$", "plasma"),
        (branch_data["trace_residual_2p1"], r"$\mathcal{T}_{\rm ref}^{(2+1)}$", "coolwarm"),
        (imex_seed["delta_chi_seed"], r"$\delta\chi_{\rm seed}$", "coolwarm"),
    ]
    for ax, (field, title, cmap) in zip(axes.reshape(-1), panels):
        im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap)
        ax.contour(xg, zg, masks["refine_mask"].astype(float), levels=[0.5], colors="white", linewidths=0.8)
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, shrink=0.86)
    fig.suptitle(f"D branch scalar-tensor reference-correction prototype, t={t:.2f}", fontsize=12)
    fig.savefig(outdir / f"d_scalar_tensor_reference_t{t:.1f}.png", dpi=180)
    plt.close(fig)


def run(
    outdir: str | Path,
    ell: float = 100.0,
    mp: float = 300.0,
    rho_floor: float = 1.0e-12,
    probe_dt: float = 1.0e-3,
    nx: int = 96,
    nz: int = 96,
) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = FlatLocalizedCrossingParams(nx=nx, nz=nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    times = [0.0, 8.0, 16.0]
    snapshots: list[dict[str, object]] = []

    for t in times:
        ref = reference_snapshot_data(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t=t,
            probe_dt=probe_dt,
            ell=ell,
            mp=mp,
            mass=params.m,
            rho_floor=rho_floor,
        )
        masks = build_refine_mask(
            rho=ref["rho"],
            r_tilde=ref["R_tilde"],
            phi_ref=ref["phi_ref"],
            residual_norm=ref["residual_norm"],
            derivative_norm=ref["derivative_norm"],
            dx=ref["dx"],
            dz=ref["dz"],
            ell=ell,
        )
        seed = imex_local_seed(
            chi_ref=ref["chi_ref"],
            r_tilde=ref["R_tilde"],
            trace_residual=ref["trace_residual_2p1"],
            derivative_norm=ref["derivative_norm"],
            ell=ell,
            dx=ref["dx"],
            dz=ref["dz"],
        )

        render_maps(out, ref, masks, seed, x, z, t)
        snapshots.append(
            {
                "time": float(t),
                "support_residual_norm": array_stats(ref["residual_norm"], masks["support_mask"]),
                "support_trace_residual_2p1": array_stats(ref["trace_residual_2p1"], masks["support_mask"]),
                "support_delta_chi_seed": array_stats(seed["delta_chi_seed"], masks["support_mask"]),
                "support_lambda_total": array_stats(seed["lambda_total"], masks["support_mask"]),
                "refine_fraction_support": float(
                    np.mean(masks["refine_mask"][masks["support_mask"]]) if np.any(masks["support_mask"]) else 0.0
                ),
                "near_transition_fraction_support": float(
                    np.mean(masks["near_transition_mask"][masks["support_mask"]]) if np.any(masks["support_mask"]) else 0.0
                ),
                "dt_global_suggested": float(seed["dt_global_suggested"]),
            }
        )

    summary = {
        "params": {
            "ell": ell,
            "mp": mp,
            "rho_floor": rho_floor,
            "probe_dt": probe_dt,
            "grid_nx": int(len(x)),
            "grid_nz": int(len(z)),
            "times": times,
            "definition_support": "rho > 1e-3 * rho_max",
            "definition_high_density": "rho > 1e-1 * rho_max",
            "definition_near_transition": "| |ell^2 R_tilde| - 1 | <= 1",
            "definition_refine_mask": "support & (near_transition or residual>p99 or |grad phi|>p99 or derivative_norm>p99)",
        },
        "snapshots": snapshots,
        "files": {
            "maps": [str((out / f"d_scalar_tensor_reference_t{t:.1f}.png").resolve()) for t in times],
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(
            Path(__file__).resolve().parent.parent
            / "visualizations"
            / "d_scalar_tensor_reference_imex_prototype"
        ),
    )
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--probe-dt", type=float, default=1.0e-3)
    parser.add_argument("--nx", type=int, default=96)
    parser.add_argument("--nz", type=int, default=96)
    args = parser.parse_args()
    summary = run(
        outdir=args.outdir,
        ell=args.ell,
        mp=args.mp,
        rho_floor=args.rho_floor,
        probe_dt=args.probe_dt,
        nx=args.nx,
        nz=args.nz,
    )
    print(json.dumps(summary, indent=2))
