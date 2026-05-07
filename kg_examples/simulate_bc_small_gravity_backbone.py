from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from calibrate_exact_kg_ps_fd_absorbing import centerline_plot, make_sponge
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map
from simulate_bc_psfd_imex import (
    apply_geometry_boundary,
    compute_matter_observables,
    edge_damp,
    rhs_branch_b,
    rhs_branch_c,
    spectral_k2,
    wave_imex_step,
)


def edge_damp_to_anchor(field: np.ndarray, sigma: np.ndarray, dt: float, anchor: float, lo=None, hi=None):
    field[:] = anchor + (field - anchor) / (1.0 + sigma * dt)
    if lo is not None or hi is not None:
        np.clip(field, lo, hi, out=field)


def build_exact_background(alpha: float):
    # Approved benchmark for the displayed A image.
    benchmark_params = Exact2p1Params(alpha=alpha)
    xb, zb, psi_bench, rho_bench = integrate_xz(alpha, benchmark_params.overlap_time, benchmark_params)
    s_bench = np.unwrap(np.unwrap(np.angle(psi_bench), axis=0), axis=1)

    # Larger exact box for driving B/C geometry without the absorber touching the support.
    solver_params = Exact2p1Params(
        alpha=alpha,
        x_half_range=20.0,
        z_half_range=20.0,
        nx=241,
        nz=241,
        nkx=81,
        nkz=81,
    )
    x, z, psi_outer, rho_outer = integrate_xz(alpha, solver_params.overlap_time, solver_params)
    s_outer = np.unwrap(np.unwrap(np.angle(psi_outer), axis=0), axis=1)
    return benchmark_params, solver_params, xb, zb, rho_bench, s_bench, x, z, rho_outer, s_outer


def run(outdir: str | Path, alpha: float = 0.5, lambda_grav: float = 1.0e-2, mp: float = 300.0, ell: float = 0.02, dt_target: float = 0.025, c_substeps: int = 1):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    benchmark_params, solver_params, xb, zb, rho_a_exact, s_a_exact, x, z, rho_a_outer, s_a_outer = build_exact_background(alpha)
    dt = dt_target
    steps = int(round(solver_params.overlap_time / dt))
    dt = solver_params.overlap_time / steps
    dx = 2.0 * solver_params.x_half_range / (solver_params.nx - 1)
    dz = 2.0 * solver_params.z_half_range / (solver_params.nz - 1)
    k2 = spectral_k2(solver_params.nx, solver_params.nz, dx, dz)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=solver_params.x_half_range, sigma_max=2.0)

    tau_b = np.zeros_like(rho_a_outer)
    beta_b = np.zeros_like(rho_a_outer)
    ptau_b = np.zeros_like(rho_a_outer)
    pbeta_b = np.zeros_like(rho_a_outer)

    tau_c = np.zeros_like(rho_a_outer)
    beta_c = np.zeros_like(rho_a_outer)
    eta_c = np.zeros_like(rho_a_outer)
    ptau_c = np.zeros_like(rho_a_outer)
    pbeta_c = np.zeros_like(rho_a_outer)
    peta_c = np.zeros_like(rho_a_outer)

    for _ in range(steps):
        src_tau_b, src_beta_b, *_ = rhs_branch_b(tau_b, ptau_b, beta_b, pbeta_b, rho_a_outer, s_a_outer, solver_params.m, mp, dx, dz, k2)
        tau_b, ptau_b = wave_imex_step(tau_b, ptau_b, lambda_grav * src_tau_b, dt, k2)
        beta_b, pbeta_b = wave_imex_step(beta_b, pbeta_b, lambda_grav * src_beta_b, dt, k2)
        edge_damp(tau_b, sigma, dt, -2.0, 2.0)
        edge_damp(beta_b, sigma, dt, -2.0, 2.0)
        edge_damp(ptau_b, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_b, sigma, dt, -10.0, 10.0)
        apply_geometry_boundary(tau_b, beta_b)

        dt_c = dt / max(c_substeps, 1)
        sigma_c = sigma / max(c_substeps, 1)
        for _sub in range(max(c_substeps, 1)):
            phi_c = 1.0 + eta_c
            src_tau_c, src_beta_c, src_phi_c, *_ = rhs_branch_c(
                tau_c, ptau_c, beta_c, pbeta_c, phi_c, peta_c, rho_a_outer, s_a_outer, solver_params.m, mp, ell, dx, dz
            )
            tau_c, ptau_c = wave_imex_step(tau_c, ptau_c, lambda_grav * src_tau_c, dt_c, k2)
            beta_c, pbeta_c = wave_imex_step(beta_c, pbeta_c, lambda_grav * src_beta_c, dt_c, k2)
            eta_c, peta_c = wave_imex_step(eta_c, peta_c, lambda_grav * src_phi_c, dt_c, k2)
            edge_damp(tau_c, sigma_c, dt_c, -2.0, 2.0)
            edge_damp(beta_c, sigma_c, dt_c, -2.0, 2.0)
            edge_damp(eta_c, sigma_c, dt_c, -0.8, 0.0)
            edge_damp(ptau_c, sigma_c, dt_c, -10.0, 10.0)
            edge_damp(pbeta_c, sigma_c, dt_c, -10.0, 10.0)
            edge_damp(peta_c, sigma_c, dt_c, -10.0, 10.0)
            phi_c = np.clip(1.0 + eta_c, 0.2, 1.0)
            eta_c = phi_c - 1.0
            apply_geometry_boundary(tau_c, beta_c, phi_c)
            eta_c[0, :] = 0.0
            eta_c[-1, :] = 0.0
            eta_c[:, 0] = 0.0
            eta_c[:, -1] = 0.0

    # Simple constraint-compatible density map for debugging:
    # when lambda_grav=0, tau=beta=0 and these reduce exactly to rho_A.
    rho_b = rho_a_outer * np.exp(-np.clip(beta_b + tau_b, -12.0, 12.0))
    rho_c = rho_a_outer * np.exp(-np.clip(beta_c + tau_c, -12.0, 12.0))

    # Crop the outer-grid B/C fields onto the exact benchmark inner grid.
    ix = np.where(np.abs(x) <= benchmark_params.x_half_range + 1e-12)[0]
    iz = np.where(np.abs(z) <= benchmark_params.z_half_range + 1e-12)[0]
    xi = x[ix]
    zi = z[iz]
    if not (np.allclose(xi, xb) and np.allclose(zi, zb)):
        raise RuntimeError("Cropped outer grid does not match approved benchmark grid.")
    A = rho_a_exact
    B = rho_b[np.ix_(ix, iz)]
    C = rho_c[np.ix_(ix, iz)]

    render_map(out / "rho_A_inner.png", xi, zi, A, "A exact benchmark")
    render_map(out / "rho_B_inner.png", xi, zi, B, "B small-gravity debug")
    render_map(out / "rho_C_inner.png", xi, zi, C, "C small-gravity debug")
    render_map(out / "rho_B_minus_A.png", xi, zi, np.abs(B - A), "|B-A|")
    render_map(out / "rho_C_minus_A.png", xi, zi, np.abs(C - A), "|C-A|")
    centerline_plot(out / "centerline_B.png", xi, A[:, len(zi) // 2], B[:, len(zi) // 2])
    centerline_plot(out / "centerline_C.png", xi, A[:, len(zi) // 2], C[:, len(zi) // 2])

    summary = {
        "alpha": alpha,
        "lambda_grav": lambda_grav,
        "mp": mp,
        "ell": ell,
        "dt": dt,
        "steps": steps,
        "c_substeps": c_substeps,
        "B_vs_A_relL1": float(np.mean(np.abs(B - A)) / np.max(A)),
        "C_vs_A_relL1": float(np.mean(np.abs(C - A)) / np.max(A)),
        "C_vs_B_relL1": float(np.mean(np.abs(C - B)) / np.max(A)),
        "B_tau_max": float(np.max(np.abs(tau_b))),
        "C_tau_max": float(np.max(np.abs(tau_c))),
        "C_phi_dev_max": float(np.max(np.abs(eta_c))),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--lambda-grav", type=float, default=1.0e-2)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--ell", type=float, default=0.02)
    parser.add_argument("--dt-target", type=float, default=0.025)
    parser.add_argument("--c-substeps", type=int, default=1)
    parser.add_argument("--outdir", type=str, default=str(Path(__file__).resolve().parent.parent / "visualizations" / "bc_small_gravity_backbone"))
    args = parser.parse_args()
    run(args.outdir, alpha=args.alpha, lambda_grav=args.lambda_grav, mp=args.mp, ell=args.ell, dt_target=args.dt_target, c_substeps=args.c_substeps)
