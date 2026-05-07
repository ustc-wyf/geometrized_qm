from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from calibrate_exact_kg_ps_fd_absorbing import centerline_plot, make_sponge
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map
from simulate_bc_psfd_imex import (
    apply_geometry_boundary,
    edge_damp,
    rhs_branch_c,
    spectral_k2,
    wave_imex_step,
)


def build_exact_background_reduced(alpha: float):
    benchmark_params = Exact2p1Params(alpha=alpha)
    xb, zb, psi_bench, rho_bench = integrate_xz(alpha, benchmark_params.overlap_time, benchmark_params)
    s_bench = np.unwrap(np.unwrap(np.angle(psi_bench), axis=0), axis=1)

    solver_params = Exact2p1Params(
        alpha=alpha,
        x_half_range=15.0,
        z_half_range=15.0,
        nx=181,
        nz=181,
        nkx=61,
        nkz=61,
    )
    x, z, psi_outer, rho_outer = integrate_xz(alpha, solver_params.overlap_time, solver_params)
    s_outer = np.unwrap(np.unwrap(np.angle(psi_outer), axis=0), axis=1)
    return benchmark_params, solver_params, xb, zb, rho_bench, s_bench, x, z, rho_outer, s_outer


def c_step(
    tau: np.ndarray,
    ptau: np.ndarray,
    beta: np.ndarray,
    pbeta: np.ndarray,
    eta: np.ndarray,
    peta: np.ndarray,
    rho_a_outer: np.ndarray,
    s_a_outer: np.ndarray,
    solver_params,
    mp: float,
    ell: float,
    dx: float,
    dz: float,
    k2: np.ndarray,
    sigma: np.ndarray,
    lambda_grav: float,
    dt: float,
):
    phi = 1.0 + eta
    src_tau, src_beta, src_phi, *_ = rhs_branch_c(
        tau, ptau, beta, pbeta, phi, peta, rho_a_outer, s_a_outer, solver_params.m, mp, ell, dx, dz
    )
    tau1, ptau1 = wave_imex_step(tau, ptau, lambda_grav * src_tau, dt, k2)
    beta1, pbeta1 = wave_imex_step(beta, pbeta, lambda_grav * src_beta, dt, k2)
    eta1, peta1 = wave_imex_step(eta, peta, lambda_grav * src_phi, dt, k2)

    edge_damp(tau1, sigma, dt, -2.0, 2.0)
    edge_damp(beta1, sigma, dt, -2.0, 2.0)
    edge_damp(eta1, sigma, dt, -0.8, 0.0)
    edge_damp(ptau1, sigma, dt, -10.0, 10.0)
    edge_damp(pbeta1, sigma, dt, -10.0, 10.0)
    edge_damp(peta1, sigma, dt, -10.0, 10.0)

    phi1 = np.clip(1.0 + eta1, 0.2, 1.0)
    eta1 = phi1 - 1.0
    apply_geometry_boundary(tau1, beta1, phi1)
    eta1[0, :] = 0.0
    eta1[-1, :] = 0.0
    eta1[:, 0] = 0.0
    eta1[:, -1] = 0.0

    return tau1, ptau1, beta1, pbeta1, eta1, peta1


def error_norm(full, refined, atol: float, rtol: float) -> float:
    errs = []
    for a, b in zip(full, refined):
        scale = atol + rtol * np.maximum(np.abs(a), np.abs(b))
        errs.append(np.max(np.abs(a - b) / scale))
    return float(max(errs))


def run(
    outdir: str | Path,
    alpha: float = 0.5,
    lambda_grav: float = 10.0,
    mp: float = 300.0,
    ell: float = 0.02,
    dt_initial: float = 0.025,
    rtol: float = 1e-3,
    atol: float = 1e-8,
):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    benchmark_params, solver_params, xb, zb, rho_a_exact, s_a_exact, x, z, rho_a_outer, s_a_outer = build_exact_background_reduced(alpha)
    total_time = solver_params.overlap_time
    dx = 2.0 * solver_params.x_half_range / (solver_params.nx - 1)
    dz = 2.0 * solver_params.z_half_range / (solver_params.nz - 1)
    k2 = spectral_k2(solver_params.nx, solver_params.nz, dx, dz)
    sigma = make_sponge(x, z, inner_half_range=12.0, outer_half_range=solver_params.x_half_range, sigma_max=2.0)

    tau = np.zeros_like(rho_a_outer)
    ptau = np.zeros_like(rho_a_outer)
    beta = np.zeros_like(rho_a_outer)
    pbeta = np.zeros_like(rho_a_outer)
    eta = np.zeros_like(rho_a_outer)
    peta = np.zeros_like(rho_a_outer)

    t = 0.0
    dt = dt_initial
    accepted = 0
    rejected = 0
    dt_min = dt
    dt_max = dt

    while t < total_time:
        dt = min(dt, total_time - t)
        state = (tau, ptau, beta, pbeta, eta, peta)

        full = c_step(*state, rho_a_outer, s_a_outer, solver_params, mp, ell, dx, dz, k2, sigma, lambda_grav, dt)
        half = c_step(*state, rho_a_outer, s_a_outer, solver_params, mp, ell, dx, dz, k2, sigma, lambda_grav, dt / 2.0)
        half2 = c_step(*half, rho_a_outer, s_a_outer, solver_params, mp, ell, dx, dz, k2, sigma, lambda_grav, dt / 2.0)

        err = error_norm(full, half2, atol=atol, rtol=rtol)
        if err <= 1.0:
            tau, ptau, beta, pbeta, eta, peta = [arr.copy() for arr in half2]
            t += dt
            accepted += 1
            dt_min = min(dt_min, dt)
            dt_max = max(dt_max, dt)
            if err == 0.0:
                growth = 2.0
            else:
                growth = min(2.0, max(1.1, 0.9 * (1.0 / err) ** 0.5))
            dt *= growth
        else:
            rejected += 1
            dt *= max(0.2, 0.9 * (1.0 / err) ** 0.5)

    phi = 1.0 + eta
    rho_c = rho_a_outer * np.exp(-np.clip(beta + tau, -12.0, 12.0))

    ix = np.where(np.abs(x) <= benchmark_params.x_half_range + 1e-12)[0]
    iz = np.where(np.abs(z) <= benchmark_params.z_half_range + 1e-12)[0]
    xi = x[ix]
    zi = z[iz]
    if not (np.allclose(xi, xb) and np.allclose(zi, zb)):
        raise RuntimeError("Cropped outer grid does not match approved benchmark grid.")

    A = rho_a_exact
    C = rho_c[np.ix_(ix, iz)]

    render_map(out / "rho_A_inner.png", xi, zi, A, "A exact benchmark")
    render_map(out / "rho_C_inner.png", xi, zi, C, "C adaptive debug")
    render_map(out / "rho_C_minus_A.png", xi, zi, np.abs(C - A), "|C-A|")
    centerline_plot(out / "centerline_C.png", xi, A[:, len(zi) // 2], C[:, len(zi) // 2])

    summary = {
        "alpha": alpha,
        "lambda_grav": lambda_grav,
        "mp": mp,
        "ell": ell,
        "rtol": rtol,
        "atol": atol,
        "accepted_steps": accepted,
        "rejected_steps": rejected,
        "dt_min": dt_min,
        "dt_max": dt_max,
        "C_vs_A_relL1": float(np.mean(np.abs(C - A)) / np.max(A)),
        "C_tau_max": float(np.max(np.abs(tau))),
        "C_phi_dev_max": float(np.max(np.abs(eta))),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--lambda-grav", type=float, default=10.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--ell", type=float, default=0.02)
    parser.add_argument("--dt-initial", type=float, default=0.025)
    parser.add_argument("--rtol", type=float, default=1e-2)
    parser.add_argument("--atol", type=float, default=1e-6)
    parser.add_argument("--outdir", type=str, default=str(Path(__file__).resolve().parent.parent / "visualizations" / "c_adaptive_backbone"))
    args = parser.parse_args()
    run(
        args.outdir,
        alpha=args.alpha,
        lambda_grav=args.lambda_grav,
        mp=args.mp,
        ell=args.ell,
        dt_initial=args.dt_initial,
        rtol=args.rtol,
        atol=args.atol,
    )
