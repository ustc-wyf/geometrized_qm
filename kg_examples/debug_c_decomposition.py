from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map
from simulate_bc_psfd_imex import (
    U_of_phi,
    apply_geometry_boundary,
    boundary_exact,
    compute_matter_observables,
    edge_damp,
    make_sponge,
    matter_rhs_conservative,
    rhs_branch_b,
    rhs_branch_c,
    spectral_k2,
    wave_imex_step,
)


def render_centerline(path: Path, x: np.ndarray, curves):
    plt.figure(figsize=(6.2, 3.8))
    for label, y in curves:
        plt.plot(x, y, label=label, linewidth=1.5)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run_case(kind: str, x, z, rho0, s0, params, steps, dt, sigma, mp, ell):
    rho = rho0.copy()
    s = s0.copy()
    tau = np.zeros_like(rho)
    beta = np.zeros_like(rho)
    phi = np.ones_like(rho)
    ptau = np.zeros_like(rho)
    pbeta = np.zeros_like(rho)
    pphi = np.zeros_like(rho)
    dx = 2.0 * params.x_half_range / (params.nx - 1)
    dz = 2.0 * params.z_half_range / (params.nz - 1)
    k2 = spectral_k2(params.nx, params.nz, dx, dz)

    for n in range(steps):
        t_next = (n + 1) * dt
        rho_bc, s_bc = boundary_exact(0.5, t_next, x, z, params)

        if kind == "Bflat":
            src_tau = np.zeros_like(rho)
            src_beta = np.zeros_like(rho)
            sx, sz, e, ncur = compute_matter_observables(rho, s, tau, beta, params.m, dx, dz)
            n_t = matter_rhs_conservative(ncur, sx, sz, e, dx, dz)
            s_t = -e
        elif kind == "Bfull":
            src_tau, src_beta, n_t, s_t, e, ncur = rhs_branch_b(tau, ptau, beta, pbeta, rho, s, params.m, mp, dx, dz, k2)
        elif kind == "Cphi1":
            # Keep phi = 1, so this should collapse to Bfull numerically
            phi[:] = 1.0
            pphi[:] = 0.0
            src_tau, src_beta, src_phi, n_t, s_t, e, ncur = rhs_branch_c(
                tau, ptau, beta, pbeta, phi, pphi, rho, s, params.m, mp, ell, dx, dz
            )
        elif kind == "Cfull":
            src_tau, src_beta, src_phi, n_t, s_t, e, ncur = rhs_branch_c(
                tau, ptau, beta, pbeta, phi, pphi, rho, s, params.m, mp, ell, dx, dz
            )
        else:
            raise ValueError(kind)

        if kind in ("Bfull", "Cphi1", "Cfull"):
            tau, ptau = wave_imex_step(tau, ptau, src_tau, dt, k2)
            beta, pbeta = wave_imex_step(beta, pbeta, src_beta, dt, k2)
            edge_damp(tau, sigma, dt, -2.0, 2.0)
            edge_damp(beta, sigma, dt, -2.0, 2.0)
            edge_damp(ptau, sigma, dt, -10.0, 10.0)
            edge_damp(pbeta, sigma, dt, -10.0, 10.0)
            apply_geometry_boundary(tau, beta)

        if kind == "Cfull":
            phi, pphi = wave_imex_step(phi, pphi, src_phi, dt, k2)
            edge_damp(phi, sigma, dt, 0.2, 1.0)
            edge_damp(pphi, sigma, dt, -10.0, 10.0)
            apply_geometry_boundary(tau, beta, phi)

        s = s + dt * s_t
        from simulate_bc_psfd_imex import apply_bottom_boundary
        apply_bottom_boundary(s, s_bc)
        ncur = ncur + dt * n_t
        _, _, e2, n2 = compute_matter_observables(rho, s, tau, beta, params.m, dx, dz)
        apply_bottom_boundary(ncur, n2)
        rho = np.maximum(ncur / np.maximum(np.exp(np.clip(beta + tau, -12.0, 12.0)) * e2, 1e-12), 1e-12)
        apply_bottom_boundary(rho, rho_bc)

    return rho, tau, beta, phi


def main():
    out = Path(__file__).resolve().parent.parent / "visualizations" / "debug_c_decomposition"
    out.mkdir(parents=True, exist_ok=True)
    params = Exact2p1Params(x_half_range=20.0, z_half_range=20.0, nx=121, nz=121, nkx=41, nkz=41)
    t_final = params.overlap_time
    steps = 220
    dt = t_final / steps
    x, z, psi0, rho0 = integrate_xz(0.5, 0.0, params)
    s0 = np.unwrap(np.unwrap(np.angle(psi0), axis=0), axis=1)
    _, _, _, rhoA = integrate_xz(0.5, t_final, params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=20.0, sigma_max=2.0)

    results = {}
    rho_fields = {}
    for kind in ("Bflat", "Bfull", "Cphi1", "Cfull"):
        rho, tau, beta, phi = run_case(kind, x, z, rho0, s0, params, steps, dt, sigma, mp=300.0, ell=0.02)
        rho_fields[kind] = rho
        results[kind] = {
            "rho_max": float(np.nanmax(rho)),
            "tau_max": float(np.nanmax(np.abs(tau))),
            "beta_max": float(np.nanmax(np.abs(beta))),
            "phi_dev_max": float(np.nanmax(np.abs(phi - 1.0))),
        }
        render_map(out / f"{kind}.png", x, z, rho, kind)

    center = len(z) // 2
    render_centerline(
        out / "centerline_compare.png",
        x,
        [
            ("A", rhoA[:, center]),
            ("Bflat", rho_fields["Bflat"][:, center]),
            ("Bfull", rho_fields["Bfull"][:, center]),
            ("Cphi1", rho_fields["Cphi1"][:, center]),
            ("Cfull", rho_fields["Cfull"][:, center]),
        ],
    )
    (out / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
