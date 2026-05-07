from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map
from simulate_bc_psfd_imex import (
    U_of_phi,
    edge_damp,
    make_sponge,
    spectral_dx,
    spectral_dz,
    spectral_k2,
    spectral_lap,
    wave_imex_step,
)


def analytic_slice(alpha: float, t: float, params: Exact2p1Params):
    x, z, psi, rho = integrate_xz(alpha, t, params)
    s = np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)
    return x, z, psi, rho, s


def compute_matter_observables(rho: np.ndarray, s: np.ndarray, tau: np.ndarray, beta: np.ndarray, m: float, dx: float, dz: float):
    sx = spectral_dx(s, dx)
    sz = spectral_dz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m * np.exp(np.clip(2.0 * tau, -12.0, 12.0)), 1e-12))
    w = np.exp(np.clip(beta + tau, -12.0, 12.0))
    n = np.maximum(rho * e * w, 1e-12)
    return sx, sz, e, n


def rhs_branch_b(tau, ptau, beta, pbeta, rho, s, m, mp, dx, dz):
    sx, sz, e, n = compute_matter_observables(rho, s, tau, beta, m, dx, dz)
    mu = (m * m * np.exp(np.clip(2.0 * tau, -12.0, 12.0)) * rho) / (mp * mp)
    taux = spectral_dx(tau, dx)
    tauz = spectral_dz(tau, dz)
    betax = spectral_dx(beta, dx)
    betaz = spectral_dz(beta, dz)
    dot_tau = ptau * ptau - taux * taux - tauz * tauz
    dot_beta = pbeta * pbeta - betax * betax - betaz * betaz
    dot_btau = pbeta * ptau - betax * taux - betaz * tauz
    src_tau = -0.5 * (dot_tau + dot_beta + mu)
    src_beta = -dot_btau - dot_beta - 0.5 * mu
    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(spectral_dx(fluxx, dx) + spectral_dz(fluxz, dz))
    s_t = -e
    return src_tau, src_beta, n_t, s_t, n


def rhs_branch_c(tau, ptau, beta, pbeta, phi, pphi, rho, s, m, mp, ell, dx, dz):
    sx, sz, e, n = compute_matter_observables(rho, s, tau, beta, m, dx, dz)
    mu = (m * m * np.exp(np.clip(2.0 * tau, -12.0, 12.0)) * rho) / (mp * mp)
    taux = spectral_dx(tau, dx)
    tauz = spectral_dz(tau, dz)
    betax = spectral_dx(beta, dx)
    betaz = spectral_dz(beta, dz)
    phix = spectral_dx(phi, dx)
    phiz = spectral_dz(phi, dz)
    dot_tau = ptau * ptau - taux * taux - tauz * tauz
    dot_beta = pbeta * pbeta - betax * betax - betaz * betaz
    dot_btau = pbeta * ptau - betax * taux - betaz * tauz
    dot_ptau = pphi * ptau - phix * taux - phiz * tauz
    dot_pbeta = pphi * pbeta - phix * betax - phiz * betaz
    u = U_of_phi(phi, ell)
    up = (1.0 - (1.0 / 3.0) * np.power(np.clip(phi, 1e-6, 1.0), -2.0 / 3.0)) / (ell * ell)
    exp2t = np.exp(np.clip(2.0 * tau, -12.0, 12.0))
    a = phi * dot_tau + dot_ptau + 0.5 * exp2t * u
    b = dot_tau + dot_btau + dot_beta + 0.5 * exp2t * up
    c = 2.0 * phi * dot_beta + 2.0 * phi * dot_btau + phi * dot_tau + 3.0 * dot_pbeta + 2.0 * dot_ptau + 1.5 * exp2t * u + 3.0 * mu
    src_tau = (c - 2.0 * a - 2.0 * phi * b) / (6.0 * np.clip(phi, 1e-6, None))
    src_beta = -b - 2.0 * src_tau
    src_phi = -a - 2.0 * phi * src_tau
    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(spectral_dx(fluxx, dx) + spectral_dz(fluxz, dz))
    s_t = -e
    return src_tau, src_beta, src_phi, n_t, s_t, n


def render_centerline(path: Path, x: np.ndarray, curves: list[tuple[str, np.ndarray]]):
    plt.figure(figsize=(6.2, 3.8))
    for label, y in curves:
        plt.plot(x, y, label=label, linewidth=1.6)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path, mp=300.0, ell=0.02):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    inner = 10.0
    outer = 20.0
    params = Exact2p1Params(x_half_range=outer, z_half_range=outer, nx=181, nz=181, nkx=61, nkz=61)
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    k2 = spectral_k2(params.nx, params.nz, dx, dz)
    t_final = params.overlap_time
    dt = 0.025
    steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, psi_prev, _ = integrate_xz(0.5, -dt, params)
    _, _, psi_a, rho0, s0 = analytic_slice(0.5, 0.0, params)
    _, _, psi_exact, rho_exact, _ = analytic_slice(0.5, t_final, params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)

    psi_now = psi_a.copy()

    rho_b = rho0.copy()
    s_b = s0.copy()
    tau_b = np.zeros_like(rho_b)
    beta_b = np.zeros_like(rho_b)
    ptau_b = np.zeros_like(rho_b)
    pbeta_b = np.zeros_like(rho_b)

    rho_c = rho0.copy()
    s_c = s0.copy()
    tau_c = np.zeros_like(rho_c)
    beta_c = np.zeros_like(rho_c)
    phi_c = np.ones_like(rho_c)
    ptau_c = np.zeros_like(rho_c)
    pbeta_c = np.zeros_like(rho_c)
    pphi_c = np.zeros_like(rho_c)

    for _ in range(steps):
        # A branch
        rhs_a = spectral_lap(psi_now, k2) - params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs_a
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        psi_prev, psi_now = psi_now, psi_next

        # B branch
        src_tau, src_beta, n_t, s_t, n_b = rhs_branch_b(tau_b, ptau_b, beta_b, pbeta_b, rho_b, s_b, params.m, mp, dx, dz)
        tau_b, ptau_b = wave_imex_step(tau_b, ptau_b, src_tau, dt, k2)
        beta_b, pbeta_b = wave_imex_step(beta_b, pbeta_b, src_beta, dt, k2)
        edge_damp(tau_b, sigma, dt, -2.0, 2.0)
        edge_damp(beta_b, sigma, dt, -2.0, 2.0)
        edge_damp(ptau_b, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_b, sigma, dt, -10.0, 10.0)
        s_b = s_b + dt * s_t
        n_b = n_b + dt * n_t
        _, _, e_b, n_b_geom = compute_matter_observables(rho_b, s_b, tau_b, beta_b, params.m, dx, dz)
        rho_b = np.maximum(n_b / np.maximum(np.exp(np.clip(beta_b + tau_b, -12.0, 12.0)) * e_b, 1e-12), 1e-12)

        # C branch
        src_tau, src_beta, src_phi, n_t, s_t, n_c = rhs_branch_c(tau_c, ptau_c, beta_c, pbeta_c, phi_c, pphi_c, rho_c, s_c, params.m, mp, ell, dx, dz)
        tau_c, ptau_c = wave_imex_step(tau_c, ptau_c, src_tau, dt, k2)
        beta_c, pbeta_c = wave_imex_step(beta_c, pbeta_c, src_beta, dt, k2)
        phi_c, pphi_c = wave_imex_step(phi_c, pphi_c, src_phi, dt, k2)
        edge_damp(tau_c, sigma, dt, -2.0, 2.0)
        edge_damp(beta_c, sigma, dt, -2.0, 2.0)
        edge_damp(phi_c, sigma, dt, 0.2, 1.0)
        edge_damp(ptau_c, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_c, sigma, dt, -10.0, 10.0)
        edge_damp(pphi_c, sigma, dt, -10.0, 10.0)
        s_c = s_c + dt * s_t
        n_c = n_c + dt * n_t
        _, _, e_c, n_c_geom = compute_matter_observables(rho_c, s_c, tau_c, beta_c, params.m, dx, dz)
        rho_c = np.maximum(n_c / np.maximum(np.exp(np.clip(beta_c + tau_c, -12.0, 12.0)) * e_c, 1e-12), 1e-12)

    rho_a_num = np.abs(psi_now) ** 2
    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    xi = x[mask_x]
    zi = z[mask_z]
    A = rho_exact[np.ix_(mask_x, mask_z)]
    A_num = rho_a_num[np.ix_(mask_x, mask_z)]
    B = rho_b[np.ix_(mask_x, mask_z)]
    C = rho_c[np.ix_(mask_x, mask_z)]

    render_map(out / "rho_A_exact.png", xi, zi, A, "A exact (inner)")
    render_map(out / "rho_A_num.png", xi, zi, A_num, "A cauchy (inner)")
    render_map(out / "rho_B_inner.png", xi, zi, B, "B IMEX cauchy (inner)")
    render_map(out / "rho_C_inner.png", xi, zi, C, "C IMEX cauchy (inner)")
    render_map(out / "rho_A_num_minus_A.png", xi, zi, np.abs(A_num - A), "|A_num-A|")
    render_map(out / "rho_B_minus_A.png", xi, zi, np.abs(B - A), "|B-A|")
    render_map(out / "rho_C_minus_A.png", xi, zi, np.abs(C - A), "|C-A|")
    render_centerline(out / "centerline_compare.png", xi, [("A", A[:, len(zi)//2]), ("A_num", A_num[:, len(zi)//2]), ("B", B[:, len(zi)//2]), ("C", C[:, len(zi)//2])])

    summary = {
        "params": {"mp": mp, "ell": ell, "dt": dt, "steps": steps},
        "errors_inner": {
            "A_num_vs_A_relL1": float(np.mean(np.abs(A_num - A)) / np.max(A)),
            "B_vs_A_relL1": float(np.mean(np.abs(B - A)) / np.max(A)),
            "C_vs_A_relL1": float(np.mean(np.abs(C - A)) / np.max(A)),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = run(Path(__file__).resolve().parent.parent / "visualizations" / "abc_cauchy_same_initial")
    print(json.dumps(summary, indent=2))
