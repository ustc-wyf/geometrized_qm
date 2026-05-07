from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz, render_map


def spectral_dx(f: np.ndarray, dx: float) -> np.ndarray:
    nx, nz = f.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    fhat = np.fft.fft2(f)
    return np.fft.ifft2((1j * kx)[:, None] * fhat).real


def spectral_dz(f: np.ndarray, dz: float) -> np.ndarray:
    nx, nz = f.shape
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    fhat = np.fft.fft2(f)
    return np.fft.ifft2((1j * kz)[None, :] * fhat).real


def spectral_lap(f: np.ndarray, dx: float, dz: float) -> np.ndarray:
    nx, nz = f.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    fhat = np.fft.fft2(f)
    return np.fft.ifft2(-(KX * KX + KZ * KZ) * fhat).real


def make_sponge(x: np.ndarray, z: np.ndarray, inner_half_range: float, outer_half_range: float, sigma_max: float) -> np.ndarray:
    X, Z = np.meshgrid(x, z, indexing="ij")
    sx = np.zeros_like(X)
    sz = np.zeros_like(Z)
    w = max(outer_half_range - inner_half_range, 1e-6)
    mx = np.abs(X) > inner_half_range
    mz = np.abs(Z) > inner_half_range
    sx[mx] = ((np.abs(X[mx]) - inner_half_range) / w) ** 2
    sz[mz] = ((np.abs(Z[mz]) - inner_half_range) / w) ** 2
    return sigma_max * (sx + sz)


def analytic_slice(alpha: float, t: float, params: Exact2p1Params):
    x, z, psi, rho = integrate_xz(alpha, t, params)
    s = np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)
    return x, z, psi, rho, s


def boundary_exact(alpha: float, t: float, x: np.ndarray, z: np.ndarray, params: Exact2p1Params):
    nx = len(x)
    nz = len(z)
    psi = np.zeros((nx, nz), dtype=np.complex128)
    psi[:, 0] = evaluate_points(alpha, x, np.full_like(x, z[0]), t, params)
    psi[:, -1] = evaluate_points(alpha, x, np.full_like(x, z[-1]), t, params)
    psi[0, :] = evaluate_points(alpha, np.full_like(z, x[0]), z, t, params)
    psi[-1, :] = evaluate_points(alpha, np.full_like(z, x[-1]), z, t, params)
    rho = np.abs(psi) ** 2 + 1e-12
    s = np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)
    return rho, s


def apply_boundary(field: np.ndarray, target: np.ndarray) -> None:
    field[0, :] = target[0, :]
    field[-1, :] = target[-1, :]
    field[:, 0] = target[:, 0]
    field[:, -1] = target[:, -1]


def apply_bottom_boundary(field: np.ndarray, target: np.ndarray) -> None:
    field[:, 0] = target[:, 0]


def apply_geometry_boundary(tau: np.ndarray, beta: np.ndarray, phi: np.ndarray | None = None) -> None:
    for f in (tau, beta):
        f[0, :] = 0.0
        f[-1, :] = 0.0
        f[:, 0] = 0.0
        f[:, -1] = 0.0
    if phi is not None:
        phi[0, :] = 1.0
        phi[-1, :] = 1.0
        phi[:, 0] = 1.0
        phi[:, -1] = 1.0


def compute_matter_observables(rho: np.ndarray, s: np.ndarray, tau: np.ndarray, beta: np.ndarray, m: float, dx: float, dz: float):
    sx = spectral_dx(s, dx)
    sz = spectral_dz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)), 1e-12))
    w = np.exp(np.clip(beta + tau, -20.0, 20.0))
    n = np.maximum(rho * e * w, 1e-12)
    return sx, sz, e, n


def rhs_branch_b(tau, ptau, beta, pbeta, rho, s, m, mp, dx, dz):
    sx = spectral_dx(s, dx)
    sz = spectral_dz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)), 1e-12))
    n = np.maximum(rho * e * np.exp(np.clip(beta + tau, -20.0, 20.0)), 1e-12)
    mu = (m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)) * rho) / (mp * mp)

    taux = spectral_dx(tau, dx)
    tauz = spectral_dz(tau, dz)
    betax = spectral_dx(beta, dx)
    betaz = spectral_dz(beta, dz)

    dot_tau = ptau * ptau - taux * taux - tauz * tauz
    dot_beta = pbeta * pbeta - betax * betax - betaz * betaz
    dot_btau = pbeta * ptau - betax * taux - betaz * tauz

    tau_tt = spectral_lap(tau, dx, dz) - 0.5 * (dot_tau + dot_beta + mu)
    beta_tt = spectral_lap(beta, dx, dz) - dot_btau - dot_beta - 0.5 * mu

    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(spectral_dx(fluxx, dx) + spectral_dz(fluxz, dz))
    s_t = -e
    return tau_tt, beta_tt, n_t, s_t, e, n


def U_of_phi(phi: np.ndarray, ell: float) -> np.ndarray:
    phi = np.clip(phi, 1e-6, 1.0)
    return -(np.power(phi, 1.0 / 3.0) - phi) / (ell * ell)


def uprime_of_phi(phi: np.ndarray, ell: float) -> np.ndarray:
    phi = np.clip(phi, 1e-6, 1.0)
    return (1.0 - (1.0 / 3.0) * np.power(phi, -2.0 / 3.0)) / (ell * ell)


def rhs_branch_c(tau, ptau, beta, pbeta, phi, pphi, rho, s, m, mp, ell, dx, dz):
    sx = spectral_dx(s, dx)
    sz = spectral_dz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)), 1e-12))
    n = np.maximum(rho * e * np.exp(np.clip(beta + tau, -20.0, 20.0)), 1e-12)
    mu = (m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)) * rho) / (mp * mp)

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
    up = uprime_of_phi(phi, ell)
    exp2t = np.exp(np.clip(2.0 * tau, -20.0, 20.0))

    a = phi * dot_tau + dot_ptau + 0.5 * exp2t * u
    b = dot_tau + dot_btau + dot_beta + 0.5 * exp2t * up
    c = 2.0 * phi * dot_beta + 2.0 * phi * dot_btau + phi * dot_tau + 3.0 * dot_pbeta + 2.0 * dot_ptau + 1.5 * exp2t * u + 3.0 * mu

    box_tau = (c - 2.0 * a - 2.0 * phi * b) / (6.0 * np.clip(phi, 1e-6, None))
    box_beta = -b - 2.0 * box_tau
    box_phi = -a - 2.0 * phi * box_tau

    tau_tt = spectral_lap(tau, dx, dz) + box_tau
    beta_tt = spectral_lap(beta, dx, dz) + box_beta
    phi_tt = spectral_lap(phi, dx, dz) + box_phi

    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(spectral_dx(fluxx, dx) + spectral_dz(fluxz, dz))
    s_t = -e
    return tau_tt, beta_tt, phi_tt, n_t, s_t, e, n


def render_centerline(path: Path, x: np.ndarray, curves: list[tuple[str, np.ndarray]]) -> None:
    plt.figure(figsize=(6.2, 3.8))
    for label, y in curves:
        plt.plot(x, y, label=label, linewidth=1.6)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path, nx=241, nz=241, nkx=81, nkz=81, steps_hint=455, mp=300.0, ell=0.02):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    inner = 10.0
    outer = 20.0
    params = Exact2p1Params(x_half_range=outer, z_half_range=outer, nx=nx, nz=nz, nkx=nkx, nkz=nkz)
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    dt = 0.025
    t_final = params.overlap_time
    steps = steps_hint
    if steps <= 0:
        steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, _, rho_a, _ = analytic_slice(0.5, t_final, params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)

    # B branch
    _, _, _, rho_b, s_b = analytic_slice(0.5, 0.0, params)
    tau_b = np.zeros_like(rho_b)
    beta_b = np.zeros_like(rho_b)
    ptau_b = np.zeros_like(rho_b)
    pbeta_b = np.zeros_like(rho_b)

    # C branch
    _, _, _, rho_c, s_c = analytic_slice(0.5, 0.0, params)
    tau_c = np.zeros_like(rho_c)
    beta_c = np.zeros_like(rho_c)
    phi_c = np.ones_like(rho_c)
    ptau_c = np.zeros_like(rho_c)
    pbeta_c = np.zeros_like(rho_c)
    pphi_c = np.zeros_like(rho_c)

    for n in range(steps):
        t_next = (n + 1) * dt
        rho_bc, s_bc = boundary_exact(0.5, t_next, x, z, params)

        tau_tt, beta_tt, n_t, s_t, e_b, n_b = rhs_branch_b(
            tau_b, ptau_b, beta_b, pbeta_b, rho_b, s_b, params.m, mp, dx, dz
        )
        ptau_b = (ptau_b + dt * tau_tt) / (1.0 + 0.5 * sigma * dt)
        pbeta_b = (pbeta_b + dt * beta_tt) / (1.0 + 0.5 * sigma * dt)
        tau_b = tau_b + dt * ptau_b
        beta_b = beta_b + dt * pbeta_b
        s_b = s_b + dt * s_t
        n_b = (n_b + dt * n_t) / (1.0 + 0.5 * sigma * dt)
        apply_geometry_boundary(tau_b, beta_b)
        apply_bottom_boundary(rho_b, rho_bc)
        apply_bottom_boundary(s_b, s_bc)
        _, _, e_b2, n_b2 = compute_matter_observables(rho_b, s_b, tau_b, beta_b, params.m, dx, dz)
        apply_bottom_boundary(n_b, n_b2)
        rho_b = np.maximum(n_b / np.maximum(np.exp(np.clip(beta_b + tau_b, -20.0, 20.0)) * e_b2, 1e-12), 1e-12)
        apply_bottom_boundary(rho_b, rho_bc)

        tau_tt, beta_tt, phi_tt, n_t, s_t, e_c, n_c = rhs_branch_c(
            tau_c, ptau_c, beta_c, pbeta_c, phi_c, pphi_c, rho_c, s_c, params.m, mp, ell, dx, dz
        )
        ptau_c = (ptau_c + dt * tau_tt) / (1.0 + 0.5 * sigma * dt)
        pbeta_c = (pbeta_c + dt * beta_tt) / (1.0 + 0.5 * sigma * dt)
        pphi_c = (pphi_c + dt * phi_tt) / (1.0 + 0.5 * sigma * dt)
        tau_c = tau_c + dt * ptau_c
        beta_c = beta_c + dt * pbeta_c
        phi_c = np.clip(phi_c + dt * pphi_c, 0.2, 1.0)
        s_c = s_c + dt * s_t
        n_c = (n_c + dt * n_t) / (1.0 + 0.5 * sigma * dt)
        apply_geometry_boundary(tau_c, beta_c, phi_c)
        apply_bottom_boundary(rho_c, rho_bc)
        apply_bottom_boundary(s_c, s_bc)
        _, _, e_c2, n_c2 = compute_matter_observables(rho_c, s_c, tau_c, beta_c, params.m, dx, dz)
        apply_bottom_boundary(n_c, n_c2)
        rho_c = np.maximum(n_c / np.maximum(np.exp(np.clip(beta_c + tau_c, -20.0, 20.0)) * e_c2, 1e-12), 1e-12)
        apply_bottom_boundary(rho_c, rho_bc)

    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    xi = x[mask_x]
    zi = z[mask_z]
    A = rho_a[np.ix_(mask_x, mask_z)]
    B = rho_b[np.ix_(mask_x, mask_z)]
    C = rho_c[np.ix_(mask_x, mask_z)]

    render_map(out / "rho_A_inner.png", xi, zi, A, "A exact (inner)")
    render_map(out / "rho_B_inner.png", xi, zi, B, "B ps+fd (inner)")
    render_map(out / "rho_C_inner.png", xi, zi, C, "C ps+fd (inner)")
    render_map(out / "rho_B_minus_A.png", xi, zi, np.abs(B - A), "|B-A|")
    render_map(out / "rho_C_minus_A.png", xi, zi, np.abs(C - A), "|C-A|")
    render_centerline(out / "centerline_compare.png", xi, [("A", A[:, len(zi)//2]), ("B", B[:, len(zi)//2]), ("C", C[:, len(zi)//2])])

    summary = {
        "params": {
            "inner_half_range": inner,
            "outer_half_range": outer,
            "nx": params.nx,
            "nz": params.nz,
            "nkx": params.nkx,
            "nkz": params.nkz,
            "dt": dt,
            "steps": steps,
            "sigma_max": 2.0,
            "mp": mp,
            "ell": ell,
        },
        "errors_inner": {
            "B_vs_A_relL1": float(np.mean(np.abs(B - A)) / np.max(A)),
            "C_vs_A_relL1": float(np.mean(np.abs(C - A)) / np.max(A)),
            "C_vs_B_relL1": float(np.mean(np.abs(C - B)) / np.max(A)),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default=str(Path(__file__).resolve().parent.parent / "visualizations" / "bc_psfd_comparison"))
    parser.add_argument("--nx", type=int, default=241)
    parser.add_argument("--nz", type=int, default=241)
    parser.add_argument("--nkx", type=int, default=81)
    parser.add_argument("--nkz", type=int, default=81)
    parser.add_argument("--steps", dest="steps_hint", type=int, default=455)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--ell", type=float, default=0.02)
    args = parser.parse_args()
    summary = run(args.outdir, nx=args.nx, nz=args.nz, nkx=args.nkx, nkz=args.nkz, steps_hint=args.steps_hint, mp=args.mp, ell=args.ell)
    print(json.dumps(summary, indent=2))
