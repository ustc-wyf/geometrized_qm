from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz, render_map
from simulate_bc_psfd_imex import (
    U_of_phi,
    uprime_of_phi,
    apply_bottom_boundary,
    apply_geometry_boundary,
    boundary_exact,
    compute_matter_observables,
    edge_damp,
    make_sponge,
    rhs_branch_c,
    spectral_k2,
    wave_imex_step,
)


def render_line(path: Path, xs, ys, xlabel, ylabel, title):
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(xs, ys, marker="o")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def render_map_fixed_scale(path: Path, x, z, field, title, vmax):
    plt.figure(figsize=(5.6, 4.8))
    plt.imshow(
        field.T,
        origin="lower",
        extent=[x.min(), x.max(), z.min(), z.max()],
        aspect="equal",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    plt.colorbar(label="rho")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run_case(ell: float, mp: float, outdir: Path):
    inner = 10.0
    outer = 20.0
    params = Exact2p1Params(x_half_range=outer, z_half_range=outer, nx=121, nz=121, nkx=41, nkz=41)
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    k2 = spectral_k2(params.nx, params.nz, dx, dz)
    t_final = params.overlap_time
    steps = 220
    dt = t_final / steps

    x, z, _, rho_a = integrate_xz(0.5, t_final, params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)

    _, _, psi0, rho_c = integrate_xz(0.5, 0.0, params)
    s_c = np.unwrap(np.unwrap(np.angle(psi0), axis=0), axis=1)
    tau_c = np.zeros_like(rho_c)
    beta_c = np.zeros_like(rho_c)
    phi_c = np.ones_like(rho_c)
    ptau_c = np.zeros_like(rho_c)
    pbeta_c = np.zeros_like(rho_c)
    pphi_c = np.zeros_like(rho_c)

    for n in range(steps):
        t_next = (n + 1) * dt
        rho_bc, s_bc = boundary_exact(0.5, t_next, x, z, params)
        src_tau, src_beta, src_phi, n_t, s_t, e_c, n_c = rhs_branch_c(
            tau_c, ptau_c, beta_c, pbeta_c, phi_c, pphi_c, rho_c, s_c, params.m, mp, ell, dx, dz
        )
        tau_c, ptau_c = wave_imex_step(tau_c, ptau_c, src_tau, dt, k2)
        beta_c, pbeta_c = wave_imex_step(beta_c, pbeta_c, src_beta, dt, k2)
        phi_c, pphi_c = wave_imex_step(phi_c, pphi_c, src_phi, dt, k2)
        edge_damp(tau_c, sigma, dt, -2.0, 2.0)
        edge_damp(beta_c, sigma, dt, -2.0, 2.0)
        edge_damp(phi_c, sigma, dt, 0.2, 1.0)
        edge_damp(ptau_c, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_c, sigma, dt, -10.0, 10.0)
        edge_damp(pphi_c, sigma, dt, -10.0, 10.0)
        apply_geometry_boundary(tau_c, beta_c, phi_c)
        s_c = s_c + dt * s_t
        apply_bottom_boundary(s_c, s_bc)
        n_c = n_c + dt * n_t
        _, _, e_c2, n_c2 = compute_matter_observables(rho_c, s_c, tau_c, beta_c, params.m, dx, dz)
        apply_bottom_boundary(n_c, n_c2)
        rho_c = np.maximum(n_c / np.maximum(np.exp(np.clip(beta_c + tau_c, -12.0, 12.0)) * e_c2, 1e-12), 1e-12)
        apply_bottom_boundary(rho_c, rho_bc)

    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    xi = x[mask_x]
    zi = z[mask_z]
    A = rho_a[np.ix_(mask_x, mask_z)]
    C = rho_c[np.ix_(mask_x, mask_z)]
    Phi = phi_c[np.ix_(mask_x, mask_z)]

    tag = f"{ell:.0e}".replace("+0", "").replace("-0", "-")
    render_map(outdir / f"rho_C_ell_{tag}.png", xi, zi, C, f"C rho, ell={ell:g}")
    render_map_fixed_scale(
        outdir / f"rho_C_ell_{tag}_matched_to_A.png",
        xi,
        zi,
        C,
        f"C rho matched scale, ell={ell:g}",
        vmax=float(np.max(A)),
    )
    render_map(outdir / f"rho_abs_diff_ell_{tag}.png", xi, zi, np.abs(C - A), f"|C-A|, ell={ell:g}")

    # low-curvature proxy: phi close to 1
    low_mask = np.abs(Phi - 1.0) < 1e-2
    if np.any(low_mask):
        low_rel_l1 = float(np.mean(np.abs(C[low_mask] - A[low_mask])) / np.max(A))
    else:
        low_rel_l1 = float("nan")
    full_rel_l1 = float(np.mean(np.abs(C - A)) / np.max(A))
    return {
        "ell": ell,
        "mp": mp,
        "full_rel_l1": full_rel_l1,
        "low_rel_l1": low_rel_l1,
        "phi_dev_max": float(np.max(np.abs(Phi - 1.0))),
        "low_mask_fraction": float(np.mean(low_mask)),
        "rho_c_max": float(np.max(C)),
    }


def main():
    outdir = Path(__file__).resolve().parent.parent / "visualizations" / "c_small_ell_local"
    outdir.mkdir(parents=True, exist_ok=True)
    mp = 300.0
    ells = [2e-2, 1e-2, 5e-3, 2e-3, 1e-3]
    summaries = [run_case(ell, mp, outdir) for ell in ells]
    (outdir / "summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")

    render_line(outdir / "full_rel_l1_vs_ell.png", ells, [s["full_rel_l1"] for s in summaries], "ell", "rel L1", "C vs A full-domain")
    render_line(outdir / "low_rel_l1_vs_ell.png", ells, [max(s["low_rel_l1"], 1e-16) for s in summaries], "ell", "rel L1", "C vs A low-curvature proxy")


if __name__ == "__main__":
    main()
