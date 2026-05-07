from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map
from simulate_bc_psfd_imex import (
    apply_geometry_boundary,
    compute_matter_observables,
    edge_damp,
    make_sponge,
    matter_rhs_conservative,
    rhs_branch_b,
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


def apply_neumann_edges(field: np.ndarray) -> None:
    field[0, :] = field[1, :]
    field[-1, :] = field[-2, :]
    field[:, 0] = field[:, 1]
    field[:, -1] = field[:, -2]


def run_case(kind: str, alpha: float, steps: int):
    params = Exact2p1Params(alpha=alpha, x_half_range=20.0, z_half_range=20.0, nx=121, nz=121, nkx=41, nkz=41)
    t_final = params.overlap_time
    dt = t_final / steps
    x, z, psi0, rho0 = integrate_xz(alpha, 0.0, params)
    _, _, _, rho_exact = integrate_xz(alpha, t_final, params)
    s0 = np.unwrap(np.unwrap(np.angle(psi0), axis=0), axis=1)

    dx = 2.0 * params.x_half_range / (params.nx - 1)
    dz = 2.0 * params.z_half_range / (params.nz - 1)
    k2 = spectral_k2(params.nx, params.nz, dx, dz)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=20.0, sigma_max=2.0)

    rho = rho0.copy()
    s = s0.copy()
    tau = np.zeros_like(rho)
    beta = np.zeros_like(rho)
    ptau = np.zeros_like(rho)
    pbeta = np.zeros_like(rho)

    for _ in range(steps):
        if kind == "Bflat":
            sx, sz, e, n = compute_matter_observables(rho, s, tau, beta, params.m, dx, dz)
            n_t = matter_rhs_conservative(n, sx, sz, e, dx, dz)
            s_t = -e
        elif kind == "Bfull":
            src_tau, src_beta, n_t, s_t, e, n = rhs_branch_b(tau, ptau, beta, pbeta, rho, s, params.m, 300.0, dx, dz, k2)
            tau, ptau = wave_imex_step(tau, ptau, src_tau, dt, k2)
            beta, pbeta = wave_imex_step(beta, pbeta, src_beta, dt, k2)
            edge_damp(tau, sigma, dt, -2.0, 2.0)
            edge_damp(beta, sigma, dt, -2.0, 2.0)
            edge_damp(ptau, sigma, dt, -10.0, 10.0)
            edge_damp(pbeta, sigma, dt, -10.0, 10.0)
            apply_geometry_boundary(tau, beta)
        else:
            raise ValueError(kind)

        s = s + dt * s_t
        apply_neumann_edges(s)
        n = n + dt * n_t
        apply_neumann_edges(n)
        _, _, e2, _ = compute_matter_observables(rho, s, tau, beta, params.m, dx, dz)
        rho = np.maximum(n / np.maximum(np.exp(np.clip(beta + tau, -12.0, 12.0)) * e2, 1e-12), 1e-12)
        apply_neumann_edges(rho)

    return {
        "x": x,
        "z": z,
        "rho": rho,
        "rho_exact": rho_exact,
        "tau": tau,
        "beta": beta,
        "t_final": t_final,
        "dt": dt,
    }


def main():
    out = Path(__file__).resolve().parent.parent / "visualizations" / "debug_b_single_beam"
    out.mkdir(parents=True, exist_ok=True)

    alpha = 1.0
    a_exact = run_case("Bflat", alpha, 1)  # cheap access to exact grids only
    x = a_exact["x"]
    z = a_exact["z"]
    rho_exact = a_exact["rho_exact"]

    render_map(out / "A_exact_alpha_1.png", x, z, rho_exact, "A exact alpha=1")

    results = {}
    for kind in ("Bflat", "Bfull"):
        coarse = run_case(kind, alpha, 220)
        fine = run_case(kind, alpha, 440)
        render_map(out / f"{kind}_alpha_1_steps220.png", x, z, coarse["rho"], f"{kind} alpha=1 steps=220")
        render_map(out / f"{kind}_alpha_1_steps440.png", x, z, fine["rho"], f"{kind} alpha=1 steps=440")
        render_map(out / f"{kind}_minus_A_steps220.png", x, z, np.abs(coarse["rho"] - rho_exact), f"|{kind}-A| steps=220")
        render_map(out / f"{kind}_minus_A_steps440.png", x, z, np.abs(fine["rho"] - rho_exact), f"|{kind}-A| steps=440")
        render_centerline(
            out / f"{kind}_centerline_compare.png",
            x,
            [
                ("A exact", rho_exact[:, len(z) // 2]),
                (f"{kind} 220", coarse["rho"][:, len(z) // 2]),
                (f"{kind} 440", fine["rho"][:, len(z) // 2]),
            ],
        )
        rel_between = float(np.mean(np.abs(fine["rho"] - coarse["rho"])) / max(np.max(rho_exact), 1e-12))
        rel_to_exact_220 = float(np.mean(np.abs(coarse["rho"] - rho_exact)) / max(np.max(rho_exact), 1e-12))
        rel_to_exact_440 = float(np.mean(np.abs(fine["rho"] - rho_exact)) / max(np.max(rho_exact), 1e-12))
        results[kind] = {
            "rho_max_220": float(np.max(coarse["rho"])),
            "rho_max_440": float(np.max(fine["rho"])),
            "tau_max_220": float(np.max(np.abs(coarse["tau"]))),
            "tau_max_440": float(np.max(np.abs(fine["tau"]))),
            "beta_max_220": float(np.max(np.abs(coarse["beta"]))),
            "beta_max_440": float(np.max(np.abs(fine["beta"]))),
            "rel_between_220_440": rel_between,
            "rel_to_exact_220": rel_to_exact_220,
            "rel_to_exact_440": rel_to_exact_440,
        }

    (out / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
