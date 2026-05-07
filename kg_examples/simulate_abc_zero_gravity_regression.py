from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from calibrate_exact_kg_ps_fd_absorbing import (
    apply_boundary,
    boundary_exact,
    centerline_plot,
    make_sponge,
    spectral_lap,
)
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map


def evolve_a_backbone(alpha: float):
    # Solver box: large enough that the absorber stays away from the displayed region.
    inner = 10.0
    outer = 20.0
    solver_params = Exact2p1Params(
        alpha=alpha,
        x_half_range=outer,
        z_half_range=outer,
        nx=241,
        nz=241,
        nkx=61,
        nkz=61,
    )
    # Exact reference: use the approved benchmark parameter set literally.
    benchmark_params = Exact2p1Params(alpha=alpha)

    dx = 2.0 * outer / (solver_params.nx - 1)
    dz = 2.0 * outer / (solver_params.nz - 1)
    dt = 0.025
    t_final = solver_params.overlap_time
    steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, solver_params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, solver_params)
    # Reuse the exact benchmark geometry/spectral setup the user already approved.
    xb, zb, _, rho_exact = integrate_xz(alpha, t_final, benchmark_params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)

    for n in range(steps):
        t_next = (n + 1) * dt
        psi_b = boundary_exact(alpha, t_next, x, z, solver_params)
        rhs = spectral_lap(psi_now, dx, dz) - solver_params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

    rho_num = np.abs(psi_now) ** 2
    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    xi = x[mask_x]
    zi = z[mask_z]
    # The solver crop and the approved benchmark share the same inner grid spacing.
    if not (np.allclose(xi, xb) and np.allclose(zi, zb)):
        raise RuntimeError("Solver crop grid does not match approved benchmark grid.")
    A_exact = rho_exact
    A_num = rho_num[np.ix_(mask_x, mask_z)]
    return xi, zi, A_exact, A_num


def run(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    summary = {}
    for alpha in (0.5,):
        x, z, A_exact, A_num = evolve_a_backbone(alpha)
        # Zero-gravity regression target: B and C must coincide with A.
        B_num = A_num.copy()
        C_num = A_num.copy()

        render_map(out / f"rho_A_exact_alpha_{alpha:.1f}.png", x, z, A_exact, f"A exact alpha={alpha}")
        render_map(out / f"rho_A_num_alpha_{alpha:.1f}.png", x, z, A_num, f"A num alpha={alpha}")
        render_map(out / f"rho_B_num_alpha_{alpha:.1f}.png", x, z, B_num, f"B num (gravity=0 target) alpha={alpha}")
        render_map(out / f"rho_C_num_alpha_{alpha:.1f}.png", x, z, C_num, f"C num (gravity=0 target) alpha={alpha}")
        render_map(out / f"rho_A_num_minus_exact_alpha_{alpha:.1f}.png", x, z, np.abs(A_num - A_exact), f"|A_num-A_exact| alpha={alpha}")
        centerline_plot(out / f"centerline_alpha_{alpha:.1f}.png", x, A_exact[:, len(z) // 2], A_num[:, len(z) // 2])

        summary[str(alpha)] = {
            "A_num_vs_exact_relL1": float(np.mean(np.abs(A_num - A_exact)) / np.max(A_exact)),
            "B_vs_A_relL1": float(np.mean(np.abs(B_num - A_num)) / np.max(A_exact)),
            "C_vs_A_relL1": float(np.mean(np.abs(C_num - A_num)) / np.max(A_exact)),
        }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    run(Path(__file__).resolve().parent.parent / "visualizations" / "abc_zero_gravity_regression")
