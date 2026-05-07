from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz, render_map
from simulate_bc_psfd_imex import analytic_slice, apply_bottom_boundary, compute_matter_observables, render_centerline, spectral_dx, spectral_dz


def run(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    inner = 10.0
    outer = 20.0
    params = Exact2p1Params(x_half_range=outer, z_half_range=outer, nx=181, nz=181, nkx=61, nkz=61)
    t_final = params.overlap_time
    steps = 260
    dt = t_final / steps

    x, z, _, rho_a, _ = analytic_slice(0.5, t_final, params)
    _, _, _, rho_b, s_b = analytic_slice(0.5, 0.0, params)

    for n in range(steps):
        t_next = (n + 1) * dt
        rho_bc, s_bc = boundary_exact(0.5, t_next, x, z, params)
        sx = spectral_dx(s_b, 2.0 * outer / (params.nx - 1))
        sz = spectral_dz(s_b, 2.0 * outer / (params.nz - 1))
        e = np.sqrt(np.maximum(sx * sx + sz * sz + params.m**2, 1e-12))
        ncur = np.maximum(rho_b * e, 1e-12)
        fluxx = (ncur / e) * sx
        fluxz = (ncur / e) * sz
        n_t = -(spectral_dx(fluxx, 2.0 * outer / (params.nx - 1)) + spectral_dz(fluxz, 2.0 * outer / (params.nz - 1)))
        s_t = -e
        s_b = s_b + dt * s_t
        ncur = ncur + dt * n_t
        apply_bottom_boundary(s_b, s_bc)
        _, _, e2, _ = compute_matter_observables(rho_b, s_b, np.zeros_like(rho_b), np.zeros_like(rho_b), params.m, 2.0 * outer / (params.nx - 1), 2.0 * outer / (params.nz - 1))
        apply_bottom_boundary(ncur, rho_bc * e2)
        rho_b = np.maximum(ncur / e2, 1e-12)
        apply_bottom_boundary(rho_b, rho_bc)

    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    xi = x[mask_x]
    zi = z[mask_z]
    A = rho_a[np.ix_(mask_x, mask_z)]
    Bflat = rho_b[np.ix_(mask_x, mask_z)]

    render_map(out / "rho_A_inner.png", xi, zi, A, "A exact")
    render_map(out / "rho_Bflat_inner.png", xi, zi, Bflat, "B flat matter-only")
    render_map(out / "rho_abs_diff.png", xi, zi, np.abs(Bflat - A), "|Bflat-A|")
    render_centerline(out / "centerline_compare.png", xi, [("A", A[:, len(zi)//2]), ("Bflat", Bflat[:, len(zi)//2])])

    summary = {
        "errors_inner": {
            "Bflat_vs_A_relL1": float(np.mean(np.abs(Bflat - A)) / np.max(A)),
        }
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


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


if __name__ == "__main__":
    summary = run(Path(__file__).resolve().parent.parent / "visualizations" / "debug_b_boundary_compatibility")
    print(json.dumps(summary, indent=2))
