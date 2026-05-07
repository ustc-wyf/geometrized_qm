from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from simulate_bc_small_gravity_backbone import build_exact_background
from simulate_bc_psfd_imex import rhs_branch_c, spectral_k2


def chi_of_phi(phi: np.ndarray, ell: float) -> np.ndarray:
    return np.sqrt(np.maximum(np.power(phi, -2.0 / 3.0) - 1.0, 0.0)) / (ell * ell)


def dchi_dphi(phi: np.ndarray, ell: float) -> np.ndarray:
    inside = np.maximum(np.power(phi, -2.0 / 3.0) - 1.0, 1e-30)
    return -(1.0 / (3.0 * ell * ell)) * np.power(phi, -5.0 / 3.0) / np.sqrt(inside)


def main():
    out = Path(__file__).resolve().parent.parent / "visualizations" / "debug_c_phi_pathology"
    out.mkdir(parents=True, exist_ok=True)

    alpha = 0.5
    ell = 0.02
    benchmark_params, solver_params, xb, zb, rho_a_exact, s_a_exact, x, z, rho_a_outer, s_a_outer = build_exact_background(alpha)
    dx = 2.0 * solver_params.x_half_range / (solver_params.nx - 1)
    dz = 2.0 * solver_params.z_half_range / (solver_params.nz - 1)
    k2 = spectral_k2(solver_params.nx, solver_params.nz, dx, dz)

    tau = np.zeros_like(rho_a_outer)
    beta = np.zeros_like(rho_a_outer)
    phi = np.ones_like(rho_a_outer)
    ptau = np.zeros_like(rho_a_outer)
    pbeta = np.zeros_like(rho_a_outer)
    pphi = np.zeros_like(rho_a_outer)

    src_tau, src_beta, src_phi, *_ = rhs_branch_c(
        tau, ptau, beta, pbeta, phi, pphi, rho_a_outer, s_a_outer, solver_params.m, 300.0, ell, dx, dz
    )

    dt = 0.025
    trial_phi = phi + dt * src_phi
    clipped_phi = np.clip(trial_phi, 0.2, 1.0)

    scan_phi = 1.0 - np.array([1e-2, 1e-4, 1e-6, 1e-8, 1e-10], dtype=float)
    scan = {
        f"{p:.10f}": {
            "chi": float(chi_of_phi(np.array([p]), ell)[0]),
            "dchi_dphi": float(dchi_dphi(np.array([p]), ell)[0]),
        }
        for p in scan_phi
    }

    mask_x = np.abs(x) <= benchmark_params.x_half_range
    mask_z = np.abs(z) <= benchmark_params.z_half_range
    xi = x[mask_x]
    zi = z[mask_z]
    src_phi_inner = src_phi[np.ix_(mask_x, mask_z)]
    trial_inner = trial_phi[np.ix_(mask_x, mask_z)]
    clipped_inner = clipped_phi[np.ix_(mask_x, mask_z)]

    for name, field in [
        ("src_phi_inner.png", src_phi_inner),
        ("trial_phi_inner.png", trial_inner),
        ("clipped_phi_inner.png", clipped_inner),
    ]:
        plt.figure(figsize=(5.6, 4.8))
        plt.imshow(field.T, origin="lower", extent=[xi.min(), xi.max(), zi.min(), zi.max()], aspect="equal", cmap="viridis")
        plt.colorbar()
        plt.xlabel("x")
        plt.ylabel("z")
        plt.tight_layout()
        plt.savefig(out / name, dpi=180)
        plt.close()

    summary = {
        "src_phi_sign": {
            "positive_fraction": float(np.mean(src_phi > 0)),
            "negative_fraction": float(np.mean(src_phi < 0)),
            "max": float(np.max(src_phi)),
            "min": float(np.min(src_phi)),
            "mean_abs": float(np.mean(np.abs(src_phi))),
        },
        "one_step_trial": {
            "fraction_above_1": float(np.mean(trial_phi > 1.0)),
            "fraction_below_0_2": float(np.mean(trial_phi < 0.2)),
            "trial_max": float(np.max(trial_phi)),
            "trial_min": float(np.min(trial_phi)),
        },
        "one_step_clip_effect": {
            "fraction_changed_by_clip": float(np.mean(np.abs(clipped_phi - trial_phi) > 0)),
            "max_clip_delta": float(np.max(np.abs(clipped_phi - trial_phi))),
        },
        "phi_to_chi_scan": scan,
    }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
