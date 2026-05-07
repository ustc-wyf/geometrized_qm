from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from simulate_bc_psfd_imex import (
    compute_matter_observables,
    edge_damp,
    make_sponge,
    rhs_branch_b,
    spectral_k2,
    wave_imex_step,
)


def render_line(path: Path, x: np.ndarray, y: np.ndarray, title: str, ylabel: str) -> None:
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(x, y, linewidth=1.6)
    plt.xlabel("x")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    nx = 181
    nz = 181
    outer = 20.0
    x = np.linspace(-outer, outer, nx)
    z = np.linspace(-outer, outer, nz)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    k2 = spectral_k2(nx, nz, dx, dz)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)
    dt = 0.05
    steps = 200
    m = 1.0

    # Test 1: freeze geometry, evolve matter only in flat limit.
    rho = np.ones((nx, nz), dtype=float)
    s = np.zeros((nx, nz), dtype=float)
    for _ in range(steps):
        sx, sz, e, n = compute_matter_observables(rho, s, np.zeros_like(rho), np.zeros_like(rho), m, dx, dz)
        fluxx = (n / e) * sx
        fluxz = (n / e) * sz
        n_t = -(np.gradient(fluxx, dx, axis=0) + np.gradient(fluxz, dz, axis=1))
        s_t = -e
        s = s + dt * s_t
        n = n + dt * n_t
        rho = np.maximum(n / e, 1e-12)

    exact_s = -m * dt * steps
    test1 = {
        "rho_deviation_linf": float(np.max(np.abs(rho - 1.0))),
        "s_center": float(s[nx // 2, nz // 2]),
        "s_expected": float(exact_s),
        "s_error_linf": float(np.max(np.abs(s - exact_s))),
    }

    # Test 2: same initial data, but turn on B-geometry with huge M_P.
    rho2 = np.ones((nx, nz), dtype=float)
    s2 = np.zeros((nx, nz), dtype=float)
    tau = np.zeros((nx, nz), dtype=float)
    beta = np.zeros((nx, nz), dtype=float)
    ptau = np.zeros((nx, nz), dtype=float)
    pbeta = np.zeros((nx, nz), dtype=float)

    mp = 1.0e6
    for _ in range(steps):
        src_tau, src_beta, n_t, s_t, e, n = rhs_branch_b(tau, ptau, beta, pbeta, rho2, s2, m, mp, dx, dz, k2)
        tau, ptau = wave_imex_step(tau, ptau, src_tau, dt, k2)
        beta, pbeta = wave_imex_step(beta, pbeta, src_beta, dt, k2)
        edge_damp(tau, sigma, dt, -2.0, 2.0)
        edge_damp(beta, sigma, dt, -2.0, 2.0)
        edge_damp(ptau, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta, sigma, dt, -10.0, 10.0)
        s2 = s2 + dt * s_t
        n = (n + dt * n_t) / (1.0 + 0.5 * sigma * dt)
        _, _, e2, n2 = compute_matter_observables(rho2, s2, tau, beta, m, dx, dz)
        rho2 = np.maximum(n / np.maximum(np.exp(np.clip(beta + tau, -12.0, 12.0)) * e2, 1e-12), 1e-12)

    test2 = {
        "rho_deviation_linf": float(np.max(np.abs(rho2 - 1.0))),
        "tau_linf": float(np.max(np.abs(tau))),
        "beta_linf": float(np.max(np.abs(beta))),
        "s_center": float(s2[nx // 2, nz // 2]),
        "s_expected": float(exact_s),
        "s_error_linf": float(np.max(np.abs(s2 - exact_s))),
    }

    # Test 3: weak gravity, but no sponge and no boundary forcing.
    rho3 = np.ones((nx, nz), dtype=float)
    s3 = np.zeros((nx, nz), dtype=float)
    tau3 = np.zeros((nx, nz), dtype=float)
    beta3 = np.zeros((nx, nz), dtype=float)
    ptau3 = np.zeros((nx, nz), dtype=float)
    pbeta3 = np.zeros((nx, nz), dtype=float)
    for _ in range(steps):
        src_tau, src_beta, n_t, s_t, e, n = rhs_branch_b(tau3, ptau3, beta3, pbeta3, rho3, s3, m, mp, dx, dz, k2)
        tau3, ptau3 = wave_imex_step(tau3, ptau3, src_tau, dt, k2)
        beta3, pbeta3 = wave_imex_step(beta3, pbeta3, src_beta, dt, k2)
        s3 = s3 + dt * s_t
        n = n + dt * n_t
        _, _, e3, n3 = compute_matter_observables(rho3, s3, tau3, beta3, m, dx, dz)
        rho3 = np.maximum(n / np.maximum(np.exp(np.clip(beta3 + tau3, -12.0, 12.0)) * e3, 1e-12), 1e-12)

    test3 = {
        "rho_deviation_linf": float(np.max(np.abs(rho3 - 1.0))),
        "tau_linf": float(np.max(np.abs(tau3))),
        "beta_linf": float(np.max(np.abs(beta3))),
        "s_center": float(s3[nx // 2, nz // 2]),
        "s_expected": float(exact_s),
        "s_error_linf": float(np.max(np.abs(s3 - exact_s))),
    }

    center = nz // 2
    render_line(out / "test1_rho_center.png", x, rho[:, center], "test1 rho centerline", "rho")
    render_line(out / "test2_rho_center.png", x, rho2[:, center], "test2 rho centerline", "rho")
    render_line(out / "test2_tau_center.png", x, tau[:, center], "test2 tau centerline", "tau")
    render_line(out / "test2_beta_center.png", x, beta[:, center], "test2 beta centerline", "beta")

    render_line(out / "test3_rho_center.png", x, rho3[:, center], "test3 rho centerline", "rho")
    render_line(out / "test3_tau_center.png", x, tau3[:, center], "test3 tau centerline", "tau")
    render_line(out / "test3_beta_center.png", x, beta3[:, center], "test3 beta centerline", "beta")

    summary = {"test1_flat_frozen": test1, "test2_flat_weak_gravity": test2, "test3_flat_weak_gravity_no_sponge": test3}
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = run(Path(__file__).resolve().parent.parent / "visualizations" / "debug_b_flat_limit")
    print(json.dumps(summary, indent=2))
