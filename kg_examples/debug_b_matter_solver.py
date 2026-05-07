from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz
from simulate_bc_psfd_imex import spectral_dx, spectral_dz


def ddx(f, dx):
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def ddz(f, dz):
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out


def rusanov_flux_1d(q_l, q_r, u_l, u_r):
    a = np.maximum(np.abs(u_l), np.abs(u_r))
    f_l = u_l * q_l
    f_r = u_r * q_r
    return 0.5 * (f_l + f_r) - 0.5 * a * (q_r - q_l)


def step_fv_rusanov(rho, s, m, dx, dz, dt):
    sx = ddx(s, dx)
    sz = ddz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m, 1e-12))
    n = np.maximum(rho * e, 1e-12)
    vx = sx / e
    vz = sz / e

    # Conservative update for n with simple transmissive edge states.
    qx_l = n[:-1, :]
    qx_r = n[1:, :]
    vx_l = vx[:-1, :]
    vx_r = vx[1:, :]
    fx = rusanov_flux_1d(qx_l, qx_r, vx_l, vx_r)

    qz_l = n[:, :-1]
    qz_r = n[:, 1:]
    vz_l = vz[:, :-1]
    vz_r = vz[:, 1:]
    fz = rusanov_flux_1d(qz_l, qz_r, vz_l, vz_r)

    n1 = n.copy()
    n1[1:-1, 1:-1] -= (dt / dx) * (fx[1:, 1:-1] - fx[:-1, 1:-1])
    n1[1:-1, 1:-1] -= (dt / dz) * (fz[1:-1, 1:] - fz[1:-1, :-1])

    # Zero-gradient edge update to avoid injecting boundary artifacts in the debug test.
    n1[0, :] = n1[1, :]
    n1[-1, :] = n1[-2, :]
    n1[:, 0] = n1[:, 1]
    n1[:, -1] = n1[:, -2]
    n1 = np.maximum(n1, 1e-12)

    s_t = -e
    s1 = s + dt * s_t
    sx1 = ddx(s1, dx)
    sz1 = ddz(s1, dz)
    e1 = np.sqrt(np.maximum(sx1 * sx1 + sz1 * sz1 + m * m, 1e-12))
    rho1 = np.maximum(n1 / e1, 1e-12)
    return rho1, s1, n, n1


def render_line(path: Path, xs, ys, xlabel, ylabel, title):
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(xs, ys, linewidth=1.5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def step_current(rho, s, m, dx, dz, dt):
    sx = spectral_dx(s, dx)
    sz = spectral_dz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m, 1e-12))
    n = np.maximum(rho * e, 1e-12)
    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(spectral_dx(fluxx, dx) + spectral_dz(fluxz, dz))
    s_t = -e
    s1 = s + dt * s_t
    n1 = n + dt * n_t
    sx1 = spectral_dx(s1, dx)
    sz1 = spectral_dz(s1, dz)
    e1 = np.sqrt(np.maximum(sx1 * sx1 + sz1 * sz1 + m * m, 1e-12))
    rho1 = np.maximum(n1 / e1, 1e-12)
    return rho1, s1, n, n1


def step_localgrad(rho, s, m, dx, dz, dt):
    sx = ddx(s, dx)
    sz = ddz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m, 1e-12))
    n = np.maximum(rho * e, 1e-12)
    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
    s_t = -e
    s1 = s + dt * s_t
    n1 = n + dt * n_t
    sx1 = ddx(s1, dx)
    sz1 = ddz(s1, dz)
    e1 = np.sqrt(np.maximum(sx1 * sx1 + sz1 * sz1 + m * m, 1e-12))
    rho1 = np.maximum(n1 / e1, 1e-12)
    return rho1, s1, n, n1


def main():
    out = Path(__file__).resolve().parent.parent / "visualizations" / "debug_b_matter_solver"
    out.mkdir(parents=True, exist_ok=True)

    params = Exact2p1Params(x_half_range=20.0, z_half_range=20.0, nx=181, nz=181, nkx=61, nkz=61)
    x, z, psi0, rho0 = integrate_xz(1.0, 0.0, params)
    s0 = np.unwrap(np.unwrap(np.angle(psi0), axis=0), axis=1)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dt = 0.01
    steps = 200

    results = {}
    for name, stepper in [("current", step_current), ("localgrad", step_localgrad), ("fv_rusanov", step_fv_rusanov)]:
        rho = rho0.copy()
        s = s0.copy()
        total_n_hist = []
        rho_max_hist = []
        rho_min_hist = []
        for _ in range(steps):
            rho, s, n, n1 = stepper(rho, s, params.m, dx, dz, dt)
            total_n_hist.append(float(np.sum(n1) * dx * dz))
            rho_max_hist.append(float(np.max(rho)))
            rho_min_hist.append(float(np.min(rho)))
        results[name] = {
            "total_n_start": total_n_hist[0],
            "total_n_end": total_n_hist[-1],
            "total_n_rel_drift": float(abs(total_n_hist[-1] - total_n_hist[0]) / max(abs(total_n_hist[0]), 1e-12)),
            "rho_max_end": rho_max_hist[-1],
            "rho_min_end": rho_min_hist[-1],
        }
        render_line(out / f"{name}_total_n.png", np.arange(steps), np.array(total_n_hist), "step", "sum n", f"{name}: total n")
        render_line(out / f"{name}_rho_max.png", np.arange(steps), np.array(rho_max_hist), "step", "max rho", f"{name}: rho max")

    (out / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
