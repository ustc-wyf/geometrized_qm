from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz, render_map


def spectral_lap(f: np.ndarray, dx: float, dz: float) -> np.ndarray:
    nx, nz = f.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    fhat = np.fft.fft2(f)
    return np.fft.ifft2(-(KX * KX + KZ * KZ) * fhat)


def make_sponge(x: np.ndarray, z: np.ndarray, inner_half_range: float, outer_half_range: float, sigma_max: float) -> np.ndarray:
    X, Z = np.meshgrid(x, z, indexing="ij")
    sx = np.zeros_like(X)
    sz = np.zeros_like(Z)
    wx = max(outer_half_range - inner_half_range, 1e-6)
    wz = max(outer_half_range - inner_half_range, 1e-6)
    mx = np.abs(X) > inner_half_range
    mz = np.abs(Z) > inner_half_range
    sx[mx] = ((np.abs(X[mx]) - inner_half_range) / wx) ** 2
    sz[mz] = ((np.abs(Z[mz]) - inner_half_range) / wz) ** 2
    return sigma_max * (sx + sz)


def boundary_exact(alpha: float, t: float, x: np.ndarray, z: np.ndarray, params: Exact2p1Params) -> np.ndarray:
    nx = len(x)
    nz = len(z)
    psi = np.zeros((nx, nz), dtype=np.complex128)
    psi[:, 0] = evaluate_points(alpha, x, np.full_like(x, z[0]), t, params)
    psi[:, -1] = evaluate_points(alpha, x, np.full_like(x, z[-1]), t, params)
    psi[0, :] = evaluate_points(alpha, np.full_like(z, x[0]), z, t, params)
    psi[-1, :] = evaluate_points(alpha, np.full_like(z, x[-1]), z, t, params)
    return psi


def apply_boundary(field: np.ndarray, target: np.ndarray) -> None:
    field[0, :] = target[0, :]
    field[-1, :] = target[-1, :]
    field[:, 0] = target[:, 0]
    field[:, -1] = target[:, -1]


def centerline_plot(path: Path, x: np.ndarray, exact: np.ndarray, numeric: np.ndarray) -> None:
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(x, exact, label="exact", linewidth=1.8)
    plt.plot(x, numeric, label="ps+fd", linewidth=1.4)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    inner = 10.0
    outer = 20.0
    params = Exact2p1Params(
        x_half_range=outer,
        z_half_range=outer,
        nx=241,
        nz=241,
        nkx=81,
        nkz=81,
    )
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    dt = 0.025
    t_final = params.overlap_time
    steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, psi_prev, _ = integrate_xz(0.5, -dt, params)
    _, _, psi_now, _ = integrate_xz(0.5, 0.0, params)
    _, _, psi_exact, rho_exact = integrate_xz(0.5, t_final, params)
    sigma = make_sponge(x, z, inner_half_range=14.0, outer_half_range=outer, sigma_max=2.0)

    for n in range(steps):
        t_next = (n + 1) * dt
        psi_b = boundary_exact(0.5, t_next, x, z, params)
        rhs = spectral_lap(psi_now, dx, dz) - params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

    rho_num = np.abs(psi_now) ** 2
    diff = rho_num - rho_exact

    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    xi = x[mask_x]
    zi = z[mask_z]
    inner_exact = rho_exact[np.ix_(mask_x, mask_z)]
    inner_num = rho_num[np.ix_(mask_x, mask_z)]
    inner_diff = inner_num - inner_exact

    render_map(out / "rho_exact_inner.png", xi, zi, inner_exact, "exact KG (inner)")
    render_map(out / "rho_numeric_inner.png", xi, zi, inner_num, "ps+fd KG (inner)")
    render_map(out / "rho_abs_diff_inner.png", xi, zi, np.abs(inner_diff), "|diff| (inner)")
    centerline_plot(out / "centerline_inner.png", xi, inner_exact[:, len(zi) // 2], inner_num[:, len(zi) // 2])

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
        },
        "errors_inner": {
            "l1_mean": float(np.mean(np.abs(inner_diff))),
            "linf": float(np.max(np.abs(inner_diff))),
            "relative_l1": float(np.mean(np.abs(inner_diff)) / np.max(inner_exact)),
            "relative_linf": float(np.max(np.abs(inner_diff)) / np.max(inner_exact)),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = run(Path(__file__).resolve().parent.parent / "visualizations" / "a_psfd_absorbing_calibration")
    print(json.dumps(summary, indent=2))
