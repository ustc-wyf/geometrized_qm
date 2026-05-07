from __future__ import annotations

import json
from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz, render_map


def lap(f: np.ndarray, dx: float, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    # Start with a conservative 2nd-order Laplacian everywhere
    out[1:-1] += (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dx * dx)
    out[:, 1:-1] += (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / (dz * dz)
    # Upgrade the well-resolved bulk to 4th-order
    out[2:-2, 2:-2] = (
        (-f[4:, 2:-2] + 16.0 * f[3:-1, 2:-2] - 30.0 * f[2:-2, 2:-2] + 16.0 * f[1:-3, 2:-2] - f[:-4, 2:-2])
        / (12.0 * dx * dx)
        + (-f[2:-2, 4:] + 16.0 * f[2:-2, 3:-1] - 30.0 * f[2:-2, 2:-2] + 16.0 * f[2:-2, 1:-3] - f[2:-2, :-4])
        / (12.0 * dz * dz)
    )
    out[0] = out[1]
    out[-1] = out[-2]
    out[:, 0] = out[:, 1]
    out[:, -1] = out[:, -2]
    return out


def make_sponge(x: np.ndarray, z: np.ndarray, inner_half_range: float, outer_half_range: float, sigma_max: float) -> np.ndarray:
    X, Z = np.meshgrid(x, z, indexing="ij")
    sx = np.zeros_like(X)
    sz = np.zeros_like(Z)

    wx = max(outer_half_range - inner_half_range, 1e-6)
    wz = max(outer_half_range - inner_half_range, 1e-6)

    mask_x = np.abs(X) > inner_half_range
    mask_z = np.abs(Z) > inner_half_range
    sx[mask_x] = ((np.abs(X[mask_x]) - inner_half_range) / wx) ** 2
    sz[mask_z] = ((np.abs(Z[mask_z]) - inner_half_range) / wz) ** 2
    return sigma_max * (sx + sz)


def centerline_plot(path: Path, x: np.ndarray, exact: np.ndarray, numeric: np.ndarray) -> None:
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(x, exact, label="exact", linewidth=1.8)
    plt.plot(x, numeric, label="fd+sponge", linewidth=1.4)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


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


def run(outdir: str | Path, inner=10.0, outer=20.0, nx=241, nz=241, nkx=81, nkz=81, dt0=0.025, sigma_max=2.0):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = Exact2p1Params(
        x_half_range=outer,
        z_half_range=outer,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )

    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    dt = dt0
    t_final = params.overlap_time
    steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, psi_prev, _ = integrate_xz(0.5, -dt, params)
    _, _, psi_now, _ = integrate_xz(0.5, 0.0, params)
    _, _, psi_exact, rho_exact = integrate_xz(0.5, t_final, params)

    sponge_inner = inner + 4.0
    sigma = make_sponge(x, z, inner_half_range=sponge_inner, outer_half_range=outer, sigma_max=sigma_max)

    for _ in range(steps):
        t_next = (_ + 1) * dt
        psi_b = boundary_exact(0.5, t_next, x, z, params)
        rhs = lap(psi_now, dx, dz) - params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        denom = 1.0 + 0.5 * sigma * dt
        psi_next = numer / denom
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

    rho_num = np.abs(psi_now) ** 2
    diff = rho_num - rho_exact

    mask_x = np.abs(x) <= inner
    mask_z = np.abs(z) <= inner
    inner_exact = rho_exact[np.ix_(mask_x, mask_z)]
    inner_num = rho_num[np.ix_(mask_x, mask_z)]
    inner_diff = inner_num - inner_exact

    render_map(out / "rho_exact_full.png", x, z, rho_exact, "exact KG (full domain)")
    render_map(out / "rho_numeric_full.png", x, z, rho_num, "fd+sponge KG (full domain)")
    render_map(out / "rho_abs_diff_full.png", x, z, np.abs(diff), "|diff| (full domain)")

    xi = x[mask_x]
    zi = z[mask_z]
    render_map(out / "rho_exact_inner.png", xi, zi, inner_exact, "exact KG (inner window)")
    render_map(out / "rho_numeric_inner.png", xi, zi, inner_num, "fd+sponge KG (inner window)")
    render_map(out / "rho_abs_diff_inner.png", xi, zi, np.abs(inner_diff), "|diff| (inner window)")

    mid = len(zi) // 2
    centerline_plot(out / "centerline_inner.png", xi, inner_exact[:, mid], inner_num[:, mid])

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
            "sigma_max": sigma_max,
            "sponge_inner": sponge_inner,
            "k0": params.k0,
            "m": params.m,
            "sigma_k_parallel": params.sigma_k_parallel,
            "sigma_k_perp": params.sigma_k_perp,
        },
        "errors_full": {
            "l1_mean": float(np.mean(np.abs(diff))),
            "linf": float(np.max(np.abs(diff))),
            "relative_l1": float(np.mean(np.abs(diff)) / np.max(rho_exact)),
            "relative_linf": float(np.max(np.abs(diff)) / np.max(rho_exact)),
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default=str(Path(__file__).resolve().parent.parent / "visualizations" / "a_fd_absorbing_calibration"))
    parser.add_argument("--inner", type=float, default=10.0)
    parser.add_argument("--outer", type=float, default=20.0)
    parser.add_argument("--nx", type=int, default=241)
    parser.add_argument("--nz", type=int, default=241)
    parser.add_argument("--nkx", type=int, default=81)
    parser.add_argument("--nkz", type=int, default=81)
    parser.add_argument("--dt", dest="dt0", type=float, default=0.025)
    parser.add_argument("--sigma-max", type=float, default=2.0)
    args = parser.parse_args()
    summary = run(
        args.outdir,
        inner=args.inner,
        outer=args.outer,
        nx=args.nx,
        nz=args.nz,
        nkx=args.nkx,
        nkz=args.nkz,
        dt0=args.dt0,
        sigma_max=args.sigma_max,
    )
    print(json.dumps(summary, indent=2))
