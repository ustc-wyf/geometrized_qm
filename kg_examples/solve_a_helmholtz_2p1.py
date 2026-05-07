from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import eye, kron, diags
from scipy.sparse.linalg import spsolve

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map


def build_helmholtz_matrix(nx: int, nz: int, dx: float, dz: float, k2: float):
    nx_i = nx - 2
    nz_i = nz - 2

    tx = diags(
        [np.ones(nx_i - 1), -2.0 * np.ones(nx_i), np.ones(nx_i - 1)],
        offsets=[-1, 0, 1],
        shape=(nx_i, nx_i),
        format="csr",
    ) / (dx * dx)
    tz = diags(
        [np.ones(nz_i - 1), -2.0 * np.ones(nz_i), np.ones(nz_i - 1)],
        offsets=[-1, 0, 1],
        shape=(nz_i, nz_i),
        format="csr",
    ) / (dz * dz)

    a = kron(eye(nz_i, format="csr"), tx) + kron(tz, eye(nx_i, format="csr")) + k2 * eye(nx_i * nz_i, format="csr")
    return a.tocsr()


def boundary_rhs(phi: np.ndarray, dx: float, dz: float) -> np.ndarray:
    nx, nz = phi.shape
    nx_i = nx - 2
    nz_i = nz - 2
    rhs = np.zeros((nx_i, nz_i), dtype=np.complex128)

    rhs[0, :] -= phi[0, 1:-1] / (dx * dx)
    rhs[-1, :] -= phi[-1, 1:-1] / (dx * dx)
    rhs[:, 0] -= phi[1:-1, 0] / (dz * dz)
    rhs[:, -1] -= phi[1:-1, -1] / (dz * dz)

    return rhs


def solve_once(alpha: float, params: Exact2p1Params):
    t = params.overlap_time
    x, z, psi_exact, rho_exact = integrate_xz(alpha, t, params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    k2 = params.omega0**2 - params.m**2

    a = build_helmholtz_matrix(params.nx, params.nz, dx, dz, k2)
    rhs = boundary_rhs(psi_exact, dx, dz).reshape(-1, order="F")
    psi_interior = spsolve(a, rhs)

    psi_num = psi_exact.copy()
    psi_num[1:-1, 1:-1] = psi_interior.reshape((params.nx - 2, params.nz - 2), order="F")
    rho_num = np.abs(psi_num) ** 2

    diff = rho_num - rho_exact
    return x, z, rho_exact, rho_num, diff


def centerline_plot(path: Path, x: np.ndarray, exact: np.ndarray, numeric: np.ndarray):
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(x, exact, label="exact", linewidth=1.8)
    plt.plot(x, numeric, label="helmholtz", linewidth=1.4)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = Exact2p1Params()
    summary = {
        "params": {
            "extent": params.x_half_range,
            "nx": params.nx,
            "nz": params.nz,
            "nkx": params.nkx,
            "nkz": params.nkz,
            "k0": params.k0,
            "m": params.m,
            "sigma_k_parallel": params.sigma_k_parallel,
            "sigma_k_perp": params.sigma_k_perp,
            "t": params.overlap_time,
        }
    }

    for alpha in (0.0, 0.5, 1.0):
        x, z, rho_exact, rho_num, diff = solve_once(alpha, params)
        render_map(out / f"alpha_{alpha:.1f}_rho_exact.png", x, z, rho_exact, f"exact A, alpha={alpha}")
        render_map(out / f"alpha_{alpha:.1f}_rho_helmholtz.png", x, z, rho_num, f"helmholtz A, alpha={alpha}")
        render_map(out / f"alpha_{alpha:.1f}_rho_abs_diff.png", x, z, np.abs(diff), f"|diff|, alpha={alpha}")

        mid = len(z) // 2
        centerline_plot(
            out / f"alpha_{alpha:.1f}_centerline.png",
            x,
            rho_exact[:, mid],
            rho_num[:, mid],
        )

        summary[str(alpha)] = {
            "l1_mean": float(np.mean(np.abs(diff))),
            "linf": float(np.max(np.abs(diff))),
            "relative_l1": float(np.mean(np.abs(diff)) / np.max(rho_exact)),
            "relative_linf": float(np.max(np.abs(diff)) / np.max(rho_exact)),
        }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = run(Path(__file__).resolve().parent.parent / "visualizations" / "a_helmholtz_solver")
    print(json.dumps(summary, indent=2))
