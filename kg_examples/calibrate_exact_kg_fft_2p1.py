from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map


def fft_propagate_positive_frequency(psi0: np.ndarray, dx: float, dz: float, m: float, dt: float):
    nx, nz = psi0.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    omega = np.sqrt(KX**2 + KZ**2 + m**2)
    psi_k = np.fft.fft2(psi0)
    psi_t = np.fft.ifft2(psi_k * np.exp(-1j * omega * dt))
    return psi_t


def centerline_plot(path: Path, x: np.ndarray, exact: np.ndarray, numeric: np.ndarray):
    plt.figure(figsize=(6.2, 3.8))
    plt.plot(x, exact, label="exact", linewidth=1.8)
    plt.plot(x, numeric, label="fft", linewidth=1.4)
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
    t_final = params.overlap_time

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
            "t_final": t_final,
        }
    }

    for alpha in (0.0, 0.5, 1.0):
        x, z, psi0, _ = integrate_xz(alpha, 0.0, params)
        _, _, psi_exact, rho_exact = integrate_xz(alpha, t_final, params)
        dx = float(x[1] - x[0])
        dz = float(z[1] - z[0])

        psi_fft = fft_propagate_positive_frequency(psi0, dx, dz, params.m, t_final)
        rho_fft = np.abs(psi_fft) ** 2
        diff = rho_fft - rho_exact

        render_map(out / f"alpha_{alpha:.1f}_rho_exact.png", x, z, rho_exact, f"exact KG, alpha={alpha}")
        render_map(out / f"alpha_{alpha:.1f}_rho_fft.png", x, z, rho_fft, f"fft KG, alpha={alpha}")
        render_map(out / f"alpha_{alpha:.1f}_rho_abs_diff.png", x, z, np.abs(diff), f"|diff|, alpha={alpha}")

        mid = len(z) // 2
        centerline_plot(out / f"alpha_{alpha:.1f}_centerline.png", x, rho_exact[:, mid], rho_fft[:, mid])

        summary[str(alpha)] = {
            "l1_mean": float(np.mean(np.abs(diff))),
            "linf": float(np.max(np.abs(diff))),
            "relative_l1": float(np.mean(np.abs(diff)) / np.max(rho_exact)),
            "relative_linf": float(np.max(np.abs(diff)) / np.max(rho_exact)),
        }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = run(Path(__file__).resolve().parent.parent / "visualizations" / "a_fft_calibration")
    print(json.dumps(summary, indent=2))
