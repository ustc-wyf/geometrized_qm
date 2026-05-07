from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class Exact2p1Params:
    m: float = 1.0
    alpha: float = 0.5
    k0: float = 10.0
    sigma_k_parallel: float = 0.1
    sigma_k_perp: float = 1.0
    phi0: float = 0.0
    x0: float = 8.0
    z0: float = 8.0
    x_half_range: float = 10.0
    z_half_range: float = 10.0
    nx: int = 121
    nz: int = 121
    nkx: int = 81
    nkz: int = 81
    k_extent_sigma: float = 5.0
    chunk_x: int = 4

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def kx0(self) -> float:
        return float(self.k0 / np.sqrt(2.0))

    @property
    def kz0(self) -> float:
        return float(self.k0 / np.sqrt(2.0))

    @property
    def overlap_time(self) -> float:
        vz = self.kz0 / self.omega0
        return float(self.z0 / vz)

    @property
    def sigma_parallel_over_k0(self) -> float:
        return float(self.sigma_k_parallel / self.k0)

    @property
    def sigma_perp_over_k0(self) -> float:
        return float(self.sigma_k_perp / self.k0)


def make_grids(params: Exact2p1Params):
    x = np.linspace(-params.x_half_range, params.x_half_range, params.nx)
    z = np.linspace(-params.z_half_range, params.z_half_range, params.nz)
    width = params.k_extent_sigma * max(params.sigma_k_parallel, params.sigma_k_perp)
    kx = np.linspace(-params.k0 - width, params.k0 + width, params.nkx)
    kz = np.linspace(-width, params.k0 + width, params.nkz)
    return x, z, kx, kz


def spectral_amplitude(KX, KZ, params: Exact2p1Params):
    # Beam A: +45 degree propagation axis
    dkx_a = KX - params.kx0
    dkz_a = KZ - params.kz0
    dpar_a = (dkx_a + dkz_a) / np.sqrt(2.0)
    dperp_a = (dkx_a - dkz_a) / np.sqrt(2.0)
    ga = np.exp(
        -(dpar_a**2) / (4.0 * params.sigma_k_parallel**2)
        -(dperp_a**2) / (4.0 * params.sigma_k_perp**2)
    )

    # Beam B: -45 degree propagation axis, still moving upward in +z
    dkx_b = KX + params.kx0
    dkz_b = KZ - params.kz0
    dpar_b = (-dkx_b + dkz_b) / np.sqrt(2.0)
    dperp_b = (dkx_b + dkz_b) / np.sqrt(2.0)
    gb = np.exp(
        -(dpar_b**2) / (4.0 * params.sigma_k_parallel**2)
        -(dperp_b**2) / (4.0 * params.sigma_k_perp**2)
    )

    # Two beams injected from z<0, meeting near the origin at t = overlap_time
    a = np.sqrt(params.alpha) * ga * np.exp(-1j * (KX * (-params.x0) + KZ * (-params.z0)))
    b = np.sqrt(max(0.0, 1.0 - params.alpha)) * gb * np.exp(
        -1j * (KX * (+params.x0) + KZ * (-params.z0)) + 1j * params.phi0
    )
    return a + b


def evaluate_points(alpha: float, x_eval, z_eval, t: float, params: Exact2p1Params):
    params = Exact2p1Params(**{**params.__dict__, "alpha": alpha})
    x_eval = np.asarray(x_eval)
    z_eval = np.asarray(z_eval)
    _, _, kx, kz = make_grids(params)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    omega = np.sqrt(KX**2 + KZ**2 + params.m**2)
    amp = spectral_amplitude(KX, KZ, params)

    psi_out = np.zeros(x_eval.shape, dtype=np.complex128)
    flat_x = x_eval.ravel()
    flat_z = z_eval.ravel()
    for i0 in range(0, len(flat_x), params.chunk_x):
        i1 = min(len(flat_x), i0 + params.chunk_x)
        xx = flat_x[i0:i1]
        zz = flat_z[i0:i1]
        phase = np.exp(
            1j
            * (
                xx[:, None, None] * KX[None, :, :]
                + zz[:, None, None] * KZ[None, :, :]
            )
            - 1j * t * omega[None, :, :]
        )
        vals = np.trapezoid(np.trapezoid(amp[None, :, :] * phase, kz, axis=-1), kx, axis=-1)
        psi_out.ravel()[i0:i1] = vals
    return psi_out


def integrate_xz(alpha: float, t: float, params: Exact2p1Params):
    x, z, _, _ = make_grids(params)
    X, Z = np.meshgrid(x, z, indexing="ij")
    psi_out = evaluate_points(alpha, X, Z, t, params)
    rho = np.abs(psi_out) ** 2
    return x, z, psi_out, rho


def render_map(path: Path, x, z, field, title):
    plt.figure(figsize=(5.6, 4.8))
    plt.imshow(field.T, origin="lower", extent=[x.min(), x.max(), z.min(), z.max()], aspect="equal", cmap="viridis")
    plt.colorbar(label="rho")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    params = Exact2p1Params()
    t = params.overlap_time

    summary = {
        "params": {
            "m": params.m,
            "k0": params.k0,
            "sigma_k_parallel": params.sigma_k_parallel,
            "sigma_k_perp": params.sigma_k_perp,
            "sigma_parallel_over_k0": params.sigma_parallel_over_k0,
            "sigma_perp_over_k0": params.sigma_perp_over_k0,
            "x0": params.x0,
            "z0": params.z0,
            "overlap_time": t,
            "nx": params.nx,
            "nz": params.nz,
            "nkx": params.nkx,
            "nkz": params.nkz,
        }
    }

    for alpha in (0.0, 0.5, 1.0):
        x, z, psi, rho = integrate_xz(alpha, t, params)
        np.savez_compressed(out / f"exact2p1_alpha_{alpha:.1f}.npz", x=x, z=z, psi=psi, rho=rho)
        render_map(out / f"exact2p1_alpha_{alpha:.1f}_xz_linear.png", x, z, rho, f"exact 2+1d KG, alpha={alpha}")
        summary[str(alpha)] = {
            "rho_max": float(rho.max()),
            "rho_boundary_over_max": float(np.max(np.concatenate([rho[0, :], rho[-1, :], rho[:, 0], rho[:, -1]])) / rho.max()),
        }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = run(Path(__file__).resolve().parent.parent / "visualizations" / "exact_kg_2p1_benchmark")
    print(json.dumps(result, indent=2))
