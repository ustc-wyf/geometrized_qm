from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np


@dataclass
class CrossingExactParams:
    m: float
    alpha: float
    k0: float
    sigma_k: float
    phi0: float
    slit_separation_a: float
    x_half_range: float
    z_half_range: float
    nx: int = 81
    nz: int = 81
    nkx: int = 25
    nky: int = 13
    nkz: int = 49
    k_extent_sigma: float = 6.0
    chunk_x: int = 2

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def sigma_over_k0(self) -> float:
        return float(self.sigma_k / self.k0)


def default_params() -> CrossingExactParams:
    return CrossingExactParams(
        m=1.0,
        alpha=0.5,
        k0=6.0,
        sigma_k=0.6,
        phi0=0.0,
        slit_separation_a=8.0,
        x_half_range=10.0,
        z_half_range=10.0,
    )


def make_grids(params: CrossingExactParams):
    x = np.linspace(-params.x_half_range, params.x_half_range, params.nx)
    z = np.linspace(-params.z_half_range, params.z_half_range, params.nz)
    width = params.k_extent_sigma * params.sigma_k
    kx = np.linspace(-params.k0 - width, params.k0 + width, params.nkx)
    ky = np.linspace(-width, width, params.nky)
    kz = np.linspace(-params.k0 - width, params.k0 + width, params.nkz)
    return x, z, kx, ky, kz


def omega_grid(KX, KY, KZ, m):
    return np.sqrt(KX**2 + KY**2 + KZ**2 + m**2)


def spectral_amplitude(KX, KY, KZ, params: CrossingExactParams):
    kx0 = params.k0 / np.sqrt(2.0)
    kz0 = params.k0 / np.sqrt(2.0)

    xA, yA, zA = -0.5 * params.slit_separation_a, 0.0, -0.5 * params.slit_separation_a
    xB, yB, zB = +0.5 * params.slit_separation_a, 0.0, +0.5 * params.slit_separation_a

    ga = np.exp(-((KX - kx0) ** 2 + KY**2 + (KZ - kz0) ** 2) / (4.0 * params.sigma_k**2))
    gb = np.exp(-((KX + kx0) ** 2 + KY**2 + (KZ + kz0) ** 2) / (4.0 * params.sigma_k**2))

    aA = params.alpha * ga * np.exp(-1j * (KX * xA + KY * yA + KZ * zA))
    aB = np.sqrt(max(0.0, 1.0 - params.alpha)) * gb * np.exp(-1j * (KX * xB + KY * yB + KZ * zB) + 1j * params.phi0)
    return aA + aB


def _integrate3(values, kx, ky, kz):
    tmp = np.trapezoid(values, kz, axis=-1)
    tmp = np.trapezoid(tmp, ky, axis=-1)
    tmp = np.trapezoid(tmp, kx, axis=-1)
    return tmp


def compute_xz_slice(alpha: float, t: float, params: CrossingExactParams):
    params = CrossingExactParams(**{**params.__dict__, "alpha": alpha})
    x, z, kx, ky, kz = make_grids(params)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    omega = omega_grid(KX, KY, KZ, params.m)
    amp = spectral_amplitude(KX, KY, KZ, params)

    rho = np.zeros((len(x), len(z)), dtype=float)
    for i0 in range(0, len(x), params.chunk_x):
        i1 = min(len(x), i0 + params.chunk_x)
        xx = x[i0:i1]
        phase = np.exp(
            1j
            * (
                xx[:, None, None, None, None] * KX[None, None, ...]
                + z[None, :, None, None, None] * KZ[None, None, ...]
            )
            - 1j * t * omega[None, None, ...]
        )
        psi = _integrate3(amp[None, None, ...] * phase, kx, ky, kz)
        rho[i0:i1] = np.abs(psi) ** 2
    return x, z, rho


def run_all(output_dir: str | Path):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    params = default_params()
    t = float(np.sqrt(2.0) * (0.5 * params.slit_separation_a) / (params.k0 / params.omega0))

    summaries = {}
    for alpha in (0.0, 0.5, 1.0):
        x, z, rho = compute_xz_slice(alpha, t, params)
        summaries[str(alpha)] = {
            "rho_max": float(np.max(rho)),
            "rho_boundary_over_max": float(
                np.max(np.concatenate([rho[0, :], rho[-1, :], rho[:, 0], rho[:, -1]])) / np.max(rho)
            ),
        }
        np.savez_compressed(out / f"crossing_exact_alpha_{alpha:.1f}.npz", x=x, z=z, rho=rho)

    meta = {
        "params": {
            "m": params.m,
            "k0": params.k0,
            "sigma_k": params.sigma_k,
            "sigma_over_k0": params.sigma_over_k0,
            "slit_separation_a": params.slit_separation_a,
            "x_half_range": params.x_half_range,
            "z_half_range": params.z_half_range,
            "t": t,
        },
        "summaries": summaries,
    }
    (out / "crossing_exact_summary.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    return meta


if __name__ == "__main__":
    result = run_all(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
