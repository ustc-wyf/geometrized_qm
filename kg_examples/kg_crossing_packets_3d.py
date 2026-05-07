from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np


@dataclass
class KGCrossingParams:
    m: float
    alpha: float
    k0: float
    sigma_k: float
    theta: float
    phi0: float
    slit_separation_a: float
    x_half_range: float
    z_half_range: float
    nkx: int = 121
    nkz: int = 241
    k_extent_sigma: float = 8.0

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def sigma_over_k0(self) -> float:
        return float(self.sigma_k / self.k0)


def default_crossing_params() -> KGCrossingParams:
    return KGCrossingParams(
        m=1.0,
        alpha=0.5,
        k0=6.0,
        sigma_k=0.6,
        theta=np.pi / 4.0,
        phi0=0.0,
        slit_separation_a=16.0,
        x_half_range=28.0,
        z_half_range=28.0,
    )


def make_k_grids(params: KGCrossingParams) -> tuple[np.ndarray, np.ndarray]:
    width = params.k_extent_sigma * params.sigma_k
    kx = np.linspace(-params.k0 - width, params.k0 + width, params.nkx)
    kz = np.linspace(-params.k0 - width, params.k0 + width, params.nkz)
    return kx, kz


def omega_2d(KX: np.ndarray, KZ: np.ndarray, m: float) -> np.ndarray:
    return np.sqrt(KX**2 + KZ**2 + m**2)


def spectral_amplitude_crossing(kx: np.ndarray, kz: np.ndarray, params: KGCrossingParams) -> np.ndarray:
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    kx0 = params.k0 * np.cos(params.theta)
    kz0 = params.k0 * np.sin(params.theta)

    xA, zA = -0.5 * params.slit_separation_a, -0.5 * params.slit_separation_a
    xB, zB = +0.5 * params.slit_separation_a, +0.5 * params.slit_separation_a

    ga = np.exp(-((KX - kx0) ** 2 + (KZ - kz0) ** 2) / (4.0 * params.sigma_k**2))
    gb = np.exp(-((KX + kx0) ** 2 + (KZ + kz0) ** 2) / (4.0 * params.sigma_k**2))

    aA = params.alpha * ga * np.exp(-1j * (KX * xA + KZ * zA))
    aB = np.sqrt(max(0.0, 1.0 - params.alpha)) * gb * np.exp(-1j * (KX * xB + KZ * zB) + 1j * params.phi0)
    return aA + aB


def _double_trapezoid(values: np.ndarray, kx: np.ndarray, kz: np.ndarray) -> np.ndarray:
    tmp = np.trapezoid(values, kz, axis=-1)
    tmp = np.trapezoid(tmp, kx, axis=-1)
    return tmp


def compute_slice(x: np.ndarray, z: np.ndarray, t: float, params: KGCrossingParams) -> dict[str, np.ndarray]:
    kx, kz = make_k_grids(params)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    omega = omega_2d(KX, KZ, params.m)
    amp = spectral_amplitude_crossing(kx, kz, params)

    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    phase = np.exp(
        1j * (Xg[..., None, None] * KX[None, None, ...] + Zg[..., None, None] * KZ[None, None, ...])
        - 1j * t * omega[None, None, ...]
    )
    base = amp[None, None, ...] * phase

    psi = _double_trapezoid(base, kx, kz)
    rho = np.abs(psi) ** 2
    return {"x": x, "z": z, "psi": psi, "rho": rho}


def boundary_ratio(rho: np.ndarray) -> float:
    top = rho[:, -1]
    bottom = rho[:, 0]
    left = rho[0, :]
    right = rho[-1, :]
    boundary = np.concatenate([top, bottom, left, right])
    return float(np.max(boundary) / np.max(rho))


def run_crossing_case(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = default_crossing_params()
    x = np.linspace(-params.x_half_range, params.x_half_range, 241)
    z = np.linspace(-params.z_half_range, params.z_half_range, 241)
    t = float(np.sqrt(2.0) * (0.5 * params.slit_separation_a) / (params.k0 / params.omega0))

    data = compute_slice(x, z, t, params)
    rho = data["rho"]

    payload = {
        "params": {
            "m": params.m,
            "alpha": params.alpha,
            "k0": params.k0,
            "sigma_k": params.sigma_k,
            "sigma_over_k0": params.sigma_over_k0,
            "theta_rad": float(params.theta),
            "theta_deg": float(np.degrees(params.theta)),
            "phi0": params.phi0,
            "slit_separation_a": params.slit_separation_a,
            "x_min": float(x.min()),
            "x_max": float(x.max()),
            "z_min": float(z.min()),
            "z_max": float(z.max()),
            "t": t,
        },
        "summary": {
            "rho_max": float(np.max(rho)),
            "rho_boundary_over_max": boundary_ratio(rho),
        },
    }

    (out / "crossing3d_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    np.savez_compressed(out / "crossing3d_fields.npz", x=x, z=z, rho=rho, psi=data["psi"])
    return payload


if __name__ == "__main__":
    result = run_crossing_case(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
