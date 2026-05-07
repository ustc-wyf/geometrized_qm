from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np


@dataclass
class CrossingGaussianBeamParams:
    m: float
    alpha: float
    k0: float
    w0: float
    theta: float
    phi0: float
    x_half_range: float
    z_half_range: float
    nx: int = 401
    nz: int = 401

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def sigma_over_k0(self) -> float:
        # w0 and k0 satisfy paraxial-like narrow-angle criterion k0*w0 >> 1
        return float(1.0 / (self.k0 * self.w0))

    @property
    def zR(self) -> float:
        return float(self.k0 * self.w0**2 / 2.0)


def default_params() -> CrossingGaussianBeamParams:
    return CrossingGaussianBeamParams(
        m=1.0,
        alpha=0.5,
        k0=10.0,
        w0=1.0,
        theta=np.pi / 4.0,
        phi0=0.0,
        x_half_range=5.0,
        z_half_range=5.0,
    )


def make_grids(params: CrossingGaussianBeamParams) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.linspace(-params.x_half_range, params.x_half_range, params.nx)
    z = np.linspace(-params.z_half_range, params.z_half_range, params.nz)
    X, Z = np.meshgrid(x, z, indexing="ij")
    return x, z, X, Z


def rotated_coords(X: np.ndarray, Z: np.ndarray, theta: float, sign: int) -> tuple[np.ndarray, np.ndarray]:
    # Two distinct beam axes crossing at the origin:
    # sign=+1  -> axis along z = +x   (for theta=pi/4)
    # sign=-1  -> axis along z = -x   (for theta=pi/4)
    c = np.cos(theta)
    s = np.sin(theta)
    if sign == +1:
        # longitudinal axis e1 = (sinθ, cosθ)
        u = s * X + c * Z
        # transverse coordinate p1 = (cosθ, -sinθ)
        v = c * X - s * Z
    else:
        # longitudinal axis e2 = (-sinθ, cosθ)
        u = -s * X + c * Z
        # transverse coordinate p2 = (cosθ, +sinθ)
        v = c * X + s * Z
    return u, v


def gaussian_beam_rotated(X: np.ndarray, Z: np.ndarray, params: CrossingGaussianBeamParams, sign: int) -> np.ndarray:
    u, v = rotated_coords(X, Z, params.theta, sign)
    zR = params.zR
    w = params.w0 * np.sqrt(1.0 + (u / zR) ** 2)
    # curvature radius with safe handling at waist
    R = np.where(np.abs(u) < 1e-14, np.inf, u * (1.0 + (zR / u) ** 2))
    gouy = 0.5 * np.arctan2(u, zR)

    amp = np.sqrt(params.w0 / w) * np.exp(-(v**2) / (2.0 * w**2))
    phase = params.k0 * u + np.where(np.isfinite(R), params.k0 * v**2 / (2.0 * R), 0.0) - gouy
    return amp * np.exp(1j * phase)


def compute_xz_slice(t: float, params: CrossingGaussianBeamParams) -> dict[str, np.ndarray]:
    x, z, X, Z = make_grids(params)

    common_t = np.exp(-1j * params.omega0 * t)
    psi_a = params.alpha * gaussian_beam_rotated(X, Z, params, +1) * common_t
    psi_b = np.sqrt(max(0.0, 1.0 - params.alpha)) * gaussian_beam_rotated(X, Z, params, -1) * np.exp(1j * params.phi0) * common_t
    psi = psi_a + psi_b
    rho = np.abs(psi) ** 2
    return {"x": x, "z": z, "X": X, "Z": Z, "psi": psi, "rho": rho}


def boundary_ratio(rho: np.ndarray) -> float:
    top = rho[:, -1]
    bottom = rho[:, 0]
    left = rho[0, :]
    right = rho[-1, :]
    boundary = np.concatenate([top, bottom, left, right])
    return float(np.max(boundary) / np.max(rho))


def run_case(output_dir: str | Path, alpha: float) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = default_params()
    params.alpha = alpha
    # monochromatic beam: intensity is stationary in time, choose t=0
    t = 0.0
    data = compute_xz_slice(t, params)
    rho = data["rho"]

    payload = {
        "params": {
            "m": params.m,
            "alpha": params.alpha,
            "k0": params.k0,
            "omega0": params.omega0,
            "w0": params.w0,
            "sigma_over_k0": params.sigma_over_k0,
            "theta_deg": float(np.degrees(params.theta)),
            "phi0": params.phi0,
            "x_min": float(data["x"].min()),
            "x_max": float(data["x"].max()),
            "z_min": float(data["z"].min()),
            "z_max": float(data["z"].max()),
            "t": t,
        },
        "summary": {
            "rho_max": float(np.max(rho)),
            "rho_boundary_over_max": boundary_ratio(rho),
        },
    }

    (out / f"double_slit_alpha_{alpha:.1f}_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        out / f"double_slit_alpha_{alpha:.1f}_fields.npz",
        x=data["x"],
        z=data["z"],
        rho=rho,
        psi=data["psi"],
    )
    return payload


if __name__ == "__main__":
    out = Path(__file__).resolve().parent / "outputs"
    all_results = {str(alpha): run_case(out, alpha) for alpha in (0.0, 0.5, 1.0)}
    print(json.dumps(all_results, ensure_ascii=False, indent=2))
