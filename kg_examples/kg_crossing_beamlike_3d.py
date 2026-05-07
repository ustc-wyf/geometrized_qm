from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np


@dataclass
class CrossingBeamlikeParams:
    m: float
    alpha: float
    k0: float
    sigma_x: float
    sigma_z: float
    theta: float
    phi0: float
    slit_separation_a: float
    x_half_range: float
    z_half_range: float
    nx: int = 301
    nz: int = 301

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def v_group(self) -> float:
        return float(self.k0 / self.omega0)

    @property
    def sigma_over_lambda_like(self) -> float:
        return float(1.0 / (self.k0 * min(self.sigma_x, self.sigma_z)))


def default_params() -> CrossingBeamlikeParams:
    return CrossingBeamlikeParams(
        m=1.0,
        alpha=0.5,
        k0=10.0,
        sigma_x=2.2,
        sigma_z=2.2,
        theta=np.pi / 4.0,
        phi0=0.0,
        slit_separation_a=16.0,
        x_half_range=28.0,
        z_half_range=28.0,
    )


def make_grids(params: CrossingBeamlikeParams) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.linspace(-params.x_half_range, params.x_half_range, params.nx)
    z = np.linspace(-params.z_half_range, params.z_half_range, params.nz)
    X, Z = np.meshgrid(x, z, indexing="ij")
    return x, z, X, Z


def packet_centers(t: float, params: CrossingBeamlikeParams) -> tuple[tuple[float, float], tuple[float, float]]:
    c = np.cos(params.theta)
    s = np.sin(params.theta)
    v = params.v_group
    xA = -0.5 * params.slit_separation_a + v * c * t
    zA = -0.5 * params.slit_separation_a + v * s * t
    xB = +0.5 * params.slit_separation_a - v * c * t
    zB = +0.5 * params.slit_separation_a - v * s * t
    return (xA, zA), (xB, zB)


def compute_slice(t: float, params: CrossingBeamlikeParams) -> dict[str, np.ndarray]:
    x, z, X, Z = make_grids(params)
    c = np.cos(params.theta)
    s = np.sin(params.theta)

    (xA, zA), (xB, zB) = packet_centers(t, params)

    envA = np.exp(-((X - xA) ** 2) / (4.0 * params.sigma_x**2) - ((Z - zA) ** 2) / (4.0 * params.sigma_z**2))
    envB = np.exp(-((X - xB) ** 2) / (4.0 * params.sigma_x**2) - ((Z - zB) ** 2) / (4.0 * params.sigma_z**2))

    phaseA = np.exp(1j * (params.k0 * (c * X + s * Z) - params.omega0 * t))
    phaseB = np.exp(1j * (-params.k0 * (c * X + s * Z) - params.omega0 * t + params.phi0))

    psiA = params.alpha * envA * phaseA
    psiB = np.sqrt(max(0.0, 1.0 - params.alpha)) * envB * phaseB
    psi = psiA + psiB
    rho = np.abs(psi) ** 2

    return {
        "x": x,
        "z": z,
        "X": X,
        "Z": Z,
        "psi": psi,
        "rho": rho,
        "center_A": np.array([xA, zA]),
        "center_B": np.array([xB, zB]),
    }


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

    # overlap time: centers meet near origin
    t_overlap = float((0.5 * params.slit_separation_a) / (params.v_group * np.cos(params.theta)))
    data = compute_slice(t_overlap, params)
    rho = data["rho"]

    payload = {
        "params": {
            "m": params.m,
            "alpha": params.alpha,
            "k0": params.k0,
            "omega0": params.omega0,
            "sigma_x": params.sigma_x,
            "sigma_z": params.sigma_z,
            "theta_deg": float(np.degrees(params.theta)),
            "phi0": params.phi0,
            "slit_separation_a": params.slit_separation_a,
            "x_min": float(data["x"].min()),
            "x_max": float(data["x"].max()),
            "z_min": float(data["z"].min()),
            "z_max": float(data["z"].max()),
            "t_overlap": t_overlap,
            "center_A": data["center_A"].tolist(),
            "center_B": data["center_B"].tolist(),
        },
        "summary": {
            "rho_max": float(np.max(rho)),
            "rho_boundary_over_max": boundary_ratio(rho),
        },
    }

    (out / f"beamlike_alpha_{alpha:.1f}_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        out / f"beamlike_alpha_{alpha:.1f}_fields.npz",
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
