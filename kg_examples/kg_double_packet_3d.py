from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np


@dataclass
class KGDoublePacket3DParams:
    m: float
    alpha: float
    sigma_kx: float
    sigma_ky: float
    sigma_kz: float
    k0: float
    phi0: float
    slit_separation_a: float
    screen_half_width_b: float
    z0: float
    z_screen: float
    k_extent_sigma: float = 8.0
    nkx: int = 61
    nky: int = 41
    nkz: int = 121
    nx_screen: int = 801

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def v_group(self) -> float:
        return float(self.k0 / self.omega0)

    @property
    def x_screen(self) -> np.ndarray:
        return np.linspace(-20.0, 20.0, self.nx_screen)


def default_params_3d() -> KGDoublePacket3DParams:
    return KGDoublePacket3DParams(
        m=1.0,
        alpha=0.6,
        sigma_kx=0.20,
        sigma_ky=0.10,
        sigma_kz=0.35,
        k0=8.0,
        phi0=0.7,
        slit_separation_a=4.0,
        screen_half_width_b=3.0,
        z0=-8.0,
        z_screen=22.0,
    )


def make_k_grids(params: KGDoublePacket3DParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kx = np.linspace(-params.k_extent_sigma * params.sigma_kx, params.k_extent_sigma * params.sigma_kx, params.nkx)
    ky = np.linspace(-params.k_extent_sigma * params.sigma_ky, params.k_extent_sigma * params.sigma_ky, params.nky)
    kz = np.linspace(
        params.k0 - params.k_extent_sigma * params.sigma_kz,
        params.k0 + params.k_extent_sigma * params.sigma_kz,
        params.nkz,
    )
    return kx, ky, kz


def omega_3d(kx: np.ndarray, ky: np.ndarray, kz: np.ndarray, m: float) -> np.ndarray:
    return np.sqrt(kx**2 + ky**2 + kz**2 + m**2)


def spectral_amplitude_3d(
    kx: np.ndarray, ky: np.ndarray, kz: np.ndarray, params: KGDoublePacket3DParams
) -> np.ndarray:
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    gaussian = np.exp(
        -(KX**2) / (4.0 * params.sigma_kx**2)
        -(KY**2) / (4.0 * params.sigma_ky**2)
        -((KZ - params.k0) ** 2) / (4.0 * params.sigma_kz**2)
    )

    x_l0 = -0.5 * params.slit_separation_a
    x_r0 = +0.5 * params.slit_separation_a
    y0 = 0.0

    a_l = params.alpha * gaussian * np.exp(-1j * (KX * x_l0 + KY * y0 + KZ * params.z0))
    a_r = (
        np.sqrt(max(0.0, 1.0 - params.alpha))
        * gaussian
        * np.exp(-1j * (KX * x_r0 + KY * y0 + KZ * params.z0) + 1j * params.phi0)
    )
    return a_l + a_r


def _triple_trapezoid(values: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray) -> np.ndarray:
    tmp = np.trapezoid(values, kz, axis=-1)
    tmp = np.trapezoid(tmp, ky, axis=-1)
    tmp = np.trapezoid(tmp, kx, axis=-1)
    return tmp


def compute_screen_line(
    x: np.ndarray,
    y0: float,
    z_screen: float,
    t: float,
    params: KGDoublePacket3DParams,
) -> dict[str, np.ndarray]:
    kx, ky, kz = make_k_grids(params)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    omega = omega_3d(KX, KY, KZ, params.m)
    amp = spectral_amplitude_3d(kx, ky, kz, params)

    phase = np.exp(
        1j
        * (
            x[:, None, None, None] * KX[None, ...]
            + y0 * KY[None, ...]
            + z_screen * KZ[None, ...]
        )
        - 1j * t * omega[None, ...]
    )
    base = amp[None, ...] * phase

    psi = _triple_trapezoid(base, kx, ky, kz)
    dt_psi = _triple_trapezoid((-1j * omega)[None, ...] * base, kx, ky, kz)
    dx_psi = _triple_trapezoid((1j * KX)[None, ...] * base, kx, ky, kz)
    dy_psi = _triple_trapezoid((1j * KY)[None, ...] * base, kx, ky, kz)
    dz_psi = _triple_trapezoid((1j * KZ)[None, ...] * base, kx, ky, kz)
    dt2_psi = _triple_trapezoid((-(omega**2))[None, ...] * base, kx, ky, kz)
    dx2_psi = _triple_trapezoid((-(KX**2))[None, ...] * base, kx, ky, kz)
    dy2_psi = _triple_trapezoid((-(KY**2))[None, ...] * base, kx, ky, kz)
    dz2_psi = _triple_trapezoid((-(KZ**2))[None, ...] * base, kx, ky, kz)

    rho = np.abs(psi) ** 2
    rho_safe = np.maximum(rho, 1e-14)
    j0 = np.imag(np.conj(psi) * dt_psi)
    jx = -np.imag(np.conj(psi) * dx_psi)
    jy = -np.imag(np.conj(psi) * dy_psi)
    jz = -np.imag(np.conj(psi) * dz_psi)
    dtS = j0 / rho_safe
    dxS = -jx / rho_safe
    dyS = -jy / rho_safe
    dzS = -jz / rho_safe
    X_field = dtS**2 - dxS**2 - dyS**2 - dzS**2

    R = np.sqrt(rho_safe)
    dt_rho = np.conj(psi) * dt_psi + psi * np.conj(dt_psi)
    dx_rho = np.conj(psi) * dx_psi + psi * np.conj(dx_psi)
    dy_rho = np.conj(psi) * dy_psi + psi * np.conj(dy_psi)
    dz_rho = np.conj(psi) * dz_psi + psi * np.conj(dz_psi)
    dt2_rho = np.conj(psi) * dt2_psi + psi * np.conj(dt2_psi) + 2.0 * np.conj(dt_psi) * dt_psi
    dx2_rho = np.conj(psi) * dx2_psi + psi * np.conj(dx2_psi) + 2.0 * np.conj(dx_psi) * dx_psi
    dy2_rho = np.conj(psi) * dy2_psi + psi * np.conj(dy2_psi) + 2.0 * np.conj(dy_psi) * dy_psi
    dz2_rho = np.conj(psi) * dz2_psi + psi * np.conj(dz2_psi) + 2.0 * np.conj(dz_psi) * dz_psi

    dt2R = dt2_rho / (2.0 * R) - (dt_rho**2) / (4.0 * R**3)
    dx2R = dx2_rho / (2.0 * R) - (dx_rho**2) / (4.0 * R**3)
    dy2R = dy2_rho / (2.0 * R) - (dy_rho**2) / (4.0 * R**3)
    dz2R = dz2_rho / (2.0 * R) - (dz_rho**2) / (4.0 * R**3)
    Q = (dt2R - dx2R - dy2R - dz2R) / R

    return {
        "x": x,
        "y0": np.array([y0]),
        "psi": psi,
        "rho": rho,
        "j0": j0,
        "jx": jx,
        "jy": jy,
        "jz": jz,
        "dtS": dtS,
        "dxS": dxS,
        "dyS": dyS,
        "dzS": dzS,
        "X": X_field,
        "Q": np.real_if_close(Q),
    }


def count_visible_fringes(x: np.ndarray, rho: np.ndarray, threshold_fraction: float = 0.1) -> dict[str, object]:
    mx = float(np.max(rho))
    peaks: list[float] = []
    for i in range(1, len(rho) - 1):
        if rho[i] > rho[i - 1] and rho[i] >= rho[i + 1] and rho[i] > threshold_fraction * mx:
            peaks.append(float(x[i]))
    return {
        "threshold_fraction": threshold_fraction,
        "count": len(peaks),
        "positions": peaks,
    }


def run_screen_case(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = default_params_3d()
    x = params.x_screen
    y0 = 0.0
    z_screen = params.z_screen
    t = float((z_screen - params.z0) / max(params.v_group, 1e-8))

    data = compute_screen_line(x, y0, z_screen, t, params)
    rho = data["rho"]
    X_field = np.real(data["X"])
    Q_field = np.real(data["Q"])
    fringes = count_visible_fringes(x, rho, threshold_fraction=0.1)

    payload = {
        "params": {
            "m": params.m,
            "alpha": params.alpha,
            "sigma_kx": params.sigma_kx,
            "sigma_ky": params.sigma_ky,
            "sigma_kz": params.sigma_kz,
            "sigma_kx_over_k0": params.sigma_kx / params.k0,
            "sigma_ky_over_k0": params.sigma_ky / params.k0,
            "sigma_kz_over_k0": params.sigma_kz / params.k0,
            "k0": params.k0,
            "omega0": params.omega0,
            "phi0": params.phi0,
            "slit_separation_a": params.slit_separation_a,
            "screen_half_width_b": params.screen_half_width_b,
            "z0": params.z0,
            "y0": y0,
            "z_screen": z_screen,
            "t": t,
            "x_min": float(x.min()),
            "x_max": float(x.max()),
        },
        "summary": {
            "rho_max_screen_line": float(np.max(rho)),
            "rho_min_screen_line": float(np.min(rho)),
            "visibility_like": float((np.max(rho) - np.min(rho)) / (np.max(rho) + np.min(rho) + 1e-14)),
            "X_min_screen_line": float(np.min(X_field)),
            "X_max_screen_line": float(np.max(X_field)),
            "Q_abs_max_screen_line": float(np.max(np.abs(Q_field))),
            "eps_abs_max_screen_line": float(np.max(np.abs(Q_field)) / (params.m**2)),
            "visible_fringe_count_estimate": fringes["count"],
            "visible_fringe_positions": fringes["positions"],
        },
    }

    (out / "screen3d_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        out / "screen3d_fields.npz",
        x=x,
        y0=np.array([y0]),
        rho=rho,
        X=X_field,
        Q=Q_field,
        psi=data["psi"],
    )
    return payload


if __name__ == "__main__":
    result = run_screen_case(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
