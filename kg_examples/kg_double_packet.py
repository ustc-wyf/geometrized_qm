from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json

import numpy as np


@dataclass
class KGDoublePacketParams:
    m: float
    alpha: float
    sigma_l: float
    sigma_r: float
    k0: float
    phi0: float
    x_l0: float
    x_r0: float
    k_extent_sigma: float = 8.0
    nk: int = 8192

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))


def make_k_grid(params: KGDoublePacketParams) -> np.ndarray:
    sigma_max = max(params.sigma_l, params.sigma_r)
    width = params.k_extent_sigma * sigma_max
    return np.linspace(params.k0 - width, params.k0 + width, params.nk)


def omega_of_k(k: np.ndarray, m: float) -> np.ndarray:
    return np.sqrt(k**2 + m**2)


def spectral_amplitude(k: np.ndarray, params: KGDoublePacketParams) -> np.ndarray:
    a_l = (
        params.alpha
        * np.exp(-((k - params.k0) ** 2) / (4.0 * params.sigma_l**2))
        * np.exp(-1j * k * params.x_l0)
    )
    a_r = (
        np.sqrt(max(0.0, 1.0 - params.alpha))
        * np.exp(-((k - params.k0) ** 2) / (4.0 * params.sigma_r**2))
        * np.exp(-1j * k * params.x_r0 + 1j * params.phi0)
    )
    return a_l + a_r


def _integrate_over_k(values: np.ndarray, k: np.ndarray) -> np.ndarray:
    return np.trapezoid(values, k, axis=-1)


def compute_field_and_derivatives(
    x: np.ndarray, t: float, params: KGDoublePacketParams
) -> dict[str, np.ndarray]:
    k = make_k_grid(params)
    omega = omega_of_k(k, params.m)
    amp = spectral_amplitude(k, params)

    phase = np.exp(1j * np.outer(x, k) - 1j * t * omega[None, :])
    base = amp[None, :] * phase

    psi = _integrate_over_k(base, k)
    dt_psi = _integrate_over_k((-1j * omega)[None, :] * base, k)
    dx_psi = _integrate_over_k((1j * k)[None, :] * base, k)
    dt2_psi = _integrate_over_k((-(omega**2))[None, :] * base, k)
    dx2_psi = _integrate_over_k((-(k**2))[None, :] * base, k)

    return {
        "x": x,
        "k": k,
        "psi": psi,
        "dt_psi": dt_psi,
        "dx_psi": dx_psi,
        "dt2_psi": dt2_psi,
        "dx2_psi": dx2_psi,
    }


def compute_bohm_quantities(field: dict[str, np.ndarray], rho_floor: float = 1e-12) -> dict[str, np.ndarray]:
    psi = field["psi"]
    dt_psi = field["dt_psi"]
    dx_psi = field["dx_psi"]
    dt2_psi = field["dt2_psi"]
    dx2_psi = field["dx2_psi"]

    rho = np.abs(psi) ** 2
    rho_safe = np.maximum(rho, rho_floor)
    R = np.sqrt(rho_safe)

    dt_rho = np.conj(psi) * dt_psi + psi * np.conj(dt_psi)
    dx_rho = np.conj(psi) * dx_psi + psi * np.conj(dx_psi)
    dt2_rho = (
        np.conj(psi) * dt2_psi
        + psi * np.conj(dt2_psi)
        + 2.0 * np.conj(dt_psi) * dt_psi
    )
    dx2_rho = (
        np.conj(psi) * dx2_psi
        + psi * np.conj(dx2_psi)
        + 2.0 * np.conj(dx_psi) * dx_psi
    )

    dtR = dt_rho / (2.0 * R)
    dxR = dx_rho / (2.0 * R)
    dt2R = dt2_rho / (2.0 * R) - (dt_rho**2) / (4.0 * R**3)
    dx2R = dx2_rho / (2.0 * R) - (dx_rho**2) / (4.0 * R**3)

    Q = (dt2R - dx2R) / R

    j0 = np.imag(np.conj(psi) * dt_psi)
    j1 = -np.imag(np.conj(psi) * dx_psi)

    dtS = j0 / rho_safe
    dxS = -j1 / rho_safe
    X = dtS**2 - dxS**2
    S = np.unwrap(np.angle(psi))

    return {
        "rho": rho,
        "R": R,
        "S": S,
        "j0": j0,
        "j1": j1,
        "dtS": dtS,
        "dxS": dxS,
        "X": X,
        "Q": np.real_if_close(Q),
        "dtR": np.real_if_close(dtR),
        "dxR": np.real_if_close(dxR),
    }


def summarize_q(x: np.ndarray, rho: np.ndarray, Q: np.ndarray, m: float, rho_quantile: float = 0.2) -> dict[str, float]:
    support_cut = np.quantile(rho, rho_quantile)
    mask = rho >= support_cut
    q_support = np.real(Q[mask])
    eps_support = q_support / (m**2)

    idx_abs = int(np.argmax(np.abs(Q)))
    idx_support = int(np.argmax(np.abs(q_support))) if q_support.size else 0

    result = {
        "Q_abs_max_all": float(np.max(np.abs(Q))),
        "Q_abs_max_support": float(np.max(np.abs(q_support))) if q_support.size else 0.0,
        "Q_mean_support": float(np.mean(q_support)) if q_support.size else 0.0,
        "Q_std_support": float(np.std(q_support)) if q_support.size else 0.0,
        "eps_abs_max_support": float(np.max(np.abs(eps_support))) if q_support.size else 0.0,
        "x_at_abs_max_all": float(x[idx_abs]),
        "support_threshold_rho": float(support_cut),
    }
    return result


def summarize_q_relative_support(
    x: np.ndarray, rho: np.ndarray, Q: np.ndarray, m: float, fractions: tuple[float, ...] = (1e-2, 1e-3, 1e-4, 1e-6)
) -> list[dict[str, float]]:
    rho_max = float(np.max(rho))
    summaries: list[dict[str, float]] = []
    for frac in fractions:
        mask = rho >= rho_max * frac
        q = np.real(Q[mask])
        xx = x[mask]
        summaries.append(
            {
                "relative_density_cut": float(frac),
                "count": int(mask.sum()),
                "Q_abs_max": float(np.max(np.abs(q))) if q.size else 0.0,
                "Q_mean": float(np.mean(q)) if q.size else 0.0,
                "Q_std": float(np.std(q)) if q.size else 0.0,
                "eps_abs_max": float(np.max(np.abs(q)) / (m**2)) if q.size else 0.0,
                "x_at_abs_max": float(xx[np.argmax(np.abs(q))]) if q.size else 0.0,
            }
        )
    return summaries


def default_baseline_params() -> KGDoublePacketParams:
    return KGDoublePacketParams(
        m=1.0,
        alpha=0.6,
        sigma_l=0.35,
        sigma_r=0.28,
        k0=2.0,
        phi0=0.7,
        x_l0=-6.0,
        x_r0=6.0,
        nk=8192,
    )


def run_baseline_case(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = default_baseline_params()
    x = np.linspace(-40.0, 40.0, 4001)

    omega0 = params.omega0
    v_group = params.k0 / omega0
    center_sep = params.x_r0 - params.x_l0
    t_overlap = center_sep / (4.0 * max(params.sigma_l, params.sigma_r) * max(v_group, 1e-8))
    t_sample = float(max(0.0, t_overlap))

    field = compute_field_and_derivatives(x, t_sample, params)
    bohm = compute_bohm_quantities(field)
    summary = summarize_q(x, bohm["rho"], bohm["Q"], params.m)
    relative_support = summarize_q_relative_support(x, bohm["rho"], bohm["Q"], params.m)

    payload = {
        "params": {
            "m": params.m,
            "alpha": params.alpha,
            "sigma_l": params.sigma_l,
            "sigma_r": params.sigma_r,
            "k0": params.k0,
            "omega0": params.omega0,
            "phi0": params.phi0,
            "x_l0": params.x_l0,
            "x_r0": params.x_r0,
            "t_sample": t_sample,
        },
        "summary": summary,
        "relative_support": relative_support,
    }

    (out / "baseline_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    np.savez_compressed(
        out / "baseline_fields.npz",
        x=x,
        psi=field["psi"],
        rho=bohm["rho"],
        S=bohm["S"],
        j0=bohm["j0"],
        j1=bohm["j1"],
        X=bohm["X"],
        Q=bohm["Q"],
    )

    return payload


if __name__ == "__main__":
    result = run_baseline_case(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
