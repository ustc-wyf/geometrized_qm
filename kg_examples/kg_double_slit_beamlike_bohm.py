from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

sys.path.append(str(Path(__file__).resolve().parent))

from kg_double_slit_beamlike import default_params, compute_xz_slice


def compute_bohm_slice(alpha: float) -> dict[str, object]:
    params = default_params()
    params.alpha = alpha
    t = 0.0
    data = compute_xz_slice(t, params)

    x = data["x"]
    z = data["z"]
    psi = data["psi"]
    rho = data["rho"]
    Xg = data["X"]
    Zg = data["Z"]

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    rho_safe = np.maximum(rho, 1e-14)
    R = np.sqrt(rho_safe)

    # approximate time derivative from monochromatic factor only
    dt_psi = -1j * params.omega0 * psi
    dt2_psi = -(params.omega0**2) * psi

    dx_psi = np.gradient(psi, dx, axis=0)
    dz_psi = np.gradient(psi, dz, axis=1)
    dx2_psi = np.gradient(dx_psi, dx, axis=0)
    dz2_psi = np.gradient(dz_psi, dz, axis=1)

    j0 = np.imag(np.conj(psi) * dt_psi)
    jx = -np.imag(np.conj(psi) * dx_psi)
    jz = -np.imag(np.conj(psi) * dz_psi)

    dtS = j0 / rho_safe
    dxS = -jx / rho_safe
    dzS = -jz / rho_safe
    X_field = dtS**2 - dxS**2 - dzS**2

    dt_rho = np.conj(psi) * dt_psi + psi * np.conj(dt_psi)
    dx_rho = np.conj(psi) * dx_psi + psi * np.conj(dx_psi)
    dz_rho = np.conj(psi) * dz_psi + psi * np.conj(dz_psi)
    dt2_rho = np.conj(psi) * dt2_psi + psi * np.conj(dt2_psi) + 2.0 * np.conj(dt_psi) * dt_psi
    dx2_rho = np.conj(psi) * dx2_psi + psi * np.conj(dx2_psi) + 2.0 * np.conj(dx_psi) * dx_psi
    dz2_rho = np.conj(psi) * dz2_psi + psi * np.conj(dz2_psi) + 2.0 * np.conj(dz_psi) * dz_psi

    dt2R = dt2_rho / (2.0 * R) - (dt_rho**2) / (4.0 * R**3)
    dx2R = dx2_rho / (2.0 * R) - (dx_rho**2) / (4.0 * R**3)
    dz2R = dz2_rho / (2.0 * R) - (dz_rho**2) / (4.0 * R**3)
    Q = (dt2R - dx2R - dz2R) / R

    mask_main = rho >= np.max(rho) * 1e-2
    mask_loose = rho >= np.max(rho) * 1e-6

    return {
        "params": {
            "alpha": alpha,
            "m": params.m,
            "k0": params.k0,
            "omega0": params.omega0,
            "w0": params.w0,
            "x_min": float(x.min()),
            "x_max": float(x.max()),
            "z_min": float(z.min()),
            "z_max": float(z.max()),
        },
        "summary": {
            "rho_max": float(np.max(rho)),
            "Q_abs_max_all": float(np.max(np.abs(Q))),
            "Q_abs_max_main": float(np.max(np.abs(Q[mask_main]))),
            "Q_abs_max_loose": float(np.max(np.abs(Q[mask_loose]))),
            "eps_abs_max_main": float(np.max(np.abs(Q[mask_main])) / (params.m**2)),
            "X_min_main": float(np.min(X_field[mask_main])),
            "X_max_main": float(np.max(X_field[mask_main])),
        },
        "fields": {
            "x": x,
            "z": z,
            "rho": rho,
            "Q": np.real_if_close(Q),
            "X": np.real_if_close(X_field),
        },
    }


def run_all(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    summary = {}
    for alpha in (0.0, 0.5, 1.0):
        result = compute_bohm_slice(alpha)
        summary[str(alpha)] = result["summary"]
        np.savez_compressed(
            out / f"beamlike_bohm_alpha_{alpha:.1f}.npz",
            x=result["fields"]["x"],
            z=result["fields"]["z"],
            rho=result["fields"]["rho"],
            Q=result["fields"]["Q"],
            X=result["fields"]["X"],
        )
    (out / "beamlike_bohm_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = run_all(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
