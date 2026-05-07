from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def second_derivative_periodic(arr: np.ndarray, dx: float) -> np.ndarray:
    k = 2.0 * np.pi * np.fft.fftfreq(arr.size, d=dx)
    return np.fft.ifft(-(k**2) * np.fft.fft(arr))


def q_transverse(A: np.ndarray, dx: float, rho_floor: float = 1e-10) -> np.ndarray:
    rho = np.abs(A) ** 2
    R = np.sqrt(np.maximum(rho, rho_floor))
    d2R = second_derivative_periodic(R, dx)
    return np.real_if_close(d2R / R)


def r_surr_from_A(A: np.ndarray, dx: float, m: float) -> np.ndarray:
    Qx = q_transverse(A, dx)
    eps = Qx / (m**2)
    return np.real_if_close(second_derivative_periodic(eps, dx))


def evolve_branch(
    A0: np.ndarray,
    x: np.ndarray,
    z_steps: int,
    dz: float,
    k0: float,
    m: float,
    branch: str,
    lam: float,
    ell: float,
) -> np.ndarray:
    A = A0.astype(complex).copy()
    dx = float(x[1] - x[0])
    kx = 2.0 * np.pi * np.fft.fftfreq(x.size, d=dx)
    kinetic_phase = np.exp(-1j * (kx**2) * dz / (2.0 * k0))

    def potential(current_A: np.ndarray):
        if branch == "g":
            return 0.0
        R_surr = r_surr_from_A(current_A, dx, m)
        if branch == "gtilde":
            return lam * R_surr
        if branch == "fR":
            chi = (ell**2) * R_surr
            F = 1.0 / np.sqrt(1.0 + chi**2)
            return lam * F * R_surr
        raise ValueError(branch)

    for _ in range(z_steps):
        V1 = potential(A)
        if np.isscalar(V1):
            A *= np.exp(-1j * V1 * dz / 2.0)
        else:
            A *= np.exp(-1j * V1 * dz / 2.0)
        A = np.fft.ifft(np.fft.fft(A) * kinetic_phase)
        V2 = potential(A)
        if np.isscalar(V2):
            A *= np.exp(-1j * V2 * dz / 2.0)
        else:
            A *= np.exp(-1j * V2 * dz / 2.0)
    return A


def initial_envelope(x: np.ndarray, alpha: float, a: float, w0: float, phi0: float) -> np.ndarray:
    left = alpha * np.exp(-((x + a / 2.0) ** 2) / (2.0 * w0**2))
    right = np.sqrt(max(0.0, 1.0 - alpha)) * np.exp(-((x - a / 2.0) ** 2) / (2.0 * w0**2)) * np.exp(1j * phi0)
    return left + right


def run_comparison(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    m = 1.0
    k0 = 10.0
    alpha = 0.5
    a = 8.0
    w0 = 1.0
    phi0 = 0.0
    x = np.linspace(-5.0, 5.0, 1201)
    z_final = 20.0
    z_steps = 400
    dz = z_final / z_steps
    lam = 0.02
    ell = 0.15

    A0 = initial_envelope(x, alpha, a, w0, phi0)

    outputs = {}
    for branch in ("g", "gtilde", "fR"):
        A = evolve_branch(A0, x, z_steps, dz, k0, m, branch, lam, ell)
        rho = np.abs(A) ** 2
        outputs[branch] = {"rho": rho}
        np.savez_compressed(out / f"{branch}_evolution_fR.npz", x=x, rho=rho)

    rho_g = outputs["g"]["rho"]
    payload = {
        "params": {
            "m": m,
            "k0": k0,
            "alpha": alpha,
            "a": a,
            "w0": w0,
            "phi0": phi0,
            "z_final": z_final,
            "z_steps": z_steps,
            "lam": lam,
            "ell": ell,
        },
        "summary": {
            "L1_g_vs_gtilde": float(np.trapezoid(np.abs(rho_g - outputs["gtilde"]["rho"]), x)),
            "L1_g_vs_fR": float(np.trapezoid(np.abs(rho_g - outputs["fR"]["rho"]), x)),
            "L1_gtilde_vs_fR": float(np.trapezoid(np.abs(outputs["gtilde"]["rho"] - outputs["fR"]["rho"]), x)),
        },
    }
    (out / "summary_fR.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload


if __name__ == "__main__":
    result = run_comparison(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
