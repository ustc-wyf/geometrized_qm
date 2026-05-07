from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz, render_map


def ddx(f, dx):
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def ddz(f, dz):
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out


def centerline(path: Path, x: np.ndarray, curves):
    plt.figure(figsize=(6.2, 3.8))
    for label, y in curves:
        plt.plot(x, y, label=label, linewidth=1.5)
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run_scheme(outdir: Path, use_jrho: bool):
    params = Exact2p1Params(x_half_range=20, z_half_range=20, nx=181, nz=181, nkx=61, nkz=61)
    x, z, psi, rho = integrate_xz(0.5, 0.0, params)
    s = np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dt = 0.01
    steps = 200
    exact_s = -params.m * dt * steps

    for _ in range(steps):
        if use_jrho:
            psix = np.gradient(psi, dx, axis=0)
            psiz = np.gradient(psi, dz, axis=1)
            rho_safe = np.maximum(rho, 1e-9)
            sx = np.imag(np.conj(psi) * psix) / rho_safe
            sz = np.imag(np.conj(psi) * psiz) / rho_safe
        else:
            sx = ddx(s, dx)
            sz = ddz(s, dz)
        e = np.sqrt(np.maximum(sx * sx + sz * sz + params.m**2, 1e-12))
        n = np.maximum(rho * e, 1e-12)
        fluxx = (n / e) * sx
        fluxz = (n / e) * sz
        n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
        s_t = -e
        s = s + dt * s_t
        n = n + dt * n_t
        rho = np.maximum(n / e, 1e-12)

    render_map(outdir / ("rho_jrho.png" if use_jrho else "rho_fd.png"), x, z, rho, "rho")
    centerline(outdir / ("center_jrho.png" if use_jrho else "center_fd.png"), x, [("rho", rho[:, len(z)//2])])
    return {
        "rho_deviation_linf": float(np.max(np.abs(rho - np.abs(psi)**2))),
        "s_center": float(s[len(x)//2, len(z)//2]),
        "s_expected": float(exact_s),
    }


def main():
    out = Path(__file__).resolve().parent.parent / "visualizations" / "debug_b_phase_gradients"
    out.mkdir(parents=True, exist_ok=True)
    summary = {
        "fd_phase_gradient": run_scheme(out, False),
        "j_over_rho_gradient": run_scheme(out, True),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
