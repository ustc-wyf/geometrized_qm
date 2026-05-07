from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from exact_kg_2p1_benchmark import Exact2p1Params, evaluate_points, integrate_xz, render_map


def lap(f, dx, dz):
    out = np.zeros_like(f)
    out[1:-1] += (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dx * dx)
    out[:, 1:-1] += (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / (dz * dz)
    out[0] = out[1]
    out[-1] = out[-2]
    out[:, 0] = out[:, 1]
    out[:, -1] = out[:, -2]
    return out


def apply_boundary(field, target):
    field[0, :] = target[0, :]
    field[-1, :] = target[-1, :]
    field[:, 0] = target[:, 0]
    field[:, -1] = target[:, -1]


def boundary_exact(alpha, t, x, z, params):
    nx = len(x)
    nz = len(z)
    psi = np.zeros((nx, nz), dtype=np.complex128)

    top_x = x
    top_z = np.full_like(x, z[-1])
    bot_x = x
    bot_z = np.full_like(x, z[0])
    left_x = np.full_like(z, x[0])
    left_z = z
    right_x = np.full_like(z, x[-1])
    right_z = z

    psi[:, -1] = evaluate_points(alpha, top_x, top_z, t, params)
    psi[:, 0] = evaluate_points(alpha, bot_x, bot_z, t, params)
    psi[0, :] = evaluate_points(alpha, left_x, left_z, t, params)
    psi[-1, :] = evaluate_points(alpha, right_x, right_z, t, params)
    return psi


def solve_once(params: Exact2p1Params, alpha: float, dt: float):
    t_final = params.overlap_time
    steps = int(round(t_final / dt))
    dt = t_final / steps

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    for n in range(steps):
        t_next = (n + 1) * dt
        psi_target = boundary_exact(alpha, t_next, x, z, params)
        psi_next = np.zeros_like(psi_now)
        psi_next[1:-1, 1:-1] = (
            2.0 * psi_now[1:-1, 1:-1]
            - psi_prev[1:-1, 1:-1]
            + dt * dt * (lap(psi_now, dx, dz)[1:-1, 1:-1] - params.m**2 * psi_now[1:-1, 1:-1])
        )
        apply_boundary(psi_next, psi_target)
        psi_prev, psi_now = psi_now, psi_next

    _, _, psi_exact, rho_exact = integrate_xz(alpha, t_final, params)
    rho_num = np.abs(psi_now) ** 2
    diff = rho_num - rho_exact
    summary = {
        "dt": dt,
        "steps": steps,
        "l1_mean": float(np.mean(np.abs(diff))),
        "linf": float(np.max(np.abs(diff))),
        "relative_l1": float(np.mean(np.abs(diff)) / np.max(rho_exact)),
        "relative_linf": float(np.max(np.abs(diff)) / np.max(rho_exact)),
    }
    return x, z, rho_num, rho_exact, diff, summary


def run_calibration(outdir: str | Path):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = Exact2p1Params(nx=61, nz=61, nkx=41, nkz=41)
    dts = [0.2, 0.1, 0.05]
    runs = []
    best = None
    best_err = None

    for dt in dts:
        x, z, rho_num, rho_exact, diff, summary = solve_once(params, 0.5, dt)
        runs.append(summary)
        if best_err is None or summary["relative_l1"] < best_err:
            best = (x, z, rho_num, rho_exact, diff, summary)
            best_err = summary["relative_l1"]

    center = params.nz // 2
    x, z, rho_num, rho_exact, diff, best_summary = best

    render_map(out / "rho_exact.png", x, z, rho_exact, "exact KG benchmark")
    render_map(out / "rho_numeric.png", x, z, rho_num, "numeric KG solution")
    render_map(out / "rho_abs_diff.png", x, z, np.abs(diff), "|rho_num-rho_exact|")

    plt.figure(figsize=(6.2, 3.8))
    plt.plot(x, rho_exact[:, center], label="exact", linewidth=1.8)
    plt.plot(x, rho_num[:, center], label="numeric", linewidth=1.4)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "centerline_compare.png", dpi=180)
    plt.close()

    summary = {
        "params": {
            "extent": params.x_half_range,
            "nx": params.nx,
            "nz": params.nz,
            "t_final": params.overlap_time,
            "k0": params.k0,
            "sigma_k_parallel": params.sigma_k_parallel,
            "sigma_k_perp": params.sigma_k_perp,
        },
        "runs": runs,
        "best_run": best_summary,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    summary = run_calibration(Path(__file__).resolve().parent.parent / "visualizations" / "a_solver_calibration")
    print(json.dumps(summary, indent=2))
