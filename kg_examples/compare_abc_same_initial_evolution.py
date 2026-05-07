from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from calibrate_exact_kg_ps_fd_absorbing import (
    apply_boundary,
    boundary_exact,
    make_sponge,
    spectral_lap,
)
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz
from simulate_bc_psfd_imex import (
    apply_geometry_boundary,
    edge_damp,
    rhs_branch_b,
    rhs_branch_c,
    spectral_k2,
    wave_imex_step,
)


def unwrap_phase(psi: np.ndarray) -> np.ndarray:
    return np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)


def crop_inner(x: np.ndarray, z: np.ndarray, inner_half_range: float):
    ix = np.where(np.abs(x) <= inner_half_range + 1e-12)[0]
    iz = np.where(np.abs(z) <= inner_half_range + 1e-12)[0]
    return ix, iz, x[ix], z[iz]


def rel_l1(a: np.ndarray, b: np.ndarray, scale: float) -> float:
    return float(np.mean(np.abs(a - b)) / max(scale, 1e-12))


def render_montage(path: Path, x: np.ndarray, z: np.ndarray, snapshots: list[dict], keys: list[str], titles: list[str], cmap: str = "viridis") -> None:
    nrows = len(snapshots)
    ncols = len(keys)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.0 * nrows), squeeze=False)

    vmax = max(float(np.max(s[key])) for s in snapshots for key in keys)
    vmax = max(vmax, 1e-12)
    for i, snap in enumerate(snapshots):
        time_label = f"t/T={snap['time_fraction']:.2f}"
        for j, (key, title) in enumerate(zip(keys, titles)):
            ax = axes[i][j]
            field = snap[key]
            im = ax.imshow(
                field.T,
                origin="lower",
                extent=[x.min(), x.max(), z.min(), z.max()],
                aspect="equal",
                cmap=cmap,
                vmin=0.0,
                vmax=vmax,
            )
            if i == 0:
                ax.set_title(title)
            if j == 0:
                ax.set_ylabel(f"{time_label}\nz")
            else:
                ax.set_ylabel("z")
            ax.set_xlabel("x")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label="rho")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def render_diff_montage(path: Path, x: np.ndarray, z: np.ndarray, snapshots: list[dict], keys: list[str], titles: list[str]) -> None:
    nrows = len(snapshots)
    ncols = len(keys)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.0 * nrows), squeeze=False)

    vmax = max(float(np.max(s[key])) for s in snapshots for key in keys)
    vmax = max(vmax, 1e-12)
    for i, snap in enumerate(snapshots):
        time_label = f"t/T={snap['time_fraction']:.2f}"
        for j, (key, title) in enumerate(zip(keys, titles)):
            ax = axes[i][j]
            field = snap[key]
            im = ax.imshow(
                field.T,
                origin="lower",
                extent=[x.min(), x.max(), z.min(), z.max()],
                aspect="equal",
                cmap="magma",
                vmin=0.0,
                vmax=vmax,
            )
            if i == 0:
                ax.set_title(title)
            if j == 0:
                ax.set_ylabel(f"{time_label}\nz")
            else:
                ax.set_ylabel("z")
            ax.set_xlabel("x")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, label="abs diff")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def render_error_curves(path: Path, ts: np.ndarray, ba: np.ndarray, ca: np.ndarray, cb: np.ndarray) -> None:
    plt.figure(figsize=(6.6, 4.2))
    plt.semilogy(ts, ba, label="B vs A", linewidth=1.8)
    plt.semilogy(ts, ca, label="C vs A", linewidth=1.8)
    plt.semilogy(ts, cb, label="C vs B", linewidth=1.8)
    plt.xlabel("t / T_overlap")
    plt.ylabel("relative L1")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(
    outdir: str | Path,
    alpha: float = 0.5,
    lambda_grav: float = 3.0,
    mp: float = 300.0,
    ell: float = 0.02,
    dt_target: float = 0.025,
    c_substeps: int = 8,
    inner: float = 10.0,
    outer: float = 15.0,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 61,
    nkz: int = 61,
):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    solver_params = Exact2p1Params(
        alpha=alpha,
        x_half_range=outer,
        z_half_range=outer,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )
    benchmark_params = Exact2p1Params(alpha=alpha)

    dx = 2.0 * outer / (solver_params.nx - 1)
    dz = 2.0 * outer / (solver_params.nz - 1)
    dt = dt_target
    total_time = solver_params.overlap_time
    steps = int(round(total_time / dt))
    dt = total_time / steps
    k2 = spectral_k2(solver_params.nx, solver_params.nz, dx, dz)

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, solver_params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, solver_params)
    sponge_inner = max(inner + 2.0, 0.8 * outer)
    sigma = make_sponge(x, z, inner_half_range=sponge_inner, outer_half_range=outer, sigma_max=2.0)

    tau_b = np.zeros((solver_params.nx, solver_params.nz))
    beta_b = np.zeros_like(tau_b)
    ptau_b = np.zeros_like(tau_b)
    pbeta_b = np.zeros_like(tau_b)

    tau_c = np.zeros_like(tau_b)
    beta_c = np.zeros_like(tau_b)
    eta_c = np.zeros_like(tau_b)
    ptau_c = np.zeros_like(tau_b)
    pbeta_c = np.zeros_like(tau_b)
    peta_c = np.zeros_like(tau_b)

    ix, iz, xi, zi = crop_inner(x, z, inner)
    if not (np.allclose(xi, np.linspace(-inner, inner, benchmark_params.nx)) and np.allclose(zi, np.linspace(-inner, inner, benchmark_params.nz))):
        raise RuntimeError("Inner crop no longer matches approved benchmark grid.")

    snapshot_steps = [0, steps // 4, steps // 2, (3 * steps) // 4, steps]
    snapshots: list[dict] = []
    times = []
    err_ba = []
    err_ca = []
    err_cb = []

    def record(step_index: int, current_time: float):
        rho_a_outer = np.abs(psi_now) ** 2
        rho_b_outer = rho_a_outer * np.exp(-np.clip(beta_b + tau_b, -12.0, 12.0))
        rho_c_outer = rho_a_outer * np.exp(-np.clip(beta_c + tau_c, -12.0, 12.0))

        A = rho_a_outer[np.ix_(ix, iz)]
        B = rho_b_outer[np.ix_(ix, iz)]
        C = rho_c_outer[np.ix_(ix, iz)]
        scale = float(np.max(A))

        times.append(current_time / total_time)
        err_ba.append(rel_l1(B, A, scale))
        err_ca.append(rel_l1(C, A, scale))
        err_cb.append(rel_l1(C, B, scale))

        if step_index in snapshot_steps:
            snapshots.append(
                {
                    "step": step_index,
                    "time": current_time,
                    "time_fraction": current_time / total_time,
                    "A": A,
                    "B": B,
                    "C": C,
                    "BA": np.abs(B - A),
                    "CA": np.abs(C - A),
                    "CB": np.abs(C - B),
                }
            )

    record(0, 0.0)

    for n in range(steps):
        rho_a_outer = np.abs(psi_now) ** 2
        s_a_outer = unwrap_phase(psi_now)

        src_tau_b, src_beta_b, *_ = rhs_branch_b(
            tau_b, ptau_b, beta_b, pbeta_b, rho_a_outer, s_a_outer, solver_params.m, mp, dx, dz, k2
        )
        tau_b, ptau_b = wave_imex_step(tau_b, ptau_b, lambda_grav * src_tau_b, dt, k2)
        beta_b, pbeta_b = wave_imex_step(beta_b, pbeta_b, lambda_grav * src_beta_b, dt, k2)
        edge_damp(tau_b, sigma, dt, -2.0, 2.0)
        edge_damp(beta_b, sigma, dt, -2.0, 2.0)
        edge_damp(ptau_b, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_b, sigma, dt, -10.0, 10.0)
        apply_geometry_boundary(tau_b, beta_b)

        dt_c = dt / max(c_substeps, 1)
        sigma_c = sigma / max(c_substeps, 1)
        for _ in range(max(c_substeps, 1)):
            phi_c = 1.0 + eta_c
            src_tau_c, src_beta_c, src_phi_c, *_ = rhs_branch_c(
                tau_c, ptau_c, beta_c, pbeta_c, phi_c, peta_c, rho_a_outer, s_a_outer, solver_params.m, mp, ell, dx, dz
            )
            tau_c, ptau_c = wave_imex_step(tau_c, ptau_c, lambda_grav * src_tau_c, dt_c, k2)
            beta_c, pbeta_c = wave_imex_step(beta_c, pbeta_c, lambda_grav * src_beta_c, dt_c, k2)
            eta_c, peta_c = wave_imex_step(eta_c, peta_c, lambda_grav * src_phi_c, dt_c, k2)
            edge_damp(tau_c, sigma_c, dt_c, -2.0, 2.0)
            edge_damp(beta_c, sigma_c, dt_c, -2.0, 2.0)
            edge_damp(eta_c, sigma_c, dt_c, -0.8, 0.0)
            edge_damp(ptau_c, sigma_c, dt_c, -10.0, 10.0)
            edge_damp(pbeta_c, sigma_c, dt_c, -10.0, 10.0)
            edge_damp(peta_c, sigma_c, dt_c, -10.0, 10.0)
            phi_c = np.clip(1.0 + eta_c, 0.2, 1.0)
            eta_c = phi_c - 1.0
            apply_geometry_boundary(tau_c, beta_c, phi_c)
            eta_c[0, :] = 0.0
            eta_c[-1, :] = 0.0
            eta_c[:, 0] = 0.0
            eta_c[:, -1] = 0.0

        t_next = (n + 1) * dt
        psi_b = boundary_exact(alpha, t_next, x, z, solver_params)
        rhs = spectral_lap(psi_now, dx, dz) - solver_params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

        record(n + 1, t_next)

    _, _, _, A_exact_final = integrate_xz(alpha, total_time, benchmark_params)
    A_final = (np.abs(psi_now) ** 2)[np.ix_(ix, iz)]
    a_exact_final = rel_l1(A_final, A_exact_final, max(float(np.max(A_exact_final)), 1e-12))

    render_montage(out / "rho_montage.png", xi, zi, snapshots, ["A", "B", "C"], ["A", "B", "C"])
    render_diff_montage(out / "diff_montage.png", xi, zi, snapshots, ["BA", "CA", "CB"], ["|B-A|", "|C-A|", "|C-B|"])
    render_error_curves(out / "relL1_vs_time.png", np.array(times), np.array(err_ba), np.array(err_ca), np.array(err_cb))

    summary = {
        "alpha": alpha,
        "lambda_grav": lambda_grav,
        "mp": mp,
        "ell": ell,
        "dt": dt,
        "steps": steps,
        "c_substeps": c_substeps,
        "inner_half_range": inner,
        "outer_half_range": outer,
        "nx": nx,
        "nz": nz,
        "nkx": nkx,
        "nkz": nkz,
        "snapshot_times_over_T": [float(s["time_fraction"]) for s in snapshots],
        "final": {
            "B_vs_A_relL1": err_ba[-1],
            "C_vs_A_relL1": err_ca[-1],
            "C_vs_B_relL1": err_cb[-1],
            "A_num_vs_exact_relL1": a_exact_final,
            "B_tau_max": float(np.max(np.abs(tau_b))),
            "C_tau_max": float(np.max(np.abs(tau_c))),
            "C_phi_dev_max": float(np.max(np.abs(eta_c))),
        },
        "time_series": {
            "t_over_T": [float(v) for v in times],
            "B_vs_A_relL1": [float(v) for v in err_ba],
            "C_vs_A_relL1": [float(v) for v in err_ca],
            "C_vs_B_relL1": [float(v) for v in err_cb],
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary["final"], indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--lambda-grav", type=float, default=3.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--ell", type=float, default=0.02)
    parser.add_argument("--dt-target", type=float, default=0.025)
    parser.add_argument("--c-substeps", type=int, default=8)
    parser.add_argument("--inner", type=float, default=10.0)
    parser.add_argument("--outer", type=float, default=15.0)
    parser.add_argument("--nx", type=int, default=181)
    parser.add_argument("--nz", type=int, default=181)
    parser.add_argument("--nkx", type=int, default=61)
    parser.add_argument("--nkz", type=int, default=61)
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(Path(__file__).resolve().parent.parent / "visualizations" / "abc_same_initial_evolution"),
    )
    args = parser.parse_args()
    run(
        args.outdir,
        alpha=args.alpha,
        lambda_grav=args.lambda_grav,
        mp=args.mp,
        ell=args.ell,
        dt_target=args.dt_target,
        c_substeps=args.c_substeps,
        inner=args.inner,
        outer=args.outer,
        nx=args.nx,
        nz=args.nz,
        nkx=args.nkx,
        nkz=args.nkz,
    )
