from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from adm_sources_abc import b_branch_sources_from_ns
from adm_ykilling_geometry import momentum_constraint_residual_ykilling
from simulate_ab_adm_synchronous import hamiltonian_residual
from simulate_bc_from_a_initial_data import (
    B_STATE_KEYS,
    advance_rk4,
    b_rhs_n,
    make_a_flat_initial_data,
    make_bc_initial_states,
)
from simulate_bc_geometry_from_a_reference import solve_linearized_conformal_omega
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    evolve_positive_frequency,
    initial_wavefunction,
    make_grid,
    packet_centers,
    phase_field,
    render_montage,
    spectral_omega,
)


def render_scalar_montage(
    x: np.ndarray,
    z: np.ndarray,
    frames: list[np.ndarray],
    times: np.ndarray,
    centers: list[dict[str, tuple[float, float]]],
    out_path: Path,
    title_prefix: str,
    label: str,
    cmap: str = "coolwarm",
) -> None:
    n = len(frames)
    fig, axes = plt.subplots(2, (n + 1) // 2, figsize=(3.8 * ((n + 1) // 2), 6.8), constrained_layout=True)
    axes = np.array(axes).reshape(-1)
    vmin = min(float(np.min(f)) for f in frames)
    vmax = max(float(np.max(f)) for f in frames)
    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    last_im = None
    for ax, field, t, center in zip(axes, frames, times, centers):
        last_im = ax.pcolormesh(Xg, Zg, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(f"{title_prefix}, t={t:.4f}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.plot(center["A"][0], center["A"][1], "wo", ms=3.5, mec="k", mew=0.4)
        ax.plot(center["B"][0], center["B"][1], "w^", ms=4.0, mec="k", mew=0.4)
        ax.set_aspect("equal")
    for ax in axes[len(frames) :]:
        ax.axis("off")
    if last_im is not None:
        fig.colorbar(last_im, ax=axes.tolist(), shrink=0.92, label=label)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_difference_curve(
    times: np.ndarray,
    rel_l1: np.ndarray,
    max_abs: np.ndarray,
    out_path: Path,
) -> None:
    fig, ax1 = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    ax1.plot(times, rel_l1, marker="o", lw=1.8, color="#005f73", label="Relative L1")
    ax1.set_xlabel("t")
    ax1.set_ylabel(r"$\|\rho_B-\rho_A\|_1/\|\rho_A\|_1$")
    ax1.grid(True, alpha=0.3)
    ax2 = ax1.twinx()
    ax2.plot(times, max_abs, marker="s", lw=1.5, color="#bb3e03", label="Max abs diff")
    ax2.set_ylabel(r"$\max|\rho_B-\rho_A|$")
    lines = ax1.get_lines() + ax2.get_lines()
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc="upper left")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out = root / "visualizations" / "ab_localized_reference_vs_b"
    out.mkdir(parents=True, exist_ok=True)

    loc = FlatLocalizedCrossingParams(
        m=1.0,
        alpha=0.5,
        k0=6.0,
        sigma_parallel=1.6,
        sigma_perp=1.2,
        phi0=0.0,
        x_a0=-6.0,
        z_a0=-6.0,
        x_b0=6.0,
        z_b0=-6.0,
        x_half_range=20.0,
        z_half_range=20.0,
        nx=256,
        nz=256,
        t_final=16.0,
        n_frames=7,
    )
    x, z, X, Z = make_grid(loc)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(X, Z, loc)
    rho0 = np.abs(psi0) ** 2 + loc.rho_floor
    s0 = phase_field(psi0)
    omega = spectral_omega(x, z, loc.m)

    mass = 1.0
    mp = 300.0
    dt = 2.5e-4
    steps = 200
    sample_steps = np.array([0, 40, 80, 120, 160, 200], dtype=int)
    sample_times = sample_steps * dt

    _, a_src0 = make_a_flat_initial_data(rho0, s0, mass, dx, dz)
    a_weak0 = solve_linearized_conformal_omega(a_src0["energy"], mp, dx, dz)
    omega_init = a_weak0["omega"]
    b_state, _ = make_bc_initial_states(
        rho0=rho0,
        s0=s0,
        omega_init=omega_init,
        mass=mass,
        mp=mp,
        ell=0.2,
        dx=dx,
        dz=dz,
        c_phi_mode="flat",
        elliptic_lambda=1.0,
    )
    params = {
        "dx": dx,
        "dz": dz,
        "m": mass,
        "mp": mp,
        "ell": 0.2,
        "alpha": loc.alpha,
        "k0": loc.k0,
        "w0": np.nan,
        "phi0": loc.phi0,
        "theta": np.pi / 4.0,
        "x_half_range": loc.x_half_range,
        "z_half_range": loc.z_half_range,
        "nx": loc.nx,
        "nz": loc.nz,
        "initial_mode": "localized",
        "localized_psi0": psi0,
        "localized_omega": omega,
    }

    b_last = b_rhs_n(b_state, params)
    rho_a_frames: list[np.ndarray] = []
    rho_b_frames: list[np.ndarray] = []
    hxx_b_frames: list[np.ndarray] = []
    diff_frames: list[np.ndarray] = []
    centers: list[dict[str, tuple[float, float]]] = []
    rel_l1_history: list[float] = []
    max_abs_history: list[float] = []

    current_step = 0
    for target_step in sample_steps:
        while current_step < target_step:
            t_next = (current_step + 1) * dt
            b_last = advance_rk4(b_state, B_STATE_KEYS, b_rhs_n, dt, params, t_next)
            current_step += 1
        t_now = target_step * dt
        psi_a = evolve_positive_frequency(psi0, omega, float(t_now))
        rho_a = np.abs(psi_a) ** 2
        rho_b = b_last["rho_current"].copy()
        diff = rho_b - rho_a
        rho_a_frames.append(rho_a)
        rho_b_frames.append(rho_b)
        hxx_b_frames.append(b_state["h_xx"] - 1.0)
        diff_frames.append(diff)
        centers.append(packet_centers(loc, float(t_now)))
        rel_l1_history.append(float(np.sum(np.abs(diff)) / max(np.sum(np.abs(rho_a)), 1.0e-30)))
        max_abs_history.append(float(np.max(np.abs(diff))))

    render_montage(x, z, rho_b_frames, sample_times, centers, out / "rho_b_montage.png", log_scale=False)
    render_montage(x, z, rho_b_frames, sample_times, centers, out / "rho_b_log10_montage.png", log_scale=True)
    render_scalar_montage(
        x,
        z,
        diff_frames,
        sample_times,
        centers,
        out / "rho_b_minus_a_montage.png",
        title_prefix=r"$\rho_B-\rho_A$",
        label=r"$\rho_B-\rho_A$",
        cmap="coolwarm",
    )
    render_scalar_montage(
        x,
        z,
        hxx_b_frames,
        sample_times,
        centers,
        out / "hxx_b_minus_1_montage.png",
        title_prefix=r"$h^{(B)}_{xx}-1$",
        label=r"$h_{xx}^{(B)}-1$",
        cmap="coolwarm",
    )
    render_difference_curve(
        sample_times,
        np.asarray(rel_l1_history),
        np.asarray(max_abs_history),
        out / "ab_difference_vs_time.png",
    )

    b_src_full = b_branch_sources_from_ns(
        rho=b_last["rho_current"],
        s=b_state["s"],
        lapse=np.ones_like(b_last["rho_current"]),
        shift_x=np.zeros_like(b_last["rho_current"]),
        shift_z=np.zeros_like(b_last["rho_current"]),
        h_xx=b_state["h_xx"],
        h_xz=b_state["h_xz"],
        h_zz=b_state["h_zz"],
        beta=b_state["beta"],
        mass=mass,
        dx=dx,
        dz=dz,
    )
    mom_x_res, mom_z_res = momentum_constraint_residual_ykilling(
        b_state["h_xx"],
        b_state["h_xz"],
        b_state["h_zz"],
        b_state["beta"],
        b_state["k_xx"],
        b_state["k_xz"],
        b_state["k_zz"],
        b_state["k_beta"],
        b_src_full["mom_x"],
        b_src_full["mom_z"],
        mp,
        dx,
        dz,
    )
    ham_res = hamiltonian_residual(
        b_state["h_xx"],
        b_state["h_xz"],
        b_state["h_zz"],
        b_state["beta"],
        b_state["k_xx"],
        b_state["k_xz"],
        b_state["k_zz"],
        b_state["k_beta"],
        b_last["energy_current"],
        mp,
        dx,
        dz,
    )

    summary = {
        "localized_params": {
            "m": loc.m,
            "alpha": loc.alpha,
            "k0": loc.k0,
            "sigma_parallel": loc.sigma_parallel,
            "sigma_perp": loc.sigma_perp,
            "phi0": loc.phi0,
            "x_a0": loc.x_a0,
            "z_a0": loc.z_a0,
            "x_b0": loc.x_b0,
            "z_b0": loc.z_b0,
            "x_half_range": loc.x_half_range,
            "z_half_range": loc.z_half_range,
            "nx": loc.nx,
            "nz": loc.nz,
        },
        "evolution": {
            "steps": steps,
            "dt": dt,
            "t_final": float(steps * dt),
            "sample_steps": sample_steps.tolist(),
            "sample_times": sample_times.tolist(),
            "mp": mp,
        },
        "ab_difference": {
            "sample_rel_l1": rel_l1_history,
            "sample_max_abs_rho_diff": max_abs_history,
            "final_rel_l1": rel_l1_history[-1],
            "final_max_abs_rho_diff": max_abs_history[-1],
            "peak_rel_l1_over_samples": float(np.max(rel_l1_history)),
            "peak_max_abs_rho_diff_over_samples": float(np.max(max_abs_history)),
        },
        "b_geometry": {
            "final_h_xx_minus_1_max": float(np.max(np.abs(b_state["h_xx"] - 1.0))),
            "final_h_xz_max": float(np.max(np.abs(b_state["h_xz"]))),
            "final_h_zz_minus_1_max": float(np.max(np.abs(b_state["h_zz"] - 1.0))),
            "final_beta_max": float(np.max(np.abs(b_state["beta"]))),
            "final_rho_max": float(np.max(b_last["rho_current"])),
            "final_hamiltonian_residual_max": float(np.max(np.abs(ham_res))),
            "final_momentum_residual_max": float(
                max(np.max(np.abs(mom_x_res)), np.max(np.abs(mom_z_res)))
            ),
        },
        "files": {
            "rho_b_montage": str((out / "rho_b_montage.png").resolve()),
            "rho_b_log10_montage": str((out / "rho_b_log10_montage.png").resolve()),
            "rho_b_minus_a_montage": str((out / "rho_b_minus_a_montage.png").resolve()),
            "hxx_b_minus_1_montage": str((out / "hxx_b_minus_1_montage.png").resolve()),
            "ab_difference_vs_time": str((out / "ab_difference_vs_time.png").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
