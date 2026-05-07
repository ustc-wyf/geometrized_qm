from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_ab_measurable_observables import crop_inner, unwrap_phase
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
    spectral_k2,
    wave_imex_step,
)


def render_pair_with_split(
    path: Path,
    ts: np.ndarray,
    a_values: np.ndarray,
    b_values: np.ndarray,
    split_values: np.ndarray,
    title: str,
    ylabel: str,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)

    axes[0].plot(ts, a_values, label="A", linewidth=1.7)
    axes[0].plot(ts, b_values, label="B", linewidth=1.7)
    axes[0].set_ylabel(ylabel)
    axes[0].set_title(title)
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(ts, split_values, color="#b24a2f", linewidth=1.7)
    axes[1].set_xlabel("t / T_overlap")
    axes[1].set_ylabel("B - A")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def simulate_count_series(
    alpha: float,
    phi0: float,
    lambda_grav: float = 1.0,
    mp: float = 300.0,
    dt_target: float = 0.0125,
    inner: float = 10.0,
    outer: float = 15.0,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 81,
    nkz: int = 81,
    central_disk_radius: float = 1.5,
    box_halfwidth: float = 2.0,
) -> dict[str, np.ndarray | float]:
    params = Exact2p1Params(
        alpha=alpha,
        phi0=phi0,
        x_half_range=outer,
        z_half_range=outer,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    total_time = params.overlap_time
    steps = int(round(total_time / dt_target))
    dt = total_time / steps
    k2 = spectral_k2(params.nx, params.nz, dx, dz)

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, params)
    sponge_inner = max(inner + 2.0, 0.8 * outer)
    sigma = make_sponge(x, z, inner_half_range=sponge_inner, outer_half_range=outer, sigma_max=2.0)

    tau_b = np.zeros((params.nx, params.nz))
    beta_b = np.zeros_like(tau_b)
    ptau_b = np.zeros_like(tau_b)
    pbeta_b = np.zeros_like(tau_b)

    ix, iz, xi, zi = crop_inner(x, z, inner)
    xi_grid, zi_grid = np.meshgrid(xi, zi, indexing="ij")
    disk_mask = (xi_grid**2 + zi_grid**2) <= central_disk_radius**2
    box_mask = (np.abs(xi_grid) <= box_halfwidth) & (np.abs(zi_grid) <= box_halfwidth)

    ts: list[float] = []
    a_disk: list[float] = []
    b_disk: list[float] = []
    a_box: list[float] = []
    b_box: list[float] = []

    def record(current_time: float) -> None:
        rho_a_outer = np.abs(psi_now) ** 2
        rho_b_outer = rho_a_outer * np.exp(-np.clip(beta_b + tau_b, -12.0, 12.0))
        a = rho_a_outer[np.ix_(ix, iz)]
        b = rho_b_outer[np.ix_(ix, iz)]

        ts.append(current_time / total_time)
        a_disk.append(float(np.sum(a[disk_mask])))
        b_disk.append(float(np.sum(b[disk_mask])))
        a_box.append(float(np.sum(a[box_mask])))
        b_box.append(float(np.sum(b[box_mask])))

    record(0.0)

    for n in range(steps):
        rho_a_outer = np.abs(psi_now) ** 2
        s_a_outer = unwrap_phase(psi_now)

        src_tau_b, src_beta_b, *_ = rhs_branch_b(
            tau_b, ptau_b, beta_b, pbeta_b, rho_a_outer, s_a_outer, params.m, mp, dx, dz, k2
        )
        tau_b, ptau_b = wave_imex_step(tau_b, ptau_b, lambda_grav * src_tau_b, dt, k2)
        beta_b, pbeta_b = wave_imex_step(beta_b, pbeta_b, lambda_grav * src_beta_b, dt, k2)
        edge_damp(tau_b, sigma, dt, -2.0, 2.0)
        edge_damp(beta_b, sigma, dt, -2.0, 2.0)
        edge_damp(ptau_b, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_b, sigma, dt, -10.0, 10.0)
        apply_geometry_boundary(tau_b, beta_b)

        t_next = (n + 1) * dt
        psi_bdry = boundary_exact(alpha, t_next, x, z, params)
        rhs = spectral_lap(psi_now, dx, dz) - params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_bdry)
        psi_prev, psi_now = psi_now, psi_next

        record(t_next)

    return {
        "t_over_T": np.array(ts),
        "A_disk": np.array(a_disk),
        "B_disk": np.array(b_disk),
        "A_box": np.array(a_box),
        "B_box": np.array(b_box),
        "dt": dt,
        "steps": steps,
    }


def normalized_difference(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    return num / np.maximum(den, 1.0e-30)


def summarize_split(a_values: np.ndarray, b_values: np.ndarray, ts: np.ndarray) -> dict[str, float]:
    split = b_values - a_values
    idx = int(np.argmax(np.abs(split)))
    rel = float(split[-1] / a_values[-1]) if abs(a_values[-1]) > 1.0e-30 else float("nan")
    return {
        "A_final": float(a_values[-1]),
        "B_final": float(b_values[-1]),
        "split_final_B_minus_A": float(split[-1]),
        "split_final_relative_to_A": rel,
        "max_abs_split": float(np.max(np.abs(split))),
        "time_of_max_abs_split": float(ts[idx]),
        "split_at_time_of_max_abs": float(split[idx]),
    }


def run(
    outdir: str | Path,
    lambda_grav: float = 1.0,
    dt_target: float = 0.0125,
    inner: float = 10.0,
    outer: float = 15.0,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 81,
    nkz: int = 81,
    central_disk_radius: float = 1.5,
    box_halfwidth: float = 2.0,
) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    common = {
        "lambda_grav": lambda_grav,
        "dt_target": dt_target,
        "inner": inner,
        "outer": outer,
        "nx": nx,
        "nz": nz,
        "nkx": nkx,
        "nkz": nkz,
        "central_disk_radius": central_disk_radius,
        "box_halfwidth": box_halfwidth,
    }

    both_phi0 = simulate_count_series(alpha=0.5, phi0=0.0, **common)
    both_phipi = simulate_count_series(alpha=0.5, phi0=np.pi, **common)
    beam_1 = simulate_count_series(alpha=1.0, phi0=0.0, **common)
    beam_2 = simulate_count_series(alpha=0.0, phi0=0.0, **common)

    ts = both_phi0["t_over_T"]

    phase_flip_A_disk = normalized_difference(
        both_phi0["A_disk"] - both_phipi["A_disk"],
        both_phi0["A_disk"] + both_phipi["A_disk"],
    )
    phase_flip_B_disk = normalized_difference(
        both_phi0["B_disk"] - both_phipi["B_disk"],
        both_phi0["B_disk"] + both_phipi["B_disk"],
    )
    phase_flip_A_box = normalized_difference(
        both_phi0["A_box"] - both_phipi["A_box"],
        both_phi0["A_box"] + both_phipi["A_box"],
    )
    phase_flip_B_box = normalized_difference(
        both_phi0["B_box"] - both_phipi["B_box"],
        both_phi0["B_box"] + both_phipi["B_box"],
    )

    interference_A_disk = normalized_difference(
        both_phi0["A_disk"] - (beam_1["A_disk"] + beam_2["A_disk"]),
        beam_1["A_disk"] + beam_2["A_disk"],
    )
    interference_B_disk = normalized_difference(
        both_phi0["B_disk"] - (beam_1["B_disk"] + beam_2["B_disk"]),
        beam_1["B_disk"] + beam_2["B_disk"],
    )
    interference_A_box = normalized_difference(
        both_phi0["A_box"] - (beam_1["A_box"] + beam_2["A_box"]),
        beam_1["A_box"] + beam_2["A_box"],
    )
    interference_B_box = normalized_difference(
        both_phi0["B_box"] - (beam_1["B_box"] + beam_2["B_box"]),
        beam_1["B_box"] + beam_2["B_box"],
    )

    phase_flip_split_disk = phase_flip_B_disk - phase_flip_A_disk
    phase_flip_split_box = phase_flip_B_box - phase_flip_A_box
    interference_split_disk = interference_B_disk - interference_A_disk
    interference_split_box = interference_B_box - interference_A_box

    render_pair_with_split(
        out / "phase_flip_central_disk.png",
        ts,
        phase_flip_A_disk,
        phase_flip_B_disk,
        phase_flip_split_disk,
        "Phase-flip differential, central disk",
        "D_phi",
    )
    render_pair_with_split(
        out / "phase_flip_interference_box.png",
        ts,
        phase_flip_A_box,
        phase_flip_B_box,
        phase_flip_split_box,
        "Phase-flip differential, interference box",
        "D_phi",
    )
    render_pair_with_split(
        out / "interference_excess_central_disk.png",
        ts,
        interference_A_disk,
        interference_B_disk,
        interference_split_disk,
        "Interference excess, central disk",
        "E_int",
    )
    render_pair_with_split(
        out / "interference_excess_interference_box.png",
        ts,
        interference_A_box,
        interference_B_box,
        interference_split_box,
        "Interference excess, interference box",
        "E_int",
    )

    np.savez_compressed(
        out / "series.npz",
        t_over_T=ts,
        phase_flip_A_disk=phase_flip_A_disk,
        phase_flip_B_disk=phase_flip_B_disk,
        phase_flip_split_disk=phase_flip_split_disk,
        phase_flip_A_box=phase_flip_A_box,
        phase_flip_B_box=phase_flip_B_box,
        phase_flip_split_box=phase_flip_split_box,
        interference_A_disk=interference_A_disk,
        interference_B_disk=interference_B_disk,
        interference_split_disk=interference_split_disk,
        interference_A_box=interference_A_box,
        interference_B_box=interference_B_box,
        interference_split_box=interference_split_box,
    )

    summary = {
        "lambda_grav": lambda_grav,
        "simulation_settings": {
            "dt_target": dt_target,
            "nx": nx,
            "nz": nz,
            "nkx": nkx,
            "nkz": nkz,
            "inner": inner,
            "outer": outer,
        },
        "detector_windows": {
            "central_disk_radius": central_disk_radius,
            "interference_box_halfwidth": box_halfwidth,
        },
        "phase_flip_differential": {
            "central_disk": summarize_split(phase_flip_A_disk, phase_flip_B_disk, ts),
            "interference_box": summarize_split(phase_flip_A_box, phase_flip_B_box, ts),
        },
        "interference_excess": {
            "central_disk": summarize_split(interference_A_disk, interference_B_disk, ts),
            "interference_box": summarize_split(interference_A_box, interference_B_box, ts),
        },
        "consistency_checks": {
            "single_beam_final_disk_A_alpha1": float(beam_1["A_disk"][-1]),
            "single_beam_final_disk_A_alpha0": float(beam_2["A_disk"][-1]),
            "single_beam_final_box_A_alpha1": float(beam_1["A_box"][-1]),
            "single_beam_final_box_A_alpha0": float(beam_2["A_box"][-1]),
        },
        "notes": {
            "phase_flip_definition": "D_phi = (N(phi=0) - N(phi=pi)) / (N(phi=0) + N(phi=pi))",
            "interference_excess_definition": "E_int = (N_both - (N_beam1 + N_beam2)) / (N_beam1 + N_beam2)",
            "theory_split_definition": "For each experimental observable, the reported AB split is B_prediction - A_prediction for the same differential readout.",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    compact = {
        "phase_flip_differential": summary["phase_flip_differential"],
        "interference_excess": summary["interference_excess"],
    }
    print(json.dumps(compact, indent=2))
    return summary


if __name__ == "__main__":
    result = run(
        Path(__file__).resolve().parent.parent / "visualizations" / "ab_differential_experiment_observables_lambda1"
    )
