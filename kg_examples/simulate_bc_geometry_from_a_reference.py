from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adm_initial_data_abc import a_complex_state_from_rho_s_time_symmetric
from adm_sources_abc import a_branch_sources_from_complex, b_branch_sources_from_ns, c_branch_effective_sources
from simulate_2p1_full_dynamics import edge_damp, initial_beamlike_fields
from simulate_ab_adm_synchronous import a_rhs, geometry_rhs, hamiltonian_residual, set_edge_constant


A_STATE_KEYS = ("phi1", "phi2", "pi1", "pi2")
B_STATE_KEYS = ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta")
C_STATE_KEYS = ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta", "phi", "pi_phi")


def solve_linearized_conformal_omega(
    energy: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
    max_iters: int = 2000,
    tol: float = 1.0e-10,
    relax: float = 0.7,
) -> dict[str, np.ndarray | int | float]:
    """
    线性化弱场估计：

        h_ij = (1 + 2 omega) delta_ij,
        Delta omega = - energy / M_P^2,

    边界条件取 omega = 0。
    """
    if not np.allclose(dx, dz, rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("当前线性化估计器假设 dx = dz。")

    h2 = dx * dx
    omega = np.zeros_like(energy)
    rhs = -energy / (mp * mp)
    delta = np.inf

    for it in range(max_iters):
        omega_new = omega.copy()
        omega_new[1:-1, 1:-1] = 0.25 * (
            omega[2:, 1:-1]
            + omega[:-2, 1:-1]
            + omega[1:-1, 2:]
            + omega[1:-1, :-2]
            - h2 * rhs[1:-1, 1:-1]
        )
        omega_new = relax * omega_new + (1.0 - relax) * omega
        omega_new[0, :] = 0.0
        omega_new[-1, :] = 0.0
        omega_new[:, 0] = 0.0
        omega_new[:, -1] = 0.0
        delta = float(np.max(np.abs(omega_new - omega)))
        omega = omega_new
        if delta < tol:
            break

    return {
        "omega": omega,
        "iterations": it + 1,
        "final_delta": delta,
    }


def apply_geometry_damping(state: dict[str, np.ndarray]) -> None:
    for key in ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta"):
        state[key] = edge_damp(state[key])
    if "phi" in state:
        state["phi"] = set_edge_constant(edge_damp(state["phi"]), 1.0)
    if "pi_phi" in state:
        state["pi_phi"] = set_edge_constant(edge_damp(state["pi_phi"]), 0.0)


def state_with_increment(
    state: dict[str, np.ndarray],
    evolve_keys: tuple[str, ...],
    deriv: dict[str, np.ndarray],
    scale: float,
) -> dict[str, np.ndarray]:
    trial = {key: value.copy() for key, value in state.items()}
    for key in evolve_keys:
        trial[key] = trial[key] + scale * deriv[f"{key}_t"]
    return trial


def rk4_update_in_place(
    state: dict[str, np.ndarray],
    evolve_keys: tuple[str, ...],
    dt: float,
    k1: dict[str, np.ndarray],
    k2: dict[str, np.ndarray],
    k3: dict[str, np.ndarray],
    k4: dict[str, np.ndarray],
) -> None:
    for key in evolve_keys:
        state[key] = state[key] + (dt / 6.0) * (
            k1[f"{key}_t"] + 2.0 * k2[f"{key}_t"] + 2.0 * k3[f"{key}_t"] + k4[f"{key}_t"]
        )


def b_reference_rhs(
    state: dict[str, np.ndarray],
    rho_ref: np.ndarray,
    s_ref: np.ndarray,
    mass: float,
    mp: float,
    dx: float,
    dz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    src = b_branch_sources_from_ns(
        rho_ref,
        s_ref,
        np.ones_like(rho_ref),
        np.zeros_like(rho_ref),
        np.zeros_like(rho_ref),
        state["h_xx"],
        state["h_xz"],
        state["h_zz"],
        state["beta"],
        mass,
        dx,
        dz,
    )
    geo = geometry_rhs(
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        k_xx=state["k_xx"],
        k_xz=state["k_xz"],
        k_zz=state["k_zz"],
        k_beta=state["k_beta"],
        energy=src["energy"],
        s_xx=src["s_xx"],
        s_xz=src["s_xz"],
        s_zz=src["s_zz"],
        s_yy_reduced=src["s_yy_reduced"],
        mp=mp,
        dx=dx,
        dz=dz,
    )
    return geo, src


def c_reference_rhs(
    state: dict[str, np.ndarray],
    rho_ref: np.ndarray,
    s_ref: np.ndarray,
    mass: float,
    mp: float,
    ell: float,
    dx: float,
    dz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    src = c_branch_effective_sources(
        rho=rho_ref,
        s=s_ref,
        phi=state["phi"],
        pi_phi=state["pi_phi"],
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        k_xx=state["k_xx"],
        k_xz=state["k_xz"],
        k_zz=state["k_zz"],
        k_beta=state["k_beta"],
        lapse=np.ones_like(rho_ref),
        shift_x=np.zeros_like(rho_ref),
        shift_z=np.zeros_like(rho_ref),
        mass=mass,
        mp=mp,
        ell=ell,
        dx=dx,
        dz=dz,
    )
    geo = geometry_rhs(
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        k_xx=state["k_xx"],
        k_xz=state["k_xz"],
        k_zz=state["k_zz"],
        k_beta=state["k_beta"],
        energy=src["energy"],
        s_xx=src["s_xx"],
        s_xz=src["s_xz"],
        s_zz=src["s_zz"],
        s_yy_reduced=src["s_yy_reduced"],
        mp=mp,
        dx=dx,
        dz=dz,
    )
    geo["phi_t"] = state["pi_phi"]
    geo["pi_phi_t"] = src["lap3_phi"] + src["k_trace"] * state["pi_phi"] + src["box_phi_corr"]
    return geo, src


def advance_rk4_b(
    state: dict[str, np.ndarray],
    rho_ref: np.ndarray,
    s_ref: np.ndarray,
    dt: float,
    mass: float,
    mp: float,
    dx: float,
    dz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    k1, _ = b_reference_rhs(state, rho_ref, s_ref, mass, mp, dx, dz)
    trial2 = state_with_increment(state, B_STATE_KEYS, k1, 0.5 * dt)
    k2, _ = b_reference_rhs(trial2, rho_ref, s_ref, mass, mp, dx, dz)
    trial3 = state_with_increment(state, B_STATE_KEYS, k2, 0.5 * dt)
    k3, _ = b_reference_rhs(trial3, rho_ref, s_ref, mass, mp, dx, dz)
    trial4 = state_with_increment(state, B_STATE_KEYS, k3, dt)
    k4, _ = b_reference_rhs(trial4, rho_ref, s_ref, mass, mp, dx, dz)
    rk4_update_in_place(state, B_STATE_KEYS, dt, k1, k2, k3, k4)
    apply_geometry_damping(state)
    return b_reference_rhs(state, rho_ref, s_ref, mass, mp, dx, dz)


def advance_rk4_c(
    state: dict[str, np.ndarray],
    rho_ref: np.ndarray,
    s_ref: np.ndarray,
    dt: float,
    mass: float,
    mp: float,
    ell: float,
    dx: float,
    dz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    k1, _ = c_reference_rhs(state, rho_ref, s_ref, mass, mp, ell, dx, dz)
    trial2 = state_with_increment(state, C_STATE_KEYS, k1, 0.5 * dt)
    k2, _ = c_reference_rhs(trial2, rho_ref, s_ref, mass, mp, ell, dx, dz)
    trial3 = state_with_increment(state, C_STATE_KEYS, k2, 0.5 * dt)
    k3, _ = c_reference_rhs(trial3, rho_ref, s_ref, mass, mp, ell, dx, dz)
    trial4 = state_with_increment(state, C_STATE_KEYS, k3, dt)
    k4, _ = c_reference_rhs(trial4, rho_ref, s_ref, mass, mp, ell, dx, dz)
    rk4_update_in_place(state, C_STATE_KEYS, dt, k1, k2, k3, k4)
    apply_geometry_damping(state)
    return c_reference_rhs(state, rho_ref, s_ref, mass, mp, ell, dx, dz)


def make_flat_bc_initial_states(
    shape: tuple[int, int],
    c_phi_mode: str = "flat",
    omega: np.ndarray | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    one = np.ones(shape)
    zero = np.zeros(shape)
    if omega is None:
        h_init = one.copy()
    else:
        h_init = 1.0 + 2.0 * omega

    b_state = {
        "h_xx": h_init.copy(),
        "h_xz": zero.copy(),
        "h_zz": h_init.copy(),
        "beta": zero.copy(),
        "k_xx": zero.copy(),
        "k_xz": zero.copy(),
        "k_zz": zero.copy(),
        "k_beta": zero.copy(),
    }

    c_state = {
        "h_xx": h_init.copy(),
        "h_xz": zero.copy(),
        "h_zz": h_init.copy(),
        "beta": zero.copy(),
        "k_xx": zero.copy(),
        "k_xz": zero.copy(),
        "k_zz": zero.copy(),
        "k_beta": zero.copy(),
        "phi": one.copy(),
        "pi_phi": zero.copy(),
    }

    if c_phi_mode != "flat":
        raise ValueError(f"当前只支持 c_phi_mode='flat'，收到 {c_phi_mode}")

    return b_state, c_state


def unwrap_2d_phase(phase: np.ndarray) -> np.ndarray:
    return np.unwrap(np.unwrap(phase, axis=0), axis=1)


def make_flat_a_reference_state(rho0: np.ndarray, s0: np.ndarray, mass: float, dx: float, dz: float) -> dict[str, np.ndarray]:
    one = np.ones_like(rho0)
    zero = np.zeros_like(rho0)
    a_init = a_complex_state_from_rho_s_time_symmetric(rho0, s0, one, zero, one, zero, mass, dx, dz)
    return {
        "h_xx": one.copy(),
        "h_xz": zero.copy(),
        "h_zz": one.copy(),
        "beta": zero.copy(),
        "k_xx": zero.copy(),
        "k_xz": zero.copy(),
        "k_zz": zero.copy(),
        "k_beta": zero.copy(),
        "phi1": a_init["phi1"].copy(),
        "phi2": a_init["phi2"].copy(),
        "pi1": a_init["pi1"].copy(),
        "pi2": a_init["pi2"].copy(),
    }


def flat_reference_rho_s(state: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    psi = state["phi1"] + 1j * state["phi2"]
    rho = np.maximum(np.abs(psi) ** 2, 1.0e-10)
    s = unwrap_2d_phase(np.angle(psi))
    return rho, s


def evolve_a_flat_reference_step(
    state: dict[str, np.ndarray],
    dt: float,
    mass: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    params = {"dx": dx, "dz": dz, "m": mass, "mp": 1.0}
    k1 = a_rhs(state, params)
    trial2 = state_with_increment(state, A_STATE_KEYS, k1, 0.5 * dt)
    k2 = a_rhs(trial2, params)
    trial3 = state_with_increment(state, A_STATE_KEYS, k2, 0.5 * dt)
    k3 = a_rhs(trial3, params)
    trial4 = state_with_increment(state, A_STATE_KEYS, k3, dt)
    k4 = a_rhs(trial4, params)
    rk4_update_in_place(state, A_STATE_KEYS, dt, k1, k2, k3, k4)
    for key in A_STATE_KEYS:
        state[key] = edge_damp(state[key], coef=0.999)
    rho, s = flat_reference_rho_s(state)
    src = a_branch_sources_from_complex(
        state["phi1"],
        state["phi2"],
        state["pi1"],
        state["pi2"],
        state["h_xx"],
        state["h_xz"],
        state["h_zz"],
        state["beta"],
        mass,
        dx,
        dz,
    )
    return rho, s, src


def run_case(
    output_dir: Path,
    steps: int = 50,
    dt: float = 5.0e-4,
    c_phi_mode: str = "flat",
) -> dict[str, object]:
    nx = nz = 161
    x_half_range = z_half_range = 5.0
    alpha = 0.5
    k0 = 10.0
    w0 = 1.0
    phi0 = 0.0
    mass = 1.0
    mp = 300.0
    ell = 0.02

    x, z, _, _, rho0, s0 = initial_beamlike_fields(x_half_range, z_half_range, nx, nz, alpha, k0, w0, phi0)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    a_ref_state = make_flat_a_reference_state(rho0, s0, mass, dx, dz)
    rho_ref, s_ref = flat_reference_rho_s(a_ref_state)

    a_src_init = a_branch_sources_from_complex(
        a_ref_state["phi1"],
        a_ref_state["phi2"],
        a_ref_state["pi1"],
        a_ref_state["pi2"],
        a_ref_state["h_xx"],
        a_ref_state["h_xz"],
        a_ref_state["h_zz"],
        a_ref_state["beta"],
        mass,
        dx,
        dz,
    )
    a_weak_init = solve_linearized_conformal_omega(a_src_init["energy"], mp, dx, dz)
    omega_init = a_weak_init["omega"]

    b_state, c_state = make_flat_bc_initial_states(rho_ref.shape, c_phi_mode=c_phi_mode, omega=omega_init)

    b_ham_history: list[float] = []
    c_ham_history: list[float] = []
    b_r3_history: list[float] = []
    c_r3_history: list[float] = []
    a_omega_max_history: list[float] = []
    a_energy_max_history: list[float] = []

    for _ in range(steps):
        a_src_flat = a_branch_sources_from_complex(
            a_ref_state["phi1"],
            a_ref_state["phi2"],
            a_ref_state["pi1"],
            a_ref_state["pi2"],
            a_ref_state["h_xx"],
            a_ref_state["h_xz"],
            a_ref_state["h_zz"],
            a_ref_state["beta"],
            mass,
            dx,
            dz,
        )
        a_weak = solve_linearized_conformal_omega(a_src_flat["energy"], mp, dx, dz)
        a_omega_max_history.append(float(np.max(np.abs(a_weak["omega"]))))
        a_energy_max_history.append(float(np.max(a_src_flat["energy"])))

        b_geo, b_src = advance_rk4_b(b_state, rho_ref, s_ref, dt, mass, mp, dx, dz)
        b_r3_history.append(float(np.max(np.abs(b_geo["r3_scalar"]))))
        b_ham = hamiltonian_residual(
            b_state["h_xx"],
            b_state["h_xz"],
            b_state["h_zz"],
            b_state["beta"],
            b_state["k_xx"],
            b_state["k_xz"],
            b_state["k_zz"],
            b_state["k_beta"],
            b_src["energy"],
            mp,
            dx,
            dz,
        )
        b_ham_history.append(float(np.max(np.abs(b_ham))))

        c_geo, c_src = advance_rk4_c(c_state, rho_ref, s_ref, dt, mass, mp, ell, dx, dz)
        c_r3_history.append(float(np.max(np.abs(c_geo["r3_scalar"]))))
        c_ham = hamiltonian_residual(
            c_state["h_xx"],
            c_state["h_xz"],
            c_state["h_zz"],
            c_state["beta"],
            c_state["k_xx"],
            c_state["k_xz"],
            c_state["k_zz"],
            c_state["k_beta"],
            c_src["energy"],
            mp,
            dx,
            dz,
        )
        c_ham_history.append(float(np.max(np.abs(c_ham))))

        rho_ref, s_ref, _ = evolve_a_flat_reference_step(a_ref_state, dt, mass, dx, dz)

    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "grid": {"nx": nx, "nz": nz, "x_half_range": x_half_range, "z_half_range": z_half_range},
        "beam_params": {
            "m": mass,
            "alpha": alpha,
            "k0": k0,
            "w0": w0,
            "phi0": phi0,
            "theta_deg": 45.0,
        },
        "evolution": {"steps": steps, "dt": dt, "mp": mp, "ell": ell, "c_phi_mode": c_phi_mode},
        "a_reference_weakfield": {
            "omega_max_history": a_omega_max_history,
            "energy_max_history": a_energy_max_history,
            "omega_max_final": a_omega_max_history[-1] if a_omega_max_history else 0.0,
        },
        "b_geometry": {
            "max_abs_r3_history": b_r3_history,
            "max_abs_hamiltonian_residual": b_ham_history,
            "final_h_xx_minus_1_max": float(np.max(np.abs(b_state["h_xx"] - 1.0))),
            "final_h_xz_max": float(np.max(np.abs(b_state["h_xz"]))),
            "final_h_zz_minus_1_max": float(np.max(np.abs(b_state["h_zz"] - 1.0))),
            "final_beta_max": float(np.max(np.abs(b_state["beta"]))),
        },
        "c_geometry": {
            "max_abs_r3_history": c_r3_history,
            "max_abs_hamiltonian_residual": c_ham_history,
            "final_h_xx_minus_1_max": float(np.max(np.abs(c_state["h_xx"] - 1.0))),
            "final_h_xz_max": float(np.max(np.abs(c_state["h_xz"]))),
            "final_h_zz_minus_1_max": float(np.max(np.abs(c_state["h_zz"] - 1.0))),
            "final_beta_max": float(np.max(np.abs(c_state["beta"]))),
            "final_phi_minus_1_max": float(np.max(np.abs(c_state["phi"] - 1.0))),
        },
        "reference_final": {
            "rho_max": float(np.max(rho_ref)),
            "rho_min": float(np.min(rho_ref)),
            "phase_max": float(np.max(s_ref)),
            "phase_min": float(np.min(s_ref)),
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=50)
    parser.add_argument("--dt", type=float, default=5.0e-4)
    parser.add_argument("--c-phi-mode", choices=["flat"], default="flat")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "visualizations" / "bc_geometry_from_a_reference",
    )
    args = parser.parse_args()
    summary = run_case(args.output, steps=args.steps, dt=args.dt, c_phi_mode=args.c_phi_mode)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
