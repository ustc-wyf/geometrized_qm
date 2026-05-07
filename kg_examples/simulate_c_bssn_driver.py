from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adm_initial_data_abc import a_complex_state_from_rho_s_time_symmetric
from adm_matter_2p1 import matter_observables, matter_rhs
from adm_sources_abc import a_branch_sources_from_complex, c_branch_effective_sources
from simulate_b_bssn_driver import (
    adm_to_bssn_state,
    apply_bssn_boundary,
    bssn_to_adm_geometry,
    conformal_connection_functions,
    physical_geometry_rhs_shifted,
    rk4_update_in_place,
    state_metrics,
)
from simulate_bc_from_a_initial_data import (
    exact_beamlike_boundary_fields,
    exact_localized_boundary_fields,
    make_bc_initial_states,
    recover_rho_from_n,
)
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    phase_field,
    spectral_omega,
)
from mixed_tilde_initial_data import localized_direct_tilde_adm_initial


def adm_c_to_bssn_state(adm_state: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    state = adm_to_bssn_state(adm_state)
    state["phi"] = adm_state["phi"].copy()
    state["pi_phi"] = adm_state["pi_phi"].copy()
    return state


def apply_c_bssn_boundary(state: dict[str, np.ndarray]) -> None:
    apply_bssn_boundary(state)
    state["phi"] = state["phi"].copy()
    state["pi_phi"] = state["pi_phi"].copy()
    state["phi"][:8, :] = 1.0
    state["phi"][-8:, :] = 1.0
    state["phi"][:, :8] = 1.0
    state["phi"][:, -8:] = 1.0
    state["pi_phi"][:8, :] = 0.0
    state["pi_phi"][-8:, :] = 0.0
    state["pi_phi"][:, :8] = 0.0
    state["pi_phi"][:, -8:] = 0.0


def scalar_advection(field: np.ndarray, shift_x: np.ndarray, shift_z: np.ndarray, dx: float, dz: float) -> np.ndarray:
    from adm_matter_2p1 import ddx, ddz

    return shift_x * ddx(field, dx) + shift_z * ddz(field, dz)


def apply_matter_characteristic_bc_c_bssn(state: dict[str, np.ndarray], t: float, params: dict[str, float]) -> None:
    geom = bssn_to_adm_geometry(state)
    if params["initial_mode"] == "localized":
        rho_ref, s_ref = exact_localized_boundary_fields(t, params)
    else:
        rho_ref, s_ref = exact_beamlike_boundary_fields(t, params)

    rho_cur, _ = recover_rho_from_n(
        n=state["n"],
        s=state["s"],
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta=geom["beta_geom"],
        mass=params["m"],
        dx=params["dx"],
        dz=params["dz"],
    )
    obs = matter_observables(
        rho=rho_cur,
        s=state["s"],
        lapse=state["lapse"],
        shift_x=state["shift_x"],
        shift_z=state["shift_z"],
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta=geom["beta_geom"],
        mass=params["m"],
        dx=params["dx"],
        dz=params["dz"],
    )
    vx = obs["flux_x"] / np.maximum(obs["density"], 1.0e-12)
    vz = obs["flux_z"] / np.maximum(obs["density"], 1.0e-12)

    rho_ref_obs = matter_observables(
        rho=rho_ref,
        s=s_ref,
        lapse=state["lapse"],
        shift_x=state["shift_x"],
        shift_z=state["shift_z"],
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta=geom["beta_geom"],
        mass=params["m"],
        dx=params["dx"],
        dz=params["dz"],
    )
    n_ref = rho_ref_obs["density"]

    inflow = vx[0, :] > 0.0
    state["n"][0, inflow] = n_ref[0, inflow]
    state["s"][0, inflow] = s_ref[0, inflow]
    state["n"][0, ~inflow] = state["n"][1, ~inflow]
    state["s"][0, ~inflow] = state["s"][1, ~inflow]

    inflow = vx[-1, :] < 0.0
    state["n"][-1, inflow] = n_ref[-1, inflow]
    state["s"][-1, inflow] = s_ref[-1, inflow]
    state["n"][-1, ~inflow] = state["n"][-2, ~inflow]
    state["s"][-1, ~inflow] = state["s"][-2, ~inflow]

    inflow = vz[:, 0] > 0.0
    state["n"][inflow, 0] = n_ref[inflow, 0]
    state["s"][inflow, 0] = s_ref[inflow, 0]
    state["n"][~inflow, 0] = state["n"][~inflow, 1]
    state["s"][~inflow, 0] = state["s"][~inflow, 1]

    inflow = vz[:, -1] < 0.0
    state["n"][inflow, -1] = n_ref[inflow, -1]
    state["s"][inflow, -1] = s_ref[inflow, -1]
    state["n"][~inflow, -1] = state["n"][~inflow, -2]
    state["s"][~inflow, -1] = state["s"][~inflow, -2]

    state["n"][0, 0] = 0.5 * (state["n"][1, 0] + state["n"][0, 1])
    state["n"][0, -1] = 0.5 * (state["n"][1, -1] + state["n"][0, -2])
    state["n"][-1, 0] = 0.5 * (state["n"][-2, 0] + state["n"][-1, 1])
    state["n"][-1, -1] = 0.5 * (state["n"][-2, -1] + state["n"][-1, -2])
    state["s"][0, 0] = 0.5 * (state["s"][1, 0] + state["s"][0, 1])
    state["s"][0, -1] = 0.5 * (state["s"][1, -1] + state["s"][0, -2])
    state["s"][-1, 0] = 0.5 * (state["s"][-2, 0] + state["s"][-1, 1])
    state["s"][-1, -1] = 0.5 * (state["s"][-2, -1] + state["s"][-1, -2])


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


def assert_finite_state(state: dict[str, np.ndarray], label: str) -> None:
    for key, value in state.items():
        if not np.all(np.isfinite(value)):
            raise FloatingPointError(f"{label}: 变量 {key} 出现非有限值")


def c_bssn_rhs(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]
    ell = params["ell"]
    mu_shift = params.get("mu_shift", 0.5)
    eta_shift = params.get("eta_shift", 0.5)

    geom = bssn_to_adm_geometry(state)
    rho, _ = recover_rho_from_n(
        n=state["n"],
        s=state["s"],
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta=geom["beta_geom"],
        mass=mass,
        dx=dx,
        dz=dz,
    )

    src = c_branch_effective_sources(
        rho=rho,
        s=state["s"],
        phi=state["phi"],
        pi_phi=state["pi_phi"],
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta=geom["beta_geom"],
        k_xx=geom["k_xx"],
        k_xz=geom["k_xz"],
        k_zz=geom["k_zz"],
        k_beta=geom["k_beta"],
        lapse=state["lapse"],
        shift_x=state["shift_x"],
        shift_z=state["shift_z"],
        mass=mass,
        mp=mp,
        ell=ell,
        dx=dx,
        dz=dz,
    )

    phys_rhs = physical_geometry_rhs_shifted(
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta_geom=geom["beta_geom"],
        k_xx=geom["k_xx"],
        k_xz=geom["k_xz"],
        k_zz=geom["k_zz"],
        k_beta=geom["k_beta"],
        lapse=state["lapse"],
        shift_x=state["shift_x"],
        shift_z=state["shift_z"],
        energy=src["energy"],
        s_xx=src["s_xx"],
        s_xz=src["s_xz"],
        s_zz=src["s_zz"],
        s_yy_reduced=src["s_yy_reduced"],
        mp=mp,
        dx=dx,
        dz=dz,
    )

    mrhs = matter_rhs(
        rho=rho,
        s=state["s"],
        lapse=state["lapse"],
        shift_x=state["shift_x"],
        shift_z=state["shift_z"],
        h_xx=geom["h_xx"],
        h_xz=geom["h_xz"],
        h_zz=geom["h_zz"],
        beta=geom["beta_geom"],
        mass=mass,
        dx=dx,
        dz=dz,
    )

    det_h = np.clip(geom["h_xx"] * geom["h_zz"] - geom["h_xz"] * geom["h_xz"], 1.0e-18, None)
    det_h_t = phys_rhs["h_xx_t"] * geom["h_zz"] + geom["h_xx"] * phys_rhs["h_zz_t"] - 2.0 * geom["h_xz"] * phys_rhs["h_xz_t"]
    chi_t = -(state["chi"] / 3.0) * (det_h_t / det_h + 2.0 * phys_rhs["beta_geom_t"])

    def inv_rhs_2x2(
        a: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        a_t: np.ndarray,
        b_t: np.ndarray,
        c_t: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        from adm_ykilling_geometry import inverse_2metric

        ai, bi, ci, _ = inverse_2metric(a, b, c)
        inv_t_xx = -(ai * a_t * ai + 2.0 * ai * b_t * bi + bi * c_t * bi)
        inv_t_xz = -(ai * a_t * bi + ai * b_t * ci + bi * b_t * bi + bi * c_t * ci)
        inv_t_zz = -(bi * a_t * bi + 2.0 * bi * b_t * ci + ci * c_t * ci)
        return inv_t_xx, inv_t_xz, inv_t_zz

    h_inv_t_xx, h_inv_t_xz, h_inv_t_zz = inv_rhs_2x2(
        geom["h_xx"],
        geom["h_xz"],
        geom["h_zz"],
        phys_rhs["h_xx_t"],
        phys_rhs["h_xz_t"],
        phys_rhs["h_zz_t"],
    )

    k_trace_t = (
        h_inv_t_xx * geom["k_xx"]
        + 2.0 * h_inv_t_xz * geom["k_xz"]
        + h_inv_t_zz * geom["k_zz"]
        + (
            geom["h_xx_inv"] * phys_rhs["k_xx_t"]
            + 2.0 * geom["h_xz_inv"] * phys_rhs["k_xz_t"]
            + geom["h_zz_inv"] * phys_rhs["k_zz_t"]
        )
        + phys_rhs["k_beta_t"]
    )

    gt_xx_t = chi_t * geom["h_xx"] + state["chi"] * phys_rhs["h_xx_t"]
    gt_xz_t = chi_t * geom["h_xz"] + state["chi"] * phys_rhs["h_xz_t"]
    gt_zz_t = chi_t * geom["h_zz"] + state["chi"] * phys_rhs["h_zz_t"]

    at_xx_t = chi_t * (geom["k_xx"] - geom["h_xx"] * state["k_trace"] / 3.0) + state["chi"] * (
        phys_rhs["k_xx_t"] - phys_rhs["h_xx_t"] * state["k_trace"] / 3.0 - geom["h_xx"] * k_trace_t / 3.0
    )
    at_xz_t = chi_t * (geom["k_xz"] - geom["h_xz"] * state["k_trace"] / 3.0) + state["chi"] * (
        phys_rhs["k_xz_t"] - phys_rhs["h_xz_t"] * state["k_trace"] / 3.0 - geom["h_xz"] * k_trace_t / 3.0
    )
    at_zz_t = chi_t * (geom["k_zz"] - geom["h_zz"] * state["k_trace"] / 3.0) + state["chi"] * (
        phys_rhs["k_zz_t"] - phys_rhs["h_zz_t"] * state["k_trace"] / 3.0 - geom["h_zz"] * k_trace_t / 3.0
    )

    gamma_x, gamma_z = conformal_connection_functions(state["gt_xx"], state["gt_xz"], state["gt_zz"], dx, dz)
    lapse_t = scalar_advection(state["lapse"], state["shift_x"], state["shift_z"], dx, dz) - 2.0 * state["lapse"] * state["k_trace"]
    shift_x_t = mu_shift * gamma_x - eta_shift * state["shift_x"]
    shift_z_t = mu_shift * gamma_z - eta_shift * state["shift_z"]
    phi_t = scalar_advection(state["phi"], state["shift_x"], state["shift_z"], dx, dz) + state["lapse"] * state["pi_phi"]
    pi_phi_t = scalar_advection(state["pi_phi"], state["shift_x"], state["shift_z"], dx, dz) + state["lapse"] * (
        src["lap3_phi"] + src["k_trace"] * state["pi_phi"] + src["box_phi_corr"]
    )

    return {
        "chi_t": chi_t,
        "gt_xx_t": gt_xx_t,
        "gt_xz_t": gt_xz_t,
        "gt_zz_t": gt_zz_t,
        "k_trace_t": k_trace_t,
        "at_xx_t": at_xx_t,
        "at_xz_t": at_xz_t,
        "at_zz_t": at_zz_t,
        "lapse_t": lapse_t,
        "shift_x_t": shift_x_t,
        "shift_z_t": shift_z_t,
        "n_t": mrhs["n_t"],
        "s_t": mrhs["s_t"],
        "phi_t": phi_t,
        "pi_phi_t": pi_phi_t,
        "rho_current": rho,
        "energy_current": src["energy"],
        "beta_geom_current": geom["beta_geom"],
        "gamma_driver_x": gamma_x,
        "gamma_driver_z": gamma_z,
        "r3_scalar": phys_rhs["r3_scalar"],
    }


def advance_rk4_c(
    state: dict[str, np.ndarray],
    evolve_keys: tuple[str, ...],
    dt: float,
    params: dict[str, float],
    t_next: float,
    do_finite_check: bool = True,
) -> dict[str, np.ndarray]:
    if do_finite_check:
        assert_finite_state(state, "步进前")
    k1 = c_bssn_rhs(state, params)
    trial2 = state_with_increment(state, evolve_keys, k1, 0.5 * dt)
    if do_finite_check:
        assert_finite_state(trial2, "trial2")
    k2 = c_bssn_rhs(trial2, params)
    trial3 = state_with_increment(state, evolve_keys, k2, 0.5 * dt)
    if do_finite_check:
        assert_finite_state(trial3, "trial3")
    k3 = c_bssn_rhs(trial3, params)
    trial4 = state_with_increment(state, evolve_keys, k3, dt)
    if do_finite_check:
        assert_finite_state(trial4, "trial4")
    k4 = c_bssn_rhs(trial4, params)
    rk4_update_in_place(state, evolve_keys, dt, k1, k2, k3, k4)
    apply_matter_characteristic_bc_c_bssn(state, t_next, params)
    apply_c_bssn_boundary(state)
    if do_finite_check:
        assert_finite_state(state, "步进后")
    return c_bssn_rhs(state, params)


def main() -> None:
    parser = argparse.ArgumentParser(description="C 支 BSSN 型强场原型")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "visualizations" / "c_bssn_driver")
    parser.add_argument("--steps", type=int, default=800)
    parser.add_argument("--dt", type=float, default=2.5e-4)
    parser.add_argument("--finite-check-every", type=int, default=20)
    parser.add_argument("--mu-shift", type=float, default=0.5)
    parser.add_argument("--eta-shift", type=float, default=0.5)
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--c-phi-mode", choices=["flat", "r3_proxy", "elliptic"], default="r3_proxy")
    parser.add_argument("--elliptic-lambda", type=float, default=1.0e-4)
    args = parser.parse_args()

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
    omega_loc = spectral_omega(x, z, loc.m)

    mass = 1.0
    mp = 300.0
    direct_init = localized_direct_tilde_adm_initial(
        psi0=psi0,
        psi0_hat=np.fft.fft2(psi0),
        omega=omega_loc,
        x=x,
        z=z,
    )
    _, c_adm = make_bc_initial_states(
        rho0, s0, None, mass, mp, args.ell, dx, dz,
        c_phi_mode=args.c_phi_mode, elliptic_lambda=args.elliptic_lambda,
        tilde_adm_init=direct_init["adm"],
    )
    c_state = adm_c_to_bssn_state(c_adm)

    params = {
        "dx": dx,
        "dz": dz,
        "m": mass,
        "mp": mp,
        "ell": args.ell,
        "mu_shift": args.mu_shift,
        "eta_shift": args.eta_shift,
        "alpha": loc.alpha,
        "k0": loc.k0,
        "phi0": loc.phi0,
        "theta": np.pi / 4.0,
        "x_half_range": loc.x_half_range,
        "z_half_range": loc.z_half_range,
        "nx": loc.nx,
        "nz": loc.nz,
        "initial_mode": "localized",
        "localized_psi0": psi0,
        "localized_omega": omega_loc,
        "localized_psi0_hat": np.fft.fft2(psi0),
    }

    evolve_keys = (
        "chi", "gt_xx", "gt_xz", "gt_zz", "k_trace", "at_xx", "at_xz", "at_zz",
        "lapse", "shift_x", "shift_z", "n", "s", "phi", "pi_phi"
    )
    history = []
    t = 0.0
    for step in range(1, args.steps + 1):
        t_next = t + args.dt
        do_check = args.finite_check_every > 0 and (step % args.finite_check_every == 0 or step == 1)
        rhs_out = advance_rk4_c(c_state, evolve_keys, args.dt, params, t_next, do_finite_check=do_check)
        t = t_next
        if step % 20 == 0 or step == 1:
            history.append({
                "step": step,
                "t": t,
                **state_metrics(c_state),
                "phi_minus_1_max": float(np.max(np.abs(c_state["phi"] - 1.0))),
                "pi_phi_abs_max": float(np.max(np.abs(c_state["pi_phi"]))),
                "r3_abs_max": float(np.max(np.abs(rhs_out["r3_scalar"]))),
                "rho_max": float(np.max(rhs_out["rho_current"])),
                "gamma_driver_x_abs_max": float(np.max(np.abs(rhs_out["gamma_driver_x"]))),
                "gamma_driver_z_abs_max": float(np.max(np.abs(rhs_out["gamma_driver_z"]))),
            })

    args.output.mkdir(parents=True, exist_ok=True)
    summary = {
        "steps": args.steps,
        "dt": args.dt,
        "mu_shift": args.mu_shift,
        "eta_shift": args.eta_shift,
        "ell": args.ell,
        "c_phi_mode": args.c_phi_mode,
        "elliptic_lambda": args.elliptic_lambda,
        "final": {
            **state_metrics(c_state),
            "phi_minus_1_max": float(np.max(np.abs(c_state["phi"] - 1.0))),
            "pi_phi_abs_max": float(np.max(np.abs(c_state["pi_phi"]))),
        },
        "history": history,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
