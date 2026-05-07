from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adm_initial_data_abc import (
    a_complex_state_from_rho_s_time_symmetric,
    solve_average_momentum_extrinsic_curvature,
)
from adm_sources_abc import (
    b_branch_sources_from_ns,
    c_branch_effective_sources,
    u_of_phi,
    u_phi_of_phi,
)
from adm_matter_2p1 import matter_observables, matter_rhs, recover_rho
from adm_sources_abc import a_branch_sources_from_complex
from adm_ykilling_geometry import beta_hessian_and_spatial_ricci, momentum_constraint_residual_ykilling, scalar_laplacian_ykilling
from adm_ykilling_geometry import inverse_2metric
from kg_double_slit_beamlike import compute_xz_slice, default_params
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    evolve_positive_frequency,
    initial_wavefunction,
    make_grid,
    phase_field,
    spectral_omega,
)
from simulate_2p1_full_dynamics import initial_beamlike_fields, unwrap_2d_phase
from simulate_ab_adm_synchronous import geometry_rhs, hamiltonian_residual, set_edge_constant
from simulate_bc_geometry_from_a_reference import solve_linearized_conformal_omega


B_STATE_KEYS = ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta", "n", "s")
C_STATE_KEYS = ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta", "n", "s", "phi", "pi_phi")


def boundary_ratio(rho: np.ndarray) -> float:
    boundary = np.concatenate([rho[0, :], rho[-1, :], rho[:, 0], rho[:, -1]])
    return float(np.max(boundary) / max(float(np.max(rho)), 1.0e-30))


def make_a_flat_initial_data(
    rho0: np.ndarray,
    s0: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    one = np.ones_like(rho0)
    zero = np.zeros_like(rho0)
    a_init = a_complex_state_from_rho_s_time_symmetric(rho0, s0, one, zero, one, zero, mass, dx, dz)
    a_src = a_branch_sources_from_complex(
        a_init["phi1"],
        a_init["phi2"],
        a_init["pi1"],
        a_init["pi2"],
        one,
        zero,
        one,
        zero,
        mass,
        dx,
        dz,
    )
    return a_init, a_src


def recover_rho_from_n(
    n: np.ndarray,
    s: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    ones = np.ones_like(n)
    zeros = np.zeros_like(n)
    obs_unit = matter_observables(
        rho=ones,
        s=s,
        lapse=ones,
        shift_x=zeros,
        shift_z=zeros,
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        mass=mass,
        dx=dx,
        dz=dz,
    )
    rho = recover_rho(n, obs_unit["energy"], beta, obs_unit["sqrt_h"])
    return rho, obs_unit


def make_bc_initial_states(
    rho0: np.ndarray,
    s0: np.ndarray,
    omega_init: np.ndarray | None,
    mass: float,
    mp: float,
    ell: float,
    dx: float,
    dz: float,
    c_phi_mode: str = "flat",
    elliptic_lambda: float = 1.0,
    tilde_adm_init: dict[str, np.ndarray] | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    zero = np.zeros_like(rho0)
    one = np.ones_like(rho0)
    if tilde_adm_init is not None:
        h_xx_init = tilde_adm_init["h_xx"].copy()
        h_xz_init = tilde_adm_init["h_xz"].copy()
        h_zz_init = tilde_adm_init["h_zz"].copy()
        beta_init = tilde_adm_init["beta"].copy()
        lapse_init = tilde_adm_init["lapse"].copy()
        shift_x_init = tilde_adm_init["shift_x"].copy()
        shift_z_init = tilde_adm_init["shift_z"].copy()
        k_xx_init = tilde_adm_init["k_xx"].copy()
        k_xz_init = tilde_adm_init["k_xz"].copy()
        k_zz_init = tilde_adm_init["k_zz"].copy()
        k_beta_init = tilde_adm_init["k_beta"].copy()
    else:
        if omega_init is None:
            raise ValueError("当未提供 tilde_adm_init 时，omega_init 不能为空。")
        h_init = 1.0 + 2.0 * omega_init
        h_xx_init = h_init.copy()
        h_xz_init = zero.copy()
        h_zz_init = h_init.copy()
        beta_init = zero.copy()
        lapse_init = one.copy()
        shift_x_init = zero.copy()
        shift_z_init = zero.copy()

        rho_tmp = rho0.copy()
        b_src_init = b_branch_sources_from_ns(
            rho=rho_tmp,
            s=s0,
            lapse=lapse_init,
            shift_x=shift_x_init,
            shift_z=shift_z_init,
            h_xx=h_xx_init,
            h_xz=h_xz_init,
            h_zz=h_zz_init,
            beta=beta_init,
            mass=mass,
            dx=dx,
            dz=dz,
        )
        ksolve = solve_average_momentum_extrinsic_curvature(
            b_src_init["mom_x"],
            b_src_init["mom_z"],
            mp,
            dx,
            dz,
        )
        k_xx_init = ksolve["k_xx"].copy()
        k_xz_init = ksolve["k_xz"].copy()
        k_zz_init = ksolve["k_zz"].copy()
        k_beta_init = zero.copy()

    obs_init = matter_observables(
        rho=rho0,
        s=s0,
        lapse=lapse_init,
        shift_x=shift_x_init,
        shift_z=shift_z_init,
        h_xx=h_xx_init,
        h_xz=h_xz_init,
        h_zz=h_zz_init,
        beta=beta_init,
        mass=mass,
        dx=dx,
        dz=dz,
    )
    rho_init = rho0.copy()
    # 由同一组 (rho, S, tilde g) 直接构造守恒密度 n，避免在 t=0 人为引入 A/B 偏差。
    n_init = np.exp(beta_init) * obs_init["sqrt_h"] * rho_init * obs_init["energy"]

    b_state = {
        "h_xx": h_xx_init.copy(),
        "h_xz": h_xz_init.copy(),
        "h_zz": h_zz_init.copy(),
        "beta": beta_init.copy(),
        "lapse": lapse_init.copy(),
        "shift_x": shift_x_init.copy(),
        "shift_z": shift_z_init.copy(),
        "k_xx": k_xx_init.copy(),
        "k_xz": k_xz_init.copy(),
        "k_zz": k_zz_init.copy(),
        "k_beta": k_beta_init.copy(),
        "n": n_init.copy(),
        "s": s0.copy(),
    }

    c_state = {
        "h_xx": h_xx_init.copy(),
        "h_xz": h_xz_init.copy(),
        "h_zz": h_zz_init.copy(),
        "beta": beta_init.copy(),
        "lapse": lapse_init.copy(),
        "shift_x": shift_x_init.copy(),
        "shift_z": shift_z_init.copy(),
        "k_xx": k_xx_init.copy(),
        "k_xz": k_xz_init.copy(),
        "k_zz": k_zz_init.copy(),
        "k_beta": k_beta_init.copy(),
        "n": n_init.copy(),
        "s": s0.copy(),
        "phi": one.copy(),
        "pi_phi": zero.copy(),
    }
    if c_phi_mode == "r3_proxy":
        geo_init = beta_hessian_and_spatial_ricci(h_init, zero, h_init, zero, dx, dz)
        r_proxy = np.maximum(geo_init["r3_scalar"], 0.0)
        c_state["phi"] = np.power(1.0 + np.power((ell * ell) * r_proxy, 2.0), -1.5)
    elif c_phi_mode == "elliptic":
        phi_raw = solve_phi_initial_elliptic(
            rho=rho_init,
            s=s0,
            h_xx=h_init,
            h_xz=zero,
            h_zz=h_init,
            beta=zero,
            mass=mass,
            mp=mp,
            ell=ell,
            dx=dx,
            dz=dz,
        )
        c_state["phi"] = np.clip(1.0 - elliptic_lambda * (1.0 - phi_raw), 1.0e-6, 1.0)
    elif c_phi_mode != "flat":
        raise ValueError(f"未知的 c_phi_mode: {c_phi_mode}")
    return b_state, c_state


def solve_phi_initial_elliptic(
    rho: np.ndarray,
    s: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    mp: float,
    ell: float,
    dx: float,
    dz: float,
    iterations: int = 4000,
    relax: float | None = None,
) -> np.ndarray:
    one = np.ones_like(rho)
    zero = np.zeros_like(rho)
    if relax is None:
        relax = 0.1 * min(dx * dx, dz * dz)

    bsrc = b_branch_sources_from_ns(
        rho=rho,
        s=s,
        lapse=one,
        shift_x=zero,
        shift_z=zero,
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        mass=mass,
        dx=dx,
        dz=dz,
    )
    trace_b = 2.0 * mass * mass * bsrc["rho"]
    phi = np.ones_like(rho)
    for _ in range(iterations):
        lap3_phi = scalar_laplacian_ykilling(h_xx, h_xz, h_zz, beta, phi, dx, dz)
        source = (2.0 * u_of_phi(phi, ell) - phi * u_phi_of_phi(phi, ell)) / 3.0 - trace_b / (3.0 * mp * mp)
        residual = lap3_phi - source
        phi = np.clip(phi - relax * residual, 1.0e-6, 1.0)
        phi = set_edge_constant(phi, 1.0)
    return phi


def apply_damping_bc(state: dict[str, np.ndarray]) -> None:
    for key, value in (
        ("h_xx", 1.0),
        ("h_xz", 0.0),
        ("h_zz", 1.0),
        ("beta", 0.0),
        ("lapse", 1.0),
        ("k_xx", 0.0),
        ("k_xz", 0.0),
        ("k_zz", 0.0),
        ("k_beta", 0.0),
    ):
        if key in state:
            state[key] = set_edge_constant(state[key], value)

    if "phi" in state:
        state["phi"] = set_edge_constant(state["phi"], 1.0)
    if "pi_phi" in state:
        state["pi_phi"] = set_edge_constant(state["pi_phi"], 0.0)


def exact_beamlike_boundary_fields(t: float, params: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    beam = default_params()
    beam.m = params["m"]
    beam.alpha = params["alpha"]
    beam.k0 = params["k0"]
    beam.w0 = params["w0"]
    beam.phi0 = params["phi0"]
    beam.theta = params["theta"]
    beam.x_half_range = params["x_half_range"]
    beam.z_half_range = params["z_half_range"]
    beam.nx = int(params["nx"])
    beam.nz = int(params["nz"])
    data = compute_xz_slice(t, beam)
    rho_ref = data["rho"] + 1.0e-8
    s_ref = unwrap_2d_phase(np.angle(data["psi"]))
    return rho_ref, s_ref


def exact_localized_boundary_fields(t: float, params: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    if params.get("localized_psi0_hat") is not None:
        psi = np.fft.ifft2(params["localized_psi0_hat"] * np.exp(-1j * params["localized_omega"] * float(t)))
    else:
        psi = evolve_positive_frequency(params["localized_psi0"], params["localized_omega"], float(t))
    rho_ref = np.abs(psi) ** 2 + 1.0e-12
    s_ref = phase_field(psi)
    return rho_ref, s_ref


def apply_matter_characteristic_bc(state: dict[str, np.ndarray], t: float, params: dict[str, float]) -> None:
    if params["initial_mode"] == "localized":
        rho_ref, s_ref = exact_localized_boundary_fields(t, params)
    else:
        rho_ref, s_ref = exact_beamlike_boundary_fields(t, params)
    rho_cur, _ = recover_rho_from_n(
        n=state["n"],
        s=state["s"],
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=params["m"],
        dx=params["dx"],
        dz=params["dz"],
    )
    lapse = state.get("lapse", np.ones_like(state["n"]))
    obs = matter_observables(
        rho=rho_cur,
        s=state["s"],
        lapse=lapse,
        shift_x=np.zeros_like(rho_cur),
        shift_z=np.zeros_like(rho_cur),
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=params["m"],
        dx=params["dx"],
        dz=params["dz"],
    )
    vx = obs["flux_x"] / np.maximum(obs["density"], 1.0e-12)
    vz = obs["flux_z"] / np.maximum(obs["density"], 1.0e-12)

    rho_ref_obs = matter_observables(
        rho=rho_ref,
        s=s_ref,
        lapse=lapse,
        shift_x=np.zeros_like(rho_ref),
        shift_z=np.zeros_like(rho_ref),
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=params["m"],
        dx=params["dx"],
        dz=params["dz"],
    )
    n_ref = rho_ref_obs["density"]

    # 左边界 x = x_min：入流当且仅当 v^x > 0
    inflow = vx[0, :] > 0.0
    state["n"][0, inflow] = n_ref[0, inflow]
    state["s"][0, inflow] = s_ref[0, inflow]
    state["n"][0, ~inflow] = state["n"][1, ~inflow]
    state["s"][0, ~inflow] = state["s"][1, ~inflow]

    # 右边界 x = x_max：入流当且仅当 v^x < 0
    inflow = vx[-1, :] < 0.0
    state["n"][-1, inflow] = n_ref[-1, inflow]
    state["s"][-1, inflow] = s_ref[-1, inflow]
    state["n"][-1, ~inflow] = state["n"][-2, ~inflow]
    state["s"][-1, ~inflow] = state["s"][-2, ~inflow]

    # 下边界 z = z_min：入流当且仅当 v^z > 0
    inflow = vz[:, 0] > 0.0
    state["n"][inflow, 0] = n_ref[inflow, 0]
    state["s"][inflow, 0] = s_ref[inflow, 0]
    state["n"][~inflow, 0] = state["n"][~inflow, 1]
    state["s"][~inflow, 0] = state["s"][~inflow, 1]

    # 上边界 z = z_max：入流当且仅当 v^z < 0
    inflow = vz[:, -1] < 0.0
    state["n"][inflow, -1] = n_ref[inflow, -1]
    state["s"][inflow, -1] = s_ref[inflow, -1]
    state["n"][~inflow, -1] = state["n"][~inflow, -2]
    state["s"][~inflow, -1] = state["s"][~inflow, -2]

    # 角点不单独定义特征法向，简单取相邻边界平均。
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


def b_rhs_n(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]

    rho, _ = recover_rho_from_n(
        n=state["n"],
        s=state["s"],
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=mass,
        dx=dx,
        dz=dz,
    )
    lapse = state.get("lapse", np.ones_like(rho))
    zeros = np.zeros_like(rho)

    src = b_branch_sources_from_ns(
        rho=rho,
        s=state["s"],
        lapse=lapse,
        shift_x=zeros,
        shift_z=zeros,
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=mass,
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
        lapse=lapse,
    )
    mrhs = matter_rhs(
        rho=rho,
        s=state["s"],
        lapse=lapse,
        shift_x=zeros,
        shift_z=zeros,
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=mass,
        dx=dx,
        dz=dz,
    )
    lapse_t = np.zeros_like(lapse)
    if params.get("gauge_mode", "synchronous") == "1plog":
        hxx_inv, hxz_inv, hzz_inv, _ = inverse_2metric(state["h_xx"], state["h_xz"], state["h_zz"])
        k_trace = hxx_inv * state["k_xx"] + 2.0 * hxz_inv * state["k_xz"] + hzz_inv * state["k_zz"] + state["k_beta"]
        lapse_t = -2.0 * lapse * k_trace
    return {
        **geo,
        "n_t": mrhs["n_t"],
        "s_t": mrhs["s_t"],
        "lapse_t": lapse_t,
        "rho_current": rho,
        "energy_current": src["energy"],
    }


def c_rhs_n(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]
    ell = params["ell"]

    rho, _ = recover_rho_from_n(
        n=state["n"],
        s=state["s"],
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=mass,
        dx=dx,
        dz=dz,
    )
    lapse = state.get("lapse", np.ones_like(rho))
    zeros = np.zeros_like(rho)
    src = c_branch_effective_sources(
        rho=rho,
        s=state["s"],
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
        lapse=lapse,
        shift_x=zeros,
        shift_z=zeros,
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
        lapse=lapse,
    )
    mrhs = matter_rhs(
        rho=rho,
        s=state["s"],
        lapse=lapse,
        shift_x=zeros,
        shift_z=zeros,
        h_xx=state["h_xx"],
        h_xz=state["h_xz"],
        h_zz=state["h_zz"],
        beta=state["beta"],
        mass=mass,
        dx=dx,
        dz=dz,
    )
    lapse_t = np.zeros_like(lapse)
    if params.get("gauge_mode", "synchronous") == "1plog":
        lapse_t = -2.0 * lapse * src["k_trace"]
    return {
        **geo,
        "n_t": mrhs["n_t"],
        "s_t": mrhs["s_t"],
        "lapse_t": lapse_t,
        "phi_t": state["pi_phi"],
        "pi_phi_t": src["lap3_phi"] + src["k_trace"] * state["pi_phi"] + src["box_phi_corr"],
        "rho_current": rho,
        "energy_current": src["energy"],
    }


def advance_rk4(
    state: dict[str, np.ndarray],
    evolve_keys: tuple[str, ...],
    rhs_fn,
    dt: float,
    params: dict[str, float],
    t_next: float,
    do_finite_check: bool = True,
) -> dict[str, np.ndarray]:
    if do_finite_check:
        assert_finite_state(state, "步进前")
    k1 = rhs_fn(state, params)
    trial2 = state_with_increment(state, evolve_keys, k1, 0.5 * dt)
    if do_finite_check:
        assert_finite_state(trial2, "trial2")
    k2 = rhs_fn(trial2, params)
    trial3 = state_with_increment(state, evolve_keys, k2, 0.5 * dt)
    if do_finite_check:
        assert_finite_state(trial3, "trial3")
    k3 = rhs_fn(trial3, params)
    trial4 = state_with_increment(state, evolve_keys, k3, dt)
    if do_finite_check:
        assert_finite_state(trial4, "trial4")
    k4 = rhs_fn(trial4, params)
    rk4_update_in_place(state, evolve_keys, dt, k1, k2, k3, k4)
    apply_matter_characteristic_bc(state, t_next, params)
    apply_damping_bc(state)
    if do_finite_check:
        assert_finite_state(state, "步进后")
    return rhs_fn(state, params)


def run_case(
    output_dir: Path,
    steps: int = 10,
    dt: float = 5.0e-4,
    c_phi_mode: str = "flat",
    ell: float = 0.02,
    elliptic_lambda: float = 1.0,
    initial_mode: str = "localized",
    finite_check_every: int = 1,
    gauge_mode: str = "synchronous",
) -> dict[str, object]:
    mass = 1.0
    mp = 300.0
    if initial_mode == "localized":
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
        psi0 = initial_wavefunction(X, Z, loc)
        rho0 = np.abs(psi0) ** 2 + loc.rho_floor
        s0 = phase_field(psi0)
        x_half_range = loc.x_half_range
        z_half_range = loc.z_half_range
        nx = loc.nx
        nz = loc.nz
        alpha = loc.alpha
        k0 = loc.k0
        phi0 = loc.phi0
        w0 = np.nan
        theta = np.pi / 4.0
        localized_omega = spectral_omega(x, z, loc.m)
        localized_psi0_hat = np.fft.fft2(psi0)
    elif initial_mode == "beamlike":
        nx = nz = 161
        x_half_range = z_half_range = 5.0
        alpha = 0.5
        k0 = 10.0
        w0 = 1.0
        phi0 = 0.0
        theta = np.pi / 4.0
        x, z, _, _, rho0, s0 = initial_beamlike_fields(x_half_range, z_half_range, nx, nz, alpha, k0, w0, phi0)
        localized_omega = None
        psi0 = None
        localized_psi0_hat = None
    else:
        raise ValueError(f"未知的 initial_mode: {initial_mode}")
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    _, a_src0 = make_a_flat_initial_data(rho0, s0, mass, dx, dz)
    a_weak0 = solve_linearized_conformal_omega(a_src0["energy"], mp, dx, dz)
    omega_init = a_weak0["omega"]

    b_state, c_state = make_bc_initial_states(
        rho0,
        s0,
        omega_init,
        mass,
        mp,
        ell,
        dx,
        dz,
        c_phi_mode=c_phi_mode,
        elliptic_lambda=elliptic_lambda,
    )
    phi_init_minus_1_max = float(np.max(np.abs(c_state["phi"] - 1.0)))
    params = {
        "dx": dx,
        "dz": dz,
        "m": mass,
        "mp": mp,
        "ell": ell,
        "gauge_mode": gauge_mode,
        "alpha": alpha,
        "k0": k0,
        "w0": w0,
        "phi0": phi0,
        "theta": theta,
        "x_half_range": x_half_range,
        "z_half_range": z_half_range,
        "nx": nx,
        "nz": nz,
        "initial_mode": initial_mode,
        "localized_psi0": psi0,
        "localized_omega": localized_omega,
        "localized_psi0_hat": localized_psi0_hat,
    }

    b_ham_history: list[float] = []
    c_ham_history: list[float] = []
    b_mom_history: list[float] = []
    c_mom_history: list[float] = []
    b_r3_history: list[float] = []
    c_r3_history: list[float] = []
    b_boundary_history: list[float] = []
    c_boundary_history: list[float] = []

    b_last = b_rhs_n(b_state, params)
    c_last = c_rhs_n(c_state, params)
    init_b_ham = float(np.max(np.abs(hamiltonian_residual(
        b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"],
        b_state["k_xx"], b_state["k_xz"], b_state["k_zz"], b_state["k_beta"],
        b_last["energy_current"], mp, dx, dz
    ))))
    init_c_ham = float(np.max(np.abs(hamiltonian_residual(
        c_state["h_xx"], c_state["h_xz"], c_state["h_zz"], c_state["beta"],
        c_state["k_xx"], c_state["k_xz"], c_state["k_zz"], c_state["k_beta"],
        c_last["energy_current"], mp, dx, dz
    ))))
    init_b_src_full = b_branch_sources_from_ns(
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
    init_b_mom_x, init_b_mom_z = momentum_constraint_residual_ykilling(
        b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"],
        b_state["k_xx"], b_state["k_xz"], b_state["k_zz"], b_state["k_beta"],
        init_b_src_full["mom_x"],
        init_b_src_full["mom_z"],
        mp,
        dx,
        dz,
    )
    init_c_src_full = c_branch_effective_sources(
        rho=c_last["rho_current"],
        s=c_state["s"],
        phi=c_state["phi"],
        pi_phi=c_state["pi_phi"],
        h_xx=c_state["h_xx"],
        h_xz=c_state["h_xz"],
        h_zz=c_state["h_zz"],
        beta=c_state["beta"],
        k_xx=c_state["k_xx"],
        k_xz=c_state["k_xz"],
        k_zz=c_state["k_zz"],
        k_beta=c_state["k_beta"],
        lapse=np.ones_like(c_last["rho_current"]),
        shift_x=np.zeros_like(c_last["rho_current"]),
        shift_z=np.zeros_like(c_last["rho_current"]),
        mass=mass,
        mp=mp,
        ell=ell,
        dx=dx,
        dz=dz,
    )
    init_c_mom_x, init_c_mom_z = momentum_constraint_residual_ykilling(
        c_state["h_xx"], c_state["h_xz"], c_state["h_zz"], c_state["beta"],
        c_state["k_xx"], c_state["k_xz"], c_state["k_zz"], c_state["k_beta"],
        init_c_src_full["mom_x"],
        init_c_src_full["mom_z"],
        mp,
        dx,
        dz,
    )
    init_b_mom = float(max(np.max(np.abs(init_b_mom_x)), np.max(np.abs(init_b_mom_z))))
    init_c_mom = float(max(np.max(np.abs(init_c_mom_x)), np.max(np.abs(init_c_mom_z))))

    t = 0.0
    b_evolve_keys = B_STATE_KEYS + (("lapse",) if gauge_mode == "1plog" else ())
    c_evolve_keys = C_STATE_KEYS + (("lapse",) if gauge_mode == "1plog" else ())
    for _ in range(steps):
        t_next = t + dt
        step_index = int(round(t / dt)) + 1
        do_finite_check = finite_check_every > 0 and (step_index % finite_check_every == 0 or step_index == 1)
        b_last = advance_rk4(b_state, b_evolve_keys, b_rhs_n, dt, params, t_next, do_finite_check=do_finite_check)
        c_last = advance_rk4(c_state, c_evolve_keys, c_rhs_n, dt, params, t_next, do_finite_check=do_finite_check)
        b_boundary_history.append(boundary_ratio(b_last["rho_current"]))
        c_boundary_history.append(boundary_ratio(c_last["rho_current"]))
        b_r3_history.append(float(np.max(np.abs(b_last["r3_scalar"]))))
        c_r3_history.append(float(np.max(np.abs(c_last["r3_scalar"]))))
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
        b_mom_x_res, b_mom_z_res = momentum_constraint_residual_ykilling(
            b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"],
            b_state["k_xx"], b_state["k_xz"], b_state["k_zz"], b_state["k_beta"],
            b_src_full["mom_x"], b_src_full["mom_z"], mp, dx, dz
        )
        b_mom_history.append(float(max(np.max(np.abs(b_mom_x_res)), np.max(np.abs(b_mom_z_res)))))
        b_ham_history.append(float(np.max(np.abs(hamiltonian_residual(
            b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"],
            b_state["k_xx"], b_state["k_xz"], b_state["k_zz"], b_state["k_beta"],
            b_last["energy_current"], mp, dx, dz
        )))))
        c_src_full = c_branch_effective_sources(
            rho=c_last["rho_current"],
            s=c_state["s"],
            phi=c_state["phi"],
            pi_phi=c_state["pi_phi"],
            h_xx=c_state["h_xx"],
            h_xz=c_state["h_xz"],
            h_zz=c_state["h_zz"],
            beta=c_state["beta"],
            k_xx=c_state["k_xx"],
            k_xz=c_state["k_xz"],
            k_zz=c_state["k_zz"],
            k_beta=c_state["k_beta"],
            lapse=np.ones_like(c_last["rho_current"]),
            shift_x=np.zeros_like(c_last["rho_current"]),
            shift_z=np.zeros_like(c_last["rho_current"]),
            mass=mass,
            mp=mp,
            ell=ell,
            dx=dx,
            dz=dz,
        )
        c_mom_x_res, c_mom_z_res = momentum_constraint_residual_ykilling(
            c_state["h_xx"], c_state["h_xz"], c_state["h_zz"], c_state["beta"],
            c_state["k_xx"], c_state["k_xz"], c_state["k_zz"], c_state["k_beta"],
            c_src_full["mom_x"], c_src_full["mom_z"], mp, dx, dz
        )
        c_mom_history.append(float(max(np.max(np.abs(c_mom_x_res)), np.max(np.abs(c_mom_z_res)))))
        c_ham_history.append(float(np.max(np.abs(hamiltonian_residual(
            c_state["h_xx"], c_state["h_xz"], c_state["h_zz"], c_state["beta"],
            c_state["k_xx"], c_state["k_xz"], c_state["k_zz"], c_state["k_beta"],
            c_last["energy_current"], mp, dx, dz
        )))))
        t = t_next

    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "grid": {"nx": nx, "nz": nz, "x_half_range": x_half_range, "z_half_range": z_half_range},
        "initial_mode": initial_mode,
        "beam_params": {"m": mass, "alpha": alpha, "k0": k0, "w0": w0, "phi0": phi0, "theta_deg": 45.0},
        "evolution": {
            "steps": steps,
            "dt": dt,
            "mp": mp,
            "ell": ell,
            "gauge_mode": gauge_mode,
            "c_phi_mode": c_phi_mode,
            "elliptic_lambda": elliptic_lambda,
            "finite_check_every": finite_check_every,
        },
        "a_initial_weakfield": {"omega_max_init": float(np.max(np.abs(omega_init)))},
        "b_geometry": {
            "initial_hamiltonian_residual_max": init_b_ham,
            "initial_momentum_residual_max": init_b_mom,
            "max_abs_r3_history": b_r3_history,
            "max_abs_hamiltonian_residual": b_ham_history,
            "max_abs_momentum_residual": b_mom_history,
            "final_h_xx_minus_1_max": float(np.max(np.abs(b_state["h_xx"] - 1.0))),
            "final_h_xz_max": float(np.max(np.abs(b_state["h_xz"]))),
            "final_h_zz_minus_1_max": float(np.max(np.abs(b_state["h_zz"] - 1.0))),
            "final_beta_max": float(np.max(np.abs(b_state["beta"]))),
            "final_lapse_min": float(np.min(b_state["lapse"])),
            "final_lapse_max": float(np.max(b_state["lapse"])),
            "final_rho_max": float(np.max(b_last["rho_current"])),
            "max_boundary_ratio": float(max(b_boundary_history) if b_boundary_history else boundary_ratio(b_last["rho_current"])),
        },
        "c_geometry": {
            "initial_hamiltonian_residual_max": init_c_ham,
            "initial_momentum_residual_max": init_c_mom,
            "initial_phi_minus_1_max": phi_init_minus_1_max,
            "max_abs_r3_history": c_r3_history,
            "max_abs_hamiltonian_residual": c_ham_history,
            "max_abs_momentum_residual": c_mom_history,
            "final_h_xx_minus_1_max": float(np.max(np.abs(c_state["h_xx"] - 1.0))),
            "final_h_xz_max": float(np.max(np.abs(c_state["h_xz"]))),
            "final_h_zz_minus_1_max": float(np.max(np.abs(c_state["h_zz"] - 1.0))),
            "final_beta_max": float(np.max(np.abs(c_state["beta"]))),
            "final_lapse_min": float(np.min(c_state["lapse"])),
            "final_lapse_max": float(np.max(c_state["lapse"])),
            "final_phi_minus_1_max": float(np.max(np.abs(c_state["phi"] - 1.0))),
            "final_rho_max": float(np.max(c_last["rho_current"])),
            "max_boundary_ratio": float(max(c_boundary_history) if c_boundary_history else boundary_ratio(c_last["rho_current"])),
        },
        "bc_difference": {
            "final_max_abs_h_xx_diff": float(np.max(np.abs(c_state["h_xx"] - b_state["h_xx"]))),
            "final_max_abs_h_xz_diff": float(np.max(np.abs(c_state["h_xz"] - b_state["h_xz"]))),
            "final_max_abs_h_zz_diff": float(np.max(np.abs(c_state["h_zz"] - b_state["h_zz"]))),
            "final_max_abs_beta_diff": float(np.max(np.abs(c_state["beta"] - b_state["beta"]))),
            "final_max_abs_k_xx_diff": float(np.max(np.abs(c_state["k_xx"] - b_state["k_xx"]))),
            "final_max_abs_k_xz_diff": float(np.max(np.abs(c_state["k_xz"] - b_state["k_xz"]))),
            "final_max_abs_k_zz_diff": float(np.max(np.abs(c_state["k_zz"] - b_state["k_zz"]))),
            "final_max_abs_k_beta_diff": float(np.max(np.abs(c_state["k_beta"] - b_state["k_beta"]))),
            "final_max_abs_n_diff": float(np.max(np.abs(c_state["n"] - b_state["n"]))),
            "final_max_abs_s_diff": float(np.max(np.abs(c_state["s"] - b_state["s"]))),
            "final_max_abs_rho_diff": float(np.max(np.abs(c_last["rho_current"] - b_last["rho_current"]))),
            "final_max_abs_phi_minus_1": float(np.max(np.abs(c_state["phi"] - 1.0))),
            "final_hamiltonian_residual_gap": float(abs(c_ham_history[-1] - b_ham_history[-1])) if b_ham_history and c_ham_history else 0.0,
            "final_momentum_residual_gap": float(abs(c_mom_history[-1] - b_mom_history[-1])) if b_mom_history and c_mom_history else 0.0,
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--dt", type=float, default=5.0e-4)
    parser.add_argument("--c-phi-mode", choices=["flat", "r3_proxy", "elliptic"], default="flat")
    parser.add_argument("--ell", type=float, default=0.02)
    parser.add_argument("--elliptic-lambda", type=float, default=1.0)
    parser.add_argument("--initial-mode", choices=["localized", "beamlike"], default="localized")
    parser.add_argument("--finite-check-every", type=int, default=1)
    parser.add_argument("--gauge-mode", choices=["synchronous", "1plog"], default="synchronous")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "visualizations" / "bc_from_a_initial_data",
    )
    args = parser.parse_args()
    summary = run_case(
        args.output,
        steps=args.steps,
        dt=args.dt,
        c_phi_mode=args.c_phi_mode,
        ell=args.ell,
        elliptic_lambda=args.elliptic_lambda,
        initial_mode=args.initial_mode,
        finite_check_every=args.finite_check_every,
        gauge_mode=args.gauge_mode,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
