from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adm_matter_2p1 import matter_observables, matter_rhs
from adm_sources_abc import (
    a_branch_sources_from_complex,
    b_branch_sources_from_ns,
    c_branch_effective_sources,
    complex_initial_data_from_rho_s,
)
from adm_initial_data_abc import (
    a_complex_state_from_rho_s_time_symmetric,
    branch_a_sources_from_rho_s_time_symmetric,
    branch_b_sources_from_rho_s,
    solve_average_momentum_extrinsic_curvature,
    solve_time_symmetric_conformal_geometry,
)
from adm_ykilling_geometry import (
    beta_hessian_and_spatial_ricci,
    ddx,
    ddz,
    inverse_2metric,
    scalar_laplacian_ykilling,
    scalar_hessian_ykilling,
)
from simulate_2p1_full_dynamics import initial_beamlike_fields


def edge_damp(field: np.ndarray, coef: float = 0.995, width: int = 8) -> np.ndarray:
    out = field.copy()
    for i in range(width):
        fac = coef ** (width - i)
        out[i, :] *= fac
        out[-1 - i, :] *= fac
        out[:, i] *= fac
        out[:, -1 - i] *= fac
    return out


def set_edge_constant(field: np.ndarray, value: float, width: int = 8) -> np.ndarray:
    out = field.copy()
    out[:width, :] = value
    out[-width:, :] = value
    out[:, :width] = value
    out[:, -width:] = value
    return out


def geometry_rhs(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    k_xx: np.ndarray,
    k_xz: np.ndarray,
    k_zz: np.ndarray,
    k_beta: np.ndarray,
    energy: np.ndarray,
    s_xx: np.ndarray,
    s_xz: np.ndarray,
    s_zz: np.ndarray,
    s_yy_reduced: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
    lapse: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    geo = beta_hessian_and_spatial_ricci(h_xx, h_xz, h_zz, beta, dx, dz)
    hxx_inv = geo["hxx_inv"]
    hxz_inv = geo["hxz_inv"]
    hzz_inv = geo["hzz_inv"]
    if lapse is None:
        lapse = np.ones_like(h_xx)
    lapse_hess = scalar_hessian_ykilling(h_xx, h_xz, h_zz, lapse, dx, dz)
    beta_x = ddx(beta, dx)
    beta_z = ddz(beta, dz)
    grad_beta_dot_grad_lapse = (
        hxx_inv * beta_x * lapse_hess["fx"]
        + hxz_inv * (beta_x * lapse_hess["fz"] + beta_z * lapse_hess["fx"])
        + hzz_inv * beta_z * lapse_hess["fz"]
    )

    k_trace_h = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz
    k_trace = k_trace_h + k_beta

    kx_mixed_x = hxx_inv * k_xx + hxz_inv * k_xz
    kx_mixed_z = hxz_inv * k_xx + hzz_inv * k_xz
    kz_mixed_x = hxx_inv * k_xz + hxz_inv * k_zz
    kz_mixed_z = hxz_inv * k_xz + hzz_inv * k_zz

    kk_xx = k_xx * kx_mixed_x + k_xz * kz_mixed_x
    kk_xz = k_xx * kx_mixed_z + k_xz * kz_mixed_z
    kk_zz = k_xz * kx_mixed_z + k_zz * kz_mixed_z

    s_trace = (
        hxx_inv * s_xx
        + 2.0 * hxz_inv * s_xz
        + hzz_inv * s_zz
        + s_yy_reduced
    )

    source_common = 0.5 * (s_trace - energy)

    h_xx_t = -2.0 * lapse * k_xx
    h_xz_t = -2.0 * lapse * k_xz
    h_zz_t = -2.0 * lapse * k_zz
    beta_t = -lapse * k_beta

    geom_xx = geo["r3_xx"] + k_trace * k_xx - 2.0 * kk_xx - (s_xx - h_xx * source_common) / (mp * mp)
    geom_xz = geo["r3_xz"] + k_trace * k_xz - 2.0 * kk_xz - (s_xz - h_xz * source_common) / (mp * mp)
    geom_zz = geo["r3_zz"] + k_trace * k_zz - 2.0 * kk_zz - (s_zz - h_zz * source_common) / (mp * mp)
    k_xx_t = -lapse_hess["hess_xx"] + lapse * geom_xx
    k_xz_t = -lapse_hess["hess_xz"] + lapse * geom_xz
    k_zz_t = -lapse_hess["hess_zz"] + lapse * geom_zz
    k_beta_t = (
        -grad_beta_dot_grad_lapse
        + lapse * (
            geo["r3_yy_reduced"]
            + k_trace * k_beta
            - (s_yy_reduced - source_common) / (mp * mp)
        )
    )

    return {
        "h_xx_t": h_xx_t,
        "h_xz_t": h_xz_t,
        "h_zz_t": h_zz_t,
        "beta_t": beta_t,
        "k_xx_t": k_xx_t,
        "k_xz_t": k_xz_t,
        "k_zz_t": k_zz_t,
        "k_beta_t": k_beta_t,
        "r3_scalar": geo["r3_scalar"],
    }


def hamiltonian_residual(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    k_xx: np.ndarray,
    k_xz: np.ndarray,
    k_zz: np.ndarray,
    k_beta: np.ndarray,
    energy: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
) -> np.ndarray:
    geo = beta_hessian_and_spatial_ricci(h_xx, h_xz, h_zz, beta, dx, dz)
    hxx_inv = geo["hxx_inv"]
    hxz_inv = geo["hxz_inv"]
    hzz_inv = geo["hzz_inv"]

    k_trace_h = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz
    k_trace = k_trace_h + k_beta
    k_sq = (
        hxx_inv * hxx_inv * k_xx * k_xx
        + 2.0 * (hxx_inv * hxz_inv + hxz_inv * hzz_inv) * k_xx * k_xz
        + 2.0 * (hxz_inv * hxz_inv + hxx_inv * hzz_inv) * k_xz * k_xz
        + 2.0 * hxz_inv * hzz_inv * k_xz * k_zz
        + hzz_inv * hzz_inv * k_zz * k_zz
    )
    k_sq = k_sq + k_beta * k_beta
    return geo["r3_scalar"] + k_trace * k_trace - k_sq - 2.0 * energy / (mp * mp)


def a_rhs(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]

    h_xx = state["h_xx"]
    h_xz = state["h_xz"]
    h_zz = state["h_zz"]
    beta = state["beta"]
    k_xx = state["k_xx"]
    k_xz = state["k_xz"]
    k_zz = state["k_zz"]
    k_beta = state["k_beta"]
    phi1 = state["phi1"]
    phi2 = state["phi2"]
    pi1 = state["pi1"]
    pi2 = state["pi2"]

    src = a_branch_sources_from_complex(phi1, phi2, pi1, pi2, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    geo_rhs = geometry_rhs(
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        k_xx=k_xx,
        k_xz=k_xz,
        k_zz=k_zz,
        k_beta=k_beta,
        energy=src["energy"],
        s_xx=src["s_xx"],
        s_xz=src["s_xz"],
        s_zz=src["s_zz"],
        s_yy_reduced=src["s_yy_reduced"],
        mp=mp,
        dx=dx,
        dz=dz,
    )

    hxx_inv, hxz_inv, hzz_inv, _ = inverse_2metric(h_xx, h_xz, h_zz)
    k_trace = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz + k_beta

    lap_phi1 = scalar_laplacian_ykilling(h_xx, h_xz, h_zz, beta, phi1, dx, dz)
    lap_phi2 = scalar_laplacian_ykilling(h_xx, h_xz, h_zz, beta, phi2, dx, dz)

    return {
        **geo_rhs,
        "phi1_t": pi1,
        "phi2_t": pi2,
        "pi1_t": lap_phi1 + k_trace * pi1 - mass * mass * phi1,
        "pi2_t": lap_phi2 + k_trace * pi2 - mass * mass * phi2,
    }


def b_rhs(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]

    h_xx = state["h_xx"]
    h_xz = state["h_xz"]
    h_zz = state["h_zz"]
    beta = state["beta"]
    k_xx = state["k_xx"]
    k_xz = state["k_xz"]
    k_zz = state["k_zz"]
    k_beta = state["k_beta"]
    rho = state["rho"]
    s = state["s"]

    lapse = state.get("lapse", np.ones_like(rho))
    zero = np.zeros_like(rho)

    src = b_branch_sources_from_ns(rho, s, lapse, zero, zero, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    geo_rhs = geometry_rhs(
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        k_xx=k_xx,
        k_xz=k_xz,
        k_zz=k_zz,
        k_beta=k_beta,
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

    mrhs = matter_rhs(rho, s, lapse, zero, zero, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    hxx_inv, hxz_inv, hzz_inv, _ = inverse_2metric(h_xx, h_xz, h_zz)
    k_trace = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz + k_beta
    lapse_t = np.zeros_like(lapse)
    if params.get("gauge_mode", "synchronous") == "1plog":
        lapse_t = -2.0 * lapse * k_trace

    return {
        **geo_rhs,
        "rho_t": np.zeros_like(rho),  # rho is reconstructed after n update
        "s_t": mrhs["s_t"],
        "n_t": mrhs["n_t"],
        "density": mrhs["density"],
        "energy_matter": mrhs["energy"],
        "sqrt_h": mrhs["sqrt_h"],
        "lapse_t": lapse_t,
    }


def c_rhs(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]
    ell = params["ell"]

    h_xx = state["h_xx"]
    h_xz = state["h_xz"]
    h_zz = state["h_zz"]
    beta = state["beta"]
    k_xx = state["k_xx"]
    k_xz = state["k_xz"]
    k_zz = state["k_zz"]
    k_beta = state["k_beta"]
    rho = state["rho"]
    s = state["s"]
    phi = state["phi"]
    pi_phi = state["pi_phi"]

    lapse = state.get("lapse", np.ones_like(rho))
    zero = np.zeros_like(rho)

    src = c_branch_effective_sources(
        rho=rho,
        s=s,
        phi=phi,
        pi_phi=pi_phi,
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        k_xx=k_xx,
        k_xz=k_xz,
        k_zz=k_zz,
        k_beta=k_beta,
        lapse=lapse,
        shift_x=zero,
        shift_z=zero,
        mass=mass,
        mp=mp,
        ell=ell,
        dx=dx,
        dz=dz,
    )
    geo_rhs = geometry_rhs(
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        k_xx=k_xx,
        k_xz=k_xz,
        k_zz=k_zz,
        k_beta=k_beta,
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

    mrhs = matter_rhs(rho, s, lapse, zero, zero, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    pi_phi_t = src["lap3_phi"] + src["k_trace"] * pi_phi + src["box_phi_corr"]
    lapse_t = np.zeros_like(lapse)
    if params.get("gauge_mode", "synchronous") == "1plog":
        lapse_t = -2.0 * lapse * src["k_trace"]

    return {
        **geo_rhs,
        "rho_t": np.zeros_like(rho),
        "s_t": mrhs["s_t"],
        "n_t": mrhs["n_t"],
        "density": mrhs["density"],
        "energy_matter": mrhs["energy"],
        "sqrt_h": mrhs["sqrt_h"],
        "phi_t": pi_phi,
        "pi_phi_t": pi_phi_t,
        "lapse_t": lapse_t,
    }


def make_initial_states(
    nx: int,
    nz: int,
    x_half_range: float,
    z_half_range: float,
    alpha: float,
    k0: float,
    w0: float,
    phi0: float,
    mass: float,
    mp: float,
    ell: float,
    c_phi_mode: str = "flat",
):
    x, z, _, _, rho0, s0 = initial_beamlike_fields(x_half_range, z_half_range, nx, nz, alpha, k0, w0, phi0)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    shape = rho0.shape
    one = np.ones(shape)
    zero = np.zeros(shape)

    a_flat = branch_a_sources_from_rho_s_time_symmetric(rho0, s0, one, zero, one, zero, mass, dx, dz)
    b_flat = branch_b_sources_from_rho_s(rho0, s0, one, zero, one, zero, mass, dx, dz)
    energy_common = 0.5 * (a_flat["energy"] + b_flat["energy"])
    geo = solve_time_symmetric_conformal_geometry(energy_common, mp, dx, dz)

    h_xx = geo["h_xx"]
    h_xz = geo["h_xz"]
    h_zz = geo["h_zz"]
    beta = geo["beta"]

    a_src_common = branch_a_sources_from_rho_s_time_symmetric(rho0, s0, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    b_src_common = branch_b_sources_from_rho_s(rho0, s0, h_xx, h_xz, h_zz, beta, mass, dx, dz)
    mom_x_common = 0.5 * (a_src_common["mom_x"] + b_src_common["mom_x"])
    mom_z_common = 0.5 * (a_src_common["mom_z"] + b_src_common["mom_z"])
    ksolve = solve_average_momentum_extrinsic_curvature(mom_x_common, mom_z_common, mp, dx, dz)

    a_init = a_src_common

    a_state = {
        "h_xx": h_xx.copy(),
        "h_xz": h_xz.copy(),
        "h_zz": h_zz.copy(),
        "beta": beta.copy(),
        "k_xx": ksolve["k_xx"].copy(),
        "k_xz": ksolve["k_xz"].copy(),
        "k_zz": ksolve["k_zz"].copy(),
        "k_beta": zero.copy(),
        "phi1": a_init["phi1"],
        "phi2": a_init["phi2"],
        "pi1": a_init["pi1"],
        "pi2": a_init["pi2"],
    }

    b_state = {
        "h_xx": h_xx.copy(),
        "h_xz": h_xz.copy(),
        "h_zz": h_zz.copy(),
        "beta": beta.copy(),
        "k_xx": ksolve["k_xx"].copy(),
        "k_xz": ksolve["k_xz"].copy(),
        "k_zz": ksolve["k_zz"].copy(),
        "k_beta": zero.copy(),
        "rho": rho0.copy(),
        "s": s0.copy(),
    }

    c_state = {
        "h_xx": h_xx.copy(),
        "h_xz": h_xz.copy(),
        "h_zz": h_zz.copy(),
        "beta": beta.copy(),
        "k_xx": ksolve["k_xx"].copy(),
        "k_xz": ksolve["k_xz"].copy(),
        "k_zz": ksolve["k_zz"].copy(),
        "k_beta": zero.copy(),
        "rho": rho0.copy(),
        "s": s0.copy(),
        "phi": one.copy(),
        "pi_phi": zero.copy(),
    }

    if c_phi_mode == "r3_proxy":
        r_proxy = np.maximum(geo["r3_scalar"], 0.0)
        c_state["phi"] = np.power(1.0 + np.power((ell * ell) * r_proxy, 2.0), -1.5)
    elif c_phi_mode != "flat":
        raise ValueError(f"未知的 c_phi_mode: {c_phi_mode}")

    meta = {
        "rho0": rho0,
        "s0": s0,
        "energy_a_flat": a_flat["energy"],
        "energy_b_flat": b_flat["energy"],
        "energy_common": energy_common,
        "omega": geo["omega"],
        "r3_scalar": geo["r3_scalar"],
        "ham_constraint_iterations": geo["iterations"],
        "ham_constraint_final_delta": geo["final_delta"],
        "ham_constraint_common_max_abs_residual": float(np.max(np.abs(geo["hamiltonian_residual"]))),
        "mom_constraint_iterations": ksolve["iterations"],
        "mom_constraint_final_delta": ksolve["final_delta"],
        "c_phi_mode": c_phi_mode,
        "c_phi_initial_max": float(np.max(c_state["phi"])),
        "c_phi_initial_min": float(np.min(c_state["phi"])),
    }

    return x, z, a_state, b_state, c_state, meta


def apply_damping(state: dict[str, np.ndarray]) -> None:
    for key, value in state.items():
        state[key] = edge_damp(value)
    if "phi" in state:
        state["phi"] = set_edge_constant(state["phi"], 1.0)
    if "pi_phi" in state:
        state["pi_phi"] = set_edge_constant(state["pi_phi"], 0.0)


def run_case(output_dir: Path, steps: int = 5, dt: float = 0.0025, c_phi_mode: str = "flat") -> dict[str, object]:
    nx = nz = 161
    x_half_range = z_half_range = 5.0
    alpha = 0.5
    k0 = 10.0
    w0 = 1.0
    phi0 = 0.0

    x, z, a_state, b_state, c_state, meta = make_initial_states(
        nx,
        nz,
        x_half_range,
        z_half_range,
        alpha,
        k0,
        w0,
        phi0,
        mass=1.0,
        mp=300.0,
        ell=0.02,
        c_phi_mode=c_phi_mode,
    )
    params = {
        "dx": float(x[1] - x[0]),
        "dz": float(z[1] - z[0]),
        "m": 1.0,
        "mp": 300.0,
        "ell": 0.02,
    }

    a_history: list[float] = []
    b_history: list[float] = []
    c_history: list[float] = []
    a_ham: list[float] = []
    b_ham: list[float] = []
    c_ham: list[float] = []

    a_src0 = a_branch_sources_from_complex(a_state["phi1"], a_state["phi2"], a_state["pi1"], a_state["pi2"], a_state["h_xx"], a_state["h_xz"], a_state["h_zz"], a_state["beta"], params["m"], params["dx"], params["dz"])
    b_src0 = b_branch_sources_from_ns(b_state["rho"], b_state["s"], np.ones_like(b_state["rho"]), np.zeros_like(b_state["rho"]), np.zeros_like(b_state["rho"]), b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"], params["m"], params["dx"], params["dz"])
    c_src0 = c_branch_effective_sources(
        rho=c_state["rho"],
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
        lapse=np.ones_like(c_state["rho"]),
        shift_x=np.zeros_like(c_state["rho"]),
        shift_z=np.zeros_like(c_state["rho"]),
        mass=params["m"],
        mp=params["mp"],
        ell=params["ell"],
        dx=params["dx"],
        dz=params["dz"],
    )
    init_a_ham = float(np.max(np.abs(hamiltonian_residual(a_state["h_xx"], a_state["h_xz"], a_state["h_zz"], a_state["beta"], a_state["k_xx"], a_state["k_xz"], a_state["k_zz"], a_state["k_beta"], a_src0["energy"], params["mp"], params["dx"], params["dz"]))))
    init_b_ham = float(np.max(np.abs(hamiltonian_residual(b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"], b_state["k_xx"], b_state["k_xz"], b_state["k_zz"], b_state["k_beta"], b_src0["energy"], params["mp"], params["dx"], params["dz"]))))
    init_c_ham = float(np.max(np.abs(hamiltonian_residual(c_state["h_xx"], c_state["h_xz"], c_state["h_zz"], c_state["beta"], c_state["k_xx"], c_state["k_xz"], c_state["k_zz"], c_state["k_beta"], c_src0["energy"], params["mp"], params["dx"], params["dz"]))))

    for _ in range(steps):
        ar = a_rhs(a_state, params)
        for key in ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta", "phi1", "phi2", "pi1", "pi2"):
            a_state[key] = a_state[key] + dt * ar[f"{key}_t"]
        apply_damping(a_state)
        a_history.append(float(np.max(np.abs(ar["r3_scalar"]))))
        a_src = a_branch_sources_from_complex(a_state["phi1"], a_state["phi2"], a_state["pi1"], a_state["pi2"], a_state["h_xx"], a_state["h_xz"], a_state["h_zz"], a_state["beta"], params["m"], params["dx"], params["dz"])
        a_ham.append(float(np.max(np.abs(hamiltonian_residual(a_state["h_xx"], a_state["h_xz"], a_state["h_zz"], a_state["beta"], a_state["k_xx"], a_state["k_xz"], a_state["k_zz"], a_state["k_beta"], a_src["energy"], params["mp"], params["dx"], params["dz"])))))

        br = b_rhs(b_state, params)
        for key in ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta", "s"):
            b_state[key] = b_state[key] + dt * br[f"{key}_t"]
        density = br["density"] + dt * br["n_t"]
        b_state["rho"] = np.maximum(density / np.maximum(np.exp(b_state["beta"]) * br["sqrt_h"] * br["energy_matter"], 1.0e-12), 1.0e-12)
        apply_damping(b_state)
        b_history.append(float(np.max(np.abs(br["r3_scalar"]))))
        b_src = b_branch_sources_from_ns(b_state["rho"], b_state["s"], np.ones_like(b_state["rho"]), np.zeros_like(b_state["rho"]), np.zeros_like(b_state["rho"]), b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"], params["m"], params["dx"], params["dz"])
        b_ham.append(float(np.max(np.abs(hamiltonian_residual(b_state["h_xx"], b_state["h_xz"], b_state["h_zz"], b_state["beta"], b_state["k_xx"], b_state["k_xz"], b_state["k_zz"], b_state["k_beta"], b_src["energy"], params["mp"], params["dx"], params["dz"])))))

        cr = c_rhs(c_state, params)
        for key in ("h_xx", "h_xz", "h_zz", "beta", "k_xx", "k_xz", "k_zz", "k_beta", "s", "phi", "pi_phi"):
            c_state[key] = c_state[key] + dt * cr[f"{key}_t"]
        density_c = cr["density"] + dt * cr["n_t"]
        c_state["rho"] = np.maximum(density_c / np.maximum(np.exp(c_state["beta"]) * cr["sqrt_h"] * cr["energy_matter"], 1.0e-12), 1.0e-12)
        apply_damping(c_state)
        c_history.append(float(np.max(np.abs(cr["r3_scalar"]))))
        c_src = c_branch_effective_sources(
            rho=c_state["rho"],
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
            lapse=np.ones_like(c_state["rho"]),
            shift_x=np.zeros_like(c_state["rho"]),
            shift_z=np.zeros_like(c_state["rho"]),
            mass=params["m"],
            mp=params["mp"],
            ell=params["ell"],
            dx=params["dx"],
            dz=params["dz"],
        )
        c_ham.append(float(np.max(np.abs(hamiltonian_residual(c_state["h_xx"], c_state["h_xz"], c_state["h_zz"], c_state["beta"], c_state["k_xx"], c_state["k_xz"], c_state["k_zz"], c_state["k_beta"], c_src["energy"], params["mp"], params["dx"], params["dz"])))))

    output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "steps": steps,
        "dt": dt,
        "initial_geometry_mode": "time_symmetric_conformal_plus_average_momentum_AB",
        "initial_c_phi_mode": meta["c_phi_mode"],
        "initial_c_phi_min": meta["c_phi_initial_min"],
        "initial_c_phi_max": meta["c_phi_initial_max"],
        "initial_ham_constraint_iterations": meta["ham_constraint_iterations"],
        "initial_ham_constraint_final_delta": meta["ham_constraint_final_delta"],
        "initial_common_ham_constraint_max_abs_residual": meta["ham_constraint_common_max_abs_residual"],
        "initial_mom_constraint_iterations": meta["mom_constraint_iterations"],
        "initial_mom_constraint_final_delta": meta["mom_constraint_final_delta"],
        "initial_a_max_abs_hamiltonian_residual": init_a_ham,
        "initial_b_max_abs_hamiltonian_residual": init_b_ham,
        "initial_c_max_abs_hamiltonian_residual": init_c_ham,
        "a_max_abs_r3_history": a_history,
        "b_max_abs_r3_history": b_history,
        "c_max_abs_r3_history": c_history,
        "a_max_abs_hamiltonian_residual": a_ham,
        "b_max_abs_hamiltonian_residual": b_ham,
        "c_max_abs_hamiltonian_residual": c_ham,
        "a_final_phi_norm": float(np.max(np.sqrt(a_state["phi1"] ** 2 + a_state["phi2"] ** 2))),
        "b_final_rho_max": float(np.max(b_state["rho"])),
        "c_final_rho_max": float(np.max(c_state["rho"])),
        "c_final_phi_max": float(np.max(c_state["phi"])),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=5)
    parser.add_argument("--dt", type=float, default=0.0025)
    parser.add_argument("--c-phi-mode", choices=["flat", "r3_proxy"], default="flat")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "visualizations" / "ab_adm_synchronous_prototype")
    args = parser.parse_args()
    summary = run_case(args.output, steps=args.steps, dt=args.dt, c_phi_mode=args.c_phi_mode)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
