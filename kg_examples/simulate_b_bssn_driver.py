from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adm_matter_2p1 import ddx, ddz, matter_observables, matter_rhs
from adm_sources_abc import b_branch_sources_from_ns
from adm_ykilling_geometry import (
    beta_hessian_and_spatial_ricci,
    inverse_2metric,
    scalar_hessian_ykilling,
)
from simulate_bc_from_a_initial_data import (
    exact_beamlike_boundary_fields,
    exact_localized_boundary_fields,
    make_bc_initial_states,
    recover_rho_from_n,
    set_edge_constant,
)
from mixed_tilde_initial_data import localized_direct_tilde_adm_initial
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    phase_field,
    spectral_omega,
)


def conformal_connection_functions(
    gt_xx: np.ndarray,
    gt_xz: np.ndarray,
    gt_zz: np.ndarray,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    gt_xx_inv, gt_xz_inv, gt_zz_inv, _ = inverse_2metric(gt_xx, gt_xz, gt_zz)
    gamma_x = -(ddx(gt_xx_inv, dx) + ddz(gt_xz_inv, dz))
    gamma_z = -(ddx(gt_xz_inv, dx) + ddz(gt_zz_inv, dz))
    return gamma_x, gamma_z


def adm_to_bssn_state(adm_state: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    h_xx = adm_state["h_xx"]
    h_xz = adm_state["h_xz"]
    h_zz = adm_state["h_zz"]
    beta_geom = adm_state["beta"]
    k_xx = adm_state["k_xx"]
    k_xz = adm_state["k_xz"]
    k_zz = adm_state["k_zz"]
    k_beta = adm_state["k_beta"]

    _, _, _, det_h = inverse_2metric(h_xx, h_xz, h_zz)
    det_gamma = np.clip(det_h * np.exp(2.0 * beta_geom), 1.0e-18, None)
    chi = np.power(det_gamma, -1.0 / 3.0)

    gt_xx = chi * h_xx
    gt_xz = chi * h_xz
    gt_zz = chi * h_zz

    hxx_inv, hxz_inv, hzz_inv, _ = inverse_2metric(h_xx, h_xz, h_zz)
    k_trace = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz + k_beta

    at_xx = chi * (k_xx - h_xx * k_trace / 3.0)
    at_xz = chi * (k_xz - h_xz * k_trace / 3.0)
    at_zz = chi * (k_zz - h_zz * k_trace / 3.0)

    return {
        "chi": chi,
        "gt_xx": gt_xx,
        "gt_xz": gt_xz,
        "gt_zz": gt_zz,
        "k_trace": k_trace,
        "at_xx": at_xx,
        "at_xz": at_xz,
        "at_zz": at_zz,
        "lapse": adm_state.get("lapse", np.ones_like(chi)),
        "shift_x": np.zeros_like(chi),
        "shift_z": np.zeros_like(chi),
        "n": adm_state["n"].copy(),
        "s": adm_state["s"].copy(),
    }


def bssn_to_adm_geometry(
    state: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    chi = np.clip(state["chi"], 1.0e-18, None)
    gt_xx = state["gt_xx"]
    gt_xz = state["gt_xz"]
    gt_zz = state["gt_zz"]

    gt_xx_inv, gt_xz_inv, gt_zz_inv, det_gt = inverse_2metric(gt_xx, gt_xz, gt_zz)
    det_gt = np.clip(det_gt, 1.0e-18, None)
    gtyy = 1.0 / det_gt

    h_xx = gt_xx / chi
    h_xz = gt_xz / chi
    h_zz = gt_zz / chi
    beta_geom = -0.5 * (np.log(chi) + np.log(det_gt))
    h_xx_inv, h_xz_inv, h_zz_inv, _ = inverse_2metric(h_xx, h_xz, h_zz)

    at_xx = state["at_xx"]
    at_xz = state["at_xz"]
    at_zz = state["at_zz"]
    k_trace = state["k_trace"]

    at_trace_2d = gt_xx_inv * at_xx + 2.0 * gt_xz_inv * at_xz + gt_zz_inv * at_zz
    at_yy = -gtyy * at_trace_2d

    k_xx = at_xx / chi + h_xx * k_trace / 3.0
    k_xz = at_xz / chi + h_xz * k_trace / 3.0
    k_zz = at_zz / chi + h_zz * k_trace / 3.0
    k_beta = at_yy / gtyy + k_trace / 3.0

    return {
        "h_xx": h_xx,
        "h_xz": h_xz,
        "h_zz": h_zz,
        "beta_geom": beta_geom,
        "k_xx": k_xx,
        "k_xz": k_xz,
        "k_zz": k_zz,
        "k_beta": k_beta,
        "h_xx_inv": h_xx_inv,
        "h_xz_inv": h_xz_inv,
        "h_zz_inv": h_zz_inv,
        "gt_xx_inv": gt_xx_inv,
        "gt_xz_inv": gt_xz_inv,
        "gt_zz_inv": gt_zz_inv,
        "det_gt": det_gt,
        "gtyy": gtyy,
        "at_yy": at_yy,
    }


def tensor_lie_rhs_2d(
    t_xx: np.ndarray,
    t_xz: np.ndarray,
    t_zz: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    shift_x: np.ndarray,
    shift_z: np.ndarray,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sx_x = ddx(shift_x, dx)
    sx_z = ddz(shift_x, dz)
    sz_x = ddx(shift_z, dx)
    sz_z = ddz(shift_z, dz)

    adv_xx = shift_x * ddx(t_xx, dx) + shift_z * ddz(t_xx, dz)
    adv_xz = shift_x * ddx(t_xz, dx) + shift_z * ddz(t_xz, dz)
    adv_zz = shift_x * ddx(t_zz, dx) + shift_z * ddz(t_zz, dz)

    lie_xx = adv_xx + 2.0 * (t_xx * sx_x + t_xz * sz_x)
    lie_xz = adv_xz + h_xx * sx_z + h_xz * sz_z + h_xz * sx_x + h_zz * sz_x
    lie_zz = adv_zz + 2.0 * (t_xz * sx_z + t_zz * sz_z)
    return lie_xx, lie_xz, lie_zz


def scalar_advection(field: np.ndarray, shift_x: np.ndarray, shift_z: np.ndarray, dx: float, dz: float) -> np.ndarray:
    return shift_x * ddx(field, dx) + shift_z * ddz(field, dz)


def physical_geometry_rhs_shifted(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta_geom: np.ndarray,
    k_xx: np.ndarray,
    k_xz: np.ndarray,
    k_zz: np.ndarray,
    k_beta: np.ndarray,
    lapse: np.ndarray,
    shift_x: np.ndarray,
    shift_z: np.ndarray,
    energy: np.ndarray,
    s_xx: np.ndarray,
    s_xz: np.ndarray,
    s_zz: np.ndarray,
    s_yy_reduced: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    geo = beta_hessian_and_spatial_ricci(h_xx, h_xz, h_zz, beta_geom, dx, dz)
    hxx_inv = geo["hxx_inv"]
    hxz_inv = geo["hxz_inv"]
    hzz_inv = geo["hzz_inv"]
    lapse_hess = scalar_hessian_ykilling(h_xx, h_xz, h_zz, lapse, dx, dz)

    beta_x = ddx(beta_geom, dx)
    beta_z = ddz(beta_geom, dz)
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

    s_trace = hxx_inv * s_xx + 2.0 * hxz_inv * s_xz + hzz_inv * s_zz + s_yy_reduced
    source_common = 0.5 * (s_trace - energy)

    l_h_xx, l_h_xz, l_h_zz = tensor_lie_rhs_2d(
        h_xx, h_xz, h_zz, h_xx, h_xz, h_zz, shift_x, shift_z, dx, dz
    )
    l_k_xx, l_k_xz, l_k_zz = tensor_lie_rhs_2d(
        k_xx, k_xz, k_zz, h_xx, h_xz, h_zz, shift_x, shift_z, dx, dz
    )
    adv_beta = scalar_advection(beta_geom, shift_x, shift_z, dx, dz)
    adv_k_beta = scalar_advection(k_beta, shift_x, shift_z, dx, dz)

    h_xx_t = -2.0 * lapse * k_xx + l_h_xx
    h_xz_t = -2.0 * lapse * k_xz + l_h_xz
    h_zz_t = -2.0 * lapse * k_zz + l_h_zz
    beta_geom_t = -lapse * k_beta + adv_beta

    geom_xx = geo["r3_xx"] + k_trace * k_xx - 2.0 * kk_xx - (s_xx - h_xx * source_common) / (mp * mp)
    geom_xz = geo["r3_xz"] + k_trace * k_xz - 2.0 * kk_xz - (s_xz - h_xz * source_common) / (mp * mp)
    geom_zz = geo["r3_zz"] + k_trace * k_zz - 2.0 * kk_zz - (s_zz - h_zz * source_common) / (mp * mp)
    k_xx_t = -lapse_hess["hess_xx"] + lapse * geom_xx + l_k_xx
    k_xz_t = -lapse_hess["hess_xz"] + lapse * geom_xz + l_k_xz
    k_zz_t = -lapse_hess["hess_zz"] + lapse * geom_zz + l_k_zz
    k_beta_t = -grad_beta_dot_grad_lapse + lapse * (
        geo["r3_yy_reduced"] + k_trace * k_beta - (s_yy_reduced - source_common) / (mp * mp)
    ) + adv_k_beta

    return {
        "h_xx_t": h_xx_t,
        "h_xz_t": h_xz_t,
        "h_zz_t": h_zz_t,
        "beta_geom_t": beta_geom_t,
        "k_xx_t": k_xx_t,
        "k_xz_t": k_xz_t,
        "k_zz_t": k_zz_t,
        "k_beta_t": k_beta_t,
        "r3_scalar": geo["r3_scalar"],
    }


def b_bssn_rhs(state: dict[str, np.ndarray], params: dict[str, float]) -> dict[str, np.ndarray]:
    dx = params["dx"]
    dz = params["dz"]
    mass = params["m"]
    mp = params["mp"]
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

    src = b_branch_sources_from_ns(
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
        "rho_current": rho,
        "energy_current": src["energy"],
        "beta_geom_current": geom["beta_geom"],
        "gamma_driver_x": gamma_x,
        "gamma_driver_z": gamma_z,
        "r3_scalar": phys_rhs["r3_scalar"],
    }


def apply_bssn_boundary(state: dict[str, np.ndarray]) -> None:
    for key, value in (
        ("chi", 1.0),
        ("gt_xx", 1.0),
        ("gt_xz", 0.0),
        ("gt_zz", 1.0),
        ("k_trace", 0.0),
        ("at_xx", 0.0),
        ("at_xz", 0.0),
        ("at_zz", 0.0),
        ("lapse", 1.0),
        ("shift_x", 0.0),
        ("shift_z", 0.0),
    ):
        state[key] = set_edge_constant(state[key], value)


def apply_matter_characteristic_bc_bssn(state: dict[str, np.ndarray], t: float, params: dict[str, float]) -> None:
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


def advance_rk4(
    state: dict[str, np.ndarray],
    evolve_keys: tuple[str, ...],
    dt: float,
    params: dict[str, float],
    t_next: float,
    do_finite_check: bool = True,
) -> dict[str, np.ndarray]:
    if do_finite_check:
        assert_finite_state(state, "步进前")
    k1 = b_bssn_rhs(state, params)
    trial2 = state_with_increment(state, evolve_keys, k1, 0.5 * dt)
    if do_finite_check:
        assert_finite_state(trial2, "trial2")
    k2 = b_bssn_rhs(trial2, params)
    trial3 = state_with_increment(state, evolve_keys, k2, 0.5 * dt)
    if do_finite_check:
        assert_finite_state(trial3, "trial3")
    k3 = b_bssn_rhs(trial3, params)
    trial4 = state_with_increment(state, evolve_keys, k3, dt)
    if do_finite_check:
        assert_finite_state(trial4, "trial4")
    k4 = b_bssn_rhs(trial4, params)
    rk4_update_in_place(state, evolve_keys, dt, k1, k2, k3, k4)
    apply_matter_characteristic_bc_bssn(state, t_next, params)
    apply_bssn_boundary(state)
    if do_finite_check:
        assert_finite_state(state, "步进后")
    return b_bssn_rhs(state, params)


def state_metrics(state: dict[str, np.ndarray]) -> dict[str, float]:
    geom = bssn_to_adm_geometry(state)
    return {
        "chi_min": float(np.nanmin(state["chi"])),
        "chi_max": float(np.nanmax(state["chi"])),
        "lapse_min": float(np.nanmin(state["lapse"])),
        "lapse_max": float(np.nanmax(state["lapse"])),
        "shift_abs_max": float(max(np.nanmax(np.abs(state["shift_x"])), np.nanmax(np.abs(state["shift_z"])))),
        "beta_geom_abs_max": float(np.nanmax(np.abs(geom["beta_geom"]))),
        "h_xx_minus_1_max": float(np.nanmax(np.abs(geom["h_xx"] - 1.0))),
        "h_xz_abs_max": float(np.nanmax(np.abs(geom["h_xz"]))),
        "h_zz_minus_1_max": float(np.nanmax(np.abs(geom["h_zz"] - 1.0))),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="B 支 BSSN 型强场原型")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "visualizations" / "b_bssn_driver")
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--dt", type=float, default=2.5e-4)
    parser.add_argument("--finite-check-every", type=int, default=20)
    parser.add_argument("--mu-shift", type=float, default=0.5)
    parser.add_argument("--eta-shift", type=float, default=0.5)
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
    b_adm, _ = make_bc_initial_states(
        rho0,
        s0,
        None,
        mass,
        mp,
        0.2,
        dx,
        dz,
        c_phi_mode="flat",
        tilde_adm_init=direct_init["adm"],
    )
    b_state = adm_to_bssn_state(b_adm)

    params = {
        "dx": dx,
        "dz": dz,
        "m": mass,
        "mp": mp,
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

    evolve_keys = ("chi", "gt_xx", "gt_xz", "gt_zz", "k_trace", "at_xx", "at_xz", "at_zz", "lapse", "shift_x", "shift_z", "n", "s")
    history = []
    t = 0.0
    for step in range(1, args.steps + 1):
        t_next = t + args.dt
        do_check = args.finite_check_every > 0 and (step % args.finite_check_every == 0 or step == 1)
        rhs_out = advance_rk4(b_state, evolve_keys, args.dt, params, t_next, do_finite_check=do_check)
        t = t_next
        if step % 20 == 0 or step == 1:
            history.append({
                "step": step,
                "t": t,
                **state_metrics(b_state),
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
        "final": state_metrics(b_state),
        "history": history,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
