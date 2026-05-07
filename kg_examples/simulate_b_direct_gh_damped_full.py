from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coordinate_matter_evolution import coordinate_matter_rhs_covector, conservative_density_from_rho_u
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_initial
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


PAIRS = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]


def d1x(f: np.ndarray, dx: float) -> np.ndarray:
    return np.gradient(f, dx, axis=0, edge_order=2)


def d1z(f: np.ndarray, dz: float) -> np.ndarray:
    return np.gradient(f, dz, axis=1, edge_order=2)


def d2xx(f: np.ndarray, dx: float) -> np.ndarray:
    return np.gradient(np.gradient(f, dx, axis=0, edge_order=2), dx, axis=0, edge_order=2)


def d2zz(f: np.ndarray, dz: float) -> np.ndarray:
    return np.gradient(np.gradient(f, dz, axis=1, edge_order=2), dz, axis=1, edge_order=2)


def d2xz(f: np.ndarray, dx: float, dz: float) -> np.ndarray:
    return np.gradient(np.gradient(f, dx, axis=0, edge_order=2), dz, axis=1, edge_order=2)


def dilate_mask(mask: np.ndarray, rounds: int = 8) -> np.ndarray:
    out = mask.copy()
    for _ in range(rounds):
        out = (
            out
            | np.roll(out, 1, axis=0)
            | np.roll(out, -1, axis=0)
            | np.roll(out, 1, axis=1)
            | np.roll(out, -1, axis=1)
        )
    return out


def metric_jets(metric_cov: np.ndarray, metric_t: np.ndarray, dx: float, dz: float) -> tuple[np.ndarray, np.ndarray]:
    nx, nz = metric_cov.shape[:2]
    dg = np.zeros((nx, nz, 3, 3, 3), dtype=float)
    d2g = np.zeros((nx, nz, 3, 3, 3, 3), dtype=float)

    dg[..., 0, :, :] = metric_t
    dg[..., 1, :, :] = d1x(metric_cov, dx)
    dg[..., 2, :, :] = d1z(metric_cov, dz)

    d2g[..., 0, 1, :, :] = d1x(metric_t, dx)
    d2g[..., 1, 0, :, :] = d2g[..., 0, 1, :, :]
    d2g[..., 0, 2, :, :] = d1z(metric_t, dz)
    d2g[..., 2, 0, :, :] = d2g[..., 0, 2, :, :]
    d2g[..., 1, 1, :, :] = d2xx(metric_cov, dx)
    d2g[..., 1, 2, :, :] = d2xz(metric_cov, dx, dz)
    d2g[..., 2, 1, :, :] = d2g[..., 1, 2, :, :]
    d2g[..., 2, 2, :, :] = d2zz(metric_cov, dz)
    return dg, d2g


def scalar_jets(beta: np.ndarray, beta_t: np.ndarray, dx: float, dz: float) -> tuple[np.ndarray, np.ndarray]:
    db = np.zeros(beta.shape + (3,), dtype=float)
    d2b = np.zeros(beta.shape + (3, 3), dtype=float)
    db[..., 0] = beta_t
    db[..., 1] = d1x(beta, dx)
    db[..., 2] = d1z(beta, dz)
    d2b[..., 0, 1] = d1x(beta_t, dx)
    d2b[..., 1, 0] = d2b[..., 0, 1]
    d2b[..., 0, 2] = d1z(beta_t, dz)
    d2b[..., 2, 0] = d2b[..., 0, 2]
    d2b[..., 1, 1] = d2xx(beta, dx)
    d2b[..., 1, 2] = d2xz(beta, dx, dz)
    d2b[..., 2, 1] = d2b[..., 1, 2]
    d2b[..., 2, 2] = d2zz(beta, dz)
    return db, d2b


def point_reduced_ricci(
    g: np.ndarray,
    dg: np.ndarray,
    d2g: np.ndarray,
    h: np.ndarray,
    dh: np.ndarray,
    damping_kappa: float = 0.0,
    damping_lambda: float = 1.0,
    n_floor: float = 1.0e-12,
) -> np.ndarray:
    ginv = np.linalg.inv(g)

    gamma1 = np.zeros((3, 3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma1[a, b, c] = 0.5 * (dg[b, a, c] + dg[c, a, b] - dg[a, b, c])
    gamma2 = np.einsum("ad,dbc->abc", ginv, gamma1, optimize=True)

    dginv = np.zeros((3, 3, 3), dtype=float)
    for e in range(3):
        dginv[e] = -ginv @ dg[e] @ ginv

    dgamma1 = np.zeros((3, 3, 3, 3), dtype=float)
    for e in range(3):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    dgamma1[e, a, b, c] = 0.5 * (
                        d2g[e, b, a, c] + d2g[e, c, a, b] - d2g[e, a, b, c]
                    )

    dgamma2 = np.zeros((3, 3, 3, 3), dtype=float)
    for e in range(3):
        dgamma2[e] = np.einsum("ad,dbc->abc", dginv[e], gamma1, optimize=True) + np.einsum(
            "ad,dbc->abc", ginv, dgamma1[e], optimize=True
        )

    ricci = np.zeros((3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            term1 = sum(dgamma2[c, c, a, b] for c in range(3))
            term2 = sum(dgamma2[b, c, a, c] for c in range(3))
            quad1 = 0.0
            quad2 = 0.0
            for c in range(3):
                for d in range(3):
                    quad1 += gamma2[c, a, b] * gamma2[d, c, d]
                    quad2 += gamma2[c, a, d] * gamma2[d, b, c]
            ricci[a, b] = term1 - term2 + quad1 - quad2

    gamma_low = np.einsum("bc,abc->a", ginv, gamma1, optimize=True)
    dgamma_low = np.zeros((3, 3), dtype=float)
    for e in range(3):
        dgamma_low[e] = np.einsum("bc,abc->a", dginv[e], gamma1, optimize=True) + np.einsum(
            "bc,abc->a", ginv, dgamma1[e], optimize=True
        )
    cvec = gamma_low - h
    dc = dgamma_low - dh

    red = np.zeros((3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            nabla_sym = 0.5 * (dc[a, b] + dc[b, a])
            for d in range(3):
                nabla_sym -= gamma2[d, a, b] * cvec[d]
            red[a, b] = ricci[a, b] - nabla_sym

    if damping_kappa > 0.0:
        g00_abs = max(abs(float(ginv[0, 0])), n_floor)
        denom = np.sqrt(g00_abs)
        n_cov = np.array([1.0 / denom, 0.0, 0.0], dtype=float)
        n_con = ginv[:, 0] / denom
        ndotc = float(np.dot(n_con, cvec))
        damp = np.zeros((3, 3), dtype=float)
        for a in range(3):
            for b in range(3):
                damp[a, b] = (
                    n_cov[a] * cvec[b]
                    + n_cov[b] * cvec[a]
                    - damping_lambda * g[a, b] * ndotc
                )
        red = red + damping_kappa * damp
    return red


def point_beta_tt(
    g: np.ndarray,
    dg: np.ndarray,
    db: np.ndarray,
    d2b: np.ndarray,
    rho: float,
    mass: float,
    mp: float,
) -> float:
    ginv = np.linalg.inv(g)
    gamma1 = np.zeros((3, 3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma1[a, b, c] = 0.5 * (dg[b, a, c] + dg[c, a, b] - dg[a, b, c])
    gamma2 = np.einsum("ad,dbc->abc", ginv, gamma1, optimize=True)

    gamma_contract = 0.0
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma_contract += ginv[a, b] * gamma2[c, a, b] * db[c]

    rest = (
        2.0 * ginv[0, 1] * d2b[0, 1]
        + 2.0 * ginv[0, 2] * d2b[0, 2]
        + ginv[1, 1] * d2b[1, 1]
        + 2.0 * ginv[1, 2] * d2b[1, 2]
        + ginv[2, 2] * d2b[2, 2]
        - gamma_contract
    )
    rhs = -(mass * mass) * np.exp(-2.0 * db[0] * 0.0) * rho / (mp * mp)
    return (rhs - rest) / ginv[0, 0]


def build_h_source(metric_cov: np.ndarray, metric_t: np.ndarray, dx: float, dz: float) -> tuple[np.ndarray, np.ndarray]:
    dg, _ = metric_jets(metric_cov, metric_t, dx, dz)
    nx, nz = metric_cov.shape[:2]
    h = np.zeros((nx, nz, 3), dtype=float)
    dh = np.zeros((nx, nz, 3, 3), dtype=float)
    for i in range(nx):
        for j in range(nz):
            g = metric_cov[i, j]
            ginv = np.linalg.pinv(g, rcond=1.0e-12, hermitian=True)
            gamma1 = np.zeros((3, 3, 3), dtype=float)
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        gamma1[a, b, c] = 0.5 * (dg[i, j, b, a, c] + dg[i, j, c, a, b] - dg[i, j, a, b, c])
            h[i, j] = np.einsum("bc,abc->a", ginv, gamma1, optimize=True)
    dh[..., 1, :] = d1x(h, dx)
    dh[..., 2, :] = d1z(h, dz)
    return h, dh


def metric_second_time_rhs(
    metric_cov: np.ndarray,
    metric_t: np.ndarray,
    beta: np.ndarray,
    beta_t: np.ndarray,
    rho: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    h_source: np.ndarray,
    dh_source: np.ndarray,
    mp: float,
    mass: float,
    dx: float,
    dz: float,
    active_mask: np.ndarray,
    damping_kappa: float = 0.0,
    damping_lambda: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    nx, nz = metric_cov.shape[:2]
    dg, d2g = metric_jets(metric_cov, metric_t, dx, dz)
    db, d2b = scalar_jets(beta, beta_t, dx, dz)

    g_tt = np.zeros_like(metric_cov)
    beta_tt = np.zeros_like(beta)
    for i in range(nx):
        for j in range(nz):
            if not active_mask[i, j]:
                continue
            g = metric_cov[i, j]
            ginv = np.linalg.inv(g)
            gradb = db[i, j]
            beta_sq = float(ginv @ gradb @ gradb)
            u = np.array([u_t[i, j], u_x[i, j], u_z[i, j]], dtype=float)
            t_beta = 2.0 * mp * mp * (np.outer(gradb, gradb) - 0.5 * g * beta_sq)
            t_m = 2.0 * rho[i, j] * np.outer(u, u)
            t_tot = t_beta + t_m
            trace_t = float(np.einsum("ab,ab->", ginv, t_tot, optimize=True))
            source = (t_tot - g * trace_t) / (mp * mp)

            base_d2g = d2g[i, j].copy()
            base_d2g[0, 0, :, :] = 0.0
            red0 = point_reduced_ricci(
                g,
                dg[i, j],
                base_d2g,
                h_source[i, j],
                dh_source[i, j],
                damping_kappa=damping_kappa,
                damping_lambda=damping_lambda,
            )
            rhs0 = source - red0

            a_mat = np.zeros((6, 6), dtype=float)
            for col, (a, b) in enumerate(PAIRS):
                trial = base_d2g.copy()
                trial[0, 0, a, b] = 1.0
                trial[0, 0, b, a] = 1.0
                red1 = point_reduced_ricci(
                    g,
                    dg[i, j],
                    trial,
                    h_source[i, j],
                    dh_source[i, j],
                    damping_kappa=damping_kappa,
                    damping_lambda=damping_lambda,
                )
                delta = red1 - red0
                for row, (c, d) in enumerate(PAIRS):
                    a_mat[row, col] = delta[c, d]

            b_vec = np.array([rhs0[a, b] for (a, b) in PAIRS], dtype=float)
            try:
                sol = np.linalg.solve(a_mat, b_vec)
            except np.linalg.LinAlgError:
                sol = np.zeros(6, dtype=float)

            for val, (a, b) in zip(sol, PAIRS):
                g_tt[i, j, a, b] = val
                g_tt[i, j, b, a] = val

            beta_tt[i, j] = point_beta_tt(g, dg[i, j], db[i, j], d2b[i, j], float(rho[i, j]), mass, mp)
    return g_tt, beta_tt


def state_zero_like(state: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {k: np.zeros_like(v) for k, v in state.items()}


def state_add_scaled(state: dict[str, np.ndarray], delta: dict[str, np.ndarray], scale: float) -> dict[str, np.ndarray]:
    out = {}
    for k, v in state.items():
        if k == "ut":
            out[k] = v.copy()
        else:
            out[k] = v + scale * delta[k]
    return out


def project_state(
    state: dict[str, np.ndarray],
    active: np.ndarray,
    eta_cov: np.ndarray,
    mass: float,
) -> dict[str, np.ndarray]:
    out = {k: np.array(v, copy=True) for k, v in state.items()}
    out["metric"] = np.where(active[..., None, None], out["metric"], eta_cov)
    out["pi"] = np.where(active[..., None, None], out["pi"], 0.0)
    out["phix"] = np.where(active[..., None, None], out["phix"], 0.0)
    out["phiz"] = np.where(active[..., None, None], out["phiz"], 0.0)
    out["beta"] = np.where(active, out["beta"], 0.0)
    out["pi_beta"] = np.where(active, out["pi_beta"], 0.0)
    out["phix_beta"] = np.where(active, out["phix_beta"], 0.0)
    out["phiz_beta"] = np.where(active, out["phiz_beta"], 0.0)
    out["n_cons"] = np.where(active, out["n_cons"], 0.0)
    out["ux"] = np.where(active, out["ux"], 0.0)
    out["uz"] = np.where(active, out["uz"], 0.0)
    out["ut"] = np.where(active, out["ut"], -mass)
    return out


def gamma_contract_scalar(ginv: np.ndarray, dg: np.ndarray, db: np.ndarray) -> float:
    gamma1 = np.zeros((3, 3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma1[a, b, c] = 0.5 * (dg[b, a, c] + dg[c, a, b] - dg[a, b, c])
    gamma2 = np.einsum("ad,dbc->abc", ginv, gamma1, optimize=True)
    gamma_contract = 0.0
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma_contract += ginv[a, b] * gamma2[c, a, b] * db[c]
    return float(gamma_contract)


def rhs_full_gh(
    state: dict[str, np.ndarray],
    h_source: np.ndarray,
    dh_source: np.ndarray,
    active: np.ndarray,
    dx: float,
    dz: float,
    mp: float,
    mass: float,
    damping_kappa: float,
    damping_lambda: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    rhs = state_zero_like(state)

    metric = state["metric"]
    pi = state["pi"]
    phix = state["phix"]
    phiz = state["phiz"]
    beta = state["beta"]
    pi_beta = state["pi_beta"]
    phix_beta = state["phix_beta"]
    phiz_beta = state["phiz_beta"]

    mass_eff = mass * np.exp(-beta)
    mrhs = coordinate_matter_rhs_covector(
        n_cons=state["n_cons"],
        u_x=state["ux"],
        u_z=state["uz"],
        metric_cov_txz=metric,
        mass=mass_eff,
        dx=dx,
        dz=dz,
        u_t_reference=state["ut"],
    )

    dx_pi = d1x(pi, dx)
    dz_pi = d1z(pi, dz)
    dxx_g = d1x(phix, dx)
    dzz_g = d1z(phiz, dz)
    dx_phiz = d1x(phiz, dx)
    dz_phix = d1z(phix, dz)
    cross_g = 0.5 * (dx_phiz + dz_phix)

    dx_pi_beta = d1x(pi_beta, dx)
    dz_pi_beta = d1z(pi_beta, dz)
    dxx_beta = d1x(phix_beta, dx)
    dzz_beta = d1z(phiz_beta, dz)
    dx_phiz_beta = d1x(phiz_beta, dx)
    dz_phix_beta = d1z(phix_beta, dz)
    cross_beta = 0.5 * (dx_phiz_beta + dz_phix_beta)

    rhs["metric"] = pi
    rhs["phix"] = dx_pi
    rhs["phiz"] = dz_pi
    rhs["beta"] = pi_beta
    rhs["phix_beta"] = dx_pi_beta
    rhs["phiz_beta"] = dz_pi_beta
    rhs["n_cons"] = np.where(active, mrhs["n_t"], 0.0)
    rhs["ux"] = np.where(active, mrhs["u_x_t"], 0.0)
    rhs["uz"] = np.where(active, mrhs["u_z_t"], 0.0)
    rhs["ut"] = np.zeros_like(state["ut"])

    zero_d2g = np.zeros((3, 3, 3, 3), dtype=float)
    for i in range(metric.shape[0]):
        for j in range(metric.shape[1]):
            if not active[i, j]:
                continue
            g = metric[i, j]
            ginv = np.linalg.pinv(g, rcond=1.0e-12, hermitian=True)
            g00 = float(ginv[0, 0])
            g00_safe = g00 if abs(g00) > 1.0e-12 else (1.0e-12 if g00 >= 0.0 else -1.0e-12)

            dg_pt = np.stack([pi[i, j], phix[i, j], phiz[i, j]], axis=0)
            db_pt = np.array([pi_beta[i, j], phix_beta[i, j], phiz_beta[i, j]], dtype=float)
            beta_sq = float(ginv @ db_pt @ db_pt)
            u = np.array([mrhs["u_t"][i, j], state["ux"][i, j], state["uz"][i, j]], dtype=float)
            t_beta = 2.0 * mp * mp * (np.outer(db_pt, db_pt) - 0.5 * g * beta_sq)
            t_m = 2.0 * mrhs["rho"][i, j] * np.outer(u, u)
            t_tot = t_beta + t_m
            trace_t = float(np.einsum("ab,ab->", ginv, t_tot, optimize=True))
            source = (t_tot - g * trace_t) / (mp * mp)

            red0 = point_reduced_ricci(
                g,
                dg_pt,
                zero_d2g,
                h_source[i, j],
                dh_source[i, j],
                damping_kappa=damping_kappa,
                damping_lambda=damping_lambda,
            )

            spatial_term = (
                ginv[1, 1] * dxx_g[i, j]
                + 2.0 * ginv[1, 2] * cross_g[i, j]
                + ginv[2, 2] * dzz_g[i, j]
            )
            mixed_term = 2.0 * ginv[0, 1] * dx_pi[i, j] + 2.0 * ginv[0, 2] * dz_pi[i, j]
            rhs["pi"][i, j] = (-2.0 * (source - red0) - mixed_term - spatial_term) / g00_safe

            gamma_contract = gamma_contract_scalar(ginv, dg_pt, db_pt)
            beta_rhs = -(mass * mass) * np.exp(-2.0 * beta[i, j]) * mrhs["rho"][i, j] / (mp * mp)
            beta_spatial = (
                2.0 * ginv[0, 1] * dx_pi_beta[i, j]
                + 2.0 * ginv[0, 2] * dz_pi_beta[i, j]
                + ginv[1, 1] * dxx_beta[i, j]
                + 2.0 * ginv[1, 2] * cross_beta[i, j]
                + ginv[2, 2] * dzz_beta[i, j]
                - gamma_contract
            )
            rhs["pi_beta"][i, j] = (beta_rhs - beta_spatial) / g00_safe

    return rhs, mrhs


def run_case(
    output_dir: str | Path,
    nx: int = 64,
    nz: int = 64,
    steps: int = 8,
    dt: float = 2.5e-4,
    mp: float = 300.0,
    active_rel_cut: float = 1.0e-2,
    active_dilate_rounds: int = 6,
    damping_kappa: float = 0.0,
    damping_lambda: float = 1.0,
) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = FlatLocalizedCrossingParams(nx=nx, nz=nz)
    x, z, X, Z = make_grid(params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    init = localized_direct_tilde_coordinate_initial(psi0, psi0_hat, omega, x, z)
    d = init["diagnostics"]

    metric = init["metric_cov_txz_0"].copy()
    pi = init["metric_cov_txz_t0"].copy()
    beta = np.zeros((nx, nz), dtype=float)
    pi_beta = np.zeros((nx, nz), dtype=float)
    ux = d["s_x0"].copy()
    uz = d["s_z0"].copy()
    ut = d["s_t0"].copy()

    rho0 = init["rho0"]
    support = rho0 > active_rel_cut * float(np.max(rho0))
    active = dilate_mask(support, rounds=active_dilate_rounds)

    eta_cov = np.zeros_like(metric)
    eta_cov[..., 0, 0] = 1.0
    eta_cov[..., 1, 1] = -1.0
    eta_cov[..., 2, 2] = -1.0
    state = {
        "metric": metric,
        "pi": pi,
        "phix": d1x(metric, dx),
        "phiz": d1z(metric, dz),
        "beta": beta,
        "pi_beta": pi_beta,
        "phix_beta": np.zeros((nx, nz), dtype=float),
        "phiz_beta": np.zeros((nx, nz), dtype=float),
        "n_cons": conservative_density_from_rho_u(rho0, ux, uz, metric, params.m, u_t_reference=ut)["n_cons"],
        "ux": ux,
        "uz": uz,
        "ut": ut,
    }
    state = project_state(state, active, eta_cov, params.m)
    h_source, dh_source = build_h_source(state["metric"], state["pi"], dx, dz)

    records: list[dict[str, float]] = []
    for step in range(steps):
        k1, aux1 = rhs_full_gh(
            state,
            h_source,
            dh_source,
            active,
            dx,
            dz,
            mp,
            params.m,
            damping_kappa,
            damping_lambda,
        )
        s2 = project_state(state_add_scaled(state, k1, 0.5 * dt), active, eta_cov, params.m)
        k2, aux2 = rhs_full_gh(
            s2,
            h_source,
            dh_source,
            active,
            dx,
            dz,
            mp,
            params.m,
            damping_kappa,
            damping_lambda,
        )
        s3 = project_state(state_add_scaled(state, k2, 0.5 * dt), active, eta_cov, params.m)
        k3, aux3 = rhs_full_gh(
            s3,
            h_source,
            dh_source,
            active,
            dx,
            dz,
            mp,
            params.m,
            damping_kappa,
            damping_lambda,
        )
        s4 = project_state(state_add_scaled(state, k3, dt), active, eta_cov, params.m)
        k4, aux4 = rhs_full_gh(
            s4,
            h_source,
            dh_source,
            active,
            dx,
            dz,
            mp,
            params.m,
            damping_kappa,
            damping_lambda,
        )

        new_state = {}
        for key, arr in state.items():
            if key == "ut":
                new_state[key] = aux4["u_t"]
            else:
                new_state[key] = arr + (dt / 6.0) * (k1[key] + 2.0 * k2[key] + 2.0 * k3[key] + k4[key])
        state = project_state(new_state, active, eta_cov, params.m)
        current = coordinate_matter_rhs_covector(
            n_cons=state["n_cons"],
            u_x=state["ux"],
            u_z=state["uz"],
            metric_cov_txz=state["metric"],
            mass=params.m * np.exp(-state["beta"]),
            dx=dx,
            dz=dz,
            u_t_reference=state["ut"],
        )
        state["ut"] = np.where(active, current["u_t"], -params.m)

        if (step + 1) % 2 == 0 or step == steps - 1:
            records.append(
                {
                    "step": int(step + 1),
                    "t": float((step + 1) * dt),
                    "rho_max": float(np.nanmax(current["rho"])),
                    "disc_min_support": float(np.nanmin(current["discriminant"][support])),
                    "metric_abs_max": float(np.nanmax(np.abs(state["metric"]))),
                    "pi_abs_max": float(np.nanmax(np.abs(state["pi"]))),
                    "beta_abs_max": float(np.nanmax(np.abs(state["beta"]))),
                    "pi_beta_abs_max": float(np.nanmax(np.abs(state["pi_beta"]))),
                    "active_count": int(np.sum(active)),
                }
            )

    summary = {
        "nx": int(nx),
        "nz": int(nz),
        "steps": int(steps),
        "dt": float(dt),
        "t_final": float(steps * dt),
        "damping_kappa": float(damping_kappa),
        "damping_lambda": float(damping_lambda),
        "records": records,
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="固定实验室坐标下的 B 支广义调和约束阻尼一阶原型")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--nz", type=int, default=64)
    parser.add_argument("--steps", type=int, default=8)
    parser.add_argument("--dt", type=float, default=2.5e-4)
    parser.add_argument("--active-rel-cut", type=float, default=1.0e-2)
    parser.add_argument("--active-dilate-rounds", type=int, default=6)
    parser.add_argument("--damping-kappa", type=float, default=0.0)
    parser.add_argument("--damping-lambda", type=float, default=1.0)
    args = parser.parse_args()
    summary = run_case(
        args.output_dir,
        nx=args.nx,
        nz=args.nz,
        steps=args.steps,
        dt=args.dt,
        active_rel_cut=args.active_rel_cut,
        active_dilate_rounds=args.active_dilate_rounds,
        damping_kappa=args.damping_kappa,
        damping_lambda=args.damping_lambda,
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
