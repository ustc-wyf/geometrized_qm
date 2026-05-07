from __future__ import annotations

from pathlib import Path
import json

import numpy as np

from kg_double_packet import (
    default_baseline_params,
    compute_field_and_derivatives,
    compute_bohm_quantities,
)


def inverse_metric_2x2(g_cov: np.ndarray) -> np.ndarray:
    det = g_cov[..., 0, 0] * g_cov[..., 1, 1] - g_cov[..., 0, 1] * g_cov[..., 1, 0]
    inv = np.empty_like(g_cov)
    inv[..., 0, 0] = g_cov[..., 1, 1] / det
    inv[..., 1, 1] = g_cov[..., 0, 0] / det
    inv[..., 0, 1] = -g_cov[..., 0, 1] / det
    inv[..., 1, 0] = -g_cov[..., 1, 0] / det
    return inv


def compute_metric_on_grid(
    x: np.ndarray, t: np.ndarray, m: float | None = None
) -> dict[str, np.ndarray]:
    params = default_baseline_params()
    if m is not None:
        params.m = m

    nt = len(t)
    nx = len(x)

    rho = np.empty((nt, nx), dtype=float)
    Q = np.empty((nt, nx), dtype=float)
    X = np.empty((nt, nx), dtype=float)
    dtS = np.empty((nt, nx), dtype=float)
    dxS = np.empty((nt, nx), dtype=float)

    for i, tt in enumerate(t):
        field = compute_field_and_derivatives(x, float(tt), params)
        bohm = compute_bohm_quantities(field)
        rho[i] = np.real(bohm["rho"])
        Q[i] = np.real(bohm["Q"])
        X[i] = np.real(bohm["X"])
        dtS[i] = np.real(bohm["dtS"])
        dxS[i] = np.real(bohm["dxS"])

    # timelike branch only, where X > 0
    eps = Q / (params.m**2)
    u_cov = np.stack([dtS, dxS], axis=-1)
    g_cov = np.zeros((nt, nx, 2, 2), dtype=float)
    g_cov[..., 0, 0] = 1.0
    g_cov[..., 1, 1] = -1.0

    safe_X = np.where(np.abs(X) > 1e-12, X, np.nan)
    coeff = Q / (safe_X * params.m**2)

    g_tilde = g_cov.copy()
    for a in range(2):
        for b in range(2):
            g_tilde[..., a, b] += coeff * u_cov[..., a] * u_cov[..., b]

    det = g_tilde[..., 0, 0] * g_tilde[..., 1, 1] - g_tilde[..., 0, 1] * g_tilde[..., 1, 0]
    inv = inverse_metric_2x2(g_tilde)

    return {
        "rho": rho,
        "Q": Q,
        "X": X,
        "eps": eps,
        "g_tilde_cov": g_tilde,
        "g_tilde_inv": inv,
        "det_g_tilde": det,
    }


def compute_ricci_scalar_2d(g_cov: np.ndarray, x: np.ndarray, t: np.ndarray) -> np.ndarray:
    g_inv = inverse_metric_2x2(g_cov)
    dt = float(t[1] - t[0])
    dx = float(x[1] - x[0])

    dg = np.empty((2, 2, 2) + g_cov.shape[:2], dtype=float)
    # dg[sigma, mu, nu, it, ix] = ∂_sigma g_{mu nu}
    for mu in range(2):
        for nu in range(2):
            dg[0, mu, nu] = np.gradient(g_cov[..., mu, nu], dt, axis=0)
            dg[1, mu, nu] = np.gradient(g_cov[..., mu, nu], dx, axis=1)

    Gamma = np.zeros(g_cov.shape[:2] + (2, 2, 2), dtype=float)
    # Gamma[..., rho, mu, nu]
    for rho in range(2):
        for mu in range(2):
            for nu in range(2):
                acc = 0.0
                for sigma in range(2):
                    acc += g_inv[..., rho, sigma] * (
                        dg[mu, sigma, nu] + dg[nu, sigma, mu] - dg[sigma, mu, nu]
                    )
                Gamma[..., rho, mu, nu] = 0.5 * acc

    dGamma = np.empty((2,) + Gamma.shape, dtype=float)
    for rho in range(2):
        dGamma[0, ..., rho, :, :] = np.gradient(Gamma[..., rho, :, :], dt, axis=0)
        dGamma[1, ..., rho, :, :] = np.gradient(Gamma[..., rho, :, :], dx, axis=1)

    Ric = np.zeros(g_cov.shape[:2] + (2, 2), dtype=float)
    for mu in range(2):
        for nu in range(2):
            term = np.zeros(g_cov.shape[:2], dtype=float)
            for rho in range(2):
                term += dGamma[rho, ..., rho, mu, nu]
                term -= dGamma[nu, ..., rho, mu, rho]
                for lam in range(2):
                    term += Gamma[..., rho, rho, lam] * Gamma[..., lam, mu, nu]
                    term -= Gamma[..., rho, nu, lam] * Gamma[..., lam, mu, rho]
            Ric[..., mu, nu] = term

    R = np.zeros(g_cov.shape[:2], dtype=float)
    for mu in range(2):
        for nu in range(2):
            R += g_inv[..., mu, nu] * Ric[..., mu, nu]
    return R


def run_curvature_case(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = default_baseline_params()
    x = np.linspace(-18.0, 18.0, 901)
    t0 = 9.583148474999101
    t = np.linspace(t0 - 2.0, t0 + 2.0, 121)

    geom = compute_metric_on_grid(x, t)
    R_tilde = compute_ricci_scalar_2d(geom["g_tilde_cov"], x, t)

    rho = geom["rho"]
    eps = geom["eps"]
    mask_support = rho >= rho.max() * 1e-2
    mask_loose = rho >= rho.max() * 1e-6

    payload = {
        "window": {
            "x_min": float(x.min()),
            "x_max": float(x.max()),
            "t_min": float(t.min()),
            "t_max": float(t.max()),
        },
        "summary": {
            "R_tilde_abs_max_all": float(np.nanmax(np.abs(R_tilde))),
            "R_tilde_abs_max_support_1e-2": float(np.nanmax(np.abs(R_tilde[mask_support]))),
            "R_tilde_abs_max_support_1e-6": float(np.nanmax(np.abs(R_tilde[mask_loose]))),
            "eps_abs_max_support_1e-2": float(np.nanmax(np.abs(eps[mask_support]))),
            "eps_abs_max_support_1e-6": float(np.nanmax(np.abs(eps[mask_loose]))),
            "det_sign_flip_count": int(np.sum(geom["det_g_tilde"] >= 0.0)),
        },
    }

    (out / "curvature_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    np.savez_compressed(
        out / "curvature_fields.npz",
        x=x,
        t=t,
        rho=rho,
        Q=geom["Q"],
        X=geom["X"],
        eps=eps,
        R_tilde=R_tilde,
        det_g_tilde=geom["det_g_tilde"],
    )
    return payload


if __name__ == "__main__":
    result = run_curvature_case(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
