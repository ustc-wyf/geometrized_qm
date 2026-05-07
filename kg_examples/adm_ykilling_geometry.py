from __future__ import annotations

import numpy as np


def ddx(field: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * dx)


def ddz(field: np.ndarray, dz: float) -> np.ndarray:
    return (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dz)


def ddxx(field: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(field, -1, axis=0) - 2.0 * field + np.roll(field, 1, axis=0)) / (dx * dx)


def ddzz(field: np.ndarray, dz: float) -> np.ndarray:
    return (np.roll(field, -1, axis=1) - 2.0 * field + np.roll(field, 1, axis=1)) / (dz * dz)


def ddxz(field: np.ndarray, dx: float, dz: float) -> np.ndarray:
    return (
        np.roll(np.roll(field, -1, axis=0), -1, axis=1)
        - np.roll(np.roll(field, -1, axis=0), 1, axis=1)
        - np.roll(np.roll(field, 1, axis=0), -1, axis=1)
        + np.roll(np.roll(field, 1, axis=0), 1, axis=1)
    ) / (4.0 * dx * dz)


def inverse_2metric(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    det_h = h_xx * h_zz - h_xz * h_xz
    hxx_inv = h_zz / det_h
    hxz_inv = -h_xz / det_h
    hzz_inv = h_xx / det_h
    return hxx_inv, hxz_inv, hzz_inv, det_h


def gamma_det_sqrt(det_h: np.ndarray, beta: np.ndarray) -> np.ndarray:
    return np.exp(beta) * np.sqrt(det_h)


def christoffel_2d(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    dx: float,
    dz: float,
) -> tuple[dict[tuple[int, int, int], np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    h_inv_xx, h_inv_xz, h_inv_zz, det_h = inverse_2metric(h_xx, h_xz, h_zz)
    g = {
        (0, 0): h_xx,
        (0, 1): h_xz,
        (1, 0): h_xz,
        (1, 1): h_zz,
    }
    g_inv = {
        (0, 0): h_inv_xx,
        (0, 1): h_inv_xz,
        (1, 0): h_inv_xz,
        (1, 1): h_inv_zz,
    }
    d = {
        0: {
            (0, 0): ddx(h_xx, dx),
            (0, 1): ddx(h_xz, dx),
            (1, 0): ddx(h_xz, dx),
            (1, 1): ddx(h_zz, dx),
        },
        1: {
            (0, 0): ddz(h_xx, dz),
            (0, 1): ddz(h_xz, dz),
            (1, 0): ddz(h_xz, dz),
            (1, 1): ddz(h_zz, dz),
        },
    }
    gamma: dict[tuple[int, int, int], np.ndarray] = {}
    for lam in (0, 1):
        for i in (0, 1):
            for j in (0, 1):
                acc = np.zeros_like(h_xx)
                for k in (0, 1):
                    acc = acc + 0.5 * g_inv[(lam, k)] * (
                        d[i][(k, j)] + d[j][(k, i)] - d[k][(i, j)]
                    )
                gamma[(lam, i, j)] = acc
    return gamma, (h_inv_xx, h_inv_xz, h_inv_zz, det_h)


def ricci_scalar_2d(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, dict[tuple[int, int, int], np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    gamma, inv_pack = christoffel_2d(h_xx, h_xz, h_zz, dx, dz)
    deriv = {
        0: {(lam, i, j): ddx(val, dx) for (lam, i, j), val in gamma.items()},
        1: {(lam, i, j): ddz(val, dz) for (lam, i, j), val in gamma.items()},
    }
    ricci = {(0, 0): np.zeros_like(h_xx), (0, 1): np.zeros_like(h_xx), (1, 1): np.zeros_like(h_xx)}
    for i, j in ((0, 0), (0, 1), (1, 1)):
        term = np.zeros_like(h_xx)
        for k in (0, 1):
            term = term + deriv[k][(k, i, j)] - deriv[j][(k, i, k)]
            for ell in (0, 1):
                term = term + gamma[(k, i, j)] * gamma[(ell, k, ell)] - gamma[(k, i, ell)] * gamma[(ell, j, k)]
        ricci[(i, j)] = term
    h_inv_xx, h_inv_xz, h_inv_zz, _ = inv_pack
    r2 = h_inv_xx * ricci[(0, 0)] + 2.0 * h_inv_xz * ricci[(0, 1)] + h_inv_zz * ricci[(1, 1)]
    return r2, gamma, inv_pack


def beta_hessian_and_spatial_ricci(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    r2, gamma, inv_pack = ricci_scalar_2d(h_xx, h_xz, h_zz, dx, dz)
    h_inv_xx, h_inv_xz, h_inv_zz, det_h = inv_pack

    beta_x = ddx(beta, dx)
    beta_z = ddz(beta, dz)
    beta_xx = ddxx(beta, dx)
    beta_zz = ddzz(beta, dz)
    beta_xz = ddxz(beta, dx, dz)

    hess_xx = beta_xx - gamma[(0, 0, 0)] * beta_x - gamma[(1, 0, 0)] * beta_z
    hess_xz = beta_xz - gamma[(0, 0, 1)] * beta_x - gamma[(1, 0, 1)] * beta_z
    hess_zz = beta_zz - gamma[(0, 1, 1)] * beta_x - gamma[(1, 1, 1)] * beta_z

    grad2 = h_inv_xx * beta_x * beta_x + 2.0 * h_inv_xz * beta_x * beta_z + h_inv_zz * beta_z * beta_z
    lap = h_inv_xx * hess_xx + 2.0 * h_inv_xz * hess_xz + h_inv_zz * hess_zz

    ricci2_xx = 0.5 * r2 * h_xx
    ricci2_xz = 0.5 * r2 * h_xz
    ricci2_zz = 0.5 * r2 * h_zz

    r3_xx = ricci2_xx - hess_xx - beta_x * beta_x
    r3_xz = ricci2_xz - hess_xz - beta_x * beta_z
    r3_zz = ricci2_zz - hess_zz - beta_z * beta_z
    r3_yy_reduced = -(lap + grad2)
    r3 = r2 - 2.0 * lap - 2.0 * grad2

    return {
        "det_h": det_h,
        "hxx_inv": h_inv_xx,
        "hxz_inv": h_inv_xz,
        "hzz_inv": h_inv_zz,
        "r2": r2,
        "lap_beta": lap,
        "grad_beta_sq": grad2,
        "hess_beta_xx": hess_xx,
        "hess_beta_xz": hess_xz,
        "hess_beta_zz": hess_zz,
        "r3_scalar": r3,
        "r3_xx": r3_xx,
        "r3_xz": r3_xz,
        "r3_zz": r3_zz,
        "r3_yy_reduced": r3_yy_reduced,
    }


def momentum_constraint_residual_ykilling(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    k_xx: np.ndarray,
    k_xz: np.ndarray,
    k_zz: np.ndarray,
    k_beta: np.ndarray,
    mom_x: np.ndarray,
    mom_z: np.ndarray,
    mp: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    gamma, inv_pack = christoffel_2d(h_xx, h_xz, h_zz, dx, dz)
    hxx_inv, hxz_inv, hzz_inv, _ = inv_pack

    kx_x = hxx_inv * k_xx + hxz_inv * k_xz
    kz_x = hxz_inv * k_xx + hzz_inv * k_xz
    kx_z = hxx_inv * k_xz + hxz_inv * k_zz
    kz_z = hxz_inv * k_xz + hzz_inv * k_zz

    k_trace = hxx_inv * k_xx + 2.0 * hxz_inv * k_xz + hzz_inv * k_zz + k_beta

    beta_x = ddx(beta, dx)
    beta_z = ddz(beta, dz)
    k_trace_x = ddx(k_trace, dx)
    k_trace_z = ddz(k_trace, dz)

    div_k_x = ddx(kx_x, dx) + ddz(kz_x, dz)
    div_k_x = div_k_x + (gamma[(0, 0, 0)] + gamma[(1, 0, 1)]) * kx_x
    div_k_x = div_k_x + (gamma[(0, 0, 1)] + gamma[(1, 1, 1)]) * kz_x
    div_k_x = div_k_x - gamma[(0, 0, 0)] * kx_x - gamma[(1, 0, 0)] * kx_z
    div_k_x = div_k_x - gamma[(0, 1, 0)] * kz_x - gamma[(1, 1, 0)] * kz_z
    div_k_x = div_k_x + beta_x * kx_x + beta_z * kz_x

    div_k_z = ddx(kx_z, dx) + ddz(kz_z, dz)
    div_k_z = div_k_z + (gamma[(0, 0, 0)] + gamma[(1, 0, 1)]) * kx_z
    div_k_z = div_k_z + (gamma[(0, 0, 1)] + gamma[(1, 1, 1)]) * kz_z
    div_k_z = div_k_z - gamma[(0, 0, 1)] * kx_x - gamma[(1, 0, 1)] * kx_z
    div_k_z = div_k_z - gamma[(0, 1, 1)] * kz_x - gamma[(1, 1, 1)] * kz_z
    div_k_z = div_k_z + beta_x * kx_z + beta_z * kz_z

    res_x = div_k_x - k_trace_x - mom_x / (mp * mp)
    res_z = div_k_z - k_trace_z - mom_z / (mp * mp)
    return res_x, res_z


def scalar_laplacian_ykilling(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    field: np.ndarray,
    dx: float,
    dz: float,
) -> np.ndarray:
    gamma, inv_pack = christoffel_2d(h_xx, h_xz, h_zz, dx, dz)
    h_inv_xx, h_inv_xz, h_inv_zz, _ = inv_pack

    fx = ddx(field, dx)
    fz = ddz(field, dz)
    fxx = ddxx(field, dx)
    fzz = ddzz(field, dz)
    fxz = ddxz(field, dx, dz)

    hess_xx = fxx - gamma[(0, 0, 0)] * fx - gamma[(1, 0, 0)] * fz
    hess_xz = fxz - gamma[(0, 0, 1)] * fx - gamma[(1, 0, 1)] * fz
    hess_zz = fzz - gamma[(0, 1, 1)] * fx - gamma[(1, 1, 1)] * fz

    beta_x = ddx(beta, dx)
    beta_z = ddz(beta, dz)
    grad_beta_dot_grad_f = h_inv_xx * beta_x * fx + h_inv_xz * (beta_x * fz + beta_z * fx) + h_inv_zz * beta_z * fz

    return h_inv_xx * hess_xx + 2.0 * h_inv_xz * hess_xz + h_inv_zz * hess_zz + grad_beta_dot_grad_f


def scalar_hessian_ykilling(
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    field: np.ndarray,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    gamma, inv_pack = christoffel_2d(h_xx, h_xz, h_zz, dx, dz)
    h_inv_xx, h_inv_xz, h_inv_zz, _ = inv_pack

    fx = ddx(field, dx)
    fz = ddz(field, dz)
    fxx = ddxx(field, dx)
    fzz = ddzz(field, dz)
    fxz = ddxz(field, dx, dz)

    hess_xx = fxx - gamma[(0, 0, 0)] * fx - gamma[(1, 0, 0)] * fz
    hess_xz = fxz - gamma[(0, 0, 1)] * fx - gamma[(1, 0, 1)] * fz
    hess_zz = fzz - gamma[(0, 1, 1)] * fx - gamma[(1, 1, 1)] * fz
    lap2 = h_inv_xx * hess_xx + 2.0 * h_inv_xz * hess_xz + h_inv_zz * hess_zz

    return {
        "fx": fx,
        "fz": fz,
        "hess_xx": hess_xx,
        "hess_xz": hess_xz,
        "hess_zz": hess_zz,
        "lap2": lap2,
    }
