from __future__ import annotations

import numpy as np


def spectral_dx(f: np.ndarray, dx: float) -> np.ndarray:
    nx, nz = f.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    return np.fft.ifft2((1j * kx)[:, None] * np.fft.fft2(f)).real


def spectral_dz(f: np.ndarray, dz: float) -> np.ndarray:
    nx, nz = f.shape
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    return np.fft.ifft2((1j * kz)[None, :] * np.fft.fft2(f)).real


def project_to_gradient(
    u_x: np.ndarray,
    u_z: np.ndarray,
    dx: float,
    dz: float,
    *,
    mean_mode: str = "zero",
) -> dict[str, np.ndarray | float]:
    """Helmholtz-style projection of a 2D vector field to a potential gradient.

    The returned phase potential S satisfies
        grad S = Proj_grad(u)
    in a periodic spectral sense.  This is a purely numerical compatibility
    step for phase-clock evolution; it does not alter the continuum PDE.
    """
    if not np.allclose(dx, dz, rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("phase_clock_projection currently assumes dx = dz.")

    ux = np.asarray(u_x, dtype=float)
    uz = np.asarray(u_z, dtype=float)
    nx, nz = ux.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    k2 = KX * KX + KZ * KZ

    ux_hat = np.fft.fft2(ux)
    uz_hat = np.fft.fft2(uz)
    div_hat = 1j * KX * ux_hat + 1j * KZ * uz_hat

    psi_hat = np.zeros_like(div_hat)
    mask = k2 > 1.0e-300
    psi_hat[mask] = -div_hat[mask] / k2[mask]
    if mean_mode == "zero":
        psi_hat[0, 0] = 0.0
    elif mean_mode != "keep":
        raise ValueError("mean_mode must be 'zero' or 'keep'.")

    s = np.fft.ifft2(psi_hat).real
    s_x = spectral_dx(s, dx)
    s_z = spectral_dz(s, dz)
    resid_x = ux - s_x
    resid_z = uz - s_z
    resid_norm = np.sqrt(resid_x * resid_x + resid_z * resid_z)
    input_norm = np.sqrt(ux * ux + uz * uz)
    denom = np.maximum(input_norm, 1.0e-300)
    rel_resid = resid_norm / denom
    return {
        "s": s,
        "s_x": s_x,
        "s_z": s_z,
        "resid_x": resid_x,
        "resid_z": resid_z,
        "resid_norm": resid_norm,
        "rel_resid": rel_resid,
        "divergence": np.real(np.fft.ifft2(div_hat)),
    }


def project_to_weighted_gradient(
    u_x: np.ndarray,
    u_z: np.ndarray,
    dx: float,
    dz: float,
    *,
    weight: np.ndarray | None = None,
    weight_floor: float = 1.0e-4,
    mean_mode: str = "weighted",
) -> dict[str, np.ndarray | float]:
    """Weighted local phase reconstruction on a finite rectangular grid.

    The solve minimizes a discrete weighted least-squares mismatch between
    ``(u_x, u_z)`` and ``grad S`` with homogeneous Dirichlet data on the outer
    boundary.  This is better suited to localized wave packets than a global
    periodic Helmholtz projection, and it avoids wrap-around coupling between
    opposite edges of the box.
    """
    ux = np.asarray(u_x, dtype=float)
    uz = np.asarray(u_z, dtype=float)
    if ux.shape != uz.shape:
        raise ValueError("u_x and u_z must have the same shape.")
    if ux.ndim != 2:
        raise ValueError("project_to_weighted_gradient expects a 2D grid.")
    nx, nz = ux.shape
    if nx < 3 or nz < 3:
        raise ValueError("grid must have at least 3x3 points.")

    if weight is None:
        w = np.ones_like(ux, dtype=float)
    else:
        w = np.asarray(weight, dtype=float)
        if w.shape != ux.shape:
            raise ValueError("weight must have the same shape as u_x/u_z.")
        w = np.maximum(w, 0.0)
    w = w + float(weight_floor)

    face_x = 0.5 * (w[1:, :] + w[:-1, :])
    face_z = 0.5 * (w[:, 1:] + w[:, :-1])

    flux_x = 0.5 * (w[1:, :] * ux[1:, :] + w[:-1, :] * ux[:-1, :])
    flux_z = 0.5 * (w[:, 1:] * uz[:, 1:] + w[:, :-1] * uz[:, :-1])
    div_wu = np.zeros_like(ux)
    div_wu[1:-1, 1:-1] = (flux_x[1:, 1:-1] - flux_x[:-1, 1:-1]) / dx + (flux_z[1:-1, 1:] - flux_z[1:-1, :-1]) / dz

    rhs = -div_wu[1:-1, 1:-1].copy()
    diag = (
        (face_x[1:, 1:-1] + face_x[:-1, 1:-1]) / (dx * dx)
        + (face_z[1:-1, 1:] + face_z[1:-1, :-1]) / (dz * dz)
    )
    diag = np.maximum(diag, weight_floor)

    interior_shape = rhs.shape

    def unpack(vec: np.ndarray) -> np.ndarray:
        out = np.zeros_like(ux)
        out[1:-1, 1:-1] = np.asarray(vec, dtype=float).reshape(interior_shape)
        return out

    def pack(field: np.ndarray) -> np.ndarray:
        return np.asarray(field[1:-1, 1:-1], dtype=float).reshape(-1)

    def apply_operator(vec: np.ndarray) -> np.ndarray:
        s_full = unpack(vec)
        grad_x = (s_full[1:, :] - s_full[:-1, :]) / dx
        grad_z = (s_full[:, 1:] - s_full[:, :-1]) / dz
        flux_x_s = face_x * grad_x
        flux_z_s = face_z * grad_z
        div = np.zeros_like(s_full)
        div[1:-1, 1:-1] = (flux_x_s[1:, 1:-1] - flux_x_s[:-1, 1:-1]) / dx + (flux_z_s[1:-1, 1:] - flux_z_s[1:-1, :-1]) / dz
        return pack(-div)

    b = rhs.reshape(-1)
    x = np.zeros_like(b)
    r = b - apply_operator(x)
    z = r / diag.reshape(-1)
    p = z.copy()
    rz_old = float(np.dot(r, z))
    b_norm = max(float(np.linalg.norm(b)), 1.0e-300)
    if not np.isfinite(b_norm):
        b_norm = 1.0
    max_iter = min(max(4 * x.size, 64), 2000)
    tol = 1.0e-8
    for _ in range(max_iter):
        Ap = apply_operator(p)
        denom = float(np.dot(p, Ap))
        if abs(denom) <= 1.0e-300:
            break
        alpha = rz_old / denom
        x = x + alpha * p
        r = r - alpha * Ap
        if float(np.linalg.norm(r)) <= tol * b_norm:
            break
        z = r / diag.reshape(-1)
        rz_new = float(np.dot(r, z))
        if abs(rz_old) <= 1.0e-300:
            break
        beta = rz_new / rz_old
        p = z + beta * p
        rz_old = rz_new

    s = unpack(x)
    if mean_mode not in {"weighted", "zero"}:
        raise ValueError("mean_mode must be 'weighted' or 'zero'.")
    if mean_mode == "weighted":
        finite = np.isfinite(s)
        if np.any(finite):
            wf = w[finite]
            sf = s[finite]
            denom_w = max(float(np.sum(wf)), 1.0e-300)
            s -= float(np.sum(wf * sf) / denom_w)

    s_x = np.gradient(s, dx, axis=0, edge_order=2)
    s_z = np.gradient(s, dz, axis=1, edge_order=2)
    resid_x = ux - s_x
    resid_z = uz - s_z
    resid_norm = np.sqrt(resid_x * resid_x + resid_z * resid_z)
    input_norm = np.sqrt(ux * ux + uz * uz)
    denom = np.maximum(input_norm, 1.0e-300)
    rel_resid = resid_norm / denom
    weighted_resid = float(np.sum(w * resid_norm) / max(np.sum(w * input_norm), 1.0e-300))
    return {
        "s": s,
        "s_x": s_x,
        "s_z": s_z,
        "resid_x": resid_x,
        "resid_z": resid_z,
        "resid_norm": resid_norm,
        "rel_resid": rel_resid,
        "weighted_rel_resid": weighted_resid,
        "weight": w,
        "divergence": div_wu,
    }
