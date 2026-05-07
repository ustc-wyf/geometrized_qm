from __future__ import annotations

import numpy as np

from adm_matter_2p1 import ddx, ddz
from simulate_flat_localized_crossing_packets import phase_field


def exact_localized_wave_derivatives(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
) -> dict[str, np.ndarray]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    kx = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(len(z), d=dz)
    kx_grid, kz_grid = np.meshgrid(kx, kz, indexing="ij")

    phase = np.exp(-1j * omega * float(t))
    psi_hat_t = psi0_hat * phase

    psi = np.fft.ifft2(psi_hat_t)
    psi_t = np.fft.ifft2((-1j * omega) * psi_hat_t)
    psi_tt = np.fft.ifft2((-(omega**2)) * psi_hat_t)
    psi_x = np.fft.ifft2((1j * kx_grid) * psi_hat_t)
    psi_z = np.fft.ifft2((1j * kz_grid) * psi_hat_t)
    psi_xx = np.fft.ifft2((-(kx_grid**2)) * psi_hat_t)
    psi_zz = np.fft.ifft2((-(kz_grid**2)) * psi_hat_t)

    return {
        "psi": psi,
        "psi_t": psi_t,
        "psi_tt": psi_tt,
        "psi_x": psi_x,
        "psi_z": psi_z,
        "psi_xx": psi_xx,
        "psi_zz": psi_zz,
    }


def bohm_fields_from_wave_derivatives(
    wave: dict[str, np.ndarray],
    rho_floor: float = 1.0e-12,
) -> dict[str, np.ndarray]:
    psi = wave["psi"]
    psi_t = wave["psi_t"]
    psi_tt = wave["psi_tt"]
    psi_x = wave["psi_x"]
    psi_z = wave["psi_z"]
    psi_xx = wave["psi_xx"]
    psi_zz = wave["psi_zz"]

    rho = np.abs(psi) ** 2
    rho_safe = np.maximum(rho, rho_floor)
    root_rho = np.sqrt(rho_safe)

    dt_rho = np.real(np.conj(psi) * psi_t + psi * np.conj(psi_t))
    dx_rho = np.real(np.conj(psi) * psi_x + psi * np.conj(psi_x))
    dz_rho = np.real(np.conj(psi) * psi_z + psi * np.conj(psi_z))

    dt2_rho = np.real(np.conj(psi) * psi_tt + psi * np.conj(psi_tt) + 2.0 * np.conj(psi_t) * psi_t)
    dx2_rho = np.real(np.conj(psi) * psi_xx + psi * np.conj(psi_xx) + 2.0 * np.conj(psi_x) * psi_x)
    dz2_rho = np.real(np.conj(psi) * psi_zz + psi * np.conj(psi_zz) + 2.0 * np.conj(psi_z) * psi_z)

    dt_root = dt_rho / (2.0 * root_rho)
    dx_root = dx_rho / (2.0 * root_rho)
    dz_root = dz_rho / (2.0 * root_rho)

    dt2_root = dt2_rho / (2.0 * root_rho) - (dt_rho**2) / (4.0 * root_rho**3)
    dx2_root = dx2_rho / (2.0 * root_rho) - (dx_rho**2) / (4.0 * root_rho**3)
    dz2_root = dz2_rho / (2.0 * root_rho) - (dz_rho**2) / (4.0 * root_rho**3)

    q_field = (dt2_root - dx2_root - dz2_root) / root_rho

    s_t = np.imag(np.conj(psi) * psi_t) / rho_safe
    s_x = np.imag(np.conj(psi) * psi_x) / rho_safe
    s_z = np.imag(np.conj(psi) * psi_z) / rho_safe

    x_field = s_t * s_t - s_x * s_x - s_z * s_z
    y_field = s_t * dt_root - s_x * dx_root - s_z * dz_root
    z_field = dt_root * dt_root - dx_root * dx_root - dz_root * dz_root
    delta_field = x_field * z_field - y_field * y_field

    s = phase_field(psi)

    return {
        "rho": rho,
        "rho_safe": rho_safe,
        "root_rho": root_rho,
        "s": s,
        "s_t": s_t,
        "s_x": s_x,
        "s_z": s_z,
        "r_t": dt_root,
        "r_x": dx_root,
        "r_z": dz_root,
        "Q": np.real_if_close(q_field),
        "X": np.real_if_close(x_field),
        "Y": np.real_if_close(y_field),
        "Z": np.real_if_close(z_field),
        "Delta": np.real_if_close(delta_field),
    }


def projector_from_u_r_pinv(
    s_t: np.ndarray,
    s_x: np.ndarray,
    s_z: np.ndarray,
    r_t: np.ndarray,
    r_x: np.ndarray,
    r_z: np.ndarray,
    x_field: np.ndarray,
    y_field: np.ndarray,
    z_field: np.ndarray,
    pinv_rcond: float = 1.0e-10,
) -> np.ndarray:
    u_up = np.stack([s_t, -s_x, -s_z], axis=-1)
    r_up = np.stack([r_t, -r_x, -r_z], axis=-1)
    basis = np.stack([u_up, r_up], axis=-1)  # (..., 3, 2)

    gram = np.empty(x_field.shape + (2, 2), dtype=float)
    gram[..., 0, 0] = x_field
    gram[..., 0, 1] = y_field
    gram[..., 1, 0] = y_field
    gram[..., 1, 1] = z_field

    gram_pinv = np.linalg.pinv(gram, rcond=pinv_rcond, hermitian=True)
    projector = np.einsum("...ma,...ab,...nb->...mn", basis, gram_pinv, basis, optimize=True)
    return np.real_if_close(projector)


def tilde_inverse_from_flat_mixed_transform(
    bohm: dict[str, np.ndarray],
    x_floor: float = 1.0e-10,
    pinv_rcond: float = 1.0e-10,
) -> dict[str, np.ndarray]:
    x_field = np.where(np.abs(bohm["X"]) > x_floor, bohm["X"], np.sign(bohm["X"]) * x_floor + (bohm["X"] == 0.0) * x_floor)
    projector = projector_from_u_r_pinv(
        bohm["s_t"],
        bohm["s_x"],
        bohm["s_z"],
        bohm["r_t"],
        bohm["r_x"],
        bohm["r_z"],
        bohm["X"],
        bohm["Y"],
        bohm["Z"],
        pinv_rcond=pinv_rcond,
    )

    eta = np.zeros(bohm["X"].shape + (3, 3), dtype=float)
    eta[..., 0, 0] = 1.0
    eta[..., 1, 1] = -1.0
    eta[..., 2, 2] = -1.0

    coeff = (bohm["Q"] / x_field)[..., None, None]
    ginv = eta - coeff * projector
    return {
        "tilde_ginv_txz": np.real_if_close(ginv),
        "projector": np.real_if_close(projector),
        "coeff_Q_over_X": np.real_if_close(bohm["Q"] / x_field),
    }


def safe_covariant_from_inverse(
    tilde_ginv_txz: np.ndarray,
    rho: np.ndarray,
    det_floor: float = 1.0e-12,
    tail_rel_cut: float = 1.0e-8,
    pinv_rcond: float = 1.0e-12,
) -> dict[str, np.ndarray]:
    mats = tilde_ginv_txz.reshape(-1, 3, 3)
    dets = np.linalg.det(mats)
    rho_flat = rho.reshape(-1)
    rho_cut = float(np.max(rho)) * tail_rel_cut
    tail_mask = rho_flat <= rho_cut
    safe_mask = np.abs(dets) > det_floor

    cov = np.linalg.pinv(mats, rcond=pinv_rcond, hermitian=True)

    singular_mask = (~safe_mask).reshape(tilde_ginv_txz.shape[:2])
    tail_mask_2d = tail_mask.reshape(tilde_ginv_txz.shape[:2])
    return {
        "cov_txz": np.real_if_close(cov.reshape(tilde_ginv_txz.shape)),
        "det_txz": np.real_if_close(dets.reshape(tilde_ginv_txz.shape[:2])),
        "singular_mask": singular_mask,
        "tail_mask": tail_mask_2d,
    }


def adm_geometry_from_tilde_inverse(
    tilde_ginv_txz: np.ndarray,
    ginv_floor: float = 1.0e-12,
) -> dict[str, np.ndarray]:
    gtt = tilde_ginv_txz[..., 0, 0]
    if np.any(gtt <= ginv_floor):
        raise ValueError("初始 tilde g^{tt} 非正，无法提取 lapse。")

    lapse = 1.0 / np.sqrt(np.maximum(gtt, ginv_floor))
    shift_x = -tilde_ginv_txz[..., 0, 1] / np.maximum(gtt, ginv_floor)
    shift_z = -tilde_ginv_txz[..., 0, 2] / np.maximum(gtt, ginv_floor)

    cov = np.linalg.inv(tilde_ginv_txz.reshape(-1, 3, 3)).reshape(tilde_ginv_txz.shape)
    h_xx = -cov[..., 1, 1]
    h_xz = -cov[..., 1, 2]
    h_zz = -cov[..., 2, 2]
    beta_geom = np.zeros_like(h_xx)
    return {
        "lapse": np.real_if_close(lapse),
        "shift_x": np.real_if_close(shift_x),
        "shift_z": np.real_if_close(shift_z),
        "h_xx": np.real_if_close(h_xx),
        "h_xz": np.real_if_close(h_xz),
        "h_zz": np.real_if_close(h_zz),
        "beta": beta_geom,
        "cov_txz": np.real_if_close(cov),
    }


def lie_metric_2d(
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

    adv_xx = shift_x * ddx(h_xx, dx) + shift_z * ddz(h_xx, dz)
    adv_xz = shift_x * ddx(h_xz, dx) + shift_z * ddz(h_xz, dz)
    adv_zz = shift_x * ddx(h_zz, dx) + shift_z * ddz(h_zz, dz)

    lie_xx = adv_xx + 2.0 * (h_xx * sx_x + h_xz * sz_x)
    lie_xz = adv_xz + h_xx * sx_z + h_xz * sz_z + h_xz * sx_x + h_zz * sz_x
    lie_zz = adv_zz + 2.0 * (h_xz * sx_z + h_zz * sz_z)
    return lie_xx, lie_xz, lie_zz


def localized_direct_tilde_adm_initial(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    probe_dt: float = 1.0e-3,
    rho_floor: float = 1.0e-12,
    x_floor: float = 1.0e-10,
    pinv_rcond: float = 1.0e-10,
) -> dict[str, object]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    wave_0 = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x, z, 0.0)
    bohm_0 = bohm_fields_from_wave_derivatives(wave_0, rho_floor=rho_floor)
    inv_0 = tilde_inverse_from_flat_mixed_transform(bohm_0, x_floor=x_floor, pinv_rcond=pinv_rcond)
    adm_0 = adm_geometry_from_tilde_inverse(inv_0["tilde_ginv_txz"])

    wave_p = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x, z, probe_dt)
    bohm_p = bohm_fields_from_wave_derivatives(wave_p, rho_floor=rho_floor)
    inv_p = tilde_inverse_from_flat_mixed_transform(bohm_p, x_floor=x_floor, pinv_rcond=pinv_rcond)
    adm_p = adm_geometry_from_tilde_inverse(inv_p["tilde_ginv_txz"])

    wave_m = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x, z, -probe_dt)
    bohm_m = bohm_fields_from_wave_derivatives(wave_m, rho_floor=rho_floor)
    inv_m = tilde_inverse_from_flat_mixed_transform(bohm_m, x_floor=x_floor, pinv_rcond=pinv_rcond)
    adm_m = adm_geometry_from_tilde_inverse(inv_m["tilde_ginv_txz"])

    h_xx_t = (adm_p["h_xx"] - adm_m["h_xx"]) / (2.0 * probe_dt)
    h_xz_t = (adm_p["h_xz"] - adm_m["h_xz"]) / (2.0 * probe_dt)
    h_zz_t = (adm_p["h_zz"] - adm_m["h_zz"]) / (2.0 * probe_dt)

    lie_xx, lie_xz, lie_zz = lie_metric_2d(
        adm_0["h_xx"], adm_0["h_xz"], adm_0["h_zz"],
        adm_0["shift_x"], adm_0["shift_z"], dx, dz
    )

    lapse0 = np.maximum(adm_0["lapse"], 1.0e-12)
    k_xx = -(h_xx_t - lie_xx) / (2.0 * lapse0)
    k_xz = -(h_xz_t - lie_xz) / (2.0 * lapse0)
    k_zz = -(h_zz_t - lie_zz) / (2.0 * lapse0)
    k_beta = np.zeros_like(k_xx)

    adm_0.update(
        {
            "k_xx": np.real_if_close(k_xx),
            "k_xz": np.real_if_close(k_xz),
            "k_zz": np.real_if_close(k_zz),
            "k_beta": k_beta,
        }
    )

    diagnostics = {
        "rho0": bohm_0["rho"],
        "s0": bohm_0["s"],
        "Q0": bohm_0["Q"],
        "X0": bohm_0["X"],
        "Y0": bohm_0["Y"],
        "Z0": bohm_0["Z"],
        "Delta0": bohm_0["Delta"],
        "projector0": inv_0["projector"],
        "coeff_Q_over_X_0": inv_0["coeff_Q_over_X"],
    }
    return {
        "adm": adm_0,
        "diagnostics": diagnostics,
    }


def localized_direct_tilde_coordinate_snapshot(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    rho_floor: float = 1.0e-12,
    x_floor: float = 1.0e-10,
    pinv_rcond: float = 1.0e-10,
) -> dict[str, object]:
    wave_t = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x, z, t)
    bohm_t = bohm_fields_from_wave_derivatives(wave_t, rho_floor=rho_floor)
    inv_t = tilde_inverse_from_flat_mixed_transform(bohm_t, x_floor=x_floor, pinv_rcond=pinv_rcond)
    safe_cov = safe_covariant_from_inverse(
        inv_t["tilde_ginv_txz"],
        bohm_t["rho"],
        det_floor=1.0e-12,
        tail_rel_cut=1.0e-8,
        pinv_rcond=pinv_rcond,
    )
    return {
        "wave": wave_t,
        "bohm": bohm_t,
        "inverse": inv_t,
        "cov_txz": safe_cov["cov_txz"],
        "det_txz": safe_cov["det_txz"],
        "singular_mask": safe_cov["singular_mask"],
        "tail_mask": safe_cov["tail_mask"],
    }


def localized_direct_tilde_coordinate_initial(
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    initial_time: float = 0.0,
    probe_dt: float = 1.0e-3,
    rho_floor: float = 1.0e-12,
    x_floor: float = 1.0e-10,
    pinv_rcond: float = 1.0e-10,
) -> dict[str, object]:
    snap_0 = localized_direct_tilde_coordinate_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=initial_time,
        rho_floor=rho_floor,
        x_floor=x_floor,
        pinv_rcond=pinv_rcond,
    )
    snap_p = localized_direct_tilde_coordinate_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=initial_time + probe_dt,
        rho_floor=rho_floor,
        x_floor=x_floor,
        pinv_rcond=pinv_rcond,
    )
    snap_m = localized_direct_tilde_coordinate_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=initial_time - probe_dt,
        rho_floor=rho_floor,
        x_floor=x_floor,
        pinv_rcond=pinv_rcond,
    )

    cov_0 = np.real_if_close(snap_0["cov_txz"])
    cov_t = np.real_if_close((snap_p["cov_txz"] - snap_m["cov_txz"]) / (2.0 * probe_dt))
    ginv_0 = np.real_if_close(snap_0["inverse"]["tilde_ginv_txz"])
    ginv_t = np.real_if_close(
        (snap_p["inverse"]["tilde_ginv_txz"] - snap_m["inverse"]["tilde_ginv_txz"]) / (2.0 * probe_dt)
    )

    bohm_0 = snap_0["bohm"]
    rho0 = bohm_0["rho"]
    s0 = bohm_0["s"]
    support = rho0 > 1.0e-3 * float(np.max(rho0))

    eigvals = np.linalg.eigvalsh(ginv_0.reshape(-1, 3, 3)).reshape(ginv_0.shape[:-2] + (3,))
    pos_count = (eigvals > 0.0).sum(axis=-1)
    neg_count = (eigvals < 0.0).sum(axis=-1)
    lorentz_mask = (pos_count == 1) & (neg_count == 2)

    diagnostics = {
        "rho0": rho0,
        "s0": s0,
        "s_t0": bohm_0["s_t"],
        "s_x0": bohm_0["s_x"],
        "s_z0": bohm_0["s_z"],
        "r_t0": bohm_0["r_t"],
        "r_x0": bohm_0["r_x"],
        "r_z0": bohm_0["r_z"],
        "Q0": bohm_0["Q"],
        "X0": bohm_0["X"],
        "Y0": bohm_0["Y"],
        "Z0": bohm_0["Z"],
        "Delta0": bohm_0["Delta"],
        "coeff_Q_over_X_0": snap_0["inverse"]["coeff_Q_over_X"],
        "projector0": snap_0["inverse"]["projector"],
        "support_mask": support,
        "lorentz_mask": lorentz_mask,
        "gtt0": ginv_0[..., 0, 0],
        "eigvals0": eigvals,
        "det0": snap_0["det_txz"],
        "singular_mask0": snap_0["singular_mask"],
        "tail_mask0": snap_0["tail_mask"],
    }
    return {
        "metric_cov_txz_0": cov_0,
        "metric_cov_txz_t0": cov_t,
        "metric_inv_txz_0": ginv_0,
        "metric_inv_txz_t0": ginv_t,
        "rho0": rho0,
        "s0": s0,
        "initial_time": float(initial_time),
        "diagnostics": diagnostics,
    }
