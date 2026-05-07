from __future__ import annotations

import numpy as np


def _mass_sq_field(mass: float | np.ndarray, shape: tuple[int, ...]) -> np.ndarray:
    arr = np.asarray(mass, dtype=float)
    if arr.ndim == 0:
        return np.full(shape, float(arr) * float(arr), dtype=float)
    if arr.shape != shape:
        raise ValueError(f"质量场形状不匹配: 期望 {shape}, 得到 {arr.shape}")
    return arr * arr


def ddx(field: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * dx)


def ddz(field: np.ndarray, dz: float) -> np.ndarray:
    return (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dz)


def inverse_metric_block(metric_cov_txz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mats = metric_cov_txz.reshape(-1, 3, 3)
    inv = np.linalg.pinv(mats, rcond=1.0e-12, hermitian=True).reshape(metric_cov_txz.shape)
    det = np.linalg.det(mats).reshape(metric_cov_txz.shape[:-2])
    return inv, det


def solve_phase_time_derivative(
    ginv_txz: np.ndarray,
    s_x: np.ndarray,
    s_z: np.ndarray,
    mass: float | np.ndarray,
    branch: str = "negative_frequency",
    s_t_reference: np.ndarray | None = None,
    disc_floor: float = 1.0e-14,
    a_floor: float = 1.0e-14,
) -> dict[str, np.ndarray]:
    a = ginv_txz[..., 0, 0]
    b = ginv_txz[..., 0, 1] * s_x + ginv_txz[..., 0, 2] * s_z
    mass_sq = _mass_sq_field(mass, s_x.shape)
    c = (
        ginv_txz[..., 1, 1] * s_x * s_x
        + 2.0 * ginv_txz[..., 1, 2] * s_x * s_z
        + ginv_txz[..., 2, 2] * s_z * s_z
        - mass_sq
    )

    disc = b * b - a * c
    disc_safe = np.maximum(disc, disc_floor)
    sqrt_disc = np.sqrt(disc_safe)

    a_safe = np.where(np.abs(a) > a_floor, a, np.sign(a) * a_floor + (a == 0.0) * a_floor)
    root_plus = (-b + sqrt_disc) / a_safe
    root_minus = (-b - sqrt_disc) / a_safe

    if s_t_reference is not None:
        choose_plus = np.abs(root_plus - s_t_reference) <= np.abs(root_minus - s_t_reference)
        s_t = np.where(choose_plus, root_plus, root_minus)
    elif branch == "negative_frequency":
        s_t = np.where(np.abs(root_minus) >= np.abs(root_plus), root_minus, root_plus)
    elif branch == "positive_frequency":
        s_t = np.where(np.abs(root_plus) >= np.abs(root_minus), root_plus, root_minus)
    else:
        raise ValueError(f"未知相位分支: {branch}")

    return {
        "s_t": s_t,
        "discriminant": disc,
        "root_plus": root_plus,
        "root_minus": root_minus,
        "a_coeff": a,
        "b_coeff": b,
        "c_coeff": c,
    }


def current_components_from_rho_s(
    rho: np.ndarray,
    s_t: np.ndarray,
    s_x: np.ndarray,
    s_z: np.ndarray,
    ginv_txz: np.ndarray,
) -> dict[str, np.ndarray]:
    j_t = rho * (
        ginv_txz[..., 0, 0] * s_t
        + ginv_txz[..., 0, 1] * s_x
        + ginv_txz[..., 0, 2] * s_z
    )
    j_x = rho * (
        ginv_txz[..., 1, 0] * s_t
        + ginv_txz[..., 1, 1] * s_x
        + ginv_txz[..., 1, 2] * s_z
    )
    j_z = rho * (
        ginv_txz[..., 2, 0] * s_t
        + ginv_txz[..., 2, 1] * s_x
        + ginv_txz[..., 2, 2] * s_z
    )
    return {
        "j_t": j_t,
        "j_x": j_x,
        "j_z": j_z,
    }


def solve_covector_time_component(
    ginv_txz: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    mass: float | np.ndarray,
    branch: str = "negative_frequency",
    u_t_reference: np.ndarray | None = None,
    disc_floor: float = 1.0e-14,
    a_floor: float = 1.0e-14,
) -> dict[str, np.ndarray]:
    a = ginv_txz[..., 0, 0]
    b = ginv_txz[..., 0, 1] * u_x + ginv_txz[..., 0, 2] * u_z
    mass_sq = _mass_sq_field(mass, u_x.shape)
    c = (
        ginv_txz[..., 1, 1] * u_x * u_x
        + 2.0 * ginv_txz[..., 1, 2] * u_x * u_z
        + ginv_txz[..., 2, 2] * u_z * u_z
        - mass_sq
    )

    disc = b * b - a * c
    disc_safe = np.maximum(disc, disc_floor)
    sqrt_disc = np.sqrt(disc_safe)
    a_safe = np.where(np.abs(a) > a_floor, a, np.sign(a) * a_floor + (a == 0.0) * a_floor)
    root_plus = (-b + sqrt_disc) / a_safe
    root_minus = (-b - sqrt_disc) / a_safe

    if u_t_reference is not None:
        choose_plus = np.abs(root_plus - u_t_reference) <= np.abs(root_minus - u_t_reference)
        u_t = np.where(choose_plus, root_plus, root_minus)
    elif branch == "negative_frequency":
        u_t = np.where(np.abs(root_minus) >= np.abs(root_plus), root_minus, root_plus)
    elif branch == "positive_frequency":
        u_t = np.where(np.abs(root_plus) >= np.abs(root_minus), root_plus, root_minus)
    else:
        raise ValueError(f"未知相位分支: {branch}")

    return {
        "u_t": u_t,
        "discriminant": disc,
        "root_plus": root_plus,
        "root_minus": root_minus,
        "a_coeff": a,
        "b_coeff": b,
        "c_coeff": c,
    }


def current_components_from_rho_u(
    rho: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    ginv_txz: np.ndarray,
) -> dict[str, np.ndarray]:
    j_t = rho * (
        ginv_txz[..., 0, 0] * u_t
        + ginv_txz[..., 0, 1] * u_x
        + ginv_txz[..., 0, 2] * u_z
    )
    j_x = rho * (
        ginv_txz[..., 1, 0] * u_t
        + ginv_txz[..., 1, 1] * u_x
        + ginv_txz[..., 1, 2] * u_z
    )
    j_z = rho * (
        ginv_txz[..., 2, 0] * u_t
        + ginv_txz[..., 2, 1] * u_x
        + ginv_txz[..., 2, 2] * u_z
    )
    return {
        "j_t": j_t,
        "j_x": j_x,
        "j_z": j_z,
    }


def conservative_density_from_rho_s(
    rho: np.ndarray,
    s: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    phase_branch: str = "negative_frequency",
    s_t_reference: np.ndarray | None = None,
    rho_floor: float = 1.0e-14,
) -> dict[str, np.ndarray]:
    ginv, det_cov = inverse_metric_block(metric_cov_txz)
    sqrt_abs_g = np.sqrt(np.abs(det_cov))

    s_x = ddx(s, dx)
    s_z = ddz(s, dz)
    phase = solve_phase_time_derivative(
        ginv,
        s_x,
        s_z,
        mass,
        branch=phase_branch,
        s_t_reference=s_t_reference,
    )
    cur = current_components_from_rho_s(np.maximum(rho, rho_floor), phase["s_t"], s_x, s_z, ginv)

    n_cons = sqrt_abs_g * cur["j_t"]
    flux_x = sqrt_abs_g * cur["j_x"]
    flux_z = sqrt_abs_g * cur["j_z"]
    return {
        "ginv": ginv,
        "det_cov": det_cov,
        "sqrt_abs_g": sqrt_abs_g,
        "s_x": s_x,
        "s_z": s_z,
        "s_t": phase["s_t"],
        "discriminant": phase["discriminant"],
        "j_t": cur["j_t"],
        "j_x": cur["j_x"],
        "j_z": cur["j_z"],
        "n_cons": n_cons,
        "flux_x": flux_x,
        "flux_z": flux_z,
    }


def recover_rho_from_conservative_density(
    n_cons: np.ndarray,
    s: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    phase_branch: str = "negative_frequency",
    s_t_reference: np.ndarray | None = None,
    jt_floor: float = 1.0e-14,
) -> dict[str, np.ndarray]:
    ginv, det_cov = inverse_metric_block(metric_cov_txz)
    sqrt_abs_g = np.sqrt(np.abs(det_cov))
    s_x = ddx(s, dx)
    s_z = ddz(s, dz)
    phase = solve_phase_time_derivative(
        ginv,
        s_x,
        s_z,
        mass,
        branch=phase_branch,
        s_t_reference=s_t_reference,
    )

    jt_unit = (
        ginv[..., 0, 0] * phase["s_t"]
        + ginv[..., 0, 1] * s_x
        + ginv[..., 0, 2] * s_z
    )
    jt_safe = np.where(np.abs(jt_unit) > jt_floor, jt_unit, np.sign(jt_unit) * jt_floor + (jt_unit == 0.0) * jt_floor)
    rho = n_cons / (sqrt_abs_g * jt_safe)
    jx = rho * (
        ginv[..., 1, 0] * phase["s_t"]
        + ginv[..., 1, 1] * s_x
        + ginv[..., 1, 2] * s_z
    )
    jz = rho * (
        ginv[..., 2, 0] * phase["s_t"]
        + ginv[..., 2, 1] * s_x
        + ginv[..., 2, 2] * s_z
    )
    return {
        "rho": rho,
        "sqrt_abs_g": sqrt_abs_g,
        "s_t": phase["s_t"],
        "s_x": s_x,
        "s_z": s_z,
        "discriminant": phase["discriminant"],
        "j_t": rho * jt_unit,
        "j_x": jx,
        "j_z": jz,
        "flux_x": sqrt_abs_g * jx,
        "flux_z": sqrt_abs_g * jz,
    }


def coordinate_matter_rhs(
    n_cons: np.ndarray,
    s: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    phase_branch: str = "negative_frequency",
    s_t_reference: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    rec = recover_rho_from_conservative_density(
        n_cons=n_cons,
        s=s,
        metric_cov_txz=metric_cov_txz,
        mass=mass,
        dx=dx,
        dz=dz,
        phase_branch=phase_branch,
        s_t_reference=s_t_reference,
    )
    n_t = -(ddx(rec["flux_x"], dx) + ddz(rec["flux_z"], dz))
    s_t = rec["s_t"]
    return {
        "n_t": n_t,
        "s_t": s_t,
        "rho": rec["rho"],
        "s_x": rec["s_x"],
        "s_z": rec["s_z"],
        "discriminant": rec["discriminant"],
        "j_t": rec["j_t"],
        "j_x": rec["j_x"],
        "j_z": rec["j_z"],
        "flux_x": rec["flux_x"],
        "flux_z": rec["flux_z"],
        "sqrt_abs_g": rec["sqrt_abs_g"],
    }


def conservative_density_from_rho_u(
    rho: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    branch: str = "negative_frequency",
    u_t_reference: np.ndarray | None = None,
    rho_floor: float = 1.0e-14,
    metric_inv_txz: np.ndarray | None = None,
    det_cov_txz: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    if metric_inv_txz is None or det_cov_txz is None:
        ginv, det_cov = inverse_metric_block(metric_cov_txz)
    else:
        ginv = metric_inv_txz
        det_cov = det_cov_txz
    sqrt_abs_g = np.sqrt(np.abs(det_cov))
    phase = solve_covector_time_component(
        ginv,
        u_x,
        u_z,
        mass,
        branch=branch,
        u_t_reference=u_t_reference,
    )
    cur = current_components_from_rho_u(np.maximum(rho, rho_floor), phase["u_t"], u_x, u_z, ginv)
    return {
        "ginv": ginv,
        "det_cov": det_cov,
        "sqrt_abs_g": sqrt_abs_g,
        "u_t": phase["u_t"],
        "u_x": u_x,
        "u_z": u_z,
        "discriminant": phase["discriminant"],
        "j_t": cur["j_t"],
        "j_x": cur["j_x"],
        "j_z": cur["j_z"],
        "n_cons": sqrt_abs_g * cur["j_t"],
        "flux_x": sqrt_abs_g * cur["j_x"],
        "flux_z": sqrt_abs_g * cur["j_z"],
    }


def recover_rho_from_conservative_density_and_u(
    n_cons: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    branch: str = "negative_frequency",
    u_t_reference: np.ndarray | None = None,
    jt_floor: float = 1.0e-14,
    metric_inv_txz: np.ndarray | None = None,
    det_cov_txz: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    if metric_inv_txz is None or det_cov_txz is None:
        ginv, det_cov = inverse_metric_block(metric_cov_txz)
    else:
        ginv = metric_inv_txz
        det_cov = det_cov_txz
    sqrt_abs_g = np.sqrt(np.abs(det_cov))
    phase = solve_covector_time_component(
        ginv,
        u_x,
        u_z,
        mass,
        branch=branch,
        u_t_reference=u_t_reference,
    )
    jt_unit = (
        ginv[..., 0, 0] * phase["u_t"]
        + ginv[..., 0, 1] * u_x
        + ginv[..., 0, 2] * u_z
    )
    jt_safe = np.where(np.abs(jt_unit) > jt_floor, jt_unit, np.sign(jt_unit) * jt_floor + (jt_unit == 0.0) * jt_floor)
    rho = n_cons / (sqrt_abs_g * jt_safe)
    jx = rho * (
        ginv[..., 1, 0] * phase["u_t"]
        + ginv[..., 1, 1] * u_x
        + ginv[..., 1, 2] * u_z
    )
    jz = rho * (
        ginv[..., 2, 0] * phase["u_t"]
        + ginv[..., 2, 1] * u_x
        + ginv[..., 2, 2] * u_z
    )
    return {
        "rho": rho,
        "sqrt_abs_g": sqrt_abs_g,
        "u_t": phase["u_t"],
        "u_x": u_x,
        "u_z": u_z,
        "discriminant": phase["discriminant"],
        "j_t": rho * jt_unit,
        "j_x": jx,
        "j_z": jz,
        "flux_x": sqrt_abs_g * jx,
        "flux_z": sqrt_abs_g * jz,
    }


def coordinate_matter_rhs_covector(
    n_cons: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    branch: str = "negative_frequency",
    u_t_reference: np.ndarray | None = None,
    metric_inv_txz: np.ndarray | None = None,
    det_cov_txz: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    rec = recover_rho_from_conservative_density_and_u(
        n_cons=n_cons,
        u_x=u_x,
        u_z=u_z,
        metric_cov_txz=metric_cov_txz,
        mass=mass,
        branch=branch,
        u_t_reference=u_t_reference,
        metric_inv_txz=metric_inv_txz,
        det_cov_txz=det_cov_txz,
    )
    n_t = -(ddx(rec["flux_x"], dx) + ddz(rec["flux_z"], dz))
    u_t = rec["u_t"]
    ux_t = ddx(u_t, dx)
    uz_t = ddz(u_t, dz)
    return {
        "n_t": n_t,
        "u_t": u_t,
        "u_x_t": ux_t,
        "u_z_t": uz_t,
        "rho": rec["rho"],
        "discriminant": rec["discriminant"],
        "j_t": rec["j_t"],
        "j_x": rec["j_x"],
        "j_z": rec["j_z"],
        "flux_x": rec["flux_x"],
        "flux_z": rec["flux_z"],
        "sqrt_abs_g": rec["sqrt_abs_g"],
    }
