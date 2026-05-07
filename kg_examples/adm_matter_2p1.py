from __future__ import annotations

import numpy as np


def ddx(f: np.ndarray, dx: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def ddz(f: np.ndarray, dz: float) -> np.ndarray:
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out


def inverse_2metric(h_xx: np.ndarray, h_xz: np.ndarray, h_zz: np.ndarray, floor: float = 1.0e-12):
    det_h = h_xx * h_zz - h_xz * h_xz
    det_h = np.where(np.abs(det_h) < floor, np.sign(det_h + floor) * floor, det_h)
    hxx_inv = h_zz / det_h
    hxz_inv = -h_xz / det_h
    hzz_inv = h_xx / det_h
    sqrt_h = np.sqrt(np.maximum(det_h, floor))
    return hxx_inv, hxz_inv, hzz_inv, det_h, sqrt_h


def matter_observables(
    rho: np.ndarray,
    s: np.ndarray,
    lapse: np.ndarray,
    shift_x: np.ndarray,
    shift_z: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    floor: float = 1.0e-12,
):
    sx = ddx(s, dx)
    sz = ddz(s, dz)
    hxx_inv, hxz_inv, hzz_inv, det_h, sqrt_h = inverse_2metric(h_xx, h_xz, h_zz, floor=floor)
    p2 = hxx_inv * sx * sx + 2.0 * hxz_inv * sx * sz + hzz_inv * sz * sz
    energy = np.sqrt(np.maximum(mass * mass + p2, floor))
    dens = np.maximum(np.exp(beta) * sqrt_h * rho * energy, floor)
    s_t = shift_x * sx + shift_z * sz - lapse * energy
    flux_x = dens * shift_x + (lapse * dens / energy) * (hxx_inv * sx + hxz_inv * sz)
    flux_z = dens * shift_z + (lapse * dens / energy) * (hxz_inv * sx + hzz_inv * sz)
    return {
        "sx": sx,
        "sz": sz,
        "energy": energy,
        "density": dens,
        "s_t": s_t,
        "flux_x": flux_x,
        "flux_z": flux_z,
        "sqrt_h": sqrt_h,
        "det_h": det_h,
        "hxx_inv": hxx_inv,
        "hxz_inv": hxz_inv,
        "hzz_inv": hzz_inv,
    }


def recover_rho(
    density: np.ndarray,
    energy: np.ndarray,
    beta: np.ndarray,
    sqrt_h: np.ndarray,
    floor: float = 1.0e-12,
):
    return np.maximum(density / np.maximum(np.exp(beta) * sqrt_h * energy, floor), floor)


def matter_rhs(
    rho: np.ndarray,
    s: np.ndarray,
    lapse: np.ndarray,
    shift_x: np.ndarray,
    shift_z: np.ndarray,
    h_xx: np.ndarray,
    h_xz: np.ndarray,
    h_zz: np.ndarray,
    beta: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    floor: float = 1.0e-12,
):
    obs = matter_observables(
        rho=rho,
        s=s,
        lapse=lapse,
        shift_x=shift_x,
        shift_z=shift_z,
        h_xx=h_xx,
        h_xz=h_xz,
        h_zz=h_zz,
        beta=beta,
        mass=mass,
        dx=dx,
        dz=dz,
        floor=floor,
    )
    n_t = -(ddx(obs["flux_x"], dx) + ddz(obs["flux_z"], dz))
    return {
        **obs,
        "n_t": n_t,
        "rho_recovered": recover_rho(obs["density"], obs["energy"], beta, obs["sqrt_h"], floor=floor),
    }
