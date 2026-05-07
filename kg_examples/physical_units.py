from __future__ import annotations

import math
from dataclasses import asdict, dataclass

from simulate_flat_localized_crossing_packets import FlatLocalizedCrossingParams


HBAR_C_EV_M = 1.973269804e-7
HBAR_EV_S = 6.582119569e-16
C_M_PER_S = 299_792_458.0
PLANCK_MASS_EV = 1.220890128e28
PLANCK_LENGTH_EV_INV = 1.0 / PLANCK_MASS_EV


@dataclass(frozen=True)
class OpticalPhysicalScale:
    wavelength_nm: float
    k0_ev: float
    omega0_ev: float
    mass_ev: float
    mass_over_omega: float
    sigma_parallel_ev_inv: float
    sigma_perp_ev_inv: float
    x_center_abs_ev_inv: float
    z_center_abs_ev_inv: float
    x_half_range_ev_inv: float
    z_half_range_ev_inv: float
    sigma_parallel_um: float
    sigma_perp_um: float
    x_center_abs_um: float
    z_center_abs_um: float
    x_half_range_um: float
    z_half_range_um: float
    t_meet_ev_inv: float
    t_meet_fs: float
    old_dimensionless_scale_ev_inv: float
    spectral_width_over_k0: float

    def to_json(self) -> dict[str, float]:
        return asdict(self)


def ev_inv_to_meter(value: float) -> float:
    return value * HBAR_C_EV_M


def meter_to_ev_inv(value: float) -> float:
    return value / HBAR_C_EV_M


def ev_inv_to_fs(value: float) -> float:
    return value * HBAR_EV_S * 1.0e15


def photon_k0_ev_from_wavelength_nm(wavelength_nm: float) -> float:
    wavelength_m = wavelength_nm * 1.0e-9
    k0_m_inv = 2.0 * math.pi / wavelength_m
    return k0_m_inv * HBAR_C_EV_M


def massive_omega_and_mass_from_ratio(k0_ev: float, mass_over_omega: float) -> tuple[float, float]:
    """Solve m = r omega0, omega0^2 = k0^2 + m^2."""
    r = float(mass_over_omega)
    if not (0.0 <= r < 1.0):
        raise ValueError("mass_over_omega must satisfy 0 <= r < 1")
    omega0 = float(k0_ev) / math.sqrt(1.0 - r * r)
    mass = r * omega0
    return omega0, mass


def optical_crossing_params_from_physical_scale(
    *,
    wavelength_nm: float = 1550.0,
    mass_over_omega: float = 0.1,
    resolution: int = 64,
    old_k0: float = 6.0,
    old_sigma_parallel: float = 1.6,
    old_sigma_perp: float = 1.2,
    old_center_abs: float = 6.0,
    old_half_range: float = 20.0,
    alpha: float = 0.5,
    phi0: float = 0.0,
    rho_floor: float = 1.0e-12,
) -> tuple[FlatLocalizedCrossingParams, OpticalPhysicalScale]:
    """Build the localized crossing packet in eV natural units.

    The old dimensionless setup had k0=6, sigma_parallel=1.6, sigma_perp=1.2,
    centers at |x|=|z|=6, and half range 20.  We preserve these ratios by
    scaling every length by old_k0/k0_physical.
    """
    k0_ev = photon_k0_ev_from_wavelength_nm(wavelength_nm)
    omega0_ev, mass_ev = massive_omega_and_mass_from_ratio(k0_ev, mass_over_omega)
    scale = old_k0 / k0_ev

    sigma_parallel = old_sigma_parallel * scale
    sigma_perp = old_sigma_perp * scale
    center_abs = old_center_abs * scale
    half_range = old_half_range * scale

    params = FlatLocalizedCrossingParams(
        m=mass_ev,
        alpha=alpha,
        k0=k0_ev,
        sigma_parallel=sigma_parallel,
        sigma_perp=sigma_perp,
        phi0=phi0,
        x_a0=-center_abs,
        z_a0=-center_abs,
        x_b0=center_abs,
        z_b0=-center_abs,
        x_half_range=half_range,
        z_half_range=half_range,
        nx=resolution,
        nz=resolution,
        rho_floor=rho_floor,
    )
    t_meet = params.t_meet
    scale_info = OpticalPhysicalScale(
        wavelength_nm=wavelength_nm,
        k0_ev=k0_ev,
        omega0_ev=omega0_ev,
        mass_ev=mass_ev,
        mass_over_omega=mass_over_omega,
        sigma_parallel_ev_inv=sigma_parallel,
        sigma_perp_ev_inv=sigma_perp,
        x_center_abs_ev_inv=center_abs,
        z_center_abs_ev_inv=center_abs,
        x_half_range_ev_inv=half_range,
        z_half_range_ev_inv=half_range,
        sigma_parallel_um=ev_inv_to_meter(sigma_parallel) * 1.0e6,
        sigma_perp_um=ev_inv_to_meter(sigma_perp) * 1.0e6,
        x_center_abs_um=ev_inv_to_meter(center_abs) * 1.0e6,
        z_center_abs_um=ev_inv_to_meter(center_abs) * 1.0e6,
        x_half_range_um=ev_inv_to_meter(half_range) * 1.0e6,
        z_half_range_um=ev_inv_to_meter(half_range) * 1.0e6,
        t_meet_ev_inv=t_meet,
        t_meet_fs=ev_inv_to_fs(t_meet),
        old_dimensionless_scale_ev_inv=scale,
        spectral_width_over_k0=1.0 / (sigma_parallel * k0_ev),
    )
    return params, scale_info
