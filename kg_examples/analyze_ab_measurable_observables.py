from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from calibrate_exact_kg_ps_fd_absorbing import (
    apply_boundary,
    boundary_exact,
    make_sponge,
    spectral_lap,
)
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz
from simulate_bc_psfd_imex import (
    apply_geometry_boundary,
    edge_damp,
    rhs_branch_b,
    spectral_k2,
    wave_imex_step,
)

G_NEWTON = 6.67430e-11
M_EARTH = 5.9722e24
C_LIGHT = 299792458.0
R_EARTH = 6.371e6


def unwrap_phase(psi: np.ndarray) -> np.ndarray:
    return np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)


def crop_inner(x: np.ndarray, z: np.ndarray, inner_half_range: float):
    ix = np.where(np.abs(x) <= inner_half_range + 1e-12)[0]
    iz = np.where(np.abs(z) <= inner_half_range + 1e-12)[0]
    return ix, iz, x[ix], z[iz]


def simulate_final_ab(
    alpha: float = 0.5,
    lambda_grav: float = 1.0,
    mp: float = 300.0,
    dt_target: float = 0.025,
    inner: float = 10.0,
    outer: float = 15.0,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 61,
    nkz: int = 61,
):
    solver_params = Exact2p1Params(
        alpha=alpha,
        x_half_range=outer,
        z_half_range=outer,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )
    dx = 2.0 * outer / (solver_params.nx - 1)
    dz = 2.0 * outer / (solver_params.nz - 1)
    dt = dt_target
    total_time = solver_params.overlap_time
    steps = int(round(total_time / dt))
    dt = total_time / steps
    k2 = spectral_k2(solver_params.nx, solver_params.nz, dx, dz)

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, solver_params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, solver_params)

    sponge_inner = max(inner + 2.0, 0.8 * outer)
    sigma = make_sponge(x, z, inner_half_range=sponge_inner, outer_half_range=outer, sigma_max=2.0)

    tau_b = np.zeros((solver_params.nx, solver_params.nz))
    beta_b = np.zeros_like(tau_b)
    ptau_b = np.zeros_like(tau_b)
    pbeta_b = np.zeros_like(tau_b)

    for n in range(steps):
        rho_a_outer = np.abs(psi_now) ** 2
        s_a_outer = unwrap_phase(psi_now)

        src_tau_b, src_beta_b, *_ = rhs_branch_b(
            tau_b, ptau_b, beta_b, pbeta_b, rho_a_outer, s_a_outer, solver_params.m, mp, dx, dz, k2
        )
        tau_b, ptau_b = wave_imex_step(tau_b, ptau_b, lambda_grav * src_tau_b, dt, k2)
        beta_b, pbeta_b = wave_imex_step(beta_b, pbeta_b, lambda_grav * src_beta_b, dt, k2)
        edge_damp(tau_b, sigma, dt, -2.0, 2.0)
        edge_damp(beta_b, sigma, dt, -2.0, 2.0)
        edge_damp(ptau_b, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_b, sigma, dt, -10.0, 10.0)
        apply_geometry_boundary(tau_b, beta_b)

        t_next = (n + 1) * dt
        psi_b = boundary_exact(alpha, t_next, x, z, solver_params)
        rhs = spectral_lap(psi_now, dx, dz) - solver_params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

    rho_a_outer = np.abs(psi_now) ** 2
    rho_b_outer = rho_a_outer * np.exp(-np.clip(beta_b + tau_b, -12.0, 12.0))

    ix, iz, xi, zi = crop_inner(x, z, inner)
    A = rho_a_outer[np.ix_(ix, iz)]
    B = rho_b_outer[np.ix_(ix, iz)]

    return {
        "x": xi,
        "z": zi,
        "A": A,
        "B": B,
        "tau_b": tau_b[np.ix_(ix, iz)],
        "beta_b": beta_b[np.ix_(ix, iz)],
        "dt": dt,
        "steps": steps,
        "t_over_T": 1.0,
    }


def moving_average(y: np.ndarray, win: int) -> np.ndarray:
    kernel = np.ones(win, dtype=float) / float(win)
    return np.convolve(y, kernel, mode="same")


def local_extrema(y: np.ndarray):
    dy = np.diff(y)
    sgn = np.sign(dy)
    peaks = np.where((sgn[:-1] > 0) & (sgn[1:] <= 0))[0] + 1
    mins = np.where((sgn[:-1] < 0) & (sgn[1:] >= 0))[0] + 1
    return peaks, mins


def extract_centerline_observables(
    x: np.ndarray,
    z: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    x_window: float = 4.0,
    fine_step: float = 0.001,
):
    mid = len(z) // 2
    Aline = A[:, mid]
    Bline = B[:, mid]
    mask = (x >= -x_window) & (x <= x_window)
    xw = x[mask]
    Aw = Aline[mask]
    Bw = Bline[mask]

    xf = np.arange(float(xw[0]), float(xw[-1]) + fine_step, fine_step)
    Af = np.interp(xf, xw, Aw)
    Bf = np.interp(xf, xw, Bw)

    smooth_win = 401
    Ad = (Af - moving_average(Af, smooth_win)) * np.hanning(len(xf))
    Bd = (Bf - moving_average(Bf, smooth_win)) * np.hanning(len(xf))
    F_A = np.fft.rfft(Ad)
    F_B = np.fft.rfft(Bd)
    ks = 2.0 * np.pi * np.fft.rfftfreq(len(xf), d=fine_step)
    power = np.abs(F_A) ** 2
    valid = ks > 2.0
    idx_dom = np.where(valid)[0][np.argmax(power[valid])]
    k_dom = float(ks[idx_dom])
    fringe_spacing_fft = float(2.0 * np.pi / k_dom)
    phase_shift = float(np.angle(F_B[idx_dom] * np.conj(F_A[idx_dom])))
    equivalent_shift = float(phase_shift / k_dom)
    fractional_shift = float(equivalent_shift / fringe_spacing_fft)

    peaks_a, mins_a = local_extrema(Af)
    peaks_b, mins_b = local_extrema(Bf)
    peaks_a = peaks_a[Af[peaks_a] > 0.2 * float(np.max(Af))]
    peaks_b = peaks_b[Bf[peaks_b] > 0.2 * float(np.max(Bf))]
    peaks_a_center = peaks_a[np.abs(xf[peaks_a]) < 3.0]
    peaks_b_center = peaks_b[np.abs(xf[peaks_b]) < 3.0]

    spacing_a = float(np.mean(np.diff(xf[peaks_a_center])))
    spacing_b = float(np.mean(np.diff(xf[peaks_b_center])))
    spacing_rel_diff = float((spacing_b - spacing_a) / spacing_a)

    pa = int(peaks_a[np.argmin(np.abs(xf[peaks_a]))])
    pb = int(peaks_b[np.argmin(np.abs(xf[peaks_b]))])
    la = int(mins_a[mins_a < pa][-1])
    ra = int(mins_a[mins_a > pa][0])
    lb = int(mins_b[mins_b < pb][-1])
    rb = int(mins_b[mins_b > pb][0])

    vis_a = float((Af[pa] - 0.5 * (Af[la] + Af[ra])) / (Af[pa] + 0.5 * (Af[la] + Af[ra])))
    vis_b = float((Bf[pb] - 0.5 * (Bf[lb] + Bf[rb])) / (Bf[pb] + 0.5 * (Bf[lb] + Bf[rb])))
    vis_rel_diff = float((vis_b - vis_a) / vis_a)

    central_peak_rel_intensity = float((Bf[pb] - Af[pa]) / Af[pa])
    mean_peak_rel_intensity = float(
        np.mean((Bf[peaks_b_center] - Af[peaks_a_center]) / np.maximum(Af[peaks_a_center], 1.0e-14))
    )
    support_rel_pointwise_max = float(
        np.max(np.abs((B - A)[A > 0.1 * np.max(A)] / np.maximum(A[A > 0.1 * np.max(A)], 1.0e-14)))
    )

    return {
        "x_window": x_window,
        "dominant_k_centerline": k_dom,
        "fringe_spacing_fft": fringe_spacing_fft,
        "fringe_spacing_A": spacing_a,
        "fringe_spacing_B": spacing_b,
        "fringe_spacing_rel_diff": spacing_rel_diff,
        "effective_phase_shift_rad": phase_shift,
        "equivalent_fringe_shift": equivalent_shift,
        "equivalent_fringe_shift_over_spacing": fractional_shift,
        "visibility_A": vis_a,
        "visibility_B": vis_b,
        "visibility_rel_diff": vis_rel_diff,
        "central_peak_rel_intensity_shift": central_peak_rel_intensity,
        "mean_peak_rel_intensity_shift": mean_peak_rel_intensity,
        "support_pointwise_rel_diff_max": support_rel_pointwise_max,
        "peak_positions_A": [float(v) for v in xf[peaks_a_center]],
        "peak_positions_B": [float(v) for v in xf[peaks_b_center]],
        "profile_x": xf,
        "profile_A": Af,
        "profile_B": Bf,
    }


def earth_curvature_scale() -> float:
    return G_NEWTON * M_EARTH / (C_LIGHT**2 * R_EARTH**3)


def schwarzschild_orientation_table(readout_span_m: float):
    kappa = earth_curvature_scale()
    rows = []
    for label, theta_deg in [
        ("radial", 0.0),
        ("tilted_45deg", 45.0),
        ("horizontal_tangential", 90.0),
        ("magic_angle", math.degrees(math.acos(1.0 / math.sqrt(3.0)))),
    ]:
        theta = math.radians(theta_deg)
        quad = abs(3.0 * math.cos(theta) ** 2 - 1.0)
        rows.append(
            {
                "orientation": label,
                "angle_to_radial_deg": theta_deg,
                "quadrupole_factor_abs": quad,
                "fractional_readout_modulation_estimate": kappa * (readout_span_m**2) * quad,
            }
        )
    return {
        "curvature_scale_m_minus_2": kappa,
        "readout_span_m": readout_span_m,
        "table": rows,
    }


def gw_orientation_table(h_gw: float):
    rows = []
    for label, theta_deg, plus, cross in [
        ("plus_aligned", 0.0, 0.5 * h_gw, 0.0),
        ("plus_45deg", 45.0, 0.0, 0.5 * h_gw),
        ("plus_orthogonal", 90.0, -0.5 * h_gw, 0.0),
        ("propagation_axis", None, 0.0, 0.0),
    ]:
        rows.append(
            {
                "orientation": label,
                "angle_to_plus_axis_deg": theta_deg,
                "plus_mode_fractional_modulation": plus,
                "cross_mode_fractional_modulation": cross,
                "max_abs_fractional_modulation": max(abs(plus), abs(cross)),
            }
        )
    return {
        "characteristic_strain_reference": h_gw,
        "max_best_case_fractional_modulation": 0.5 * h_gw,
        "table": rows,
    }


def apply_fractional_modulation(flat_value: float, frac: float) -> float:
    return float(flat_value * (1.0 + frac))


def build_background_estimates(obs: dict, d_phys_m: float, readout_span_m: float, h_gw: float):
    flat = {
        "fringe_spacing_fractional_shift": obs["fringe_spacing_rel_diff"],
        "effective_phase_shift_rad": obs["effective_phase_shift_rad"],
        "equivalent_fringe_shift_over_spacing": obs["equivalent_fringe_shift_over_spacing"],
        "equivalent_fringe_shift_for_input_spacing_m": obs["equivalent_fringe_shift_over_spacing"] * d_phys_m,
        "central_peak_rel_intensity_shift": obs["central_peak_rel_intensity_shift"],
    }

    schw = schwarzschild_orientation_table(readout_span_m)
    for row in schw["table"]:
        frac = row["fractional_readout_modulation_estimate"]
        row["estimated_total_AB_spacing_shift"] = apply_fractional_modulation(flat["fringe_spacing_fractional_shift"], frac)
        row["estimated_total_AB_phase_shift_rad"] = apply_fractional_modulation(flat["effective_phase_shift_rad"], frac)
        row["estimated_total_AB_equivalent_shift_m"] = apply_fractional_modulation(
            flat["equivalent_fringe_shift_for_input_spacing_m"], frac
        )
        row["estimated_total_AB_central_peak_shift"] = apply_fractional_modulation(
            flat["central_peak_rel_intensity_shift"], frac
        )

    gw = gw_orientation_table(h_gw)
    for row in gw["table"]:
        frac = row["max_abs_fractional_modulation"]
        row["estimated_total_AB_spacing_shift"] = apply_fractional_modulation(flat["fringe_spacing_fractional_shift"], frac)
        row["estimated_total_AB_phase_shift_rad"] = apply_fractional_modulation(flat["effective_phase_shift_rad"], frac)
        row["estimated_total_AB_equivalent_shift_m"] = apply_fractional_modulation(
            flat["equivalent_fringe_shift_for_input_spacing_m"], frac
        )
        row["estimated_total_AB_central_peak_shift"] = apply_fractional_modulation(
            flat["central_peak_rel_intensity_shift"], frac
        )

    return {
        "flat_reference": flat,
        "schwarzschild_earth_surface": schw,
        "nanohertz_gw_background": gw,
    }


def render_centerline(path: Path, x: np.ndarray, a: np.ndarray, b: np.ndarray):
    plt.figure(figsize=(6.8, 4.2))
    plt.plot(x, a, label="A", linewidth=1.8)
    plt.plot(x, b, label="B", linewidth=1.8)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def render_relative_diff(path: Path, x: np.ndarray, a: np.ndarray, b: np.ndarray):
    plt.figure(figsize=(6.8, 4.2))
    rel = (b - a) / np.maximum(a, 1.0e-14)
    plt.plot(x, rel, linewidth=1.8)
    plt.xlabel("x")
    plt.ylabel("(B-A)/A")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(
    outdir: str | Path,
    lambda_grav: float = 1.0,
    fringe_spacing_phys_m: float = 1.0e-6,
    readout_span_m: float = 1.0,
    h_gw: float = 2.4e-15,
    dt_target: float = 0.0125,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 81,
    nkz: int = 81,
):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    sim = simulate_final_ab(
        lambda_grav=lambda_grav,
        dt_target=dt_target,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )
    obs = extract_centerline_observables(sim["x"], sim["z"], sim["A"], sim["B"])
    bgs = build_background_estimates(obs, fringe_spacing_phys_m, readout_span_m, h_gw)

    render_centerline(out / "centerline_AB.png", obs["profile_x"], obs["profile_A"], obs["profile_B"])
    render_relative_diff(out / "centerline_AB_rel_diff.png", obs["profile_x"], obs["profile_A"], obs["profile_B"])

    summary = {
        "lambda_grav": lambda_grav,
        "measurement_assumptions": {
            "fringe_spacing_phys_m_for_absolute_shift_example": fringe_spacing_phys_m,
            "readout_span_m_for_schwarzschild_estimate": readout_span_m,
            "gw_characteristic_strain": h_gw,
            "dt_target": dt_target,
            "nx": nx,
            "nz": nz,
            "nkx": nkx,
            "nkz": nkz,
        },
        "centerline_observables": {
            k: v
            for k, v in obs.items()
            if k not in {"profile_x", "profile_A", "profile_B"}
        },
        "background_estimates": bgs,
        "notes": {
            "phase_caveat": "Current shared-baseline scaffold evolves A dynamically and reconstructs B via geometry on the same matter backbone, so the quoted phase shift is an effective fringe-pattern shift, not an independently propagated B phase variable.",
            "curved_background_caveat": "Earth Schwarzschild and GW entries are weak-field estimates on top of the current flat-background AB observable split, not a full curved-background self-consistent PDE evolution.",
            "vacuum_background_note": "Both exterior Schwarzschild and a vacuum plane GW have G_mu_nu=0, so the first-order bulk Einstein-Hilbert correction from note 024 vanishes; the listed background effects are therefore readout-scale weak-field modulations.",
        },
    }

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default=str(Path(__file__).resolve().parent.parent / "visualizations" / "ab_measurable_observables_lambda1"))
    parser.add_argument("--lambda-grav", type=float, default=1.0)
    parser.add_argument("--fringe-spacing-phys-m", type=float, default=1.0e-6)
    parser.add_argument("--readout-span-m", type=float, default=1.0)
    parser.add_argument("--h-gw", type=float, default=2.4e-15)
    parser.add_argument("--dt-target", type=float, default=0.0125)
    parser.add_argument("--nx", type=int, default=181)
    parser.add_argument("--nz", type=int, default=181)
    parser.add_argument("--nkx", type=int, default=81)
    parser.add_argument("--nkz", type=int, default=81)
    args = parser.parse_args()
    run(
        args.outdir,
        lambda_grav=args.lambda_grav,
        fringe_spacing_phys_m=args.fringe_spacing_phys_m,
        readout_span_m=args.readout_span_m,
        h_gw=args.h_gw,
        dt_target=args.dt_target,
        nx=args.nx,
        nz=args.nz,
        nkx=args.nkx,
        nkz=args.nkz,
    )
