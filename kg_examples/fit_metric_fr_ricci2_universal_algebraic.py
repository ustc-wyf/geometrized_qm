from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial import chebyshev as cheb

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from diagnose_metric_fr_direction_matching import stats_abs, tensor_norm
from diagnose_metric_fr_ricci2_direction_matching import invariant_contract2
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference, make_cropped_snapshot, real_array
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det


@dataclass
class CaseData:
    tau_old: float
    x: np.ndarray
    z: np.ndarray
    rho: np.ndarray
    metric_cov: np.ndarray
    ricci: np.ndarray
    ricci_square_tensor: np.ndarray
    r_scalar: np.ndarray
    p_scalar: np.ndarray
    target: np.ndarray
    rhs_norm: np.ndarray
    masks: dict[str, np.ndarray]


def make_snapshot_pack(args: argparse.Namespace, time_ev_inv: float, ix: np.ndarray, iz: np.ndarray, ref: dict[str, object]) -> dict[str, object]:
    return make_cropped_snapshot(
        psi0=np.asarray(ref["psi0"]),
        psi0_hat=np.asarray(ref["psi0_hat"]),
        omega=np.asarray(ref["omega"]),
        x_full=np.asarray(ref["x_full"]),
        z_full=np.asarray(ref["z_full"]),
        ix=ix,
        iz=iz,
        t=time_ev_inv,
        rho_floor=float(args.rho_floor),
        x_floor=float(args.x_floor),
        pinv_rcond=float(args.pinv_rcond),
        with_metric=True,
    )


def build_case(args: argparse.Namespace, ref: dict[str, object], tau_old: float) -> CaseData:
    params = ref["params"]
    scale = ref["scale"]
    x_full = np.asarray(ref["x_full"], dtype=float)
    z_full = np.asarray(ref["z_full"], dtype=float)
    old_scale = float(scale.old_dimensionless_scale_ev_inv)
    t_meet_old = float(params.t_meet / old_scale)
    old_time = t_meet_old + float(tau_old)
    time_ev_inv = old_time * old_scale
    dt = float(args.dt_old) * old_scale

    window_ev_inv = float(args.window_um) * 1.0e-6 / HBAR_C_EV_M
    ix = np.where(np.abs(x_full) <= window_ev_inv)[0]
    iz = np.where(np.abs(z_full) <= window_ev_inv)[0]
    x = x_full[ix]
    z = z_full[iz]
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    snaps = {
        -1: make_snapshot_pack(args, time_ev_inv - dt, ix, iz, ref),
        0: make_snapshot_pack(args, time_ev_inv, ix, iz, ref),
        1: make_snapshot_pack(args, time_ev_inv + dt, ix, iz, ref),
    }

    metric_m = real_array(snaps[-1]["cov_txz"])
    metric_0 = real_array(snaps[0]["cov_txz"])
    metric_p = real_array(snaps[1]["cov_txz"])
    dg_0, d2g_0 = metric_jets_full(metric_m, metric_0, metric_p, dt, dx, dz)
    metric_inv, _, ricci, r_scalar = geometry_data(metric_0, dg_0, d2g_0)

    ricci_up = np.einsum("...ac,...cb->...ab", metric_inv, ricci, optimize=True)
    ricci_square_tensor = np.einsum("...ac,...cb->...ab", ricci, ricci_up, optimize=True)
    p_scalar = invariant_contract2(metric_inv, ricci, ricci)

    bohm = snaps[0]["bohm"]
    rho = real_array(bohm["rho"])
    x_field = real_array(bohm["X"])
    measure = np.abs(x_field) * rho / (float(params.m) * float(params.m))
    sqrt_abs_g = safe_sqrt_abs_det(metric_0)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)
    stress, _ = stress_tensor_tilde(
        metric_0,
        metric_inv,
        rho_tilde,
        real_array(bohm["s_t"]),
        real_array(bohm["s_x"]),
        real_array(bohm["s_z"]),
        m=float(params.m),
    )
    target = real_array(stress) / (float(args.mp) * float(args.mp))
    rhs_norm = tensor_norm(target)

    support = (rho > float(args.support_rho_frac) * float(np.max(rho))) & (
        measure > float(args.support_measure_frac) * float(np.max(measure))
    )
    trusted = erode_mask_8(support, int(args.trusted_erosion))
    if not np.any(trusted):
        trusted = support.copy()
    masks = {
        "support": support,
        "trusted": trusted,
        "core1": support & (rho > 1.0e-2 * float(np.max(rho))),
        "core10": support & (rho > 1.0e-1 * float(np.max(rho))),
    }
    return CaseData(
        tau_old=float(tau_old),
        x=x,
        z=z,
        rho=rho,
        metric_cov=metric_0,
        ricci=ricci,
        ricci_square_tensor=ricci_square_tensor,
        r_scalar=r_scalar,
        p_scalar=p_scalar,
        target=target,
        rhs_norm=rhs_norm,
        masks=masks,
    )


def cheb_derivative_vandermonde(t: np.ndarray, degree: int) -> np.ndarray:
    out = np.empty((t.size, degree + 1), dtype=float)
    for n in range(degree + 1):
        out[:, n] = cheb.Chebyshev.basis(n).deriv(1)(t)
    return out


def transform_values(values: np.ndarray, scale0: float, s_min: float, s_max: float) -> tuple[np.ndarray, np.ndarray]:
    vals = np.asarray(values, dtype=float)
    s = np.arcsinh(vals / scale0)
    t_scale = 2.0 / max(s_max - s_min, 1.0e-300)
    t = np.clip(t_scale * (s - s_min) - 1.0, -1.0, 1.0)
    ds_dv = 1.0 / np.sqrt(vals * vals + scale0 * scale0)
    dt_dv = t_scale * ds_dv
    return t, dt_dv


def basis_values_2d(
    r: np.ndarray,
    p: np.ndarray,
    degree_r: int,
    degree_p: int,
    r0: float,
    p0: float,
    r_s_min: float,
    r_s_max: float,
    p_s_min: float,
    p_s_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.asarray(r, dtype=float)
    p = np.asarray(p, dtype=float)
    tr, dtr_dr = transform_values(r, r0, r_s_min, r_s_max)
    tp, dtp_dp = transform_values(p, p0, p_s_min, p_s_max)

    vr = cheb.chebvander(tr, degree_r)
    vp = cheb.chebvander(tp, degree_p)
    dvr = cheb_derivative_vandermonde(tr, degree_r) * dtr_dr[:, None]
    dvp = cheb_derivative_vandermonde(tp, degree_p) * dtp_dp[:, None]

    f = np.einsum("in,im->inm", vr, vp, optimize=True).reshape(r.size, -1)
    f_r = np.einsum("in,im->inm", dvr, vp, optimize=True).reshape(r.size, -1)
    f_p = np.einsum("in,im->inm", vr, dvp, optimize=True).reshape(r.size, -1)
    return f, f_r, f_p


def make_fit_scales(cases: list[CaseData], fit_region: str) -> dict[str, float]:
    r_fit = np.concatenate([case.r_scalar[case.masks[fit_region]] for case in cases])
    p_fit = np.concatenate([case.p_scalar[case.masks[fit_region]] for case in cases])
    r_abs = np.abs(r_fit[np.isfinite(r_fit)])
    p_abs = np.abs(p_fit[np.isfinite(p_fit)])
    r0 = max(float(np.percentile(r_abs, 50.0)), 1.0e-30)
    p0 = max(float(np.percentile(p_abs, 50.0)), 1.0e-30)
    sr = np.arcsinh(r_fit / r0)
    sp = np.arcsinh(p_fit / p0)
    return {
        "r0": r0,
        "p0": p0,
        "r_s_min": float(np.percentile(sr[np.isfinite(sr)], 0.5)),
        "r_s_max": float(np.percentile(sr[np.isfinite(sr)], 99.5)),
        "p_s_min": float(np.percentile(sp[np.isfinite(sp)], 0.5)),
        "p_s_max": float(np.percentile(sp[np.isfinite(sp)], 99.5)),
    }


def design_rows_for_case(
    case: CaseData,
    mask: np.ndarray,
    degree_r: int,
    degree_p: int,
    scales: dict[str, float],
    target_floor_frac: float,
) -> tuple[np.ndarray, np.ndarray]:
    idx = np.where(mask.reshape(-1))[0]
    r = case.r_scalar.reshape(-1)[idx]
    p = case.p_scalar.reshape(-1)[idx]
    rho = case.rho.reshape(-1)[idx]
    rhs_norm = case.rhs_norm.reshape(-1)[idx]
    target = case.target.reshape(-1, 9)[idx]
    metric = case.metric_cov.reshape(-1, 9)[idx]
    ricci = case.ricci.reshape(-1, 9)[idx]
    q_tensor = case.ricci_square_tensor.reshape(-1, 9)[idx]

    f, f_r, f_p = basis_values_2d(
        r,
        p,
        degree_r,
        degree_p,
        scales["r0"],
        scales["p0"],
        scales["r_s_min"],
        scales["r_s_max"],
        scales["p_s_min"],
        scales["p_s_max"],
    )
    n = idx.size
    k = f.shape[1]
    rows = np.zeros((n * 9, k), dtype=float)
    target_rows = target.reshape(n * 9)
    for comp in range(9):
        block = f_r * ricci[:, comp, None] + 2.0 * f_p * q_tensor[:, comp, None] - 0.5 * f * metric[:, comp, None]
        rows[comp::9, :] = block

    rho_w = np.sqrt(np.maximum(rho, 0.0) / max(float(np.max(case.rho)), 1.0e-300))
    floor = float(target_floor_frac) * max(float(np.percentile(rhs_norm[np.isfinite(rhs_norm)], 95.0)), 1.0e-300)
    rel_w = 1.0 / np.maximum(rhs_norm, floor)
    w = rho_w * rel_w
    w_rows = np.repeat(w, 9)
    return rows * w_rows[:, None], target_rows * w_rows


def fit_universal(
    cases: list[CaseData],
    fit_region: str,
    degree_r: int,
    degree_p: int,
    target_floor_frac: float,
) -> dict[str, object]:
    scales = make_fit_scales(cases, fit_region)
    matrices = [
        design_rows_for_case(case, case.masks[fit_region], degree_r, degree_p, scales, target_floor_frac)
        for case in cases
    ]
    a = np.vstack([m[0] for m in matrices])
    b = np.concatenate([m[1] for m in matrices])
    column_scale = np.linalg.norm(a, axis=0)
    column_scale = np.where(column_scale > 0.0, column_scale, 1.0)
    a_scaled = a / column_scale[None, :]
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_scaled, b, rcond=1.0e-12)
    coeff = coeff_scaled / column_scale
    residual_norm = float(np.linalg.norm(a @ coeff - b))
    target_norm = float(np.linalg.norm(b))
    return {
        "degree_r": int(degree_r),
        "degree_p": int(degree_p),
        "fit_region": fit_region,
        "coeff": coeff,
        "column_scale": column_scale,
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "scaled_condition": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
        "weighted_residual_norm": residual_norm,
        "weighted_target_norm": target_norm,
        "weighted_relative_residual": residual_norm / max(target_norm, 1.0e-300),
        **scales,
    }


def evaluate_case(case: CaseData, fit: dict[str, object], region: str) -> dict[str, object]:
    mask = case.masks[region]
    r = case.r_scalar
    p = case.p_scalar
    f, f_r, f_p = basis_values_2d(
        r.reshape(-1),
        p.reshape(-1),
        int(fit["degree_r"]),
        int(fit["degree_p"]),
        float(fit["r0"]),
        float(fit["p0"]),
        float(fit["r_s_min"]),
        float(fit["r_s_max"]),
        float(fit["p_s_min"]),
        float(fit["p_s_max"]),
    )
    coeff = np.asarray(fit["coeff"], dtype=float)
    f = (f @ coeff).reshape(r.shape)
    f_r = (f_r @ coeff).reshape(r.shape)
    f_p = (f_p @ coeff).reshape(r.shape)
    lhs = (
        f_r[..., None, None] * case.ricci
        + 2.0 * f_p[..., None, None] * case.ricci_square_tensor
        - 0.5 * f[..., None, None] * case.metric_cov
    )
    residual = tensor_norm(lhs - case.target) / np.maximum(case.rhs_norm, 1.0e-300)
    lhs_norm = tensor_norm(lhs)
    rel_lhs_rhs = tensor_norm(lhs - case.target) / np.maximum(lhs_norm + case.rhs_norm, 1.0e-300)
    return {
        "count": int(np.count_nonzero(mask)),
        "residual_to_rhs": stats_abs(residual[mask]),
        "residual_to_lhs_plus_rhs": stats_abs(rel_lhs_rhs[mask]),
        "lhs_norm": stats_abs(lhs_norm[mask]),
        "rhs_norm": stats_abs(case.rhs_norm[mask]),
        "f": stats_abs(f[mask]),
        "f_R": stats_abs(f_r[mask]),
        "f_I2": stats_abs(f_p[mask]),
    }


def render_residual_plot(out_path: Path, cases: list[CaseData], fit: dict[str, object], region: str) -> None:
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.8), constrained_layout=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        r = case.r_scalar
        p = case.p_scalar
        f, f_r, f_p = basis_values_2d(
            r.reshape(-1),
            p.reshape(-1),
            int(fit["degree_r"]),
            int(fit["degree_p"]),
            float(fit["r0"]),
            float(fit["p0"]),
            float(fit["r_s_min"]),
            float(fit["r_s_max"]),
            float(fit["p_s_min"]),
            float(fit["p_s_max"]),
        )
        coeff = np.asarray(fit["coeff"], dtype=float)
        f = (f @ coeff).reshape(r.shape)
        f_r = (f_r @ coeff).reshape(r.shape)
        f_p = (f_p @ coeff).reshape(r.shape)
        lhs = (
            f_r[..., None, None] * case.ricci
            + 2.0 * f_p[..., None, None] * case.ricci_square_tensor
            - 0.5 * f[..., None, None] * case.metric_cov
        )
        residual = tensor_norm(lhs - case.target) / np.maximum(case.rhs_norm, 1.0e-300)
        extent = [
            float(case.x[0] * HBAR_C_EV_M * 1.0e6),
            float(case.x[-1] * HBAR_C_EV_M * 1.0e6),
            float(case.z[0] * HBAR_C_EV_M * 1.0e6),
            float(case.z[-1] * HBAR_C_EV_M * 1.0e6),
        ]
        im = ax.imshow(np.log10(1.0e-6 + np.clip(residual, 0.0, 1.0e6)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        ax.contour(case.x * HBAR_C_EV_M * 1.0e6, case.z * HBAR_C_EV_M * 1.0e6, case.masks[region].T.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        ax.set_title(f"tau={case.tau_old:g}; log10 residual")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_coeff_heatmap(out_path: Path, fits: list[dict[str, object]]) -> None:
    cols = len(fits)
    fig, axes = plt.subplots(1, cols, figsize=(4.6 * cols, 4.2), constrained_layout=True)
    if cols == 1:
        axes = [axes]
    for ax, fit in zip(axes, fits):
        coeff = np.asarray(fit["coeff"], dtype=float).reshape(int(fit["degree_r"]) + 1, int(fit["degree_p"]) + 1)
        im = ax.imshow(np.sign(coeff) * np.log10(1.0 + np.abs(coeff)), origin="lower", aspect="auto", cmap="coolwarm")
        ax.set_title(f"deg R={fit['degree_r']}, I2={fit['degree_p']}")
        ax.set_xlabel("I2 Cheb index")
        ax.set_ylabel("R Cheb index")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    if args.ell is None:
        args.ell = float(args.ell_over_planck) * PLANCK_LENGTH_EV_INV

    ref = build_physical_reference(args)
    cases = [build_case(args, ref, tau) for tau in [float(x) for x in args.taus.split(",") if x.strip()]]
    degrees = [tuple(int(part) for part in item.split(":")) for item in args.degrees.split(",") if item.strip()]
    fits = [
        fit_universal(cases, args.fit_region, degree_r, degree_p, float(args.target_floor_frac))
        for degree_r, degree_p in degrees
    ]

    report_fits: list[dict[str, object]] = []
    for fit in fits:
        entry: dict[str, object] = {k: v for k, v in fit.items() if k not in {"coeff", "column_scale"}}
        entry["coeff_abs"] = stats_abs(np.asarray(fit["coeff"]))
        entry["column_scale_abs"] = stats_abs(np.asarray(fit["column_scale"]))
        entry["by_case"] = {}
        for case in cases:
            entry["by_case"][f"tau={case.tau_old:g}"] = {
                region: evaluate_case(case, fit, region)
                for region in ["trusted", "core1", "core10"]
            }
        report_fits.append(entry)

    best_fit = fits[-1]
    residual_plot = out / "universal_f_of_R_I2_algebraic_residual_maps.png"
    coeff_plot = out / "universal_f_of_R_I2_algebraic_coefficients.png"
    render_residual_plot(residual_plot, cases, best_fit, args.fit_region)
    render_coeff_heatmap(coeff_plot, fits)

    np.savez_compressed(
        out / "universal_f_of_R_I2_algebraic_fit_coefficients.npz",
        **{f"degree_{fit['degree_r']}_{fit['degree_p']}_coeff": np.asarray(fit["coeff"]) for fit in fits},
        **{f"degree_{fit['degree_r']}_{fit['degree_p']}_column_scale": np.asarray(fit["column_scale"]) for fit in fits},
        degrees=np.asarray(degrees),
    )
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "dt_old": float(args.dt_old),
            "taus": [case.tau_old for case in cases],
            "fit_region": args.fit_region,
            "degrees": [list(d) for d in degrees],
            "mp": float(args.mp),
            "ell": float(args.ell),
            "ell_over_planck_length": float(args.ell) / PLANCK_LENGTH_EV_INV,
            "target_floor_frac": float(args.target_floor_frac),
        },
        "definition": {
            "f_model": "f(R,I2)=sum c_nm T_n(t(asinh(R/R0))) T_m(t(asinh(I2/I20))). f_R and f_I2 are analytic derivatives of the same fitted function.",
            "I2": "I2=R_mn R^mn, with contractions using gtilde inverse.",
            "algebraic_lhs": "f_R Ricci_mn + 2 f_I2 Ricci_malpha Ricci^alpha_n - 0.5 f g_mn.",
            "scope": "This script fits only the algebraic part of metric f(R,I2). The full metric variation also contains derivative terms of f_R and f_I2 R_mn.",
        },
        "fits": report_fits,
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "coefficients_npz": str((out / "universal_f_of_R_I2_algebraic_fit_coefficients.npz").resolve()),
            "residual_plot": str(residual_plot.resolve()),
            "coefficients_plot": str(coeff_plot.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--fit-region", choices=["trusted", "core1", "core10"], default="core10")
    parser.add_argument("--degrees", type=str, default="4:4,6:6,8:8,10:10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=640)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
