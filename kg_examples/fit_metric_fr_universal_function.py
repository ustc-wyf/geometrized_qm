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

from analyze_bcd_residuals_from_a_reference import (
    covariant_hessian_from_jets,
    geometry_data,
    stress_tensor_tilde,
)
from analyze_c_terms_from_a_reference import metric_jets_full, scalar_jets_full
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import (
    build_physical_reference,
    make_cropped_snapshot,
    real_array,
)
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det


@dataclass
class CaseData:
    tau_old: float
    x: np.ndarray
    z: np.ndarray
    rho: np.ndarray
    metric_cov: np.ndarray
    ricci: np.ndarray
    r_scalar: np.ndarray
    b2: np.ndarray
    b3: np.ndarray
    target: np.ndarray
    rhs_norm: np.ndarray
    masks: dict[str, np.ndarray]


def stats_abs(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = np.abs(vals[np.isfinite(vals)])
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p90": 0.0, "p95": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p90": float(np.percentile(vals, 90.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }


def tensor_norm(tensor: np.ndarray) -> np.ndarray:
    return np.sqrt(np.einsum("...ab,...ab->...", tensor, tensor, optimize=True))


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
        -2: make_snapshot_pack(args, time_ev_inv - 2.0 * dt, ix, iz, ref),
        -1: make_snapshot_pack(args, time_ev_inv - dt, ix, iz, ref),
        0: make_snapshot_pack(args, time_ev_inv, ix, iz, ref),
        1: make_snapshot_pack(args, time_ev_inv + dt, ix, iz, ref),
        2: make_snapshot_pack(args, time_ev_inv + 2.0 * dt, ix, iz, ref),
    }

    metric_m2 = real_array(snaps[-2]["cov_txz"])
    metric_m = real_array(snaps[-1]["cov_txz"])
    metric_0 = real_array(snaps[0]["cov_txz"])
    metric_p = real_array(snaps[1]["cov_txz"])
    metric_p2 = real_array(snaps[2]["cov_txz"])

    dg_0, d2g_0 = metric_jets_full(metric_m, metric_0, metric_p, dt, dx, dz)
    metric_inv, gamma2, ricci, r_scalar = geometry_data(metric_0, dg_0, d2g_0)

    dg_m, d2g_m = metric_jets_full(metric_m2, metric_m, metric_0, dt, dx, dz)
    _, _, _, r_scalar_m = geometry_data(metric_m, dg_m, d2g_m)
    dg_p, d2g_p = metric_jets_full(metric_0, metric_p, metric_p2, dt, dx, dz)
    _, _, _, r_scalar_p = geometry_data(metric_p, dg_p, d2g_p)

    dR, d2R = scalar_jets_full(r_scalar_m, r_scalar, r_scalar_p, dt, dx, dz)
    cov_hess_R = covariant_hessian_from_jets(dR, d2R, gamma2)
    box_R = np.einsum("...ab,...ab->...", metric_inv, cov_hess_R, optimize=True)
    grad_R_sq = np.einsum("...ab,...a,...b->...", metric_inv, dR, dR, optimize=True)
    b2 = metric_0 * box_R[..., None, None] - cov_hess_R
    b3 = metric_0 * grad_R_sq[..., None, None] - np.einsum("...a,...b->...ab", dR, dR, optimize=True)

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
        r_scalar=r_scalar,
        b2=b2,
        b3=b3,
        target=target,
        rhs_norm=rhs_norm,
        masks=masks,
    )


def cheb_derivative_vandermonde(t: np.ndarray, degree: int, deriv: int) -> np.ndarray:
    out = np.empty((t.size, degree + 1), dtype=float)
    for n in range(degree + 1):
        poly = cheb.Chebyshev.basis(n).deriv(deriv)
        out[:, n] = poly(t)
    return out


def basis_values_for_r(
    r: np.ndarray,
    degree: int,
    r0: float,
    s_min: float,
    s_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r = np.asarray(r, dtype=float)
    s = np.arcsinh(r / r0)
    scale = 2.0 / max(s_max - s_min, 1.0e-300)
    t = scale * (s - s_min) - 1.0
    t = np.clip(t, -1.0, 1.0)

    v0 = cheb.chebvander(t, degree)
    vt1 = cheb_derivative_vandermonde(t, degree, 1)
    vt2 = cheb_derivative_vandermonde(t, degree, 2)
    vt3 = cheb_derivative_vandermonde(t, degree, 3)

    vs1 = vt1 * scale
    vs2 = vt2 * scale * scale
    vs3 = vt3 * scale * scale * scale

    denom2 = r * r + r0 * r0
    j = 1.0 / np.sqrt(denom2)
    j_r = -r / np.power(denom2, 1.5)
    j_rr = (2.0 * r * r - r0 * r0) / np.power(denom2, 2.5)

    f = v0
    f_r = j[:, None] * vs1
    f_rr = (j * j)[:, None] * vs2 + j_r[:, None] * vs1
    f_rrr = (j**3)[:, None] * vs3 + (3.0 * j * j_r)[:, None] * vs2 + j_rr[:, None] * vs1
    return f, f_r, f_rr, f_rrr


def design_rows_for_case(
    case: CaseData,
    mask: np.ndarray,
    degree: int,
    r0: float,
    s_min: float,
    s_max: float,
    target_floor_frac: float,
) -> tuple[np.ndarray, np.ndarray]:
    idx = np.where(mask.reshape(-1))[0]
    shape = case.r_scalar.shape
    r = case.r_scalar.reshape(-1)[idx]
    rho = case.rho.reshape(-1)[idx]
    rhs_norm = case.rhs_norm.reshape(-1)[idx]
    target = case.target.reshape(-1, 9)[idx]
    metric = case.metric_cov.reshape(-1, 9)[idx]
    ricci = case.ricci.reshape(-1, 9)[idx]
    b2 = case.b2.reshape(-1, 9)[idx]
    b3 = case.b3.reshape(-1, 9)[idx]

    f, f_r, f_rr, f_rrr = basis_values_for_r(r, degree, r0, s_min, s_max)
    n = idx.size
    k = degree + 1
    rows = np.zeros((n * 9, k), dtype=float)
    target_rows = target.reshape(n * 9)
    for comp in range(9):
        block = (
            f_r * ricci[:, comp, None]
            - 0.5 * f * metric[:, comp, None]
            + f_rr * b2[:, comp, None]
            + f_rrr * b3[:, comp, None]
        )
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
    degree: int,
    target_floor_frac: float,
) -> dict[str, object]:
    r_fit = np.concatenate([case.r_scalar[case.masks[fit_region]] for case in cases])
    r_abs = np.abs(r_fit[np.isfinite(r_fit)])
    r0 = max(float(np.percentile(r_abs, 50.0)), 1.0e-30)
    s_all = np.arcsinh(r_fit / r0)
    s_min = float(np.percentile(s_all[np.isfinite(s_all)], 0.5))
    s_max = float(np.percentile(s_all[np.isfinite(s_all)], 99.5))
    if s_max <= s_min:
        s_min = float(np.min(s_all))
        s_max = float(np.max(s_all))

    matrices = [
        design_rows_for_case(case, case.masks[fit_region], degree, r0, s_min, s_max, target_floor_frac)
        for case in cases
    ]
    a = np.vstack([m[0] for m in matrices])
    b = np.concatenate([m[1] for m in matrices])
    column_scale = np.linalg.norm(a, axis=0)
    column_scale = np.where(column_scale > 0.0, column_scale, 1.0)
    a_scaled = a / column_scale[None, :]
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_scaled, b, rcond=1.0e-12)
    coeff = coeff_scaled / column_scale
    weighted_residual_norm = float(np.linalg.norm(a @ coeff - b))
    weighted_target_norm = float(np.linalg.norm(b))
    return {
        "degree": int(degree),
        "fit_region": fit_region,
        "coeff": coeff,
        "column_scale": column_scale,
        "r0": float(r0),
        "s_min": float(s_min),
        "s_max": float(s_max),
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "scaled_condition": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
        "weighted_residual_norm": weighted_residual_norm,
        "weighted_target_norm": weighted_target_norm,
        "weighted_relative_residual": weighted_residual_norm / max(weighted_target_norm, 1.0e-300),
    }


def evaluate_case(case: CaseData, fit: dict[str, object], region: str) -> dict[str, object]:
    mask = case.masks[region]
    r = case.r_scalar
    f, f_r, f_rr, f_rrr = basis_values_for_r(
        r.reshape(-1),
        int(fit["degree"]),
        float(fit["r0"]),
        float(fit["s_min"]),
        float(fit["s_max"]),
    )
    coeff = np.asarray(fit["coeff"], dtype=float)
    f = (f @ coeff).reshape(r.shape)
    f_r = (f_r @ coeff).reshape(r.shape)
    f_rr = (f_rr @ coeff).reshape(r.shape)
    f_rrr = (f_rrr @ coeff).reshape(r.shape)
    lhs = (
        f_r[..., None, None] * case.ricci
        - 0.5 * f[..., None, None] * case.metric_cov
        + f_rr[..., None, None] * case.b2
        + f_rrr[..., None, None] * case.b3
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
        "f_RR": stats_abs(f_rr[mask]),
        "f_RRR": stats_abs(f_rrr[mask]),
    }


def render_function_plot(out_path: Path, fits: list[dict[str, object]], r_values: np.ndarray) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.4), constrained_layout=True)
    labels = ["f", "f_R", "f_RR", "f_RRR"]
    for fit in fits:
        r = np.linspace(float(np.percentile(r_values, 1.0)), float(np.percentile(r_values, 99.0)), 1000)
        vals = basis_values_for_r(r, int(fit["degree"]), float(fit["r0"]), float(fit["s_min"]), float(fit["s_max"]))
        coeff = np.asarray(fit["coeff"], dtype=float)
        for ax, label, val in zip(axes.ravel(), labels, vals):
            y = val @ coeff
            ax.plot(r, y, label=f"deg {fit['degree']}", linewidth=1.2)
            ax.set_xscale("symlog", linthresh=1.0e-20)
            ax.set_yscale("symlog", linthresh=1.0e-80)
            ax.set_xlabel("R_tilde")
            ax.set_ylabel(label)
            ax.set_title(label + "(R_tilde)")
    axes[0, 0].legend(fontsize=8)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_residual_plot(out_path: Path, cases: list[CaseData], fit: dict[str, object], region: str) -> None:
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.8), constrained_layout=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        r = case.r_scalar
        f, f_r, f_rr, f_rrr = basis_values_for_r(
            r.reshape(-1),
            int(fit["degree"]),
            float(fit["r0"]),
            float(fit["s_min"]),
            float(fit["s_max"]),
        )
        coeff = np.asarray(fit["coeff"], dtype=float)
        f = (f @ coeff).reshape(r.shape)
        f_r = (f_r @ coeff).reshape(r.shape)
        f_rr = (f_rr @ coeff).reshape(r.shape)
        f_rrr = (f_rrr @ coeff).reshape(r.shape)
        lhs = (
            f_r[..., None, None] * case.ricci
            - 0.5 * f[..., None, None] * case.metric_cov
            + f_rr[..., None, None] * case.b2
            + f_rrr[..., None, None] * case.b3
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


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    if args.ell is None:
        args.ell = float(args.ell_over_planck) * PLANCK_LENGTH_EV_INV

    ref = build_physical_reference(args)
    cases = [build_case(args, ref, tau) for tau in [float(x) for x in args.taus.split(",") if x.strip()]]
    degrees = [int(x) for x in args.degrees.split(",") if x.strip()]
    fits = [fit_universal(cases, args.fit_region, degree, float(args.target_floor_frac)) for degree in degrees]

    report_fits: list[dict[str, object]] = []
    for fit in fits:
        fit_entry: dict[str, object] = {k: v for k, v in fit.items() if k != "coeff"}
        fit_entry.pop("column_scale", None)
        fit_entry["coeff_abs"] = stats_abs(np.asarray(fit["coeff"]))
        fit_entry["column_scale_abs"] = stats_abs(np.asarray(fit["column_scale"]))
        fit_entry["by_case"] = {}
        for case in cases:
            fit_entry["by_case"][f"tau={case.tau_old:g}"] = {
                region: evaluate_case(case, fit, region)
                for region in ["trusted", "core1", "core10"]
            }
        report_fits.append(fit_entry)

    best_fit = fits[-1]
    r_all = np.concatenate([case.r_scalar[case.masks[args.fit_region]] for case in cases])
    function_plot = out / "universal_f_of_R_fit_functions.png"
    residual_plot = out / "universal_f_of_R_residual_maps.png"
    render_function_plot(function_plot, fits, r_all)
    render_residual_plot(residual_plot, cases, best_fit, args.fit_region)

    np.savez_compressed(
        out / "universal_f_of_R_fit_coefficients.npz",
        **{f"degree_{fit['degree']}_coeff": np.asarray(fit["coeff"]) for fit in fits},
        **{f"degree_{fit['degree']}_column_scale": np.asarray(fit["column_scale"]) for fit in fits},
        degrees=np.asarray(degrees),
        r0=np.asarray([fit["r0"] for fit in fits]),
        s_min=np.asarray([fit["s_min"] for fit in fits]),
        s_max=np.asarray([fit["s_max"] for fit in fits]),
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
            "degrees": degrees,
            "mp": float(args.mp),
            "ell": float(args.ell),
            "ell_over_planck_length": float(args.ell) / PLANCK_LENGTH_EV_INV,
            "target_floor_frac": float(args.target_floor_frac),
        },
        "definition": {
            "f_model": "f(R)=sum c_n T_n(t(asinh(R/R0))). f_R, f_RR, f_RRR are analytic R-derivatives of the same fitted function.",
            "field_equation_lhs": "f_R Ricci_mn - 0.5 f g_mn + f_RR*(g BoxR-HessR)_mn + f_RRR*(g gradR^2-dR dR)_mn",
            "fit": "single universal f(R) is fitted simultaneously across all requested tau slices on fit_region.",
        },
        "fits": report_fits,
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "coefficients_npz": str((out / "universal_f_of_R_fit_coefficients.npz").resolve()),
            "function_plot": str(function_plot.resolve()),
            "residual_plot": str(residual_plot.resolve()),
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
    parser.add_argument("--degrees", type=str, default="6,10,14,18")
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
