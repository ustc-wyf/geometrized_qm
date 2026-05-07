from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fit_metric_fr_ricci2_universal_full import (
    CaseFullData,
    build_case,
    full_lhs_from_fields,
)
from physical_units import HBAR_C_EV_M, PLANCK_LENGTH_EV_INV, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from diagnose_metric_fr_direction_matching import stats_abs, tensor_norm


# Minimal local finite-jet basis whose highest-order terms are 6-derivative
# scalar invariants in 2+1d.
BASIS_SPECS: list[tuple[str, int, int]] = [
    ("1", 0, 0),
    ("R", 1, 0),
    ("P", 0, 1),
    ("R2", 2, 0),
    ("RP", 1, 1),
    ("R3", 3, 0),
]


def monomial_fields(r: np.ndarray, p: np.ndarray, a: int, b: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = np.asarray(r, dtype=float)
    p = np.asarray(p, dtype=float)
    f = np.power(r, a) * np.power(p, b)
    if a > 0:
        fr = a * np.power(r, a - 1) * np.power(p, b)
    else:
        fr = np.zeros_like(f)
    if b > 0:
        fy = b * np.power(r, a) * np.power(p, b - 1)
    else:
        fy = np.zeros_like(f)
    return f, fr, fy


def basis_lhs_tensors(case: CaseFullData, dt: float, dx: float, dz: float) -> list[np.ndarray]:
    tensors: list[np.ndarray] = []
    for _, a, b in BASIS_SPECS:
        f_m, fr_m, fy_m = monomial_fields(case.geom_m.r_scalar, case.geom_m.p_scalar, a, b)
        f_0, fr_0, fy_0 = monomial_fields(case.geom_0.r_scalar, case.geom_0.p_scalar, a, b)
        f_p, fr_p, fy_p = monomial_fields(case.geom_p.r_scalar, case.geom_p.p_scalar, a, b)
        parts = full_lhs_from_fields(
            case.geom_m,
            case.geom_0,
            case.geom_p,
            case.dgamma_0,
            f_m,
            f_0,
            f_p,
            fr_m,
            fr_0,
            fr_p,
            fy_m,
            fy_0,
            fy_p,
            dt,
            dx,
            dz,
        )
        tensors.append(parts["lhs"])
    return tensors


def design_rows_for_case(
    case: CaseFullData,
    mask: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    target_floor_frac: float,
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    idx = np.where(mask.reshape(-1))[0]
    target = case.target.reshape(-1, 9)[idx]
    rhs_norm = case.rhs_norm.reshape(-1)[idx]
    rho = case.rho.reshape(-1)[idx]
    rho_w = np.sqrt(np.maximum(rho, 0.0) / max(float(np.max(case.rho)), 1.0e-300))
    floor = float(target_floor_frac) * max(float(np.percentile(rhs_norm[np.isfinite(rhs_norm)], 95.0)), 1.0e-300)
    rel_w = 1.0 / np.maximum(rhs_norm, floor)
    w = rho_w * rel_w
    w_rows = np.repeat(w, 9)

    tensors = basis_lhs_tensors(case, dt, dx, dz)
    basis_rows = [tensor.reshape(-1, 9)[idx] for tensor in tensors]
    a_point = np.stack(basis_rows, axis=-1)  # (npts, 9, k)
    a = a_point.reshape(idx.size * 9, len(BASIS_SPECS))
    b = target.reshape(idx.size * 9)
    return a * w_rows[:, None], b * w_rows, tensors


def fit_finitejet6(
    cases: list[CaseFullData],
    fit_region: str,
    target_floor_frac: float,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, object]:
    matrices: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    basis_cache: list[list[np.ndarray]] = []
    for case in cases:
        a, b, tensors = design_rows_for_case(case, case.masks[fit_region], dt, dx, dz, target_floor_frac)
        matrices.append(a)
        targets.append(b)
        basis_cache.append(tensors)

    a = np.vstack(matrices)
    b = np.concatenate(targets)
    column_scale = np.linalg.norm(a, axis=0)
    column_scale = np.where(column_scale > 0.0, column_scale, 1.0)
    a_scaled = a / column_scale[None, :]
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_scaled, b, rcond=1.0e-12)
    coeff = coeff_scaled / column_scale
    residual_norm = float(np.linalg.norm(a @ coeff - b))
    target_norm = float(np.linalg.norm(b))
    return {
        "basis_specs": BASIS_SPECS,
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
        "basis_cache": basis_cache,
    }


def evaluate_case(
    case: CaseFullData,
    fit: dict[str, object],
    region: str,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, object]:
    coeff = np.asarray(fit["coeff"], dtype=float)
    tensors = basis_lhs_tensors(case, dt, dx, dz)
    lhs = np.zeros_like(case.target, dtype=float)
    for value, tensor in zip(coeff, tensors):
        lhs += float(value) * tensor
    residual = tensor_norm(lhs - case.target) / np.maximum(case.rhs_norm, 1.0e-300)
    lhs_norm = tensor_norm(lhs)
    rel_lhs_rhs = tensor_norm(lhs - case.target) / np.maximum(lhs_norm + case.rhs_norm, 1.0e-300)
    mask = case.masks[region]
    return {
        "count": int(np.count_nonzero(mask)),
        "residual_to_rhs": stats_abs(residual[mask]),
        "residual_to_lhs_plus_rhs": stats_abs(rel_lhs_rhs[mask]),
        "lhs_norm": stats_abs(lhs_norm[mask]),
        "rhs_norm": stats_abs(case.rhs_norm[mask]),
    }


def render_residual_plot(
    out_path: Path,
    cases: list[CaseFullData],
    fit: dict[str, object],
    region: str,
    dt: float,
    dx: float,
    dz: float,
) -> None:
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.8), constrained_layout=True)
    if len(cases) == 1:
        axes = [axes]
    coeff = np.asarray(fit["coeff"], dtype=float)
    for ax, case in zip(axes, cases):
        tensors = basis_lhs_tensors(case, dt, dx, dz)
        lhs = np.zeros_like(case.target, dtype=float)
        for value, tensor in zip(coeff, tensors):
            lhs += float(value) * tensor
        residual = tensor_norm(lhs - case.target) / np.maximum(case.rhs_norm, 1.0e-300)
        extent = [
            float(case.x[0] * HBAR_C_EV_M * 1.0e6),
            float(case.x[-1] * HBAR_C_EV_M * 1.0e6),
            float(case.z[0] * HBAR_C_EV_M * 1.0e6),
            float(case.z[-1] * HBAR_C_EV_M * 1.0e6),
        ]
        im = ax.imshow(np.log10(1.0e-6 + np.clip(residual, 0.0, 1.0e6)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        ax.contour(
            case.x * HBAR_C_EV_M * 1.0e6,
            case.z * HBAR_C_EV_M * 1.0e6,
            case.masks[region].T.astype(float),
            levels=[0.5],
            colors="white",
            linewidths=0.6,
        )
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
    scale = ref["scale"]
    dt = float(args.dt_old) * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])

    fit = fit_finitejet6(cases, args.fit_region, float(args.target_floor_frac), dt, dx, dz)
    report_fit = {k: v for k, v in fit.items() if k not in {"coeff", "column_scale", "basis_cache"}}
    report_fit["coeff_abs"] = stats_abs(np.asarray(fit["coeff"]))
    report_fit["column_scale_abs"] = stats_abs(np.asarray(fit["column_scale"]))
    report_fit["by_case"] = {}
    for case in cases:
        report_fit["by_case"][f"tau={case.tau_old:g}"] = {
            region: evaluate_case(case, fit, region, dt, dx, dz)
            for region in ["trusted", "core1", "core10"]
        }

    residual_plot = out / "finitejet6_local_residual_maps.png"
    render_residual_plot(residual_plot, cases, fit, args.fit_region, dt, dx, dz)

    np.savez_compressed(
        out / "finitejet6_local_fit_coefficients.npz",
        coeff=np.asarray(fit["coeff"]),
        column_scale=np.asarray(fit["column_scale"]),
        basis_specs=np.asarray(BASIS_SPECS, dtype=object),
    )
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(dt),
            "taus": [case.tau_old for case in cases],
            "fit_region": args.fit_region,
            "mp": float(args.mp),
            "ell": float(args.ell),
            "ell_over_planck_length": float(args.ell) / PLANCK_LENGTH_EV_INV,
            "target_floor_frac": float(args.target_floor_frac),
        },
        "definition": {
            "basis": [
                "1",
                "R",
                "P=R_mn R^mn",
                "R^2",
                "R*P",
                "R^3",
            ],
            "interpretation": "Minimal local finite-jet basis whose highest-order terms are 6-derivative scalar invariants in 2+1d. This is action-level, not a higher-degree single-variable f(R) fit.",
            "diagnostic_scope": "Fixed-background回代/拟合诊断；尚不是完整D支动力学求解器。",
        },
        "fit": report_fit,
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "coefficients_npz": str((out / "finitejet6_local_fit_coefficients.npz").resolve()),
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
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=320)
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
