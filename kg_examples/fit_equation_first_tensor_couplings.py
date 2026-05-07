from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import (
    CaseData,
    build_case,
    lower_mixed_square,
    symmetric_rows,
    tensor_norm,
)
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det


SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }
    if weights is not None:
        w = np.asarray(weights, dtype=float)[mask]
        w = np.maximum(w, 0.0)
        if float(np.sum(w)) > 0.0:
            out["weighted_mean"] = float(np.sum(w * vals) / np.sum(w))
    return out


def sym_outer(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return 0.5 * (
        np.einsum("...a,...b->...ab", a, b, optimize=True)
        + np.einsum("...a,...b->...ab", b, a, optimize=True)
    )


def normal_covector(case: CaseData) -> np.ndarray:
    # 2+1d covariant pseudo-normal to the u-r plane:
    # n_mu = sqrt(|g|) eps_{mu nu lambda} u^nu r^lambda.
    u_up = np.einsum("...ab,...b->...a", case.metric_inv, case.u_cov, optimize=True)
    r_up = np.einsum("...ab,...b->...a", case.metric_inv, case.r_cov, optimize=True)
    eps = np.zeros((3, 3, 3), dtype=float)
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1.0
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1.0
    return safe_sqrt_abs_det(case.metric_cov)[..., None] * np.einsum("abc,...b,...c->...a", eps, u_up, r_up, optimize=True)


def invariant_i3(case: CaseData) -> np.ndarray:
    mixed = np.einsum("...ac,...cb->...ab", case.metric_inv, case.ricci, optimize=True)
    mixed2 = np.einsum("...ac,...cb->...ab", mixed, mixed, optimize=True)
    mixed3 = np.einsum("...ac,...cb->...ab", mixed2, mixed, optimize=True)
    return np.trace(mixed3, axis1=-2, axis2=-1)


def scalar_features(case: CaseData, family: str) -> dict[str, np.ndarray]:
    r = case.r_scalar
    ricci_sq = lower_mixed_square(case.ricci, case.metric_inv)
    i2 = np.einsum("...ab,...ab->...", case.metric_inv, ricci_sq, optimize=True)
    i3 = invariant_i3(case)
    if family == "constant":
        return {"1": np.ones_like(r)}
    if family == "linear_RI2":
        return {"1": np.ones_like(r), "R": r, "I2": i2}
    if family == "quadratic_RI2":
        return {
            "1": np.ones_like(r),
            "R": r,
            "I2": i2,
            "R2": r * r,
            "R_I2": r * i2,
            "I2_2": i2 * i2,
        }
    if family == "cubic_RI2I3":
        return {
            "1": np.ones_like(r),
            "R": r,
            "I2": i2,
            "I3": i3,
            "R2": r * r,
            "R_I2": r * i2,
            "I2_2": i2 * i2,
            "R3": r * r * r,
            "R2_I2": r * r * i2,
            "R_I2_2": r * i2 * i2,
            "I2_3": i2 * i2 * i2,
        }
    raise ValueError(f"unknown scalar feature family: {family}")


def tensor_atoms(case: CaseData, family: str) -> dict[str, np.ndarray]:
    u = case.u_cov
    r = case.r_cov
    n = normal_covector(case)
    atoms = {
        "g": case.metric_cov,
        "uu": sym_outer(u, u),
        "rr": sym_outer(r, r),
        "ur": sym_outer(u, r),
    }
    if family == "matter4":
        return atoms
    if family == "matter6_normal":
        return {
            "uu": sym_outer(u, u),
            "rr": sym_outer(r, r),
            "nn": sym_outer(n, n),
            "ur": sym_outer(u, r),
            "un": sym_outer(u, n),
            "rn": sym_outer(r, n),
        }
    if family == "matter4_plus_normal":
        atoms["nn"] = sym_outer(n, n)
        return atoms
    raise ValueError(f"unknown tensor atom family: {family}")


def expanded_basis(case: CaseData, atom_family: str, feature_family: str) -> dict[str, np.ndarray]:
    atoms = tensor_atoms(case, atom_family)
    features = scalar_features(case, feature_family)
    basis: dict[str, np.ndarray] = {}
    for fname, fval in features.items():
        for aname, aval in atoms.items():
            basis[f"{fname}*{aname}"] = fval[..., None, None] * aval
    return basis


def weights_for(case: CaseData, idx: np.ndarray, target_rows: np.ndarray, target_floor_frac: float) -> np.ndarray:
    target_norm = np.sqrt(np.sum(target_rows**2, axis=1))
    finite = target_norm[np.isfinite(target_norm)]
    floor = float(target_floor_frac) * max(float(np.percentile(finite, 95.0)) if finite.size else 0.0, 1.0e-300)
    rho_w = np.sqrt(np.maximum(case.rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    return rho_w / np.maximum(target_norm, floor)


def fit_global(
    cases: list[CaseData],
    atom_family: str,
    feature_family: str,
    region: str,
    target_floor_frac: float,
) -> dict[str, object]:
    names: list[str] | None = None
    matrices: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    for case in cases:
        mask = getattr(case, region)
        idx = np.where(mask.reshape(-1))[0]
        target_rows = symmetric_rows(case.r_need, idx)
        w = weights_for(case, idx, target_rows, target_floor_frac)
        basis = expanded_basis(case, atom_family, feature_family)
        if names is None:
            names = list(basis.keys())
        rows = np.stack([symmetric_rows(basis[name], idx) for name in names], axis=-1)
        a = rows.reshape(idx.size * len(SYMMETRIC_COMPONENTS), len(names))
        b = target_rows.reshape(idx.size * len(SYMMETRIC_COMPONENTS))
        w_rows = np.repeat(w, len(SYMMETRIC_COMPONENTS))
        matrices.append(a * w_rows[:, None])
        targets.append(b * w_rows)
    if names is None:
        raise RuntimeError("no fit rows")
    a_all = np.vstack(matrices)
    b_all = np.concatenate(targets)
    col_scale = np.linalg.norm(a_all, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_all / col_scale[None, :], b_all, rcond=1.0e-12)
    coeff = coeff_scaled / col_scale
    residual = a_all @ coeff - b_all
    return {
        "model": f"{feature_family}:{atom_family}",
        "basis_names": names,
        "coeff": coeff,
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "scaled_condition": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
        "weighted_relative_residual": float(np.linalg.norm(residual) / max(np.linalg.norm(b_all), 1.0e-300)),
    }


def predict_global(case: CaseData, fit: dict[str, object], atom_family: str, feature_family: str) -> np.ndarray:
    basis = expanded_basis(case, atom_family, feature_family)
    coeff = np.asarray(fit["coeff"], dtype=float)
    out = np.zeros_like(case.r_need)
    for c, name in zip(coeff, fit["basis_names"]):
        out += float(c) * basis[str(name)]
    return out


def evaluate_global(case: CaseData, fit: dict[str, object], atom_family: str, feature_family: str, region: str) -> dict[str, object]:
    pred = predict_global(case, fit, atom_family, feature_family)
    rel = tensor_norm(pred - case.r_need) / np.maximum(tensor_norm(case.r_need), 1.0e-300)
    mask = getattr(case, region)
    w = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    return {
        "count": int(np.count_nonzero(mask)),
        "relative_residual": stats(rel[mask], w),
        "target_norm": stats(tensor_norm(case.r_need)[mask], w),
        "prediction_norm": stats(tensor_norm(pred)[mask], w),
    }


def fit_pointwise(case: CaseData, atom_family: str, region: str) -> dict[str, object]:
    mask = getattr(case, region)
    idx = np.where(mask.reshape(-1))[0]
    basis = tensor_atoms(case, atom_family)
    names = list(basis.keys())
    rows = np.stack([symmetric_rows(basis[name], idx) for name in names], axis=-1)
    target = symmetric_rows(case.r_need, idx)
    residual = np.full(idx.size, np.nan, dtype=float)
    cond = np.full(idx.size, np.nan, dtype=float)
    for p in range(idx.size):
        a = rows[p]
        b = target[p]
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        sol, *_ = np.linalg.lstsq(a, b, rcond=1.0e-12)
        fitted = a @ sol
        residual[p] = float(np.linalg.norm(fitted - b) / max(np.linalg.norm(b), 1.0e-300))
        svals = np.linalg.svd(a, compute_uv=False)
        cond[p] = float(np.max(svals) / max(np.min(svals), 1.0e-300))
    w = np.sqrt(np.maximum(case.rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    return {
        "atom_family": atom_family,
        "region": region,
        "count": int(idx.size),
        "pointwise_relative_residual": stats(residual, w),
        "pointwise_condition": stats(cond, w),
    }


def render_summary_plot(out_path: Path, cases: list[CaseData], fit: dict[str, object], atom_family: str, feature_family: str, region: str) -> None:
    fig, axes = plt.subplots(1, len(cases), figsize=(5.2 * len(cases), 4.8), constrained_layout=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, case in zip(axes, cases):
        pred = predict_global(case, fit, atom_family, feature_family)
        rel = tensor_norm(pred - case.r_need) / np.maximum(tensor_norm(case.r_need), 1.0e-300)
        extent = [
            float(case.x[0] * HBAR_C_EV_M * 1.0e6),
            float(case.x[-1] * HBAR_C_EV_M * 1.0e6),
            float(case.z[0] * HBAR_C_EV_M * 1.0e6),
            float(case.z[-1] * HBAR_C_EV_M * 1.0e6),
        ]
        im = ax.imshow(np.log10(1.0e-6 + np.clip(rel, 0.0, 1.0e6)).T, origin="lower", extent=extent, aspect="equal", cmap="viridis")
        ax.contour(case.x * HBAR_C_EV_M * 1.0e6, case.z * HBAR_C_EV_M * 1.0e6, getattr(case, region).T.astype(float), levels=[0.5], colors="white", linewidths=0.6)
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
    ref = build_physical_reference(args)
    cases = [build_case(args, ref, tau) for tau in [float(x) for x in args.taus.split(",") if x.strip()]]

    specs = [
        ("matter4", "constant"),
        ("matter4", "linear_RI2"),
        ("matter4", "quadratic_RI2"),
        ("matter4", "cubic_RI2I3"),
        ("matter4_plus_normal", "quadratic_RI2"),
        ("matter4_plus_normal", "cubic_RI2I3"),
    ]
    global_fits: list[dict[str, object]] = []
    for atom_family, feature_family in specs:
        fit = fit_global(cases, atom_family, feature_family, args.fit_region, float(args.target_floor_frac))
        entry = {k: v for k, v in fit.items() if k != "coeff"}
        entry["coeff_abs"] = stats(np.abs(np.asarray(fit["coeff"])))
        entry["by_case"] = {
            f"tau={case.tau_old:g}": evaluate_global(case, fit, atom_family, feature_family, args.fit_region)
            for case in cases
        }
        global_fits.append(entry)
        np.savez_compressed(out / f"coeff_{feature_family}_{atom_family}.npz", coeff=np.asarray(fit["coeff"]), basis_names=np.asarray(fit["basis_names"], dtype=object))

    pointwise = []
    for atom_family in ["matter4", "matter4_plus_normal", "matter6_normal"]:
        pointwise.append(
            {
                f"tau={case.tau_old:g}": fit_pointwise(case, atom_family, args.fit_region)
                for case in cases
            }
        )

    best_index = int(np.argmin([float(item["weighted_relative_residual"]) for item in global_fits]))
    best_atom, best_feature = specs[best_index]
    best_fit = fit_global(cases, best_atom, best_feature, args.fit_region, float(args.target_floor_frac))
    plot_path = out / "equation_first_best_residual_maps.png"
    render_summary_plot(plot_path, cases, best_fit, best_atom, best_feature, args.fit_region)

    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "dt_old": float(args.dt_old),
            "taus": [case.tau_old for case in cases],
            "fit_region": args.fit_region,
            "mp": float(args.mp),
            "target_floor_frac": float(args.target_floor_frac),
        },
        "definition": {
            "equation_form": "Gtilde_mn - Ttilde_mn/Mp^2 = C_mn[gtilde,u,r].",
            "target": "r_need_mn = Gtilde_mn - Ttilde_mn/Mp^2, evaluated on A-reference-generated gtilde.",
            "global_models": "C_mn is expanded in tensor atoms {g, uu, rr, ur, optionally nn} with scalar coefficients polynomial in local curvature invariants {R,I2,I3}. Coefficients are universal constants fitted across all requested time slices.",
            "pointwise_models": "Diagnostic lower bounds where coefficients are allowed to vary pointwise; these are not final closed equations unless the coefficients are later shown to be scalar functions of local invariants.",
        },
        "global_fits": global_fits,
        "pointwise_diagnostics": pointwise,
        "outputs": {
            "report_json": str((out / "summary.json").resolve()),
            "best_residual_plot": str(plot_path.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=160)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
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
