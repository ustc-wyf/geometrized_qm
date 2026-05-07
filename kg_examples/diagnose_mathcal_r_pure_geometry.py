from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import (
    build_physical_reference,
    make_cropped_snapshot,
    real_array,
)
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det, tilde_measure_from_bohm


SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


@dataclass
class CaseData:
    tau_old: float
    x: np.ndarray
    z: np.ndarray
    rho_a: np.ndarray
    measure_tilde: np.ndarray
    rho_tilde: np.ndarray
    support: np.ndarray
    trusted: np.ndarray
    core10: np.ndarray
    metric_cov: np.ndarray
    metric_inv: np.ndarray
    ricci: np.ndarray
    r_scalar: np.ndarray
    einstein: np.ndarray
    stress_tilde: np.ndarray
    source_tilde: np.ndarray
    r_need: np.ndarray
    u_cov: np.ndarray
    r_cov: np.ndarray


def tensor_norm(tensor: np.ndarray) -> np.ndarray:
    return np.sqrt(np.sum(np.asarray(tensor, dtype=float) ** 2, axis=(-2, -1)))


def stats(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    abs_vals = np.abs(vals)
    return {
        "count": int(vals.size),
        "min": float(np.min(abs_vals)),
        "p50": float(np.percentile(abs_vals, 50.0)),
        "p95": float(np.percentile(abs_vals, 95.0)),
        "max": float(np.max(abs_vals)),
    }


def symmetric_rows(tensor: np.ndarray, idx: np.ndarray) -> np.ndarray:
    flat = tensor.reshape(-1, 3, 3)
    rows = [flat[idx, a, b] for a, b in SYMMETRIC_COMPONENTS]
    return np.stack(rows, axis=1)


def lower_mixed_square(cov2: np.ndarray, metric_inv: np.ndarray) -> np.ndarray:
    mixed = np.einsum("...ac,...cb->...ab", metric_inv, cov2, optimize=True)
    return np.einsum("...ac,...cb->...ab", cov2, mixed, optimize=True)


def make_bases(case: CaseData) -> dict[str, dict[str, np.ndarray]]:
    metric = case.metric_cov
    ricci = case.ricci
    scalar = case.r_scalar
    einstein = case.einstein
    ricci_sq = lower_mixed_square(ricci, case.metric_inv)
    i2 = np.einsum("...ab,...ab->...", case.metric_inv, ricci_sq, optimize=True)

    u_outer = np.einsum("...a,...b->...ab", case.u_cov, case.u_cov, optimize=True)
    r_outer = np.einsum("...a,...b->...ab", case.r_cov, case.r_cov, optimize=True)
    ur_sym = 0.5 * (
        np.einsum("...a,...b->...ab", case.u_cov, case.r_cov, optimize=True)
        + np.einsum("...a,...b->...ab", case.r_cov, case.u_cov, optimize=True)
    )

    return {
        "pure_geometry_without_eh": {
            "gtilde_mn": metric,
            "Rtilde_mn": ricci,
            "Rtilde_gtilde_mn": scalar[..., None, None] * metric,
            "Rtilde2_mn": ricci_sq,
            "I2_gtilde_mn": i2[..., None, None] * metric,
        },
        "pure_geometry_with_eh": {
            "Gtilde_mn": einstein,
            "gtilde_mn": metric,
            "Rtilde2_mn": ricci_sq,
            "I2_gtilde_mn": i2[..., None, None] * metric,
        },
        "matter_directions": {
            "gtilde_mn": metric,
            "u_m_u_n": u_outer,
            "r_m_r_n": r_outer,
            "u_(m_r_n)": ur_sym,
        },
        "combined_geom_matter": {
            "Gtilde_mn": einstein,
            "gtilde_mn": metric,
            "Rtilde2_mn": ricci_sq,
            "I2_gtilde_mn": i2[..., None, None] * metric,
            "u_m_u_n": u_outer,
            "r_m_r_n": r_outer,
            "u_(m_r_n)": ur_sym,
        },
    }


def weighted_fit(
    *,
    cases: list[CaseData],
    target_name: str,
    basis_name: str,
    region: str,
    target_floor_frac: float,
) -> dict[str, object]:
    matrices: list[np.ndarray] = []
    targets: list[np.ndarray] = []
    point_records: list[tuple[int, np.ndarray, np.ndarray]] = []
    basis_names: list[str] | None = None

    for case_index, case in enumerate(cases):
        mask = getattr(case, region)
        idx = np.where(mask.reshape(-1))[0]
        target_tensor = getattr(case, target_name)
        target_rows = symmetric_rows(target_tensor, idx)
        target_norm = np.sqrt(np.sum(target_rows**2, axis=1))
        finite_norm = target_norm[np.isfinite(target_norm)]
        floor = float(target_floor_frac) * max(
            float(np.percentile(finite_norm, 95.0)) if finite_norm.size else 0.0,
            1.0e-300,
        )
        rho_w = np.sqrt(np.maximum(case.rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
        rel_w = 1.0 / np.maximum(target_norm, floor)
        point_w = rho_w * rel_w

        bases = make_bases(case)[basis_name]
        if basis_names is None:
            basis_names = list(bases.keys())
        elif basis_names != list(bases.keys()):
            raise RuntimeError("basis order changed across cases")
        a_point = np.stack([symmetric_rows(bases[name], idx) for name in basis_names], axis=-1)
        a = a_point.reshape(idx.size * len(SYMMETRIC_COMPONENTS), len(basis_names))
        b = target_rows.reshape(idx.size * len(SYMMETRIC_COMPONENTS))
        row_w = np.repeat(point_w, len(SYMMETRIC_COMPONENTS))
        matrices.append(a * row_w[:, None])
        targets.append(b * row_w)
        point_records.append((case_index, idx, target_rows))

    if basis_names is None:
        raise RuntimeError("no fit data")
    a_all = np.vstack(matrices)
    b_all = np.concatenate(targets)
    column_scale = np.linalg.norm(a_all, axis=0)
    column_scale = np.where(column_scale > 0.0, column_scale, 1.0)
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_all / column_scale[None, :], b_all, rcond=1.0e-12)
    coeff = coeff_scaled / column_scale
    residual = a_all @ coeff - b_all
    rel = float(np.linalg.norm(residual) / max(np.linalg.norm(b_all), 1.0e-300))
    return {
        "target": target_name,
        "basis": basis_name,
        "region": region,
        "basis_names": basis_names,
        "coeff": coeff,
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "scaled_condition": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
        "weighted_relative_residual": rel,
    }


def predict(case: CaseData, fit: dict[str, object]) -> np.ndarray:
    bases = make_bases(case)[str(fit["basis"])]
    coeff = np.asarray(fit["coeff"], dtype=float)
    out = np.zeros_like(case.metric_cov)
    for value, name in zip(coeff, fit["basis_names"]):
        out += float(value) * bases[str(name)]
    return out


def evaluate_fit(case: CaseData, fit: dict[str, object], region: str) -> dict[str, object]:
    target = getattr(case, str(fit["target"]))
    model = predict(case, fit)
    rel = tensor_norm(model - target) / np.maximum(tensor_norm(target), 1.0e-300)
    mask = getattr(case, region)
    return {
        "count": int(np.count_nonzero(mask)),
        "point_relative_residual": stats(rel[mask]),
        "target_norm": stats(tensor_norm(target)[mask]),
        "model_norm": stats(tensor_norm(model)[mask]),
    }


def build_case(args: argparse.Namespace, ref: dict[str, object], tau_old: float) -> CaseData:
    params = ref["params"]
    scale = ref["scale"]
    x_full = np.asarray(ref["x_full"], dtype=float)
    z_full = np.asarray(ref["z_full"], dtype=float)
    psi0 = np.asarray(ref["psi0"])
    psi0_hat = np.asarray(ref["psi0_hat"])
    omega = np.asarray(ref["omega"])
    old_scale = float(scale.old_dimensionless_scale_ev_inv)
    t_meet_old = float(params.t_meet / old_scale)
    t = (t_meet_old + float(tau_old)) * old_scale
    dt = float(args.dt_old) * old_scale

    window_ev_inv = float(args.window_um) * 1.0e-6 / HBAR_C_EV_M
    ix = np.where(np.abs(x_full) <= window_ev_inv)[0]
    iz = np.where(np.abs(z_full) <= window_ev_inv)[0]
    if ix.size < 8 or iz.size < 8:
        raise ValueError("window contains too few grid points")
    x = x_full[ix]
    z = z_full[iz]

    snap_m = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=t - dt,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=True,
    )
    snap_0 = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=t,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=True,
    )
    snap_p = make_cropped_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x_full=x_full,
        z_full=z_full,
        ix=ix,
        iz=iz,
        t=t + dt,
        rho_floor=args.rho_floor,
        x_floor=args.x_floor,
        pinv_rcond=args.pinv_rcond,
        with_metric=True,
    )

    metric_m = real_array(snap_m["cov_txz"])
    metric_0 = real_array(snap_0["cov_txz"])
    metric_p = real_array(snap_p["cov_txz"])
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dg, d2g = metric_jets_full(metric_m, metric_0, metric_p, dt, dx, dz)
    metric_inv, _, ricci, r_scalar = geometry_data(metric_0, dg, d2g)
    einstein = ricci - 0.5 * metric_0 * r_scalar[..., None, None]

    bohm = snap_0["bohm"]
    rho_a = real_array(bohm["rho"])
    measure_tilde = tilde_measure_from_bohm(bohm, float(params.m))
    sqrt_abs_g = safe_sqrt_abs_det(metric_0)
    rho_tilde = np.maximum(measure_tilde / np.maximum(sqrt_abs_g, 1.0e-300), float(args.rho_floor))
    stress_tilde, _ = stress_tensor_tilde(
        metric_0,
        metric_inv,
        rho_tilde,
        real_array(bohm["s_t"]),
        real_array(bohm["s_x"]),
        real_array(bohm["s_z"]),
        m=float(params.m),
    )
    source_tilde = real_array(stress_tilde) / (float(args.mp) * float(args.mp))
    r_need = einstein - source_tilde

    support = (rho_a > float(args.support_rho_frac) * float(np.max(rho_a))) & (
        measure_tilde > float(args.support_measure_frac) * float(np.max(measure_tilde))
    )
    trusted = erode_mask_8(support, int(args.trusted_erosion))
    if not np.any(trusted):
        trusted = support.copy()
    core10 = support & (rho_a > 1.0e-1 * float(np.max(rho_a)))
    if not np.any(core10):
        core10 = trusted.copy()

    u_cov = np.stack([real_array(bohm["s_t"]), real_array(bohm["s_x"]), real_array(bohm["s_z"])], axis=-1)
    r_cov = np.stack([real_array(bohm["r_t"]), real_array(bohm["r_x"]), real_array(bohm["r_z"])], axis=-1)
    return CaseData(
        tau_old=float(tau_old),
        x=x,
        z=z,
        rho_a=rho_a,
        measure_tilde=measure_tilde,
        rho_tilde=rho_tilde,
        support=support,
        trusted=trusted,
        core10=core10,
        metric_cov=metric_0,
        metric_inv=metric_inv,
        ricci=real_array(ricci),
        r_scalar=real_array(r_scalar),
        einstein=real_array(einstein),
        stress_tilde=real_array(stress_tilde),
        source_tilde=source_tilde,
        r_need=real_array(r_need),
        u_cov=u_cov,
        r_cov=r_cov,
    )


def render_maps(
    out_path: Path,
    cases: list[CaseData],
    fits: dict[tuple[str, str], dict[str, object]],
    region: str,
) -> None:
    fig, axes = plt.subplots(len(cases), 4, figsize=(18.0, 4.6 * len(cases)), constrained_layout=True)
    if len(cases) == 1:
        axes = axes[None, :]
    for row, case in enumerate(cases):
        x_um = case.x * HBAR_C_EV_M * 1.0e6
        z_um = case.z * HBAR_C_EV_M * 1.0e6
        xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
        mask = getattr(case, region)
        fields = [
            ("rho_A", case.rho_a, "viridis"),
            ("|R_need|", tensor_norm(case.r_need), "magma"),
            (
                "geom-only residual for source",
                tensor_norm(predict(case, fits[("source_tilde", "pure_geometry_without_eh")]) - case.source_tilde)
                / np.maximum(tensor_norm(case.source_tilde), 1.0e-300),
                "inferno",
            ),
            (
                "matter-dir residual for source",
                tensor_norm(predict(case, fits[("source_tilde", "matter_directions")]) - case.source_tilde)
                / np.maximum(tensor_norm(case.source_tilde), 1.0e-300),
                "inferno",
            ),
        ]
        for col, (title, field, cmap) in enumerate(fields):
            data = np.asarray(field, dtype=float)
            finite = data[np.isfinite(data)]
            if title.startswith("geom") or title.startswith("matter"):
                vmax = max(float(np.percentile(data[mask & np.isfinite(data)], 95.0)) if np.any(mask) else 1.0, 1.0e-12)
                vmin = 0.0
                plot_data = np.clip(data, 0.0, vmax)
            else:
                vmax = max(float(np.percentile(np.abs(finite), 99.0)) if finite.size else 1.0, 1.0e-300)
                vmin = 0.0
                plot_data = np.clip(np.abs(data), vmin, vmax)
            im = axes[row, col].pcolormesh(xg, zg, plot_data, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            axes[row, col].contour(xg, zg, mask.astype(float), levels=[0.5], colors="white", linewidths=0.6)
            axes[row, col].set_title(f"tau={case.tau_old:g}: {title}")
            axes[row, col].set_xlabel("x [um]")
            axes[row, col].set_ylabel("z [um]")
            axes[row, col].set_aspect("equal")
            fig.colorbar(im, ax=axes[row, col], fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    cases = [build_case(args, ref, tau) for tau in [float(x) for x in args.taus.split(",") if x.strip()]]

    fit_specs = [
        ("r_need", "pure_geometry_with_eh"),
        ("r_need", "pure_geometry_without_eh"),
        ("r_need", "matter_directions"),
        ("source_tilde", "pure_geometry_without_eh"),
        ("source_tilde", "matter_directions"),
        ("source_tilde", "combined_geom_matter"),
    ]
    fits: dict[tuple[str, str], dict[str, object]] = {}
    fit_reports: list[dict[str, object]] = []
    for target_name, basis_name in fit_specs:
        fit = weighted_fit(
            cases=cases,
            target_name=target_name,
            basis_name=basis_name,
            region=args.fit_region,
            target_floor_frac=float(args.target_floor_frac),
        )
        fits[(target_name, basis_name)] = fit
        report = {key: value for key, value in fit.items() if key != "coeff"}
        report["coefficients"] = {
            name: float(value) for name, value in zip(fit["basis_names"], np.asarray(fit["coeff"], dtype=float))
        }
        report["by_case"] = {
            f"tau={case.tau_old:g}": evaluate_fit(case, fit, args.fit_region) for case in cases
        }
        fit_reports.append(report)

    map_path = args.output / "mathcal_r_pure_geometry_maps.png"
    render_maps(map_path, cases, fits, args.fit_region)

    np.savez_compressed(
        args.output / "mathcal_r_fit_coefficients.npz",
        **{
            f"{target}_{basis}_coeff": np.asarray(fit["coeff"], dtype=float)
            for (target, basis), fit in fits.items()
        },
    )
    scale = ref["scale"]
    params = ref["params"]
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "m_ev": float(params.m),
            "full_resolution": int(args.full_resolution),
            "cropped_shape": [int(cases[0].rho_a.shape[0]), int(cases[0].rho_a.shape[1])],
            "window_um": float(args.window_um),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(args.dt_old) * float(scale.old_dimensionless_scale_ev_inv),
            "taus": [case.tau_old for case in cases],
            "fit_region": args.fit_region,
            "mp_ev": float(args.mp),
            "target_floor_frac": float(args.target_floor_frac),
            "kg_norm_after": float(ref["kg_norm_after"]),
        },
        "definitions": {
            "r_need": "R_need_mn/Mp^2 = Gtilde_mn - Ttilde_mn/Mp^2. This is the directly computable extra source needed by EH[gtilde] plus the corrected tilde-matter source on the flat-A reference.",
            "source_tilde": "Ttilde_mn/Mp^2, using sqrt(|gtilde|) rho_tilde = |X| rho_A/m^2.",
            "pure_geometry_with_eh": "span{Gtilde_mn, gtilde_mn, R_mixed_square_mn, I2 gtilde_mn}; this can trivially absorb the EH piece of r_need.",
            "pure_geometry_without_eh": "span{gtilde_mn, Rtilde_mn, Rtilde gtilde_mn, R_mixed_square_mn, I2 gtilde_mn}.",
            "matter_directions": "span{gtilde_mn, u_m u_n, r_m r_n, u_(m r_n)}.",
            "caveat": "The flat Gaussian KG snapshot is not an exact finite-Mp Einstein-KG A solution. Therefore this is a structural diagnostic of the reference-generated gtilde field, not a proof for the exact transformed-A Einstein equation.",
        },
        "case_stats": {
            f"tau={case.tau_old:g}": {
                "support_count": int(np.count_nonzero(case.support)),
                "trusted_count": int(np.count_nonzero(case.trusted)),
                "core10_count": int(np.count_nonzero(case.core10)),
                "rho_a": stats(case.rho_a[case.support]),
                "Rtilde": stats(case.r_scalar[case.support]),
                "|Gtilde|": stats(tensor_norm(case.einstein)[case.support]),
                "|Ttilde|/Mp^2": stats(tensor_norm(case.source_tilde)[case.support]),
                "|R_need|": stats(tensor_norm(case.r_need)[case.support]),
            }
            for case in cases
        },
        "fits": fit_reports,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "coefficients_npz": str((args.output / "mathcal_r_fit_coefficients.npz").resolve()),
            "maps_png": str(map_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=320)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
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
    run(parser.parse_args())


if __name__ == "__main__":
    main()
