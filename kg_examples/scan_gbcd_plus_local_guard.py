from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from coordinate_matter_evolution import conservative_density_from_rho_u, coordinate_matter_rhs_covector
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_reduced_dynamic_same_initial import rk4_matter_step
from simulate_d_tridomain_full_dynamics import erode_mask_8, safe_sqrt_abs_det
from solve_gbcd_full_linear_metric_update_sparse import load_lambda_cov
from test_gbcd_v0_one_step_lambda import weighted_stats


def relative_l1(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(np.abs(candidate[valid] - reference[valid])) / denom)


def optional_array(data: np.lib.npyio.NpzFile, key: str, shape_like: np.ndarray) -> np.ndarray:
    if key in data.files:
        return np.asarray(data[key], dtype=float)
    return np.zeros_like(shape_like)


def _norm_sym(field: np.ndarray) -> np.ndarray:
    return tensor_norm(field)


def build_scale_maps(
    *,
    case,
    full_metric_plus: np.ndarray,
    delta_plus: np.ndarray,
    fit_mask: np.ndarray,
    trust_limits: tuple[float, ...],
    rho_tapers: tuple[float, ...],
) -> dict[str, np.ndarray]:
    shape = fit_mask.shape
    maps: dict[str, np.ndarray] = {}
    maps["full"] = np.ones(shape, dtype=float)
    maps["fit_only"] = fit_mask.astype(float)
    maps["fit_erode1"] = erode_mask_8(fit_mask, 1).astype(float)
    maps["fit_erode2"] = erode_mask_8(fit_mask, 2).astype(float)
    maps["trusted_only"] = np.asarray(case.trusted, dtype=float)
    maps["core10_only"] = np.asarray(case.core10, dtype=float)

    rho_frac = np.maximum(case.rho_a, 0.0) / max(float(np.max(case.rho_a)), 1.0e-300)
    for cutoff in rho_tapers:
        key = f"rho_taper_{cutoff:g}"
        maps[key] = np.clip(rho_frac / max(float(cutoff), 1.0e-300), 0.0, 1.0)

    metric_norm = _norm_sym(full_metric_plus)
    delta_norm = _norm_sym(delta_plus)
    for limit in trust_limits:
        local = float(limit) * np.maximum(metric_norm, 1.0e-300) / np.maximum(delta_norm, 1.0e-300)
        maps[f"rel_clamp_{limit:g}"] = np.minimum(1.0, local)

    for limit in trust_limits:
        local = float(limit) * np.maximum(metric_norm, 1.0e-300) / np.maximum(delta_norm, 1.0e-300)
        maps[f"fit_rel_clamp_{limit:g}"] = fit_mask.astype(float) * np.minimum(1.0, local)

    return maps


def evaluate_candidate(
    *,
    name: str,
    scale: np.ndarray,
    case,
    case_next,
    full_case,
    cov0: np.ndarray,
    fit_mask: np.ndarray,
    delta_m: np.ndarray,
    delta_0: np.ndarray,
    delta_p: np.ndarray,
    n_next: np.ndarray,
    ux_next: np.ndarray,
    uz_next: np.ndarray,
    u_t_reference: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    mass: float,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    scale4 = scale[..., None, None]
    metric_m = full_case.geom_m.metric_cov + scale4 * delta_m
    metric_0 = full_case.geom_0.metric_cov + scale4 * delta_0
    metric_p = full_case.geom_p.metric_cov + scale4 * delta_p
    geom, _ = make_time_geometry(metric_0, metric_m, metric_p, dt, dx, dz, with_dgamma=False)
    ein = geom.ricci - 0.5 * metric_0 * geom.r_scalar[..., None, None]
    r_need = ein - case.source_tilde
    rel = tensor_norm(cov0 - r_need) / np.maximum(tensor_norm(r_need), 1.0e-300)
    weights = np.sqrt(np.maximum(case.rho_a[fit_mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))

    current_plus = coordinate_matter_rhs_covector(
        n_cons=n_next,
        u_x=ux_next,
        u_z=uz_next,
        metric_cov_txz=metric_p,
        mass=mass,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=u_t_reference,
    )
    measure_plus = safe_sqrt_abs_det(metric_p) * current_plus["rho"]
    det_plus = np.linalg.det(metric_p.reshape(-1, 3, 3)).reshape(metric_p.shape[:2])
    support = np.asarray(case.support, dtype=bool)
    bad_disc = current_plus["discriminant"] < 0.0
    bad_det = det_plus <= 0.0
    rho_mass = np.maximum(case.rho_a, 0.0)
    support_mass = max(float(np.sum(rho_mass[support])), 1.0e-300)
    delta_norm = _norm_sym(scale4 * delta_p)
    metric_norm = _norm_sym(full_case.geom_p.metric_cov)
    ratio = delta_norm / np.maximum(metric_norm, 1.0e-300)
    row = {
        "name": name,
        "field_residual_weighted_mean": float(weighted_stats(rel[fit_mask], weights)["weighted_mean"]),
        "field_residual_p50": float(weighted_stats(rel[fit_mask], weights)["p50"]),
        "field_residual_p95": float(weighted_stats(rel[fit_mask], weights)["p95"]),
        "one_step_plus_vs_A_relative_l1_support": relative_l1(case_next.measure_tilde, measure_plus, support),
        "negative_discriminant_fraction_support": float(np.count_nonzero(bad_disc[support]) / max(np.count_nonzero(support), 1)),
        "negative_discriminant_rho_mass_fraction_support": float(np.sum(rho_mass[support & bad_disc]) / support_mass),
        "nonpositive_det_fraction_support": float(np.count_nonzero(bad_det[support]) / max(np.count_nonzero(support), 1)),
        "nonpositive_det_rho_mass_fraction_support": float(np.sum(rho_mass[support & bad_det]) / support_mass),
        "bad_disc_core10_fraction": float(np.count_nonzero(bad_disc[case.core10]) / max(np.count_nonzero(case.core10), 1)),
        "bad_det_core10_fraction": float(np.count_nonzero(bad_det[case.core10]) / max(np.count_nonzero(case.core10), 1)),
        "det_plus_min_support": float(np.min(det_plus[support])),
        "disc_plus_min_support": float(np.min(current_plus["discriminant"][support])),
        "relative_delta_plus_p95_support": float(np.percentile(ratio[support], 95.0)),
        "relative_delta_plus_max_support": float(np.max(ratio[support])),
        "scale_mean_support": float(np.mean(scale[support])),
        "scale_p50_support": float(np.percentile(scale[support], 50.0)),
        "scale_p95_support": float(np.percentile(scale[support], 95.0)),
    }
    fields = {
        "field_residual": rel,
        "measure_plus": measure_plus,
        "bad_disc": bad_disc.astype(float),
        "bad_det": bad_det.astype(float),
        "scale": scale,
    }
    return row, fields


def render_scatter(out_path: Path, rows: list[dict[str, float]]) -> None:
    residual = np.asarray([row["field_residual_weighted_mean"] for row in rows], dtype=float)
    measure = np.asarray([row["one_step_plus_vs_A_relative_l1_support"] for row in rows], dtype=float)
    neg = np.asarray([row["negative_discriminant_fraction_support"] for row in rows], dtype=float)
    labels = [row["name"] for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), constrained_layout=True)
    axes[0].scatter(residual, np.maximum(measure, 1.0e-16), c=neg, cmap="magma", s=65)
    axes[0].set_yscale("log")
    axes[0].set_xlabel("D field-equation residual, weighted mean")
    axes[0].set_ylabel("one-step transformed-measure deviation from A")
    axes[0].set_title("Lower-left is better; color = bad mass-shell fraction")
    axes[0].grid(alpha=0.25)
    for x, y, label in zip(residual, np.maximum(measure, 1.0e-16), labels):
        axes[0].annotate(label, (x, y), fontsize=7, alpha=0.75)

    order = np.argsort(residual)
    axes[1].barh(np.arange(len(rows)), residual[order], color="tab:blue", alpha=0.72)
    axes[1].set_yticks(np.arange(len(rows)), [labels[i] for i in order], fontsize=7)
    axes[1].set_xlabel("D field-equation residual, weighted mean")
    axes[1].set_title("Residual ranking")
    axes[1].grid(axis="x", alpha=0.25)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_maps(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho_a: np.ndarray,
    support: np.ndarray,
    fit_mask: np.ndarray,
    name: str,
    fields: dict[str, np.ndarray],
) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 10.0), constrained_layout=True)
    panels = [
        (rho_a, "A rho", "viridis", False),
        (fields["scale"], f"local scale: {name}", "viridis", False),
        (fields["field_residual"], "log10(1 + D residual)", "magma", True),
        (fields["bad_disc"] + 2.0 * fields["bad_det"], "bad cells: 1=disc, 2=det, 3=both", "coolwarm", False),
    ]
    for ax, (data, title, cmap, log_scale) in zip(axes.flat, panels):
        arr = np.asarray(data, dtype=float)
        if log_scale:
            arr = np.log10(1.0 + np.maximum(arr, 0.0))
        im = ax.pcolormesh(xg, zg, arr, shading="auto", cmap=cmap)
        ax.contour(xg, zg, support.astype(float), levels=[0.5], colors=["white"], linewidths=0.6)
        ax.contour(xg, zg, fit_mask.astype(float), levels=[0.5], colors=["cyan"], linewidths=0.6)
        ax.set_title(title + "; white=support, cyan=fit mask")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def parse_tuple(text: str) -> tuple[float, ...]:
    values = tuple(float(part.strip()) for part in text.split(",") if part.strip())
    if not values:
        raise ValueError("empty numeric list")
    return values


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    case = build_case(args, ref, float(args.tau))
    case_next = build_case(args, ref, float(args.tau) + float(args.force_probe_dt_old))
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(case.x[1] - case.x[0])
    dz = float(case.z[1] - case.z[0])
    cov0, coeff_mask = load_lambda_cov({0: case}, args.coefficients, str(args.atom_family))
    data = np.load(args.data)
    fit_mask = np.asarray(data["mask"], dtype=bool)
    if not np.array_equal(fit_mask, coeff_mask):
        raise ValueError("data mask and coefficient mask differ")
    delta_m = optional_array(data, "delta_metric_minus", full_case.geom_m.metric_cov)
    delta_0 = optional_array(data, "delta_metric_center", full_case.geom_0.metric_cov)
    delta_p = optional_array(data, "delta_metric_plus", full_case.geom_p.metric_cov)

    rho_tilde0 = np.maximum(case.measure_tilde / np.maximum(safe_sqrt_abs_det(full_case.geom_0.metric_cov), 1.0e-300), float(args.rho_floor))
    cons0 = conservative_density_from_rho_u(
        rho=np.where(case.support, rho_tilde0, 0.0),
        u_x=case.u_cov[..., 1],
        u_z=case.u_cov[..., 2],
        metric_cov_txz=full_case.geom_0.metric_cov,
        mass=float(ref["params"].m),
        branch="negative_frequency",
        u_t_reference=case.u_cov[..., 0],
        rho_floor=float(args.rho_floor),
    )
    n_next, ux_next, uz_next, current_center = rk4_matter_step(
        n_cons=cons0["n_cons"],
        u_x=case.u_cov[..., 1],
        u_z=case.u_cov[..., 2],
        metric_cov=full_case.geom_0.metric_cov,
        mass=float(ref["params"].m),
        dx=dx,
        dz=dz,
        dt=dt,
        u_t_reference=case.u_cov[..., 0],
        active_mask=case.support,
    )
    scale_maps = build_scale_maps(
        case=case,
        full_metric_plus=full_case.geom_p.metric_cov,
        delta_plus=delta_p,
        fit_mask=fit_mask,
        trust_limits=parse_tuple(args.trust_limits),
        rho_tapers=parse_tuple(args.rho_tapers),
    )
    rows: list[dict[str, float]] = []
    fields_by_name: dict[str, dict[str, np.ndarray]] = {}
    for name, scale in scale_maps.items():
        row, fields = evaluate_candidate(
            name=name,
            scale=scale,
            case=case,
            case_next=case_next,
            full_case=full_case,
            cov0=cov0,
            fit_mask=fit_mask,
            delta_m=delta_m,
            delta_0=delta_0,
            delta_p=delta_p,
            n_next=n_next,
            ux_next=ux_next,
            uz_next=uz_next,
            u_t_reference=current_center["u_t"],
            dt=dt,
            dx=dx,
            dz=dz,
            mass=float(ref["params"].m),
        )
        rows.append(row)
        fields_by_name[name] = fields

    full_fields = fields_by_name.get("full")
    if full_fields is not None:
        bad_full = (full_fields["bad_disc"] > 0.5) | (full_fields["bad_det"] > 0.5)
        auto_maps = {
            "auto_bad_zero": np.where(bad_full, 0.0, 1.0),
            "auto_bad_dilate1_zero": np.where(base.dilate_mask_8(bad_full, 1), 0.0, 1.0),
            "auto_bad_dilate2_zero": np.where(base.dilate_mask_8(bad_full, 2), 0.0, 1.0),
        }
        for name, scale in auto_maps.items():
            row, fields = evaluate_candidate(
                name=name,
                scale=scale,
                case=case,
                case_next=case_next,
                full_case=full_case,
                cov0=cov0,
                fit_mask=fit_mask,
                delta_m=delta_m,
                delta_0=delta_0,
                delta_p=delta_p,
                n_next=n_next,
                ux_next=ux_next,
                uz_next=uz_next,
                u_t_reference=current_center["u_t"],
                dt=dt,
                dx=dx,
                dz=dz,
                mass=float(ref["params"].m),
            )
            rows.append(row)
            fields_by_name[name] = fields

    feasible = [
        row
        for row in rows
        if row["negative_discriminant_fraction_support"] <= float(args.max_negative_discriminant_fraction)
        and row["nonpositive_det_fraction_support"] <= float(args.max_nonpositive_det_fraction)
        and row["one_step_plus_vs_A_relative_l1_support"] <= float(args.max_measure_deviation)
    ]
    best_residual = min(rows, key=lambda row: row["field_residual_weighted_mean"])
    best_feasible = min(feasible, key=lambda row: row["field_residual_weighted_mean"]) if feasible else None
    map_names = ["full", best_residual["name"]]
    if best_feasible is not None and best_feasible["name"] not in map_names:
        map_names.append(best_feasible["name"])
    scatter_path = args.output / "gbcd_plus_local_guard_scan.png"
    render_scatter(scatter_path, rows)
    map_outputs: dict[str, str] = {}
    for name in map_names:
        path = args.output / f"gbcd_plus_local_guard_{name}.png"
        render_maps(path, case.x, case.z, case.rho_a, case.support, fit_mask, name, fields_by_name[name])
        map_outputs[name] = str(path.resolve())

    report = {
        "parameters": {
            "tau": float(args.tau),
            "dt_old": float(args.force_probe_dt_old),
            "data": str(args.data.resolve()),
            "coefficients": str(args.coefficients.resolve()),
            "trust_limits": parse_tuple(args.trust_limits),
            "rho_tapers": parse_tuple(args.rho_tapers),
            "max_measure_deviation": float(args.max_measure_deviation),
            "max_negative_discriminant_fraction": float(args.max_negative_discriminant_fraction),
            "max_nonpositive_det_fraction": float(args.max_nonpositive_det_fraction),
        },
        "meaning": {
            "local_scale": "Pointwise multiplier on the already solved future-slice metric update delta g_+. This does not change the field equation; it probes what kind of trust-region guard the next solver must impose.",
            "field_residual": "Relative mismatch in Gtilde - Ttilde/Mp^2 = C on the initial fitting mask. Smaller is better.",
            "admissibility": "Matter can be reconstructed on candidate g_+ only if the mass-shell discriminant stays nonnegative and det(g_+) keeps the intended Lorentzian branch.",
        },
        "best_residual": best_residual,
        "best_feasible": best_feasible,
        "rows": rows,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "scatter_png": str(scatter_path.resolve()),
            "map_pngs": map_outputs,
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=-3.5)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
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
    parser.add_argument("--atom-family", choices=["auto", "matter4", "matter4_plus_normal", "matter6_normal"], default="auto")
    parser.add_argument("--trust-limits", default="1,3,10,30,100,300,1000")
    parser.add_argument("--rho-tapers", default="0.001,0.003,0.01,0.03,0.1")
    parser.add_argument("--max-measure-deviation", type=float, default=0.05)
    parser.add_argument("--max-negative-discriminant-fraction", type=float, default=0.0)
    parser.add_argument("--max-nonpositive-det-fraction", type=float, default=0.0)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
