from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from analyze_c_terms_from_a_reference import metric_jets_full
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from fit_metric_fr_ricci2_universal_full import make_time_geometry
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from solve_gbcd_full_linear_metric_update_sparse import load_lambda_cov
from test_gbcd_v0_one_step_lambda import ATOM_NAMES, weighted_stats


COMPONENT_LABELS = ["tt", "tx", "tz", "xx", "xz", "zz"]


def erode8(mask: np.ndarray) -> np.ndarray:
    out = mask.copy()
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            shifted = np.zeros_like(mask, dtype=bool)
            src_i0 = max(0, -di)
            src_i1 = mask.shape[0] - max(0, di)
            src_j0 = max(0, -dj)
            src_j1 = mask.shape[1] - max(0, dj)
            dst_i0 = max(0, di)
            dst_i1 = mask.shape[0] - max(0, -di)
            dst_j0 = max(0, dj)
            dst_j1 = mask.shape[1] - max(0, -dj)
            shifted[dst_i0:dst_i1, dst_j0:dst_j1] = mask[src_i0:src_i1, src_j0:src_j1]
            out &= shifted
    return out


def boundary_layers(mask: np.ndarray) -> np.ndarray:
    dist = -np.ones(mask.shape, dtype=int)
    remaining = mask.copy()
    layer = 0
    while np.any(remaining):
        eroded = erode8(remaining)
        shell = remaining & ~eroded
        dist[shell] = layer
        remaining = eroded
        layer += 1
    return dist


def finite_stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float | int]:
    vals = np.asarray(values, dtype=float)
    good = np.isfinite(vals)
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        good &= np.isfinite(w) & (w >= 0.0)
        w = w[good]
    else:
        w = None
    vals = vals[good]
    if vals.size == 0:
        return {"count": 0}
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }
    if w is not None:
        out["weighted_mean"] = float(np.sum(vals * w) / max(float(np.sum(w)), 1.0e-300))
    return out


def region_report(name: str, region: np.ndarray, residual: np.ndarray, weights: np.ndarray) -> dict[str, object]:
    if not np.any(region):
        return {"name": name, "count": 0}
    return {"name": name, "stats": weighted_stats(residual[region], weights[region])}


def quantile_regions(field: np.ndarray, mask: np.ndarray, labels: list[str]) -> list[tuple[str, np.ndarray]]:
    vals = field[mask & np.isfinite(field)]
    if vals.size == 0:
        return []
    cuts = np.percentile(vals, [25.0, 50.0, 75.0])
    regions = []
    lo = -np.inf
    for label, hi in zip(labels, [cuts[0], cuts[1], cuts[2], np.inf]):
        regions.append((label, mask & np.isfinite(field) & (field >= lo) & (field < hi)))
        lo = hi
    return regions


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    base.ATOM_NAMES = ATOM_NAMES

    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    cases = {key: build_case(args, ref, float(args.tau) + key * float(args.force_probe_dt_old)) for key in base.TIME_KEYS}
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])

    cov0, coeff_mask = load_lambda_cov(cases, args.coefficients)
    data = np.load(args.data)
    corrected_metric_minus = np.asarray(
        data["corrected_metric_minus"] if "corrected_metric_minus" in data.files else full_case.geom_m.metric_cov,
        dtype=float,
    )
    corrected_metric_center = np.asarray(
        data["corrected_metric_center"] if "corrected_metric_center" in data.files else full_case.geom_0.metric_cov,
        dtype=float,
    )
    corrected_metric_plus = np.asarray(data["corrected_metric_plus"], dtype=float)
    delta_metric = np.asarray(data["delta_metric_plus"], dtype=float)
    mask = np.asarray(data["mask"], dtype=bool)
    if not np.array_equal(mask, coeff_mask):
        raise ValueError("data mask and coefficient mask differ")

    geom_corr, _ = make_time_geometry(corrected_metric_center, corrected_metric_minus, corrected_metric_plus, dt, dx, dz, with_dgamma=False)
    ein_corr = geom_corr.ricci - 0.5 * corrected_metric_center * geom_corr.r_scalar[..., None, None]
    exact_r_need = ein_corr - cases[0].source_tilde
    err = cov0 - exact_r_need
    denom = np.maximum(tensor_norm(exact_r_need), 1.0e-300)
    rel = tensor_norm(err) / denom
    rho = np.asarray(cases[0].rho_a, dtype=float)
    rho_max = max(float(np.max(rho)), 1.0e-300)
    weights_field = np.sqrt(np.maximum(rho, 0.0) / rho_max)
    metric_norm = tensor_norm(full_case.geom_p.metric_cov)
    delta_rel = tensor_norm(delta_metric) / np.maximum(metric_norm, 1.0e-300)
    det = np.linalg.det(corrected_metric_plus)
    dist = boundary_layers(mask)

    component_rel = np.zeros(mask.shape + (len(base.SYMMETRIC_COMPONENTS),), dtype=float)
    for cpos, (a, b) in enumerate(base.SYMMETRIC_COMPONENTS):
        component_rel[..., cpos] = np.abs(err[..., a, b]) / denom
    dominant_component = np.argmax(component_rel, axis=-1)

    regions: list[dict[str, object]] = []
    regions.append(region_report("core50: rho >= 0.5 rho_max", mask & (rho >= 0.5 * rho_max), rel, weights_field))
    regions.append(region_report("core10_only: 0.1 rho_max <= rho < 0.5 rho_max", mask & (rho >= 0.1 * rho_max) & (rho < 0.5 * rho_max), rel, weights_field))
    for label, region in [
        ("boundary_layer_0", mask & (dist == 0)),
        ("boundary_layer_1", mask & (dist == 1)),
        ("boundary_layer_2_to_3", mask & (dist >= 2) & (dist <= 3)),
        ("boundary_layer_ge4", mask & (dist >= 4)),
    ]:
        regions.append(region_report(label, region, rel, weights_field))
    for label, region in quantile_regions(np.abs(det), mask, ["det_abs_q0_25", "det_abs_q25_50", "det_abs_q50_75", "det_abs_q75_100"]):
        regions.append(region_report(label, region, rel, weights_field))
    for label, region in quantile_regions(delta_rel, mask, ["delta_rel_q0_25", "delta_rel_q25_50", "delta_rel_q50_75", "delta_rel_q75_100"]):
        regions.append(region_report(label, region, rel, weights_field))
    regions.append(region_report("det_nonpositive", mask & (det <= 0.0), rel, weights_field))

    comp_stats = []
    for cpos, label in enumerate(COMPONENT_LABELS):
        comp_stats.append(
            {
                "component": label,
                "relative_abs_component_over_tensor_denominator": weighted_stats(component_rel[..., cpos][mask], weights_field[mask]),
                "dominant_count": int(np.count_nonzero(mask & (dominant_component == cpos))),
            }
        )

    x_um = cases[0].x * HBAR_C_EV_M * 1.0e6
    z_um = cases[0].z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.3), constrained_layout=True)
    finite_rel = rel[mask & np.isfinite(rel)]
    vmax_rel = float(np.percentile(finite_rel, 95.0)) if finite_rel.size else 1.0
    panels = [
        ("rho", rho, "viridis", 0.0, float(np.percentile(rho[mask], 99.5))),
        ("exact residual", rel, "magma", 0.0, vmax_rel),
        ("boundary layer index", np.where(mask, dist, np.nan), "plasma", 0.0, max(float(np.nanmax(np.where(mask, dist, np.nan))), 1.0)),
        ("log10 |det(corrected g+)|", np.log10(np.maximum(np.abs(det), 1.0e-300)), "cividis", None, None),
        ("log10 |delta g+|/|g+|", np.log10(np.maximum(delta_rel, 1.0e-300)), "cividis", None, None),
        ("dominant residual component", np.where(mask, dominant_component, np.nan), "tab10", 0.0, len(COMPONENT_LABELS) - 1),
    ]
    for ax, (title, field, cmap, vmin, vmax) in zip(axes.ravel(), panels):
        arr = np.asarray(field, dtype=float)
        if vmin is None:
            vals = arr[mask & np.isfinite(arr)]
            vmin = float(np.percentile(vals, 2.0)) if vals.size else 0.0
        if vmax is None:
            vals = arr[mask & np.isfinite(arr)]
            vmax = float(np.percentile(vals, 98.0)) if vals.size else 1.0
        if vmax <= vmin:
            vmax = vmin + 1.0
        im = ax.pcolormesh(xg, zg, np.clip(arr, vmin, vmax), shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour(xg, zg, mask.astype(float), levels=[0.5], colors="white", linewidths=0.8)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle("White contour = core10 fitting mask boundary; boundary layer 0 is the mask edge.", fontsize=11)
    plot = args.output / "gbcd_residual_source_diagnostics.png"
    fig.savefig(plot, dpi=180)
    plt.close(fig)

    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "data": str(args.data.resolve()),
            "coefficients": str(args.coefficients.resolve()),
        },
        "definitions": {
            "boundary_layer_0": "points in the fitting mask with at least one 8-neighbor outside the mask",
            "component_rel": "|tensor residual component| / ||exact target tensor||",
            "dominant_component": "component with largest component_rel at that point",
        },
        "overall": {
            "exact_relative_residual_mask": weighted_stats(rel[mask], weights_field[mask]),
            "delta_rel_mask": weighted_stats(delta_rel[mask], weights_field[mask]),
            "det_nonpositive_count": int(np.count_nonzero(mask & (det <= 0.0))),
            "mask_count": int(np.count_nonzero(mask)),
        },
        "region_breakdown": regions,
        "component_breakdown": comp_stats,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot.resolve()),
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
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="core10")
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
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
