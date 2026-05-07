from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from test_gbcd_v0_one_step_lambda import weighted_stats


def masked_percentile(field: np.ndarray, mask: np.ndarray, q: float, default: float = 1.0) -> float:
    vals = field[mask & np.isfinite(field)]
    if vals.size == 0:
        return default
    return float(np.percentile(vals, q))


def rho_bin_stats(rho: np.ndarray, residual: np.ndarray, mask: np.ndarray) -> list[dict[str, object]]:
    rho_max = max(float(np.max(rho)), 1.0e-300)
    bins = [
        ("core50", 0.5, 1.0),
        ("core10", 0.1, 1.0),
        ("rho_1pct_to_10pct", 0.01, 0.1),
        ("rho_0p1pct_to_1pct", 0.001, 0.01),
    ]
    out = []
    for name, lo, hi in bins:
        region = mask & (rho >= lo * rho_max) & (rho < hi * rho_max if hi < 1.0 else rho <= hi * rho_max)
        weights = np.sqrt(np.maximum(rho[region], 0.0) / rho_max)
        out.append(
            {
                "name": name,
                "definition": f"mask and {lo:g} <= rho/rho_max <={hi:g}" if hi == 1.0 else f"mask and {lo:g} <= rho/rho_max < {hi:g}",
                "stats": weighted_stats(residual[region], weights) if np.any(region) else {"count": 0},
            }
        )
    return out


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV

    ref = build_physical_reference(args)
    case = build_case(args, ref, float(args.tau))
    full_case = build_full_case(args, ref, float(args.tau))
    data = np.load(args.data)
    mask = np.asarray(data["mask"], dtype=bool)
    exact_after = np.asarray(data["exact_after_rel"], dtype=float)
    linear_after = np.asarray(data["linear_after_rel"], dtype=float)
    before = np.asarray(data["before_rel"], dtype=float)
    corrected_metric_plus = np.asarray(data["corrected_metric_plus"], dtype=float)
    delta_metric_plus = np.asarray(data["delta_metric_plus"], dtype=float)

    rho = np.asarray(case.rho_a, dtype=float)
    det_corrected = np.linalg.det(corrected_metric_plus)
    metric_norm = tensor_norm(full_case.geom_p.metric_cov)
    delta_rel = tensor_norm(delta_metric_plus) / np.maximum(metric_norm, 1.0e-300)

    x_um = case.x * HBAR_C_EV_M * 1.0e6
    z_um = case.z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")

    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.4), constrained_layout=True)
    panels = [
        ("A-branch density rho", rho, "viridis", 0.0, masked_percentile(rho, mask, 99.5)),
        ("D equation residual before metric update", before, "magma", 0.0, masked_percentile(before, mask, 95.0)),
        ("linearized residual after update", linear_after, "magma", 0.0, masked_percentile(linear_after, mask, 95.0)),
        ("exact nonlinear residual after update", exact_after, "magma", 0.0, masked_percentile(exact_after, mask, 95.0)),
        ("relative next-slice metric correction |delta g+|/|g+|", delta_rel, "cividis", 0.0, masked_percentile(delta_rel, mask, 95.0)),
        ("det(corrected g+)", det_corrected, "coolwarm", masked_percentile(det_corrected, mask, 2.0), masked_percentile(det_corrected, mask, 98.0)),
    ]
    for ax, (title, field, cmap, vmin, vmax) in zip(axes.ravel(), panels):
        if vmax <= vmin:
            vmax = vmin + 1.0
        im = ax.pcolormesh(xg, zg, np.clip(field, vmin, vmax), shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour(xg, zg, mask.astype(float), levels=[0.5], colors="white", linewidths=0.8)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(
        "White contour = fitting mask boundary; core10 means rho >= 0.1 rho_max within this mask.",
        fontsize=11,
    )
    plot = args.output / "gbcd_metric_update_diagnostic_maps.png"
    fig.savefig(plot, dpi=180)
    plt.close(fig)

    weights = np.sqrt(np.maximum(rho[mask], 0.0) / max(float(np.max(rho)), 1.0e-300))
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "data": str(args.data.resolve()),
        },
        "definitions": {
            "white_contour": "boundary of the mask used in the full-linear metric update",
            "exact_after_relative_residual": "||C_lambda - (G[corrected_g+] - T_tilde)|| / ||G[corrected_g+] - T_tilde|| pointwise tensor norm",
            "relative_next_slice_metric_correction": "||delta g_+|| / ||g_+|| pointwise tensor norm",
            "core10": "rho >= 0.1 * max(rho) inside the fitting mask",
        },
        "stats": {
            "exact_after_relative_residual_mask": weighted_stats(exact_after[mask], weights),
            "linear_after_relative_residual_mask": weighted_stats(linear_after[mask], weights),
            "before_relative_residual_mask": weighted_stats(before[mask], weights),
            "relative_delta_metric_mask": weighted_stats(delta_rel[mask], weights),
            "det_corrected_metric_plus_mask": {
                "min": float(np.min(det_corrected[mask])),
                "p50": float(np.percentile(det_corrected[mask], 50.0)),
                "p95_abs": float(np.percentile(np.abs(det_corrected[mask]), 95.0)),
                "max": float(np.max(det_corrected[mask])),
                "nonpositive_count": int(np.count_nonzero(det_corrected[mask] <= 0.0)),
            },
            "rho_binned_exact_after": rho_bin_stats(rho, exact_after, mask),
        },
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
    parser.add_argument("--tau", type=float, default=0.0)
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
