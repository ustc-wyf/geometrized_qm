from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from coordinate_matter_evolution import coordinate_matter_rhs_covector
from diagnose_mathcal_r_pure_geometry import build_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_reduced_dynamic_same_initial import rk4_matter_step
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det


def relative_l1(reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(np.abs(candidate[valid] - reference[valid])) / denom)


def weighted_relative_l1(reference: np.ndarray, candidate: np.ndarray, weight: np.ndarray) -> float:
    valid = np.isfinite(reference) & np.isfinite(candidate) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(weight[valid] * np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(weight[valid] * np.abs(candidate[valid] - reference[valid])) / denom)


def stats(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }


def render(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    support: np.ndarray,
    reference_measure: np.ndarray,
    initial_measure: np.ndarray,
    next_measure: np.ndarray,
    discriminant_next: np.ndarray,
) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    vmax = max(
        float(np.percentile(reference_measure[support], 99.5)) if np.any(support) else 0.0,
        float(np.percentile(initial_measure[support], 99.5)) if np.any(support) else 0.0,
        float(np.percentile(next_measure[support], 99.5)) if np.any(support) else 0.0,
        1.0e-300,
    )
    diff = next_measure - reference_measure
    diff_v = max(float(np.percentile(np.abs(diff[support]), 99.0)) if np.any(support) else 0.0, 1.0e-300)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 10.0), constrained_layout=True)
    panels = [
        (initial_measure, "initial D transformed measure", "viridis", 0.0, vmax),
        (next_measure, "one-step D transformed measure", "viridis", 0.0, vmax),
        (diff, "D one-step - A reference measure", "coolwarm", -diff_v, diff_v),
        (discriminant_next, "mass-shell discriminant after one step", "magma", None, None),
    ]
    for ax, (field, title, cmap, vmin, vmax_panel) in zip(axes.flat, panels):
        im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax_panel)
        ax.contour(xg, zg, support.astype(float), levels=[0.5], colors=["white"], linewidths=0.6)
        ax.set_title(title + "; white=support")
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    package = np.load(args.package)
    x = np.asarray(package["x"], dtype=float)
    z = np.asarray(package["z"], dtype=float)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    support = np.asarray(package["support"], dtype=bool)
    trusted = np.asarray(package["trusted"], dtype=bool)
    metric_center = np.asarray(package["metric_center"], dtype=float)
    metric_plus = np.asarray(package["metric_plus"], dtype=float)
    n_cons = np.asarray(package["n_cons"], dtype=float)
    u_x = np.asarray(package["u_x"], dtype=float)
    u_z = np.asarray(package["u_z"], dtype=float)
    u_t = np.asarray(package["u_t"], dtype=float)
    initial_measure = np.asarray(package["measure_tilde"], dtype=float)

    ref = build_physical_reference(args)
    mass = float(ref["params"].m)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.dt_old) * old_scale
    reference_next = build_case(args, ref, float(args.tau) + float(args.dt_old))

    # First microstep: matter sees the initial metric g_0.  Then we reconstruct
    # diagnostic densities both on g_0 and on the already solved next slice g_+.
    n_next, ux_next, uz_next, current_center = rk4_matter_step(
        n_cons=n_cons,
        u_x=u_x,
        u_z=u_z,
        metric_cov=metric_center,
        mass=mass,
        dx=dx,
        dz=dz,
        dt=dt,
        u_t_reference=u_t,
        active_mask=support,
    )
    current_plus = coordinate_matter_rhs_covector(
        n_cons=n_next,
        u_x=ux_next,
        u_z=uz_next,
        metric_cov_txz=metric_plus,
        mass=mass,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=current_center["u_t"],
    )
    measure_next_center = safe_sqrt_abs_det(metric_center) * current_center["rho"]
    measure_next_plus = safe_sqrt_abs_det(metric_plus) * current_plus["rho"]
    weight = np.where(support, np.maximum(reference_next.rho_a, 0.0), 0.0)

    plot_path = args.output / "gbcd_plus_initial_one_step.png"
    render(plot_path, x, z, support, reference_next.measure_tilde, initial_measure, measure_next_plus, current_plus["discriminant"])
    report = {
        "parameters": {
            "package": str(args.package.resolve()),
            "tau_initial": float(args.tau),
            "dt_old": float(args.dt_old),
            "dt_ev_inv": float(dt),
            "mass_ev": mass,
        },
        "meaning": {
            "one_step": "Advance only the D matter variables for one tiny time step from the plus-only initial package.",
            "metric_center_measure": "Diagnostic density reconstructed on the initial metric g_0.",
            "metric_plus_measure": "Diagnostic density reconstructed on the solved next metric g_+.",
            "reference": "A/KG reference at tau_initial + dt_old.  It is only a comparison target, not imposed during the step.",
        },
        "diagnostics": {
            "support_points": int(np.count_nonzero(support)),
            "trusted_points": int(np.count_nonzero(trusted)),
            "initial_vs_A_next_relative_l1_support": relative_l1(reference_next.measure_tilde, initial_measure, support),
            "one_step_center_vs_A_next_relative_l1_support": relative_l1(reference_next.measure_tilde, measure_next_center, support),
            "one_step_plus_vs_A_next_relative_l1_support": relative_l1(reference_next.measure_tilde, measure_next_plus, support),
            "one_step_plus_vs_A_next_weighted_l1_support": weighted_relative_l1(reference_next.measure_tilde, measure_next_plus, weight),
            "disc_center_after_step": stats(current_center["discriminant"][support]),
            "disc_plus_after_step": stats(current_plus["discriminant"][support]),
            "negative_discriminant_fraction_center": float(np.count_nonzero(current_center["discriminant"][support] < 0.0) / max(np.count_nonzero(support), 1)),
            "negative_discriminant_fraction_plus": float(np.count_nonzero(current_plus["discriminant"][support] < 0.0) / max(np.count_nonzero(support), 1)),
            "n_cons_relative_change_support": relative_l1(n_cons, n_next, support),
            "u_x_relative_change_support": relative_l1(u_x, ux_next, support),
            "u_z_relative_change_support": relative_l1(u_z, uz_next, support),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=-3.5)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
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
