from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_equation_first_transverse_force import (
    covector_norm_euclidean,
    divergence_cov2,
    derivative_tensor,
    transverse_covector,
)
from diagnose_gbcd_trace_local_closure import (
    ATOM_NAMES,
    qe_trace_coefficients,
    solve_local,
    trace_coefficients,
    weighted_stats,
)
from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


TIME_KEYS = (-1, 0, 1)
SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def tensor_field_from_coeff(case, coeff: np.ndarray) -> np.ndarray:
    atoms = tensor_atoms(case, "matter4")
    out = np.zeros_like(case.r_need)
    for apos, name in enumerate(ATOM_NAMES):
        out += coeff[..., apos, None, None] * atoms[name]
    return out


def closure_rows(case, closure: str) -> np.ndarray | None:
    if closure == "unconstrained":
        return None
    if closure == "trace0":
        return trace_coefficients(case)
    if closure == "qe0":
        return qe_trace_coefficients(case)
    raise ValueError(f"unknown closure: {closure}")


def make_active_mask(cases: dict[int, object], region: str, force_region: str, active_dilation: int) -> np.ndarray:
    base.ATOM_NAMES = ATOM_NAMES
    return base.make_active_mask(cases, region, force_region, active_dilation)


def local_coefficients(case, closure: str, active_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    atoms = tensor_atoms(case, "matter4")
    coeff = np.zeros(case.rho_a.shape + (len(ATOM_NAMES),), dtype=float)
    rel_res = np.full(case.rho_a.shape, np.nan, dtype=float)
    c_rows = closure_rows(case, closure)
    points = np.argwhere(active_mask)
    target_norm = tensor_norm(case.r_need)
    finite_target = target_norm[active_mask & np.isfinite(target_norm)]
    floor = 0.05 * max(float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0, 1.0e-300)
    for i_raw, j_raw in points:
        i = int(i_raw)
        j = int(j_raw)
        a = np.zeros((len(SYMMETRIC_COMPONENTS), len(ATOM_NAMES)), dtype=float)
        for apos, atom_name in enumerate(ATOM_NAMES):
            for cpos, (mu, nu) in enumerate(SYMMETRIC_COMPONENTS):
                a[cpos, apos] = atoms[atom_name][i, j, mu, nu]
        b = np.asarray([case.r_need[i, j, mu, nu] for mu, nu in SYMMETRIC_COMPONENTS], dtype=float)
        constraint = None if c_rows is None else c_rows[i, j]
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        c, res = solve_local(a, b, constraint)
        if not np.all(np.isfinite(c)) or not np.isfinite(res):
            continue
        coeff[i, j, :] = c
        rel_res[i, j] = res / max(float(target_norm[i, j]), floor, 1.0e-300)
    return coeff, rel_res


def evaluate_center(
    args: argparse.Namespace,
    ref: dict[str, object],
    tau: float,
    closure: str,
) -> dict[str, object]:
    probe = float(args.force_probe_dt_old)
    cases = {
        -1: build_case(args, ref, tau - probe),
        0: build_case(args, ref, tau),
        1: build_case(args, ref, tau + probe),
    }
    full = build_full_case(args, ref, tau)
    scale = ref["scale"]
    dt = probe * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    active = make_active_mask(cases, args.fit_region, args.force_region, int(args.active_dilation))
    force_mask = base.make_force_mask(cases[0], args.force_region, int(args.force_mask_erosion))

    coeffs: dict[int, np.ndarray] = {}
    rels: dict[int, np.ndarray] = {}
    cov: dict[int, np.ndarray] = {}
    for key in TIME_KEYS:
        coeffs[key], rels[key] = local_coefficients(cases[key], closure, active)
        cov[key] = tensor_field_from_coeff(cases[key], coeffs[key])

    dcov = derivative_tensor(cov[-1], cov[0], cov[1], dt, dx, dz)
    j_cov = divergence_cov2(cov[0], dcov, full.geom_0.metric_inv, full.geom_0.gamma)
    j_perp = transverse_covector(j_cov, cases[0].u_cov, full.geom_0.metric_inv)
    total = covector_norm_euclidean(j_cov)
    perp = covector_norm_euclidean(j_perp)
    c_norm = tensor_norm(cov[0])
    h = min(abs(dt), abs(dx), abs(dz))
    weights_force = np.sqrt(
        np.maximum(cases[0].rho_a[force_mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300)
    )

    by_time = {}
    for key in TIME_KEYS:
        mask = getattr(cases[key], args.fit_region)
        weights = np.sqrt(np.maximum(cases[key].rho_a[mask], 0.0) / max(float(np.max(cases[key].rho_a)), 1.0e-300))
        by_time[str(key)] = {
            "tau": float(cases[key].tau_old),
            "algebraic_relative_residual": weighted_stats(rels[key][mask], weights),
        }

    return {
        "tau": float(tau),
        "closure": closure,
        "fit_region": args.fit_region,
        "force_region": args.force_region,
        "active_count": int(np.count_nonzero(active)),
        "force_count": int(np.count_nonzero(force_mask)),
        "by_time": by_time,
        "central_divergence": {
            "full_divergence_norm": weighted_stats(total[force_mask], weights_force),
            "transverse_divergence_norm": weighted_stats(perp[force_mask], weights_force),
            "full_over_derivative_scale": weighted_stats(
                total[force_mask] / np.maximum(c_norm[force_mask] / max(h, 1.0e-300), 1.0e-300),
                weights_force,
            ),
            "transverse_over_derivative_scale": weighted_stats(
                perp[force_mask] / np.maximum(c_norm[force_mask] / max(h, 1.0e-300), 1.0e-300),
                weights_force,
            ),
            "transverse_fraction": weighted_stats(perp[force_mask] / np.maximum(total[force_mask], 1.0e-300), weights_force),
        },
    }


def render(path: Path, results: list[dict[str, object]]) -> None:
    closures = sorted({str(r["closure"]) for r in results})
    taus = sorted({float(r["tau"]) for r in results})
    x = np.arange(len(taus))
    width = 0.8 / max(len(closures), 1)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    for cpos, closure in enumerate(closures):
        selected = [r for tau in taus for r in results if float(r["tau"]) == tau and str(r["closure"]) == closure]
        alg = [
            float(r["by_time"]["0"]["algebraic_relative_residual"]["weighted_mean"])  # type: ignore[index]
            for r in selected
        ]
        div = [
            float(r["central_divergence"]["full_over_derivative_scale"]["weighted_mean"])  # type: ignore[index]
            for r in selected
        ]
        offset = (cpos - (len(closures) - 1) / 2.0) * width
        axes[0].bar(x + offset, alg, width=width, label=closure)
        axes[1].bar(x + offset, div, width=width, label=closure)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels([f"tau={tau:g}" for tau in taus])
        ax.set_yscale("log")
        ax.legend(fontsize=8)
    axes[0].set_title("local algebraic residual")
    axes[0].set_ylabel("weighted mean")
    axes[1].set_title("natural divergence")
    axes[1].set_ylabel(r"weighted mean |nabla C| / (|C|/h)")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ATOM_NAMES
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    taus = [float(x.strip()) for x in args.taus.split(",") if x.strip()]
    closures = [x.strip() for x in args.closures.split(",") if x.strip()]
    results = []
    for tau in taus:
        for closure in closures:
            print(f"[local-conservation] tau={tau:g} closure={closure}", flush=True)
            results.append(evaluate_center(args, ref, tau, closure))
    plot_path = args.output / "trace_local_conservation_summary.png"
    render(plot_path, results)
    report = {
        "parameters": {
            "taus": taus,
            "closures": closures,
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "force_probe_dt_old": float(args.force_probe_dt_old),
        },
        "definition": {
            "purpose": "Check whether pointwise local closure fields naturally satisfy full conservation. This is not a hard-conservation solve.",
            "algebraic_residual": "Pointwise local fit residual of C_mn to R_need in span{g,uu,rr,ur}.",
            "natural_divergence": "|tilde nabla^mu C_mu nu| / (|C|/h), evaluated at the central time slice from tau±force_probe_dt_old.",
            "interpretation": "Small natural_divergence means local closure is already compatible with conservation; large natural_divergence means a genuine global conservation solve is required.",
        },
        "results": results,
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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--closures", type=str, default="unconstrained,trace0,qe0")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
    parser.add_argument("--active-dilation", type=int, default=1)
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
