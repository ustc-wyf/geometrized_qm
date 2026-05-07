from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coordinate_matter_evolution import inverse_metric_block
from physical_units import (
    PLANCK_LENGTH_EV_INV,
    PLANCK_MASS_EV,
    massive_omega_and_mass_from_ratio,
    photon_k0_ev_from_wavelength_nm,
)
from simulate_d_reduced_dynamic_same_initial import relative_l1
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det, stress_tensor_tilde, tensor_norm


def stats(values: np.ndarray | list[float]) -> dict[str, float]:
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


def weighted_relative_l1(reference: np.ndarray, candidate: np.ndarray, weight: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return 0.0
    denom = max(float(np.sum(weight[valid] * np.abs(reference[valid]))), 1.0e-300)
    return float(np.sum(weight[valid] * np.abs(candidate[valid] - reference[valid])) / denom)


def load_rows(path: Path) -> list[dict[str, float]]:
    if not path.exists():
        return []
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def infer_mass(summary: dict[str, object], explicit_mass: float | None) -> float:
    if explicit_mass is not None:
        return float(explicit_mass)
    params = summary.get("params", {})
    if not isinstance(params, dict):
        raise ValueError("summary params missing; pass --mass explicitly")
    wavelength_nm = float(params["wavelength_nm"])
    mass_over_omega = float(params["mass_over_omega"])
    k0 = photon_k0_ev_from_wavelength_nm(wavelength_nm)
    _, mass = massive_omega_and_mass_from_ratio(k0, mass_over_omega)
    return float(mass)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    summary_path = args.case_dir / "summary.json"
    fields_path = args.case_dir / "fields_final.npz"
    rows_path = args.case_dir / "tensor_interface_rows.jsonl"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    fields = np.load(fields_path)

    params = summary.get("params", {})
    if not isinstance(params, dict):
        params = {}
    mass = infer_mass(summary, args.mass)
    mp = float(args.mp if args.mp is not None else params.get("mp", PLANCK_MASS_EV))
    ell = float(args.ell if args.ell is not None else params.get("ell", args.ell_over_planck * PLANCK_LENGTH_EV_INV))

    metric_cov = np.asarray(fields["metric_cov_D"], dtype=float)
    metric_inv, _ = inverse_metric_block(metric_cov)
    measure = np.asarray(fields["ntilde_measure_D"], dtype=float)
    rho_a = np.asarray(fields["rho_A"], dtype=float)
    rho_pull = np.asarray(fields["rho_D_to_A"], dtype=float)
    support = np.asarray(fields["support"], dtype=bool)
    trusted = np.asarray(fields["trusted"], dtype=bool)
    u_cov = np.stack(
        [
            np.asarray(fields["u_t_D"], dtype=float),
            np.asarray(fields["u_x_D"], dtype=float),
            np.asarray(fields["u_z_D"], dtype=float),
        ],
        axis=-1,
    )
    tilde_x = np.einsum("...ab,...a,...b->...", metric_inv, u_cov, u_cov, optimize=True)
    mass_shell_defect = tilde_x - mass * mass
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)

    rows = load_rows(rows_path)
    row_residuals = np.array([float(row.get("direct_tensor_residual_relative", np.inf)) for row in rows], dtype=float)
    row_solved = np.array([float(row.get("direct_tensor_solved", 0.0)) > 0.5 for row in rows], dtype=bool)
    finite_rows = np.isfinite(row_residuals)
    tensor_pass = bool(
        rows
        and np.all(finite_rows)
        and float(np.count_nonzero(row_solved) / max(len(rows), 1)) >= float(args.require_tensor_solved_fraction)
        and float(np.nanmax(row_residuals)) <= float(args.tensor_rel_tol)
    )

    raw_y = np.asarray(fields["tensor_raw_y"], dtype=float) if "tensor_raw_y" in fields else None
    plateau_summary: dict[str, object] = {"available": raw_y is not None}
    plateau_pass = True
    if raw_y is not None:
        clipped = np.clip(np.abs(raw_y), 0.0, 350.0)
        phi = 1.0 / (np.cosh(clipped) ** 2)
        saturated = support & (phi <= float(args.saturation_phi))
        stress, tilde_x_check = stress_tensor_tilde(
            metric_cov,
            metric_inv,
            rho_tilde,
            u_cov[..., 0],
            u_cov[..., 1],
            u_cov[..., 2],
            m=mass,
        )
        rhs = np.asarray(stress, dtype=float) / (mp * mp)
        f = np.tanh(raw_y) / (ell * ell)
        lhs = -0.5 * f[..., None, None] * metric_cov
        residual = lhs - rhs
        lhs_norm = tensor_norm(lhs)
        rhs_norm = tensor_norm(rhs)
        residual_norm = tensor_norm(residual)
        rel = residual_norm / np.maximum(lhs_norm + rhs_norm, 1.0e-300)
        plateau_summary = {
            "available": True,
            "mask_count": int(np.count_nonzero(saturated)),
            "lhs_norm": stats(lhs_norm[saturated]),
            "rhs_norm": stats(rhs_norm[saturated]),
            "residual_norm": stats(residual_norm[saturated]),
            "relative_residual": stats(rel[saturated]),
            "tilde_x_minus_m2": stats((tilde_x_check - mass * mass)[saturated]),
            "definition": "Saturated-bulk algebraic D equation check: -0.5*f(Rtilde)*gtilde_mn ~= T_mn/Mp^2 where f=tanh(ell^2 Rtilde)/ell^2 and f_R is negligible.",
        }
        plateau_pass = bool(
            np.count_nonzero(saturated) == 0
            or float(plateau_summary["relative_residual"]["p95"]) <= float(args.plateau_rel_p95_tol)
        )

    mass_shell_pass = bool(stats(mass_shell_defect[trusted])["p95"] <= float(args.mass_shell_p95_tol))
    rho_positive_pass = bool(np.nanmin(rho_tilde[trusted]) >= -abs(float(args.rho_tilde_min_tol)))
    pullback_weight = np.where(support, rho_a, 0.0)
    pullback_weighted_l1 = weighted_relative_l1(rho_a, rho_pull, pullback_weight, support)
    pullback_pass = bool(pullback_weighted_l1 <= float(args.pullback_weighted_l1_tol))

    checks = {
        "mass_shell": {
            "pass": mass_shell_pass,
            "tilde_x_minus_m2_stats_trusted": stats(mass_shell_defect[trusted]),
            "tolerance_p95": float(args.mass_shell_p95_tol),
        },
        "positive_tilde_density": {
            "pass": rho_positive_pass,
            "rho_tilde_stats_trusted": stats(rho_tilde[trusted]),
            "minimum_tolerance": -abs(float(args.rho_tilde_min_tol)),
        },
        "pullback_rho_objective": {
            "pass": pullback_pass,
            "weighted_l1_support": pullback_weighted_l1,
            "relative_l1_support": float(relative_l1(rho_a, rho_pull, support)),
            "relative_l1_trusted": float(relative_l1(rho_a, rho_pull, trusted)),
            "tolerance": float(args.pullback_weighted_l1_tol),
            "note": "This is an optimization target, not itself a D field equation.",
        },
        "full_tensor_interface": {
            "pass": tensor_pass,
            "row_count": int(len(rows)),
            "solved_count": int(np.count_nonzero(row_solved)),
            "solved_fraction": float(np.count_nonzero(row_solved) / max(len(rows), 1)),
            "residual_stats": stats(row_residuals[finite_rows]),
            "required_solved_fraction": float(args.require_tensor_solved_fraction),
            "tensor_rel_tol": float(args.tensor_rel_tol),
        },
        "saturated_bulk_algebraic": {
            "pass": plateau_pass,
            **plateau_summary,
            "relative_p95_tolerance": float(args.plateau_rel_p95_tol),
        },
    }
    certified = bool(all(item.get("pass", False) for item in checks.values()))
    result = {
        "case_dir": str(args.case_dir.resolve()),
        "certified": certified,
        "certification_scope": (
            "Residual certificate for the D equations currently implemented in this code path. "
            "It is not a mathematical proof of continuum convergence and does not replace a full implicit D field-equation solver."
        ),
        "parameters": {
            "mass": mass,
            "mp": mp,
            "ell": ell,
            "ell_over_planck": ell / PLANCK_LENGTH_EV_INV,
        },
        "checks": checks,
        "failure_summary": [name for name, item in checks.items() if not item.get("pass", False)],
    }
    output_path = out / "d_branch_equation_certificate.json"
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mass", type=float, default=None)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--mass-shell-p95-tol", type=float, default=1.0e-10)
    parser.add_argument("--rho-tilde-min-tol", type=float, default=1.0e-14)
    parser.add_argument("--pullback-weighted-l1-tol", type=float, default=1.0e-3)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument("--require-tensor-solved-fraction", type=float, default=1.0)
    parser.add_argument("--saturation-phi", type=float, default=1.0e-3)
    parser.add_argument("--plateau-rel-p95-tol", type=float, default=1.0e-2)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
