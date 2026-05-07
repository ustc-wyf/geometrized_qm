from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from coordinate_matter_evolution import inverse_metric_block
from physical_units import (
    HBAR_C_EV_M,
    PLANCK_LENGTH_EV_INV,
    PLANCK_MASS_EV,
    massive_omega_and_mass_from_ratio,
    photon_k0_ev_from_wavelength_nm,
)
from simulate_d_reduced_dynamic_same_initial import relative_l1
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det


def stats(values: np.ndarray | list[float]) -> dict[str, float]:
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


def abs_stats(values: np.ndarray | list[float]) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = np.abs(vals[np.isfinite(vals)])
    if vals.size == 0:
        return {"count": 0, "min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    return {
        "count": int(vals.size),
        "min": float(np.min(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
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


def matrix_rank_from_svd(singular_values: np.ndarray, rel_tol: float) -> np.ndarray:
    leading = singular_values[..., :1]
    threshold = np.maximum(leading * rel_tol, 1.0e-300)
    return np.sum(singular_values > threshold, axis=-1)


def lorentz_signature_counts(metric_cov: np.ndarray, mask: np.ndarray) -> dict[str, int]:
    if not np.any(mask):
        return {"lorentz_1pos_2neg": 0, "lorentz_2pos_1neg": 0, "degenerate_or_other": 0}
    mats = metric_cov[mask].reshape(-1, 3, 3)
    eig = np.linalg.eigvalsh(mats)
    pos = np.sum(eig > 1.0e-12 * np.maximum(np.max(np.abs(eig), axis=-1, keepdims=True), 1.0), axis=-1)
    neg = np.sum(eig < -1.0e-12 * np.maximum(np.max(np.abs(eig), axis=-1, keepdims=True), 1.0), axis=-1)
    return {
        "lorentz_1pos_2neg": int(np.count_nonzero((pos == 1) & (neg == 2))),
        "lorentz_2pos_1neg": int(np.count_nonzero((pos == 2) & (neg == 1))),
        "degenerate_or_other": int(np.count_nonzero(~(((pos == 1) & (neg == 2)) | ((pos == 2) & (neg == 1))))),
    }


def row_summary(rows: list[dict[str, float]], tensor_rel_tol: float) -> dict[str, object]:
    if not rows:
        return {
            "available": False,
            "row_count": 0,
            "solved_count": 0,
            "solved_fraction": 0.0,
            "residual_stats": abs_stats([]),
        }
    residuals = np.array([float(row.get("direct_tensor_residual_relative", np.nan)) for row in rows], dtype=float)
    solved = np.array(
        [
            (float(row.get("direct_tensor_solved", 0.0)) > 0.5)
            and np.isfinite(float(row.get("direct_tensor_residual_relative", np.nan)))
            and abs(float(row.get("direct_tensor_residual_relative", np.nan))) <= tensor_rel_tol
            for row in rows
        ],
        dtype=bool,
    )
    return {
        "available": True,
        "row_count": int(len(rows)),
        "solved_count": int(np.count_nonzero(solved)),
        "solved_fraction": float(np.count_nonzero(solved) / max(len(rows), 1)),
        "residual_stats": abs_stats(residuals),
        "tensor_rel_tol": float(tensor_rel_tol),
    }


def render_infeasibility_plot(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho_a: np.ndarray,
    raw_y: np.ndarray,
    saturated: np.ndarray,
    matter_saturated: np.ndarray,
    infeasible: np.ndarray,
) -> None:
    extent = [
        float(x[0] * HBAR_C_EV_M * 1.0e6),
        float(x[-1] * HBAR_C_EV_M * 1.0e6),
        float(z[0] * HBAR_C_EV_M * 1.0e6),
        float(z[-1] * HBAR_C_EV_M * 1.0e6),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4), constrained_layout=True)
    im0 = axes[0].imshow(rho_a.T, origin="lower", extent=extent, aspect="equal", cmap="magma")
    axes[0].set_title("A branch rho")
    fig.colorbar(im0, ax=axes[0], fraction=0.046)

    im1 = axes[1].imshow(np.sign(raw_y).T * np.log10(1.0 + np.abs(raw_y)).T, origin="lower", extent=extent, aspect="equal", cmap="coolwarm")
    axes[1].set_title("sign(raw_y) log10(1+|raw_y|)")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)

    axes[2].imshow(saturated.T.astype(float), origin="lower", extent=extent, aspect="equal", cmap="Greys", vmin=0.0, vmax=1.0)
    axes[2].contour(x * HBAR_C_EV_M * 1.0e6, z * HBAR_C_EV_M * 1.0e6, matter_saturated.T.astype(float), levels=[0.5], colors=["tab:orange"], linewidths=0.8)
    axes[2].set_title("saturated bulk; orange=matter support")

    axes[3].imshow(infeasible.T.astype(float), origin="lower", extent=extent, aspect="equal", cmap="Reds", vmin=0.0, vmax=1.0)
    axes[3].contour(x * HBAR_C_EV_M * 1.0e6, z * HBAR_C_EV_M * 1.0e6, rho_a.T, levels=[float(np.nanmax(rho_a) * 1.0e-3)], colors=["black"], linewidths=0.6)
    axes[3].set_title("hard-constraint infeasible cells")

    for ax in axes:
        ax.set_xlabel("x (um)")
        ax.set_ylabel("z (um)")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


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

    x = np.asarray(fields["x"], dtype=float)
    z = np.asarray(fields["z"], dtype=float)
    rho_a = np.asarray(fields["rho_A"], dtype=float)
    rho_pull = np.asarray(fields["rho_D_to_A"], dtype=float)
    measure = np.asarray(fields["ntilde_measure_D"], dtype=float)
    metric_cov = np.asarray(fields["metric_cov_D"], dtype=float)
    support = np.asarray(fields["support"], dtype=bool)
    trusted = np.asarray(fields["trusted"], dtype=bool)
    raw_y = np.asarray(fields["tensor_raw_y"], dtype=float) if "tensor_raw_y" in fields else np.zeros_like(rho_a)
    u_cov = np.stack(
        [
            np.asarray(fields["u_t_D"], dtype=float),
            np.asarray(fields["u_x_D"], dtype=float),
            np.asarray(fields["u_z_D"], dtype=float),
        ],
        axis=-1,
    )

    metric_inv, _ = inverse_metric_block(metric_cov)
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)
    tilde_x = np.einsum("...ab,...a,...b->...", metric_inv, u_cov, u_cov, optimize=True)
    mass_shell_defect = tilde_x - mass * mass

    clipped = np.clip(np.abs(raw_y), 0.0, 350.0)
    f_r = 1.0 / (np.cosh(clipped) ** 2)
    f = np.tanh(raw_y) / (ell * ell)
    saturated = support & (f_r <= float(args.saturation_phi))
    rho_threshold = float(np.nanmax(rho_a[support])) * float(args.matter_rho_frac) if np.any(support) else 0.0
    matter_support = support & (rho_a >= rho_threshold)
    matter_saturated = saturated & matter_support
    mass_shell_ok = np.abs(mass_shell_defect) <= float(args.mass_shell_abs_tol)
    algebraic_f_nonzero = np.abs(f) >= float(args.f_abs_floor)
    rank_obstruction_mask = matter_saturated & mass_shell_ok & algebraic_f_nonzero

    # In the saturated bulk, the exact metric f(R) equation drops to
    #   -0.5 M_P^2 f g_mn = T_mn.
    # The matter equation enforces g^mn u_m u_n=m^2, so T_mn is proportional
    # to u_m u_n and has rank one wherever rho_tilde is nonzero.
    stress_on_shell = 2.0 * rho_tilde[..., None, None] * np.einsum("...a,...b->...ab", u_cov, u_cov, optimize=True)
    stress_svd = np.linalg.svd(stress_on_shell.reshape(-1, 3, 3), compute_uv=False).reshape(rho_a.shape + (3,))
    stress_rank = matrix_rank_from_svd(stress_svd, float(args.rank_rel_tol))

    with np.errstate(divide="ignore", invalid="ignore"):
        required_metric = -2.0 * stress_on_shell / np.maximum((mp * mp) * np.abs(f)[..., None, None], 1.0e-300)
        required_metric *= np.sign(f)[..., None, None]
    required_det = np.linalg.det(required_metric.reshape(-1, 3, 3)).reshape(rho_a.shape)
    required_svd = np.linalg.svd(required_metric.reshape(-1, 3, 3), compute_uv=False).reshape(rho_a.shape + (3,))
    required_rank = matrix_rank_from_svd(required_svd, float(args.rank_rel_tol))
    degenerate_required_metric = rank_obstruction_mask & (required_rank < 3)
    stress_rank_obstruction = rank_obstruction_mask & (stress_rank < 3)
    infeasible = rank_obstruction_mask & (degenerate_required_metric | stress_rank_obstruction)

    rows = load_rows(rows_path)
    tensor_rows = row_summary(rows, float(args.tensor_rel_tol))
    pullback_weight = np.where(support, rho_a, 0.0)
    pullback_weighted_l1 = weighted_relative_l1(rho_a, rho_pull, pullback_weight, support)
    strict_plateau_blocked = bool(np.count_nonzero(infeasible) > 0)
    tensor_rows_all_pass = bool(tensor_rows.get("solved_fraction", 0.0) == 1.0)
    candidate_can_evolve = bool(args.strict_plateau_algebraic_hard and (not strict_plateau_blocked) and tensor_rows_all_pass)
    if candidate_can_evolve:
        status = "solved_under_strict_plateau_gate"
    elif args.strict_plateau_algebraic_hard and strict_plateau_blocked:
        status = "strict_plateau_algebraic_infeasible"
    else:
        status = "requires_full_fR_equation_solve"

    plot_path = out / "d_initial_solver_constraint_map.png"
    render_infeasibility_plot(plot_path, x, z, rho_a, raw_y, saturated, matter_saturated, infeasible)

    report = {
        "case_dir": str(args.case_dir.resolve()),
        "status": status,
        "solver_scope": (
            "Hard-constrained D initial-data gate for the three-domain metric f(R) solver developed in this project. "
            "It enforces the saturated-bulk algebraic tensor equation as a necessary exact D-equation condition, "
            "then checks whether the existing full tensor interface rows can be accepted everywhere."
        ),
        "parameters": {
            "mass": mass,
            "mp": mp,
            "ell": ell,
            "ell_over_planck_length": ell / PLANCK_LENGTH_EV_INV,
            "saturation_phi": float(args.saturation_phi),
            "matter_rho_frac": float(args.matter_rho_frac),
            "rho_threshold": rho_threshold,
            "mass_shell_abs_tol": float(args.mass_shell_abs_tol),
            "rank_rel_tol": float(args.rank_rel_tol),
        },
        "pullback_objective_before_projection": {
            "weighted_relative_l1_support": pullback_weighted_l1,
            "relative_l1_support": float(relative_l1(rho_a, rho_pull, support)),
            "relative_l1_trusted": float(relative_l1(rho_a, rho_pull, trusted)),
            "meaning": "This is how close the A-derived candidate already is in the rho objective before enforcing hard D equations.",
        },
        "mass_shell": {
            "trusted_abs_stats": abs_stats(mass_shell_defect[trusted]),
            "support_abs_stats": abs_stats(mass_shell_defect[support]),
            "ok_count_in_matter_saturated": int(np.count_nonzero(rank_obstruction_mask)),
            "matter_saturated_count": int(np.count_nonzero(matter_saturated)),
        },
        "saturated_bulk_rank_condition": {
            "equation": "-0.5*M_P^2*f*gtilde_mn = Ttilde_mn in saturated bulk where f_R and its derivatives are negligible",
            "on_shell_source": "Ttilde_mn = 2*rho_tilde*u_m*u_n when gtilde^mn u_m u_n=m^2",
            "necessary_condition": "For nonzero rho_tilde and nonzero f, the RHS is rank one while a Lorentzian metric must be rank three in this 2+1d simulation.",
            "support_count": int(np.count_nonzero(support)),
            "trusted_count": int(np.count_nonzero(trusted)),
            "saturated_count": int(np.count_nonzero(saturated)),
            "matter_support_count": int(np.count_nonzero(matter_support)),
            "matter_saturated_count": int(np.count_nonzero(matter_saturated)),
            "rank_obstruction_candidate_count": int(np.count_nonzero(rank_obstruction_mask)),
            "infeasible_count": int(np.count_nonzero(infeasible)),
            "infeasible_fraction_of_support": float(np.count_nonzero(infeasible) / max(np.count_nonzero(support), 1)),
            "infeasible_rho_mass_fraction": float(np.sum(rho_a[infeasible]) / max(np.sum(rho_a[support]), 1.0e-300)),
            "stress_on_shell_rank_stats": stats(stress_rank[rank_obstruction_mask]),
            "required_metric_rank_stats": stats(required_rank[rank_obstruction_mask]),
            "required_metric_det_abs_stats": abs_stats(required_det[rank_obstruction_mask]),
            "current_metric_signature_on_support": lorentz_signature_counts(metric_cov, support),
            "verdict": "strict_plateau_gate_fails" if np.count_nonzero(infeasible) else "strict_plateau_gate_passes",
            "caveat": (
                "This is only the strict smooth-plateau algebraic limit. It must not be used as a proof that the full D equation is infeasible, "
                "because the full metric f(R) equation also contains phi*Ricci and (g Box - nabla nabla) phi terms."
            ),
        },
        "full_tensor_interface_rows": tensor_rows,
        "evolution_test": {
            "run": candidate_can_evolve,
            "reason_not_run": None
            if candidate_can_evolve
            else (
                "This script no longer treats the strict plateau algebraic rank test as a full-D infeasibility proof. "
                "A full f(R) equation residual/implicit solve is required before evolving a certified D candidate."
            ),
        },
        "outputs": {
            "report_json": str((out / "d_initial_solver_report.json").resolve()),
            "constraint_map_png": str(plot_path.resolve()),
        },
        "interpretation": (
            "The strict smooth-plateau algebraic reduction fails on this A-like optical snapshot. "
            "This does not prove the full D branch impossible; it proves that dropping the phi*Ricci and derivative phi terms is too strong here."
        )
        if np.count_nonzero(infeasible)
        else (
            "No strict plateau rank obstruction was detected. The next gate is still the full f(R) equation residual/implicit solve."
        ),
    }

    report_path = out / "d_initial_solver_report.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mass", type=float, default=None)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--saturation-phi", type=float, default=1.0e-3)
    parser.add_argument("--matter-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--mass-shell-abs-tol", type=float, default=1.0e-10)
    parser.add_argument("--rank-rel-tol", type=float, default=1.0e-8)
    parser.add_argument("--f-abs-floor", type=float, default=1.0e-300)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument(
        "--strict-plateau-algebraic-hard",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "If true, treat the smooth saturated-bulk algebraic limit as a hard rejection gate. "
            "Default false because the full D equation includes phi*Ricci and derivative phi terms."
        ),
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
