from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_gbcd_conservation_principal_symbol import principal_matrix, stats as stats_plain
from diagnose_mathcal_r_pure_geometry import tensor_norm
from fit_gbcd_auxiliary_gauge import solve_hard_quadratic
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


TIME_KEYS = (-1, 0, 1)


def stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = np.abs(vals[mask])
    if weights is not None:
        w = np.asarray(weights, dtype=float)[mask]
        w = np.where(np.isfinite(w) & (w > 0.0), w, 0.0)
    else:
        w = None
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "weighted_mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    if w is not None and np.sum(w) > 0.0:
        weighted = float(np.sum(w * vals) / np.sum(w))
    else:
        weighted = float(np.mean(vals))
    return {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "weighted_mean": weighted,
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }


def flat_q_from_case(case, mass: float) -> np.ndarray:
    u = case.u_cov
    x_flat = u[..., 0] ** 2 - u[..., 1] ** 2 - u[..., 2] ** 2
    return x_flat - float(mass) ** 2


def trace_coefficients(case) -> np.ndarray:
    atoms = base.tensor_atoms(case, "matter4")
    out = np.zeros(case.rho_a.shape + (len(base.ATOM_NAMES),), dtype=float)
    for apos, atom_name in enumerate(base.ATOM_NAMES):
        out[..., apos] = np.einsum("...ab,...ab->...", case.metric_inv, atoms[atom_name], optimize=True)
    return out


def make_closure_constraints(
    cases: dict[int, object],
    meta: dict[str, object],
    *,
    region: str,
    closure: str,
    q_scale: float,
) -> tuple[np.ndarray, list[str]]:
    col = np.asarray(meta["column_index"], dtype=int)
    n_base = int(meta["ncols"])
    extra_names: list[str] = []
    if closure in {"zero", "shear0"}:
        n_extra = 0
    elif closure == "q1":
        extra_names = ["trace_q1_alpha"]
        n_extra = 1
    elif closure == "q2":
        extra_names = ["trace_q1_alpha", "trace_q2_beta"]
        n_extra = 2
    else:
        raise ValueError(f"unknown closure: {closure}")

    rows: list[np.ndarray] = []
    for tpos, key in enumerate(TIME_KEYS):
        case = cases[key]
        coeffs = trace_coefficients(case)
        q_hat = flat_q_from_case(case, cases["mass"]) / max(abs(q_scale), 1.0e-300)
        mask = getattr(case, region)
        rho_w = np.sqrt(np.maximum(case.rho_a, 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
        for i_raw, j_raw in np.argwhere(mask):
            i = int(i_raw)
            j = int(j_raw)
            row = np.zeros(n_base + n_extra, dtype=float)
            if closure == "shear0":
                apos = base.ATOM_NAMES.index("ur")
                c = int(col[tpos, apos, i, j])
                if c >= 0:
                    row[c] = 1.0
            else:
                for apos in range(len(base.ATOM_NAMES)):
                    c = int(col[tpos, apos, i, j])
                    if c >= 0:
                        row[c] = coeffs[i, j, apos]
                if closure in {"q1", "q2"}:
                    row[n_base] = -q_hat[i, j]
                if closure == "q2":
                    row[n_base + 1] = -(q_hat[i, j] ** 2)
            norm = np.linalg.norm(row)
            if norm > 0.0:
                # The row is a hard equality.  This normalization improves the
                # nullspace SVD without changing the constraint surface.
                row *= float(rho_w[i, j]) / norm
            rows.append(row)
    if not rows:
        raise RuntimeError("trace closure produced no rows")
    return np.vstack(rows), extra_names


def solve_trace_closure(
    a_alg: np.ndarray,
    b_alg: np.ndarray,
    equality: np.ndarray,
    *,
    ridge: float,
) -> tuple[np.ndarray, dict[str, object]]:
    ncols = equality.shape[1]
    if a_alg.shape[1] < ncols:
        pad = np.zeros((a_alg.shape[0], ncols - a_alg.shape[1]), dtype=float)
        a_obj = np.hstack([a_alg, pad])
    else:
        a_obj = a_alg
    col_scale = np.linalg.norm(a_obj, axis=0) + np.linalg.norm(equality, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_s = a_obj / col_scale[None, :]
    f_s = equality / col_scale[None, :]
    h_reg = float(ridge) * np.eye(ncols)
    coeff_scaled, info = solve_hard_quadratic(a_s, b_alg, f_s, h_reg, solver="nullspace")
    coeff = coeff_scaled / col_scale
    return coeff, {
        **info,
        "ridge": float(ridge),
        "objective_residual_norm": float(np.linalg.norm(a_obj @ coeff - b_alg)),
        "equality_residual_norm": float(np.linalg.norm(equality @ coeff)),
        "equality_residual_max_abs": float(np.max(np.abs(equality @ coeff))) if equality.size else 0.0,
    }


def evaluate_closure(cases: dict[int, object], coeff: np.ndarray, meta: dict[str, object], closure: str, q_scale: float, region: str) -> dict[str, object]:
    n_base = int(meta["ncols"])
    coeff_base = coeff[:n_base]
    fields = base.coeff_to_fields(coeff_base, meta, cases[0].rho_a.shape)
    out: dict[str, object] = {}
    for key in TIME_KEYS:
        case = cases[key]
        coeffs = trace_coefficients(case)
        if closure == "shear0":
            value = fields[key][..., base.ATOM_NAMES.index("ur")]
            target = np.zeros_like(value)
        else:
            value = np.einsum("...a,...a->...", fields[key], coeffs, optimize=True)
            q_hat = flat_q_from_case(case, cases["mass"]) / max(abs(q_scale), 1.0e-300)
            target = np.zeros_like(value)
            if closure in {"q1", "q2"}:
                target += coeff[n_base] * q_hat
            if closure == "q2":
                target += coeff[n_base + 1] * q_hat**2
        mask = getattr(case, region)
        weights = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
        denom = max(float(np.percentile(np.abs(value[mask]), 95.0)), 1.0e-300)
        out[str(key)] = {
            "value": stats(value[mask], weights),
            "target": stats(target[mask], weights),
            "absolute_error": stats((value - target)[mask], weights),
            "relative_to_value_p95": stats((value - target)[mask] / denom, weights),
        }
    return out


def trace_nullspace_basis(trace_row: np.ndarray) -> np.ndarray:
    row = np.asarray(trace_row, dtype=float).reshape(1, 4)
    _, s, vh = np.linalg.svd(row, full_matrices=True)
    tol = 4 * np.finfo(float).eps * (float(s[0]) if s.size else 0.0)
    rank = int(np.count_nonzero(s > tol))
    return vh[rank:].T


def closure_row_at(case, closure: str, i: int, j: int) -> np.ndarray:
    if closure == "shear0":
        row = np.zeros(4, dtype=float)
        row[base.ATOM_NAMES.index("ur")] = 1.0
        return row
    return trace_coefficients(case)[i, j]


def analyze_reduced_symbol(cases: dict[int, object], region: str, closure: str) -> dict[str, object]:
    k_list = {
        "lab_t": np.array([1.0, 0.0, 0.0]),
        "lab_x": np.array([0.0, 1.0, 0.0]),
        "lab_z": np.array([0.0, 0.0, 1.0]),
    }
    records = []
    for key in TIME_KEYS:
        case = cases[key]
        mask = getattr(case, region)
        for k_name, k_cov in k_list.items():
            p = principal_matrix(case, k_cov)
            ranks = []
            conds = []
            for i, j in np.argwhere(mask):
                n = trace_nullspace_basis(closure_row_at(case, closure, int(i), int(j)))
                reduced = p[int(i), int(j)] @ n
                s = np.linalg.svd(reduced, compute_uv=False)
                tol = max(reduced.shape) * np.finfo(float).eps * (float(s[0]) if s.size else 0.0)
                rank = int(np.count_nonzero(s > tol))
                ranks.append(rank)
                conds.append(float(s[0] / max(s[-1], 1.0e-300)) if s.size else 0.0)
            records.append(
                {
                    "time_key": int(key),
                    "tau": float(case.tau_old),
                    "k_name": k_name,
                    "rank": stats_plain(np.asarray(ranks, dtype=float)),
                    "condition": stats_plain(np.asarray(conds, dtype=float)),
                }
            )
    return {"closure": closure, "records": records}


def render_plot(path: Path, records: list[dict[str, object]]) -> None:
    labels = [r["closure"] for r in records]
    center = [r["evaluation"]["by_time"]["0"]["algebraic_relative_residual"]["weighted_mean"] for r in records]
    minus = [r["evaluation"]["by_time"]["-1"]["algebraic_relative_residual"]["weighted_mean"] for r in records]
    plus = [r["evaluation"]["by_time"]["1"]["algebraic_relative_residual"]["weighted_mean"] for r in records]
    eq = [r["solve"]["equality_residual_max_abs"] for r in records]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    axes[0].bar(x - 0.25, minus, width=0.25, label="tau-probe")
    axes[0].bar(x, center, width=0.25, label="tau")
    axes[0].bar(x + 0.25, plus, width=0.25, label="tau+probe")
    axes[0].set_yscale("log")
    axes[0].set_title("algebraic residual weighted mean")
    axes[0].legend(fontsize=8)
    axes[1].bar(x, eq)
    axes[1].set_yscale("log")
    axes[1].set_title("hard equality max residual")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ("g", "uu", "rr", "ur")
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    probe = float(args.force_probe_dt_old)
    cases = {
        -1: build_full_reference_case(args, ref, float(args.tau) - probe),
        0: build_full_reference_case(args, ref, float(args.tau)),
        1: build_full_reference_case(args, ref, float(args.tau) + probe),
    }
    cases["mass"] = float(ref["params"].m)
    full = build_full_case(args, ref, float(args.tau))
    scale = ref["scale"]
    dt = probe * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    a, b, meta = base.assemble_system(
        cases,
        full,
        region=args.fit_region,
        force_region=args.force_region,
        target_floor_frac=float(args.target_floor_frac),
        force_weight=1.0,
        force_mode="full",
        force_erosion=int(args.force_mask_erosion),
        active_dilation=int(args.active_dilation),
        dt=dt,
        dx=dx,
        dz=dz,
    )
    n_force = 3 * int(meta["force_count"])
    a_alg = a[:-n_force]
    b_alg = b[:-n_force]
    force_rows = a[-n_force:]

    q_vals = []
    for key in TIME_KEYS:
        case = cases[key]
        mask = getattr(case, args.fit_region)
        q_vals.append(np.abs(flat_q_from_case(case, cases["mass"])[mask]))
    q_all = np.concatenate(q_vals)
    q_scale = float(np.percentile(q_all[np.isfinite(q_all)], 95.0)) if np.any(np.isfinite(q_all)) else 1.0
    q_scale = max(q_scale, 1.0e-300)

    records = []
    coeff_arrays = {}
    for closure in [x.strip() for x in args.closures.split(",") if x.strip()]:
        trace_rows, extra_names = make_closure_constraints(cases, meta, region=args.fit_region, closure=closure, q_scale=q_scale)
        if force_rows.shape[1] < trace_rows.shape[1]:
            pad = np.zeros((force_rows.shape[0], trace_rows.shape[1] - force_rows.shape[1]), dtype=float)
            force_eq = np.hstack([force_rows, pad])
        else:
            force_eq = force_rows
        equality = np.vstack([force_eq, trace_rows])
        coeff, solve = solve_trace_closure(a_alg, b_alg, equality, ridge=float(args.ridge))
        coeff_base = coeff[: int(meta["ncols"])]
        evaluation = base.evaluate_solution(
            cases,
            full,
            coeff_base,
            meta,
            region=args.fit_region,
            force_region=args.force_region,
            force_erosion=int(args.force_mask_erosion),
            dt=dt,
            dx=dx,
            dz=dz,
        )
        closure_eval = evaluate_closure(cases, coeff, meta, closure, q_scale, args.fit_region)
        record = {
            "closure": closure,
            "extra_names": extra_names,
            "extra_coefficients": {name: float(coeff[int(meta["ncols"]) + i]) for i, name in enumerate(extra_names)},
            "solve": solve,
            "evaluation": evaluation,
            "closure_diagnostic": closure_eval,
        }
        records.append(record)
        fields = base.coeff_to_fields(coeff_base, meta, cases[0].rho_a.shape)
        for key, suffix in [(-1, "m"), (0, "0"), (1, "p")]:
            for apos, atom in enumerate(base.ATOM_NAMES):
                coeff_arrays[f"{closure}_{atom}_{suffix}"] = fields[key][..., apos]

    symbols = {
        closure: analyze_reduced_symbol(cases, args.fit_region, closure)
        for closure in [x.strip() for x in args.closures.split(",") if x.strip()]
    }
    coeff_path = args.output / "gbcd_trace_closure_coefficients.npz"
    np.savez_compressed(coeff_path, **coeff_arrays)
    plot_path = args.output / "gbcd_trace_closure_summary.png"
    render_plot(plot_path, records)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "closures": [r["closure"] for r in records],
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "full_resolution": int(args.full_resolution),
            "q_definition": "Q_flat = X_flat - m^2, with X_flat=s_t^2-s_x^2-s_z^2 under current Q convention.",
            "q_scale_p95_abs": q_scale,
        },
        "definition": {
            "trace_closure": "C^mu_mu = F(Q), with F(0)=0. zero: F=0; q1: F=alpha Q/q_scale; q2: F=alpha Q/q_scale + beta (Q/q_scale)^2.",
            "shear0_closure": "D=0 in C_mn=A g_mn+B uu+C rr+D u_(m r_n), i.e. no u-r shear in this chosen matter basis.",
            "hard_constraints": "Full conservation rows plus trace-closure rows are imposed through a nullspace solve.",
            "residual": "Algebraic residual compares C_mn to Gtilde_mn-Ttilde_mn/Mp^2 on the fixed A-reference generated gtilde.",
        },
        "records": records,
        "reduced_principal_symbol": symbols,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "coefficients_npz": str(coeff_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def build_full_reference_case(args: argparse.Namespace, ref: dict[str, object], tau: float):
    from diagnose_mathcal_r_pure_geometry import build_case

    return build_case(args, ref, tau)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--closures", type=str, default="zero,q1,q2")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--ridge", type=float, default=1.0e-10)
    parser.add_argument("--full-resolution", type=int, default=96)
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
