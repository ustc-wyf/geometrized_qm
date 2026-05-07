from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from diagnose_mathcal_r_pure_geometry import build_case, tensor_norm
from fit_gbcd_auxiliary_gauge import stats
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


TIME_KEYS = (-1, 0, 1)
ATOM_NAMES = ("g", "uu", "rr", "ur")


def weighted_stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    base_stats = stats(values)
    vals = np.asarray(values, dtype=float)
    finite = np.isfinite(vals)
    if weights is None:
        base_stats["weighted_mean"] = base_stats["mean"]
        return base_stats
    w = np.asarray(weights, dtype=float)
    finite &= np.isfinite(w) & (w >= 0.0)
    if not np.any(finite):
        base_stats["weighted_mean"] = 0.0
        return base_stats
    abs_vals = np.abs(vals[finite])
    wf = w[finite]
    denom = float(np.sum(wf))
    base_stats["weighted_mean"] = float(np.sum(wf * abs_vals) / denom) if denom > 0.0 else float(np.mean(abs_vals))
    return base_stats


def split_system(a: np.ndarray, b: np.ndarray, meta: dict[str, object], region_count: int) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    n_per = int(region_count) * len(base.SYMMETRIC_COMPONENTS)
    n_force = 3 * int(meta["force_count"])
    return {
        "alg_m": (a[:n_per], b[:n_per]),
        "alg_0": (a[n_per : 2 * n_per], b[n_per : 2 * n_per]),
        "alg_p": (a[2 * n_per : 3 * n_per], b[2 * n_per : 3 * n_per]),
        "force": (a[-n_force:], b[-n_force:]),
    }


def load_full_coeff(path: Path, meta: dict[str, object], case_name: str, shape: tuple[int, int]) -> np.ndarray:
    data = np.load(path)
    col = np.asarray(meta["column_index"], dtype=int)
    coeff = np.zeros(int(meta["ncols"]), dtype=float)
    suffix_by_tpos = {0: "m", 1: "0", 2: "p"}
    for tpos in range(3):
        suffix = suffix_by_tpos[tpos]
        for apos, atom in enumerate(ATOM_NAMES):
            field = np.asarray(data[f"{case_name}_{atom}_{suffix}"], dtype=float)
            if field.shape != shape:
                raise ValueError(f"shape mismatch for {case_name}_{atom}_{suffix}: {field.shape} != {shape}")
            ids = col[tpos, apos]
            mask = ids >= 0
            coeff[ids[mask]] = field[mask]
    return coeff


def plus_columns(meta: dict[str, object]) -> tuple[np.ndarray, dict[int, int]]:
    col = np.asarray(meta["column_index"], dtype=int)
    ids = np.unique(col[2][col[2] >= 0])
    ids = ids.astype(int)
    mapping = {int(c): n for n, c in enumerate(ids)}
    return ids, mapping


def compact_rows(a_full: np.ndarray, unknown_ids: np.ndarray) -> np.ndarray:
    return a_full[:, unknown_ids]


def known_subtracted_rhs(a_full: np.ndarray, b_full: np.ndarray, coeff_known: np.ndarray, unknown_ids: np.ndarray) -> np.ndarray:
    known_coeff = coeff_known.copy()
    known_coeff[unknown_ids] = 0.0
    return b_full - a_full @ known_coeff


def atom_scales_from_known(coeff_full: np.ndarray, meta: dict[str, object]) -> dict[str, float]:
    col = np.asarray(meta["column_index"], dtype=int)
    out = {}
    for apos, atom in enumerate(ATOM_NAMES):
        vals = []
        for tpos in range(3):
            ids = col[tpos, apos]
            mask = ids >= 0
            vals.append(coeff_full[ids[mask]])
        arr = np.concatenate(vals)
        finite = arr[np.isfinite(arr)]
        out[atom] = max(float(np.percentile(np.abs(finite), 95.0)) if finite.size else 1.0, 1.0e-30)
    return out


def add_regularization_rows(
    meta: dict[str, object],
    coeff_known: np.ndarray,
    unknown_ids: np.ndarray,
    mapping: dict[int, int],
    *,
    norm_weight: float,
    time_weight: float,
    target_mode: str,
    atom_scales: dict[str, float],
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    col = np.asarray(meta["column_index"], dtype=int)
    active = np.asarray(meta["active_mask"], dtype=bool)
    rows: list[np.ndarray] = []
    rhs: list[float] = []
    counts = {"norm_rows": 0, "time_rows": 0}
    n = len(unknown_ids)
    sqrt_norm = float(np.sqrt(max(norm_weight, 0.0)))
    sqrt_time = float(np.sqrt(max(time_weight, 0.0)))
    for apos, atom in enumerate(ATOM_NAMES):
        scale = atom_scales[atom]
        for i, j in np.argwhere(active):
            c_p = int(col[2, apos, int(i), int(j)])
            if c_p < 0:
                continue
            compact = mapping[c_p]
            if sqrt_norm > 0.0:
                row = np.zeros(n, dtype=float)
                row[compact] = sqrt_norm / scale
                rows.append(row)
                rhs.append(0.0)
                counts["norm_rows"] += 1
            if sqrt_time > 0.0:
                c_0 = int(col[1, apos, int(i), int(j)])
                c_m = int(col[0, apos, int(i), int(j)])
                if c_0 >= 0:
                    if target_mode == "hold":
                        target = coeff_known[c_0]
                    elif target_mode == "inertial" and c_m >= 0:
                        target = 2.0 * coeff_known[c_0] - coeff_known[c_m]
                    else:
                        target = coeff_known[c_0]
                    row = np.zeros(n, dtype=float)
                    row[compact] = sqrt_time / scale
                    rows.append(row)
                    rhs.append(float(sqrt_time * target / scale))
                    counts["time_rows"] += 1
    if not rows:
        return np.zeros((0, n), dtype=float), np.zeros(0, dtype=float), counts
    return np.vstack(rows), np.asarray(rhs, dtype=float), counts


def solve_nonhomogeneous_hard(a: np.ndarray, b: np.ndarray, f: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    col_scale = np.linalg.norm(a, axis=0) + np.linalg.norm(f, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_s = a / col_scale[None, :]
    f_s = f / col_scale[None, :]
    u, s, vh = np.linalg.svd(f_s, full_matrices=True)
    if s.size:
        tol = max(f_s.shape) * np.finfo(float).eps * float(s[0])
        rank = int(np.count_nonzero(s > tol))
    else:
        tol = 0.0
        rank = 0
    if rank:
        y0 = vh[:rank].T @ ((u[:, :rank].T @ d) / s[:rank])
    else:
        y0 = np.zeros(a.shape[1], dtype=float)
    null = vh[rank:].T
    if null.shape[1] > 0:
        lhs = a_s @ null
        rhs = b - a_s @ y0
        z = np.linalg.lstsq(lhs, rhs, rcond=1.0e-12)[0]
        y = y0 + null @ z
    else:
        y = y0
    x = y / col_scale
    eq = f @ x - d
    obj = a @ x - b
    return x, {
        "rank": rank,
        "nullity": int(null.shape[1]),
        "singular_min": float(s[-1]) if s.size else 0.0,
        "singular_max": float(s[0]) if s.size else 0.0,
        "singular_tol": float(tol),
        "equality_norm": float(np.linalg.norm(eq)),
        "equality_max_abs": float(np.max(np.abs(eq))) if eq.size else 0.0,
        "objective_norm": float(np.linalg.norm(obj)),
    }


def coeff_full_with_prediction(coeff_known: np.ndarray, unknown_ids: np.ndarray, x_plus: np.ndarray) -> np.ndarray:
    out = coeff_known.copy()
    out[unknown_ids] = x_plus
    return out


def evaluate_plus(cases: dict[int, object], coeff_full: np.ndarray, meta: dict[str, object], region: str) -> dict[str, object]:
    fields = base.coeff_to_fields(coeff_full, meta, cases[0].rho_a.shape)
    cov = base.coeff_fields_to_cov(cases, fields)
    case = cases[1]
    mask = getattr(case, region)
    rel = tensor_norm(cov[1] - case.r_need) / np.maximum(tensor_norm(case.r_need), 1.0e-300)
    weights = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    return {
        "plus_algebraic_relative_residual": weighted_stats(rel[mask], weights),
        "plus_model_norm": weighted_stats(tensor_norm(cov[1])[mask], weights),
        "plus_target_norm": weighted_stats(tensor_norm(case.r_need)[mask], weights),
    }


def coefficient_difference(coeff_pred_full: np.ndarray, coeff_ref_full: np.ndarray, meta: dict[str, object]) -> dict[str, object]:
    fields_p = base.coeff_to_fields(coeff_pred_full, meta, np.asarray(meta["active_mask"], dtype=bool).shape)
    fields_r = base.coeff_to_fields(coeff_ref_full, meta, np.asarray(meta["active_mask"], dtype=bool).shape)
    active = np.asarray(meta["active_mask"], dtype=bool)
    out = {}
    for apos, atom in enumerate(ATOM_NAMES):
        ref = fields_r[1][..., apos][active]
        pred = fields_p[1][..., apos][active]
        scale = max(float(np.percentile(np.abs(ref[np.isfinite(ref)]), 95.0)) if np.any(np.isfinite(ref)) else 1.0, 1.0e-300)
        out[atom] = stats((pred - ref) / scale)
    return out


def render_plot(path: Path, records: list[dict[str, object]]) -> None:
    labels = [r["target_mode"] for r in records]
    alg = [r["evaluation"]["plus_algebraic_relative_residual"]["weighted_mean"] for r in records]
    diff = []
    for r in records:
        diff.append(float(np.mean([v["p95"] for v in r["plus_coeff_difference"].values()])))
    eq = [r["solve"]["equality_max_abs"] for r in records]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), constrained_layout=True)
    axes[0].bar(x, alg)
    axes[0].set_yscale("log")
    axes[0].set_title("one-step plus residual")
    axes[1].bar(x, diff)
    axes[1].set_yscale("log")
    axes[1].set_title("plus coeff diff p95 avg")
    axes[2].bar(x, eq)
    axes[2].set_yscale("log")
    axes[2].set_title("conservation equality max")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=20, ha="right")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def default_coeff_path(tau: float, time_weight: float) -> Path:
    if abs(time_weight - 1.0e-5) > 1.0e-16:
        raise ValueError("default_coeff_path currently only knows time_weight=1e-5 outputs")
    if tau == 0.0:
        return Path("visualizations/equation_first_gbcd_aux_noq_time_scan_n96_tau0_tw_1em5/gbcd_auxiliary_gauge_coefficients.npz")
    tag = ("taum" if tau < 0 else "taup") + str(abs(tau)).replace(".", "p")
    return Path(f"visualizations/equation_first_gbcd_aux_noq_time1em5_n96_{tag}/gbcd_auxiliary_gauge_coefficients.npz")


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ATOM_NAMES
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    probe = float(args.force_probe_dt_old)
    cases = {
        -1: build_case(args, ref, float(args.tau) - probe),
        0: build_case(args, ref, float(args.tau)),
        1: build_case(args, ref, float(args.tau) + probe),
    }
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
    region_count = int(np.count_nonzero(getattr(cases[0], args.fit_region)))
    blocks = split_system(a, b, meta, region_count)
    unknown_ids, mapping = plus_columns(meta)
    coeff_path = args.coefficients or default_coeff_path(float(args.tau), float(args.time_weight))
    coeff_ref = load_full_coeff(coeff_path, meta, "time", cases[0].rho_a.shape)
    coeff_known = coeff_ref.copy()
    coeff_known[unknown_ids] = 0.0

    a_plus_full, b_plus_full = blocks["alg_p"]
    a_plus = compact_rows(a_plus_full, unknown_ids)
    b_plus = known_subtracted_rhs(a_plus_full, b_plus_full, coeff_known, unknown_ids)
    f_full, d_full = blocks["force"]
    f_plus = compact_rows(f_full, unknown_ids)
    d_plus = known_subtracted_rhs(f_full, d_full, coeff_known, unknown_ids)

    atom_scales = atom_scales_from_known(coeff_ref, meta)
    records = []
    coeff_out = {}
    for target_mode in [x.strip() for x in args.target_modes.split(",") if x.strip()]:
        a_reg, b_reg, reg_counts = add_regularization_rows(
            meta,
            coeff_ref,
            unknown_ids,
            mapping,
            norm_weight=float(args.norm_weight),
            time_weight=float(args.time_weight),
            target_mode=target_mode,
            atom_scales=atom_scales,
        )
        a_obj = np.vstack([a_plus, a_reg])
        b_obj = np.concatenate([b_plus, b_reg])
        x_plus, solve = solve_nonhomogeneous_hard(a_obj, b_obj, f_plus, d_plus)
        coeff_pred = coeff_full_with_prediction(coeff_known, unknown_ids, x_plus)
        evaluation = evaluate_plus(cases, coeff_pred, meta, args.fit_region)
        diff = coefficient_difference(coeff_pred, coeff_ref, meta)
        records.append(
            {
                "target_mode": target_mode,
                "regularizer_counts": reg_counts,
                "solve": solve,
                "evaluation": evaluation,
                "plus_coeff_difference": diff,
            }
        )
        fields = base.coeff_to_fields(coeff_pred, meta, cases[0].rho_a.shape)
        for apos, atom in enumerate(ATOM_NAMES):
            coeff_out[f"{target_mode}_{atom}_p"] = fields[1][..., apos]
    coeff_out["active_mask"] = np.asarray(meta["active_mask"], dtype=bool)
    coeff_out_path = args.output / "gbcd_v0_one_step_coefficients.npz"
    np.savez_compressed(coeff_out_path, **coeff_out)
    plot_path = args.output / "gbcd_v0_one_step_summary.png"
    render_plot(plot_path, records)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "norm_weight": float(args.norm_weight),
            "time_weight": float(args.time_weight),
            "coefficients": str(coeff_path.resolve()),
            "target_modes": [r["target_mode"] for r in records],
        },
        "definition": {
            "goal": "Given lambda_- and lambda_0 from the v0 three-layer projection, solve only lambda_+ with plus algebraic residual, center full-conservation hard constraint, and v0 norm/time regularization.",
            "hold": "time regularizer pulls lambda_+ toward lambda_0.",
            "inertial": "time regularizer pulls lambda_+ toward 2 lambda_0-lambda_-.",
        },
        "records": records,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "coefficients_npz": str(coeff_out_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--coefficients", type=Path, default=None)
    parser.add_argument("--target-modes", type=str, default="hold,inertial")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--norm-weight", type=float, default=1.0e-8)
    parser.add_argument("--time-weight", type=float, default=1.0e-5)
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
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
