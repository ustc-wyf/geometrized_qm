from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import fit_equation_first_constrained_bcd as base
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from diagnose_mathcal_r_pure_geometry import build_case
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case


TIME_KEYS = (-1, 0, 1)


def stats(values: np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    abs_vals = np.abs(vals)
    return {
        "count": int(vals.size),
        "mean": float(np.mean(abs_vals)),
        "p50": float(np.percentile(abs_vals, 50.0)),
        "p95": float(np.percentile(abs_vals, 95.0)),
        "max": float(np.max(abs_vals)),
    }


def split_system(a: np.ndarray, b: np.ndarray, meta: dict[str, object]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_force = 3 * int(meta["force_count"])
    if n_force <= 0:
        raise ValueError("full/transverse hard constraint requires force rows")
    return a[:-n_force], b[:-n_force], a[-n_force:]


def preliminary_solution(a_alg: np.ndarray, b_alg: np.ndarray, f_rows: np.ndarray, ridge_scaled: float) -> tuple[np.ndarray, np.ndarray]:
    col_scale = np.linalg.norm(a_alg, axis=0) + np.linalg.norm(f_rows, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_s = a_alg / col_scale[None, :]
    f_s = f_rows / col_scale[None, :]
    kkt = np.block(
        [
            [a_s.T @ a_s + float(ridge_scaled) * np.eye(a_s.shape[1]), f_s.T],
            [f_s, np.zeros((f_s.shape[0], f_s.shape[0]))],
        ]
    )
    rhs = np.concatenate([a_s.T @ b_alg, np.zeros(f_s.shape[0], dtype=float)])
    y = np.linalg.lstsq(kkt, rhs, rcond=1.0e-12)[0][: a_s.shape[1]]
    return y / col_scale, col_scale


def atom_scales_from_coeff(coeff: np.ndarray, meta: dict[str, object], atoms: tuple[str, ...]) -> dict[str, float]:
    col = np.asarray(meta["column_index"], dtype=int)
    scales = {}
    for apos, atom in enumerate(atoms):
        vals = []
        for tpos in range(len(TIME_KEYS)):
            ids = col[tpos, apos]
            mask = ids >= 0
            vals.append(coeff[ids[mask]])
        arr = np.concatenate(vals) if vals else np.asarray([1.0])
        scale = float(np.percentile(np.abs(arr[np.isfinite(arr)]), 95.0)) if np.any(np.isfinite(arr)) else 1.0
        scales[atom] = max(scale, 1.0e-30)
    return scales


def flat_q_from_case(case, mass: float) -> np.ndarray:
    u = case.u_cov
    x_flat = u[..., 0] ** 2 - u[..., 1] ** 2 - u[..., 2] ** 2
    return x_flat - float(mass) ** 2


def build_q_gate_weights(
    cases: dict[int, object],
    active: np.ndarray,
    *,
    mass: float,
    eps: float,
    power: float,
    max_weight: float,
) -> tuple[dict[int, np.ndarray], dict[str, float]]:
    vals = []
    for key in TIME_KEYS:
        vals.append(np.abs(flat_q_from_case(cases[key], mass)[active]))
    all_vals = np.concatenate(vals) if vals else np.asarray([1.0])
    finite = all_vals[np.isfinite(all_vals)]
    q_scale = float(np.percentile(finite, 95.0)) if finite.size else 1.0
    q_scale = max(q_scale, 1.0e-300)
    weights: dict[int, np.ndarray] = {}
    all_weights = []
    for key in TIME_KEYS:
        q_hat = np.abs(flat_q_from_case(cases[key], mass)) / q_scale
        w = 1.0 / np.maximum(float(eps) + q_hat, 1.0e-300) ** float(power)
        w = np.minimum(w, float(max_weight))
        weights[key] = w
        all_weights.append(w[active])
    w_all = np.concatenate(all_weights) if all_weights else np.asarray([1.0])
    return weights, {
        "q_scale_p95_abs": q_scale,
        "eps": float(eps),
        "power": float(power),
        "max_weight": float(max_weight),
        "weight_mean_active": float(np.mean(w_all[np.isfinite(w_all)])) if np.any(np.isfinite(w_all)) else 0.0,
        "weight_p95_active": float(np.percentile(w_all[np.isfinite(w_all)], 95.0)) if np.any(np.isfinite(w_all)) else 0.0,
        "weight_max_active": float(np.max(w_all[np.isfinite(w_all)])) if np.any(np.isfinite(w_all)) else 0.0,
    }


def add_reg_row_to_h(h: np.ndarray, entries: list[tuple[int, float]], weight: float) -> None:
    if weight <= 0.0:
        return
    for i, vi in entries:
        h[i, i] += weight * vi * vi
    for a in range(len(entries)):
        i, vi = entries[a]
        for b in range(a + 1, len(entries)):
            j, vj = entries[b]
            val = weight * vi * vj
            h[i, j] += val
            h[j, i] += val


def build_regularizer_h(
    meta: dict[str, object],
    col_scale: np.ndarray,
    atoms: tuple[str, ...],
    atom_scales: dict[str, float],
    *,
    w_norm: float,
    w_space: float,
    w_time: float,
    norm_point_weights: dict[int, np.ndarray] | None = None,
) -> tuple[np.ndarray, dict[str, int]]:
    col = np.asarray(meta["column_index"], dtype=int)
    active = np.asarray(meta["active_mask"], dtype=bool)
    ncols = int(meta["ncols"])
    h = np.zeros((ncols, ncols), dtype=float)
    counts = {"norm_rows": 0, "space_rows": 0, "time_rows": 0}

    if w_norm > 0.0:
        for tpos in range(len(TIME_KEYS)):
            key = TIME_KEYS[tpos]
            for apos, atom in enumerate(atoms):
                ids = col[tpos, apos]
                mask = ids >= 0
                coeff = 1.0 / (atom_scales[atom] * col_scale[ids[mask]])
                if norm_point_weights is not None:
                    point_w = np.asarray(norm_point_weights[key], dtype=float)[mask]
                else:
                    point_w = 1.0
                h[ids[mask], ids[mask]] += w_norm * point_w * coeff * coeff
                counts["norm_rows"] += int(np.count_nonzero(mask))

    if w_space > 0.0:
        nx, nz = active.shape
        for tpos in range(len(TIME_KEYS)):
            for apos, atom in enumerate(atoms):
                denom = atom_scales[atom]
                for i in range(nx):
                    for j in range(nz):
                        c0 = int(col[tpos, apos, i, j])
                        if c0 < 0:
                            continue
                        if i + 1 < nx:
                            c1 = int(col[tpos, apos, i + 1, j])
                            if c1 >= 0:
                                add_reg_row_to_h(
                                    h,
                                    [(c1, 1.0 / (denom * col_scale[c1])), (c0, -1.0 / (denom * col_scale[c0]))],
                                    w_space,
                                )
                                counts["space_rows"] += 1
                        if j + 1 < nz:
                            c1 = int(col[tpos, apos, i, j + 1])
                            if c1 >= 0:
                                add_reg_row_to_h(
                                    h,
                                    [(c1, 1.0 / (denom * col_scale[c1])), (c0, -1.0 / (denom * col_scale[c0]))],
                                    w_space,
                                )
                                counts["space_rows"] += 1

    if w_time > 0.0:
        for apos, atom in enumerate(atoms):
            denom = atom_scales[atom]
            for i, j in np.argwhere(active):
                for t0, t1 in [(0, 1), (1, 2)]:
                    c0 = int(col[t0, apos, i, j])
                    c1 = int(col[t1, apos, i, j])
                    if c0 >= 0 and c1 >= 0:
                        add_reg_row_to_h(
                            h,
                            [(c1, 1.0 / (denom * col_scale[c1])), (c0, -1.0 / (denom * col_scale[c0]))],
                            w_time,
                        )
                        counts["time_rows"] += 1
    return h, counts


def build_nullspace_cache(a_s: np.ndarray, b_alg: np.ndarray, f_s: np.ndarray) -> dict[str, object]:
    u, s, vh = np.linalg.svd(f_s, full_matrices=True)
    del u
    if s.size:
        tol = max(f_s.shape) * np.finfo(float).eps * float(s[0])
        rank = int(np.count_nonzero(s > tol))
    else:
        tol = 0.0
        rank = 0
    null = vh[rank:].T
    if null.shape[1] == 0:
        a_n = np.zeros((a_s.shape[0], 0), dtype=float)
        base_lhs = np.zeros((0, 0), dtype=float)
        base_rhs = np.zeros(0, dtype=float)
    else:
        a_n = a_s @ null
        base_lhs = a_n.T @ a_n
        base_rhs = a_n.T @ b_alg
    return {
        "null": null,
        "base_lhs": base_lhs,
        "base_rhs": base_rhs,
        "force_rank": rank,
        "force_singular_tol": float(tol),
        "force_singular_min": float(s[-1]) if s.size else 0.0,
        "force_singular_max": float(s[0]) if s.size else 0.0,
    }


def solve_hard_quadratic(
    a_s: np.ndarray,
    b_alg: np.ndarray,
    f_s: np.ndarray,
    h_reg: np.ndarray,
    *,
    solver: str,
    nullspace_cache: dict[str, object] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    h = a_s.T @ a_s + h_reg
    rhs_h = a_s.T @ b_alg
    if solver == "kkt":
        kkt = np.block([[h, f_s.T], [f_s, np.zeros((f_s.shape[0], f_s.shape[0]))]])
        rhs = np.concatenate([rhs_h, np.zeros(f_s.shape[0], dtype=float)])
        y = np.linalg.lstsq(kkt, rhs, rcond=1.0e-12)[0][: a_s.shape[1]]
        return y, {
            "hard_solver": "kkt",
            "scaled_force_residual_norm": float(np.linalg.norm(f_s @ y)),
            "scaled_force_residual_max_abs": float(np.max(np.abs(f_s @ y))) if f_s.size else 0.0,
        }

    if solver != "nullspace":
        raise ValueError(f"unknown hard solver: {solver}")

    # The constraint is homogeneous: f_s @ y = 0.  Solving in the numerical
    # nullspace keeps the conservation law hard even when the KKT matrix is
    # ill-conditioned.
    cache = nullspace_cache or build_nullspace_cache(a_s, b_alg, f_s)
    null = np.asarray(cache["null"], dtype=float)
    if null.shape[1] == 0:
        y = np.zeros(a_s.shape[1], dtype=float)
    else:
        lhs = np.asarray(cache["base_lhs"], dtype=float).copy()
        if np.any(h_reg):
            hn = h_reg @ null
            lhs = lhs + null.T @ hn
        rhs = np.asarray(cache["base_rhs"], dtype=float)
        z = np.linalg.lstsq(lhs, rhs, rcond=1.0e-12)[0]
        y = null @ z
    scaled_force = f_s @ y
    return y, {
        "hard_solver": "nullspace",
        "force_rank": int(cache["force_rank"]),
        "force_singular_tol": float(cache["force_singular_tol"]),
        "force_singular_min": float(cache["force_singular_min"]),
        "force_singular_max": float(cache["force_singular_max"]),
        "nullity": int(null.shape[1]),
        "scaled_force_residual_norm": float(np.linalg.norm(scaled_force)),
        "scaled_force_residual_max_abs": float(np.max(np.abs(scaled_force))) if scaled_force.size else 0.0,
    }


def solve_with_gauge(
    a_alg: np.ndarray,
    b_alg: np.ndarray,
    f_rows: np.ndarray,
    meta: dict[str, object],
    atoms: tuple[str, ...],
    *,
    w_norm: float,
    w_space: float,
    w_time: float,
    prelim_ridge: float,
    hard_solver: str,
    col_scale: np.ndarray | None = None,
    atom_scales: dict[str, float] | None = None,
    nullspace_cache: dict[str, object] | None = None,
    norm_point_weights: dict[int, np.ndarray] | None = None,
    norm_weight_info: dict[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    if col_scale is None or atom_scales is None:
        prelim, col_scale = preliminary_solution(a_alg, b_alg, f_rows, prelim_ridge)
        atom_scales = atom_scales_from_coeff(prelim, meta, atoms)
    a_s = a_alg / col_scale[None, :]
    f_s = f_rows / col_scale[None, :]
    h_reg, reg_counts = build_regularizer_h(
        meta,
        col_scale,
        atoms,
        atom_scales,
        w_norm=w_norm,
        w_space=w_space,
        w_time=w_time,
        norm_point_weights=norm_point_weights,
    )
    y, solver_info = solve_hard_quadratic(
        a_s,
        b_alg,
        f_s,
        h_reg,
        solver=hard_solver,
        nullspace_cache=nullspace_cache,
    )
    coeff = y / col_scale
    alg_res = a_alg @ coeff - b_alg
    force_res = f_rows @ coeff
    info = {
        "w_norm": float(w_norm),
        "w_space": float(w_space),
        "w_time": float(w_time),
        "prelim_ridge": float(prelim_ridge),
        "atom_scales": atom_scales,
        "norm_weight_info": norm_weight_info or {},
        "regularizer_counts": reg_counts,
        **solver_info,
        "algebraic_weighted_residual_norm": float(np.linalg.norm(alg_res)),
        "force_residual_norm": float(np.linalg.norm(force_res)),
        "force_residual_max_abs": float(np.max(np.abs(force_res))) if force_res.size else 0.0,
    }
    return coeff, info


def coefficient_smoothness(coeff: np.ndarray, meta: dict[str, object], atoms: tuple[str, ...], atom_scales: dict[str, float]) -> dict[str, object]:
    fields = base.coeff_to_fields(coeff, meta, np.asarray(meta["active_mask"], dtype=bool).shape)
    active = np.asarray(meta["active_mask"], dtype=bool)
    out = {}
    for apos, atom in enumerate(atoms):
        scale = atom_scales[atom]
        vals0 = fields[0][..., apos][active] / scale
        spatial_diffs = []
        for key in TIME_KEYS:
            f = fields[key][..., apos] / scale
            mask = active
            spatial_diffs.append((f[1:, :] - f[:-1, :])[mask[1:, :] & mask[:-1, :]])
            spatial_diffs.append((f[:, 1:] - f[:, :-1])[mask[:, 1:] & mask[:, :-1]])
        spatial = np.concatenate([x for x in spatial_diffs if x.size]) if any(x.size for x in spatial_diffs) else np.asarray([])
        time_diffs = []
        for key0, key1 in [(-1, 0), (0, 1)]:
            time_diffs.append((fields[key1][..., apos][active] - fields[key0][..., apos][active]) / scale)
        temporal = np.concatenate(time_diffs) if time_diffs else np.asarray([])
        out[atom] = {
            "central_coeff_normalized": stats(vals0),
            "spatial_neighbor_diff_normalized": stats(spatial),
            "time_neighbor_diff_normalized": stats(temporal),
        }
    return out


def render_summary(out_path: Path, records: list[dict[str, object]]) -> None:
    labels = [r["name"] for r in records]
    alg = [r["evaluation"]["by_time"]["0"]["algebraic_relative_residual"].get("weighted_mean", 0.0) for r in records]
    force = [r["evaluation"]["central_force"]["transverse_over_derivative_scale"].get("weighted_mean", 0.0) for r in records]
    smooth = []
    temporal = []
    for r in records:
        vals = []
        tv = []
        for atom_report in r["smoothness"].values():
            vals.append(atom_report["spatial_neighbor_diff_normalized"]["p95"])
            tv.append(atom_report["time_neighbor_diff_normalized"]["p95"])
        smooth.append(float(np.nanmean(vals)))
        temporal.append(float(np.nanmean(tv)))
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4), constrained_layout=True)
    axes[0].bar(x, alg)
    axes[0].set_yscale("log")
    axes[0].set_title("central algebraic residual")
    axes[1].bar(x, force)
    axes[1].set_yscale("log")
    axes[1].set_title("transverse force scale")
    axes[2].bar(x - 0.18, smooth, width=0.36, label="space p95")
    axes[2].bar(x + 0.18, temporal, width=0.36, label="time p95")
    axes[2].set_yscale("log")
    axes[2].set_title("normalized coefficient roughness")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = tuple(x.strip() for x in args.atoms.split(",") if x.strip())
    atoms = base.ATOM_NAMES
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
        force_mode=args.force_mode,
        force_erosion=int(args.force_mask_erosion),
        active_dilation=int(args.active_dilation),
        dt=dt,
        dx=dx,
        dz=dz,
    )
    a_alg, b_alg, f_rows = split_system(a, b, meta)
    prelim, col_scale = preliminary_solution(a_alg, b_alg, f_rows, float(args.prelim_ridge))
    atom_scales = atom_scales_from_coeff(prelim, meta, atoms)
    a_s = a_alg / col_scale[None, :]
    f_s = f_rows / col_scale[None, :]
    nullspace_cache = build_nullspace_cache(a_s, b_alg, f_s) if args.hard_solver == "nullspace" else None
    norm_point_weights = None
    norm_weight_info: dict[str, float] = {}
    if args.q_gated_norm:
        norm_point_weights, norm_weight_info = build_q_gate_weights(
            cases,
            np.asarray(meta["active_mask"], dtype=bool),
            mass=float(ref["params"].m),
            eps=float(args.q_gate_eps),
            power=float(args.q_gate_power),
            max_weight=float(args.q_gate_max),
        )
    all_specs = {
        "norm": (float(args.norm_weight), 0.0, 0.0),
        "space": (float(args.norm_weight), float(args.space_weight), 0.0),
        "time": (float(args.norm_weight), 0.0, float(args.time_weight)),
        "space_time": (float(args.norm_weight), float(args.space_weight), float(args.time_weight)),
    }
    selected = [x.strip() for x in args.gauge_cases.split(",") if x.strip()]
    specs = {name: all_specs[name] for name in selected}
    records = []
    coeff_arrays = {}
    for name, (w_norm, w_space, w_time) in specs.items():
        coeff, info = solve_with_gauge(
            a_alg,
            b_alg,
            f_rows,
            meta,
            atoms,
            w_norm=w_norm,
            w_space=w_space,
            w_time=w_time,
            prelim_ridge=float(args.prelim_ridge),
            hard_solver=args.hard_solver,
            col_scale=col_scale,
            atom_scales=atom_scales,
            nullspace_cache=nullspace_cache,
            norm_point_weights=norm_point_weights,
            norm_weight_info=norm_weight_info,
        )
        evaluation = base.evaluate_solution(
            cases,
            full,
            coeff,
            meta,
            region=args.fit_region,
            force_region=args.force_region,
            force_erosion=int(args.force_mask_erosion),
            dt=dt,
            dx=dx,
            dz=dz,
        )
        smoothness = coefficient_smoothness(coeff, meta, atoms, info["atom_scales"])
        records.append({"name": name, "gauge": info, "evaluation": evaluation, "smoothness": smoothness})
        fields = base.coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
        for key, suffix in [(-1, "m"), (0, "0"), (1, "p")]:
            for apos, atom in enumerate(atoms):
                coeff_arrays[f"{name}_{atom}_{suffix}"] = fields[key][..., apos]
    coeff_arrays["active_mask"] = np.asarray(meta["active_mask"], dtype=bool)
    coeff_arrays["force_mask"] = np.asarray(meta["force_mask"], dtype=bool)
    coeff_path = args.output / "gbcd_auxiliary_gauge_coefficients.npz"
    np.savez_compressed(coeff_path, **coeff_arrays)
    plot_path = args.output / "gbcd_auxiliary_gauge_summary.png"
    render_summary(plot_path, records)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "atoms": list(atoms),
            "force_mode": args.force_mode,
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "norm_weight": float(args.norm_weight),
            "space_weight": float(args.space_weight),
            "time_weight": float(args.time_weight),
            "prelim_ridge": float(args.prelim_ridge),
            "hard_solver": args.hard_solver,
            "gauge_cases": selected,
            "q_gated_norm": bool(args.q_gated_norm),
            "q_gate": norm_weight_info,
        },
        "definition": {
            "goal": "Choose a representative gauge for auxiliary anisotropic stress coefficients A,B,C,D under the same hard conservation constraint.",
            "norm": "minimize normalized coefficient magnitude.",
            "space": "minimize normalized spatial neighbor differences.",
            "time": "minimize normalized time-layer differences across tau-probe, tau, tau+probe.",
            "space_time": "combine space and time smoothing.",
        },
        "records": records,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "coefficients_npz": str(coeff_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--atoms", type=str, default="g,uu,rr,ur")
    parser.add_argument("--force-mode", choices=["transverse", "full"], default="full")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--norm-weight", type=float, default=1.0e-8)
    parser.add_argument("--space-weight", type=float, default=1.0e-6)
    parser.add_argument("--time-weight", type=float, default=1.0e-6)
    parser.add_argument("--prelim-ridge", type=float, default=1.0e-8)
    parser.add_argument("--hard-solver", choices=["nullspace", "kkt"], default="nullspace")
    parser.add_argument("--gauge-cases", type=str, default="norm,space,time,space_time")
    parser.add_argument("--q-gated-norm", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--q-gate-eps", type=float, default=0.03)
    parser.add_argument("--q-gate-power", type=float, default=2.0)
    parser.add_argument("--q-gate-max", type=float, default=1.0e4)
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
