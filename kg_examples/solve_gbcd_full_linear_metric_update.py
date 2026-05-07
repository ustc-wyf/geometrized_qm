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
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from test_gbcd_v0_one_step_lambda import ATOM_NAMES, weighted_stats


COMPONENTS = base.SYMMETRIC_COMPONENTS


def sym_to_vec(tensor: np.ndarray) -> np.ndarray:
    return np.asarray([tensor[a, b] for a, b in COMPONENTS], dtype=float)


def vec_to_sym(values: np.ndarray) -> np.ndarray:
    out = np.zeros((3, 3), dtype=float)
    for value, (a, b) in zip(values, COMPONENTS):
        out[a, b] = float(value)
        out[b, a] = float(value)
    return out


def add_symmetric_component(arr: np.ndarray, a: int, b: int, value: float) -> None:
    arr[..., a, b] += value
    if a != b:
        arr[..., b, a] += value


def einstein_point(metric_cov: np.ndarray, dg: np.ndarray, d2g: np.ndarray) -> np.ndarray:
    ginv = np.linalg.pinv(metric_cov, rcond=1.0e-12, hermitian=True)
    gamma1 = np.zeros((3, 3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gamma1[a, b, c] = 0.5 * (dg[b, a, c] + dg[c, a, b] - dg[a, b, c])
    gamma2 = np.einsum("ad,dbc->abc", ginv, gamma1, optimize=True)

    dginv = np.zeros((3, 3, 3), dtype=float)
    for e in range(3):
        dginv[e] = -ginv @ dg[e] @ ginv

    dgamma1 = np.zeros((3, 3, 3, 3), dtype=float)
    for e in range(3):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    dgamma1[e, a, b, c] = 0.5 * (
                        d2g[e, b, a, c] + d2g[e, c, a, b] - d2g[e, a, b, c]
                    )

    dgamma2 = np.zeros((3, 3, 3, 3), dtype=float)
    for e in range(3):
        dgamma2[e] = np.einsum("ad,dbc->abc", dginv[e], gamma1, optimize=True) + np.einsum(
            "ad,dbc->abc", ginv, dgamma1[e], optimize=True
        )

    ricci = np.zeros((3, 3), dtype=float)
    for a in range(3):
        for b in range(3):
            term1 = sum(dgamma2[c, c, a, b] for c in range(3))
            term2 = sum(dgamma2[b, c, a, c] for c in range(3))
            quad1 = 0.0
            quad2 = 0.0
            for c in range(3):
                for d in range(3):
                    quad1 += gamma2[c, a, b] * gamma2[d, c, d]
                    quad2 += gamma2[c, a, d] * gamma2[d, b, c]
            ricci[a, b] = term1 - term2 + quad1 - quad2
    scalar = float(np.einsum("ab,ab->", ginv, ricci))
    return ricci - 0.5 * metric_cov * scalar


def metric_variable_index(mask: np.ndarray) -> tuple[list[tuple[int, int, int, int]], np.ndarray]:
    variables: list[tuple[int, int, int, int]] = []
    col = -np.ones(mask.shape + (3, 3), dtype=int)
    for i, j in np.argwhere(mask):
        for a, b in COMPONENTS:
            n = len(variables)
            variables.append((int(i), int(j), int(a), int(b)))
            col[int(i), int(j), int(a), int(b)] = n
            col[int(i), int(j), int(b), int(a)] = n
    return variables, col


def affected_points(i: int, j: int, shape: tuple[int, int]) -> list[tuple[int, int]]:
    pts = {(i, j)}
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ii = i + di
        jj = j + dj
        if 0 <= ii < shape[0] and 0 <= jj < shape[1]:
            pts.add((ii, jj))
    return sorted(pts)


def perturb_jets_for_metric_plus_variable(
    point: tuple[int, int],
    var: tuple[int, int, int, int],
    *,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray]:
    pi, pj = point
    qi, qj, a, b = var
    ddg = np.zeros((3, 3, 3), dtype=float)
    dd2g = np.zeros((3, 3, 3, 3), dtype=float)
    if pi == qi and pj == qj:
        add_symmetric_component(ddg[0], a, b, 1.0 / (2.0 * dt))
        add_symmetric_component(dd2g[0, 0], a, b, 1.0 / (dt * dt))
    if pj == qj and abs(pi - qi) == 1:
        # d/dx metric_t at row pi uses central difference
        # (metric_t[pi+1]-metric_t[pi-1])/(2 dx).
        sign = 1.0 if qi == pi + 1 else -1.0
        coeff = sign / (4.0 * dx * dt)
        add_symmetric_component(dd2g[0, 1], a, b, coeff)
        add_symmetric_component(dd2g[1, 0], a, b, coeff)
    if pi == qi and abs(pj - qj) == 1:
        sign = 1.0 if qj == pj + 1 else -1.0
        coeff = sign / (4.0 * dz * dt)
        add_symmetric_component(dd2g[0, 2], a, b, coeff)
        add_symmetric_component(dd2g[2, 0], a, b, coeff)
    return ddg, dd2g


def build_full_linear_operator(
    *,
    metric_cov: np.ndarray,
    dg_base: np.ndarray,
    d2g_base: np.ndarray,
    variables: list[tuple[int, int, int, int]],
    row_points: list[tuple[int, int]],
    row_weights: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    eps: float,
) -> np.ndarray:
    point_to_row = {pt: n for n, pt in enumerate(row_points)}
    nrows = len(row_points) * len(COMPONENTS)
    mat = np.zeros((nrows, len(variables)), dtype=float)
    for col, var in enumerate(variables):
        for pt in affected_points(var[0], var[1], metric_cov.shape[:2]):
            row_index = point_to_row.get(pt)
            if row_index is None:
                continue
            ddg, dd2g = perturb_jets_for_metric_plus_variable(pt, var, dt=dt, dx=dx, dz=dz)
            i, j = pt
            e_plus = einstein_point(metric_cov[i, j], dg_base[i, j] + eps * ddg, d2g_base[i, j] + eps * dd2g)
            e_minus = einstein_point(metric_cov[i, j], dg_base[i, j] - eps * ddg, d2g_base[i, j] - eps * dd2g)
            deriv = (e_plus - e_minus) / (2.0 * eps)
            mat[row_index * len(COMPONENTS) : (row_index + 1) * len(COMPONENTS), col] = (
                row_weights[row_index] * sym_to_vec(deriv)
            )
    return mat


def solve_scaled_lstsq(a: np.ndarray, b: np.ndarray, ridge: float) -> tuple[np.ndarray, dict[str, object]]:
    col_scale = np.linalg.norm(a, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_s = a / col_scale[None, :]
    if ridge > 0.0:
        a_solve = np.vstack([a_s, float(ridge) * np.eye(a_s.shape[1])])
        b_solve = np.concatenate([b, np.zeros(a_s.shape[1], dtype=float)])
    else:
        a_solve = a_s
        b_solve = b
    y, residuals, rank, singular = np.linalg.lstsq(a_solve, b_solve, rcond=1.0e-12)
    x = y / col_scale
    residual = a @ x - b
    return x, {
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "condition_scaled": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_relative_residual": float(np.linalg.norm(residual) / max(np.linalg.norm(b), 1.0e-300)),
        "weighted_residual_norm": float(np.linalg.norm(residual)),
        "weighted_rhs_norm": float(np.linalg.norm(b)),
        "ridge": float(ridge),
        "lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
    }


def render(path: Path, x: np.ndarray, z: np.ndarray, mask: np.ndarray, before: np.ndarray, linear: np.ndarray, exact: np.ndarray) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(1, 3, figsize=(14.7, 4.5), constrained_layout=True)
    fields = [
        ("before full-linear update", before, "magma"),
        ("linearized residual after update", linear, "magma"),
        ("exact nonlinear residual after update", exact, "magma"),
    ]
    for ax, (title, field, cmap) in zip(axes, fields):
        finite = field[mask & np.isfinite(field)]
        vmax = max(float(np.percentile(finite, 95.0)) if finite.size else 1.0, 1.0e-16)
        im = ax.pcolormesh(xg, zg, np.clip(field, 0.0, vmax), shading="auto", cmap=cmap, vmin=0.0, vmax=vmax)
        ax.contour(xg, zg, mask.astype(float), levels=[0.5], colors="white", linewidths=0.7)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    base.ATOM_NAMES = ATOM_NAMES
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    old_scale = float(ref["scale"].old_dimensionless_scale_ev_inv)
    dt = float(args.force_probe_dt_old) * old_scale
    cases = {
        key: build_case(args, ref, float(args.tau) + key * float(args.force_probe_dt_old))
        for key in base.TIME_KEYS
    }
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])

    coeff_path = args.coefficients
    if coeff_path is None:
        tag = "tau0" if abs(float(args.tau)) < 1.0e-12 else ("taum3p5" if float(args.tau) < 0.0 else "taup3p5")
        coeff_path = Path(f"visualizations/equation_first_gbcd_principal_constraint_projection_n96_{tag}/principal_constraint_coefficients.npz")
    coeff_data = np.load(coeff_path)
    fields = {}
    for key, suffix in [(-1, "m"), (0, "0"), (1, "p")]:
        arr = np.zeros(cases[0].rho_a.shape + (len(ATOM_NAMES),), dtype=float)
        for apos, atom in enumerate(ATOM_NAMES):
            arr[..., apos] = coeff_data[f"{atom}_{suffix}"]
        fields[key] = arr
    atoms = tensor_atoms(cases[0], "matter4")
    cov0 = np.zeros_like(cases[0].r_need)
    for apos, atom in enumerate(ATOM_NAMES):
        cov0 += fields[0][..., apos, None, None] * atoms[atom]

    residual = cov0 - cases[0].r_need
    fit_mask = np.asarray(coeff_data["mask"], dtype=bool)
    active_metric = base.dilate_mask_8(fit_mask, int(args.metric_active_dilation))
    active_metric[0, :] = False
    active_metric[-1, :] = False
    active_metric[:, 0] = False
    active_metric[:, -1] = False
    variables, _ = metric_variable_index(active_metric)
    row_points = [(int(i), int(j)) for i, j in np.argwhere(fit_mask)]
    target_norm = np.asarray([max(float(tensor_norm(cases[0].r_need)[i, j]), 1.0e-300) for i, j in row_points])
    finite_target = target_norm[np.isfinite(target_norm)]
    floor = float(args.target_floor_frac) * max(
        float(np.percentile(finite_target, 95.0)) if finite_target.size else 0.0,
        1.0e-300,
    )
    row_weights = np.asarray(
        [
            np.sqrt(max(cases[0].rho_a[i, j], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
            / max(target_norm[n], floor)
            for n, (i, j) in enumerate(row_points)
        ],
        dtype=float,
    )
    rhs = np.concatenate([row_weights[n] * sym_to_vec(residual[i, j]) for n, (i, j) in enumerate(row_points)])
    dg_base, d2g_base = metric_jets_full(
        full_case.geom_m.metric_cov,
        full_case.geom_0.metric_cov,
        full_case.geom_p.metric_cov,
        dt,
        dx,
        dz,
    )
    operator = build_full_linear_operator(
        metric_cov=full_case.geom_0.metric_cov,
        dg_base=dg_base,
        d2g_base=d2g_base,
        variables=variables,
        row_points=row_points,
        row_weights=row_weights,
        dt=dt,
        dx=dx,
        dz=dz,
        eps=float(args.linear_eps),
    )
    delta_vec, solve = solve_scaled_lstsq(operator, rhs, float(args.ridge))
    delta_metric = np.zeros_like(full_case.geom_p.metric_cov)
    for value, (i, j, a, b) in zip(delta_vec, variables):
        delta_metric[i, j, a, b] += value
        if a != b:
            delta_metric[i, j, b, a] += value
    predicted = np.zeros_like(residual)
    weighted_pred = operator @ delta_vec
    for n, (i, j) in enumerate(row_points):
        predicted[i, j] = vec_to_sym(weighted_pred[n * len(COMPONENTS) : (n + 1) * len(COMPONENTS)] / max(row_weights[n], 1.0e-300))

    corrected_metric_plus = full_case.geom_p.metric_cov + delta_metric
    from fit_metric_fr_ricci2_universal_full import make_time_geometry

    geom_corr, _ = make_time_geometry(
        full_case.geom_0.metric_cov,
        full_case.geom_m.metric_cov,
        corrected_metric_plus,
        dt,
        dx,
        dz,
        with_dgamma=False,
    )
    ein_corr = geom_corr.ricci - 0.5 * full_case.geom_0.metric_cov * geom_corr.r_scalar[..., None, None]
    exact_r_need = ein_corr - cases[0].source_tilde

    before_rel = tensor_norm(residual) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    linear_after = residual - predicted
    linear_after_rel = tensor_norm(linear_after) / np.maximum(tensor_norm(cases[0].r_need + predicted), 1.0e-300)
    exact_after_rel = tensor_norm(cov0 - exact_r_need) / np.maximum(tensor_norm(exact_r_need), 1.0e-300)
    metric_norm = tensor_norm(full_case.geom_p.metric_cov)
    delta_norm = tensor_norm(delta_metric)
    weights_eval = np.sqrt(np.maximum(cases[0].rho_a[fit_mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))

    plot_path = args.output / "gbcd_full_linear_metric_update.png"
    render(plot_path, cases[0].x, cases[0].z, fit_mask, before_rel, linear_after_rel, exact_after_rel)
    npz_path = args.output / "full_linear_metric_update.npz"
    np.savez_compressed(
        npz_path,
        delta_metric_plus=delta_metric,
        corrected_metric_plus=corrected_metric_plus,
        before_rel=before_rel,
        linear_after_rel=linear_after_rel,
        exact_after_rel=exact_after_rel,
        mask=fit_mask,
    )
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "metric_active_dilation": int(args.metric_active_dilation),
            "dt_ev_inv": float(dt),
            "dx_ev_inv": float(dx),
            "dz_ev_inv": float(dz),
            "linear_eps": float(args.linear_eps),
            "ridge": float(args.ridge),
            "coefficients": str(coeff_path.resolve()),
        },
        "definition": {
            "goal": "Global dense prototype for full-linear metric update: solve one coupled linearized Einstein update for delta gtilde_+ instead of pointwise principal correction.",
            "included_terms": "Full pointwise linearization of the existing Einstein-tensor stencil with respect to gtilde_+, including local tt terms, mixed tx/tz stencil terms, and connection-quadratic linear terms.",
            "current_limit": "Dense prototype for low-resolution diagnostics only; high-resolution runs need sparse or matrix-free solve.",
        },
        "counts": {
            "fit_points": int(np.count_nonzero(fit_mask)),
            "metric_active_points": int(np.count_nonzero(active_metric)),
            "variables": int(len(variables)),
            "rows": int(operator.shape[0]),
        },
        "solve": solve,
        "evaluation": {
            "before_relative_residual": weighted_stats(before_rel[fit_mask], weights_eval),
            "linear_after_relative_residual": weighted_stats(linear_after_rel[fit_mask], weights_eval),
            "exact_after_relative_residual": weighted_stats(exact_after_rel[fit_mask], weights_eval),
            "relative_delta_metric_plus_norm": weighted_stats(
                delta_norm[active_metric] / np.maximum(metric_norm[active_metric], 1.0e-300)
            ),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "data_npz": str(npz_path.resolve()),
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
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--metric-active-dilation", type=int, default=1)
    parser.add_argument("--linear-eps", type=float, default=1.0e-8)
    parser.add_argument("--ridge", type=float, default=1.0e-10)
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
