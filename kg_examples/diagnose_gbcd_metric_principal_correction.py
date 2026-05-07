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
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from test_gbcd_v0_one_step_lambda import ATOM_NAMES, load_full_coeff, weighted_stats


COMPONENTS = base.SYMMETRIC_COMPONENTS


def sym_to_vec(tensor: np.ndarray) -> np.ndarray:
    return np.asarray([tensor[a, b] for a, b in COMPONENTS], dtype=float)


def vec_to_sym(vec: np.ndarray) -> np.ndarray:
    out = np.zeros((3, 3), dtype=float)
    for value, (a, b) in zip(vec, COMPONENTS):
        out[a, b] = float(value)
        out[b, a] = float(value)
    return out


def principal_delta_einstein(metric_cov: np.ndarray, metric_inv: np.ndarray, hdd: np.ndarray) -> np.ndarray:
    """Time-time principal part of delta G_mn from hdd_mn = partial_t^2 h_mn.

    This keeps only terms containing two lab-time derivatives.  In the current
    1550 nm setup dt << dx,dz, so this is the dominant local part of the
    second-order metric evolution operator.
    """
    out_ric = np.zeros((3, 3), dtype=float)
    trace_hdd = float(np.einsum("ab,ab->", metric_inv, hdd))
    for mu in range(3):
        for nu in range(3):
            term = -metric_inv[0, 0] * hdd[mu, nu]
            if mu == 0:
                term += float(np.einsum("b,b->", metric_inv[0], hdd[:, nu]))
            if nu == 0:
                term += float(np.einsum("b,b->", metric_inv[0], hdd[:, mu]))
            if mu == 0 and nu == 0:
                term -= trace_hdd
            out_ric[mu, nu] = 0.5 * term
    delta_r = float(np.einsum("ab,ab->", metric_inv, out_ric))
    return out_ric - 0.5 * metric_cov * delta_r


def local_principal_matrix(metric_cov: np.ndarray, metric_inv: np.ndarray, dt: float) -> np.ndarray:
    mat = np.zeros((len(COMPONENTS), len(COMPONENTS)), dtype=float)
    for col, (a, b) in enumerate(COMPONENTS):
        h = np.zeros((3, 3), dtype=float)
        h[a, b] = 1.0 / (dt * dt)
        h[b, a] = 1.0 / (dt * dt)
        dg = principal_delta_einstein(metric_cov, metric_inv, h)
        mat[:, col] = sym_to_vec(dg)
    return mat


def render(path: Path, x: np.ndarray, z: np.ndarray, mask: np.ndarray, before: np.ndarray, after: np.ndarray, hnorm: np.ndarray) -> None:
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), constrained_layout=True)
    fields = [
        ("before residual", before, "magma"),
        ("after principal metric correction", after, "magma"),
        (r"$|\delta \tilde g_+|$ principal estimate", hnorm, "viridis"),
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
    dt = float(args.dt_old) * old_scale
    probe = float(args.force_probe_dt_old)
    cases = {key: build_case(args, ref, float(args.tau) + key * probe) for key in base.TIME_KEYS}
    full_case = build_full_case(args, ref, float(args.tau))
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    _, _, meta = base.assemble_system(
        cases,
        full_case,
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
    coeff_path = args.coefficients
    if coeff_path is None:
        if abs(float(args.tau)) < 1.0e-12:
            coeff_path = Path("visualizations/equation_first_gbcd_aux_noq_time_scan_n96_tau0_tw_1em5/gbcd_auxiliary_gauge_coefficients.npz")
        elif float(args.tau) < 0:
            coeff_path = Path("visualizations/equation_first_gbcd_aux_noq_time1em5_n96_taum3p5/gbcd_auxiliary_gauge_coefficients.npz")
        else:
            coeff_path = Path("visualizations/equation_first_gbcd_aux_noq_time1em5_n96_taup3p5/gbcd_auxiliary_gauge_coefficients.npz")
    coeff = load_full_coeff(coeff_path, meta, "time", cases[0].rho_a.shape)
    fields = base.coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
    cov = base.coeff_fields_to_cov(cases, fields)

    residual = cov[0] - cases[0].r_need
    delta_g_plus = np.zeros_like(cases[0].metric_cov)
    predicted_delta_g = np.zeros_like(residual)
    ranks = []
    conds = []
    mask = getattr(cases[0], args.fit_region)
    for i, j in np.argwhere(mask):
        i = int(i)
        j = int(j)
        mat = local_principal_matrix(cases[0].metric_cov[i, j], cases[0].metric_inv[i, j], dt)
        rhs = sym_to_vec(residual[i, j])
        sol, *_ = np.linalg.lstsq(mat, rhs, rcond=1.0e-12)
        delta_g_plus[i, j] = vec_to_sym(sol)
        predicted_delta_g[i, j] = vec_to_sym(mat @ sol)
        singular = np.linalg.svd(mat, compute_uv=False)
        ranks.append(int(np.count_nonzero(singular > max(mat.shape) * np.finfo(float).eps * singular[0])) if singular.size else 0)
        conds.append(float(singular[0] / max(singular[-1], 1.0e-300)) if singular.size else 0.0)

    before_rel = tensor_norm(residual) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    after_residual = residual - predicted_delta_g
    after_rel = tensor_norm(after_residual) / np.maximum(tensor_norm(cases[0].r_need + predicted_delta_g), 1.0e-300)
    hnorm = tensor_norm(delta_g_plus)
    metric_norm = tensor_norm(cases[0].metric_cov)
    weights = np.sqrt(np.maximum(cases[0].rho_a[mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    plot_path = args.output / "gbcd_metric_principal_correction.png"
    render(plot_path, cases[0].x, cases[0].z, mask, before_rel, after_rel, hnorm)
    data_path = args.output / "metric_principal_correction.npz"
    np.savez_compressed(
        data_path,
        delta_g_plus=delta_g_plus,
        predicted_delta_g=predicted_delta_g,
        before_rel=before_rel,
        after_rel=after_rel,
        mask=mask,
    )
    report = {
        "parameters": {
            "tau": float(args.tau),
            "full_resolution": int(args.full_resolution),
            "fit_region": args.fit_region,
            "dt_ev_inv": float(dt),
            "dx_ev_inv": float(dx),
            "dz_ev_inv": float(dz),
            "dt_over_dx": float(abs(dt / dx)),
            "dt_over_dz": float(abs(dt / dz)),
            "coefficients": str(coeff_path.resolve()),
        },
        "definition": {
            "goal": "Pointwise principal-part estimate of the next metric slice correction delta gtilde_+ needed to satisfy the center algebraic gBCD field equation.",
            "kept_terms": "Only terms with two lab-time derivatives in delta G_mn are kept.",
            "why_valid_as_first_estimate": "In this setup dt/dx and dt/dz are about 1e-6, so mixed/spatial derivative terms are much smaller in the explicit time-step operator.",
        },
        "matrix": {
            "rank_counts": {str(k): int(ranks.count(k)) for k in sorted(set(ranks))},
            "condition_p50": float(np.percentile(conds, 50.0)) if conds else 0.0,
            "condition_p95": float(np.percentile(conds, 95.0)) if conds else 0.0,
            "condition_max": float(np.max(conds)) if conds else 0.0,
        },
        "evaluation": {
            "before_relative_residual": weighted_stats(before_rel[mask], weights),
            "after_relative_residual": weighted_stats(after_rel[mask], weights),
            "delta_metric_plus_norm": weighted_stats(hnorm[mask], weights),
            "relative_delta_metric_plus_norm": weighted_stats(hnorm[mask] / np.maximum(metric_norm[mask], 1.0e-300), weights),
        },
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot_path.resolve()),
            "data_npz": str(data_path.resolve()),
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
