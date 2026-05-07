from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_c_terms_from_a_reference import d1x, d1z
from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_tridomain_full_dynamics import erode_mask_8


SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


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


def pointwise_fit_tensor(case, atom_names: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    atoms = tensor_atoms(case, "matter4")
    idx = np.arange(case.r_need.shape[0] * case.r_need.shape[1])
    rows = np.stack([symmetric_rows(atoms[name], idx) for name in atom_names], axis=-1)
    target = symmetric_rows(case.r_need, idx)
    coeff = np.full((idx.size, len(atom_names)), np.nan, dtype=float)
    residual = np.full(idx.size, np.nan, dtype=float)
    fitted_rows = np.zeros_like(target)
    for p in range(idx.size):
        a = rows[p]
        b = target[p]
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        sol, *_ = np.linalg.lstsq(a, b, rcond=1.0e-12)
        fitted = a @ sol
        coeff[p] = sol
        fitted_rows[p] = fitted
        residual[p] = float(np.linalg.norm(fitted - b) / max(np.linalg.norm(b), 1.0e-300))
    shape = case.r_need.shape[:2]
    fitted_cov = np.zeros_like(case.r_need)
    for k, (a, b) in enumerate(SYMMETRIC_COMPONENTS):
        fitted_cov[..., a, b] = fitted_rows[:, k].reshape(shape)
        fitted_cov[..., b, a] = fitted_rows[:, k].reshape(shape)
    return coeff.reshape(shape + (len(atom_names),)), residual.reshape(shape), fitted_cov


def derivative_tensor(c_m: np.ndarray, c_0: np.ndarray, c_p: np.ndarray, dt: float, dx: float, dz: float) -> np.ndarray:
    dc = np.zeros(c_0.shape[:2] + (3,) + c_0.shape[2:], dtype=float)
    dc[..., 0, :, :] = (c_p - c_m) / (2.0 * dt)
    dc[..., 1, :, :] = d1x(c_0, dx)
    dc[..., 2, :, :] = d1z(c_0, dz)
    return dc


def divergence_cov2(cov: np.ndarray, dcov: np.ndarray, metric_inv: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    # J_nu = nabla^mu C_{mu nu}.
    j = np.zeros(cov.shape[:2] + (3,), dtype=float)
    for nu in range(3):
        val = np.zeros(cov.shape[:2], dtype=float)
        for mu in range(3):
            for alpha in range(3):
                term = dcov[..., alpha, mu, nu].copy()
                for lam in range(3):
                    term -= gamma[..., lam, alpha, mu] * cov[..., lam, nu]
                    term -= gamma[..., lam, alpha, nu] * cov[..., mu, lam]
                val += metric_inv[..., mu, alpha] * term
        j[..., nu] = val
    return j


def transverse_covector(j_cov: np.ndarray, u_cov: np.ndarray, metric_inv: np.ndarray) -> np.ndarray:
    u_up = np.einsum("...ab,...b->...a", metric_inv, u_cov, optimize=True)
    u2 = np.einsum("...a,...a->...", u_up, u_cov, optimize=True)
    uj = np.einsum("...a,...a->...", u_up, j_cov, optimize=True)
    return j_cov - u_cov * (uj / np.where(np.abs(u2) > 1.0e-300, u2, np.nan))[..., None]


def covector_norm_euclidean(v: np.ndarray) -> np.ndarray:
    return np.sqrt(np.sum(np.asarray(v, dtype=float) ** 2, axis=-1))


def analyze_one(args: argparse.Namespace, ref: dict[str, object], tau: float, atom_names: list[str]) -> dict[str, object]:
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

    coeff_m, residual_m, cov_m = pointwise_fit_tensor(cases[-1], atom_names)
    coeff_0, residual_0, cov_0 = pointwise_fit_tensor(cases[0], atom_names)
    coeff_p, residual_p, cov_p = pointwise_fit_tensor(cases[1], atom_names)
    _ = (coeff_m, coeff_0, coeff_p, residual_m, residual_p)

    dcov = derivative_tensor(cov_m, cov_0, cov_p, dt, dx, dz)
    j_cov = divergence_cov2(cov_0, dcov, full.geom_0.metric_inv, full.geom_0.gamma)
    j_perp = transverse_covector(j_cov, cases[0].u_cov, full.geom_0.metric_inv)

    total = covector_norm_euclidean(j_cov)
    perp = covector_norm_euclidean(j_perp)
    c_norm = tensor_norm(cov_0)
    deriv_step = min(abs(dt), abs(dx), abs(dz))
    scale_norm = c_norm / max(deriv_step, 1.0e-300)

    mask = erode_mask_8(cases[0].core10, int(args.force_mask_erosion))
    if not np.any(mask):
        mask = cases[0].core10
    return {
        "tau": float(tau),
        "atom_names": atom_names,
        "mask_count": int(np.count_nonzero(mask)),
        "algebraic_residual": stats(residual_0[mask]),
        "divergence_norm": stats(total[mask]),
        "transverse_force_norm": stats(perp[mask]),
        "transverse_fraction": stats(perp[mask] / np.maximum(total[mask], 1.0e-300)),
        "transverse_over_derivative_scale": stats(perp[mask] / np.maximum(scale_norm[mask], 1.0e-300)),
    }


def render_bar_plot(out_path: Path, report: dict[str, object]) -> None:
    cases = report["models"]
    labels = [item["model"] for item in cases]
    residual = [item["overall"]["algebraic_residual"]["p50"] for item in cases]
    force = [item["overall"]["transverse_fraction"]["p50"] for item in cases]
    scaled = [item["overall"]["transverse_over_derivative_scale"]["p50"] for item in cases]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9.8, 5.2), constrained_layout=True)
    ax.bar(x - 0.25, residual, width=0.25, label="algebraic residual p50")
    ax.bar(x, force, width=0.25, label="transverse fraction p50")
    ax.bar(x + 0.25, scaled, width=0.25, label="transverse / (|C|/h) p50")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_yscale("log")
    ax.set_title("Minimal equation candidates: residual vs geodesic-force diagnostic")
    ax.legend(fontsize=8)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def aggregate(items: list[dict[str, object]]) -> dict[str, object]:
    keys = [
        "algebraic_residual",
        "divergence_norm",
        "transverse_force_norm",
        "transverse_fraction",
        "transverse_over_derivative_scale",
    ]
    out: dict[str, object] = {}
    for key in keys:
        vals = []
        # We aggregate only summary p50/p95 across slices here to avoid
        # retaining large grids in the report.
        for item in items:
            vals.append(float(item[key]["p50"]))
        out[key] = stats(np.asarray(vals))
    return out


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    taus = [float(x) for x in args.taus.split(",") if x.strip()]
    model_specs = {
        "uu": ["uu"],
        "uu_rr": ["uu", "rr"],
        "uu_rr_ur": ["uu", "rr", "ur"],
        "g_uu_rr_ur": ["g", "uu", "rr", "ur"],
    }
    models = []
    for model, atoms in model_specs.items():
        by_tau = [analyze_one(args, ref, tau, atoms) for tau in taus]
        models.append({"model": model, "atom_names": atoms, "overall": aggregate(by_tau), "by_tau": by_tau})
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "taus": taus,
            "fit_region": "core10",
            "force_probe_dt_old": float(args.force_probe_dt_old),
            "force_mask_erosion": int(args.force_mask_erosion),
        },
        "definition": {
            "candidate": "C_mn = sum atoms with pointwise least-squares coefficients; this diagnoses the best algebraic fit before a closed coefficient law is imposed.",
            "geodesic_condition": "h^nu_alpha nabla^mu C_mu^alpha = 0. We report transverse_fraction=|J_perp|/|J| and transverse_over_derivative_scale=|J_perp|/(|C|/h).",
            "caveat": "Pointwise coefficients use independent fits at tau-probe, tau, tau+probe; this is a diagnostic, not a closed evolution law.",
        },
        "models": models,
        "outputs": {
            "summary_json": str((out / "summary.json").resolve()),
            "bar_plot": str((out / "residual_vs_transverse_force.png").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    render_bar_plot(out / "residual_vs_transverse_force.png", report)
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--full-resolution", type=int, default=160)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
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
