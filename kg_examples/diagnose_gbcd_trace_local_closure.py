from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import SYMMETRIC_COMPONENTS, build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


ATOM_NAMES = ("g", "uu", "rr", "ur")


def weighted_stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
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
    if w is not None and float(np.sum(w)) > 0.0:
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


def trace_coefficients(case) -> np.ndarray:
    atoms = tensor_atoms(case, "matter4")
    out = np.zeros(case.rho_a.shape + (len(ATOM_NAMES),), dtype=float)
    for apos, atom_name in enumerate(ATOM_NAMES):
        out[..., apos] = np.einsum("...ab,...ab->...", case.metric_inv, atoms[atom_name], optimize=True)
    return out


def qe_trace_coefficients(case) -> np.ndarray:
    atoms = tensor_atoms(case, "matter4")
    g_inv = case.metric_inv
    u = case.u_cov
    r = case.r_cov
    u_up = np.einsum("...ab,...b->...a", g_inv, u, optimize=True)
    r_up = np.einsum("...ab,...b->...a", g_inv, r, optimize=True)
    u2 = np.einsum("...a,...a->...", u_up, u, optimize=True)
    ur = np.einsum("...a,...a->...", u_up, r, optimize=True)
    safe_u2 = np.where(np.abs(u2) > 1.0e-300, u2, np.nan)
    r_perp = r - (ur / safe_u2)[..., None] * u
    r_perp_up = np.einsum("...ab,...b->...a", g_inv, r_perp, optimize=True)
    r_perp2 = np.einsum("...a,...a->...", r_perp_up, r_perp, optimize=True)

    e0_cov = u / np.sqrt(np.maximum(u2, 1.0e-300))[..., None]
    e1_cov = r_perp / np.sqrt(np.maximum(-r_perp2, 1.0e-300))[..., None]
    e0_up = np.einsum("...ab,...b->...a", g_inv, e0_cov, optimize=True)
    e1_up = np.einsum("...ab,...b->...a", g_inv, e1_cov, optimize=True)
    qe = (
        np.einsum("...a,...b->...ab", e0_up, e0_up, optimize=True)
        + np.einsum("...a,...b->...ab", e1_up, e1_up, optimize=True)
    )

    out = np.zeros(case.rho_a.shape + (len(ATOM_NAMES),), dtype=float)
    for apos, atom_name in enumerate(ATOM_NAMES):
        out[..., apos] = np.einsum("...ab,...ab->...", qe, atoms[atom_name], optimize=True)
    bad = (~np.isfinite(u2)) | (~np.isfinite(r_perp2)) | (u2 <= 0.0) | (r_perp2 >= 0.0)
    out[bad, :] = np.nan
    return out


def nullspace_basis(row: np.ndarray) -> np.ndarray | None:
    row = np.asarray(row, dtype=float)
    if not np.all(np.isfinite(row)):
        return None
    norm = float(np.linalg.norm(row))
    if norm <= 1.0e-300:
        return None
    _, s, vh = np.linalg.svd(row.reshape(1, -1), full_matrices=True)
    tol = row.size * np.finfo(float).eps * (float(s[0]) if s.size else 0.0)
    rank = int(np.count_nonzero(s > tol))
    return vh[rank:].T


def solve_local(a: np.ndarray, b: np.ndarray, constraint: np.ndarray | None) -> tuple[np.ndarray, float]:
    if constraint is None:
        coeff = np.linalg.lstsq(a, b, rcond=1.0e-12)[0]
    else:
        n = nullspace_basis(constraint)
        if n is None or n.shape[1] == 0:
            return np.full(a.shape[1], np.nan), np.nan
        z = np.linalg.lstsq(a @ n, b, rcond=1.0e-12)[0]
        coeff = n @ z
    residual = float(np.linalg.norm(a @ coeff - b))
    return coeff, residual


def evaluate_case(case, region: str, target_floor_frac: float) -> dict[str, object]:
    mask = getattr(case, region)
    idx = np.where(mask.reshape(-1))[0]
    atoms = tensor_atoms(case, "matter4")
    atom_rows = np.stack([symmetric_rows(atoms[name], idx) for name in ATOM_NAMES], axis=-1)
    target = symmetric_rows(case.r_need, idx)
    target_norm = np.sqrt(np.sum(target**2, axis=1))
    finite_norm = target_norm[np.isfinite(target_norm)]
    floor = float(target_floor_frac) * max(
        float(np.percentile(finite_norm, 95.0)) if finite_norm.size else 0.0,
        1.0e-300,
    )
    denom = np.maximum(target_norm, floor)
    weights = np.sqrt(
        np.maximum(case.rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300)
    )
    trace_rows = trace_coefficients(case).reshape(-1, len(ATOM_NAMES))[idx]
    qe_rows = qe_trace_coefficients(case).reshape(-1, len(ATOM_NAMES))[idx]

    residuals: dict[str, list[float]] = {"unconstrained": [], "trace0": [], "qe0": []}
    coeff_norms: dict[str, list[float]] = {"unconstrained": [], "trace0": [], "qe0": []}
    valid_weights: dict[str, list[float]] = {"unconstrained": [], "trace0": [], "qe0": []}
    skipped = {"trace0": 0, "qe0": 0}

    for p in range(idx.size):
        a = atom_rows[p]
        b = target[p]
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        for name, constraint in (
            ("unconstrained", None),
            ("trace0", trace_rows[p]),
            ("qe0", qe_rows[p]),
        ):
            coeff, res = solve_local(a, b, constraint)
            if not np.isfinite(res):
                if name in skipped:
                    skipped[name] += 1
                continue
            residuals[name].append(res / denom[p])
            coeff_norms[name].append(float(np.linalg.norm(coeff)))
            valid_weights[name].append(float(weights[p]))

    by_closure = {}
    for name in residuals:
        r = np.asarray(residuals[name], dtype=float)
        w = np.asarray(valid_weights[name], dtype=float)
        by_closure[name] = {
            "relative_residual": weighted_stats(r, w),
            "coeff_norm": weighted_stats(np.asarray(coeff_norms[name], dtype=float), w),
        }
    return {
        "tau": float(case.tau_old),
        "region": region,
        "mask_count": int(np.count_nonzero(mask)),
        "target_norm": weighted_stats(target_norm, weights),
        "closures": by_closure,
        "skipped": skipped,
        "fields": {
            "rho_a": weighted_stats(case.rho_a[mask]),
            "r_need_norm": weighted_stats(tensor_norm(case.r_need)[mask], weights),
        },
    }


def render_summary(path: Path, results: list[dict[str, object]]) -> None:
    closures = ["unconstrained", "trace0", "qe0"]
    labels = [f"tau={r['tau']:g}" for r in results]
    x = np.arange(len(labels))
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.5, 4.2), constrained_layout=True)
    for k, name in enumerate(closures):
        vals = [
            float(r["closures"][name]["relative_residual"]["weighted_mean"])  # type: ignore[index]
            for r in results
        ]
        ax.bar(x + (k - 1) * width, vals, width=width, label=name)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_ylabel("weighted mean local relative residual")
    ax.set_title("n=384 local algebraic closure check")
    ax.legend()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    taus = [float(x.strip()) for x in args.taus.split(",") if x.strip()]
    results = []
    for tau in taus:
        print(f"[local-closure] building case tau={tau:g}", flush=True)
        case = build_case(args, ref, tau)
        print(f"[local-closure] evaluating tau={tau:g} region={args.fit_region}", flush=True)
        results.append(evaluate_case(case, args.fit_region, float(args.target_floor_frac)))
    plot_path = args.output / "local_trace_closure_summary.png"
    render_summary(plot_path, results)
    report = {
        "parameters": {
            "taus": taus,
            "fit_region": args.fit_region,
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "target_floor_frac": float(args.target_floor_frac),
        },
        "definition": {
            "purpose": "Lightweight n=384 local algebraic sanity check. It does not impose the full conservation rows or time-coupled hard nullspace solve.",
            "unconstrained": "Pointwise least-squares fit of R_need in span{g,uu,rr,ur}.",
            "trace0": "Pointwise least-squares fit in span{g,uu,rr,ur} with g^{mn} C_mn=0.",
            "qe0": "Pointwise least-squares fit with q_E^{mn} C_mn=0, where q_E is the positive u-r plane contraction.",
            "residual": "||A lambda - R_need|| / max(||R_need||, target_floor_frac*p95(||R_need||)). Weighted means use sqrt(rho_a/max rho_a).",
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
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
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
