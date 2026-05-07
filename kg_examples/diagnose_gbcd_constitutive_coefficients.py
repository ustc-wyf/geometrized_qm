from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import build_case, lower_mixed_square
from physical_units import PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


COEFF_NAMES = ("g", "uu", "rr", "ur")


def stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(np.abs(vals))),
        "p50": float(np.percentile(np.abs(vals), 50.0)),
        "p95": float(np.percentile(np.abs(vals), 95.0)),
        "max": float(np.max(np.abs(vals))),
    }
    if weights is not None:
        w = np.asarray(weights, dtype=float)[mask]
        denom = float(np.sum(np.maximum(w, 0.0)))
        if denom > 0.0:
            out["weighted_mean"] = float(np.sum(np.maximum(w, 0.0) * np.abs(vals)) / denom)
    return out


def coefficient_path(root: Path, tau: float, mode: str) -> Path:
    if tau < 0:
        tag = f"m{abs(tau):g}".replace(".", "p")
    elif tau > 0:
        tag = f"p{tau:g}".replace(".", "p")
    else:
        tag = "0"
    if mode == "full":
        if tag == "0":
            return root / "equation_first_constrained_gbcd_hard_full_n96_tau0_trustedfit" / "constrained_bcd_best_coefficients.npz"
        return root / f"equation_first_constrained_gbcd_hard_full_n96_tau_{tag}_trustedfit" / "constrained_bcd_best_coefficients.npz"
    if tag == "0":
        return root / "equation_first_constrained_gbcd_hard_n96_tau0_trustedfit_cli_recheck" / "constrained_bcd_best_coefficients.npz"
    return root / f"equation_first_constrained_gbcd_hard_n96_tau_{tag}_trustedfit_cli" / "constrained_bcd_best_coefficients.npz"


def coefficient_tag(tau: float) -> str:
    if tau < 0:
        return f"m{abs(tau):g}".replace(".", "p")
    if tau > 0:
        return f"p{tau:g}".replace(".", "p")
    return "0"


def invariant_features(case) -> dict[str, np.ndarray]:
    u = case.u_cov
    r = case.r_cov
    ginv = case.metric_inv
    ricci_sq = lower_mixed_square(case.ricci, ginv)
    i2 = np.einsum("...ab,...ab->...", ginv, ricci_sq, optimize=True)
    u2 = np.einsum("...a,...ab,...b->...", u, ginv, u, optimize=True)
    r2 = np.einsum("...a,...ab,...b->...", r, ginv, r, optimize=True)
    ur = np.einsum("...a,...ab,...b->...", u, ginv, r, optimize=True)
    rho_rel = case.rho_a / max(float(np.max(case.rho_a)), 1.0e-300)
    return {
        "log_rho_rel": np.log10(np.maximum(rho_rel, 1.0e-300)),
        "R": case.r_scalar,
        "I2": i2,
        "u2": u2,
        "r2": r2,
        "ur": ur,
    }


def make_basis(features: np.ndarray, degree: int) -> tuple[np.ndarray, list[str]]:
    n, k = features.shape
    cols = [np.ones(n, dtype=float)]
    names = ["1"]
    base_names = [f"f{i}" for i in range(k)]
    for i in range(k):
        cols.append(features[:, i])
        names.append(base_names[i])
    if degree >= 2:
        for i in range(k):
            for j in range(i, k):
                cols.append(features[:, i] * features[:, j])
                names.append(f"{base_names[i]}*{base_names[j]}")
    return np.stack(cols, axis=1), names


def weighted_lstsq_predict(
    train_x: np.ndarray,
    train_y: np.ndarray,
    train_w: np.ndarray,
    test_x: np.ndarray,
    degree: int,
) -> np.ndarray:
    a, _ = make_basis(train_x, degree)
    b, _ = make_basis(test_x, degree)
    scale = np.linalg.norm(a * train_w[:, None], axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    coeff, *_ = np.linalg.lstsq((a / scale[None, :]) * train_w[:, None], train_y * train_w, rcond=1.0e-12)
    coeff = coeff / scale
    return b @ coeff


def knn_predict(train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray, k: int) -> np.ndarray:
    out = np.zeros(test_x.shape[0], dtype=float)
    kk = max(1, min(int(k), train_x.shape[0]))
    for start in range(0, test_x.shape[0], 256):
        chunk = test_x[start : start + 256]
        dist2 = np.sum((chunk[:, None, :] - train_x[None, :, :]) ** 2, axis=2)
        idx = np.argpartition(dist2, kk - 1, axis=1)[:, :kk]
        w = 1.0 / np.maximum(np.take_along_axis(dist2, idx, axis=1), 1.0e-30)
        vals = train_y[idx]
        out[start : start + chunk.shape[0]] = np.sum(w * vals, axis=1) / np.sum(w, axis=1)
    return out


def relative_error(pred: np.ndarray, target: np.ndarray, weights: np.ndarray, floor_frac: float) -> dict[str, float]:
    scale = max(float(np.percentile(np.abs(target[np.isfinite(target)]), 95.0)), 1.0e-300)
    denom = np.maximum(np.abs(target), float(floor_frac) * scale)
    rel = np.abs(pred - target) / denom
    return stats(rel, weights)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    taus = [float(x) for x in args.taus.split(",") if x.strip()]
    cases = {tau: build_case(args, ref, tau) for tau in taus}

    rows = []
    for tau in taus:
        if args.coefficient_template:
            coeff_file = Path(args.coefficient_template.format(tag=coefficient_tag(tau), tau=f"{tau:g}"))
        else:
            coeff_file = coefficient_path(args.visualization_root, tau, args.force_mode)
        data = np.load(coeff_file)
        case = cases[tau]
        mask = getattr(case, args.region)
        feats = invariant_features(case)
        feature_names = list(feats.keys())
        x = np.stack([feats[name][mask] for name in feature_names], axis=1)
        weights = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
        coeffs = {name: np.asarray(data[f"{name}_0"], dtype=float)[mask] for name in COEFF_NAMES}
        rows.append({"tau": tau, "x": x, "weights": weights, "coeffs": coeffs})

    all_x = np.vstack([row["x"] for row in rows])
    feature_scale = np.percentile(np.abs(all_x), 95.0, axis=0)
    feature_scale = np.where(feature_scale > 0.0, feature_scale, 1.0)
    for row in rows:
        row["x_scaled"] = row["x"] / feature_scale[None, :]

    reports = []
    for coeff_name in COEFF_NAMES:
        y_all = np.concatenate([row["coeffs"][coeff_name] for row in rows])
        w_all = np.concatenate([row["weights"] for row in rows])
        x_all = np.vstack([row["x_scaled"] for row in rows])
        in_sample = {}
        for degree in (1, 2):
            pred = weighted_lstsq_predict(x_all, y_all, w_all, x_all, degree)
            in_sample[f"poly_degree_{degree}"] = relative_error(pred, y_all, w_all, args.coeff_floor_frac)
        by_holdout = {}
        for test_i, test_row in enumerate(rows):
            train = [row for i, row in enumerate(rows) if i != test_i]
            tx = np.vstack([row["x_scaled"] for row in train])
            ty = np.concatenate([row["coeffs"][coeff_name] for row in train])
            tw = np.concatenate([row["weights"] for row in train])
            test_x = test_row["x_scaled"]
            test_y = test_row["coeffs"][coeff_name]
            test_w = test_row["weights"]
            by_holdout[f"tau={test_row['tau']:g}"] = {
                "poly_degree_1": relative_error(
                    weighted_lstsq_predict(tx, ty, tw, test_x, 1),
                    test_y,
                    test_w,
                    args.coeff_floor_frac,
                ),
                "poly_degree_2": relative_error(
                    weighted_lstsq_predict(tx, ty, tw, test_x, 2),
                    test_y,
                    test_w,
                    args.coeff_floor_frac,
                ),
                "knn_k5": relative_error(knn_predict(tx, ty, test_x, 5), test_y, test_w, args.coeff_floor_frac),
            }
        reports.append({"coefficient": coeff_name, "in_sample": in_sample, "leave_one_tau_out": by_holdout})

    plot_path = args.output / "gbcd_constitutive_coefficients_residuals.png"
    labels = [r["coefficient"] for r in reports]
    in_poly2 = [r["in_sample"]["poly_degree_2"].get("weighted_mean", 0.0) for r in reports]
    hold_poly2 = []
    hold_knn = []
    for r in reports:
        vals_p = [v["poly_degree_2"].get("weighted_mean", 0.0) for v in r["leave_one_tau_out"].values()]
        vals_k = [v["knn_k5"].get("weighted_mean", 0.0) for v in r["leave_one_tau_out"].values()]
        hold_poly2.append(float(np.mean(vals_p)))
        hold_knn.append(float(np.mean(vals_k)))
    x_axis = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.bar(x_axis - 0.25, in_poly2, width=0.25, label="in-sample quadratic")
    ax.bar(x_axis, hold_poly2, width=0.25, label="leave-one-tau quadratic")
    ax.bar(x_axis + 0.25, hold_knn, width=0.25, label="leave-one-tau kNN")
    ax.set_xticks(x_axis)
    ax.set_xticklabels(labels)
    ax.set_yscale("log")
    ax.set_title("Can gBCD coefficients be simple local scalar functions?")
    ax.set_ylabel("weighted relative error")
    ax.legend(fontsize=8)
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)

    report = {
        "parameters": {
            "taus": taus,
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "region": args.region,
            "force_mode": args.force_mode,
            "feature_names": feature_names,
            "feature_scale_p95": {name: float(scale) for name, scale in zip(feature_names, feature_scale)},
            "coefficient_template": args.coefficient_template,
            "coeff_floor_frac": float(args.coeff_floor_frac),
        },
        "definition": {
            "goal": "Test whether hard-constrained gBCD coefficient fields A,B,C,D can be represented by simple local scalar functions.",
            "caveat": "The KKT coefficients can contain null-space/gauge choices; failure here means no simple constitutive law in this feature set, not a final no-go for gBCD.",
            "features": "log rho_rel, Rtilde, I2=RmnR^mn, u2, r2, ur, scaled by global p95.",
        },
        "reports": reports,
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
    parser.add_argument("--visualization-root", type=Path, default=Path("visualizations"))
    parser.add_argument("--coefficient-template", type=str, default="")
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--force-mode", choices=["transverse", "full"], default="full")
    parser.add_argument("--region", choices=["support", "trusted", "core10"], default="trusted")
    parser.add_argument("--coeff-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--mp", type=float, default=None)
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
