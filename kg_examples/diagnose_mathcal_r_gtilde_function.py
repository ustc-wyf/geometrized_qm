from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import (
    CaseData,
    build_case,
    symmetric_rows,
    tensor_norm,
)
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


def stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    result = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }
    if weights is not None:
        w = np.asarray(weights, dtype=float)[mask]
        w = np.maximum(w, 0.0)
        if float(np.sum(w)) > 0.0:
            result["weighted_mean"] = float(np.sum(w * vals) / np.sum(w))
    return result


def ricci_square_and_i2(case: CaseData) -> tuple[np.ndarray, np.ndarray]:
    mixed = np.einsum("...ac,...cb->...ab", case.metric_inv, case.ricci, optimize=True)
    ricci_sq = np.einsum("...ac,...cb->...ab", case.ricci, mixed, optimize=True)
    i2 = np.einsum("...ab,...ab->...", case.metric_inv, ricci_sq, optimize=True)
    return ricci_sq, i2


def ricci_i3(case: CaseData) -> np.ndarray:
    mixed = np.einsum("...ac,...cb->...ab", case.metric_inv, case.ricci, optimize=True)
    mixed2 = np.einsum("...ac,...cb->...ab", mixed, mixed, optimize=True)
    mixed3 = np.einsum("...ac,...cb->...ab", mixed2, mixed, optimize=True)
    return np.trace(mixed3, axis1=-2, axis2=-1)


def feature_matrix(case: CaseData, feature_set: str, region: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mask = getattr(case, region)
    idx = np.where(mask.reshape(-1))[0]
    ricci_sq, i2 = ricci_square_and_i2(case)
    i3 = ricci_i3(case)

    pieces: list[np.ndarray] = []
    if feature_set == "curvature_scalars":
        scalars = np.stack(
            [
                case.r_scalar.reshape(-1)[idx],
                i2.reshape(-1)[idx],
                i3.reshape(-1)[idx],
                np.log10(1.0e-300 + np.abs(np.linalg.det(case.metric_cov.reshape(-1, 3, 3))[idx])),
            ],
            axis=1,
        )
        pieces.append(scalars)
    elif feature_set in {"local_geometry_tensors", "local_geometry_plus_matter"}:
        tensors = [
            case.metric_cov,
            case.ricci,
            case.einstein,
            ricci_sq,
            i2[..., None, None] * case.metric_cov,
        ]
        pieces.extend(symmetric_rows(tensor, idx) for tensor in tensors)
        scalars = np.stack([case.r_scalar.reshape(-1)[idx], i2.reshape(-1)[idx], i3.reshape(-1)[idx]], axis=1)
        pieces.append(scalars)
        if feature_set == "local_geometry_plus_matter":
            u_outer = np.einsum("...a,...b->...ab", case.u_cov, case.u_cov, optimize=True)
            r_outer = np.einsum("...a,...b->...ab", case.r_cov, case.r_cov, optimize=True)
            ur_sym = 0.5 * (
                np.einsum("...a,...b->...ab", case.u_cov, case.r_cov, optimize=True)
                + np.einsum("...a,...b->...ab", case.r_cov, case.u_cov, optimize=True)
            )
            pieces.extend(symmetric_rows(tensor, idx) for tensor in [u_outer, r_outer, ur_sym])
    else:
        raise ValueError(f"unknown feature_set: {feature_set}")

    features = np.concatenate(pieces, axis=1)
    # Signed log compression preserves signs while preventing a few curvature
    # spikes from defining the whole nearest-neighbour geometry.
    features = np.sign(features) * np.log1p(np.abs(features))
    target = symmetric_rows(case.source_tilde, idx)
    weights = np.sqrt(np.maximum(case.rho_a.reshape(-1)[idx], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
    return features, target, weights


def knn_predict_leave_one_time(
    cases: list[CaseData],
    *,
    feature_set: str,
    region: str,
    k: int,
    chunk_size: int,
) -> dict[str, object]:
    per_case: dict[str, object] = {}
    all_rel: list[np.ndarray] = []
    all_weights: list[np.ndarray] = []
    for held_index, held in enumerate(cases):
        train_features: list[np.ndarray] = []
        train_target: list[np.ndarray] = []
        for index, case in enumerate(cases):
            if index == held_index:
                continue
            feat, target, _ = feature_matrix(case, feature_set, region)
            train_features.append(feat)
            train_target.append(target)
        x_train = np.vstack(train_features)
        y_train = np.vstack(train_target)
        x_test, y_test, weights = feature_matrix(held, feature_set, region)

        mean = np.mean(x_train, axis=0)
        std = np.std(x_train, axis=0)
        std = np.where(std > 1.0e-12, std, 1.0)
        xt = (x_train - mean) / std
        xs = (x_test - mean) / std

        pred = np.zeros_like(y_test)
        nn_dist = np.zeros(y_test.shape[0], dtype=float)
        k_eff = min(int(k), y_train.shape[0])
        for start in range(0, xs.shape[0], int(chunk_size)):
            stop = min(start + int(chunk_size), xs.shape[0])
            block = xs[start:stop]
            dist2 = np.sum((block[:, None, :] - xt[None, :, :]) ** 2, axis=2)
            idx = np.argpartition(dist2, kth=k_eff - 1, axis=1)[:, :k_eff]
            dsel = np.take_along_axis(dist2, idx, axis=1)
            # Inverse-distance weights; exact matches become dominant but finite.
            w = 1.0 / np.maximum(dsel, 1.0e-18)
            w = w / np.sum(w, axis=1, keepdims=True)
            pred[start:stop] = np.einsum("bk,bkc->bc", w, y_train[idx], optimize=True)
            nn_dist[start:stop] = np.sqrt(np.min(dsel, axis=1))

        rel = np.linalg.norm(pred - y_test, axis=1) / np.maximum(np.linalg.norm(y_test, axis=1), 1.0e-300)
        all_rel.append(rel)
        all_weights.append(weights)
        per_case[f"tau={held.tau_old:g}"] = {
            "count": int(y_test.shape[0]),
            "relative_residual": stats(rel, weights),
            "nearest_distance": stats(nn_dist, weights),
            "target_norm": stats(np.linalg.norm(y_test, axis=1), weights),
            "prediction_norm": stats(np.linalg.norm(pred, axis=1), weights),
        }
    rel_all = np.concatenate(all_rel)
    weights_all = np.concatenate(all_weights)
    return {
        "feature_set": feature_set,
        "region": region,
        "k": int(k),
        "overall": stats(rel_all, weights_all),
        "by_heldout_time": per_case,
    }


def render_knn_bars(out_path: Path, diagnostics: list[dict[str, object]]) -> None:
    labels = [str(item["feature_set"]) for item in diagnostics]
    p50 = [float(item["overall"]["p50"]) for item in diagnostics]
    p95 = [float(item["overall"]["p95"]) for item in diagnostics]
    wmean = [float(item["overall"].get("weighted_mean", item["overall"]["mean"])) for item in diagnostics]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    ax.bar(x - 0.25, p50, width=0.25, label="p50")
    ax.bar(x, wmean, width=0.25, label="rho-weighted mean")
    ax.bar(x + 0.25, p95, width=0.25, label="p95")
    ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("leave-one-time kNN relative residual for Ttilde/Mp^2")
    ax.set_title("Can the nontrivial part be predicted from local gtilde features?")
    ax.legend()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    cases = [build_case(args, ref, tau) for tau in [float(x) for x in args.taus.split(",") if x.strip()]]
    feature_sets = [item.strip() for item in args.feature_sets.split(",") if item.strip()]
    diagnostics = [
        knn_predict_leave_one_time(
            cases,
            feature_set=feature_set,
            region=args.region,
            k=args.k,
            chunk_size=args.chunk_size,
        )
        for feature_set in feature_sets
    ]
    plot = args.output / "mathcal_r_gtilde_function_knn.png"
    render_knn_bars(plot, diagnostics)
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "full_resolution": int(args.full_resolution),
            "cropped_shape": [int(cases[0].rho_a.shape[0]), int(cases[0].rho_a.shape[1])],
            "window_um": float(args.window_um),
            "dt_old": float(args.dt_old),
            "taus": [case.tau_old for case in cases],
            "region": args.region,
            "k": int(args.k),
            "mp_ev": float(args.mp),
        },
        "definition": {
            "target": "The nontrivial part source_tilde = Ttilde_mn/Mp^2. Full R_need is not used as the target because it is trivially dominated by Gtilde_mn.",
            "method": "Leave-one-time-slice k-nearest-neighbour prediction. Train on two tau slices, predict the held-out tau slice from local gtilde-derived features. Low residual would support single-valued local dependence on gtilde features.",
            "curvature_scalars": "R, tr(R^2), tr(R^3), log|det(gtilde)|.",
            "local_geometry_tensors": "Symmetric components of gtilde, Ricci, Einstein, Ricci^2, I2*gtilde, plus R/I2/I3.",
            "local_geometry_plus_matter": "local_geometry_tensors plus u_mu u_nu, r_mu r_nu, u_(mu r_nu); this is a contrast, not pure gtilde.",
        },
        "diagnostics": diagnostics,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "plot_png": str(plot.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument(
        "--feature-sets",
        type=str,
        default="curvature_scalars,local_geometry_tensors,local_geometry_plus_matter",
    )
    parser.add_argument("--k", type=int, default=16)
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--full-resolution", type=int, default=320)
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
    run(parser.parse_args())


if __name__ == "__main__":
    main()
