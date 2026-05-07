from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from analyze_bcd_residuals_from_a_reference import stress_tensor_tilde
from coordinate_matter_evolution import inverse_metric_block
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_tridomain_full_dynamics import safe_sqrt_abs_det


def load_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


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


def bilinear_sample(field: np.ndarray, x: np.ndarray, z: np.ndarray, xp: np.ndarray, zp: np.ndarray) -> np.ndarray:
    ix = np.searchsorted(x, xp) - 1
    iz = np.searchsorted(z, zp) - 1
    ix = np.clip(ix, 0, len(x) - 2)
    iz = np.clip(iz, 0, len(z) - 2)
    tx = (xp - x[ix]) / np.maximum(x[ix + 1] - x[ix], 1.0e-300)
    tz = (zp - z[iz]) / np.maximum(z[iz + 1] - z[iz], 1.0e-300)
    f00 = field[ix, iz]
    f10 = field[ix + 1, iz]
    f01 = field[ix, iz + 1]
    f11 = field[ix + 1, iz + 1]
    return (1.0 - tx) * (1.0 - tz) * f00 + tx * (1.0 - tz) * f10 + (1.0 - tx) * tz * f01 + tx * tz * f11


def weighted_rel_l2(reference: np.ndarray, candidate: np.ndarray, weight: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return 0.0
    rel = (candidate[valid] - reference[valid]) / np.maximum(np.abs(reference[valid]), 1.0e-300)
    return float(np.sqrt(np.sum(weight[valid] * rel * rel) / np.sum(weight[valid])))


def weighted_rel_l1(reference: np.ndarray, candidate: np.ndarray, weight: np.ndarray, mask: np.ndarray) -> float:
    valid = mask & np.isfinite(reference) & np.isfinite(candidate) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(valid):
        return 0.0
    rel = np.abs(candidate[valid] - reference[valid]) / np.maximum(np.abs(reference[valid]), 1.0e-300)
    return float(np.sum(weight[valid] * rel) / np.sum(weight[valid]))


def row_segments_um(rows: list[dict[str, float]], max_segments: int | None = None) -> list[list[tuple[float, float]]]:
    selected = rows
    if max_segments is not None and len(rows) > max_segments:
        idx = np.linspace(0, len(rows) - 1, max_segments, dtype=int)
        selected = [rows[int(i)] for i in idx]
    scale = HBAR_C_EV_M * 1.0e6
    return [[(r["x0"] * scale, r["z0"] * scale), (r["x1"] * scale, r["z1"] * scale)] for r in selected]


def run(args: argparse.Namespace) -> dict[str, object]:
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    fields_path = args.case_dir / "fields_final.npz"
    rows_path = args.case_dir / "tensor_interface_rows.jsonl"
    fields = np.load(fields_path)
    rows = load_rows(rows_path)

    x = np.asarray(fields["x"], dtype=float)
    z = np.asarray(fields["z"], dtype=float)
    rho_a = np.asarray(fields["rho_A"], dtype=float)
    rho_pull = np.asarray(fields["rho_D_to_A"], dtype=float)
    support = np.asarray(fields["support"], dtype=bool)
    trusted = np.asarray(fields["trusted"], dtype=bool)
    measure = np.asarray(fields["ntilde_measure_D"], dtype=float)
    metric_cov = np.asarray(fields["metric_cov_D"], dtype=float)
    u_t = np.asarray(fields["u_t_D"], dtype=float)
    u_x = np.asarray(fields["u_x_D"], dtype=float)
    u_z = np.asarray(fields["u_z_D"], dtype=float)
    raw_y = np.asarray(fields["tensor_raw_y"], dtype=float) if "tensor_raw_y" in fields else np.zeros_like(rho_a)

    metric_inv, _ = inverse_metric_block(metric_cov)
    sqrt_abs_g = safe_sqrt_abs_det(metric_cov)
    rho_tilde = np.where(sqrt_abs_g > 0.0, measure / np.maximum(sqrt_abs_g, 1.0e-300), 0.0)
    stress, _ = stress_tensor_tilde(metric_cov, metric_inv, rho_tilde, u_t, u_x, u_z, m=args.mass)
    rhs = np.asarray(stress, dtype=float) / (float(args.mp) * float(args.mp))
    rhs_norm = np.sqrt(np.einsum("...ab,...ab->...", rhs, rhs, optimize=True))

    xp = np.array([row["x"] for row in rows], dtype=float)
    zp = np.array([row["z"] for row in rows], dtype=float)
    row_rho = bilinear_sample(rho_a, x, z, xp, zp)
    row_rhs_norm = bilinear_sample(rhs_norm, x, z, xp, zp)
    row_alg_norm = np.array([abs(float(row.get("algebraic_tensor_norm", np.nan))) for row in rows], dtype=float)
    row_residual = np.array([abs(float(row.get("direct_tensor_residual_relative", np.nan))) for row in rows], dtype=float)
    solved = np.array([float(row.get("direct_tensor_solved", 0.0)) > 0.5 for row in rows], dtype=bool)
    source_to_geometry = row_rhs_norm / np.maximum(row_alg_norm, 1.0e-300)
    scale_needed = 1.0 / np.maximum(source_to_geometry, 1.0e-300)

    rho_weight = np.where(support, rho_a, 0.0)
    uniform_support_weight = support.astype(float)
    alphas = np.concatenate(
        [
            np.linspace(0.8, 1.2, 161),
            np.array([1.0e2, 1.0e10, 1.0e30, 1.0e60, 1.0e63]),
        ]
    )
    scan = []
    for alpha in alphas:
        candidate = alpha * rho_pull
        scan.append(
            {
                "alpha": float(alpha),
                "rho_weighted_rel_l1": weighted_rel_l1(rho_a, candidate, rho_weight, support),
                "rho_weighted_rel_l2": weighted_rel_l2(rho_a, candidate, rho_weight, support),
                "uniform_support_rel_l1": weighted_rel_l1(rho_a, candidate, uniform_support_weight, support),
                "uniform_support_rel_l2": weighted_rel_l2(rho_a, candidate, uniform_support_weight, support),
            }
        )
    best_l2 = min(scan, key=lambda item: item["rho_weighted_rel_l2"])
    best_l1 = min(scan, key=lambda item: item["rho_weighted_rel_l1"])

    accepted_rows = [row for row, ok in zip(rows, solved) if ok]
    rejected_rows = [row for row, ok in zip(rows, solved) if not ok]
    summary = {
        "case_dir": str(args.case_dir.resolve()),
        "scope": (
            "Feasibility check for a D-initial-data projection. It keeps the current A-derived geometry fixed "
            "and asks whether changing the transformed matter measure can both minimize the pulled-back rho error "
            "and repair D full-tensor interface residuals."
        ),
        "objective": {
            "rho_pullback_definition": "rho_pullback = m^2 * (sqrt(|gtilde|) rho_tilde) / |X_g|; with fixed metric and u, scaling rho_tilde scales rho_pullback linearly.",
            "best_alpha_weighted_l1": best_l1,
            "best_alpha_weighted_l2": best_l2,
        },
        "full_tensor_rows": {
            "count": int(len(rows)),
            "accepted_count": int(np.count_nonzero(solved)),
            "accepted_fraction": float(np.count_nonzero(solved) / max(len(rows), 1)),
            "rejected_count": int(len(rows) - np.count_nonzero(solved)),
            "residual_all": stats(row_residual),
            "residual_accepted": stats(row_residual[solved]),
            "residual_rejected": stats(row_residual[~solved]),
        },
        "density_location_of_rows": {
            "rho_at_all_rows": stats(row_rho),
            "rho_at_accepted_rows": stats(row_rho[solved]),
            "rho_at_rejected_rows": stats(row_rho[~solved]),
            "rho_weighted_accepted_fraction": float(np.sum(row_rho[solved]) / max(float(np.sum(row_rho)), 1.0e-300)),
        },
        "source_vs_geometry": {
            "rhs_tensor_norm_at_rows": stats(row_rhs_norm),
            "algebraic_tensor_norm_at_rows": stats(row_alg_norm),
            "rhs_over_algebraic_norm_at_rows": stats(source_to_geometry),
            "matter_scale_needed_for_order_one_geometric_effect": stats(scale_needed),
            "interpretation": (
                "If only rho_tilde is changed while metric jets are fixed, the D tensor condition changes through T/Mp^2. "
                "The required matter scaling for an order-one effect is the inverse of rhs/algebraic norm."
            ),
        },
        "conclusion": (
            "The rho-pullback objective is minimized at alpha≈1, i.e. the A-derived transformed density. "
            "But T/Mp^2 is tens of orders smaller than the geometric tensor-jump scale, so changing rho_tilde alone "
            "cannot make the rejected full-tensor rows satisfy the D equation without destroying the rho objective. "
            "A true D initial-data projection must adjust gtilde metric jets/interface geometry, not merely rho_tilde."
        ),
    }

    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.5), constrained_layout=True)
    im = axes[0, 0].pcolormesh(xg, zg, rho_a, shading="auto", cmap="viridis")
    axes[0, 0].add_collection(LineCollection(row_segments_um(rejected_rows, 3500), colors="tomato", linewidths=0.45, alpha=0.55))
    axes[0, 0].add_collection(LineCollection(row_segments_um(accepted_rows, 3500), colors="lime", linewidths=0.75, alpha=0.85))
    axes[0, 0].set_title("rho_A with tensor rows: green accepted, red rejected")
    fig.colorbar(im, ax=axes[0, 0])

    axes[0, 1].hist(np.log10(np.maximum(row_rho[~solved], 1.0e-300)), bins=80, alpha=0.65, label="rejected")
    axes[0, 1].hist(np.log10(np.maximum(row_rho[solved], 1.0e-300)), bins=80, alpha=0.65, label="accepted")
    axes[0, 1].set_xlabel("log10 rho_A at interface row center")
    axes[0, 1].set_ylabel("row count")
    axes[0, 1].legend()
    axes[0, 1].set_title("where accepted/rejected rows sit in rho")

    near = [item for item in scan if 0.8 <= item["alpha"] <= 1.2]
    axes[1, 0].plot([item["alpha"] for item in near], [item["rho_weighted_rel_l1"] for item in near], label="rho-weighted L1")
    axes[1, 0].plot([item["alpha"] for item in near], [item["rho_weighted_rel_l2"] for item in near], label="rho-weighted L2")
    axes[1, 0].axvline(1.0, color="black", linewidth=0.9, linestyle="--")
    axes[1, 0].set_xlabel("global scale alpha applied to rho_tilde measure")
    axes[1, 0].set_ylabel("weighted pointwise relative error")
    axes[1, 0].legend()
    axes[1, 0].set_title("rho objective: fixed geometry optimum is alpha=1")

    axes[1, 1].hist(np.log10(np.maximum(source_to_geometry, 1.0e-300)), bins=90, color="slateblue", alpha=0.8)
    axes[1, 1].set_xlabel("log10(||T/Mp^2|| / ||interface geometric tensor||)")
    axes[1, 1].set_ylabel("row count")
    axes[1, 1].set_title("matter source is negligible for tensor jump repair")

    for ax in [axes[0, 0]]:
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
    figure_path = out / "d_initial_projection_feasibility.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    summary["artifacts"] = {
        "figure": str(figure_path.resolve()),
        "summary": str((out / "summary.json").resolve()),
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mass", type=float, default=0.08061535492817485)
    parser.add_argument("--mp", type=float, default=PLANCK_MASS_EV)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
