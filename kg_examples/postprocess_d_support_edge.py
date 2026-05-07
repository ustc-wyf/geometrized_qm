from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from physical_units import HBAR_C_EV_M
from simulate_d_tridomain_full_dynamics import (
    dilate_mask_4,
    erode_mask_8,
    support_taper_weights,
    weighted_relative_l1,
)
from simulate_d_reduced_dynamic_same_initial import relative_l1


def load_case(path: Path, taper_decades: float, rho_floor: float) -> dict[str, object]:
    summary = json.loads((path / "summary.json").read_text(encoding="utf-8"))
    fields = np.load(path / "fields_final.npz")
    support = np.asarray(fields["support"], dtype=bool)
    trusted = np.asarray(fields["trusted"], dtype=bool)
    active = np.asarray(fields["active"], dtype=bool)
    if not np.any(trusted):
        trusted = erode_mask_8(support, int(summary["params"].get("trusted_erosion", 0)))
    support_edge = support & (~trusted)
    if not np.any(support_edge):
        support_edge = support.copy()

    rho_a = np.asarray(fields["rho_A"], dtype=float)
    measure_a = np.asarray(fields["ntilde_measure_A"], dtype=float)
    measure_d = np.asarray(fields["ntilde_measure_D"], dtype=float)
    taper = support_taper_weights(
        rho=rho_a,
        measure=measure_a,
        rho_frac=float(summary["params"].get("support_rho_frac", 1.0e-3)),
        measure_frac=float(summary["params"].get("support_measure_frac", 1.0e-3)),
        decades=taper_decades,
        floor=rho_floor,
    )
    taper_active = np.where(active, taper, 0.0)
    active_from_support = dilate_mask_4(support, int(summary["params"].get("active_dilation", 0)))
    active_diff_count = int(np.count_nonzero(active_from_support ^ active))

    diff = measure_d - measure_a
    denom = np.maximum(np.abs(measure_a), 1.0e-300)
    point_rel = np.abs(diff) / denom
    result = {
        "label": path.name,
        "path": str(path.resolve()),
        "summary": str((path / "summary.json").resolve()),
        "fields": str((path / "fields_final.npz").resolve()),
        "figure": None,
        "resolution": int(summary["params"]["resolution"]),
        "dt": float(summary["params"]["dt"]),
        "dt_fs": summary["params"].get("dt_fs"),
        "initial_time_old_units": summary["params"].get("initial_time_old_units"),
        "t_completed_fs": summary["params"].get("t_completed_fs"),
        "support_count": int(np.count_nonzero(support)),
        "trusted_count": int(np.count_nonzero(trusted)),
        "support_edge_count": int(np.count_nonzero(support_edge)),
        "active_count": int(np.count_nonzero(active)),
        "active_wrap_diff_count": active_diff_count,
        "taper_weight_sum": float(np.sum(taper_active)),
        "taper_transition_count": int(np.count_nonzero((taper > 0.0) & (taper < 1.0))),
        "taper_outside_active_fraction": float(np.sum(taper[~active]) / max(float(np.sum(taper)), 1.0e-300)),
        "measure_rel_l1_support": relative_l1(measure_a, measure_d, support),
        "measure_rel_l1_trusted": relative_l1(measure_a, measure_d, trusted),
        "measure_rel_l1_support_edge": relative_l1(measure_a, measure_d, support_edge),
        "measure_rel_l1_tapered": weighted_relative_l1(measure_a, measure_d, taper_active),
        "point_rel_p95_support": float(np.percentile(point_rel[support], 95.0)),
        "point_rel_p95_trusted": float(np.percentile(point_rel[trusted], 95.0)),
        "point_rel_p95_support_edge": float(np.percentile(point_rel[support_edge], 95.0)),
    }
    return {
        "summary": summary,
        "fields": fields,
        "support": support,
        "trusted": trusted,
        "active": active,
        "support_edge": support_edge,
        "taper": taper,
        "taper_active": taper_active,
        "point_rel": point_rel,
        "result": result,
    }


def render_case(out_path: Path, case: dict[str, object], axis_unit: str) -> None:
    fields = case["fields"]
    support = case["support"]
    trusted = case["trusted"]
    active = case["active"]
    support_edge = case["support_edge"]
    taper = case["taper"]
    point_rel = case["point_rel"]
    result = case["result"]

    x = np.asarray(fields["x"], dtype=float)
    z = np.asarray(fields["z"], dtype=float)
    if axis_unit == "um":
        factor = HBAR_C_EV_M * 1.0e6
        x = x * factor
        z = z * factor
    xg, zg = np.meshgrid(x, z, indexing="ij")
    measure_a = np.asarray(fields["ntilde_measure_A"], dtype=float)
    measure_d = np.asarray(fields["ntilde_measure_D"], dtype=float)
    diff = measure_d - measure_a

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.2), constrained_layout=True)
    vmax = max(float(np.percentile(measure_a[support], 99.5)), float(np.percentile(measure_d[support], 99.5)), 1.0e-16)
    im0 = axes[0, 0].pcolormesh(xg, zg, measure_a, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.8)
    axes[0, 0].contour(xg, zg, trusted.astype(float), levels=[0.5], colors="cyan", linewidths=0.8)
    axes[0, 0].set_title("A measure: white support, cyan trusted")
    fig.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].pcolormesh(xg, zg, taper, shading="auto", cmap="magma", vmin=0.0, vmax=1.0)
    axes[0, 1].contour(xg, zg, active.astype(float), levels=[0.5], colors="white", linewidths=0.8)
    axes[0, 1].set_title("diagnostic tapered-support weight")
    fig.colorbar(im1, ax=axes[0, 1])

    dmax = max(float(np.percentile(np.abs(diff[support]), 99.0)), 1.0e-16)
    im2 = axes[0, 2].pcolormesh(xg, zg, diff, shading="auto", cmap="coolwarm", vmin=-dmax, vmax=dmax)
    axes[0, 2].contour(xg, zg, support_edge.astype(float), levels=[0.5], colors="black", linewidths=0.8)
    axes[0, 2].set_title("D-A measure; black support edge")
    fig.colorbar(im2, ax=axes[0, 2])

    category = np.zeros_like(taper)
    category[active] = 1.0
    category[support] = 2.0
    category[trusted] = 3.0
    im3 = axes[1, 0].pcolormesh(xg, zg, category, shading="auto", cmap="cividis", vmin=0.0, vmax=3.0)
    axes[1, 0].set_title("mask categories: 1 active, 2 support, 3 trusted")
    fig.colorbar(im3, ax=axes[1, 0])

    rel_clip = max(float(np.percentile(point_rel[support], 99.0)), 1.0e-16)
    im4 = axes[1, 1].pcolormesh(
        xg,
        zg,
        np.clip(point_rel, 0.0, rel_clip),
        shading="auto",
        cmap="inferno",
        vmin=0.0,
        vmax=rel_clip,
    )
    axes[1, 1].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.6)
    axes[1, 1].set_title("|D-A|/|A| pointwise, clipped at support p99")
    fig.colorbar(im4, ax=axes[1, 1])

    labels = ["support", "trusted", "edge", "taper"]
    vals = [
        result["measure_rel_l1_support"],
        result["measure_rel_l1_trusted"],
        result["measure_rel_l1_support_edge"],
        result["measure_rel_l1_tapered"],
    ]
    axes[1, 2].bar(labels, vals, color=["#4c78a8", "#54a24b", "#f58518", "#b279a2"])
    axes[1, 2].set_yscale("log")
    axes[1, 2].set_title("final measure rel L1 by support expression")
    axes[1, 2].grid(axis="y", alpha=0.25)

    for ax in axes.ravel()[:5]:
        ax.set_xlabel(f"x [{axis_unit}]")
        ax.set_ylabel(f"z [{axis_unit}]")
        ax.set_aspect("equal")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", type=Path, action="append", required=True)
    parser.add_argument("--taper-decades", type=float, default=1.0)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--axis-unit", choices=["code", "um"], default="um")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    cases = [load_case(path, args.taper_decades, args.rho_floor) for path in args.case]
    results = []
    for case in cases:
        result = case["result"]
        fig = args.output / f"{result['label']}_support_edge_taper.png"
        render_case(fig, case, args.axis_unit)
        result["figure"] = str(fig.resolve())
        results.append(result)

    fig, ax = plt.subplots(figsize=(9.4, 5.4), constrained_layout=True)
    labels = [r["label"].replace("d_physical_1550nm_kg_", "") for r in results]
    x = np.arange(len(results))
    width = 0.2
    series = [
        ("support", [r["measure_rel_l1_support"] for r in results]),
        ("trusted", [r["measure_rel_l1_trusted"] for r in results]),
        ("edge", [r["measure_rel_l1_support_edge"] for r in results]),
        ("taper", [r["measure_rel_l1_tapered"] for r in results]),
    ]
    for idx, (name, vals) in enumerate(series):
        ax.bar(x + (idx - 1.5) * width, vals, width=width, label=name)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_yscale("log")
    ax.set_ylabel("final measure relative L1")
    ax.set_title("Support-edge expression postprocess: no equation/source changes")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    overview = args.output / "support_edge_taper_overview.png"
    fig.savefig(overview, dpi=180)
    plt.close(fig)

    payload = {
        "scope": "Postprocess existing D-branch fields only.  This changes no evolution equation, source, support mask used during the run, or physical boundary condition.",
        "taper_decades": args.taper_decades,
        "results": results,
        "overview_figure": str(overview.resolve()),
        "summary": str((args.output / "support_edge_taper_summary.json").resolve()),
    }
    (args.output / "support_edge_taper_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
