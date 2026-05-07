from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap

from integrate_d_interface_jumps import scalar_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot


MODE_ORDER = [
    "accepted",
    "no_real_covector",
    "rank_tail_only",
    "negative_only",
    "tensor_residual_only",
    "rank_tail+negative",
    "rank_tail+tensor_residual",
    "negative+tensor_residual",
    "rank_tail+negative+tensor_residual",
]

MODE_COLORS = [
    "#111111",
    "#8d99ae",
    "#0077b6",
    "#f77f00",
    "#d62828",
    "#2a9d8f",
    "#9d4edd",
    "#e76f51",
    "#6d597a",
]


def read_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def classify(row: dict[str, float], *, rank1_tol: float, tensor_rel_tol: float) -> str:
    if row.get("direct_tensor_solved", 0.0) > 0.5:
        return "accepted"
    if row.get("candidate", 0.0) <= 0.5:
        return "no_real_covector"
    failures: list[str] = []
    if row.get("rank1_tail_relative", np.inf) > rank1_tol:
        failures.append("rank_tail")
    if row.get("rank1_negative_relative", np.inf) > rank1_tol:
        failures.append("negative")
    if row.get("direct_tensor_residual_relative", np.inf) > tensor_rel_tol:
        failures.append("tensor_residual")
    return "+".join(failures) if failures else "unclassified_numeric"


def add_lines(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    *,
    colors: str | list[str] = "black",
    linewidth: float = 0.85,
    linestyle: str = "solid",
) -> LineCollection | None:
    if not rows:
        return None
    collection = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth, linestyles=linestyle)
    ax.add_collection(collection)
    return collection


def mode_stats(rows: list[dict[str, float]], modes: list[str]) -> dict[str, object]:
    grouped: dict[str, list[dict[str, float]]] = defaultdict(list)
    for row, mode in zip(rows, modes):
        grouped[mode].append(row)
    out: dict[str, object] = {}
    for mode in MODE_ORDER + sorted(set(grouped) - set(MODE_ORDER)):
        group = grouped.get(mode, [])
        if not group:
            continue
        out[mode] = {
            "count": len(group),
            "fraction": float(len(group) / max(len(rows), 1)),
            "rank1_tail_relative": scalar_stats([row.get("rank1_tail_relative", np.nan) for row in group]),
            "rank1_negative_relative": scalar_stats([row.get("rank1_negative_relative", np.nan) for row in group]),
            "direct_tensor_residual_relative": scalar_stats(
                [row.get("direct_tensor_residual_relative", np.nan) for row in group]
            ),
            "spatial_alignment_with_current": scalar_stats(
                [row.get("spatial_alignment_with_current", np.nan) for row in group]
            ),
            "direct_speed": scalar_stats([row.get("direct_speed", np.nan) for row in group]),
        }
    return out


def reference_rho(x: np.ndarray, z: np.ndarray, *, time: float, rho_floor: float) -> np.ndarray:
    params = FlatLocalizedCrossingParams(nx=x.size, nz=z.size)
    X, Z = np.meshgrid(x, z, indexing="ij")
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    snap = localized_direct_tilde_coordinate_snapshot(psi0, psi0_hat, omega, x, z, time, rho_floor=rho_floor)
    return np.real_if_close(snap["bohm"]["rho"]).astype(float)


def render(
    out_path: Path,
    *,
    rows: list[dict[str, float]],
    modes: list[str],
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.8), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))

    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_lines(axes[0, 0], rows, colors=["#1b9e77" if row["level"] < 0.0 else "#d95f02" for row in rows], linewidth=0.7)
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0, 0].set_title("linear rho + current raw_y=+/-1")
    fig.colorbar(im0, ax=axes[0, 0])

    mode_to_idx = {mode: i for i, mode in enumerate(MODE_ORDER)}
    colors = [MODE_COLORS[mode_to_idx.get(mode, len(MODE_COLORS) - 1)] for mode in modes]
    add_lines(axes[0, 1], rows, colors=colors, linewidth=0.9)
    axes[0, 1].set_title("direct tensor rejection modes")

    accepted = [row for row, mode in zip(rows, modes) if mode == "accepted"]
    unresolved = [row for row, mode in zip(rows, modes) if mode != "accepted"]
    add_lines(axes[0, 2], unresolved, colors="0.75", linewidth=0.5, linestyle="dotted")
    add_lines(axes[0, 2], accepted, colors="black", linewidth=0.9)
    axes[0, 2].set_title("accepted black; unresolved dotted")

    candidate_rows = [row for row in rows if row.get("candidate", 0.0) > 0.5]
    for ax, key, title, vmax_default in [
        (axes[1, 0], "rank1_tail_relative", "rank-one tail mismatch", None),
        (axes[1, 1], "rank1_negative_relative", "negative eigenweight / leading", None),
        (axes[1, 2], "direct_tensor_residual_relative", "direct tensor residual / algebraic", None),
    ]:
        if not candidate_rows:
            ax.set_title(title)
            continue
        vals = np.asarray([row.get(key, np.nan) for row in candidate_rows], dtype=float)
        finite = vals[np.isfinite(vals)]
        vmax = max(float(np.percentile(finite, 95.0)), 1.0e-12) if finite.size else (vmax_default or 1.0)
        collection = LineCollection(row_segments(candidate_rows), cmap="magma", linewidths=0.95)
        collection.set_array(vals)
        collection.set_clim(0.0, vmax)
        ax.add_collection(collection)
        fig.colorbar(collection, ax=ax)
        ax.set_title(title)

    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
        ax.set_aspect("equal")

    handles = [
        plt.Line2D([0], [0], color=MODE_COLORS[i], lw=3, label=mode)
        for i, mode in enumerate(MODE_ORDER)
        if mode in set(modes)
    ]
    if handles:
        axes[0, 1].legend(handles=handles, loc="upper right", fontsize=7)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    rows = read_rows(args.rows_jsonl)
    modes = [classify(row, rank1_tol=args.rank1_tol, tensor_rel_tol=args.tensor_rel_tol) for row in rows]

    params = FlatLocalizedCrossingParams(nx=args.resolution, nz=args.resolution)
    x, z, _, _ = make_grid(params)
    rho = reference_rho(x, z, time=args.time, rho_floor=args.rho_floor)

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    fig_path = out / f"d_direct_tensor_rejection_modes_n{args.resolution}_t{args.time:g}.png"
    render(fig_path, rows=rows, modes=modes, x=x, z=z, rho=rho, rho_min_frac=args.rho_min_frac)

    counts = Counter(modes)
    summary = {
        "params": {
            "rows_jsonl": str(args.rows_jsonl.resolve()),
            "rank1_tol": args.rank1_tol,
            "tensor_rel_tol": args.tensor_rel_tol,
            "time": args.time,
            "resolution": args.resolution,
            "rho_floor": args.rho_floor,
        },
        "definitions": {
            "accepted": "Passed the direct tensor rank-one and full residual tests.",
            "no_real_covector": "The tensor target could not produce a positive real rank-one covector candidate.",
            "rank_tail": "g_ab F^2-H_ab is not close enough to rank one.",
            "negative": "The rank-one target has too much negative eigenweight for a real covector outer product.",
            "tensor_residual": "The reconstructed covector fails the full tensor residual tolerance.",
        },
        "summary": {
            "count": len(rows),
            "mode_counts": dict(counts),
            "mode_fractions": {key: float(value / max(len(rows), 1)) for key, value in counts.items()},
            "mode_stats": mode_stats(rows, modes),
        },
        "files": {"figure": str(fig_path.resolve())},
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows-jsonl", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rank1-tol", type=float, default=0.1)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolution", type=int, default=128)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
