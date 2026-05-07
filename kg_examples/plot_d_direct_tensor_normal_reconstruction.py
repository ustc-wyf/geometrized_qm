from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from diagnose_d_direct_tensor_interface import direct_tensor_rows
from diagnose_d_trace_interface_weak_form import parse_positive_levels, real_array
from diagnose_d_tensor_interface_jump import reference_tensor_data
from integrate_d_interface_jumps import scalar_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def linelet_segments(
    rows: list[dict[str, float]],
    half_length: float,
    *,
    tangent_from_target: bool,
) -> list[list[tuple[float, float]]]:
    segments: list[list[tuple[float, float]]] = []
    for row in rows:
        fx = row.get("target_covector_x", np.nan)
        fz = row.get("target_covector_z", np.nan)
        spatial_norm = float(np.hypot(fx, fz))
        if not np.isfinite(spatial_norm) or spatial_norm <= 1.0e-300:
            continue
        if tangent_from_target:
            tx, tz = -fz / spatial_norm, fx / spatial_norm
        else:
            tx, tz = fx / spatial_norm, fz / spatial_norm
        xc = row["x"]
        zc = row["z"]
        segments.append(
            [
                (xc - half_length * tx, zc - half_length * tz),
                (xc + half_length * tx, zc + half_length * tz),
            ]
        )
    return segments


def add_collection(
    ax: plt.Axes,
    segments: list[list[tuple[float, float]]],
    *,
    colors: str | list[str] = "black",
    values: np.ndarray | None = None,
    cmap: str = "viridis",
    linewidth: float = 0.8,
    linestyle: str = "solid",
    vmin: float | None = None,
    vmax: float | None = None,
) -> LineCollection | None:
    if not segments:
        return None
    if values is None:
        collection = LineCollection(segments, colors=colors, linewidths=linewidth, linestyles=linestyle)
    else:
        collection = LineCollection(segments, cmap=cmap, linewidths=linewidth, linestyles=linestyle)
        collection.set_array(values)
        if vmin is not None and vmax is not None:
            collection.set_clim(vmin, vmax)
    ax.add_collection(collection)
    return collection


def finite_values(rows: list[dict[str, float]], key: str) -> np.ndarray:
    vals = np.asarray([row.get(key, np.nan) for row in rows], dtype=float)
    return vals[np.isfinite(vals)]


def summarize(rows: list[dict[str, float]]) -> dict[str, object]:
    candidates = [row for row in rows if row.get("candidate", 0.0) > 0.5]
    accepted = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    rejected = [row for row in rows if row.get("direct_tensor_solved", 0.0) <= 0.5]
    return {
        "count": len(rows),
        "candidate_count": len(candidates),
        "candidate_fraction": float(len(candidates) / max(len(rows), 1)),
        "accepted_count": len(accepted),
        "accepted_fraction": float(len(accepted) / max(len(rows), 1)),
        "rejected_count": len(rejected),
        "accepted_direct_tensor_residual_relative": scalar_stats(
            [row["direct_tensor_residual_relative"] for row in accepted]
        ),
        "accepted_rank1_tail_relative": scalar_stats([row["rank1_tail_relative"] for row in accepted]),
        "accepted_target_speed": scalar_stats([row["direct_speed"] for row in accepted]),
        "accepted_target_current_alignment": scalar_stats(
            [row["spatial_alignment_with_current"] for row in accepted]
        ),
        "all_candidate_direct_tensor_residual_relative": scalar_stats(
            [row["direct_tensor_residual_relative"] for row in candidates]
        ),
        "all_candidate_rank1_tail_relative": scalar_stats([row["rank1_tail_relative"] for row in candidates]),
        "all_candidate_target_current_alignment": scalar_stats(
            [row["spatial_alignment_with_current"] for row in candidates]
        ),
    }


def render(
    out_path: Path,
    *,
    x: np.ndarray,
    z: np.ndarray,
    rho: np.ndarray,
    rows: list[dict[str, float]],
    rank1_tol: float,
    tensor_rel_tol: float,
    rho_min_frac: float,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    support = rho > rho_min_frac * float(np.max(rho))
    accepted = [row for row in rows if row.get("direct_tensor_solved", 0.0) > 0.5]
    rejected = [row for row in rows if row.get("direct_tensor_solved", 0.0) <= 0.5]
    candidates = [row for row in rows if row.get("candidate", 0.0) > 0.5]

    dx = float(abs(x[1] - x[0]))
    dz = float(abs(z[1] - z[0]))
    half_length = 0.45 * float(np.hypot(dx, dz))
    target_tangents = linelet_segments(accepted, half_length, tangent_from_target=True)
    target_normals = linelet_segments(accepted, 0.32 * half_length, tangent_from_target=False)

    fig, axes = plt.subplots(2, 3, figsize=(16.6, 9.8), constrained_layout=True)
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))

    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_collection(axes[0, 0], row_segments(rows), colors=["#1b9e77" if r["level"] < 0.0 else "#d95f02" for r in rows])
    axes[0, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7, linestyles="dashed")
    axes[0, 0].set_title("linear rho + current raw_y=+/-1")
    fig.colorbar(im0, ax=axes[0, 0])

    add_collection(axes[0, 1], row_segments(rejected), colors="0.55", linewidth=0.6, linestyle="dotted")
    add_collection(axes[0, 1], row_segments(accepted), colors="black", linewidth=0.9)
    axes[0, 1].set_title(f"direct tensor accepted (rank/res <= {rank1_tol:g}/{tensor_rel_tol:g})")

    add_collection(axes[0, 2], row_segments(rejected), colors="0.75", linewidth=0.45, linestyle="dotted")
    add_collection(axes[0, 2], target_tangents, colors="#005f73", linewidth=0.9)
    add_collection(axes[0, 2], target_normals, colors="#ae2012", linewidth=0.7)
    axes[0, 2].set_title("target linelets: blue tangent, red normal")

    if candidates:
        vals = np.asarray([row["direct_tensor_residual_relative"] for row in candidates], dtype=float)
        vmax = max(float(np.percentile(vals[np.isfinite(vals)], 95.0)), 1.0e-12)
        lc = add_collection(
            axes[1, 0],
            row_segments(candidates),
            values=vals,
            cmap="magma",
            linewidth=1.0,
            vmin=0.0,
            vmax=vmax,
        )
        fig.colorbar(lc, ax=axes[1, 0])
    axes[1, 0].set_title("direct tensor residual / algebraic")

    if candidates:
        vals = np.asarray([row["spatial_alignment_with_current"] for row in candidates], dtype=float)
        lc = add_collection(
            axes[1, 1],
            row_segments(candidates),
            values=vals,
            cmap="viridis",
            linewidth=1.0,
            vmin=0.0,
            vmax=1.0,
        )
        fig.colorbar(lc, ax=axes[1, 1])
    axes[1, 1].set_title("target-current normal alignment")

    if candidates:
        vals = np.asarray([row["direct_speed"] for row in candidates], dtype=float)
        finite = vals[np.isfinite(vals)]
        vmax = max(float(np.percentile(np.abs(finite), 95.0)), 1.0e-12) if finite.size else 1.0
        lc = add_collection(
            axes[1, 2],
            row_segments(candidates),
            values=vals,
            cmap="coolwarm",
            linewidth=1.0,
            vmin=-vmax,
            vmax=vmax,
        )
        fig.colorbar(lc, ax=axes[1, 2])
    axes[1, 2].set_title("target coordinate speed -F_t/|F_spatial|")

    for ax in axes.ravel():
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_xlim(float(x[0]), float(x[-1]))
        ax.set_ylim(float(z[0]), float(z[-1]))
        ax.set_aspect("equal")

    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    params = FlatLocalizedCrossingParams(nx=args.resolution, nz=args.resolution)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    data = reference_tensor_data(
        x=x,
        z=z,
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        t=args.time,
        probe_dt=args.probe_dt,
        ell=args.ell,
        mp=args.mp,
        mass=params.m,
        rho_floor=args.rho_floor,
    )

    levels = parse_positive_levels(args.levels)
    rows: list[dict[str, float]] = []
    for level in [sign * lev for lev in levels for sign in (-1.0, 1.0)]:
        rows.extend(
            direct_tensor_rows(
                raw_y=real_array(data["raw_y"]),
                metric_cov=real_array(data["metric_cov"]),
                metric_inv=real_array(data["metric_inv"]),
                ricci=real_array(data["ricci"]),
                rho=real_array(data["rho"]),
                rhs_tensor=real_array(data["rhs_tensor"]),
                x=x,
                z=z,
                level=level,
                ell=args.ell,
                half_width_y=args.half_width_y,
                samples=args.samples,
                rank1_tol=args.rank1_tol,
                tensor_rel_tol=args.tensor_rel_tol,
            )
        )

    fig_path = out / f"d_direct_tensor_normal_reconstruction_ell{args.ell:g}_n{args.resolution}_t{args.time:g}.png"
    render(
        fig_path,
        x=x,
        z=z,
        rho=real_array(data["rho"]),
        rows=rows,
        rank1_tol=args.rank1_tol,
        tensor_rel_tol=args.tensor_rel_tol,
        rho_min_frac=args.rho_min_frac,
    )

    rows_path = out / "direct_tensor_interface_rows.jsonl"
    with rows_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")

    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolution": args.resolution,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "levels_abs_raw_y": levels,
            "half_width_y": args.half_width_y,
            "samples": args.samples,
            "rho_min_frac": args.rho_min_frac,
            "rank1_tol": args.rank1_tol,
            "tensor_rel_tol": args.tensor_rel_tol,
        },
        "definitions": {
            "direct_tensor_condition": "The acceptance condition is the full leading tensor equation H_ab=-F_aF_b+g_ab F^2. No trace-speed solve is used.",
            "accepted": "A segment is accepted only when the tensor target is a real rank-one covector and its full tensor residual passes the configured numerical tolerances.",
            "target_linelets": "Blue linelets are local tangents perpendicular to the direct-tensor target spatial covector; red linelets show that target normal. They are a reconstruction diagnostic, not a new physics rule.",
        },
        "summary": summarize(rows),
        "files": {
            "figure": str(fig_path.resolve()),
            "rows_jsonl": str(rows_path.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolution", type=int, default=128)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument("--rho-min-frac", type=float, default=1.0e-3)
    parser.add_argument("--rank1-tol", type=float, default=0.1)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
