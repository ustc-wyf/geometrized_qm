from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

from diagnose_d_direct_tensor_failure_features import add_features
from diagnose_d_tensor_interface_jump import reference_tensor_data
from diagnose_d_trace_interface_weak_form import parse_positive_levels, real_array
from integrate_d_interface_jumps import scalar_stats
from prototype_d_tensor_reconstructed_interface import least_squares_tensor_rows
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


@dataclass(frozen=True)
class Window:
    name: str
    xmin: float
    xmax: float
    zmin: float
    zmax: float

    def contains(self, row: dict[str, float]) -> bool:
        return self.xmin <= row["x"] <= self.xmax and self.zmin <= row["z"] <= self.zmax


def parse_ints(text: str) -> list[int]:
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_lines(
    ax: plt.Axes,
    rows: list[dict[str, float]],
    *,
    colors: str | list[str] = "black",
    values: np.ndarray | None = None,
    cmap: str = "magma",
    linewidth: float = 0.8,
    linestyle: str = "solid",
    vmin: float | None = None,
    vmax: float | None = None,
) -> LineCollection | None:
    if not rows:
        return None
    if values is None:
        collection = LineCollection(row_segments(rows), colors=colors, linewidths=linewidth, linestyles=linestyle)
    else:
        collection = LineCollection(row_segments(rows), cmap=cmap, linewidths=linewidth, linestyles=linestyle)
        collection.set_array(values)
        if vmin is not None and vmax is not None:
            collection.set_clim(vmin, vmax)
    ax.add_collection(collection)
    return collection


def compute_rows(args: argparse.Namespace, resolution: int) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], list[dict[str, float]]]:
    params = FlatLocalizedCrossingParams(nx=resolution, nz=resolution)
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
            least_squares_tensor_rows(
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
                tensor_rel_tol=args.tensor_rel_tol,
                max_iter=args.max_iter,
                multistart=args.multistart,
                max_seeds=args.max_seeds,
            )
        )
    rows = add_features(rows, data=data, x=x, z=z)
    return x, z, data, rows


def row_mode(row: dict[str, float], *, tensor_rel_tol: float) -> str:
    if row.get("direct_tensor_solved", 0.0) > 0.5:
        return "accepted"
    if row.get("candidate", 0.0) <= 0.5:
        return "no_candidate"
    if row.get("direct_tensor_residual_relative", np.inf) > tensor_rel_tol:
        return "ls_residual_fail"
    return "unclassified_unresolved"


def summarize_rows(rows: list[dict[str, float]], *, tensor_rel_tol: float) -> dict[str, object]:
    modes = [row_mode(row, tensor_rel_tol=tensor_rel_tol) for row in rows]
    accepted = [row for row, mode in zip(rows, modes) if mode == "accepted"]
    unresolved = [row for row, mode in zip(rows, modes) if mode != "accepted"]
    candidates = [row for row in rows if row.get("candidate", 0.0) > 0.5]
    return {
        "row_count": len(rows),
        "candidate_count": len(candidates),
        "candidate_fraction": float(len(candidates) / max(len(rows), 1)),
        "accepted_count": len(accepted),
        "accepted_fraction": float(len(accepted) / max(len(rows), 1)),
        "unresolved_count": len(unresolved),
        "unresolved_fraction": float(len(unresolved) / max(len(rows), 1)),
        "mode_counts": {mode: int(modes.count(mode)) for mode in sorted(set(modes))},
        "direct_tensor_residual_relative_all_candidates": scalar_stats(
            [row.get("direct_tensor_residual_relative", np.nan) for row in candidates]
        ),
        "direct_tensor_residual_relative_accepted": scalar_stats(
            [row.get("direct_tensor_residual_relative", np.nan) for row in accepted]
        ),
        "direct_tensor_residual_relative_unresolved": scalar_stats(
            [row.get("direct_tensor_residual_relative", np.nan) for row in unresolved]
        ),
        "rho_relative_unresolved": scalar_stats([row.get("rho_relative", np.nan) for row in unresolved]),
        "raw_y_grad_norm_unresolved": scalar_stats([row.get("raw_y_grad_norm", np.nan) for row in unresolved]),
        "metric_condition_unresolved": scalar_stats([row.get("metric_condition_abs_eig", np.nan) for row in unresolved]),
    }


def select_windows(
    rows: list[dict[str, float]],
    *,
    tensor_rel_tol: float,
    focus_rho_min_frac: float,
    width_x: float,
    width_z: float,
    top_count: int,
    xmin: float,
    xmax: float,
    zmin: float,
    zmax: float,
) -> list[Window]:
    unresolved = [
        row
        for row in rows
        if row_mode(row, tensor_rel_tol=tensor_rel_tol) != "accepted"
        and row.get("rho_relative", 0.0) >= focus_rho_min_frac
    ]
    if not unresolved:
        return []

    bins: dict[tuple[int, int], int] = {}
    for row in unresolved:
        ix = int(np.floor((row["x"] - xmin) / width_x))
        iz = int(np.floor((row["z"] - zmin) / width_z))
        bins[(ix, iz)] = bins.get((ix, iz), 0) + 1

    selected: list[Window] = []
    for (ix, iz), _count in sorted(bins.items(), key=lambda item: item[1], reverse=True):
        wx0 = xmin + ix * width_x
        wz0 = zmin + iz * width_z
        win = Window(
            name=f"W{len(selected)+1}",
            xmin=max(xmin, wx0),
            xmax=min(xmax, wx0 + width_x),
            zmin=max(zmin, wz0),
            zmax=min(zmax, wz0 + width_z),
        )
        overlaps = False
        for old in selected:
            x_overlap = max(0.0, min(win.xmax, old.xmax) - max(win.xmin, old.xmin))
            z_overlap = max(0.0, min(win.zmax, old.zmax) - max(win.zmin, old.zmin))
            if x_overlap * z_overlap > 0.25 * width_x * width_z:
                overlaps = True
                break
        if overlaps:
            continue
        selected.append(win)
        if len(selected) >= top_count:
            break
    return selected


def summarize_windows(
    rows: list[dict[str, float]],
    windows: list[Window],
    *,
    tensor_rel_tol: float,
) -> dict[str, object]:
    output: dict[str, object] = {}
    for window in windows:
        subset = [row for row in rows if window.contains(row)]
        output[window.name] = {
            "bbox": {
                "xmin": window.xmin,
                "xmax": window.xmax,
                "zmin": window.zmin,
                "zmax": window.zmax,
            },
            "summary": summarize_rows(subset, tensor_rel_tol=tensor_rel_tol),
        }
    return output


def render(
    out_path: Path,
    *,
    results: dict[int, dict[str, object]],
    windows: list[Window],
    highest_resolution: int,
    tensor_rel_tol: float,
) -> None:
    high = results[highest_resolution]
    x = high["x"]
    z = high["z"]
    data = high["data"]
    rows = high["rows"]
    xg, zg = np.meshgrid(x, z, indexing="ij")
    rho = real_array(data["rho"])
    accepted = [row for row in rows if row_mode(row, tensor_rel_tol=tensor_rel_tol) == "accepted"]
    unresolved = [row for row in rows if row_mode(row, tensor_rel_tol=tensor_rel_tol) != "accepted"]

    fig, axes = plt.subplots(2, 2, figsize=(14.2, 11.2), constrained_layout=True)
    support = rho > 1.0e-3 * float(np.max(rho))
    rho_vmax = max(float(np.percentile(rho[support], 99.5)), 1.0e-16) if np.any(support) else float(np.max(rho))
    im0 = axes[0, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=rho_vmax)
    add_lines(axes[0, 0], unresolved, colors="#d62828", linewidth=0.85)
    add_lines(axes[0, 0], accepted, colors="black", linewidth=0.55)
    for window in windows:
        rect = Rectangle(
            (window.xmin, window.zmin),
            window.xmax - window.xmin,
            window.zmax - window.zmin,
            fill=False,
            edgecolor="white",
            linewidth=1.6,
        )
        axes[0, 0].add_patch(rect)
        axes[0, 0].text(window.xmin, window.zmax, window.name, color="white", fontsize=9, va="bottom")
    axes[0, 0].set_title(f"highest n={highest_resolution}: rho + accepted black / unresolved red")
    fig.colorbar(im0, ax=axes[0, 0])

    if unresolved:
        vals = np.asarray([row.get("direct_tensor_residual_relative", np.nan) for row in unresolved], dtype=float)
        finite = vals[np.isfinite(vals)]
        vmax = max(float(np.percentile(finite, 95.0)), 1.0e-12) if finite.size else 1.0
        lc = add_lines(axes[0, 1], unresolved, values=vals, cmap="magma", linewidth=1.0, vmin=0.0, vmax=vmax)
        fig.colorbar(lc, ax=axes[0, 1])
    add_lines(axes[0, 1], accepted, colors="0.72", linewidth=0.45)
    for window in windows:
        axes[0, 1].add_patch(
            Rectangle(
                (window.xmin, window.zmin),
                window.xmax - window.xmin,
                window.zmax - window.zmin,
                fill=False,
                edgecolor="#1d3557",
                linewidth=1.2,
            )
        )
    axes[0, 1].set_title("unresolved residual after least-squares")

    resolutions = sorted(results)
    global_accept = [results[n]["summary"]["accepted_fraction"] for n in resolutions]
    axes[1, 0].plot(resolutions, global_accept, marker="o", linewidth=2.0, label="global")
    for window in windows:
        vals = [
            results[n]["windows"][window.name]["summary"]["accepted_fraction"]
            if results[n]["windows"][window.name]["summary"]["row_count"] > 0
            else np.nan
            for n in resolutions
        ]
        axes[1, 0].plot(resolutions, vals, marker="o", linewidth=1.5, label=window.name)
    axes[1, 0].set_ylim(0.0, 1.02)
    axes[1, 0].set_xlabel("resolution")
    axes[1, 0].set_ylabel("accepted fraction")
    axes[1, 0].set_title("local convergence of accepted fraction")
    axes[1, 0].legend(loc="best", fontsize=8)
    axes[1, 0].grid(alpha=0.25)

    global_unres_p95 = [
        results[n]["summary"]["direct_tensor_residual_relative_unresolved"]["p95"] for n in resolutions
    ]
    axes[1, 1].plot(resolutions, global_unres_p95, marker="o", linewidth=2.0, label="global unresolved p95")
    for window in windows:
        vals = [
            results[n]["windows"][window.name]["summary"]["direct_tensor_residual_relative_unresolved"]["p95"]
            if results[n]["windows"][window.name]["summary"]["unresolved_count"] > 0
            else np.nan
            for n in resolutions
        ]
        axes[1, 1].plot(resolutions, vals, marker="o", linewidth=1.5, label=f"{window.name} unresolved p95")
    axes[1, 1].set_xlabel("resolution")
    axes[1, 1].set_ylabel("relative tensor residual p95")
    axes[1, 1].set_title("least-squares unresolved residual")
    axes[1, 1].legend(loc="best", fontsize=8)
    axes[1, 1].grid(alpha=0.25)

    for ax in axes[:1, :].ravel():
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
    resolutions = parse_ints(args.resolutions)
    if not resolutions:
        raise ValueError("At least one resolution is required.")

    results: dict[int, dict[str, object]] = {}
    for resolution in resolutions:
        x, z, data, rows = compute_rows(args, resolution)
        rows_path = out / f"least_squares_rows_n{resolution}.jsonl"
        with rows_path.open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
        results[resolution] = {
            "x": x,
            "z": z,
            "data": data,
            "rows": rows,
            "rows_jsonl": str(rows_path.resolve()),
            "summary": summarize_rows(rows, tensor_rel_tol=args.tensor_rel_tol),
        }

    base_resolution = args.base_resolution if args.base_resolution in results else sorted(results)[len(results) // 2]
    bx = results[base_resolution]["x"]
    bz = results[base_resolution]["z"]
    windows = select_windows(
        results[base_resolution]["rows"],
        tensor_rel_tol=args.tensor_rel_tol,
        focus_rho_min_frac=args.focus_rho_min_frac,
        width_x=args.window_width_x,
        width_z=args.window_width_z,
        top_count=args.top_windows,
        xmin=float(bx[0]),
        xmax=float(bx[-1]),
        zmin=float(bz[0]),
        zmax=float(bz[-1]),
    )

    for result in results.values():
        result["windows"] = summarize_windows(result["rows"], windows, tensor_rel_tol=args.tensor_rel_tol)

    highest_resolution = max(resolutions)
    fig_path = out / f"d_ls_local_resolution_convergence_ell{args.ell:g}_t{args.time:g}.png"
    render(
        fig_path,
        results=results,
        windows=windows,
        highest_resolution=highest_resolution,
        tensor_rel_tol=args.tensor_rel_tol,
    )

    serializable_results: dict[str, object] = {}
    for resolution, result in results.items():
        serializable_results[str(resolution)] = {
            "rows_jsonl": result["rows_jsonl"],
            "summary": result["summary"],
            "windows": result["windows"],
        }
    summary = {
        "params": {
            "branch": "D",
            "ell": args.ell,
            "time": args.time,
            "resolutions": resolutions,
            "base_resolution": base_resolution,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "levels_abs_raw_y": parse_positive_levels(args.levels),
            "half_width_y": args.half_width_y,
            "samples": args.samples,
            "tensor_rel_tol": args.tensor_rel_tol,
            "max_iter": args.max_iter,
            "multistart": args.multistart,
            "max_seeds": args.max_seeds,
            "focus_rho_min_frac": args.focus_rho_min_frac,
            "window_width_x": args.window_width_x,
            "window_width_z": args.window_width_z,
        },
        "definitions": {
            "purpose": "Checks whether least-squares unresolved direct-tensor segments disappear or converge under local/high resolution comparison.",
            "window_selection": "Windows are selected from the base resolution by binning unresolved rows above focus_rho_min_frac; they are diagnostic windows, not physical boundary conditions.",
            "accepted": "A row is accepted only when the full direct-tensor residual after least-squares is <= tensor_rel_tol.",
            "multistart": "If enabled, several covector seeds are tried for the same tensor residual. This is a solver diagnostic, not a new closure rule.",
            "no_closure_added": "No trace speed, damping, clipping, interpolation, or promotion of failed rows is used.",
        },
        "windows": [window.__dict__ for window in windows],
        "results": serializable_results,
        "files": {"figure": str(fig_path.resolve())},
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--resolutions", type=str, default="128,160,192")
    parser.add_argument("--base-resolution", type=int, default=160)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--levels", type=str, default="1")
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument("--tensor-rel-tol", type=float, default=0.1)
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--multistart", action="store_true")
    parser.add_argument("--max-seeds", type=int, default=24)
    parser.add_argument("--focus-rho-min-frac", type=float, default=1.0e-4)
    parser.add_argument("--window-width-x", type=float, default=4.0)
    parser.add_argument("--window-width-z", type=float, default=4.0)
    parser.add_argument("--top-windows", type=int, default=4)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
