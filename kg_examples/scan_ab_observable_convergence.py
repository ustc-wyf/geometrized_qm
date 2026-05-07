from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_ab_measurable_observables import extract_centerline_observables, simulate_final_ab


OBS_KEYS = [
    "central_peak_rel_intensity_shift",
    "mean_peak_rel_intensity_shift",
    "support_pointwise_rel_diff_max",
    "visibility_rel_diff",
    "effective_phase_shift_rad",
    "equivalent_fringe_shift_over_spacing",
]

OBS_LABELS = {
    "central_peak_rel_intensity_shift": "central peak",
    "mean_peak_rel_intensity_shift": "mean peak",
    "support_pointwise_rel_diff_max": "support max",
    "visibility_rel_diff": "visibility",
    "effective_phase_shift_rad": "phase (rad)",
    "equivalent_fringe_shift_over_spacing": "shift / spacing",
}


def collect_case(**sim_kwargs) -> dict[str, object]:
    sim = simulate_final_ab(**sim_kwargs)
    obs = extract_centerline_observables(sim["x"], sim["z"], sim["A"], sim["B"])
    a = sim["A"]
    b = sim["B"]
    return {
        "simulation": {
            "dt": sim["dt"],
            "steps": sim["steps"],
            "dt_target": sim_kwargs.get("dt_target", 0.025),
            "nx_outer": sim_kwargs.get("nx", 181),
            "nz_outer": sim_kwargs.get("nz", 181),
            "nx_inner": len(sim["x"]),
            "nz_inner": len(sim["z"]),
            "nkx": sim_kwargs.get("nkx", 61),
            "nkz": sim_kwargs.get("nkz", 61),
        },
        "field_metrics": {
            "B_vs_A_relL1": float(np.mean(np.abs(b - a)) / max(float(np.max(a)), 1.0e-14)),
            "A_max": float(np.max(a)),
            "B_max": float(np.max(b)),
        },
        "observables": {
            key: float(obs[key])
            for key in OBS_KEYS
        },
    }


def attach_deltas(cases: list[dict[str, object]]) -> None:
    ref_obs = cases[-1]["observables"]
    for case in cases:
        obs = case["observables"]
        deltas = {}
        signal_to_floor = {}
        for key in OBS_KEYS:
            delta = abs(float(obs[key]) - float(ref_obs[key]))
            deltas[key] = float(delta)
            if delta == 0.0:
                signal_to_floor[key] = None
            else:
                signal_to_floor[key] = float(abs(float(ref_obs[key])) / delta)
        case["delta_to_finest"] = deltas
        case["signal_to_delta_to_finest"] = signal_to_floor


def estimate_floor(cases: list[dict[str, object]]) -> dict[str, float]:
    if len(cases) < 2:
        return {key: 0.0 for key in OBS_KEYS}
    finest = cases[-1]["observables"]
    previous = cases[-2]["observables"]
    return {
        key: float(abs(float(finest[key]) - float(previous[key])))
        for key in OBS_KEYS
    }


def render_observable_scan(
    path: Path,
    xs: list[float],
    cases: list[dict[str, object]],
    xlabel: str,
    xscale: str = "log",
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), squeeze=False)
    plot_keys = [
        "central_peak_rel_intensity_shift",
        "mean_peak_rel_intensity_shift",
        "support_pointwise_rel_diff_max",
        "B_vs_A_relL1",
    ]

    for ax, key in zip(axes.ravel(), plot_keys):
        if key == "B_vs_A_relL1":
            ys = [float(case["field_metrics"][key]) for case in cases]
            label = "field relL1"
        else:
            ys = [float(case["observables"][key]) for case in cases]
            label = OBS_LABELS[key]
        ax.plot(xs, ys, marker="o", linewidth=1.7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.set_xscale(xscale)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def render_delta_scan(
    path: Path,
    xs: list[float],
    cases: list[dict[str, object]],
    xlabel: str,
    xscale: str = "log",
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4), squeeze=False)
    plot_keys = [
        "central_peak_rel_intensity_shift",
        "mean_peak_rel_intensity_shift",
        "support_pointwise_rel_diff_max",
        "visibility_rel_diff",
    ]

    for ax, key in zip(axes.ravel(), plot_keys):
        ys = [max(float(case["delta_to_finest"][key]), 1.0e-30) for case in cases]
        ax.plot(xs, ys, marker="o", linewidth=1.7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"|delta| to finest ({OBS_LABELS[key]})")
        ax.grid(alpha=0.3)
        ax.set_xscale(xscale)
        ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(outdir: str | Path) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    common = {
        "alpha": 0.5,
        "lambda_grav": 1.0,
        "mp": 300.0,
        "inner": 10.0,
        "outer": 15.0,
        "nkx": 81,
        "nkz": 81,
    }

    dt_cases: list[dict[str, object]] = []
    for dt_target in [0.05, 0.025, 0.0125, 0.00625]:
        dt_cases.append(
            collect_case(
                **common,
                dt_target=dt_target,
                nx=181,
                nz=181,
            )
        )
    attach_deltas(dt_cases)

    grid_cases: list[dict[str, object]] = []
    for n in [141, 181, 221]:
        grid_cases.append(
            collect_case(
                **common,
                dt_target=0.0125,
                nx=n,
                nz=n,
            )
        )
    attach_deltas(grid_cases)

    dt_xs = [float(case["simulation"]["dt"]) for case in dt_cases]
    dx_xs = [
        2.0 * common["outer"] / (int(case["simulation"]["nx_outer"]) - 1)
        for case in grid_cases
    ]

    render_observable_scan(out / "dt_scan_observables.png", dt_xs, dt_cases, "dt")
    render_delta_scan(out / "dt_scan_deltas.png", dt_xs, dt_cases, "dt")
    render_observable_scan(out / "dx_scan_observables.png", dx_xs, grid_cases, "dx")
    render_delta_scan(out / "dx_scan_deltas.png", dx_xs, grid_cases, "dx")

    summary = {
        "lambda_grav": 1.0,
        "dt_scan": dt_cases,
        "grid_scan": grid_cases,
        "estimated_numerical_floors": {
            "dt_scan": estimate_floor(dt_cases),
            "grid_scan": estimate_floor(grid_cases),
        },
        "notes": {
            "purpose": "Direct convergence scan of the lambda=1 AB observables to test whether the nonzero signal survives dt and grid refinement.",
            "main_signal_candidate": "central_peak_rel_intensity_shift",
            "reference_convention": "Each scan uses its finest entry as the local reference for delta-to-finest estimates.",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = run(Path(__file__).resolve().parent.parent / "visualizations" / "ab_observable_convergence")
    print(json.dumps(result, indent=2))
