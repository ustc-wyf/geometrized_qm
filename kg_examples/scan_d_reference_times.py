from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from prototype_d_reference_imex import (
    build_refine_mask,
    imex_local_seed,
    reference_snapshot_data,
)
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def run(
    outdir: str | Path,
    ell: float = 100.0,
    mp: float = 300.0,
    rho_floor: float = 1.0e-12,
    probe_dt: float = 1.0e-3,
    nx: int = 64,
    nz: int = 64,
    n_times: int = 9,
    t_final: float = 16.0,
) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = FlatLocalizedCrossingParams(nx=nx, nz=nz, t_final=t_final)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    times = np.linspace(0.0, t_final, n_times)
    rows: list[dict[str, float]] = []

    for t in times:
        ref = reference_snapshot_data(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t=float(t),
            probe_dt=probe_dt,
            ell=ell,
            mp=mp,
            mass=params.m,
            rho_floor=rho_floor,
        )
        masks = build_refine_mask(
            rho=ref["rho"],
            r_tilde=ref["R_tilde"],
            phi_ref=ref["phi_ref"],
            residual_norm=ref["residual_norm"],
            derivative_norm=ref["derivative_norm"],
            dx=ref["dx"],
            dz=ref["dz"],
            ell=ell,
        )
        seed = imex_local_seed(
            chi_ref=ref["chi_ref"],
            r_tilde=ref["R_tilde"],
            trace_residual=ref["trace_residual_2p1"],
            derivative_norm=ref["derivative_norm"],
            ell=ell,
            dx=ref["dx"],
            dz=ref["dz"],
        )
        support = masks["support_mask"]
        support_count = int(np.count_nonzero(support))
        refine_count = int(np.count_nonzero(masks["refine_mask"] & support))
        near_count = int(np.count_nonzero(masks["near_transition_mask"] & support))
        row = {
            "time": float(t),
            "support_residual_p95": float(masks["residual_p95_support"]),
            "support_residual_p99": float(masks["residual_p99_support"]),
            "support_residual_abs_max": float(np.max(ref["residual_norm"][support])) if support_count else 0.0,
            "support_trace_p95": float(np.percentile(np.abs(ref["trace_residual_2p1"][support]), 95.0)) if support_count else 0.0,
            "support_trace_abs_max": float(np.max(np.abs(ref["trace_residual_2p1"][support]))) if support_count else 0.0,
            "support_dt_global_suggested": float(seed["dt_global_suggested"]),
            "support_refine_fraction": float(refine_count / support_count) if support_count else 0.0,
            "support_near_transition_fraction": float(near_count / support_count) if support_count else 0.0,
        }
        rows.append(row)

    residual_p95 = np.array([row["support_residual_p95"] for row in rows])
    trace_p95 = np.array([row["support_trace_p95"] for row in rows])
    refine_frac = np.array([row["support_refine_fraction"] for row in rows])
    dt_suggest = np.array([row["support_dt_global_suggested"] for row in rows])
    dt_log = -np.log10(np.maximum(dt_suggest, 1.0e-30))

    def zscore(arr: np.ndarray) -> np.ndarray:
        std = float(np.std(arr))
        if std < 1.0e-30:
            return np.zeros_like(arr)
        return (arr - float(np.mean(arr))) / std

    severity = zscore(residual_p95) + 0.7 * zscore(trace_p95) + 0.5 * zscore(refine_frac) + 0.5 * zscore(dt_log)
    for row, score in zip(rows, severity):
        row["severity_score"] = float(score)

    rows_sorted = sorted(rows, key=lambda r: r["severity_score"], reverse=True)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), constrained_layout=True)
    axes = axes.reshape(-1)
    axes[0].plot(times, residual_p95, marker="o", lw=1.6)
    axes[0].set_title("support residual p95")
    axes[1].plot(times, trace_p95, marker="o", lw=1.6)
    axes[1].set_title("support trace residual p95")
    axes[2].plot(times, refine_frac, marker="o", lw=1.6)
    axes[2].set_title("support refine fraction")
    axes[3].plot(times, dt_suggest, marker="o", lw=1.6)
    axes[3].set_yscale("log")
    axes[3].set_title("suggested global dt")
    for ax in axes:
        ax.set_xlabel("t")
        ax.grid(True, alpha=0.3)
    fig.savefig(out / "d_time_scan_metrics.png", dpi=180)
    plt.close(fig)

    summary = {
        "params": {
            "ell": ell,
            "mp": mp,
            "rho_floor": rho_floor,
            "probe_dt": probe_dt,
            "grid_nx": nx,
            "grid_nz": nz,
            "n_times": n_times,
            "t_final": t_final,
            "definition_support": "rho > 1e-3 * rho_max",
            "definition_refine_mask": "support & (near_transition or residual>p99 or |grad phi|>p99 or derivative_norm>p99)",
            "severity_score": "z(residual_p95) + 0.7 z(trace_p95) + 0.5 z(refine_fraction) + 0.5 z(-log10(dt)))",
        },
        "rows": rows,
        "worst_times_by_severity": rows_sorted[: min(5, len(rows_sorted))],
        "files": {
            "metrics_plot": str((out / "d_time_scan_metrics.png").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(
            Path(__file__).resolve().parent.parent
            / "visualizations"
            / "d_time_scan_reference_64"
        ),
    )
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--probe-dt", type=float, default=1.0e-3)
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--nz", type=int, default=64)
    parser.add_argument("--n-times", type=int, default=9)
    parser.add_argument("--t-final", type=float, default=16.0)
    args = parser.parse_args()
    summary = run(
        outdir=args.outdir,
        ell=args.ell,
        mp=args.mp,
        rho_floor=args.rho_floor,
        probe_dt=args.probe_dt,
        nx=args.nx,
        nz=args.nz,
        n_times=args.n_times,
        t_final=args.t_final,
    )
    print(json.dumps(summary, indent=2))
