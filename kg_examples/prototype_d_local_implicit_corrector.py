from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np

from prototype_d_reference_imex import (
    array_stats,
    build_refine_mask,
    imex_local_seed,
    reference_snapshot_data,
    trace_residual_2p1_from_chi,
)
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def local_implicit_correct_chi(
    ref: dict[str, np.ndarray],
    refine_mask: np.ndarray,
    ell: float,
    probe_dt: float,
    n_iter: int = 3,
    damping: float = 0.8,
) -> dict[str, np.ndarray | list[dict[str, float]]]:
    chi_m = np.array(ref["chi_m"], copy=True)
    chi_0 = np.array(ref["chi_ref"], copy=True)
    chi_p = np.array(ref["chi_p"], copy=True)
    metric_m = ref["metric_cov_m"]
    metric_0 = ref["metric_cov"]
    metric_p = ref["metric_cov_p"]
    stress_trace = ref["stress_trace"]
    dx = float(ref["dx"])
    dz = float(ref["dz"])

    refine_indices = np.argwhere(refine_mask)
    history: list[dict[str, float]] = []
    for it in range(n_iter):
        residual_data = trace_residual_2p1_from_chi(
            chi_m=chi_m,
            chi_0=chi_0,
            chi_p=chi_p,
            metric_m=metric_m,
            metric_0=metric_0,
            metric_p=metric_p,
            stress_trace=stress_trace,
            ell=ell,
            probe_dt=probe_dt,
            dx=dx,
            dz=dz,
        )
        current_residual = residual_data["trace_residual_2p1"]
        update = np.zeros_like(chi_0)
        for i, j in refine_indices:
            base = float(chi_0[i, j])
            f0 = float(current_residual[i, j])
            eps = max(1.0e-4, 1.0e-3 * max(abs(base), 1.0))

            chi_plus = np.array(chi_0, copy=True)
            chi_minus = np.array(chi_0, copy=True)
            chi_plus[i, j] = base + eps
            chi_minus[i, j] = base - eps

            f_plus = float(
                trace_residual_2p1_from_chi(
                    chi_m=chi_m,
                    chi_0=chi_plus,
                    chi_p=chi_p,
                    metric_m=metric_m,
                    metric_0=metric_0,
                    metric_p=metric_p,
                    stress_trace=stress_trace,
                    ell=ell,
                    probe_dt=probe_dt,
                    dx=dx,
                    dz=dz,
                )["trace_residual_2p1"][i, j]
            )
            f_minus = float(
                trace_residual_2p1_from_chi(
                    chi_m=chi_m,
                    chi_0=chi_minus,
                    chi_p=chi_p,
                    metric_m=metric_m,
                    metric_0=metric_0,
                    metric_p=metric_p,
                    stress_trace=stress_trace,
                    ell=ell,
                    probe_dt=probe_dt,
                    dx=dx,
                    dz=dz,
                )["trace_residual_2p1"][i, j]
            )
            jac = (f_plus - f_minus) / (2.0 * eps)
            if abs(jac) < 1.0e-10:
                continue
            raw_delta = -damping * f0 / jac
            raw_delta = float(np.clip(raw_delta, -0.25, 0.25))

            best_delta = 0.0
            best_abs = abs(f0)
            for scale in (1.0, 0.5, 0.25, 0.1):
                trial_delta = raw_delta * scale
                chi_trial = np.array(chi_0, copy=True)
                chi_trial[i, j] = base + trial_delta
                f_trial = float(
                    trace_residual_2p1_from_chi(
                        chi_m=chi_m,
                        chi_0=chi_trial,
                        chi_p=chi_p,
                        metric_m=metric_m,
                        metric_0=metric_0,
                        metric_p=metric_p,
                        stress_trace=stress_trace,
                        ell=ell,
                        probe_dt=probe_dt,
                        dx=dx,
                        dz=dz,
                    )["trace_residual_2p1"][i, j]
                )
                if abs(f_trial) < best_abs:
                    best_abs = abs(f_trial)
                    best_delta = trial_delta
            update[i, j] = best_delta

        chi_0 = chi_0 + update
        history.append(
            {
                "iter": float(it),
                "residual_support_p95": array_stats(current_residual, refine_mask)["p95"],
                "residual_support_abs_max": array_stats(current_residual, refine_mask)["abs_max"],
                "delta_support_p95": array_stats(update, refine_mask)["p95"],
                "delta_support_abs_max": array_stats(update, refine_mask)["abs_max"],
            }
        )

    final_residual = trace_residual_2p1_from_chi(
        chi_m=chi_m,
        chi_0=chi_0,
        chi_p=chi_p,
        metric_m=metric_m,
        metric_0=metric_0,
        metric_p=metric_p,
        stress_trace=stress_trace,
        ell=ell,
        probe_dt=probe_dt,
        dx=dx,
        dz=dz,
    )
    return {
        "chi_corrected": np.real_if_close(chi_0),
        "trace_residual_corrected": np.real_if_close(final_residual["trace_residual_2p1"]),
        "phi_corrected": np.real_if_close(final_residual["phi"]),
        "history": history,
    }


def render_comparison(
    outdir: Path,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    ref: dict[str, np.ndarray],
    masks: dict[str, np.ndarray | float],
    corrected: dict[str, np.ndarray | list[dict[str, float]]],
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    panels = [
        (
            ref["chi_ref"],
            r"$\chi_{\rm ref} = \tilde R_{\rm ref}$",
            "Reference signed curvature scalar built from the reference metric.",
            "coolwarm",
        ),
        (
            corrected["chi_corrected"],
            r"$\chi_{\rm corr}$",
            "Signed curvature proxy after local implicit correction.",
            "coolwarm",
        ),
        (
            corrected["chi_corrected"] - ref["chi_ref"],
            r"$\delta\chi = \chi_{\rm corr}-\chi_{\rm ref}$",
            "Local implicit correction applied to the signed curvature proxy.",
            "coolwarm",
        ),
        (
            ref["trace_residual_2p1"],
            r"$\mathcal{T}_{\rm ref}^{(2+1)}$",
            "Reference 2+1d trace residual of the D-branch scalar-tensor field equation.",
            "coolwarm",
        ),
        (
            corrected["trace_residual_corrected"],
            r"$\mathcal{T}_{\rm corr}^{(2+1)}$",
            "Trace residual after the local implicit correction.",
            "coolwarm",
        ),
        (
            np.abs(corrected["trace_residual_corrected"]) - np.abs(ref["trace_residual_2p1"]),
            r"$|\mathcal{T}_{\rm corr}|-|\mathcal{T}_{\rm ref}|$",
            "Change in residual magnitude; negative values mean local improvement.",
            "coolwarm",
        ),
    ]
    description = (
        "White contour = refine mask = support & "
        "(near_transition or residual>p99 or |grad phi|>p99 or derivative_norm>p99).\n"
        "support means rho > 1e-3 * rho_max. "
        "chi is the signed curvature proxy and T is the 2+1d trace residual."
    )

    def draw_figure(use_symlog: bool, suffix: str) -> None:
        fig, axes = plt.subplots(2, 3, figsize=(15.0, 9.6), constrained_layout=True)
        for ax, (field, title, subtitle, cmap) in zip(axes.reshape(-1), panels):
            if use_symlog:
                max_abs = float(np.nanmax(np.abs(field)))
                max_abs = max(max_abs, 1.0e-12)
                linthresh = max(1.0e-8, 1.0e-3 * max_abs)
                norm = SymLogNorm(linthresh=linthresh, vmin=-max_abs, vmax=max_abs)
                im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap, norm=norm)
            else:
                im = ax.pcolormesh(xg, zg, field, shading="auto", cmap=cmap)
            ax.contour(
                xg,
                zg,
                masks["refine_mask"].astype(float),
                levels=[0.5],
                colors="white",
                linewidths=0.9,
            )
            ax.set_title(title, fontsize=10)
            ax.text(
                0.02,
                0.02,
                subtitle,
                transform=ax.transAxes,
                fontsize=7.5,
                color="black",
                bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 1.8},
                va="bottom",
                ha="left",
            )
            ax.set_xlabel("x")
            ax.set_ylabel("z")
            ax.set_aspect("equal")
            fig.colorbar(im, ax=ax, shrink=0.86)
        scale_label = "symlog color scale" if use_symlog else "linear color scale"
        fig.suptitle(
            f"D branch local implicit chi corrector, t={t:.2f} ({scale_label})",
            fontsize=13,
        )
        fig.text(0.02, 0.985, description, fontsize=8.5, va="top")
        fig.savefig(outdir / f"d_local_corrector_t{t:.1f}{suffix}.png", dpi=180)
        plt.close(fig)

    draw_figure(use_symlog=False, suffix="_linear")
    draw_figure(use_symlog=True, suffix="_symlog")


def run(
    outdir: str | Path,
    ell: float = 100.0,
    mp: float = 300.0,
    rho_floor: float = 1.0e-12,
    probe_dt: float = 1.0e-3,
    nx: int = 96,
    nz: int = 96,
    times: list[float] | None = None,
) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = FlatLocalizedCrossingParams(nx=nx, nz=nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    times = [8.0] if times is None else times
    snapshots: list[dict[str, object]] = []
    for t in times:
        ref = reference_snapshot_data(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t=t,
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
        corrected = local_implicit_correct_chi(
            ref=ref,
            refine_mask=masks["refine_mask"],
            ell=ell,
            probe_dt=probe_dt,
        )
        render_comparison(out, x, z, t, ref, masks, corrected)

        before_stats = array_stats(ref["trace_residual_2p1"], masks["support_mask"])
        after_stats = array_stats(corrected["trace_residual_corrected"], masks["support_mask"])
        before_refine_stats = array_stats(ref["trace_residual_2p1"], masks["refine_mask"])
        after_refine_stats = array_stats(corrected["trace_residual_corrected"], masks["refine_mask"])
        snapshots.append(
            {
                "time": float(t),
                "support_trace_before": before_stats,
                "support_trace_after": after_stats,
                "refine_trace_before": before_refine_stats,
                "refine_trace_after": after_refine_stats,
                "refine_fraction_support": float(
                    np.mean(masks["refine_mask"][masks["support_mask"]]) if np.any(masks["support_mask"]) else 0.0
                ),
                "history": corrected["history"],
            }
        )

    summary = {
        "params": {
            "ell": ell,
            "mp": mp,
            "rho_floor": rho_floor,
            "probe_dt": probe_dt,
            "grid_nx": int(len(x)),
            "grid_nz": int(len(z)),
            "times": times,
            "definition_support": "rho > 1e-3 * rho_max",
            "definition_refine_mask": "support & (near_transition or residual>p99 or |grad phi|>p99 or derivative_norm>p99)",
            "corrector": "frozen-geometry local implicit chi corrector on refine mask only",
        },
        "snapshots": snapshots,
        "files": {
            "maps_linear": [str((out / f"d_local_corrector_t{t:.1f}_linear.png").resolve()) for t in times],
            "maps_symlog": [str((out / f"d_local_corrector_t{t:.1f}_symlog.png").resolve()) for t in times],
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
            / "d_local_implicit_corrector_96"
        ),
    )
    parser.add_argument("--ell", type=float, default=100.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--probe-dt", type=float, default=1.0e-3)
    parser.add_argument("--nx", type=int, default=96)
    parser.add_argument("--nz", type=int, default=96)
    parser.add_argument("--times", type=float, nargs="*", default=None)
    args = parser.parse_args()
    summary = run(
        outdir=args.outdir,
        ell=args.ell,
        mp=args.mp,
        rho_floor=args.rho_floor,
        probe_dt=args.probe_dt,
        nx=args.nx,
        nz=args.nz,
        times=args.times,
    )
    print(json.dumps(summary, indent=2))
