from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from physical_units import HBAR_C_EV_M
from simulate_d_local_window_from_a_snapshot import build_physical_reference, make_cropped_snapshot, real_array


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    ref = build_physical_reference(args)
    params = ref["params"]
    scale = ref["scale"]
    old_scale = float(scale.old_dimensionless_scale_ev_inv)
    t_meet_old = float(params.t_meet / old_scale)
    x_full = np.asarray(ref["x_full"], dtype=float)
    z_full = np.asarray(ref["z_full"], dtype=float)
    window_ev_inv = float(args.window_um) * 1.0e-6 / HBAR_C_EV_M
    ix = np.where(np.abs(x_full) <= window_ev_inv)[0]
    iz = np.where(np.abs(z_full) <= window_ev_inv)[0]
    x = x_full[ix]
    z = z_full[iz]
    x_um = x * HBAR_C_EV_M * 1.0e6
    z_um = z * HBAR_C_EV_M * 1.0e6
    dx_um = float((x[1] - x[0]) * HBAR_C_EV_M * 1.0e6)
    wavelength_um = float(args.wavelength_nm) * 1.0e-3
    taus = [float(t.strip()) for t in args.taus.split(",") if t.strip()]
    rhos = []
    stats = []
    for tau in taus:
        t = (t_meet_old + tau) * old_scale
        snap = make_cropped_snapshot(
            psi0=np.asarray(ref["psi0"]),
            psi0_hat=np.asarray(ref["psi0_hat"]),
            omega=np.asarray(ref["omega"]),
            x_full=x_full,
            z_full=z_full,
            ix=ix,
            iz=iz,
            t=t,
            rho_floor=float(args.rho_floor),
            x_floor=float(args.x_floor),
            pinv_rcond=float(args.pinv_rcond),
            with_metric=False,
        )
        rho = real_array(snap["bohm"]["rho"])
        rhos.append(rho)
        edge = np.concatenate([rho[0, :], rho[-1, :], rho[:, 0], rho[:, -1]])
        stats.append(
            {
                "tau": tau,
                "rho_max": float(np.max(rho)),
                "rho_edge_over_max": float(np.max(edge) / max(float(np.max(rho)), 1.0e-300)),
            }
        )

    vmax = max(float(np.percentile(rho, 99.5)) for rho in rhos)
    fig, axes = plt.subplots(1, len(taus), figsize=(5.1 * len(taus), 4.7), constrained_layout=True)
    if len(taus) == 1:
        axes = [axes]
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    for ax, tau, rho in zip(axes, taus, rhos):
        im = ax.pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
        ax.set_title(f"A branch rho, tau={tau:g}")
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    plot = args.output / "a_branch_rho_resolution_check.png"
    fig.savefig(plot, dpi=180)
    plt.close(fig)

    report = {
        "parameters": {
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "local_shape": [int(len(ix)), int(len(iz))],
            "dx_um": dx_um,
            "wavelength_um": wavelength_um,
            "points_per_wavelength": float(wavelength_um / dx_um),
            "points_per_half_wavelength": float(0.5 * wavelength_um / dx_um),
            "taus": taus,
            "kg_norm_after": float(ref["kg_norm_after"]),
        },
        "stats": stats,
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
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
