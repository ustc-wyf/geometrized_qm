from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from coordinate_matter_evolution import (
    conservative_density_from_rho_u,
    coordinate_matter_rhs_covector,
)
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_initial
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def rk4_step(
    n_cons: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov_txz: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    dt: float,
    u_t_reference: np.ndarray,
    active_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    def rhs(n0: np.ndarray, ux0: np.ndarray, uz0: np.ndarray, ut_ref: np.ndarray):
        n0_eff = np.where(active_mask, n0, 0.0)
        out = coordinate_matter_rhs_covector(
            n_cons=n0_eff,
            u_x=ux0,
            u_z=uz0,
            metric_cov_txz=metric_cov_txz,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=ut_ref,
        )
        return out["n_t"], out["u_x_t"], out["u_z_t"], out

    k1_n, k1_ux, k1_uz, o1 = rhs(n_cons, u_x, u_z, u_t_reference)
    k2_n, k2_ux, k2_uz, _ = rhs(
        n_cons + 0.5 * dt * k1_n,
        u_x + 0.5 * dt * k1_ux,
        u_z + 0.5 * dt * k1_uz,
        o1["u_t"],
    )
    k3_n, k3_ux, k3_uz, _ = rhs(
        n_cons + 0.5 * dt * k2_n,
        u_x + 0.5 * dt * k2_ux,
        u_z + 0.5 * dt * k2_uz,
        o1["u_t"],
    )
    k4_n, k4_ux, k4_uz, _ = rhs(
        n_cons + dt * k3_n,
        u_x + dt * k3_ux,
        u_z + dt * k3_uz,
        o1["u_t"],
    )

    n_next = n_cons + (dt / 6.0) * (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n)
    ux_next = u_x + (dt / 6.0) * (k1_ux + 2.0 * k2_ux + 2.0 * k3_ux + k4_ux)
    uz_next = u_z + (dt / 6.0) * (k1_uz + 2.0 * k2_uz + 2.0 * k3_uz + k4_uz)
    n_next = np.where(active_mask, n_next, 0.0)

    final = coordinate_matter_rhs_covector(
        n_cons=n_next,
        u_x=ux_next,
        u_z=uz_next,
        metric_cov_txz=metric_cov_txz,
        mass=mass,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=o1["u_t"],
    )
    return n_next, ux_next, uz_next, final


def run_case(
    output_dir: str | Path,
    steps: int = 80,
    dt: float = 2.5e-4,
    rho_floor: float = 1.0e-12,
) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    params = FlatLocalizedCrossingParams()
    x, z, X, Z = make_grid(params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    init = localized_direct_tilde_coordinate_initial(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        rho_floor=rho_floor,
    )
    metric_cov = init["metric_cov_txz_0"]
    d = init["diagnostics"]
    rho0 = init["rho0"]
    ux0 = d["s_x0"]
    uz0 = d["s_z0"]
    ut0 = d["s_t0"]

    cons0 = conservative_density_from_rho_u(
        rho=rho0,
        u_x=ux0,
        u_z=uz0,
        metric_cov_txz=metric_cov,
        mass=params.m,
        branch="negative_frequency",
        u_t_reference=ut0,
        rho_floor=rho_floor,
    )

    n = cons0["n_cons"].copy()
    ux = ux0.copy()
    uz = uz0.copy()
    current = coordinate_matter_rhs_covector(
        n_cons=n,
        u_x=ux,
        u_z=uz,
        metric_cov_txz=metric_cov,
        mass=params.m,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=ut0,
    )
    records: list[dict[str, float]] = []
    support = d["support_mask"]
    active_mask = support.copy()
    for _ in range(8):
        active_mask = (
            active_mask
            | np.roll(active_mask, 1, axis=0)
            | np.roll(active_mask, -1, axis=0)
            | np.roll(active_mask, 1, axis=1)
            | np.roll(active_mask, -1, axis=1)
        )
    n = np.where(active_mask, n, 0.0)

    for step in range(steps):
        n, ux, uz, current = rk4_step(
            n_cons=n,
            u_x=ux,
            u_z=uz,
            metric_cov_txz=metric_cov,
            mass=params.m,
            dx=dx,
            dz=dz,
            dt=dt,
            u_t_reference=current["u_t"],
            active_mask=active_mask,
        )
        if (step + 1) % 20 == 0 or step == steps - 1:
            records.append(
                {
                    "step": int(step + 1),
                    "t": float((step + 1) * dt),
                    "rho_min": float(np.nanmin(current["rho"])),
                    "rho_max": float(np.nanmax(current["rho"])),
                    "disc_min_support": float(np.nanmin(current["discriminant"][support])),
                    "disc_min_all": float(np.nanmin(current["discriminant"])),
                    "n_cons_max_abs": float(np.nanmax(np.abs(n))),
                    "ux_max_abs": float(np.nanmax(np.abs(ux))),
                    "uz_max_abs": float(np.nanmax(np.abs(uz))),
                }
            )

    rho_final = current["rho"]
    summary = {
        "steps": int(steps),
        "dt": float(dt),
        "t_final": float(steps * dt),
        "rho_rel_l1_vs_initial": float(np.nanmean(np.abs(rho_final - rho0)) / max(np.nanmean(np.abs(rho0)), 1.0e-30)),
        "rho_max_abs_vs_initial": float(np.nanmax(np.abs(rho_final - rho0))),
        "disc_min_all_final": float(np.nanmin(current["discriminant"])),
        "disc_min_support_final": float(np.nanmin(current["discriminant"][support])),
        "records": records,
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="固定实验室坐标、冻结初始 tilde g 的 B 支物质原型（推进 n_cons,u_x,u_z）")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--dt", type=float, default=2.5e-4)
    args = parser.parse_args()
    summary = run_case(args.output_dir, steps=args.steps, dt=args.dt)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
