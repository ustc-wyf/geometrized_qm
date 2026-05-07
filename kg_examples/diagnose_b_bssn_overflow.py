from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from simulate_b_bssn_driver import (
    adm_to_bssn_state,
    advance_rk4,
    b_bssn_rhs,
    state_metrics,
)
from simulate_bc_from_a_initial_data import make_bc_initial_states
from simulate_bc_geometry_from_a_reference import solve_linearized_conformal_omega
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    phase_field,
    spectral_omega,
)
from adm_initial_data_abc import a_complex_state_from_rho_s_time_symmetric
from adm_sources_abc import a_branch_sources_from_complex


def main() -> None:
    parser = argparse.ArgumentParser(description="B 支 BSSN 型长时间溢出诊断")
    parser.add_argument("--output", type=str, default="visualizations/b_bssn_overflow_diagnosis")
    parser.add_argument("--checkpoint-every", type=int, default=1000)
    parser.add_argument("--steps", type=int, default=64000)
    parser.add_argument("--dt", type=float, default=2.5e-4)
    parser.add_argument("--finite-check-every", type=int, default=100)
    parser.add_argument("--mu-shift", type=float, default=0.5)
    parser.add_argument("--eta-shift", type=float, default=0.5)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    out_dir = root / args.output
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "checkpoints.jsonl"
    result_path = out_dir / "result.json"

    loc = FlatLocalizedCrossingParams(
        m=1.0,
        alpha=0.5,
        k0=6.0,
        sigma_parallel=1.6,
        sigma_perp=1.2,
        phi0=0.0,
        x_a0=-6.0,
        z_a0=-6.0,
        x_b0=6.0,
        z_b0=-6.0,
        x_half_range=20.0,
        z_half_range=20.0,
        nx=256,
        nz=256,
        t_final=16.0,
        n_frames=7,
    )
    x, z, X, Z = make_grid(loc)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(X, Z, loc)
    rho0 = np.abs(psi0) ** 2 + loc.rho_floor
    s0 = phase_field(psi0)
    omega_loc = spectral_omega(x, z, loc.m)

    one = np.ones_like(rho0)
    zero = np.zeros_like(rho0)
    a_init = a_complex_state_from_rho_s_time_symmetric(rho0, s0, one, zero, one, zero, 1.0, dx, dz)
    a_src = a_branch_sources_from_complex(
        a_init["phi1"], a_init["phi2"], a_init["pi1"], a_init["pi2"], one, zero, one, zero, 1.0, dx, dz
    )
    omega_init = solve_linearized_conformal_omega(a_src["energy"], 300.0, dx, dz)["omega"]
    b_adm, _ = make_bc_initial_states(rho0, s0, omega_init, 1.0, 300.0, 0.2, dx, dz, c_phi_mode="flat")
    state = adm_to_bssn_state(b_adm)

    params = {
        "dx": dx,
        "dz": dz,
        "m": 1.0,
        "mp": 300.0,
        "mu_shift": args.mu_shift,
        "eta_shift": args.eta_shift,
        "alpha": loc.alpha,
        "k0": loc.k0,
        "phi0": loc.phi0,
        "theta": np.pi / 4.0,
        "x_half_range": loc.x_half_range,
        "z_half_range": loc.z_half_range,
        "nx": loc.nx,
        "nz": loc.nz,
        "initial_mode": "localized",
        "localized_psi0": psi0,
        "localized_omega": omega_loc,
        "localized_psi0_hat": np.fft.fft2(psi0),
    }
    evolve_keys = (
        "chi", "gt_xx", "gt_xz", "gt_zz", "k_trace",
        "at_xx", "at_xz", "at_zz", "lapse", "shift_x", "shift_z", "n", "s",
    )

    np.seterr(over="raise", invalid="raise")

    with log_path.open("w", encoding="utf-8") as f:
        t = 0.0
        try:
            for step in range(1, args.steps + 1):
                t_next = step * args.dt
                do_finite_check = args.finite_check_every > 0 and (
                    step % args.finite_check_every == 0 or step == 1
                )
                rhs = advance_rk4(
                    state=state,
                    evolve_keys=evolve_keys,
                    dt=args.dt,
                    params=params,
                    t_next=t_next,
                    do_finite_check=do_finite_check,
                )
                t = t_next
                if step % args.checkpoint_every == 0 or step == 1:
                    rec = {
                        "step": step,
                        "t": t,
                        **state_metrics(state),
                        "r3_abs_max": float(np.max(np.abs(rhs["r3_scalar"]))),
                        "rho_max": float(np.max(rhs["rho_current"])),
                        "gamma_driver_x_abs_max": float(np.max(np.abs(rhs["gamma_driver_x"]))),
                        "gamma_driver_z_abs_max": float(np.max(np.abs(rhs["gamma_driver_z"]))),
                    }
                    f.write(json.dumps(rec, ensure_ascii=False) + "\n")
                    f.flush()
            result = {"status": "completed", "step": args.steps, "t": t, **state_metrics(state)}
        except Exception as exc:  # noqa: BLE001
            result = {
                "status": "failed",
                "step": step,
                "t": t,
                "t_next": t_next,
                "exception_type": type(exc).__name__,
                "exception_message": str(exc),
                **state_metrics(state),
            }

    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
