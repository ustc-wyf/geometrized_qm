from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from adm_initial_data_abc import a_complex_state_from_rho_s_time_symmetric
from simulate_bc_from_a_initial_data import (
    B_STATE_KEYS,
    C_STATE_KEYS,
    advance_rk4,
    b_rhs_n,
    c_rhs_n,
    make_bc_initial_states,
)
from simulate_bc_geometry_from_a_reference import solve_linearized_conformal_omega
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    phase_field,
    spectral_omega,
)
from adm_sources_abc import a_branch_sources_from_complex


def make_a_flat_initial_data(
    rho0: np.ndarray,
    s0: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    one = np.ones_like(rho0)
    zero = np.zeros_like(rho0)
    a_init = a_complex_state_from_rho_s_time_symmetric(rho0, s0, one, zero, one, zero, mass, dx, dz)
    a_src = a_branch_sources_from_complex(
        a_init["phi1"],
        a_init["phi2"],
        a_init["pi1"],
        a_init["pi2"],
        one,
        zero,
        one,
        zero,
        mass,
        dx,
        dz,
    )
    return a_init, a_src


def state_metrics(state: dict[str, np.ndarray]) -> dict[str, float]:
    beta = state["beta"]
    h_xx = state["h_xx"]
    h_xz = state["h_xz"]
    h_zz = state["h_zz"]
    out = {
        "beta_min": float(np.nanmin(beta)),
        "beta_max": float(np.nanmax(beta)),
        "beta_abs_max": float(np.nanmax(np.abs(beta))),
        "lapse_min": float(np.nanmin(state["lapse"])) if "lapse" in state else 1.0,
        "lapse_max": float(np.nanmax(state["lapse"])) if "lapse" in state else 1.0,
        "h_xx_abs_max": float(np.nanmax(np.abs(h_xx))),
        "h_xz_abs_max": float(np.nanmax(np.abs(h_xz))),
        "h_zz_abs_max": float(np.nanmax(np.abs(h_zz))),
        "h_xx_minus_1_max": float(np.nanmax(np.abs(h_xx - 1.0))),
        "h_zz_minus_1_max": float(np.nanmax(np.abs(h_zz - 1.0))),
    }
    if "phi" in state:
        out["phi_minus_1_abs_max"] = float(np.nanmax(np.abs(state["phi"] - 1.0)))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="定位 B/C 长时间演化中的首次溢出或非有限值")
    parser.add_argument("--output", type=str, default="visualizations/bc_overflow_diagnosis", help="输出目录")
    parser.add_argument("--checkpoint-every", type=int, default=1000, help="检查点写入步数间隔")
    parser.add_argument("--gauge-mode", choices=["synchronous", "1plog"], default="synchronous", help="几何切片规范")
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
    omega = spectral_omega(x, z, loc.m)
    psi0_hat = np.fft.fft2(psi0)

    mass = 1.0
    mp = 300.0
    ell = 0.2
    elliptic_lambda = 1.0e-4
    dt = 2.5e-4
    steps = 64000
    finite_check_every = 100
    checkpoint_every = args.checkpoint_every

    _, a_src0 = make_a_flat_initial_data(rho0, s0, mass, dx, dz)
    a_weak0 = solve_linearized_conformal_omega(a_src0["energy"], mp, dx, dz)
    omega_init = a_weak0["omega"]

    b_state, c_state = make_bc_initial_states(
        rho0=rho0,
        s0=s0,
        omega_init=omega_init,
        mass=mass,
        mp=mp,
        ell=ell,
        dx=dx,
        dz=dz,
        c_phi_mode="elliptic",
        elliptic_lambda=elliptic_lambda,
    )

    params = {
        "dx": dx,
        "dz": dz,
        "m": mass,
        "mp": mp,
        "ell": ell,
        "gauge_mode": args.gauge_mode,
        "alpha": loc.alpha,
        "k0": loc.k0,
        "w0": np.nan,
        "phi0": loc.phi0,
        "theta": np.pi / 4.0,
        "x_half_range": loc.x_half_range,
        "z_half_range": loc.z_half_range,
        "nx": loc.nx,
        "nz": loc.nz,
        "initial_mode": "localized",
        "localized_psi0": psi0,
        "localized_omega": omega,
        "localized_psi0_hat": psi0_hat,
    }

    np.seterr(over="raise", invalid="raise")
    b_evolve_keys = B_STATE_KEYS + (("lapse",) if args.gauge_mode == "1plog" else ())
    c_evolve_keys = C_STATE_KEYS + (("lapse",) if args.gauge_mode == "1plog" else ())

    with log_path.open("w", encoding="utf-8") as f:
        t = 0.0
        try:
            for step in range(1, steps + 1):
                t_next = step * dt
                do_finite_check = finite_check_every > 0 and (step % finite_check_every == 0 or step == 1)
                _ = advance_rk4(b_state, b_evolve_keys, b_rhs_n, dt, params, t_next, do_finite_check=do_finite_check)
                _ = advance_rk4(c_state, c_evolve_keys, c_rhs_n, dt, params, t_next, do_finite_check=do_finite_check)
                t = t_next
                if step % checkpoint_every == 0 or step == 1:
                    rec = {
                        "step": step,
                        "t": t,
                        "b": state_metrics(b_state),
                        "c": state_metrics(c_state),
                    }
                    f.write(json.dumps(rec, ensure_ascii=False) + "\n")
                    f.flush()
            result = {"status": "completed", "step": steps, "t": t}
        except Exception as exc:  # noqa: BLE001
            result = {
                "status": "failed",
                "step": step,
                "t": t,
                "t_next": t_next,
                "exception_type": type(exc).__name__,
                "exception_message": str(exc),
                "b": state_metrics(b_state),
                "c": state_metrics(c_state),
            }

    result_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
