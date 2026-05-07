from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from plot_bcd_transition_layers import make_fields
from scan_cd_ell_window import branch_f_derivatives
from scan_cd_resolution_left_arc import LEFT_ARC_BBOX, support_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def fraction(mask: np.ndarray, within: np.ndarray) -> float:
    denom = int(np.sum(within))
    if denom == 0:
        return 0.0
    return float(np.sum(mask & within) / denom)


def safe_ratio_stats(numerator: np.ndarray, denominator: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    denom = np.maximum(np.abs(denominator), 1.0e-300)
    return support_stats(np.abs(numerator) / denom, mask)


def signed_stats(field: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    vals = np.asarray(field[mask], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"mean": 0.0, "median": 0.0, "p05": 0.0, "p95": 0.0, "sum": 0.0}
    return {
        "mean": float(np.mean(vals)),
        "median": float(np.median(vals)),
        "p05": float(np.percentile(vals, 5.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "sum": float(np.sum(vals)),
    }


def signed_to_abs_sum_ratio(field: np.ndarray, mask: np.ndarray) -> float:
    vals = np.asarray(field[mask], dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return 0.0
    return float(np.abs(np.sum(vals)) / max(np.sum(np.abs(vals)), 1.0e-300))


def spatial_derivatives(field: np.ndarray, x: np.ndarray, z: np.ndarray) -> dict[str, np.ndarray]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    fx, fz = np.gradient(field, dx, dz, edge_order=2)
    fxx, fxz = np.gradient(fx, dx, dz, edge_order=2)
    fzx, fzz = np.gradient(fz, dx, dz, edge_order=2)
    return {
        "fx": fx,
        "fz": fz,
        "fxx": fxx,
        "fxz": 0.5 * (fxz + fzx),
        "fzz": fzz,
    }


def normal_layer_quantities(phi: np.ndarray, y: np.ndarray, x: np.ndarray, z: np.ndarray) -> dict[str, np.ndarray]:
    y_der = spatial_derivatives(y, x, z)
    phi_der = spatial_derivatives(phi, x, z)

    grad_y_norm = np.sqrt(y_der["fx"] ** 2 + y_der["fz"] ** 2)
    nx = y_der["fx"] / np.maximum(grad_y_norm, 1.0e-300)
    nz = y_der["fz"] / np.maximum(grad_y_norm, 1.0e-300)

    dphi_dn = phi_der["fx"] * nx + phi_der["fz"] * nz
    d2phi_dnn = (
        nx * nx * phi_der["fxx"]
        + 2.0 * nx * nz * phi_der["fxz"]
        + nz * nz * phi_der["fzz"]
    )

    # Local proxy for integrating across an O(1) interval in y=ell^2 R.
    # It is not a replacement for a resolved line integral; it flags whether
    # a thin layer is likely to have finite signed jump or only large absolute strength.
    width_proxy = 1.0 / np.maximum(grad_y_norm, 1.0e-300)
    signed_jump_proxy = d2phi_dnn * width_proxy
    abs_layer_strength_proxy = np.abs(d2phi_dnn) * width_proxy

    return {
        "grad_y_norm": grad_y_norm,
        "dphi_dn": dphi_dn,
        "d2phi_dnn": d2phi_dnn,
        "width_proxy": width_proxy,
        "signed_jump_proxy": signed_jump_proxy,
        "abs_layer_strength_proxy": abs_layer_strength_proxy,
    }


def summarize(branch: str, ell: float, fields: dict[str, np.ndarray], x: np.ndarray, z: np.ndarray, mp: float) -> dict[str, object]:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    xmin, xmax, zmin, zmax = LEFT_ARC_BBOX
    window = (xg >= xmin) & (xg <= xmax) & (zg >= zmin) & (zg <= zmax)

    rho = np.real_if_close(fields["rho"])
    r_tilde = np.real_if_close(fields["R_tilde"])
    rhs_norm_scaled = np.real_if_close(fields["rhs_norm"])
    algebraic_norm_scaled = np.real_if_close(fields["algebraic_norm"])
    derivative_norm_scaled = np.real_if_close(fields["derivative_norm"])
    residual_norm_scaled = np.real_if_close(fields["residual_norm"])

    support = rho > 1.0e-3 * float(np.max(rho))
    dense = rho > 1.0e-1 * float(np.max(rho))
    win_support = window & support
    win_dense = window & dense

    derivs = branch_f_derivatives(branch, r_tilde, ell)
    y = derivs["y"]
    f = np.real_if_close(derivs["f"])
    f_r = np.real_if_close(derivs["fR"])

    strict = (np.abs(y) >= 0.5) & (np.abs(y) <= 2.0)
    broad = (np.abs(y) >= 1.0) & (np.abs(y) <= 100.0)

    layer = normal_layer_quantities(f_r, y, x, z)

    mp2 = mp * mp
    t_norm = mp2 * rhs_norm_scaled
    algebraic_norm = mp2 * algebraic_norm_scaled
    derivative_norm = mp2 * derivative_norm_scaled
    residual_norm = mp2 * residual_norm_scaled
    effective_planck_m2 = mp2 * f_r
    saturated_cosmological_scale = 0.5 * mp2 * np.abs(f)

    masks = {
        "window_support": win_support,
        "window_dense": win_dense,
        "window_strict_support": window & support & strict,
        "window_broad_support": window & support & broad,
    }

    records: dict[str, object] = {}
    for name, mask in masks.items():
        records[name] = {
            "count": int(np.sum(mask)),
            "fR": support_stats(f_r, mask),
            "M_P2_fR": support_stats(effective_planck_m2, mask),
            "M_P2_abs_f_over_2": support_stats(saturated_cosmological_scale, mask),
            "T_norm": support_stats(t_norm, mask),
            "M_P2_algebraic_norm": support_stats(algebraic_norm, mask),
            "M_P2_derivative_norm": support_stats(derivative_norm, mask),
            "M_P2_residual_norm": support_stats(residual_norm, mask),
            "algebraic_to_T_ratio": safe_ratio_stats(algebraic_norm, t_norm, mask),
            "derivative_to_T_ratio": safe_ratio_stats(derivative_norm, t_norm, mask),
            "cosmological_to_T_ratio": safe_ratio_stats(saturated_cosmological_scale, t_norm, mask),
            "dphi_dn": support_stats(layer["dphi_dn"], mask),
            "d2phi_dnn": support_stats(layer["d2phi_dnn"], mask),
            "width_proxy": support_stats(layer["width_proxy"], mask),
            "signed_jump_proxy": support_stats(layer["signed_jump_proxy"], mask),
            "signed_jump_proxy_signed": signed_stats(layer["signed_jump_proxy"], mask),
            "signed_to_abs_jump_proxy_sum_ratio": signed_to_abs_sum_ratio(layer["signed_jump_proxy"], mask),
            "abs_layer_strength_proxy": support_stats(layer["abs_layer_strength_proxy"], mask),
        }

    return {
        "branch": branch,
        "ell": float(ell),
        "bbox": list(LEFT_ARC_BBOX),
        "counts": {
            "support": int(np.sum(support)),
            "window_support": int(np.sum(win_support)),
            "window_strict_support": int(np.sum(window & support & strict)),
            "window_broad_support": int(np.sum(window & support & broad)),
        },
        "fractions_window_support": {
            "strict_transition": fraction(strict, win_support),
            "broad_transition": fraction(broad, win_support),
            "fR_lt_1e-4": fraction(f_r < 1.0e-4, win_support),
            "M_P2_fR_lt_1": fraction(effective_planck_m2 < 1.0, win_support),
            "M_P2_fR_lt_Tnorm_p95_global": fraction(
                effective_planck_m2 < np.percentile(t_norm[support], 95.0) if np.any(support) else False,
                win_support,
            ),
        },
        "statistics": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ells", type=str, default="10,30,100,300")
    parser.add_argument("--branches", type=str, default="D")
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--nx", type=int, default=160)
    parser.add_argument("--nz", type=int, default=160)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    args = parser.parse_args()

    ells = [float(tok) for tok in args.ells.split(",") if tok.strip()]
    branches = [tok.strip() for tok in args.branches.split(",") if tok.strip()]

    params = FlatLocalizedCrossingParams(nx=args.nx, nz=args.nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    records: list[dict[str, object]] = []
    for branch in branches:
        for ell in ells:
            fields = make_fields(
                branch=branch,
                psi0=psi0,
                psi0_hat=psi0_hat,
                omega=omega,
                x=x,
                z=z,
                t=args.time,
                probe_dt=args.probe_dt,
                ell=ell,
                mp=args.mp,
                mass=params.m,
                rho_floor=args.rho_floor,
            )
            records.append(summarize(branch=branch, ell=ell, fields=fields, x=x, z=z, mp=args.mp))

    summary = {
        "params": {
            "ells": ells,
            "branches": branches,
            "time": args.time,
            "nx": args.nx,
            "nz": args.nz,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
        },
        "definitions": {
            "support": "rho > 1e-3 * rho_max",
            "window": "left arc window: -7.4 <= x <= -2.2, 2.6 <= z <= 6.2",
            "strict_transition": "0.5 <= |ell^2 R_tilde| <= 2.0",
            "broad_transition": "1 <= |ell^2 R_tilde| <= 100",
            "dphi_dn": "spatial normal derivative of phi=f_R, with normal along grad(ell^2 R_tilde)",
            "d2phi_dnn": "spatial second derivative of phi=f_R along the same normal; ignores derivative-of-normal terms",
            "width_proxy": "1 / |grad(ell^2 R_tilde)|",
            "signed_jump_proxy": "d2phi_dnn * width_proxy; local proxy, not a resolved line integral",
            "abs_layer_strength_proxy": "|d2phi_dnn| * width_proxy; local absolute thin-layer strength proxy",
        },
        "records": records,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    compact = {
        f"{r['branch']}_ell_{r['ell']:g}": {
            "window_strict_support": r["counts"]["window_strict_support"],
            "window_broad_support": r["counts"]["window_broad_support"],
            "fR_p95_window": r["statistics"]["window_support"]["fR"]["p95"],
            "MP2_fR_p95_window": r["statistics"]["window_support"]["M_P2_fR"]["p95"],
            "Tnorm_p95_window": r["statistics"]["window_support"]["T_norm"]["p95"],
            "derivative_to_T_p95_broad": r["statistics"]["window_broad_support"]["derivative_to_T_ratio"]["p95"],
            "abs_layer_strength_p95_strict": r["statistics"]["window_strict_support"]["abs_layer_strength_proxy"]["p95"],
            "signed_jump_proxy_p95_strict": r["statistics"]["window_strict_support"]["signed_jump_proxy"]["p95"],
            "signed_to_abs_jump_ratio_strict": r["statistics"]["window_strict_support"]["signed_to_abs_jump_proxy_sum_ratio"],
        }
        for r in records
    }
    print(json.dumps(compact, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
