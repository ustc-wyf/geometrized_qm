from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from plot_bcd_transition_layers import make_fields
from scan_cd_resolution_left_arc import support_stats
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def branch_f_derivatives(branch: str, r_tilde: np.ndarray, ell: float) -> dict[str, np.ndarray]:
    l2 = ell * ell
    y = l2 * r_tilde
    if branch == "C":
        denom = 1.0 + y * y
        f = r_tilde / np.sqrt(denom)
        fr = np.power(denom, -1.5)
        frr = -3.0 * l2 * y * np.power(denom, -2.5)
        return {"y": y, "f": f, "fR": fr, "fRR": frr}
    if branch == "D":
        tanh_y = np.tanh(y)
        fr = 1.0 - tanh_y * tanh_y
        f = tanh_y / l2
        frr = -2.0 * l2 * fr * tanh_y
        return {"y": y, "f": f, "fR": fr, "fRR": frr}
    raise ValueError(f"Unsupported branch: {branch}")


def signed_scalaron_mass_squared(
    r_tilde: np.ndarray,
    fr: np.ndarray,
    frr: np.ndarray,
    frr_floor: float,
) -> np.ndarray:
    out = np.full_like(r_tilde, np.nan, dtype=float)
    mask = np.abs(frr) > frr_floor
    out[mask] = (fr[mask] - r_tilde[mask] * frr[mask]) / (3.0 * frr[mask])
    return out


def finite_fraction(mask: np.ndarray, within: np.ndarray) -> float:
    denom = int(np.sum(within))
    if denom == 0:
        return 0.0
    return float(np.sum(mask & within) / denom)


def finite_sign_fractions(field: np.ndarray, within: np.ndarray) -> dict[str, float]:
    finite = within & np.isfinite(field)
    denom = int(np.sum(finite))
    if denom == 0:
        return {"positive": 0.0, "negative": 0.0, "zero": 0.0, "undefined": 1.0}
    return {
        "positive": float(np.sum(finite & (field > 0.0)) / denom),
        "negative": float(np.sum(finite & (field < 0.0)) / denom),
        "zero": float(np.sum(finite & (field == 0.0)) / denom),
        "undefined": float(1.0 - denom / max(int(np.sum(within)), 1)),
    }


def masked_abs_percentiles(field: np.ndarray, mask: np.ndarray, percentiles: tuple[float, ...]) -> dict[str, float]:
    vals = np.abs(field[mask])
    if vals.size == 0:
        return {f"p{p:g}": 0.0 for p in percentiles}
    return {f"p{p:g}": float(np.percentile(vals, p)) for p in percentiles}


def transition_width_proxy(
    r_tilde: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    ell: float,
) -> np.ndarray:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    dr_dx, dr_dz = np.gradient(r_tilde, dx, dz, edge_order=2)
    grad_norm = np.sqrt(dr_dx * dr_dx + dr_dz * dr_dz)
    return 1.0 / np.maximum((ell * ell) * grad_norm, 1.0e-300)


def summarize_branch_time(
    branch: str,
    ell: float,
    t: float,
    fields: dict[str, np.ndarray],
    x: np.ndarray,
    z: np.ndarray,
    frr_floor: float,
) -> dict[str, object]:
    rho = np.real_if_close(fields["rho"])
    r_tilde = np.real_if_close(fields["R_tilde"])
    residual = np.real_if_close(fields["residual_norm"])
    derivative_norm = np.real_if_close(fields["derivative_norm"])
    algebraic_norm = np.real_if_close(fields["algebraic_norm"])

    support = rho > 1.0e-3 * float(np.max(rho))
    dense = rho > 1.0e-1 * float(np.max(rho))

    derivs = branch_f_derivatives(branch, r_tilde, ell)
    y = derivs["y"]
    fr = derivs["fR"]
    frr = derivs["fRR"]
    ms2 = signed_scalaron_mass_squared(r_tilde, fr, frr, frr_floor=frr_floor)

    strict = (np.abs(y) >= 0.5) & (np.abs(y) <= 2.0)
    broad = (np.abs(y) >= 1.0) & (np.abs(y) <= 100.0)
    hotspot_thr = float(np.percentile(residual[support], 99.0)) if np.any(support) else float(np.max(residual))
    hotspot = support & (residual >= hotspot_thr)

    width_proxy = transition_width_proxy(r_tilde, x, z, ell)

    return {
        "branch": branch,
        "ell": float(ell),
        "time": float(t),
        "counts": {
            "support": int(np.sum(support)),
            "dense": int(np.sum(dense)),
            "strict_support": int(np.sum(strict & support)),
            "broad_support": int(np.sum(broad & support)),
            "hotspot_support": int(np.sum(hotspot)),
        },
        "fractions_support": {
            "R_tilde_positive": finite_fraction(r_tilde > 0.0, support),
            "fR_lt_1e-2": finite_fraction(fr < 1.0e-2, support),
            "fR_lt_1e-4": finite_fraction(fr < 1.0e-4, support),
            "fR_lt_1e-8": finite_fraction(fr < 1.0e-8, support),
            "fRR_negative": finite_fraction(frr < 0.0, support),
            "fRR_positive": finite_fraction(frr > 0.0, support),
            "strict_transition": finite_fraction(strict, support),
            "broad_transition": finite_fraction(broad, support),
            "hotspot_inside_strict": finite_fraction(strict, hotspot),
            "hotspot_inside_broad": finite_fraction(broad, hotspot),
        },
        "scalaron_m2_sign_support": finite_sign_fractions(ms2, support),
        "scalaron_m2_abs_support": masked_abs_percentiles(ms2, support & np.isfinite(ms2), (50.0, 95.0, 99.0)),
        "transition_width_proxy_support": masked_abs_percentiles(width_proxy, strict & support, (50.0, 95.0, 99.0)),
        "R_tilde_support": support_stats(r_tilde, support),
        "fR_support": support_stats(fr, support),
        "fRR_support": support_stats(frr, support),
        "residual_norm_support": support_stats(residual, support),
        "residual_norm_dense": support_stats(residual, dense),
        "derivative_norm_support": support_stats(derivative_norm, support),
        "algebraic_norm_support": support_stats(algebraic_norm, support),
        "hotspot_threshold_support_p99": hotspot_thr,
    }


def aggregate(records: list[dict[str, object]]) -> dict[str, object]:
    by_key: dict[tuple[str, float], list[dict[str, object]]] = {}
    for record in records:
        by_key.setdefault((str(record["branch"]), float(record["ell"])), []).append(record)

    out: dict[str, object] = {}
    for (branch, ell), rows in sorted(by_key.items(), key=lambda item: (item[0][0], item[0][1])):
        key = f"{branch}_ell_{ell:g}"
        out[key] = {
            "branch": branch,
            "ell": ell,
            "max_support_residual_p95_over_times": max(float(r["residual_norm_support"]["p95"]) for r in rows),
            "max_support_residual_abs_max_over_times": max(float(r["residual_norm_support"]["abs_max"]) for r in rows),
            "max_support_derivative_abs_max_over_times": max(float(r["derivative_norm_support"]["abs_max"]) for r in rows),
            "max_fR_lt_1e-4_fraction_support": max(float(r["fractions_support"]["fR_lt_1e-4"]) for r in rows),
            "max_fRR_negative_fraction_support": max(float(r["fractions_support"]["fRR_negative"]) for r in rows),
            "max_scalaron_m2_negative_fraction_support": max(float(r["scalaron_m2_sign_support"]["negative"]) for r in rows),
            "max_strict_transition_fraction_support": max(float(r["fractions_support"]["strict_transition"]) for r in rows),
            "max_hotspot_inside_broad_fraction": max(float(r["fractions_support"]["hotspot_inside_broad"]) for r in rows),
        }
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ells", type=str, default="1,3,10,30,100,300")
    parser.add_argument("--branches", type=str, default="C,D")
    parser.add_argument("--times", type=str, default="0,8,16")
    parser.add_argument("--nx", type=int, default=64)
    parser.add_argument("--nz", type=int, default=64)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--frr-floor", type=float, default=1.0e-14)
    args = parser.parse_args()

    ells = [float(tok) for tok in args.ells.split(",") if tok.strip()]
    branches = [tok.strip() for tok in args.branches.split(",") if tok.strip()]
    times = [float(tok) for tok in args.times.split(",") if tok.strip()]

    params = FlatLocalizedCrossingParams(nx=args.nx, nz=args.nz)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)

    records: list[dict[str, object]] = []
    for ell in ells:
        for branch in branches:
            for t in times:
                fields = make_fields(
                    branch=branch,
                    psi0=psi0,
                    psi0_hat=psi0_hat,
                    omega=omega,
                    x=x,
                    z=z,
                    t=t,
                    probe_dt=args.probe_dt,
                    ell=ell,
                    mp=args.mp,
                    mass=params.m,
                    rho_floor=args.rho_floor,
                )
                records.append(
                    summarize_branch_time(
                        branch=branch,
                        ell=ell,
                        t=t,
                        fields=fields,
                        x=x,
                        z=z,
                        frr_floor=args.frr_floor,
                    )
                )

    summary = {
        "params": {
            "ells": ells,
            "branches": branches,
            "times": times,
            "nx": args.nx,
            "nz": args.nz,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "frr_floor": args.frr_floor,
        },
        "definitions": {
            "support": "rho > 1e-3 * rho_max",
            "dense": "rho > 1e-1 * rho_max",
            "strict_transition": "0.5 <= |ell^2 R_tilde| <= 2.0",
            "broad_transition": "1 <= |ell^2 R_tilde| <= 100",
            "transition_width_proxy": "1 / (ell^2 * |spatial grad R_tilde|), evaluated statistically on strict_transition & support",
            "scalaron_m2": "(f_R - R_tilde f_RR) / (3 f_RR), undefined where |f_RR| <= frr_floor",
        },
        "aggregate": aggregate(records),
        "records": records,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps(summary["aggregate"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
