from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from diagnose_cd_transition_jump_mp2 import signed_to_abs_sum_ratio
from plot_bcd_transition_layers import make_fields
from scan_cd_ell_window import branch_f_derivatives, masked_abs_percentiles
from scan_cd_resolution_left_arc import LEFT_ARC_BBOX
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def scalar_stats(values: list[float] | np.ndarray) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"count": 0, "abs_max": 0.0, "abs_mean": 0.0, "abs_median": 0.0, "p95": 0.0}
    abs_vals = np.abs(vals)
    return {
        "count": int(vals.size),
        "abs_max": float(np.max(abs_vals)),
        "abs_mean": float(np.mean(abs_vals)),
        "abs_median": float(np.median(abs_vals)),
        "p95": float(np.percentile(abs_vals, 95.0)),
    }


def bilinear_sample(field: np.ndarray, x: np.ndarray, z: np.ndarray, xp: np.ndarray, zp: np.ndarray) -> np.ndarray:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    u = (xp - float(x[0])) / dx
    v = (zp - float(z[0])) / dz
    i0 = np.floor(u).astype(int)
    j0 = np.floor(v).astype(int)
    i0 = np.clip(i0, 0, len(x) - 2)
    j0 = np.clip(j0, 0, len(z) - 2)
    tx = np.clip(u - i0, 0.0, 1.0)
    tz = np.clip(v - j0, 0.0, 1.0)
    f00 = field[i0, j0]
    f10 = field[i0 + 1, j0]
    f01 = field[i0, j0 + 1]
    f11 = field[i0 + 1, j0 + 1]
    return (
        (1.0 - tx) * (1.0 - tz) * f00
        + tx * (1.0 - tz) * f10
        + (1.0 - tx) * tz * f01
        + tx * tz * f11
    )


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


def find_interface_points(
    y: np.ndarray,
    rho: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    level: float,
) -> list[tuple[float, float]]:
    xmin, xmax, zmin, zmax = LEFT_ARC_BBOX
    rho_threshold = 1.0e-3 * float(np.max(rho))
    y_der = spatial_derivatives(y, x, z)
    points: list[tuple[float, float]] = []

    for i in range(len(x) - 1):
        for j in range(len(z) - 1):
            cell_y = y[i : i + 2, j : j + 2]
            if not (np.nanmin(cell_y) <= level <= np.nanmax(cell_y)):
                continue
            corners = [
                (float(x[i]), float(z[j]), float(y[i, j])),
                (float(x[i + 1]), float(z[j]), float(y[i + 1, j])),
                (float(x[i + 1]), float(z[j + 1]), float(y[i + 1, j + 1])),
                (float(x[i]), float(z[j + 1]), float(y[i, j + 1])),
            ]
            edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
            edge_points: list[tuple[float, float]] = []
            for a, b in edges:
                x0, z0, y0 = corners[a]
                x1, z1, y1 = corners[b]
                dy = y1 - y0
                if abs(dy) <= 1.0e-300:
                    continue
                tau = (level - y0) / dy
                if 0.0 <= tau <= 1.0:
                    xp = x0 + tau * (x1 - x0)
                    zp = z0 + tau * (z1 - z0)
                    if xmin <= xp <= xmax and zmin <= zp <= zmax:
                        edge_points.append((float(xp), float(zp)))
            if not edge_points:
                continue
            # A level set usually crosses a cell in two points. The midpoint is a
            # more stable representative than either edge intersection.
            xp = float(np.mean([p[0] for p in edge_points]))
            zp = float(np.mean([p[1] for p in edge_points]))
            rho_p = float(bilinear_sample(rho, x, z, np.array([xp]), np.array([zp]))[0])
            if rho_p <= rho_threshold:
                continue
            gx = float(bilinear_sample(y_der["fx"], x, z, np.array([xp]), np.array([zp]))[0])
            gz = float(bilinear_sample(y_der["fz"], x, z, np.array([xp]), np.array([zp]))[0])
            if gx * gx + gz * gz <= 1.0e-300:
                continue
            points.append((xp, zp))
    return points


def integrate_lines(
    y: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    points: list[tuple[float, float]],
    half_width_y: float,
    samples: int,
) -> list[dict[str, float]]:
    y_der = spatial_derivatives(y, x, z)
    xmin, xmax = float(x[0]), float(x[-1])
    zmin, zmax = float(z[0]), float(z[-1])
    rows: list[dict[str, float]] = []

    for xp, zp in points:
        gx = float(bilinear_sample(y_der["fx"], x, z, np.array([xp]), np.array([zp]))[0])
        gz = float(bilinear_sample(y_der["fz"], x, z, np.array([xp]), np.array([zp]))[0])
        grad = float(np.hypot(gx, gz))
        if grad <= 1.0e-300:
            continue
        nx = gx / grad
        nz = gz / grad
        half_width = half_width_y / grad
        s = np.linspace(-half_width, half_width, samples)
        xs = xp + s * nx
        zs = zp + s * nz
        valid = (xs >= xmin) & (xs <= xmax) & (zs >= zmin) & (zs <= zmax)
        if int(np.sum(valid)) < max(5, samples // 2):
            continue
        xs = xs[valid]
        zs = zs[valid]
        s_valid = s[valid]

        y_vals = bilinear_sample(y, x, z, xs, zs)
        tanh_y = np.tanh(y_vals)
        phi_vals = 1.0 - tanh_y * tanh_y
        dphi_dn = np.gradient(phi_vals, s_valid, edge_order=2)
        d2phi_dnn = np.gradient(dphi_dn, s_valid, edge_order=2)

        signed_integral = float(np.trapezoid(d2phi_dnn, s_valid))
        abs_integral = float(np.trapezoid(np.abs(d2phi_dnn), s_valid))
        endpoint_jump = float(dphi_dn[-1] - dphi_dn[0])
        rows.append(
            {
                "x": float(xp),
                "z": float(zp),
                "grad_y": grad,
                "half_width": float(half_width),
                "y_min": float(np.min(y_vals)),
                "y_max": float(np.max(y_vals)),
                "endpoint_jump": endpoint_jump,
                "signed_integral_d2phi": signed_integral,
                "abs_integral_d2phi": abs_integral,
                "jump_integral_mismatch": float(endpoint_jump - signed_integral),
                "abs_signed_over_abs_integral": float(abs(signed_integral) / max(abs_integral, 1.0e-300)),
            }
        )
    return rows


def summarize_rows(rows: list[dict[str, float]]) -> dict[str, object]:
    if not rows:
        return {
            "line_count": 0,
            "endpoint_jump": scalar_stats([]),
            "signed_integral_d2phi": scalar_stats([]),
            "abs_integral_d2phi": scalar_stats([]),
            "jump_integral_mismatch": scalar_stats([]),
            "abs_signed_over_abs_integral": scalar_stats([]),
            "global_signed_to_abs_integral_ratio": 0.0,
        }
    signed = np.array([r["signed_integral_d2phi"] for r in rows], dtype=float)
    abs_int = np.array([r["abs_integral_d2phi"] for r in rows], dtype=float)
    return {
        "line_count": len(rows),
        "endpoint_jump": scalar_stats([r["endpoint_jump"] for r in rows]),
        "signed_integral_d2phi": scalar_stats(signed),
        "abs_integral_d2phi": scalar_stats(abs_int),
        "jump_integral_mismatch": scalar_stats([r["jump_integral_mismatch"] for r in rows]),
        "abs_signed_over_abs_integral": scalar_stats([r["abs_signed_over_abs_integral"] for r in rows]),
        "global_signed_to_abs_integral_ratio": float(abs(np.sum(signed)) / max(np.sum(abs_int), 1.0e-300)),
    }


def run_case(
    ell: float,
    resolution: int,
    t: float,
    half_width_y: float,
    samples: int,
    mp: float,
    probe_dt: float,
    rho_floor: float,
) -> dict[str, object]:
    params = FlatLocalizedCrossingParams(nx=resolution, nz=resolution)
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    fields = make_fields(
        branch="D",
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
    rho = np.real_if_close(fields["rho"])
    r_tilde = np.real_if_close(fields["R_tilde"])
    derivs = branch_f_derivatives("D", r_tilde, ell)
    y = np.real_if_close(derivs["y"])

    level_summaries: dict[str, object] = {}
    all_rows: list[dict[str, float]] = []
    for level in (-1.0, 1.0):
        points = find_interface_points(y, rho, x, z, level=level)
        rows = integrate_lines(y, x, z, points, half_width_y=half_width_y, samples=samples)
        for row in rows:
            row["level"] = level
        all_rows.extend(rows)
        level_summaries[f"level_{level:+.0f}"] = {
            "candidate_points": len(points),
            **summarize_rows(rows),
        }

    return {
        "ell": float(ell),
        "resolution": int(resolution),
        "time": float(t),
        "dx": float(x[1] - x[0]),
        "dz": float(z[1] - z[0]),
        "levels": level_summaries,
        "all_levels": summarize_rows(all_rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ells", type=str, default="30,100")
    parser.add_argument("--resolutions", type=str, default="160,240,320")
    parser.add_argument("--time", type=float, default=16.0)
    parser.add_argument("--half-width-y", type=float, default=2.0)
    parser.add_argument("--samples", type=int, default=81)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--probe-dt", type=float, default=5.0e-4)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    args = parser.parse_args()

    ells = [float(tok) for tok in args.ells.split(",") if tok.strip()]
    resolutions = [int(tok) for tok in args.resolutions.split(",") if tok.strip()]

    cases: list[dict[str, object]] = []
    for ell in ells:
        for resolution in resolutions:
            cases.append(
                run_case(
                    ell=ell,
                    resolution=resolution,
                    t=args.time,
                    half_width_y=args.half_width_y,
                    samples=args.samples,
                    mp=args.mp,
                    probe_dt=args.probe_dt,
                    rho_floor=args.rho_floor,
                )
            )

    summary = {
        "params": {
            "branch": "D",
            "ells": ells,
            "resolutions": resolutions,
            "time": args.time,
            "half_width_y": args.half_width_y,
            "samples": args.samples,
            "mp": args.mp,
            "probe_dt": args.probe_dt,
            "rho_floor": args.rho_floor,
            "bbox": list(LEFT_ARC_BBOX),
        },
        "definitions": {
            "interface": "dynamic level sets y=ell^2 R_tilde = +/-1 on the selected time slice",
            "normal": "n = grad(y) / |grad(y)|",
            "line_parameter": "physical coordinate distance along n",
            "half_width_y": "line half-width chosen so linearized y changes by this amount",
            "endpoint_jump": "d_n f_R at line end minus d_n f_R at line start",
            "signed_integral_d2phi": "integral d_n^2 f_R dn along the short line",
            "abs_integral_d2phi": "integral |d_n^2 f_R| dn along the short line",
        },
        "cases": cases,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    compact = {
        f"ell_{case['ell']:g}_n_{case['resolution']}": {
            "line_count": case["all_levels"]["line_count"],
            "endpoint_jump_p95": case["all_levels"]["endpoint_jump"]["p95"],
            "signed_integral_p95": case["all_levels"]["signed_integral_d2phi"]["p95"],
            "abs_integral_p95": case["all_levels"]["abs_integral_d2phi"]["p95"],
            "global_signed_to_abs_ratio": case["all_levels"]["global_signed_to_abs_integral_ratio"],
            "mismatch_p95": case["all_levels"]["jump_integral_mismatch"]["p95"],
        }
        for case in cases
    }
    print(json.dumps(compact, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
