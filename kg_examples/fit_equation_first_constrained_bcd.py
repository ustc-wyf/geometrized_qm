from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyze_c_terms_from_a_reference import d1x, d1z
from diagnose_equation_first_transverse_force import (
    covector_norm_euclidean,
    divergence_cov2,
    derivative_tensor,
    transverse_covector,
)
from diagnose_mathcal_r_pure_geometry import CaseData, build_case, symmetric_rows, tensor_norm
from fit_equation_first_tensor_couplings import tensor_atoms
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference
from simulate_d_tridomain_full_dynamics import erode_mask_8


SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
TIME_KEYS = (-1, 0, 1)
ATOM_NAMES = ("uu", "rr", "ur")
ATOM_FAMILY = "matter4"


def stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0}
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "max": float(np.max(vals)),
    }
    if weights is not None:
        w = np.asarray(weights, dtype=float)[mask]
        w = np.maximum(w, 0.0)
        denom = float(np.sum(w))
        if denom > 0.0:
            out["weighted_mean"] = float(np.sum(w * vals) / denom)
    return out


def dilate_mask_8(mask: np.ndarray, steps: int) -> np.ndarray:
    out = np.asarray(mask, dtype=bool).copy()
    for _ in range(max(int(steps), 0)):
        new = out.copy()
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                src_i0 = max(0, -di)
                src_i1 = out.shape[0] - max(0, di)
                src_j0 = max(0, -dj)
                src_j1 = out.shape[1] - max(0, dj)
                dst_i0 = max(0, di)
                dst_i1 = out.shape[0] - max(0, -di)
                dst_j0 = max(0, dj)
                dst_j1 = out.shape[1] - max(0, -dj)
                new[dst_i0:dst_i1, dst_j0:dst_j1] |= out[src_i0:src_i1, src_j0:src_j1]
        out = new
    return out


def relative_row_weight(case: CaseData, idx: np.ndarray, target_rows: np.ndarray, floor_frac: float) -> np.ndarray:
    target_norm = np.sqrt(np.sum(target_rows**2, axis=1))
    finite = target_norm[np.isfinite(target_norm)]
    floor = float(floor_frac) * max(float(np.percentile(finite, 95.0)) if finite.size else 0.0, 1.0e-300)
    rho = np.maximum(case.rho_a.reshape(-1)[idx], 0.0)
    rho_w = np.sqrt(rho / max(float(np.max(case.rho_a)), 1.0e-300))
    return rho_w / np.maximum(target_norm, floor)


def force_row_weight(case: CaseData, mask: np.ndarray, h: float, floor_frac: float, force_weight: float) -> np.ndarray:
    idx = np.where(mask.reshape(-1))[0]
    target_norm = tensor_norm(case.r_need).reshape(-1)[idx]
    finite = target_norm[np.isfinite(target_norm)]
    floor = float(floor_frac) * max(float(np.percentile(finite, 95.0)) if finite.size else 0.0, 1.0e-300)
    rho = np.maximum(case.rho_a.reshape(-1)[idx], 0.0)
    rho_w = np.sqrt(rho / max(float(np.max(case.rho_a)), 1.0e-300))
    # The force has one extra derivative.  Dividing by |target|/h makes the
    # penalty dimensionless and matches the diagnostic |J_perp|/(|C|/h).
    return float(force_weight) * rho_w / np.maximum(target_norm / max(h, 1.0e-300), floor / max(h, 1.0e-300))


def make_active_mask(cases: dict[int, CaseData], region: str, force_region: str, active_dilation: int) -> np.ndarray:
    mask = np.zeros_like(cases[0].rho_a, dtype=bool)
    for key in TIME_KEYS:
        mask |= getattr(cases[key], region)
    mask |= dilate_mask_8(getattr(cases[0], force_region), active_dilation)
    return mask


def make_force_mask(case: CaseData, force_region: str, erosion: int) -> np.ndarray:
    mask = erode_mask_8(getattr(case, force_region), int(erosion))
    if not np.any(mask):
        mask = getattr(case, force_region).copy()
    mask = mask.copy()
    # Avoid one-sided finite differences in this first prototype.
    mask[0, :] = False
    mask[-1, :] = False
    mask[:, 0] = False
    mask[:, -1] = False
    return mask


def build_column_index(active_mask: np.ndarray) -> tuple[np.ndarray, int]:
    col = -np.ones((len(TIME_KEYS), len(ATOM_NAMES)) + active_mask.shape, dtype=int)
    n = 0
    active_points = np.argwhere(active_mask)
    for tpos, _ in enumerate(TIME_KEYS):
        for apos, _ in enumerate(ATOM_NAMES):
            for i, j in active_points:
                col[tpos, apos, int(i), int(j)] = n
                n += 1
    return col, n


def add_coeff(row: np.ndarray, col_index: np.ndarray, tpos: int, apos: int, i: int, j: int, value: float) -> None:
    c = int(col_index[tpos, apos, i, j])
    if c >= 0 and np.isfinite(value):
        row[c] += float(value)


def assemble_system(
    cases: dict[int, CaseData],
    full_gamma_case,
    *,
    region: str,
    force_region: str,
    target_floor_frac: float,
    force_weight: float,
    force_mode: str,
    force_erosion: int,
    active_dilation: int,
    dt: float,
    dx: float,
    dz: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    active_mask = make_active_mask(cases, region, force_region, active_dilation)
    force_mask = make_force_mask(cases[0], force_region, force_erosion)
    col_index, ncols = build_column_index(active_mask)
    rows: list[np.ndarray] = []
    rhs: list[float] = []

    atom_cache = {key: tensor_atoms(cases[key], "matter4") for key in TIME_KEYS}
    for tpos, key in enumerate(TIME_KEYS):
        case = cases[key]
        mask = getattr(case, region)
        idx = np.where(mask.reshape(-1))[0]
        target_rows = symmetric_rows(case.r_need, idx)
        weights = relative_row_weight(case, idx, target_rows, target_floor_frac)
        points = np.argwhere(mask)
        for p, (i_raw, j_raw) in enumerate(points):
            i = int(i_raw)
            j = int(j_raw)
            for comp_pos, (mu, nu) in enumerate(SYMMETRIC_COMPONENTS):
                row = np.zeros(ncols, dtype=float)
                for apos, atom_name in enumerate(ATOM_NAMES):
                    add_coeff(row, col_index, tpos, apos, i, j, weights[p] * atom_cache[key][atom_name][i, j, mu, nu])
                rows.append(row)
                rhs.append(float(weights[p] * target_rows[p, comp_pos]))

    force_idx = np.where(force_mask.reshape(-1))[0]
    force_weights = force_row_weight(
        cases[0],
        force_mask,
        min(abs(dt), abs(dx), abs(dz)),
        target_floor_frac,
        force_weight,
    )
    metric_inv = full_gamma_case.geom_0.metric_inv
    gamma = full_gamma_case.geom_0.gamma
    u_cov = cases[0].u_cov
    u_up = np.einsum("...ab,...b->...a", metric_inv, u_cov, optimize=True)
    u2 = np.einsum("...a,...a->...", u_up, u_cov, optimize=True)
    projector = np.zeros(force_mask.shape + (3, 3), dtype=float)
    for nu in range(3):
        for sig in range(3):
            projector[..., nu, sig] = (1.0 if nu == sig else 0.0) - u_cov[..., nu] * u_up[..., sig] / np.where(
                np.abs(u2) > 1.0e-300, u2, np.nan
            )

    points_force = np.argwhere(force_mask)
    for p, (i_raw, j_raw) in enumerate(points_force):
        i = int(i_raw)
        j = int(j_raw)
        row_j_sigma = np.zeros((3, ncols), dtype=float)
        for sigma in range(3):
            for mu in range(3):
                for alpha in range(3):
                    pref = metric_inv[i, j, mu, alpha]
                    if not np.isfinite(pref):
                        continue
                    # partial_alpha C_{mu sigma}
                    if alpha == 0:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            val_p = atom_cache[1][atom_name][i, j, mu, sigma] / (2.0 * dt)
                            val_m = -atom_cache[-1][atom_name][i, j, mu, sigma] / (2.0 * dt)
                            add_coeff(row_j_sigma[sigma], col_index, 2, apos, i, j, pref * val_p)
                            add_coeff(row_j_sigma[sigma], col_index, 0, apos, i, j, pref * val_m)
                    elif alpha == 1:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            val_p = atom_cache[0][atom_name][i + 1, j, mu, sigma] / (2.0 * dx)
                            val_m = -atom_cache[0][atom_name][i - 1, j, mu, sigma] / (2.0 * dx)
                            add_coeff(row_j_sigma[sigma], col_index, 1, apos, i + 1, j, pref * val_p)
                            add_coeff(row_j_sigma[sigma], col_index, 1, apos, i - 1, j, pref * val_m)
                    else:
                        for apos, atom_name in enumerate(ATOM_NAMES):
                            val_p = atom_cache[0][atom_name][i, j + 1, mu, sigma] / (2.0 * dz)
                            val_m = -atom_cache[0][atom_name][i, j - 1, mu, sigma] / (2.0 * dz)
                            add_coeff(row_j_sigma[sigma], col_index, 1, apos, i, j + 1, pref * val_p)
                            add_coeff(row_j_sigma[sigma], col_index, 1, apos, i, j - 1, pref * val_m)

                    # -Gamma^lambda_{alpha mu} C_{lambda sigma}
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, mu]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom_name in enumerate(ATOM_NAMES):
                                add_coeff(
                                    row_j_sigma[sigma],
                                    col_index,
                                    1,
                                    apos,
                                    i,
                                    j,
                                    conn * atom_cache[0][atom_name][i, j, lam, sigma],
                                )
                    # -Gamma^lambda_{alpha sigma} C_{mu lambda}
                    for lam in range(3):
                        conn = -pref * gamma[i, j, lam, alpha, sigma]
                        if np.isfinite(conn) and conn != 0.0:
                            for apos, atom_name in enumerate(ATOM_NAMES):
                                add_coeff(
                                    row_j_sigma[sigma],
                                    col_index,
                                    1,
                                    apos,
                                    i,
                                    j,
                                    conn * atom_cache[0][atom_name][i, j, mu, lam],
                                )

        if force_mode == "transverse":
            for nu in range(3):
                row = np.zeros(ncols, dtype=float)
                for sigma in range(3):
                    row += projector[i, j, nu, sigma] * row_j_sigma[sigma]
                rows.append(force_weights[p] * row)
                rhs.append(0.0)
        elif force_mode == "full":
            for sigma in range(3):
                rows.append(force_weights[p] * row_j_sigma[sigma])
                rhs.append(0.0)
        else:
            raise ValueError(f"unknown force_mode: {force_mode}")

    a = np.vstack(rows)
    b = np.asarray(rhs, dtype=float)
    meta = {
        "active_count": int(np.count_nonzero(active_mask)),
        "force_count": int(np.count_nonzero(force_mask)),
        "force_mode": force_mode,
        "ncols": int(ncols),
        "nrows": int(a.shape[0]),
        "active_mask": active_mask,
        "force_mask": force_mask,
        "column_index": col_index,
    }
    return a, b, meta


def solve_scaled_lstsq(a: np.ndarray, b: np.ndarray, ridge: float) -> tuple[np.ndarray, dict[str, float]]:
    col_scale = np.linalg.norm(a, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_scaled = a / col_scale[None, :]
    if ridge > 0.0:
        lam = float(ridge)
        a_solve = np.vstack([a_scaled, lam * np.eye(a_scaled.shape[1])])
        b_solve = np.concatenate([b, np.zeros(a_scaled.shape[1], dtype=float)])
    else:
        a_solve = a_scaled
        b_solve = b
    coeff_scaled, residuals, rank, singular = np.linalg.lstsq(a_solve, b_solve, rcond=1.0e-12)
    coeff = coeff_scaled / col_scale
    residual = a @ coeff - b
    return coeff, {
        "rank": int(rank),
        "singular_min": float(np.min(singular)) if singular.size else 0.0,
        "singular_max": float(np.max(singular)) if singular.size else 0.0,
        "condition_scaled": float(np.max(singular) / max(np.min(singular), 1.0e-300)) if singular.size else 0.0,
        "weighted_system_relative_residual": float(np.linalg.norm(residual) / max(np.linalg.norm(b), 1.0e-300)),
        "weighted_system_residual_norm": float(np.linalg.norm(residual)),
        "weighted_system_rhs_norm": float(np.linalg.norm(b)),
        "ridge": float(ridge),
        "lstsq_residual_sum": float(residuals[0]) if residuals.size else 0.0,
    }


def solve_hard_force_constrained(a: np.ndarray, b: np.ndarray, meta: dict[str, object], ridge: float) -> tuple[np.ndarray, dict[str, object]]:
    n_force = 3 * int(meta["force_count"])
    if n_force <= 0:
        raise ValueError("hard force constraint requires nonzero force rows")
    a_alg = a[:-n_force]
    b_alg = b[:-n_force]
    f_rows = a[-n_force:]

    col_scale = np.linalg.norm(a_alg, axis=0) + np.linalg.norm(f_rows, axis=0)
    col_scale = np.where(col_scale > 0.0, col_scale, 1.0)
    a_scaled = a_alg / col_scale[None, :]
    f_scaled = f_rows / col_scale[None, :]
    eps = float(ridge)
    kkt = np.block(
        [
            [a_scaled.T @ a_scaled + eps * np.eye(a_scaled.shape[1]), f_scaled.T],
            [f_scaled, np.zeros((f_scaled.shape[0], f_scaled.shape[0]))],
        ]
    )
    rhs = np.concatenate([a_scaled.T @ b_alg, np.zeros(f_scaled.shape[0], dtype=float)])
    sol = np.linalg.lstsq(kkt, rhs, rcond=1.0e-12)[0]
    coeff = sol[: a_scaled.shape[1]] / col_scale
    alg_residual = a_alg @ coeff - b_alg
    force_residual = f_rows @ coeff
    return coeff, {
        "hard_constraint": True,
        "force_rows": int(n_force),
        "kkt_shape": [int(kkt.shape[0]), int(kkt.shape[1])],
        "ridge": eps,
        "force_residual_norm": float(np.linalg.norm(force_residual)),
        "force_residual_max_abs": float(np.max(np.abs(force_residual))) if force_residual.size else 0.0,
        "algebraic_weighted_residual_norm": float(np.linalg.norm(alg_residual)),
        "force_over_algebraic_residual_norm": float(
            np.linalg.norm(force_residual) / max(np.linalg.norm(alg_residual), 1.0e-300)
        ),
    }


def coeff_to_fields(coeff: np.ndarray, meta: dict[str, object], shape: tuple[int, int]) -> dict[int, np.ndarray]:
    col_index = np.asarray(meta["column_index"], dtype=int)
    fields: dict[int, np.ndarray] = {}
    for tpos, key in enumerate(TIME_KEYS):
        arr = np.zeros(shape + (len(ATOM_NAMES),), dtype=float)
        for apos, _ in enumerate(ATOM_NAMES):
            ids = col_index[tpos, apos]
            mask = ids >= 0
            arr[..., apos][mask] = coeff[ids[mask]]
        fields[key] = arr
    return fields


def coeff_fields_to_cov(cases: dict[int, CaseData], coeff_fields: dict[int, np.ndarray]) -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    for key in TIME_KEYS:
        atoms = tensor_atoms(cases[key], ATOM_FAMILY)
        cov = np.zeros_like(cases[key].r_need)
        for apos, atom_name in enumerate(ATOM_NAMES):
            cov += coeff_fields[key][..., apos, None, None] * atoms[atom_name]
        out[key] = cov
    return out


def evaluate_solution(
    cases: dict[int, CaseData],
    full_gamma_case,
    coeff: np.ndarray,
    meta: dict[str, object],
    *,
    region: str,
    force_region: str,
    force_erosion: int,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, object]:
    coeff_fields = coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
    cov = coeff_fields_to_cov(cases, coeff_fields)
    by_time: dict[str, object] = {}
    for key in TIME_KEYS:
        case = cases[key]
        rel = tensor_norm(cov[key] - case.r_need) / np.maximum(tensor_norm(case.r_need), 1.0e-300)
        mask = getattr(case, region)
        weights = np.sqrt(np.maximum(case.rho_a[mask], 0.0) / max(float(np.max(case.rho_a)), 1.0e-300))
        by_time[str(key)] = {
            "tau": float(case.tau_old),
            "algebraic_relative_residual": stats(rel[mask], weights),
            "model_norm": stats(tensor_norm(cov[key])[mask], weights),
            "target_norm": stats(tensor_norm(case.r_need)[mask], weights),
        }

    dcov = derivative_tensor(cov[-1], cov[0], cov[1], dt, dx, dz)
    j_cov = divergence_cov2(cov[0], dcov, full_gamma_case.geom_0.metric_inv, full_gamma_case.geom_0.gamma)
    j_perp = transverse_covector(j_cov, cases[0].u_cov, full_gamma_case.geom_0.metric_inv)
    total = covector_norm_euclidean(j_cov)
    perp = covector_norm_euclidean(j_perp)
    c_norm = tensor_norm(cov[0])
    h = min(abs(dt), abs(dx), abs(dz))
    force_mask = make_force_mask(cases[0], force_region, force_erosion)
    weights_force = np.sqrt(np.maximum(cases[0].rho_a[force_mask], 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    return {
        "by_time": by_time,
        "central_force": {
            "force_region": force_region,
            "force_count": int(np.count_nonzero(force_mask)),
            "divergence_norm": stats(total[force_mask], weights_force),
            "transverse_force_norm": stats(perp[force_mask], weights_force),
            "transverse_fraction": stats(perp[force_mask] / np.maximum(total[force_mask], 1.0e-300), weights_force),
            "transverse_over_derivative_scale": stats(
                perp[force_mask] / np.maximum(c_norm[force_mask] / max(h, 1.0e-300), 1.0e-300),
                weights_force,
            ),
        },
        "coefficient_stats": {
            name: stats(coeff_fields[0][..., apos][np.asarray(meta["active_mask"], dtype=bool)])
            for apos, name in enumerate(ATOM_NAMES)
        },
    }


def render_tradeoff(out_path: Path, records: list[dict[str, object]]) -> None:
    lambdas = np.asarray([float(r["force_weight"]) for r in records], dtype=float)
    x = np.arange(len(lambdas))
    alg = np.asarray(
        [float(r["evaluation"]["by_time"]["0"]["algebraic_relative_residual"].get("weighted_mean", 0.0)) for r in records]
    )
    force = np.asarray(
        [float(r["evaluation"]["central_force"]["transverse_over_derivative_scale"].get("weighted_mean", 0.0)) for r in records]
    )
    frac = np.asarray([float(r["evaluation"]["central_force"]["transverse_fraction"].get("p50", 0.0)) for r in records])
    fig, ax = plt.subplots(figsize=(9.5, 5.0), constrained_layout=True)
    ax.plot(x, alg, "o-", label="algebraic residual weighted mean")
    ax.plot(x, force, "o-", label="transverse/(|C|/h) weighted mean")
    ax.plot(x, frac, "o-", label="transverse fraction p50")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:g}" for v in lambdas], rotation=25, ha="right")
    ax.set_xlabel("force penalty weight lambda")
    ax.set_title("Constrained B,C,D equation-first fit: residual/force tradeoff")
    ax.legend(fontsize=8)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_maps(
    out_path: Path,
    cases: dict[int, CaseData],
    full_gamma_case,
    coeff: np.ndarray,
    meta: dict[str, object],
    region: str,
    force_region: str,
    force_erosion: int,
    dt: float,
    dx: float,
    dz: float,
) -> None:
    coeff_fields = coeff_to_fields(coeff, meta, cases[0].rho_a.shape)
    cov = coeff_fields_to_cov(cases, coeff_fields)
    dcov = derivative_tensor(cov[-1], cov[0], cov[1], dt, dx, dz)
    j_cov = divergence_cov2(cov[0], dcov, full_gamma_case.geom_0.metric_inv, full_gamma_case.geom_0.gamma)
    j_perp = transverse_covector(j_cov, cases[0].u_cov, full_gamma_case.geom_0.metric_inv)
    perp = covector_norm_euclidean(j_perp)
    c_norm = tensor_norm(cov[0])
    h = min(abs(dt), abs(dx), abs(dz))
    force_scaled = perp / np.maximum(c_norm / max(h, 1.0e-300), 1.0e-300)
    rel = tensor_norm(cov[0] - cases[0].r_need) / np.maximum(tensor_norm(cases[0].r_need), 1.0e-300)
    x_um = cases[0].x * HBAR_C_EV_M * 1.0e6
    z_um = cases[0].z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    active = np.asarray(meta["active_mask"], dtype=bool)
    force_mask = make_force_mask(cases[0], force_region, force_erosion)
    fields = [
        ("rho_A", cases[0].rho_a, "viridis", False),
        ("algebraic relative residual", rel, "inferno", True),
        ("transverse force / (|C|/h)", force_scaled, "magma", True),
        ("active coefficient region", active.astype(float), "gray", False),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(18.0, 4.5), constrained_layout=True)
    for ax, (title, data, cmap, clip_masked) in zip(axes, fields):
        values = np.asarray(data, dtype=float)
        if clip_masked:
            mask = force_mask if title.startswith("transverse") else getattr(cases[0], region)
            finite = values[mask & np.isfinite(values)]
            vmax = max(float(np.percentile(finite, 95.0)) if finite.size else 1.0, 1.0e-14)
            plot_data = np.clip(values, 0.0, vmax)
        else:
            finite = values[np.isfinite(values)]
            vmax = max(float(np.percentile(np.abs(finite), 99.0)) if finite.size else 1.0, 1.0e-300)
            plot_data = np.clip(values, 0.0, vmax)
        im = ax.pcolormesh(xg, zg, plot_data, shading="auto", cmap=cmap, vmin=0.0, vmax=vmax)
        ax.contour(xg, zg, getattr(cases[0], region).astype(float), levels=[0.5], colors="white", linewidths=0.7)
        ax.contour(xg, zg, force_mask.astype(float), levels=[0.5], colors="cyan", linewidths=0.6)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    global ATOM_NAMES
    ATOM_NAMES = tuple(x.strip() for x in str(args.atoms).split(",") if x.strip())
    allowed_atoms = {"g", "uu", "rr", "ur"}
    unknown_atoms = [name for name in ATOM_NAMES if name not in allowed_atoms]
    if unknown_atoms:
        raise ValueError(f"unknown atom names {unknown_atoms}; allowed atoms are {sorted(allowed_atoms)}")
    if not ATOM_NAMES:
        raise ValueError("at least one atom is required")

    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    probe = float(args.force_probe_dt_old)
    cases = {
        -1: build_case(args, ref, float(args.tau) - probe),
        0: build_case(args, ref, float(args.tau)),
        1: build_case(args, ref, float(args.tau) + probe),
    }
    full_gamma_case = build_full_case(args, ref, float(args.tau))
    scale = ref["scale"]
    dt = probe * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    force_weights = [1.0] if args.hard_constraint else [float(x) for x in args.force_weights.split(",") if x.strip()]
    records: list[dict[str, object]] = []
    best_record: dict[str, object] | None = None
    best_coeff: np.ndarray | None = None
    best_meta: dict[str, object] | None = None

    for force_weight in force_weights:
        a, b, meta = assemble_system(
            cases,
            full_gamma_case,
            region=args.fit_region,
            force_region=args.force_region,
            target_floor_frac=float(args.target_floor_frac),
            force_weight=force_weight,
            force_mode=args.force_mode,
            force_erosion=int(args.force_mask_erosion),
            active_dilation=int(args.active_dilation),
            dt=dt,
            dx=dx,
            dz=dz,
        )
        if args.hard_constraint:
            coeff, solve_info = solve_hard_force_constrained(a, b, meta, float(args.ridge))
        else:
            coeff, solve_info = solve_scaled_lstsq(a, b, float(args.ridge))
        evaluation = evaluate_solution(
            cases,
            full_gamma_case,
            coeff,
            meta,
            region=args.fit_region,
            force_region=args.force_region,
            force_erosion=int(args.force_mask_erosion),
            dt=dt,
            dx=dx,
            dz=dz,
        )
        record = {
            "force_weight": float(force_weight),
            "atoms": list(ATOM_NAMES),
            "system": {k: v for k, v in meta.items() if k not in {"active_mask", "force_mask", "column_index"}},
            "solve": solve_info,
            "evaluation": evaluation,
        }
        records.append(record)
        alg = float(evaluation["by_time"]["0"]["algebraic_relative_residual"].get("weighted_mean", 1.0e300))
        force = float(evaluation["central_force"]["transverse_over_derivative_scale"].get("weighted_mean", 1.0e300))
        score = alg + force
        if best_record is None or score < float(best_record["selection_score"]):
            record["selection_score"] = score
            best_record = record
            best_coeff = coeff
            best_meta = meta
        else:
            record["selection_score"] = score

    if best_record is None or best_coeff is None or best_meta is None:
        raise RuntimeError("no constrained fits were produced")

    tradeoff_path = args.output / "constrained_bcd_tradeoff.png"
    maps_path = args.output / "constrained_bcd_best_maps.png"
    coeff_path = args.output / "constrained_bcd_best_coefficients.npz"
    render_tradeoff(tradeoff_path, records)
    render_maps(
        maps_path,
        cases,
        full_gamma_case,
        best_coeff,
        best_meta,
        args.fit_region,
        args.force_region,
        int(args.force_mask_erosion),
        dt,
        dx,
        dz,
    )
    coeff_fields = coeff_to_fields(best_coeff, best_meta, cases[0].rho_a.shape)
    coeff_arrays = {}
    suffix_by_time = {-1: "m", 0: "0", 1: "p"}
    for key in TIME_KEYS:
        for apos, atom_name in enumerate(ATOM_NAMES):
            coeff_arrays[f"{atom_name}_{suffix_by_time[key]}"] = coeff_fields[key][..., apos]
    coeff_arrays["active_mask"] = np.asarray(best_meta["active_mask"], dtype=bool)
    coeff_arrays["force_mask"] = np.asarray(best_meta["force_mask"], dtype=bool)
    np.savez_compressed(coeff_path, **coeff_arrays)
    report = {
        "parameters": {
            "tau": float(args.tau),
            "force_probe_dt_old": float(args.force_probe_dt_old),
            "full_resolution": int(args.full_resolution),
            "cropped_shape": [int(cases[0].rho_a.shape[0]), int(cases[0].rho_a.shape[1])],
            "window_um": float(args.window_um),
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "fit_region": args.fit_region,
            "force_region": args.force_region,
            "force_mode": args.force_mode,
            "force_mask_erosion": int(args.force_mask_erosion),
            "atoms": list(ATOM_NAMES),
            "hard_constraint": bool(args.hard_constraint),
            "target_floor_frac": float(args.target_floor_frac),
            "ridge": float(args.ridge),
            "force_weights": force_weights,
            "mp_ev": float(args.mp),
        },
        "definition": {
            "equation": "Gtilde_mn = Ttilde_mn/Mp^2 + C_mn, with C_mn expanded in the listed atoms.",
            "algebraic_residual": "|C_mn - (Gtilde_mn - Ttilde_mn/Mp^2)| / |Gtilde_mn - Ttilde_mn/Mp^2|.",
            "geodesic_constraint": "transverse mode enforces h^nu_alpha nabla^mu C_mu^alpha = 0; full mode enforces nabla^mu C_mu_nu = 0.",
            "force_scale": "transverse_over_derivative_scale = |J_perp|/(|C|/h), where h is the smallest grid/time finite-difference step.",
            "status": "Prototype constrained fit on a fixed A-reference generated tilde geometry. It tests closure feasibility; it is not yet a full D evolution solver.",
        },
        "records": records,
        "best_force_weight": float(best_record["force_weight"]),
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "tradeoff_png": str(tradeoff_path.resolve()),
            "best_maps_png": str(maps_path.resolve()),
            "coefficients_npz": str(coeff_path.resolve()),
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tau", type=float, default=0.0)
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-mode", choices=["transverse", "full"], default="transverse")
    parser.add_argument("--atoms", type=str, default="uu,rr,ur")
    parser.add_argument("--hard-constraint", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--force-weights", type=str, default="0,0.03,0.1,0.3,1,3,10,30")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--ridge", type=float, default=0.0)
    parser.add_argument("--full-resolution", type=int, default=96)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-mask-erosion", type=int, default=1)
    parser.add_argument("--active-dilation", type=int, default=1)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--mp", type=float, default=None)
    parser.add_argument("--ell", type=float, default=None)
    parser.add_argument("--ell-over-planck", type=float, default=1.0e60)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--x-floor", type=float, default=1.0e-10)
    parser.add_argument("--pinv-rcond", type=float, default=1.0e-10)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
