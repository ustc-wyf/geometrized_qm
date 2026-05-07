from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_d_jt_degeneracy import (
    parse_float_list,
    real_array,
    scalar_stats,
)
from mixed_tilde_initial_data import localized_direct_tilde_coordinate_snapshot
from physical_units import HBAR_C_EV_M, ev_inv_to_fs
from simulate_d_tridomain_full_dynamics import (
    build_crossing_params,
    erode_mask_8,
    kg_positive_frequency_norm,
    tilde_measure_from_bohm,
)
from simulate_flat_localized_crossing_packets import initial_wavefunction, make_grid, spectral_omega


@dataclass(frozen=True)
class CandidateCovector:
    a: float
    b: float

    @property
    def label(self) -> str:
        return f"dt{self.a:+.3g}dx{self.b:+.3g}dz"

    @property
    def flat_norm(self) -> float:
        return 1.0 - self.a * self.a - self.b * self.b


def parse_candidate_list(text: str) -> list[CandidateCovector]:
    candidates: list[CandidateCovector] = []
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        parts = [float(x.strip()) for x in item.split(",") if x.strip()]
        if len(parts) != 2:
            raise ValueError(f"Candidate '{item}' must be 'a,b'.")
        candidates.append(CandidateCovector(parts[0], parts[1]))
    return unique_candidates(candidates)


def unique_candidates(candidates: list[CandidateCovector]) -> list[CandidateCovector]:
    seen: set[tuple[float, float]] = set()
    unique: list[CandidateCovector] = []
    for candidate in candidates:
        key = (round(candidate.a, 12), round(candidate.b, 12))
        if key in seen:
            continue
        seen.add(key)
        unique.append(candidate)
    return unique


def build_candidate_grid(max_tilt: float, step: float, manual: str) -> list[CandidateCovector]:
    if step <= 0.0:
        raise ValueError("--tilt-step must be positive.")
    values = np.arange(-max_tilt, max_tilt + 0.5 * step, step)
    candidates = [CandidateCovector(0.0, 0.0)]
    candidates.extend(CandidateCovector(float(a), float(b)) for a in values for b in values)
    if manual:
        candidates.extend(parse_candidate_list(manual))
    return unique_candidates(candidates)


def sign_oriented_boundary_count(labels: np.ndarray, signs: np.ndarray, mask: np.ndarray) -> int:
    oriented = 2 * np.asarray(labels, dtype=int) + (np.asarray(signs, dtype=float) >= 0.0).astype(int)
    m = np.asarray(mask, dtype=bool)
    x_edges = m[:-1, :] & m[1:, :] & (oriented[:-1, :] != oriented[1:, :])
    z_edges = m[:, :-1] & m[:, 1:] & (oriented[:, :-1] != oriented[:, 1:])
    return int(np.count_nonzero(x_edges) + np.count_nonzero(z_edges))


def measure_fraction(mask: np.ndarray, measure: np.ndarray, base: np.ndarray) -> float:
    denom = max(float(np.sum(measure[base])), 1.0e-300)
    return float(np.sum(measure[base & mask]) / denom)


def summarize_cover(
    *,
    name: str,
    abs_values: np.ndarray,
    measure: np.ndarray,
    mask: np.ndarray,
    thresholds: list[float],
) -> dict[str, object]:
    rows = {}
    for threshold in thresholds:
        bad = abs_values < threshold
        rows[f"abs_jtau_lt_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask & bad)),
            "count_fraction": float(np.count_nonzero(mask & bad) / max(int(np.count_nonzero(mask)), 1)),
            "measure_fraction": measure_fraction(bad, measure, mask),
        }
    return {
        "name": name,
        "count": int(np.count_nonzero(mask)),
        "measure_sum": float(np.sum(measure[mask])),
        "abs_jtau": scalar_stats(abs_values[mask]),
        "thresholds": rows,
    }


def compute_slice_data(
    *,
    psi0: np.ndarray,
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    t: float,
    t_old: float,
    mass: float,
    rho_floor: float,
    support_rho_frac: float,
    support_measure_frac: float,
    trusted_erosion: int,
) -> dict[str, object]:
    snap = localized_direct_tilde_coordinate_snapshot(
        psi0=psi0,
        psi0_hat=psi0_hat,
        omega=omega,
        x=x,
        z=z,
        t=t,
        rho_floor=rho_floor,
    )
    bohm = snap["bohm"]
    rho = real_array(bohm["rho"])
    measure = tilde_measure_from_bohm(bohm, mass)
    metric_inv = real_array(snap["inverse"]["tilde_ginv_txz"])
    u_cov = np.stack(
        [
            real_array(bohm["s_t"]),
            real_array(bohm["s_x"]),
            real_array(bohm["s_z"]),
        ],
        axis=-1,
    )
    j_vec = np.einsum("...mn,...n->...m", metric_inv, u_cov)
    support = (rho > support_rho_frac * float(np.max(rho))) & (
        measure > support_measure_frac * float(np.max(measure))
    )
    trusted = erode_mask_8(support, trusted_erosion)
    if not np.any(trusted):
        trusted = support.copy()
    return {
        "t_old": float(t_old),
        "t": float(t),
        "t_fs": ev_inv_to_fs(t),
        "rho": rho,
        "measure": measure,
        "j_vec": j_vec,
        "metric_inv": metric_inv,
        "support": support,
        "trusted": trusted,
    }


def evaluate_candidates(j_vec: np.ndarray, candidates: list[CandidateCovector]) -> np.ndarray:
    a = np.asarray([candidate.a for candidate in candidates], dtype=float)
    b = np.asarray([candidate.b for candidate in candidates], dtype=float)
    return j_vec[..., 0, None] + j_vec[..., 1, None] * a[None, None, :] + j_vec[..., 2, None] * b[None, None, :]


def evaluate_candidate_norms(metric_inv: np.ndarray, candidates: list[CandidateCovector]) -> np.ndarray:
    covectors = np.asarray([[1.0, candidate.a, candidate.b] for candidate in candidates], dtype=float)
    return np.einsum("...mn,am,an->...a", metric_inv, covectors, covectors)


def flatten_training_samples(
    slices: list[dict[str, object]],
    candidate_values: list[np.ndarray],
    mask_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    arrays: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for item, values in zip(slices, candidate_values):
        mask = np.asarray(item[mask_name], dtype=bool)
        measure = np.asarray(item["measure"], dtype=float)
        if not np.any(mask):
            continue
        arrays.append(np.abs(values[mask, :]))
        weights.append(measure[mask])
    if not arrays:
        return np.zeros((0, 0), dtype=float), np.zeros(0, dtype=float)
    return np.concatenate(arrays, axis=0), np.concatenate(weights, axis=0)


def greedy_select_candidates(
    abs_values: np.ndarray,
    weights: np.ndarray,
    candidates: list[CandidateCovector],
    *,
    max_patches: int,
    cover_threshold: float,
    force_lab: bool,
) -> list[int]:
    if abs_values.size == 0:
        return [0]
    selected: list[int] = []
    current_best = np.zeros(abs_values.shape[0], dtype=float)
    if force_lab:
        selected.append(0)
        current_best = np.maximum(current_best, abs_values[:, 0])

    total_weight = max(float(np.sum(weights)), 1.0e-300)
    for _ in range(max(0, max_patches - len(selected))):
        current_covered = current_best >= cover_threshold
        best_score = -np.inf
        best_idx = None
        for idx in range(len(candidates)):
            if idx in selected:
                continue
            proposed_best = np.maximum(current_best, abs_values[:, idx])
            proposed_covered = proposed_best >= cover_threshold
            new_coverage = float(np.sum(weights[proposed_covered]) / total_weight)
            old_coverage = float(np.sum(weights[current_covered]) / total_weight)
            p01_gain = float(np.percentile(proposed_best, 1.0) - np.percentile(current_best, 1.0))
            p05_gain = float(np.percentile(proposed_best, 5.0) - np.percentile(current_best, 5.0))
            score = (new_coverage - old_coverage) + 1.0e-3 * p05_gain + 1.0e-4 * p01_gain
            if score > best_score:
                best_score = score
                best_idx = idx
        if best_idx is None:
            break
        selected.append(best_idx)
        current_best = np.maximum(current_best, abs_values[:, best_idx])
        if float(np.sum(weights[current_best >= cover_threshold]) / total_weight) >= 0.999999:
            break
    return selected


def label_stats(
    labels: np.ndarray,
    signs: np.ndarray,
    measure: np.ndarray,
    mask: np.ndarray,
    selected_candidates: list[CandidateCovector],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    total_measure = max(float(np.sum(measure[mask])), 1.0e-300)
    for local_idx, candidate in enumerate(selected_candidates):
        selected = mask & (labels == local_idx)
        if not np.any(selected):
            continue
        positive = selected & (signs >= 0.0)
        negative = selected & (signs < 0.0)
        rows.append(
            {
                "local_index": int(local_idx),
                "candidate": {"a": candidate.a, "b": candidate.b, "flat_norm": candidate.flat_norm},
                "count": int(np.count_nonzero(selected)),
                "measure_fraction": float(np.sum(measure[selected]) / total_measure),
                "positive_measure_fraction_within_label": float(
                    np.sum(measure[positive]) / max(float(np.sum(measure[selected])), 1.0e-300)
                ),
                "negative_measure_fraction_within_label": float(
                    np.sum(measure[negative]) / max(float(np.sum(measure[selected])), 1.0e-300)
                ),
            }
        )
    return rows


def analyze_slice_cover(
    item: dict[str, object],
    values_all: np.ndarray,
    norms_all: np.ndarray,
    selected_indices: list[int],
    candidates: list[CandidateCovector],
    thresholds: list[float],
    norm_floor: float,
) -> dict[str, object]:
    measure = np.asarray(item["measure"], dtype=float)
    support = np.asarray(item["support"], dtype=bool)
    trusted = np.asarray(item["trusted"], dtype=bool)

    lab = values_all[..., 0]
    abs_lab = np.abs(lab)
    admissible_all = norms_all > norm_floor
    values_admissible = np.where(admissible_all, values_all, 0.0)
    abs_all = np.abs(values_admissible)
    all_label = np.argmax(abs_all, axis=-1)
    all_best = np.take_along_axis(abs_all, all_label[..., None], axis=-1)[..., 0]
    all_sign = np.take_along_axis(values_admissible, all_label[..., None], axis=-1)[..., 0]
    admissible_count = np.sum(admissible_all, axis=-1)

    selected_values = values_admissible[..., selected_indices]
    abs_selected = np.abs(selected_values)
    selected_label = np.argmax(abs_selected, axis=-1)
    selected_best = np.take_along_axis(abs_selected, selected_label[..., None], axis=-1)[..., 0]
    selected_sign = np.take_along_axis(selected_values, selected_label[..., None], axis=-1)[..., 0]

    selected_candidates = [candidates[idx] for idx in selected_indices]
    return {
        "t_old": item["t_old"],
        "t": item["t"],
        "t_fs": item["t_fs"],
        "support": {
            "lab": summarize_cover(
                name="lab_time_dt",
                abs_values=abs_lab,
                measure=measure,
                mask=support,
                thresholds=thresholds,
            ),
            "all_grid_best": summarize_cover(
                name="admissible_grid_best",
                abs_values=all_best,
                measure=measure,
                mask=support,
                thresholds=thresholds,
            ),
            "selected_best": summarize_cover(
                name="selected_best",
                abs_values=selected_best,
                measure=measure,
                mask=support,
                thresholds=thresholds,
            ),
            "selected_oriented_boundary_edges": sign_oriented_boundary_count(selected_label, selected_sign, support),
            "no_admissible_covector_measure_fraction": measure_fraction(admissible_count == 0, measure, support),
            "admissible_covector_count": scalar_stats(admissible_count[support]),
            "selected_label_stats": label_stats(
                selected_label, selected_sign, measure, support, selected_candidates
            ),
        },
        "trusted": {
            "lab": summarize_cover(
                name="lab_time_dt",
                abs_values=abs_lab,
                measure=measure,
                mask=trusted,
                thresholds=thresholds,
            ),
            "all_grid_best": summarize_cover(
                name="admissible_grid_best",
                abs_values=all_best,
                measure=measure,
                mask=trusted,
                thresholds=thresholds,
            ),
            "selected_best": summarize_cover(
                name="selected_best",
                abs_values=selected_best,
                measure=measure,
                mask=trusted,
                thresholds=thresholds,
            ),
            "selected_oriented_boundary_edges": sign_oriented_boundary_count(selected_label, selected_sign, trusted),
            "no_admissible_covector_measure_fraction": measure_fraction(admissible_count == 0, measure, trusted),
            "admissible_covector_count": scalar_stats(admissible_count[trusted]),
            "selected_label_stats": label_stats(
                selected_label, selected_sign, measure, trusted, selected_candidates
            ),
        },
        "fields": {
            "abs_lab": abs_lab,
            "all_best": all_best,
            "selected_best": selected_best,
            "selected_label": selected_label,
            "selected_sign": selected_sign,
        },
    }


def choose_plot_indices(times_old: list[float], plot_times_old: list[float]) -> list[int]:
    if not times_old:
        return []
    indices: list[int] = []
    arr = np.asarray(times_old, dtype=float)
    for requested in plot_times_old:
        idx = int(np.argmin(np.abs(arr - requested)))
        if idx not in indices:
            indices.append(idx)
    return indices


def render_cover_plot(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    slices: list[dict[str, object]],
    analyses: list[dict[str, object]],
    plot_indices: list[int],
    selected_candidates: list[CandidateCovector],
) -> None:
    plot_factor = HBAR_C_EV_M * 1.0e6
    x_plot = x * plot_factor
    z_plot = z * plot_factor
    xg, zg = np.meshgrid(x_plot, z_plot, indexing="ij")
    n = len(plot_indices)
    fig, axes = plt.subplots(n, 4, figsize=(20.0, 4.6 * n), constrained_layout=True)
    if n == 1:
        axes = axes[None, :]

    for row, idx in enumerate(plot_indices):
        item = slices[idx]
        analysis = analyses[idx]
        fields = analysis["fields"]
        rho = np.asarray(item["rho"], dtype=float)
        support = np.asarray(item["support"], dtype=bool)
        trusted = np.asarray(item["trusted"], dtype=bool)
        title = f"t_old={float(item['t_old']):g}, t={float(item['t_fs']):.3g} fs"

        vmax_rho = max(float(np.percentile(rho[support], 99.5)) if np.any(support) else float(np.max(rho)), 1.0e-16)
        im0 = axes[row, 0].pcolormesh(xg, zg, rho, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax_rho)
        axes[row, 0].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.7)
        axes[row, 0].contour(xg, zg, trusted.astype(float), levels=[0.5], colors="cyan", linewidths=0.7)
        axes[row, 0].set_title(f"{title}: rho_A")
        fig.colorbar(im0, ax=axes[row, 0])

        lab_log = np.log10(np.maximum(fields["abs_lab"], 1.0e-16))
        im1 = axes[row, 1].pcolormesh(xg, zg, lab_log, shading="auto", cmap="magma")
        axes[row, 1].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        axes[row, 1].set_title("lab time: log10 |j^t|")
        fig.colorbar(im1, ax=axes[row, 1])

        selected_log = np.log10(np.maximum(fields["selected_best"], 1.0e-16))
        im2 = axes[row, 2].pcolormesh(xg, zg, selected_log, shading="auto", cmap="magma")
        axes[row, 2].contour(xg, zg, support.astype(float), levels=[0.5], colors="white", linewidths=0.6)
        axes[row, 2].set_title("local-time cover: log10 max |j^tau_A|")
        fig.colorbar(im2, ax=axes[row, 2])

        oriented = 2 * fields["selected_label"].astype(int) + (fields["selected_sign"] >= 0.0).astype(int)
        masked_oriented = np.where(support, oriented, np.nan)
        im3 = axes[row, 3].pcolormesh(xg, zg, masked_oriented, shading="auto", cmap="tab20")
        axes[row, 3].contour(xg, zg, support.astype(float), levels=[0.5], colors="black", linewidths=0.6)
        axes[row, 3].set_title("selected local-time patch/sign label")
        fig.colorbar(im3, ax=axes[row, 3])

        for ax in axes[row, :]:
            ax.set_xlabel("x [um]")
            ax.set_ylabel("z [um]")
            ax.set_aspect("equal")

    subtitle = ", ".join(f"{i}:dt{c.a:+.2g}dx{c.b:+.2g}dz" for i, c in enumerate(selected_candidates[:12]))
    fig.suptitle(f"Selected covectors: {subtitle}", fontsize=11)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    params, physical_scale = build_crossing_params(args)
    if not args.physical_optical or physical_scale is None:
        raise ValueError("This diagnostic expects --physical-optical for old-time labels.")

    x, z, x_grid, z_grid = make_grid(params)
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    psi0 = initial_wavefunction(x_grid, z_grid, params)
    psi0_hat = np.fft.fft2(psi0)
    omega = spectral_omega(x, z, params.m)
    psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_before = kg_positive_frequency_norm(psi0, psi0_t, dx, dz)
    if args.normalize_kg:
        psi0 = psi0 / np.sqrt(max(kg_norm_before, 1.0e-300))
        psi0_hat = np.fft.fft2(psi0)
        psi0_t = np.fft.ifft2((-1j * omega) * psi0_hat)
    kg_norm_after = kg_positive_frequency_norm(psi0, psi0_t, dx, dz)

    old_scale = physical_scale["old_dimensionless_scale_ev_inv"]
    times_old = parse_float_list(args.times_old)
    thresholds = parse_float_list(args.jtau_thresholds)
    candidates = build_candidate_grid(args.tilt_max, args.tilt_step, args.manual_candidates)

    slices = [
        compute_slice_data(
            psi0=psi0,
            psi0_hat=psi0_hat,
            omega=omega,
            x=x,
            z=z,
            t=t_old * old_scale,
            t_old=t_old,
            mass=params.m,
            rho_floor=args.rho_floor,
            support_rho_frac=args.support_rho_frac,
            support_measure_frac=args.support_measure_frac,
            trusted_erosion=args.trusted_erosion,
        )
        for t_old in times_old
    ]
    values_all = [evaluate_candidates(np.asarray(item["j_vec"], dtype=float), candidates) for item in slices]
    norms_all = [
        evaluate_candidate_norms(np.asarray(item["metric_inv"], dtype=float), candidates) for item in slices
    ]
    if args.require_tilde_timelike:
        training_values = [
            np.where(norms > args.time_covector_norm_floor, values, 0.0)
            for values, norms in zip(values_all, norms_all)
        ]
    else:
        training_values = values_all
        norms_all = [np.ones_like(values, dtype=float) for values in values_all]
    train_abs, train_weights = flatten_training_samples(slices, training_values, args.training_mask)
    selected_indices = greedy_select_candidates(
        train_abs,
        train_weights,
        candidates,
        max_patches=args.max_patches,
        cover_threshold=args.cover_threshold,
        force_lab=args.force_lab,
    )
    selected_candidates = [candidates[idx] for idx in selected_indices]

    analyses = [
        analyze_slice_cover(
            item,
            values,
            norms,
            selected_indices,
            candidates,
            thresholds,
            args.time_covector_norm_floor if args.require_tilde_timelike else -np.inf,
        )
        for item, values, norms in zip(slices, values_all, norms_all)
    ]

    plot_times = parse_float_list(args.plot_times_old)
    plot_indices = choose_plot_indices(times_old, plot_times)
    fig_path = out / "d_local_time_cover_summary.png"
    render_cover_plot(fig_path, x, z, slices, analyses, plot_indices, selected_candidates)

    summary = {
        "params": {
            "resolution": args.resolution,
            "physical_optical": bool(args.physical_optical),
            "physical_scale": physical_scale,
            "normalize_kg": bool(args.normalize_kg),
            "kg_norm_before": kg_norm_before,
            "kg_norm_after": kg_norm_after,
            "support_rho_frac": args.support_rho_frac,
            "support_measure_frac": args.support_measure_frac,
            "trusted_erosion": args.trusted_erosion,
            "times_old": times_old,
            "plot_times_old": plot_times,
            "jtau_thresholds": thresholds,
            "candidate_grid": {
                "tilt_max": args.tilt_max,
                "tilt_step": args.tilt_step,
                "candidate_count": len(candidates),
                "manual_candidates": args.manual_candidates,
            },
            "greedy": {
                "training_mask": args.training_mask,
                "cover_threshold": args.cover_threshold,
                "max_patches": args.max_patches,
                "force_lab": bool(args.force_lab),
                "require_tilde_timelike": bool(args.require_tilde_timelike),
                "time_covector_norm_floor": args.time_covector_norm_floor,
            },
        },
        "definitions": {
            "j_vec": "j^mu/rho_tilde = gtilde^{mu nu} u_nu, with u_nu=partial_nu S.",
            "candidate_time_covector": "d tau_A = dt + a_A dx + b_A dz.",
            "candidate_timelike_test": "A candidate is admissible at a point when gtilde^{mu nu}(d tau_A)_mu(d tau_A)_nu is above time_covector_norm_floor.  This prevents the diagnostic from using arbitrary strongly tilted non-time covectors as fake time coordinates.",
            "jtau": "j^{tau_A}/rho_tilde = (d tau_A)_mu gtilde^{mu nu} u_nu.  A local-time patch is non-characteristic where |jtau| is bounded away from zero.",
            "selected_best": "For each point, choose the selected candidate A that maximizes |jtau_A|.  This is a coordinate-patch diagnostic, not a physical damping or clipping rule.",
            "oriented_patch_label": "The plotted label includes both selected covector and sign of jtau.  Opposite signs mean opposite local time orientation, not negative density.",
            "uncovered_fraction": "threshold rows report the transformed-measure fraction in the mask where |jtau| is still below the threshold.",
        },
        "selected_candidates": [
            {"global_index": int(idx), "a": candidates[idx].a, "b": candidates[idx].b, "flat_norm": candidates[idx].flat_norm}
            for idx in selected_indices
        ],
        "slices": [
            {
                "t_old": analysis["t_old"],
                "t": analysis["t"],
                "t_fs": analysis["t_fs"],
                "support": analysis["support"],
                "trusted": analysis["trusted"],
            }
            for analysis in analyses
        ],
        "files": {
            "figure": str(fig_path.resolve()),
            "summary": str((out / "summary.json").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--physical-optical", action="store_true", default=True)
    parser.add_argument("--wavelength-nm", type=float, default=1550.0)
    parser.add_argument("--mass-over-omega", type=float, default=0.1)
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--normalize-kg", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--times-old", default="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16")
    parser.add_argument("--plot-times-old", default="0,8,16")
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--support-measure-frac", type=float, default=1.0e-3)
    parser.add_argument("--trusted-erosion", type=int, default=1)
    parser.add_argument("--jtau-thresholds", default="1e-6,1e-4,1e-3,1e-2")
    parser.add_argument("--tilt-max", type=float, default=4.0)
    parser.add_argument("--tilt-step", type=float, default=0.5)
    parser.add_argument("--manual-candidates", default="")
    parser.add_argument("--training-mask", choices=("support", "trusted"), default="trusted")
    parser.add_argument("--cover-threshold", type=float, default=1.0e-3)
    parser.add_argument("--max-patches", type=int, default=12)
    parser.add_argument("--force-lab", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--require-tilde-timelike", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--time-covector-norm-floor", type=float, default=1.0e-8)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
