from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_ab_measurable_observables import crop_inner, extract_centerline_observables, unwrap_phase
from calibrate_exact_kg_ps_fd_absorbing import (
    apply_boundary,
    boundary_exact,
    make_sponge,
    spectral_lap,
)
from exact_kg_2p1_benchmark import Exact2p1Params, integrate_xz
from simulate_bc_psfd_imex import (
    apply_geometry_boundary,
    edge_damp,
    rhs_branch_b,
    spectral_k2,
    wave_imex_step,
)


def rel_l1(a: np.ndarray, b: np.ndarray, scale: float) -> float:
    return float(np.mean(np.abs(a - b)) / max(scale, 1.0e-14))


def first_crossing(times: list[float], values: list[float], threshold: float) -> float | None:
    for t, value in zip(times, values):
        if abs(value) >= threshold:
            return float(t)
    return None


def safe_centerline_summary(x: np.ndarray, z: np.ndarray, a: np.ndarray, b: np.ndarray) -> dict[str, float | None]:
    try:
        obs = extract_centerline_observables(x, z, a, b)
    except Exception:
        return {
            "central_peak_rel_intensity_shift": None,
            "mean_peak_rel_intensity_shift": None,
            "visibility_rel_diff": None,
            "effective_phase_shift_rad": None,
            "equivalent_fringe_shift_over_spacing": None,
        }
    return {
        "central_peak_rel_intensity_shift": float(obs["central_peak_rel_intensity_shift"]),
        "mean_peak_rel_intensity_shift": float(obs["mean_peak_rel_intensity_shift"]),
        "visibility_rel_diff": float(obs["visibility_rel_diff"]),
        "effective_phase_shift_rad": float(obs["effective_phase_shift_rad"]),
        "equivalent_fringe_shift_over_spacing": float(obs["equivalent_fringe_shift_over_spacing"]),
    }


def render_signed_curves(path: Path, ts: np.ndarray, curves: list[tuple[str, np.ndarray]], ylabel: str) -> None:
    plt.figure(figsize=(7.0, 4.4))
    for label, values in curves:
        plt.plot(ts, values, label=label, linewidth=1.7)
    plt.xlabel("t / T_overlap")
    plt.ylabel(ylabel)
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def render_abs_curves(path: Path, ts: np.ndarray, curves: list[tuple[str, np.ndarray]]) -> None:
    plt.figure(figsize=(7.0, 4.4))
    for label, values in curves:
        plt.semilogy(ts, np.maximum(np.abs(values), 1.0e-16), label=label, linewidth=1.7)
    plt.xlabel("t / T_overlap")
    plt.ylabel("absolute magnitude")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def render_rel_diff_montage(path: Path, x: np.ndarray, z: np.ndarray, snapshots: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(len(snapshots), 1, figsize=(5.0, 3.2 * len(snapshots)), squeeze=False)
    vmax = max(float(np.max(np.abs(s["rel_diff"]))) for s in snapshots)
    vmax = max(vmax, 1.0e-12)
    for ax, snap in zip(axes.ravel(), snapshots):
        im = ax.imshow(
            snap["rel_diff"].T,
            origin="lower",
            extent=[x.min(), x.max(), z.min(), z.max()],
            aspect="equal",
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
        )
        ax.set_title(f"t/T={snap['time_fraction']:.2f}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.9, label="(B-A)/A")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(
    outdir: str | Path,
    alpha: float = 0.5,
    lambda_grav: float = 1.0,
    mp: float = 300.0,
    dt_target: float = 0.0125,
    inner: float = 10.0,
    outer: float = 15.0,
    nx: int = 181,
    nz: int = 181,
    nkx: int = 81,
    nkz: int = 81,
    central_disk_radius: float = 1.5,
    bright_fraction: float = 0.5,
    interference_box_halfwidth: float = 2.0,
) -> dict[str, object]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    params = Exact2p1Params(
        alpha=alpha,
        x_half_range=outer,
        z_half_range=outer,
        nx=nx,
        nz=nz,
        nkx=nkx,
        nkz=nkz,
    )
    dx = 2.0 * outer / (params.nx - 1)
    dz = 2.0 * outer / (params.nz - 1)
    total_time = params.overlap_time
    steps = int(round(total_time / dt_target))
    dt = total_time / steps
    k2 = spectral_k2(params.nx, params.nz, dx, dz)

    x, z, psi_prev, _ = integrate_xz(alpha, -dt, params)
    _, _, psi_now, _ = integrate_xz(alpha, 0.0, params)
    sponge_inner = max(inner + 2.0, 0.8 * outer)
    sigma = make_sponge(x, z, inner_half_range=sponge_inner, outer_half_range=outer, sigma_max=2.0)

    tau_b = np.zeros((params.nx, params.nz))
    beta_b = np.zeros_like(tau_b)
    ptau_b = np.zeros_like(tau_b)
    pbeta_b = np.zeros_like(tau_b)

    ix, iz, xi, zi = crop_inner(x, z, inner)
    xi_grid, zi_grid = np.meshgrid(xi, zi, indexing="ij")
    disk_mask = (xi_grid**2 + zi_grid**2) <= central_disk_radius**2
    box_mask = (np.abs(xi_grid) <= interference_box_halfwidth) & (np.abs(zi_grid) <= interference_box_halfwidth)
    midx = len(xi) // 2
    midz = len(zi) // 2

    ts: list[float] = []
    rel_l1_values: list[float] = []
    center_point_shift: list[float] = []
    central_disk_shift: list[float] = []
    bright_region_shift: list[float] = []
    box_shift: list[float] = []
    center_peak_shift: list[float] = []
    mean_peak_shift: list[float] = []
    support_max_shift: list[float] = []
    visibility_shift: list[float] = []
    phase_shift: list[float] = []
    frac_shift: list[float] = []

    snapshot_steps = [0, steps // 4, steps // 2, (3 * steps) // 4, steps]
    snapshots: list[dict[str, object]] = []

    def record(step_index: int, current_time: float) -> None:
        rho_a_outer = np.abs(psi_now) ** 2
        rho_b_outer = rho_a_outer * np.exp(-np.clip(beta_b + tau_b, -12.0, 12.0))
        a = rho_a_outer[np.ix_(ix, iz)]
        b = rho_b_outer[np.ix_(ix, iz)]
        scale = float(np.max(a))
        bright_mask = a >= bright_fraction * scale
        centerline = safe_centerline_summary(xi, zi, a, b)

        ts.append(current_time / total_time)
        rel_l1_values.append(rel_l1(a, b, scale))
        center_point_shift.append(float((b[midx, midz] - a[midx, midz]) / max(a[midx, midz], 1.0e-14)))
        central_disk_shift.append(float((np.sum(b[disk_mask]) - np.sum(a[disk_mask])) / max(np.sum(a[disk_mask]), 1.0e-14)))
        bright_region_shift.append(
            float((np.sum(b[bright_mask]) - np.sum(a[bright_mask])) / max(np.sum(a[bright_mask]), 1.0e-14))
        )
        box_shift.append(float((np.sum(b[box_mask]) - np.sum(a[box_mask])) / max(np.sum(a[box_mask]), 1.0e-14)))

        center_peak_value = centerline["central_peak_rel_intensity_shift"]
        mean_peak_value = centerline["mean_peak_rel_intensity_shift"]
        visibility_value = centerline["visibility_rel_diff"]
        phase_value = centerline["effective_phase_shift_rad"]
        frac_value = centerline["equivalent_fringe_shift_over_spacing"]

        center_peak_shift.append(float(center_peak_value) if center_peak_value is not None else np.nan)
        mean_peak_shift.append(float(mean_peak_value) if mean_peak_value is not None else np.nan)
        visibility_shift.append(float(visibility_value) if visibility_value is not None else np.nan)
        phase_shift.append(float(phase_value) if phase_value is not None else np.nan)
        frac_shift.append(float(frac_value) if frac_value is not None else np.nan)
        support_max_shift.append(
            float(np.max(np.abs((b - a)[a > 0.1 * scale] / np.maximum(a[a > 0.1 * scale], 1.0e-14))))
        )

        if step_index in snapshot_steps:
            snapshots.append(
                {
                    "time_fraction": current_time / total_time,
                    "rel_diff": (b - a) / np.maximum(a, 1.0e-14),
                }
            )

    record(0, 0.0)

    for n in range(steps):
        rho_a_outer = np.abs(psi_now) ** 2
        s_a_outer = unwrap_phase(psi_now)

        src_tau_b, src_beta_b, *_ = rhs_branch_b(
            tau_b, ptau_b, beta_b, pbeta_b, rho_a_outer, s_a_outer, params.m, mp, dx, dz, k2
        )
        tau_b, ptau_b = wave_imex_step(tau_b, ptau_b, lambda_grav * src_tau_b, dt, k2)
        beta_b, pbeta_b = wave_imex_step(beta_b, pbeta_b, lambda_grav * src_beta_b, dt, k2)
        edge_damp(tau_b, sigma, dt, -2.0, 2.0)
        edge_damp(beta_b, sigma, dt, -2.0, 2.0)
        edge_damp(ptau_b, sigma, dt, -10.0, 10.0)
        edge_damp(pbeta_b, sigma, dt, -10.0, 10.0)
        apply_geometry_boundary(tau_b, beta_b)

        t_next = (n + 1) * dt
        psi_b = boundary_exact(alpha, t_next, x, z, params)
        rhs = spectral_lap(psi_now, dx, dz) - params.m**2 * psi_now
        numer = 2.0 * psi_now - (1.0 - 0.5 * sigma * dt) * psi_prev + dt * dt * rhs
        psi_next = numer / (1.0 + 0.5 * sigma * dt)
        apply_boundary(psi_next, psi_b)
        psi_prev, psi_now = psi_now, psi_next

        record(n + 1, t_next)

    ts_array = np.array(ts)
    render_signed_curves(
        out / "signed_intensity_shifts_vs_time.png",
        ts_array,
        [
            ("center point", np.array(center_point_shift)),
            ("central disk", np.array(central_disk_shift)),
            ("bright region", np.array(bright_region_shift)),
            ("interference box", np.array(box_shift)),
        ],
        "relative shift",
    )
    render_abs_curves(
        out / "observable_magnitudes_vs_time.png",
        ts_array,
        [
            ("field relL1", np.array(rel_l1_values)),
            ("center peak", np.array(center_peak_shift)),
            ("mean peak", np.array(mean_peak_shift)),
            ("support max", np.array(support_max_shift)),
        ],
    )
    render_signed_curves(
        out / "phase_visibility_vs_time.png",
        ts_array,
        [
            ("visibility", np.array(visibility_shift)),
            ("phase (rad)", np.array(phase_shift)),
            ("shift / spacing", np.array(frac_shift)),
        ],
        "signed value",
    )
    render_rel_diff_montage(out / "rel_diff_montage.png", xi, zi, snapshots)

    summary = {
        "lambda_grav": lambda_grav,
        "dt": dt,
        "steps": steps,
        "time_series": {
            "t_over_T": [float(v) for v in ts],
            "B_vs_A_relL1": [float(v) for v in rel_l1_values],
            "center_point_rel_shift": [float(v) for v in center_point_shift],
            "central_disk_rel_shift": [float(v) for v in central_disk_shift],
            "bright_region_rel_shift": [float(v) for v in bright_region_shift],
            "interference_box_rel_shift": [float(v) for v in box_shift],
            "central_peak_rel_intensity_shift": [None if np.isnan(v) else float(v) for v in center_peak_shift],
            "mean_peak_rel_intensity_shift": [None if np.isnan(v) else float(v) for v in mean_peak_shift],
            "support_pointwise_rel_diff_max": [float(v) for v in support_max_shift],
            "visibility_rel_diff": [None if np.isnan(v) else float(v) for v in visibility_shift],
            "effective_phase_shift_rad": [None if np.isnan(v) else float(v) for v in phase_shift],
            "equivalent_fringe_shift_over_spacing": [None if np.isnan(v) else float(v) for v in frac_shift],
        },
        "milestones": {
            "B_vs_A_relL1_ge_1e-6": first_crossing(ts, rel_l1_values, 1.0e-6),
            "center_peak_rel_shift_ge_1e-5": first_crossing(ts, center_peak_shift, 1.0e-5),
            "central_disk_rel_shift_ge_1e-5": first_crossing(ts, central_disk_shift, 1.0e-5),
            "bright_region_rel_shift_ge_1e-5": first_crossing(ts, bright_region_shift, 1.0e-5),
        },
        "extrema": {
            "max_abs_B_vs_A_relL1": float(np.max(np.abs(rel_l1_values))),
            "max_abs_center_point_rel_shift": float(np.max(np.abs(center_point_shift))),
            "max_abs_central_disk_rel_shift": float(np.max(np.abs(central_disk_shift))),
            "max_abs_bright_region_rel_shift": float(np.max(np.abs(bright_region_shift))),
            "max_abs_central_peak_rel_shift": float(np.nanmax(np.abs(center_peak_shift))),
            "max_abs_mean_peak_rel_shift": float(np.nanmax(np.abs(mean_peak_shift))),
            "max_abs_support_pointwise_rel_diff_max": float(np.max(np.abs(support_max_shift))),
        },
        "final_values": {
            "B_vs_A_relL1": float(rel_l1_values[-1]),
            "center_point_rel_shift": float(center_point_shift[-1]),
            "central_disk_rel_shift": float(central_disk_shift[-1]),
            "bright_region_rel_shift": float(bright_region_shift[-1]),
            "interference_box_rel_shift": float(box_shift[-1]),
            "central_peak_rel_intensity_shift": None if np.isnan(center_peak_shift[-1]) else float(center_peak_shift[-1]),
            "mean_peak_rel_intensity_shift": None if np.isnan(mean_peak_shift[-1]) else float(mean_peak_shift[-1]),
            "support_pointwise_rel_diff_max": float(support_max_shift[-1]),
            "visibility_rel_diff": None if np.isnan(visibility_shift[-1]) else float(visibility_shift[-1]),
            "effective_phase_shift_rad": None if np.isnan(phase_shift[-1]) else float(phase_shift[-1]),
            "equivalent_fringe_shift_over_spacing": None if np.isnan(frac_shift[-1]) else float(frac_shift[-1]),
        },
        "notes": {
            "purpose": "Track when the AB split turns on during the shared-backbone evolution and whether integrated readouts behave more cleanly than final-time centerline quantities alone.",
            "readout_masks": {
                "central_disk_radius": central_disk_radius,
                "bright_fraction_of_Amax": bright_fraction,
                "interference_box_halfwidth": interference_box_halfwidth,
            },
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = run(Path(__file__).resolve().parent.parent / "visualizations" / "ab_time_observables_lambda1")
    print(json.dumps(result, indent=2))
