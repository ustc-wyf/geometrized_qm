from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from analyze_bcd_residuals_from_a_reference import geometry_data, stress_tensor_tilde
from analyze_c_terms_from_a_reference import metric_jets_full
from coordinate_matter_evolution import (
    conservative_density_from_rho_u,
    coordinate_matter_rhs_covector,
)
from diagnose_d_trace_interface_weak_form import find_level_segments, real_array
from integrate_d_interface_jumps import bilinear_sample, scalar_stats
from mixed_tilde_initial_data import (
    bohm_fields_from_wave_derivatives,
    exact_localized_wave_derivatives,
    localized_direct_tilde_coordinate_initial,
    safe_covariant_from_inverse,
    tilde_inverse_from_flat_mixed_transform,
)
from prototype_d_signed_distance_moving_interface import compute_distance_interface_rows
from simulate_flat_localized_crossing_packets import (
    FlatLocalizedCrossingParams,
    initial_wavefunction,
    make_grid,
    spectral_omega,
)


def ddx(field: np.ndarray, dx: float) -> np.ndarray:
    return (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * dx)


def ddz(field: np.ndarray, dz: float) -> np.ndarray:
    return (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dz)


def lap(field: np.ndarray, dx: float, dz: float) -> np.ndarray:
    return (
        (np.roll(field, -1, axis=0) - 2.0 * field + np.roll(field, 1, axis=0)) / (dx * dx)
        + (np.roll(field, -1, axis=1) - 2.0 * field + np.roll(field, 1, axis=1)) / (dz * dz)
    )


def exact_a_fields(
    psi0_hat: np.ndarray,
    omega: np.ndarray,
    t: float,
    rho_floor: float,
) -> dict[str, np.ndarray]:
    psi = np.fft.ifft2(psi0_hat * np.exp(-1j * omega * float(t)))
    wave = {
        "psi": psi,
        "psi_t": np.fft.ifft2((-1j * omega) * psi0_hat * np.exp(-1j * omega * float(t))),
        "psi_tt": np.fft.ifft2((-(omega**2)) * psi0_hat * np.exp(-1j * omega * float(t))),
        # The exact derivative helper is used elsewhere for full diagnostics;
        # these entries are replaced by the caller when needed.
    }
    rho = np.abs(psi) ** 2
    phase = np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)
    return {"psi": psi, "rho": rho, "S": phase, "rho_safe": np.maximum(rho, rho_floor)}


def bohm_like_from_dynamic_state(
    rho: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    root_pprev: np.ndarray,
    root_prev: np.ndarray,
    root_cur: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    rho_floor: float,
) -> dict[str, np.ndarray]:
    root_safe = np.sqrt(np.maximum(rho, rho_floor))
    r_t = (root_cur - root_prev) / max(dt, 1.0e-300)
    r_tt = (root_cur - 2.0 * root_prev + root_pprev) / max(dt * dt, 1.0e-300)
    r_x = ddx(root_safe, dx)
    r_z = ddz(root_safe, dz)
    q_field = (r_tt - lap(root_safe, dx, dz)) / np.maximum(root_safe, np.sqrt(rho_floor))

    x_field = u_t * u_t - u_x * u_x - u_z * u_z
    y_field = u_t * r_t - u_x * r_x - u_z * r_z
    z_field = r_t * r_t - r_x * r_x - r_z * r_z
    delta_field = x_field * z_field - y_field * y_field
    return {
        "rho": rho,
        "rho_safe": np.maximum(rho, rho_floor),
        "root_rho": root_safe,
        "s_t": u_t,
        "s_x": u_x,
        "s_z": u_z,
        "r_t": r_t,
        "r_x": r_x,
        "r_z": r_z,
        "Q": np.real_if_close(q_field),
        "X": np.real_if_close(x_field),
        "Y": np.real_if_close(y_field),
        "Z": np.real_if_close(z_field),
        "Delta": np.real_if_close(delta_field),
    }


def reconstruct_metric_from_state(
    rho: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    root_pprev: np.ndarray,
    root_prev: np.ndarray,
    root_cur: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
    rho_floor: float,
) -> dict[str, np.ndarray]:
    bohm = bohm_like_from_dynamic_state(
        rho=rho,
        u_t=u_t,
        u_x=u_x,
        u_z=u_z,
        root_pprev=root_pprev,
        root_prev=root_prev,
        root_cur=root_cur,
        dt=dt,
        dx=dx,
        dz=dz,
        rho_floor=rho_floor,
    )
    inv = tilde_inverse_from_flat_mixed_transform(bohm, x_floor=1.0e-10, pinv_rcond=1.0e-10)
    cov = safe_covariant_from_inverse(
        inv["tilde_ginv_txz"],
        rho,
        det_floor=1.0e-12,
        tail_rel_cut=1.0e-8,
        pinv_rcond=1.0e-12,
    )
    return {
        "bohm": bohm,
        "metric_inv": real_array(inv["tilde_ginv_txz"]),
        "metric_cov": real_array(cov["cov_txz"]),
        "det": real_array(cov["det_txz"]),
        "singular_mask": cov["singular_mask"],
        "tail_mask": cov["tail_mask"],
    }


def sanitize_metric_tail(metric_cov: np.ndarray, active_mask: np.ndarray) -> np.ndarray:
    """Replace low-density/tail metric data by flat metric.

    This is a numerical domain restriction: the D comparison is only trusted
    on the matter support.  Outside the active mask the localized wave packet
    carries negligible density, while the disformal reconstruction is often
    ill-conditioned.
    """
    out = np.asarray(metric_cov, dtype=float).copy()
    eta_cov = np.zeros((3, 3), dtype=float)
    eta_cov[0, 0] = 1.0
    eta_cov[1, 1] = -1.0
    eta_cov[2, 2] = -1.0
    bad = (~active_mask) | (~np.all(np.isfinite(out), axis=(-1, -2)))
    out[bad] = eta_cov
    return out


def rk4_matter_step(
    n_cons: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    metric_cov: np.ndarray,
    mass: float,
    dx: float,
    dz: float,
    dt: float,
    u_t_reference: np.ndarray,
    active_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    def rusanov_density_rhs(rec: dict[str, np.ndarray], n0: np.ndarray) -> np.ndarray:
        n_safe = np.where(np.abs(n0) > 1.0e-300, n0, np.sign(n0) * 1.0e-300 + (n0 == 0.0) * 1.0e-300)
        vx = rec["flux_x"] / n_safe
        vz = rec["flux_z"] / n_safe
        vx = np.where(active_mask, vx, 0.0)
        vz = np.where(active_mask, vz, 0.0)

        def flux_1d(q_l: np.ndarray, q_r: np.ndarray, v_l: np.ndarray, v_r: np.ndarray) -> np.ndarray:
            speed = np.maximum(np.abs(v_l), np.abs(v_r))
            return 0.5 * (q_l * v_l + q_r * v_r) - 0.5 * speed * (q_r - q_l)

        fx = flux_1d(n0[:-1, :], n0[1:, :], vx[:-1, :], vx[1:, :])
        fz = flux_1d(n0[:, :-1], n0[:, 1:], vz[:, :-1], vz[:, 1:])
        n_t = np.zeros_like(n0)
        n_t[1:-1, 1:-1] -= (fx[1:, 1:-1] - fx[:-1, 1:-1]) / dx
        n_t[1:-1, 1:-1] -= (fz[1:-1, 1:] - fz[1:-1, :-1]) / dz
        return np.where(active_mask, n_t, 0.0)

    def rhs(n0: np.ndarray, ux0: np.ndarray, uz0: np.ndarray, ut_ref: np.ndarray):
        n_eff = np.where(active_mask, n0, 0.0)
        out = coordinate_matter_rhs_covector(
            n_cons=n_eff,
            u_x=ux0,
            u_z=uz0,
            metric_cov_txz=metric_cov,
            mass=mass,
            dx=dx,
            dz=dz,
            branch="negative_frequency",
            u_t_reference=ut_ref,
        )
        n_t = rusanov_density_rhs(out, n_eff)
        return n_t, out["u_x_t"], out["u_z_t"], out

    k1_n, k1_x, k1_z, o1 = rhs(n_cons, u_x, u_z, u_t_reference)
    k2_n, k2_x, k2_z, _ = rhs(n_cons + 0.5 * dt * k1_n, u_x + 0.5 * dt * k1_x, u_z + 0.5 * dt * k1_z, o1["u_t"])
    k3_n, k3_x, k3_z, _ = rhs(n_cons + 0.5 * dt * k2_n, u_x + 0.5 * dt * k2_x, u_z + 0.5 * dt * k2_z, o1["u_t"])
    k4_n, k4_x, k4_z, _ = rhs(n_cons + dt * k3_n, u_x + dt * k3_x, u_z + dt * k3_z, o1["u_t"])

    n_next = n_cons + (dt / 6.0) * (k1_n + 2.0 * k2_n + 2.0 * k3_n + k4_n)
    ux_next = u_x + (dt / 6.0) * (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x)
    uz_next = u_z + (dt / 6.0) * (k1_z + 2.0 * k2_z + 2.0 * k3_z + k4_z)
    n_next = np.where(active_mask, n_next, 0.0)
    final = coordinate_matter_rhs_covector(
        n_cons=n_next,
        u_x=ux_next,
        u_z=uz_next,
        metric_cov_txz=metric_cov,
        mass=mass,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=o1["u_t"],
    )
    n_next = np.where(active_mask, n_next, 0.0)
    ux_next = np.where(active_mask, ux_next, 0.0)
    uz_next = np.where(active_mask, uz_next, 0.0)
    for key in ("rho", "u_t", "discriminant", "j_t", "j_x", "j_z"):
        if key in final:
            final[key] = np.where(active_mask, final[key], 0.0)
    if not (
        np.all(np.isfinite(n_next[active_mask]))
        and np.all(np.isfinite(ux_next[active_mask]))
        and np.all(np.isfinite(uz_next[active_mask]))
        and np.all(np.isfinite(final["rho"][active_mask]))
        and np.all(np.isfinite(final["u_t"][active_mask]))
    ):
        raise FloatingPointError("non-finite value produced by D matter RK4 step")
    return n_next, ux_next, uz_next, final


def metric_curvature_from_history(
    metric_pprev: np.ndarray,
    metric_prev: np.ndarray,
    metric_cur: np.ndarray,
    dt: float,
    dx: float,
    dz: float,
) -> dict[str, np.ndarray]:
    # geometry_data expects m,0,p central data.  With available history, this
    # evaluates curvature on metric_prev.  It is used as a lagged diagnostic and
    # interface extraction field, not as an implicit solve.
    dg, d2g = metric_jets_full(metric_pprev, metric_prev, metric_cur, dt, dx, dz)
    ginv, _, ricci, r_scalar = geometry_data(metric_prev, dg, d2g)
    return {
        "metric_cov": metric_prev,
        "metric_inv": real_array(ginv),
        "ricci": real_array(ricci),
        "R_tilde": real_array(r_scalar),
    }


def stress_trace_from_state(
    metric_cov: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    u_t: np.ndarray,
    u_x: np.ndarray,
    u_z: np.ndarray,
    mass: float,
    mp: float,
) -> np.ndarray:
    stress, _ = stress_tensor_tilde(metric_cov, metric_inv, rho, u_t, u_x, u_z, m=mass)
    rhs = real_array(stress) / (mp * mp)
    return np.einsum("...ab,...ab->...", metric_inv, rhs, optimize=True)


def interface_jump_diagnostics(
    raw_y: np.ndarray,
    metric_inv: np.ndarray,
    rho: np.ndarray,
    stress_trace: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    ell: float,
    half_width_y: float,
    samples: int,
) -> dict[str, object]:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    raw_x, raw_z = np.gradient(raw_y, dx, dz, edge_order=2)
    alpha = np.sqrt(raw_x * raw_x + raw_z * raw_z)
    rows: list[dict[str, float]] = []
    zero_dt = np.zeros_like(raw_y)
    for level in (-1.0, 1.0):
        segments = find_level_segments(raw_y, rho, x, z, level)
        if not segments:
            continue
        # Use signed distance initialized from the current raw_y-level segments.
        # This converts the trace jump law into a normal speed diagnostic.
        distance = signed_distance_from_segments(raw_y - level, x, z, segments)
        rows.extend(
            compute_distance_interface_rows(
                distance_field=distance,
                raw_y_ref_grad_norm=alpha,
                metric_inv=metric_inv,
                rho=rho,
                stress_trace=stress_trace,
                x=x,
                z=z,
                branch_level=level,
                ell=ell,
                half_width_y=half_width_y,
                samples=samples,
                dt_reference_grid=zero_dt,
            )
        )
    solved = [row for row in rows if row.get("solved", 0.0) > 0.5]
    return {
        "rows": rows,
        "segment_count": len(rows),
        "solved_count": len(solved),
        "solved_fraction": float(len(solved) / max(len(rows), 1)),
        "speed": scalar_stats([row.get("coordinate_normal_speed_corrected", np.nan) for row in solved]),
    }


def point_segment_distance_squared(px: np.ndarray, pz: np.ndarray, seg: dict[str, float]) -> np.ndarray:
    vx = seg["x1"] - seg["x0"]
    vz = seg["z1"] - seg["z0"]
    denom = max(vx * vx + vz * vz, 1.0e-300)
    tau = ((px - seg["x0"]) * vx + (pz - seg["z0"]) * vz) / denom
    tau = np.clip(tau, 0.0, 1.0)
    cx = seg["x0"] + tau * vx
    cz = seg["z0"] + tau * vz
    return (px - cx) ** 2 + (pz - cz) ** 2


def signed_distance_from_segments(sign_source: np.ndarray, x: np.ndarray, z: np.ndarray, segments: list[dict[str, float]]) -> np.ndarray:
    if not segments:
        return sign_source.copy()
    xg, zg = np.meshgrid(x, z, indexing="ij")
    dist2 = np.full_like(sign_source, np.inf, dtype=float)
    for seg in segments:
        dist2 = np.minimum(dist2, point_segment_distance_squared(xg, zg, seg))
    return np.where(sign_source >= 0.0, 1.0, -1.0) * np.sqrt(dist2)


def row_segments(rows: list[dict[str, float]]) -> list[list[tuple[float, float]]]:
    return [[(row["x0"], row["z0"]), (row["x1"], row["z1"])] for row in rows]


def add_lines(ax: plt.Axes, rows: list[dict[str, float]], *, colors: str = "white", linewidth: float = 0.8) -> None:
    if not rows:
        return
    ax.add_collection(LineCollection(row_segments(rows), colors=colors, linewidths=linewidth))


def relative_l1(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    denom = max(float(np.mean(np.abs(a[mask]))), 1.0e-300)
    return float(np.mean(np.abs(a[mask] - b[mask])) / denom)


def render_summary(
    out_path: Path,
    x: np.ndarray,
    z: np.ndarray,
    rho_a: np.ndarray,
    rho_d: np.ndarray,
    raw_y: np.ndarray | None,
    interface_rows: list[dict[str, float]],
    records: list[dict[str, float]],
    support: np.ndarray,
) -> None:
    xg, zg = np.meshgrid(x, z, indexing="ij")
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 9.2), constrained_layout=True)
    vmax = max(float(np.percentile(rho_a[support], 99.5)), float(np.percentile(rho_d[support], 99.5)), 1.0e-16)
    im0 = axes[0, 0].pcolormesh(xg, zg, rho_a, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
    axes[0, 0].set_title("A exact rho")
    fig.colorbar(im0, ax=axes[0, 0])
    im1 = axes[0, 1].pcolormesh(xg, zg, rho_d, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
    add_lines(axes[0, 1], interface_rows, colors="white", linewidth=0.65)
    axes[0, 1].set_title("D reduced rho + current interfaces")
    fig.colorbar(im1, ax=axes[0, 1])
    diff = rho_d - rho_a
    dmax = max(float(np.percentile(np.abs(diff[support]), 99.0)), 1.0e-16)
    im2 = axes[0, 2].pcolormesh(xg, zg, diff, shading="auto", cmap="coolwarm", vmin=-dmax, vmax=dmax)
    axes[0, 2].set_title("rho_D - rho_A")
    fig.colorbar(im2, ax=axes[0, 2])

    if raw_y is not None:
        clip = max(float(np.percentile(np.abs(raw_y[support]), 95.0)), 1.0)
        im3 = axes[1, 0].pcolormesh(xg, zg, np.clip(raw_y, -clip, clip), shading="auto", cmap="coolwarm", vmin=-clip, vmax=clip)
        add_lines(axes[1, 0], interface_rows, colors="black", linewidth=0.55)
        axes[1, 0].set_title("lagged raw_y=ell^2 R_tilde")
        fig.colorbar(im3, ax=axes[1, 0])
    times = np.asarray([r["t"] for r in records], dtype=float)
    axes[1, 1].plot(times, [r["rho_rel_l1_support"] for r in records], marker="o", label="rho rel L1")
    axes[1, 1].plot(times, [r["phase_grad_rel_l1_support"] for r in records], marker="o", label="grad S rel L1")
    axes[1, 1].set_title("D vs A observables")
    axes[1, 1].set_xlabel("t")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.25)
    axes[1, 2].plot(times, [r["interface_solved_fraction"] for r in records], marker="o", label="trace jump solved fraction")
    axes[1, 2].plot(times, [r["disc_min_support"] for r in records], marker="o", label="mass-shell discriminant min")
    axes[1, 2].set_title("D diagnostic")
    axes[1, 2].set_xlabel("t")
    axes[1, 2].legend()
    axes[1, 2].grid(alpha=0.25)
    for ax in axes[:1, :].ravel().tolist() + [axes[1, 0]]:
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    params = FlatLocalizedCrossingParams(nx=args.resolution, nz=args.resolution)
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
        probe_dt=args.dt,
        rho_floor=args.rho_floor,
    )
    rho0 = real_array(init["rho0"])
    diag = init["diagnostics"]
    ux = real_array(diag["s_x0"])
    uz = real_array(diag["s_z0"])
    ut = real_array(diag["s_t0"])
    support = rho0 > args.support_rho_frac * float(np.max(rho0))
    active = support.copy()
    for _ in range(args.active_dilation):
        active = active | np.roll(active, 1, 0) | np.roll(active, -1, 0) | np.roll(active, 1, 1) | np.roll(active, -1, 1)
    metric_cov = sanitize_metric_tail(real_array(init["metric_cov_txz_0"]), active)
    metric_prev = metric_cov.copy()
    metric_pprev = sanitize_metric_tail(real_array(init["metric_cov_txz_0"] - args.dt * init["metric_cov_txz_t0"]), active)

    cons = conservative_density_from_rho_u(
        rho=rho0,
        u_x=ux,
        u_z=uz,
        metric_cov_txz=metric_cov,
        mass=params.m,
        branch="negative_frequency",
        u_t_reference=ut,
        rho_floor=args.rho_floor,
    )
    n_cons = np.where(active, cons["n_cons"], 0.0)
    current = coordinate_matter_rhs_covector(
        n_cons=n_cons,
        u_x=ux,
        u_z=uz,
        metric_cov_txz=metric_cov,
        mass=params.m,
        dx=dx,
        dz=dz,
        branch="negative_frequency",
        u_t_reference=ut,
    )

    wave_m = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x, z, -args.dt)
    root_pprev = bohm_fields_from_wave_derivatives(wave_m, rho_floor=args.rho_floor)["root_rho"]
    root_prev = np.sqrt(np.maximum(rho0, args.rho_floor))
    root_cur = root_prev.copy()

    records: list[dict[str, float]] = []
    last_raw_y: np.ndarray | None = None
    last_rows: list[dict[str, float]] = []

    stopped_reason = "completed"
    completed_step = 0
    for step in range(args.steps + 1):
        completed_step = step
        t = step * args.dt
        a_wave = exact_localized_wave_derivatives(psi0, psi0_hat, omega, x, z, t)
        a_bohm = bohm_fields_from_wave_derivatives(a_wave, rho_floor=args.rho_floor)
        rho_a = real_array(a_bohm["rho"])
        rho_d = real_array(current["rho"])
        sx_a = real_array(a_bohm["s_x"])
        sz_a = real_array(a_bohm["s_z"])
        sx_d = real_array(ux)
        sz_d = real_array(uz)

        raw_y = None
        interface = {"rows": [], "solved_fraction": 0.0, "segment_count": 0, "solved_count": 0}
        do_interface_diagnostic = step >= 2 and (
            args.interface_every > 0
            and (step % args.interface_every == 0 or step == args.steps)
        )
        if do_interface_diagnostic:
            geom = metric_curvature_from_history(metric_pprev, metric_prev, metric_cov, args.dt, dx, dz)
            raw_y = args.ell * args.ell * geom["R_tilde"]
            stress_trace = stress_trace_from_state(
                geom["metric_cov"],
                geom["metric_inv"],
                rho_d,
                current["u_t"],
                ux,
                uz,
                params.m,
                args.mp,
            )
            interface = interface_jump_diagnostics(
                raw_y=raw_y,
                metric_inv=geom["metric_inv"],
                rho=rho_d,
                stress_trace=stress_trace,
                x=x,
                z=z,
                ell=args.ell,
                half_width_y=args.half_width_y,
                samples=args.samples,
            )
            last_raw_y = raw_y
            last_rows = interface["rows"]

        grad_a = np.sqrt(sx_a * sx_a + sz_a * sz_a)
        grad_d = np.sqrt(sx_d * sx_d + sz_d * sz_d)
        records.append(
            {
                "step": int(step),
                "t": float(t),
                "rho_rel_l1_support": relative_l1(rho_a, rho_d, support),
                "rho_max_abs_support": float(np.max(np.abs((rho_d - rho_a)[support]))),
                "phase_grad_rel_l1_support": relative_l1(grad_a, grad_d, support),
                "disc_min_support": float(np.nanmin(current["discriminant"][support])),
                "rho_d_min_support": float(np.nanmin(rho_d[support])),
                "rho_d_max_support": float(np.nanmax(rho_d[support])),
                "interface_segment_count": int(interface["segment_count"]),
                "interface_solved_count": int(interface["solved_count"]),
                "interface_solved_fraction": float(interface["solved_fraction"]),
            }
        )
        disc_min_support = float(np.nanmin(current["discriminant"][support]))
        if args.stop_on_negative_discriminant and disc_min_support < -abs(args.disc_tolerance):
            stopped_reason = f"negative mass-shell discriminant on support: {disc_min_support:.6e}"
            break
        if step == args.steps:
            break

        try:
            n_cons, ux, uz, current = rk4_matter_step(
                n_cons=n_cons,
                u_x=ux,
                u_z=uz,
                metric_cov=metric_cov,
                mass=params.m,
                dx=dx,
                dz=dz,
                dt=args.dt,
                u_t_reference=current["u_t"],
                active_mask=active,
            )
        except FloatingPointError as exc:
            stopped_reason = str(exc)
            break
        rho_new = real_array(current["rho"])
        root_next = np.sqrt(np.maximum(rho_new, args.rho_floor))

        if args.metric_mode == "instant_transform":
            try:
                reconstructed = reconstruct_metric_from_state(
                    rho=rho_new,
                    u_t=current["u_t"],
                    u_x=ux,
                    u_z=uz,
                    root_pprev=root_pprev,
                    root_prev=root_prev,
                    root_cur=root_next,
                    dt=args.dt,
                    dx=dx,
                    dz=dz,
                    rho_floor=args.rho_floor,
                )
            except Exception as exc:
                stopped_reason = f"metric reconstruction failed: {type(exc).__name__}: {exc}"
                break
            metric_pprev = metric_prev
            metric_prev = metric_cov
            metric_cov = sanitize_metric_tail(reconstructed["metric_cov"], active)
            if not np.all(np.isfinite(metric_cov)):
                stopped_reason = "metric reconstruction produced non-finite metric_cov"
                break
        elif args.metric_mode == "frozen_initial":
            metric_pprev = metric_prev
            metric_prev = metric_cov
            metric_cov = metric_cov
        else:
            raise ValueError(args.metric_mode)
        root_pprev, root_prev, root_cur = root_prev, root_next, root_next

    t_completed = completed_step * args.dt
    final_a = exact_a_fields(psi0_hat, omega, t_completed, args.rho_floor)
    fig_path = out / f"d_reduced_dynamic_vs_a_{args.metric_mode}_ell{args.ell:g}_n{args.resolution}_t{t_completed:g}.png"
    render_summary(
        fig_path,
        x=x,
        z=z,
        rho_a=final_a["rho"],
        rho_d=real_array(current["rho"]),
        raw_y=last_raw_y,
        interface_rows=last_rows,
        records=records,
        support=support,
    )
    npz_path = out / "fields_final.npz"
    np.savez_compressed(
        npz_path,
        x=x,
        z=z,
        rho_A=final_a["rho"],
        rho_D=real_array(current["rho"]),
        u_x=ux,
        u_z=uz,
        metric_cov=metric_cov,
        raw_y=np.zeros_like(rho0) if last_raw_y is None else last_raw_y,
    )

    summary = {
        "params": {
            "branch": "D",
            "scope": "Reduced dynamic D prototype: same A initial matter data; D matter equations evolve rho/S covector on a tilde metric; metric is either frozen initial or instantaneously reconstructed from current matter. Interface jump law is evaluated dynamically as a diagnostic, not yet a full independent tilde-g bulk solve.",
            "metric_mode": args.metric_mode,
            "ell": args.ell,
            "mp": args.mp,
            "resolution": args.resolution,
            "dt": args.dt,
            "steps": args.steps,
            "steps_completed": completed_step,
            "stopped_reason": stopped_reason,
            "t_requested": args.steps * args.dt,
            "t_completed": t_completed,
            "rho_floor": args.rho_floor,
            "support_rho_frac": args.support_rho_frac,
            "active_dilation": args.active_dilation,
            "interface_every": args.interface_every,
            "stop_on_negative_discriminant": args.stop_on_negative_discriminant,
            "disc_tolerance": args.disc_tolerance,
        },
        "definitions": {
            "rho_rel_l1_support": "Mean |rho_D-rho_A| divided by mean |rho_A| on the initial support mask.",
            "phase_grad_rel_l1_support": "Mean difference of |grad S| divided by A value on the initial support mask.",
            "interface_solved_fraction": "Trace jump-law speed solvability fraction on the dynamically extracted lagged raw_y=+/-1 interfaces. This is not yet full tensor jump feedback.",
            "important_caveat": "This is the first executable dynamic step beyond A-reference residual tests, not the final D branch with independent tilde-g variational evolution.",
        },
        "final": records[-1],
        "records": records,
        "files": {
            "figure": str(fig_path.resolve()),
            "fields_npz": str(npz_path.resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metric-mode", choices=["frozen_initial", "instant_transform"], default="instant_transform")
    parser.add_argument("--ell", type=float, default=30.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--dt", type=float, default=2.5e-4)
    parser.add_argument("--steps", type=int, default=40)
    parser.add_argument("--rho-floor", type=float, default=1.0e-12)
    parser.add_argument("--support-rho-frac", type=float, default=1.0e-3)
    parser.add_argument("--active-dilation", type=int, default=0)
    parser.add_argument(
        "--interface-every",
        type=int,
        default=5,
        help="Compute interface/jump diagnostics every N steps. Matter evolution is still performed every step.",
    )
    parser.add_argument("--disc-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--stop-on-negative-discriminant", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--half-width-y", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=61)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
