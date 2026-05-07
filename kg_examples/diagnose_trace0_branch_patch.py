from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diagnose_mathcal_r_pure_geometry import build_case, symmetric_rows, tensor_norm
from fit_gbcd_trace0_sparse_conservation import (
    TIME_KEYS,
    build_basis,
    build_column_index,
    coeff_to_lambda_fields,
    make_active_mask,
    tensor_fields,
)
from fit_metric_fr_ricci2_universal_full import build_case as build_full_case
from physical_units import HBAR_C_EV_M, PLANCK_MASS_EV
from simulate_d_local_window_from_a_snapshot import build_physical_reference


SYMMETRIC_COMPONENTS = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


def tau_tag(tau: float) -> str:
    if abs(tau) < 5.0e-13:
        return "tau0"
    sign = "p" if tau > 0.0 else "m"
    return "tau" + sign + f"{abs(tau):g}".replace(".", "p")


def coeff_path_for(args: argparse.Namespace, tau: float) -> Path:
    if args.coeff_template:
        return Path(str(args.coeff_template).format(tau=tau, tag=tau_tag(tau)))
    return (
        args.coeff_root
        / f"equation_first_gbcd_trace0_hard_alm_n{args.full_resolution}_{tau_tag(tau)}_{args.fit_region}"
        / "trace0_sparse_conservation_coefficients.npz"
    )


def finite_stats(values: np.ndarray, weights: np.ndarray | None = None) -> dict[str, float]:
    vals = np.asarray(values, dtype=float)
    mask = np.isfinite(vals)
    vals = vals[mask]
    if vals.size == 0:
        return {"count": 0, "mean": 0.0, "p50": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0}
    out = {
        "count": int(vals.size),
        "mean": float(np.mean(vals)),
        "p50": float(np.percentile(vals, 50.0)),
        "p95": float(np.percentile(vals, 95.0)),
        "p99": float(np.percentile(vals, 99.0)),
        "max": float(np.max(vals)),
    }
    if weights is not None:
        w_all = np.asarray(weights, dtype=float)
        w = np.maximum(w_all[mask], 0.0)
        w_sum = float(np.sum(w))
        if w_sum > 0.0:
            order = np.argsort(vals)
            sorted_vals = vals[order]
            sorted_w = w[order]
            cdf = np.cumsum(sorted_w) / w_sum
            out["weighted_mean"] = float(np.sum(vals * w) / w_sum)
            out["weighted_p50"] = float(sorted_vals[np.searchsorted(cdf, 0.50, side="left")])
            out["weighted_p95"] = float(sorted_vals[np.searchsorted(cdf, 0.95, side="left")])
    return out


def raise_covector(metric_inv: np.ndarray, cov: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...b->...a", metric_inv, cov, optimize=True)


def scalar_contract(metric_inv: np.ndarray, a_cov: np.ndarray, b_cov: np.ndarray) -> np.ndarray:
    return np.einsum("...ab,...a,...b->...", metric_inv, a_cov, b_cov, optimize=True)


def qtilde_rho(cases: dict[int, object], gamma: np.ndarray, dt: float, dx: float, dz: float) -> np.ndarray:
    r_up = {key: raise_covector(cases[key].metric_inv, cases[key].r_cov) for key in TIME_KEYS}
    div = (r_up[1][..., 0] - r_up[-1][..., 0]) / (2.0 * dt)
    div = div + np.gradient(r_up[0][..., 1], dx, axis=0, edge_order=2)
    div = div + np.gradient(r_up[0][..., 2], dz, axis=1, edge_order=2)
    conn = np.zeros_like(div)
    for mu in range(3):
        for lam in range(3):
            conn += gamma[..., mu, mu, lam] * r_up[0][..., lam]
    r2 = scalar_contract(cases[0].metric_inv, cases[0].r_cov, cases[0].r_cov)
    return div + conn + r2


def symmetric_vector(tensor: np.ndarray) -> np.ndarray:
    return np.asarray([tensor[a, b] for a, b in SYMMETRIC_COMPONENTS], dtype=float)


def pointwise_span_residual(
    target: np.ndarray,
    bases: list[np.ndarray],
    mask: np.ndarray,
    *,
    target_floor: float,
) -> np.ndarray:
    rel = np.full(mask.shape, np.nan, dtype=float)
    for i_raw, j_raw in np.argwhere(mask):
        i = int(i_raw)
        j = int(j_raw)
        b = symmetric_vector(target[i, j])
        a = np.stack([symmetric_vector(base[i, j]) for base in bases], axis=1)
        if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            continue
        scale = np.linalg.norm(a, axis=0)
        keep = scale > 1.0e-300
        if not np.any(keep):
            continue
        a = a[:, keep]
        scale = scale[keep]
        try:
            coeff_scaled, *_ = np.linalg.lstsq(a / scale[None, :], b, rcond=1.0e-12)
        except np.linalg.LinAlgError:
            continue
        pred = a @ (coeff_scaled / scale)
        denom = max(float(np.linalg.norm(b)), float(target_floor), 1.0e-300)
        rel[i, j] = float(np.linalg.norm(pred - b) / denom)
    return rel


def reconstruct_c_field(args: argparse.Namespace, ref: dict[str, object], tau: float) -> tuple[dict[int, object], np.ndarray, np.ndarray]:
    probe = float(args.force_probe_dt_old)
    cases = {
        -1: build_case(args, ref, tau - probe),
        0: build_case(args, ref, tau),
        1: build_case(args, ref, tau + probe),
    }
    active = make_active_mask(cases, args.fit_region, args.force_region, int(args.active_dilation))
    col_index, ncols = build_column_index(active)
    basis = build_basis(cases, active)
    path = coeff_path_for(args, tau)
    if not path.exists():
        raise FileNotFoundError(f"missing coefficient file for tau={tau:g}: {path}")
    z = np.load(path)
    if args.coeff_key not in z.files:
        raise KeyError(f"{path} does not contain {args.coeff_key}; available={z.files}")
    coeff = np.asarray(z[args.coeff_key], dtype=float)
    if coeff.size != ncols:
        raise ValueError(f"coefficient size mismatch for tau={tau:g}: file has {coeff.size}, rebuilt ncols={ncols}")
    lambdas = coeff_to_lambda_fields(coeff, col_index, basis, cases[0].rho_a.shape)
    cov = tensor_fields(cases, lambdas)
    return cases, cov[0], active


def analyze_one(args: argparse.Namespace, ref: dict[str, object], tau: float) -> dict[str, object]:
    cases, c_cov, active = reconstruct_c_field(args, ref, tau)
    full = build_full_case(args, ref, tau)
    scale = ref["scale"]
    params = ref["params"]
    dt = float(args.force_probe_dt_old) * float(scale.old_dimensionless_scale_ev_inv)
    dx = float(cases[0].x[1] - cases[0].x[0])
    dz = float(cases[0].z[1] - cases[0].z[0])
    q_raw = qtilde_rho(cases, full.geom_0.gamma, dt, dx, dz)
    q_dimless = q_raw / max(float(params.m) ** 2, 1.0e-300)
    region = getattr(cases[0], args.diagnostic_region)
    weight = np.sqrt(np.maximum(cases[0].rho_a, 0.0) / max(float(np.max(cases[0].rho_a)), 1.0e-300))
    q_abs_region = np.abs(q_dimless[region & np.isfinite(q_dimless)])
    q95 = max(float(np.percentile(q_abs_region, 95.0)) if q_abs_region.size else 0.0, 1.0e-300)
    q_rel = np.abs(q_dimless) / q95
    chi = q_rel * q_rel / (q_rel * q_rel + float(args.chi_q0) ** 2)
    c_norm = tensor_norm(c_cov)
    chat_norm = c_norm / np.maximum(chi, float(args.chi_floor))

    u2 = scalar_contract(cases[0].metric_inv, cases[0].u_cov, cases[0].u_cov)
    r2 = scalar_contract(cases[0].metric_inv, cases[0].r_cov, cases[0].r_cov)
    ur = scalar_contract(cases[0].metric_inv, cases[0].u_cov, cases[0].r_cov)
    delta = u2 * r2 - ur * ur
    delta_rel = np.abs(delta) / np.maximum(np.abs(u2 * r2) + np.abs(ur * ur), 1.0e-300)
    r_parallel = (ur / np.maximum(u2, 1.0e-300))[..., None] * cases[0].u_cov
    r_perp = cases[0].r_cov - r_parallel
    rperp2 = scalar_contract(cases[0].metric_inv, r_perp, r_perp)
    rperp_rel = np.sqrt(np.abs(rperp2)) / np.maximum(np.sqrt(np.abs(r2)), 1.0e-300)

    target_norm_region = tensor_norm(cases[0].r_need)[region & np.isfinite(tensor_norm(cases[0].r_need))]
    target_floor = float(args.target_floor_frac) * max(
        float(np.percentile(target_norm_region, 95.0)) if target_norm_region.size else 0.0,
        1.0e-300,
    )
    uu = np.einsum("...a,...b->...ab", cases[0].u_cov, cases[0].u_cov, optimize=True)
    h_trace0 = uu - (u2 / 3.0)[..., None, None] * cases[0].metric_cov
    eu2_rel = pointwise_span_residual(cases[0].r_need, [cases[0].metric_cov, uu], region, target_floor=target_floor)
    eu_trace0_rel = pointwise_span_residual(cases[0].r_need, [h_trace0], region, target_floor=target_floor)

    branch_bins: dict[str, object] = {}
    for threshold in [float(x) for x in args.qrel_thresholds.split(",") if x.strip()]:
        mask = region & (q_rel <= threshold)
        branch_bins[f"qrel_le_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask)),
            "q_abs": finite_stats(np.abs(q_dimless[mask]), weight[mask]),
            "C_norm": finite_stats(c_norm[mask], weight[mask]),
            "C_over_chi_norm": finite_stats(chat_norm[mask], weight[mask]),
        }
    for threshold in [float(x) for x in args.qabs_thresholds.split(",") if x.strip()]:
        mask = region & (np.abs(q_dimless) <= threshold)
        branch_bins[f"qabs_le_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask)),
            "q_abs": finite_stats(np.abs(q_dimless[mask]), weight[mask]),
            "C_norm": finite_stats(c_norm[mask], weight[mask]),
            "C_over_chi_norm": finite_stats(chat_norm[mask], weight[mask]),
        }

    patch_bins: dict[str, object] = {}
    for threshold in [float(x) for x in args.delta_thresholds.split(",") if x.strip()]:
        mask = region & (delta_rel <= threshold)
        patch_bins[f"delta_rel_le_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask)),
            "delta_rel": finite_stats(delta_rel[mask], weight[mask]),
            "rperp_rel": finite_stats(rperp_rel[mask], weight[mask]),
            "Eu2_residual": finite_stats(eu2_rel[mask], weight[mask]),
            "Eu_trace0_residual": finite_stats(eu_trace0_rel[mask], weight[mask]),
        }
    for threshold in [float(x) for x in args.rperp_thresholds.split(",") if x.strip()]:
        mask = region & (rperp_rel <= threshold)
        patch_bins[f"rperp_rel_le_{threshold:g}"] = {
            "count": int(np.count_nonzero(mask)),
            "delta_rel": finite_stats(delta_rel[mask], weight[mask]),
            "rperp_rel": finite_stats(rperp_rel[mask], weight[mask]),
            "Eu2_residual": finite_stats(eu2_rel[mask], weight[mask]),
            "Eu_trace0_residual": finite_stats(eu_trace0_rel[mask], weight[mask]),
        }

    map_path = args.output / f"branch_patch_{tau_tag(tau)}.png"
    render_maps(
        map_path,
        cases[0],
        region,
        active,
        q_rel,
        c_norm,
        chat_norm,
        delta_rel,
        rperp_rel,
        eu2_rel,
        eu_trace0_rel,
    )

    return {
        "tau": float(tau),
        "coefficient_path": str(coeff_path_for(args, tau).resolve()),
        "map_png": str(map_path.resolve()),
        "shape": [int(cases[0].rho_a.shape[0]), int(cases[0].rho_a.shape[1])],
        "counts": {
            "support": int(np.count_nonzero(cases[0].support)),
            "trusted": int(np.count_nonzero(cases[0].trusted)),
            "core10": int(np.count_nonzero(cases[0].core10)),
            "active": int(np.count_nonzero(active)),
            "diagnostic_region": int(np.count_nonzero(region)),
        },
        "q": {
            "q95_abs_in_region": float(q95),
            "q_abs_stats": finite_stats(np.abs(q_dimless[region]), weight[region]),
            "qrel_stats": finite_stats(q_rel[region], weight[region]),
        },
        "branch_bins": branch_bins,
        "patch": {
            "delta_rel_stats": finite_stats(delta_rel[region], weight[region]),
            "rperp_rel_stats": finite_stats(rperp_rel[region], weight[region]),
            "Eu2_residual_stats": finite_stats(eu2_rel[region], weight[region]),
            "Eu_trace0_residual_stats": finite_stats(eu_trace0_rel[region], weight[region]),
            "near_degenerate_bins": patch_bins,
        },
    }


def render_maps(
    out_path: Path,
    case,
    region: np.ndarray,
    active: np.ndarray,
    q_rel: np.ndarray,
    c_norm: np.ndarray,
    chat_norm: np.ndarray,
    delta_rel: np.ndarray,
    rperp_rel: np.ndarray,
    eu2_rel: np.ndarray,
    eu_trace0_rel: np.ndarray,
) -> None:
    x_um = case.x * HBAR_C_EV_M * 1.0e6
    z_um = case.z * HBAR_C_EV_M * 1.0e6
    xg, zg = np.meshgrid(x_um, z_um, indexing="ij")
    fields = [
        ("rho_A", case.rho_a, "linear", "viridis"),
        ("log10 q_rel", q_rel, "log", "magma"),
        ("log10 |C|", c_norm, "log", "cividis"),
        ("log10 |C|/chi", chat_norm, "log", "inferno"),
        ("log10 Delta_rel", delta_rel, "log", "plasma"),
        ("log10 r_perp_rel", rperp_rel, "log", "plasma"),
        ("Eu span{g,uu} residual", eu2_rel, "linear", "rocket" if "rocket" in plt.colormaps() else "magma"),
        ("Eu trace0 residual", eu_trace0_rel, "linear", "rocket" if "rocket" in plt.colormaps() else "magma"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(18.0, 8.8), constrained_layout=True)
    for ax, (title, field, scale, cmap) in zip(axes.flat, fields):
        data = np.asarray(field, dtype=float)
        if scale == "log":
            plot = np.log10(np.maximum(np.abs(data), 1.0e-300))
            valid = plot[region & np.isfinite(plot)]
            if valid.size:
                vmin = float(np.percentile(valid, 5.0))
                vmax = float(np.percentile(valid, 95.0))
                if vmax <= vmin:
                    vmax = vmin + 1.0
            else:
                vmin, vmax = -1.0, 1.0
        else:
            plot = data
            valid = plot[region & np.isfinite(plot)]
            vmin = 0.0
            vmax = max(float(np.percentile(valid, 95.0)) if valid.size else 1.0, 1.0e-300)
        im = ax.pcolormesh(xg, zg, np.clip(plot, vmin, vmax), shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.contour(xg, zg, region.astype(float), levels=[0.5], colors="white", linewidths=0.5)
        if np.any(active):
            ax.contour(xg, zg, active.astype(float), levels=[0.5], colors="cyan", linewidths=0.35)
        ax.set_title(title)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("z [um]")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_summary(out_path: Path, records: list[dict[str, object]]) -> None:
    labels = [f"tau={r['tau']:g}" for r in records]
    q_small = []
    chat = []
    eu2 = []
    eutr = []
    for r in records:
        branch = r["branch_bins"].get("qrel_le_0.1") or next(iter(r["branch_bins"].values()))
        q_small.append(int(branch["count"]))
        chat.append(float(branch["C_over_chi_norm"]["weighted_p95"] if "weighted_p95" in branch["C_over_chi_norm"] else branch["C_over_chi_norm"]["p95"]))
        patch = r["patch"]
        eu2.append(float(patch["Eu2_residual_stats"].get("weighted_p95", patch["Eu2_residual_stats"]["p95"])))
        eutr.append(float(patch["Eu_trace0_residual_stats"].get("weighted_p95", patch["Eu_trace0_residual_stats"]["p95"])))
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), constrained_layout=True)
    x = np.arange(len(labels))
    axes[0].bar(x, q_small)
    axes[0].set_title("q_rel <= 0.1 count")
    axes[0].set_xticks(x, labels, rotation=20, ha="right")
    axes[1].plot(x, chat, "o-")
    axes[1].set_yscale("log")
    axes[1].set_title("p95 |C|/chi in q_rel<=0.1")
    axes[1].set_xticks(x, labels, rotation=20, ha="right")
    axes[2].plot(x, eu2, "o-", label="span{g,uu}")
    axes[2].plot(x, eutr, "s-", label="trace0 only")
    axes[2].set_yscale("log")
    axes[2].set_title("p95 low-dim patch residual")
    axes[2].set_xticks(x, labels, rotation=20, ha="right")
    axes[2].legend()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict[str, object]:
    args.output.mkdir(parents=True, exist_ok=True)
    if args.mp is None:
        args.mp = PLANCK_MASS_EV
    ref = build_physical_reference(args)
    records = [analyze_one(args, ref, float(t.strip())) for t in args.taus.split(",") if t.strip()]
    summary_path = args.output / "branch_patch_summary.png"
    render_summary(summary_path, records)
    params = ref["params"]
    scale = ref["scale"]
    report = {
        "parameters": {
            "wavelength_nm": float(args.wavelength_nm),
            "mass_over_omega": float(args.mass_over_omega),
            "m_ev": float(params.m),
            "mp_ev": float(args.mp),
            "full_resolution": int(args.full_resolution),
            "window_um": float(args.window_um),
            "fit_region": str(args.fit_region),
            "force_region": str(args.force_region),
            "diagnostic_region": str(args.diagnostic_region),
            "force_probe_dt_old": float(args.force_probe_dt_old),
            "force_probe_dt_ev_inv": float(args.force_probe_dt_old) * float(scale.old_dimensionless_scale_ev_inv),
            "chi_definition": f"chi=(q_rel^2)/(q_rel^2+{args.chi_q0:g}^2), q_rel=|q|/p95(|q| in diagnostic region)",
            "q_definition": "q=(tilde nabla_mu r^mu + r_mu r^mu)/m^2 with r_mu=tilde nabla_mu ln sqrt(rho_tilde)",
            "Eu_patch_definition": "span{tilde g_mn,u_m u_n}; trace0 low patch uses u_m u_n-(u^2/3)tilde g_mn because this numerical test is 2+1d.",
        },
        "records": records,
        "outputs": {
            "summary_json": str((args.output / "summary.json").resolve()),
            "summary_png": str(summary_path.resolve()),
            "map_pngs": [r["map_png"] for r in records],
        },
    }
    (args.output / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--taus", type=str, default="-3.5,0,3.5")
    parser.add_argument("--fit-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--force-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--diagnostic-region", choices=["support", "trusted", "core10"], default="core10")
    parser.add_argument("--target-floor-frac", type=float, default=0.05)
    parser.add_argument("--full-resolution", type=int, default=384)
    parser.add_argument("--window-um", type=float, default=9.0)
    parser.add_argument("--dt-old", type=float, default=2.5e-7)
    parser.add_argument("--force-probe-dt-old", type=float, default=2.5e-7)
    parser.add_argument("--active-dilation", type=int, default=1)
    parser.add_argument("--coeff-root", type=Path, default=Path("visualizations"))
    parser.add_argument("--coeff-template", type=str, default="")
    parser.add_argument("--coeff-key", type=str, default="coeff_hard_alm")
    parser.add_argument("--chi-q0", type=float, default=0.1)
    parser.add_argument("--chi-floor", type=float, default=1.0e-12)
    parser.add_argument("--qrel-thresholds", type=str, default="0.01,0.03,0.1,0.3")
    parser.add_argument("--qabs-thresholds", type=str, default="0.01,0.1,1")
    parser.add_argument("--delta-thresholds", type=str, default="1e-4,1e-3,1e-2,1e-1")
    parser.add_argument("--rperp-thresholds", type=str, default="1e-3,1e-2,1e-1")
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
    run(parser.parse_args())


if __name__ == "__main__":
    main()
