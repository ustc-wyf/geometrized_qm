from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from analyze_ab_measurable_observables import run as run_flat_observables

G_NEWTON = 6.67430e-11
C_LIGHT = 299792458.0

M_EARTH = 5.9722e24
R_EARTH = 6.371e6
M_SUN = 1.98847e30
AU = 1.495978707e11


def point_mass_curvature_scale(mass_kg: float, radius_m: float) -> float:
    return G_NEWTON * mass_kg / (C_LIGHT**2 * radius_m**3)


def orientation_rows_for_axisymmetric_background(
    kappa_total: float,
    readout_span_m: float,
    axis_label: str,
    flat: dict[str, float],
    numerical_floor: float | None,
) -> list[dict[str, float | str | None]]:
    rows = []
    for label, theta_deg in [
        ("axial", 0.0),
        ("tilted_45deg", 45.0),
        ("transverse", 90.0),
        ("magic_angle", math.degrees(math.acos(1.0 / math.sqrt(3.0)))),
    ]:
        theta = math.radians(theta_deg)
        quad = 3.0 * math.cos(theta) ** 2 - 1.0
        frac = float(kappa_total * (readout_span_m**2) * quad)
        delta_peak = float(flat["central_peak_rel_intensity_shift"] * frac)
        if numerical_floor is None or numerical_floor == 0.0:
            ratio = None
        else:
            ratio = float(abs(delta_peak) / numerical_floor)
        rows.append(
            {
                "orientation": label,
                "angle_to_axis_deg": theta_deg,
                "axis_reference": axis_label,
                "quadrupole_factor_signed": float(quad),
                "fractional_modulation_estimate": frac,
                "estimated_total_AB_central_peak_shift": float(flat["central_peak_rel_intensity_shift"] * (1.0 + frac)),
                "incremental_background_change_in_central_peak_shift": delta_peak,
                "increment_vs_numerical_floor": ratio,
            }
        )
    return rows


def gw_rows(
    h_gw: float,
    flat: dict[str, float],
    numerical_floor: float | None,
) -> list[dict[str, float | str | None]]:
    rows = []
    for label, theta_deg, frac in [
        ("plus_aligned", 0.0, 0.5 * h_gw),
        ("plus_45deg", 45.0, 0.5 * h_gw),
        ("plus_orthogonal", 90.0, -0.5 * h_gw),
        ("propagation_axis", None, 0.0),
    ]:
        delta_peak = float(flat["central_peak_rel_intensity_shift"] * frac)
        if numerical_floor is None or numerical_floor == 0.0:
            ratio = None
        else:
            ratio = float(abs(delta_peak) / numerical_floor)
        rows.append(
            {
                "orientation": label,
                "angle_to_plus_axis_deg": theta_deg,
                "fractional_modulation_estimate": float(frac),
                "estimated_total_AB_central_peak_shift": float(flat["central_peak_rel_intensity_shift"] * (1.0 + frac)),
                "incremental_background_change_in_central_peak_shift": delta_peak,
                "increment_vs_numerical_floor": ratio,
            }
        )
    return rows


def max_abs_fraction(rows: list[dict[str, object]]) -> float:
    return float(max(abs(float(row["fractional_modulation_estimate"])) for row in rows))


def max_abs_increment(rows: list[dict[str, object]]) -> float:
    return float(max(abs(float(row["incremental_background_change_in_central_peak_shift"])) for row in rows))


def load_flat_reference(root: Path) -> dict[str, object]:
    conv_path = root / "visualizations" / "ab_observable_convergence" / "summary.json"
    if conv_path.exists():
        conv = json.loads(conv_path.read_text(encoding="utf-8"))
        finest_dt = conv["dt_scan"][-1]["observables"]
        return {
            "centerline_observables": {
                "central_peak_rel_intensity_shift": float(finest_dt["central_peak_rel_intensity_shift"]),
                "mean_peak_rel_intensity_shift": float(finest_dt["mean_peak_rel_intensity_shift"]),
                "support_pointwise_rel_diff_max": float(finest_dt["support_pointwise_rel_diff_max"]),
                "effective_phase_shift_rad": float(finest_dt["effective_phase_shift_rad"]),
                "equivalent_fringe_shift_over_spacing": float(finest_dt["equivalent_fringe_shift_over_spacing"]),
            }
        }
    path = root / "visualizations" / "ab_measurable_observables_lambda1" / "summary.json"
    if not path.exists():
        run_flat_observables(root / "visualizations" / "ab_measurable_observables_lambda1")
    return json.loads(path.read_text(encoding="utf-8"))


def load_numerical_floor(root: Path) -> dict[str, float] | None:
    path = root / "visualizations" / "ab_observable_convergence" / "summary.json"
    if not path.exists():
        return None
    summary = json.loads(path.read_text(encoding="utf-8"))
    floors = summary.get("estimated_numerical_floors", {})
    dt_floor = floors.get("dt_scan", {})
    dx_floor = floors.get("grid_scan", {})
    merged = {}
    keys = set(dt_floor) | set(dx_floor)
    for key in keys:
        merged[key] = float(max(float(dt_floor.get(key, 0.0)), float(dx_floor.get(key, 0.0))))
    return merged


def render_bar(path: Path, labels: list[str], values: list[float], ylabel: str) -> None:
    plt.figure(figsize=(9.4, 4.8))
    plt.bar(labels, values, color="#3b6fb6")
    plt.yscale("log")
    plt.ylabel(ylabel)
    plt.xticks(rotation=30, ha="right")
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(
    outdir: str | Path,
    earth_readout_span_m: float = 1.0,
    l1_readout_span_m: float = 1.0,
    cavendish_readout_span_m: float = 0.01,
    gw_strain: float = 2.4e-15,
    cavendish_mass_kg: float = 1000.0,
    cavendish_standoff_m: float = 0.1,
) -> dict[str, object]:
    root = Path(__file__).resolve().parent.parent
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    flat_summary = load_flat_reference(root)
    flat = flat_summary["centerline_observables"]
    numerical_floor = load_numerical_floor(root)
    peak_floor = None if numerical_floor is None else numerical_floor.get("central_peak_rel_intensity_shift")

    earth_rows = orientation_rows_for_axisymmetric_background(
        point_mass_curvature_scale(M_EARTH, R_EARTH),
        earth_readout_span_m,
        "local radial to Earth center",
        flat,
        peak_floor,
    )
    leo_rows = orientation_rows_for_axisymmetric_background(
        point_mass_curvature_scale(M_EARTH, R_EARTH + 4.0e5),
        earth_readout_span_m,
        "local radial to Earth center",
        flat,
        peak_floor,
    )
    geo_rows = orientation_rows_for_axisymmetric_background(
        point_mass_curvature_scale(M_EARTH, 4.2164e7),
        earth_readout_span_m,
        "local radial to Earth center",
        flat,
        peak_floor,
    )
    l1_kappa = point_mass_curvature_scale(M_EARTH, 1.5e9) + point_mass_curvature_scale(M_SUN, AU - 1.5e9)
    l1_rows = orientation_rows_for_axisymmetric_background(
        l1_kappa,
        l1_readout_span_m,
        "Sun-Earth line",
        flat,
        peak_floor,
    )
    cav_rows = orientation_rows_for_axisymmetric_background(
        point_mass_curvature_scale(cavendish_mass_kg, cavendish_standoff_m),
        cavendish_readout_span_m,
        "line to nearby source mass",
        flat,
        peak_floor,
    )
    gw_background_rows = gw_rows(gw_strain, flat, peak_floor)

    scenarios = {
        "earth_surface": {
            "type": "vacuum_tidal",
            "curvature_scale_m_minus_2": point_mass_curvature_scale(M_EARTH, R_EARTH),
            "readout_span_m": earth_readout_span_m,
            "rows": earth_rows,
        },
        "leo_400km": {
            "type": "vacuum_tidal",
            "curvature_scale_m_minus_2": point_mass_curvature_scale(M_EARTH, R_EARTH + 4.0e5),
            "readout_span_m": earth_readout_span_m,
            "rows": leo_rows,
        },
        "geo": {
            "type": "vacuum_tidal",
            "curvature_scale_m_minus_2": point_mass_curvature_scale(M_EARTH, 4.2164e7),
            "readout_span_m": earth_readout_span_m,
            "rows": geo_rows,
        },
        "sun_earth_l1": {
            "type": "vacuum_tidal",
            "curvature_scale_m_minus_2": l1_kappa,
            "readout_span_m": l1_readout_span_m,
            "rows": l1_rows,
        },
        "cavendish_1tonne_10cm": {
            "type": "vacuum_tidal",
            "curvature_scale_m_minus_2": point_mass_curvature_scale(cavendish_mass_kg, cavendish_standoff_m),
            "readout_span_m": cavendish_readout_span_m,
            "rows": cav_rows,
        },
        "nanohertz_gw_background": {
            "type": "vacuum_wave",
            "characteristic_strain": gw_strain,
            "rows": gw_background_rows,
        },
    }

    labels = []
    mod_values = []
    inc_values = []
    for name, info in scenarios.items():
        labels.append(name)
        mod_values.append(max_abs_fraction(info["rows"]))
        inc_values.append(max_abs_increment(info["rows"]))

    render_bar(out / "scenario_max_fractional_modulation.png", labels, mod_values, "max |fractional modulation|")
    render_bar(
        out / "scenario_central_peak_increment.png",
        labels,
        inc_values,
        "max |background-induced change in central peak shift|",
    )

    summary = {
        "flat_reference": {
            "central_peak_rel_intensity_shift": float(flat["central_peak_rel_intensity_shift"]),
            "mean_peak_rel_intensity_shift": float(flat["mean_peak_rel_intensity_shift"]),
            "support_pointwise_rel_diff_max": float(flat["support_pointwise_rel_diff_max"]),
            "effective_phase_shift_rad": float(flat["effective_phase_shift_rad"]),
            "equivalent_fringe_shift_over_spacing": float(flat["equivalent_fringe_shift_over_spacing"]),
        },
        "numerical_floor_reference": {
            "central_peak_rel_intensity_shift": peak_floor,
        },
        "scenarios": scenarios,
        "notes": {
            "scope": "This script upgrades the previous single-point weak-field estimate into a scenario ladder with explicit orientation scans and a direct comparison to the current numerical floor of the main AB intensity signal.",
            "vacuum_caveat": "Earth, LEO, GEO, Sun-Earth L1, nearby source mass, and a plane GW are all vacuum backgrounds in the beam region, so these entries remain tidal/readout modulations rather than first-order bulk Einstein-Hilbert differences.",
            "new_candidate": "The modulated nearby source mass is included because it is experimentally switchable and can exceed Earth-surface tidal curvature in a small apparatus even though it remains a vacuum tidal test in the current approximation.",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = run(Path(__file__).resolve().parent.parent / "visualizations" / "ab_background_scenarios_lambda1")
    print(json.dumps(result, indent=2))
