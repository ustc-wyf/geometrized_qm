from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def raw_required_counts(signal: float, sigma_level: float) -> float:
    return (sigma_level / abs(signal)) ** 2


def phase_required_counts(signal: float, d_ref: float, sigma_level: float) -> float:
    coeff = max(1.0 - d_ref * d_ref, 1.0e-30)
    return sigma_level * sigma_level * coeff / (abs(signal) ** 2)


def interference_required_counts(signal: float, e_ref: float, sigma_level: float) -> tuple[float, float]:
    coeff = 2.0 + 3.0 * e_ref + e_ref * e_ref
    denom_counts = sigma_level * sigma_level * coeff / (abs(signal) ** 2)
    total_counts = denom_counts * (2.0 + e_ref)
    return denom_counts, total_counts


def phase_pair_imbalance_tolerance(signal: float, d_ref: float, budget_fraction: float = 1.0) -> float:
    coeff = max(1.0 - d_ref * d_ref, 1.0e-30)
    return budget_fraction * abs(signal) / coeff


def interference_norm_tolerance(signal: float, e_ref: float, budget_fraction: float = 1.0) -> float:
    coeff = max(abs(1.0 + e_ref), 1.0e-30)
    return budget_fraction * abs(signal) / coeff


def render_required_counts(path: Path, channel_names: list[str], count_sets: dict[str, list[float]], ylabel: str) -> None:
    x = np.arange(len(channel_names))
    width = 0.22
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    for i, (label, values) in enumerate(count_sets.items()):
        ax.bar(x + (i - 1) * width, values, width=width, label=label)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(channel_names, rotation=18, ha="right")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def render_snr_curves(path: Path, count_grid: np.ndarray, curves: list[tuple[str, np.ndarray]]) -> None:
    plt.figure(figsize=(7.2, 4.8))
    for label, values in curves:
        plt.semilogx(count_grid, values, label=label, linewidth=1.7)
    plt.xlabel("counts")
    plt.ylabel("shot-noise SNR")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def run(outdir: str | Path) -> dict[str, object]:
    root = Path(__file__).resolve().parent.parent
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    raw = load_json(root / "visualizations" / "ab_time_observables_lambda1" / "summary.json")
    diff = load_json(root / "visualizations" / "ab_differential_experiment_observables_lambda1" / "summary.json")
    diff_series = np.load(root / "visualizations" / "ab_differential_experiment_observables_lambda1" / "series.npz")

    raw_peak_signal = float(raw["final_values"]["center_point_rel_shift"])
    raw_box_signal = float(raw["final_values"]["interference_box_rel_shift"])

    phase_box = diff["phase_flip_differential"]["interference_box"]
    phase_box_signal = abs(float(phase_box["max_abs_split"]))
    phase_box_ref = float(diff_series["phase_flip_A_box"][int(np.argmax(np.abs(diff_series["phase_flip_B_box"] - diff_series["phase_flip_A_box"])))])

    eint_box = diff["interference_excess"]["interference_box"]
    eint_box_signal = abs(float(eint_box["max_abs_split"]))
    eint_box_ref = float(diff_series["interference_A_box"][int(np.argmax(np.abs(diff_series["interference_B_box"] - diff_series["interference_A_box"])))])

    channels = {
        "raw_center_peak": {
            "type": "raw_fractional",
            "signal": raw_peak_signal,
            "description": "raw intensity shift at center peak",
        },
        "raw_interference_box": {
            "type": "raw_fractional",
            "signal": raw_box_signal,
            "description": "raw integrated intensity shift in interference box",
        },
        "phase_flip_box": {
            "type": "phase_flip",
            "signal": phase_box_signal,
            "reference_observable_A": phase_box_ref,
            "description": "phase-flip differential in interference box",
        },
        "interference_excess_box": {
            "type": "interference_excess",
            "signal": eint_box_signal,
            "reference_observable_A": eint_box_ref,
            "description": "interference excess differential in interference box",
        },
    }

    sigma_levels = [1.0, 3.0, 5.0]
    count_requirements: dict[str, dict[str, float]] = {}
    tolerances: dict[str, dict[str, float]] = {}
    snr_curves: list[tuple[str, np.ndarray]] = []
    count_grid = np.logspace(6, 14, 300)

    for name, info in channels.items():
        signal = float(info["signal"])
        entry_counts: dict[str, float] = {}
        entry_tol: dict[str, float] = {
            "signal_size_budget": abs(signal),
            "ten_percent_budget": 0.1 * abs(signal),
        }

        if info["type"] == "raw_fractional":
            for level in sigma_levels:
                entry_counts[f"{int(level)}sigma_window_counts"] = raw_required_counts(signal, level)
            snr_curves.append((name, abs(signal) * np.sqrt(count_grid)))

        elif info["type"] == "phase_flip":
            d_ref = float(info["reference_observable_A"])
            for level in sigma_levels:
                entry_counts[f"{int(level)}sigma_total_counts_two_phase_states"] = phase_required_counts(signal, d_ref, level)
            entry_tol["pair_imbalance_signal_size_budget"] = phase_pair_imbalance_tolerance(signal, d_ref, 1.0)
            entry_tol["pair_imbalance_ten_percent_budget"] = phase_pair_imbalance_tolerance(signal, d_ref, 0.1)
            snr_curves.append((name, abs(signal) * np.sqrt(count_grid / max(1.0 - d_ref * d_ref, 1.0e-30))))

        elif info["type"] == "interference_excess":
            e_ref = float(info["reference_observable_A"])
            for level in sigma_levels:
                denom_counts, total_counts = interference_required_counts(signal, e_ref, level)
                entry_counts[f"{int(level)}sigma_denominator_counts_single_beam_sum"] = denom_counts
                entry_counts[f"{int(level)}sigma_total_counts_all_three_runs"] = total_counts
            entry_tol["single_beam_normalization_signal_size_budget"] = interference_norm_tolerance(signal, e_ref, 1.0)
            entry_tol["single_beam_normalization_ten_percent_budget"] = interference_norm_tolerance(signal, e_ref, 0.1)
            coeff = 2.0 + 3.0 * e_ref + e_ref * e_ref
            snr_curves.append((name, abs(signal) * np.sqrt(count_grid / coeff)))

        count_requirements[name] = entry_counts
        tolerances[name] = entry_tol

    count_sets_5sigma = {
        "5sigma": [
            count_requirements["raw_center_peak"]["5sigma_window_counts"],
            count_requirements["raw_interference_box"]["5sigma_window_counts"],
            count_requirements["phase_flip_box"]["5sigma_total_counts_two_phase_states"],
            count_requirements["interference_excess_box"]["5sigma_total_counts_all_three_runs"],
        ],
        "3sigma": [
            count_requirements["raw_center_peak"]["3sigma_window_counts"],
            count_requirements["raw_interference_box"]["3sigma_window_counts"],
            count_requirements["phase_flip_box"]["3sigma_total_counts_two_phase_states"],
            count_requirements["interference_excess_box"]["3sigma_total_counts_all_three_runs"],
        ],
        "1sigma": [
            count_requirements["raw_center_peak"]["1sigma_window_counts"],
            count_requirements["raw_interference_box"]["1sigma_window_counts"],
            count_requirements["phase_flip_box"]["1sigma_total_counts_two_phase_states"],
            count_requirements["interference_excess_box"]["1sigma_total_counts_all_three_runs"],
        ],
    }

    render_required_counts(
        out / "required_counts_by_channel.png",
        ["raw peak", "raw box", "phase flip", "interference excess"],
        count_sets_5sigma,
        "required counts",
    )
    render_snr_curves(out / "shot_noise_snr_vs_counts.png", count_grid, snr_curves)

    summary = {
        "channels": channels,
        "count_requirements": count_requirements,
        "systematic_tolerances": tolerances,
        "notes": {
            "raw_fractional_model": "For a raw relative intensity shift delta, Poisson shot noise gives sigma_rel ~ 1/sqrt(N_window).",
            "phase_flip_model": "For D_phi=(N0-Npi)/(N0+Npi), Poisson shot noise gives Var(D_phi) ~ (1-D_phi^2)/N_total across the two phase states.",
            "interference_excess_model": "For E_int=(N_both-(N1+N2))/(N1+N2), Poisson shot noise gives Var(E_int) ~ (2+3E+E^2)/B where B=(N1+N2) is the denominator count budget.",
            "systematics": "signal_size_budget means the systematic alone equals the predicted AB split; ten_percent_budget means it consumes only 10% of that split.",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    run(Path(__file__).resolve().parent.parent / "visualizations" / "ab_noise_baselines_lambda1")
