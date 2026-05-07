from __future__ import annotations

import json
from pathlib import Path

from simulate_bc_small_gravity_backbone import run


def main():
    out = Path(__file__).resolve().parent.parent / "visualizations" / "c_dt_convergence_scan"
    out.mkdir(parents=True, exist_ok=True)

    cases = {
        "lambda_1": {
            "lambda_grav": 1.0,
            "dt_targets": [0.025, 0.0125, 0.00625, 0.003125],
        },
        "lambda_10": {
            "lambda_grav": 10.0,
            "dt_targets": [0.025, 0.0125, 0.00625],
        },
    }

    summary: dict[str, dict[str, object]] = {}

    for label, cfg in cases.items():
        case_summary: dict[str, object] = {}
        for dt_target in cfg["dt_targets"]:
            tag = f"dt_{str(dt_target).replace('.', 'p')}"
            case_out = out / f"{label}_{tag}"
            result = run(
                case_out,
                alpha=0.5,
                lambda_grav=cfg["lambda_grav"],
                mp=300.0,
                ell=0.02,
                dt_target=dt_target,
            )
            case_summary[str(dt_target)] = result
        summary[label] = case_summary

    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
