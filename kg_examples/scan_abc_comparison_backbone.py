from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

from simulate_bc_small_gravity_backbone import run


def choose_c_substeps(lambda_grav: float) -> int:
    if lambda_grav <= 1.0e-1:
        return 1
    if lambda_grav <= 3.0e-1:
        return 2
    if lambda_grav <= 1.0:
        return 4
    return 8


def first_crossing(results: list[dict[str, float]], key: str, threshold: float) -> float | None:
    for row in results:
        if row["lambda_grav"] <= 0.0:
            continue
        if row[key] >= threshold:
            return float(row["lambda_grav"])
    return None


def main(rerun: bool = False):
    out = Path(__file__).resolve().parent.parent / "visualizations" / "abc_comparison_backbone_scan"
    out.mkdir(parents=True, exist_ok=True)

    # Denser log-spaced continuation scan, with more C subcycling only once
    # the shared-baseline runs enter the known stiff regime.
    lambda_grid = [
        0.0,
        1.0e-4,
        3.0e-4,
        1.0e-3,
        3.0e-3,
        1.0e-2,
        3.0e-2,
        1.0e-1,
        3.0e-1,
        1.0,
        3.0,
        10.0,
    ]
    grid = [{"lambda_grav": lam, "c_substeps": choose_c_substeps(lam)} for lam in lambda_grid]

    results = []
    for cfg in grid:
        lam = cfg["lambda_grav"]
        tag = str(lam).replace(".", "p").replace("-", "m")
        case_out = out / f"lambda_{tag}"
        summary_path = case_out / "summary.json"
        if summary_path.exists() and not rerun:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        else:
            summary = run(
                case_out,
                alpha=0.5,
                lambda_grav=lam,
                mp=300.0,
                ell=0.02,
                dt_target=0.025,
                c_substeps=cfg["c_substeps"],
            )
        summary["c_substeps"] = cfg["c_substeps"]
        results.append(summary)

    results.sort(key=lambda row: row["lambda_grav"])

    (out / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")

    thresholds = {
        "heuristic_first_crossings": {
            "B_vs_A_relL1_ge_1e-6": first_crossing(results, "B_vs_A_relL1", 1.0e-6),
            "B_vs_A_relL1_ge_1e-5": first_crossing(results, "B_vs_A_relL1", 1.0e-5),
            "C_vs_B_relL1_ge_1e-5": first_crossing(results, "C_vs_B_relL1", 1.0e-5),
            "C_vs_B_relL1_ge_1e-4": first_crossing(results, "C_vs_B_relL1", 1.0e-4),
            "C_vs_A_relL1_ge_1e-4": first_crossing(results, "C_vs_A_relL1", 1.0e-4),
        }
    }
    (out / "thresholds.json").write_text(json.dumps(thresholds, indent=2), encoding="utf-8")

    xs = [r["lambda_grav"] for r in results]
    b = [r["B_vs_A_relL1"] for r in results]
    c = [r["C_vs_A_relL1"] for r in results]
    c_b = [r["C_vs_B_relL1"] for r in results]

    plt.figure(figsize=(6.2, 4.2))
    plt.loglog(xs[1:], b[1:], marker="o", label="B vs A")
    plt.loglog(xs[1:], c[1:], marker="o", label="C vs A")
    plt.xlabel("lambda_grav")
    plt.ylabel("relative L1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "BA_CA_vs_lambda.png", dpi=180)
    plt.close()

    plt.figure(figsize=(6.2, 4.2))
    plt.loglog(xs[1:], c_b[1:], marker="o", label="C vs B")
    plt.xlabel("lambda_grav")
    plt.ylabel("relative L1")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out / "CB_vs_lambda.png", dpi=180)
    plt.close()

    print(
        json.dumps(
            {
                "results": results,
                **thresholds,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--rerun", action="store_true", help="Recompute cases even if their summary.json already exists.")
    args = parser.parse_args()
    main(rerun=args.rerun)
