from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt

from simulate_bc_psfd_imex import run


def main():
    outdir = Path(__file__).resolve().parent.parent / "visualizations" / "bc_weak_coupling_scan"
    outdir.mkdir(parents=True, exist_ok=True)

    cases = [
        {"name": "case_1", "mp": 300.0, "ell": 0.02},
        {"name": "case_2", "mp": 1000.0, "ell": 0.01},
        {"name": "case_3", "mp": 3000.0, "ell": 0.005},
        {"name": "case_4", "mp": 10000.0, "ell": 0.001},
    ]

    summaries = []
    for case in cases:
        case_dir = outdir / case["name"]
        summary = run(
            case_dir,
            nx=121,
            nz=121,
            nkx=41,
            nkz=41,
            steps_hint=220,
            mp=case["mp"],
            ell=case["ell"],
        )
        summary["case"] = case
        summaries.append(summary)

    # Plot trend lines
    mp_values = [s["case"]["mp"] for s in summaries]
    ell_values = [s["case"]["ell"] for s in summaries]
    b_err = [s["errors_inner"]["B_vs_A_relL1"] for s in summaries]
    c_err = [s["errors_inner"]["C_vs_A_relL1"] for s in summaries]

    plt.figure(figsize=(6.4, 4.2))
    plt.plot(mp_values, b_err, marker="o", label="B vs A")
    plt.plot(mp_values, c_err, marker="o", label="C vs A")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("M_P")
    plt.ylabel("relative L1 error")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "error_vs_mp.png", dpi=180)
    plt.close()

    plt.figure(figsize=(6.4, 4.2))
    plt.plot(ell_values, c_err, marker="o", label="C vs A")
    plt.gca().invert_xaxis()
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("ell")
    plt.ylabel("relative L1 error")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "c_error_vs_ell.png", dpi=180)
    plt.close()

    (outdir / "summary.json").write_text(json.dumps(summaries, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
