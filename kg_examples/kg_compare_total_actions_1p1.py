from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

sys.path.append(str(Path(__file__).resolve().parent))

from kg_double_packet_curvature import compute_metric_on_grid, compute_ricci_scalar_2d


def run_comparison(output_dir: str | Path) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    x = np.linspace(-12.0, 12.0, 241)
    t0 = 9.583148474999101
    t = np.linspace(t0 - 1.0, t0 + 1.0, 41)

    geom = compute_metric_on_grid(x, t)
    R_tilde = compute_ricci_scalar_2d(geom["g_tilde_cov"], x, t)

    det = geom["det_g_tilde"]
    sqrt_minus_g_tilde = np.sqrt(np.maximum(-det, 0.0))
    L_eh_tilde = sqrt_minus_g_tilde * R_tilde

    rho = geom["rho"]
    eps = geom["eps"]
    mask_support = rho >= rho.max() * 1e-2
    mask_loose = rho >= rho.max() * 1e-6

    dx = float(x[1] - x[0])
    dt = float(t[1] - t[0])

    result = {
        "statement": {
            "eh_g_bulk_density": 0.0,
            "eh_gtilde_bulk_density_nonzero": True,
            "bulk_prediction_difference_for_pure_1p1_EH": "absent_in_bulk_topological_limit",
        },
        "summary": {
            "R_tilde_abs_max_all": float(np.nanmax(np.abs(R_tilde))),
            "R_tilde_abs_max_support_1e-2": float(np.nanmax(np.abs(R_tilde[mask_support]))),
            "R_tilde_abs_max_support_1e-6": float(np.nanmax(np.abs(R_tilde[mask_loose]))),
            "L_eh_tilde_abs_max_all": float(np.nanmax(np.abs(L_eh_tilde))),
            "L_eh_tilde_abs_max_support_1e-2": float(np.nanmax(np.abs(L_eh_tilde[mask_support]))),
            "L_eh_tilde_abs_max_support_1e-6": float(np.nanmax(np.abs(L_eh_tilde[mask_loose]))),
            "eh_tilde_box_integral": float(np.nansum(L_eh_tilde) * dx * dt),
            "eps_abs_max_support_1e-2": float(np.nanmax(np.abs(eps[mask_support]))),
            "eps_abs_max_support_1e-6": float(np.nanmax(np.abs(eps[mask_loose]))),
            "det_nonnegative_count": int(np.sum(det >= 0.0)),
        },
    }

    (out / "compare_total_actions_1p1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    np.savez_compressed(
        out / "compare_total_actions_1p1_fields.npz",
        x=x,
        t=t,
        rho=rho,
        eps=eps,
        R_tilde=R_tilde,
        L_eh_tilde=L_eh_tilde,
    )
    return result


if __name__ == "__main__":
    result = run_comparison(Path(__file__).resolve().parent / "outputs")
    print(json.dumps(result, ensure_ascii=False, indent=2))
