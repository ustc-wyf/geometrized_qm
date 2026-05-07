import json
import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parent))

from kg_double_slit_beamlike import compute_xz_slice, default_params


def main():
    params = default_params()
    params.alpha = 0.5
    params.k0 = 10.0
    params.w0 = 1.0
    params.x_half_range = 5.0
    params.z_half_range = 5.0
    params.nx = 401
    params.nz = 401

    data = compute_xz_slice(0.0, params)
    psi = data["psi"]
    x = data["x"]
    z = data["z"]
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    omega = params.omega0
    m = params.m

    psix = np.gradient(psi, dx, axis=0)
    psiz = np.gradient(psi, dz, axis=1)
    psixx = np.gradient(psix, dx, axis=0)
    psizz = np.gradient(psiz, dz, axis=1)

    residual = psixx + psizz + (omega**2 - m**2) * psi
    mask = np.abs(psi) > 1e-4 * np.abs(psi).max()
    rel = np.abs(residual[mask]) / np.maximum(np.abs(((omega**2 - m**2) * psi)[mask]), 1e-12)

    summary = {
        "omega": float(omega),
        "abs_psi_max": float(np.abs(psi).max()),
        "abs_residual_max_main": float(np.abs(residual[mask]).max()),
        "abs_residual_mean_main": float(np.abs(residual[mask]).mean()),
        "relative_residual_max_main": float(rel.max()),
        "relative_residual_mean_main": float(rel.mean()),
    }

    outdir = Path(__file__).resolve().parent / "outputs"
    outdir.mkdir(exist_ok=True)
    (outdir / "benchmark_residual_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
