from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent))

from kg_double_slit_beamlike_bohm import compute_bohm_slice


def gradient_magnitude(arr: np.ndarray, dx: float, dz: float) -> np.ndarray:
    ax = np.gradient(arr, dx, axis=0)
    az = np.gradient(arr, dz, axis=1)
    return np.sqrt(ax**2 + az**2)


def render_maps(alpha: float, kappa: float, out_dir: Path) -> dict[str, str]:
    result = compute_bohm_slice(alpha)
    x = result["fields"]["x"]
    z = result["fields"]["z"]
    rho = np.real(result["fields"]["rho"])
    Q = np.real(result["fields"]["Q"])
    X = np.real(result["fields"]["X"])
    Xg, Zg = np.meshgrid(x, z, indexing="ij")

    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])

    eps = Q / result["params"]["m"] ** 2
    F = 1.0 / (1.0 + (kappa * Q) ** 2)
    gradF = gradient_magnitude(F, dx, dz)

    plots = {
        "rho": (rho, "viridis"),
        "Q_clipped": (np.clip(Q, -50, 50), "coolwarm"),
        "X_clipped": (np.clip(X, -50, 50), "coolwarm"),
        "eps_clipped": (np.clip(eps, -50, 50), "coolwarm"),
        "F": (F, "viridis"),
        "gradF": (gradF, "magma"),
    }

    out_paths: dict[str, str] = {}
    for name, (arr, cmap) in plots.items():
        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
        im = ax.pcolormesh(Xg, Zg, arr, shading="auto", cmap=cmap)
        ax.set_title(f"{name}(x,z), alpha={alpha}, kappa={kappa}")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, label=name)
        fp = out_dir / f"alpha_{alpha:.1f}_{name}.png"
        fig.savefig(fp, dpi=180)
        plt.close(fig)
        out_paths[name] = str(fp)

    return out_paths


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out = root / "visualizations" / "action_entry_maps"
    out.mkdir(parents=True, exist_ok=True)

    kappa = 0.5
    summary = {}
    for alpha in (0.0, 0.5, 1.0):
        summary[str(alpha)] = render_maps(alpha, kappa, out)

    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
