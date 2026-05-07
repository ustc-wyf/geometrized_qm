from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from kg_crossing_packets_3d import run_crossing_case


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out = root / "visualizations" / "crossing3d_case"
    out.mkdir(parents=True, exist_ok=True)

    # ensure latest data
    result = run_crossing_case(root / "kg_examples" / "outputs")

    data = np.load(root / "kg_examples" / "outputs" / "crossing3d_fields.npz")
    x = data["x"]
    z = data["z"]
    rho = data["rho"]

    Xg, Zg = np.meshgrid(x, z, indexing="ij")

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    im = ax.pcolormesh(Xg, Zg, rho, shading="auto", cmap="viridis")
    ax.set_title("3+1d crossing packets: density rho(x,z) on y=0 slice")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label="rho")
    fig.savefig(out / "rho_linear.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    im = ax.pcolormesh(Xg, Zg, np.log10(np.maximum(rho, 1e-16)), shading="auto", cmap="magma")
    ax.set_title("3+1d crossing packets: log10 rho(x,z) on y=0 slice")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label="log10 rho")
    fig.savefig(out / "rho_log10.png", dpi=180)
    plt.close(fig)

    # boundary diagnostics
    top = rho[:, -1]
    bottom = rho[:, 0]
    left = rho[0, :]
    right = rho[-1, :]

    fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ax[0, 0].plot(x, top)
    ax[0, 0].set_title("top boundary z=z_max")
    ax[0, 1].plot(x, bottom)
    ax[0, 1].set_title("bottom boundary z=z_min")
    ax[1, 0].plot(z, left)
    ax[1, 0].set_title("left boundary x=x_min")
    ax[1, 1].plot(z, right)
    ax[1, 1].set_title("right boundary x=x_max")
    for row in ax:
        for a in row:
            a.grid(True, alpha=0.3)
    fig.savefig(out / "boundary_cuts.png", dpi=180)
    plt.close(fig)

    summary = {
        **result,
        "files": {
            "rho_linear": str((out / "rho_linear.png").resolve()),
            "rho_log10": str((out / "rho_log10.png").resolve()),
            "boundary_cuts": str((out / "boundary_cuts.png").resolve()),
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
