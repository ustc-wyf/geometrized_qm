from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import imageio.v2 as imageio
except Exception:  # pragma: no cover
    imageio = None


@dataclass
class FlatLocalizedCrossingParams:
    m: float = 1.0
    alpha: float = 0.5
    k0: float = 6.0
    sigma_parallel: float = 1.6
    sigma_perp: float = 1.2
    phi0: float = 0.0
    x_a0: float = -6.0
    z_a0: float = -6.0
    x_b0: float = 6.0
    z_b0: float = -6.0
    x_half_range: float = 20.0
    z_half_range: float = 20.0
    nx: int = 256
    nz: int = 256
    t_final: float = 16.0
    n_frames: int = 7
    rho_floor: float = 1e-12

    @property
    def omega0(self) -> float:
        return float(np.sqrt(self.k0**2 + self.m**2))

    @property
    def v_group(self) -> float:
        return float(self.k0 / self.omega0)

    @property
    def t_meet(self) -> float:
        vz = self.v_group / np.sqrt(2.0)
        return float(-self.z_a0 / max(vz, 1e-12))


def make_grid(params: FlatLocalizedCrossingParams) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x = np.linspace(-params.x_half_range, params.x_half_range, params.nx, endpoint=False)
    z = np.linspace(-params.z_half_range, params.z_half_range, params.nz, endpoint=False)
    X, Z = np.meshgrid(x, z, indexing="ij")
    return x, z, X, Z


def packet_a_local_coords(X: np.ndarray, Z: np.ndarray, params: FlatLocalizedCrossingParams) -> tuple[np.ndarray, np.ndarray]:
    dx = X - params.x_a0
    dz = Z - params.z_a0
    q = (dx + dz) / np.sqrt(2.0)
    p = (dx - dz) / np.sqrt(2.0)
    return q, p


def packet_b_local_coords(X: np.ndarray, Z: np.ndarray, params: FlatLocalizedCrossingParams) -> tuple[np.ndarray, np.ndarray]:
    dx = X - params.x_b0
    dz = Z - params.z_b0
    q = (-dx + dz) / np.sqrt(2.0)
    p = (dx + dz) / np.sqrt(2.0)
    return q, p


def initial_wavefunction(X: np.ndarray, Z: np.ndarray, params: FlatLocalizedCrossingParams) -> np.ndarray:
    q_a, p_a = packet_a_local_coords(X, Z, params)
    q_b, p_b = packet_b_local_coords(X, Z, params)

    env_a = np.exp(
        -0.5 * ((q_a / params.sigma_parallel) ** 2 + (p_a / params.sigma_perp) ** 2)
    )
    env_b = np.exp(
        -0.5 * ((q_b / params.sigma_parallel) ** 2 + (p_b / params.sigma_perp) ** 2)
    )

    psi_a = params.alpha * env_a * np.exp(1j * params.k0 * q_a)
    psi_b = np.sqrt(max(0.0, 1.0 - params.alpha)) * env_b * np.exp(1j * (params.k0 * q_b + params.phi0))
    return psi_a + psi_b


def spectral_omega(x: np.ndarray, z: np.ndarray, m: float) -> np.ndarray:
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    kx = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(len(z), d=dz)
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    return np.sqrt(KX**2 + KZ**2 + m**2)


def evolve_positive_frequency(psi0: np.ndarray, omega: np.ndarray, t: float) -> np.ndarray:
    psi0_hat = np.fft.fft2(psi0)
    return np.fft.ifft2(psi0_hat * np.exp(-1j * omega * t))


def boundary_ratio(rho: np.ndarray) -> float:
    boundary = np.concatenate([rho[0, :], rho[-1, :], rho[:, 0], rho[:, -1]])
    return float(np.max(boundary) / max(float(np.max(rho)), 1e-30))


def packet_centers(params: FlatLocalizedCrossingParams, t: float) -> dict[str, tuple[float, float]]:
    v = params.v_group / np.sqrt(2.0)
    return {
        "A": (params.x_a0 + v * t, params.z_a0 + v * t),
        "B": (params.x_b0 - v * t, params.z_b0 + v * t),
    }


def phase_field(psi: np.ndarray) -> np.ndarray:
    return np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)


def render_montage(
    x: np.ndarray,
    z: np.ndarray,
    rho_frames: list[np.ndarray],
    times: np.ndarray,
    centers: list[dict[str, tuple[float, float]]],
    out_path: Path,
    log_scale: bool = False,
) -> None:
    n = len(rho_frames)
    fig, axes = plt.subplots(2, (n + 1) // 2, figsize=(3.8 * ((n + 1) // 2), 6.8), constrained_layout=True)
    axes = np.array(axes).reshape(-1)

    if log_scale:
        fields = [np.log10(np.maximum(rho, 1e-16)) for rho in rho_frames]
        cmap = "magma"
        cbar_label = r"$\log_{10}\rho$"
    else:
        fields = rho_frames
        cmap = "viridis"
        cbar_label = r"$\rho$"

    vmin = min(float(np.min(f)) for f in fields)
    vmax = max(float(np.max(f)) for f in fields)

    Xg, Zg = np.meshgrid(x, z, indexing="ij")
    last_im = None
    for ax, field, t, center in zip(axes, fields, times, centers):
        last_im = ax.pcolormesh(Xg, Zg, field, shading="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(fr"$t={t:.2f}$")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.plot(center["A"][0], center["A"][1], "wo", ms=3.5, mec="k", mew=0.4)
        ax.plot(center["B"][0], center["B"][1], "w^", ms=4.0, mec="k", mew=0.4)
        ax.set_aspect("equal")
    for ax in axes[len(rho_frames) :]:
        ax.axis("off")
    if last_im is not None:
        fig.colorbar(last_im, ax=axes.tolist(), shrink=0.92, label=cbar_label)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def render_boundary_plot(times: np.ndarray, ratios: np.ndarray, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.2), constrained_layout=True)
    ax.plot(times, ratios, marker="o", lw=1.8)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\rho_{\partial\Omega,\max}/\rho_{\max}$")
    ax.set_title("Boundary density ratio")
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)


def maybe_render_gif(frame_paths: list[Path], out_path: Path, duration_ms: int = 220) -> str | None:
    if imageio is None:
        return None
    images = [imageio.imread(path) for path in frame_paths]
    imageio.mimsave(out_path, images, duration=duration_ms / 1000.0, loop=0)
    return str(out_path)


def run_case(output_dir: str | Path, params: FlatLocalizedCrossingParams | None = None) -> dict[str, object]:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    frames_dir = out / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    params = params or FlatLocalizedCrossingParams()
    x, z, X, Z = make_grid(params)
    psi0 = initial_wavefunction(X, Z, params)
    omega = spectral_omega(x, z, params.m)

    times = np.linspace(0.0, params.t_final, params.n_frames)
    rho_frames: list[np.ndarray] = []
    S_frames: list[np.ndarray] = []
    centers: list[dict[str, tuple[float, float]]] = []
    boundary_ratios: list[float] = []
    frame_paths: list[Path] = []

    for idx, t in enumerate(times):
        psi = evolve_positive_frequency(psi0, omega, float(t))
        rho = np.abs(psi) ** 2
        S = phase_field(psi)
        rho_frames.append(rho)
        S_frames.append(S)
        centers_now = packet_centers(params, float(t))
        centers.append(centers_now)
        boundary_ratios.append(boundary_ratio(rho))

        fig, ax = plt.subplots(figsize=(5.2, 4.8), constrained_layout=True)
        Xg, Zg = np.meshgrid(x, z, indexing="ij")
        im = ax.pcolormesh(Xg, Zg, rho, shading="auto", cmap="viridis")
        ax.plot(centers_now["A"][0], centers_now["A"][1], "wo", ms=4, mec="k", mew=0.4)
        ax.plot(centers_now["B"][0], centers_now["B"][1], "w^", ms=4.5, mec="k", mew=0.4)
        ax.set_title(fr"Localized packet density, $t={t:.2f}$")
        ax.set_xlabel("x")
        ax.set_ylabel("z")
        ax.set_aspect("equal")
        fig.colorbar(im, ax=ax, label=r"$\rho$")
        frame_path = frames_dir / f"rho_{idx:03d}.png"
        fig.savefig(frame_path, dpi=180)
        plt.close(fig)
        frame_paths.append(frame_path)

    render_montage(x, z, rho_frames, times, centers, out / "rho_montage.png", log_scale=False)
    render_montage(x, z, rho_frames, times, centers, out / "rho_log10_montage.png", log_scale=True)
    render_boundary_plot(times, np.asarray(boundary_ratios), out / "boundary_ratio_vs_time.png")
    gif_path = maybe_render_gif(frame_paths, out / "rho_evolution.gif")

    np.savez_compressed(
        out / "flat_localized_crossing_packets_fields.npz",
        x=x,
        z=z,
        times=times,
        rho=np.stack(rho_frames, axis=0),
        S=np.stack(S_frames, axis=0),
        psi0=psi0,
    )

    summary = {
        "params": {
            **asdict(params),
            "omega0": params.omega0,
            "v_group": params.v_group,
            "t_meet": params.t_meet,
        },
        "summary": {
            "max_boundary_ratio_over_time": float(np.max(boundary_ratios)),
            "boundary_ratios": [float(v) for v in boundary_ratios],
            "times": [float(t) for t in times],
        },
        "files": {
            "rho_montage": str((out / "rho_montage.png").resolve()),
            "rho_log10_montage": str((out / "rho_log10_montage.png").resolve()),
            "boundary_ratio_vs_time": str((out / "boundary_ratio_vs_time.png").resolve()),
            "fields_npz": str((out / "flat_localized_crossing_packets_fields.npz").resolve()),
            "gif": gif_path,
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    out = root / "visualizations" / "flat_localized_crossing_packets"
    result = run_case(out)
    print(json.dumps(result, ensure_ascii=False, indent=2))
