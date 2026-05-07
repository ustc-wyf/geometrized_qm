import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent))

from kg_double_slit_beamlike import default_params, compute_xz_slice


def ddx(f, dx):
    out = np.zeros_like(f)
    out[1:-1] = (f[2:] - f[:-2]) / (2.0 * dx)
    out[0] = (f[1] - f[0]) / dx
    out[-1] = (f[-1] - f[-2]) / dx
    return out


def ddz(f, dz):
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    out[:, 0] = (f[:, 1] - f[:, 0]) / dz
    out[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return out


def lap(f, dx, dz):
    out = np.zeros_like(f)
    out[1:-1] += (f[2:] - 2.0 * f[1:-1] + f[:-2]) / (dx * dx)
    out[:, 1:-1] += (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / (dz * dz)
    out[0] = out[1]
    out[-1] = out[-2]
    out[:, 0] = out[:, 1]
    out[:, -1] = out[:, -2]
    return out


def edge_damp(f, coef=0.995, width=6):
    g = f.copy()
    for i in range(width):
        fac = coef ** (width - i)
        g[i, :] *= fac
        g[-1 - i, :] *= fac
        g[:, i] *= fac
        g[:, -1 - i] *= fac
    return g


def unwrap_2d_phase(phase):
    return np.unwrap(np.unwrap(phase, axis=0), axis=1)


def initial_beamlike_fields(x_half_range, z_half_range, nx, nz, alpha, k0, w0, phi0):
    params = default_params()
    params.alpha = alpha
    params.k0 = k0
    params.w0 = w0
    params.phi0 = phi0
    params.x_half_range = x_half_range
    params.z_half_range = z_half_range
    params.nx = nx
    params.nz = nz
    data = compute_xz_slice(0.0, params)
    rho = data["rho"] + 1e-8
    S = unwrap_2d_phase(np.angle(data["psi"]))
    return data["x"], data["z"], data["X"], data["Z"], rho, S


def compute_matter_observables(rho, S, tau, beta, m, dx, dz):
    Sx = ddx(S, dx)
    Sz = ddz(S, dz)
    E = np.sqrt(np.maximum(Sx * Sx + Sz * Sz + m * m * np.exp(2.0 * tau), 1e-10))
    W = np.exp(beta + tau)
    n = np.maximum(rho * E * W, 1e-10)
    return Sx, Sz, E, W, n


def rhs_reference(state, params):
    rho, S = state
    dx = params["dx"]
    dz = params["dz"]
    m = params["m"]
    Sx = ddx(S, dx)
    Sz = ddz(S, dz)
    E = np.sqrt(np.maximum(Sx * Sx + Sz * Sz + m * m, 1e-10))
    n = np.maximum(rho * E, 1e-10)
    fluxx = (n / E) * Sx
    fluxz = (n / E) * Sz
    n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
    S_t = -E
    return {"S_t": S_t, "n_t": n_t, "E": E, "n": n}


def rhs_branch_b(state, params):
    tau, ptau, beta, pbeta, rho, S = state
    dx = params["dx"]
    dz = params["dz"]
    m = params["m"]
    mp = params["mp"]

    Sx, Sz, E, W, n = compute_matter_observables(rho, S, tau, beta, m, dx, dz)
    mu = (m * m * np.exp(2.0 * tau) * rho) / (mp * mp)

    taux = ddx(tau, dx)
    tauz = ddz(tau, dz)
    betax = ddx(beta, dx)
    betaz = ddz(beta, dz)

    dot_tau = ptau * ptau - taux * taux - tauz * tauz
    dot_beta = pbeta * pbeta - betax * betax - betaz * betaz
    dot_btau = pbeta * ptau - betax * taux - betaz * tauz

    tau_tt = lap(tau, dx, dz) - 0.5 * (dot_tau + dot_beta + mu)
    beta_tt = lap(beta, dx, dz) - dot_btau - dot_beta - 0.5 * mu

    fluxx = (n / E) * Sx
    fluxz = (n / E) * Sz
    n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
    S_t = -E

    return {
        "tau_t": ptau,
        "ptau_t": tau_tt,
        "beta_t": pbeta,
        "pbeta_t": beta_tt,
        "S_t": S_t,
        "n_t": n_t,
        "E": E,
        "n": n,
    }


def U_of_phi(phi, ell):
    phi = np.clip(phi, 1e-6, 1.0)
    return -(np.power(phi, 1.0 / 3.0) - phi) / (ell * ell)


def Uprime_of_phi(phi, ell):
    phi = np.clip(phi, 1e-6, 1.0)
    return (1.0 - (1.0 / 3.0) * np.power(phi, -2.0 / 3.0)) / (ell * ell)


def rhs_branch_c(state, params):
    tau, ptau, beta, pbeta, phi, pphi, rho, S = state
    dx = params["dx"]
    dz = params["dz"]
    m = params["m"]
    mp = params["mp"]
    ell = params["ell"]

    Sx, Sz, E, W, n = compute_matter_observables(rho, S, tau, beta, m, dx, dz)
    mu = (m * m * np.exp(2.0 * tau) * rho) / (mp * mp)

    taux = ddx(tau, dx)
    tauz = ddz(tau, dz)
    betax = ddx(beta, dx)
    betaz = ddz(beta, dz)
    phix = ddx(phi, dx)
    phiz = ddz(phi, dz)

    dot_tau = ptau * ptau - taux * taux - tauz * tauz
    dot_beta = pbeta * pbeta - betax * betax - betaz * betaz
    dot_btau = pbeta * ptau - betax * taux - betaz * tauz
    dot_ptau = pphi * ptau - phix * taux - phiz * tauz
    dot_pbeta = pphi * pbeta - phix * betax - phiz * betaz

    U = U_of_phi(phi, ell)
    Up = Uprime_of_phi(phi, ell)

    A = phi * dot_tau + dot_ptau + 0.5 * np.exp(2.0 * tau) * U
    B = dot_tau + dot_btau + dot_beta + 0.5 * np.exp(2.0 * tau) * Up
    C = (
        2.0 * phi * dot_beta
        + 2.0 * phi * dot_btau
        + phi * dot_tau
        + 3.0 * dot_pbeta
        + 2.0 * dot_ptau
        + 1.5 * np.exp(2.0 * tau) * U
        + 3.0 * mu
    )

    box_tau = (C - 2.0 * A - 2.0 * phi * B) / (6.0 * np.clip(phi, 1e-6, None))
    box_beta = -B - 2.0 * box_tau
    box_phi = -A - 2.0 * phi * box_tau

    tau_tt = lap(tau, dx, dz) + box_tau
    beta_tt = lap(beta, dx, dz) + box_beta
    phi_tt = lap(phi, dx, dz) + box_phi

    fluxx = (n / E) * Sx
    fluxz = (n / E) * Sz
    n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
    S_t = -E

    return {
        "tau_t": ptau,
        "ptau_t": tau_tt,
        "beta_t": pbeta,
        "pbeta_t": beta_tt,
        "phi_t": pphi,
        "pphi_t": phi_tt,
        "S_t": S_t,
        "n_t": n_t,
        "E": E,
        "n": n,
    }


def integrate_reference(params):
    rho = params["rho0"].copy()
    S = params["S0"].copy()
    dt = params["dt"]
    steps = params["steps"]

    for _ in range(steps):
        rhs = rhs_reference((rho, S), params)
        S += dt * rhs["S_t"]
        n = rhs["n"] + dt * rhs["n_t"]
        n = edge_damp(n, coef=0.999)
        E = np.sqrt(
            np.maximum(ddx(S, params["dx"]) ** 2 + ddz(S, params["dz"]) ** 2 + params["m"] ** 2, 1e-10)
        )
        rho = np.maximum(n / E, 1e-10)
    return {"rho": rho, "S": S}


def integrate_branch_b(params):
    rho = params["rho0"].copy()
    S = params["S0"].copy()
    tau = np.zeros_like(rho)
    beta = np.zeros_like(rho)
    ptau = np.zeros_like(rho)
    pbeta = np.zeros_like(rho)
    dt = params["dt"]
    steps = params["steps"]

    for _ in range(steps):
        rhs = rhs_branch_b((tau, ptau, beta, pbeta, rho, S), params)
        ptau += dt * rhs["ptau_t"]
        pbeta += dt * rhs["pbeta_t"]
        tau += dt * ptau
        beta += dt * pbeta
        S += dt * rhs["S_t"]
        n = rhs["n"] + dt * rhs["n_t"]

        ptau = edge_damp(ptau)
        pbeta = edge_damp(pbeta)
        n = edge_damp(n, coef=0.999)

        E = np.sqrt(
            np.maximum(
                ddx(S, params["dx"]) ** 2
                + ddz(S, params["dz"]) ** 2
                + params["m"] ** 2 * np.exp(2.0 * tau),
                1e-10,
            )
        )
        rho = np.maximum(n / (np.exp(beta + tau) * E), 1e-10)
    return {"rho": rho, "S": S, "tau": tau, "beta": beta}


def integrate_branch_c(params):
    rho = params["rho0"].copy()
    S = params["S0"].copy()
    tau = np.zeros_like(rho)
    beta = np.zeros_like(rho)
    phi = np.ones_like(rho)
    ptau = np.zeros_like(rho)
    pbeta = np.zeros_like(rho)
    pphi = np.zeros_like(rho)
    dt = params["dt"]
    steps = params["steps"]

    for _ in range(steps):
        rhs = rhs_branch_c((tau, ptau, beta, pbeta, phi, pphi, rho, S), params)
        ptau += dt * rhs["ptau_t"]
        pbeta += dt * rhs["pbeta_t"]
        pphi += dt * rhs["pphi_t"]
        tau += dt * ptau
        beta += dt * pbeta
        phi += dt * pphi
        phi = np.clip(phi, 1e-4, 1.0)
        S += dt * rhs["S_t"]
        n = rhs["n"] + dt * rhs["n_t"]

        ptau = edge_damp(ptau)
        pbeta = edge_damp(pbeta)
        pphi = edge_damp(pphi)
        n = edge_damp(n, coef=0.999)

        E = np.sqrt(
            np.maximum(
                ddx(S, params["dx"]) ** 2
                + ddz(S, params["dz"]) ** 2
                + params["m"] ** 2 * np.exp(2.0 * tau),
                1e-10,
            )
        )
        rho = np.maximum(n / (np.exp(beta + tau) * E), 1e-10)
    return {"rho": rho, "S": S, "tau": tau, "beta": beta, "phi": phi}


def render_map(path, x, z, field, title):
    plt.figure(figsize=(5.6, 4.8))
    plt.imshow(
        field.T,
        origin="lower",
        extent=[x.min(), x.max(), z.min(), z.max()],
        aspect="equal",
        cmap="magma",
    )
    plt.colorbar()
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def render_cut(path, x, curves):
    plt.figure(figsize=(6.0, 3.8))
    for label, y in curves:
        plt.plot(x, y, label=label, linewidth=1.8)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, default=72)
    parser.add_argument("--nz", type=int, default=72)
    parser.add_argument("--range", dest="extent", type=float, default=5.0)
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--dt", type=float, default=0.003)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--k0", type=float, default=3.0)
    parser.add_argument("--w0", type=float, default=1.0)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--m", type=float, default=1.0)
    parser.add_argument("--mp", type=float, default=12.0)
    parser.add_argument("--ell", type=float, default=0.2)
    parser.add_argument("--outdir", type=str, default="E:\\量子势几何化_codex\\visualizations\\full_dynamics_2p1")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    x = np.linspace(-args.extent, args.extent, args.nx)
    z = np.linspace(-args.extent, args.extent, args.nz)
    dx = x[1] - x[0]
    dz = z[1] - z[0]
    x, z, _, _, rho0, S0 = initial_beamlike_fields(
        args.extent, args.extent, args.nx, args.nz, args.alpha, args.k0, args.w0, args.phi0
    )

    params = {
        "dx": dx,
        "dz": dz,
        "dt": args.dt,
        "steps": args.steps,
        "m": args.m,
        "mp": args.mp,
        "ell": args.ell,
        "rho0": rho0,
        "S0": S0,
    }

    ref = integrate_reference(params)
    b = integrate_branch_b(params)
    c = integrate_branch_c(params)

    render_map(outdir / "rho_reference.png", x, z, ref["rho"], "flat reference")
    render_map(outdir / "rho_R_tilde.png", x, z, b["rho"], "EH[g~]")
    render_map(outdir / "rho_fR_tilde.png", x, z, c["rho"], "F(R~)R~")
    render_map(outdir / "tau_R_tilde.png", x, z, b["tau"], "tau: EH[g~]")
    render_map(outdir / "tau_fR_tilde.png", x, z, c["tau"], "tau: F(R~)R~")
    render_map(outdir / "beta_R_tilde.png", x, z, b["beta"], "beta: EH[g~]")
    render_map(outdir / "beta_fR_tilde.png", x, z, c["beta"], "beta: F(R~)R~")
    render_map(outdir / "phi_fR_tilde.png", x, z, c["phi"], "Phi: F(R~)R~")

    mid = args.nz // 2
    render_cut(
        outdir / "rho_centerline_compare.png",
        x,
        [
            ("reference", ref["rho"][:, mid]),
            ("R_tilde", b["rho"][:, mid]),
            ("fR_tilde", c["rho"][:, mid]),
        ],
    )

    summary = {
        "grid": {"nx": args.nx, "nz": args.nz, "extent": args.extent, "dt": args.dt, "steps": args.steps},
        "params": {
            "alpha": args.alpha,
            "k0": args.k0,
            "w0": args.w0,
            "phi0": args.phi0,
            "m": args.m,
            "mp": args.mp,
            "ell": args.ell,
        },
        "rho_max": {
            "reference": float(ref["rho"].max()),
            "R_tilde": float(b["rho"].max()),
            "fR_tilde": float(c["rho"].max()),
        },
        "l1_differences": {
            "ref_vs_R_tilde": float(np.mean(np.abs(ref["rho"] - b["rho"]))),
            "ref_vs_fR_tilde": float(np.mean(np.abs(ref["rho"] - c["rho"]))),
            "R_tilde_vs_fR_tilde": float(np.mean(np.abs(b["rho"] - c["rho"]))),
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
