import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent))

from kg_double_slit_beamlike import compute_xz_slice, default_params


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


def edge_damp_interior(f, coef=0.998, width=8):
    g = f.copy()
    nx, nz = f.shape
    for i in range(1, width + 1):
        fac = coef ** (width + 1 - i)
        g[i, :] *= fac
        g[nx - 1 - i, :] *= fac
        g[:, 1 + i] *= fac
        g[:, nz - 1 - i] *= fac
    return g


def apply_boundary_from_target(field, target):
    field[0, :] = target[0, :]
    field[-1, :] = target[-1, :]
    field[:, 0] = target[:, 0]
    field[:, -1] = target[:, -1]


def apply_geometry_boundary(tau, beta, phi=None):
    for f in (tau, beta):
        f[0, :] = 0.0
        f[-1, :] = 0.0
        f[:, 0] = 0.0
        f[:, -1] = 0.0
    if phi is not None:
        phi[0, :] = 1.0
        phi[-1, :] = 1.0
        phi[:, 0] = 1.0
        phi[:, -1] = 1.0


def clip_field(f, lo, hi):
    np.clip(f, lo, hi, out=f)


def analytic_slice(extent, nx, nz, alpha, k0, w0, phi0, t):
    params = default_params()
    params.alpha = alpha
    params.k0 = k0
    params.w0 = w0
    params.phi0 = phi0
    params.x_half_range = extent
    params.z_half_range = extent
    params.nx = nx
    params.nz = nz
    data = compute_xz_slice(t, params)
    psi = data["psi"]
    rho = np.abs(psi) ** 2 + 1e-10
    s = np.unwrap(np.unwrap(np.angle(psi), axis=0), axis=1)
    return data["x"], data["z"], psi, rho, s


def render_map(path, x, z, field, title):
    plt.figure(figsize=(5.6, 4.8))
    plt.imshow(field.T, origin="lower", extent=[x.min(), x.max(), z.min(), z.max()], aspect="equal", cmap="magma")
    plt.colorbar(label="rho")
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def render_cut(path, x, curves):
    plt.figure(figsize=(6.2, 3.8))
    for label, y in curves:
        plt.plot(x, y, label=label, linewidth=1.8)
    plt.xlabel("x")
    plt.ylabel("rho(x, z=0)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def compute_matter_observables(rho, s, tau, beta, m, dx, dz):
    sx = ddx(s, dx)
    sz = ddz(s, dz)
    e = np.sqrt(np.maximum(sx * sx + sz * sz + m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)), 1e-12))
    w = np.exp(np.clip(beta + tau, -20.0, 20.0))
    n = np.maximum(rho * e * w, 1e-12)
    return sx, sz, e, w, n


def rhs_branch_b(state, params):
    tau, ptau, beta, pbeta, rho, s = state
    dx = params["dx"]
    dz = params["dz"]
    m = params["m"]
    mp = params["mp"]

    sx, sz, e, _, n = compute_matter_observables(rho, s, tau, beta, m, dx, dz)
    mu = (m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)) * rho) / (mp * mp)

    taux = ddx(tau, dx)
    tauz = ddz(tau, dz)
    betax = ddx(beta, dx)
    betaz = ddz(beta, dz)

    dot_tau = ptau * ptau - taux * taux - tauz * tauz
    dot_beta = pbeta * pbeta - betax * betax - betaz * betaz
    dot_btau = pbeta * ptau - betax * taux - betaz * tauz

    tau_tt = lap(tau, dx, dz) - 0.5 * (dot_tau + dot_beta + mu)
    beta_tt = lap(beta, dx, dz) - dot_btau - dot_beta - 0.5 * mu

    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
    s_t = -e
    return {"tau_tt": tau_tt, "beta_tt": beta_tt, "n_t": n_t, "s_t": s_t, "n": n}


def U_of_phi(phi, ell):
    phi = np.clip(phi, 1e-6, 1.0)
    return -(np.power(phi, 1.0 / 3.0) - phi) / (ell * ell)


def uprime_of_phi(phi, ell):
    phi = np.clip(phi, 1e-6, 1.0)
    return (1.0 - (1.0 / 3.0) * np.power(phi, -2.0 / 3.0)) / (ell * ell)


def rhs_branch_c(state, params):
    tau, ptau, beta, pbeta, phi, pphi, rho, s = state
    dx = params["dx"]
    dz = params["dz"]
    m = params["m"]
    mp = params["mp"]
    ell = params["ell"]

    sx, sz, e, _, n = compute_matter_observables(rho, s, tau, beta, m, dx, dz)
    mu = (m * m * np.exp(np.clip(2.0 * tau, -20.0, 20.0)) * rho) / (mp * mp)

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

    u = U_of_phi(phi, ell)
    up = uprime_of_phi(phi, ell)
    exp2t = np.exp(np.clip(2.0 * tau, -20.0, 20.0))

    a = phi * dot_tau + dot_ptau + 0.5 * exp2t * u
    b = dot_tau + dot_btau + dot_beta + 0.5 * exp2t * up
    c = (
        2.0 * phi * dot_beta
        + 2.0 * phi * dot_btau
        + phi * dot_tau
        + 3.0 * dot_pbeta
        + 2.0 * dot_ptau
        + 1.5 * exp2t * u
        + 3.0 * mu
    )

    box_tau = (c - 2.0 * a - 2.0 * phi * b) / (6.0 * np.clip(phi, 1e-6, None))
    box_beta = -b - 2.0 * box_tau
    box_phi = -a - 2.0 * phi * box_tau

    tau_tt = lap(tau, dx, dz) + box_tau
    beta_tt = lap(beta, dx, dz) + box_beta
    phi_tt = lap(phi, dx, dz) + box_phi

    fluxx = (n / e) * sx
    fluxz = (n / e) * sz
    n_t = -(ddx(fluxx, dx) + ddz(fluxz, dz))
    s_t = -e
    return {"tau_tt": tau_tt, "beta_tt": beta_tt, "phi_tt": phi_tt, "n_t": n_t, "s_t": s_t, "n": n}


def benchmark_branch_a(params):
    # A branch is the "no Bekenstein transform" reference branch.
    # For now we treat the user-approved beamlike solution as the benchmark
    # preparation instead of mixing in a second, inconsistent numerical IBVP.
    _, _, psi, rho, _ = analytic_slice(
        params["extent"], params["nx"], params["nz"], params["alpha"], params["k0"], params["w0"], params["phi0"], 0.0
    )
    return {"rho": rho, "psi": psi}


def integrate_branch_b(params):
    _, _, _, rho, s = analytic_slice(
        params["extent"], params["nx"], params["nz"], params["alpha"], params["k0"], params["w0"], params["phi0"], 0.0
    )
    tau = np.zeros_like(rho)
    beta = np.zeros_like(rho)
    ptau = np.zeros_like(rho)
    pbeta = np.zeros_like(rho)

    for step in range(params["steps"]):
        t_next = (step + 1) * params["dt"]
        _, _, _, rho_b, s_b = analytic_slice(
            params["extent"], params["nx"], params["nz"], params["alpha"], params["k0"], params["w0"], params["phi0"], t_next
        )
        rhs = rhs_branch_b((tau, ptau, beta, pbeta, rho, s), params)
        ptau += params["dt"] * rhs["tau_tt"]
        pbeta += params["dt"] * rhs["beta_tt"]
        tau += params["dt"] * ptau
        beta += params["dt"] * pbeta
        s += params["dt"] * rhs["s_t"]
        n = rhs["n"] + params["dt"] * rhs["n_t"]

        ptau = edge_damp_interior(ptau)
        pbeta = edge_damp_interior(pbeta)
        n = edge_damp_interior(n, coef=0.999)
        clip_field(ptau, -10.0, 10.0)
        clip_field(pbeta, -10.0, 10.0)
        clip_field(tau, -2.0, 2.0)
        clip_field(beta, -2.0, 2.0)

        apply_geometry_boundary(tau, beta)
        apply_boundary_from_target(s, s_b)

        _, _, e, _, _ = compute_matter_observables(rho, s, tau, beta, params["m"], params["dx"], params["dz"])
        n_target = rho_b * e * np.exp(np.clip(beta + tau, -20.0, 20.0))
        apply_boundary_from_target(n, n_target)
        rho = np.maximum(n / (np.exp(np.clip(beta + tau, -20.0, 20.0)) * e), 1e-10)
        apply_boundary_from_target(rho, rho_b)
    return {"rho": rho, "tau": tau, "beta": beta}


def integrate_branch_c(params):
    _, _, _, rho, s = analytic_slice(
        params["extent"], params["nx"], params["nz"], params["alpha"], params["k0"], params["w0"], params["phi0"], 0.0
    )
    tau = np.zeros_like(rho)
    beta = np.zeros_like(rho)
    phi = np.ones_like(rho)
    ptau = np.zeros_like(rho)
    pbeta = np.zeros_like(rho)
    pphi = np.zeros_like(rho)

    for step in range(params["steps"]):
        t_next = (step + 1) * params["dt"]
        _, _, _, rho_b, s_b = analytic_slice(
            params["extent"], params["nx"], params["nz"], params["alpha"], params["k0"], params["w0"], params["phi0"], t_next
        )
        rhs = rhs_branch_c((tau, ptau, beta, pbeta, phi, pphi, rho, s), params)
        ptau += params["dt"] * rhs["tau_tt"]
        pbeta += params["dt"] * rhs["beta_tt"]
        pphi += params["dt"] * rhs["phi_tt"]
        tau += params["dt"] * ptau
        beta += params["dt"] * pbeta
        phi += params["dt"] * pphi
        s += params["dt"] * rhs["s_t"]
        n = rhs["n"] + params["dt"] * rhs["n_t"]

        ptau = edge_damp_interior(ptau)
        pbeta = edge_damp_interior(pbeta)
        pphi = edge_damp_interior(pphi)
        n = edge_damp_interior(n, coef=0.999)
        clip_field(ptau, -10.0, 10.0)
        clip_field(pbeta, -10.0, 10.0)
        clip_field(pphi, -10.0, 10.0)
        clip_field(tau, -2.0, 2.0)
        clip_field(beta, -2.0, 2.0)
        clip_field(phi, 0.2, 1.0)

        apply_geometry_boundary(tau, beta, phi)
        apply_boundary_from_target(s, s_b)

        _, _, e, _, _ = compute_matter_observables(rho, s, tau, beta, params["m"], params["dx"], params["dz"])
        n_target = rho_b * e * np.exp(np.clip(beta + tau, -20.0, 20.0))
        apply_boundary_from_target(n, n_target)
        rho = np.maximum(n / (np.exp(np.clip(beta + tau, -20.0, 20.0)) * e), 1e-10)
        apply_boundary_from_target(rho, rho_b)
    return {"rho": rho, "tau": tau, "beta": beta, "phi": phi}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, default=100)
    parser.add_argument("--nz", type=int, default=100)
    parser.add_argument("--extent", type=float, default=5.0)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--dt", type=float, default=0.002)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--k0", type=float, default=10.0)
    parser.add_argument("--w0", type=float, default=1.0)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--m", type=float, default=1.0)
    parser.add_argument("--mp", type=float, default=300.0)
    parser.add_argument("--ell", type=float, default=0.02)
    parser.add_argument("--outdir", type=str, default=str(Path.cwd() / "visualizations" / "boundary_driven_2p1"))
    args = parser.parse_args()

    x = np.linspace(-args.extent, args.extent, args.nx)
    z = np.linspace(-args.extent, args.extent, args.nz)
    params = {
        "nx": args.nx,
        "nz": args.nz,
        "extent": args.extent,
        "dx": x[1] - x[0],
        "dz": z[1] - z[0],
        "steps": args.steps,
        "dt": args.dt,
        "alpha": args.alpha,
        "k0": args.k0,
        "w0": args.w0,
        "phi0": args.phi0,
        "m": args.m,
        "mp": args.mp,
        "ell": args.ell,
    }

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    a = benchmark_branch_a(params)
    b = integrate_branch_b(params)
    c = integrate_branch_c(params)

    render_map(outdir / "rho_A_flatKG.png", x, z, a["rho"], "A: analytic flat KG benchmark")
    render_map(outdir / "rho_B_Rtilde.png", x, z, b["rho"], "B: EH[g~]")
    render_map(outdir / "rho_C_fRtilde.png", x, z, c["rho"], "C: F(R~)R~")

    mid = args.nz // 2
    render_cut(
        outdir / "rho_centerline_compare.png",
        x,
        [("A", a["rho"][:, mid]), ("B", b["rho"][:, mid]), ("C", c["rho"][:, mid])],
    )

    summary = {
        "rho_max": {"A": float(np.nanmax(a["rho"])), "B": float(np.nanmax(b["rho"])), "C": float(np.nanmax(c["rho"]))},
        "l1": {
            "A_vs_B": float(np.nanmean(np.abs(a["rho"] - b["rho"]))),
            "A_vs_C": float(np.nanmean(np.abs(a["rho"] - c["rho"]))),
            "B_vs_C": float(np.nanmean(np.abs(b["rho"] - c["rho"]))),
        },
        "params": {
            "nx": args.nx,
            "nz": args.nz,
            "extent": args.extent,
            "dx": params["dx"],
            "dz": params["dz"],
            "steps": args.steps,
            "dt": args.dt,
            "alpha": args.alpha,
            "k0": args.k0,
            "w0": args.w0,
            "phi0": args.phi0,
            "m": args.m,
            "mp": args.mp,
            "ell": args.ell,
        },
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
