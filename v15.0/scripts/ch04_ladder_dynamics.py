#!/usr/bin/env python3
"""
ch04_ladder_dynamics.py — band, Bessel spreading, geometry superpositions (Ch. 4).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.special import jv

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def c_closed(k, n, s):
    if abs(k - n * s) < 1e-12:
        return 1 / np.sqrt(s)
    return (2 * s**1.5 / np.pi) * ((-1)**n) * n * np.sin(k * np.pi / s) / (k**2 - n**2 * s**2)

def main():
    # --- fig 4.1: band
    q = np.linspace(0, np.pi, 300)
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(q, -2 * np.cos(q))
    ax.plot(q, 2 * np.sin(q), "--", lw=1.0, label=r"group velocity $v_g = 2J\sin q$")
    ax.axvline(np.pi / 2, color="k", ls=":", lw=0.8)
    ax.text(np.pi / 2 + 0.05, -1.7, "band center:\nfastest size change", fontsize=8)
    ax.set_xlabel(r"quasimomentum $q$ (conjugate to ladder position)")
    ax.set_ylabel(r"$E(q)/J$")
    ax.set_title(r"The geometry band: $E(q) = -2J\cos q$")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch04_fig1_band_structure.png")

    # --- fig 4.2: Bessel carpet (exact evolution vs propagator)
    N, J, n0 = 240, 1.0, 120
    H = np.zeros((N, N))
    for i in range(N - 1):
        H[i, i + 1] = H[i + 1, i] = J
    w, U = np.linalg.eigh(H)
    psi0 = np.zeros(N); psi0[n0] = 1.0
    c0 = U.T @ psi0
    ts = np.linspace(0, 40, 160)
    P = np.array([np.abs(U @ (np.exp(-1j * w * t) * c0))**2 for t in ts])
    # check against Bessel propagator at one time
    t = 25.0
    psi = U @ (np.exp(-1j * w * t) * c0)
    d = np.arange(N) - n0
    bess = np.abs(jv(d, 2 * J * t))**2
    err = np.max(np.abs(np.abs(psi)**2 - bess)[20:-20])
    print(f"Bessel propagator check at t={t}: max |P_exact - J_d^2| = {err:.2e}")

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    pc = ax.pcolormesh(np.arange(N) - n0, ts, np.log10(P + 1e-12), shading="nearest",
                       cmap="magma", vmin=-6, vmax=0)
    ax.plot(2 * J * ts, ts, "w--", lw=0.9); ax.plot(-2 * J * ts, ts, "w--", lw=0.9)
    fig.colorbar(pc, label=r"$\log_{10}|\psi_n(t)|^2$")
    ax.set_xlabel(r"$n - n_0$"), ax.set_ylabel(r"$t\,J$")
    ax.set_title("A definite geometry delocalizing: ballistic fronts at $\\pm 2Jt$")
    fig.savefig(FIGS / "ch04_fig2_bessel_spreading.png")

    # --- fig 4.3: superpositions + participation ratio
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2))
    for ax, n in ((axes[0], 1), (axes[1], 100)):
        s = (n + 1) / n
        ks = np.arange(1, max(12, n + 8))
        cs = np.array([c_closed(k, n, s)**2 for k in ks])
        ax.bar(ks - n, cs, width=0.8)
        ax.set_title(fr"$n = {n}$"), ax.set_xlabel(r"$k - n$")
        ax.set_xlim(-6, 8)
    axes[0].set_ylabel(r"$|c_k|^2$")
    ns = np.arange(1, 120)
    R = []
    for n in ns:
        s = (n + 1) / n
        cs = np.array([c_closed(k, n, s)**2 for k in range(1, n + 400)])
        R.append(1.0 / np.sum(cs**2))
    axes[2].semilogx(ns, R)
    axes[2].axhline(1, color="k", ls=":", lw=0.8)
    axes[2].set_xlabel(r"$n$"), axes[2].set_ylabel(r"participation ratio $\mathcal{R}(n)$")
    axes[2].set_title("quantum $\\to$ classical")
    print(f"participation ratio R(1) = {R[0]:.3f}, R(5) = {R[4]:.3f}, R(100) = {R[99]:.3f}")
    fig.savefig(FIGS / "ch04_fig3_geometry_superposition.png")
    print("figures -> ch04_fig1/2/3")

if __name__ == "__main__":
    main()
