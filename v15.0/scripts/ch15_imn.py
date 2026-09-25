#!/usr/bin/env python3
"""
ch15_imn.py — mode-volume normalization diagnostics (Ch. 15, Figs. 15.1-15.2).

Prints: 1+1D closed form pi/8L vs partial sums; 3+1D plateau; the Casimir
negative result (zero-point sums NOT regulated).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def main():
    # --- 1+1D at x = L/2 (odd n only contribute)
    N = 5000
    n = np.arange(1, N + 1)
    odd = n % 2 == 1
    std = np.cumsum(np.where(odd, 1 / (np.pi * n), 0))
    mvn = np.cumsum(np.where(odd, 1 / (np.pi * n**2), 0))
    print(f"1+1D <phi^2>(L/2): partial sum to N={N}: {mvn[-1]:.6f}  vs  pi/8 = {np.pi/8:.6f}")
    print(f"1+1D standard partial sum (log-divergent): {std[-1]:.4f} and climbing")

    # --- per-mode decay (fig 15.1)
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    nn = n[odd][:200]
    ax.loglog(nn, 1 / (np.pi * nn), "o", ms=2.6, label=r"standard: $\sim n^{-1}$")
    ax.loglog(nn, 1 / (np.pi * nn**2), "s", ms=2.6, label=r"mode-volume: $\sim n^{-2}$")
    ax.loglog(nn, 1 / (np.pi * nn), "-", lw=0.7, color="#7f7f7f")
    ax.loglog(nn, 1 / (np.pi * nn**2), ":", lw=0.7, color="#7f7f7f")
    ax.set_xlabel("$n$"), ax.set_ylabel("per-mode contribution to $\\langle\\phi^2\\rangle$")
    ax.set_title("One extra power per mode: the entire mechanism")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch15_fig1_per_mode_decay.png")

    # --- plateaus (fig 15.2): 1+1D and 3+1D side by side
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    ax = axes[0]
    ax.semilogx(n, std, label="standard (log-divergent)")
    ax.semilogx(n, mvn, label="mode-volume")
    ax.axhline(np.pi / 8, color="k", ls=":", lw=0.9)
    ax.text(4, np.pi / 8 + 0.04, r"$\pi/8$", fontsize=9)
    ax.set_xlabel("modes kept $N$"), ax.set_ylabel(r"$\langle\phi^2(L/2)\rangle\,L$")
    ax.set_title("1+1D")
    ax.legend(fontsize=8)
    # 3+1D: all-odd lattice, cumulative in |n|
    M = 31
    g = np.arange(1, M + 1, 2)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2).ravel()
    r.sort()
    cum_std = np.cumsum(1 / r)
    cum_mvn = np.cumsum(1 / r**4)
    ax = axes[1]
    ax.plot(r, cum_std / cum_std[-1], label="standard (quadratic, normalized)")
    ax.plot(r, cum_mvn / cum_mvn[-1], label="mode-volume (plateaus)")
    ax.set_xlabel(r"$|\vec n|$"), ax.set_title("3+1D (cumulative, all-odd lattice)")
    ax.legend(fontsize=8)
    print(f"3+1D mode-volume sum saturated at |n| ~ {r[np.searchsorted(cum_mvn, 0.99*cum_mvn[-1])]:.0f} of {r[-1]:.0f}")
    fig.savefig(FIGS / "ch15_fig2_phi2_plateau.png")

    # --- Casimir negative result
    zp_std = np.cumsum(n.astype(float))            # ~ N^2
    zp_mvn = np.cumsum(np.ones_like(n, dtype=float))  # omega_n / V_n ~ const => ~ N
    print(f"zero-point sums at N={N}: standard ~ N^2 ({zp_std[-1]:.2e}), "
          f"mode-volume ~ N ({zp_mvn[-1]:.2e}) — both divergent: NOT regulated")
    print("figures -> ch15_fig1_per_mode_decay.png, ch15_fig2_phi2_plateau.png")

if __name__ == "__main__":
    main()
