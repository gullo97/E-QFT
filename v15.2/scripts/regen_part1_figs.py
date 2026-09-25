#!/usr/bin/env python3
"""
regen_part1_figs.py — fresh v15-styled versions of three formerly-adapted figures:
ch03_fig3 (overlap anatomy), ch06_fig1 (spectator projection), ch08_fig1 (bag spectrum).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def c_closed(k, n, s):
    if abs(k - n * s) < 1e-12:
        return 1 / np.sqrt(s)
    return (2 * s**1.5 / np.pi) * ((-1)**n) * n * np.sin(k * np.pi / s) / (k**2 - n**2 * s**2)

def ch03_fig3():
    n, L = 3, 1.0
    s = (n + 1) / n
    Lp = s * L
    x = np.linspace(0, Lp, 600)
    old = np.where(x <= L, np.sqrt(2 / L) * np.sin(n * np.pi * np.clip(x, 0, L) / L), 0.0)
    new = np.sqrt(2 / Lp) * np.sin((n + 1) * np.pi * x / Lp)
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.4), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    ax.plot(x, old, label=fr"old mode $\varphi_{n}(x; L)$, extended by zero", lw=1.6)
    ax.plot(x, new, "--", label=fr"promoted mode $\varphi_{n+1}(x; L')$ — same wavelength", lw=1.6)
    ax.axvspan(L, Lp, color="#dddddd", alpha=0.5)
    ax.text(L + (Lp - L) / 2, 1.25, "fresh\nregion", fontsize=8, ha="center")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("$x$"), ax.legend(fontsize=7.5, loc="lower left")
    ax.set_title(fr"iso-energy step $n = {n} \to {n+1}$: wavelength preserved, coverage lost")
    ax = axes[1]
    ks = np.arange(1, 13)
    cs = np.array([c_closed(k, n, s)**2 for k in ks])
    ax.bar(ks, cs, width=0.7)
    ax.bar([n + 1], [cs[n]], width=0.7, color="#d95f00")
    ax.set_xlabel("$k$"), ax.set_ylabel(r"$|c_k|^2$")
    ax.set_title(fr"promoted peak $= n/(n{{+}}1) = {n/(n+1):.3f}$; sinc tail", fontsize=9)
    fig.savefig(FIGS / "ch03_fig3_overlap_wavefunctions.png"); plt.close(fig)
    print("  -> ch03_fig3_overlap_wavefunctions.png")

def ch06_fig1():
    n2 = 3
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), sharey=True)
    # expansion: s = (n1+1)/n1 with n1 = 5
    n1 = 5
    s = (n1 + 1) / n1
    ks = np.arange(1, 25)
    cs_e = np.array([c_closed(k, n2, s)**2 for k in ks])
    ax = axes[0]
    ax.bar(ks, cs_e, width=0.7)
    ax.set_title(fr"expansion ($s = {n1+1}/{n1}$): $\sum_k |c_k|^2 = {np.sum([c_closed(k, n2, s)**2 for k in range(1, 4001)]):.5f}$", fontsize=9)
    ax.set_xlabel("$k$"), ax.set_ylabel(r"$|c_k|^2$"), ax.set_xlim(0, 14)
    # contraction: L'' = (n1-1)/n1, overlaps of phi_{n2}(L) truncated onto basis of L''
    Lpp = (n1 - 1) / n1
    def c_con(k):
        kk, nn = k * np.pi / Lpp, n2 * np.pi
        if abs(kk - nn) < 1e-12:
            return np.sqrt(4 / Lpp) * (Lpp / 2 - np.sin(2 * kk * Lpp) / (4 * kk))
        return np.sqrt(4 / Lpp) * (np.sin((kk - nn) * Lpp) / (2 * (kk - nn)) - np.sin((kk + nn) * Lpp) / (2 * (kk + nn)))
    cs_c = np.array([c_con(k)**2 for k in ks])
    tot = np.sum([c_con(k)**2 for k in range(1, 4001)])
    ax = axes[1]
    ax.bar(ks, cs_c, width=0.7, color="#d95f00")
    Ncon = (n1 - 1) / n1 + np.sin(2 * np.pi * n2 / n1) / (2 * np.pi * n2)
    ax.set_title(fr"contraction: $\sum_k |c_k|^2 = {tot:.5f} = \mathcal{{N}}_{{\rm con}}$ ($= {Ncon:.5f}$)", fontsize=9)
    ax.set_xlabel("$k$"), ax.set_xlim(0, 14)
    ax.text(8, 0.5, fr"deficit $= {1-Ncon:.3f}$:" + "\nweight in the\nexcised region", fontsize=8)
    fig.suptitle(fr"A spectator ($n_2 = {n2}$) rides a geometry change driven by $n_1 = {n1}$", fontsize=10)
    fig.savefig(FIGS / "ch06_fig1_spectator_projection.png"); plt.close(fig)
    print(f"  -> ch06_fig1_spectator_projection.png (deficit check {1-tot:.5f} vs {1-Ncon:.5f})")

def ch08_fig1():
    L = 1.0
    mLs = np.logspace(-2, 2, 120)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    for j in range(1, 6):
        ps = []
        for mL in mLs:
            m = mL / L
            lo, hi = (j - 0.4999) * np.pi / L + 1e-9, (j + 0.4999) * np.pi / L
            try:
                ps.append(brentq(lambda p: np.tan(p * L) + p / m, lo, hi, xtol=1e-12) * L / np.pi)
            except ValueError:
                ps.append(np.nan)
        ax.semilogx(mLs, ps, lw=1.4)
        ax.axhline(j - 0.5, color="#7f7f7f", lw=0.5, ls=":")
        ax.axhline(j, color="#7f7f7f", lw=0.5, ls="--")
    ax.text(0.012, 0.62, "massless limit:\nhalf-integer modes", fontsize=8)
    ax.text(20, 0.78, "non-relativistic limit:\nDirichlet modes", fontsize=8)
    ax.set_xlabel(r"$mL$")
    ax.set_ylabel(r"$p_j L/\pi$")
    ax.set_title(r"MIT bag momenta: $\tan(pL) = -p/m$, between its two limits")
    fig.savefig(FIGS / "ch08_fig1_bag_spectrum.png"); plt.close(fig)
    print("  -> ch08_fig1_bag_spectrum.png")

if __name__ == "__main__":
    ch03_fig3(); ch06_fig1(); ch08_fig1()
