#!/usr/bin/env python3
"""
ch10_mirror_spectrum.py — Spectral Mirror Theorem data (Ch. 10, Fig. 10.3).

Bulk-phase MIT bag spectra vs delta: every level moves (m -> m cos d), levels move
in mirror pairs, eta stays at machine zero; overlay = roots of tan(pL) = -p/(m cos d).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq
from ch08_bag_spectrum import find_spectrum, eta_heat

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def tan_roots(mR, L, nmax=6):
    ps = []
    for j in range(1, nmax + 1):
        lo, hi = (j - 0.4999) * np.pi / L + 1e-9, (j + 0.4999) * np.pi / L
        try:
            ps.append(brentq(lambda p: np.tan(p * L) + p / mR, lo, hi, xtol=1e-13))
        except ValueError:
            pass
    return np.array(ps)

def main():
    m, L = 1.0, 1.0
    deltas = np.linspace(0, 1.35, 16)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for d in deltas:
        mR, mI = m * np.cos(d), m * np.sin(d)
        sp = find_spectrum(mR, mI, L, 0.0, 0.0, Emax=14.0, nE=3001)
        ax.plot([d] * len(sp), sp, ".", ms=2.6, color="#1f77b4")
        ps = tan_roots(mR, L)
        Ep = np.concatenate([np.sqrt(ps**2 + m**2), -np.sqrt(ps**2 + m**2)])
        ax.plot([d] * len(Ep), Ep, "o", ms=5.0, mfc="none", mec="#d95f00", mew=0.7)
    for d in (0.3, np.pi / 4):
        mR, mI = m * np.cos(d), m * np.sin(d)
        sp = find_spectrum(mR, mI, L, 0.0, 0.0, Emax=60.0, nE=9001)
        print(f"delta = {d:.3f}: levels = {len(sp)}, eta(t=0.03) = {eta_heat(sp, 0.03):+.2e}")
    ax.axhspan(-m, m, color="#dddddd", alpha=0.35)
    ax.set_xlabel(r"CP phase $\delta$")
    ax.set_ylabel(r"$E$")
    ax.set_ylim(-9, 9)
    ax.set_title(r"The mirror at work: spectra vs $\delta$ (dots: solver; circles: $m \to m\cos\delta$ bag)")
    fig.savefig(FIGS / "ch10_fig3_mirror_spectrum.png")
    print("figure -> ch10_fig3_mirror_spectrum.png")

if __name__ == "__main__":
    main()
