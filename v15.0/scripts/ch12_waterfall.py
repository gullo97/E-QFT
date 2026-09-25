#!/usr/bin/env python3
"""
ch12_waterfall.py — spectral-flow waterfall (Ch. 12, Fig. 12.3, new).

Full spectrum of the two-wall chiral bag theta = (-2, +2), m = 1, versus box size L.
Mirror pairs everywhere except the single gap level executing the crossing at the
analytically predicted L* = artanh(-cos(Delta/2)/cos(Sigma))/m = 0.44302.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from ch08_bag_spectrum import find_spectrum

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def main():
    m, th0, thL = 1.0, -2.0, 2.0
    D, S = thL - th0, 0.5 * (th0 + thL)
    Lstar = np.arctanh(-np.cos(D / 2) / np.cos(S)) / m
    print(f"predicted crossing: mL* = {Lstar:.5f}")

    Ls = np.linspace(0.12, 1.0, 60)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    gap_track = []
    for L in Ls:
        spec = find_spectrum(m, 0.0, L, th0, thL, Emax=25.0, nE=3001)
        ax.plot([L] * len(spec), spec, ".", ms=1.8, color="#1f77b4", alpha=0.6)
        ingap = spec[np.abs(spec) < m]
        if len(ingap):
            gap_track.append((L, ingap[np.argmin(np.abs(ingap))]))
    if gap_track:
        gl = np.array(gap_track)
        ax.plot(gl[:, 0], gl[:, 1], "-", color="#d95f00", lw=1.8, label="gap level (spectral flow)")
        # print the crossing seen by the track
        sgn = np.sign(gl[:, 1])
        flips = np.where(sgn[:-1] * sgn[1:] < 0)[0]
        if len(flips):
            i = flips[0]
            Lc = gl[i, 0] + (gl[i + 1, 0] - gl[i, 0]) * abs(gl[i, 1]) / (abs(gl[i, 1]) + abs(gl[i + 1, 1]))
            print(f"observed crossing in track: L = {Lc:.5f}  (predicted {Lstar:.5f})")
    ax.axvline(Lstar, color="k", ls=":", lw=1.0, label=fr"$L^* = {Lstar:.5f}$ (analytic)")
    ax.axhspan(-m, m, color="#dddddd", alpha=0.35, label="mass gap")
    ax.set_xlabel(r"box size $L$")
    ax.set_ylabel(r"$E$")
    ax.set_ylim(-8, 8)
    ax.set_title(r"Spectral flow waterfall, $\theta = (-2, +2)$: one level carries the charge")
    ax.legend(fontsize=8, loc="upper right")
    fig.savefig(FIGS / "ch12_fig3_spectral_flow_waterfall.png")
    print(f"figure -> {FIGS / 'ch12_fig3_spectral_flow_waterfall.png'}")

if __name__ == "__main__":
    main()
