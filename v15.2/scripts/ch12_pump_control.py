#!/usr/bin/env python3
"""
ch12_pump_control.py — the decisive experiment, spectral side (Ch. 12, Fig. 12.4).

Q_vac(L) = -eta/2 for the crossing configuration theta = (-2, +2) and the control
(0.7, 0.7), over the quench window 0.30 -> 0.70. Pump = Q_vac(0.30) - Q_vac(0.70).
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

def qvac(L, th0, thL, m=1.0, Emax=120.0):
    spec = find_spectrum(m, 0.0, L, th0, thL, Emax=Emax, nE=6001)
    ts = np.array([0.08, 0.04, 0.02])
    etas = np.array([np.sum(np.sign(spec) * np.exp(-t * np.abs(spec))) for t in ts])
    # Richardson (linear in t) on the last two
    eta0 = etas[-1] + (etas[-1] - etas[-2])
    return -0.5 * eta0

def main():
    Ls = np.linspace(0.26, 0.78, 14)
    q_cross = [qvac(L, -2.0, 2.0) for L in Ls]
    q_ctrl = [qvac(L, 0.7, 0.7) for L in Ls]
    pump = qvac(0.30, -2.0, 2.0) - qvac(0.70, -2.0, 2.0)
    ctrl = qvac(0.30, 0.7, 0.7) - qvac(0.70, 0.7, 0.7)
    print(f"pump   (theta = (-2,+2), L: 0.30 -> 0.70): DQ_net = {pump:+.4f}   (exact: +1)")
    print(f"control(theta = (0.7,0.7), same quench):   DQ_net = {ctrl:+.4f}   (exact: 0)")

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(Ls, q_cross, "o-", ms=4, label=r"$\theta = (-2, +2)$: crossing at $L^* = 0.443$")
    ax.plot(Ls, q_ctrl, "s-", ms=4, label=r"$\theta = (0.7, 0.7)$: control (no mismatch)")
    ax.axvline(0.44302, color="k", ls=":", lw=0.9)
    for L0 in (0.30, 0.70):
        ax.axvline(L0, color="#7f7f7f", lw=0.7, alpha=0.7)
    ax.annotate("", xy=(0.70, q_cross[0] + 0.04), xytext=(0.30, q_cross[0] + 0.04),
                arrowprops=dict(arrowstyle="->", lw=1.2, color="#2ca02c"))
    ax.text(0.5, q_cross[0] + 0.10, fr"quench: $\Delta Q_{{\rm net}} = {pump:+.3f}$",
            ha="center", fontsize=9, color="#2ca02c")
    ax.set_xlabel(r"box size $L$")
    ax.set_ylabel(r"$Q_{\mathrm{vac}}(L)$")
    ax.set_title("The pump and its control: one quantized unit vs exactly zero")
    ax.legend(fontsize=8, loc="center right")
    fig.savefig(FIGS / "ch12_fig4_pump_vs_control.png")
    print("figure -> ch12_fig4_pump_vs_control.png")

if __name__ == "__main__":
    main()
