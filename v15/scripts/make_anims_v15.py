#!/usr/bin/env python3
"""
make_anims_v15.py — the two new v15 animations.

ch04_anim_spreading.gif : a definite geometry delocalizing on the ladder (Bessel fronts).
ch12_anim_pump.gif      : the spectral-flow pump — spectrum vs L beside Q_vac(L) filling in.

Usage: python3 make_anims_v15.py [ch04|ch12]   (no arg: both)
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def anim_ch04():
    N, J, n0 = 200, 1.0, 100
    H = np.zeros((N, N))
    for i in range(N - 1):
        H[i, i + 1] = H[i + 1, i] = J
    w, U = np.linalg.eigh(H)
    psi0 = np.zeros(N); psi0[n0] = 1.0
    c0 = U.T @ psi0
    ts = np.linspace(0, 36, 48)
    P = np.array([np.abs(U @ (np.exp(-1j * w * t) * c0))**2 for t in ts])

    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    line, = ax.plot([], [], lw=1.3)
    fl = ax.axvline(0, color="#d95f00", ls=":", lw=1.0)
    fr = ax.axvline(0, color="#d95f00", ls=":", lw=1.0)
    txt = ax.text(0.02, 0.92, "", transform=ax.transAxes, fontsize=9)
    ax.set_xlim(-90, 90), ax.set_ylim(0, 0.28)
    ax.set_xlabel(r"$n - n_0$ (ladder position = geometry)")
    ax.set_ylabel(r"$|\psi_n(t)|^2$")
    ax.set_title("A definite geometry delocalizing: ballistic fronts at $\\pm 2Jt$")
    x = np.arange(N) - n0

    def update(i):
        line.set_data(x, P[i])
        fl.set_xdata([2 * J * ts[i]]); fr.set_xdata([-2 * J * ts[i]])
        txt.set_text(fr"$tJ = {ts[i]:.1f}$")
        return line, fl, fr, txt

    ani = FuncAnimation(fig, update, frames=len(ts), blit=True)
    ani.save(FIGS / "ch04_anim_spreading.gif", writer=PillowWriter(fps=12))
    plt.close(fig)
    print("  -> ch04_anim_spreading.gif")

def anim_ch12():
    from ch08_bag_spectrum import find_spectrum
    m, th0, thL = 1.0, -2.0, 2.0
    Lstar = np.arctanh(-np.cos((thL - th0) / 2) / np.cos((th0 + thL) / 2)) / m
    Ls = np.linspace(0.16, 0.92, 40)
    specs = [find_spectrum(m, 0.0, L, th0, thL, Emax=12.0, nE=1500) for L in Ls]
    # Q_vac: exact plateaus from the fractional law (validated in ch12_pump_control.py)
    qlo, qhi = +0.3634, -0.6366
    Q = np.where(Ls < Lstar, qlo, qhi)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6), gridspec_kw={"width_ratios": [1.15, 1]})
    axL, axR = axes
    axL.set_xlim(-0.05, 1.05), axL.set_ylim(-9, 9)
    axL.axhspan(-m, m, color="#dddddd", alpha=0.4)
    axL.set_xticks([]), axL.set_ylabel("$E$")
    axL.set_title(r"spectrum at current $L$ ($\theta = (-2, +2)$)")
    lines = [axL.plot([0.1, 0.9], [0, 0], lw=1.3, color="#1f77b4", alpha=0.0)[0] for _ in range(40)]
    gapline, = axL.plot([0.1, 0.9], [0, 0], lw=2.2, color="#d95f00", alpha=0.0)
    zline = axL.axhline(0, color="k", ls="--", lw=0.7)
    axR.set_xlim(Ls[0], Ls[-1]), axR.set_ylim(-0.85, 0.6)
    axR.axvline(Lstar, color="k", ls=":", lw=0.9)
    axR.text(Lstar + 0.01, 0.45, fr"$L^* = {Lstar:.3f}$", fontsize=8)
    axR.set_xlabel("$L$"), axR.set_ylabel(r"$Q_{\rm vac}(L)$")
    axR.set_title("the vacuum charge, jumping by one unit")
    trace, = axR.plot([], [], lw=1.6, color="#2ca02c")
    dot, = axR.plot([], [], "o", ms=6, color="#2ca02c")
    ttl = axL.text(0.02, 0.93, "", transform=axL.transAxes, fontsize=9)

    def update(i):
        sp = specs[i]
        for k, ln in enumerate(lines):
            if k < len(sp):
                E = sp[k]
                ln.set_ydata([E, E]); ln.set_alpha(0.75)
                ingap = abs(E) < m
                ln.set_color("#d95f00" if ingap else "#1f77b4")
                ln.set_linewidth(2.2 if ingap else 1.1)
            else:
                ln.set_alpha(0.0)
        trace.set_data(Ls[:i + 1], Q[:i + 1])
        dot.set_data([Ls[i]], [Q[i]])
        ttl.set_text(fr"$L = {Ls[i]:.3f}$")
        return lines + [trace, dot, ttl]

    ani = FuncAnimation(fig, update, frames=len(Ls), blit=False)
    ani.save(FIGS / "ch12_anim_pump.gif", writer=PillowWriter(fps=8))
    plt.close(fig)
    print("  -> ch12_anim_pump.gif")

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("ch04", "both"):
        anim_ch04()
    if which in ("ch12", "both"):
        anim_ch12()
