#!/usr/bin/env python3
"""
ch02_iso_energy.py — iso-energy curves and the selection rule (Ch. 2, Figs. 2.1–2.2).
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
    L = np.linspace(0.4, 6, 400)
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for n0 in range(1, 7):
        ax.plot(L, n0 * L / 1.0, color="#1f77b4", lw=1.0, alpha=0.75)
    for n in range(1, 13):
        for k in range(1, 13):
            if 0.4 < k / n * n <= 6:
                pass
    # physical states: integer n on each curve n/L = const (take L0 = 1)
    for n0 in range(1, 7):
        ks = np.arange(1, 13)
        Ls = ks / n0
        m = (Ls >= 0.4) & (Ls <= 6)
        ax.plot(Ls[m], ks[m], "o", ms=3.4, color="#d95f00")
    ax.set_xlabel(r"box size $L$"), ax.set_ylabel(r"quantum number $n$")
    ax.set_title(r"Iso-energy curves $n/L = $ const; dots: physical (integer-$n$) states")
    fig.savefig(FIGS / "ch02_fig1_iso_energy_curves.png")

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    n0 = 3
    Lc = np.linspace(0.5, 2.4, 200)
    ax.plot(Lc, n0 * Lc, color="#1f77b4", lw=1.4, label=fr"$E = $ const curve through $(L_0, n_0 = {n0})$")
    pts = [(2/3, 2, r"$\hat T_-$"), (1.0, 3, ""), (4/3, 4, r"$\hat T_+$")]
    for Lp, np_, lbl in pts:
        ax.plot(Lp, np_, "o", ms=6, color="#d95f00")
    ax.annotate("", xy=(4/3 - 0.03, 4 - 0.06), xytext=(1.03, 3.06),
                arrowprops=dict(arrowstyle="->", lw=1.3, color="#2ca02c"))
    ax.annotate("", xy=(2/3 + 0.03, 2 + 0.06), xytext=(0.97, 2.94),
                arrowprops=dict(arrowstyle="->", lw=1.3, color="#7f7f7f"))
    ax.text(1.18, 3.62, r"$n \to n+1,\ L \to \frac{n+1}{n}L$", fontsize=9, color="#2ca02c")
    ax.text(0.55, 2.45, r"$n \to n-1,\ L \to \frac{n-1}{n}L$", fontsize=9, color="#7f7f7f")
    ax.set_xlabel(r"$L/L_0$"), ax.set_ylabel(r"$n$")
    ax.set_title("The Iso-Energy Selection Rule: the only nearest-neighbour moves")
    ax.legend(fontsize=8, loc="upper left")
    fig.savefig(FIGS / "ch02_fig2_selection_rule.png")
    print("figures -> ch02_fig1_iso_energy_curves.png, ch02_fig2_selection_rule.png")

if __name__ == "__main__":
    main()
