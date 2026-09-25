#!/usr/bin/env python3
"""
ch03_overlaps.py — Geometric Overlap Theorem + Dilation-Generator Lemma (Ch. 3).

Prints: promoted-mode overlap table sqrt(n/(n+m)); closed form (3.8) vs quadrature;
completeness; <n+1| x d/dx |n> vs 2n(n+1)/(2n+1) at n = 10, 40, 160.
Figures: ch03_fig1_sector_ladder (schematic), ch03_fig2_telescoping not regenerated
here (reuse), overlap anatomy left to the reused figure.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def c_closed(k, n, s):
    if abs(k - n * s) < 1e-12:
        return 1 / np.sqrt(s)
    return (2 * s**1.5 / np.pi) * ((-1)**n) * n * np.sin(k * np.pi / s) / (k**2 - n**2 * s**2)

def c_quad(k, n, s):
    f = lambda u: 2 / np.sqrt(s) * np.sin(k * np.pi * u / s) * np.sin(n * np.pi * u)
    return quad(f, 0, 1, limit=200)[0]

def main():
    print("=== promoted-mode overlaps sqrt(n/(n+m)) ===")
    for n, m in [(1, 1), (5, 1), (20, 1), (5, 3)]:
        s = (n + m) / n
        print(f"  n={n:3d}, m={m}: closed {np.sqrt(n/(n+m)):.6f}   formula-limit {c_closed(n+m, n, s):.6f}")

    print("=== closed form (3.8) vs quadrature ===")
    worst = 0
    for n, s in [(1, 2.0), (3, 1.5), (7, 1.2)]:
        for k in range(1, 12):
            worst = max(worst, abs(c_closed(k, n, s) - c_quad(k, n, s)))
    print(f"  max |closed - quadrature| = {worst:.2e}")

    print("=== completeness sum_k |c_k|^2 = 1 ===")
    for n, s in [(1, 2.0), (3, 1.5), (7, 1.2)]:
        tot = sum(c_closed(k, n, s)**2 for k in range(1, 4001))
        print(f"  n={n}, s={s}: sum = {tot:.8f}")

    print("=== Dilation-Generator Lemma <n+1|x d/dx|n> = 2n(n+1)/(2n+1) ===")
    for n in (10, 40, 160):
        # numeric: <m| x psi_n' > on [0, 1]
        m = n + 1
        f = lambda x: 2 * np.sin(m * np.pi * x) * x * n * np.pi * np.cos(n * np.pi * x)
        num = quad(f, 0, 1, limit=400)[0]
        print(f"  n={n:4d}: numeric {num:.4f}   formula {2*n*(n+1)/(2*n+1):.4f}")

    # ---- Fig 3.1: the sector ladder (signature schematic)
    fig, ax = plt.subplots(figsize=(7.4, 3.4))
    n0 = 3
    for k in range(6):
        L = (n0 + k) / n0
        x0 = sum((n0 + j) / n0 for j in range(k)) * 0.32 + 0.2 * k
        ax.add_patch(plt.Rectangle((x0, 0), L * 0.3, 1.15, fill=False, lw=1.6,
                                   edgecolor="#1f77b4"))
        xs = np.linspace(0, L * 0.3, 120)
        ax.plot(x0 + xs, 0.55 + 0.42 * np.sin((n0 + k) * np.pi * xs / (L * 0.3)),
                color="#d95f00", lw=1.1)
        ax.text(x0 + L * 0.15, -0.18, fr"$|{n0+k};\,\frac{{{n0+k}}}{{{n0}}}L_0\rangle$",
                ha="center", fontsize=8)
        if k < 5:
            ax.annotate("", xy=(x0 + L * 0.3 + 0.17, 0.6), xytext=(x0 + L * 0.3 + 0.02, 0.6),
                        arrowprops=dict(arrowstyle="->", lw=1.2))
            ax.text(x0 + L * 0.3 + 0.095, 0.72, r"$\iota$", fontsize=9, ha="center")
    ax.set_xlim(0, 8.3), ax.set_ylim(-0.45, 1.5)
    ax.axis("off")
    ax.set_title(r"The geometry-sector ladder: $\mathcal{H} = \bigoplus_L \mathcal{H}(L)$, iso-energy rungs, embeddings as arrows")
    fig.savefig(FIGS / "ch03_fig1_sector_ladder.png")
    print(f"figure -> ch03_fig1_sector_ladder.png")

if __name__ == "__main__":
    main()
