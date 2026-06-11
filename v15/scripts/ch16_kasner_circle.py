#!/usr/bin/env python3
"""ch16_kasner_circle.py — the Kasner circle and the BKL cobweb (Ch. 16, Figs. 16.2-16.3)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def kasner(u):
    d = 1 + u + u**2
    return -u / d, (1 + u) / d, u * (1 + u) / d

def main():
    # the circle: parametrize all (p1,p2) with sum p =1, sum p^2 = 1
    t = np.linspace(0, 2 * np.pi, 400)
    # solve: p1 + p2 + p3 = 1, p1^2+p2^2+p3^2 = 1 -> ellipse in (p1,p2)
    # param: p1 = 1/3 + r cos, p2 = 1/3 + r cos(t - 2pi/3) with r = 2/3 works:
    r = 2 / 3
    p1 = 1 / 3 + r * np.cos(t)
    p2 = 1 / 3 + r * np.cos(t - 2 * np.pi / 3)
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ax.plot(p1, p2, lw=1.6, label=r"Kasner circle ($\Sigma_K = 1$)")
    ax.scatter([1/3], [1/3], marker="o", s=40, color="#2ca02c", zorder=5, label=r"FRW ($\Sigma_K = 1/3$)")
    u2 = kasner(2.0)
    ax.scatter([u2[0]], [u2[1]], marker="*", s=130, color="#d95f00", zorder=5,
               label=fr"$u = 2$: $({u2[0]:.3f}, {u2[1]:.3f}, {u2[2]:.3f})$")
    for u in (1.0, 1.5, 3.0, 5.0, 10.0):
        p = kasner(u)
        ax.scatter([p[0]], [p[1]], s=12, color="k", zorder=4)
        ax.annotate(fr"$u={u:g}$", (p[0], p[1]), textcoords="offset points",
                    xytext=(6, 4), fontsize=7)
    th = np.linspace(0, 2*np.pi, 100)
    for rr in (0.2, 0.4):
        ax.plot(1/3 + rr*np.cos(th), 1/3 + rr*np.cos(th - 2*np.pi/3), color="#7f7f7f",
                lw=0.6, ls=":")
    ax.set_xlabel(r"$p_1$"), ax.set_ylabel(r"$p_2$")
    ax.set_title("The Kasner circle, its interior, and the $u$-parametrization")
    ax.legend(fontsize=8, loc="lower left")
    ax.set_aspect("equal")
    fig.savefig(FIGS / "ch16_fig2_kasner_circle.png")

    # BKL cobweb
    def bkl(u):
        return u - 1 if u >= 2 else 1 / (u - 1)
    us = [4.3]
    for _ in range(12):
        us.append(bkl(us[-1]))
    print("BKL orbit from u0=4.3:", np.round(us[:8], 4))
    x = np.linspace(1.001, 6, 500)
    fmap = np.where(x >= 2, x - 1, 1 / (x - 1))
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    ax.plot(x, fmap, lw=1.4, label=r"BKL map")
    ax.plot(x, x, "k:", lw=0.8)
    cx, cy = [], []
    u = 4.3
    for _ in range(10):
        v = bkl(u)
        cx += [u, v]; cy += [v, v]
        u = v
    ax.plot(cx, cy, color="#d95f00", lw=1.0, alpha=0.9, label=r"orbit of $u_0 = 4.3$")
    ax.set_xlabel(r"$u_k$"), ax.set_ylabel(r"$u_{k+1}$")
    ax.set_xlim(1, 6), ax.set_ylim(1, 6)
    ax.set_title("Kasner eras and reversals: chaos as a continued fraction")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch16_fig3_bkl_map.png")
    print("figures -> ch16_fig2_kasner_circle.png, ch16_fig3_bkl_map.png")

if __name__ == "__main__":
    main()
