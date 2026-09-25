#!/usr/bin/env python3
"""
ch06_norms.py — norm theorems and the Expansion-Preference Theorem (Ch. 6).

Prints: N_con closed form vs quadrature; expansion norm = 1 and KE preservation;
sudden-contraction KE divergence rate; asymmetry vs N with sigma and N_crit.
Figures: ch06_fig2_norm_deficit, ch06_fig4_asymmetry_vs_N (regenerated).
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

def N_con(n1, n2):
    return (n1 - 1) / n1 + np.sin(2 * np.pi * n2 / n1) / (2 * np.pi * n2)

def main():
    print("=== Contraction Norm Theorem: closed form vs quadrature ===")
    for n1, n2 in [(5, 3), (8, 5), (12, 7)]:
        Lpp = (n1 - 1) / n1
        num = quad(lambda x: 2 * np.sin(n2 * np.pi * x)**2, 0, Lpp, limit=200)[0]
        print(f"  (n1,n2)=({n1},{n2}): closed {N_con(n1,n2):.5f}  quadrature {num:.5f}")

    print("=== Expansion: norm = 1 and KE preservation (4000-mode sums) ===")
    for n2, s in [(3, 1.25), (5, 1.2)]:
        ks = np.arange(1, 4001)
        def c(k):
            if abs(k - n2 * s) < 1e-12:
                return 1 / np.sqrt(s)
            return (2 * s**1.5 / np.pi) * ((-1)**n2) * n2 * np.sin(k * np.pi / s) / (k**2 - n2**2 * s**2)
        cs = np.array([c(k) for k in ks])
        norm = np.sum(cs**2)
        KE_new = np.sum(cs**2 * ks**2 / s**2)        # in units pi^2/2mL^2
        print(f"  n2={n2}, s={s}: norm {norm:.6f}  KE(new)/KE(old) = {KE_new/n2**2:.6f}")

    print("=== Sudden contraction: KE partial sums diverge linearly ===")
    n1, n2 = 5, 3
    Lpp = (n1 - 1) / n1
    ks = np.arange(1, 2001)
    cs = np.array([quad(lambda x: np.sqrt(2/Lpp)*np.sin(k*np.pi*x/Lpp)*np.sqrt(2)*np.sin(n2*np.pi*x),
                        0, Lpp, limit=300)[0] for k in ks[:200]])
    KE_partial = np.cumsum(cs**2 * (ks[:200] * np.pi / Lpp)**2)
    slope = (KE_partial[-1] - KE_partial[100]) / (200 - 101)
    print(f"  KE partial-sum slope per mode at large k: {slope:.3f} (nonzero => linear divergence)")

    # ---- Fig 6.2: deficit map
    n1s = np.arange(2, 31)
    n2s = np.arange(1, 31)
    Z = np.array([[1 - N_con(a, b) for a in n1s] for b in n2s])
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    pc = ax.pcolormesh(n1s, n2s, Z, shading="nearest", cmap="viridis")
    fig.colorbar(pc, label=r"deficit $1 - \mathcal{N}_{\mathrm{con}}$")
    ax.set_xlabel(r"$n_1$ (active)"), ax.set_ylabel(r"$n_2$ (spectator)")
    ax.set_title("Contraction is never free")
    fig.savefig(FIGS / "ch06_fig2_norm_deficit.png")

    # ---- Fig 6.4: asymmetry vs N
    rng = np.random.default_rng(2)
    print("=== Expansion-Preference: sigma and N_crit ===")
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for n1, color in [(5, "#1f77b4"), (12, "#d95f00")]:
        Ns = np.arange(0, 41)
        trials = []
        for _ in range(200):
            n2 = rng.integers(1, 4 * n1, size=Ns.max())
            logs = np.log(np.array([N_con(n1, b) for b in n2]))
            trials.append(-np.concatenate([[0], np.cumsum(logs)]))
        m = np.mean(trials, axis=0)
        sigma = (m[-1] - m[0]) / Ns.max()
        Ncrit = 6 * np.log(10) / sigma
        ax.plot(Ns, np.log(n1**2/(n1**2-1)) + m, color=color,
                label=fr"$n_1 = {n1}$: $\sigma = {sigma:.3f}$, $N_{{crit}}(10^6) \approx {Ncrit:.0f}$")
        print(f"  n1={n1}: sigma = {sigma:.3f} per spectator, N_crit(1e6) ~ {Ncrit:.0f}, N_crit(1e10) ~ {10*np.log(10)/sigma:.0f}")
    ax.axhline(6 * np.log(10), color="k", ls=":", lw=0.8)
    ax.text(1, 6 * np.log(10) + 0.3, r"suppression $10^{6}$", fontsize=8)
    ax.set_xlabel(r"spectator number $N$")
    ax.set_ylabel(r"$-\ln\,(P_{\mathrm{con}}/P_{\mathrm{exp}})$")
    ax.set_title("The arrow: exponential preference for expansion")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch06_fig4_asymmetry_vs_N.png")
    print("figures -> ch06_fig2_norm_deficit.png, ch06_fig4_asymmetry_vs_N.png")

if __name__ == "__main__":
    main()
