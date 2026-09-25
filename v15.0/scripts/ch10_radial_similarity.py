#!/usr/bin/env python3
"""
ch10_radial_similarity.py — Bulk-Phase Dressing Theorem, 3+1D radial check (Ch. 10, Fig. 10.2).

Integrates the same-sign-m_I radial system (10.9) directly (kappa = -1, s-wave) and
compares: (i) eigenvalues against the standard MIT radial bag with m -> m cos(delta);
(ii) the wavefunction ratio g_delta / g_0 against e^{m_I r}.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")
KAPPA = -1

def rhs(r, y, E, mR, mI):
    g, f = y
    dg = (E + mR) * f + mI * g - (1 + KAPPA) / r * g
    df = -(E - mR) * g + mI * f - (1 - KAPPA) / r * f
    return [dg, df]

def shoot(E, mR, mI, R, r0=1e-6):
    # regular series at the origin for kappa = -1: g ~ 1, f ~ -(E - mR) r / 3
    y0 = [1.0, -(E - mR) * r0 / 3.0]
    sol = solve_ivp(rhs, (r0, R), y0, args=(E, mR, mI), rtol=1e-10, atol=1e-12,
                    dense_output=True)
    g, f = sol.y[0, -1], sol.y[1, -1]
    return f + g, sol      # MIT condition f(R) = -g(R)

def eigenvalues(mR, mI, R, Emax=6.5):
    Es = np.linspace(mR + 1e-3, Emax, 240)
    vals = np.array([shoot(E, mR, mI, R)[0] for E in Es])
    roots = []
    for i in range(len(Es) - 1):
        if vals[i] * vals[i + 1] < 0:
            roots.append(brentq(lambda E: shoot(E, mR, mI, R)[0], Es[i], Es[i + 1],
                                xtol=1e-12))
    return np.array(roots)

def main():
    m, R = 1.7, 1.0
    print("=== Dressing Theorem: dressed system vs standard bag with m cos(delta) ===")
    worst = 0.0
    for d in (0.3, 0.6, 1.0):
        mR, mI = m * np.cos(d), m * np.sin(d)
        E_dressed = eigenvalues(mR, mI, R)
        E_standard = eigenvalues(mR, 0.0, R)
        n = min(len(E_dressed), len(E_standard))
        diff = np.max(np.abs(E_dressed[:n] - E_standard[:n]))
        worst = max(worst, diff)
        print(f"  delta={d:.1f}: levels {n}, max |E_dressed - E_std(m cos d)| = {diff:.2e}")
    print(f"  worst over deltas: {worst:.2e}")

    print("=== wavefunction ratio g_delta/g_0 vs e^(m_I r) ===")
    d = 0.6
    mR, mI = m * np.cos(d), m * np.sin(d)
    E1 = eigenvalues(mR, mI, R)[0]
    _, sol_d = shoot(E1, mR, mI, R)
    _, sol_0 = shoot(E1, mR, 0.0, R)
    rs = np.linspace(0.05, R, 60)
    ratio = sol_d.sol(rs)[0] / sol_0.sol(rs)[0]
    ratio /= ratio[0] / np.exp(mI * rs[0])      # fix overall normalization
    err = np.max(np.abs(ratio - np.exp(mI * rs)))
    print(f"  max |g_d/g_0 - e^(m_I r)| = {err:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    ax = axes[0]
    for d, c in ((0.3, "#1f77b4"), (0.6, "#d95f00"), (1.0, "#2ca02c")):
        mR, mI = m * np.cos(d), m * np.sin(d)
        Ed = eigenvalues(mR, mI, R)
        Es = eigenvalues(mR, 0.0, R)
        n = min(len(Ed), len(Es))
        ax.semilogy(Ed[:n], np.abs(Ed[:n] - Es[:n]) + 1e-16, "o-", ms=3.5, color=c,
                    label=fr"$\delta = {d}$")
    ax.set_xlabel(r"$E$"), ax.set_ylabel(r"$|E_{\rm dressed} - E_{m\cos\delta}|$")
    ax.set_title("eigenvalue identity (solver-tolerance floor)")
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.plot(rs, ratio, "o", ms=3.0, label=r"$g_\delta/g_0$ (numerical)")
    ax.plot(rs, np.exp(mI * rs), "-", lw=1.2, label=r"$e^{m_I r}$")
    ax.set_xlabel(r"$r$"), ax.set_title(r"the dressing, pointwise ($\delta = 0.6$)")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch10_fig2_radial_similarity.png")
    print("figure -> ch10_fig2_radial_similarity.png")

if __name__ == "__main__":
    main()
