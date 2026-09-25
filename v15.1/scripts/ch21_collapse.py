#!/usr/bin/env python3
"""
ch21_collapse.py — separate-universe cell chain: growth of structure (Ch. 21).

A chain of cells, each a separate-universe FRW patch (matter-dominated), with a
Gaussian overdensity delta_i:

    addot_i = -(H0^2/2) Omega_i / a_i^2,    Omega_i = 1 + delta_i(0),

synchronized initial a and constraint-consistent adot (curvature offset).
Linear theory: delta ~ a; overdense cells turn around and collapse while the
background expands.  Prints the growth exponent, the turnaround and collapse
epochs, and the linearly extrapolated contrast at both epochs (spherical-
collapse thresholds 1.062 / 1.686 for comparison).

Adapted from scripts/ports/newton_sim.py section (C).
Writes ../data/ch21_collapse_traj.npz; Figure 21.3.
"""
import pathlib

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
DATA = HERE.parent / "data"
DATA.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

print("# cell chain with overdensity: linear growth, turnaround, collapse")
Mc = 41
ic = np.arange(Mc) - Mc // 2
delta0 = 0.08 * np.exp(-(ic / 4.0) ** 2)
Om = 1 + delta0
H0, a0 = 1.0, 0.02
amin = 0.25 * a0


def rhs(t, y):
    a, ad = y[:Mc], y[Mc:]
    acc = -0.5 * H0**2 * Om / np.maximum(a, amin) ** 2
    frozen = a <= amin
    ad = np.where(frozen, 0.0, ad)
    acc = np.where(frozen, 0.0, acc)
    return np.concatenate([ad, acc])


ad0 = H0 * np.sqrt(Om / a0 + (1 - Om))        # synchronized constraint-consistent ICs
y0 = np.concatenate([np.full(Mc, a0), ad0])
tend = 160.0
tev = np.geomspace(2e-3, tend, 420)
solc = solve_ivp(rhs, [0, tend], y0, t_eval=tev, method="LSODA", rtol=1e-9, atol=1e-12)
A = np.maximum(solc.y[:Mc].T, amin)            # (T, Mc)
abg = A[:, 0]                                  # far cell ~ background (delta ~ 1e-8)
dcen = (1 + delta0[Mc // 2]) * (abg / A[:, Mc // 2]) ** 3 - 1
kta = int(np.argmax(A[:, Mc // 2]))
hit = A[:, Mc // 2] <= amin * 1.001
kcol = int(np.argmax(hit)) if hit.any() else len(tev) - 1
print(f"  central delta(0) = {delta0[Mc//2]:.3f}; turnaround at a_bg/a0 = {abg[kta]/a0:.1f}, "
      f"collapse at a_bg/a0 = {abg[kcol]/a0:.1f}")
w = (abg / a0 > 30) & (abg / a0 < 300)         # late linear window: pure growing mode
slope = np.polyfit(np.log(abg[w]), np.log(dcen[w]), 1)[0]
print(f"  linear-growth exponent d ln delta / d ln a = {slope:.3f}   (theory: 1)")
Cfit = np.exp(np.mean(np.log(dcen[w]) - np.log(abg[w])))
print(f"  linearly extrapolated contrast delta_lin = C+ a_bg: at turnaround "
      f"{Cfit*abg[kta]:.3f} (spherical collapse: 1.062), at collapse "
      f"{Cfit*abg[kcol]:.3f} (spherical collapse: 1.686)")

np.savez_compressed(DATA / "ch21_collapse_traj.npz", t=tev, A=A, abg=abg, dcen=dcen,
                    delta0=delta0, ic=ic, amin=amin)

fig, (c1, c2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
c1.loglog(abg / abg[0], dcen, "crimson", lw=2, label=r"central cell $\delta(t)$ (simulation)")
c1.loglog(abg / abg[0], Cfit * abg, "k--", label=r"linear theory $\delta \propto a$")
c1.axhline(1.686, color="gray", ls=":", label=r"$\delta_{\rm lin}=1.686$ at collapse (theory)")
c1.axvline(abg[kta] / abg[0], color="gray", lw=0.8)
c1.text(abg[kta] / abg[0] * 1.05, 3e-2, "turnaround", rotation=90, fontsize=8, color="gray")
c1.set_xlabel(r"background expansion $a_{\rm bg}/a_0$"); c1.set_ylabel(r"density contrast $\delta$")
c1.set_title("gravitational instability from the cell chain"); c1.legend(fontsize=8)
ks = [0, kta // 2, kta, (kta + kcol) // 2, min(kcol + 5, len(tev) - 1)]
for k, col in zip(ks, plt.cm.viridis(np.linspace(0, 0.9, 5))):
    c2.semilogy(ic, A[k] / abg[k], color=col, label=f"$a_{{bg}}/a_0={abg[k]/a0:.0f}$")
c2.set_xlabel("cell index"); c2.set_ylabel(r"$a_i / a_{\rm bg}$")
c2.set_title("the well forms: local scale factor lags,\nturns around and collapses")
c2.legend(fontsize=8)
fig.tight_layout(); fig.savefig(FIGS / "ch21_fig3_collapse.png", dpi=150); plt.close(fig)
print("figure -> ch21_fig3_collapse.png; data -> ch21_collapse_traj.npz")
