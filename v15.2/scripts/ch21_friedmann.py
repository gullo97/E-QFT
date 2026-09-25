#!/usr/bin/env python3
"""
ch21_friedmann.py — emergent Friedmann expansion from the constrained ladder (Ch. 21).

Generator-coupled hopping t(n) = g (n+1/2) sqrt(n/(n+1)) gives alphadot = 2g
(Ch. 20).  Closing the Hamiltonian constraint 3 kappa H^2 = rho slaves the
coupling to the matter density:

    g(alpha) = g0 exp(-3(1+w)(alpha - alpha0)/2),    rho ~ a^{-3(1+w)}.

Prediction: a(t) = (1 + 3(1+w) g0 t)^{2/(3(1+w))} -> t^{2/3} (matter),
t^{1/2} (radiation).  The genuine quantum wavepacket is evolved with the
self-consistently updated hopping (240 constraint updates over the run).

Adapted from scripts/ports/newton_sim.py section (B).  Figure 21.2.
"""
import pathlib

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ladder import gaussian_packet

HERE = pathlib.Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

print("# constraint-closed quantum lattice -> Friedmann expansion")


def friedmann_lattice(w, g0=0.5, n0=420, N=3600, T=12.0, nstep=240):
    n = np.arange(1, N + 1)
    psi = gaussian_packet(n, n0, 0.05 * n0, -np.pi / 2)
    alpha0 = np.log(n0)
    ts = np.linspace(0, T, nstep + 1)
    alphas = [alpha0]
    g = g0
    for k in range(nstep):
        t_n = g * (n[:-1] + 0.5) * np.sqrt(n[:-1] / (n[:-1] + 1))
        H = diags([t_n, t_n], [-1, 1], format="csc")
        psi = expm_multiply(-1j * (ts[k + 1] - ts[k]) * H, psi)
        mu = float((np.abs(psi) ** 2) @ n)
        al = np.log(mu)
        alphas.append(al)
        g = g0 * np.exp(-1.5 * (1 + w) * (al - alpha0))   # the Hamiltonian constraint
    return ts, np.array(alphas)


fig, (b1, b2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
for w, col, lab in [(0.0, "crimson", "matter $w=0$"), (1 / 3, "steelblue", "radiation $w=1/3$")]:
    ts, al = friedmann_lattice(w)
    g0, al0 = 0.5, al[0]
    al_an = al0 + 2 / (3 * (1 + w)) * np.log(1 + 3 * (1 + w) * g0 * ts)
    err = np.abs(al - al_an).max()
    print(f"  w={w:.2f}: max |alpha_quantum - alpha_Friedmann| = {err:.4f} "
          f"over Delta alpha = {al[-1]-al[0]:.2f}")
    b1.plot(ts, al - al0, color=col, lw=2, label=lab + " (quantum lattice)")
    b1.plot(ts, al_an - al0, "k--", lw=1)
    teff = ts + 1 / (3 * (1 + w) * g0)
    slope = np.gradient(al, np.log(teff))
    b2.plot(ts, slope, color=col, lw=2)
    b2.axhline(2 / (3 * (1 + w)), color=col, ls=":", lw=1)
    print(f"  w={w:.2f}: late-time exponent d ln a / d ln t = {slope[-1]:.4f} "
          f"(Friedmann: {2/(3*(1+w)):.4f})")
b1.plot([], [], "k--", label="Friedmann analytic")
b1.set_xlabel("model time"); b1.set_ylabel(r"$\ln a(t) - \ln a_0$")
b1.set_title("emergent Friedmann expansion\n(quantum ladder + Hamiltonian constraint)")
b1.legend(fontsize=9)
b2.set_xlabel("model time"); b2.set_ylabel(r"local exponent $d\ln a / d\ln t$")
b2.set_title(r"expansion exponents $\to$ 2/3 (matter), 1/2 (radiation)")
fig.tight_layout(); fig.savefig(FIGS / "ch21_fig2_friedmann.png", dpi=150); plt.close(fig)

# update-interval convergence check (App. C.3): halving the constraint-update
# interval must leave alpha(t) essentially unchanged
ts240, al240 = friedmann_lattice(0.0, nstep=240)
ts480, al480 = friedmann_lattice(0.0, nstep=480)
dal = np.abs(al480[::2] - al240).max()          # 480-grid contains the 240-grid
print(f"  update-interval check (w=0): max |alpha_480 - alpha_240| = {dal:.1e}")
print("figure -> ch21_fig2_friedmann.png")
