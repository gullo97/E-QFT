#!/usr/bin/env python3
"""
ch17_dictionary_check.py — the Dictionary Theorem, numerically (Ch. 17).

Moving-wall description: particle in [0, L], H = -(1/2m) d^2/dx^2.
Metric description: unit interval, H = -(1/2m) g^{xixi} d^2/dxi^2 with g = L^2.
Finite-difference spectra of both must agree (they are the same operator).
Also: comoving-mode redshift along an L(t) trajectory (Fig. 17.2 data).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def fd_spectrum(L, n=2000, m=1.0, k=6):
    h = L / (n + 1)
    main = np.full(n, 2.0)
    off = np.full(n - 1, -1.0)
    # eigenvalues of tridiagonal via analytic formula for the FD Laplacian
    j = np.arange(1, k + 1)
    lam = (2 - 2 * np.cos(j * np.pi / (n + 1))) / h**2
    return lam / (2 * m)

def main():
    m = 1.0
    print("=== moving-wall vs metric spectra (finite differences, n = 2000) ===")
    for L in (1.0, 1.7, 3.2):
        Ea = fd_spectrum(L)                                # wall at L, flat metric
        # metric picture: unit interval, inverse metric 1/L^2 multiplies the FD Laplacian
        Eb = fd_spectrum(1.0) / L**2
        print(f"  L = {L}: max |E_wall - E_metric| / E = {np.max(np.abs(Ea - Eb)/Ea):.2e}")
    print("  (identical by construction: the dictionary is a change of variables)")

    # redshift figure: iso-energy curves re-read as comoving modes
    t = np.linspace(0, 3, 200)
    Lt = 1.0 + 0.8 * t                                     # linear expansion history
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for n in range(1, 7):
        ax.plot(t, n * np.pi / Lt, color="#1f77b4", lw=1.2)
    ax.set_xlabel("$t$")
    ax.set_ylabel(r"physical frequency $\omega_n = n\pi/L(t)$")
    ax.set_title("Ch. 2's iso-energy curves, re-read: comoving modes redshifting as $1/a(t)$")
    fig.savefig(FIGS / "ch17_fig2_redshift_curves.png")
    print("figure -> ch17_fig2_redshift_curves.png")

if __name__ == "__main__":
    main()
