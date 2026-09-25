#!/usr/bin/env python3
"""
ch17_dictionary_check.py — the Dictionary Theorem, numerically (Ch. 17).

Moving-wall description: particle in [0, L], H = -(1/2m) d^2/dx^2.
Metric description: unit interval, H = -(1/2m) g^{xixi} d^2/dxi^2 with g_xixi = L^2.

The two operators are the same operator in two coordinate systems, so their spectra
must coincide. To make this a check that could fail, the two pictures are assembled
as INDEPENDENT finite-difference matrices on DIFFERENT grids and both are
diagonalized (eigh_tridiagonal): they can then agree only through the continuum
spectrum they share, down to the O(h^2) discretization floor.
Also: comoving-mode redshift along an L(t) trajectory (Fig. 17.2 data).
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def wall_spectrum(L, n, m=1.0, k=6):
    """Assemble -(1/2m) d^2/dx^2 on [0, L] (n interior points, h = L/(n+1)) and diagonalize."""
    h = L / (n + 1)
    diag = np.full(n, 2.0) / (2 * m * h**2)
    off = np.full(n - 1, -1.0) / (2 * m * h**2)
    E, _ = eigh_tridiagonal(diag, off, select="i", select_range=(0, k - 1))
    return E

def metric_spectrum(L, n, m=1.0, k=6):
    """Assemble -(1/2m) g^{xixi} d^2/dxi^2 on [0, 1], g^{xixi} = 1/L^2 (independent grid)."""
    h = 1.0 / (n + 1)
    diag = np.full(n, 2.0) / (2 * m * h**2 * L**2)
    off = np.full(n - 1, -1.0) / (2 * m * h**2 * L**2)
    E, _ = eigh_tridiagonal(diag, off, select="i", select_range=(0, k - 1))
    return E

def main():
    m, k = 1.0, 6
    j = np.arange(1, k + 1)
    print("=== moving-wall vs metric spectra: independent matrices, independent grids ===")
    for L in (1.0, 1.7, 3.2):
        Ea = wall_spectrum(L, n=2000)
        Eb = metric_spectrum(L, n=3000)
        Ec = np.pi**2 * j**2 / (2 * m * L**2)        # the shared continuum spectrum
        print(f"  L = {L}: max |E_wall - E_metric| / E = {np.max(np.abs(Ea - Eb)/Ec):.2e}"
              f"   (wall vs continuum {np.max(np.abs(Ea - Ec)/Ec):.2e},"
              f" metric vs continuum {np.max(np.abs(Eb - Ec)/Ec):.2e})")
    # the residual disagreement is pure O(h^2) discretization error: halve h, expect /4
    L = 1.7
    Ec = np.pi**2 * j**2 / (2 * m * L**2)
    r1 = np.max(np.abs(wall_spectrum(L, 1000) - metric_spectrum(L, 1500)) / Ec)
    r2 = np.max(np.abs(wall_spectrum(L, 2000) - metric_spectrum(L, 3000)) / Ec)
    print(f"  grid doubling shrinks the deviation by {r1/r2:.2f} (O(h^2): factor 4 expected)")
    print("  -> the two pictures agree through the continuum spectrum they share: the dictionary.")

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
