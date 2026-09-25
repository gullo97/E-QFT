#!/usr/bin/env python3
"""
ch11_wall_density.py — vacuum charge density wall layers (Ch. 11, Fig. 11.1, new).

rho_vac(x) = -(1/2) sum_k sgn(E_k) |psi_k(x)|^2 e^{-t|E_k|}, complete spectra,
bulk-phase MIT bag (delta != 0). Shows the equal-and-opposite wall layers, net ~ 0,
and the regional integral over the old box after expansion (the Polarization
Identity integrand).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from ch08_bag_spectrum import find_spectrum, eigenfunction

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def rho_vac(L, delta, m=1.0, Emax=120.0, t=0.06, nx=401):
    mR, mI = m * np.cos(delta), m * np.sin(delta)
    spec = find_spectrum(mR, mI, L, 0.0, 0.0, Emax=Emax, nE=6001)
    x = np.linspace(0, L, nx)
    rho = np.zeros(nx)
    for E in spec:
        psi = eigenfunction(E, mR, mI, L, 0.0, x)
        rho += -0.5 * np.sign(E) * (np.abs(psi)**2).sum(axis=1) * np.exp(-t * abs(E))
    return x, rho, spec

def main():
    m, delta = 1.0, np.pi / 4
    L, Lp = 1.0, 1.2

    xo, rho_o, spec_o = rho_vac(L, delta)
    xn, rho_n, spec_n = rho_vac(Lp, delta)
    net_o = np.trapezoid(rho_o, xo)
    net_n = np.trapezoid(rho_n, xn)
    # new vacuum's charge inside the old region (Polarization Identity integrand)
    mask = xn <= L
    Q_old_region = np.trapezoid(rho_n[mask], xn[mask])
    print(f"old box  L={L}:  net vacuum charge = {net_o:+.4f} (exact: 0)")
    print(f"new box  L'={Lp}: net vacuum charge = {net_n:+.4f} (exact: 0)")
    print(f"new vacuum charge inside old region [0, {L}]: {Q_old_region:+.4f}")
    print(f"  -> -DeltaQ_cross(inf); thesis benchmark value of the polarization: 0.0672")
    print(f"     (regulator t=0.06, Emax=120; tighten both to converge)")

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(xo, rho_o, label=fr"old box vacuum, $L = {L}$", lw=1.6)
    ax.plot(xn, rho_n, "--", label=fr"new box vacuum, $L' = {Lp}$", lw=1.6)
    ax.axvspan(0, L, color="#dddddd", alpha=0.45, label="old region $[0, L]$")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\rho_{\mathrm{vac}}(x)$  (regulated)")
    ax.set_title(fr"Wall charge layers at $\delta = \pi/4$: equal and opposite, net zero")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch11_fig1_wall_charge_density.png")
    print(f"figure -> {FIGS / 'ch11_fig1_wall_charge_density.png'}")

if __name__ == "__main__":
    main()
