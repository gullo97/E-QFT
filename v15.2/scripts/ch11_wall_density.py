#!/usr/bin/env python3
"""
ch11_wall_density.py — vacuum charge density wall layers (Ch. 11, Fig. 11.1).

rho_vac(x) = -(1/2) sum_k sgn(E_k) |psi_k(x)|^2 e^{-t|E_k|}, complete analytic
spectra of the bulk-phase MIT bag (exact roots of tan(pL) = -p/(m cos d), so
the +-E pairing -- and hence the vanishing of every whole-box vacuum charge --
holds BY CONSTRUCTION here; the honest solver-level mirror test is
ch10_mirror_spectrum.py). Shows the equal-and-opposite wall layers and prints
the new vacuum's charge inside the old region, heat-kernel regulated and
Richardson-extrapolated in t: the Polarization Identity value -DQ_cross(inf).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from ch11_polarization import modes, _trapz

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
try:
    plt.style.use(HERE / "style.mplstyle")
except OSError:
    pass


def rho_vac(L, delta, m=1.0, Emax=120.0, t=0.06, nx=401):
    x = np.linspace(0, L, nx)
    E, psi = modes(m, delta, L, Emax, x)
    dens = (np.abs(psi) ** 2).sum(axis=1)                      # (nmode, nx)
    rho = -0.5 * ((np.sign(E) * np.exp(-t * np.abs(E)))[:, None] * dens).sum(axis=0)
    return x, rho


def regional_integral(L, Lp, delta, m=1.0, Emax=300.0, ts=(0.04, 0.02, 0.01),
                      nx=3001):
    """int_0^L rho'_vac dx, heat-kernel regulated, Richardson-extrapolated."""
    x = np.linspace(0, Lp, nx)
    E, psi = modes(m, delta, Lp, Emax, x)
    mask = x <= L + 1e-12
    w = np.array([_trapz((np.abs(p) ** 2).sum(axis=0)[mask], x[mask]) for p in psi])
    vals = [float(np.sum(-0.5 * np.sign(E) * np.exp(-t * np.abs(E)) * w)) for t in ts]
    t1, t2 = ts[-2], ts[-1]
    return vals, (vals[-1] * t1 - vals[-2] * t2) / (t1 - t2)


def main():
    m, delta = 1.0, np.pi / 4
    L, Lp = 1.0, 1.2

    xo, rho_o = rho_vac(L, delta)
    xn, rho_n = rho_vac(Lp, delta)
    print(f"old box  L={L}:  net vacuum charge = {_trapz(rho_o, xo):+.2e}  "
          "(0 by construction: exact +-E pairing)")
    print(f"new box  L'={Lp}: net vacuum charge = {_trapz(rho_n, xn):+.2e}")
    vals, extrap = regional_integral(L, Lp, delta)
    print(f"new vacuum charge inside old region [0, {L}], heat kernel "
          f"t=(0.04, 0.02, 0.01): {np.round(vals, 5)}")
    print(f"  -> Richardson t->0: {extrap:+.5f}  =  -DQ_cross(inf)")
    print("  (the Polarization Identity value +0.0672 printed by "
          "ch11_polarization.py, with the sign unwound)")

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(xo, rho_o, label=fr"old box vacuum, $L = {L}$", lw=1.6)
    ax.plot(xn, rho_n, "--", label=fr"new box vacuum, $L' = {Lp}$", lw=1.6)
    ax.axvspan(0, L, color="#dddddd", alpha=0.45, label="old region $[0, L]$")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\rho_{\mathrm{vac}}(x)$  (regulated)")
    ax.set_title(r"Wall charge layers at $\delta = \pi/4$: equal and opposite, net zero")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch11_fig1_wall_charge_density.png")
    print(f"figure -> {FIGS / 'ch11_fig1_wall_charge_density.png'}")


if __name__ == "__main__":
    main()
