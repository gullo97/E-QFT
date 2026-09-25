#!/usr/bin/env python3
"""
ch09_eta_toy.py — spectral asymmetry on solvable spectra (Ch. 9, Figs. 9.1–9.2).

Toy A: massless MIT spectrum (mirror) -> eta(t) = 0 identically.
Toy B: same with the lowest negative level moved across zero -> eta -> 2 (Q_vac -> -1).
Regulator independence: heat kernel vs Gaussian vs normalized smooth cutoff.
Corollary 9.2, jump part: a parametric family dragging one level through zero.
Corollary 9.2, smooth part: the massless two-wall chiral bag tower -> eta = Delta/pi,
  continuous in the wall mismatch, with the gap open (no crossing anywhere).
Corollary 9.2, integrality hypothesis: fixed-angle plateau flatness of Q_vac(L) for
  the two-wall bag of Ch. 12 (walls (-2,+2), m=1) on both sides of L*.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def spectrum_mirror(L=1.0, N=400):
    p = (np.arange(1, N + 1) - 0.5) * np.pi / L
    return np.concatenate([p, -p])

def eta_reg(spec, t, kind="heat"):
    a = np.abs(spec)
    if kind == "heat":
        w = np.exp(-t * a)
    elif kind == "gauss":
        w = np.exp(-(t * a)**2)
    else:  # normalized smooth cutoff: f(0) = 1 exactly, as (9.1)'s hypothesis requires
        w = (1 - np.tanh(4 * (t * a - 1))) / (1 - np.tanh(-4.0))
    return np.sum(np.sign(spec) * w)

def main():
    sp = spectrum_mirror()
    shifted = sp.copy()
    i = np.argmax(shifted[shifted < 0])      # least-negative level
    idx = np.where(shifted == shifted[shifted < 0].max())[0][0]
    shifted[idx] = -shifted[idx]

    ts = np.geomspace(0.002, 0.5, 40)
    print("=== mirror spectrum: eta(t) ===")
    print("  max |eta(t)| =", max(abs(eta_reg(sp, t)) for t in ts), "(exact: 0)")
    print("=== shifted spectrum: eta(t->0) under three regulators ===")
    for kind in ("heat", "gauss", "cutoff"):
        vals = [eta_reg(shifted, t, kind) for t in np.geomspace(0.002, 0.02, 5)]
        # Richardson-flavored linear extrapolation in t
        t5 = np.geomspace(0.002, 0.02, 5)
        coef = np.polyfit(t5, vals, 1)
        print(f"  {kind:6s}: eta(0) ~ {coef[1]:+.6f}  (exact: +2)")

    # --- Cor. 9.2, smooth part in isolation: massless two-wall chiral bag ---
    # spectrum from the master equation (12.7) at m = 0 (L = 1):
    #   E_n^+ = (n + 1/2)pi - Delta/2,   E_n^- = -[(n + 1/2)pi + Delta/2]
    print("=== massless chiral bag: eta = Delta/pi, continuous, no crossing ===")
    nlev = np.arange(200000)
    ts2 = np.geomspace(1e-4, 1e-3, 6)
    for Delta in (0.5, 1.0, 2.0, 3.0):
        Ep = (nlev + 0.5) * np.pi - Delta / 2          # positive tower
        Em = (nlev + 0.5) * np.pi + Delta / 2          # |negative tower|
        eta0 = {}
        for kind in ("heat", "gauss"):
            if kind == "heat":
                vals = [np.sum(np.exp(-t * Ep)) - np.sum(np.exp(-t * Em)) for t in ts2]
            else:
                vals = [np.sum(np.exp(-(t * Ep)**2)) - np.sum(np.exp(-(t * Em)**2)) for t in ts2]
            eta0[kind] = np.polyfit(ts2, vals, 2)[-1]
        print(f"  Delta={Delta:3.1f}: eta_heat={eta0['heat']:+.8f}  eta_gauss={eta0['gauss']:+.8f}"
              f"  Delta/pi={Delta/np.pi:+.8f}  min|E|={(np.pi - Delta)/2:.4f} > 0")

    # --- Cor. 9.2, integrality hypothesis: fixed-angle plateaus of the Ch. 12 bag ---
    # walls (-2,+2): Sigma = 0, Delta = 4, m = 1; master equation (12.7) root-solved.
    print("=== two-wall bag (-2,+2), m=1: Q_vac(L) plateau flatness ===")
    c2, s2v = np.cos(2.0), np.sin(2.0)                 # cos(Delta/2), sin(Delta/2)

    def bag_spectrum(L, Emax=800.0):
        roots = []
        for sgn in (+1.0, -1.0):                       # above-gap branches E = sgn*sqrt(p^2+1)
            def f(p, sgn=sgn):
                return p * c2 * np.cos(p * L) - np.sin(p * L) * (sgn * np.hypot(p, 1.0) * s2v - 1.0)
            pmax = np.sqrt(Emax**2 - 1.0)
            grid = np.linspace(1e-9, pmax, max(2000, int(pmax * L / 0.25) + 100))
            fg = f(grid)
            for i in np.nonzero(np.sign(fg[:-1]) * np.sign(fg[1:]) < 0)[0]:
                p = brentq(f, grid[i], grid[i + 1], xtol=1e-13)
                roots.append(sgn * np.hypot(p, 1.0))

        def g(E):                                      # below-gap branch, q = sqrt(1 - E^2)
            q = np.sqrt(1.0 - E * E)
            return q * c2 * np.cosh(q * L) - np.sinh(q * L) * (E * s2v - 1.0)
        grid = np.linspace(-1 + 1e-9, 1 - 1e-9, 4000)
        gg = np.array([g(E) for E in grid])
        for i in np.nonzero(np.sign(gg[:-1]) * np.sign(gg[1:]) < 0)[0]:
            roots.append(brentq(g, grid[i], grid[i + 1], xtol=1e-14))
        return np.array(sorted(roots))

    def qvac(L):
        E = bag_spectrum(L)
        tl = np.geomspace(0.012, 0.09, 9)
        vals = [np.sum(np.sign(E) * np.exp(-t * np.abs(E))) for t in tl]
        return -0.5 * np.polyfit(tl, vals, 2)[-1]

    below = {L: qvac(L) for L in (0.30, 0.35, 0.40)}
    above = {L: qvac(L) for L in (0.50, 0.60, 0.70)}
    for L, q in {**below, **above}.items():
        print(f"  Q_vac(L={L:.2f}) = {q:+.6f}")
    spread_b = max(below.values()) - min(below.values())
    spread_a = max(above.values()) - min(above.values())
    print(f"  L* = artanh({-c2:.5f}) = {np.arctanh(-c2):.5f}")
    print(f"  plateau spread below L*: {spread_b:.1e}; above L*: {spread_a:.1e}"
          f"  (flat: the smooth part is L-independent at fixed angles)")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    ax = axes[0]
    for kind, lbl in (("heat", "heat kernel"), ("gauss", "Gaussian"), ("cutoff", "normalized smooth cutoff")):
        ax.semilogx(ts, [eta_reg(shifted, t, kind) for t in ts], label=lbl)
    ax.axhline(2, color="k", ls=":", lw=0.8)
    ax.set_xlabel("regulator scale $t$"), ax.set_ylabel(r"$\eta(t)$")
    ax.set_title("shifted spectrum: three regulators, one limit")
    ax.legend(fontsize=8)

    # Corollary 9.2: drag one level through zero
    ax = axes[1]
    lam = np.linspace(-1, 1, 201)
    base = spectrum_mirror(N=60)
    etas = []
    for l in lam:
        spec = np.concatenate([base, [l]])   # one extra level crossing zero at l = 0
        etas.append(eta_reg(spec, 0.01))
    ax.plot(lam, etas)
    ax.set_xlabel(r"deformation $\lambda$ (level at $E = \lambda$)")
    ax.set_ylabel(r"$\eta$")
    ax.set_title(r"$\eta$ jumps by 2 exactly at the zero crossing")
    fig.savefig(FIGS / "ch09_fig2_regulator_convergence.png")

    # Fig 9.1: schematic level diagrams
    fig2, axes2 = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
    for axx, spc, ttl in ((axes2[0], sp[:12].tolist() + sp[400:412].tolist(), "mirror: $\\eta = 0$"),
                          (axes2[1], shifted[:12].tolist() + shifted[400:412].tolist(),
                           "one level across zero: $\\eta = 2$")):
        for E in spc:
            axx.hlines(E, 0.15, 0.85, color=("#1f77b4" if E > 0 else "#d95f00"), lw=1.4)
        axx.axhline(0, color="k", lw=0.7, ls="--")
        axx.set_title(ttl), axx.set_xticks([])
        axx.set_ylim(-12, 12)
    axes2[0].set_ylabel("$E$")
    fig2.savefig(FIGS / "ch09_fig1_eta_schematic.png")
    print(f"figures -> ch09_fig1_eta_schematic.png, ch09_fig2_regulator_convergence.png")

if __name__ == "__main__":
    main()
