#!/usr/bin/env python3
"""
ch09_eta_toy.py — spectral asymmetry on solvable spectra (Ch. 9, Figs. 9.1–9.2).

Toy A: massless MIT spectrum (mirror) -> eta(t) = 0 identically.
Toy B: same with the lowest negative level moved across zero -> eta -> 2 (Q_vac -> -1).
Regulator independence: heat kernel vs Gaussian vs smoothed sharp cutoff.
Corollary 9.2 illustration: a parametric family dragging one level through zero.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

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
    else:  # smoothed sharp cutoff
        w = 0.5 * (1 - np.tanh((a - 1.0 / t) * t * 4))
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
    for kind in ("heat", "gauss", "sharp"):
        vals = [eta_reg(shifted, t, kind) for t in np.geomspace(0.002, 0.02, 5)]
        # Richardson-flavored linear extrapolation in t
        t5 = np.geomspace(0.002, 0.02, 5)
        coef = np.polyfit(t5, vals, 1)
        print(f"  {kind:6s}: eta(0) ~ {coef[1]:+.6f}  (exact: +2)")

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6))
    ax = axes[0]
    for kind, lbl in (("heat", "heat kernel"), ("gauss", "Gaussian"), ("sharp", "smoothed cutoff")):
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
