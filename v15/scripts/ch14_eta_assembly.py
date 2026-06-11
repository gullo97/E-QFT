#!/usr/bin/env python3
"""
ch14_eta_assembly.py — eta_B assembly and budget (Ch. 14, eqs. 14.8–14.9, Fig. 14.3).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def main():
    C_sph = 28 / 79
    g_star = 106.75
    dilution = 45 / (2 * np.pi**2 * g_star)
    N_dof = 3.0
    f_neq = 1.0
    pref = C_sph * dilution * N_dof * f_neq
    print(f"sphaleron conversion  C_sph     = 28/79 = {C_sph:.4f}")
    print(f"entropy dilution      45/2pi^2g* = {dilution:.5f}   (g* = {g_star})")
    print(f"dof factor            N_dof     = {N_dof}")
    print(f"prefactor product               = {pref:.3e}")
    print(f"=> eta_B = {pref:.2e} * P_x * eps_CP        (eq. 14.8: 2.3e-2 * ...)")
    obs = 6.143e-10
    req = obs / pref
    print(f"observed eta_B = {obs:.3e}  =>  required P_x * eps_CP = {req:.2e}  (eq. 14.9)")
    for Px in (0.38, 0.15, 0.07):
        print(f"   at P_x = {Px:.2f}: required eps_CP = {req/Px:.2e}")

    # budget bars (Fig. 14.3)
    scenarios = [
        ("perturbative pessimum\n(sub-threshold)", 0.0, 0.0, "#7f7f7f"),
        ("broad tails, GIM-crushed\n(P=0.07, eps=1e-10)", 0.07, 1e-10, "#1f77b4"),
        ("strong walls, GIM-crushed\n(P=0.38, eps=1e-10)", 0.38, 1e-10, "#1f77b4"),
        ("strong walls, partial evasion\n(P=0.38, eps=1e-8)", 0.38, 1e-8, "#2ca02c"),
        ("strong walls, loop-level CP\n(P=0.38, eps=1e-6)", 0.38, 1e-6, "#2ca02c"),
        ("maximal\n(P=1, eps=1)", 1.0, 1.0, "#d95f00"),
    ]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ys, labels = [], []
    for i, (name, Px, eps, color) in enumerate(scenarios):
        val = pref * Px * eps
        ys.append(i)
        labels.append(name)
        if val > 0:
            ax.barh(i, np.log10(val) + 25, left=-25, color=color, alpha=0.8, height=0.6)
            ax.text(np.log10(val) + 0.2, i, f"{val:.0e}", va="center", fontsize=7)
        else:
            ax.text(-24.8, i, "exactly 0 (no crossing)", va="center", fontsize=7, color="#7f7f7f")
    ax.axvline(np.log10(obs), color="k", lw=1.4, ls="--")
    ax.text(np.log10(obs) + 0.15, len(scenarios) - 0.4, "observed", fontsize=8)
    ax.set_yticks(ys, labels, fontsize=7)
    ax.set_xlabel(r"$\log_{10}\eta_B$")
    ax.set_xlim(-25, 0)
    ax.set_title("The honest budget: no bar is labeled 'prediction'")
    fig.savefig(FIGS / "ch14_fig3_eta_budget.png")
    print("figure -> ch14_fig3_eta_budget.png")

if __name__ == "__main__":
    main()
