#!/usr/bin/env python3
"""
ch14_eta_assembly.py — eta_B assembly and budget (Ch. 14, eqs. 14.8–14.9, Fig. 14.3).

Chain: per-bag pumped QUARK number -> baryon number (1/3 per quark) -> density
per T^3 (n_bag ~ T^3 [Postulate]) -> n_B/s (45/2pi^2 g*) -> sphaleron exchange
(28/79) -> today's photon normalization ((s/n_gamma)_0 = pi^4 g_s0 / 45 zeta3).

The eps_CP reference ladder is built from Ch. 13.5's budget lines
(J_CP ~ 3e-5, S_GIM ~ 1e-7, 1/16pi^2 per loop), each rung a NAMED assumption —
no rung is tuned to land on the observation.
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
    gs_today = 43 / 11                 # 2 + (7/8)*6*(4/11): photons + 3 nu species after e+e-
    zeta3 = 1.2020569031595943
    s_over_ngamma = np.pi**4 * gs_today / (45 * zeta3)   # today: n_B/s -> n_B/n_gamma
    dilution = 45 / (2 * np.pi**2 * g_star)              # production: n_B/T^3 -> n_B/s
    N_dof = 3.0                                          # colors: three quarks pumped per bag crossing
    B_per_quark = 1.0 / 3.0                              # eta_B counts BARYONS: n_B = n_q/3
    f_neq = 1.0
    pref = s_over_ngamma * dilution * C_sph * N_dof * B_per_quark * f_neq
    print(f"today conversion      (s/n_g)_0  = pi^4 g_s0/45 zeta3 = {s_over_ngamma:.3f}   (g_s0 = 43/11)")
    print(f"entropy conversion    45/2pi^2g* = {dilution:.5f}   (g* = {g_star}; n_B/T^3 -> n_B/s)")
    print(f"sphaleron conversion  C_sph     = 28/79 = {C_sph:.4f}")
    print(f"dof factor            N_dof     = {N_dof:.0f}   (colors; per-bag -> per-T^3 assumes n_bag ~ T^3 [Postulate])")
    print(f"baryon per quark      1/N_c     = {B_per_quark:.4f}   (three pumped quarks = one baryon; N_dof/N_c = 1)")
    print(f"prefactor product               = {pref:.3e}")
    print(f"=> eta_B = {pref:.2e} * P_x * eps_CP        (eq. 14.8: 5.3e-2 * ...)")
    obs = 6.143e-10
    req = obs / pref
    print(f"observed eta_B = {obs:.3e}  =>  required P_x * eps_CP = {req:.2e}  (eq. 14.9)")
    for Px in (1.0, 0.38, 0.15, 0.07):
        print(f"   at P_x = {Px:.2f}: required eps_CP = {req/Px:.2e}")

    # eps_CP reference ladder from Ch. 13.5's budget lines (assumption-labeled)
    J_CP = 3e-5
    S_GIM = 1e-7
    loop = 1 / (16 * np.pi**2)
    ladder = [
        ("GIM evaded, one loop (ceiling): J/16pi^2", J_CP * loop),
        ("GIM evaded, two loops: J/(16pi^2)^2", J_CP * loop**2),
        ("bulk-perturbative: J*S_GIM*(1/16pi^2)^2", J_CP * S_GIM * loop**2),
    ]
    print("eps_CP reference ladder (Ch. 13.5 budget lines):")
    for name, eps in ladder:
        print(f"   {name:45s} = {eps:.1e}   required/this (P_x=0.38) = {req/0.38/eps:.1e}")

    # budget bars (Fig. 14.3): each bar = a named assumption combo, none tuned
    scenarios = [
        ("perturbative pessimum\n(sub-threshold: P=0)", 0.0, 0.0, "#7f7f7f"),
        ("strong walls, bulk-perturbative bias\n(P=0.38, eps=J*S_GIM/(16pi^2)^2)", 0.38, J_CP * S_GIM * loop**2, "#1f77b4"),
        ("strong walls, GIM evaded at two loops\n(P=0.38, eps=J/(16pi^2)^2)", 0.38, J_CP * loop**2, "#2ca02c"),
        ("strong walls, GIM evaded at one loop\n(P=0.38, eps=J/16pi^2)", 0.38, J_CP * loop, "#2ca02c"),
        ("maximal\n(P=1, eps=1)", 1.0, 1.0, "#d95f00"),
    ]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ys, labels = [], []
    for i, (name, Px, eps, color) in enumerate(scenarios):
        val = pref * Px * eps
        ys.append(i)
        labels.append(name)
        if val > 0:
            ax.barh(i, np.log10(val) + 25, left=-25, color=color, alpha=0.8, height=0.6)
            ax.text(np.log10(val) + 0.2, i, f"{val:.0e}", va="center", fontsize=7)
            print(f"   bar: {name.splitlines()[0]:40s} eta_B = {val:.1e}")
        else:
            ax.text(-24.8, i, "exactly 0 (no crossing)", va="center", fontsize=7, color="#7f7f7f")
    ax.axvline(np.log10(obs), color="k", lw=1.4, ls="--")
    ax.text(np.log10(obs) + 0.15, len(scenarios) - 0.4, "observed", fontsize=8)
    ax.set_yticks(ys, labels, fontsize=7)
    ax.set_xlabel(r"$\log_{10}\eta_B$")
    ax.set_xlim(-25, 0)
    ax.set_title("The honest budget: assumption-labeled rungs, none tuned to the answer")
    fig.savefig(FIGS / "ch14_fig3_eta_budget.png")
    print("figure -> ch14_fig3_eta_budget.png")

if __name__ == "__main__":
    main()
