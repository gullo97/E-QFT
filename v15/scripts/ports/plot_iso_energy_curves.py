"""
Plot the iso-energy curves in the (n, L) plane.

For a given energy E, the constraint E = n²π²ℏ²/(2mL²) defines:
    L_n = ℓ₀ · n,   where ℓ₀ = πℏ / √(2mE)

This script visualises several iso-energy curves at different energies,
highlighting the discrete lattice points (n, L_n) that the iso-energy
axiom connects via transitions n → n+1.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# ── style ────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "font.family":      "serif",
    "mathtext.fontset":  "cm",
    "font.size":         12,
    "axes.labelsize":    14,
    "axes.titlesize":    15,
    "legend.fontsize":   11,
})

# ── parameters (natural units: ℏ = 1, m = 1) ────────────────────────
hbar = 1.0
m    = 1.0

energies = [0.5, 1.0, 2.0, 4.0]          # four representative energies
colours  = ["#2274A5", "#D64933", "#3EA660", "#9B5DE5"]
labels   = [rf"$E = {E:.1f}$" for E in energies]

n_max     = 15                             # max quantum number to plot
n_discrete = np.arange(1, n_max + 1)       # discrete n values (integers ≥ 1)
n_cont     = np.linspace(0.8, n_max + 0.5, 300)  # continuous curve for visual

# ── figure ───────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5.5))

for E, col, lab in zip(energies, colours, labels):
    ell0 = np.pi * hbar / np.sqrt(2 * m * E)   # fundamental spacing

    # continuous guide line  L = ℓ₀ · n
    L_cont = ell0 * n_cont
    ax.plot(n_cont, L_cont, color=col, lw=1.8, alpha=0.45, zorder=1)

    # discrete lattice points
    L_pts = ell0 * n_discrete
    ax.scatter(n_discrete, L_pts, s=48, color=col, edgecolors="k",
               linewidths=0.5, zorder=3, label=lab)

# ── annotate an example transition arrow  n→n+1 on the E=1 curve ───
E_ref  = 1.0
ell0   = np.pi * hbar / np.sqrt(2 * m * E_ref)
n_arr  = 4
L_start = ell0 * n_arr
L_end   = ell0 * (n_arr + 1)

ax.annotate("",
            xy=(n_arr + 1, L_end), xytext=(n_arr, L_start),
            arrowprops=dict(arrowstyle="-|>", color="#D64933",
                            lw=2.0, mutation_scale=16),
            zorder=4)
ax.text(n_arr + 0.5, (L_start + L_end) / 2 + 0.35,
        r"$n \!\to\! n\!+\!1$" "\n" r"$L \!\to\! \frac{n+1}{n}L$",
        ha="center", va="bottom", fontsize=10.5, color="#D64933",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#D64933",
                  alpha=0.85))

# ── axis labels, legend, grid ────────────────────────────────────────
ax.set_xlabel(r"Quantum number $n$")
ax.set_ylabel(r"Box size $L$  (units of $\pi\hbar/\sqrt{2m}$)")
ax.set_title("Iso-energy curves in the $(n,\\, L)$ plane")
ax.set_xlim(0, n_max + 1)
ax.set_ylim(0, None)
ax.set_xticks(range(0, n_max + 2, 2))
ax.legend(loc="upper left", framealpha=0.9)
ax.grid(True, ls="--", alpha=0.3)

fig.tight_layout()

# ── save ─────────────────────────────────────────────────────────────
out = Path(__file__).resolve().parent.parent / "figures"
out.mkdir(exist_ok=True)
fig.savefig(out / "iso_energy_curves.png", dpi=200, bbox_inches="tight")
fig.savefig(out / "iso_energy_curves.pdf", bbox_inches="tight")
print(f"Saved to {out / 'iso_energy_curves.png'}")
plt.show()
