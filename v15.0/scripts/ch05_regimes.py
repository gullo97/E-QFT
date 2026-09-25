#!/usr/bin/env python3
"""ch05_regimes.py — the sudden/adiabatic regime diagram (Ch. 5, Fig. 5.1)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

w = np.logspace(-2, 2, 200)
fig, ax = plt.subplots(figsize=(5.8, 4.0))
ax.fill_between(w, 1e-3 / w, 1e-2 * np.ones_like(w) * 0 + 1e-3, color="#aac8e8", alpha=0.0)
ax.loglog(w, 1 / w, "k-", lw=1.5)
ax.fill_between(w, 1e-3, 1 / w, color="#aac8e8", alpha=0.55)
ax.fill_between(w, 1 / w, 1e3, color="#f2c89b", alpha=0.55)
ax.text(0.05, 0.05, "QUENCH\n(state frozen;\nre-expand in new basis)", fontsize=9, ha="left")
ax.text(8, 30, "ADIABATIC\n(state tracks;\nquantum numbers ride)", fontsize=9, ha="left")
ax.text(1.3, 1.3, r"$\omega_n\tau = 1$", rotation=-38, fontsize=9)
ax.set_xlabel(r"mode frequency $\omega_n$")
ax.set_ylabel(r"change duration $\tau$")
ax.set_title("One change, two treatments — and the physical cutoff between them")
ax.set_xlim(1e-2, 1e2), ax.set_ylim(1e-3, 1e3)
fig.savefig(FIGS / "ch05_fig1_regimes.png")
print("figure -> ch05_fig1_regimes.png")
