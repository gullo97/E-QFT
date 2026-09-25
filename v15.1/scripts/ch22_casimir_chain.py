#!/usr/bin/env python3
"""
ch22_casimir_chain.py — cell energetics and derived wall dynamics (Ch. 22).

The matter vacuum energy of a cell is the 1D Casimir energy E_C(L) = -pi/(24 L)
(massless Dirichlet field): walls ATTRACT, F = -dE_C/dL < 0.  An occupant in
mode n (Part I nonrelativistic spectrum, m = 1) pushes back with
E_m = pi^2 n^2 / (2 L^2).  A chain of cells relaxes under these derived forces:
empty cells collapse to the short-range core scale, occupied cells survive and
share the freed length.  Endpoint identities asserted programmatically (App. C.5).

Adapted from scripts/ports/newton_sim.py section (A).  Figure 22.2.
"""
import pathlib

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

print("# induced stiffness: vacuum (Casimir) wall attraction vs matter pressure")


def E_cell(L, n):
    return -np.pi / (24 * L) + np.pi**2 * n**2 / (2 * L**2)


M = 12
n_occ = np.array([2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0])
Ltot = 24.0
rng = np.random.default_rng(3)
L0 = np.full(M, Ltot / M) * (1 + 0.04 * rng.standard_normal(M))
L0 *= Ltot / L0.sum()
x0 = np.cumsum(L0)[:-1]                       # interior wall positions

Lfloor = 0.10                                 # short-range core = minimum cell size


def force(x):
    L = np.diff(np.concatenate([[0.0], x, [Ltot]]))
    Lc = np.maximum(L, 1e-3)
    K = np.pi / (12 * Lfloor)                 # core balancing Casimir at L ~ sqrt(2) Lfloor
    dEdL = np.pi / (24 * Lc**2) - np.pi**2 * n_occ**2 / Lc**3 - K * Lfloor**3 / Lc**4
    return dEdL[1:] - dEdL[:-1]               # F_k = E'(L_{k+1}) - E'(L_k)


sol = solve_ivp(lambda t, x: 4.0 * force(x), [0, 1200], x0, method="LSODA",
                t_eval=np.linspace(0, 1200, 400), rtol=1e-8, atol=1e-10)
Lhist = np.array([np.diff(np.concatenate([[0.0], sol.y[:, k], [Ltot]]))
                  for k in range(sol.y.shape[1])])
Lend = Lhist[-1]
Lempty = Lend[n_occ == 0]
Locc = Lend[n_occ > 0]
print(f"  empty-cell sizes at end:    {np.round(Lempty, 3)}  (collapse)")
print(f"  occupied-cell sizes at end: {np.round(Locc, 3)}  (survive & share length)")
print(f"  endpoints: empty {Lempty.mean():.3f}, occupied {Locc.mean():.3f}, "
      f"total {Lend.sum():.3f}")
print(f"  core scale check: empty endpoint / (sqrt(2) Lfloor) = "
      f"{Lempty.mean()/(np.sqrt(2)*Lfloor):.3f}  "
      f"(isolated-cell balance predicts 1; the chain's inter-cell forces shift it)")

# endpoint identities, asserted programmatically (App. C.5)
assert abs(Lend.sum() - Ltot) < 1e-6, "total length not conserved"
assert np.allclose(Lempty, Lempty.mean(), rtol=2e-2), "empty endpoints not uniform"
assert np.allclose(Locc, Locc.mean(), rtol=2e-2), "occupied endpoints not uniform"
assert abs(Lempty.mean() / (np.sqrt(2) * Lfloor) - 1) < 0.1, \
    "empty-cell endpoint far from the core balance scale"

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
Ls = np.linspace(0.3, 6, 300)
a1.plot(Ls, -np.pi / (24 * Ls), "k--", label=r"vacuum: $E_C=-\pi/24L$ (attractive)")
a1.plot(Ls, np.pi**2 * 4 / (2 * Ls**2), "b:", label=r"matter $n=2$: $\pi^2 n^2/2L^2$ (repulsive)")
a1.plot(Ls, E_cell(Ls, 2), "crimson", label="total, occupied cell")
a1.set_xlabel("$L$"); a1.set_ylabel("cell energy"); a1.set_ylim(-1.5, 3)
a1.set_title("derived cell energetics (1D Sakharov)"); a1.legend(fontsize=8)
for i in range(M):
    a2.plot(sol.t, Lhist[:, i], color=("crimson" if n_occ[i] else "0.4"),
            lw=1.6 if n_occ[i] else 1.0)
a2.set_xlabel("relaxation time"); a2.set_ylabel("$L_i(t)$")
a2.set_title("chain relaxation: empty cells collapse (gray),\noccupied cells survive (red)")
fig.tight_layout(); fig.savefig(FIGS / "ch22_fig2_cell_energetics.png", dpi=150); plt.close(fig)
print("figure -> ch22_fig2_cell_energetics.png")
