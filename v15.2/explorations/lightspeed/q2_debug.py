import numpy as np
from q2_walls import simulate
rng = np.random.default_rng(3)
Nc, k = 3, 2
Xw0 = np.arange(Nc+1, dtype=float); cell = np.repeat(np.arange(Nc, dtype=np.int64), k)
xp0 = cell + 0.02 + 0.96*rng.random(cell.size); vp0 = rng.choice([-1.0, 1.0], cell.size)
for T in (0.5, 1, 2, 4, 8):
    snaps, x, v, eps, Xw, Vw, Pw, Ew, tp, tw = simulate(1, 5.0, Xw0, xp0, vp0, cell, k, T, np.array([T]), -1, 0.0)
    Xnow = snaps[0]
    print(f"T={T}: walls {np.round(Xnow,3)}  photon E sum {eps.sum():.4f}  wall E-M sum {(Ew[1:-1]-5).sum():.4f}  total-initial {eps.sum()+(Ew[1:-1]-5).sum()-2*Nc:+.2e}")
    xnow = x + v*(T-tp)
    print("   photons", np.round(xnow,3), "cells", cell)
