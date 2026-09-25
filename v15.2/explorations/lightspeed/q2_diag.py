import numpy as np
from q2_walls import simulate
rng = np.random.default_rng(1)
Nc, k = 101, 4
Xw0 = np.arange(Nc+1, dtype=float)
cell = np.repeat(np.arange(Nc, dtype=np.int64), k)
xp0 = cell + 0.02 + 0.96*rng.random(cell.size); vp0 = rng.choice([-1.0, 1.0], cell.size)
ts = np.linspace(0, 20, 41)
for kind, M in ((0, 5.0), (1, 5.0), (1, 50.0)):
    s = simulate(kind, M, Xw0, xp0, vp0, cell, k, 20.0, ts, -1, 0.0)
    L = np.diff(s, axis=1)
    dX = s - Xw0
    print(kind, M, "min cell length over time", L.min().round(4), " max |wall displacement|", np.abs(dX).max().round(3),
          " rms displacement at t=20", dX[-1].std().round(4))
