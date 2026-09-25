import numpy as np
from q1_nr_wall import run_quantum, run_classical
us = np.geomspace(0.03, 30, 25)
out = {}
for n, N in ((3, 1200), (10, 1600), (30, 2400)):
    vp = n*np.pi
    rows = []
    for r in us:
        E, P1, Pn, nm = run_quantum(n, r*vp, N=N)
        Ec = run_classical(n, r*vp, M=40000)
        rows.append((r, E, P1, Pn, Ec))
    out[n] = np.array(rows)
    a = out[n]
    for r in (0.1, 0.5, 1.0, 2.0, 5.0):
        i = np.argmin(abs(a[:, 0]-r))
        print(f"n={n} u/vp={a[i,0]:.3f}: E/E_n quantum {a[i,1]:.4f} classical {a[i,4]:.4f}; P(n+1) {a[i,2]:.4f}; P(n) {a[i,3]:.4f}")
np.save("figs/q1_scan.npy", out, allow_pickle=True)
