import numpy as np
from sim1_chain_maps import surv, M, occ, modes0, run
import sim1_chain_maps as S

# t = 0 acceptance weight per cell, averaged over its occupants (rule B: push to cheaper side)
def weights0(modes, Ls):
    wB, wC = np.zeros(M), np.zeros(M)
    for i in range(M):
        if not modes[i]: wB[i] = wC[i] = np.nan; continue
        b, c = [], []
        for n in modes[i]:
            d = Ls[i]/n; O2 = n/(n+1)
            PR = np.prod([surv(np.array(modes[j]), d/Ls[j]).prod() for j in range(i+1, M) if modes[j]])
            PL = np.prod([surv(np.array(modes[j]), d/Ls[j]).prod() for j in range(0, i) if modes[j]])
            b.append(O2*(PR+PL)/2)
            cs = []
            for j in ((i+1) % M, (i-1) % M):
                cs.append(surv(np.array(modes[j]), d/Ls[j]).prod() if modes[j] else 1.0)
            c.append(O2*np.mean(cs))
        wB[i], wC[i] = np.mean(b), np.mean(c)
    return wB, wC
wB, wC = weights0(modes0, np.ones(M))
print("rule B t=0 mean weights: edges", np.round(wB[[0,1,M-2,M-1]],3), " bulk min", np.nanmin(wB), " cell 20", wB[20])
print("rule C t=0 weights: min", np.nanmin(wC), "max", np.nanmax(wC))
np.save("figs/sim1_w0.npy", np.vstack([wB, wC]))

d = np.load("figs/sim1.npz")
HC = d["HC"]; fin = HC[-1]
for k in range(M): print(k, occ[k], f"{fin[k]:.3g}", end=" | ")
print()
print("rule C: corr(log L_final, occupancy) =", np.corrcoef(np.log(fin), occ)[0,1])
big = np.argsort(fin)[-6:]; print("largest cells", big, occ[big], np.round(fin[big],2))
print("cells < 1e-2:", np.where(fin < 1e-2)[0], occ[fin < 1e-2])
