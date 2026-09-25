"""
Q2 - exact event-driven dynamics of a 1D chain of cells: free walls of mass M, k occupants per
cell bouncing elastically between the two walls of their cell (occupants of one cell do not
interact with each other).  Two matter types:
  'nr'  : non-relativistic particles, mass m=1, speed |v|=1, Newtonian elastic collisions
  'ph'  : massless quanta, energy eps=1, speed c=1; walls relativistic with rest mass M;
          exact recoil:  eps' = eps (E_w - s P_w)/(E_w + s P_w + 2 eps)  (-> Doppler for M -> inf),
          wall momentum P' = P + s (eps + eps')          (exact 1D kinematics)
A realisation is run twice from identical initial data: unperturbed, and with the occupants
of the centre cell given 10% more energy.  The ensemble-mean wall displacement difference
gives the pulse (sound); the first nonzero difference gives the causal front.
"""
import numpy as np
from numba import njit

@njit(cache=True)
def hit_time(x, s, vel, X, V):
    # particle at x moving with velocity vel (sign s) ; wall at X moving with V
    rel = vel - V
    if rel == 0.0:
        return 1e300
    t = (X - x)/rel
    return t if t > 1e-13 else 1e300

@njit(cache=True)
def tree_build(vals):
    size = 1
    while size < vals.size:
        size *= 2
    tr = np.full(2*size, 1e300); idx = np.full(2*size, -1, np.int64)
    for i in range(vals.size):
        tr[size+i] = vals[i]; idx[size+i] = i
    for p in range(size-1, 0, -1):
        if tr[2*p] <= tr[2*p+1]:
            tr[p] = tr[2*p]; idx[p] = idx[2*p]
        else:
            tr[p] = tr[2*p+1]; idx[p] = idx[2*p+1]
    return tr, idx, size

@njit(cache=True)
def tree_set(tr, idx, size, i, val):
    p = size + i; tr[p] = val
    p //= 2
    while p >= 1:
        if tr[2*p] <= tr[2*p+1]:
            tr[p] = tr[2*p]; idx[p] = idx[2*p]
        else:
            tr[p] = tr[2*p+1]; idx[p] = idx[2*p+1]
        p //= 2

@njit(cache=True)
def simulate(kind, M, Xw0, xp0, vp0, cell, k, T, tsnap, kick_cell, kick):
    Nw = Xw0.size; P = xp0.size
    Xw = Xw0.copy(); Vw = np.zeros(Nw); Pw = np.zeros(Nw); Ew = np.full(Nw, M)   # relativistic wall state
    tw = np.zeros(Nw)           # time at which wall position was last updated
    x = xp0.copy(); v = vp0.copy(); tp = np.zeros(P)
    eps = np.abs(vp0).copy()    # photon energies (kind 1) ; for nr, speed magnitude
    if kind == 1:
        for i in range(P):
            eps[i] = 1.0
            v[i] = np.sign(vp0[i])
    for i in range(P):
        if cell[i] == kick_cell:
            if kind == 0:
                v[i] *= np.sqrt(1.0 + kick)
            else:
                eps[i] *= (1.0 + kick)
    nxt = np.empty(P); tgt = np.empty(P, np.int64)
    def_dummy = 0
    for i in range(P):
        c = cell[i]
        best = 1e300; w = -1
        for j in range(c, c+2):
            th = hit_time(x[i], 1.0, v[i], Xw[j], Vw[j])
            if th < best:
                best = th; w = j
        nxt[i] = best; tgt[i] = w
    snaps = np.empty((tsnap.size, Nw)); si = 0
    tr, tix, tsize = tree_build(nxt)
    t = 0.0
    while True:
        i = tix[1]; tn = tr[1]
        while si < tsnap.size and tsnap[si] <= tn:
            for j in range(Nw):
                snaps[si, j] = Xw[j] + Vw[j]*(tsnap[si]-tw[j])
            si += 1
        if tn > T or si >= tsnap.size:
            break
        j = tgt[i]
        # advance particle and wall to tn
        x[i] += v[i]*(tn - tp[i]); tp[i] = tn
        Xw[j] += Vw[j]*(tn - tw[j]); tw[j] = tn
        x[i] = Xw[j]
        if j == 0 or j == Nw-1:
            v[i] = -v[i]        # outer walls fixed and infinitely heavy
        elif kind == 0:
            m = 1.0
            vi, Vj = v[i], Vw[j]
            v[i] = ((m-M)*vi + 2*M*Vj)/(m+M)
            Vw[j] = ((M-m)*Vj + 2*m*vi)/(m+M)
        else:
            s = v[i]; b = Vw[j]
            e2 = eps[i]*(Ew[j] - s*Pw[j])/(Ew[j] + s*Pw[j] + 2.0*eps[i])
            Pw[j] += s*(eps[i] + e2); Ew[j] += eps[i] - e2
            eps[i] = e2; v[i] = -s
            Vw[j] = Pw[j]/Ew[j]
        # recompute next events for particles in the two cells adjacent to wall j
        for c in range(j-1, j+1):
            if c < 0 or c >= Nw-1:
                continue
            for q in range(c*k, c*k+k):
                xq = x[q] + v[q]*(tn - tp[q]); tp[q] = tn; x[q] = xq
                best = 1e300; w = -1
                for jj in range(c, c+2):
                    Xj = Xw[jj] + Vw[jj]*(tn - tw[jj])
                    th = hit_time(xq, 1.0, v[q], Xj, Vw[jj])
                    if th < best:
                        best = th; w = jj
                nxt[q] = tn + best; tgt[q] = w
                tree_set(tr, tix, tsize, q, nxt[q])
    return snaps, x, v, eps, Xw, Vw, Pw, Ew, tp, tw

def ensemble(kind, M, Nc=301, k=4, T=80.0, R=40, kick=0.1, seed=0):
    rng = np.random.default_rng(seed)
    Xw0 = np.arange(Nc+1, dtype=float)
    tsnap = np.linspace(0, T, 161)
    c0 = Nc//2
    diffs = np.zeros((R, tsnap.size, Nc+1))
    for r in range(R):
        cell = np.repeat(np.arange(Nc, dtype=np.int64), k)
        xp0 = cell + 0.02 + 0.96*rng.random(cell.size)
        vp0 = rng.choice([-1.0, 1.0], cell.size)
        a = simulate(kind, M, Xw0, xp0, vp0, cell, k, T, tsnap, -1, 0.0)[0]
        b = simulate(kind, M, Xw0, xp0, vp0, cell, k, T, tsnap, c0, kick)[0]
        diffs[r] = b - a
    return tsnap, diffs, c0

if __name__ == "__main__":
    import time
    t0 = time.time()
    ts, d, c0 = ensemble(0, 5.0, Nc=101, T=20.0, R=2)
    print("smoke test nr ok", d.shape, f"{time.time()-t0:.1f}s", "max |diff| at end", np.abs(d[:, -1]).max())
    ts, d, c0 = ensemble(1, 5.0, Nc=101, T=20.0, R=2)
    print("smoke test ph ok", f"{time.time()-t0:.1f}s", "max |diff| at end", np.abs(d[:, -1]).max())
