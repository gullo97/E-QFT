"""
Sim 1 - a 1D chain of cells with shared walls: where can an iso-energy step happen?

Three gluing/identification rules for the step 'active particle n -> n+1 in cell i,
L_i -> L_i (n+1)/n' (Part I rule, applied locally):
  A  comoving identification: neighbouring cells are relabelled, not deformed (free).
  B  physical-space sudden embedding, open chain: all cells on one side are rigidly
     translated by delta; each translated occupant keeps the surviving norm
     N(x) = 1 - x + sin(2 pi n x)/(2 pi n), x = delta/L_j (identical to Thm 6.3).
  C  physical-space sudden embedding, closed ring (fixed total length): the wall is
     pushed into one neighbour, which contracts by delta; its occupants pay N(x).
Short-time golden-rule weights (Ch. 6.5): w = O^2 * product of spectator norms,
summed over the two push directions. Stochastic (rejection) unravelling.
"""
import numpy as np
from vstyle import *

rng = np.random.default_rng(7)

def surv(n, x):
    """surviving norm of mode n after cutting a fraction x off one side (sudden)."""
    x = np.minimum(x, 1.0)
    return 1 - x + np.sin(2*np.pi*n*x)/(2*np.pi*n)

# --- analytic check of the translation norm against quadrature
L = 1.0
for n, d in [(3, 0.2), (7, 0.13), (20, 0.05)]:
    y = np.linspace(d, L, 200001)
    q = np.trapz(2/L*np.sin(n*np.pi*y/L)**2, y)
    print(f"translation norm n={n} delta={d}: quadrature {q:.8f}  formula {surv(n, d/L):.8f}")

M = 40
occ = rng.poisson(5, M); occ[[6, 7, 18, 29, 30, 31]] = 0          # a few empty cells
modes0 = [list(rng.integers(8, 31, k)) for k in occ]
print("occupancy:", occ.tolist(), " total particles", occ.sum())

def run(rule, proposals=250_000, snaps=200):
    Ls = np.ones(M); modes = [list(m) for m in modes0]
    owners = np.repeat(np.arange(M), occ); slot = np.concatenate([np.arange(k) for k in occ])
    hist, where = [], np.zeros(M)
    for p in range(proposals):
        if p % (proposals // snaps) == 0: hist.append(Ls.copy())
        k = rng.integers(len(owners)); i, a = owners[k], slot[k]
        n = modes[i][a]; s = (n+1)/n; d = Ls[i]/n
        w_act = n/(n+1)
        if rule == "A":
            opts = [(1.0, None)]
        elif rule == "B":
            PR = np.prod([surv(np.array(modes[j]), d/Ls[j]).prod() if modes[j] else 1.0 for j in range(i+1, M)]) if i < M-1 else 1.0
            PL = np.prod([surv(np.array(modes[j]), d/Ls[j]).prod() if modes[j] else 1.0 for j in range(0, i)]) if i > 0 else 1.0
            opts = [(PR, "R"), (PL, "L")]
        else:  # C: ring, contract a neighbour
            opts = []
            for side, j in (("R", (i+1) % M), ("L", (i-1) % M)):
                if d >= Ls[j]: opts.append((0.0, (side, j))); continue
                P = surv(np.array(modes[j]), d/Ls[j]).prod() if modes[j] else 1.0
                opts.append((P, (side, j)))
        tot = sum(o[0] for o in opts) / len(opts)
        if rng.random() >= w_act*tot: continue
        # accepted: choose branch
        pr = np.array([o[0] for o in opts]); choice = opts[rng.choice(len(opts), p=pr/pr.sum())][1]
        where[i] += 1
        # active cell expands; its spectators keep their wavelength (dominant component)
        modes[i] = [max(1, int(round(m*s))) for m in modes[i]]; modes[i][a] = n+1
        Ls[i] *= s
        if rule == "C":
            side, j = choice; f = (Ls[j]-d)/Ls[j]
            modes[j] = [max(1, int(round(m*f))) for m in modes[j]]; Ls[j] -= d
    hist.append(Ls.copy())
    return np.array(hist), where

if __name__ == "__main__":
    res = {}
    for r in "ABC":
        H, W = run(r); res[r] = (H, W)
        print(f"rule {r}: accepted steps {int(W.sum())}; total length {H[-1].sum():.2f}; "
              f"min/max cell {H[-1].min():.3g}/{H[-1].max():.3g}; "
              f"frac of steps in outer 3 cells each side {W[[0,1,2,M-3,M-2,M-1]].sum()/W.sum():.3f}")
    np.savez("figs/sim1.npz", occ=occ, **{f"H{r}": res[r][0] for r in "ABC"}, **{f"W{r}": res[r][1] for r in "ABC"})
