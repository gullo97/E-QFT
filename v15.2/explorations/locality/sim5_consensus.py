"""
Sim 5 - consensus inside one cell with N occupants (embedding picture, thesis rules).
Gillespie: each occupant is 'active' at rate Gamma * n/(n+1) (expansion weight; spectators free, Thm 6.1).
Step: L -> L (n+1)/n, active n -> n+1; each spectator resampled from |c_k(n_b, s)|^2 (eq. 3.8/3.9).
"""
import numpy as np
rng = np.random.default_rng(5)

def ck2(n, s, K):
    k = np.arange(1, K+1, dtype=float)
    den = k**2 - (n*s)**2
    with np.errstate(divide="ignore", invalid="ignore"):
        c = 2*s**1.5/np.pi*((-1)**n)*n*np.sin(k*np.pi/s)/den
    hit = np.isclose(den, 0, atol=1e-9); c[hit] = 1/np.sqrt(s)
    return c**2

# (1) completeness and energy moments of one sudden expansion step
for n, s in ((5, 6/5), (12, 13/12), (7, 1.37)):
    out = []
    for K in (500, 2000, 8000, 32000):
        p = ck2(n, s, K); k = np.arange(1, K+1)
        E = (k/s)**2; En = n**2                         # energies in units pi^2/(2L^2)
        out.append((K, 1-p.sum(), (p*E).sum()/En, (p*E**2).sum()/En**2))
    print(f"n={n}, s={s:.4f}: " + "; ".join(f"K={K}: 1-norm {a:.1e}, <E>/E {b:.5f}, <E^2>/E^2 {c:.1f}" for K, a, b, c in out))

def run(modes, T, Gamma=1.0, K_fac=30, rec=200):
    n = np.array(modes, int); L = 1.0; t = 0.0
    k_phys0 = np.pi*n/L; E0 = (np.pi**2*n**2/(2*L**2)).sum()
    ts, lnL, Es = [0.0], [0.0], [1.0]
    while t < T:
        w = Gamma*n/(n+1.0); R = w.sum()
        t += rng.exponential(1/R)
        a = rng.choice(len(n), p=w/R); s = (n[a]+1)/n[a]
        for b in range(len(n)):
            if b == a: continue
            K = int(K_fac*n[b]*s) + 50
            p = ck2(n[b], s, K); p /= p.sum()
            n[b] = rng.choice(K, p=p) + 1
        n[a] += 1; L *= s
        ts.append(t); lnL.append(np.log(L)); Es.append((np.pi**2*n**2/(2*L**2)).sum()/E0)
    return np.array(ts), np.array(lnL), np.array(Es), n, L

if __name__ == "__main__":
    # (2) growth law and 1/sqrt(N) consensus
    res = {}
    for N in (4, 16, 64):
        R = 40 if N < 64 else 24
        T = 3.0/N                      # same expected number of steps per particle
        fin, Ef, drift_pred = [], [], []
        for r in range(R):
            modes = rng.integers(10, 31, N)
            lam_half = (1.0/modes)*(modes/(modes+1.0))
            drift_pred.append(lam_half.sum())      # dL/dt at t=0 (Gamma=1, L=1)
            ts, lnL, Es, n, L = run(modes, T)
            fin.append(lnL[-1]); Ef.append(Es[-1])
            if r == 0: res[f"traj{N}"] = (ts, lnL, Es)
        fin = np.array(fin); Ef = np.array(Ef)
        print(f"N={N}: <ln L(T)> = {fin.mean():.4f} +- {fin.std():.4f}  -> relative spread {fin.std()/fin.mean():.4f};"
              f"  E(T)/E(0): mean {Ef.mean():.4f}, median {np.median(Ef):.4f}, max {Ef.max():.3f}")
        res[N] = (fin, Ef)
    np.save("figs/sim5.npy", res, allow_pickle=True)
