import numpy as np
from sim5_consensus import ck2
rng = np.random.default_rng(21)

def run_steps(modes, S, Gamma=1.0, K_fac=12):
    n = np.array(modes, int); L = 1.0; t = 0.0
    E0 = (n**2).sum()*np.pi**2/2
    for _ in range(S):
        w = Gamma*n/(n+1.0); R = w.sum(); t += rng.exponential(1/R)
        a = rng.choice(len(n), p=w/R); s = (n[a]+1)/n[a]
        for b in range(len(n)):
            if b == a: continue
            K = int(K_fac*n[b]*s)+50; p = ck2(n[b], s, K); p /= p.sum(); n[b] = rng.choice(K, p=p)+1
        n[a] += 1; L *= s
    return np.log(L), t, (n**2).sum()*np.pi**2/(2*L**2)/E0

out = {}
for N in (4, 16, 64):
    for S in (4, 16, 64):
        R = 60 if N*S <= 1024 else 30
        g, tt, ee = zip(*[run_steps(rng.integers(10, 31, N), S) for _ in range(R)])
        g = np.array(g); out[(N, S)] = (g.mean(), g.std(), np.mean(tt), np.array(ee))
        print(f"N={N:3d} S={S:3d}: growth G=<lnL>={g.mean():.4f}, std={g.std():.4f}, std*sqrt(G*20)={g.std()*np.sqrt(g.mean()*20):.3f}, "
              f"mean time {np.mean(tt):.4f} (x N = {np.mean(tt)*N:.3f}); E/E0 median {np.median(ee):.4f} mean {np.mean(ee):.4f} max {np.max(ee):.3f}")
np.save("figs/sim5b.npy", out, allow_pickle=True)

# long run: L(t) linear?  (physical momenta preserved on average -> dL/dt = Gamma * sum 1/(n0+1) * L0... in units L0=1)
modes = rng.integers(10, 31, 8); n0 = modes.copy()
n = modes.copy(); L = 1.0; t = 0.0; T, LL = [0.0], [1.0]
while L < 6.0:
    w = n/(n+1.0); R = w.sum(); t += rng.exponential(1/R); a = rng.choice(len(n), p=w/R); s = (n[a]+1)/n[a]
    for b in range(len(n)):
        if b == a: continue
        K = int(8*n[b]*s)+50; p = ck2(n[b], s, K); p /= p.sum(); n[b] = rng.choice(K, p=p)+1
    n[a] += 1; L *= s; T.append(t); LL.append(L)
T, LL = np.array(T), np.array(LL)
slope = np.polyfit(T, LL, 1)[0]
print(f"long run: {len(T)-1} steps, fitted dL/dt = {slope:.4f}; prediction sum_a (lambda_a/2) with lambda/2 = L0/n0: {np.sum(1.0/n0):.4f}")
print("final physical wavenumbers / initial (sorted):", np.round(np.sort(np.pi*n/L)/np.sort(np.pi*n0), 3))
np.savez("figs/sim5_long.npz", T=T, LL=LL, pred=np.sum(1.0/n0))
