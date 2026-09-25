import numpy as np, time
from concurrent.futures import ProcessPoolExecutor
from q2_walls import simulate

Nc, k, kick = 301, 16, 1.0
mus = [0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0]

def pred(kind, mu):
    return np.sqrt(3/(1+mu)) if kind == 0 else np.sqrt(1/(1+mu))

def one(args):
    kind, M, T, seed = args
    rng = np.random.default_rng(seed)
    Xw0 = np.arange(Nc+1, dtype=float); cell = np.repeat(np.arange(Nc, dtype=np.int64), k)
    xp0 = cell + 0.02 + 0.96*rng.random(cell.size); vp0 = rng.choice([-1.0, 1.0], cell.size)
    ts = np.linspace(0, T, 201)
    a = simulate(kind, M, Xw0, xp0, vp0, cell, k, T, ts, -1, 0.0)[0]
    b = simulate(kind, M, Xw0, xp0, vp0, cell, k, T, ts, Nc//2, kick)[0]
    return b - a

if __name__ == "__main__":
    res = {}
    R = 240
    with ProcessPoolExecutor(max_workers=30) as ex:
        for kind in (0, 1):
            for mu in mus:
                t0 = time.time()
                M = mu*k if kind == 0 else mu*2*k
                T = float(min(220.0, 110.0/pred(kind, mu)))
                if kind == 0 and mu < 0.3: T = 40.0          # light walls: fast signals, keep inside chain
                diffs = np.array(list(ex.map(one, [(kind, M, T, 1000*kind + s) for s in range(R)])))
                c0 = Nc//2
                right = diffs[:, :, c0+1:]; left = -diffs[:, :, c0::-1][:, :, :right.shape[2]]
                sym = np.concatenate([right, left], axis=0)                 # 2R samples, j = 1..150
                mean = sym.mean(0); sem = sym.std(0)/np.sqrt(sym.shape[0])
                nz = np.abs(sym) > 1e-9
                first = np.where(nz.any(1), np.argmax(nz, axis=1), -1)      # per sample per wall: first snapshot index
                ts = np.linspace(0, T, 201)
                arrive = np.array([ts[first[:, j][first[:, j] >= 0].min()] if (first[:, j] >= 0).any() else np.nan for j in range(first.shape[1])])
                res[(kind, mu)] = dict(T=T, mean=mean, sem=sem, arrive=arrive, M=M)
                print(f"kind={kind} mu={mu}: M={M}, T={T:.0f}, {time.time()-t0:.0f}s, plateau j=1..3 {mean[-1,:3].round(3)}, sem {sem[-1,:20].mean():.3f}", flush=True)
    np.save("figs/q2_production.npy", res, allow_pickle=True)
