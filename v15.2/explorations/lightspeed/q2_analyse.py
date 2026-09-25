import numpy as np
res = np.load("figs/q2_production.npy", allow_pickle=True).item()
def speeds(r):
    T, mean, sem, arrive = r["T"], r["mean"], r["sem"], r["arrive"]
    ts = np.linspace(0, T, mean.shape[0]); J = np.arange(1, mean.shape[1]+1)
    # half-arrival of the mean displacement at each wall
    th = np.full(J.size, np.nan)
    for j in range(J.size):
        prof = mean[:, j]; mx = prof.max()
        if mx > 4*sem[:, j].max() and mx > 0.02:
            th[j] = ts[np.argmax(prof > 0.5*mx)]
    ok = np.isfinite(th) & (J >= 5)
    vs = np.polyfit(th[ok], J[ok], 1)[0] if ok.sum() > 5 else np.nan
    ok2 = np.isfinite(arrive) & (J >= 5)
    vc = np.polyfit(arrive[ok2], J[ok2], 1)[0] if ok2.sum() > 5 else np.nan
    return vs, vc, th, ok.sum()
if __name__ == "__main__":
    for key in sorted(res):
        kind, mu = key; r = res[key]
        vs, vc, th, nok = speeds(r)
        p = np.sqrt(3/(1+mu)) if kind == 0 else np.sqrt(1/(1+mu))
        print(f"kind={kind} mu={mu:6}: pulse speed {vs:.3f} (adiabatic/enthalpy prediction {p:.3f}), causal front {vc:.3f}, walls used {nok}, "
              f"max mean dX {r['mean'].max():.3f}")
