"""
Q1b - massless field (c=1) between a fixed mirror at x=0 and a mirror moving from L0=1 to
L1=(n+1)/n at constant speed u<=1, then at rest.  Exact solution by Moore's method:
phi(x,t) = f(t-x) - f(t+x);  mirror condition f(t+L(t)) = f(t-L(t)).
Initial standing wave mode n at rest:  f(s) = -1/2 sin(n pi s) on [-1,1].
Energy E(t) = \int_{t-L}^{t+L} f'(s)^2 ds.  Mode energies of the final cavity from phi, phi_t.
"""
import numpy as np

def build_f(n, u, smax, ds):
    L0, L1 = 1.0, (n+1)/n; tau = (L1-L0)/u
    s = np.arange(-1.0, smax+ds, ds)
    f = np.where(s <= 1.0, -0.5*np.sin(n*np.pi*s), np.nan)
    s_switch = L0 + (1+u)*tau                     # characteristic hitting the mirror when it stops
    for i in np.where(s > 1.0)[0]:
        si = s[i]
        if si <= s_switch:
            t = (si-L0)/(1+u); arg = (1-u)*t - L0
        else:
            arg = si - 2*L1
        f[i] = np.interp(arg, s, f)                # arg < si, already known
    return s, f, tau, L1

def analyse(n, u, ds=2e-5):
    s, f, tau, L1 = build_f(n, u, (1.0/n)/u + (n+1)/n + 1.0, ds)
    fp = np.gradient(f, ds)
    E0 = np.trapz(fp[(s >= -1) & (s <= 1)]**2, s[(s >= -1) & (s <= 1)])
    w = (s >= tau-L1) & (s <= tau+L1)
    Ef = np.trapz(fp[w]**2, s[w])
    # mode decomposition at t = tau
    x = np.linspace(0, L1, 20001)
    phi = np.interp(tau-x, s, f) - np.interp(tau+x, s, f)
    phit = np.interp(tau-x, s, fp) - np.interp(tau+x, s, fp)
    ks = np.arange(1, 4*n+40)
    modes = np.sqrt(2/L1)*np.sin(np.outer(ks, x)*np.pi/L1)
    q = np.trapz(modes*phi, x, axis=1); p = np.trapz(modes*phit, x, axis=1)
    om = ks*np.pi/L1; em = 0.5*(p**2 + om**2*q**2)
    return Ef/E0, em[n]/Ef, em.sum()/Ef

if __name__ == "__main__":
    for n in (3, 10):
        for u in (1.0, 0.999, 0.9, 0.5, 0.2):
            r = analyse(n, u)
            pred = 1 - (1-u)/(n*(1+u)) if u >= 1/(n+1) else np.nan
            print(f"n={n} u={u}: E_f/E0={r[0]:.5f} (single-reflection formula {pred:.5f}); energy fraction in mode n+1 = {r[1]:.5f} "
                  f"(sudden n/(n+1)={n/(n+1):.5f}); modal sum check {r[2]:.4f}")
