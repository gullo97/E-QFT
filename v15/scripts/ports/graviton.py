"""
Milestone 4: the graviton test (covariant proper-time regulator).

Method: massless field on a periodic lattice with metric g(x) (exact parametrization,
A = g^{-1} per site).  All perturbations depend on x only -> block decomposition over
exact transverse momenta; spectral derivatives in x.  Vacuum energy with proper-time
smoothing E(s) = (1/2) sum omega exp(-s omega): built from the operator itself, hence
diffeomorphism-invariant in the continuum at every s.  The raw lattice sum (s=0) breaks
diffeo at the cutoff scale (measured: exact-diffeo configs cost Pi ~ 0.4); smoothing at
s ~ 1.5-2 lattice units restores covariance (exact-diffeo Pi drops below 1% of signals)
while remaining cutoff-dominated: this is the covariant induced action of the geometry.

Observables: Pi(q) = 4 dE / (eps^2 V) per polarization; gradient couplings kappa from
Pi = Pi0 + kappa q^2.  Verdict sought: kappa_TT+ = kappa_TTx > 0 (the graviton),
kappa_conf opposite sign (DeWitt signature), diffeo control ~ 0.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

def spectral_D(N):
    F = np.fft.fft(np.eye(N), axis=0)
    kk = 2j*np.pi*np.fft.fftfreq(N)
    D = np.real(np.fft.ifft(kk[:, None]*F, axis=0))
    return 0.5*(D - D.T)

def invert_metric(g, d, N):
    """g: dict (a,b)->array over x. Return A = g^{-1} per site and sqrt(det g)."""
    G = np.zeros((N, d, d))
    for a in range(d):
        for b in range(d):
            G[:, a, b] = g[(min(a,b), max(a,b))]
    A = np.linalg.inv(G)
    sg = np.sqrt(np.linalg.det(G))
    out = {}
    for a in range(d):
        for b in range(a, d):
            out[(a, b)] = A[:, a, b]
    return out, sg

def vac_energy(N, d, gfun, s):
    x = np.arange(N)
    A, sg = invert_metric(gfun(x), d, N)
    w = {k: A[k]*sg for k in A}
    rs = 1.0/np.sqrt(sg)
    Dx = spectral_D(N)
    kts = 2*np.pi*np.fft.fftfreq(N)
    Kxx = Dx.T @ (w[(0,0)][:, None]*Dx)
    E = 0.0
    if d == 2:
        Cxy = 1j*(Dx.T @ np.diag(w[(0,1)]) - np.diag(w[(0,1)]) @ Dx)
        for ky in kts:
            K = Kxx + np.diag(ky**2*w[(1,1)]) + ky*Cxy
            H = (K*rs[:, None])*rs[None, :]
            om = np.sqrt(np.clip(np.linalg.eigvalsh(H), 0, None))
            E += 0.5*np.sum(om*np.exp(-s*om))
    else:
        Cxy = 1j*(Dx.T @ np.diag(w[(0,1)]) - np.diag(w[(0,1)]) @ Dx)
        Cxz = 1j*(Dx.T @ np.diag(w[(0,2)]) - np.diag(w[(0,2)]) @ Dx)
        for ky in kts:
            for kz in kts:
                K = (Kxx + np.diag(ky**2*w[(1,1)] + kz**2*w[(2,2)] + 2*ky*kz*w[(1,2)])
                     + ky*Cxy + kz*Cxz)
                H = (K*rs[:, None])*rs[None, :]
                om = np.sqrt(np.clip(np.linalg.eigvalsh(H), 0, None))
                E += 0.5*np.sum(om*np.exp(-s*om))
    return E

def gmode(d, mode, eps, q):
    def f(x):
        c = np.cos(q*x); one = np.ones_like(x, float); z = 0*one
        g = {(a, b): (one.copy() if a == b else z.copy())
             for a in range(d) for b in range(a, d)}
        if mode == "conf":
            for a in range(d): g[(a, a)] = 1 + eps*c
        elif mode == "TT+":
            g[(1, 1)] = 1 + eps*c; g[(2, 2)] = 1 - eps*c
        elif mode == "TTx":
            g[(1, 2)] = eps*c
        elif mode == "plus2d":
            g[(0, 0)] = 1 + eps*c; g[(1, 1)] = 1 - eps*c
        elif mode == "diffeo":           # exact pullback of flat by x->x, y->y+(eps/q)sin(qx)
            g[(0, 0)] = 1 + (eps*c)**2; g[(0, 1)] = eps*c
        return g
    return f

eps = 0.08
results = {}
for s in (2.0, 1.5):
    # ---- 2D validation ----
    N2 = 63
    E0 = vac_energy(N2, 2, gmode(2, "base", 0, 0), s)
    r2 = {m: [] for m in ("conf", "plus2d", "diffeo")}
    for k in (1, 2, 3, 4):
        q = 2*np.pi*k/N2
        for m in r2:
            dE = vac_energy(N2, 2, gmode(2, m, eps, q), s) - E0
            r2[m].append((k, float(q*q), float(4*dE/(eps**2*N2**2))))
    print(f"# 2D s={s}: diffeo control Pi = " +
          ", ".join(f"{r2['diffeo'][i][2]:+.1e}" for i in range(4)))
    # ---- 3D test ----
    N3 = 45
    E0_3 = vac_energy(N3, 3, gmode(3, "base", 0, 0), s)
    r3 = {m: [] for m in ("conf", "TT+", "TTx", "diffeo")}
    for k in (1, 2, 3, 4):
        q = 2*np.pi*k/N3
        for m in r3:
            dE = vac_energy(N3, 3, gmode(3, m, eps, q), s) - E0_3
            r3[m].append((k, float(q*q), float(4*dE/(eps**2*N3**3))))
        rr = {m: r3[m][-1][2] for m in r3}
        print(f"  3D s={s} k={k}: conf={rr['conf']:+.6f}  TT+={rr['TT+']:+.6f}"
              f"  TTx={rr['TTx']:+.6f}  diffeo={rr['diffeo']:+.1e}")
    fits = {}
    for tag, res in (("2d", r2), ("3d", r3)):
        for m, R in res.items():
            R = np.array(R)
            kap, Pi0 = np.polyfit(R[:, 1], R[:, 2], 1)
            fits[f"{tag}:{m}"] = [float(Pi0), float(kap)]
    print(f"# s={s} gradient couplings kappa:")
    for key, (p0, kp) in fits.items():
        print(f"    {key:10s}: Pi0 = {p0:+.6f}   kappa = {kp:+.6f}")
    results[str(s)] = {"r2": r2, "r3": r3, "fits": fits}

json.dump(results, open("graviton_results.json", "w"), indent=1)

# figure at s=2.0
s = "2.0"
r3 = results[s]["r3"]; fits = results[s]["fits"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
for m, mk, col in (("conf", "o", "crimson"), ("plus2d", "s", "seagreen"), ("diffeo", "^", "steelblue")):
    R = np.array(results[s]["r2"][m])
    a1.plot(R[:, 1], R[:, 2], mk+"-", color=col, label=m)
a1.axhline(0, color="k", lw=0.8)
a1.set_xlabel("$q^2$"); a1.set_ylabel(r"$\Pi(q)$")
a1.set_title("2 spatial dims, covariant window $s=2$:\ndiffeo control at zero")
a1.legend(fontsize=8)
for m, mk, col in (("conf", "o", "crimson"), ("TT+", "s", "seagreen"),
                   ("TTx", "D", "olive"), ("diffeo", "^", "steelblue")):
    R = np.array(r3[m])
    a2.plot(R[:, 1], R[:, 2] - fits[f"3d:{m}"][0], mk+"-", color=col,
            label=fr"{m}: $\kappa={fits[f'3d:{m}'][1]:+.5f}$")
a2.axhline(0, color="k", lw=0.8)
a2.set_xlabel("$q^2$"); a2.set_ylabel(r"$\Pi(q)-\Pi_0$")
a2.set_title("3 spatial dims, $s=2$: the graviton test\n(TT gradient couplings vs conformal vs diffeo)")
a2.legend(fontsize=8)
fig.tight_layout(); fig.savefig("figs/fig14_graviton.png", dpi=150); plt.close(fig)
print("saved fig14 and graviton_results.json")
