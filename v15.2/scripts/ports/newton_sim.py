"""
Milestone 2: from cosmology to attraction.

(A) Induced stiffness in 1D (Sakharov in the box model): the matter vacuum energy of a
    cell is the Casimir energy E_C(L) = -pi/(24 L)  (massless Dirichlet field, hbar=c=1).
    Walls therefore ATTRACT: F = -dE_C/dL = -pi/(24 L^2) < 0.  Matter (particle in mode n,
    E_m = pi^2 n^2 / (2 L^2), m=1) pushes back.  A chain of cells relaxes: empty cells
    collapse (vacuum attraction), occupied cells survive -- derived wall dynamics.

(B) Emergent Friedmann: close the Hamiltonian constraint on the quantum geometry ladder.
    Generator-coupled hopping t(n) = g (n+1/2) sqrt(n/(n+1)) gives alphadot = 2g (milestone 1).
    The constraint  3 kappa H^2 = rho  slaves the coupling to the matter density:
        g(alpha) = g0 exp(-3(1+w)(alpha - alpha0)/2),   rho ~ a^{-3(1+w)}.
    Prediction: a(t) = (1 + 3(1+w) g0 t)^{2/(3(1+w))}  -> t^{2/3} (matter), t^{1/2} (radiation).
    We evolve the genuine quantum wavepacket with the self-consistently updated hopping.

(C) Gravitational instability ("attraction"): a chain of cells, each a separate-universe
    FRW patch (matter-dominated), with a Gaussian overdensity delta_i.  Per cell:
        addot_i = -(H0^2/2) Omega_i / a_i^2,   Omega_i = 1 + delta_i(0),
    synchronized initial a and constraint-consistent adot.  Linear theory: delta ~ a;
    overdense cells turn around and collapse while the background expands.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

# ============ (A) Casimir wall dynamics: derived attraction ============
print("# (A) induced stiffness: vacuum (Casimir) wall attraction vs matter pressure")
def E_cell(L, n):
    return -np.pi/(24*L) + np.pi**2*n**2/(2*L**2)

M = 12
n_occ = np.array([2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0])
Ltot = 24.0
rng = np.random.default_rng(3)
L0 = np.full(M, Ltot/M)*(1 + 0.04*rng.standard_normal(M))
L0 *= Ltot/L0.sum()
x0 = np.cumsum(L0)[:-1]                       # interior wall positions

Lfloor = 0.10                                      # short-range core = minimum cell size (UV cutoff)
def force(x):
    L = np.diff(np.concatenate([[0.0], x, [Ltot]]))
    Lc = np.maximum(L, 1e-3)
    K = np.pi/(12*Lfloor)                          # core balancing Casimir at L ~ sqrt(2) Lfloor
    dEdL = np.pi/(24*Lc**2) - np.pi**2*n_occ**2/Lc**3 - K*Lfloor**3/Lc**4
    return dEdL[1:] - dEdL[:-1]                    # F_k = E'(L_{k+1}) - E'(L_k)

sol = solve_ivp(lambda t, x: 4.0*force(x), [0, 1200], x0, method="LSODA",
                t_eval=np.linspace(0, 1200, 400), rtol=1e-8, atol=1e-10)
Lhist = np.array([np.diff(np.concatenate([[0.0], sol.y[:, k], [Ltot]]))
                  for k in range(sol.y.shape[1])])
print(f"  empty-cell sizes at end:    {np.round(Lhist[-1][n_occ==0], 3)}  (collapse)")
print(f"  occupied-cell sizes at end: {np.round(Lhist[-1][n_occ>0], 3)}  (survive & share length)")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
Ls = np.linspace(0.3, 6, 300)
a1.plot(Ls, -np.pi/(24*Ls), "k--", label=r"vacuum: $E_C=-\pi/24L$ (attractive)")
a1.plot(Ls, np.pi**2*4/(2*Ls**2), "b:", label=r"matter $n=2$: $\pi^2 n^2/2L^2$ (repulsive)")
a1.plot(Ls, E_cell(Ls, 2), "crimson", label="total, occupied cell")
a1.set_xlabel("$L$"); a1.set_ylabel("cell energy"); a1.set_ylim(-1.5, 3)
a1.set_title("derived cell energetics (1D Sakharov)"); a1.legend(fontsize=8)
for i in range(M):
    a2.plot(sol.t, Lhist[:, i], color=("crimson" if n_occ[i] else "0.4"),
            lw=1.6 if n_occ[i] else 1.0)
a2.set_xlabel("relaxation time"); a2.set_ylabel("$L_i(t)$")
a2.set_title("chain relaxation: empty cells collapse (gray),\noccupied cells survive (red)")
fig.tight_layout(); fig.savefig("figs/fig8_stiffness.png", dpi=150); plt.close(fig)

# ============ (B) emergent Friedmann from the constrained quantum ladder ============
print("# (B) constraint-closed quantum lattice -> Friedmann expansion")
def friedmann_lattice(w, g0=0.5, n0=420, N=3600, T=12.0, nstep=240):
    n = np.arange(1, N+1)
    psi = np.exp(-(n - n0)**2/(4*(0.05*n0)**2) - 1j*(np.pi/2)*n)
    psi /= np.linalg.norm(psi)
    alpha0 = np.log(n0)
    ts = np.linspace(0, T, nstep+1)
    alphas = [alpha0]
    g = g0
    for k in range(nstep):
        t_n = g*(n[:-1] + 0.5)*np.sqrt(n[:-1]/(n[:-1]+1))
        H = diags([t_n, t_n], [-1, 1], format="csc")
        psi = expm_multiply(-1j*(ts[k+1]-ts[k])*H, psi)
        mu = float((np.abs(psi)**2) @ n)
        al = np.log(mu)
        alphas.append(al)
        g = g0*np.exp(-1.5*(1+w)*(al - alpha0))      # the Hamiltonian constraint
    return ts, np.array(alphas)

fig, (b1, b2) = plt.subplots(1, 2, figsize=(10.5, 4.2))
for w, col, lab in [(0.0, "crimson", "matter $w=0$"), (1/3, "steelblue", "radiation $w=1/3$")]:
    ts, al = friedmann_lattice(w)
    g0, al0 = 0.5, al[0]
    al_an = al0 + 2/(3*(1+w))*np.log(1 + 3*(1+w)*g0*ts)
    err = np.abs(al - al_an).max()
    print(f"  w={w:.2f}: max |alpha_quantum - alpha_Friedmann| = {err:.4f} over Delta alpha = {al[-1]-al[0]:.2f}")
    b1.plot(ts, al - al0, color=col, lw=2, label=lab+" (quantum lattice)")
    b1.plot(ts, al_an - al0, "k--", lw=1)
    teff = ts + 1/(3*(1+w)*g0)
    slope = np.gradient(al, np.log(teff))
    b2.plot(ts, slope, color=col, lw=2)
    b2.axhline(2/(3*(1+w)), color=col, ls=":", lw=1)
b1.plot([], [], "k--", label="Friedmann analytic")
b1.set_xlabel("model time"); b1.set_ylabel(r"$\ln a(t) - \ln a_0$")
b1.set_title("emergent Friedmann expansion\n(quantum ladder + Hamiltonian constraint)")
b1.legend(fontsize=9)
b2.set_xlabel("model time"); b2.set_ylabel(r"local exponent $d\ln a / d\ln t$")
b2.set_title(r"expansion exponents $\to$ 2/3 (matter), 1/2 (radiation)")
fig.tight_layout(); fig.savefig("figs/fig9_friedmann.png", dpi=150); plt.close(fig)

# ============ (C) separate-universe chain: growth of structure & collapse ============
print("# (C) cell chain with overdensity: linear growth, turnaround, collapse")
Mc = 41
ic = np.arange(Mc) - Mc//2
delta0 = 0.08*np.exp(-(ic/4.0)**2)
Om = 1 + delta0
H0, a0 = 1.0, 0.02
amin = 0.25*a0

def rhs(t, y):
    a, ad = y[:Mc], y[Mc:]
    acc = -0.5*H0**2*Om/np.maximum(a, amin)**2
    frozen = a <= amin
    ad = np.where(frozen, 0.0, ad)
    acc = np.where(frozen, 0.0, acc)
    return np.concatenate([ad, acc])

ad0 = H0*np.sqrt(Om/a0 + (1 - Om))            # synchronized constraint-consistent ICs
y0 = np.concatenate([np.full(Mc, a0), ad0])
tend = 160.0
tev = np.geomspace(2e-3, tend, 420)
solc = solve_ivp(rhs, [0, tend], y0, t_eval=tev, method="LSODA", rtol=1e-9, atol=1e-12)
A = np.maximum(solc.y[:Mc].T, amin)            # (T, Mc)
abg = A[:, 0]                                  # far cell ~ background (delta ~ 1e-8)
dcen = (1 + delta0[Mc//2])*(abg/A[:, Mc//2])**3 - 1
kta = int(np.argmax(A[:, Mc//2]))
hit = A[:, Mc//2] <= amin*1.001
kcol = int(np.argmax(hit)) if hit.any() else len(tev)-1
print(f"  central delta(0) = {delta0[Mc//2]:.3f}; turnaround at a_bg/a0 = {abg[kta]/a0:.1f}, "
      f"collapse at a_bg/a0 = {abg[kcol]/a0:.1f}")
w = (abg/a0 > 30) & (abg/a0 < 300)             # late linear window: pure growing mode
slope = np.polyfit(np.log(abg[w]), np.log(dcen[w]), 1)[0]
print(f"  linear-growth exponent d ln delta / d ln a = {slope:.3f}   (theory: 1)")

fig, (c1, c2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
c1.loglog(abg/abg[0], dcen, "crimson", lw=2, label=r"central cell $\delta(t)$ (simulation)")
Cfit = np.exp(np.mean(np.log(dcen[w]) - np.log(abg[w])))
c1.loglog(abg/abg[0], Cfit*abg, "k--", label=r"linear theory $\delta \propto a$")
c1.axhline(1.686, color="gray", ls=":", label=r"$\delta_{\rm lin}=1.686$ at collapse (theory)")
c1.axvline(abg[kta]/abg[0], color="gray", lw=0.8)
c1.text(abg[kta]/abg[0]*1.05, 3e-2, "turnaround", rotation=90, fontsize=8, color="gray")
c1.set_xlabel(r"background expansion $a_{\rm bg}/a_0$"); c1.set_ylabel(r"density contrast $\delta$")
c1.set_title("gravitational instability from the cell chain"); c1.legend(fontsize=8)
ks = [0, kta//2, kta, (kta+kcol)//2, min(kcol+5, len(tev)-1)]
for k, col in zip(ks, plt.cm.viridis(np.linspace(0, 0.9, 5))):
    c2.semilogy(ic, A[k]/abg[k], color=col, label=f"$a_{{bg}}/a_0={abg[k]/a0:.0f}$")
c2.set_xlabel("cell index"); c2.set_ylabel(r"$a_i / a_{\rm bg}$")
c2.set_title("the well forms: local scale factor lags,\nturns around and collapses")
c2.legend(fontsize=8)
fig.tight_layout(); fig.savefig("figs/fig10_collapse.png", dpi=150); plt.close(fig)

np.savez_compressed("collapse_traj.npz", t=tev, A=A, abg=abg, dcen=dcen,
                    delta0=delta0, ic=ic, amin=amin)
print("saved figs 8-10 and collapse_traj.npz")
