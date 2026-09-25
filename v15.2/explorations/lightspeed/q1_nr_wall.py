"""
Q1a - non-relativistic particle (hbar=m=1) in mode n of [0,L0=1]; the wall moves from L0 to
L1 = L0 (n+1)/n at constant speed u (step time t0 = dL/u), then stops.
Exact TDSE in comoving coordinates xi = x/L (Ch. 17):
   i chi_t = -(1/2L^2) chi_xixi + i (Ldot/L) (xi d_xi + 1/2) chi
Crank-Nicolson on a uniform xi-grid (unitary: the dilation term is discretised antisymmetrically).
Outputs: <E>_final / E_n and P(promoted mode n+1).  Classical ensemble for comparison.
"""
import numpy as np
from scipy.linalg import solve_banded

def run_quantum(n, u, N=1600, steps_per_unit=None):
    L0, L1 = 1.0, (n+1)/n; dL = L1-L0; tau = dL/u
    h = 1.0/(N+1); xi = np.arange(1, N+1)*h
    chi = (np.sqrt(2)*np.sin(n*np.pi*xi)*np.sqrt(h)).astype(complex)        # discrete normalisation sum|chi|^2 = 1
    E_n = 0.5*(n*np.pi)**2
    nsteps = int(max(400, 40*E_n*tau, 40*u*tau/h*0 + 400))
    # accuracy: dt << 1/E_n  and wall moves < ~0.2 grid cells per step in xi-space is irrelevant (comoving)
    nsteps = int(max(nsteps, 60*E_n*tau))
    dt = tau/nsteps
    # operators: T = -(1/2) D2 (tridiagonal), A = 1/2 (xi D1 + D1 xi) antisymmetric tridiagonal
    d2_main = -2.0/h**2*np.ones(N); d2_off = 1.0/h**2*np.ones(N-1)
    xm = 0.5*(xi[:-1]+xi[1:])
    a_up = xm/(2*h)        # A[i,i+1]
    a_lo = -xm/(2*h)       # A[i+1,i]
    def apply(Lc, Ld, v):
        a = -0.5/Lc**2; b = 1j*Ld/Lc
        out = a*d2_main*v
        out[:-1] += a*d2_off*v[1:] + b*a_up*v[1:]
        out[1:] += a*d2_off*v[:-1] + b*a_lo*v[:-1]
        return out
    t = 0.0
    for k in range(nsteps):
        Lc = L0 + u*(t+0.5*dt); Ld = u
        a = -0.5/Lc**2; b = 1j*Ld/Lc
        # (1 + i dt/2 H) chi_new = (1 - i dt/2 H) chi
        rhs = chi - 0.5j*dt*apply(Lc, Ld, chi)
        ab = np.zeros((3, N), complex)
        ab[1] = 1 + 0.5j*dt*a*d2_main
        ab[0, 1:] = 0.5j*dt*(a*d2_off + b*a_up)       # upper diagonal
        ab[2, :-1] = 0.5j*dt*(a*d2_off + b*a_lo)      # lower diagonal
        chi = solve_banded((1, 1), ab, rhs); t += dt
    # final box L1: energy and mode populations (discrete sine basis = exact eigvecs of D2)
    k = np.arange(1, N+1)
    modes = np.sqrt(2*h)*np.sin(np.outer(k, xi)*np.pi)
    c = modes @ chi
    P = np.abs(c)**2
    Ek = (2/h**2)*(1-np.cos(k*np.pi*h))/2/L1**2      # discrete -D2/2 eigenvalues
    E_fin = (P*Ek).sum()
    En_disc = (2/h**2)*(1-np.cos(n*np.pi*h))/2
    return E_fin/En_disc, P[n], P[n-1], P.sum()      # P[n] is mode n+1

def run_classical(n, u, M=200000, rng=np.random.default_rng(0)):
    v = n*np.pi; L0, L1 = 1.0, (n+1)/n; tau = (L1-L0)/u
    x = rng.random(M); s = rng.choice([-1.0, 1.0], M); vel = s*v
    t = np.zeros(M); E0 = 0.5*v*v
    active = np.ones(M, bool)
    for it in range(100000):
        if not active.any(): break
        # time to next wall hit
        wall_now = L0 + u*t
        tr = np.where(vel > u, (wall_now - x)/(vel - u), np.inf)
        tl = np.where(vel < 0, -x/vel, np.inf)
        th = np.minimum(tr, tl); tnext = t + th
        done = tnext >= tau
        # finished particles: drift to tau
        active &= ~done
        idx = np.where(active)[0]
        if idx.size == 0: break
        x[idx] += vel[idx]*th[idx]; t[idx] = tnext[idx]
        hitr = idx[tr[idx] <= tl[idx]]; hitl = idx[tr[idx] > tl[idx]]
        vel[hitr] = 2*u - vel[hitr]; vel[hitl] = -vel[hitl]
    return np.mean(0.5*vel**2)/E0

if __name__ == "__main__":
    # convergence / limit checks
    for n in (3, 10):
        for N in (800, 1600):
            r = run_quantum(n, 200*n*np.pi, N=N)
            print(f"sudden check n={n} N={N}: E/E_n={r[0]:.5f} (exact 1), P(n+1)={r[1]:.5f} (exact n/(n+1)={n/(n+1):.5f}), norm {r[3]:.8f}")
        r = run_quantum(n, 0.003*n*np.pi, N=800)
        print(f"adiabatic check n={n}: E/E_n={r[0]:.5f} (exact (n/(n+1))^2={(n/(n+1))**2:.5f}), P(n)={r[2]:.5f}")
