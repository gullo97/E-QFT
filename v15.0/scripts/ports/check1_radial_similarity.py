"""
Check 1 (paper 2 foundation):
Claim (new analytic result): the 'bulk-CP radial prescription', Eq. (2) of paper 2,

    g' = (E+mR) f + mI g - (1+k)/r g
    f' = -(E-mR) g + mI f - (1-k)/r f          (mR = m cos d, mI = m sin d)

is mapped by  (g,f) = e^{mI r} (G,F)  onto the STANDARD (delta=0) MIT-bag radial
system with mass mR = m cos(delta).  Since the MIT boundary condition
f(R) = -g(R) is invariant under a common positive rescaling, the whole spectrum
at CP phase delta equals the delta=0 spectrum with mass m -> m cos(delta).

We verify:
 (a) eigenvalues of the delta>0 system  ==  eigenvalues of delta=0 system with m->m cos d
 (b) the wavefunction ratio g_delta(r)/G_0(r) == e^{mI r}  pointwise.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import spherical_jn

m, R = 1.7, 1.0

def rhs(r, y, E, mR, mI, kappa):
    g, f = y
    dg = (E + mR)*f + mI*g - (1+kappa)/r*g
    df = -(E - mR)*g + mI*f - (1-kappa)/r*f
    return [dg, df]

def residual(E, mR, mI, kappa):
    """Integrate from r0 with the regular solution and return MIT residual f(R)+g(R)."""
    r0 = 1e-6
    lup  = kappa if kappa > 0 else -kappa - 1
    # regular behaviour g ~ r^lup ; get f from the g' equation at r0 (series consistent)
    g0 = r0**lup
    # from dg eq: f = (g' - mI g + (1+k)/r g)/(E+mR); g' = lup r^{lup-1}
    f0 = (lup*r0**(lup-1) - mI*g0 + (1+kappa)/r0*g0)/(E + mR)
    sol = solve_ivp(rhs, (r0, R), [g0, f0], args=(E, mR, mI, kappa),
                    rtol=1e-9, atol=1e-300, dense_output=True, method='DOP853')
    g, f = sol.y[0, -1], sol.y[1, -1]
    nrm = np.hypot(g, f)
    return (f + g)/nrm, sol

def spectrum(mR, mI, kappa, Emax=12.0, nE=240):
    Es = np.linspace(mR*0.0 + 0.05, Emax, nE)   # particle branch E>0
    res = np.array([residual(E, mR, mI, kappa)[0] for E in Es])
    roots = []
    for i in range(len(Es)-1):
        if np.isfinite(res[i]) and np.isfinite(res[i+1]) and res[i]*res[i+1] < 0:
            roots.append(brentq(lambda E: residual(E, mR, mI, kappa)[0], Es[i], Es[i+1], xtol=1e-12))
    return np.array(roots)

print("=== (a) spectral identity: delta-system  vs  delta=0 with m->m cos d ===")
for kappa in (+1, -1, +2):
    for delta in (0.3, 0.6, 1.0):
        mR, mI = m*np.cos(delta), m*np.sin(delta)
        E_d  = spectrum(mR, mI, kappa)[:3]
        E_0  = spectrum(mR, 0.0, kappa)[:3]      # standard system, mass mR, no mI
        diff = np.max(np.abs(E_d - E_0))
        print(f" kappa={kappa:+d} d={delta:.1f}:  E(delta)={np.round(E_d,8)}  max|dE| = {diff:.2e}")

print()
print("=== (b) wavefunction dressing g_delta(r) / g_0(r) = e^{mI r} ===")
kappa, delta = +1, 0.6
mR, mI = m*np.cos(delta), m*np.sin(delta)
E1 = spectrum(mR, mI, kappa)[0]
_, sol_d = residual(E1, mR, mI, kappa)
_, sol_0 = residual(E1, mR, 0.0, kappa)
rs = np.linspace(0.05, R, 12)
gd = sol_d.sol(rs)[0]; g0 = sol_0.sol(rs)[0]
ratio = gd/g0
pred  = np.exp(mI*rs)
# both solutions share the same r->0 normalisation, so the ratio should be exp(mI r) exactly
print(" r       g_d/g_0        exp(mI r)     rel.err")
for r, a, b in zip(rs, ratio, pred):
    print(f" {r:.3f}  {a: .9f}  {b: .9f}  {abs(a/b-1):.2e}")
