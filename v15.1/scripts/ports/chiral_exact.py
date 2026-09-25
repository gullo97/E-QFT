"""
Exact spectrum of the 1+1D chiral bag, real mass m, wall angles th0, thL.

Master eigenvalue equation (derived analytically from the transfer matrix):

    p cos(D/2) cos(pL) = sin(pL) [ E sin(D/2) - m cos(S) ]
    D = thL - th0,   S = (th0 + thL)/2,   p^2 = E^2 - m^2  (p -> iq below gap)

Checks encoded here:
  * D=0, S=0      ->  tan(pL) = -p/m            (MIT)
  * D=0, S=delta  ->  tan(pL) = -p/(m cos delta) (bulk CP phase rotated to walls!)
  * E -> -E maps the equation to itself with D -> -D: asymmetry iff sin(D/2)!=0.
  * E=0 solution exists iff tanh(mL) = -cos(D/2)/cos(S): level crossing condition.
"""
import numpy as np
from scipy.optimize import brentq

def _CS(E, m, L):
    """C = cos(pL), S = sin(pL)/p as entire real functions of E (works below gap)."""
    s2 = E*E - m*m
    if s2 >= 0:
        p = np.sqrt(s2)
        return np.cos(p*L), (L if p == 0 else np.sin(p*L)/p)
    q = np.sqrt(-s2)
    return np.cosh(q*L), np.sinh(q*L)/q

def phi(E, m, L, th0, thL):
    D, S = thL - th0, 0.5*(th0 + thL)
    C, Sf = _CS(E, m, L)
    return -np.cos(D/2)*C + Sf*(E*np.sin(D/2) - m*np.cos(S))

def chiral_spectrum(m, L, th0, thL, Emax):
    """Complete spectrum in [-Emax, Emax] via sign changes of phi on a fine grid."""
    # grid spacing well below level spacing pi/L (asymptotically)
    n = int(40 * (2*Emax) * L / np.pi) + 4000
    Es = np.linspace(-Emax, Emax, n)
    vals = np.array([phi(E, m, L, th0, thL) for E in Es])
    roots = []
    for i in range(len(Es)-1):
        if vals[i] == 0.0:
            roots.append(Es[i])
        elif vals[i]*vals[i+1] < 0:
            roots.append(brentq(phi, Es[i], Es[i+1], args=(m, L, th0, thL), xtol=1e-13))
    return np.array(roots)

def crossing_L(m, th0, thL):
    """Box size at which a level sits exactly at E=0 (if it exists)."""
    D, S = thL - th0, 0.5*(th0 + thL)
    r = -np.cos(D/2)/np.cos(S)
    return np.arctanh(r)/m if 0 < r < 1 else None

if __name__ == "__main__":
    from bag1d import find_spectrum
    m, L = 1.0, 1.0
    print("=== verify master equation against transfer-matrix solver ===")
    for (t0, tL) in ((0.0, 0.0), (0.0, 0.5), (0.0, 1.0), (0.0, 1.5), (-0.7, 1.1)):
        ex = chiral_spectrum(m, L, t0, tL, Emax=12.0)
        tm = find_spectrum(m, 0.0, L, t0, tL, Emax=12.0, nE=8000)
        k = min(len(ex), len(tm))
        # align by value
        d = max(abs(a - b) for a, b in zip(np.sort(ex)[:k], np.sort(tm)[:k]))
        print(f" th=({t0:+.1f},{tL:+.1f}): exact {len(ex)} levels, TM {len(tm)} levels, "
              f"max|diff| = {d:.2e}")
    print()
    print("=== bulk CP phase == common wall shift:  tan(pL) = -p/(m cos d) ===")
    for d in (0.3, 0.8):
        ex = chiral_spectrum(m, L, d, d, Emax=12.0)       # walls shifted by d, real mass
        tm = find_spectrum(m*np.cos(d), m*np.sin(d), L, 0.0, 0.0, Emax=12.0, nE=8000)
        k = min(len(ex), len(tm))
        diff = max(abs(a-b) for a, b in zip(np.sort(ex)[:k], np.sort(tm)[:k]))
        print(f" delta={d}: max|exact - bulkCP| = {diff:.2e}   ({len(ex)} vs {len(tm)} levels)")
