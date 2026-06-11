#!/usr/bin/env python3
"""
ch08_bag_spectrum.py — the `bag1d` module (Ch. 8/10/12 workhorse) + self-tests.

1+1D Dirac bag, chiral-wall family, CP mass:
    H = -i s1 d/dx + mR s3 - mI s2,   mR = m cos(delta), mI = m sin(delta)
    psi' = K psi,  K = [[-mI, i(E+mR)], [i(E-mR), mI]],  K^2 = -p^2 I, p^2 = E^2 - m^2
    walls: -i n.gamma e^{i theta gamma5} psi = psi
      x=0 ray: v0 = (cos a, -i sin a), a = pi/4 + th0/2
      x=L ray: vL = (cos b, +i sin b), b = pi/4 - thL/2

Self-tests (printed when run as main):
  A. MIT spectrum vs tan(pL) = -p/m roots            [Ch. 8.3]
  B. spectral mirror symmetry for bulk phases        [Ch. 10.4]  (eta ~ machine zero)
  C. chiral-wall master equation vs transfer matrix  [Ch. 12.2]
  D. (delta,delta)-walls real mass == bulk complex mass spectrum [Ch. 10.6/12.2(b)]
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

# ---------------------------------------------------------------- core
def rays(th0, thL):
    a = np.pi / 4 + th0 / 2
    b = np.pi / 4 - thL / 2
    v0 = np.array([np.cos(a), -1j * np.sin(a)], dtype=complex)
    vL = np.array([np.cos(b), +1j * np.sin(b)], dtype=complex)
    return v0, vL

def transfer(E, mR, mI, L):
    """e^{KL} in closed form using K^2 = -p^2 I."""
    p2 = E**2 - (mR**2 + mI**2)
    p = np.sqrt(complex(p2))
    K = np.array([[-mI, 1j * (E + mR)], [1j * (E - mR), mI]], dtype=complex)
    if abs(p) < 1e-14:
        return np.eye(2, dtype=complex) + K * L
    return np.cos(p * L) * np.eye(2, dtype=complex) + (np.sin(p * L) / p) * K

def residual(E, mR, mI, L, th0, thL):
    """Vanishes iff E is an eigenvalue: det test for psi(L) parallel to the x=L ray."""
    v0, vL = rays(th0, thL)
    psiL = transfer(E, mR, mI, L) @ v0
    return (vL[0] * psiL[1] - vL[1] * psiL[0]) / np.linalg.norm(psiL)

def find_spectrum(mR, mI, L, th0, thL, Emax, nE=4001, tol=1e-9):
    Es = np.linspace(-Emax, Emax, nE)
    r = np.abs([residual(E, mR, mI, L, th0, thL) for E in Es])
    roots = []
    for i in range(1, nE - 1):
        if r[i] < r[i - 1] and r[i] < r[i + 1] and r[i] < 0.5:
            res = minimize_scalar(
                lambda E: np.abs(residual(E, mR, mI, L, th0, thL))**2,
                bounds=(Es[i - 1], Es[i + 1]), method="bounded",
                options={"xatol": 1e-14})
            if np.abs(residual(res.x, mR, mI, L, th0, thL)) < 1e-6:
                roots.append(res.x)
    roots = sorted(roots)
    keep = []
    for x in roots:
        if not keep or x - keep[-1] > 1e-8 * max(1.0, abs(x)):
            keep.append(x)
    return np.array(keep)

def eigenfunction(E, mR, mI, L, th0, x):
    v0, _ = rays(th0, 0.0)
    psi = np.array([transfer(E, mR, mI, xi) @ v0 for xi in x])
    # normalize on [0, L]
    norm = np.trapezoid((np.abs(psi)**2).sum(axis=1), x)
    return psi / np.sqrt(norm)

def master_eq(E, m, L, th0, thL):
    """Chiral-Wall Master Equation residual (12.7), valid both branches."""
    D, S = thL - th0, 0.5 * (th0 + thL)
    p = np.sqrt(complex(E**2 - m**2))
    if abs(p) < 1e-14:
        return np.cos(D / 2) * L - (E * np.sin(D / 2) - m * np.cos(S)) * L**2 / 2  # limit form
    return (p * np.cos(D / 2) * np.cos(p * L)
            - np.sin(p * L) * (E * np.sin(D / 2) - m * np.cos(S)))

def eta_heat(spec, t):
    return np.sum(np.sign(spec) * np.exp(-t * np.abs(spec)))

# ---------------------------------------------------------------- tests
if __name__ == "__main__":
    m, L = 1.0, 1.0
    print("=== A. MIT bag: solver vs tan(pL) = -p/m ===")
    sp = find_spectrum(m, 0.0, L, 0.0, 0.0, Emax=12.0)
    ps = []
    for j in range(1, 5):
        lo, hi = (j - 0.4999) * np.pi / L + 1e-9, (j + 0.4999) * np.pi / L
        try:
            ps.append(brentq(lambda p: np.tan(p * L) + p / m, lo, hi, xtol=1e-13))
        except ValueError:
            pass
    Ep = np.sort(np.concatenate([+np.sqrt(np.array(ps)**2 + m**2),
                                 -np.sqrt(np.array(ps)**2 + m**2)]))
    got = sp[np.abs(sp) < np.max(np.abs(Ep)) + 0.3]
    print(f"  max |solver - tan-roots| = {np.max(np.abs(np.sort(got) - Ep)):.2e}"
          if len(got) == len(Ep) else f"  LEVEL COUNT MISMATCH {len(got)} vs {len(Ep)}")

    print("=== B. bulk phase: mirror symmetry and eta ===")
    for delta in (0.3, np.pi / 4):
        mR, mI = m * np.cos(delta), m * np.sin(delta)
        sp = find_spectrum(mR, mI, L, 0.0, 0.0, Emax=40.0, nE=8001)
        asym = np.max(np.abs(np.sort(sp) + np.sort(sp)[::-1]))
        print(f"  delta={delta:.3f}: levels={len(sp)}, E<->-E asymmetry={asym:.2e}, "
              f"eta(t=0.05)={eta_heat(sp, 0.05):+.2e}")

    print("=== C. chiral master equation vs transfer matrix ===")
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(5):
        th0, thL = rng.uniform(-2, 2, 2)
        sp = find_spectrum(m, 0.0, L, th0, thL, Emax=10.0)
        for E in sp:
            # master-equation residual should vanish at solver eigenvalues
            worst = max(worst, abs(master_eq(E, m, L, th0, thL)))
    print(f"  max |master-eq residual| over 5 angle pairs = {worst:.2e}")

    print("=== D. (delta,delta) walls + real mass == bulk complex mass ===")
    delta = 0.6
    spA = find_spectrum(m, 0.0, L, delta, delta, Emax=10.0)          # walls carry delta
    spB = find_spectrum(m * np.cos(delta), m * np.sin(delta), L, 0.0, 0.0, Emax=10.0)
    n = min(len(spA), len(spB))
    print(f"  max |spec diff| = {np.max(np.abs(np.sort(spA)[:n] - np.sort(spB)[:n])):.2e}")
