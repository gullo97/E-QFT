"""
Genuine 1+1D MIT/chiral bag with CP-violating mass, NO reduction prescription.

H = -i s1 d/dx + mR s3 - mI s2        (Hermitian; mR=m cos d, mI=m sin d)
exactly the second-quantised single-particle Hamiltonian of L = psibar(i dslash - m e^{i d g5})psi
with gamma0=s3, gamma1=i s2, gamma5=s1.

Boundary conditions (chiral-bag family, theta=0 is MIT):
  at x=L  (outward +x):  -i gamma1 psi = e^{i theta_L gamma5} psi
  at x=0  (outward -x):  +i gamma1 psi = e^{i theta_0 gamma5} psi
For the plain CP-mass bag, theta_0 = theta_L = 0.

Solver: psi' = K(E) psi with K = i s1 (E - M);  psi(0) fixed by the x=0 BC ray;
eigenvalue when psi(L) lies in the x=L BC ray.
"""
import numpy as np
from scipy.linalg import expm, null_space
from scipy.optimize import brentq

s1 = np.array([[0,1],[1,0]], dtype=complex)
s2 = np.array([[0,-1j],[1j,0]], dtype=complex)
s3 = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
g1 = 1j*s2
g5 = s1

def bc_vector(theta, side):
    """Return the 1-dim ray (unit vector) allowed by the chiral BC on the given wall."""
    sgn = +1 if side == 'L' else -1          # -i n g1 psi = e^{i th g5} psi, n=+1 at L, -1 at 0
    P = expm(-1j*theta*g5) @ (sgn * (-1j) * g1)   # P psi = psi, P^2 = 1
    w, v = np.linalg.eig(P)
    idx = np.argmin(np.abs(w - 1))
    assert abs(w[idx]-1) < 1e-12, w
    return v[:, idx]

def K(E, mR, mI):
    M = mR*s3 - mI*s2
    return 1j * s1 @ (E*I2 - M)

def psi_of_x(E, mR, mI, L, th0, x):
    v0 = bc_vector(th0, '0')
    return expm(K(E, mR, mI)*x) @ v0

def boundary_residual(E, mR, mI, L, th0, thL):
    vL = bc_vector(thL, 'L')
    # orthogonal complement of vL
    w = null_space(vL.conj().reshape(1, 2))[:, 0]
    psiL = psi_of_x(E, mR, mI, L, th0, L)
    return (w.conj() @ psiL) / np.linalg.norm(psiL)

def real_residual(E, mR, mI, L, th0, thL):
    """For root finding: the residual has a smooth complex phase; use a phase-fixed real part."""
    r = boundary_residual(E, mR, mI, L, th0, thL)
    return r

def find_spectrum(mR, mI, L, th0, thL, Emax, nE=6000):
    Es = np.linspace(-Emax, Emax, nE)
    rs = np.array([boundary_residual(E, mR, mI, L, th0, thL) for E in Es])
    # the residual is complex; eigenvalues are where |r| dips to ~0.  Use sign changes of
    # the dominant real linear combination: rotate by local phase.
    roots = []
    absr = np.abs(rs)
    for i in range(1, len(Es)-1):
        if absr[i] < absr[i-1] and absr[i] < absr[i+1] and absr[i] < 0.5:
            # refine by minimising |r|^2 with golden / brent on the bracket
            from scipy.optimize import minimize_scalar
            res = minimize_scalar(lambda E: np.abs(boundary_residual(E, mR, mI, L, th0, thL))**2,
                                  bracket=None, bounds=(Es[i-1], Es[i+1]), method='bounded',
                                  options={'xatol':1e-13})
            if np.abs(boundary_residual(res.x, mR, mI, L, th0, thL)) < 1e-7:
                roots.append(res.x)
    roots = np.array(sorted(roots))
    # deduplicate
    if len(roots):
        keep = [roots[0]]
        for r in roots[1:]:
            if r - keep[-1] > 1e-8*max(1, abs(r)):
                keep.append(r)
        roots = np.array(keep)
    return roots

if __name__ == "__main__":
    m, L = 1.0, 1.0

    print("=== A. plain CP-mass bag (theta walls = 0): spectrum vs tan(pL) = -p/(m cos d) ===")
    for delta in (0.0, 0.3, 0.8):
        mR, mI = m*np.cos(delta), m*np.sin(delta)
        sp = find_spectrum(mR, mI, L, 0.0, 0.0, Emax=12.0)
        # predicted: roots of tan(pL) = -p/mR, p>0 ; E = +-sqrt(p^2+m^2)
        ps = []
        for j in range(1, 5):
            lo, hi = (j-0.4999)*np.pi/L + 1e-9, (j+0.4999)*np.pi/L
            f = lambda p: np.tan(p*L) + p/mR
            try: ps.append(brentq(f, lo, hi, xtol=1e-13))
            except ValueError: pass
        Epred = np.sort(np.concatenate([np.sqrt(np.array(ps)**2+m**2),
                                        -np.sqrt(np.array(ps)**2+m**2)]))
        got = sp[np.abs(sp) < np.max(np.abs(Epred))+0.5]
        print(f" d={delta:.1f}: solver E = {np.round(got,8)}")
        print(f"        tan-eq E = {np.round(Epred,8)}")
        print(f"        max |diff| = {np.max(np.abs(np.sort(got)-Epred)):.2e}"
              if len(got)==len(Epred) else "        LEVEL-COUNT MISMATCH")
        asym = np.max(np.abs(np.sort(got) + np.sort(got)[::-1]))
        print(f"        E -> -E asymmetry of spectrum: {asym:.2e}")

    print()
    print("=== B. chiral walls theta_0 != theta_L: spectral asymmetry appears ===")
    mR, mI = m, 0.0
    for th in (0.0, 0.5, 1.0, 1.5):
        sp = find_spectrum(mR, mI, L, 0.0, th, Emax=12.0)
        asym = np.max(np.abs(np.sort(sp) + np.sort(sp)[::-1])) if len(sp)%2==0 else np.nan
        print(f" theta_L={th:.1f}: n_levels(|E|<12)={len(sp)}, E in (-4,4): {np.round(sp[np.abs(sp)<4],6)}")
