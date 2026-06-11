"""
Sudden expansion L -> sL of the genuine 1+1D CP-mass bag.

Key analytic statements being tested:
 (T1) The spectrum at phase delta is E-symmetric for all |delta|<pi/2  => vacuum
      charge of every box vanishes => NET charge production is exactly zero.
 (T2) The infinite-basis limit of the cross-branch Bogoliubov sum is NOT zero but
      equals a charge-POLARISATION integral:
          DQ_cross(inf) = (1/2) Tr_{[0,L]} sgn(H_{L'})
                        = (1/2) sum_k sgn(E'_k) w_k ,   w_k = int_0^L |psi'_k|^2
      i.e. minus the new-bag vacuum-charge that resides outside [0,L].  The same
      quantity integrated over the whole new box vanishes (eta = 0): what the
      quench produces is a redistribution of vacuum polarisation between the two
      walls, not net charge.
 (T3) Paper 3's A(s=2, mL=1.085) = 0.0764 at nmax=18 (same-branch form) is
      reproduced by the genuine solver -- their 1D machinery is correct -- but the
      number is a truncation of the polarisation integral, not pair production.
"""
import numpy as np
from bag1d import bc_vector, K, find_spectrum, s1, s2, s3, I2
np.set_printoptions(precision=6, suppress=False)

def plain_spectrum(mR, m, Lbox, Emax):
    """Exact spectrum of the plain CP bag: tan(pL) = -p/mR, E = +-sqrt(p^2+m^2).  Complete."""
    from scipy.optimize import brentq
    ps, j = [], 1
    while True:
        lo = (j - 0.5)*np.pi/Lbox + 1e-12
        hi = j*np.pi/Lbox - 1e-12
        f = lambda p: np.tan(p*Lbox) + p/mR
        if f(lo) * f(hi) < 0:
            ps.append(brentq(f, lo, hi, xtol=1e-14))
        if np.sqrt((j*np.pi/Lbox)**2 + m**2) > Emax:
            break
        j += 1
    ps = np.array(ps)
    E = np.sqrt(ps**2 + m**2)
    return np.sort(np.concatenate([-E[E <= Emax], E[E <= Emax]]))

def modes(mR, mI, Lbox, th0, thL, Emax, ngrid=3001):
    """Return energies (sorted), and wavefunctions sampled on x-grid of [0,Lbox], unit-normalised."""
    if th0 == 0 and thL == 0:
        m = np.hypot(mR, mI)
        Es = plain_spectrum(mR, m, Lbox, Emax)
    else:
        Es = find_spectrum(mR, mI, Lbox, th0, thL, Emax=Emax,
                           nE=max(4000, int(60*Emax*Lbox/np.pi)))
    x = np.linspace(0, Lbox, ngrid)
    psis = []
    v0 = bc_vector(th0, '0')
    for E in Es:
        lam, V = np.linalg.eig(K(E, mR, mI))
        c = np.linalg.solve(V, v0)
        ph = np.exp(np.outer(x, lam))            # (ngrid,2)
        psi = (V @ (c[:, None] * ph.T))          # (2, ngrid)
        nrm = np.sqrt(np.trapezoid(np.abs(psi[0])**2 + np.abs(psi[1])**2, x))
        psis.append(psi / nrm)
    return np.array(Es), np.array(psis), x      # psis: (nmode, 2, ngrid)

def overlaps(E_old, psi_old, x_old, E_new, psi_new, x_new):
    """alpha_{k n} = int_0^L psi_new_k^dag psi_old_n dx  (old modes vanish beyond L)."""
    # interpolate new modes onto old grid
    a = np.zeros((len(E_new), len(E_old)), dtype=complex)
    for k in range(len(E_new)):
        pk = np.vstack([np.interp(x_old, x_new, psi_new[k, c].real) +
                        1j*np.interp(x_old, x_new, psi_new[k, c].imag) for c in range(2)])
        for n in range(len(E_old)):
            integ = (pk.conj() * psi_old[n]).sum(axis=0)
            a[k, n] = np.trapezoid(integ, x_old)
    return a

def dq_forms(E_old, E_new, a, Ncut):
    """Both Delta-Q forms with N lowest |E| modes per branch on each side."""
    def pick(E, N, sign):
        idx = np.where(np.sign(E) == sign)[0]
        idx = idx[np.argsort(np.abs(E[idx]))][:N]
        return idx
    op, om = pick(E_old, Ncut, +1), pick(E_old, Ncut, -1)
    np_, nm = pick(E_new, Ncut, +1), pick(E_new, Ncut, -1)
    A_pp = a[np.ix_(np_, op)]; A_pm = a[np.ix_(np_, om)]
    A_mp = a[np.ix_(nm, op)]; A_mm = a[np.ix_(nm, om)]
    dq_cross = np.sum(np.abs(A_pm)**2) - np.sum(np.abs(A_mp)**2)
    dq_same  = np.sum(np.abs(A_pp)**2) - np.sum(np.abs(A_mm)**2)
    return dq_cross, dq_same

def polarisation_prediction(mR, mI, L, Lp, th0, thL, Emax=300.0, ts=(0.04, 0.02, 0.01)):
    """ (1/2) sum_k sgn(E'_k) e^{-t|E'_k|} w_k  with w_k = int_0^L |psi'_k|^2, extrapolated t->0,
        plus the whole-box check (should -> 0)."""
    E, psi, x = modes(mR, mI, Lp, th0, thL, Emax, ngrid=4001)
    mask = x <= L + 1e-12
    w = np.array([np.trapezoid((np.abs(p[0])**2 + np.abs(p[1])**2)[mask], x[mask]) for p in psi])
    out, tot = [], []
    for t in ts:
        damp = np.exp(-t*np.abs(E))
        out.append(0.5*np.sum(np.sign(E)*damp*w))
        tot.append(0.5*np.sum(np.sign(E)*damp))
    # Richardson on the last two (assume linear in t)
    t1, t2 = ts[-2], ts[-1]
    extrap = (out[-1]*t1 - out[-2]*t2)/(t1 - t2)
    return out, tot, extrap

if __name__ == "__main__":
    import sys, json
    job = sys.argv[1] if len(sys.argv) > 1 else "all"

    if job in ("repro", "all"):
        # (T3) reproduce paper 3: s=2, mL=1.085, nmax=18, same-branch form, A = dQ/sin 2d
        m, L, s = 1.085, 1.0, 2.0
        delta = 0.015
        mR, mI = m*np.cos(delta), m*np.sin(delta)
        Emax_old = np.sqrt((18*np.pi/L)**2 + m**2) + 1
        Emax_new = np.sqrt((18*np.pi/(s*L))**2 + m**2) + 1
        Eo, po, xo = modes(mR, mI, L, 0, 0, Emax_old)
        En, pn, xn = modes(mR, mI, s*L, 0, 0, Emax_new, ngrid=4001)
        a = overlaps(Eo, po, xo, En, pn, xn)
        dqc, dqs = dq_forms(Eo, En, a, 18)
        print(f"[repro] s=2 mL=1.085 d=0.015 nmax=18:  dQ_same/sin2d = {dqs/np.sin(2*delta):.4f}"
              f"   dQ_cross/sin2d = {dqc/np.sin(2*delta):.4f}   (paper: A=0.0764, ~10% lower for cross)")

    if job in ("conv", "all"):
        # (T2) convergence study at delta=pi/4, s=1.2, mL=1  (paper 3 Fig.4: dQ=+0.069)
        m, L, s, delta = 1.0, 1.0, 1.2, np.pi/4
        mR, mI = m*np.cos(delta), m*np.sin(delta)
        Nmax = 90
        Emax_old = np.sqrt((Nmax*np.pi/L)**2 + m**2) + 1
        Emax_new = np.sqrt((Nmax*np.pi/(s*L))**2 + m**2) + 1
        Eo, po, xo = modes(mR, mI, L, 0, 0, Emax_old, ngrid=6001)
        En, pn, xn = modes(mR, mI, s*L, 0, 0, Emax_new, ngrid=7001)
        a = overlaps(Eo, po, xo, En, pn, xn)
        rows = []
        for N in (6, 10, 14, 18, 26, 38, 54, 70, 85):
            dqc, dqs = dq_forms(Eo, En, a, N)
            rows.append((N, dqc, dqs))
            print(f"[conv] N={N:3d}:  dQ_cross = {dqc:+.5f}   dQ_same = {dqs:+.5f}")
        out, tot, extrap = polarisation_prediction(mR, mI, L, s*L, 0, 0)
        print(f"[conv] polarisation prediction (heat-kernel t={0.04,0.02,0.01}): {np.round(out,5)}"
              f" -> extrap {extrap:+.5f}")
        print(f"[conv] whole-box eta check (must -> 0): {np.round(tot,6)}")
        json.dump({"rows": rows, "pred": extrap}, open("conv_results.json", "w"))

    if job in ("sym", "all"):
        # (T1) eta of every box vanishes for the plain CP bag
        m, delta = 1.0, np.pi/4
        mR, mI = m*np.cos(delta), m*np.sin(delta)
        for Lb in (1.0, 1.2):
            E, psi, x = modes(mR, mI, Lb, 0, 0, Emax=200.0, ngrid=801)
            for t in (0.05, 0.02):
                eta = np.sum(np.sign(E)*np.exp(-t*np.abs(E)))
                print(f"[sym] L={Lb}: eta(t={t}) = {eta:.2e}   (levels: {len(E)})")
