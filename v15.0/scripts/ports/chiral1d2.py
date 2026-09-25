"""
F4 numerics, exact-spectrum version.

(C1) Qvac = -eta/2 vs wall mismatch D (symmetric split th0=-D/2, thL=+D/2, S=0):
     law of fractional vacuum charge; check linearity in D at large mL.
(C2) th0=-2, thL=+2: level crossing at tanh(mL*) = -cos(2)  =>  L* = 0.44290 (m=1).
     E_min(L) crosses zero, Qvac(L) jumps by exactly 1.
(C3) Sudden expansion L=0.30 -> 0.70 across L*: spectral prediction
     DQ = Qvac(L') - Qvac(L) ~ -/+1; Bogoliubov forms must converge to it.
"""
import sys
import numpy as np
from chiral_exact import chiral_spectrum, crossing_L

def eta_extrap(E, ts=(0.06, 0.04, 0.025, 0.015)):
    v = [float(np.sum(np.sign(E)*np.exp(-t*np.abs(E)))) for t in ts]
    # Richardson: fit v = a + b t + c t^2 on last three points
    A = np.vstack([np.ones(3), ts[-3:], np.square(ts[-3:])]).T
    a = np.linalg.lstsq(A, v[-3:], rcond=None)[0][0]
    return v, a

def chiral_modes(m, L, th0, thL, Emax, ngrid):
    """Analytic eigenfunctions on a grid; returns (E, psi[nmode,2,ngrid], x)."""
    Es = chiral_spectrum(m, L, th0, thL, Emax)
    x = np.linspace(0, L, ngrid)
    a = np.pi/4 + th0/2
    psis = []
    for E in Es:
        s2 = E*E - m*m
        if s2 >= 0:
            p = np.sqrt(s2); C = np.cos(p*x); S = np.sin(p*x)/p if p else x
        else:
            q = np.sqrt(-s2); C = np.cosh(q*x); S = np.sinh(q*x)/q
        p1 = np.cos(a)*C + (E+m)*np.sin(a)*S
        p2 = -1j*(np.sin(a)*C - np.cos(a)*s2*S/(E+m))
        psi = np.vstack([p1.astype(complex), p2])
        nrm = np.sqrt(np.trapezoid(np.abs(psi[0])**2 + np.abs(psi[1])**2, x))
        psis.append(psi/nrm)
    return Es, np.array(psis), x

job = sys.argv[1] if len(sys.argv) > 1 else "all"
m = 1.0

if job in ("qvac", "all"):
    print("# (C1) Qvac vs D, symmetric walls th=(-D/2,+D/2), mL=6 (decoupled walls), Emax=350")
    L = 6.0
    for D in (0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8):
        E = chiral_spectrum(m, L, -D/2, D/2, Emax=350.0)
        v, ext = eta_extrap(E)
        print(f"  D={D:3.1f}: levels={len(E):4d}  Qvac={-0.5*ext:+.5f}   D/(2 pi)={D/(2*np.pi):+.5f}")

if job in ("cross", "all"):
    t0, tL = -2.0, 2.0
    Lstar = crossing_L(m, t0, tL)
    print(f"# (C2) th=({t0},{tL}): predicted crossing at L* = {Lstar:.5f}")
    for L in (0.25, 0.35, 0.42, 0.4429, 0.47, 0.55, 0.70, 0.90):
        E = chiral_spectrum(m, L, t0, tL, Emax=400.0)
        e0 = E[np.argmin(np.abs(E))]
        v, ext = eta_extrap(E)
        print(f"  L={L:6.4f}: E_min={e0:+.5f}   Qvac={-0.5*ext:+.5f}")

if job in ("quench", "all"):
    from quench1d import overlaps, dq_forms
    t0, tL, L, s = -2.0, 2.0, 0.30, 0.70/0.30
    print(f"# (C3) quench L={L} -> {s*L:.2f} across the crossing, th=({t0},{tL})")
    Eo, po, xo = chiral_modes(m, L,  t0, tL, Emax=420.0, ngrid=2501)
    En, pn, xn = chiral_modes(m, s*L, t0, tL, Emax=420.0, ngrid=5001)
    a = overlaps(Eo, po, xo, En, pn, xn)
    # completeness sanity: sum_k |a_{k n}|^2 = 1 for low old modes
    comp = np.sum(np.abs(a)**2, axis=0)
    print(f"  completeness of lowest 6 old modes: {np.round(comp[np.argsort(np.abs(Eo))[:6]],5)}")
    for N in (10, 20, 35, 55, 75):
        dqc, dqs = dq_forms(Eo, En, a, N)
        print(f"  N={N:3d}: dQ_cross={dqc:+.5f}  dQ_same={dqs:+.5f}")
    _, exo = eta_extrap(Eo); _, exn = eta_extrap(En)
    print(f"  spectral prediction DQ = Qvac(L)-Qvac(L') = {-0.5*(exo-exn):+.5f}"
          f"   [Qvac(L)={-0.5*exo:+.5f}, Qvac(L')={-0.5*exn:+.5f}]")
    # control: same quench but NO crossing (equal small walls -> D=0)
    print("# (C3-control) same quench with D=0 walls th=(0.7,0.7) (pure mass renorm.)")
    Eo, po, xo = chiral_modes(m, L,  0.7, 0.7, Emax=420.0, ngrid=2501)
    En, pn, xn = chiral_modes(m, s*L, 0.7, 0.7, Emax=420.0, ngrid=5001)
    a = overlaps(Eo, po, xo, En, pn, xn)
    for N in (20, 55, 75):
        dqc, dqs = dq_forms(Eo, En, a, N)
        print(f"  N={N:3d}: dQ_cross={dqc:+.5f}  dQ_same={dqs:+.5f}")
    _, exo = eta_extrap(Eo); _, exn = eta_extrap(En)
    print(f"  spectral prediction DQ = {-0.5*(exo-exn):+.5f}")
