#!/usr/bin/env python3
"""
ch12_chiral_bag.py — the two-wall chiral bag, end to end (Ch. 12).

Sections (each prints the numbers quoted in the chapter):
  A. boundary rays (12.3)-(12.4) vs direct eigendecomposition of -i n.g e^{i th g5}
  B. Chiral-Wall Master Equation (12.7) vs the independent transfer-matrix solver
  C. (delta,delta) walls + real mass == bulk complex-mass spectrum (closes Ch. 10)
  D. Fractional Charge Law: Q_vac vs -Delta/2pi at Sigma=0, mL=6   -> Fig. 12.1
  E. the crossing at theta=(-2,+2): plateaus and the unit jump      -> Fig. 12.2
  F. truncated Bogoliubov forms for the pump and control quenches (Ch. 11 vs 12)

The complete-spectrum engine is the *entire-function* form of (12.7),
    phi(E) = -cos(D/2) C(E) + S(E) [E sin(D/2) - m cos(Sig)],
    C = cos(pL), S = sin(pL)/p  (analytic through the gap, p -> iq),
root-refined by brentq; heat-kernel eta with Richardson extrapolation in t.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq
from ch08_bag_spectrum import rays, find_spectrum

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

_trapezoid = getattr(np, "trapezoid", None) or np.trapz


# ---------------------------------------------------------------- spectrum
def phi_vec(Es, m, L, th0, thL):
    """Entire-function master-equation residual, vectorized over E."""
    D, S = thL - th0, 0.5 * (th0 + thL)
    s2 = Es * Es - m * m
    r = np.sqrt(np.abs(s2))
    rg = np.where(s2 < 0, r, 0.0)          # gap argument only (avoids cosh overflow)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.where(s2 >= 0, np.cos(r * L), np.cosh(rg * L))
        Sf = np.where(s2 > 0, np.sin(r * L) / r,
                      np.where(s2 < 0, np.sinh(rg * L) / np.where(rg > 0, rg, 1.0), L))
    return -np.cos(D / 2) * C + Sf * (Es * np.sin(D / 2) - m * np.cos(S))


def chiral_spectrum(m, L, th0, thL, Emax):
    """Complete spectrum in [-Emax, Emax] from sign changes of phi (brentq-refined)."""
    n = int(40 * (2 * Emax) * L / np.pi) + 4000
    Es = np.linspace(-Emax, Emax, n)
    vals = phi_vec(Es, m, L, th0, thL)
    idx = np.where(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    f = lambda E: float(phi_vec(np.array([E]), m, L, th0, thL)[0])
    return np.array([brentq(f, Es[i], Es[i + 1], xtol=1e-13) for i in idx])


def qvac(m, L, th0, thL, Emax=350.0, ts=(0.08, 0.04, 0.02)):
    """Q_vac = -eta/2, heat-kernel regulated, Richardson-extrapolated (linear in t)."""
    spec = chiral_spectrum(m, L, th0, thL, Emax)
    etas = [np.sum(np.sign(spec) * np.exp(-t * np.abs(spec))) for t in ts]
    eta0 = etas[-1] + (etas[-1] - etas[-2])
    return -0.5 * eta0, len(spec)


# ---------------------------------------------------------------- modes/overlaps
def mode_functions(m, L, th0, thL, Emax, ngrid=3001):
    """Energies + normalized spinor wavefunctions on a grid (closed form, no loop:
    psi(x) = cos(px) v0 + sin(px)/p K v0, K the transfer generator)."""
    Es = chiral_spectrum(m, L, th0, thL, Emax)
    x = np.linspace(0.0, L, ngrid)
    v0, _ = rays(th0, thL)
    K = lambda E: np.array([[0.0, 1j * (E + m)], [1j * (E - m), 0.0]])
    psis = np.empty((len(Es), 2, ngrid), dtype=complex)
    for i, E in enumerate(Es):
        p = np.sqrt(complex(E * E - m * m))
        Kv = K(E) @ v0
        if abs(p) < 1e-14:
            psi = v0[:, None] + np.outer(Kv, x)
        else:
            psi = np.outer(v0, np.cos(p * x)) + np.outer(Kv, np.sin(p * x) / p)
        nrm = np.sqrt(_trapezoid((np.abs(psi) ** 2).sum(axis=0), x).real)
        psis[i] = psi / nrm
    return Es, psis, x


def overlap_matrix(E_old, psi_old, x_old, E_new, psi_new, x_new):
    """a_{kn} = int_0^{L_old} psi'_k^dag psi_n dx (old modes vanish beyond L_old)."""
    a = np.zeros((len(E_new), len(E_old)), dtype=complex)
    for k in range(len(E_new)):
        pk = np.vstack([np.interp(x_old, x_new, psi_new[k, c].real) +
                        1j * np.interp(x_old, x_new, psi_new[k, c].imag)
                        for c in range(2)])
        integ = np.einsum("cx,ncx->nx", pk.conj(), psi_old)
        a[k] = _trapezoid(integ, x_old, axis=1)
    return a


def dq_truncated(E_old, E_new, a, N):
    """Same-branch and cross-branch truncated forms, N lowest |E| per branch."""
    def pick(E, sign):
        idx = np.where(np.sign(E) == sign)[0]
        return idx[np.argsort(np.abs(E[idx]))][:N]
    op, om = pick(E_old, +1), pick(E_old, -1)
    npos, nneg = pick(E_new, +1), pick(E_new, -1)
    dq_cross = (np.abs(a[np.ix_(npos, om)]) ** 2).sum() - (np.abs(a[np.ix_(nneg, op)]) ** 2).sum()
    dq_same = (np.abs(a[np.ix_(npos, op)]) ** 2).sum() - (np.abs(a[np.ix_(nneg, om)]) ** 2).sum()
    return dq_cross, dq_same


# ---------------------------------------------------------------- main
def main():
    m = 1.0

    # ---- A. boundary rays vs eigendecomposition
    print("=== A. boundary rays (12.3)-(12.4) vs eigendecomposition ===")
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2m = np.array([[0, -1j], [1j, 0]], dtype=complex)
    worst = 0.0
    for th in (-2.0, -0.7, 0.0, 0.5, 1.3, 2.0):
        eth = np.cos(th) * np.eye(2) + 1j * np.sin(th) * s1     # e^{i th g5}
        for sgn, which in ((+1, "L"), (-1, "0")):
            Mop = sgn * (s2m @ eth)                             # -i n.g e^{i th g5}
            w, V = np.linalg.eig(Mop)
            v = V[:, np.argmin(np.abs(w - 1))]
            ray = rays(th, th)[0 if which == "0" else 1]
            # compare up to phase
            ph = v @ ray.conj()
            worst = max(worst, np.linalg.norm(v * (ph.conj() / abs(ph)) - ray))
    print(f"  max |ray - eigvec| over 6 angles x 2 walls = {worst:.2e}")

    # ---- B. master equation vs transfer-matrix solver
    print("=== B. master equation (12.7) vs transfer-matrix solver ===")
    worst = 0.0
    for (t0, tL) in ((0.0, 0.0), (0.0, 0.5), (0.0, 1.5), (-0.7, 1.1), (-2.0, 2.0)):
        ex = chiral_spectrum(m, 1.0, t0, tL, Emax=12.0)
        tm = find_spectrum(m, 0.0, 1.0, t0, tL, Emax=12.0, nE=8001)
        k = min(len(ex), len(tm))
        worst = max(worst, max(abs(a - b) for a, b in
                               zip(np.sort(ex)[:k], np.sort(tm)[:k])))
    print(f"  max |master-eq roots - TM solver| over 5 angle pairs = {worst:.2e}")

    # ---- C. (delta,delta) real mass == bulk complex mass
    print("=== C. (delta,delta) walls + real mass == bulk complex mass ===")
    worst = 0.0
    for d in (0.3, 0.8):
        ex = chiral_spectrum(m, 1.0, d, d, Emax=12.0)
        tm = find_spectrum(m * np.cos(d), m * np.sin(d), 1.0, 0.0, 0.0,
                           Emax=12.0, nE=8001)
        k = min(len(ex), len(tm))
        worst = max(worst, max(abs(a - b) for a, b in
                               zip(np.sort(ex)[:k], np.sort(tm)[:k])))
    print(f"  max |spec diff| = {worst:.2e}")

    # ---- D. fractional charge law at Sigma = 0, mL = 6
    print("=== D. Fractional Charge Law: Q_vac vs -Delta/2pi (Sigma=0, mL=6) ===")
    L6 = 6.0
    Deltas = np.array([0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8])
    qs, nlev = [], 0
    for D in Deltas:
        q, nl = qvac(m, L6, -D / 2, +D / 2)
        qs.append(q)
        nlev = max(nlev, nl)
    qs = np.array(qs)
    for D, q in zip(Deltas, qs):
        print(f"  Delta={D:.1f}: Q_vac = {q:+.4f}   -Delta/2pi = {-D / (2 * np.pi):+.4f}")
    print(f"  (complete spectra: ~{nlev} levels below E_max = 350)")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    Dfine = np.linspace(0, 3.0, 200)
    ax.plot(Dfine, -Dfine / (2 * np.pi), "-", color="#555", lw=1.2,
            label=r"$-\Delta/2\pi$")
    ax.plot(Deltas, qs, "o", ms=5, label=r"complete-spectrum $Q_{\rm vac}$")
    ax.set_xlabel(r"wall mismatch $\Delta$")
    ax.set_ylabel(r"$Q_{\rm vac}$")
    ax.set_title(r"The fractional law: $Q_{\rm vac} = -\Delta/2\pi$ ($\Sigma=0$, $mL=6$)")
    ax.legend()
    fig.savefig(FIGS / "ch12_fig1_qvac_law.png")
    print(f"figure -> {FIGS / 'ch12_fig1_qvac_law.png'}")

    # ---- E. the crossing: plateaus and the unit jump
    print("=== E. crossing at theta=(-2,+2): plateaus and jump ===")
    th = (-2.0, 2.0)
    Lstar = np.arctanh(-np.cos(2.0) / np.cos(0.0)) / m
    print(f"  predicted mL* = {Lstar:.5f}")
    Ls = np.concatenate([np.linspace(0.30, 0.435, 5), np.linspace(0.452, 0.60, 5)])
    qL = np.array([qvac(m, L, *th)[0] for L in Ls])
    below, above = qL[Ls < Lstar], qL[Ls > Lstar]
    print(f"  plateau below L*: {below.mean():+.4f}   "
          f"(law: -4/2pi mod 1 = {(-4 / (2 * np.pi)) % 1:+.4f})")
    print(f"  plateau above L*: {above.mean():+.4f}")
    print(f"  jump across L*  : {above.mean() - below.mean():+.5f}   (exact: -1)")

    # gap-level track for the figure
    Lt = np.linspace(0.30, 0.60, 61)
    gap = []
    for L in Lt:
        sp = chiral_spectrum(m, L, *th, Emax=6.0)
        sp = sp[np.abs(sp) < 0.999 * m]
        gap.append(sp[np.argmin(np.abs(sp))] if len(sp) else np.nan)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.2, 5.6), sharex=True)
    ax1.plot(Lt, gap, "o-", ms=3)
    ax1.axhline(0, color="#888", lw=0.7)
    ax1.axvline(Lstar, color="k", ls=":", lw=0.9)
    ax1.set_ylabel(r"gap level $E(L)$")
    ax1.set_title(r"The pump's escapement: the crossing at $mL^* = 0.44302$")
    Lq = np.sort(np.concatenate([Ls, [Lstar - 0.004, Lstar + 0.004]]))
    qq = np.array([qvac(m, L, *th)[0] for L in Lq])
    ax2.plot(Lq, qq, "o-", ms=3)
    ax2.axvline(Lstar, color="k", ls=":", lw=0.9)
    ax2.set_xlabel(r"box size $L$")
    ax2.set_ylabel(r"$Q_{\rm vac}(L)$")
    fig.savefig(FIGS / "ch12_fig2_crossing.png")
    print(f"figure -> {FIGS / 'ch12_fig2_crossing.png'}")

    # ---- F. truncated Bogoliubov forms: pump vs control
    print("=== F. truncated forms, pump (theta=(-2,+2)) and control (0.7,0.7) ===")
    Lo, Ln, Ncut = 0.30, 0.70, 75
    for label, t0, tL in (("pump   ", -2.0, 2.0), ("control", 0.7, 0.7)):
        Emax_old = np.sqrt((Ncut * np.pi / Lo) ** 2 + m * m) + 5
        Emax_new = np.sqrt((Ncut * np.pi / Ln) ** 2 + m * m) + 5
        Eo, po, xo = mode_functions(m, Lo, t0, tL, Emax_old, ngrid=3001)
        En, pn, xn = mode_functions(m, Ln, t0, tL, Emax_new, ngrid=5001)
        a = overlap_matrix(Eo, po, xo, En, pn, xn)
        for N in (25, 50, 75):
            dqc, dqs = dq_truncated(Eo, En, a, N)
            print(f"  {label} N={N:3d}: dQ_cross = {dqc:+.4f}   dQ_same = {dqs:+.4f}")
    print("  [spectral-flow answers: pump = +1 exactly, control = 0 exactly;")
    print("   see ch12_pump_control.py — the truncated forms measure neither]")


if __name__ == "__main__":
    main()
