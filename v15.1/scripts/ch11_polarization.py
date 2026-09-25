#!/usr/bin/env python3
"""
ch11_polarization.py — what truncated pair-creation sums measure (Ch. 11).

Sudden expansion L -> sL of the 1+1D bulk-CP-mass MIT bag. All spectra are the
exact roots of tan(pL) = -p/(m cos delta) with E = +-sqrt(p^2 + m^2), so the
+-E pairing is analytic here -- BY CONSTRUCTION, not a test; the honest
solver-level mirror check is ch10_mirror_spectrum.py. Eigenfunctions in closed
form: psi(x) = [cos(px) + sin(px) K(E)/p] v0 with K^2 = -p^2, v0 = (1, -i)/sqrt2.

default   : --repro + --conv
--repro   : historical point s = 2, mL = 1.085, delta = 0.015, N = 18 per
            branch; prints A = dQ/sin(2 delta) for both truncated forms (11.1)
--conv    : benchmark point delta = pi/4, s = 1.2, mL = 1; both truncated
            forms at N = 6..85, the Polarization Identity value
            (1/2) sum_k sgn(E'_k) w_k (heat kernel, Richardson t -> 0),
            whole-box eta; Fig. 11.2
--oddness : delta -> -delta antisymmetry of both truncated forms (Fig. 11.3)
--defect  : App. A (A.9) truncation defect: the old-indexed defect decays as
            1/N; the wrong-way-round (new-indexed) defect saturates at
            1 - w_k, exactly as (A.8) demands
--all     : everything
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
try:
    plt.style.use(HERE / "style.mplstyle")
except OSError:
    pass

_trapz = getattr(np, "trapezoid", np.trapz)


# ---------------------------------------------------------------- modes
def tan_roots(m, delta, L, Emax):
    """All p > 0 with tan(pL) = -p/(m cos delta) and sqrt(p^2 + m^2) <= Emax."""
    mR = m * np.cos(delta)
    f = lambda p: np.tan(p * L) + p / mR
    ps, j = [], 1
    while np.sqrt(((j - 1) * np.pi / L) ** 2 + m ** 2) <= Emax:
        lo = (j - 0.5) * np.pi / L + 1e-12
        hi = j * np.pi / L - 1e-12
        if f(lo) * f(hi) < 0:
            ps.append(brentq(f, lo, hi, xtol=1e-14))
        j += 1
    ps = np.array(ps)
    return ps[np.sqrt(ps ** 2 + m ** 2) <= Emax]


def modes(m, delta, L, Emax, x, x_eval=None):
    """Complete spectrum (exact +-E pairs) and closed-form eigenfunctions.

    Normalized on the grid x (which must span [0, L]); sampled on x_eval
    (default: x). Returns (E ascending, psi[nmode, 2, len(x_eval)]).
    """
    mR, mI = m * np.cos(delta), m * np.sin(delta)
    ps = tan_roots(m, delta, L, Emax)
    E = np.concatenate([np.sqrt(ps ** 2 + m ** 2), -np.sqrt(ps ** 2 + m ** 2)])
    pp = np.concatenate([ps, ps])
    order = np.argsort(E)
    E, pp = E[order], pp[order]
    xe = x if x_eval is None else x_eval
    v0 = np.array([1.0, -1.0j]) / np.sqrt(2.0)
    psis = np.empty((len(E), 2, len(xe)), dtype=complex)
    for i in range(len(E)):
        K = np.array([[-mI, 1j * (E[i] + mR)], [1j * (E[i] - mR), mI]])
        Kv = K @ v0

        def raw(xx):
            return (np.cos(pp[i] * xx)[None, :] * v0[:, None]
                    + (np.sin(pp[i] * xx) / pp[i])[None, :] * Kv[:, None])

        nrm = np.sqrt(_trapz((np.abs(raw(x)) ** 2).sum(axis=0), x))
        psis[i] = raw(xe) / nrm
    return E, psis


# ---------------------------------------------------------------- overlaps
def overlap_matrix(psi_new, psi_old, x):
    """alpha_{kn} = int_0^L psi'_k^dag psi_n dx (both sampled on the old grid)."""
    w = np.full(len(x), x[1] - x[0])
    w[0] *= 0.5
    w[-1] *= 0.5
    PN = (psi_new.conj() * w[None, None, :]).reshape(len(psi_new), -1)
    PO = psi_old.reshape(len(psi_old), -1)
    return PN @ PO.T


def dq_forms(E_old, E_new, a, N):
    """Both truncated charge forms of (11.1), N lowest-|E| modes per branch."""
    def pick(E, sign):
        idx = np.where(np.sign(E) == sign)[0]
        return idx[np.argsort(np.abs(E[idx]))][:N]
    op, om = pick(E_old, +1), pick(E_old, -1)
    kp, km = pick(E_new, +1), pick(E_new, -1)
    cross = ((np.abs(a[np.ix_(kp, om)]) ** 2).sum()
             - (np.abs(a[np.ix_(km, op)]) ** 2).sum())
    same = ((np.abs(a[np.ix_(kp, op)]) ** 2).sum()
            - (np.abs(a[np.ix_(km, om)]) ** 2).sum())
    return cross, same


def build(m, delta, L, s, Nmax, ngrid_old=6001, ngrid_new=7001):
    """Old/new spectra, overlap matrix, and old-region weights w_k."""
    Emax_old = np.sqrt((Nmax * np.pi / L) ** 2 + m ** 2) + 1
    Emax_new = np.sqrt((Nmax * np.pi / (s * L)) ** 2 + m ** 2) + 1
    xo = np.linspace(0, L, ngrid_old)
    xn = np.linspace(0, s * L, ngrid_new)
    Eo, po = modes(m, delta, L, Emax_old, xo)
    En, pn = modes(m, delta, s * L, Emax_new, xn, x_eval=xo)   # normed on [0, L']
    a = overlap_matrix(pn, po, xo)
    w = np.array([_trapz((np.abs(p) ** 2).sum(axis=0), xo) for p in pn])
    return Eo, En, a, w


# ---------------------------------------------------------------- jobs
def repro():
    m, L, s, delta, N = 1.085, 1.0, 2.0, 0.015, 18
    Eo, En, a, _ = build(m, delta, L, s, Nmax=N, ngrid_old=4001, ngrid_new=4001)
    dqc, dqs = dq_forms(Eo, En, a, N)
    A_s, A_c = dqs / np.sin(2 * delta), dqc / np.sin(2 * delta)
    print("[repro] s=2, mL=1.085, delta=0.015, N=18 per branch:")
    print(f"        A_same  = dQ_same /sin(2d) = {A_s:+.4f}")
    print(f"        A_cross = dQ_cross/sin(2d) = {A_c:+.4f}"
          f"   (inter-form split {100 * (A_s - A_c) / A_s:.0f}%)")


def conv():
    m, L, s, delta = 1.0, 1.0, 1.2, np.pi / 4
    Eo, En, a, w = build(m, delta, L, s, Nmax=90)
    print(f"[conv] benchmark delta=pi/4, s=1.2, mL=1 "
          f"({len(Eo)} old / {len(En)} new modes, complete analytic spectra)")
    Ns = (6, 10, 14, 18, 26, 38, 54, 70, 85)
    rows = []
    for N in Ns:
        dqc, dqs = dq_forms(Eo, En, a, N)
        rows.append((N, dqc, dqs))
        print(f"[conv] N={N:3d}:  dQ_cross = {dqc:+.5f}   dQ_same = {dqs:+.5f}")
    ts = (0.04, 0.02, 0.01)
    pol = [0.5 * float(np.sum(np.sign(En) * np.exp(-t * np.abs(En)) * w)) for t in ts]
    t1, t2 = ts[-2], ts[-1]
    extrap = (pol[-1] * t1 - pol[-2] * t2) / (t1 - t2)
    print(f"[conv] polarization (1/2) sum_k sgn(E'_k) w_k at t={ts}: "
          f"{np.round(pol, 5)}  ->  Richardson t->0: {extrap:+.5f}")
    eta = [float(np.sum(np.sign(En) * np.exp(-t * np.abs(En)))) for t in ts]
    print(f"[conv] whole-box eta(L'): {np.round(eta, 12)}  "
          "(zero BY CONSTRUCTION: analytic +-E pairing; honest solver test: "
          "ch10_mirror_spectrum.py)")

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    Narr = np.array([r[0] for r in rows])
    ax.plot(Narr, [r[1] for r in rows], "o-", ms=4, label=r"$\Delta Q_{\rm cross}$")
    ax.plot(Narr, [r[2] for r in rows], "s-", ms=4, label=r"$\Delta Q_{\rm same}$")
    ax.axhline(extrap, color="k", lw=0.9, ls="-",
               label=fr"polarization integral ${extrap:+.4f}$")
    ax.axhline(0, color="#888888", lw=0.9, ls=":", label="net charge (theorem): 0")
    ax.set_xlabel(r"modes per branch $N$")
    ax.set_ylabel(r"truncated charge forms")
    ax.set_title("Two truncated forms closing on the polarization integral")
    ax.legend(fontsize=8)
    axi = ax.inset_axes([0.55, 0.18, 0.4, 0.32])
    gap = np.abs(np.array([r[2] - r[1] for r in rows]))
    axi.loglog(Narr, gap, "o-", ms=3)
    axi.loglog(Narr, gap[0] * Narr[0] / Narr, "--", lw=0.8)
    axi.set_title(r"$|\Delta Q_{\rm same} - \Delta Q_{\rm cross}|$ vs $1/N$", fontsize=7)
    axi.tick_params(labelsize=6)
    fig.savefig(FIGS / "ch11_fig2_convergence.png")
    print("figure -> ch11_fig2_convergence.png")


def oddness():
    m, L, s, N = 1.0, 1.0, 1.2, 18
    deltas = np.linspace(0.1, 1.2, 7)
    sums_c, sums_s = [], []
    for d in deltas:
        r = []
        for sgn in (+1, -1):
            Eo, En, a, _ = build(m, sgn * d, L, s, Nmax=N,
                                 ngrid_old=3001, ngrid_new=3001)
            r.append(dq_forms(Eo, En, a, N))
        sums_c.append(r[0][0] + r[1][0])
        sums_s.append(r[0][1] + r[1][1])
    print(f"[oddness] max |dQ(+d) + dQ(-d)| over delta in [0.1, 1.2], N=18: "
          f"cross {np.max(np.abs(sums_c)):.1e}, same {np.max(np.abs(sums_s)):.1e}")
    print("[oddness] (exact oddness at any symmetric truncation: C maps the "
          "+delta problem to the -delta problem with E -> -E)")
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    ax.plot(deltas, sums_c, "o-", ms=4, label=r"$\Delta Q_{\rm cross}(+\delta) + \Delta Q_{\rm cross}(-\delta)$")
    ax.plot(deltas, sums_s, "s-", ms=4, label=r"$\Delta Q_{\rm same}(+\delta) + \Delta Q_{\rm same}(-\delta)$")
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel(r"CP phase $\delta$")
    ax.set_ylabel("oddness residual")
    ax.set_title("ECPT oddness of the truncated forms: zero at machine precision")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch11_fig3_ecpt_oddness.png")
    print("figure -> ch11_fig3_ecpt_oddness.png")


def defect():
    m, L, s, delta = 1.0, 1.0, 1.2, np.pi / 4
    Eo, En, a, w = build(m, delta, L, s, Nmax=60, ngrid_old=4001, ngrid_new=5001)
    nlow = np.argsort(np.abs(Eo))[:6]
    print("[defect] old-indexed defect (A.9): max diag of "
          "delta_nn' - sum_{|k|<=N} a_kn a*_kn' over the 6 lowest old modes")
    for N in (10, 20, 40, 58):
        kidx = np.argsort(np.abs(En))[:2 * N]
        G = a[np.ix_(kidx, nlow)]
        D = np.eye(6) - (G.conj().T @ G).real
        print(f"[defect]   N={N:3d}: {np.max(np.abs(np.diag(D))):.4f}")
    klow = np.argsort(np.abs(En))[:6]
    nidx = np.argsort(np.abs(Eo))[:2 * 58]
    G = a[np.ix_(klow, nidx)]
    D = np.eye(6) - (G @ G.conj().T).real
    print(f"[defect] wrong-direction (new-indexed) defect at N=58, diag: "
          f"{np.round(np.diag(D), 4)}")
    print(f"[defect] saturation 1 - w_k demanded by (A.8):          "
          f"{np.round(1 - w[klow], 4)}")


if __name__ == "__main__":
    args = set(sys.argv[1:])
    jobs = args if args else {"--repro", "--conv"}
    if "--all" in jobs:
        jobs = {"--repro", "--conv", "--oddness", "--defect"}
    if "--repro" in jobs:
        repro()
    if "--conv" in jobs:
        conv()
    if "--oddness" in jobs:
        oddness()
    if "--defect" in jobs:
        defect()
