#!/usr/bin/env python3
"""
ch25_graviton.py — the graviton test (Ch. 25).

Massless scalar on a periodic 63^3 lattice with static metric g(x) depending on
one coordinate: the Laplace-Beltrami operator block-diagonalizes over the 63^2
exact transverse momenta (3969 blocks of dimension 63; spectral derivatives on
an odd grid, no Nyquist mode).  Vacuum energy with the covariant proper-time
regulator E(s) = (1/2) sum_n omega_n exp(-s omega_n).  Each configuration is
diagonalized ONCE; E(s) is then evaluated at every s from the stored spectrum.

Observables (Ch. 25.1): Pi(q) = 4 dE / (eps^2 V) per sector; quartic fits
Pi = Pi0 + kappa q^2 + beta q^4 over k = 1..5; the induced gradient stiffness
kappa compared to the parameter-free one-loop coefficient C(s) = 1/(48 pi^2 s^2)
(eq. 25.5).  Prediction (25.6): kappa_TT+ = kappa_TTx = -C, kappa_conf = +C,
kappa_diffeo = 0.

Unit tests (gate, App. C.4): flat-space block spectrum vs omega^2 = k.k exact to
machine; uniform-shear spectrum vs the exactly solvable omega^2 = g^{ij} k_i k_j.
Controls (App. C.1.5): the exact-diffeomorphism pullback g = J^T J of flat space
runs in every sweep, including the raw (s = 0) sum that exhibits the covariance
window of Ch. 25.2.

SIGN AUDIT (Ch. 25.5): the same spectra are also evaluated in the 4D-covariant
proper-time scheme -- the exact T=0 mode-sum reduction of
Gamma_E = -1/2 Int ds/s Tr exp(-s(-box)):

    E_4D(s0) = 1/2 sum_n w_n erfc(sqrt(s0) w_n)
               - (1/(2 sqrt(pi s0))) sum_n exp(-s0 w_n^2),

whose Weyl-expansion prediction is E_4D contains C_4D(s0) Int sqrt(g) R^(3) with
C_4D(s0) = -1/(192 pi^2 s0): the OPPOSITE overall sign to the exp(-s w) scheme.
The quadratically divergent EH coefficient is regulator-scheme-dependent, sign
included; only the DeWitt pattern (TT vs conformal opposite, diffeo null, pure
q^2, polarization equality) is scheme-robust.  In this scheme kappa_TT = -C_4D
> 0 (stable waves) and kappa_conf < 0 (compression unstable: attraction), the
signs that match S_EH -- and the diffeo control lands at machine zero.

Writes ../data/ch25_graviton_results.json (table of record for the chapter)
and regenerates Figures 25.2 (covariance window) and 25.3 (kappa convergence).
"""
import json
import pathlib
import time

import numpy as np
from scipy.special import erfc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
DATA = HERE.parent / "data"
DATA.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

N = 63                     # lattice side (odd: no Nyquist ambiguity)
EPS = 0.08                 # perturbation amplitude
KS = (1, 2, 3, 4, 5)       # wavenumbers along x
S_LIST = (2.0, 2.5, 3.0)   # proper-time regulator values (lattice units)
S_WINDOW = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # covariance-window scan (diffeo)
S0_LIST = (4.0, 9.0)       # 4D proper-time values for the sign audit (Ch. 25.5)


def C_of_s(s):
    """Universal one-loop EH coefficient in the exp(-s w) scheme, eq. (25.5)."""
    return 1.0 / (48.0 * np.pi**2 * s**2)


def C_4d(s0):
    """EH coefficient in the 4D proper-time scheme (sign audit, Ch. 25.5)."""
    return -1.0 / (192.0 * np.pi**2 * s0)


def spectral_D(n):
    """Antisymmetrized spectral derivative on n points (n odd -> exact)."""
    F = np.fft.fft(np.eye(n), axis=0)
    kk = 2j * np.pi * np.fft.fftfreq(n)
    D = np.real(np.fft.ifft(kk[:, None] * F, axis=0))
    return 0.5 * (D - D.T)


def metric_arrays(g, n):
    """g: dict (a,b)->array over x (upper triangle). Return inverse and sqrt(det)."""
    G = np.zeros((n, 3, 3))
    for a in range(3):
        for b in range(3):
            G[:, a, b] = g[(min(a, b), max(a, b))]
    A = np.linalg.inv(G)
    sg = np.sqrt(np.linalg.det(G))
    inv = {}
    for a in range(3):
        for b in range(a, 3):
            inv[(a, b)] = A[:, a, b]
    return inv, sg


def spectrum(gfun):
    """All 63^3 frequencies of sqrt(Delta_g) for a metric depending on x only.

    Delta_g = (1/sqrt g) d_i (sqrt g g^{ij} d_j): discretized symmetrically
    (similarity-transformed by sqrt(sqrt g) to a Hermitian block), so the
    operator is self-adjoint per block by construction (App. C.4).
    """
    x = np.arange(N)
    A, sg = metric_arrays(gfun(x), N)
    w = {k: A[k] * sg for k in A}
    rs = 1.0 / np.sqrt(sg)
    Dx = spectral_D(N)
    kts = 2 * np.pi * np.fft.fftfreq(N)
    Kxx = Dx.T @ (w[(0, 0)][:, None] * Dx)
    Cxy = 1j * (Dx.T @ np.diag(w[(0, 1)]) - np.diag(w[(0, 1)]) @ Dx)
    Cxz = 1j * (Dx.T @ np.diag(w[(0, 2)]) - np.diag(w[(0, 2)]) @ Dx)
    om_all = np.empty((N, N, N))
    for iy, ky in enumerate(kts):
        for iz, kz in enumerate(kts):
            K = (Kxx
                 + np.diag(ky**2 * w[(1, 1)] + kz**2 * w[(2, 2)]
                           + 2 * ky * kz * w[(1, 2)])
                 + ky * Cxy + kz * Cxz)
            H = (K * rs[:, None]) * rs[None, :]
            om_all[iy, iz] = np.sqrt(np.clip(np.linalg.eigvalsh(H), 0.0, None))
    return om_all.ravel()


def E_of_s(om, s):
    return 0.5 * np.sum(om * np.exp(-s * om))


def E_4d(om, s0):
    """Exact T=0 mode-sum reduction of the 4D proper-time effective action."""
    return (0.5 * np.sum(om * erfc(np.sqrt(s0) * om))
            - np.sum(np.exp(-s0 * om**2)) / (2.0 * np.sqrt(np.pi * s0)))


def gmode(mode, eps, q):
    """Test configurations of Ch. 25.1 (all waves along x)."""
    def f(x):
        c = np.cos(q * x)
        one = np.ones_like(x, float)
        z = np.zeros_like(x, float)
        g = {(a, b): (one.copy() if a == b else z.copy())
             for a in range(3) for b in range(a, 3)}
        if mode == "conf":                     # P = identity
            for a in range(3):
                g[(a, a)] = 1 + eps * c
        elif mode == "TT+":                    # P = diag(0, 1, -1)
            g[(1, 1)] = 1 + eps * c
            g[(2, 2)] = 1 - eps * c
        elif mode == "TTx":                    # P = E_yz + E_zy
            g[(1, 2)] = eps * c
        elif mode == "diffeo":                 # exact pullback of flat space by
            g[(0, 0)] = 1 + (eps * c)**2       #   y -> y + (eps/q) sin(qx)
            g[(0, 1)] = eps * c
        return g
    return f


# ----------------------------------------------------------------------------
# Unit tests (gate everything downstream, App. C.4)
# ----------------------------------------------------------------------------
def unit_tests():
    kts = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY, KZ = np.meshgrid(kts, kts, kts, indexing="ij")

    # flat space: omega^2 = k.k, exact
    om_flat = np.sort(spectrum(gmode("base", 0, 0)))
    om_exact = np.sort(np.sqrt(KX**2 + KY**2 + KZ**2).ravel())
    err_flat = np.max(np.abs(om_flat - om_exact)) / np.max(om_exact)

    # uniform shear: constant g with g_yz = gamma; exactly solvable by
    # coordinate change: omega^2 = g^{ij} k_i k_j over the same lattice momenta
    gamma = 0.3

    def gshear(x):
        one = np.ones_like(x, float)
        z = np.zeros_like(x, float)
        return {(0, 0): one, (1, 1): one, (2, 2): one.copy(),
                (0, 1): z, (0, 2): z.copy(), (1, 2): gamma * one}

    om_shear = np.sort(spectrum(gshear))
    det = 1 - gamma**2       # inverse of [[1,gamma],[gamma,1]] in the yz block
    om2 = KX**2 + (KY**2 + KZ**2 - 2 * gamma * KY * KZ) / det
    om_shear_exact = np.sort(np.sqrt(om2).ravel())
    err_shear = np.max(np.abs(om_shear - om_shear_exact)) / np.max(om_shear_exact)

    print(f"# unit test  flat-space spectrum : max rel err = {err_flat:.1e}")
    print(f"# unit test  uniform-shear       : max rel err = {err_shear:.1e}")
    assert err_flat < 1e-12 and err_shear < 1e-7, "unit tests FAILED — abort"
    return err_flat, err_shear


# ----------------------------------------------------------------------------
# Production sweep
# ----------------------------------------------------------------------------
def main():
    t0 = time.time()
    err_flat, err_shear = unit_tests()
    V = N**3

    om_base = spectrum(gmode("base", 0, 0))
    E0 = {s: E_of_s(om_base, s) for s in set(S_LIST) | set(S_WINDOW)}
    E0_4d = {s0: E_4d(om_base, s0) for s0 in S0_LIST}

    sectors = ("conf", "TT+", "TTx", "diffeo")
    Pi = {s: {m: [] for m in sectors} for s in S_LIST}       # (q^2, Pi) pairs
    Pi4 = {s0: {m: [] for m in sectors} for s0 in S0_LIST}   # 4D-PT scheme
    window = {m: {} for m in ("diffeo", "TT+")}              # Pi(s) at k=1

    for m in sectors:
        for k in KS:
            q = 2 * np.pi * k / N
            om = spectrum(gmode(m, EPS, q))
            for s in S_LIST:
                dE = E_of_s(om, s) - E0[s]
                Pi[s][m].append((float(q * q), float(4 * dE / (EPS**2 * V))))
            for s0 in S0_LIST:
                dE = E_4d(om, s0) - E0_4d[s0]
                Pi4[s0][m].append((float(q * q), float(4 * dE / (EPS**2 * V))))
            if m in window and k == 1:
                for s in S_WINDOW:
                    dE = E_of_s(om, s) - E0[s]
                    window[m][s] = float(4 * dE / (EPS**2 * V))
            print(f"#   sector {m:6s} k={k}  done ({time.time()-t0:6.1f}s)")

    # quartic fits Pi = Pi0 + kappa q^2 + beta q^4 (degree 2 in q^2); the fit is
    # over 5 deterministic points, so the quality measure is the rms residual of
    # the fit (not a statistical covariance), reported alongside kappa/C
    results = {"N": N, "eps": EPS, "ks": list(KS),
               "unit_tests": {"flat": err_flat, "shear": err_shear},
               "covariance_window_k1": {m: {str(s): window[m][s] for s in S_WINDOW}
                                        for m in window},
               "table": {}}
    print("\n#  s     C(s)       sector   kappa        kappa/C     fit rms     "
          "Pi0         beta")
    for s in S_LIST:
        C = C_of_s(s)
        row = {}
        for m in sectors:
            R = np.array(Pi[s][m])
            coef = np.polyfit(R[:, 0], R[:, 1], 2)
            beta, kap, Pi0 = coef
            rms = float(np.sqrt(np.mean((np.polyval(coef, R[:, 0]) - R[:, 1]) ** 2)))
            row[m] = {"kappa": float(kap), "kappa_over_C": float(kap / C),
                      "fit_rms": rms, "Pi0": float(Pi0),
                      "beta": float(beta), "points": R.tolist()}
            print(f"  {s:.1f}  {C:.3e}  {m:6s}  {kap:+.3e}  {kap/C:+9.4f}  "
                  f"{rms:9.1e}  {Pi0:+.3e}  {beta:+.3e}")
        results["table"][str(s)] = row

    # ---- sign audit: the same spectra in the 4D proper-time scheme ----
    results["table_4dpt"] = {}
    print("\n# SIGN AUDIT (Ch. 25.5) — 4D proper-time scheme, C_4D(s0) = "
          "-1/(192 pi^2 s0):")
    print("#  s0    C_4D(s0)    sector   kappa        kappa/C_4D")
    for s0 in S0_LIST:
        C4 = C_4d(s0)
        row = {}
        for m in sectors:
            R = np.array(Pi4[s0][m])
            coef = np.polyfit(R[:, 0], R[:, 1], 2)
            beta4, kap4, Pi04 = coef
            row[m] = {"kappa": float(kap4), "kappa_over_C4d": float(kap4 / C4),
                      "Pi0": float(Pi04), "beta": float(beta4),
                      "points": R.tolist()}
            print(f"  {s0:.1f}  {C4:+.3e}  {m:6s}  {kap4:+.3e}  {kap4/C4:+9.4f}")
        results["table_4dpt"][str(s0)] = row
    print("# scheme-robust: DeWitt pattern (TT vs conf opposite, diffeo null).")
    print("# scheme-dependent: the overall sign — kappa_TT > 0 (stable waves),")
    print("# kappa_conf < 0 (attraction) in the covariant-matching scheme.")

    print("\n# covariance window (k=1): Pi(s), diffeo control vs TT+ signal")
    for s in S_WINDOW:
        print(f"    s = {s:.1f} : Pi_diffeo = {window['diffeo'][s]:+.3e}   "
              f"Pi_TT+ = {window['TT+'][s]:+.3e}")

    out = DATA / "ch25_graviton_results.json"
    json.dump(results, open(out, "w"), indent=1)
    print(f"\n# wrote {out}  ({time.time()-t0:.1f}s total)")

    # ---- Figure 25.2: the covariance window ----
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    sw = np.array(S_WINDOW[1:])                     # skip s=0 on the log plot
    ax.semilogy(sw, [abs(window["diffeo"][s]) for s in sw], "^-", color="steelblue",
                label=r"$|\Pi_{\rm diffeo}|$ (control $\to$ 0)")
    ax.semilogy(sw, [abs(window["TT+"][s]) for s in sw], "s-", color="seagreen",
                label=r"$|\Pi_{TT+}|$ (signal $\propto 1/s^2$)")
    ax.semilogy(sw, [4 * C_of_s(s) * (2 * np.pi / N) ** 2 / 4 for s in sw], "k:",
                label=r"$C(s)\,q_1^2$ scale")
    ax.axhline(abs(window["diffeo"][0.0]), color="steelblue", ls="--", lw=0.8)
    ax.text(2.35, abs(window["diffeo"][0.0]) * 1.3,
            f"raw lattice sum: {window['diffeo'][0.0]:+.2f}", fontsize=8, color="steelblue")
    ax.set_xlabel("proper-time regulator $s$ (lattice units)")
    ax.set_ylabel(r"$|\Pi(q_1)|$")
    ax.set_title("the covariance window: diffeo control collapses,\nphysical signal persists")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(FIGS / "ch25_fig2_covariance_window.png", dpi=150)
    plt.close(fig)

    # ---- Figure 25.3: convergence of kappa/C onto the DeWitt pattern ----
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    style = {"TT+": ("o-", "seagreen"), "TTx": ("D-", "olive"),
             "conf": ("s-", "crimson"), "diffeo": ("^-", "steelblue")}
    for m in sectors:
        vals = [results["table"][str(s)][m]["kappa_over_C"] for s in S_LIST]
        ax.plot(S_LIST, vals, style[m][0], color=style[m][1], label=m)
    ax.axhline(-1, color="k", ls="--", lw=0.9)
    ax.axhline(+1, color="k", ls=":", lw=0.9)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_yscale("symlog", linthresh=2)
    ax.set_xlabel("proper-time regulator $s$ (lattice units)")
    ax.set_ylabel(r"$\kappa / C(s)$")
    ax.set_title("convergence to the universal coefficient:\n"
                 r"TT $\to -1$, conformal opposite sign, diffeo $\to$ 0")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(FIGS / "ch25_fig3_kappa_convergence.png", dpi=150)
    plt.close(fig)
    print("# figures -> ch25_fig2_covariance_window.png, ch25_fig3_kappa_convergence.png")


if __name__ == "__main__":
    main()
