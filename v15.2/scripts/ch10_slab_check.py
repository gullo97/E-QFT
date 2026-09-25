#!/usr/bin/env python3
"""
ch10_slab_check.py — Slab Eigenvalue Theorem at nonzero transverse momenta
(Ch. 10, App. 10.A).

Four-component Dirac transfer-matrix solver for the bulk-CP-mass MIT slab:
confined direction x in [0, L], transverse plane waves e^{i(k_y y + k_z z)}.
Checks, at randomized (k_y, k_z):
  * every level sits at E = sqrt(p_j^2 + k_perp^2 + m^2) with p_j a root of
    tan(pL) = -p/(m cos delta) -- the confined quantization is (10.8),
    independent of the transverse momenta;
  * each level is exactly twofold (spin) degenerate (two vanishing singular
    values of the boundary-matching matrix);
  * the negative branch passes the same matching;
  * no midgap solutions, including the special point E = -m cos(delta) where
    the 1+1D derivation of (10.8) divides by zero.
"""
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar

I2 = np.eye(2)
Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
g1 = np.block([[Z2, sx], [-sx, Z2]])
g2 = np.block([[Z2, sy], [-sy, Z2]])
g3 = np.block([[Z2, sz], [-sz, Z2]])
g5 = np.block([[Z2, I2], [I2, Z2]]).astype(complex)


def M(E, ky, kz, mR, mI):
    """psi' = M psi for psi(x) e^{i(k_y y + k_z z) - iEt}."""
    return 1j * g1 @ (-E * g0 + ky * g2 + kz * g3 + mR * np.eye(4) + 1j * mI * g5)


def _plus_subspace(op):
    w, V = np.linalg.eig(op)
    cols = np.array([V[:, i] for i in range(4) if abs(w[i] - 1) < 1e-10]).T
    q, _ = np.linalg.qr(cols)
    return q


V0 = _plus_subspace(+1j * g1)     # x = 0 wall: +i gamma^1 psi = psi
VL = _plus_subspace(-1j * g1)     # x = L wall: -i gamma^1 psi = psi
P_perp = np.eye(4) - VL @ VL.conj().T


def svals(E, ky, kz, mR, mI, L):
    """Singular values of the boundary-matching matrix (zero <=> eigenvalue)."""
    A = P_perp @ expm(M(E, ky, kz, mR, mI) * L) @ V0
    return np.linalg.svd(A, compute_uv=False)


def find_levels(ky, kz, mR, mI, L, Emax, nE=3000):
    m2 = mR ** 2 + mI ** 2 + ky ** 2 + kz ** 2
    Es = np.linspace(np.sqrt(m2) + 1e-4, Emax, nE)
    r = np.array([svals(E, ky, kz, mR, mI, L)[-1] for E in Es])
    roots = []
    for i in range(1, nE - 1):
        if r[i] < r[i - 1] and r[i] < r[i + 1] and r[i] < 0.2:
            res = minimize_scalar(lambda E: svals(E, ky, kz, mR, mI, L)[-1],
                                  bounds=(Es[i - 1], Es[i + 1]), method="bounded",
                                  options={"xatol": 1e-13})
            if svals(res.x, ky, kz, mR, mI, L)[-1] < 1e-6:
                roots.append(res.x)
    return np.array(sorted(roots))


def tan_roots(mcosd, L, nmax=8):
    ps = []
    for j in range(1, nmax + 1):
        lo, hi = (j - 0.4999) * np.pi / L + 1e-9, (j + 0.4999) * np.pi / L
        try:
            ps.append(brentq(lambda p: np.tan(p * L) + p / mcosd, lo, hi,
                             xtol=1e-14))
        except ValueError:
            pass
    return np.array(ps)


def main():
    m, L, d = 1.0, 1.0, 0.7
    mR, mI = m * np.cos(d), m * np.sin(d)
    ps = tan_roots(mR, L)
    print(f"=== Slab Eigenvalue Theorem check: m={m}, L={L}, delta={d} ===")
    print(f"confined roots p_j of tan(pL) = -p/(m cos d): {np.round(ps, 6)}")
    rng = np.random.default_rng(7)
    worst, worst_deg = 0.0, 0.0
    for _ in range(4):
        ky, kz = rng.uniform(-3, 3, 2)
        Epred = np.sort(np.sqrt(ps ** 2 + ky ** 2 + kz ** 2 + m ** 2))
        lv = find_levels(ky, kz, mR, mI, L, Emax=Epred.max() + 0.5)
        n = min(len(lv), len(Epred))
        err = np.max(np.abs(lv[:n] - Epred[:n]))
        sv = svals(lv[0], ky, kz, mR, mI, L)
        worst = max(worst, err)
        worst_deg = max(worst_deg, sv[-2])
        print(f"  ky={ky:+.3f} kz={kz:+.3f}: {len(lv)} levels, "
              f"max|E - E_pred| = {err:.1e}, "
              f"two smallest svals at lowest level: {sv[-1]:.1e}, {sv[-2]:.1e}")
    print(f"worst |E - E_pred| over the four momenta: {worst:.1e}")
    print(f"worst second singular value (twofold spin degeneracy): {worst_deg:.1e}")

    ky, kz = 1.3, -0.8
    Eneg = np.sqrt(ps[0] ** 2 + ky ** 2 + kz ** 2 + m ** 2)
    print(f"negative branch residual at E = -{Eneg:.6f}: "
          f"{svals(-Eneg, ky, kz, mR, mI, L)[-1]:.1e} (should be ~0)")

    Es = np.linspace(-0.999 * m, 0.999 * m, 400)
    rmin = min(svals(E, 0.0, 0.0, mR, mI, L)[-1] for E in Es)
    print(f"midgap scan (|E| < m, k_perp = 0): min residual = {rmin:.2f} "
          "(O(1): no midgap level)")
    print(f"residual at the special point E = -m cos d = {-mR:.4f}: "
          f"{svals(-mR, 0.0, 0.0, mR, mI, L)[-1]:.2f} (nonzero: no eigenvalue)")


if __name__ == "__main__":
    main()
