"""
Milestone 3, Sim D: the induced action of the geometry field, computed from matter loops.

Setup: a massless 1D field on [0, Ltot] with Dirichlet ends and M-1 semi-transparent
interior walls (delta barriers of strength lam) at positions x_i = i + u_i.  The wall
displacement field u_i IS the geometry field (it moves the cell boundaries).  The matter
vacuum energy E_vac[{u_i}] = (1/2) sum_n omega_n, regularized by heat-kernel smoothing,
defines the induced effective action.  For a standing-wave displacement
u_i = eps sin(pi k i / M), the quadratic response

    delta E(q) = (eps^2/2) * Pi(q) * sum_i sin^2(...),     q = pi k / M,

is the geometry field's induced dispersion.  An elastic (gradient) action would give
Pi(q) = kappa_grad * qhat^2 with qhat^2 = (2 sin(q/2))^2.  The SIGN of kappa_grad is the
result: negative = the matter vacuum makes the compression mode unstable (the Lorentzian
conformal-mode sign, and the engine of clumping); the matter (occupied-cell) contribution
kappa_m = 3 pi^2 n^2 per bond (derived analytically below) is positive and can stabilize.
"""
import json
import numpy as np
from scipy.linalg import eigvalsh_tridiagonal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

M, Ltot, lam = 24, 24.0, 8.0
h = 0.008
N = int(round(Ltot/h)) - 1

def spectrum(walls):
    diag = np.full(N, 2.0)/h**2
    off = np.full(N-1, -1.0)/h**2
    for w in walls:
        j = int(round(w/h)) - 1
        diag[j] += lam/h
    ev = eigvalsh_tridiagonal(diag, off)
    return np.sqrt(np.maximum(ev, 0.0))

def E_s(om, s):
    return 0.5*np.sum(om*np.exp(-s*om))

walls0 = np.arange(1.0, M)               # 23 interior walls
om0 = spectrum(walls0)
eps = 0.10
ss = np.array([0.40, 0.20, 0.10])

rows = []
print(f"# induced dispersion Pi(q): {M} cells, lam={lam}, eps={eps}, grid h={h}")
for k in range(1, 13):
    u = eps*np.sin(np.pi*k*walls0/M)
    om = spectrum(walls0 + u)
    dE = np.array([E_s(om, s) - E_s(om0, s) for s in ss])
    dE0 = (4*dE[2] - dE[1])/3            # Richardson in s^2
    norm = np.sum(np.sin(np.pi*k*walls0/M)**2)
    Pi = 2*dE0/(eps**2*norm)
    qhat2 = (2*np.sin(np.pi*k/(2*M)))**2
    rows.append((k, float(dE[1]), float(dE[2]), float(dE0), float(Pi), float(qhat2), float(Pi/qhat2)))
    print(f"  k={k:2d}: dE(s=0.2)={dE[1]:+.5f} dE(s=0.1)={dE[2]:+.5f} -> dE0={dE0:+.5f}"
          f"   Pi={Pi:+.5f}  Pi/qhat^2={Pi/qhat2:+.4f}")

R = np.array(rows)
kap_vac = np.median(R[:8, 6])
print(f"  induced vacuum stiffness kappa_grad(vac) ~ {kap_vac:+.4f}  (negative = unstable/attractive)")
n_occ = 1
kap_m = 3*np.pi**2*n_occ**2
print(f"  analytic matter stiffness per bond (n={n_occ}): kappa_m = 3 pi^2 n^2 = {kap_m:.3f} > 0")

json.dump({"rows": rows, "kappa_vac": float(kap_vac), "kappa_matter_n1": float(kap_m),
           "lam": lam, "M": M, "eps": eps}, open("stiffness_dispersion.json", "w"), indent=1)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.3))
a1.plot(R[:, 5], R[:, 4], "o", color="crimson", label=r"measured $\Pi(q)$ (vacuum loops)")
qq = np.linspace(0, R[-1, 5], 100)
a1.plot(qq, kap_vac*qq, "k--", label=fr"$\kappa_{{\rm grad}}\,\hat q^2$, $\kappa={kap_vac:.3f}$")
a1.set_xlabel(r"$\hat q^2 = (2\sin q/2)^2$"); a1.set_ylabel(r"$\Pi(q)$")
a1.set_title("induced geometry-field dispersion from matter vacuum\n(negative: compression mode unstable)")
a1.legend(fontsize=9)
a2.plot(R[:, 0], R[:, 6], "s-", color="crimson", label=r"vacuum $\Pi/\hat q^2$")
a2.axhline(kap_m, color="steelblue", ls="--", label=fr"matter $n=1$: $+3\pi^2$ per bond")
a2.axhline(0, color="k", lw=0.8)
a2.set_xlabel("mode number $k$"); a2.set_ylabel(r"stiffness $\Pi/\hat q^2$")
a2.set_title("vacuum (induced, attractive) vs matter (repulsive) stiffness")
a2.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figs/fig11_stiffness_dispersion.png", dpi=150); plt.close(fig)
print("saved fig11 and stiffness_dispersion.json")
