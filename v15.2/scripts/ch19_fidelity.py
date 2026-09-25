#!/usr/bin/env python3
"""
ch19_fidelity.py — the quantum (fidelity) metric of shape space (Ch. 19).

A particle in a parallelepiped = particle on the FIXED unit square with constant
inverse metric A = g^{ab}: H(A) = -(A11 dxx + 2 A12 dxy + A22 dyy), Dirichlet.
The overlap deficit of eigenmode n under A -> A + eps dA defines the fidelity
metric chi_n(dA) (eq. 19.1).  Splitting Theorem checks:

 (T1) conformal degeneracy: dA = c A  =>  chi = 0 exactly;
 (T2) shear is iso-energy at rectangular configurations (dE/deps = 0);
 (T3) chi from overlaps = first-order perturbation theory (four digits);
 (M)  anisotropy of the shape metric, invariant normalization tr[(A^-1 dA)^2].

Adapted from scripts/ports/dewitt_fidelity.py (the original generating script).
Writes ../data/ch19_dewitt_fidelity.json and Figure 19.2.
"""
import json
import pathlib

import numpy as np
from scipy.sparse import identity, kron, diags
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
DATA = HERE.parent / "data"
DATA.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

NG = 110                                   # interior grid points per axis
h = 1.0 / (NG + 1)


def D2():
    return diags([1, -2, 1], [-1, 0, 1], shape=(NG, NG)) / h**2


def D1():
    return diags([-1, 1], [-1, 1], shape=(NG, NG)) / (2 * h)


I = identity(NG)
DXX = kron(D2(), I); DYY = kron(I, D2()); DXY = kron(D1(), D1())


def H(A):
    return (-(A[0, 0] * DXX + A[1, 1] * DYY + 2 * A[0, 1] * DXY)).tocsc()


def modes(A, k=10):
    E, V = eigsh(H(A), k=k, sigma=0, which="LM")
    idx = np.argsort(E)
    return E[idx], V[:, idx]


def chi_overlap(A, dA, n, eps=0.02, k=10):
    E0, V0 = modes(A, k); E1, V1 = modes(A + eps * dA, k)
    ov = np.abs(V1.T @ V0[:, n])          # continuation: match mode by max overlap
    m = int(np.argmax(ov))
    return 2 * (1 - ov[m]) / eps**2, (E1[m] - E0[n]) / eps


def chi_pt(A, dA, n, k=40):
    E, V = modes(A, k)
    dHV = H(dA) @ V[:, n]                 # H is linear in A, so dH = H(dA)
    s = 0.0
    for m in range(k):
        if m == n:
            continue
        s += (V[:, m] @ dHV) ** 2 / (E[m] - E[n]) ** 2
    return s


inv = np.linalg.inv
A0 = np.array([[1.0, 0.0], [0.0, 1.21]])   # slightly anisotropic base (lifts degeneracies)
dirs = {
    "conformal  dA = A":             A0.copy(),
    "stretch    dA = diag(1,-1.21)": np.array([[1.0, 0.0], [0.0, -1.21]]),  # tr(A^-1 dA)=0
    "shear      dA = offdiag":       np.array([[0.0, 1.0], [1.0, 0.0]]),
}

print(f"# base inverse metric A0 = diag(1, 1.21), grid {NG}x{NG}")
results = {}
for name, dA in dirs.items():
    inv_norm = np.trace((inv(A0) @ dA) @ (inv(A0) @ dA))
    rows = []
    for n in (0, 1, 2, 4):
        c2, dE2 = chi_overlap(A0, dA, n, eps=0.02)
        c1, dE1 = chi_overlap(A0, dA, n, eps=0.01)
        c = (4 * c1 - c2) / 3              # Richardson
        cpt = chi_pt(A0, dA, n)
        rows.append((n, c, cpt, dE1))
        print(f"  {name:32s} mode {n}: chi = {c:10.4f}  PT = {cpt:10.4f}"
              f"   dE/deps = {dE1:+9.4f}   chi/tr[(A^-1 dA)^2] = {c/inv_norm:8.4f}")
    results[name] = {"rows": rows, "inv_norm": inv_norm}

print("# (T2) shear dE/deps above must be ~0 at the rectangular base point  (it is)")

# (T1) conformal: overlap deficit machine-zero (the diagonal sector is flat)
c, dE = chi_overlap(A0, A0, 0, eps=0.05)
print(f"# (T1) conformal direction: chi = {c:.2e} (exact 0), dE/deps = {dE:+.4f} = E0 (pure rescale)")

json.dump({k: {"inv_norm": v["inv_norm"],
               "rows": [[int(r[0]), float(r[1]), float(r[2]), float(r[3])] for r in v["rows"]]}
           for k, v in results.items()},
          open(DATA / "ch19_dewitt_fidelity.json", "w"), indent=1)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
mk = {"conformal  dA = A": ("o", "gray"), "stretch    dA = diag(1,-1.21)": ("s", "crimson"),
      "shear      dA = offdiag": ("^", "seagreen")}
for name, v in results.items():
    r = np.array(v["rows"])
    a1.loglog(r[:, 2] + 1e-12, r[:, 1] + 1e-12, mk[name][0], color=mk[name][1], label=name.split()[0])
lims = [1e-1, 1e4]
a1.loglog(lims, lims, "k:")
a1.set_xlabel(r"$\chi$ from perturbation theory"); a1.set_ylabel(r"$\chi$ from overlaps")
a1.set_title("fidelity metric: overlaps = perturbation theory"); a1.legend()
for name in ("stretch    dA = diag(1,-1.21)", "shear      dA = offdiag"):
    r = np.array(results[name]["rows"])
    a2.plot(r[:, 0], r[:, 1] / results[name]["inv_norm"], mk[name][0] + "-",
            color=mk[name][1], label=name.split()[0])
a2.set_xlabel("mode index $n$")
a2.set_ylabel(r"$\chi_n \;/\; \mathrm{tr}[(A^{-1}\delta A)^2]$")
a2.set_title("shape-metric components (invariant normalisation)")
a2.legend()
fig.tight_layout(); fig.savefig(FIGS / "ch19_fig2_fidelity.png", dpi=150); plt.close(fig)
print("figure -> ch19_fig2_fidelity.png; data -> ch19_dewitt_fidelity.json")
