#!/usr/bin/env python3
"""
ch18_metric_updates.py — kinematic checks for Ch. 18 (prints every quoted number).

(1) generator algebra: the 9 edge-space generators (3 directional dilations E_aa,
    3 symmetric shears E_ab + E_ba, 3 rotations E_ab - E_ba) span gl(3), and all
    commutators lie in that span (closure residual at machine zero).
(2) Gram update law (18.2): an iso-energy directional transition e_a -> ((n+1)/n) e_a
    acts on g = E^T E as a congruence; over random transition sequences the update
    matches the recomputed Gram matrix at machine zero; the update is rank-1 at
    rectangular configurations and rank-2 off them.
(3) transitivity: random SPD metrics are realized exactly as Gram matrices of an
    edge triple (Cholesky) — the symmetric sector acts transitively on metrics.
(4) NEGATIVE control (the lesson box of Sec. 18.2): the covariant shortcut
    E = pi^2 n_a n_b g^{ab} read at non-diagonal g^{ab} is falsified: the exact
    Dirichlet spectrum of a sheared cell is EVEN in the shear (reflection x -> 1-x),
    while the shortcut is linear in it; at g^{12} = 0.30 the shortcut misses the
    ground state by ~29% (n=(1,1)) / ~26% (n=(1,-1)), where the diagonal formula at
    the same stretch entries misses only the genuine O(shear^2) shift (~1.5%).
(5) POSITIVE control: at zero shear the FD spectrum matches pi^2(m^2 A11 + n^2 A22).
"""
import numpy as np
from scipy.sparse import identity, kron, diags
from scipy.sparse.linalg import eigsh

rng = np.random.default_rng(7)

# ------------------------- (1) gl(3) span and closure -------------------------
def basis9():
    Ms = []
    for a in range(3):                        # directional dilations
        M = np.zeros((3, 3)); M[a, a] = 1.0; Ms.append(M)
    for (a, b) in ((0, 1), (0, 2), (1, 2)):   # symmetric shears
        M = np.zeros((3, 3)); M[a, b] = M[b, a] = 1.0; Ms.append(M)
    for (a, b) in ((0, 1), (0, 2), (1, 2)):   # frame rotations
        M = np.zeros((3, 3)); M[a, b] = 1.0; M[b, a] = -1.0; Ms.append(M)
    return Ms

B = basis9()
V = np.array([M.flatten() for M in B]).T            # 9 x 9
rank = np.linalg.matrix_rank(V)
res = 0.0
for X in B:
    for Y in B:
        C = (X @ Y - Y @ X).flatten()
        coef, *_ = np.linalg.lstsq(V, C, rcond=None)
        res = max(res, float(np.max(np.abs(V @ coef - C))))
print(f"(1) gl(3): span rank = {rank}/9; max commutator closure residual = {res:.1e}")

# ------------------------- (2) Gram congruence updates ------------------------
def gram(E):                                        # edges as columns of E
    return E.T @ E

def rank_eps(M, eps=1e-10):
    return int(np.sum(np.linalg.svd(M, compute_uv=False) > eps))

Erun = np.diag([1.0, 1.3, 0.8])                     # rectangular start
nrun = np.array([7.0, 11.0, 5.0])
worst = 0.0
for _ in range(200):
    a = int(rng.integers(0, 3))
    lam = (nrun[a] + 1) / nrun[a]
    D = np.eye(3); D[a, a] = lam
    g_pred = D @ gram(Erun) @ D                     # the congruence update (18.2)
    Erun = Erun @ D                                 # dilate edge a
    nrun[a] += 1
    worst = max(worst, float(np.max(np.abs(gram(Erun) - g_pred))))
print(f"(2) congruence update (18.2) over 200 random transitions: max residual = {worst:.1e}")
Erect = np.diag([1.0, 1.3, 0.8])
Eshear = Erect.copy(); Eshear[0, 1] = 0.45          # edge 2 tilted toward edge 1
for tag, E0 in (("rectangular", Erect), ("sheared", Eshear)):
    D = np.eye(3); D[0, 0] = 1 + 1.0/7.0
    dg = D @ gram(E0) @ D - gram(E0)
    print(f"    rank of the update at a {tag} configuration: {rank_eps(dg)}")

# ------------------------------ (3) transitivity ------------------------------
w = 0.0
for _ in range(50):
    A = rng.standard_normal((3, 3))
    g_t = A @ A.T + 3*np.eye(3)                     # random SPD metric
    Echol = np.linalg.cholesky(g_t).T               # edge triple with g = E^T E
    w = max(w, float(np.max(np.abs(gram(Echol) - g_t))))
print(f"(3) transitivity: 50 random SPD metrics realized as Gram matrices, residual = {w:.1e}")

# ------------- (4)/(5) sheared-cell spectrum vs the covariant shortcut --------
NG = 120
h = 1.0/(NG + 1)
D2 = diags([1, -2, 1], [-1, 0, 1], shape=(NG, NG))/h**2
D1 = diags([-1, 1], [-1, 1], shape=(NG, NG))/(2*h)
I = identity(NG)
DXX = kron(D2, I); DYY = kron(I, D2); DXY = kron(D1, D1)

def spec(A11, A12, A22, k=4):
    H = (-(A11*DXX + A22*DYY + 2*A12*DXY)).tocsc()
    E, _ = eigsh(H, k=k, sigma=0, which="LM")
    return np.sort(E)

A11, A22 = 1.0, 1.21
E0 = spec(A11, 0.0, A22)
Erect0 = np.pi**2*(A11 + A22)
print(f"(5) rectangular control: FD ground state {E0[0]:.4f} vs pi^2(A11+A22) = {Erect0:.4f}"
      f"   (rel {abs(E0[0]-Erect0)/E0[0]:.1e}, discretization)")
s = 0.30
Ep, Em = spec(A11, +s, A22), spec(A11, -s, A22)
print(f"(4) sheared cell, g^12 = {s}: even-in-shear check max|E(+s)-E(-s)| = {np.max(np.abs(Ep-Em)):.1e}")
Eex = Ep[0]
for n1, n2 in ((1, 1), (1, -1)):
    Ef = np.pi**2*(n1*n1*A11 + 2*n1*n2*s + n2*n2*A22)
    print(f"    covariant shortcut, n=({n1:+d},{n2:+d}): {Ef:.4f} vs exact {Eex:.4f}"
          f"   (rel err {abs(Ef-Eex)/Eex:.1%})")
print(f"    diagonal formula at the same stretch entries: {Erect0:.4f}"
      f"   (rel err {abs(Erect0-Eex)/Eex:.1%} = the genuine O(shear^2) shift)")
print("    verdict: the shortcut's linear-in-shear cross term is spurious (Ch. 19, Splitting (ii));")
print("    (18.3) is exact on the rectangular sector and first-order-stationary off it.")
