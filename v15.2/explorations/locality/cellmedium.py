"""
Elastic medium made of quantum cells (Sim 3).  A cell is the unit square mapped by a
constant deformation gradient F (parallelogram); its metric g = F^T F, A = g^{-1}.
Occupants: a closed-shell Fermi sea of the Dirichlet problem H = -1/2 A^{ab} d_a d_b.
W(F) = sum of the lowest N_f eigenvalues (exact Galerkin in the product sine basis).
"""
import numpy as np
from itertools import product

def D1(K):
    m = np.arange(1, K+1)[:, None]; n = np.arange(1, K+1)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        D = 2*m*n*(1-(-1.0)**(m+n))/(m**2-n**2)
    D[m.repeat(K, 1) == n.repeat(K, 0)] = 0.0
    return D                                   # <m| d/dx |n> on [0,1]

class Cell:
    def __init__(self, R2=100, K=34):
        self.K = K
        idx = np.array(list(product(range(1, K+1), repeat=2)))
        self.m1, self.m2 = idx[:, 0], idx[:, 1]
        D = D1(K); self.V = np.kron(D, D)      # <m1 m2| dx dy |n1 n2>
        occ = self.m1**2 + self.m2**2 <= R2
        self.Nf = int(occ.sum())
        e = np.sort(self.m1**2 + self.m2**2)
        self.gap = e[self.Nf] - e[self.Nf-1]
    def H(self, A):
        d = 0.5*np.pi**2*(A[0, 0]*self.m1**2 + A[1, 1]*self.m2**2)
        return np.diag(d) - A[0, 1]*self.V
    def W(self, F):
        A = np.linalg.inv(F.T @ F)
        ev = np.linalg.eigvalsh(self.H(A))
        return ev[:self.Nf].sum()

def hessian_F(cell, h=2e-3):
    I = np.eye(2); f0 = cell.W(I)
    grad = np.zeros((2, 2)); C = np.zeros((2, 2, 2, 2))
    E = lambda i, j: np.outer(np.eye(2)[i], np.eye(2)[j])
    for i, j in product(range(2), repeat=2):
        grad[i, j] = (cell.W(I+h*E(i, j)) - cell.W(I-h*E(i, j)))/(2*h)
    for (i, j), (k, l) in product(product(range(2), repeat=2), repeat=2):
        C[i, j, k, l] = (cell.W(I+h*E(i, j)+h*E(k, l)) - cell.W(I+h*E(i, j)-h*E(k, l))
                         - cell.W(I-h*E(i, j)+h*E(k, l)) + cell.W(I-h*E(i, j)-h*E(k, l)))/(4*h*h)
    return f0, grad, C

if __name__ == "__main__":
    for K in (24, 34):
        c = Cell(K=K); f0, g, C = hessian_F(c)
        print(f"K={K}: N_f={c.Nf}, shell gap={c.gap}, W={f0:.4f}, first Piola P=\n{np.round(g,4)}")
        print("  C1111 C2222 C1122 C1212 C1221 C2121:", np.round([C[0,0,0,0], C[1,1,1,1], C[0,0,1,1], C[0,1,0,1], C[0,1,1,0], C[1,0,1,0]], 3))
    # rotation invariance check
    th = 0.3; Rm = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    F = np.array([[1.05, 0.02], [-0.01, 0.97]])
    print("rotation invariance W(RF)-W(F):", c.W(Rm @ F) - c.W(F))
    # Polya check: area-preserving shear raises energy
    for e in (0.02, 0.05):
        F = np.array([[1, e], [0, 1]]); print(f"area-preserving shear {e}: dW = {c.W(F)-c.W(np.eye(2)):+.5f}")
    np.save("figs/C_tensor.npy", C); np.save("figs/P0.npy", g)
