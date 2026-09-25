"""
Sim 6 - the 'vote' of many occupants on a LOCAL metric change: dE = -1/2 \int T^{ab} h_ab.
Check against exact finite-difference Laplace-Beltrami eigenvalues on the unit square with
g = diag(1 + h(x,y), 1), h = eps * Gaussian bump at X0 (width w).
"""
import numpy as np
from scipy.sparse import diags, kron, identity
from scipy.sparse.linalg import eigsh

Nx = 90; hgrid = 1.0/(Nx+1); xs = np.arange(1, Nx+1)*hgrid
def bump(X0, Y0, w):
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    return np.exp(-((X-X0)**2+(Y-Y0)**2)/(2*w*w))

def spectrum(hfun, k=60):
    # flux form: -(1/sqrtg)[d_x(sqrtg g^xx d_x) + d_y(sqrtg d_y)], g=diag(1+h,1): sqrtg=sqrt(1+h), sqrtg g^xx = 1/sqrt(1+h)
    X = np.concatenate([[0], xs, [1]])
    Xf, Yf = np.meshgrid(X, X, indexing="ij")
    hf = hfun(Xf, Yf); sg = np.sqrt(1+hf); cx = 1/sg; cy = sg
    n = Nx; idx = lambda i, j: (i-1)*n + (j-1)
    rows, cols, vals = [], [], []
    for i in range(1, n+1):
        for j in range(1, n+1):
            p = idx(i, j); diag = 0.0
            for di, dj, cf in ((1, 0, 0.5*(cx[i, j]+cx[i+1, j])), (-1, 0, 0.5*(cx[i, j]+cx[i-1, j])),
                               (0, 1, 0.5*(cy[i, j]+cy[i, j+1])), (0, -1, 0.5*(cy[i, j]+cy[i, j-1]))):
                diag += cf
                ii, jj = i+di, j+dj
                if 1 <= ii <= n and 1 <= jj <= n:
                    rows.append(p); cols.append(idx(ii, jj)); vals.append(-cf)
            rows.append(p); cols.append(p); vals.append(diag)
    from scipy.sparse import coo_matrix
    A = coo_matrix((np.array(vals)/hgrid**2, (rows, cols)), shape=(n*n, n*n)).tocsc()
    Mv = sg[1:-1, 1:-1].ravel()
    Mh = diags(1/np.sqrt(Mv))                         # symmetric form M^-1/2 A M^-1/2
    ev = eigsh(Mh @ A @ Mh, k=k, sigma=0, which="LM", return_eigenvectors=False)
    return 0.5*np.sort(ev)                            # E = 1/2 * eigenvalue of -Laplace-Beltrami

def T_xx(modes):
    X, Y = np.meshgrid(xs, xs, indexing="ij"); T = np.zeros_like(X)
    for m, n in modes:
        psi = 2*np.sin(m*np.pi*X)*np.sin(n*np.pi*Y)
        px = 2*m*np.pi*np.cos(m*np.pi*X)*np.sin(n*np.pi*Y); py = 2*n*np.pi*np.sin(m*np.pi*X)*np.cos(n*np.pi*Y)
        E = 0.5*np.pi**2*(m*m+n*n)
        T += px**2 + E*psi**2 - 0.5*(px**2+py**2)
    return T

if __name__ == "__main__":
    E0 = spectrum(lambda X, Y: 0*X)
    # closed shell: modes with m^2+n^2 <= 50
    shell = [(m, n) for m in range(1, 9) for n in range(1, 9) if m*m+n*n <= 50]
    Nf = len(shell); print("Fermi sea N_f =", Nf, " FD gap at Fermi level:", E0[Nf]-E0[Nf-1])
    eps, w = 0.02, 0.07
    for (X0, Y0) in ((0.5, 0.5), (0.3, 0.62), (0.2, 0.2), (0.5, 0.25)):
        hp = lambda X, Y: eps*np.exp(-((X-X0)**2+(Y-Y0)**2)/(2*w*w))
        hm = lambda X, Y: -eps*np.exp(-((X-X0)**2+(Y-Y0)**2)/(2*w*w))
        Ep, Em = spectrum(hp), spectrum(hm)
        dE_exact = (Ep[:Nf].sum() - Em[:Nf].sum())/2
        dE_T = -0.5*np.sum(T_xx(shell)*eps*bump(X0, Y0, w))*hgrid**2
        dE1_exact = (Ep[0]-Em[0])/2; dE1_T = -0.5*np.sum(T_xx([(1, 1)])*eps*bump(X0, Y0, w))*hgrid**2
        print(f"bump at ({X0},{Y0}): Fermi sea dE exact {dE_exact:+.5e}, stress formula {dE_T:+.5e};"
              f"  ground state exact {dE1_exact:+.4e}, formula {dE1_T:+.4e}")
