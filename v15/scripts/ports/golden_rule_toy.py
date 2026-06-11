"""
F3: paper 1's projection-norm probability rule DERIVED from unitary dynamics.

Enlarged Hilbert space H = H_box(L) (+) H_box(L'); coupling V = g * embedding
(truncation for contraction, zero-extension for expansion).

Paper-1 convention: n1 = ACTIVE particle mode (iso-energy step => L' = (n1-1)L/n1
for one contraction step); n2 = SPECTATOR mode.  First-order PT for the spectator:

  P_transfer(t)/(g t)^2  ->  sum_f |<f|V|n2>|^2 / g^2  =  || psi_{n2} restricted ||^2
                          =  (n1-1)/n1 + sin(2 pi n2/n1)/(2 pi n2)  =  Ncon(n1,n2)

(Parseval over the complete new-box basis), and = 1 for expansion (zero-extension
is an isometry).  Factorisation over N spectators: V acts as a tensor product.
"""
import numpy as np

def box_modes(L, nmax, x):
    return np.array([np.sqrt(2/L)*np.sin(k*np.pi*x/L)*(x <= L) for k in range(1, nmax+1)])

def weights(spec_mode, Lp, L=1.0, nmax=220, ngrid=24001, g=1e-3, ts=(0.02, 0.008)):
    x = np.linspace(0, max(L, Lp), ngrid)
    A = box_modes(L,  nmax, x); B = box_modes(Lp, nmax, x)
    E1 = (np.arange(1, nmax+1)*np.pi/L)**2
    E2 = (np.arange(1, nmax+1)*np.pi/Lp)**2
    wts = np.gradient(x)
    Vmat = g * ((B * wts) @ A.T)
    W_parseval = np.sum(np.abs(Vmat[:, spec_mode-1])**2) / g**2     # golden-rule sum
    H = np.block([[np.diag(E1), Vmat.T], [Vmat, np.diag(E2)]])
    w, U = np.linalg.eigh(H)
    psi0 = np.zeros(2*nmax); psi0[spec_mode-1] = 1.0
    Ps = []
    for t in ts:
        psit = U @ (np.exp(-1j*w*t) * (U.T @ psi0))
        Ps.append(np.sum(np.abs(psit[nmax:])**2) / (g*t)**2)
    return W_parseval, Ps

def Ncon(n1, n2):
    return (n1-1)/n1 + np.sin(2*np.pi*n2/n1)/(2*np.pi*n2)

print("# contraction L -> (n1-1)L/n1, spectator mode n2; golden-rule weight vs Ncon(n1,n2)")
for (n1, n2) in ((5, 3), (5, 4), (8, 5), (12, 7)):
    W, Ps = weights(n2, (n1-1)/n1)
    print(f"  n1={n1:2d}, n2={n2:2d}:  sum|<f|V|i>|^2 = {W:.5f}   "
          f"P/(gt)^2 (t->0): {Ps[0]:.5f}, {Ps[1]:.5f}   Ncon = {Ncon(n1,n2):.5f}")
print("# expansion L -> 2L (isometry): weight must be 1, any spectator")
for n2 in (3, 7):
    W, Ps = weights(n2, 2.0)
    print(f"  n2={n2:2d}:  sum|<f|V|i>|^2 = {W:.5f}   P/(gt)^2: {Ps[0]:.5f}, {Ps[1]:.5f}")
