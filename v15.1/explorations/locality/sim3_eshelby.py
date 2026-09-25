"""
Sim 3 - embedded (conforming) cells: static response of the neighbours to a step in one cell.
Linear elasticity of the quantum-cell medium (moduli from cellmedium.py) on a periodic lattice,
solved exactly in Fourier space.  Source: change of first Piola stress in the central cell when one
fermion is promoted (n1,n2)->(n1+1,n2) (directional) or (n,n)->(n+1,n+1) (isotropic).
"""
import numpy as np
C = np.load("figs/C_tensor.npy")
W0 = 17854.11436

def solve(dP, Ccell=C, Ncell=128, p=4, sig=0.12):
    N = Ncell*p; h = 1.0/p
    k = 2*np.pi*np.fft.fftfreq(N, d=h)
    K1, K2 = np.meshgrid(k, k, indexing="ij")           # axis0 = x, axis1 = y
    kv = np.stack([K1, K2])
    chi = np.zeros((N, N)); c0 = N//2 - p//2; chi[c0:c0+p, c0:c0+p] = 1.0
    chih = np.fft.fft2(chi)*np.exp(-0.5*sig**2*(K1**2+K2**2))
    G = np.einsum("ijkl,jab,lab->ikab", Ccell, kv, kv)   # acoustic tensor per k
    rhs = 1j*np.einsum("ij,jab->iab", dP, kv)*chih       # k_j C k_l u_k = i k_j dP_ij chi
    G = np.moveaxis(G, (0, 1), (2, 3)); rhs = np.moveaxis(rhs, 0, 2)
    G[0, 0] = np.eye(2); rhs[0, 0] = 0
    uh = np.linalg.solve(G, rhs[..., None])[..., 0]
    uh = np.moveaxis(uh, 2, 0)
    grad = np.real(np.fft.ifft2(1j*kv[None, :, :, :]*uh[:, None, :, :], axes=(-2, -1)))  # grad[i,j] = d_j u_i
    u = np.real(np.fft.ifft2(uh, axes=(-2, -1)))
    x = (np.arange(N) - (c0 + p/2 - 0.5))*h
    return x, u, grad, chi

def analyse(label, dP, n_target):
    x, u, grad, chi = solve(dP)
    eps = 0.5*(grad + grad.transpose(1, 0, 2, 3))
    th = eps[0, 0] + eps[1, 1]; dev = np.sqrt(((eps[0, 0]-eps[1, 1])/2)**2 + eps[0, 1]**2)
    X, Y = np.meshgrid(x, x, indexing="ij"); R = np.hypot(X, Y)
    inc = chi > 0
    e11c = eps[0, 0][inc].mean(); e22c = eps[1, 1][inc].mean()
    rb = np.geomspace(2.5, 30, 12)
    prof = np.array([np.sqrt((eps**2).sum(axis=(0, 1)))[(R > r*0.93) & (R < r*1.07)].mean() for r in rb])
    slope = np.polyfit(np.log(rb), np.log(prof), 1)[0]
    ring = (R > 4) & (R < 20)
    frac_vol = np.abs(th[ring]).mean()/np.abs(dev[ring]).mean()
    print(f"[{label}] centre cell strain e11={e11c:.3e} e22={e22c:.3e}; iso-energy target 1/n = {1/n_target:.3e}"
          f" -> realised fraction {e11c*n_target:.3e}")
    print(f"   |strain| decay exponent (r=2.5..30 cells) = {slope:.3f} (2D point dipole: -2)")
    print(f"   outside (4<r<20): <|dilation|>/<|shear|> = {frac_vol:.3f}")
    return dict(x=x, u=u, eps=eps, th=th, dev=dev, rb=rb, prof=prof, slope=slope, e11c=e11c, frac=frac_vol)

if __name__ == "__main__":
    n1 = 7
    dP_dir = np.array([[-np.pi**2*(2*n1+1), 0], [0, 0]])
    dP_iso = np.array([[-np.pi**2*(2*n1+1), 0], [0, -np.pi**2*(2*n1+1)]])
    print("promotion energy at fixed shape:", 0.5*np.pi**2*(2*n1+1), " cell energy W0:", W0)
    Rd = analyse("directional", dP_dir, n1); Ri = analyse("isotropic", dP_iso, n1)
    # isotropic control medium: same bulk, isotropic shear -> dilation outside should vanish
    lam, mu = 0.0, 866.6
    Ciso = np.zeros((2, 2, 2, 2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    Ciso[i, j, k, l] = lam*(i == j)*(k == l) + mu*((i == k)*(j == l) + (i == l)*(j == k))
    x, u, grad, chi = solve(dP_iso, Ciso)
    eps = 0.5*(grad + grad.transpose(1, 0, 2, 3)); th = eps[0, 0]+eps[1, 1]
    dev = np.sqrt(((eps[0, 0]-eps[1, 1])/2)**2 + eps[0, 1]**2)
    X, Y = np.meshgrid(x, x, indexing="ij"); R = np.hypot(X, Y); ring = (R > 4) & (R < 20)
    print(f"[control: isotropic medium, isotropic step] outside <|dilation|>/<|shear|> = {np.abs(th[ring]).mean()/np.abs(dev[ring]).mean():.2e}")
    np.savez("figs/sim3.npz", **{f"d_{k}": v for k, v in Rd.items() if k != "u"}, **{f"i_{k}": v for k, v in Ri.items() if k != "u"})
