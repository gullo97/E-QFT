"""
Milestone 3, Sim E: the static weak-field limit on the 3D cell lattice.

Linearizing the constraint-closed cell dynamics about a uniform background, with the
gradient stiffness of Sim D coupling neighbouring cells, the static sector reduces to
the lattice Poisson equation for the local log-scale-factor perturbation:

    kappa * (lattice Laplacian) delta_alpha(x) = - c * delta_rho(x),

i.e. delta_alpha = -Phi_Newton (units kappa = c = 1 here; the structure is the result).
We solve it on a 64^3 periodic cell lattice by FFT and test:
  (1) point mass  ->  delta_alpha(r) ~ 1/r  (the Newtonian potential profile);
  (2) two masses  ->  interaction energy E_int(d) ~ -1/d, force ~ 1/d^2 (inverse square).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
Ng = 128
k = 2*np.pi*np.fft.fftfreq(Ng)
KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
khat2 = (2*np.sin(KX/2))**2 + (2*np.sin(KY/2))**2 + (2*np.sin(KZ/2))**2
khat2[0, 0, 0] = 1.0

def poisson(rho):
    rho = rho - rho.mean()                  # periodic solvability
    F = np.fft.fftn(rho)/(-khat2)
    F[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(F))         # solves  lap phi = rho  ->  alpha = -phi

# (1) point mass
rho = np.zeros((Ng,)*3); rho[0, 0, 0] = 1.0
alpha = -poisson(rho)                        # delta_alpha = -Phi > 0 near the mass
ii0 = np.minimum(np.arange(Ng), Ng - np.arange(Ng))
RX0, RY0, RZ0 = np.meshgrid(ii0, ii0, ii0, indexing="ij")
rr = np.sqrt(RX0**2 + RY0**2 + RZ0**2)
alpha = alpha - alpha[rr > 0.46*Ng].mean()   # remove periodic far-field offset
ii = np.minimum(np.arange(Ng), Ng - np.arange(Ng))
RX, RY, RZ = np.meshgrid(ii, ii, ii, indexing="ij")
r = np.sqrt(RX**2 + RY**2 + RZ**2).ravel()
av = alpha.ravel()
rb = np.linspace(3, 16, 14)
prof = np.array([av[(r > x-0.5) & (r <= x+0.5)].mean() for x in rb])
X = np.vstack([1/rb, np.ones_like(rb)]).T          # alpha = A/r + c  (c = zero-mode constant)
(Acoef, coff), *_ = np.linalg.lstsq(X, prof, rcond=None)
prof_c = prof - coff
slope = np.polyfit(np.log(rb), np.log(prof_c), 1)[0]
resid = np.abs(prof - (Acoef/rb + coff)).max()/prof.max()
print(f"# point mass: A/r+c fit: A={Acoef:.4f}, zero-mode c={coff:+.5f}, max residual {resid:.1e}")
print(f"# point mass: exponent after zero-mode removal = {slope:.3f}   (Newton: -1)")

# (2) two masses: interaction energy and force
ds = np.arange(4, 19)
Eint = []
for d in ds:
    Eint.append(-alpha[d % Ng, 0, 0])        # E_int = -int Phi_1 rho_2 = -(+alpha at d)... sign: bound state
Eint = np.array(Eint)
XE = np.vstack([1/ds, np.ones_like(ds, float)]).T
(AE, cE), *_ = np.linalg.lstsq(XE, Eint, rcond=None)
Eint = Eint - cE
sE = np.polyfit(np.log(ds), np.log(-Eint), 1)[0]
F = -np.gradient(Eint, ds.astype(float))
sF = np.polyfit(np.log(ds[2:-2]), np.log(np.abs(F[2:-2])), 1)[0]
print(f"# two masses: slope of |E_int(d)| = {sE:.3f} (Newton: -1);  slope of |F(d)| = {sF:.3f} (Newton: -2)")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.3))
a1.loglog(rb, prof_c, "o", color="crimson", label=r"lattice $\delta\alpha(r)$")
a1.loglog(rb, prof_c[4]*rb[4]/rb, "k--", label=r"$\propto 1/r$  (Newtonian potential)")
a1.set_xlabel("cell distance $r$"); a1.set_ylabel(r"$\delta\alpha = -\Phi$")
a1.set_title(f"point mass on the cell lattice: slope {slope:.2f}")
a1.legend(fontsize=9)
a2.loglog(ds, -Eint, "s", color="crimson", label=r"$|E_{\rm int}(d)|$, slope %.2f" % sE)
a2.loglog(ds[2:-2], np.abs(F[2:-2]), "^", color="steelblue", label=r"$|F(d)|$, slope %.2f" % sF)
a2.loglog(ds, (-Eint[2])*ds[2]/ds, "k--", lw=1)
a2.loglog(ds, np.abs(F[4])*(ds[6]/ds)**2, "k:", lw=1)
a2.set_xlabel("separation $d$ (cells)"); a2.set_ylabel("magnitude")
a2.set_title("two masses: $1/d$ binding energy, $1/d^2$ force")
a2.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figs/fig12_newton.png", dpi=150); plt.close(fig)
print("saved fig12")
