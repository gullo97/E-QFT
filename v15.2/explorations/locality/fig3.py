import numpy as np
from vstyle import *
import sim3_eshelby as S
n1 = 7; dP = np.array([[-np.pi**2*(2*n1+1), 0], [0, 0]])
x, u, grad, chi = S.solve(dP, Ncell=256, p=4)
eps = 0.5*(grad + grad.transpose(1, 0, 2, 3)); th = eps[0, 0]+eps[1, 1]
X, Y = np.meshgrid(x, x, indexing="ij"); R = np.hypot(X, Y)
xi, ui, gi, _ = S.solve(np.diag([dP[0, 0]]*2), Ncell=256, p=4)
epsi = 0.5*(gi + gi.transpose(1, 0, 2, 3))
p = 4; N = len(x)

fig = plt.figure(figsize=(11.4, 7.6))
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.36)
# (a) deformed lattice
ax = fig.add_subplot(gs[0, 0]); amp = 60
c = np.argmin(abs(x)); half = 6
idx = np.arange(-half, half+2)*p + int(np.argmin(abs(x+0.5)))
for a_ in idx:
    ax.plot(x[a_] + amp*u[0][a_, idx[0]:idx[-1]+1], x[idx[0]:idx[-1]+1] + amp*u[1][a_, idx[0]:idx[-1]+1], color=INK2, lw=0.8)
    ax.plot(x[idx[0]:idx[-1]+1] + amp*u[0][idx[0]:idx[-1]+1, a_], x[a_] + amp*u[1][idx[0]:idx[-1]+1, a_], color=INK2, lw=0.8)
ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, fc="#f6c7b3", ec="none", zorder=0))
ax.set_aspect("equal"); ax.set_xlim(-6, 6); ax.set_ylim(-6, 6); ax.grid(False)
ax.set_title("(a) cell lattice after the step\n" + r"(displacements $\times$" + f"{amp}; shaded = stepped cell)")
ax.set_xlabel("$x$ (cells)"); ax.set_ylabel("$y$ (cells)")
# (b) shear map, (c) dilation map
for k, (F, lab, tag) in enumerate([(eps[0, 1], r"shear strain $\varepsilon_{xy}$", "b"), (th, r"area change $\varepsilon_{xx}+\varepsilon_{yy}$", "c")]):
    ax = fig.add_subplot(gs[0, 1+k]); sl = (R < 9)
    lim = np.percentile(abs(F[(R > 1.2) & (R < 9)]), 97)
    im = ax.imshow(F.T, origin="lower", extent=[x[0], x[-1], x[0], x[-1]], cmap=DIV, vmin=-lim, vmax=lim)
    ax.add_patch(plt.Rectangle((-0.5, -0.5), 1, 1, fc="none", ec=INK, lw=0.8))
    ax.set_xlim(-8, 8); ax.set_ylim(-8, 8); ax.grid(False)
    ax.set_title(f"({tag}) {lab}\n(directional step)"); ax.set_xlabel("$x$ (cells)")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03); cb.formatter.set_powerlimits((0, 0)); cb.ax.yaxis.get_offset_text().set_x(3.2)
# (d) decay
ax = fig.add_subplot(gs[1, 0])
rb = np.geomspace(1.5, 60, 18)
for E_, col, lab in ((eps, C2, "directional step"), (epsi, C1, "isotropic step")):
    mag = np.sqrt((E_**2).sum(axis=(0, 1)))
    prof = np.array([mag[(R > r*0.93) & (R < r*1.07)].mean() for r in rb])
    w = (rb > 5) & (rb < 40); sl = np.polyfit(np.log(rb[w]), np.log(prof[w]), 1)[0]
    ax.loglog(rb, prof, "o", color=col, ms=4.5, label=f"{lab}: slope {sl:.2f} (fit 5–40)")
    print(lab, "slope", sl)
ax.loglog(rb, prof[3]*(rb[3]/rb)**2, "--", color=MUTED, lw=1, label=r"$r^{-2}$ (2D force dipole)")
ax.set_xlabel("distance $r$ from stepped cell (cells)"); ax.set_ylabel(r"$|\varepsilon|$ (ring average)")
ax.set_title("(d) how far the step is felt"); ax.legend(fontsize=7.5)
# (e) dilation/shear outside vs box size
ax = fig.add_subplot(gs[1, 1])
lam, mu = 0.0, 866.6
Ciso = np.zeros((2, 2, 2, 2))
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                Ciso[i, j, k, l] = lam*(i == j)*(k == l) + mu*((i == k)*(j == l) + (i == l)*(j == k))
Ns = [32, 64, 128, 256]; res = {"q": [], "c": []}
for Nc in Ns:
    for Cc, key in ((S.C, "q"), (Ciso, "c")):
        xx, uu, gg, ch = S.solve(np.diag([dP[0, 0]]*2), Cc, Ncell=Nc, p=4)
        ee = 0.5*(gg + gg.transpose(1, 0, 2, 3)); t = ee[0, 0]+ee[1, 1]; dv = np.sqrt(((ee[0, 0]-ee[1, 1])/2)**2+ee[0, 1]**2)
        XX, YY = np.meshgrid(xx, xx, indexing="ij"); RR = np.hypot(XX, YY); ring = (RR > 4) & (RR < 12)
        res[key].append(np.abs(t[ring]).mean()/np.abs(dv[ring]).mean())
print(res)
ax.loglog(Ns, res["q"], "o-", color=C2, label="quantum-cell medium (this model)")
ax.loglog(Ns, res["c"], "s-", color=C1, label="isotropic control medium")
ax.loglog(Ns, res["c"][0]*(Ns[0]/np.array(Ns))**2, "--", color=MUTED, lw=1, label=r"$\propto N^{-2}$ (closed-box compensation)")
ax.set_xlabel("lattice size $N$ (cells per side)"); ax.set_ylabel("neighbours' |area change| / |shear|")
ax.set_xticks(Ns); ax.set_xticklabels([str(v) for v in Ns]); ax.minorticks_off()
ax.set_title("(e) isotropic step, ring $4<r<12$:\nneighbours are mostly sheared"); ax.legend(fontsize=7.5)
# (f) realised stretch of the stepped cell
ax = fig.add_subplot(gs[1, 2])
e11 = eps[0, 0][chi > 0].mean()
vals = [1/n1, e11]; labs = ["iso-energy rule\n$1/n_1$", "realised in\nlattice"]
b = ax.bar(labs, vals, color=[MUTED, C2], width=0.55)
ax.set_yscale("log"); ax.set_ylabel(r"stretch $\varepsilon_{xx}$ of the stepped cell")
for bb, v in zip(b, vals): ax.text(bb.get_x()+bb.get_width()/2, v*1.15, f"{v:.2e}", ha="center", fontsize=8.5, color=INK)
ax.set_title(r"(f) stepped cell, fermion $(7,7)\to(8,7)$" + "\namong $N_f=69$ occupants"); ax.set_ylim(1e-3, 0.5)
print("realised e11", e11, "ratio", e11*n1)
fig.savefig("figs/fig3_embedded_response.png"); print("saved")
