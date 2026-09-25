import numpy as np
from vstyle import *
from sim5_consensus import ck2
import sim6_stress as S6
fig = plt.figure(figsize=(11.4, 7.8)); gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.48)
# (a) energy moments of one sudden expansion step
ax = fig.add_subplot(gs[0, 0]); Ks = np.array([250, 500, 1000, 2000, 4000, 8000, 16000, 32000])
for (n, s), col, mk in (((5, 6/5), C1, "o"), ((12, 13/12), C2, "s"), ((7, 1.37), C3, "^")):
    m1, m2 = [], []
    for K in Ks:
        p = ck2(n, s, K); k = np.arange(1, K+1); E = (k/s)**2
        m1.append((p*E).sum()/n**2); m2.append((p*E**2).sum()/n**4 - ((p*E).sum()/n**2)**2)
    ax.loglog(Ks, m2, mk+"-", color=col, ms=4, lw=1.3, label=f"$n={n}$, $s={s:.3f}$")
    print(n, s, "mean", np.round(m1[-3:], 5))
ax.loglog(Ks, 4.4*Ks/500*0.5, "--", color=MUTED, lw=1, label=r"$\propto K$")
ax.set_xlabel("mode cutoff $K$"); ax.set_ylabel(r"Var$(E)/E^2$ of one spectator")
ax.set_title("(a) one sudden expansion step: mean energy\nexact, variance diverges with the cutoff"); ax.legend(fontsize=7.5)
# (b) + (c) from sim5b
out = np.load("figs/sim5b.npy", allow_pickle=True).item()
ax = fig.add_subplot(gs[0, 1])
for N, col, mk in ((4, C1, "o"), (16, C2, "s"), (64, C3, "^")):
    Ss = [4, 16, 64]; sd = [out[(N, S)][1] for S in Ss]; G = [out[(N, S)][0] for S in Ss]
    ax.plot(G, sd, mk+"-", color=col, ms=6, lw=1.4, label=f"$N={N}$ occupants")
ax.set_xlabel(r"mean growth $\langle \ln L\rangle$ (after $S$ = 4, 16, 64 steps)"); ax.set_ylabel(r"spread std$(\ln L)$ across runs")
ax.set_title("(b) spread vs growth: no " + r"$1/\sqrt{N}$" + " shrinkage\nwith the number of occupants"); ax.legend(fontsize=8)
ax = fig.add_subplot(gs[0, 2])
Eall = np.concatenate([out[(N, 64)][3] for N in (4, 16, 64)])
ax.hist(Eall, bins=np.linspace(0.4, 2.5, 43), color=C1, edgecolor="white", lw=0.5)
ax.axvline(np.median(Eall), color=INK, ls=":", lw=1.2); ax.axvline(Eall.mean(), color=C2, ls="--", lw=1.2)
ax.text(np.median(Eall)-0.03, ax.get_ylim()[1]*0.9, f"median {np.median(Eall):.2f}", ha="right", fontsize=8.5)
ax.text(Eall.mean()+0.03, ax.get_ylim()[1]*0.78, f"mean {Eall.mean():.2f}", ha="left", fontsize=8.5, color=INK)
ax.set_xlabel(r"cell energy $E/E_0$ after 64 steps"); ax.set_ylabel("runs")
ax.set_title(f"(c) iso-energy only on average: typical cells\nlose energy, rare ones are kicked ({len(Eall)} runs)")
from matplotlib.ticker import MaxNLocator
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
print("E after 64 steps: median", np.median(Eall), "mean", Eall.mean(), "max", Eall.max())
# (d) T^xx maps
X, Y = np.meshgrid(S6.xs, S6.xs, indexing="ij")
shell = [(m, n) for m in range(1, 9) for n in range(1, 9) if m*m+n*n <= 50]
for k, (modes, ttl) in enumerate([([(3, 3)], "one occupant, mode (3,3)"), (shell, f"Fermi sea, $N_f={len(shell)}$")]):
    ax = fig.add_subplot(gs[1, k]); T = S6.T_xx(modes); T = T/T.mean()
    lim = np.abs(T).max()
    im = ax.imshow(T.T, origin="lower", extent=[0, 1, 0, 1], cmap=DIV, vmin=-lim, vmax=lim); ax.grid(False)
    ax.set_title(f"({'de'[k]}) local vote " + r"$T^{xx}/\langle T^{xx}\rangle$" + f":\n{ttl}");ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    print(ttl, "T/<T> range", T.min(), T.max())
# (f) exact vs formula along diagonal
ax = fig.add_subplot(gs[1, 2]); eps, w = 0.02, 0.07
line = np.linspace(0.12, 0.88, 40)
for modes, col, lab in (([(3, 3)], C2, "mode (3,3)"), (shell, C1, "Fermi sea")):
    T = S6.T_xx(modes); nrm = None
    f = np.array([-0.5*np.sum(T*eps*S6.bump(t, t, w))*S6.hgrid**2 for t in line])
    nrm = np.abs(f).mean(); ax.plot(line, f/nrm, color=col, lw=1.6, label=f"{lab}: stress formula")
    pts = [0.2, 0.35, 0.5, 0.65, 0.8]; ex = []
    idx = list(range(len(shell))) if len(modes) > 1 else None
    for t in pts:
        hp = lambda X, Y: eps*np.exp(-((X-t)**2+(Y-t)**2)/(2*w*w)); hm = lambda X, Y: -hp(X, Y)
        Ep, Em = S6.spectrum(hp), S6.spectrum(hm); E0 = S6.spectrum(lambda X, Y: 0*X)
        if len(modes) > 1: ex.append((Ep[:len(shell)].sum()-Em[:len(shell)].sum())/2)
        else:
            j = np.argmin(abs(E0 - 0.5*np.pi**2*18)); ex.append((Ep[j]-Em[j])/2)
    ax.plot(pts, np.array(ex)/nrm, "o", color=col, ms=6, mec="white", label=f"{lab}: exact FD eigenvalues")
    print(lab, "exact", np.round(ex, 5), "formula at pts", np.round(np.interp(pts, line, f), 5))
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlabel("position $x=y$ of a small local stretch $h_{xx}$"); ax.set_ylabel(r"$\delta E$ / mean $|\delta E|$")
ax.set_title("(f) energy response to a local stretch:\nsingle occupant vs many (exact check)"); ax.set_ylim(-1.85, 0.2); ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, columnspacing=0.8)
fig.savefig("figs/fig5_consensus.png"); print("saved")
