import numpy as np
from vstyle import *
d = np.load("figs/sim4.npz"); ts, H1, H2, c, cs, n0 = d["ts"], d["H1"], d["H2"], int(d["c"]), float(d["cs"]), int(d["n0"])
fig = plt.figure(figsize=(11.2, 7.2)); gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.28)
ax = fig.add_subplot(gs[0, 0]); win = 240
lim = 0.02
im = ax.imshow((H2[:, c-win:c+win+1]-1).T, aspect="auto", origin="lower", cmap=DIV, vmin=-lim, vmax=lim,
               extent=[ts[0], ts[-1], -win-0.5, win+0.5])
ax.plot(ts, cs*ts, color=INK, lw=0.8, ls="--"); ax.plot(ts, -cs*ts, color=INK, lw=0.8, ls="--")
ax.text(40, cs*40+14, r"$\pm c_s t$", fontsize=9, color=INK)
ax.set_ylim(-win, win); ax.grid(False)
ax.set_xlabel("time"); ax.set_ylabel("cell index relative to stepped cell")
ax.set_title("(a) $L_i(t)-L$ after the thesis-literal step (P2)\nblue = compressed, red = stretched; clipped at $\pm 0.02$")
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
ax = fig.add_subplot(gs[0, 1])
idx = np.arange(-win, win+1)
for tsel, col in ((5, C1), (20, C2), (45, C3)):
    k = np.argmin(abs(ts-tsel)); ax.plot(idx, H2[k, c-win:c+win+1]-1, color=col, lw=1.2, label=f"t = {ts[k]:.0f}")
ax.set_ylim(-0.012, 0.012); ax.set_xlabel("cell index relative to stepped cell"); ax.set_ylabel("$L_i - L$")
ax.set_title("(b) snapshots: compression fronts at " + r"$\pm c_s t$" + "\nplus a cell-scale ringing wake inside the cone"); ax.legend(loc="lower left")
ax = fig.add_subplot(gs[1, 0])
ax.plot(ts, H1[:, c], color=C1, lw=1.4, label="P1: promote at fixed walls, release")
ax.plot(ts, H2[:, c], color=C2, lw=1.4, label=r"P2: thesis step $L\to L(n{+}1)/n$, release")
ax.axhline((n0+1)/n0, color=MUTED, ls="--", lw=1); ax.text(59, (n0+1)/n0+0.002, "iso-energy rule $(n{+}1)/n$", ha="right", fontsize=8.5, color=INK2)
ax.axhline(((n0+1)/n0)**(2/3), color=INK, ls=":", lw=1.2); ax.text(59, ((n0+1)/n0)**(2/3)-0.006, "pressure balance $((n{+}1)/n)^{2/3}$", ha="right", fontsize=8.5, color=INK2)
ax.set_xlabel("time"); ax.set_ylabel("stepped cell size $L_c/L$"); ax.set_title(f"(c) the stepped cell settles by pressure balance ($n={n0}$)")
ax.legend(loc="lower right", fontsize=8); ax.set_ylim(0.99, 1.115)
ax = fig.add_subplot(gs[1, 1])
k = np.argmin(abs(ts-45))
disp = np.cumsum(H2[k]-1)            # wall displacement = cumulative size change from the left end
ax.plot(np.arange(len(disp))-c, disp, color=C2, lw=1.5)
ax.axhline(0.5*(((n0+1)/n0)**(2/3)-1), color=INK, ls=":", lw=1.2)
ax.text(-290, 0.5*(((n0+1)/n0)**(2/3)-1)+0.001, "half of the settled excess length", fontsize=8.5, color=INK2)
ax.set_xlim(-300, 300); ax.set_xlabel("wall index relative to stepped cell"); ax.set_ylabel("wall displacement at t = 45")
ax.set_title("(d) what X+dX keeps: a rigid translation (not a size\nchange) inside the light-cone, nothing outside it")
fig.savefig("figs/fig4_pulse.png"); print("saved")
