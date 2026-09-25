import numpy as np
from vstyle import *
from q2_analyse import speeds
res = np.load("figs/q2_production.npy", allow_pickle=True).item()
fig = plt.figure(figsize=(11.4, 7.8)); gs = fig.add_gridspec(2, 3, hspace=0.5, wspace=0.34)
maps = [((1, 0.01), "(a) massless occupants, light walls ($\mu=0.01$)"), ((1, 5.0), "(b) massless occupants, heavy walls ($\mu=5$)"),
        ((0, 0.05), "(c) massive occupants, light walls ($\mu=0.05$)")]
for k, (key, ttl) in enumerate(maps):
    r = res[key]; ax = fig.add_subplot(gs[0, k]); T = r["T"]; m = r["mean"]
    lim = 0.25
    im = ax.imshow(m.T, origin="lower", aspect="auto", extent=[0, T, 0.5, m.shape[1]+0.5], cmap=DIV, vmin=-lim, vmax=lim)
    tt = np.linspace(0, T, 50)
    ax.plot(tt, tt, color=INK, lw=0.9, ls="--"); ax.text(T*0.3, min(T*0.3, 140)+8, "$j=ct$" if key[0] == 1 else "$j=v_p t$", fontsize=8.5)
    mu = key[1]; cs = np.sqrt(1/(1+mu)) if key[0] == 1 else np.sqrt(3/(1+mu))
    if abs(cs-1) > 0.05: ax.plot(tt, cs*tt, color=INK, lw=0.9, ls=":"); ax.text(T*0.45, cs*T*0.45+6, "$c_s t$ (prediction)", fontsize=8.5)
    ax.set_ylim(0.5, 150.5); ax.set_xlabel("time (cell lengths / particle speed)"); ax.set_ylabel("wall index $j$ from kicked cell")
    ax.set_title(ttl, fontsize=9.5); ax.grid(False)
cb = fig.colorbar(im, ax=fig.axes[:3], fraction=0.02, pad=0.01); cb.set_label(r"mean wall displacement $\langle\delta X_j\rangle$")
ax = fig.add_subplot(gs[1, 0:2])
mus = np.array(sorted({k[1] for k in res})); mm = np.geomspace(0.008, 130, 200)
for kind, col, lab, f in ((1, C1, "massless occupants (energy $E$ per cell)", lambda x: np.sqrt(1/(1+x))),
                          (0, C2, "massive occupants (mass $km$ per cell)", lambda x: np.sqrt(3/(1+x)))):
    vs = np.array([speeds(res[(kind, mu)])[0] for mu in mus]); vc = np.array([speeds(res[(kind, mu)])[1] for mu in mus])
    ax.loglog(mus, vs, "o", color=col, ms=6, mec="white", label=f"{lab}: pulse speed")
    ax.loglog(mus, vc, "^", color=col, ms=5, mfc="none", label="  leading edge ($|\delta X|>10^{-9}$)")
    ax.loglog(mm, f(mm), color=col, lw=1.4)
ax.loglog([1.0], [0.722], "s", color=C1, ms=6, mec=INK, label="small kick (30%)"); ax.loglog([2.0], [1.020], "s", color=C2, ms=6, mec=INK)
ax.axhline(1, color=INK, lw=0.8, ls="--")
ax.text(0.009, 1.06, "$c$ (massless) = $v_p$ (massive) = 1", fontsize=8.5, color=INK2)
ax.set_xlabel(r"wall mass $\mu = M_w c^2/2E_\mathrm{cell}$ (massless)  or  $M_w/km$ (massive)")
ax.set_ylabel("speed (units of $c$, resp. $v_p$)")
ax.set_title("(d) disturbance speed vs wall mass; lines: $c/\sqrt{1+\mu}$ (enthalpy) and $v_p\sqrt{3/(1+\mu)}$ (adiabatic)")
ax.legend(fontsize=7.5, loc="lower left", ncol=2)
# (e) causal step rate bound
import q1_rate_bound as RB
ax = fig.add_subplot(gs[1, 2]); xs = np.geomspace(0.01, 64, 14)
b = [RB.run(16, x/16) for x in xs]; nv = [RB.run(16, x/16, rule="add") for x in xs]
ax.loglog(xs, nv, "s", color=MUTED, ms=5, label="naive: steps superpose")
ax.loglog(xs, b, "o", color=C1, ms=5, mec="white", label="causal: one step at a time")
xx = np.geomspace(0.01, 64, 100); ax.loglog(xx, xx/(1+xx), color=C1, lw=1.2); ax.axhline(1, color=INK, lw=0.8, ls="--")
ax.set_xlabel(r"total step rate $\times\,t_0$  ($N\Gamma t_0$)"); ax.set_ylabel("mean wall speed / $c$")
ax.set_title("(e) steps of duration $t_0=\Delta L/c$ cap\ncell growth at $c$ (line: $x/(1+x)$)"); ax.legend(fontsize=7.5)
fig.savefig("figs/figB_wall_mass.png"); print("saved")
