"""Animation: emergent gravitational collapse in the separate-universe cell chain."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

d = np.load("collapse_traj.npz")
t, A, abg, dcen, delta0, ic, amin = (d[k] for k in
                                     ("t", "A", "abg", "dcen", "delta0", "ic", "amin"))
a0 = abg[0]
Mc = A.shape[1]
frames = list(range(0, len(t), 5))
cmap = plt.cm.inferno_r

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.3),
                               gridspec_kw={"width_ratios": [1.35, 1]})
ax2.set_xscale("log"); ax2.set_yscale("log")
hit = A[:, Mc//2] <= amin*1.001
kcol = int(np.argmax(hit)) if hit.any() else len(t)-1
dplot = np.where(np.arange(len(t)) <= kcol, dcen, np.nan)
ax2.plot(abg/a0, dplot, color="0.6", lw=1.2)
wfit = (abg/a0 > 30) & (abg/a0 < 300)
Cfit = np.exp(np.mean(np.log(dcen[wfit]) - np.log(abg[wfit])))
ax2.plot(abg/a0, Cfit*abg, "k--", lw=1, label=r"$\delta\propto a$")
ax2.axhline(1.686, color="gray", ls=":", lw=1)
ax2.set_ylim(2e-2, 5e2)
dot, = ax2.plot([], [], "o", color="crimson", ms=8)
ax2.set_xlabel(r"$a_{\rm bg}/a_0$"); ax2.set_ylabel(r"central $\delta$")
ax2.set_title("density contrast growth"); ax2.legend(fontsize=8, loc="upper left")

def update(k):
    ax1.clear()
    w = A[k]/abg[k]                                  # comoving cell widths
    delta = (1 + delta0)*(abg[k]/A[k])**3 - 1
    edges = np.concatenate([[0.0], np.cumsum(w)])
    edges -= edges[-1]/2
    lv = np.log10(1 + np.maximum(delta, 0))
    cols = cmap(0.08 + 0.9*(lv - lv.min())/(np.ptp(lv) + 1e-12))
    for i in range(Mc):
        ax1.fill_between([edges[i], edges[i+1]], 0, 1, color=cols[i],
                         edgecolor="white", lw=0.4)
    ax1.set_xlim(-Mc/2 - 1, Mc/2 + 1); ax1.set_ylim(0, 1)
    ax1.set_yticks([])
    ax1.set_xlabel("comoving position (cell widths $a_i/a_{\\rm bg}$)")
    ax1.set_title("emergent attraction: walls fall toward the overdensity\n"
                  f"$a_{{\\rm bg}}/a_0 = {abg[k]/a0:7.1f}$,   "
                  f"central $\\delta = {min(delta[Mc//2], 9999):8.2f}$   (colour: relative overdensity)")
    dot.set_data([abg[k]/a0], [float(np.clip(dcen[k], 2.5e-2, 4.5e2))])
    return dot,

anim = FuncAnimation(fig, update, frames=frames, blit=False)
anim.save("figs/anim_collapse.gif", writer=PillowWriter(fps=11), dpi=88)
print(f"saved figs/anim_collapse.gif ({len(frames)} frames)")
