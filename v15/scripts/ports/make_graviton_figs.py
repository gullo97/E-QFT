"""Figure 14: the graviton measurement vs the one-loop coefficient, and anim of the wave."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
out = json.load(open("graviton_refined.json"))

# ---------- panel A: gradient part of Pi(q) at s=3.0 with EH prediction lines
s = "3.0"
C = out[s]["C"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6))
style = {"conf": ("o", "crimson"), "TT+": ("s", "seagreen"),
         "TTx": ("D", "olive"), "diffeo": ("^", "steelblue")}
for m, (mk, col) in style.items():
    R = np.array(out[s]["modes"][m]["rows"])
    X = np.vstack([np.ones(len(R)), R[:, 1], R[:, 1]**2]).T
    coef, *_ = np.linalg.lstsq(X, R[:, 2], rcond=None)
    a1.plot(R[:, 1], (R[:, 2] - coef[0])/C, mk, color=col, ms=7,
            label=fr"{m}:  $\kappa/C = {out[s]['modes'][m]['kappa_over_C']:+.2f}$")
qq = np.linspace(0, 0.26, 50)
a1.plot(qq, -qq, "k--", lw=1.2)
a1.plot(qq, +qq, "k:", lw=1.2)
a1.text(0.165, -0.225, r"$-C(s)\,q^2$  (EH, TT)", fontsize=9, rotation=-13)
a1.text(0.13, +0.165, r"$+C(s)\,q^2$  (EH, conf)", fontsize=9, rotation=13)
a1.set_xlabel("$q^2$"); a1.set_ylabel(r"$[\Pi(q)-\Pi_0]\;/\;C(s)$")
a1.set_title("induced gradient couplings at $s=3$ vs the\none-loop Einstein–Hilbert prediction $C(s)=1/48\\pi^2 s^2$")
a1.legend(fontsize=8, loc="upper left")
a1.set_ylim(-0.55, 1.0); a1.set_xlim(0, 0.26)

# ---------- panel B: convergence of kappa/C with the regulator scale
ss = np.array([2.0, 2.5, 3.0])
for m, (mk, col) in style.items():
    vals = [out[str(x)]["modes"][m]["kappa_over_C"] for x in ss]
    a2.plot(ss, vals, mk+"-", color=col, ms=8, label=m)
a2.axhline(-1, color="k", ls="--", lw=1.2)
a2.axhline(0, color="0.6", lw=0.8)
a2.text(2.32, -0.93, "graviton: $\\kappa_{TT} \\to -C(s)$", fontsize=10)
confvals = [out[str(x)]["modes"]["conf"]["kappa_over_C"] for x in ss]
a2.annotate("conf (off scale, opposite sign):\n$+14.2 \\to +7.7 \\to +3.8$\ndescending toward $+1$",
            xy=(2.03, -0.30), fontsize=9, color="crimson")
a2.set_xlabel("regulator proper time $s$  (covariance improves $\\to$)")
a2.set_ylabel(r"$\kappa\,/\,C(s)$")
a2.set_title("convergence to the universal coefficient:\nTT$+$, TT$\\times$ $\\to$ $-1$;  diffeo control $\\to$ 0")
a2.set_ylim(-1.7, 0.45)
a2.legend(fontsize=8, loc="lower right")
fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.12, wspace=0.25); fig.savefig("figs/fig14_graviton.png", dpi=150); plt.close(fig)
print("saved fig14")

# ---------- animation: a + polarized wave crossing the cell lattice
NX, NY = 22, 12
amp, q = 0.28, 2*np.pi*3/NX           # exaggerated amplitude for visibility
frames = 56
fig, ax = plt.subplots(figsize=(8.2, 4.6))
ax.set_xlim(-0.8, NX - 0.2); ax.set_ylim(-0.8, NY - 0.2)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("a transverse-traceless wave on the cell lattice:  $h_+=\\varepsilon\\cos(qx-\\omega t)$\n"
             "cells $1{+}h$ by $1{-}h$ — stiffness $\\kappa_{TT}=-C(s)$ measured from matter loops",
             fontsize=10)
from matplotlib.patches import Rectangle
rects = []
for i in range(NX):
    for j in range(NY):
        r = Rectangle((i, j), 0.86, 0.86, facecolor="#f5e6c8",
                      edgecolor="#8a5a2b", lw=1.0)
        ax.add_patch(r); rects.append((i, j, r))
cmap = plt.cm.RdBu_r

def update(f):
    t = f/frames*2*np.pi
    for i, j, r in rects:
        h = amp*np.cos(q*i - 2*t)
        w, hh = 0.86*(1 + h), 0.86*(1 - h)
        r.set_bounds(i + 0.43 - w/2, j + 0.43 - hh/2, w, hh)
        r.set_facecolor(cmap(0.5 + 0.45*h/amp*0.9))
    return [r for _, _, r in rects]

anim = FuncAnimation(fig, update, frames=frames, blit=False)
anim.save("figs/anim_graviton.gif", writer=PillowWriter(fps=14), dpi=85)
print("saved figs/anim_graviton.gif")
