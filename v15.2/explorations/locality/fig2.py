import numpy as np
from matplotlib.patches import Rectangle
from vstyle import *
from sim2_curvature import fields
s = 4/3
fig = plt.figure(figsize=(11.5, 8.2))
gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1.15, 1.15, 0.07], height_ratios=[1, 0.95], hspace=0.5, wspace=0.45)

def grid(ax, stretch_col=None, stretch_cell=None, title=""):
    for i in range(5):
        for j in range(5):
            x0 = i + (s-1 if (stretch_col is not None and i > stretch_col) else 0)
            wdt = s if (i == stretch_col) else 1
            fc = "#cde2fb" if i == stretch_col else "#f0efec"
            ax.add_patch(Rectangle((x0, j), wdt, 1, fc=fc, ec=INK2, lw=0.8))
    if stretch_cell:
        i, j = stretch_cell
        ax.add_patch(Rectangle((i-(s-1)/2, j), s, 1, fc="#f6c7b3", ec=C2, lw=1.8, alpha=0.9, zorder=3))
        for yy in (j, j+1):
            for xa, xb in ((i-(s-1)/2, i), (i+1, i+1+(s-1)/2)):
                ax.plot([xa, xb], [yy, yy], color="#b3261e", lw=5, zorder=4, solid_capstyle="butt")
    ax.set_xlim(-0.4, 5.8); ax.set_ylim(-0.2, 5.2); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=9.5)

grid(fig.add_subplot(gs[0, 0]), stretch_cell=(2, 2), title="(a) one cell stretched by\n$s=(n{+}1)/n$ along $x$\n(red: face with no partner)")
grid(fig.add_subplot(gs[0, 1]), stretch_col=2, title="(b) conforming gluing:\nthe whole column must\nfollow (a slab)")
x, h, a, Kd, Ki = fields(s, 0.04); lim = 25
for k, (K, tag, ttl) in enumerate([(Kd, "c", "directional step"), (Ki, "d", "isotropic step")]):
    ax = fig.add_subplot(gs[0, 2+k])
    im = ax.imshow(np.clip(K, -lim, lim), extent=[x[0], x[-1], x[0], x[-1]], origin="lower", cmap=DIV, vmin=-lim, vmax=lim)
    ax.plot([-.5, .5, .5, -.5, -.5], [-.5, -.5, .5, .5, -.5], color=INK2, lw=0.6, ls=":")
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.grid(False)
    ax.set_title(f"({tag}) $K$, {ttl}\n(abstract gluing)", fontsize=9.5)
    ax.set_xlabel("$x$ (cell units)"); ax.set_ylabel("$y$")
cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, 4])); cb.set_label(r"$K$ (clipped at $\pm 25$), $s=4/3$, $w=0.04$")

ax = fig.add_subplot(gs[1, 0:2])
for w, col in ((0.02, C1), (0.04, C2), (0.08, C3)):
    x, h, a, Kd, Ki = fields(s, w); j0 = np.argmin(abs(x))
    ax.plot(x, Kd[:, j0]*w**2, color=col, lw=1.6, label=f"$w={w}$")
ax.set_xlim(0.25, 0.75); ax.set_xlabel("$y$ at $x=0$ (across the top face at $y=0.5$)"); ax.set_ylabel(r"$K\,w^2$")
ax.set_title("(e) each face is a curvature double layer:\n$K\,w^2$ collapses, so peak $|K|\propto 1/w^2$", fontsize=9.5); ax.legend()

ax = fig.add_subplot(gs[1, 2:4])
ws = np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.11])
for ss, col, mk in ((4/3, C1, "o"), (6/5, C2, "s"), (11/10, C3, "^")):
    dips = []
    for w in ws:
        x, h, a, Kd, Ki = fields(ss, w, N=1201); j0 = np.argmin(abs(x)); up = x > 0
        dips.append(-np.sum((x[up]-0.5)*Kd[up, j0]*np.exp(a[up, j0]))*h)
    ax.plot(ws, dips, mk+"-", color=col, ms=5, lw=1.4, label=f"$s={ss:.3f}$")
    ax.axhline(ss-1, color=col, lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_xlabel("smoothing width $w$ of the cell face")
ax.set_ylabel(r"$-\int (y-y_f)\,K\sqrt{g}\,dy$ per unit length")
ax.set_title("(f) width-independent content of a face: dipole moment $= s-1$\n(dashed), the edge-length mismatch per unit length", fontsize=9.5)
ax.legend(title="step ratio", loc="upper right", ncol=3); ax.set_ylim(0, 0.42)
fig.savefig("figs/fig2_curvature.png"); print("saved")
