"""Animation: the Mixmaster universe -- Kasner epochs hopping around the circle."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

d = np.load("mixmaster_traj.npz")
t, A, ptau, u = d["t"], d["A"], d["ptau"], d["u"]
# restrict to the well-resolved era (through the 5 epochs)
kmax = np.searchsorted(t, 1250.0)
t, A, ptau, u = t[:kmax], A[:, :kmax], ptau[:, :kmax], u[:kmax]
frames = list(range(0, kmax, kmax//110))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.4),
                               gridspec_kw={"width_ratios": [1.35, 1]})
cols = ["#c0392b", "#2980b9", "#27ae60"]
for i in range(3):
    ax1.plot(t, A[i], color=cols[i], lw=1.4, label=fr"$\alpha_{i+1}$")
cursor = ax1.axvline(0, color="k", lw=1)
ax1.set_xlabel(r"BKL time $\tau$"); ax1.set_ylabel(r"$\alpha_i=\ln a_i$")
ax1.set_title("Mixmaster: free Kasner flights + curvature bounces")
ax1.legend(fontsize=8, loc="lower left")
P1, P2 = np.meshgrid(np.linspace(-0.6, 1.2, 300), np.linspace(-0.6, 1.2, 300))
ax2.contour(P1, P2, P1**2 + P2**2 + (1 - P1 - P2)**2, levels=[1.0],
            colors="k", linestyles="--")
ax2.plot([1/3], [1/3], "k*", ms=8)
trail, = ax2.plot([], [], "-", color="0.7", lw=0.8)
dot, = ax2.plot([], [], "o", color="crimson", ms=9)
txt = ax2.text(0.03, 0.95, "", transform=ax2.transAxes, fontsize=9)
ax2.set_xlabel("$p_1$"); ax2.set_ylabel("$p_2$"); ax2.set_aspect("equal")
ax2.set_title("exponents hop around the Kasner circle")

def update(k):
    cursor.set_xdata([t[k]])
    trail.set_data(ptau[0, :k], ptau[1, :k])
    dot.set_data([ptau[0, k]], [ptau[1, k]])
    txt.set_text(fr"$u = {u[k]:.3f}$")
    return cursor, trail, dot, txt

anim = FuncAnimation(fig, update, frames=frames, blit=False)
anim.save("figs/anim_mixmaster.gif", writer=PillowWriter(fps=12), dpi=88)
print(f"saved figs/anim_mixmaster.gif ({len(frames)} frames)")
