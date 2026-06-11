"""Animations of the emergent gravitational dynamics (from kasner_traj.npz)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Line3DCollection

d = np.load("kasner_traj.npz")

def cuboid_edges(sx, sy, sz):
    v = np.array([[x, y, z] for x in (-sx/2, sx/2) for y in (-sy/2, sy/2) for z in (-sz/2, sz/2)])
    E = [(0,1),(2,3),(4,5),(6,7),(0,2),(1,3),(4,6),(5,7),(0,4),(1,5),(2,6),(3,7)]
    return [(v[i], v[j]) for i, j in E]

def make_anim(name, al, times, S2, p, title, kasner_target=None, shape_norm=False,
              S2_panel=True, frames_step=2):
    al0 = al[:, :1]
    F = range(0, al.shape[1], frames_step)
    fig = plt.figure(figsize=(10.2, 4.2))
    ax3 = fig.add_subplot(1, 3, (1, 2), projection="3d")
    ax2 = fig.add_subplot(1, 3, 3)
    if S2_panel:
        ax2.plot(times, S2, color="gray", lw=1)
        ax2.axhline(1, color="k", ls="--", lw=0.9, label=r"vacuum Kasner $\Sigma p^2=1$")
        ax2.axhline(1/3, color="k", ls=":", lw=0.9, label="isotropic $1/3$")
        dot2, = ax2.plot([], [], "o", color="crimson", ms=8)
        ax2.set_xlabel("model time"); ax2.set_ylabel(r"$\Sigma\, p_a^2$")
        ax2.legend(fontsize=8)
    else:
        P1, P2 = np.meshgrid(np.linspace(-0.6, 1.2, 300), np.linspace(-0.6, 1.2, 300))
        ax2.contour(P1, P2, P1**2 + P2**2 + (1-P1-P2)**2, levels=[1.0],
                    colors="k", linestyles="--")
        ax2.plot([1/3], [1/3], "k*", ms=10)
        if kasner_target is not None:
            ax2.plot([kasner_target[0]], [kasner_target[1]], "s", mfc="none",
                     color="seagreen", ms=10)
        dot2, = ax2.plot([], [], "o", color="crimson", ms=8)
        ax2.plot(p[0], p[1], "-", color="crimson", lw=0.8, alpha=0.5)
        ax2.set_xlabel(r"$p_1$"); ax2.set_ylabel(r"$p_2$"); ax2.set_aspect("equal")
        ax2.set_title("Kasner plane")
    lc = Line3DCollection(cuboid_edges(1.0, 1.0, 1.0), colors="tab:blue", lw=2)
    ax3.add_collection3d(lc)
    R = 3.2
    ax3.set_xlim(-R, R); ax3.set_ylim(-R, R); ax3.set_zlim(-R, R)
    ax3.set_box_aspect((1, 1, 1)); ax3.set_title(title, fontsize=11)
    ax3.set_xticks([]); ax3.set_yticks([]); ax3.set_zticks([])
    txt = ax3.text2D(0.02, 0.92, "", transform=ax3.transAxes, fontsize=9)

    def update(fi):
        if shape_norm:                       # show the expansion-RATE box (Hubble ellipsoid)
            s = np.clip(p[:, fi], 0.02, None)
            s = s/np.prod(s)**(1/3)
        else:
            s = np.exp(al[:, fi] - al0[:, 0])
        s = np.clip(s, 0.02, 2*R)
        lc.set_segments(cuboid_edges(*s))
        if S2_panel:
            dot2.set_data([times[fi]], [S2[fi]])
        else:
            dot2.set_data([p[0, fi]], [p[1, fi]])
        txt.set_text(f"t = {times[fi]:.2f}   sides = ({s[0]:.2f}, {s[1]:.2f}, {s[2]:.2f})")
        return lc, dot2, txt

    anim = FuncAnimation(fig, update, frames=list(F), blit=False)
    anim.save(f"figs/{name}.gif", writer=PillowWriter(fps=12), dpi=88)
    plt.close(fig)
    print(f"saved figs/{name}.gif  ({len(list(F))} frames)")

# A: generator coupling, vacuum-Kasner point (-2/7, 3/7, 6/7): axis 1 contracts
make_anim("anim_kasner_epoch", d["al_g"], d["times_g"], d["S2_g"], d["p_g"],
          r"Kasner epoch from the box model: $p=(-\frac{2}{7},\frac{3}{7},\frac{6}{7})$"
          "\n(one axis contracts while two expand)",
          kasner_target=(-2/7, 3/7), S2_panel=False, frames_step=2)

# B: v13.3 hopping: shape isotropization (volume normalised out)
make_anim("anim_isotropization", d["al_i"], d["times"], d["S2_i"], d["p_i"],
          "v13.3 hopping: expansion-rate anisotropy decays\n(box sides = relational rates $p_a$; cube = isotropic FRW)",
          shape_norm=True, S2_panel=True, frames_step=2)
