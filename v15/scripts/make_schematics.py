#!/usr/bin/env python3
"""
make_schematics.py — all drawn (concept) figures of the thesis, generated
programmatically so that even schematics are reproducible (P7).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

BLUE, ORANGE, GREEN, GRAY, PURPLE = "#1f77b4", "#d95f00", "#2ca02c", "#7f7f7f", "#9467bd"

def box(ax, x, y, w, h, text, fc="#eef4fb", ec=BLUE, fs=8, lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)

def arrow(ax, x0, y0, x1, y1, color="k", lw=1.2, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=12, color=color, lw=lw, linestyle=ls))

def blank(figsize=(7.2, 4.4)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    return fig, ax

def save(fig, name):
    fig.savefig(FIGS / name)
    plt.close(fig)
    print("  ->", name)

# ---------------------------------------------------------------- ch01 fig1: dependency map
def ch01_fig1():
    fig, ax = blank((8.6, 5.2))
    pos = {1: (0.5, 4.6), 2: (0.5, 3.8), 3: (1.6, 3.8), 4: (2.7, 3.8), 6: (3.8, 3.8), 7: (4.9, 3.8),
           5: (2.7, 4.6), 8: (0.5, 2.6), 9: (1.6, 2.6), 10: (2.7, 2.6), 11: (3.8, 2.6),
           12: (4.9, 2.6), 13: (6.0, 2.6), 14: (7.1, 2.6), 15: (6.0, 3.8),
           16: (0.5, 1.4), 17: (1.6, 1.4), 18: (2.7, 1.4), 19: (3.8, 1.4), 20: (4.9, 1.4),
           21: (6.0, 1.4), 22: (7.1, 1.4), 23: (2.7, 0.5), 24: (3.8, 0.5), 25: (4.9, 0.5),
           26: (6.0, 0.5), 27: (7.1, 0.5)}
    inter = {5, 8, 13, 16}
    chain = [(2, 3), (3, 4), (4, 6), (6, 7), (3, 9), (8, 9), (9, 10), (10, 11), (11, 12),
             (12, 14), (13, 14), (12, 15), (3, 17), (16, 17), (17, 18), (18, 19), (19, 20),
             (20, 21), (21, 22), (22, 23), (22, 24), (22, 25), (25, 26), (14, 26), (26, 27)]
    for a, b in chain:
        xa, ya = pos[a]; xb, yb = pos[b]
        arrow(ax, xa + 0.42, ya + 0.13, xb - 0.03, yb + 0.13, color=GRAY, lw=0.9)
    for n, (x, y) in pos.items():
        fc = "#fdf3e7" if n in inter else "#eef4fb"
        ec = ORANGE if n in inter else BLUE
        box(ax, x - 0.03, y, 0.45 + 0.0, 0.28, f"{n}", fc=fc, ec=ec, fs=8)
    ax.text(0.5, 5.15, "Part I: the expanding box", fontsize=9, color=BLUE)
    ax.text(0.5, 3.25, "Part II: charge from geometry", fontsize=9, color=BLUE)
    ax.text(0.5, 2.0, "Part III: emergent gravity", fontsize=9, color=BLUE)
    ax.text(6.4, 4.6, "orange = interludes\n(defer until cited)", fontsize=8, color=ORANGE)
    ax.set_xlim(0, 8.2), ax.set_ylim(0.2, 5.4)
    ax.set_title("Chapter dependency map")
    save(fig, "ch01_fig1_dependency_map.png")

# ---------------------------------------------------------------- ch01 fig2 / ch26 fig1: three faces
def three_faces(name, labeled):
    fig, ax = blank((7.6, 5.4))
    box(ax, 2.45, 2.3, 2.9, 1.0,
        "matter on a dynamical,\ndiscrete geometry whose\nboundaries carry the CP physics",
        fc="#fdf3e7", ec=ORANGE, fs=9)
    box(ax, 0.3, 4.3, 2.2, 0.7, "arrow of expansion\n(Part I)", fs=9)
    box(ax, 5.3, 4.3, 2.2, 0.7, "baryon asymmetry\n(Part II)", fs=9)
    box(ax, 2.85, 0.4, 2.2, 0.7, "emergent gravity\n(Part III)", fs=9)
    arrow(ax, 3.3, 3.35, 1.6, 4.25, color=GRAY)
    arrow(ax, 4.5, 3.35, 6.2, 4.25, color=GRAY)
    arrow(ax, 3.95, 2.25, 3.95, 1.15, color=GRAY)
    if labeled:
        ax.text(1.9, 3.95, "Expansion-Preference Thm (Ch. 6)\nmany-body overlap measure", fontsize=7.5, ha="center")
        ax.text(6.0, 3.95, "Spectral-Flow Master Thm (Ch. 9)\n+ quantized pumping (Ch. 12, 14)", fontsize=7.5, ha="center")
        ax.text(4.05, 1.6, "Dictionary Thm (Ch. 17), Kasner selection (Ch. 20),\ninduced stiffness + graviton (Ch. 22, 25)", fontsize=7.5, ha="left")
        ax.text(3.9, 5.25, "the lock between the arrows: ECPT [Conjecture] (Ch. 7)", fontsize=8, ha="center", color=PURPLE)
        arrow(ax, 2.5, 4.65, 5.3, 4.65, color=PURPLE, lw=1.0, ls=":")
    ax.set_xlim(0, 7.8), ax.set_ylim(0.1, 5.5)
    ax.set_title("Three faces of one structure" + (" — closed" if labeled else ""))
    save(fig, name)

# ---------------------------------------------------------------- ch07 figs
def ch07_fig1():
    fig, ax = blank((7.0, 4.6))
    items = [("$E$", 0.6, "broken (Ch. 6)", ORANGE), ("$C$", 1.8, "broken (SM)", ORANGE),
             ("$P$", 3.0, "broken (SM)", ORANGE), ("$T$", 4.2, "broken (SM)", ORANGE),
             ("$CP$", 5.4, "broken, weak (SM)", ORANGE)]
    for lbl, x, st, c in items:
        box(ax, x, 2.9, 1.0, 0.6, f"{lbl}\n{st}", fc="#fdf3e7", ec=c, fs=8)
    box(ax, 1.8, 1.6, 1.6, 0.6, "$CPT$\nexact (QFT theorem)", fc="#eaf6ea", ec=GREEN, fs=8)
    box(ax, 4.0, 1.6, 1.6, 0.6, "$ECPT$\nexact? [Conjecture]", fc="#f3eefb", ec=PURPLE, fs=8)
    box(ax, 2.9, 0.3, 1.8, 0.6, "$\\mathbb{Z}_2^4$, order 16\n(the full group)", fc="#eef4fb", fs=8)
    arrow(ax, 3.8, 2.25, 3.6, 2.85, color=GRAY); arrow(ax, 4.8, 2.25, 4.8, 2.85, color=GRAY)
    arrow(ax, 3.5, 0.95, 2.8, 1.55, color=GRAY); arrow(ax, 4.3, 0.95, 4.8, 1.55, color=GRAY)
    ax.set_xlim(0, 7), ax.set_ylim(0, 4.2)
    ax.set_title("The ECPT lattice: what is broken, what is exact, what is conjectured")
    save(fig, "ch07_fig1_ecpt_lattice.png")

def ch07_fig2():
    fig, ax = blank((6.8, 4.2))
    arrow(ax, 1.0, 1.0, 5.6, 3.4, color=PURPLE, lw=2.0)
    ax.text(3.3, 2.5, "$ECPT$ exact (diagonal)", rotation=27, fontsize=9, color=PURPLE)
    arrow(ax, 1.0, 1.0, 5.6, 1.0, color=ORANGE, lw=1.6)
    ax.text(3.3, 0.7, "$E$ broken: expansion arrow", fontsize=9, color=ORANGE)
    arrow(ax, 1.0, 1.0, 1.0, 3.4, color=BLUE, lw=1.6)
    ax.text(0.75, 2.2, "$CP$ broken: matter arrow", rotation=90, fontsize=9, color=BLUE)
    ax.text(3.4, 3.8, r"$\Delta Q(+\delta) = -\Delta Q(-\delta)$ (verified to $10^{-15}$, Ch. 11)",
            fontsize=8, ha="center")
    ax.set_xlim(0, 6.8), ax.set_ylim(0, 4.2)
    ax.set_title("One exact diagonal, two broken projections")
    save(fig, "ch07_fig2_two_arrows.png")

# ---------------------------------------------------------------- ch08 figs
def ch08_fig2():
    fig, ax = blank((7.0, 3.6))
    ax.add_patch(Rectangle((1.2, 0.7), 4.4, 2.0, fill=False, lw=2.0, ec=BLUE))
    arrow(ax, 1.2, 1.7, 0.4, 1.7, color=ORANGE, lw=1.6)
    arrow(ax, 5.6, 1.7, 6.4, 1.7, color=ORANGE, lw=1.6)
    ax.text(0.55, 2.0, r"$n = -\hat x$", fontsize=9)
    ax.text(5.8, 2.0, r"$n = +\hat x$", fontsize=9)
    ax.text(1.35, 0.9, r"$-\sigma_2\psi = \psi$" + "\nray $(1, -i)$", fontsize=9)
    ax.text(4.1, 0.9, r"$\sigma_2\psi = \psi$" + "\nray $(1, +i)$", fontsize=9)
    ax.text(3.4, 2.35, r"first-order system $\Rightarrow$ one condition per wall", fontsize=9, ha="center")
    ax.text(3.4, 1.65, r"$n_\mu\bar\psi\gamma^\mu\psi = 0$: nothing leaks", fontsize=9, ha="center", color=GRAY)
    ax.set_xlim(0, 7), ax.set_ylim(0.4, 3.0)
    ax.set_title("Walls as rays: the MIT boundary geometry")
    save(fig, "ch08_fig2_boundary_geometry.png")

def ch08_fig3():
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.8), sharey=True)
    for k, (ax, ttl) in enumerate(zip(axes, ["the vacuum $|0\\rangle$: sea filled",
                                             "one particle + one antiparticle (hole)"])):
        for E in np.arange(-4.6, -0.4, 0.55):
            ax.hlines(E, 0.2, 0.8, color=ORANGE, lw=2.2)
        for E in np.arange(0.6, 4.8, 0.55):
            ax.hlines(E, 0.2, 0.8, color=BLUE, lw=1.0)
        ax.axhline(0, color="k", ls="--", lw=0.7)
        ax.axhspan(-0.45, 0.45, color="#dddddd", alpha=0.4)
        ax.set_xticks([]); ax.set_title(ttl, fontsize=9)
        if k == 1:
            ax.plot([0.5], [1.7], "o", ms=9, color=BLUE)
            ax.plot([0.5], [-1.5], "o", ms=11, mfc="white", mec=ORANGE, mew=1.6)
            ax.text(0.62, 1.65, "particle", fontsize=8)
            ax.text(0.62, -1.6, "hole = antiparticle", fontsize=8)
            ax.hlines(-1.5, 0.2, 0.8, color="#f5c9a8", lw=2.2, zorder=0)
    axes[0].set_ylabel("$E$")
    fig.savefig(FIGS / "ch08_fig3_dirac_sea.png"); plt.close(fig)
    print("  -> ch08_fig3_dirac_sea.png")

# ---------------------------------------------------------------- ch10 fig1
def ch10_fig1():
    fig, ax = blank((7.6, 3.8))
    box(ax, 0.3, 1.4, 2.6, 1.2, "bulk phase $\\delta$\nmass $m\\,e^{i\\delta\\gamma_5}$\nMIT walls $(0, 0)$", fs=9)
    box(ax, 4.7, 1.4, 2.6, 1.2, "real mass $m$\nchiral walls\n$(\\theta_0 + \\delta,\\ \\theta_L + \\delta)$", fs=9)
    arrow(ax, 2.95, 2.0, 4.65, 2.0, color=GREEN, lw=1.6)
    ax.text(3.8, 2.25, r"chiral rotation $\psi \to e^{-i\delta\gamma_5/2}\psi$", fontsize=9, ha="center", color=GREEN)
    ax.text(3.8, 0.9, r"both walls shift EQUALLY: a bulk phase can never make the walls differ",
            fontsize=9, ha="center")
    ax.text(3.8, 0.45, r"$\Rightarrow$ everything CP-physical lives in $\Delta = \theta_L - \theta_0$ (Ch. 12)",
            fontsize=9, ha="center", color=ORANGE)
    ax.set_xlim(0, 7.6), ax.set_ylim(0.2, 3.2)
    ax.set_title("Where the phase can live: a change of address, not of physics")
    save(fig, "ch10_fig1_chiral_rotation_flow.png")

# ---------------------------------------------------------------- ch13 figs
def ch13_fig1():
    fig, ax = blank((7.8, 4.4))
    rows = [("B violation", "sphalerons convert\nany charge excess\n($C_{\\rm sph} = 28/79$)",
             "spectral-flow pumping:\nquantized $\\Delta Q$ per crossing\n(Ch. 9, 12)"),
            ("C & CP violation", "CKM phase $J_{CP}$,\nGIM-crushed in bulk",
             "wall-angle mismatch $\\Delta$;\nmagnitude = $\\varepsilon_{CP}$ [Open]\n(Ch. 12, 14)"),
            ("non-equilibrium", "none (crossover;\nSM-EWBG fails here)",
             "the expansion arrow,\npermanent [Theorem]\n(Ch. 6)")]
    ax.text(2.05, 4.1, "standard physics", fontsize=9, ha="center", color=GRAY)
    ax.text(5.45, 4.1, "this framework", fontsize=9, ha="center", color=BLUE)
    for i, (cond, std, fw) in enumerate(rows):
        y = 2.9 - 1.25 * i
        box(ax, 0.05, y, 1.25, 0.95, cond, fc="#fdf3e7", ec=ORANGE, fs=8)
        box(ax, 1.45, y, 2.5, 0.95, std, fc="#f4f4f4", ec=GRAY, fs=7.5)
        box(ax, 4.2, y, 2.5, 0.95, fw, fc="#eef4fb", fs=7.5)
    ax.set_xlim(0, 7.0), ax.set_ylim(0.3, 4.3)
    ax.set_title("Sakharov's checklist, and who pays each line")
    save(fig, "ch13_fig1_sakharov_map.png")

def ch13_fig2():
    fig, ax = blank((7.4, 3.8))
    box(ax, 0.3, 2.2, 3.0, 0.9, "Failure 1: equilibrium\ncrossover, no bubbles\n(lattice result)", fc="#fdecec", ec="#c0392b", fs=8.5)
    box(ax, 4.0, 2.2, 3.0, 0.9, "Failure 2: CP magnitude\nGIM-crushed CKM\n($\\eta_B \\sim 10^{-19}$)", fc="#fdecec", ec="#c0392b", fs=8.5)
    box(ax, 0.3, 0.5, 3.0, 0.9, "repaired structurally:\nexpansion arrow = permanent\nnon-equilibrium [Theorem]", fc="#eaf6ea", ec=GREEN, fs=8.5)
    box(ax, 4.0, 0.5, 3.0, 0.9, "inherited as [Open]:\nboundary-GIM evasion\n= the named calculation", fc="#fdf3e7", ec=ORANGE, fs=8.5)
    arrow(ax, 1.8, 2.15, 1.8, 1.45, color=GREEN, lw=1.5)
    arrow(ax, 5.5, 2.15, 5.5, 1.45, color=ORANGE, lw=1.5)
    ax.set_xlim(0, 7.4), ax.set_ylim(0.2, 3.5)
    ax.set_title("SM electroweak baryogenesis fails twice; the framework's honest ledger")
    save(fig, "ch13_fig2_ewbg_failure_modes.png")

# ---------------------------------------------------------------- ch14 fig1
def ch14_fig1():
    fig, ax = blank((8.6, 3.4))
    steps = [("microphysics:\nwall angles\n$P(\\theta_0, \\theta_L)$", "#eef4fb", BLUE),
             ("crossing window\n$\\mathcal{P}_\\times$ [Theorem]\n(geometric, $\\leq 1$)", "#eaf6ea", GREEN),
             ("sign bias\n$\\varepsilon_{CP}$ [Open]\n(the CP cost)", "#fdf3e7", ORANGE),
             ("yield/crossing\n$\\pm 1$ exactly\n[anomaly-protected]", "#eaf6ea", GREEN),
             ("sphalerons\n$28/79$\n[Standard]", "#f4f4f4", GRAY),
             ("dilution\n$45/2\\pi^2 g_*$\n[Standard]", "#f4f4f4", GRAY),
             ("$\\eta_B$", "#f3eefb", PURPLE)]
    x = 0.15
    for i, (t, fc, ec) in enumerate(steps):
        w = 1.05 if i < 6 else 0.7
        box(ax, x, 1.2, w, 1.0, t, fc=fc, ec=ec, fs=7.5)
        if i < 6:
            arrow(ax, x + w + 0.02, 1.7, x + w + 0.16, 1.7, color="k")
        x += w + 0.18
    ax.text(4.3, 0.6, r"bold-bordered boxes are protected (theorem/standard); the single open factor is $\varepsilon_{CP}$",
            fontsize=8.5, ha="center")
    ax.set_xlim(0, 8.6), ax.set_ylim(0.3, 2.8)
    ax.set_title(r"The assembly pipeline: $\eta_B = 2.3\times10^{-2}\,\mathcal{P}_\times\,\varepsilon_{CP}$")
    save(fig, "ch14_fig1_pipeline.png")

# ---------------------------------------------------------------- ch15 fig3
def ch15_fig3():
    fig, ax = blank((6.6, 4.0))
    ax.add_patch(Rectangle((0.6, 0.6), 5.2, 2.8, fill=False, lw=1.8, ec=BLUE))
    ax.add_patch(Rectangle((1.1, 1.0), 2.2, 1.6, fill=False, lw=1.5, ec=ORANGE))
    ax.add_patch(Rectangle((2.6, 1.5), 2.6, 1.5, fill=False, lw=1.5, ec=GREEN))
    ax.text(0.75, 3.15, "system A: which $L$?", fontsize=9, color=BLUE)
    ax.text(1.2, 2.4, "subsystem B:\nits own $L$?", fontsize=8.5, color=ORANGE)
    ax.text(3.9, 2.7, "overlapping C:\nboth?", fontsize=8.5, color=GREEN)
    ax.text(3.3, 0.25, "fixed-$L$ damping cannot serve nested/overlapping systems consistently;\n"
                       "in the sector arena $L$ is a quantum label and the question evaporates (Ch. 3)",
            fontsize=8.5, ha="center")
    ax.set_xlim(0.2, 6.4), ax.set_ylim(-0.4, 3.6)
    ax.set_title('The unanswerable question: "which $L$?"')
    save(fig, "ch15_fig3_which_L.png")

# ---------------------------------------------------------------- ch16 fig1
def ch16_fig1():
    fig, ax = blank((6.8, 4.2))
    x = np.linspace(0.5, 6.3, 100)
    for i, (y0, amp) in enumerate([(0.8, 0.06), (1.9, 0.10), (3.0, 0.05)]):
        ax.plot(x, y0 + amp * np.sin(2.2 * x), color=BLUE, lw=1.6)
        ax.text(6.35, y0, fr"$t_{i}$", fontsize=9)
    arrow(ax, 2.2, 0.95, 2.2, 1.85, color=ORANGE, lw=1.5)
    ax.text(2.32, 1.35, r"lapse $N$: proper time\nbetween slices", fontsize=8, color=ORANGE)
    arrow(ax, 4.0, 1.0, 4.7, 1.95, color=GREEN, lw=1.5)
    ax.text(4.6, 1.35, r"shift $N^i$: coordinates\nslide along slices", fontsize=8, color=GREEN)
    ax.text(1.0, 3.45, r"each slice carries the dynamical $g_{ij}(\mathbf{x})$", fontsize=9)
    ax.set_xlim(0.3, 7.2), ax.set_ylim(0.4, 3.8)
    ax.set_title("The ADM split: space-at-an-instant, stacked")
    save(fig, "ch16_fig1_adm_foliation.png")

# ---------------------------------------------------------------- ch18 figs
def ch18_fig2():
    fig, ax = blank((7.8, 4.2))
    box(ax, 0.3, 2.7, 3.2, 1.0, "diagonal / Cartan (3)\nconformal + 2 stretches\nPart I's $\\hat T^{(a)}_\\pm$: discrete ladder", fs=8.5)
    box(ax, 4.2, 2.7, 3.2, 1.0, "off-diagonal shears (3)\nedge angles\n$\\hat S^{(ab)}$: continuous, iso-energy", fc="#eaf6ea", ec=GREEN, fs=8.5)
    box(ax, 0.3, 1.0, 3.2, 1.0, "frame rotations (3)\n$SO(3)$: pure gauge\n(edge relabeling)", fc="#f4f4f4", ec=GRAY, fs=8.5)
    box(ax, 4.2, 1.0, 3.2, 1.0, "closure: $\\mathfrak{gl}(3)$\nsymmetric 6 = metric\nantisymmetric 3 = gauge", fc="#f3eefb", ec=PURPLE, fs=8.5)
    ax.text(3.85, 0.45, "light bending + gravitational waves live in the shear (Weyl) sector:\n"
                        "dilations alone $\\Rightarrow$ conformally flat $\\Rightarrow$ no GR (Obstruction Theorem)",
            fontsize=8.5, ha="center", color=ORANGE)
    ax.set_xlim(0, 7.8), ax.set_ylim(-0.2, 4.1)
    ax.set_title("The generator audit: 6 metric directions, 3 gauge rotations")
    save(fig, "ch18_fig2_gl3_decomposition.png")

def ch18_fig3():
    fig, ax = blank((7.2, 3.8))
    th = np.linspace(0, 2 * np.pi, 100)
    for r, x0 in ((0.55, 1.4), (0.8, 1.4), (1.05, 1.4)):
        ax.plot(x0 + r * np.cos(th), 1.9 + r * np.sin(th), color=BLUE, lw=1.0)
    ax.text(1.4, 0.45, "dilations only:\nconcentric (conformally flat)\nWeyl tensor $\\equiv 0$", fontsize=8.5, ha="center")
    for ecc, x0 in ((0.3, 5.3),):
        ax.plot(x0 + 1.0 * np.cos(th), 1.9 + 0.62 * np.sin(th), color=GREEN, lw=1.2)
        ax.plot(x0 + 0.78 * np.cos(th + 0.5), 1.9 + 0.78 * np.sin(th + 0.5) * 0.85, color=GREEN, lw=1.0)
        ax.plot(x0 + 0.5 * np.cos(th - 0.4), 1.9 + 0.62 * np.sin(th - 0.4), color=GREEN, lw=0.9)
    ax.text(5.3, 0.45, "with shears:\nshapes tilt and squeeze —\nWeyl curvature available", fontsize=8.5, ha="center", color=GREEN)
    ax.text(3.35, 3.3, "1919 eclipse / LIGO live here $\\rightarrow$", fontsize=8.5, color=ORANGE)
    ax.set_xlim(0, 7.2), ax.set_ylim(0.1, 3.7)
    ax.set_title("The conformal obstruction, drawn")
    save(fig, "ch18_fig3_conformal_obstruction.png")

# ---------------------------------------------------------------- ch19 fig1
def ch19_fig1():
    fig, ax = blank((6.8, 4.4))
    th = np.linspace(0, 2 * np.pi, 100)
    for r in (0.5, 1.0, 1.5):
        ax.plot(3.4 + r * np.cos(th), 2.2 + 0.62 * r * np.sin(th), color=GRAY, lw=0.7, ls=":")
    arrow(ax, 3.4, 2.2, 5.2, 2.2, color=BLUE, lw=1.6)
    ax.text(4.9, 2.4, "diagonal: energy flows,\nstates frozen ($\\chi = 0$)\n$\\Rightarrow$ must JUMP", fontsize=8, color=BLUE)
    arrow(ax, 3.4, 2.2, 2.4, 3.4, color=GREEN, lw=1.6)
    ax.text(1.3, 3.5, "shear: iso-energy,\nstates move ($\\chi > 0$)\n$\\Rightarrow$ smooth FLOW", fontsize=8, color=GREEN)
    ax.text(3.4, 0.6, "dotted: energy level sets in the 6-dim space of metrics", fontsize=8.5, ha="center", color=GRAY)
    ax.set_xlim(0.3, 6.6), ax.set_ylim(0.3, 4.2)
    ax.set_title("Shape space split: the clock direction vs the propagating directions")
    save(fig, "ch19_fig1_shape_space.png")

# ---------------------------------------------------------------- ch20 fig1
def ch20_fig1():
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    t = np.linspace(0.2, 4, 100)
    ax = axes[0]
    for v, c in ((0.4, BLUE), (0.9, ORANGE), (1.6, GREEN)):
        ax.plot(t, 1 + v * t, color=c, lw=1.4)
    ax.set_title("embedding coupling: straight in $L$\n(isotropizes; Milne class)", fontsize=9)
    ax.set_xlabel("$t$"), ax.set_ylabel("$L_a(t)$")
    ax = axes[1]
    for p, c in ((-2/7, BLUE), (3/7, ORANGE), (6/7, GREEN)):
        ax.plot(t, np.exp(p * t), color=c, lw=1.4, label=fr"$p = {p:+.2f}$")
    ax.set_yscale("log")
    ax.set_title("generator coupling: straight in $\\ln L$\n(Kasner class — one axis contracts)", fontsize=9)
    ax.set_xlabel("$t$"), ax.set_ylabel("$L_a(t)$ (log)")
    ax.legend(fontsize=7)
    fig.savefig(FIGS / "ch20_fig1_two_couplings.png"); plt.close(fig)
    print("  -> ch20_fig1_two_couplings.png")

# ---------------------------------------------------------------- ch21 fig1
def ch21_fig1():
    fig, ax = blank((6.6, 3.8))
    ax.add_patch(Rectangle((1.0, 1.9), 2.2, 1.0, fc="#eef4fb", ec=BLUE, lw=1.4))
    ax.text(2.1, 2.4, "matter\n$+\\rho V$", ha="center", fontsize=10)
    ax.add_patch(Rectangle((3.6, 1.9), 2.2, 1.0, fc="#fdf3e7", ec=ORANGE, lw=1.4))
    ax.text(4.7, 2.4, "geometry\n$-3\\kappa V H^2$", ha="center", fontsize=10)
    ax.text(3.35, 1.55, "+", fontsize=14, ha="center")
    ax.text(3.35, 0.9, "= 0   (the constraint = the Friedmann equation)", fontsize=10, ha="center")
    arrow(ax, 2.1, 1.85, 4.0, 1.0, color=GRAY); arrow(ax, 4.7, 1.85, 4.0, 1.0, color=GRAY)
    ax.text(3.35, 3.3, r"messenger: $g(\alpha) = \frac{1}{2}\sqrt{\rho(\alpha)/3\kappa}$ — the ladder coupling slaved to matter",
            fontsize=8.5, ha="center", color=GREEN)
    ax.set_xlim(0.5, 6.3), ax.set_ylim(0.5, 3.6)
    ax.set_title("The energy budget of a closed universe")
    save(fig, "ch21_fig1_constraint_budget.png")

# ---------------------------------------------------------------- ch22 figs
def ch22_fig1():
    fig, ax = blank((7.4, 3.6))
    box(ax, 0.3, 1.4, 2.0, 1.0, "matter spectra\non geometry $g$\n(zero-point sums)", fs=8.5)
    box(ax, 2.9, 1.4, 1.9, 1.0, "integrate out\n(one loop;\nheat kernel)", fc="#f4f4f4", ec=GRAY, fs=8.5)
    box(ax, 5.4, 1.4, 1.8, 1.0, "induced action\n$c_0\\Lambda^4 + c_1\\Lambda^2 R$\n= vacuum elasticity", fc="#eaf6ea", ec=GREEN, fs=8.5)
    arrow(ax, 2.35, 1.9, 2.85, 1.9); arrow(ax, 4.85, 1.9, 5.35, 1.9)
    ax.text(3.85, 0.8, r"no gravitational action postulated: $G^{-1} = N_{\rm eff}\Lambda^2_{\rm UV}/12\pi$ (Sakharov)",
            fontsize=9, ha="center")
    ax.set_xlim(0, 7.4), ax.set_ylim(0.4, 3.0)
    ax.set_title("Induced gravity: geometry's stiffness from matter loops")
    save(fig, "ch22_fig1_induced_logic.png")

def ch22_fig4():
    fig, ax = blank((6.4, 4.6))
    items = ["matter spectra\n(cell content)", "induced stiffness\n$G^{-1} = N_{\\rm eff}\\Lambda^2/12\\pi$",
             "Planck length\n$\\ell_P = \\sqrt{G}$", "cell scale\n$\\bar L = 1.6\\,c\\,\\ell_P$",
             "cutoff\n$\\Lambda = c/\\bar L$"]
    R0, cx, cy = 1.62, 3.2, 2.3
    for i, t in enumerate(items):
        ang = np.pi / 2 - 2 * np.pi * i / len(items)
        x, y = cx + R0 * np.cos(ang), cy + R0 * np.sin(ang)
        box(ax, x - 0.75, y - 0.3, 1.5, 0.6, t, fs=7.5)
        ang2 = np.pi / 2 - 2 * np.pi * (i + 0.5) / len(items)
        arrow(ax, cx + R0 * 1.12 * np.cos(ang2 + 0.45), cy + R0 * 1.12 * np.sin(ang2 + 0.45),
              cx + R0 * 1.12 * np.cos(ang2 - 0.1), cy + R0 * 1.12 * np.sin(ang2 - 0.1), color=GRAY, lw=1.0)
    ax.text(cx, cy, "self-consistent:\nthe loop closes", ha="center", fontsize=9, color=GREEN)
    ax.set_xlim(0.4, 6.2), ax.set_ylim(-0.1, 4.6)
    ax.set_title("The Planck loop")
    save(fig, "ch22_fig4_planck_loop.png")

# ---------------------------------------------------------------- ch23 figs
def ch23_fig1():
    fig, ax = blank((8.2, 3.4))
    steps = [("global Weyl\n(box spectrum,\nCh. 2/17) ", "#eef4fb", BLUE),
             ("local Weyl\n+ $W_\\mu$\n(gauge principle)", "#f4f4f4", GRAY),
             ("obstruction:\nEH not invariant\n(4-derivative dead end)", "#fdecec", "#c0392b"),
             ("compensator $\\sigma$\n$= c_\\sigma/L_{\\rm cell}$\n(identified!)", "#fdf3e7", ORANGE),
             ("gauge-fix $\\sigma_0$:\nEH + massive $W$\n$M_P^2 = \\sigma_0^2/6$", "#eaf6ea", GREEN)]
    x = 0.2
    for i, (t, fc, ec) in enumerate(steps):
        box(ax, x, 1.1, 1.4, 1.1, t, fc=fc, ec=ec, fs=7.5)
        if i < 4:
            arrow(ax, x + 1.43, 1.65, x + 1.57, 1.65, color="k")
        x += 1.6
    ax.text(4.1, 0.55, "framework-supplied inputs highlighted (orange/blue); the rest is standard machinery",
            fontsize=8.5, ha="center")
    ax.set_xlim(0, 8.2), ax.set_ylim(0.2, 2.8)
    ax.set_title("The gauging ladder, with its two framework inputs")
    save(fig, "ch23_fig1_gauging_ladder.png")

def ch23_fig2():
    fig, ax = blank((7.0, 4.0))
    box(ax, 2.5, 3.0, 2.0, 0.7, "spectrum homogeneity\n$(n, L) \\to (\\lambda n, \\lambda L)$", fs=8)
    box(ax, 0.4, 1.6, 2.4, 0.8, "Weyl gauging route:\n$\\sigma = c_\\sigma/L_{\\rm cell}$,\n$M_P^2 = \\sigma_0^2/6$", fs=8)
    box(ax, 4.2, 1.6, 2.4, 0.8, "Sakharov route:\n$M_P^2 \\sim N_{\\rm eff}/\\bar L^2$\n(matter loops)", fs=8)
    box(ax, 2.3, 0.3, 2.4, 0.7, "$M_P \\sim 1/\\bar L$:\ntwo routes, one value", fc="#eaf6ea", ec=GREEN, fs=8.5)
    arrow(ax, 3.1, 2.95, 1.8, 2.45, color=GRAY); arrow(ax, 3.9, 2.95, 5.2, 2.45, color=GRAY)
    arrow(ax, 1.7, 1.55, 3.0, 1.05, color=GRAY); arrow(ax, 5.3, 1.55, 4.0, 1.05, color=GRAY)
    ax.set_xlim(0, 7.0), ax.set_ylim(0, 4.0)
    ax.set_title("The compensator identification: overdetermined, and consistent")
    save(fig, "ch23_fig2_compensator_id.png")

# ---------------------------------------------------------------- ch24 fig2
def ch24_fig2():
    fig, ax = blank((6.6, 3.8))
    rng = np.random.default_rng(0)
    for i in range(9):
        for j in range(5):
            d = np.hypot(i - 4, j - 2)
            s = 0.42 + 0.16 * np.exp(-d**2 / 3.0)
            ax.add_patch(Rectangle((i * 0.7 + (0.7 - s * 0.7)/2, j * 0.7 + (0.7 - s * 0.7)/2),
                                   0.7 * s, 0.7 * s, fill=False, ec=BLUE, lw=1.0))
    ax.plot([4 * 0.7 + 0.35], [2 * 0.7 + 0.35], "o", ms=10, color=ORANGE)
    ax.text(4 * 0.7 + 0.35, -0.35, "mass here: cells dilated around it\n$\\delta\\alpha = -\\Phi_{\\rm Newton} > 0$",
            fontsize=8.5, ha="center")
    ax.set_xlim(-0.3, 6.6), ax.set_ylim(-0.8, 3.8)
    ax.set_title("The weak-field metric, cell-model edition")
    save(fig, "ch24_fig2_well_schematic.png")

# ---------------------------------------------------------------- ch25 figs
def ch25_fig1():
    fig, ax = blank((7.6, 3.6))
    box(ax, 0.3, 1.5, 2.0, 1.1, "conformal (1)\npredicted $+C(s)$\nsign VERIFIED,\nmagnitude trending", fc="#fdf3e7", ec=ORANGE, fs=7.5)
    box(ax, 2.8, 1.5, 2.0, 1.1, "TT $+$ and $\\times$ (2)\npredicted $-C(s)$\nVERIFIED 1–2%\n= the graviton", fc="#eaf6ea", ec=GREEN, fs=7.5)
    box(ax, 5.3, 1.5, 2.0, 1.1, "diffeo (3)\npredicted $0$\nNULLED in the\ncovariance window", fc="#f4f4f4", ec=GRAY, fs=7.5)
    ax.text(3.8, 0.85, r"$6 = 1_{\rm conf} + 2_{TT} + 3_{\rm diffeo}$: the DeWitt pattern, sector by sector",
            fontsize=9, ha="center")
    ax.set_xlim(0, 7.6), ax.set_ylim(0.4, 3.0)
    ax.set_title("The six metric directions and their measured verdicts (Ch. 25 table)")
    save(fig, "ch25_fig1_sector_split.png")

def ch25_fig2():
    # drawn from the measured table (s, |kappa/C|): diffeo falls, TT converges
    s = np.array([2.0, 2.5, 3.0])
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    ax.plot(s, [0.690, 0.893, 0.987], "o-", label=r"TT$+$  $|\kappa|/C$")
    ax.plot(s, [0.779, 0.940, 1.008], "s-", label=r"TT$\times$  $|\kappa|/C$")
    ax.plot(s, [3.5, 1.34, 0.453], "^--", color="#7f7f7f", label="diffeo control (spurious)")
    ax.axhline(1.0, color="k", ls=":", lw=0.9)
    ax.set_xlabel(r"regulator proper time $s$ (lattice units)")
    ax.set_ylabel(r"$|\kappa|/C(s)$")
    ax.set_yscale("log")
    ax.set_title("The covariance window: control collapses, signal converges to 1")
    ax.legend(fontsize=8)
    fig.savefig(FIGS / "ch25_fig2_covariance_window.png"); plt.close(fig)
    print("  -> ch25_fig2_covariance_window.png")

# ---------------------------------------------------------------- ch26 fig2
def ch26_fig2():
    fig, ax = blank((7.8, 4.6))
    rows = [("arrow of expansion", "Theorem + Computed", GREEN),
            ("master theorem $\\Delta Q = -\\Delta\\eta/2$", "Theorem", GREEN),
            ("bulk CP inert / polarization diagnosed", "Theorem + $10^{-16}$", GREEN),
            ("quantized pump $+1.0000$ / control $0$", "Computed", GREEN),
            ("$\\eta_B$ budget", "$\\varepsilon_{CP}$ Open", ORANGE),
            ("Kasner $\\Sigma_K = 0.9997$ / Friedmann sub-%", "Computed", GREEN),
            ("stiffness, Newton, BKL", "Computed/measured", GREEN),
            ("graviton $\\kappa_{TT} = -C(s)$ to 1–2%", "Computed, parameter-free", GREEN),
            ("$\\Lambda$, signature, lapse/shift", "Open, posed", ORANGE)]
    for i, (t, st, c) in enumerate(rows):
        y = 4.0 - 0.45 * i
        ax.text(0.2, y, t, fontsize=8.5, va="center")
        box(ax, 5.1, y - 0.16, 2.4, 0.33, st, fc=("#eaf6ea" if c == GREEN else "#fdf3e7"), ec=c, fs=7.5)
    ax.set_xlim(0, 7.8), ax.set_ylim(-0.3, 4.5)
    ax.set_title("The chain, graded (full table: Ch. 26.2)")
    save(fig, "ch26_fig2_chain_table.png")

# ---------------------------------------------------------------- ch27 fig1
def ch27_fig1():
    fig, ax = blank((8.2, 4.6))
    box(ax, 0.3, 3.3, 2.2, 0.8, "1. microscopic\nconstraint test", fc="#eaf6ea", ec=GREEN, fs=8)
    box(ax, 2.9, 3.3, 2.2, 0.8, "2. $N_{\\rm eff}$ and $c$\n(physical $G$)", fc="#eaf6ea", ec=GREEN, fs=8)
    box(ax, 5.5, 3.3, 2.2, 0.8, "3. 3+1D chiral bag\n(crossing density)", fc="#eaf6ea", ec=GREEN, fs=8)
    box(ax, 1.6, 1.9, 2.2, 0.8, "5/6. signature, lapse;\n$\\Pi(\\omega, q)$ wave speed", fc="#fdf3e7", ec=ORANGE, fs=8)
    box(ax, 4.4, 1.9, 2.2, 0.8, "4. boundary-GIM\n$\\Rightarrow \\varepsilon_{CP}$", fc="#fdf3e7", ec=ORANGE, fs=8)
    box(ax, 0.6, 0.5, 2.2, 0.7, "7–8. Bianchi IX walls;\nconformal cleanup", fc="#f4f4f4", ec=GRAY, fs=8)
    box(ax, 3.2, 0.5, 2.2, 0.7, "9. $\\Lambda$ (untouched,\nstated every time)", fc="#f4f4f4", ec=GRAY, fs=8)
    box(ax, 5.8, 0.5, 2.2, 0.7, "10. quantum-geometry\nphenomenology", fc="#f4f4f4", ec=GRAY, fs=8)
    arrow(ax, 1.4, 3.25, 2.4, 2.75, color=GRAY); arrow(ax, 4.0, 3.25, 2.9, 2.75, color=GRAY)
    arrow(ax, 6.6, 3.25, 5.7, 2.75, color=GRAY)
    ax.text(4.1, 4.35, "Tier I: finite calculations (months) — each hardens or breaks a load-bearing link",
            fontsize=8.5, ha="center", color=GREEN)
    ax.text(4.1, 1.55, "Tier II: the two decision points (gravity / baryogenesis)", fontsize=8.5, ha="center", color=ORANGE)
    ax.set_xlim(0, 8.2), ax.set_ylim(0.1, 4.6)
    ax.set_title("The research program as a dependency graph")
    save(fig, "ch27_fig1_program_map.png")

# ---------------------------------------------------------------- ch17 fig1
def ch17_fig1():
    fig, ax = blank((7.4, 3.6))
    ax.add_patch(Rectangle((0.6, 1.2), 1.6, 1.2, fill=False, lw=1.8, ec=BLUE))
    ax.add_patch(Rectangle((0.6, 1.2), 2.4, 1.2, fill=False, lw=1.2, ec=BLUE, ls="--"))
    ax.text(1.8, 0.8, "walls move:\n$L \\to L'$, metric trivial", fontsize=8.5, ha="center")
    ax.add_patch(Rectangle((4.6, 1.2), 1.6, 1.2, fill=False, lw=1.8, ec=GREEN))
    for xg in np.linspace(4.6, 6.2, 7):
        ax.plot([xg, xg], [1.2, 2.4], color=GREEN, lw=0.5, alpha=0.6)
    ax.text(5.4, 0.8, "walls fixed:\nmetric grows, $g_{\\xi\\xi} = L^2$", fontsize=8.5, ha="center")
    ax.text(3.7, 1.85, "$=$", fontsize=18, ha="center")
    ax.text(3.7, 2.9, "one system, two coordinate descriptions\n(a diffeomorphism — the Dictionary Theorem)",
            fontsize=9, ha="center")
    ax.set_xlim(0.2, 7.2), ax.set_ylim(0.3, 3.4)
    ax.set_title("The dictionary")
    save(fig, "ch17_fig1_dictionary.png")

if __name__ == "__main__":
    print("generating schematics:")
    ch01_fig1()
    three_faces("ch01_fig2_three_faces.png", labeled=False)
    three_faces("ch26_fig1_three_faces_closed.png", labeled=True)
    ch07_fig1(); ch07_fig2()
    ch08_fig2(); ch08_fig3()
    ch10_fig1()
    ch13_fig1(); ch13_fig2()
    ch14_fig1()
    ch15_fig3()
    ch16_fig1()
    ch17_fig1()
    ch18_fig2(); ch18_fig3()
    ch19_fig1()
    ch20_fig1()
    ch21_fig1()
    ch22_fig1(); ch22_fig4()
    ch23_fig1(); ch23_fig2()
    ch24_fig2()
    ch25_fig1(); ch25_fig2()
    ch26_fig2()
    ch27_fig1()
    print("done.")
