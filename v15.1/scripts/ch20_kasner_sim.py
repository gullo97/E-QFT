#!/usr/bin/env python3
"""
ch20_kasner_sim.py — minisuperspace dynamics on the geometry ladder (Ch. 20).

Per axis a, the geometry ladder |n_a> is a tight-binding chain
H_a = sum_n t_a(n)(|n+1><n| + h.c.).  Two coupling laws compared:

  'flat' : t(n) = lam sqrt(n/(n+1))              -> ballistic in n: isotropization
  'gen'  : t(n) = g (n + 1/2) sqrt(n/(n+1))      (dilation-generator coupling,
           |<n+1| x d/dx |n>| = 2n(n+1)/(2n+1) ~ n + 1/2) -> ballistic in
           alpha = ln n: straight orbits in alpha-space, the Kasner class.

Semiclassics: ndot = 2 t(n) sin(theta), thetadot = -2 t'(n) cos(theta); at the
fixed points theta = +-pi/2, alphadot = -+2g exactly.  Per-axis sign choice
sigma_a = -sign(theta_a) realizes p_a = sigma_a g_a / sum_b sigma_b g_b (the
sign vector of Thm 20.1(iii)); target here: Kasner point (-2/7, 3/7, 6/7).

Adapted from scripts/ports/kasner_sim.py (the original generating script).
Writes ../data/ch20_kasner_traj.npz; Figures 20.2-20.4.
"""
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ladder import evolve_chain, relational

HERE = pathlib.Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
DATA = HERE.parent / "data"
DATA.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


# ----- generator-coupling matrix element  <m| x d/dx |n> = -2mn(-1)^{m+n}/(m^2-n^2)
def check_generator_element():
    L, ngrid = 1.0, 20001
    x = np.linspace(0, L, ngrid)
    for n in (10, 40, 160):
        m = n + 1
        fm = np.sqrt(2 / L) * np.sin(m * np.pi * x / L)
        dfn = np.sqrt(2 / L) * (n * np.pi / L) * np.cos(n * np.pi * x / L)
        val = trapz(fm * x * dfn, x)
        pred = 2 * n * (n + 1) / (2 * n + 1)
        print(f"  <{m}|x d/dx|{n}> = {val:+.4f}   analytic 2n(n+1)/(2n+1) = {pred:.4f}")


print("# generator-coupling check (theory: lambda_n ~ n + 1/2)")
check_generator_element()

# ================= Variant (i): flat hopping -> isotropization =================
lam = 1.0
T, NT = 800.0, 161
times = np.linspace(0, T, NT)
n0s, sig = (30, 600, 2500), (8, 20, 40)
N = 4600
flat = lambda n: lam * np.sqrt(n / (n + 1))
means_i, snaps_i = [], []
for n0, sg in zip(n0s, sig):
    mu, var, sn, ngrid_i = evolve_chain(N, flat, n0, sg, -np.pi / 2, times)
    means_i.append(mu); snaps_i.append(sn)
    print(f"[flat] n0={n0}: n(T)={mu[-1]:.0f}  spread={np.sqrt(var[-1]):.0f}")
al_i, p_i, S2_i = relational(times, means_i)
print(f"[flat] Sum p^2: start {S2_i[2]:.3f} -> end {S2_i[-1]:.3f}  (vacuum Kasner=1, isotropic=1/3)")

# ============ Variant (ii): generator coupling -> Kasner orbits ============
# Kasner-circle target p = (-2/7, 3/7, 6/7):  Sum p = Sum p^2 = 1.
g0 = 0.30
p_t = np.array([-2 / 7, 3 / 7, 6 / 7])
Tg, NTg = 3.2, 161
times_g = np.linspace(0, Tg, NTg)
N_g = 5200
means_g, vars_g, snaps_g = [], [], []
for a, p in enumerate(p_t):
    ga = g0 * abs(p)
    th = -np.pi / 2 * np.sign(p)        # per-axis sign choice sigma_a (Thm 20.1 iii)
    genf = (lambda gg: (lambda n: gg * (n + 0.5) * np.sqrt(n / (n + 1))))(ga)
    n0 = 700 if p < 0 else 500
    mu, var, sn, ngrid_g = evolve_chain(N_g, genf, n0, 0.06 * n0, th, times_g)
    means_g.append(mu); vars_g.append(var); snaps_g.append(sn)
    rate = np.gradient(np.log(mu), times_g)
    print(f"[gen]  p_target={p:+.3f}: measured alphadot = {rate[NTg//2]:+.4f}"
          f"  predicted 2 g_a sgn = {2*ga*np.sign(p):+.4f}")
al_g, p_g, S2_g = relational(times_g, means_g)
print(f"[gen]  exponents at mid-time: {np.round(p_g[:, NTg//2], 4)}  target {np.round(p_t, 4)}")
print(f"[gen]  Sum p^2 mid-time: {S2_g[NTg//2]:.4f}   (target 1.0000, vacuum Kasner)")

# packet-width tracking (relative spread stays small: semiclassical validity)
for a, p in enumerate(p_t):
    mu_end, sd_end = means_g[a][-1], np.sqrt(vars_g[a][-1])
    print(f"[gen]  axis {a+1}: <n>(T) = {mu_end:.0f}, sigma_n = {sd_end:.0f}, "
          f"sigma_n/<n> = {sd_end/mu_end:.4f}")

# off-circle run: isotropic stiff-matter point p=(1/3,1/3,1/3) -> Sum p^2 = 1/3
means_c = []
for a in range(3):
    genf = (lambda gg: (lambda n: gg * (n + 0.5) * np.sqrt(n / (n + 1))))(g0 / 3)
    mu, var, sn, _ = evolve_chain(N_g, genf, 500, 30, -np.pi / 2, times_g)
    means_c.append(mu)
al_c, p_c, S2_c = relational(times_g, means_c)
print(f"[gen-iso] Sum p^2 mid-time: {S2_c[NTg//2]:.4f}   (target 1/3 = {1/3:.4f})")

np.savez_compressed(DATA / "ch20_kasner_traj.npz",
                    times=times, al_i=al_i, p_i=p_i, S2_i=S2_i,
                    snaps_i=np.array(snaps_i), ngrid_i=ngrid_i,
                    times_g=times_g, al_g=al_g, p_g=p_g, S2_g=S2_g,
                    snaps_g=np.array(snaps_g), ngrid_g=ngrid_g, p_target=p_t)

# ============================== figures ==============================
fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.0))
for a in range(3):
    axs[0, 0].plot(times, p_i[a], label=fr"$p_{a+1}$")
axs[0, 0].axhline(1 / 3, color="k", ls=":"); axs[0, 0].set_title("flat hopping: exponents isotropize")
axs[0, 0].set_xlabel("model time"); axs[0, 0].set_ylabel(r"$p_a(t)$"); axs[0, 0].legend()
axs[0, 1].plot(times, S2_i, color="crimson")
axs[0, 1].axhline(1, color="k", ls="--", label=r"vacuum Kasner $\Sigma p^2=1$")
axs[0, 1].axhline(1 / 3, color="k", ls=":", label=r"isotropic $\Sigma p^2=1/3$")
axs[0, 1].set_title(r"flat hopping: $\Sigma p_a^2$ flows to isotropy")
axs[0, 1].set_xlabel("model time"); axs[0, 1].legend()
for a in range(3):
    axs[1, 0].plot(times_g, p_g[a], label=fr"$p_{a+1}$")
    axs[1, 0].axhline(p_t[a], color="gray", ls="--", lw=0.8)
axs[1, 0].set_title(r"generator coupling: exponents lock to Kasner $(-\frac{2}{7},\frac{3}{7},\frac{6}{7})$")
axs[1, 0].set_xlabel("model time"); axs[1, 0].set_ylabel(r"$p_a(t)$"); axs[1, 0].legend()
axs[1, 1].plot(times_g, S2_g, color="seagreen", label="Kasner-circle run")
axs[1, 1].plot(times_g, S2_c, color="steelblue", label="isotropic run")
axs[1, 1].axhline(1, color="k", ls="--"); axs[1, 1].axhline(1 / 3, color="k", ls=":")
axs[1, 1].set_title(r"generator coupling: $\Sigma p_a^2$ is constant (Kasner class)")
axs[1, 1].set_xlabel("model time"); axs[1, 1].legend()
fig.tight_layout(); fig.savefig(FIGS / "ch20_fig2_relational_exponents.png", dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.4, 6.0))
P1, P2 = np.meshgrid(np.linspace(-0.6, 1.2, 400), np.linspace(-0.6, 1.2, 400))
F = P1**2 + P2**2 + (1 - P1 - P2)**2
ax.contour(P1, P2, F, levels=[1.0], colors="k", linestyles="--")
ax.plot([1 / 3], [1 / 3], "k*", ms=14, label=r"isotropic FRW $(\Sigma p^2=\frac{1}{3})$")
ax.plot(p_i[0], p_i[1], "-", color="crimson", lw=2, label="flat hopping (flows to isotropy)")
ax.plot(p_i[0, 2], p_i[1, 2], "o", color="crimson")
ax.plot(p_g[0], p_g[1], ".", color="seagreen", ms=4, label="generator coupling (sits on circle)")
ax.plot([p_t[0]], [p_t[1]], "s", color="seagreen", ms=9, mfc="none", label=r"target $(-\frac{2}{7},\frac{3}{7})$")
ax.set_xlabel(r"$p_1$"); ax.set_ylabel(r"$p_2$")
ax.set_title("Kasner plane: dashed = vacuum Kasner circle")
ax.legend(fontsize=9); ax.set_aspect("equal")
fig.tight_layout(); fig.savefig(FIGS / "ch20_fig3_kasner_plane.png", dpi=150); plt.close(fig)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
L1, L2 = np.exp(al_i[0]), np.exp(al_i[1])
a1.plot(L1, L2, "crimson"); c = np.polyfit(L1, L2, 1)
a1.plot(L1, np.polyval(c, L1), "k:",
        label=f"linear fit, residual {np.abs(L2-np.polyval(c,L1)).max()/L2.max():.1e}")
a1.set_xlabel(r"$L_1$"); a1.set_ylabel(r"$L_2$")
a1.set_title("flat hopping: orbits are straight in $L$-space")
a1.legend(fontsize=9)
a2.plot(al_g[1], al_g[2], "seagreen"); c2 = np.polyfit(al_g[1], al_g[2], 1)
a2.plot(al_g[1], np.polyval(c2, al_g[1]), "k:",
        label=f"slope {c2[0]:.3f} vs Kasner $p_3/p_2$ = {p_t[2]/p_t[1]:.3f}")
a2.set_xlabel(r"$\alpha_2=\ln L_2$"); a2.set_ylabel(r"$\alpha_3=\ln L_3$")
a2.set_title(r"generator coupling: orbits straight in $\alpha$-space (Kasner class)")
a2.legend(fontsize=9)
fig.tight_layout(); fig.savefig(FIGS / "ch20_fig4_orbit_geometry.png", dpi=150); plt.close(fig)
print("figures -> ch20_fig2/fig3/fig4; data -> ch20_kasner_traj.npz")
