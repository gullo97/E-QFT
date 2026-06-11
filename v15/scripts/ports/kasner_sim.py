"""
Gravitational dynamics from the expanding-box model: Bianchi-I / Kasner milestone.

Per axis a, the geometry ladder |n_a> (with L_a = n_a L0/n0) is a tight-binding chain:

    H_a = sum_n t_a(n) ( |n+1><n| + |n><n+1| )

Two coupling laws are compared:

  'flat' : t(n) = lam * sqrt(n/(n+1))            (v13.3 Sec.9 hopping; -> ballistic in n)
  'gen'  : t(n) = g * (n + 1/2) * sqrt(n/(n+1))  (coupling through the DILATION GENERATOR;
                                                  |<n+1| x d/dx |n>| = 2n(n+1)/(2n+1) ~ n+1/2
                                                  -> ballistic in alpha = ln n)

Theory (derived in the report):
  semiclassical eqs  ndot = 2 t(n) sin(theta),  thetadot = -2 t'(n) cos(theta).
  flat: n(t) linear  -> orbits linear in L-space -> relational exponents isotropize
        (Bianchi I + matter behaviour: K = sum q^2 flows 9 -> 3).
  gen : at theta = +-pi/2 (fixed point), alphadot = -+2g exactly -> straight orbits in
        alpha-space = the Kasner/Jacobs orbit class; exponents p_a tunable via g_a.

Outputs: trajectories npz + figures fig4..fig6.
"""
import json
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})


# ----- check the generator-coupling matrix element  <m| x d/dx |n> = -2mn(-1)^{m+n}/(m^2-n^2)
def check_generator_element():
    L, ngrid = 1.0, 20001
    x = np.linspace(0, L, ngrid)
    out = []
    for n in (10, 40, 160):
        m = n + 1
        fn = np.sqrt(2/L)*np.sin(n*np.pi*x/L)
        fm = np.sqrt(2/L)*np.sin(m*np.pi*x/L)
        dfn = np.sqrt(2/L)*(n*np.pi/L)*np.cos(n*np.pi*x/L)
        val = np.trapezoid(fm * x * dfn, x)
        pred = 2*n*(n+1)/(2*n+1)
        out.append((n, val, pred))
        print(f"  <{m}|x d/dx|{n}> = {val:+.4f}   analytic 2n(n+1)/(2n+1) = {pred:.4f}")
    return out


def evolve_chain(N, tfun, n0, sigma, theta, times):
    """Evolve a Gaussian wavepacket on the chain; return <n>(t), Var(t), snapshots."""
    n = np.arange(1, N+1)
    t_n = tfun(n[:-1])                       # bond n -> n+1
    H = diags([t_n, t_n], [-1, 1], format="csc")
    psi = np.exp(-(n - n0)**2/(4*sigma**2) + 1j*theta*n)
    psi /= np.linalg.norm(psi)
    means, vars_, snaps = [], [], []
    tprev = 0.0
    for tt in times:
        if tt > tprev:
            psi = expm_multiply((-1j*(tt - tprev))*H, psi)
            tprev = tt
        P = np.abs(psi)**2
        mu = float(P @ n)
        means.append(mu)
        vars_.append(float(P @ n**2) - mu**2)
        snaps.append(P[::4].astype(np.float32))
    return np.array(means), np.array(vars_), np.array(snaps), n[::4]


def relational(times, means_list):
    """alpha_a, relational exponents p_a = (dalpha_a/dt)/(sum_b dalpha_b/dt), and Sum p^2."""
    al = np.array([np.log(m) for m in means_list])         # (3, T)
    dal = np.array([np.gradient(a, times) for a in al])
    tot = dal.sum(axis=0)
    p = dal / tot
    return al, p, (p**2).sum(axis=0)


print("# generator-coupling check (theory: lambda_n ~ n + 1/2)")
check_generator_element()

# ================= Variant (i): v13.3 hopping ('flat') -> isotropization =================
lam = 1.0
T, NT = 800.0, 161
times = np.linspace(0, T, NT)
n0s, sig = (30, 600, 2500), (8, 20, 40)
N = 4600
flat = lambda n: lam*np.sqrt(n/(n+1))
means_i, snaps_i = [], []
for n0, sg in zip(n0s, sig):
    mu, var, sn, ngrid_i = evolve_chain(N, flat, n0, sg, -np.pi/2, times)
    means_i.append(mu); snaps_i.append(sn)
    print(f"[flat] n0={n0}: n(T)={mu[-1]:.0f}  spread={np.sqrt(var[-1]):.0f}")
al_i, p_i, S2_i = relational(times, means_i)
print(f"[flat] Sum p^2: start {S2_i[2]:.3f} -> end {S2_i[-1]:.3f}  (vacuum Kasner=1, isotropic=1/3)")

# ================= Variant (ii): generator coupling ('gen') -> Kasner orbits ==============
# Kasner-circle target p = (-2/7, 3/7, 6/7):  Sum p = Sum p^2 = 1.
g0 = 0.30
p_t = np.array([-2/7, 3/7, 6/7])
Tg, NTg = 3.2, 161
times_g = np.linspace(0, Tg, NTg)
N_g = 5200
means_g, snaps_g = [], []
for a, p in enumerate(p_t):
    ga = g0*abs(p)
    th = -np.pi/2*np.sign(p)
    genf = (lambda gg: (lambda n: gg*(n + 0.5)*np.sqrt(n/(n+1))))(ga)
    n0 = 700 if p < 0 else 500
    mu, var, sn, ngrid_g = evolve_chain(N_g, genf, n0, 0.06*n0, th, times_g)
    means_g.append(mu); snaps_g.append(sn)
    rate = np.gradient(np.log(mu), times_g)
    print(f"[gen]  p_target={p:+.3f}: measured alphadot = {rate[NTg//2]:+.4f}"
          f"  predicted 2 g_a sgn = {2*ga*np.sign(p):+.4f}")
al_g, p_g, S2_g = relational(times_g, means_g)
print(f"[gen]  exponents at mid-time: {np.round(p_g[:, NTg//2], 4)}  target {np.round(p_t,4)}")
print(f"[gen]  Sum p^2 mid-time: {S2_g[NTg//2]:.4f}   (target 1.0000, vacuum Kasner)")

# off-circle run: isotropic stiff-matter point p=(1/3,1/3,1/3) -> Sum p^2 = 1/3
means_c = []
for a in range(3):
    genf = (lambda gg: (lambda n: gg*(n + 0.5)*np.sqrt(n/(n+1))))(g0/3)
    mu, var, sn, _ = evolve_chain(N_g, genf, 500, 30, -np.pi/2, times_g)
    means_c.append(mu)
al_c, p_c, S2_c = relational(times_g, means_c)
print(f"[gen-iso] Sum p^2 mid-time: {S2_c[NTg//2]:.4f}   (target 1/3 = {1/3:.4f})")

np.savez_compressed("kasner_traj.npz",
                    times=times, al_i=al_i, p_i=p_i, S2_i=S2_i,
                    snaps_i=np.array(snaps_i), ngrid_i=ngrid_i,
                    times_g=times_g, al_g=al_g, p_g=p_g, S2_g=S2_g,
                    snaps_g=np.array(snaps_g), ngrid_g=ngrid_g, p_target=p_t)

# ============================== figures ==============================
# Fig 4: relational exponents and Sum p^2 for both variants
fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.0))
for a in range(3):
    axs[0,0].plot(times, p_i[a], label=fr"$p_{a+1}$")
axs[0,0].axhline(1/3, color="k", ls=":"); axs[0,0].set_title("v13.3 hopping: exponents isotropize")
axs[0,0].set_xlabel("model time"); axs[0,0].set_ylabel(r"$p_a(t)$"); axs[0,0].legend()
axs[0,1].plot(times, S2_i, color="crimson")
axs[0,1].axhline(1, color="k", ls="--", label=r"vacuum Kasner $\Sigma p^2=1$")
axs[0,1].axhline(1/3, color="k", ls=":", label=r"isotropic $\Sigma p^2=1/3$")
axs[0,1].set_title(r"v13.3 hopping: $\Sigma p_a^2$ flows to isotropy")
axs[0,1].set_xlabel("model time"); axs[0,1].legend()
for a in range(3):
    axs[1,0].plot(times_g, p_g[a], label=fr"$p_{a+1}$")
    axs[1,0].axhline(p_t[a], color="gray", ls="--", lw=0.8)
axs[1,0].set_title(r"generator coupling: exponents lock to Kasner $(-\frac{2}{7},\frac{3}{7},\frac{6}{7})$")
axs[1,0].set_xlabel("model time"); axs[1,0].set_ylabel(r"$p_a(t)$"); axs[1,0].legend()
axs[1,1].plot(times_g, S2_g, color="seagreen", label="Kasner-circle run")
axs[1,1].plot(times_g, S2_c, color="steelblue", label="isotropic run")
axs[1,1].axhline(1, color="k", ls="--"); axs[1,1].axhline(1/3, color="k", ls=":")
axs[1,1].set_title(r"generator coupling: $\Sigma p_a^2$ is constant (Kasner class)")
axs[1,1].set_xlabel("model time"); axs[1,1].legend()
fig.tight_layout(); fig.savefig("figs/fig4_relational_exponents.png", dpi=150); plt.close(fig)

# Fig 5: the Kasner plane (p1, p2) with the vacuum circle Sum p^2 = 1, Sum p = 1
fig, ax = plt.subplots(figsize=(6.4, 6.0))
th = np.linspace(0, 2*np.pi, 400)
# parametrize the Kasner ellipse: p = (1/3,1/3) + r(cos,sin) with constraint -> solve numerically
P1, P2 = np.meshgrid(np.linspace(-0.6, 1.2, 400), np.linspace(-0.6, 1.2, 400))
F = P1**2 + P2**2 + (1 - P1 - P2)**2
ax.contour(P1, P2, F, levels=[1.0], colors="k", linestyles="--")
ax.plot([1/3], [1/3], "k*", ms=14, label=r"isotropic FRW $(\Sigma p^2=\frac{1}{3})$")
ax.plot(p_i[0], p_i[1], "-", color="crimson", lw=2, label="v13.3 hopping (flows to isotropy)")
ax.plot(p_i[0,2], p_i[1,2], "o", color="crimson")
ax.plot(p_g[0], p_g[1], ".", color="seagreen", ms=4, label="generator coupling (sits on circle)")
ax.plot([p_t[0]], [p_t[1]], "s", color="seagreen", ms=9, mfc="none", label=r"target $(-\frac{2}{7},\frac{3}{7})$")
ax.set_xlabel(r"$p_1$"); ax.set_ylabel(r"$p_2$")
ax.set_title("Kasner plane: dashed = vacuum Kasner circle")
ax.legend(fontsize=9); ax.set_aspect("equal")
fig.tight_layout(); fig.savefig("figs/fig5_kasner_plane.png", dpi=150); plt.close(fig)

# Fig 6: orbit geometry: L2 vs L1 (linear for flat), alpha2 vs alpha1 (linear for gen)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.4))
L1, L2 = np.exp(al_i[0]), np.exp(al_i[1])
a1.plot(L1, L2, "crimson"); c = np.polyfit(L1, L2, 1)
a1.plot(L1, np.polyval(c, L1), "k:", label=f"linear fit, residual {np.abs(L2-np.polyval(c,L1)).max()/L2.max():.1e}")
a1.set_xlabel(r"$L_1$"); a1.set_ylabel(r"$L_2$")
a1.set_title("v13.3 hopping: orbits are straight in $L$-space")
a1.legend(fontsize=9)
a2.plot(al_g[1], al_g[2], "seagreen"); c2 = np.polyfit(al_g[1], al_g[2], 1)
a2.plot(al_g[1], np.polyval(c2, al_g[1]), "k:",
        label=f"slope {c2[0]:.3f} vs Kasner $p_3/p_2$ = {p_t[2]/p_t[1]:.3f}")
a2.set_xlabel(r"$\alpha_2=\ln L_2$"); a2.set_ylabel(r"$\alpha_3=\ln L_3$")
a2.set_title(r"generator coupling: orbits straight in $\alpha$-space (Kasner class)")
a2.legend(fontsize=9)
fig.tight_layout(); fig.savefig("figs/fig6_orbits.png", dpi=150); plt.close(fig)
print("saved figs 4-6 and kasner_traj.npz")
