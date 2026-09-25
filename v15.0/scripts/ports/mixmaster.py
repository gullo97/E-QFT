"""
Milestone 3, Sim F: Bianchi IX (Mixmaster) -- Kasner epochs bouncing off curvature walls.

Vacuum Bianchi IX in BKL time tau (dt = abc dtau), with A_i = ln a_i:

   2 A1'' = (a2^2 - a3^2)^2 - a1^4    (+ cyclic),
   constraint  C = A1'A2' + A2'A3' + A3'A1'
                 - (1/4)[a1^4+a2^4+a3^4 - 2a1^2a2^2 - 2a2^2a3^2 - 2a3^2a1^2] = 0.

Between wall activations the dynamics is exactly the milestone-1 Kasner flight
(straight lines in alpha-space, the generator-coupled ladder's orbit class); each wall
bounce maps Kasner parameter u -> u-1 (u>=2) or 1/(u-1) (1<u<2): the BKL map, the
engine of the chaotic approach to the singularity.  We integrate toward the singularity,
extract u per epoch from the relational exponents, and test the map.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

def rhs(tau, y):
    A, P = y[:3], y[3:]
    a2 = np.exp(np.clip(2*A, -700, 60))
    w = np.array([(a2[1]-a2[2])**2 - a2[0]**2,
                  (a2[2]-a2[0])**2 - a2[1]**2,
                  (a2[0]-a2[1])**2 - a2[2]**2])
    return np.concatenate([P, 0.5*w])

def constraint(y):
    A, P = y[:3], y[3:]
    a2 = np.exp(2*A)
    W = 0.25*(a2[0]**2 + a2[1]**2 + a2[2]**2 - 2*(a2[0]*a2[1] + a2[1]*a2[2] + a2[2]*a2[0]))
    return P[0]*P[1] + P[1]*P[2] + P[2]*P[0] - W

def kasner_p(u):
    D = 1 + u + u**2
    return np.array([-u/D, (1+u)/D, u*(1+u)/D])

# ---- initial data: deep Kasner regime, contracting toward the singularity
u0 = 4.3
p = kasner_p(u0)
A0 = np.array([-4.0, -7.0, -16.0])           # hierarchical: all wall terms tiny
Pk = -p                                       # volume decreasing (sum P < 0)
y0 = np.concatenate([A0, Pk])
lam = brentq(lambda L: constraint(np.concatenate([A0, Pk + L])), -0.5, 0.5)
y0[3:] += lam
print(f"# ICs: u0={u0}, constraint residual = {constraint(y0):.2e}, lambda shift = {lam:.2e}")

T = 1400.0
sol = solve_ivp(rhs, [0, T], y0, method="DOP853", rtol=1e-11, atol=1e-12,
                dense_output=False, t_eval=np.linspace(0, T, 30000))
A, P = sol.y[:3], sol.y[3:]
Cdrift = np.max(np.abs([constraint(sol.y[:, k]) for k in range(0, 30000, 300)]))
print(f"# constraint drift over run: {Cdrift:.2e}")

# relational exponents and u(tau)
ptau = P/np.sum(P, axis=0)
psort = np.sort(ptau, axis=0)
u_tau = psort[2]/psort[1]                     # p_max / p_mid

# ---- epoch detection: plateaus of u(tau)
du = np.abs(np.gradient(u_tau, sol.t))
flat = du < 2e-3
epochs, k = [], 0
while k < len(flat):
    if flat[k]:
        j = k
        while j < len(flat) and flat[j]:
            j += 1
        if sol.t[j-1] - sol.t[k] > 6.0:
            seg = slice(k + (j-k)//4, j - (j-k)//4)
            epochs.append((sol.t[k], sol.t[j-1], float(np.median(u_tau[seg]))))
        k = j
    else:
        k += 1

def bkl(u):
    return u - 1 if u >= 2 else 1/(u - 1)

print("# epochs and the BKL map:")
print(f"  epoch 1: u = {epochs[0][2]:.4f}   (initial condition u0 = {u0})")
upred = epochs[0][2]
for i in range(1, len(epochs)):
    upred = bkl(upred)
    print(f"  epoch {i+1}: u_measured = {epochs[i][2]:.4f}   BKL prediction = {upred:.4f}"
          f"   diff = {epochs[i][2]-upred:+.4f}")

np.savez_compressed("mixmaster_traj.npz", t=sol.t, A=A, P=P, ptau=ptau, u=u_tau,
                    epochs=np.array([(a, b, c) for a, b, c in epochs]))

# ---- figure
fig = plt.figure(figsize=(11, 7.2))
ax1 = fig.add_subplot(2, 2, (1, 2))
for i in range(3):
    ax1.plot(sol.t, A[i], label=fr"$\alpha_{i+1}$")
ax1.set_xlabel(r"BKL time $\tau$"); ax1.set_ylabel(r"$\alpha_i = \ln a_i$")
ax1.set_title("Mixmaster: Kasner epochs (straight flights) and curvature bounces")
for (t0, t1, uu) in epochs:
    ax1.axvspan(t0, t1, color="0.92", zorder=0)
    ax1.text((t0+t1)/2, ax1.get_ylim()[1]*0.0 + A.max()*0.9, f"u={uu:.2f}",
             ha="center", fontsize=8, color="0.3")
ax1.legend(fontsize=9, loc="lower left")
ax2 = fig.add_subplot(2, 2, 3)
ax2.plot(sol.t, u_tau, color="crimson", lw=1)
ax2.set_ylim(0.8, 5); ax2.set_xlabel(r"$\tau$"); ax2.set_ylabel("Kasner parameter $u$")
ax2.set_title("the BKL map in action: $u \\to u-1 \\to \\dots \\to 1/(u-1)$")
ax3 = fig.add_subplot(2, 2, 4)
th = np.linspace(0, 2*np.pi, 400)
P1, P2 = np.meshgrid(np.linspace(-0.6, 1.2, 300), np.linspace(-0.6, 1.2, 300))
ax3.contour(P1, P2, P1**2 + P2**2 + (1-P1-P2)**2, levels=[1.0], colors="k", linestyles="--")
ax3.plot(ptau[0], ptau[1], ".", color="crimson", ms=2)
ax3.plot([1/3], [1/3], "k*")
ax3.set_xlabel("$p_1$"); ax3.set_ylabel("$p_2$"); ax3.set_aspect("equal")
ax3.set_title("exponents hop around the Kasner circle")
fig.tight_layout(); fig.savefig("figs/fig13_mixmaster.png", dpi=150); plt.close(fig)
print("saved fig13 and mixmaster_traj.npz")
