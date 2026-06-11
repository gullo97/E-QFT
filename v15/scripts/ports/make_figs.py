import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from chiral_exact import chiral_spectrum, crossing_L

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

# ---- Fig 1: convergence of the two DQ forms vs polarisation prediction
d = json.load(open("conv_results.json"))
rows = np.array(d["rows"]); pred = d["pred"]
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(rows[:,0], rows[:,1], "o-", label=r"$\Delta Q$ cross-branch form")
ax.plot(rows[:,0], rows[:,2], "s-", label=r"$\Delta Q$ same-branch form (paper 3)")
ax.axhline(pred, color="k", ls="--", lw=1.2,
           label=rf"polarisation $\frac{{1}}{{2}}\,\mathrm{{Tr}}_{{[0,L]}}\,\mathrm{{sgn}}\,H' = {pred:+.4f}$")
ax.axhline(0, color="crimson", ls=":", lw=1.6, label="net charge (exact theorem) = 0")
ax.set_xlabel(r"basis truncation $N$ (modes per branch)")
ax.set_ylabel(r"$\Delta Q$")
ax.set_title(r"CP-mass bag quench $L\to1.2L$, $\delta=\pi/4$, $mL=1$ (paper 3 Fig. 4 point)")
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig("figs/fig1_convergence.png", dpi=160); plt.close(fig)

# ---- Fig 2: fractional vacuum charge law
D = np.array([0.0,0.4,0.8,1.2,1.6,2.0,2.4,2.8])
Q = np.array([0.0,-0.07253,-0.13489,-0.19725,-0.25962,-0.32199,-0.38436,-0.44673])
fig, ax = plt.subplots(figsize=(6.0, 4.0))
ax.plot(D, Q, "o", ms=7, label=r"numerical $Q_{\rm vac}$ ($mL=6$, heat-kernel $\eta$)")
Dd = np.linspace(0, 3, 50)
ax.plot(Dd, -Dd/(2*np.pi), "k--", lw=1.2, label=r"$-\Delta\theta/2\pi$ (chiral-bag law)")
ax.set_xlabel(r"wall mismatch $\Delta\theta=\theta_L-\theta_0$")
ax.set_ylabel(r"$Q_{\rm vac}=-\eta/2$")
ax.set_title("Fractional vacuum charge of the 1+1D chiral bag")
ax.legend(); fig.tight_layout()
fig.savefig("figs/fig2_qvac_law.png", dpi=160); plt.close(fig)

# ---- Fig 3: level crossing and quantized Qvac jump
t0, tL, m = -2.0, 2.0, 1.0
Lstar = crossing_L(m, t0, tL)
Ls = np.linspace(0.2, 0.95, 31)
emin, qv = [], []
for L in Ls:
    E = chiral_spectrum(m, L, t0, tL, Emax=300.0)
    emin.append(E[np.argmin(np.abs(E))])
    ts = (0.06, 0.04, 0.025, 0.015)
    v = [np.sum(np.sign(E)*np.exp(-t*np.abs(E))) for t in ts]
    A = np.vstack([np.ones(3), ts[-3:], np.square(ts[-3:])]).T
    qv.append(-0.5*np.linalg.lstsq(A, v[-3:], rcond=None)[0][0])
fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.2, 6.4), sharex=True)
a1.plot(Ls, emin, "o-", ms=4); a1.axhline(0, color="k", lw=0.8)
a1.axvline(Lstar, color="crimson", ls="--", label=rf"$L^*={Lstar:.4f}$ (analytic)")
a1.set_ylabel(r"level nearest zero  $E_0(L)$"); a1.legend()
a1.set_title(r"Chiral bag $\theta=(-2,+2)$: level crossing pumps one unit of charge")
a2.plot(Ls, qv, "s-", ms=4, color="seagreen"); a2.axvline(Lstar, color="crimson", ls="--")
a2.set_ylabel(r"$Q_{\rm vac}(L)$"); a2.set_xlabel(r"box size $L$ ($m=1$)")
fig.tight_layout(); fig.savefig("figs/fig3_crossing.png", dpi=160); plt.close(fig)
print("figures written:", end=" ")
import os; print(sorted(os.listdir("figs")))
