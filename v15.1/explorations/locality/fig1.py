import numpy as np
from vstyle import *
import sim1_chain_maps as S
M, occ = S.M, S.occ
wB, wC = np.load("figs/sim1_w0.npy")
wA = np.array([np.mean([n/(n+1) for n in m]) if m else np.nan for m in S.modes0])
edge = np.array([i for i in list(range(0, 6)) + list(range(34, 40)) if occ[i] > 0])
bulk = np.array([i for i in range(12, 28) if occ[i] > 0])
runs = {}
for r in "AB":
    acc = []
    for seed in range(6):
        S.rng = np.random.default_rng(100+seed)
        H, W = S.run(r, proposals=12000, snaps=60); acc.append(np.log(H/H[0]))
        if seed == 0: runs[r+"H"] = H
    runs[r] = np.array(acc)
tt = np.linspace(0, 12000/occ.sum(), runs["A"].shape[1])
for r in "AB":
    e = runs[r][:, :, edge].mean(axis=(0, 2)); b = runs[r][:, :, bulk].mean(axis=(0, 2))
    k = np.argmin(abs(tt-10)); print(f"rule {r}: at 10 proposals/particle  <lnL> edge {e[k]:.3f} bulk {b[k]:.3f}; at end edge {e[-1]:.3f} bulk {b[-1]:.3f}")
HC = np.load("figs/sim1.npz")["HC"]

fig = plt.figure(figsize=(11, 7.6))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 0.32, 1.1], hspace=0.55, wspace=0.34)
ax = fig.add_subplot(gs[0, :2]); x = np.arange(M)
ax.semilogy(x, wA, "o-", color=C1, ms=4, lw=1.4, label="A  comoving labels (neighbours relabelled)")
ax.semilogy(x, wC, "^-", color=C3, ms=4, lw=1.4, label="C  physical-space embedding, closed ring")
ax.semilogy(x, wB, "s-", color=C2, ms=4, lw=1.4, label="B  physical-space embedding, open chain")
ax.set_ylabel("mean step weight at $t=0$"); ax.set_xlim(-0.8, M-0.2)
ax.set_title("(a) Short-time weight of an iso-energy step in cell $i$ (Ch. 6 rule, applied locally)")
ax.legend(loc="lower left", fontsize=8); ax.tick_params(labelbottom=False)
axo = fig.add_subplot(gs[1, :2], sharex=ax)
axo.bar(x, occ, color="#b9b8b2", width=0.8); axo.set_ylabel("occupants"); axo.set_xlabel("cell index $i$")
axo.grid(axis="x", visible=False)

ax = fig.add_subplot(gs[0:2, 2])
ax.scatter(occ, HC[-1], s=24, color=C3, edgecolor="white", lw=0.6)
ax.set_yscale("log"); ax.set_xlabel("occupants in cell"); ax.set_ylabel("final cell size $L_i$")
ax.set_title(r"(b) Rule C after $2.5\times10^{5}$ proposals:" "\n" "full cells grow, sparse cells are crushed")
ax.axhline(1, color=MUTED, lw=0.8, ls="--"); ax.text(10, 1.3, "initial size", color=MUTED, fontsize=8, ha="right")

vmax = max(runs["AH"].max(), runs["BH"].max()); vmax = np.log(vmax)
for k, (H, lab, tag) in enumerate([(runs["AH"], "A  comoving labels", "c"), (runs["BH"], "B  physical-space, open", "d")]):
    ax = fig.add_subplot(gs[2, k])
    im = ax.imshow(np.log(H.T/H[0][:, None]), aspect="auto", origin="lower", cmap=SEQ,
                   extent=[0, tt[-1], -0.5, M-0.5], vmin=0, vmax=vmax)
    ax.set_xlabel("proposals per particle"); ax.set_ylabel("cell index $i$"); ax.grid(False)
    ax.set_title(f"({tag}) rule {lab}")
cb = fig.colorbar(im, ax=fig.axes[-2:], fraction=0.03, pad=0.02); cb.set_label("$\ln L_i/L_i(0)$")
ax = fig.add_subplot(gs[2, 2])
for r, col in (("A", C1), ("B", C2)):
    ax.plot(tt, runs[r][:, :, edge].mean(axis=(0, 2)), color=col, lw=2, label=f"rule {r}, edge cells")
    ax.plot(tt, runs[r][:, :, bulk].mean(axis=(0, 2)), color=col, lw=2, ls="--", label=f"rule {r}, bulk cells")
ax.set_xlabel("proposals per particle"); ax.set_ylabel("mean $\ln L_i/L_i(0)$")
ax.set_title("(e) Edge vs bulk growth\n(6 seeds averaged)"); ax.legend(fontsize=7.5)
fig.savefig("figs/fig1_chain_maps.png"); print("saved")
