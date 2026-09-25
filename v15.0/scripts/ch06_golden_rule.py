#!/usr/bin/env python3
"""
ch06_golden_rule.py — the Golden-Rule Embedding Theorem, numerically (Ch. 6, Fig. 6.3).

Two-sector model: sector A = box L (modes 1..N), sector B = contracted box
L'' = (n1-1)L/n1 (modes 1..N). Coupling V = g * (truncation overlap matrix).
Checks: (i) Parseval sums -> N_con (residual ~ 1/N from sinc tails);
        (ii) exact unitary evolution P(t)/(gt)^2 -> embedding norm as t -> 0.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIGS = HERE.parent / "thesis" / "figures"
plt.style.use(HERE / "style.mplstyle")

def overlap_matrix(n1, N, L=1.0):
    """O[k, n] = <phi_k(L'') | phi_n(L)> restricted to [0, L'']  (analytic)."""
    Lpp = (n1 - 1) / n1 * L
    k = np.arange(1, N + 1)[:, None] * np.pi / Lpp
    n = np.arange(1, N + 1)[None, :] * np.pi / L
    pref = 2 / np.sqrt(L * Lpp)
    with np.errstate(divide="ignore", invalid="ignore"):
        O = pref * (np.sin((k - n) * Lpp) / (2 * (k - n)) - np.sin((k + n) * Lpp) / (2 * (k + n)))
    res = np.where(np.abs(k - n) < 1e-12, pref * (Lpp / 2 - np.sin(2 * k * Lpp) / (4 * k)), O)
    return res

def N_con(n1, n2):
    return (n1 - 1) / n1 + np.sin(2 * np.pi * n2 / n1) / (2 * np.pi * n2)

def main():
    print("=== (i) Parseval sums vs closed-form N_con (220 modes/sector) ===")
    N = 220
    for n1, n2 in [(5, 3), (8, 5), (12, 7)]:
        O = overlap_matrix(n1, N)
        parseval = np.sum(O[:, n2 - 1]**2)
        print(f"  (n1,n2)=({n1},{n2}): Parseval {parseval:.5f}   N_con {N_con(n1, n2):.5f}"
              f"   residual {N_con(n1,n2)-parseval:+.5f}")

    print("=== (ii) exact evolution P(t)/(gt)^2 -> embedding norm ===")
    n1, n2, N, g = 5, 3, 160, 0.05
    O = overlap_matrix(n1, N)
    L, Lpp = 1.0, (n1 - 1) / n1
    EA = (np.arange(1, N + 1) * np.pi / L)**2 / 2
    EB = (np.arange(1, N + 1) * np.pi / Lpp)**2 / 2
    H = np.zeros((2 * N, 2 * N))
    H[:N, :N] = np.diag(EA)
    H[N:, N:] = np.diag(EB)
    H[N:, :N] = g * O
    H[:N, N:] = g * O.T
    w, U = np.linalg.eigh(H)
    psi0 = np.zeros(2 * N); psi0[n2 - 1] = 1.0
    c0 = U.T @ psi0
    ts, Ps = np.geomspace(2e-4, 2e-2, 12), []
    for t in ts:
        psi = U @ (np.exp(-1j * w * t) * c0)
        Ps.append(np.sum(np.abs(psi[N:])**2) / (g * t)**2)
    target = np.sum(O[:, n2 - 1]**2)
    print(f"  P(t)/(gt)^2 at t = {ts[0]:.1e} .. {ts[-1]:.1e}: {Ps[0]:.4f} .. {Ps[-1]:.4f}")
    print(f"  embedding norm (Parseval, this basis): {target:.4f}; N_con = {N_con(n1, n2):.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6))
    ax = axes[0]
    Ns = np.array([20, 40, 80, 160, 320, 640])
    for (a, b), color in zip([(5, 3), (8, 5), (12, 7)], ["#1f77b4", "#d95f00", "#2ca02c"]):
        vals = [np.sum(overlap_matrix(a, n)[:, b - 1]**2) for n in Ns]
        ax.semilogx(Ns, vals, "o-", color=color, ms=3.5, label=fr"$(n_1, n_2) = ({a}, {b})$")
        ax.axhline(N_con(a, b), color=color, ls=":", lw=0.8)
    ax.set_xlabel("basis size per sector"), ax.set_ylabel("Parseval sum")
    ax.set_title(r"postulate $\to$ theorem: sums $\to \mathcal{N}_{\mathrm{con}}$")
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.semilogx(ts, Ps, "o-", ms=3.5)
    ax.axhline(target, color="k", ls=":", lw=0.9)
    ax.set_xlabel("$t$"), ax.set_ylabel(r"$P(t)/(gt)^2$")
    ax.set_title("exact evolution $\\to$ embedding norm as $t \\to 0$")
    fig.savefig(FIGS / "ch06_fig3_golden_rule_convergence.png")
    print("figure -> ch06_fig3_golden_rule_convergence.png")

if __name__ == "__main__":
    main()
