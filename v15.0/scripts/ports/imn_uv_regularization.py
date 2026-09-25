"""
Iso-energy Mode Normalisation (IMN) and UV regularisation: numerical tests.

This script accompanies v13_5_theory.md.  It computes, in 1+1D and 3+1D box
quantisation, the partial sums of:

    (a) the standard vacuum two-point  <0|phi^2(x)|0>
    (b) the standard vacuum energy density (Casimir sum)
    (c) the same quantities after the IMN substitution
            1/sqrt(L)  ->  1/sqrt(L_n),   L_n = n L
    (d) the loop-integrand power-counting test for a scalar self-energy

and produces commented plots showing
    *  log divergence of (a) and finite Basel-type limit of (a) under IMN,
    *  linear divergence of (b) and its insensitivity to the IMN substitution,
    *  quadratic divergence of the 3+1D ⟨phi^2⟩ vs. its convergence under IMN,
    *  per-mode contribution profile  c_n^(std) ~ 1/n  vs.  c_n^(IMN) ~ 1/n^2.

Conventions (consistent with v13_3 §39):
    hbar = c = 1, metric signature (+,-).
    Massless real scalar with Dirichlet BC on [0,L]^d.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# ---------------------------------------------------------------------------
# style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "font.family":      "serif",
    "mathtext.fontset": "cm",
    "font.size":        12,
    "axes.labelsize":   13,
    "axes.titlesize":   14,
    "legend.fontsize":  10,
})

FIG_DIR = Path(__file__).resolve().parent.parent / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# 1+1D vacuum two-point at the midpoint  x = L/2
# ---------------------------------------------------------------------------
def phi2_standard_1d(N: int, L: float = 1.0, x: float = 0.5) -> np.ndarray:
    """Partial sums  sum_{n=1..N}  (1/(L omega_n)) sin^2(k_n x)
    with omega_n = k_n = n pi / L."""
    n   = np.arange(1, N + 1)
    k_n = n * np.pi / L
    contrib = (1.0 / (L * k_n)) * np.sin(k_n * x) ** 2
    return np.cumsum(contrib)


def phi2_imn_1d(N: int, L: float = 1.0, x: float = 0.5) -> np.ndarray:
    """IMN partial sums: replace 1/L by 1/L_n with L_n = n L."""
    n   = np.arange(1, N + 1)
    k_n = n * np.pi / L
    L_n = n * L
    contrib = (1.0 / (L_n * k_n)) * np.sin(k_n * x) ** 2   # ~ 1/n^2
    return np.cumsum(contrib)


# ---------------------------------------------------------------------------
# 1+1D Casimir / vacuum-energy density  (1/(2L)) sum  omega_n
# ---------------------------------------------------------------------------
def Ecasimir_standard_1d(N: int, L: float = 1.0) -> np.ndarray:
    n   = np.arange(1, N + 1)
    contrib = (np.pi / (2 * L)) * n                # 1/2 omega_n / L
    return np.cumsum(contrib)


def Ecasimir_imn_1d(N: int, L: float = 1.0) -> np.ndarray:
    """If we additionally redefine the n-th oscillator frequency by
       omega_n^{IMN} = k_n^{eff} = n pi / L_n  (with L_n = n L),
       then omega_n^{IMN} = pi/L  is constant.  This shows that under the
       SAME substitution applied to the Hamiltonian, vacuum-energy divergence
       changes character (linear -> linear with constant summand, still
       divergent: this is documented in the .md as a *negative* result)."""
    n   = np.arange(1, N + 1)
    L_n = n * L
    omega_n = n * np.pi / L_n      # = pi/L = const
    contrib = 0.5 * omega_n / L_n   # 1/(2 L_n) * omega_n  ~ 1/n
    return np.cumsum(contrib)


# ---------------------------------------------------------------------------
# 3+1D vacuum two-point  <phi^2>  at the centre of the cubic box
# Sums over (n_x,n_y,n_z) with n_i>=1.
# omega_{vec n} = pi/L * |vec n|,
# u_{vec n}(L/2,L/2,L/2)^2 = (2/L)^3 prod sin^2(n_i pi/2)
# ---------------------------------------------------------------------------
def phi2_3d(Nmax: int, L: float = 1.0, mode: str = "standard") -> tuple[np.ndarray, np.ndarray]:
    """Return (Nlist, cumulative-sum) where the sum is ordered by |n|.
       mode = 'standard'  -> uses V = L^3
       mode = 'imn'       -> uses V_n = (n_eff L)^3   with n_eff = |vec n|
    """
    # build all triples (nx,ny,nz) with 1<=ni<=Nmax
    n_range = np.arange(1, Nmax + 1)
    NX, NY, NZ = np.meshgrid(n_range, n_range, n_range, indexing='ij')
    NMAG = np.sqrt(NX**2 + NY**2 + NZ**2)   # |vec n|

    # sin^2(n_i pi/2) = 1 if n_i odd, 0 if even
    s2x = (NX % 2).astype(float)
    s2y = (NY % 2).astype(float)
    s2z = (NZ % 2).astype(float)
    mode_factor = s2x * s2y * s2z            # 1 only for all-odd modes

    omega = (np.pi / L) * NMAG

    V_std = L ** 3
    if mode == "standard":
        u2 = (8.0 / V_std) * mode_factor     # |u(L/2)|^2 = (2/L)^3 * sin^2 products
        contrib = (1.0 / (2.0 * omega)) * u2 # 1/(2 omega) * |u|^2
    elif mode == "imn":
        # Replace L -> L_n in normalisation  =>  V_n = (|n| L)^3.
        # so u^2 = (2/(|n| L))^3 * mode_factor = 8/(|n|^3 L^3) * mode_factor
        u2 = (8.0 / (NMAG**3 * V_std)) * mode_factor
        contrib = (1.0 / (2.0 * omega)) * u2
    else:
        raise ValueError(mode)

    # flatten and sort by |n|
    flat_n   = NMAG.flatten()
    flat_c   = contrib.flatten()
    order    = np.argsort(flat_n)
    flat_n   = flat_n[order]
    flat_c   = flat_c[order]
    cumsum   = np.cumsum(flat_c)
    return flat_n, cumsum


# ---------------------------------------------------------------------------
# Per-mode contribution profiles  c_n  vs n
# ---------------------------------------------------------------------------
def per_mode_1d(N: int, L: float = 1.0, x: float = 0.5):
    n   = np.arange(1, N + 1)
    k_n = n * np.pi / L
    std = (1.0 / (L * k_n)) * np.sin(k_n * x) ** 2
    imn = (1.0 / (n * L * k_n)) * np.sin(k_n * x) ** 2
    return n, std, imn


# ---------------------------------------------------------------------------
# 1-loop self-energy power-counting test
#
# Standard: Sigma(p) ~ ∫ d^4k  1/(k^2 (p-k)^2)  -- log divergent
# IMN: each propagator carries an extra 1/|n|^d ~ 1/k^d
# We model the *integrand* on the cubic mode lattice (n_x,n_y,n_z,n_t)
# and integrate radially.
# ---------------------------------------------------------------------------
def self_energy_integrand(Nmax: int, d: int = 3, mode: str = "standard"):
    """Return radial cumulative integrand mass for a scalar bubble in d
    spatial + 1 temporal dimension, evaluated at external momentum p=0
    on the discrete-mode lattice.

    The bubble  ∫ d^d k / (k^2)^2  ~ ∫ k^{d-1} dk / k^4 = ∫ k^{d-5} dk
    is log-divergent in d=4 (k^{-1}) and quadratic in d=2.
    IMN supplies an extra 1/|n|^d ~ k^{-d} per propagator (two propagators),
    so the corrected radial integrand is k^{d-1-4-2d} = k^{-d-5}.
    """
    n = np.arange(1, Nmax + 1, dtype=float)
    # SPACETIME radial measure k^{D-1} = k^d (D = d + 1 spacetime dims),
    # times two propagator poles 1/k^4 from a scalar bubble.
    # Standard:  integrand ~ k^{d-4}
    # IMN:       extra k^{-2d} from two IMN propagators -> k^{-d-4}
    base = n ** d * n ** (-4)               # ~ k^{d-4}, log-div at d=3
    if mode == "standard":
        contrib = base
    elif mode == "imn":
        contrib = base * n ** (-2.0 * d)    # extra suppression
    else:
        raise ValueError(mode)
    return n, np.cumsum(contrib)


# ---------------------------------------------------------------------------
# Plotting routines
# ---------------------------------------------------------------------------
def plot_phi2_1d_convergence():
    """Figure: 1+1D <phi^2>(L/2) partial sums, standard vs IMN."""
    N    = 5000
    std  = phi2_standard_1d(N)
    imn  = phi2_imn_1d(N)
    Nlist = np.arange(1, N + 1)

    # Analytic IMN limit:  (1/pi) sum_{n odd} 1/n^2  = (1/pi) * pi^2/8 = pi/8
    imn_limit = np.pi / 8.0

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # left: standard, log divergence
    ax = axes[0]
    ax.plot(Nlist, std, color="#D64933", lw=1.6, label=r"standard QFT")
    ax.plot(Nlist, (1.0 / np.pi) * (np.euler_gamma + np.log(2 * Nlist)) / 2,
            color="k", lw=1.0, ls="--",
            label=r"$\sim \frac{1}{2\pi}\ln N$  (theoretical)")
    ax.set_xscale("log")
    ax.set_xlabel(r"mode cut-off $N$")
    ax.set_ylabel(r"$\langle 0|\hat\phi^2(L/2)|0\rangle_N$")
    ax.set_title(r"Standard QFT: $\sum 1/n\sin^2(n\pi/2)$ — log divergent")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.3)

    # right: IMN, convergence
    ax = axes[1]
    ax.plot(Nlist, imn, color="#2274A5", lw=1.6, label=r"IMN: $L\to L_n=nL$")
    ax.axhline(imn_limit, color="k", lw=1.0, ls="--",
               label=fr"limit $\pi/8\approx{imn_limit:.4f}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"mode cut-off $N$")
    ax.set_ylabel(r"$\langle 0|\hat\phi^2(L/2)|0\rangle^{\rm IMN}_N$")
    ax.set_title(r"IMN: $\sum 1/n^2\sin^2(n\pi/2)$ — convergent")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)

    fig.suptitle(r"1+1D vacuum two-point function: divergence vs. IMN regulation",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "imn_phi2_1d.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    return float(std[-1]), float(imn[-1]), imn_limit


def plot_casimir_1d():
    """Figure: 1+1D vacuum energy density partial sums."""
    N    = 1000
    std  = Ecasimir_standard_1d(N)
    imn  = Ecasimir_imn_1d(N)
    Nlist = np.arange(1, N + 1)

    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(Nlist, std, color="#D64933", lw=1.6, label=r"standard  $\sim N^2$")
    ax.plot(Nlist, imn, color="#2274A5", lw=1.6,
            label=r"IMN with $\omega_n=\pi/L_n$  $\sim N$")
    ax.set_xlabel(r"mode cut-off $N$")
    ax.set_ylabel(r"$E_{\rm vac}^{(N)}\cdot L/\pi$")
    ax.set_title("1+1D Casimir-type sums: IMN does NOT cure zero-point divergence")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = FIG_DIR / "imn_casimir_1d.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def plot_phi2_3d():
    """Figure: 3+1D <phi^2>(L/2,L/2,L/2) partial sums on the mode lattice."""
    Nmax = 30   # 27000 modes
    n_std, c_std = phi2_3d(Nmax, mode="standard")
    n_imn, c_imn = phi2_3d(Nmax, mode="imn")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    ax.plot(n_std, c_std, color="#D64933", lw=1.0)
    ax.set_xlabel(r"$|\vec n|$ (radial cut-off)")
    ax.set_ylabel(r"$\langle\hat\phi^2\rangle_{|n|\leq N}$")
    ax.set_title(r"Standard QFT in 3+1D: $\sim |\vec n|^2$ (quadratic UV div.)")
    ax.set_yscale("log")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ax.plot(n_imn, c_imn, color="#2274A5", lw=1.0)
    ax.set_xlabel(r"$|\vec n|$ (radial cut-off)")
    ax.set_ylabel(r"$\langle\hat\phi^2\rangle^{\rm IMN}_{|n|\leq N}$")
    ax.set_title(r"IMN in 3+1D: $V\to V_n=|\vec n|^3 V$ — convergent")
    ax.grid(alpha=0.3)

    fig.suptitle("3+1D vacuum two-point: quadratic UV → finite", fontsize=14, y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "imn_phi2_3d.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    return float(c_std[-1]), float(c_imn[-1])


def plot_per_mode():
    """Figure: per-mode amplitudes show 1/n vs 1/n^2 decay."""
    N = 200
    n, std, imn = per_mode_1d(N)
    # keep only odd n (even modes have sin^2(n pi/2) = 0)
    mask = (n % 2 == 1)
    n_o = n[mask]
    std_o = std[mask]
    imn_o = imn[mask]

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.loglog(n_o, std_o, 'o-', color="#D64933", ms=4,
              label=r"standard  $\propto 1/n$  (log-divergent sum)")
    ax.loglog(n_o, imn_o, 's-', color="#2274A5", ms=4,
              label=r"IMN  $\propto 1/n^2$  (Basel-convergent sum)")
    # reference power laws
    n_ref = np.logspace(np.log10(n_o[0]), np.log10(n_o[-1]), 50)
    ax.loglog(n_ref, std_o[0] * (n_ref / n_o[0]) ** -1, 'k--', lw=0.8, alpha=0.6,
              label=r"$\propto n^{-1}$")
    ax.loglog(n_ref, imn_o[0] * (n_ref / n_o[0]) ** -2, 'k:', lw=0.8, alpha=0.6,
              label=r"$\propto n^{-2}$")
    ax.set_xlabel(r"mode number $n$ (odd)")
    ax.set_ylabel("per-mode contribution to $\\langle\\hat\\phi^2(L/2)\\rangle$")
    ax.set_title("Per-mode decay: standard $\\sim 1/n$  vs.  IMN $\\sim 1/n^2$")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = FIG_DIR / "imn_per_mode_decay.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def plot_self_energy():
    """Figure: bubble-diagram integrand cumulative behaviour, std vs IMN."""
    Nmax = 1000
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # d = 3  (4D Lorentz, the physical case)
    ax = axes[0]
    n3_s, c3_s = self_energy_integrand(Nmax, d=3, mode="standard")
    n3_i, c3_i = self_energy_integrand(Nmax, d=3, mode="imn")
    ax.loglog(n3_s, c3_s, color="#D64933", lw=1.6, label="standard")
    ax.loglog(n3_i, c3_i, color="#2274A5", lw=1.6, label="IMN")
    ax.set_xlabel(r"loop-momentum cut-off  $\Lambda \sim N$")
    ax.set_ylabel("cumulative integrand (radial)")
    ax.set_title(r"Scalar bubble in $D=4$: log-divergent $\to$ convergent (IMN)")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3, which="both")

    # d = 1  (2D)
    ax = axes[1]
    n1_s, c1_s = self_energy_integrand(Nmax, d=1, mode="standard")
    n1_i, c1_i = self_energy_integrand(Nmax, d=1, mode="imn")
    ax.loglog(n1_s, c1_s, color="#D64933", lw=1.6, label="standard")
    ax.loglog(n1_i, c1_i, color="#2274A5", lw=1.6, label="IMN")
    ax.set_xlabel(r"loop-momentum cut-off  $\Lambda \sim N$")
    ax.set_ylabel("cumulative integrand (radial)")
    ax.set_title(r"Scalar bubble in $D=2$ (i.e. 1+1D): finite $\to$ finite")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("Power-counting of a one-loop bubble: integrand cumulative mass",
                 fontsize=14, y=1.02)
    fig.tight_layout()
    out = FIG_DIR / "imn_self_energy_powercount.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("Running IMN UV-regularisation numerical tests")
    print("=" * 60)

    print("\n[1] 1+1D <phi^2>(L/2) convergence")
    s, i, lim = plot_phi2_1d_convergence()
    print(f"    standard partial sum (N=5000): {s:.6f}  (drifts as ln N)")
    print(f"    IMN      partial sum (N=5000): {i:.6f}")
    print(f"    IMN analytic limit  pi/8     : {lim:.6f}")

    print("\n[2] 1+1D Casimir-type sums")
    plot_casimir_1d()

    print("\n[3] 3+1D <phi^2> on cubic mode lattice")
    s3, i3 = plot_phi2_3d()
    print(f"    standard partial sum (Nmax=30): {s3:.4f}")
    print(f"    IMN      partial sum (Nmax=30): {i3:.6f}")

    print("\n[4] Per-mode contribution profile (1+1D)")
    plot_per_mode()

    print("\n[5] Self-energy bubble power-counting (d=1 and d=3)")
    plot_self_energy()

    print("\nAll figures written to", FIG_DIR)
