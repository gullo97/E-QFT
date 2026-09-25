r"""
================================================================================
 GALAXY COLLAPSE IN 3D  --  a cosmological simulation on a discrete grid of cells
 built *exactly* from the model of "The Expanding Quantum Box" (v15 thesis).
================================================================================

This script is a faithful 3-D extension of the thesis's own gravitational-collapse
simulation (Ch. 21, Animation 21.A, `ch21_collapse.py` / `newton_sim.py`).  Nothing
in the v15 folder is modified or imported; the physics below is transcribed from
the thesis text and reproduces its validated numbers before being lifted to 3-D.

--------------------------------------------------------------------------------
THE MODEL, IN ONE PAGE  (chapter references are to v15/thesis)
--------------------------------------------------------------------------------
Ch. 17  THE DICTIONARY.  A box of size L is *identically* a unit cell with spatial
        metric g = L^2; a time-dependent size L(t) lives in the FRW spacetime
        ds^2 = dt^2 - L(t)^2 dxi^2 with scale factor a(t) = L(t).  Each cell of
        space therefore carries its own local scale factor.

Ch. 18  THE CELL LATTICE.  Tile space with cells; cell shape is the emergent metric
        g_ij = sum_a e_a^i e_a^j (the triad).  The local log-cell-size is the
        gravitational variable; delta_alpha = ln a.

Ch. 21  CONSTRAINT CLOSURE  ->  FRIEDMANN + ATTRACTION.  Imposing the Hamiltonian
        constraint 3 kappa H^2 = rho slaves the ladder coupling to the local matter
        density, turning each cell into a Friedmann universe a ~ t^{2/3(1+w)}.
        Applied cell-by-cell to an inhomogeneous lattice (the "separate-universe"
        picture, exact in the long-wavelength limit) it yields EMERGENT ATTRACTION:

           * each cell x is its own FRW patch with comoving density Omega(x)=1+delta0(x)
           * matter-dominated acceleration (the closed-constraint EOM, eq. 21 of the
             thesis script):
                       a_ddot(x)  =  - (H0^2 / 2) * Omega(x) / a(x)^2
           * synchronized, constraint-consistent initial data:
                       a(x,0) = a0 ,   a_dot(x,0) = H0 * sqrt( Omega/a0 + (1-Omega) )
           * an over-dense cell (Omega>1) is a slightly CLOSED universe: it expands
             ever slower than the background, reaches TURNAROUND (a_dot=0) and
             RECOLLAPSES -- gravitational instability, with no force law postulated.
        Density contrast follows from mass conservation in the comoving cell:
                       delta(x,t) = (1+delta0(x)) (a_bg(t)/a(x,t))^3 - 1 .
        Linear theory (Cell-Growth Theorem, eq. 21.4): delta ~ a (growing mode);
        spherical-collapse thresholds: turnaround at delta_lin = 1.062, collapse at
        delta_lin = 1.686.

Ch. 24  NEWTON ON THE LATTICE.  The static weak field obeys the lattice Poisson
        equation  kappa * lap(delta_alpha) = -c * rho * delta,  delta_alpha = -Phi.
        We solve it (FFT, exactly as `newton_lattice.py`) to display the Newtonian
        potential WELL forming around the collapsing proto-galaxy.

--------------------------------------------------------------------------------
WHAT THIS SCRIPT DOES
--------------------------------------------------------------------------------
1. Lays down a spherically-symmetric Gaussian proto-galaxy over-density on a 3-D
   grid of cells (the thesis's exact profile delta0(r) = delta_c * exp(-(r/r0)^2)).
2. Evolves every cell as a separate-universe FRW patch by integrating the model's
   own acceleration equation (above) -- the genuine "simulation over a discrete
   grid cell".
3. Prints the validation numbers and checks them against the thesis (linear-growth
   exponent ~ 1, spherical-collapse thresholds, turnaround/collapse epochs).
4. Renders a 3-D animation of the galaxy collapsing: tracer matter (seeded on the
   comoving lattice, displaced by the model's own cell-collapse map -- the 3-D
   analogue of the thesis's 1-D "walls falling toward the over-density") forms a
   dense halo while the surrounding universe keeps expanding, beside the live
   density-contrast growth arc and the radial density profile.
5. Saves a static summary figure (including the Ch. 24 Newtonian potential well)
   and the trajectory archive.

Standalone: numpy / scipy / matplotlib only.   hbar = c = 1.
================================================================================
"""
import argparse
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.animation import FuncAnimation, PillowWriter

# galaxy colormap: dim outskirts visible on black -> white-hot collapsed core
GAL_CMAP = LinearSegmentedColormap.from_list(
    "galaxy", ["#1b2740", "#2f5e95", "#54a8d0", "#f4c430", "#fff4d6", "#ffffff"])
# fixed colour scale: saturate at delta ~ 40 so the collapsed core glows white-hot
# while moderate over-densities are orange and the smooth cloud stays dim-blue
COL_VMAX = np.log10(1.0 + 40.0)

# ----------------------------------------------------------------------------- #
#  PARAMETERS  (the dynamical ones are the thesis's; the rendering ones are free)
# ----------------------------------------------------------------------------- #
H0          = 1.0      # background Hubble rate at a0           (thesis: 1.0)
A0          = 0.02     # initial scale factor                   (thesis: 0.02)
AMIN_FRAC   = 0.25     # collapse floor a_min = AMIN_FRAC*a0    (thesis: 0.25)
DELTA_C     = 0.08     # central over-density contrast          (thesis: 0.08)
R0          = 4.0      # Gaussian radius in cells               (thesis: 4.0)
N_GRID      = 64       # cells per side of the 3-D lattice
A_BG_END    = 6000.0   # run until background has grown a_bg/a0 = this
N_FRAMES    = 160      # animation frames (uniform in log a_bg)
N_TRACERS   = 17000    # tracer-matter particles for the 3-D render
R_VIR       = 1.5      # virial-core display radius (collapsed cells form a finite
                       # halo, not a point: a standard visualization of virialization)
SEED        = 7

rng = np.random.default_rng(SEED)
AMIN = AMIN_FRAC * A0


# ----------------------------------------------------------------------------- #
#  BACKGROUND  --  the unperturbed (Omega=1) Friedmann patch, a_bg ~ t^{2/3}
# ----------------------------------------------------------------------------- #
def a_background(t):
    """Exact flat matter-dominated solution of a_ddot = -(H0^2/2)/a^2 from a0."""
    return (A0**1.5 + 1.5 * H0 * t) ** (2.0 / 3.0)


def t_of_abg(abg):
    """Invert a_background."""
    return (abg ** 1.5 - A0 ** 1.5) / (1.5 * H0)


T_END = t_of_abg(A_BG_END * A0)


# ----------------------------------------------------------------------------- #
#  THE CELL EQUATION  --  one separate-universe FRW patch (Ch. 21).  The thesis's
#  cell equation  a_ddot = -(H0^2/2) Omega / a^2  with the constraint-consistent IC
#  a(0)=a0, a_dot(0)=H0*sqrt(Omega/a0+(1-Omega))  has the EXACT closed-universe
#  cycloid as its solution when Omega>1 (an over-dense cell is a closed universe):
#
#       a(eta) = A (1 - cos eta),   t(eta) = B (eta - sin eta),
#       A = Omega / (2(Omega-1)),   B = Omega / (2 H0 (Omega-1)^{3/2}),
#
#  with a -> a_min freeze at recollapse.  We use this exact solution (cross-checked
#  below against a direct ODE integration) -- it is the model solved exactly, and
#  avoids the stiffness of integrating the recollapse numerically cell by cell.
# ----------------------------------------------------------------------------- #
def evolve_cell(Omega, t):
    """Exact a(t) for one cell of comoving density Omega over the time array t."""
    if Omega <= 1.0 + 1e-7:                 # background/under-dense: flat matter patch
        return np.maximum((A0**1.5 + 1.5 * H0 * t) ** (2.0 / 3.0), AMIN)
    A = Omega / (2.0 * (Omega - 1.0))
    B = Omega / (2.0 * H0 * (Omega - 1.0) ** 1.5)
    eta0 = np.arccos(np.clip(1.0 - A0 / A, -1.0, 1.0))      # a(0)=a0  ->  eta0
    t0 = B * (eta0 - np.sin(eta0))                          # sim t=0 maps to eta0
    eta_fr = 2 * np.pi - np.arccos(np.clip(1.0 - AMIN / A, -1.0, 1.0))   # a=a_min freeze
    t_fr = B * (eta_fr - np.sin(eta_fr))
    T = t0 + t
    eta = np.clip(np.cbrt(6.0 * np.maximum(T, 0.0) / B), 1e-6, eta_fr)   # small-eta guess
    for _ in range(80):                                    # vectorized Newton on t(eta)=T
        f = B * (eta - np.sin(eta)) - T
        fp = np.maximum(B * (1.0 - np.cos(eta)), 1e-12)
        eta = np.clip(eta - f / fp, 1e-7, eta_fr)
    a = np.where(T >= t_fr, AMIN, A * (1.0 - np.cos(eta)))
    return np.maximum(a, AMIN)


def crosscheck_central(Omega, t_max):
    """Direct ODE integration of the central cell up to ~turnaround, to confirm the
    exact cycloid solution. Returns (t, a) from solve_ivp (loose, non-stiff window)."""
    def rhs(tt, y):
        a, ad = y
        return [ad, -0.5 * H0**2 * Omega / max(a, 1e-6)**2]
    ad0 = H0 * np.sqrt(max(Omega / A0 + (1.0 - Omega), 1e-30))
    te = np.linspace(1e-3, t_max, 400)
    sol = solve_ivp(rhs, [1e-3, t_max], [A0, ad0], t_eval=te, method="RK45",
                    rtol=1e-8, atol=1e-11)
    return sol.t, sol.y[0]


# ----------------------------------------------------------------------------- #
#  SIMULATION  --  evolve the radial profile of cells, build derived fields
# ----------------------------------------------------------------------------- #
def run_simulation():
    print("=" * 78)
    print(" GALAXY COLLAPSE 3D  --  separate-universe cell lattice (v15 Ch. 21)")
    print("=" * 78)
    print(f"  grid {N_GRID}^3 cells | delta_c={DELTA_C} | r0={R0} cells | "
          f"a0={A0} | a_min={AMIN:.4f}")
    print(f"  integrating to a_bg/a0 = {A_BG_END:.0f}  (t_end = {T_END:.1f})")

    # time grid: uniform in log(a_bg/a0) so every expansion decade -- and the
    # turnaround/collapse phase -- is evenly sampled in the animation
    abg = np.geomspace(1.0, A_BG_END, N_FRAMES) * A0
    t = t_of_abg(abg)

    # radial bins covering the dynamic region (delta0 negligible beyond ~5 r0)
    r_dyn_max = 5.5 * R0
    r_bins = np.arange(0.0, r_dyn_max + 0.2, 0.2)
    delta0_bins = DELTA_C * np.exp(-(r_bins / R0) ** 2)     # thesis profile, in 3-D r
    Omega_bins = 1.0 + delta0_bins

    print(f"  evolving {len(r_bins)} radial cell-shells as separate-universe FRW "
          f"patches ...")
    a_bins = np.array([evolve_cell(Om, t) for Om in Omega_bins])   # (nbins, nt)

    # derived fields on the shells -------------------------------------------------
    # density contrast from comoving mass conservation (thesis eq.)
    delta_bins = (1.0 + delta0_bins)[:, None] * (abg[None, :] / a_bins) ** 3 - 1.0
    # log-cell-size perturbation = local gravitational potential variable (Ch. 18/24)
    dalpha_bins = np.log(a_bins / abg[None, :])

    central = delta_bins[0]                                 # central cell contrast
    # ---- validation against the thesis ------------------------------------------
    kta = int(np.argmax(a_bins[0]))                         # central turnaround
    a_ta = abg[kta] / A0
    hit = a_bins[0] <= AMIN * 1.0001
    kcol = int(np.argmax(hit)) if hit.any() else len(t) - 1
    a_col = abg[kcol] / A0
    win = (abg / A0 > 30) & (abg / A0 < 300)                # linear window
    growth_exp = np.polyfit(np.log(abg[win]), np.log(central[win]), 1)[0]
    Cfit = np.exp(np.mean(np.log(central[win]) - np.log(abg[win])))
    delta_lin_col = Cfit * abg[kcol]                        # linearly-extrapolated delta

    print("-" * 78)
    print("  VALIDATION  (compare to thesis Ch. 21)")
    print(f"    linear growth exponent  d ln delta / d ln a = {growth_exp:.3f}"
          f"      [thesis 1.08, theory 1]")
    print(f"    central turnaround at    a_bg/a0 = {a_ta:7.1f}"
          f"            [thesis 1168 for delta_c=0.08]")
    print(f"    central collapse at      a_bg/a0 = {a_col:7.1f}"
          f"            [thesis 1863 for delta_c=0.08]")
    print(f"    delta_lin extrapolated to collapse  = {delta_lin_col:.3f}"
          f"        [spherical-collapse threshold 1.686]")
    # cross-check the exact cycloid against a direct ODE integration of the central cell
    tcc, acc = crosscheck_central(Omega_bins[0], 1.35 * t[kta])
    a_ta_cyc = a_bins[0].max()
    a_ta_ode = acc.max()
    print(f"    cross-check (central cell): turnaround a  cycloid={a_ta_cyc:.4f} "
          f"vs direct ODE={a_ta_ode:.4f}  (rel.diff {abs(a_ta_cyc-a_ta_ode)/a_ta_cyc:.1e})")
    print(f"    shells collapsed by end of run: "
          f"{int(np.sum(a_bins[:, -1] <= AMIN*1.0001))}/{len(r_bins)} "
          f"(galaxy core radius ~ "
          f"{r_bins[np.max(np.where(a_bins[:, -1] <= AMIN*1.0001)[0]) if (a_bins[:,-1]<=AMIN*1.0001).any() else 0]:.1f} cells)")
    print("-" * 78)

    return dict(t=t, abg=abg, r_bins=r_bins, delta0_bins=delta0_bins,
                a_bins=a_bins, delta_bins=delta_bins, dalpha_bins=dalpha_bins,
                central=central, kta=kta, kcol=kcol, Cfit=Cfit,
                growth_exp=growth_exp, a_ta=a_ta, a_col=a_col)


# ----------------------------------------------------------------------------- #
#  TRACER MATTER  --  seed particles on the comoving lattice and displace them by
#  the model's own cell-collapse map.  For the spherically-symmetric proto-galaxy
#  the displacement is the 3-D radial analogue of the thesis's 1-D wall cumsum:
#  the background-normalized physical radius of a comoving shell R is
#       r_tilde(R,t) = (1/a_bg) * Integral_0^R a(r,t) dr
#  -- inner shells whose cells collapse (a -> a_min) fall toward the centre while
#  the outskirts keep riding the Hubble flow.
# ----------------------------------------------------------------------------- #
def seed_tracers(sim):
    r_bins = sim["r_bins"]
    # comoving radii: dense inside the over-density, sparser background halo
    R_core = rng.uniform(0, 1, int(N_TRACERS * 0.80)) ** (1 / 3) * (3.0 * R0)
    R_halo = rng.uniform(0, 1, N_TRACERS - R_core.size) ** (1 / 3) * (5.0 * R0)
    R = np.concatenate([R_core, R_halo])
    # isotropic directions on the sphere
    u = rng.uniform(-1, 1, R.size)
    phi = rng.uniform(0, 2 * np.pi, R.size)
    s = np.sqrt(1 - u**2)
    dirs = np.stack([s * np.cos(phi), s * np.sin(phi), u], axis=1)
    # per-particle virial radius (uniform in a ball of radius R_VIR): the size of the
    # finite halo a collapsed cell virializes into, instead of a point singularity
    rv = R_VIR * rng.uniform(0, 1, R.size) ** (1 / 3)
    return R, dirs, rv


def display_radius(rtilde, delta, rv):
    """Background-normalized radius, floored to the virial-core size as a cell
    collapses (smooth ramp in delta) so the galaxy is a finite glowing clump."""
    L = np.log10(1.0 + np.clip(delta, 0, None))
    g = np.clip((L - np.log10(9.0)) / (np.log10(151.0) - np.log10(9.0)), 0, 1)
    return np.maximum(rtilde, rv * g)


def embedding_radius(sim):
    """Precompute r_tilde(R,t) on the radial-bin grid: cumulative integral of a / a_bg."""
    r_bins, a_bins, abg = sim["r_bins"], sim["a_bins"], sim["abg"]
    # cumulative proper radius P(R,t) = int_0^R a dr  (trapezoid along r-axis)
    dr = np.gradient(r_bins)
    # trapezoidal cumulative sum over the radial axis for every time column
    P = np.cumsum(a_bins * dr[:, None], axis=0)            # (nbins, nt)
    P -= 0.5 * a_bins * dr[:, None]                        # midpoint correction
    rt = P / abg[None, :]                                   # (nbins, nt) background-normalized
    return rt


# ----------------------------------------------------------------------------- #
#  Ch. 24 NEWTONIAN WELL  --  FFT lattice-Poisson solve (delta_alpha = -Phi)
# ----------------------------------------------------------------------------- #
def potential_well(sim, frame):
    """Radial Newtonian potential delta_alpha=-Phi from the density field at `frame`."""
    N = N_GRID
    c = (N - 1) / 2.0
    ax = np.arange(N) - c
    RX, RY, RZ = np.meshgrid(ax, ax, ax, indexing="ij")
    r = np.sqrt(RX**2 + RY**2 + RZ**2)
    delta_field = np.interp(r.ravel(), sim["r_bins"], sim["delta_bins"][:, frame],
                            right=0.0).reshape(N, N, N)
    rho = delta_field - delta_field.mean()                 # periodic solvability
    k = 2 * np.pi * np.fft.fftfreq(N)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    khat2 = (2*np.sin(KX/2))**2 + (2*np.sin(KY/2))**2 + (2*np.sin(KZ/2))**2
    khat2[0, 0, 0] = 1.0
    F = np.fft.fftn(rho) / (-khat2)
    F[0, 0, 0] = 0.0
    phi = np.real(np.fft.ifftn(F))                         # lap phi = -rho  ->  Phi
    dalpha = -phi                                          # delta_alpha = -Phi > 0 near mass
    # radial average
    rb = np.arange(1, 0.45 * N, 1.0)
    rr = r.ravel(); dd = dalpha.ravel()
    prof = np.array([dd[(rr > x - 0.5) & (rr <= x + 0.5)].mean()
                     if np.any((rr > x - 0.5) & (rr <= x + 0.5)) else np.nan
                     for x in rb])
    prof -= dd[rr > 0.40 * N].mean()                       # remove far-field offset
    return rb, prof


# ----------------------------------------------------------------------------- #
#  ANIMATION
# ----------------------------------------------------------------------------- #
def make_animation(sim, out_gif, quick=False):
    t, abg = sim["t"], sim["abg"]
    central, Cfit = sim["central"], sim["Cfit"]
    nt = len(t)
    R, dirs, rv = seed_tracers(sim)
    rt = embedding_radius(sim)                              # (nbins, nt)

    # per-particle: background-normalized radius and local contrast vs time
    # interpolate r_tilde(R,t) and delta(R,t) at each particle's comoving radius R
    from scipy.interpolate import interp1d
    rtilde_of_R = np.array([np.interp(R, sim["r_bins"], rt[:, k]) for k in range(nt)])  # (nt, Np)
    delta_of_R = np.array([np.interp(R, sim["r_bins"], sim["delta_bins"][:, k],
                                     right=0.0) for k in range(nt)])                    # (nt, Np)

    frames = list(range(nt))
    if quick:
        frames = frames[::6]

    fig = plt.figure(figsize=(13.2, 6.6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.65, 1.0], height_ratios=[1, 1],
                           wspace=0.22, hspace=0.32,
                           left=0.02, right=0.965, top=0.91, bottom=0.10)
    ax3d = fig.add_subplot(gs[:, 0], projection="3d")
    axg = fig.add_subplot(gs[0, 1])
    axp = fig.add_subplot(gs[1, 1])

    # --- static decoration of the growth panel (Ch. 21 Fig 21.3) ----------------
    axg.set_xscale("log"); axg.set_yscale("log")
    axg.plot(abg / A0, np.where(t <= t[sim["kcol"]], central, np.nan),
             color="0.55", lw=1.4, zorder=2)
    axg.plot(abg / A0, Cfit * abg, "k--", lw=1.0, label=r"linear theory $\delta\propto a$")
    axg.axhline(1.686, color="crimson", ls=":", lw=1.2, label=r"collapse $\delta_{\rm lin}=1.686$")
    axg.axhline(1.062, color="darkorange", ls=":", lw=1.0, label=r"turnaround $1.062$")
    axg.set_xlim(1, A_BG_END); axg.set_ylim(2e-2, 6e2)
    axg.set_xlabel(r"background expansion  $a_{\rm bg}/a_0$")
    axg.set_ylabel(r"central density contrast  $\delta$")
    axg.set_title("gravitational growth  (Ch. 21)", fontsize=10)
    axg.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
    axg.grid(alpha=0.25)
    gdot, = axg.plot([], [], "o", color="crimson", ms=9, zorder=5)

    rmax_disp = float(rtilde_of_R[0].max()) * 1.03

    def update(k):
        ax3d.clear()
        dlt = delta_of_R[k]
        rr = display_radius(rtilde_of_R[k], dlt, rv)        # virial-floored radius
        x = rr[:, None] * dirs                              # (Np,3) display positions
        col = np.log10(1.0 + np.clip(dlt, 0, None))
        frac = np.clip(col / COL_VMAX, 0, 1)
        order = np.argsort(dlt)                             # draw dense core last (on top)
        sz = 1.4 + 26.0 * frac ** 2
        ax3d.scatter(x[order, 0], x[order, 1], x[order, 2],
                     c=col[order], cmap=GAL_CMAP, vmin=0, vmax=COL_VMAX,
                     s=sz[order], alpha=0.72, edgecolors="none", depthshade=True)
        L = rmax_disp
        ax3d.set_xlim(-L, L); ax3d.set_ylim(-L, L); ax3d.set_zlim(-L, L)
        ax3d.set_box_aspect((1, 1, 1))
        ax3d.set_xticks([]); ax3d.set_yticks([]); ax3d.set_zticks([])
        ax3d.set_facecolor("black")
        ax3d.xaxis.set_pane_color((0, 0, 0, 1)); ax3d.yaxis.set_pane_color((0, 0, 0, 1))
        ax3d.zaxis.set_pane_color((0, 0, 0, 1))
        ax3d.grid(False)
        ax3d.view_init(elev=22, azim=-65 + 0.18 * k)
        dc = central[k]
        ax3d.set_title("emergent gravitational collapse of a proto-galaxy\n"
                       "(matter falling toward the over-density on the cell lattice)\n"
                       fr"$a_{{\rm bg}}/a_0={abg[k]/A0:6.0f}$   "
                       fr"central $\delta={min(dc,9999):7.2f}$   "
                       "(colour: local over-density)", fontsize=10, color="0.15")

        # growth panel marker
        gdot.set_data([abg[k] / A0], [float(np.clip(central[k], 2.2e-2, 5.5e2))])

        # radial density-profile panel
        axp.clear()
        prof = sim["delta_bins"][:, k]
        axp.plot(sim["r_bins"], np.clip(prof, 1e-3, None), color="crimson", lw=2)
        axp.fill_between(sim["r_bins"], 1e-3, np.clip(prof, 1e-3, None),
                         color="crimson", alpha=0.18)
        axp.set_yscale("log")
        axp.set_xlim(0, 4.2 * R0); axp.set_ylim(1e-2, 6e2)
        axp.set_xlabel("comoving radius  $r$  (cells)")
        axp.set_ylabel(r"density contrast  $\delta(r)$")
        axp.set_title("the galaxy forming: radial profile", fontsize=10)
        axp.grid(alpha=0.25)
        return ()

    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    fps = 12
    print(f"  rendering {len(frames)} frames -> {out_gif} ...")
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=(70 if quick else 96))
    plt.close(fig)
    print(f"  saved {out_gif}")


# ----------------------------------------------------------------------------- #
#  STATIC SUMMARY FIGURE
# ----------------------------------------------------------------------------- #
def make_summary(sim, out_png):
    t, abg = sim["t"], sim["abg"]
    fig = plt.figure(figsize=(13.5, 8.2))
    gs = gridspec.GridSpec(2, 3, wspace=0.30, hspace=0.34,
                           left=0.07, right=0.97, top=0.92, bottom=0.08)

    # (a) growth arc
    a = fig.add_subplot(gs[0, 0])
    a.loglog(abg / A0, np.where(t <= t[sim["kcol"]], sim["central"], np.nan),
             color="crimson", lw=2, label=r"central $\delta(t)$ (sim)")
    a.loglog(abg / A0, sim["Cfit"] * abg, "k--", lw=1, label=r"linear $\delta\propto a$")
    a.axhline(1.686, color="gray", ls=":", label=r"$\delta_{\rm lin}=1.686$")
    a.axvline(sim["a_ta"], color="darkorange", lw=0.8); a.text(sim["a_ta"]*1.05, 3e-2,
              "turnaround", rotation=90, fontsize=8, color="darkorange")
    a.set_xlabel(r"$a_{\rm bg}/a_0$"); a.set_ylabel(r"density contrast $\delta$")
    a.set_title(f"(a) growth, turnaround, collapse\n"
                f"exponent {sim['growth_exp']:.2f}  (theory 1)", fontsize=10)
    a.legend(fontsize=8); a.grid(alpha=0.3)

    # (b) radial profiles at epochs
    b = fig.add_subplot(gs[0, 1])
    ks = [0, sim["kta"]//2, sim["kta"], (sim["kta"]+sim["kcol"])//2, sim["kcol"]]
    for k, col in zip(ks, plt.cm.plasma(np.linspace(0.05, 0.92, len(ks)))):
        b.semilogy(sim["r_bins"], np.clip(sim["delta_bins"][:, k], 1e-3, None),
                   color=col, lw=1.8, label=fr"$a_{{\rm bg}}/a_0={abg[k]/A0:.0f}$")
    b.set_xlim(0, 4.2*R0); b.set_ylim(1e-2, 6e2)
    b.set_xlabel("comoving radius $r$ (cells)"); b.set_ylabel(r"$\delta(r)$")
    b.set_title("(b) the over-density sharpening into a galaxy", fontsize=10)
    b.legend(fontsize=7.5); b.grid(alpha=0.3)

    # (c) cell scale-factor profile a_i/a_bg (the well, Ch. 21 Fig 21.3 right)
    c = fig.add_subplot(gs[0, 2])
    for k, col in zip(ks, plt.cm.viridis(np.linspace(0.05, 0.92, len(ks)))):
        c.semilogy(sim["r_bins"], sim["a_bins"][:, k] / abg[k], color=col, lw=1.8,
                   label=fr"$a_{{\rm bg}}/a_0={abg[k]/A0:.0f}$")
    c.set_xlim(0, 4.2*R0)
    c.set_xlabel("comoving radius $r$ (cells)"); c.set_ylabel(r"$a_i/a_{\rm bg}$")
    c.set_title("(c) local scale factor lags &\ncollapses (the well deepens)", fontsize=10)
    c.legend(fontsize=7.5); c.grid(alpha=0.3)

    # (d) Ch. 24 Newtonian potential well at a representative epoch
    d = fig.add_subplot(gs[1, 0])
    rb, prof = potential_well(sim, sim["kta"])
    good = prof > 0
    d.loglog(rb[good], prof[good], "o", color="crimson", ms=5,
             label=r"lattice $\delta\alpha=-\Phi$")
    if good.sum() > 3:
        ref = rb[good]
        d.loglog(ref, prof[good][2]*ref[2]/ref, "k--", lw=1, label=r"$\propto 1/r$")
    d.set_xlabel("cell distance $r$"); d.set_ylabel(r"$\delta\alpha=-\Phi$")
    d.set_title("(d) Newtonian potential well\n(Ch. 24 lattice Poisson, FFT)", fontsize=10)
    d.legend(fontsize=8); d.grid(alpha=0.3, which="both")

    # (e) 3-D snapshot near collapse
    e = fig.add_subplot(gs[1, 1], projection="3d")
    R, dirs, rv = seed_tracers(sim)
    rt = embedding_radius(sim)
    ksnap = min(sim["kcol"] + 4, len(t) - 1)
    rtil = np.interp(R, sim["r_bins"], rt[:, ksnap])
    dl = np.interp(R, sim["r_bins"], sim["delta_bins"][:, ksnap], right=0.0)
    rr = display_radius(rtil, dl, rv)
    x = rr[:, None] * dirs
    col = np.log10(1 + np.clip(dl, 0, None)); frac = np.clip(col / COL_VMAX, 0, 1)
    order = np.argsort(dl)
    e.scatter(x[order, 0], x[order, 1], x[order, 2], c=col[order], cmap=GAL_CMAP,
              vmin=0, vmax=COL_VMAX, s=1.4 + 26*frac[order]**2,
              alpha=0.72, edgecolors="none")
    L = float(np.interp(R, sim["r_bins"], rt[:, 0]).max())
    e.set_xlim(-L, L); e.set_ylim(-L, L); e.set_zlim(-L, L); e.set_box_aspect((1, 1, 1))
    e.set_xticks([]); e.set_yticks([]); e.set_zticks([])
    e.set_facecolor("black")
    for pane in (e.xaxis, e.yaxis, e.zaxis): pane.set_pane_color((0, 0, 0, 1))
    e.set_title(f"(e) collapsed proto-galaxy\n$a_{{\\rm bg}}/a_0={abg[ksnap]/A0:.0f}$",
                fontsize=10)

    # (f) initial vs final 3-D (contrast)
    f = fig.add_subplot(gs[1, 2], projection="3d")
    rr0 = np.interp(R, sim["r_bins"], rt[:, 0])
    x0 = rr0[:, None] * dirs
    f.scatter(x0[:, 0], x0[:, 1], x0[:, 2], c="0.6", s=2, alpha=0.5, edgecolors="none")
    f.set_xlim(-L, L); f.set_ylim(-L, L); f.set_zlim(-L, L); f.set_box_aspect((1, 1, 1))
    f.set_xticks([]); f.set_yticks([]); f.set_zticks([])
    f.set_facecolor("black")
    for pane in (f.xaxis, f.yaxis, f.zaxis): pane.set_pane_color((0, 0, 0, 1))
    f.set_title("(f) initial smooth cloud\n(comoving, $a_{\\rm bg}/a_0=1$)", fontsize=10)

    fig.suptitle("Emergent gravitational collapse of a galaxy on a discrete cell lattice "
                 "— The Expanding Quantum Box (v15, Ch. 17–24)", fontsize=12.5)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    print(f"  saved {out_png}")


# ----------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="fast low-res test render")
    ap.add_argument("--no-anim", action="store_true", help="skip the GIF")
    args = ap.parse_args()

    sim = run_simulation()

    out_dir = "."
    np.savez_compressed(f"{out_dir}/galaxy_collapse_traj.npz",
                        t=sim["t"], abg=sim["abg"], r_bins=sim["r_bins"],
                        a_bins=sim["a_bins"], delta_bins=sim["delta_bins"],
                        dalpha_bins=sim["dalpha_bins"])
    print("  saved galaxy_collapse_traj.npz")

    make_summary(sim, f"{out_dir}/galaxy_collapse_summary.png")
    if not args.no_anim:
        make_animation(sim, f"{out_dir}/galaxy_collapse_3d.gif", quick=args.quick)
    print("done.")


if __name__ == "__main__":
    main()
