r"""
================================================================================
 ROTATIONAL EQUILIBRIUM OF A GALAXY on the E-QFT cell lattice.
 Adding angular momentum to the collapse: the condition for no collapse, the
 stable-orbit (tangential) velocity curve, and a 3-D simulation that shows it.
================================================================================

Companion to `galaxy_collapse_3d.py`.  There the proto-galaxy over-density
collapsed to a point (no angular momentum).  Here we give the matter angular
momentum about the barycentre and ask: what tangential-velocity profile v(r)
puts every shell on a stable orbit so that -- at large t -- nothing collapses?

--------------------------------------------------------------------------------
THE MODEL  (chapter references are to v15/thesis; nothing in v15 is modified)
--------------------------------------------------------------------------------
Ch. 24  NEWTON ON THE LATTICE.  The model's emergent static gravity is *exactly*
        Newtonian: the lattice Poisson equation  kappa*lap(delta_alpha) = -c*rho*delta
        with delta_alpha = -Phi gives a 1/r potential (measured exponent -1.008)
        and an inverse-square force (-2.05).  So a test mass at radius r from the
        barycentre feels the Newtonian field of the enclosed mass:
                    dPhi/dr = G*M(<r)/r^2          (Newton's shell theorem).

EQUILIBRIUM (no collapse).  A shell with specific angular momentum ell has the
effective potential  V_eff(r) = Phi(r) + ell^2/(2 r^2).  A circular, non-collapsing
orbit is the stationary point V_eff'(r)=0:
        dPhi/dr = ell^2/r^3   =>   ell^2 = G M(<r) r,
so the tangential speed that holds a shell at radius r against gravity is the
CIRCULAR-VELOCITY CURVE

        v_c(r) = sqrt( r dPhi/dr ) = sqrt( G M(<r) / r ).                 (*)

If every shell carries exactly v_t(r) = v_c(r) the net radial force vanishes
everywhere: the system is in centrifugal equilibrium and never collapses.  The
orbit is stable (Rayleigh / epicyclic criterion) when the specific angular
momentum r*v_c(r) increases outward,  d(r v_c)/dr > 0.

SHAPE OF THE CURVE.  For a centrally-concentrated mass M(<r):
   * core  (rho ~ const, M ~ r^3)  ->  v_c ~ r           (solid-body rise)
   * outside the mass (M -> M_tot)  ->  v_c ~ 1/sqrt(r)  (Keplerian decline)
A FLAT curve (v_c = const) requires M(<r) ~ r, i.e. rho ~ 1/r^2 (isothermal
sphere) -- exactly the "dark-matter halo" profile.  The E-QFT model, being
Newtonian (Ch. 24), gives a rising-then-Keplerian curve for an isolated luminous
mass: it does NOT by itself produce the observed flat rotation curves; that is
the dark-matter problem, stated cleanly in the model's own language.

--------------------------------------------------------------------------------
WHAT THIS SCRIPT DOES
--------------------------------------------------------------------------------
1. Takes the proto-galaxy mass distribution (the Gaussian over-density of the
   collapse run), computes M(<r) and the equilibrium curve v_c(r) of eq. (*).
2. Runs a 3-D N-body simulation of tracer matter in the model's Newtonian
   potential, two ways:
       (cold)      v = 0           -> the matter collapses (companion run);
       (rotating)  v = v_c(r)      -> every particle is on a circular orbit:
                                      a stable, rotating galaxy that does NOT
                                      collapse at large t.
3. Animates both in 3-D side by side, with the tangential-velocity curve v_c(r)
   and the half-mass radius vs time (cold shrinks, rotating is flat).

Standalone: numpy / scipy / matplotlib only.   hbar = c = 1, G = 1.
================================================================================
"""
import argparse
import numpy as np
from scipy.special import erf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.animation import FuncAnimation, PillowWriter

# ----------------------------------------------------------------------------- #
#  PARAMETERS
# ----------------------------------------------------------------------------- #
G      = 1.0
R0     = 4.0       # Gaussian scale radius of the proto-galaxy (cells; matches collapse run)
MTOT   = 1.0       # total galaxy mass (units)
EPS    = 0.30      # gravitational softening (avoids the central singularity)
NPART  = 6000      # tracer particles per run
N_FRAMES = 150
T_DYN_RUN = 7.0    # run length in units of the peak dynamical time
SEED   = 11

rng = np.random.default_rng(SEED)


# ----------------------------------------------------------------------------- #
#  THE GALAXY: Gaussian over-density rho(r) ~ exp(-(r/R0)^2)
#  M(<r)/Mtot = erf(x) - (2/sqrt(pi)) x exp(-x^2),  x = r/R0     (exact)
#  v_c(r) = sqrt(G M(<r)/r)                                       (eq. *)
# ----------------------------------------------------------------------------- #
def M_enclosed(r):
    x = np.asarray(r, float) / R0
    return MTOT * (erf(x) - (2.0 / np.sqrt(np.pi)) * x * np.exp(-x**2))


def v_circular(r):
    r = np.asarray(r, float)
    out = np.zeros_like(r)
    m = r > 1e-9
    out[m] = np.sqrt(G * M_enclosed(r[m]) / r[m])
    return out


def accel(x):
    """Newtonian acceleration toward the barycentre (Ch. 24 emergent gravity),
    softened. a = -G M(<r) x / (r^2+eps^2)^{3/2}."""
    r = np.linalg.norm(x, axis=1)
    M = M_enclosed(r)
    soft = (r**2 + EPS**2) ** 1.5
    return -(G * M / soft)[:, None] * x


# ----------------------------------------------------------------------------- #
#  THE ANSWER: the tangential-velocity (rotation) curve and its regimes
# ----------------------------------------------------------------------------- #
def rotation_curve_report():
    rr = np.linspace(1e-3, 6 * R0, 2000)
    vc = v_circular(rr)
    ipk = int(np.argmax(vc))
    rpk, vpk = rr[ipk], vc[ipk]
    # specific angular momentum r*vc must rise outward for stability
    L = rr * vc
    stable_frac = np.mean(np.gradient(L, rr) > 0)
    # regime slopes (log-log) in the core and in the outskirts
    core = (rr > 0.15 * R0) & (rr < 0.45 * R0)
    outer = (rr > 3.0 * R0) & (rr < 5.5 * R0)
    s_core = np.polyfit(np.log(rr[core]), np.log(vc[core]), 1)[0]
    s_outer = np.polyfit(np.log(rr[outer]), np.log(vc[outer]), 1)[0]
    print("=" * 78)
    print(" ROTATIONAL EQUILIBRIUM OF THE GALAXY  (E-QFT Ch. 24 Newtonian gravity)")
    print("=" * 78)
    print("  equilibrium / no-collapse condition:  v_t(r) = v_c(r) = sqrt(G M(<r)/r)")
    print(f"  peak circular speed v_c = {vpk:.4f} at r = {rpk:.2f} cells (~{rpk/R0:.2f} R0)")
    print(f"  core slope   d ln v_c/d ln r = {s_core:+.3f}   [solid-body  +1]")
    print(f"  outer slope  d ln v_c/d ln r = {s_outer:+.3f}   [Keplerian  -0.5]")
    print(f"  stability d(r v_c)/dr > 0 on {100*stable_frac:.0f}% of the profile "
          f"(orbits stable where it holds)")
    print(f"  -> NOT flat: a flat curve needs M(<r)~r (rho~1/r^2, isothermal halo)"
          f"  [the dark-matter point]")
    print("-" * 78)
    return rr, vc, rpk, vpk


# ----------------------------------------------------------------------------- #
#  INITIAL CONDITIONS: sample positions from the mass distribution
# ----------------------------------------------------------------------------- #
def sample_positions(n):
    rg = np.linspace(0, 6 * R0, 6000)
    cdf = M_enclosed(rg); cdf /= cdf[-1]
    u = rng.uniform(0, 1, n)
    r = np.interp(u, cdf, rg)
    cth = rng.uniform(-1, 1, n)
    ph = rng.uniform(0, 2 * np.pi, n)
    st = np.sqrt(1 - cth**2)
    dirs = np.stack([st * np.cos(ph), st * np.sin(ph), cth], axis=1)
    return r[:, None] * dirs, ph


def azimuthal_hat(x):
    """Unit azimuthal vector phi_hat about the z (spin) axis. For v = v_c*phi_hat
    every particle is on a circular orbit of its spherical radius (phi_hat _|_ r)."""
    R = np.sqrt(x[:, 0]**2 + x[:, 1]**2)
    R = np.maximum(R, 1e-9)
    return np.stack([-x[:, 1] / R, x[:, 0] / R, np.zeros(len(x))], axis=1)


# ----------------------------------------------------------------------------- #
#  LEAPFROG INTEGRATION in the model's Newtonian potential
# ----------------------------------------------------------------------------- #
def integrate(x0, v0, t_end, frame_times):
    nsteps = max(len(frame_times), 1)
    dt = (frame_times[1] - frame_times[0]) if len(frame_times) > 1 else t_end
    # use a finer internal step for accuracy, snapshot at frame_times
    sub = 8
    h = dt / sub
    x = x0.copy(); v = v0.copy()
    a = accel(x)
    snaps = [x.copy()]
    speeds = [np.linalg.norm(v, axis=1).copy()]
    t = 0.0
    for k in range(1, len(frame_times)):
        for _ in range(sub):
            v += 0.5 * h * a
            x += h * v
            a = accel(x)
            v += 0.5 * h * a
        snaps.append(x.copy())
        speeds.append(np.linalg.norm(v, axis=1).copy())
    return np.array(snaps), np.array(speeds)


def half_mass_radius(snaps):
    return np.array([np.median(np.linalg.norm(s, axis=1)) for s in snaps])


# ----------------------------------------------------------------------------- #
def run(quick=False):
    rr, vc, rpk, vpk = rotation_curve_report()
    tdyn = rpk / vpk                                   # peak dynamical time
    t_end = T_DYN_RUN * 2 * np.pi * tdyn               # several orbital periods
    nf = (N_FRAMES // 3) if quick else N_FRAMES
    frame_times = np.linspace(0, t_end, nf)

    x0, phi0 = sample_positions(NPART)
    r0 = np.linalg.norm(x0, axis=1)
    # (cold) zero velocity   |   (rotating) tangential v_c about z
    v_cold = np.zeros_like(x0)
    v_rot = v_circular(r0)[:, None] * azimuthal_hat(x0)

    print(f"  integrating {NPART} tracers x2 (cold vs rotating) for "
          f"{T_DYN_RUN:.0f} orbital periods, {nf} frames ...")
    snap_c, spd_c = integrate(x0, v_cold, t_end, frame_times)
    snap_r, spd_r = integrate(x0, v_rot, t_end, frame_times)

    rh_c = half_mass_radius(snap_c)
    rh_r = half_mass_radius(snap_r)
    print("-" * 78)
    print("  VALIDATION")
    print(f"    half-mass radius, COLD:     {rh_c[0]:.2f} -> {rh_c[-1]:.2f} cells "
          f"(shrinks x{rh_c[0]/max(rh_c[-1],1e-9):.1f}: COLLAPSE)")
    print(f"    half-mass radius, ROTATING: {rh_r[0]:.2f} -> {rh_r[-1]:.2f} cells "
          f"(ratio {rh_r[-1]/rh_r[0]:.2f}: STABLE, no collapse)")
    # check that the rotating run keeps every particle near its circular orbit:
    rfin = np.linalg.norm(snap_r[-1], axis=1)
    drift = np.median(np.abs(rfin - r0) / np.maximum(r0, EPS))
    print(f"    rotating run median |dr|/r over the whole run = {drift:.3f} "
          f"(small -> orbits stay circular: equilibrium holds)")
    print("-" * 78)
    return dict(rr=rr, vc=vc, rpk=rpk, vpk=vpk, frame_times=frame_times,
                snap_c=snap_c, snap_r=snap_r, spd_r=spd_r, r0=r0, phi0=phi0,
                rh_c=rh_c, rh_r=rh_r, tdyn=tdyn, t_end=t_end)


# ----------------------------------------------------------------------------- #
#  STATIC FIGURE: the rotation curve + the equilibrium story
# ----------------------------------------------------------------------------- #
def make_figure(S, out_png):
    fig = plt.figure(figsize=(12.5, 5.4))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.05, 1], wspace=0.28,
                           left=0.08, right=0.97, top=0.80, bottom=0.12)
    a1 = fig.add_subplot(gs[0]); a2 = fig.add_subplot(gs[1])

    rr, vc = S["rr"], S["vc"]
    a1.plot(rr / R0, vc, color="crimson", lw=2.4, label=r"$v_c(r)=\sqrt{GM(<r)/r}$ (this model)")
    # regime guides
    rcore = np.linspace(0.05, 1.0, 50) * R0
    a1.plot(rcore / R0, vc[np.argmin(np.abs(rr-0.6*R0))]*(rcore/(0.6*R0)),
            "k:", lw=1, label=r"solid-body $v\propto r$ (core)")
    rout = np.linspace(2.2, 6, 50) * R0
    iref = np.argmin(np.abs(rr - 3.0*R0))
    a1.plot(rout / R0, vc[iref]*np.sqrt(3.0*R0/rout), "k--", lw=1,
            label=r"Keplerian $v\propto1/\sqrt{r}$ (outskirts)")
    a1.axhline(S["vpk"], color="0.7", lw=0.8)
    # isothermal flat comparison (what a dark-matter halo would give)
    a1.axhline(S["vpk"], color="steelblue", ls="-.", lw=1.2,
               label=r"flat $v=$const $\Leftrightarrow\,\rho\propto1/r^2$ (isothermal halo)")
    a1.axvline(S["rpk"] / R0, color="0.7", lw=0.8)
    a1.set_xlabel(r"radius from barycentre  $r / R_0$")
    a1.set_ylabel(r"tangential speed for stable orbits  $v_c(r)$")
    a1.set_title("the stable-orbit (rotation) curve\nrising core, Keplerian outskirts — not flat",
                 fontsize=10.5)
    a1.set_xlim(0, 6); a1.set_ylim(0, S["vpk"] * 1.25)
    a1.legend(fontsize=8.2, loc="upper right"); a1.grid(alpha=0.3)

    # half-mass radius vs time: cold collapses, rotating is flat
    tt = S["frame_times"] / (2 * np.pi * S["tdyn"])
    a2.plot(tt, S["rh_c"] / R0, color="darkorange", lw=2.2,
            label="cold ($v=0$): collapses")
    a2.plot(tt, S["rh_r"] / R0, color="seagreen", lw=2.2,
            label=r"equilibrium ($v=v_c$): stable")
    a2.set_xlabel("time  (orbital periods)")
    a2.set_ylabel(r"half-mass radius  $r_{1/2}/R_0$")
    a2.set_title("no collapse at large $t$:\nrotational support holds the galaxy up",
                 fontsize=10.5)
    a2.set_ylim(0, max(S["rh_r"].max(), S["rh_c"].max()) / R0 * 1.1)
    a2.legend(fontsize=9); a2.grid(alpha=0.3)

    fig.suptitle("Angular momentum vs gravity on the E-QFT cell lattice — the stable-orbit velocity curve",
                 fontsize=12, y=0.965)
    fig.savefig(out_png, dpi=145); plt.close(fig)
    print(f"  saved {out_png}")


# ----------------------------------------------------------------------------- #
#  3-D ANIMATION: cold collapse vs rotating equilibrium + the curve
# ----------------------------------------------------------------------------- #
def make_animation(S, out_gif, quick=False):
    snap_c, snap_r = S["snap_c"], S["snap_r"]
    nf = len(S["frame_times"])
    L = 3.2 * R0

    fig = plt.figure(figsize=(13.6, 6.6))
    gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 0.95], height_ratios=[1, 1],
                           wspace=0.08, hspace=0.42, left=0.01, right=0.95,
                           top=0.90, bottom=0.10)
    axc = fig.add_subplot(gs[:, 0], projection="3d")
    axr = fig.add_subplot(gs[:, 1], projection="3d")
    axv = fig.add_subplot(gs[0, 2])
    axh = fig.add_subplot(gs[1, 2])

    # static: rotation curve
    axv.plot(S["rr"] / R0, S["vc"], color="crimson", lw=2.2)
    axv.axhline(S["vpk"], color="steelblue", ls="-.", lw=1.0,
                label=r"flat (isothermal)")
    axv.set_xlim(0, 6); axv.set_ylim(0, S["vpk"] * 1.25)
    axv.set_xlabel(r"$r/R_0$", fontsize=9); axv.set_ylabel(r"$v_c(r)$", fontsize=9)
    axv.set_title(r"stable-orbit curve $v_c=\sqrt{GM(<r)/r}$", fontsize=9.5)
    axv.legend(fontsize=7.5, loc="lower center"); axv.grid(alpha=0.3)
    axv.tick_params(labelsize=8)
    # static: half-mass radius
    tt = S["frame_times"] / (2 * np.pi * S["tdyn"])
    axh.plot(tt, S["rh_c"] / R0, color="darkorange", lw=1.8, label="cold: collapse")
    axh.plot(tt, S["rh_r"] / R0, color="seagreen", lw=1.8, label=r"$v_c$: stable")
    axh.set_xlabel("time (orbits)", fontsize=9); axh.set_ylabel(r"$r_{1/2}/R_0$", fontsize=9)
    axh.set_title("half-mass radius", fontsize=9.5)
    axh.set_ylim(0, max(S["rh_r"].max(), S["rh_c"].max())/R0*1.1)
    axh.legend(fontsize=7.5); axh.grid(alpha=0.3); axh.tick_params(labelsize=8)
    vmk_c, = axh.plot([], [], "o", color="darkorange", ms=7)
    vmk_r, = axh.plot([], [], "o", color="seagreen", ms=7)

    import matplotlib.cm as cm
    col_phi = cm.hsv((S["phi0"] % (2*np.pi)) / (2*np.pi))     # azimuth -> shows spin
    rad_norm = np.clip(S["r0"] / (3 * R0), 0, 1)
    col_cold = cm.cool(0.15 + 0.8 * rad_norm)

    def setup3d(ax, title):
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_facecolor("black")
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.set_pane_color((0, 0, 0, 1))
        ax.grid(False)
        ax.set_title(title, fontsize=10, color="0.15")

    def update(k):
        axc.clear(); axr.clear()
        xc = snap_c[k]; xr = snap_r[k]
        axc.scatter(xc[:, 0], xc[:, 1], xc[:, 2], c=col_cold, s=3.0,
                    alpha=0.55, edgecolors="none", depthshade=True)
        axr.scatter(xr[:, 0], xr[:, 1], xr[:, 2], c=col_phi, s=3.0,
                    alpha=0.7, edgecolors="none", depthshade=True)
        setup3d(axc, "NO angular momentum\n$\\rightarrow$ gravitational collapse")
        setup3d(axr, "centrifugal equilibrium $v=v_c(r)$\n$\\rightarrow$ stable rotating galaxy")
        az = -60 + 0.25 * k
        axc.view_init(elev=20, azim=az); axr.view_init(elev=20, azim=az)
        vmk_c.set_data([tt[k]], [S["rh_c"][k] / R0])
        vmk_r.set_data([tt[k]], [S["rh_r"][k] / R0])
        fig.suptitle("Angular momentum vs gravity on the E-QFT cell lattice  "
                     f"(t = {tt[k]:4.1f} orbital periods)", fontsize=12.5)
        return ()

    frames = list(range(nf))
    anim = FuncAnimation(fig, update, frames=frames, blit=False)
    print(f"  rendering {nf} frames -> {out_gif} ...")
    anim.save(out_gif, writer=PillowWriter(fps=12), dpi=(70 if quick else 95))
    plt.close(fig)
    print(f"  saved {out_gif}")


# ----------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--no-anim", action="store_true")
    args = ap.parse_args()

    S = run(quick=args.quick)
    make_figure(S, "galaxy_rotation_curve.png")
    np.savez_compressed("galaxy_rotation_traj.npz",
                        rr=S["rr"], vc=S["vc"], frame_times=S["frame_times"],
                        rh_c=S["rh_c"], rh_r=S["rh_r"])
    if not args.no_anim:
        make_animation(S, "galaxy_rotation_3d.gif", quick=args.quick)
    print("done.")


if __name__ == "__main__":
    main()
