# Chapter 24 — Newton's law and BKL chaos on the cell lattice

---

With the stiffness computed (Ch. 22) and the constraint closed (Ch. 21), the model owes us the two most recognizable behaviours of gravity it has not yet exhibited: the **inverse-square law** — the phenomenon everyone means by "gravity" — and the **chaotic approach to a singularity** that general relativity predicts in the strong field (BKL), arguably the deepest known statement about classical GR. This chapter collects both debts on the same cell lattice that did Kasner and Friedmann.

## 24.1 The lattice Poisson equation

Linearize the constraint-closed cell dynamics about a static uniform background: cell $\mathbf x$ carries $\alpha(\mathbf x) = \bar\alpha + \delta\alpha(\mathbf x)$ and density $\bar\rho(1 + \delta(\mathbf x))$ ($\delta$ = density contrast throughout this section). Two derived objects enter: the induced gradient stiffness of the geometry field (Ch. 22) and the matter coupling from the constraint (Ch. 21). A word on which stiffness, since its sign was the point of Ch. 22: $\delta\alpha$ is the *conformal* direction, whose induced energy is *negative* — an energy-minimization reading of the static sector would be unbounded and meaningless. That is not a defect but the GR structure itself: the conformal mode is not a dynamical direction to be minimized over; it is *eliminated by the Hamiltonian constraint* (Ch. 16.4, 21), and in canonical GR the weak-field limit of that constraint — not any energy principle — is what produces $\nabla^2\Phi = 4\pi G\rho$. The same elimination here, with $\kappa \equiv |\kappa_{\text{grad}}|$ the magnitude of the induced stiffness setting the units, gives the static sector as

> **Lattice Poisson Result [Theorem, given Ch. 21–22].**
>
> $$\kappa\,\hat\nabla^2\,\delta\alpha(\mathbf x) \;=\; -\,c_m\,\bar\rho\,\delta(\mathbf x), \qquad \delta\alpha \;=\; -\,\Phi_{\text{Newton}}, \tag{24.1}$$
>
> with $\hat\nabla^2$ the discrete Laplacian: **the local log-cell-size is (minus) the Newtonian potential.** Cells are *larger* near mass — which is precisely the weak-field spatial metric of GR, $g_{ij} = (1 - 2\Phi)\delta_{ij}$ with $\Phi < 0$ near a source (Ch. 16.1), read through the dictionary ($g \sim L^2$, $\delta\ln g = 2\,\delta\alpha = -2\Phi$). For two sources, the interaction energy is $E_{\text{int}}(d) = \int \Phi_1\rho_2 \propto -1/d$ (negative — binding — since $\Phi_1 < 0$ where $\rho_2$ sits), whose gradient is the inverse-square force. (The matter coupling is written $c_m$ — not $c$, which this Part already uses for the lattice constant $\Lambda = c/\bar L$; App. D.5.)

*Derivation sketch (full steps in the script's docstring and App. C).* Linearize the per-cell Hamiltonian constraint (geometry gradient energy + matter energy = 0, Ch. 21) in $\delta\alpha$ and $\delta$; the static sector is the finite-difference relation $\kappa\,\hat\nabla^2\delta\alpha = -c_m\,\bar\rho\,\delta$, i.e. (24.1) — a constraint solved for $\delta\alpha$, not a functional minimized. The *sign* of the outcome — cells dilate where matter sits, $\delta\alpha = -\Phi > 0$ near a source, attraction — traces to the measured negative compression stiffness (Ch. 22.2, regulator-independent; the 3+1D sign audit is Ch. 25.5), and the identification $\delta\alpha = -\Phi$ then carries the standard weak-field reading. $\blacksquare$

## 24.2 The measurement: $1/r$, $1/d$, $1/d^2$

**[Computed]** (`ch24_newton_lattice.py`; $128^3$ periodic cell lattice, FFT solution of (24.1); the finite-volume zero-mode constant fitted and removed — standard periodic-box bookkeeping):

| observable | measured exponent / value | Newton |
|---|---|---|
| point-source profile $\delta\alpha(r)$ | $r^{-1.008}$ | $r^{-1}$ |
| profile amplitude | $0.0766$ | $1/4\pi = 0.0796$ ($96\%$; lattice discreteness at small $r$) |
| two-mass binding energy | $d^{-0.994}$ | $d^{-1}$ |
| force | $d^{-2.05}$ | $d^{-2}$ |

**Newton's inverse-square law, emergent on the lattice of cells** — the same lattice whose homogeneous modes did Kasner (Ch. 20) and Friedmann (Ch. 21). The chain is complete at leading order with no link assumed: matter loops give the stiffness (Ch. 22), stiffness plus constraint give Poisson (24.1), Poisson gives $1/r^2$.

![Newton on the lattice](figures/ch24_fig1_newton.png)

*Figure 24.1 — The $1/r$ potential profile (log–log, fitted exponent $-1.008$), the $1/d$ binding energy, and the $1/d^2$ force on the $128^3$ lattice.*

![Cells dilated around a mass](figures/ch24_fig2_well_schematic.png)

*Figure 24.2 — The weak-field metric, cell-model edition: cells dilated near the source ($\delta\alpha = -\Phi > 0$), the potential well of Ch. 21's collapse now static and quantitative.*

### 24.2.1 Are the residuals physics? What is numerical, what is falsifiable

A careful reader will ask the right question about the table above: the measured exponent is $-1.008$, not $-1$; the amplitude is $96\%$, not $100\%$; the force runs as $d^{-2.05}$. *Is the model predicting a measurable deviation from Newton?* The answer must be given in two sharply separated parts, because conflating them would be the worst kind of dishonesty — manufacturing a "testable prediction" out of a fit residual.

**Part one: the residuals in the table are numerical artifacts, and provably so.** Equation (24.1) is exactly the Poisson equation with a discrete Laplacian; its continuum limit is exactly $1/r$. The deviations in the table have three mundane, checkable sources: (i) *periodic images* — the $128^3$ torus superimposes copies of the source, and the zero-mode subtraction leaves a residual quadrupole of the image lattice in any finite fitting window; (ii) *fit-window contamination* — the log–log fit includes radii within a few lattice spacings of the source, where the discrete Green's function genuinely differs from $1/4\pi r$; (iii) *cubic anisotropy* of $\hat\nabla^2$, which enters the lattice Green's function only at relative order $(a/r)^2$ with a quadrupolar angular structure **[Standard]** (the standard large-$r$ expansion of the cubic-lattice Green's function). All three shrink under lattice refinement at fixed physical window; none survives the continuum limit. The honest statement of the table is: *the lattice measurement is consistent with Newton to the accuracy the box permits* — the $0.8\%$ on the exponent is the quality of the simulation, not a property of nature. Anyone quoting "$r^{-1.008}$" as this model's predicted gravity law has mistaken error bars for physics, and this section exists so that no reader — friendly or hostile — can make that mistake in either direction.

**Part two: the model does predict real deviations — elsewhere, and each with its scale attached.** Because the continuum limit is exact, deviations from GR+Newton live at the *edges* of the construction, and the framework is specific about where:

1. **Planck-suppressed corrections to Newton.** The genuine lattice correction to the potential is the $(a/r)^2$ term above with $a = \bar L \approx 1.6\,c\,\ell_P$ (Ch. 22): fractional deviations of order $(\ell_P/r)^2 \sim 10^{-70}$ at micron scales. Real, absolutely predicted, and hopelessly unobservable — stated so that no one wastes a torsion balance on it.
2. **Gravitational-wave birefringence at the cell scale.** The two TT polarizations agree at $2\%$ at accessible momenta (Ch. 25) because rotational invariance is *emergent*; at wavelengths approaching $\bar L$ the cubic cell structure must split $+$ from $\times$. This is the model's cleanest structural difference from fundamental GR — a definite prediction that Lorentz/rotation invariance fails at the cell scale — but the frequencies are Planckian. Its accessible shadow is item 3.
3. **Modified dispersion at $\mathcal O(\omega\bar L)^n$.** Granular geometry generically modifies propagation, $\omega^2 = k^2[1 + \xi(k\bar L)^n + \cdots]$. Time-of-flight limits (GRB photons; GW170817's $|c_{gw}/c_\gamma - 1| < 10^{-15}$) already exclude $n = 1$ with $\xi = \mathcal O(1)$ — so the model is *falsifiable today* on this axis: it survives only if its dispersion starts at $n = 2$ (as naive $(a/r)^2$-type lattice corrections suggest) or if the coefficient is dynamically suppressed. Computing $n$ and $\xi$ from the cell dynamics is Ch. 27, item 10 — a finite calculation with a kill condition attached.
4. **The massive Weyl vector.** Gauge-fixed Weyl theory leaves a vector of mass $M_W = \sqrt6\,e_W M_P$ (Ch. 23): a Yukawa-range fifth force at $\sim\ell_P$ — invisible — *unless* $e_W$ is small, in which case fifth-force searches bound $e_W M_P$ from below. The model prefers but does not yet fix $e_W = \mathcal O(1)$; stated as the constraint it is.
5. **The baryogenesis threshold.** The sharpest near-term falsifiable structure in the thesis is not gravitational: Ch. 14's threshold — perturbatively narrow wall-angle ensembles pump *exactly zero* — plus the required $\varepsilon_{CP} \sim 10^{-8}$–$10^{-7}$ window make the boundary-GIM calculation (Ch. 27, item 4) a genuine kill test. If boundary CP violation is GIM-crushed like bulk CP violation, the framework's $\eta_B$ is zero and the mechanism is dead; no gravitational observation is needed.

The pattern is deliberate and worth stating: *where the model overlaps GR and standard cosmology it is constructed to agree exactly in the continuum limit, and the numerics above verify that construction; where it must differ — granularity signatures and the baryogenesis engine — the differences are either Planck-suppressed with known exponents or concentrated in one computable CP quantity.* The model's falsifiability is real but lives at these joints, not in the third decimal of a lattice fit.

## 24.3 Bianchi IX: Kasner flights between curvature walls

The second debt is strong-field. General relativity's deepest classical prediction concerns the approach to a generic spacelike singularity, and it is *not* the smooth crunch naive extrapolation suggests. The arena is **Bianchi IX** — the anisotropic universe with the 3-sphere's homogeneous curvature (Ch. 16.7) — whose dynamics in BKL time $\tau$ ($dt = a_1a_2a_3\,d\tau$, $\alpha_i = \ln a_i$) reads

$$2\,\alpha_1'' \;=\; \big(a_2^2 - a_3^2\big)^2 \;-\; a_1^4 \quad (+\ \text{cyclic}), \qquad \sum_{i<j}\alpha_i'\alpha_j' \;=\; \frac14\Big[\sum_i a_i^4 - 2\sum_{i<j}a_i^2 a_j^2\Big], \tag{24.2}$$

evolution plus constraint. The mechanics is a two-act loop. **Act one:** while the right-hand "wall" terms are negligible, (24.2) is free motion — straight lines in $\alpha$-space, *exactly the Kasner flights the generator-coupled ladder produces* (Ch. 20): the model owns this act outright. **Act two:** every Kasner flight carries one negative-exponent axis, and toward the singularity that axis's scale factor **grows** ($a_1 = t^{p_1}$ with $p_1 < 0$ as $t \to 0$); its quartic wall term $a_1^4$ ignites, the trajectory reflects — a **bounce** — and the universe exits onto a new Kasner flight with the axes' roles reshuffled. The reflection law is exact (the bounce is an integrable Bianchi II interlude), and in the $u$-parametrization of the Kasner circle (16.14) it is the **BKL map**:

$$u \;\to\; u - 1 \quad (u \ge 2), \qquad u \;\to\; \frac{1}{u - 1} \quad (1 < u < 2). \tag{24.3}$$

Iterated, (24.3) is conjugate to the Gauss continued-fraction map — ergodic, mixing, chaotic: the collapse is an unending, never-repeating tumble of epochs. The honesty sentence, as the production spec mandates: *the epochs are the model's own derived dynamics (Ch. 20); the wall terms' Bianchi-IX functional form is fixed by the $S^3$ structure constants and is the one imported ingredient of this section* — its derivation from the model's inter-cell coupling is posed as Ch. 27, item 7.

## 24.4 The simulation: five epochs, four decimal places

**[Computed]** (`ch24_mixmaster.py`, re-run verified in this build): constraint-consistent initial data deep in a Kasner epoch at $u_0 = 4.3$; the constraint residual is $7\times10^{-18}$ at $\tau = 0$ and drifts to only $9\times10^{-13}$ across the entire integration — tenth-digit conservation, the run's quality certificate — through five epochs toward the singularity. The measured Kasner parameter per epoch, against the map (24.3):

$$4.3001 \;\longrightarrow\; 3.3000 \;\longrightarrow\; 2.3000 \;\longrightarrow\; 1.3000 \;\longrightarrow\; 3.3333,$$

each step agreeing with the BKL prediction to $\le 6\times10^{-4}$ — **including the final step, the era reversal** $u = 1.3 \to 1/(1.3 - 1) = 10/3 = 3.333\ldots$, the continued-fraction move that makes the dynamics chaotic rather than a countdown. The figure shows the signature picture: piecewise-straight $\alpha_i(\tau)$ flights with the contracting role hopping between axes at each bounce, beside the exponent point leaping around the Kasner circle; the animation plays it:

![Mixmaster, live](figures/ch24_anim_mixmaster.gif)

*Animation 24.A — The tumble toward the singularity: Kasner flights, wall bounces, and the exponent point hopping around the circle by the BKL map — era reversal included.* A collapsing cell-model universe does not approach its singularity smoothly; it tumbles, chaotically, exactly as Belinskii, Khalatnikov and Lifshitz said a general-relativistic universe must.

![Mixmaster](figures/ch24_fig3_mixmaster.png)

*Figure 24.3 — Five Kasner epochs and four bounces. Top: the $\alpha_i(\tau)$ flights with role-swapping reflections. Bottom: the measured $u$-sequence against the BKL map, era reversal included.*

## 24.5 Summary

Two debts paid on one lattice: the static weak field prints the Newtonian $1/r$ potential with $96\%$ of the continuum amplitude and an inverse-square force (exponents $-1.008$, $-2.05$ — simulation residuals, not predicted deviations; §24.2.1 separates the numerical from the falsifiable), and the strong-field collapse executes BKL chaos to $6\times10^{-4}$ through an era reversal, with the constraint conserved to the tenth digit. Between them, these close the classical phenomenology available to the model short of radiation — and radiation is exactly what remains: the propagating, transverse-traceless sector, where the model must either produce a graviton with the right kinetic structure or be falsified at general relativity's door. That computation is next.

---

**Validation.** `ch24_newton_lattice.py`: the $128^3$ Poisson solution, profile/binding/force fits, regenerates Fig. 24.1. `ch24_mixmaster.py`: the five-epoch run, $u$-sequence table, constraint drift with the hard $10^{-10}$ discard gate; writes `data/ch24_mixmaster_traj.npz`, regenerates Fig. 24.3. (Both adapted from the archived `scripts/ports/` originals; animation reused from the validated archive, rendered by `scripts/ports/make_mixmaster_anim.py` from the saved trajectory.) All quoted numbers printed by the scripts.
