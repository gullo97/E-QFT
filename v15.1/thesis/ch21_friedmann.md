# Chapter 21 — Closing the constraint: Friedmann universes and emergent attraction

---

Chapter 20 put the model on the right *orbit family* — but an orbit family is a kinematic statement. Which orbit a universe takes — how fast it expands for a given matter content, whether it sits on the vacuum circle or in the matter interior — is fixed in general relativity by the **Hamiltonian constraint**: the statement that the total energy of geometry plus matter vanishes, the geometry's kinetic energy entering with a negative sign (Ch. 16.3). This chapter imposes that constraint on the generator-coupled ladder — imported as structural **[Postulate]**s (the balance itself and, for inhomogeneous chains, its per-cell curvature extension), their microscopic derivation deferred with a precise attack plan (Ch. 27, item 1) — and shows that the quantum ladder then becomes a *Friedmann universe* to sub-percent accuracy: expansion decelerating because matter dilutes, at the rate matter dictates. The same closure, applied cell by cell to an inhomogeneous chain, then produces the thing Part I never had: **attraction** — overdensities growing at exactly the Newtonian rate, turning around, and collapsing at the classic spherical-collapse threshold. Matter pulling matter together, inside the model — with the honest attribution stated up front: the per-cell acceleration law that the imported constraint structure dictates *is* the force law; what the model adds is everything else.

## 21.1 What the constraint adds — and what it costs

For isotropic FRW (Ch. 16.3), GR's constraint reads

$$\mathcal H_{\text{grav}} + \mathcal H_{\text{matter}} = 0, \qquad \mathcal H_{\text{grav}} = -\,3\kappa\,V H^2, \quad \kappa = \frac{1}{8\pi G}, \tag{21.1}$$

equivalent to Friedmann's equation $3\kappa H^2 = \rho$. Five ingredients, separately accounted:

- the **sign of the form** (geometry's energy negative along the volume direction): *measured* — the compression stiffness is negative (Ch. 22, $\kappa_{\text{grad}}^{\text{vac}} < 0$, the regulator-independent 1D anchor; the 3+1D sign audit is Ch. 25.5), so this piece is the model's own;
- the **exact-balance statement** (total energy exactly zero, and this specific quadratic form): imported **[Postulate]** — the balance is not derived from the ladder, and Ch. 27 item 1's conservation test is its named discharge;
- the **dynamics** it governs: derived (the ladder of Ch. 20, $\dot\alpha = 2g$);
- the **stiffness** $\kappa$ setting the units: computable from matter loops (Ch. 22), used here as given;
- the **matter side itself**: a continuum perfect fluid with dilution law $\rho = \rho_0 e^{-3(1+w)(\alpha - \alpha_0)}$ and a *free* equation-of-state parameter $w$ — an external graft **[Postulate]**, and the least advertised of the five, so it is stated loudly: nothing in this chapter connects $w$ to the model's actual matter content. A mode-$n$ occupant has a perfectly definite, derivable energy law (per cell, $E \propto n^2/L^2$ in the Part I spectrum), i.e. the model *owns* an equation of state that is never computed here; deriving it and checking which $w$ the cells actually realize belongs to the same conservation test (Ch. 27, item 1);
- the **per-cell offset** (spatial curvature), needed from §21.3 on: a perturbed cell must be allowed to miss the flat balance by a constant — the $-k/a^2$ term of a closed patch. This is not a nicety but the whole ballgame for structure: under the strict flat closure (21.2) every cell has $\dot\alpha = 2g \ge 0$ and *can never turn around*; growth, turnaround and collapse live entirely in this offset. Imported **[Postulate]** alongside the balance, same attack plan (Ch. 27, item 1).

(That is six bullets for five ingredients because the old "form" bullet is now honestly split: its *sign* is measured, its *exact balance* is postulated.) The honest framing matters for the defense: this chapter does not derive Friedmann from nothing — it shows that the standard structural statements above, grafted onto the model's own derived dynamics, yield the *entire* quantitative phenomenology of homogeneous and (at leading order) inhomogeneous cosmology. The graft is finite and itemized — dynamics and sign are the model's; balance, fluid, and offset are imports — and the yield is large; Ch. 27 narrows what remains to a finite conservation test.

## 21.2 Constraint closure on the ladder

The generator-coupled ladder expands at $\dot\alpha = 2g$ (Ch. 20, the fixed-point law) — with *constant* $g$, a de Sitter-like exponential: nothing yet tells the coupling about matter. The constraint is precisely that missing instruction. Impose $3\kappa\,\dot\alpha^2 = \rho$ with matter diluting as $\rho = \rho_0\,e^{-3(1+w)(\alpha - \alpha_0)}$:

$$3\kappa\,(2g)^2 = \rho(\alpha) \qquad\Longrightarrow\qquad g(\alpha) \;=\; g_0\,e^{-\frac{3(1+w)}{2}(\alpha - \alpha_0)}, \qquad g_0 = \tfrac12\sqrt{\rho_0/3\kappa}. \tag{21.2}$$

The coupling is **slaved to the local matter density** — the ladder's hopping becomes a self-consistent functional of what the cells contain. The ladder's own equation of motion $\dot\alpha = 2g(\alpha)$ then integrates in closed form:

$$\alpha(t) = \alpha_0 + \frac{2}{3(1+w)}\ln\!\Big(1 + \tfrac{3(1+w)}{2}H_0 t\Big) \qquad\Longleftrightarrow\qquad a(t) \propto t^{\frac{2}{3(1+w)}}, \quad H_0 = 2g_0, \tag{21.3}$$

— exactly the Friedmann power laws: $a \propto t^{2/3}$ for matter ($w = 0$), $t^{1/2}$ for radiation ($w = \tfrac13$). **[Theorem, given (21.1)]**

**The quantum verification.** A genuine wavepacket on a 3600-rung ladder, hopping $t(n) = g\,(n + \tfrac12)\sqrt{n/(n+1)}$ with $g$ updated self-consistently from the constraint at every step (240 updates across the run): **[Computed]** (`ch21_friedmann.py`) maximum deviation $|\alpha_{\text{quantum}} - \alpha_{\text{Friedmann}}| = 0.009$ over $\Delta\alpha \approx 2.0$ for matter and $0.009$ over $1.6$ for radiation — sub-percent agreement across two e-folds — with the local exponent $d\ln a/d\ln t$ converging to $2/3$ and $1/2$ respectively. The expanding quantum box, constraint-closed, *is* a Friedmann universe.

![Friedmann from the ladder](figures/ch21_fig2_friedmann.png)

*Figure 21.2 — Quantum ladder vs analytic Friedmann for matter and radiation: $\alpha(t)$ overlay (deviation $\le 0.009$) and the running exponent locking onto $2/3$ and $1/2$.*

![The energy budget](figures/ch21_fig1_constraint_budget.png)

*Figure 21.1 — The constraint as bookkeeping. Geometry's (negative) kinetic energy balancing matter's (positive) density at every instant; the coupling $g(\alpha)$ as the messenger between them.*

## 21.3 Inhomogeneity at leading order: the separate-universe chain

The model's universe is a *lattice* of cells, each carrying its own geometry — which makes it natively suited to the leading approximation of inhomogeneous cosmology, the **separate-universe picture**: each region evolves as its own FRW patch with its own density, gradients entering at next order (this is exact in the long-wavelength limit and standard in structure-formation theory; Toolbox below). The honest first inhomogeneous simulation is therefore a chain of cells, each evolved as its own FRW patch under the *curvature-extended* closure of §21.1 — not the flat law (21.2), which forbids turnaround — with an initial overdensity across the middle.

> **Toolbox: the separate-universe approximation.** For perturbations of wavelength far exceeding all causal kernels, a perturbed region cannot be distinguished, locally, from an unperturbed universe with slightly different parameters; evolving each region with its own Friedmann equations *is* perturbation theory at leading order in gradients. The model's cell structure implements this without further approximation — the inter-cell coupling that generates the next order is exactly the gradient stiffness measured in Ch. 22.

> **Cell-Growth Theorem [Theorem].** Let the background obey $\ddot a = -\tfrac{4\pi G}{3}\bar\rho\,a$ and a perturbed cell $a_p = a(1 - \lambda)$, $\lambda \ll 1$. Mass conservation per comoving volume gives the density contrast $\delta = 3\lambda$. Subtracting the two acceleration equations and linearizing:
>
> $$\ddot\lambda + 2H\dot\lambda = \tfrac{4\pi G}{3}\bar\rho\,\delta \qquad\Longrightarrow\qquad \boxed{\;\ddot\delta + 2H\dot\delta - 4\pi G\,\bar\rho\,\delta = 0\;} \tag{21.4}$$
>
> — the textbook growth equation of cosmological structure formation, obtained by comparing neighbouring cells *given* that each obeys the per-cell acceleration law. The attribution deserves one plain sentence: that law is exactly the imported force law — the comparison supplies the *structure* of the growth equation; the postulate supplies the $4\pi G$. In matter domination ($a \propto t^{2/3}$, $4\pi G\bar\rho = \tfrac{2}{3t^2}$) the growing mode is
>
> $$\delta \;\propto\; t^{2/3} \;\propto\; a$$
>
> — density contrast grows linearly with the scale factor: matter falls toward overdensities at exactly the Newtonian rate. $\blacksquare$

Beyond linear order, the cell picture keeps paying — and here the per-cell curvature offset of §21.1 earns its keep, since without it $\dot a_p > 0$ always. An overdense patch is a slightly closed universe: its exact parametric solution $a_p \propto (1 - \cos\vartheta)$, $t \propto (\vartheta - \sin\vartheta)$ has the linear-theory shadow

$$\delta_{\text{lin}}(t) \;=\; \frac{3}{5}\Big[\tfrac{3}{4}\big(\vartheta - \sin\vartheta\big)\Big]^{2/3},$$

so **turnaround** ($\dot a_p = 0$, $\vartheta = \pi$) occurs when $\delta_{\text{lin}} = \tfrac35\big(\tfrac{3\pi}{4}\big)^{2/3} = 1.062$ and **recollapse** ($\vartheta = 2\pi$) when $\delta_{\text{lin}} = \tfrac35\big(\tfrac{3\pi}{2}\big)^{2/3} = 1.686$ — the spherical-collapse thresholds underlying all of structure formation, derived in the two lines above rather than quoted from folklore. **[Standard]**

## 21.4 The collapse simulation: watching a potential well form

**[Computed]** (`ch21_collapse.py`): forty-one cells evolved as *classical* per-cell FRW patches, $\ddot a_i = -\tfrac{H_0^2}{2}\,\Omega_i/a_i^2$ with $\Omega_i = 1 + \delta_i(0)$, from synchronized initial data on the per-cell curvature-extended constraint $\dot a_i^2 = H_0^2\big(\Omega_i/a_i + 1 - \Omega_i\big)$ — the $(1 - \Omega_i)$ offset is §21.1's fourth ingredient doing its work. (The quantum ladder stands behind the homogeneous background of §21.2; it does not appear in this chain — the chain tests the constraint structure, not the ladder.) Central contrast $\delta(0) = 0.08$ with a Gaussian profile, starting at $a_0 = 0.02$ in the constraint's units, run over a factor $\sim 2000$ of background expansion:

- the synchronized start is almost entirely decaying mode: linearizing the initial data gives growing amplitude $C_+ = \tfrac35\,a_0\,\delta_0 = 9.6\times10^{-4}$ per unit $a/a_0$, so the thresholds *predict* turnaround at $a_{\text{bg}}/a_0 \approx 1.062/C_+ \approx 1.1\times10^3$ and collapse at $\approx 1.8\times10^3$ — the arc below is checkable, not folkloric;
- the central contrast grows with $d\ln\delta/d\ln a = 1.08$ across the linear window — theory: $1$, the excess being the decaying-mode transient plus incipient nonlinearity;
- it peels upward off the linear law, **turns around** at $a_{\text{bg}}/a_0 = 1168$, and **collapses** at $1863$ — matching the $C_+$-based predictions as the thresholds dictate;
- every cell in the Gaussian tail lags proportionally less: the profile of $a_i/a_{\text{bg}}$ dips smoothly toward the center — **a potential well forming**, deepest where the matter is (the quantitative weak-field identification $\delta\alpha = -\Phi_{\text{Newton}}$ is Ch. 24's first section).

The animation shows the unmistakable picture: in comoving coordinates the cell walls drift toward the overdensity from both sides — *matter attracting matter* — until the central cells crash while the outskirts keep expanding. This is gravitational instability — the mechanism that turned a smooth early universe into galaxies — running on the constraint-closed cell chain. The honest attribution, once more: the microscopic model supplies the cell ontology and the quantum background law (§21.2); the instability arc itself runs on the imported constraint structure, curvature offset included.

![Collapse, live](figures/ch21_anim_collapse.gif)

*Animation 21.A — Emergent attraction in motion: comoving cell walls falling toward the central overdensity beside the growing density contrast $\delta(a)$, through turnaround and collapse.*

![Growth, turnaround, collapse](figures/ch21_fig3_collapse.png)

*Figure 21.3 — The full arc. Left: central contrast vs scale factor — linear growth ($\delta \propto a$), the peel-off, turnaround and collapse markers at the threshold values. Right: the $a_i/a_{\text{bg}}$ profile at successive epochs — the well deepening around the matter.*

## 21.5 Summary and the price tag

Constraint closure converts the Kasner-capable ladder into a quantitative Friedmann cosmology (sub-percent over two e-folds, correct matter and radiation exponents) and, applied cell-wise, yields emergent attraction with the full arc — $\delta \propto a$, turnaround at $1.062$, collapse at $1.686$ — of gravitational instability. The bill for all of it: the imported constraint structure — the balance (21.1) *and* its per-cell curvature offset, without which nothing ever turns around — and one number, the stiffness $\kappa$, which entered only as a unit. The next chapter computes that number from the model's own matter content — sign first, magnitude after — and in doing so converts the constraint's "geometry energy is negative" from postulate to measurement.

---

**Validation.** `ch21_friedmann.py`: the constraint-closed quantum ladder vs (21.3) for $w = 0, \tfrac13$ (Fig. 21.2; deviations $0.0090/0.0091$ and exponents as quoted). `ch21_collapse.py`: the 41-cell chain (Fig. 21.3; growth exponent $1.082$, turnaround $1168$, collapse $1863$, linearly extrapolated contrast $1.17/1.87$ at the two epochs vs the spherical thresholds $1.062/1.686$), trajectory archive `data/ch21_collapse_traj.npz`. (Both adapted from the archived `scripts/ports/newton_sim.py`; `scripts/ports/make_collapse_anim.py` renders the animation.) All quoted numbers printed by the scripts.
