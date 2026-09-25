# Chapter 19 — Shape space: the Splitting Theorem

> **Author: Fable 5 · Complexity ◇4–◇5 · Status: draft v1 (pending global review pass)**
> Depends on: Ch. 16.4, 18. Feeds: Ch. 20 (why diagonal motion is a discrete ladder), Ch. 25 (the sectors whose stiffnesses get measured).

---

Chapter 18 completed the kinematics: six metric directions per cell, three reached by Part I's discrete dilation ladder, three by the continuous shears. That asymmetry of *origin* — jumps versus flows — looks like an accident of how the theory was built. This chapter proves it is not. The two sectors are distinguished by the *quantum geometry of the matter itself*: against the natural metric on shape space (how much the matter's eigenstates change per unit deformation), the diagonal directions are exactly **flat but energy-carrying**, while the shears are exactly **energy-flat but geometry-carrying**, at every rectangular configuration. The split is a theorem, and it lands the model squarely on one of general relativity's deepest structural features: the DeWitt supermetric's segregation of the volume direction (the "clock") from the shape directions (the "propagating" modes). The model does not imitate that split; it *derives* one — with an honest caveat about the Lorentzian sign that later chapters discharge dynamically.

## 19.1 The quantum metric on shape space

How much does a quantum state care about a deformation of its geometry? The canonical answer is the **fidelity (quantum) metric**. Let $H(A)$ depend on parameters $A$ (here: the inverse cell metric), with normalized eigenstate $|\psi_n(A)\rangle$. For a deformation direction $\delta A$,

$$1 - \big|\big\langle \psi_n(A)\,\big|\,\psi_n(A + \varepsilon\,\delta A)\big\rangle\big| \;=\; \frac{\varepsilon^2}{2}\,\chi_n(\delta A) \;+\; \mathcal O(\varepsilon^4), \tag{19.1}$$

and standard first-order perturbation theory gives the equivalent spectral form (derivation: insert $|\psi_n(A + \varepsilon\delta A)\rangle = |n\rangle + \varepsilon\sum_{m\neq n}|m\rangle\frac{\langle m|H_1|n\rangle}{E_n - E_m} + \cdots$ with $H_1 = \partial_\varepsilon H$, and expand the overlap):

$$\chi_n(\delta A) \;=\; \sum_{m \neq n} \frac{\big|\langle m|\,H_1(\delta A)\,|n\rangle\big|^2}{(E_m - E_n)^2}. \tag{19.2}$$

$\chi_n \ge 0$ always (a sum of squares): the fidelity metric is positive-definite *by construction* — remember this for §19.4. Both (19.1) and (19.2) are computed independently in the numerics as a cross-check of implementation, not of mathematics.

**The arena.** Following the dictionary, a particle in a general parallelepiped is a particle on the *fixed* unit square (we work in $d = 2$, where the numerics are clean and nothing conceptual changes in $d = 3$) with constant inverse metric $A = g^{ab}$:

$$H(A) \;=\; -\Big(A^{11}\,\partial_x^2 \;+\; 2A^{12}\,\partial_x\partial_y \;+\; A^{22}\,\partial_y^2\Big), \qquad \psi\big|_{\partial} = 0. \tag{19.3}$$

Deformations of the cell are deformations of $A$: diagonal directions (conformal $\delta A \propto A$, traceless stretch $\delta A \propto \mathrm{diag}(1,-1)$-type) and the shear $\delta A = E_{12} + E_{21}$.

## 19.2 The Splitting Theorem

> **Shape-Space Splitting Theorem [Theorem].** At any rectangular configuration ($A^{12} = 0$) of (19.3):
>
> **(i) Diagonal sector — overlap-flat, energy-carrying.** For every diagonal deformation $\delta A = \mathrm{diag}(\delta A^{11}, \delta A^{22})$ (including the conformal direction):
>
> $$\chi_n(\delta A_{\text{diag}}) \;=\; 0 \quad \text{exactly}, \qquad \frac{dE_{mn}}{d\varepsilon} \;=\; \pi^2\big(m^2\,\delta A^{11} + n^2\,\delta A^{22}\big) \;\neq\; 0 . $$
>
> **(ii) Shear sector — iso-energy, geometry-carrying.** For the shear deformation $\delta A = E_{12} + E_{21}$:
>
> $$\frac{dE_{mn}}{d\varepsilon} \;=\; 0 \quad \text{exactly}, \qquad \chi_{mn}(\delta A_{\text{shear}}) \;>\; 0 \;\;\text{finite}.$$

*Proof.* **(i)** At $A^{12} = 0$ the eigenfunctions of (19.3) are the product modes $\psi_{mn}(x, y) = 2\sin(m\pi x)\sin(n\pi y)$ — and these *contain no reference to $A^{11}, A^{22}$ whatsoever*: the diagonal entries multiply the separated second-derivative terms and enter only the eigenvalues $E_{mn} = \pi^2(m^2 A^{11} + n^2 A^{22})$. A deformation that leaves every eigenfunction identically unchanged has $\chi_n = 0$ by (19.1) directly (the overlap is exactly 1 at every $\varepsilon$ within the rectangular family); the energy slope is read off $E_{mn}$. **(ii)** First-order energy shift under the shear: $\langle\psi_{mn}|(-2\partial_x\partial_y)|\psi_{mn}\rangle \propto \int_0^1\sin(m\pi x)\cos(m\pi x)\,dx \cdot \int_0^1 \sin(n\pi y)\cos(n\pi y)\,dy = 0$ — each factor integrates $\tfrac12\sin(2\pi k u)$ over a full period. The shear is *automatically* iso-energy at rectangular points: no quantum-number jump is needed to stay on the energy level set. Its fidelity metric, by contrast, is generically nonzero: $\langle m'n'|(-2\partial_x\partial_y)|mn\rangle \ne 0$ for $(m', n')$ of opposite parities in both factors, and (19.2) sums their squares with no cancellation (every term $\geq 0$). $\blacksquare$

**[Computed]** (`ch19_fidelity.py`, base configuration $A = \mathrm{diag}(1, 1.21)$): diagonal-sector fidelity at machine zero — $\chi = 1.8\times10^{-13}$ (conformal direction), $\lesssim 10^{-13}$ (traceless stretch), across all modes tested; shear-sector fidelity finite and mode-structured — $\chi = 0.049,\ 5.32,\ 5.32,\ 3.73$ for the first four modes — agreeing with the perturbation-theory form (19.2) to four significant digits; shear energy slopes consistent with zero at first order (residuals $\sim \varepsilon\,\times$ curvature).

![Fidelity vs perturbation theory](figures/ch19_fig2_fidelity.png)

*Figure 19.2 — The split, measured. Left: fidelity metric from direct overlaps (19.1) vs the perturbative form (19.2) — four-digit agreement. Right: $\chi$ by deformation direction and mode — diagonal directions at machine zero, shear directions finite.*

## 19.3 What the split means

Put the two halves side by side with Part I and with general relativity.

**Why diagonal motion is a discrete ladder.** Half (i) says diagonal deformations change energies *without touching wavefunctions*. Iso-energy motion in the diagonal sector therefore cannot be continuous — to hold $E$ fixed while the spectrum slides, the state must **jump** quantum numbers. That is precisely Part I's dilation ladder $\hat T^{(a)}_\pm$, with its $\sqrt{n/(n+1)}$ overlap tolls: not a modeling choice after all, but the *only* way the diagonal sector admits iso-energy motion. Conversely, half (ii) says the shear sector needs no ladder: it flows continuously along the energy level set, accumulating quantum geometry ($\chi > 0$) as it goes. The two sectors of Ch. 18's audit are dynamically heterogeneous because the *matter* makes them so.

**The DeWitt parallel.** In canonical GR (Ch. 16.4), the kinetic term $\tfrac{1}{2N}G^{abcd}\dot g_{ab}\dot g_{cd}$ built on the DeWitt supermetric segregates the conformal direction (opposite signature; the "clock" against which the rest evolves) from the five shape directions (the propagating sector — gravitational waves live here). The box model has now *derived* a split with the same membership: the volume/diagonal sector is the special, clock-like one — first-order energy flow, discrete, and (by Ch. 6) irreversibly directed: the expansion arrow itself is the ticking — while the shear sector is the smooth, propagating one. The membership matches; the roles match; the arrow even explains *why* the clock runs one way.

**The honest caveat — and its promissory note.** The DeWitt metric's conformal direction carries a *negative sign*; the fidelity metric (19.2) is positive-definite by construction, so the Lorentzian signature **cannot** be kinematic here — at this stage the conformal direction is distinguished by flatness ($\chi = 0$) and directedness, not by sign. The sign must come from *dynamics*: the energy functional of geometry, not the geometry of states. The thesis discharges this note twice, quantitatively: the induced compression stiffness is measured **negative** in 1D (Ch. 22, $\kappa^{\text{vac}}_{\text{grad}} \approx -0.13$) and the full DeWitt pattern — TT positive against conformal negative, with the universal coefficient — is measured in 3D (Ch. 25). Flagging the gap here and closing it there is the difference between an analogy and a derivation, and the thesis insists on the second.

## 19.4 The identification-map subtlety, resolved

Half (i) contains an apparent paradox worth defusing in print, because it is the same trap as Ch. 10's, wearing geometric clothes. *"Diagonal deformations leave eigenfunctions unchanged"* (this chapter) versus *"expansion redistributes the state over new modes with $\sqrt{n/(n+1)}$ overlaps"* (Ch. 3) — both about the same physical stretch. Which is true? **Both**, relative to different identification maps between the Hilbert spaces before and after:

- the **fixed-coordinates (metric) picture**: the domain stays $[0,1]^2$, the metric changes; states are compared as functions on the fixed domain — diagonal stretches then change nothing about the functions (this chapter);
- the **embedding (moving-walls) picture**: the domain physically grows; states are compared through extension-by-zero on physical space — overlaps are then the $\sqrt{n/(n+m)}$ structure of Ch. 3.

Neither is "the" answer; *declaring the map is part of the physics* (the principle minted in Ch. 10). The deformation-versus-state bookkeeping freedom is, in fact, an old acquaintance in new dress: it is the lattice ancestor of the **lapse-and-shift freedom** of canonical gravity — how much of "the geometry changed" is attributed to coordinates versus to physical motion. The two natural maps bound the physics from two sides, and the dynamical chapters are scrupulous about which is in force: Ch. 20's ladder dynamics lives in the embedding picture (transitions are physical traffic); Ch. 25's static response lives in the metric picture (deformations of a background). Their agreement on overlapping questions is one more overdetermination check the framework passes.

## 19.5 Summary

The fidelity metric splits shape space exactly: **diagonal = energy-carrying, overlap-flat, hence discrete, clock-like, arrow-equipped; shear = iso-energy, geometry-carrying, hence continuous and propagating** — the DeWitt division, derived from matter rather than postulated, with the Lorentzian sign honestly deferred to its dynamical origin (Ch. 22, 25). With kinematics (Ch. 18) and quantum geometry (this chapter) in place, one ingredient remains before the model can be run as a universe: the *dynamics* — which operator couples neighbouring geometry sectors, and with what hopping law. That choice, Chapter 20 proves, is not a detail: it selects between a wrong universality class and exactly general relativity's.

---

**Validation.** `ch19_fidelity.py` (port): direct-overlap vs perturbative fidelity (four digits); the machine-zero diagonal sector and finite shear sector values quoted above; shear energy-slope nulls. Figure 19.1 (shape-space schematic: `ch19_fig1_shape_space.png`, drawn) and Fig. 19.2 (data).
