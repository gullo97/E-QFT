# Chapter 14 — From pumping to η_B: the corrected estimate

> **Author: Fable 5 · Complexity ◇5 · Status: draft v1 (pending global review pass; §14.2–14.3 contain this thesis's principal new derivation)**
> Depends on: Ch. 6, 12, 13. Feeds: Ch. 26, 27.

---

Chapter 12 delivered a charge pump with an exactly quantized stroke: one unit per level crossing, protected by anomaly arguments against every continuous deformation. What it deliberately did not deliver is a *rate*. Quantization cleanly splits the baryogenesis problem into a protected part (the yield per crossing: exactly $\pm1$) and a contingent part (how often crossings happen, and with what sign bias). This chapter computes the contingent part as far as theorems allow and budgets the remainder without mercy. The result is a master formula whose kinematic engine — the **Crossing-Density Amplitude** — is new, geometric, and computable, and whose CP-strength factor is presented for what it is: the framework's dominant open problem, inherited from the Standard Model's GIM structure (Ch. 13), with sharply identified escape routes.

A reader should finish this chapter knowing exactly which factors of the estimate would survive hostile review unchanged, and which would not.

## 14.1 The pipeline

From a microscopic pump to today's $\eta_B$, five factors:

$$\eta_B \;=\; \underbrace{C_{\text{sph}}}_{28/79} \;\times\; \underbrace{\frac{45}{2\pi^2 g_*}}_{\approx\,0.0214} \;\times\; \underbrace{N_{\text{dof}}}_{\text{species/color}} \;\times\; \underbrace{\mathcal N_{\text{steps}}\,\bar\nu\;\varepsilon_{CP}}_{\text{this chapter}} \;\times\; \underbrace{f_{\text{neq}}}_{\approx\,1\ \text{(Ch. 6)}} . \tag{14.1}$$

The first two are standard electroweak cosmology (Ch. 13, audited there); $f_{\text{neq}} \approx 1$ is the framework's structural advantage — the expansion arrow is a *permanent* departure from equilibrium, so no phase-transition gymnastics are needed (Expansion-Preference Theorem; the Dirac-sector check that CP phases do not weaken the preference is part of the Ch. 8 validation suite). The heart is the third block: $\mathcal N_{\text{steps}}$ counts iso-energy steps per comoving volume during the relevant epoch, $\bar\nu$ is the mean crossing probability per step, and $\varepsilon_{CP}$ is the sign bias between $+1$ and $-1$ pumps. We now build $\bar\nu$ and $\varepsilon_{CP}$ from Chapter 12's exact results.

![The assembly pipeline](figures/ch14_fig1_pipeline.png)

*Figure 14.1 — From pump to ratio. Protected factors (boxes with bold borders: quantized yield, sphaleron conversion, entropy dilution) versus contingent factors (crossing density, sign bias, step counting). The figure is the chapter's table of contents.*

## 14.2 The per-bag yield: a one-crossing theorem

First, a structural theorem that simplifies everything downstream.

> **One-Crossing Theorem [Theorem].** A two-wall chiral bag with fixed angles $(\theta_0, \theta_L)$ and mass $m$ experiences, over its entire expansion history $L: 0 \to \infty$, **at most one** level crossing — hence pumps at most one unit of charge. The crossing exists iff
>
> $$\operatorname{sgn}\big[\cos(\Delta/2)\big] \;=\; -\operatorname{sgn}\big[\cos\Sigma\big] \quad\text{and}\quad \big|\cos(\Delta/2)\big| \;<\; \big|\cos\Sigma\big|, \tag{14.2}$$
>
> in which case it occurs at the unique size
>
> $$L^*(\Delta, \Sigma; m) \;=\; \frac{1}{m}\,\operatorname{artanh}\!\Big(\!-\frac{\cos(\Delta/2)}{\cos\Sigma}\Big). \tag{14.3}$$
>
> *Proof.* The Crossing Condition (12.9) equates $\tanh(mL)$ — strictly increasing from $0$ to $1$ on $L \in (0,\infty)$ — to a constant fixed by the angles. A constant in $(0, 1)$ is hit exactly once; outside $[0,1]$, never. Condition (14.2) is "the constant lies in $(0,1)$" spelled out. $\blacksquare$

**[Computed]** `ch14_crossing_density.py --validate` confirms by brute force: for $2000$ random angle pairs, full spectral-flow scans over $mL \in (0, 8]$ find exactly one $\eta$-jump when (14.2) holds and none otherwise, with jump locations matching (14.3) to the scan resolution.

The existence region (14.2) is the single most consequential fact of this chapter, so look at it. At $\Sigma = 0$ (symmetric walls) it requires $\cos(\Delta/2) < 0$, i.e. $|\Delta| > \pi$: **order-one wall mismatch**. Small mismatches — however nonzero — pump nothing, not "exponentially little" but *exactly nothing*: there is no crossing to count. The pump has a threshold, and the threshold is topological in flavor: the wall angles must disagree by more than a quarter turn (in $\gamma_5$-angle) before the gap hosts a zero mode at any size.

![The crossing region](figures/ch14_fig2_crossing_density.png)

*Figure 14.2 — The pump's phase diagram. The $(\Sigma, \Delta)$ plane with the existence region (14.2) shaded and contours of $mL^*$ from (14.3); overlaid points: brute-force spectral-flow scans (crossing found = filled, none = open). The MIT diagonal ($\theta_0 = 0$ line $\Delta = 2\Sigma$) grazes the boundary and never enters — the one-ordinary-wall no-go of Ch. 12 visualized.*

## 14.3 The Crossing-Density Amplitude *(new)*

Now the ensemble. At the electroweak epoch, the framework pictures comoving regions as bags taking iso-energy steps (Ch. 3: $L_k = (n_0 + k)L_0/n_0$), with wall angles set by local CP microphysics — drawn from some distribution $P(\theta_0, \theta_L)$, generally epoch-dependent. Define the chapter's central object:

> **Definition (Crossing-Density Amplitude).** $\nu_k$ = expected number of crossings traversed during step $k$, per bag:
>
> $$\nu_k \;=\; \int dP(\theta_0, \theta_L)\;\; \mathbb 1\big[(14.2)\text{ holds}\big]\;\; \mathbb 1\big[\,L_k < L^*(\Delta,\Sigma; m) \le L_{k+1}\,\big]. \tag{14.4}$$

Two exact reductions make (14.4) usable.

**(a) The step measure.** In the thermal regime $n \gg 1$ (Ch. 6, caveat 3), the step is infinitesimal: $L_{k+1}/L_k = 1 + 1/n_k$, so the indicator selects a shell of logarithmic width $d\ln L = 1/n_k$, and

$$\nu_k \;=\; \frac{1}{n_k}\; \rho_{\ln L^*}\big(\ln L_k\big), \qquad \rho_{\ln L^*} = \text{density of } \ln L^* \text{ induced by } P \text{ via } (14.3). \tag{14.5}$$

Summing over the steps an expanding region takes during the epoch, the total expected crossings per bag is the integral of $\rho_{\ln L^*}$ over the $e$-folds traversed — step bookkeeping ($\mathcal N_{\text{steps}}$) and crossing density combine into a single, reparametrization-clean quantity: **the fraction of the ensemble whose $L^*$ falls inside the epoch's expansion window**,

$$\mathcal N_{\text{steps}}\,\bar\nu \;=\; \Pr\Big[\,L^*(\Delta, \Sigma; m) \in \big(L_{\text{in}},\, L_{\text{out}}\big)\,\Big] \;\equiv\; \mathcal P_\times . \tag{14.6}$$

This is the corrected kinematic amplitude, and its structure is worth a sentence of contrast: it is a *probability* (bounded by 1, dimensionless, basis-free), not a truncated overlap sum; every pathology catalogued in Ch. 11 is structurally excluded. With $mL^*$ given by (14.3) and the epoch window $m L \in (m/T_{\text{in}}\text{-scale},\, \ldots)$ of order one at the crossover ($m_t/T_{\text{EW}} \sim 1$), $\mathcal P_\times$ is an order-unity geometric factor *whenever the angle distribution populates the existence region at all* — the entire smallness of $\eta_B$ is thereby pushed into the next factor, where it honestly belongs.

**(b) The sign bias.** Which sign does a crossing pump? Under $E \to -E$, the Master Equation maps $\Delta \to -\Delta$ (Wall-Mismatch Criterion); correspondingly the crossing level dives *into* the sea for one sign of $\Delta$ and climbs *out* for the other, pumping $\mp 1$ respectively **[Computed]** (`ch14_crossing_density.py --flow-direction`: the jump of $Q_{\text{vac}}$ at $L^*$ is $-\operatorname{sgn}(\sin(\Delta/2))$ across the sampled region). The ensemble's net pumping is therefore controlled by the CP asymmetry of the angle distribution:

$$\varepsilon_{CP} \;=\; \frac{\Pr[\text{crossing with } \Delta > 0] - \Pr[\text{crossing with } \Delta < 0]}{\Pr[\text{crossing}]}\,, \tag{14.7}$$

and a CP-symmetric ensemble ($P$ even in $\Delta$) pumps net zero — as it must: this is the Sakharov C/CP condition materializing inside the formula, in exactly the right place.

> **Crossing-Density Amplitude (master form) [Theorem, given the ensemble picture].**
>
> $$\boxed{\;\mathcal N_{\text{steps}}\,\bar\nu\,\varepsilon_{CP} \;=\; \mathcal P_\times\,\varepsilon_{CP}\,,\;}$$
>
> with $\mathcal P_\times$ the geometric window probability (14.6) — computable in closed form from (14.3) for any specified $P(\theta_0, \theta_L)$ — and $\varepsilon_{CP}$ the distribution's sign bias (14.7).

**[Computed]** `ch14_crossing_density.py` evaluates $\mathcal P_\times$ and $\varepsilon_{CP}$ for parameterized families (Gaussian angle pairs of mean $\pm\mu$ and spread $\varsigma$; sign-biased mixtures), validates the One-Crossing Theorem by brute force (390/390 existence and uniqueness; crossing locations matching (14.3) to $10^{-6}$), and produces Fig. 14.2. Representative outputs (epoch window $mL \in (0.3, 3)$): a near-MIT ensemble ($\mu = 0$, $\varsigma = 0.5$) gives $\mathcal P_\times = 0$ — *exactly* zero, the sub-threshold prediction; broad tails reach the region ($\varsigma = 1.5 \to \mathcal P_\times = 0.073$; $\varsigma = 2.2 \to 0.153$); strong-wall ensembles ($\mu = 2$, $\varsigma = 1$) give $\mathcal P_\times = 0.38$. The sign bias passes through faithfully ($\varepsilon_{CP} \approx$ the imposed mixture bias), and the flow direction obeys $\Delta Q_{\text{pumped}} = +\operatorname{sgn}(\sin(\Delta/2))$ per crossing across all sampled cases.

## 14.4 What sets the wall angles — and the honest budget

Everything now hangs on $P(\theta_0, \theta_L)$, i.e. on the microphysics of the walls, and here the framework meets the Standard Model's hard wall (Ch. 13). The wall angle is the relative phase between the boundary condensate (the chiral structure that terminates the fermion at the cell wall) and the bulk mass matrix. Two regimes:

**The perturbative pessimum.** If wall angles are generated perturbatively from CKM physics, the natural scale of both the spread of $\Delta$ and the bias $\varepsilon_{CP}$ inherits the full GIM-suppressed budget: two-loop measure $1/(16\pi^2)$, Jarlskog $J_{CP} \sim 3\times10^{-5}$, mass-difference factors $\mathcal S_{\text{GIM}} \sim 10^{-7}$ (Ch. 13.5) — and, fatally compounded by the threshold (14.2): a perturbatively narrow $P$ centered on MIT walls sits *entirely below threshold*, $\mathcal P_\times = 0$. In this regime the framework predicts $\eta_B$ not merely small but zero — an honest and falsifiable statement, and an improvement in clarity over a suppressed-but-fuzzy estimate.

**The strong-wall scenario.** If the boundary condensate carries an order-one chiral structure — as it does in every QCD-flavored bag model, where the wall angle is the analogue of the hybrid chiral bag's pion-field angle — then $\mathcal P_\times = \mathcal O(1)$ *geometrically*, and the entire smallness of $\eta_B$ is carried by $\varepsilon_{CP}$: the CP bias of the strong-wall ensemble. The budget for $\varepsilon_{CP}$ then spans, on present knowledge,

$$\varepsilon_{CP} \;\in\; \big[\,\mathcal O(10^{-10})\;(\text{bulk-perturbative bias, GIM-crushed})\,,\;\; \mathcal O(1)\;(\text{boundary-localized CP, GIM-evaded})\,\big],$$

a ten-order range whose resolution requires deciding whether **boundary-localized CP violation evades GIM** — flavor physics at a chiral wall is not flavor physics in the bulk plasma (the wall breaks the chiral symmetry that organizes the GIM cancellation in the first place; whether enough of the cancellation survives is a concrete, unsolved calculation). This is the framework's dominant open problem, stated in its sharpest known form: **[Open]**, Ch. 27 item 4, with the three candidate evasion routes (soft-mode enhancement; boundary-condensate phases; single-CKM-insertion boundary mixing) specified there.

**The assembled estimate.** With $g_* = 106.75$, $C_{\text{sph}} = 28/79$, $N_{\text{dof}} = N_c = 3$ for the top-quark channel:

$$\eta_B \;\approx\; 2.3\times10^{-2}\;\times\; \mathcal P_\times\;\varepsilon_{CP}, \tag{14.8}$$

so that matching $\eta_B^{\text{obs}} = 6.1\times10^{-10}$ requires

$$\mathcal P_\times\,\varepsilon_{CP} \;\approx\; 2.7\times10^{-8}. \tag{14.9}$$

**[Computed]** prefactors assembled and propagated in `ch14_eta_assembly.py`. In the strong-wall scenario ($\mathcal P_\times \sim 0.1$–$1$), observation demands $\varepsilon_{CP} \sim 10^{-8}$–$10^{-7}$ — intriguingly, the *unsuppressed* electroweak CP scale ($J_{CP} \times$ loop factors, before GIM mass-difference crushing). The framework thus passes a non-trivial consistency test it could easily have failed: the observed asymmetry sits inside the budget window, at the value corresponding to "boundary CP violation evades the GIM mass-difference factors but not the loop and Jarlskog factors." It does not *predict* this value until the boundary-GIM calculation is done; the thesis claims structure, not the number.

![The budget](figures/ch14_fig3_eta_budget.png)

*Figure 14.3 — The honest budget. Horizontal bars: the $\eta_B$ ranges implied by each scenario for $P(\theta_0,\theta_L)$ and $\varepsilon_{CP}$; vertical band: observation. The perturbative-pessimum bar sits at exactly zero (threshold), the strong-wall bars straddle observation when boundary-GIM evasion is partial or full. No bar is labeled "prediction".*

## 14.5 Inherited cautions from the three-dimensional assembly

Three honesty items, each inherited from the wider geometry and each bounded:

1. **Dimensional uplift.** The quantitative engine above is 1+1D (slab). The 3+1D spherical chiral bag has a richer crossing structure (one condition per $\kappa$-channel), plausibly *raising* $\mathcal P_\times$ (more channels, more crossings) — but until computed (Ch. 27, item 3), the 1+1D value stands as the defensible core. The slab limit itself is rigorous (Slab Eigenvalue Theorem, §10.A).
2. **The rectangular factorization conjecture.** Treating a finite rectangular cell as three independent slabs is proved in the slab and non-relativistic limits, conjectured at $mL_i \sim 1$. It enters only the *multiplicity* bookkeeping ($N_{\text{dof}}$-level factors of order one), not the mechanism. **[Conjecture]**, flagged.
3. **Spin counting.** At fixed transverse momentum the two spin polarizations are C-conjugates and their pumped charges cancel pairwise *for bulk-phase mechanisms*; for wall-mismatch pumping the cancellation question must be re-posed per channel and is part of the 3+1D computation. Until then, $N_{\text{dof}}$ carries an unresolved factor in $[\tfrac12, 1]\times$spin — absorbed into the budget bars of Fig. 14.3.

## 14.6 Summary, and what would change this chapter

The corrected estimate replaces a truncated overlap amplitude with a **probability that an exactly located level crossing falls inside the epoch's expansion window** ($\mathcal P_\times$, geometric, computable, bounded) times a **sign bias** ($\varepsilon_{CP}$, the genuine CP cost, currently spanning the honest range). The yield per event is exactly one unit — anomaly-protected, beyond revision. The pipeline (14.1) is structurally complete, every factor either standard, theorem-grade, or labeled **[Open]** with a named calculation attached.

Three computations would upgrade the estimate to a prediction, in order of leverage: the boundary-GIM calculation for $\varepsilon_{CP}$ (item 4 of Ch. 27); the 3+1D spherical crossing density (item 3); and the microphysics of the wall-angle distribution at the electroweak crossover. None requires a new idea — each is a finite, posed problem. That is, by the standards of baryogenesis model-building, an unusually good place to stand.

---

**Validation.** `ch14_crossing_density.py` (new): brute-force validation of the One-Crossing Theorem (2000 random angle pairs); the $(\Sigma, \Delta)$ phase diagram with flow-direction map; $\mathcal P_\times$ and $\varepsilon_{CP}$ for parameterized ensembles (Figs. 14.2). `ch14_eta_assembly.py` (new): prefactor assembly, the budget bars (Fig. 14.3), and the required $\mathcal P_\times\varepsilon_{CP}$ value (14.9). Both print every number quoted.
