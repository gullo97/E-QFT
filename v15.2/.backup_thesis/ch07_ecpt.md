# Chapter 7 — ECPT: one symmetry, two arrows

> **Author: Fable 5 · Complexity ◇4 · Status: draft v1 (pending global review pass)**
> Depends on: Ch. 3, 6. Feeds: Part II framing (Ch. 11–12, 14).

---

Chapter 6 ended with a broken symmetry: the dynamics of confined many-body matter distinguishes expansion from contraction, exponentially. In particle physics, the discovery that a discrete operation is *not* a symmetry of nature has never been the end of a story — it has always been the beginning of one. Parity violation led to charge–parity ($CP$); $CP$ violation led to the realization that only the full product $CPT$ stands on a theorem. The pattern is: when a discrete operation breaks, ask *what completion of it survives*.

This chapter applies that pattern to expansion-reversal. We define the operation $\hat E$ that swaps expansion with contraction, establish precisely in what senses it is broken, and then state the framework's organizing conjecture: that the *product* $\hat E\hat C\hat P\hat T$ is exact. The conjecture's bite is immediate: if the product is exact and $E$ is broken, then $CPT$ — exact in quantum field theory by theorem — cannot absorb the breaking, so $CP$ must break *in correlation with* the expansion arrow. One symmetry, two arrows: the universe grows, and it grows *matter*.

We are scrupulous about epistemic status throughout: the definitions and breaking statements are theorems within the framework; the ECPT completion is a **[Conjecture]**, and §7.6 states exactly what survives without it (answer: every mechanism of Part II — they need the Sakharov ingredients, not the conjecture).

## 7.1 The expansion-reversal operator

> **Definition ($\hat E$) [Definition].** On the geometry-sector arena, $\hat E$ is the antiunitary-free, purely structural operation that exchanges the transition operators:
>
> $$\hat E\,\hat T_+\,\hat E^{-1} = \hat T_-\,, \qquad \hat E\,\hat T_-\,\hat E^{-1} = \hat T_+\,, \qquad \hat E^2 = \mathbb 1. \tag{7.1}$$
>
> On the iso-energy ladder it is realized as the reflection of the ladder about the current rung: configurations are untouched, *directions of motion* are reversed.

Three properties pin down what $\hat E$ is and is not.

**(a) $\hat E$ is not time reversal.** Time reversal $\hat T$ is antiunitary, conjugates phases, and reverses *all* motion. $\hat E$ is unitary and reverses only the *geometry* directionality: a state moving up the ladder maps to one moving down it, with matter content untouched. In a world with no geometry dynamics, $\hat T$ is as meaningful as ever while $\hat E$ is empty. The two can be composed independently — which is exactly why $E$ earns a separate slot in the group of §7.3.

**(b) $\hat E$ commutes with the single-sector physics.** Within a fixed sector $\mathcal H(L)$, $\hat E$ acts as the identity; it has nontrivial action only on the inter-sector dynamics. Hence any breaking of $E$ must come from the *coupling between geometry sectors* — it is a property of geometry change, not of matter.

**(c) $\hat E$ is the lattice shadow of $a \to 1/a$.** In the cosmological reading (Ch. 17: $L$ = scale factor), $\hat E$ inverts the expansion arrow, the discrete ancestor of the scale-factor inversions familiar from bounce cosmologies. We note the connection and do not lean on it.

## 7.2 E-breaking: two mechanisms, one origin

> **Proposition (E-breaking) [Theorem].** The framework breaks $E$-symmetry twice over:
>
> **(i) Kinematically.** The iso-energy ladder has a floor and no ceiling: $\hat T_-$ annihilates the $n = 1$ rung (Natural Boundary, Ch. 2/3), while $\hat T_+$ never annihilates. No relabeling can restore the symmetry of a half-infinite ladder.
>
> **(ii) Dynamically.** With $N$ spectators, the Expansion-Preference Theorem gives $P_{\text{con}}/P_{\text{exp}} = (1 - n_1^{-2})\prod_\alpha \mathcal N_{\text{con}} \to 0$ exponentially in $N$ (Ch. 6): even far from the floor, the measure of quantum mechanics is one-sided.

Both mechanisms have the same root, and it is worth saying plainly: *boundaries*. The floor at $n = 1$ is the boundary of the ladder; the contraction deficit is the boundary of the box cutting into wavefunctions (Lemma 3.1(b)). The slogan — boundaries break the symmetries that the bulk would keep — is about to become the load-bearing theme of this thesis: in Part II it is a *boundary* (the bag wall) that makes CP phases physical and a boundary mismatch that pumps charge, and in Ch. 26 the two appearances are identified as one structure. The reader should file the slogan now.

## 7.3 The ECPT group

Adjoin $\hat E$ to the standard discrete operations. Each of $E, C, P, T$ squares to the identity (on the relevant ray space) and they mutually commute as abstract labels, so the group generated is

$$G_{ECPT} \;=\; \mathbb Z_2^E \times \mathbb Z_2^C \times \mathbb Z_2^P \times \mathbb Z_2^T \;\cong\; \mathbb Z_2^4, \qquad |G_{ECPT}| = 16, \tag{7.2}$$

with sixteen elements $\{\,\mathbb 1, E, C, P, T, EC, EP, \ldots, ECPT\,\}$ and a lattice of $67$ subgroups, of which the chain that matters here is

$$\langle ECPT \rangle \;\subset\; \langle CPT,\, E \rangle \;\subset\; G_{ECPT}.$$

![The ECPT group and its breaking pattern](figures/ch07_fig1_ecpt_lattice.png)

*Figure 7.1 — The corner of the $\mathbb Z_2^4$ subgroup lattice relevant to the conjecture. Solid: operations exact in the framework + Standard Model ($CPT$, conjecturally $ECPT$). Broken: $E$ (Ch. 6), $P$, $C$, $CP$, $T$ (Standard Model). The conjecture is that exactness propagates down the highlighted edge: $CPT$ exact and $ECPT$ exact, with all proper "partial products" containing $E$ broken.*

The Standard Model supplies the known pattern on the $C, P, T$ factors: $P$ and $C$ broken maximally, $CP$ and $T$ broken weakly, $CPT$ exact **[Standard]** (the CPT theorem: locality + Lorentz invariance + hermiticity). The framework supplies $E$ broken (§7.2). The open slot is the product.

## 7.4 The conjecture

> **ECPT Conjecture [Conjecture].** The combined operation $\hat E\hat C\hat P\hat T$ — reverse the expansion arrow, conjugate charges, reflect space, reverse time — is an exact symmetry of the full dynamics on the geometry-sector arena.

*Motivation, honestly weighed.* (i) **Pattern**: every previously discovered breaking of a discrete spacetime-flavored operation has been compensated within a larger product; ECPT is the minimal completion available once $E$ exists. (ii) **Structural analogy**: the CPT theorem rests on the impossibility of separating its factors using only local, Lorentz-invariant dynamics; in the sector arena, expansion-reversal composed with motion-reversal ($T$) relates the two directions of ladder traffic, and the matter operations ($C, P$) relate the two charge branches that ladder traffic populates (Part II) — the product is the operation with a fighting chance of exactness. (iii) **Numerical oddness check**: in the relativistic computations of Part II, the charge observable obeys $\Delta Q(+\delta) = -\Delta Q(-\delta)$ to machine precision ($\sim 10^{-15}$; Ch. 11, `ch11_polarization.py`) — exactly the parity in the CP-phase that ECPT-invariance of the underlying dynamics requires of an $E$-odd observable.

*What would constitute proof.* A derivation from the (yet-to-be-specified) microscopic Lagrangian of the coupled matter–geometry system, parallel to the CPT theorem's derivation — currently out of reach because the time/lapse sector of the geometry dynamics is itself open (Ch. 27, item 5). The conjecture is a target, not an assumption silently used: we now state its consequences, and then prove that Part II does not depend on it.

## 7.5 The physical content: two correlated arrows

Suppose the conjecture holds. Then in any epoch and ensemble:

$$\underbrace{ECPT \text{ exact}}_{\text{conjecture}} \;\wedge\; \underbrace{E \text{ broken}}_{\text{Ch. 6 theorem}} \;\wedge\; \underbrace{CPT \text{ exact}}_{\text{QFT theorem}} \;\;\Longrightarrow\;\; CP \text{ broken, correlated with } E. \tag{7.3}$$

The logic: $ECPT = E \cdot (CPT)$ as group elements; if both $ECPT$ and $CPT$ were realized exactly *on the same dynamics in which $E$ is broken*, then $E = (ECPT)(CPT)^{-1}$ would be exact — contradiction. The escape is that $CPT$ exactness (a fixed-geometry QFT statement) and $ECPT$ exactness (a statement including geometry traffic) refer to the two factors of the dynamics, and consistency forces the *matter-sector* compensation: the $E$-breaking bias of geometry traffic must be accompanied by a $CP$-breaking bias of the matter transported by that traffic, with the two biases locked in sign.

In plain terms: **the same dynamics that prefers growing boxes must prefer one charge sign in what the growing boxes create.** The cosmological arrow (expansion) and the baryonic arrow (matter over antimatter) become two projections of a single broken generator — not two coincidences requiring two explanations. Part II is the quantitative cash-out: *which* CP violation, *what* it produces under expansion, and *how much*.

![Two arrows from one breaking](figures/ch07_fig2_two_arrows.png)

*Figure 7.2 — The conjecture's content. The exact diagonal ($ECPT$) with the two broken projections: geometry direction ($E$) and matter charge ($CP$). Flipping the sign of the CP phase maps the universe to its expansion-reversed conjugate — the origin of the $\Delta Q(\delta) = -\Delta Q(-\delta)$ oddness verified in Ch. 11.*

## 7.6 The Sakharov conditions, previewed

Any mechanism for a baryon asymmetry must satisfy Sakharov's three conditions (derived and referenced properly in Ch. 13; the table is a preview and a promise):

| Sakharov condition | Where this framework supplies it | Status |
|---|---|---|
| Baryon-number violation | charge production by geometry change: spectral flow under expansion (Ch. 9, 12), converted by sphalerons (Ch. 13) | mechanism **[Theorem/Computed]**; conversion **[Standard]** |
| C and CP violation | wall-localized CP violation: mismatched chiral boundary angles (Ch. 12); magnitude budget (Ch. 14) | mechanism **[Theorem]**; magnitude **[Open]** |
| Departure from equilibrium | the expansion arrow itself: Expansion-Preference Theorem (Ch. 6) | **[Theorem]** within the framework |

The third row is the structural novelty worth advertising: standard electroweak baryogenesis fails partly for lack of a strong first-order phase transition to push the plasma out of equilibrium; here the out-of-equilibrium push is *permanent* — it is the arrow of Ch. 6, operating at every epoch.

## 7.7 What survives if the conjecture fails

Everything quantitative. The Part II results — the Spectral-Flow Master Theorem, the triviality of bulk CP phases, the polarization diagnosis, quantized pumping at wall mismatches, the $\eta_B$ assembly — are theorems and computations about Dirac fields in changing boxes; none of their proofs invokes $\hat E$ or the conjecture. What would be lost is the *organizing claim*: that the sign correlation between the two arrows is law rather than accident, i.e. that a universe which expands must (given the conjecture) tilt the same way in charge. ECPT is the frame of the picture, not the canvas. We flag its derivation as Ch. 27 item 5-adjacent, and proceed — for the next eight chapters — entirely on the canvas.

---

**Validation.** No numerics native to this chapter. The ECPT-oddness check $\Delta Q(+\delta) + \Delta Q(-\delta) = 0$ (machine precision) is produced in Ch. 11 by `ch11_polarization.py --oddness` and cross-referenced here. Figures 7.1–7.2 are drawn schematics (no script).
