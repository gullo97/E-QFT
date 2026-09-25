# Chapter 1 — Introduction: three puzzles, one box

> **Author: Fable 5 · Complexity ◇4 · Status: draft v1 (written after Parts I–IV, per the production plan)**

---

## 1.1 Three facts in search of a mechanism

Cosmology rests on three structural facts that textbooks treat as inputs. **The universe expands** — and while dynamics tells us how expansion proceeds, the *direction* is an initial condition: run Einstein's equations backward and they comply without complaint. **Matter exists** — the baryon-to-photon ratio $\eta_B \approx 6\times10^{-10}$ is small but emphatically nonzero, and the Standard Model, run honestly through an early universe, predicts a value many orders of magnitude smaller still. **Geometry gravitates** — spacetime obeys Einstein's equations, whose form general relativity takes as its starting axiom and whose strength ($G$) it cannot address at all.

Three facts, conventionally three separate research programs. This thesis develops a framework in which they are projections of a single mechanism, and the mechanism begins — this is either its charm or its scandal, and the reader will decide — with the first system every quantum mechanics student meets: a particle in a box.

## 1.2 The observation, and the wager

The energy levels of a particle in a box of size $L$ satisfy $E_{\lambda n}(\lambda L) = E_n(L)$: stretch the box and promote the quantum number in proportion, and the energy does not change. On this homogeneity the thesis places one structural wager (stated as a postulate, defended by its consequences): *the box size is not a parameter but a quantum variable*, and the elementary dynamics of geometry is the discrete, energy-preserving step the homogeneity singles out — $n \to n \pm 1$ with $L \to \frac{n\pm1}{n}L$. The arena this forces is a direct sum of Hilbert spaces, one per geometry, with transitions between neighbouring sectors as the elementary events (Ch. 3).

Everything else in the thesis is consequences, and the wager pays out three times:

**The arrow (Part I).** In a many-body box, every geometry change must carry all the spectators along. Carrying them into a *larger* box costs nothing (an isometry); into a *smaller* box, each spectator pays a computable overlap toll. The tolls multiply: contraction is suppressed *exponentially in the particle number* (Ch. 6 — including the conversion of the historic probability-rule and factorization assumptions into theorems). For any macroscopic content, expansion is not preferred by initial conditions; it is preferred by the measure of quantum mechanics. Why the universe grows: because it is full.

**The matter (Part II).** Can the growing box create charge? The thesis's central negative result says: not the obvious way. Net charge production in any sudden geometry change is a purely *spectral* quantity — the change in the vacuum's own charge, $\Delta Q = -\tfrac12\Delta\eta$ (the Spectral-Flow Master Theorem, Ch. 9) — and a bulk CP-violating phase provably cannot move it (Ch. 10): the reproducible "asymmetries" such setups generate are diagnosed, exactly, as wall polarization (Ch. 11). The positive result is where the CP physics must then live: on the *boundaries*. Mismatched chiral wall angles charge the vacuum fractionally ($-\Delta/2\pi$), host level crossings at analytically known box sizes, and pump **exactly one unit of charge per crossing** under expansion — anomaly-protected, verified end to end with its control experiment (Ch. 12). The corrected cosmological estimate (Ch. 14) is then a *probability* (that a bag's crossing falls inside the epoch's expansion window) times a *sign bias* (the genuine CP cost), with the budget stated without mercy and the deciding calculation named.

**The geometry (Part III).** A box of size $L$ *is* a unit box with metric $L^2$ — a change of coordinates, not an analogy (Ch. 17) — so the sector ladder is a quantum scale factor, and a lattice of cells is a quantum spatial geometry, triad included (Ch. 18). Coupling the sectors through the framework's own symmetry generator selects, exactly, the Kasner orbit family of general relativity (Ch. 20: a quantum wavepacket holds the vacuum Kasner circle to $\Sigma_K = 0.9997$, one axis contracting). Closing the energy budget yields Friedmann cosmology, structure growth, and collapse at the classic thresholds (Ch. 21). The strength of gravity stops being an input: matter loops *induce* the geometry's stiffness — attractive, with the conformal mode's notorious wrong sign as an output (Ch. 22) — predicting the cell scale at the Planck length and printing Newton's $1/r^2$ (Ch. 24) before the finale: the transverse-traceless stiffness of the cell lattice, measured, converges on the parameter-free one-loop coefficient $1/48\pi^2s^2$ for both polarizations, massless, with the DeWitt signature — the graviton (Ch. 25).

The closing synthesis (Ch. 26) reads the three payouts as one sentence: *matter on a dynamical, discrete geometry whose boundaries carry the CP physics* — the expansion arrow is the motion, the baryons are the cargo, the emergent spacetime is the arena, and the pump that links them is anomaly inflow on a moving boundary.

## 1.3 How this thesis argues

Three commitments distinguish the document's method, and knowing them will make it a faster read.

**Labels, everywhere.** Every claim of consequence carries one of six tags: **[Theorem]** (proved here, in full — no derivation is "left to the reader"), **[Computed]** (converged numerics; every number printed by a named script, regenerable via one command), **[Standard]** (established physics, rederived where load-bearing), **[Postulate]**, **[Conjecture]**, **[Open]** (with an attack plan in Ch. 27). The summary table of §1.4 is written in these tags, and Ch. 26.4 lists what the thesis does *not* claim.

**Dead ends are taught, not buried.** Part II's path runs through a seductive wrong answer — a convergent, symmetry-respecting, nonzero "charge asymmetry" that is provably not charge. Rather than excising the episode, the thesis dissects it (Ch. 11): the diagnosis (the Polarization Identity) is itself a theorem, the methodological yield (truncated pair-creation sums never measure net charge; spectral flow always does) governs every later computation, and the control-experiment discipline it forced is now built into the validation suite. A framework that can catch its own errors at theorem strength is worth more, not less, for having made one.

**Nothing before its introduction.** Concepts are announced when needed and built before use — small tools in marked Toolbox blocks, large ones in dedicated interludes (Ch. 5: quenches; Ch. 8: relativistic fermions in finite volume; Ch. 13: baryogenesis background; Ch. 16: the general-relativity toolbox). A reader with quantum mechanics and field-theory basics needs no prior acquaintance with bag models, spectral asymmetry, ADM variables, or BKL chaos; each arrives with its motivation attached.

## 1.4 Results at a glance

| Result | Status | Where |
|---|---|---|
| Iso-energy selection rule; sector arena; overlap closed forms | [Theorem] (2 kinematic postulates) | Ch. 2–3 |
| Golden-rule upgrade: probability rule + factorization derived | [Theorem + Computed] | Ch. 6 |
| Arrow of expansion: $P_{\text{con}}/P_{\text{exp}} \sim e^{-\sigma N}$, $\sigma \approx 1/n_1$ | [Theorem + Computed] | Ch. 6 |
| ECPT: one symmetry, two arrows | [Conjecture] (results independent of it) | Ch. 7 |
| Spectral-Flow Master Theorem $\Delta Q = -\tfrac12\Delta\eta$ | [Theorem] | Ch. 9 |
| Bulk CP phases exactly inert (dressing; mirror; zero charge) | [Theorem + Computed $10^{-16}$] | Ch. 10 |
| Truncated sums = wall polarization ($\to +0.0672$ at benchmark) | [Theorem + Computed] | Ch. 11 |
| Quantized pumping at wall mismatch; control at zero | [Theorem + Computed] | Ch. 12 |
| $\eta_B$: window probability × sign bias; budget and decider named | [Theorem + Computed + Open] | Ch. 14 |
| Box = metric (exact); expansion = FRW; quench = Parker | [Theorem] | Ch. 17 |
| Triad derived; conformal obstruction; $\mathfrak{gl}(3)$ completion | [Theorem] | Ch. 18 |
| Shape-space splitting (DeWitt-type, derived) | [Theorem + Computed] | Ch. 19 |
| Generator coupling selects Kasner class; $\Sigma_K = 0.9997$ | [Theorem + Computed] | Ch. 20 |
| Friedmann sub-percent; growth/collapse at $1.062/1.686$ | [Computed; constraint Postulate] | Ch. 21 |
| Stiffness induced: sign measured; $\bar L \approx 1.6c\,\ell_P$ | [Computed + Theorem] | Ch. 22 |
| Compensator $= 1/L_{\text{cell}}$; Kibble–Sciama completion | [Framework + Standard] | Ch. 23 |
| Newton $r^{-1.008}$, $d^{-2.05}$; BKL map to $6\times10^{-4}$ | [Computed] | Ch. 24 |
| Graviton: $\kappa_{TT} \to -1/48\pi^2s^2$ ($1$–$2\%$), 2 polarizations, massless, DeWitt sign | [Computed, parameter-free] | Ch. 25 |
| Not claimed: $\Lambda$, signature derivation, $\eta_B$ value, ECPT proof | — | Ch. 26.4 |

## 1.5 Reading routes

*The hour:* Ch. 1 → 9 → 12 → 20 → 25 → 26. *The skeptic's route:* Ch. 9–11 (the theorems and the dissected artifact), then Ch. 25 (the falsifiable measurement), then Ch. 27 (what would kill it). *The full read:* in order; interludes (5, 8, 13, 16) may be deferred until first cited. The dependency map below is the document's circuit diagram.

![Dependency map](figures/ch01_fig1_dependency_map.png)

*Figure 1.1 — Chapters as nodes, logical dependence as arrows; interludes dashed. (drawn)*

![Three faces](figures/ch01_fig2_three_faces.png)

*Figure 1.2 — The thesis in one picture: one structure (center), three projections. Chapter 26 returns to this figure with every arrow labeled by its theorem. (drawn)*
