# The Expanding Quantum Box

## Iso-Energetic Geometry, the Arrow of Expansion, Spectral-Flow Baryogenesis, and Emergent Gravity

*Master document — title page, abstract, reading guide, table of contents.*

---

## Abstract

This thesis develops, from first principles, a single framework in which three of cosmology's structural facts — that the universe expands rather than contracts, that matter outnumbers antimatter, and that geometry obeys Einstein-like dynamics — arise from one underlying mechanism: quantum matter confined on a *dynamical, discrete geometry* whose elementary changes are constrained by an iso-energy rule.

Part I builds the kinematics. The spectrum of a particle in a box is invariant under the simultaneous rescaling of quantum number and box size; promoting the box size to a quantum variable turns this homogeneity into a selection rule for elementary geometry transitions, a Hilbert space that is a direct sum of geometry sectors, and a derived tight-binding dynamics on the resulting ladder. In a many-body system the transition probabilities are governed by elementary overlap integrals, and contraction is exponentially suppressed relative to expansion: the arrow of expansion is a theorem about confined quantum matter, not an initial condition.

Part II makes the matter relativistic and asks whether expansion can create charge. The central result is that net charge production in any sudden quench is a purely spectral quantity — the change in the regulated spectral asymmetry of the Dirac operator. Bulk CP-violating phases are shown to be exactly inert (they dress a CP-conserving theory), and the pair-creation sums that seem to say otherwise are diagnosed, exactly, as a redistribution of vacuum polarization between the bag walls. The genuine mechanism is wall-localized CP violation: mismatched chiral boundary angles produce fractional vacuum charge, level crossings at analytically known box sizes, and quantized, anomaly-protected charge pumping under expansion. A corrected baryogenesis estimate is assembled on this engine, with an honest uncertainty budget.

Part III turns the same structure into gravity. A box of size $L$ is exactly a unit box with metric $L^2$; a lattice of cells supplies a triad (vielbein); the dilation Noether generator selects the Kasner orbit class of general relativity; the Hamiltonian constraint closes the dynamics into Friedmann cosmology, structure growth, and collapse; matter vacuum loops *compute* the geometry's stiffness — its attractive sign, Newton's inverse-square law, BKL chaos at the singularity, the Planck scale of the cells, and finally the graviton: a massless, two-polarization transverse-traceless wave whose kinetic coefficient equals the universal one-loop value $1/48\pi^2 s^2$, with the conformal mode carrying the opposite (DeWitt) sign.

The closing synthesis identifies the three threads as faces of one object — matter on a dynamical geometry whose boundaries carry the CP physics — and states precisely what is proved, what is computed, what is conjectured, and what remains open.

---

## How to read this thesis

**If you have one hour:** Ch. 1 (the program), Ch. 9 (the master theorem of Part II), Ch. 12 (the pumping mechanism), Ch. 20 (Kasner from the ladder), Ch. 25 (the graviton), Ch. 26 (synthesis).

**If you are checking the physics:** every major claim carries an honesty label — **[Theorem]**, **[Computed]**, **[Standard]**, **[Postulate]**, **[Conjecture]**, **[Open]** — and every quoted number is printed by a named script in `v15/scripts/` (App. C describes the standards; `run_all.py` regenerates everything).

**If you are new to a tool:** large background topics have purpose-built interlude chapters — Ch. 5 (quenches), Ch. 8 (relativistic fermions in finite volume), Ch. 13 (baryogenesis background), Ch. 16 (general-relativity toolbox) — and small ones appear as `> **Toolbox:**` blocks where first needed. Nothing is used before it is introduced.

**Dependency map.** Parts must be read in order, but interludes can be deferred until referenced: see `figures/ch01_fig1_dependency_map.png`.

---

## Table of contents

All chapters are **written (draft v1, Fable 5)**; the document is complete end to end and awaits the global review pass (gates G1–G4, skeleton §6).

| # | Chapter |
|---|---------|
| — | **Part I — The expanding quantum box** |
| 1 | [Introduction: three puzzles, one box](ch01_introduction.md) |
| 2 | [The particle in a box and the iso-energy structure](ch02_iso_energy.md) |
| 3 | [Geometry as a quantum variable: sectors, transitions, overlaps](ch03_geometry_sectors.md) |
| 4 | [Dynamics on the geometry ladder](ch04_ladder_dynamics.md) |
| 5 | [Interlude: sudden, adiabatic, and the quench concept](ch05_quench_interlude.md) |
| 6 | [Many bodies and the arrow of expansion](ch06_arrow_of_expansion.md) |
| 7 | [ECPT: one symmetry, two arrows](ch07_ecpt.md) |
| — | **Part II — Charge from geometry** |
| 8 | [Interlude: relativistic fermions in a finite box](ch08_dirac_interlude.md) |
| 9 | [Vacuum charge and the Spectral-Flow Master Theorem](ch09_master_theorem.md) |
| 10 | [CP phases in a box: what is physical and what is frame](ch10_cp_phases.md) |
| 11 | [The polarization diagnostic](ch11_polarization.md) |
| 12 | [Wall-localized CP violation and quantized charge pumping](ch12_pumping.md) |
| 13 | [Interlude: baryogenesis background](ch13_baryogenesis_interlude.md) |
| 14 | [From pumping to η_B: the corrected estimate](ch14_eta_b.md) |
| 15 | [Why geometry must be dynamical: the normalization diagnostic](ch15_imn.md) |
| — | **Part III — Emergent gravity** |
| 16 | [Interlude: the general-relativity toolbox](ch16_gr_toolbox.md) |
| 17 | [The dictionary: a box is a metric, expansion is cosmology](ch17_dictionary.md) |
| 18 | [A lattice of cells: the triad, the shears, the conformal obstruction](ch18_triad.md) |
| 19 | [Shape space: the Splitting Theorem](ch19_splitting.md) |
| 20 | [Minisuperspace dynamics: the Kasner Selection Theorem](ch20_kasner.md) |
| 21 | [Closing the constraint: Friedmann universes and emergent attraction](ch21_friedmann.md) |
| 22 | [The stiffness from matter loops: induced gravity and the Planck cell](ch22_stiffness.md) |
| 23 | [Weyl gauging, the compensator identified, and the road to GR](ch23_weyl.md) |
| 24 | [Newton's law and BKL chaos on the cell lattice](ch24_newton_bkl.md) |
| 25 | [The graviton, measured](ch25_graviton.md) |
| 26 | [Synthesis: three faces of one structure](ch26_synthesis.md) |
| — | **Part IV — Outlook** |
| 27 | [Open problems, with attack plans](ch27_open_problems.md) |
| — | **Appendices** |
| A | [Bogoliubov transformations in finite volume](appA_bogoliubov.md) |
| B | [Heat kernel, zeta function, Seeley–DeWitt](appB_heat_kernel.md) |
| C | [Numerical methods and reproducibility](appC_numerics.md) |
| D | [Notation and glossary](appD_notation.md) |

---

## Honesty labels (used throughout)

| Label | Meaning |
|---|---|
| **[Theorem]** | Proved in this thesis, derivation shown in full. |
| **[Computed]** | Converged numerical result; script named, number printed by `run_all.py`. |
| **[Standard]** | Established literature; restated (and derived where load-bearing). |
| **[Postulate]** | Structural choice of the framework, stated as such. |
| **[Conjecture]** | Motivated but unproven. |
| **[Open]** | Known gap; attack plan in Ch. 27. |

## Global conventions (quick reference; full table in App. D)

$\hbar = c = 1$ (except Ch. 2, where $\hbar$ is kept explicit). Metric signature $(+,-,-,-)$. 1+1D Dirac representation: $\gamma^0 = \sigma_3$, $\gamma^1 = i\sigma_2$, $\gamma_5 = \gamma^0\gamma^1 = \sigma_1$. MIT bag boundary condition $-i\,n\!\cdot\!\gamma\,\psi = \psi$ with $n$ the outward normal; chiral-bag condition $-i\,n\!\cdot\!\gamma\,e^{i\theta\gamma_5}\psi = \psi$. Reduced Planck mass $M_P^2 = 1/8\pi G$. Box modes $\varphi_n(x) = \sqrt{2/L}\,\sin(n\pi x/L)$. Spectral asymmetry $\eta = \sum_k \operatorname{sgn}(E_k)$ (heat-kernel regulated). Wall-angle combinations $\Delta = \theta_L - \theta_0$, $\Sigma = (\theta_0 + \theta_L)/2$. Logarithmic scale factor $\alpha_a = \ln a_a$; relational exponents $p_a = \dot\alpha_a / \sum_b \dot\alpha_b$; anisotropy invariant $\Sigma_K = \sum_a p_a^2$ (subscript $K$ distinguishes it from the wall-angle mean).
