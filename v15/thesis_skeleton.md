# v15 Thesis Skeleton

## *The Expanding Quantum Box: Iso-Energetic Geometry, the Arrow of Expansion, Spectral-Flow Baryogenesis, and Emergent Gravity*

**Document role.** This is the construction blueprint for the v15 master thesis. It fixes the chapter structure, the named theorems, the pedagogical strategy, the figure plan, the validation-script plan, the per-chapter complexity rating, and the model assignment (Fable 5 vs Opus 4.8) for the writing phase. The thesis itself will live in `v15/thesis/` as one file per chapter plus a master document; nothing in the final thesis will reference v12/v13/v14 documents or call results "the v13.x result" — every result is derived, named, and validated in place.

**Thesis statement (the one-sentence version).** A single structure — quantum matter confined on a dynamical, discrete geometry whose changes are constrained by an iso-energy rule — produces, from its own mechanics: the arrow of cosmological expansion (exponential many-body preference), a matter–antimatter asymmetry (anomaly-protected spectral flow driven by boundary-localized CP violation), and gravity itself (Kasner/Friedmann/Newton/BKL dynamics and a massless two-polarization graviton induced by matter loops), with the Planck scale emerging self-consistently as the cell size.

---

## 0. Design principles (binding for every chapter)

These rules implement the requirements: precise but pedagogical, fully explicit derivations, no unexplained jargon, narrative built around figures.

**P1 — No concept before its introduction.** Every technical concept (Bogoliubov transformation, spectral asymmetry, DeWitt supermetric, anomaly inflow, vielbein, …) is *announced* when first needed ("to make this precise we need X, which the next section/chapter develops"), then *explained in a purpose-made section or interlude chapter* with motivation tied to why the narrative needs it, and only then *used*. Large topics get dedicated interlude chapters (Ch. 8, 13, 16); small ones get `> **Toolbox:**` blocks inside chapters.

**P2 — Every important result is a named theorem.** Stated in a quote block with name, hypotheses, statement, then a full derivation (every algebraic step that a careful reader would otherwise have to reproduce on paper), then a numerical validation pointing at a figure and a script. The canonical theorem names are fixed in §2 of this skeleton and must be used consistently.

**P3 — Figure-led narrative.** Wherever possible a section opens with (or quickly reaches) a figure that shows the phenomenon, and the derivation explains the figure. Each results chapter carries 2–6 figures: at least one *concept figure* (drawn schematic) and one *validation figure* (produced by a script, with the script named in the caption). Figure files: `thesis/figures/chNN_figM_shortname.png`. Captions are self-contained paragraphs.

**P4 — Honesty labels.** Every claim carries one of: **[Theorem]** (proved here), **[Computed]** (converged numerics, script cited), **[Standard]** (established literature, restated with full derivation where load-bearing), **[Postulate]** (structural choice of the framework, stated as such), **[Conjecture]** (motivated, unproven), **[Open]** (known gap, with attack plan in Ch. 27). The thesis defends itself by never blurring these.

**P5 — Dead ends are taught as diagnostics, not history.** The bulk-CP "charge production", the per-channel asymmetries, the IMN regulator: each appears inside the corrected narrative as a *named diagnostic* (e.g. the Polarization Identity explains exactly what the truncated sums measure), never as "an earlier version claimed…". No version numbers anywhere in the thesis.

**P6 — Markdown-safe mathematics.** KaTeX-compatible LaTeX only: `$…$` inline, `$$…$$` display with a blank line before and after; `aligned` (not `align`) inside display math; no custom macros; no `\labels` — equations referenced by name or by local tags like *(6.3)* written manually. Unicode for occasional inline symbols (η, κ, δ) is allowed in prose but formulas always in math mode.

**P7 — Reproducibility.** Every number quoted in the thesis is printed by a script in `v15/scripts/`, runnable standalone (`numpy/scipy/matplotlib` only), with a `run_all.py` driver and a shared `style.mplstyle`. Scripts are ported/adapted from the existing validated codebase where possible (provenance tracked in `v15/coverage_matrix.md`, an internal file that is *not* part of the thesis).

**P8 — Conventions fixed once.** $\hbar = c = 1$ except in Ch. 2 where $\hbar$ is kept explicitly for pedagogy; metric signature $(+,-,-,-)$; 1+1D Dirac representation $\gamma^0=\sigma_3,\ \gamma^1=i\sigma_2,\ \gamma_5=\sigma_1$; reduced Planck mass $M_P^2 = 1/8\pi G$. A notation appendix (App. D) holds the global symbol table; chapters may not silently redefine symbols.

---

## 1. Complexity scale and model assignment policy

**Complexity ratings** (per chapter, also flagging the hardest section):

| Rating | Meaning |
|---|---|
| ◇1 | Textbook material, low risk, exposition only |
| ◇2 | Standard graduate physics, careful but routine derivations |
| ◇3 | Involved derivations of established framework results; errors possible but checkable against existing scripts |
| ◇4 | Critical synthesis or load-bearing derivation; the thesis's credibility rests on these being exactly right |
| ◇5 | Frontier: new formulation/derivation or subtle corrected results; requires the strongest model and adversarial self-checking |

**Model policy** (as requested):

- **Fable 5**: all ◇4–◇5 chapters; all novel formulations and derivations; *every* validation script (including ports of old scripts — each port gets re-run and its numbers re-verified before being cited); final review pass of every chapter regardless of author.
- **Opus 4.8**: ◇1–◇3 chapters that expose established physics (Dirac/MIT bag review, GR toolbox, Sakharov/EWBG background, Weyl gauging history, appendices A/B/D), writing *to this skeleton's specification* with theorem names and notation fixed here.
- **Split chapters** are marked section-by-section.

---

## 2. Canonical theorem and definition names (global registry)

These names are fixed now so all chapters cite them identically. (Chapter of proof in parentheses.)

**Foundations**
- **Iso-Energy Selection Rule** — $n \to n\pm1$, $L \to \frac{n\pm1}{n}L$ as the unique nearest-neighbour energy-preserving step (Ch. 2).
- **Geometry-Sector Decomposition** — the arena $\mathcal H = \bigoplus_L \mathcal H_{\text{matter}}(L)$ with embedding maps $\iota$ (Ch. 3). *(Definition, elevated to the framework's primary postulate.)*
- **Telescoping Composition Theorem** — $\hat T_+^k|n_0;L_0\rangle = |n_0{+}k;\frac{n_0+k}{n_0}L_0\rangle$ (Ch. 3).
- **Geometric Overlap Theorem** — promoted-mode overlap $\mathcal O_{n\to n+m} = \sqrt{n/(n+m)}$; sinc-type spectator coefficients (Ch. 3).
- **Dilation-Generator Matrix-Element Lemma** — $\langle n{+}1|\,x\partial_x\,|n\rangle = \frac{2n(n+1)}{2n+1} \approx n+\tfrac12$ (Ch. 3; the seed of the Kasner Selection Theorem in Ch. 20).
- **Natural Boundary Theorem** — the iso-energy ladder terminates at $n=1$ (self-adjointness; the kinematic seed of E-breaking) (Ch. 2, used in Ch. 7).
- **Expansion Norm Theorem** — extension-by-zero is an isometry: spectator weight $=1$ under expansion (Ch. 6).
- **Kinetic-Energy Preservation Theorem** — the spectator's expected kinetic energy is preserved under sudden expansion (sum rule) (Ch. 6).
- **Contraction Norm Theorem** — $\mathcal N_{\text{con}}(n_1,n_2) = \frac{n_1-1}{n_1} + \frac{\sin(2\pi n_2/n_1)}{2\pi n_2}$ (Ch. 6).
- **Sudden-Contraction Divergence Theorem** — spectator kinetic energy diverges under strict sudden contraction (Ch. 6).
- **Golden-Rule Embedding Theorem** — the probability rule "transition weight = embedding norm" *and* the $N$-spectator factorization are first-order perturbation theory in the geometry-sector arena: $P(t) = (gt)^2\|\iota\psi\|^2$, products over spectators (Ch. 6). *(This converts the framework's two historic postulates into theorems.)*
- **Expansion-Preference Theorem** — $P_{\text{con}}/P_{\text{exp}} = (1-1/n_1^2)\prod_\alpha \mathcal N_{\text{con}}(n_1, n_{2,\alpha})$, exponentially small in $N$ (Ch. 6).
- **ECPT Conjecture** — exactness of $\hat E\hat C\hat P\hat T$; E-breaking ⇒ correlated CP-breaking (Ch. 7).

**Charge from geometry**
- **Vacuum-Charge Formula** — $Q_{\text{vac}} = -\tfrac12\eta$, $\eta = \sum_k \operatorname{sgn}E_k$ (regulated) (Ch. 9).
- **Spectral-Flow Master Theorem** — for any bilinear charge-conserving sudden quench, $\Delta Q_{\text{net}} = Q_{\text{vac}}(L) - Q_{\text{vac}}(L') = -\tfrac12[\eta(L)-\eta(L')]$ (Ch. 9). *(The centerpiece of Part II.)*
- **Bulk-Phase Dressing Theorem** — the bulk-CP radial system is the ordinary MIT bag with mass $m\cos\delta$ under the non-unitary map $(g,f) \to e^{m_I r}(G,F)$; flat-measure per-channel asymmetries are frame artifacts (Ch. 10).
- **Spectral Mirror Theorem** — the 1+1D bulk-phase MIT spectrum is exactly $E \leftrightarrow -E$ symmetric for all $|\delta|<\pi/2$, all $L$ (Ch. 10).
- **Bulk-Phase Zero-Charge Theorem** — corollary of the two above + Master Theorem: sudden expansion with a bulk CP phase produces exactly zero net charge (Ch. 10).
- **Polarization Identity** — $\Delta Q_{\text{cross}}(\infty) = \tfrac12\operatorname{Tr}_{[0,L]}\operatorname{sgn}(H_{L'}) = \tfrac12\sum_k \operatorname{sgn}(E_k')\,w_k$: truncated Bogoliubov sums measure wall-polarization redistribution, never net charge (Ch. 11).
- **Chiral-Wall Master Equation** — $p\cos\!\frac{\Delta}{2}\cos(pL) = \sin(pL)\,[\,E\sin\!\frac{\Delta}{2} - m\cos\Sigma\,]$, $\Delta=\theta_L-\theta_0$, $\Sigma=\frac{\theta_0+\theta_L}{2}$ (Ch. 12).
- **Wall-Mismatch Criterion** — spectral asymmetry (hence vacuum charge, hence net production) exists iff $\sin(\Delta/2)\neq 0$ (Ch. 12).
- **Fractional Charge Law** — $Q_{\text{vac}} \approx -\Delta/2\pi \pmod 1$ (Goldstone–Wilczek/chiral-bag charge in 1+1D) (Ch. 12).
- **Crossing Condition** — zero mode at $\tanh(mL^*) = -\cos(\Delta/2)/\cos\Sigma$; **Quantized Pumping Theorem** — quench across a crossing pumps exactly one unit of charge (Ch. 12).
- **Crossing-Density Amplitude** — the corrected kinematic factor of the baryogenesis estimate: density of level crossings per iso-energy expansion step, in closed form from the Crossing Condition (Ch. 14; *new derivation*).
- **IMN No-Continuum Theorem** — mode-dependent normalization at fixed $L$ has no thermodynamic limit and violates microcausality; the sound instinct is implemented by the Geometry-Sector Decomposition instead (Ch. 15).

**Emergent gravity**
- **Dictionary Theorem** — "box of size $L$" = "unit box with spatial metric $g_{\xi\xi}=L^2$"; $L(t)$ = FRW scale factor; sudden expansion = step-function cosmological particle production (Ch. 17).
- **Metric Emergence Theorem** — cell edge vectors define $g_{ij} = \sum_a e_a^i e_a^j$; directional transitions act as computable rank-1 metric updates (Ch. 18).
- **Flatness (Commutativity) Theorem** — $[\hat T^{(a)}_+,\hat T^{(b)}_+]=0$: single-cell transitions are a flat connection (Ch. 18).
- **Conformal Obstruction Theorem** — gauging dilations alone generates only conformally flat metrics (zero Weyl tensor): no light bending, no gravitational waves; shears $\hat S^{(ab)}$ are the necessary completion, closing $\mathfrak{gl}(3)$ (Ch. 18).
- **Shape-Space Splitting Theorem** — diagonal deformations are overlap-flat and energy-carrying (discrete ladder); shear deformations are iso-energy and metric-carrying (continuous flow): the model's derived version of the DeWitt clock/propagating split (Ch. 19).
- **Kasner Selection Theorem** — embedding coupling ⇒ Milne class (isotropization); dilation-generator coupling $t(n) = g(n+\tfrac12)$ ⇒ exact Kasner/Jacobs orbit class with $p_a = g_a/\sum_b g_b$ (Ch. 20).
- **Constraint-Closure Result** — imposing $3\kappa\dot\alpha^2=\rho$ on the generator ladder yields the exact Friedmann power laws; **Cell-Growth Theorem** — neighbouring constraint-closed cells obey $\ddot\delta + 2H\dot\delta - 4\pi G\bar\rho\,\delta = 0$ (Ch. 21).
- **Induced-Stiffness Results** — 1D: $E_C = -\pi/24L$, vacuum Jeans instability, measured $\kappa^{\text{vac}}_{\text{grad}}<0$; 3+1D: $G_{\text{ind}}^{-1} = N_{\text{eff}}\Lambda_{\text{UV}}^2/12\pi$ with $\Lambda_{\text{UV}} = c/\bar L$ ⇒ **Planck-Cell Prediction** $\bar L \approx 1.6\,c\,\ell_P$ (Ch. 22).
- **Compensator Identification** — the Weyl compensator is the conformal part of the cell geometry, $\sigma \sim 1/L_{\text{cell}}$, making $M_P \sim 1/\bar L$ a consistency relation, not a postulate (Ch. 23).
- **Lattice Poisson Result** — $\kappa\hat\nabla^2\delta\alpha = -c\bar\rho\,\delta$, measured exponents $-1.008$ (potential) and $-2.05$ (force) (Ch. 24).
- **BKL Recurrence Verification** — the cell-model Kasner epochs obey $u\to u-1$, $u\to 1/(u-1)$ to $6\times10^{-4}$ including era reversal (Ch. 24).
- **DeWitt-Signature Measurement** — covariantly regulated matter loops induce $\kappa_{TT}=-C(s)$, $C(s)=1/48\pi^2 s^2$, both polarizations, conformal mode opposite sign, diffeo sector null: the graviton with the universal one-loop coefficient (Ch. 25).

---

## 3. Folder layout

```
v15/
├── thesis_skeleton.md            # this document
├── coverage_matrix.md            # internal provenance map (legacy → v15); not part of the thesis
├── thesis/
│   ├── 00_master.md              # title page, abstract, reading guide, linked TOC, notation quick-ref
│   ├── ch01_introduction.md … ch27_open_problems.md
│   ├── appA_bogoliubov.md  appB_heat_kernel.md  appC_numerics.md  appD_notation.md
│   └── figures/                  # chNN_figM_name.png (+ a few .gif animations)
├── scripts/
│   ├── style.mplstyle            # unified plot style (port of papers/style)
│   ├── chNN_*.py                 # one or more validation/figure scripts per results chapter
│   └── run_all.py                # regenerates every figure and prints every quoted number
└── data/                         # cached .json/.npz outputs (committed, small)
```

---

## 4. The narrative arc (how a reader travels through the thesis)

**Part I (Ch. 1–7), Foundations:** an elementary observation about the particle-in-a-box spectrum is taken seriously; geometry becomes a quantum variable living in a sum of sectors; transitions between sectors are governed by overlap integrals; many-body statistics makes expansion an arrow. *Reader needs only undergraduate QM; everything else is built.*

**Part II (Ch. 8–15), Charge from geometry:** to ask whether expansion can create matter, the box is made relativistic (interlude Ch. 8). The central discovery is that net charge creation is a purely *spectral* quantity (Spectral-Flow Master Theorem); bulk CP phases are shown to be secretly trivial, the popular truncated pair-creation sums are diagnosed as wall polarization, and the genuine mechanism — wall-localized CP violation, level crossings, quantized anomaly-protected pumping — is derived end-to-end and assembled into an honest baryogenesis estimate.

**Part III (Ch. 16–26), Emergent gravity:** the dictionary "box = metric" turns Part I's machinery into cosmology (interlude Ch. 16 supplies GR background). The cell lattice supplies a triad; the generator coupling selects the Kasner orbit class; the constraint closes the ladder into Friedmann universes with growth and collapse; matter loops *compute* the stiffness, Newton's constant, the Planck cell size, the conformal-factor sign, BKL chaos, and finally the graviton's kinetic structure with a parameter-free coefficient.

**Part IV (Ch. 27), Outlook:** the honest remainder with concrete attack plans.

The two big payoffs are deliberately placed as the closing chapters of Parts II and III, and Ch. 26 fuses them: *charge pumping is anomaly inflow on a moving boundary of the emergent spacetime* — the baryon asymmetry, the expansion arrow, and gravity are three faces of one structure.

---

# PART I — THE EXPANDING QUANTUM BOX

---

### Ch. 1 — Introduction: three puzzles, one box

**Goal.** State the three target phenomena (why the universe expands; why matter outnumbers antimatter; why geometry obeys Einstein-like dynamics), preview the unified mechanism, and give the reading map (figure: thesis dependency graph). Set expectations: what is theorem, what is computed, what is conjectured (P4 labels introduced here).

- Contents: the elevator pitch from the iso-energy observation up; the thesis statement; a "results at a glance" table with honesty labels; chapter map; conventions pointer.
- **Figures:** `ch01_fig1_dependency_map` (new schematic: chapters as nodes, arrows = logical dependence); `ch01_fig2_three_faces` (new schematic: expansion arrow / baryon asymmetry / emergent gravity as faces of the box+boundary structure).
- **Scripts:** none.
- **Complexity:** ◇4 (framing is load-bearing). **Model: Fable 5.** Length ≈ 5 pp.

---

### Ch. 2 — The particle in a box and the iso-energy structure

**Goal.** From $E_n(L) = n^2\pi^2\hbar^2/2mL^2$ to the homogeneity property $(n,L)\to(\lambda n,\lambda L)$ and the **Iso-Energy Selection Rule**. Honest framing from the start: this is a homogeneity of the spectrum function, promoted to a dynamical principle by an explicit **[Postulate]** (nearest-neighbour, energy-preserving steps).

- Contents: spectrum review; iso-energy hyperbolae; the selection rule derived from integrality; what the rule does *not* determine (timescales, dynamics) — flagged for Ch. 3–4; relation to Noether (the dilation generator $\hat D = \tfrac12(\hat x\hat p + \hat p\hat x)$ announced, full treatment deferred to Ch. 3).
- **Pedagogy:** keep $\hbar$ explicit here; one Toolbox block on "what counts as a symmetry vs a homogeneity".
- **Figures:** `ch02_fig1_iso_energy_curves` (rebuild of the iso-energy hyperbola plot, new style, no legacy labels); `ch02_fig2_selection_rule` (new: zoom on one hyperbola with the two allowed integer steps highlighted).
- **Scripts:** `ch02_iso_energy.py` (port of `plot_iso_energy_curves.py`).
- **Complexity:** ◇2. **Model: Opus 4.8** (Fable review). Length ≈ 6 pp.

---

### Ch. 3 — Geometry as a quantum variable: sectors, transitions, overlaps

**Goal.** Install the framework's true arena — the **Geometry-Sector Decomposition** $\mathcal H = \bigoplus_L \mathcal H(L)$ with embedding maps — *as the primary definition* (this is the corrected, theorem-grade formulation; historically it arrived late, here it comes first so that everything downstream is derivable). Define transition operators, prove **Telescoping** and the **Geometric Overlap Theorem**, introduce the dilation generator and prove the **Dilation-Generator Matrix-Element Lemma** with an explicit forward pointer: *"this innocuous $n+\tfrac12$ will select general relativity in Ch. 20."*

- Contents: sectors and embeddings ($\iota$ = extension-by-zero / restriction — the unique maps that "do nothing to the wavefunction", i.e. the sudden approximation made structural); $\hat T_\pm$ defined on sectors; telescoping; overlap integrals in closed form ($\sqrt{n/(n+m)}$ promoted modes, sinc-type otherwise — explicitly *not* Bessel); $\hat D$, its Ward identity, its matrix elements between box modes (full elementary integral shown).
- **Pedagogy:** Toolbox block "direct sums of Hilbert spaces"; Toolbox block "the sudden approximation" (short version; full interlude in Ch. 5).
- **Figures:** `ch03_fig1_sector_ladder` (new schematic: the $\bigoplus_L$ ladder with embeddings as arrows — *this becomes the thesis's signature diagram, reused in Ch. 20*); `ch03_fig2_overlap_wavefunctions` (rebuild); `ch03_fig3_telescoping` (rebuild).
- **Scripts:** `ch03_overlaps.py` (port/merge of overlap + telescoping plot scripts; prints $\mathcal O_{n\to n+m}$ table and the $\langle n{+}1|x\partial_x|n\rangle$ check $10.4762, 40.4938, 160.4983$).
- **Complexity:** ◇4 (novel ordering of the formulation; correctness of the embedding framing is load-bearing). **Model: Fable 5.** Length ≈ 9 pp.

---

### Ch. 4 — Dynamics on the geometry ladder

**Goal.** The tight-binding Hamiltonian on the iso-energy ladder; band structure $E(q) = -2J\cos q$; Bessel-function spreading of wavepackets; geometry superpositions at small $n$ and the quantum→classical crossover.

- Contents: construction of $H_{\text{TB}}$ from the derived hoppings; Bloch waves and the band; time evolution; low-$n$ regime as "quantum geometry" (superpositions of box shapes), large-$n$ classicality; discrete breaking of dilation symmetry at low quantum numbers. Forward pointer: *which* hopping law nature uses is a dynamical choice with dramatic consequences (Ch. 20).
- **Figures:** `ch04_fig1_band_structure` (rebuild); `ch04_fig2_bessel_spreading` (new: wavepacket $|c_n(t)|^2$ carpet plot); `ch04_fig3_geometry_superposition` (rebuild: $|c_k|^2$ distributions at $n=1$ vs $n=100$).
- **Scripts:** `ch04_ladder_dynamics.py`.
- **Complexity:** ◇3. **Model: Opus 4.8** (Fable review). Length ≈ 6 pp.

---

### Ch. 5 — Interlude: sudden, adiabatic, and the quench concept

**Goal.** Purpose-made short chapter on the sudden approximation, the adiabatic theorem, the adiabatic parameter of the framework, and the notion of a *quench* — the word that powers Part II. Includes the physically principled mode cutoff: a quench of duration $\tau$ is sudden for $\omega_n\tau \ll 1$, adiabatic for $\omega_n\tau \gg 1$ (this later justifies cutting mode sums physically rather than numerically).

- **Figures:** `ch05_fig1_regimes` (new: regime diagram in the $(\omega_n,\tau)$ plane).
- **Scripts:** none (conceptual).
- **Complexity:** ◇2. **Model: Opus 4.8.** Length ≈ 4 pp.

---

### Ch. 6 — Many bodies and the arrow of expansion

**Goal.** The first physics payoff. Spectator projection; the four norm theorems; then the chapter's centerpiece: the **Golden-Rule Embedding Theorem**, which converts the historic probability-rule and factorization *postulates* into first-order perturbation theory in the sector arena; finally the **Expansion-Preference Theorem** and the critical spectator number.

- Contents: active/spectator decomposition; extension-by-zero vs truncation; **Expansion Norm Theorem** (isometry, weight 1); **Contraction Norm Theorem** with the full elementary integral $\int_0^{\frac{n_1-1}{n_1}L}\frac2L\sin^2\frac{n_2\pi x}{L}dx = \frac{n_1-1}{n_1}+\frac{\sin(2\pi n_2/n_1)}{2\pi n_2}$; **Sudden-Contraction Divergence Theorem** (kinetic energy diverges ⇒ strict sudden contraction forbidden; resonant exceptions noted); the Golden-Rule derivation $P(t) = (gt)^2\|\iota\psi\|^2$ by Parseval, factorization as tensor-product structure, unitarity bookkeeping (missing weight = the transition *fails*); compounding over $N$ spectators, $P_{\text{con}}/P_{\text{exp}}$ exponentially small, $N_{\text{crit}}\approx 15$–$20$; interpretation as the arrow of expansion, with the honest caveats (iso-energy label misleading under contraction; representative $s=2$ point vs thermal regime $n_1\gg1$; rate-interpretation is a modelling choice now *derived* within the sector dynamics).
- **Pedagogy:** Toolbox block "Fermi's golden rule and first-order perturbation theory" (compact, complete).
- **Figures:** `ch06_fig1_spectator_projection` (rebuild); `ch06_fig2_norm_deficit` (rebuild: $1-\mathcal N_{\text{con}}$ vs $n_2$ for several $n_1$); `ch06_fig3_golden_rule_convergence` (new from existing toy: Parseval sum vs exact short-time evolution, the 0.965→0.992→1 convergence); `ch06_fig4_asymmetry_vs_N` (rebuild: exponential preference; annotate $N_{\text{crit}}$).
- **Scripts:** `ch06_norms.py`, `ch06_golden_rule.py` (port of `golden_rule_toy.py`; re-verify 0.76749 vs 0.76882 etc.).
- **Complexity:** ◇4–◇5 (the theorem upgrade is the part the supervisor will probe). **Model: Fable 5.** Length ≈ 11 pp.

---

### Ch. 7 — ECPT: one symmetry, two arrows

**Goal.** Define the expansion-reversal operator $\hat E$, establish E-breaking (kinematic: natural boundary at $n=1$; dynamical: spectator suppression), construct the ECPT group ($\mathbb Z_2^4$), state the **ECPT Conjecture** with scrupulous honesty, and derive its physical content: if $\hat E\hat C\hat P\hat T$ is exact and $E$ is broken, CP-breaking must be correlated with the expansion arrow. Map onto the three Sakharov conditions (preview table; full treatment Ch. 13–14).

- Contents: $\hat E$ definition and properties; distinction from time reversal; E-symmetry breaking (both mechanisms share the boundary origin); ECPT group and subgroup lattice; breaking pattern; the two-arrows claim; what survives if the conjecture fails (the mechanism chapters stand on their own — important defensive point).
- **Figures:** `ch07_fig1_ecpt_lattice` (new: $\mathbb Z_2^4$ subgroup/breaking diagram); `ch07_fig2_two_arrows` (new schematic).
- **Scripts:** none here (the ECPT-oddness numerical check lives in Ch. 11–12 where $\Delta Q(\delta)$ is computed).
- **Complexity:** ◇4 (conjecture framing must be exactly right). **Model: Fable 5.** Length ≈ 6 pp.

---

# PART II — CHARGE FROM GEOMETRY

---

### Ch. 8 — Interlude: relativistic fermions in a finite box

**Goal.** Everything Dirac the thesis needs, self-contained: the Dirac equation, why negative energies force the sea/Fock construction, MIT bag boundary conditions, chiral-bag boundary conditions, second quantization in finite volume, charge conjugation, and — written with special care because Ch. 9 stands on it — the symmetric ordering of the charge operator.

- Contents: 1+1D Dirac with fixed conventions (P8); plane-wave and box solutions; MIT condition $-i\,n\!\cdot\!\gamma\,\psi=\psi$, current confinement, spectrum $\tan(pL)=-p/m$, massless limit half-integer modes, $mL\to\infty$ Dirichlet limit; chiral-bag condition $-i\,n\!\cdot\!\gamma\,e^{i\theta\gamma_5}\psi=\psi$ introduced *here* as the natural one-parameter family (so Ch. 12 is no surprise); 3+1D spherical bag essentials ($\kappa$ quantum number, spinor harmonics — compact, only what Ch. 10 needs); second quantization, Fock space, Dirac sea; charge conjugation; $\hat Q = \tfrac12\int[\psi^\dagger\psi - \psi\psi^\dagger]$ forced by C-covariance — flagged as the hinge of Part II.
- **Pedagogy:** this is the model interlude chapter: each tool is introduced with "why the narrative needs it".
- **Figures:** `ch08_fig1_bag_spectrum` (rebuild: MIT eigenvalues vs $mL$, both branches); `ch08_fig2_boundary_geometry` (new schematic: walls, normals, boundary rays $(1,\pm i)$); `ch08_fig3_dirac_sea` (new schematic: sea filling, particle/hole bookkeeping).
- **Scripts:** `ch08_bag_spectrum.py` (transfer-matrix solver port — this becomes the workhorse `bag1d` module for Ch. 10–12).
- **Complexity:** ◇2–◇3 but foundational. **Model: Opus 4.8**, with **Fable 5 writing §8.x on the symmetrized charge operator** (load-bearing). Length ≈ 13 pp.

---

### Ch. 9 — Vacuum charge and the Spectral-Flow Master Theorem

**Goal.** The conceptual centerpiece of Part II. Define the spectral asymmetry $\eta$ with honest regulator discussion; derive the **Vacuum-Charge Formula** $Q_{\text{vac}}=-\tfrac12\eta$; connect to fractional fermion number (Jackiw–Rebbi; Goldstone–Jaffe) as **[Standard]** anchors; then prove the **Spectral-Flow Master Theorem** $\Delta Q_{\text{net}} = -\tfrac12[\eta(L)-\eta(L')]$ in full detail (charge conservation through bilinear quenches + change of normal-ordering reference).

- Contents: symmetric ordering recap → $Q_{\text{vac}}$; heat-kernel regulator $\eta(t)=\sum_k\operatorname{sgn}(E_k)e^{-t|E_k|}$, $t\to0$, with a worked toy example; fractional charge context; the two-line master proof written out in ten lines with every step justified; immediate corollaries (symmetric spectrum ⇒ zero production *to all orders, at every truncation*); what this implies for any claimed "pair-creation asymmetry" (sets up Ch. 10–11).
- **Pedagogy:** Toolbox "regulating divergent sums: heat kernel and zeta" (full machinery in App. B).
- **Figures:** `ch09_fig1_eta_schematic` (new: lopsided vs mirror spectra, $\eta$ as level-counting); `ch09_fig2_regulator_convergence` (new: $\eta(t)$ extrapolation on an exactly solvable example).
- **Scripts:** `ch09_eta_toy.py` (new, small).
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 8 pp.

---

### Ch. 10 — CP phases in a box: what is physical and what is frame

**Goal.** Full treatment of the CP-violating mass $m\,e^{i\delta\gamma_5}$ in finite volume, organized around a single question: *where can $\delta$ hide?* Bulk phases are eliminated twice — by the **Bulk-Phase Dressing Theorem** (3+1D radial) and by the **Spectral Mirror Theorem** (1+1D exact) — yielding the **Bulk-Phase Zero-Charge Theorem**. Taught as diagnostics per P5: the per-channel asymmetries and the κ-cancellation are dissected as artifacts of an inconsistent identification map (non-orthogonal dressed modes under the flat measure), which is *the same lesson* as the identification-map subtlety that returns in Ch. 19.

- Contents: chiral rotations and why free space erases $\delta$; why walls break chiral symmetry (the bag as finite-volume analogue of the θ-angle story); the honest 1+1D Hamiltonian $H = -i\sigma_1\partial_x + m_R\sigma_3 - m_I\sigma_2$ derived from the Lagrangian step by step (Hermiticity emphasized); eigenvalue equation $\tan(pL) = -p/(m\cos\delta)$ derived in full (the boundary algebra including the $W = E+m_R$ cancellation); both energy branches per root + no-midgap-solution analysis ⇒ Spectral Mirror Theorem ⇒ $\eta\equiv0$ ⇒ zero charge **[Theorem + Computed $10^{-16}$]**; the 3+1D radial system and the dressing map $e^{m_I r}$ with the same-sign observation made transparent; why the dressed modes are non-orthonormal and what the consistent inner product is; the spherical κ-cancellation as an automatic consequence of secret CP-conservation; the lesson box: *every overlap statement is relative to an identification map* — declared a recurring principle of the thesis.
- **Figures:** `ch10_fig1_chiral_rotation_flow` (new schematic: where $\delta$ can live — bulk ↔ walls); `ch10_fig2_radial_similarity` (new from existing check: eigenvalue agreement to $5\times10^{-10}$, wavefunction ratio $=e^{m_I r}$); `ch10_fig3_mirror_spectrum` (new: spectrum vs $\delta$ collapsing onto the $m\cos\delta$ bag).
- **Scripts:** `ch10_radial_similarity.py` (port of `check1_radial_similarity.py`), `ch10_mirror_spectrum.py` (uses `bag1d` module; prints the $\eta = 4\times10^{-16}$ checks).
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 11 pp.

---

### Ch. 11 — The polarization diagnostic: what truncated pair-creation sums measure

**Goal.** Resolve the apparent paradox: a reproducible nonzero "asymmetry" ($A \approx 0.0764$-type numbers, $\propto\sin2\delta$) from Bogoliubov sums versus the Zero-Charge Theorem. Prove the **Polarization Identity** and show, with converged numerics, that the truncated sums chase a *wall-polarization redistribution integral* — finite, odd in $\delta$, and integrating to zero over the whole box. Include the ECPT-oddness check ($\Delta Q(\delta) = -\Delta Q(-\delta)$ to machine precision) — real, but a property of the polarization observable.

- Contents: Bogoliubov coefficients between sectors (recap from Ch. 8/App. A); the same-branch vs cross-branch truncated forms and their honest 20% mutual discrepancy as the smoking gun of non-unitarity at finite cutoff; the Polarization Identity derived via completeness/projector algebra (full steps: $\rho_\pm$, $P'_\pm$, $\operatorname{sgn}(H')$); wall charge layers $\pm q(\delta)$; convergence study (forms close on each other and on the independently computed polarization integral $+0.0672$ while whole-box $\eta\equiv0$); the un-definability of "empty fresh region" for a Dirac field (linearly divergent charge; symmetric cutoff implicitly half-fills) — a short, sharp section; the takeaway rule for the field: *truncated Bogoliubov charge sums never measure net charge; spectral flow always does.*
- **Figures:** `ch11_fig1_wall_charge_density` (**new script**: vacuum charge density $\rho(x)$ profile in the bag at $\delta\neq0$ showing the equal-and-opposite wall layers — the chapter's concept figure, currently missing from the corpus); `ch11_fig2_convergence` (rebuild of the convergence study with thesis styling); `ch11_fig3_ecpt_oddness` (rebuild: $\Delta Q(\delta)+\Delta Q(-\delta)$ at machine zero).
- **Scripts:** `ch11_polarization.py` (port of `quench1d.py` conv study), `ch11_wall_density.py` (**new**: spatial charge density from complete spectra).
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 9 pp.

---

### Ch. 12 — The mechanism: wall-localized CP violation and quantized charge pumping

**Goal.** The constructive payoff of Part II. The two-wall chiral bag solved exactly: **Chiral-Wall Master Equation**, **Wall-Mismatch Criterion**, **Fractional Charge Law**, **Crossing Condition**, **Quantized Pumping Theorem** — each with complete derivation and converged numerics, then the physics: this is anomaly inflow / spectral flow (Callan–Harvey), structurally the electroweak-baryogenesis configuration (CP on bubble walls) reached from the expanding-box side.

- Contents: boundary rays for general wall angles (the projector algebra worked in full, with the MIT case as check); the master equation derived (including the trigonometric collapse via $a\pm b$); reading off everything: $E\to-E \Leftrightarrow \Delta\to-\Delta$ (Wall-Mismatch Criterion), $\Delta=0$ ⇒ pure mass renormalization $m\to m\cos\Sigma$ (one line containing all of Ch. 10); below-gap branch; fractional charge numerics vs $-\Delta/2\pi$; the Goldstone–Wilczek current derivation of the coefficient (Toolbox block on the GW current — short derivation included, not just cited); crossing condition with the side payoff (one MIT wall ⇒ no crossings ever — explains why naive scans find nothing); the pump: quench across $L^*$ gives $\Delta Q = +0.996$, control quench gives exactly 0 while truncated forms still read $-0.06$ (the final nail from Ch. 11); physical reading: CP domain wall through the bag ≡ wall mismatch after chiral rotation; Callan–Harvey anomaly-inflow framing.
- **Pedagogy:** Toolbox "spectral flow and the index" ; Toolbox "Goldstone–Wilczek current in 1+1D".
- **Figures:** `ch12_fig1_qvac_law` (adapt existing: $Q_{\text{vac}}$ vs $\Delta$ against $-\Delta/2\pi$); `ch12_fig2_crossing` (adapt: level diving through zero at $L^*=0.44302$, unit jump of $Q_{\text{vac}}$); `ch12_fig3_spectral_flow_waterfall` (**new**: full spectrum vs $L$ waterfall with the crossing highlighted — the chapter's signature visual); `ch12_fig4_pump_vs_control` (new: quench results, crossing vs control, spectral vs truncated forms side by side).
- **Scripts:** `ch12_chiral_bag.py` (port/merge of `chiral_exact.py` + `chiral1d2.py`; re-verify $5\times10^{-8}$ master-equation check, plateaus $+0.3647/-0.6402$, jump $-0.99994$), `ch12_waterfall.py` (**new**).
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 13 pp.

---

### Ch. 13 — Interlude: baryogenesis background

**Goal.** Everything cosmological the estimate needs: Sakharov's three conditions; sphalerons and $C_{\text{sph}} = 28/79$; entropy dilution $45/(2\pi^2 g_*)$; the electroweak crossover; CKM matrix, Jarlskog invariant, GIM cancellation; why Standard-Model electroweak baryogenesis fails (its two failure modes — equilibrium and CP-magnitude — stated precisely, because the framework repairs exactly one of them and that honesty is strategic).

- **Figures:** `ch13_fig1_sakharov_map` (new schematic: the three conditions and where this framework supplies each); `ch13_fig2_ewbg_failure_modes` (new schematic/table-figure).
- **Scripts:** none.
- **Complexity:** ◇3. **Model: Opus 4.8** (Fable review of the GIM section — it feeds Ch. 14's budget). Length ≈ 9 pp.

---

### Ch. 14 — From pumping to η_B: the corrected estimate

**Goal.** Rebuild the baryon-asymmetry estimate on the correct engine. The kinematic amplitude is no longer a polarization integral but the **Crossing-Density Amplitude**: the expected number of level crossings per iso-energy expansion step, derived in closed form from the Crossing Condition (*new derivation, the thesis's main genuinely new calculation*). Assemble $\eta_B$ with the standard cosmological prefactors, and give the honest CP-coefficient budget (two-loop measure $1/16\pi^2$, GIM suppression, NDA range spanning $\sim$10 orders, enhancement mechanisms E1–E3 as **[Open]** with attack plans), including the spin-counting caution (the two spin polarizations at fixed transverse momentum are C-conjugates; naive multiplicities collapse).

- Contents: from one pump to a cosmological history (step counting per Hubble volume — stated as **[Postulate]** with sensitivity analysis); derivation of the crossing density $\nu(\Delta,\Sigma, m/T, s)$ from $\tanh(mL^*) = -\cos(\Delta/2)/\cos\Sigma$ along the iso-energy trajectory $L_k = (n_0{+}k)L_0/n_0$ (**new**; validated by `ch14_crossing_density.py` against brute-force spectral flow counts); what sets $\Delta\theta$ microscopically (relative phase between boundary condensate and bulk mass matrix — where $\mathcal C_2^{\text{eff}}$ now anchors); isotropic-vs-anisotropic reassessment (anisotropy not *required* for the mechanism — crossings are — but may control magnitude; the rectangular-bag factorisation conjecture stated with status); assembly $\eta_B = (\text{prefactors})\times\nu\times\sin\Delta_{\text{eff}}$-type master formula; uncertainty budget table; comparison with observation $(6.143\pm0.019)\times10^{-10}$; falsifiability discussion.
- **Figures:** `ch14_fig1_pipeline` (new: factor-flow diagram from microphysics to $\eta_B$ — the narrative anchor); `ch14_fig2_crossing_density` (**new**: heatmap of $\nu$ over $(mL,\Delta)$ with the analytic boundary overlaid); `ch14_fig3_eta_budget` (new: NDA range waterfall, observation band).
- **Scripts:** `ch14_crossing_density.py` (**new**, with brute-force validation mode), `ch14_eta_assembly.py` (**new**; supersedes `compute_eta_b*.py` with the corrected kinematics).
- **Complexity:** ◇5 (new derivation + the chapter most likely to be attacked). **Model: Fable 5.** Length ≈ 11 pp.

---

### Ch. 15 — Why geometry must be dynamical: the normalization diagnostic

**Goal.** Close Part II by answering a natural objection/idea in one sweep: could the iso-energy structure act *inside* a fixed-$L$ field theory, e.g. as a mode-dependent normalization that regulates UV divergences? Present the prescription (IMN, unnamed or renamed "mode-volume normalization"), its striking successes ($\langle\phi^2\rangle$ finite in every $d$, loop power-counting $\omega \to \omega - dI$), and then the **IMN No-Continuum Theorem**: correlations die in the thermodynamic limit, the field is a generalized free field violating microcausality, zero-point sums stay divergent. Constructive reading: the instinct "normalization should know the geometry" is correct and is implemented *consistently* by the Geometry-Sector Decomposition of Ch. 3 — high modes decohere into neighbouring sectors dynamically instead of being damped kinematically. The cell-as-cutoff idea returns, correctly housed, in Ch. 22.

- **Figures:** `ch15_fig1_per_mode_decay` (rebuild: $n^{-1}$ vs $n^{-2}$ per-mode contributions); `ch15_fig2_phi2_plateau` (rebuild: 1+1D and 3+1D saturation); `ch15_fig3_which_L` (new schematic: nested boxes and the "which $L$?" ambiguity).
- **Scripts:** `ch15_imn.py` (port of `imn_uv_regularization.py`, trimmed).
- **Complexity:** ◇3 with a ◇4 closing section. **Model: Opus 4.8 draft + Fable 5 writes the No-Continuum Theorem and the constructive-reading section.** Length ≈ 7 pp.

# PART III — EMERGENT GRAVITY

---

### Ch. 16 — Interlude: the general-relativity toolbox

**Goal.** A self-contained, purpose-built GR background chapter — only what Part III consumes, but each item complete. The reader who has never seen ADM or Kasner must be able to follow every later chapter from this one.

- Contents: metric, curvature, Einstein–Hilbert action **[Standard]**; the 3+1 (ADM) split: lapse, shift, extrinsic curvature, and *the Hamiltonian constraint as the statement that the total energy of geometry+matter vanishes* (the framing Ch. 21 needs); the **DeWitt supermetric** and its indefinite signature (conformal mode opposite sign) — flagged: *Part III will derive this sign twice, dynamically*; FRW cosmology and the Friedmann equations, matter/radiation power laws; Bianchi I, the Kasner solution derived in four lines, the Kasner circle $\sum p_a = \sum p_a^2 = 1$ and the $u$-parametrization; Jacobs stiff-matter interior; relational (clock-free) exponents $p_a = \dot\alpha_a/\sum\dot\alpha_b$ and the diagnostic $\Sigma = \sum p_a^2$; Bianchi IX and the BKL map $u\to u-1$, $u\to1/(u-1)$ (statement + meaning; the model's verification is Ch. 24); conformal flatness of 2D metrics; cosmological particle production (Parker) and moving mirrors (Davies–Fulling) — two paragraphs each, anchoring Ch. 17's dictionary; vielbein/triad formalism $g_{ij}=e^a_i e^a_j$ and local frame rotations; one page on Kibble–Sciama: GR as the gauge theory of the frame group (full use in Ch. 23).
- **Pedagogy:** this chapter is deliberately *modular*: each section is labeled with "used in Ch. …" so the reader can defer.
- **Figures:** `ch16_fig1_adm_foliation` (new schematic); `ch16_fig2_kasner_circle` (new: the circle in the $(p_1,p_2)$ plane with $u$-parametrization and special points); `ch16_fig3_bkl_map` (new: cobweb/staircase plot of the $u$-map, chaos visualized).
- **Scripts:** `ch16_kasner_circle.py` (tiny, figure only).
- **Complexity:** ◇3 (volume, not depth). **Model: Opus 4.8** (Fable review of the constraint and DeWitt sections). Length ≈ 15 pp.

---

### Ch. 17 — The dictionary: a box is a metric, expansion is cosmology

**Goal.** The exact change-of-variables that powers all of Part III. **Dictionary Theorem**: pulling the box $[0,L]$ back to the unit interval gives $H = -\partial_\xi^2/2mL^2$ = Laplace–Beltrami with $g_{\xi\xi}=L^2$; "size" vs "metric on fixed coordinates" is pure diffeomorphism. Time-dependent $L(t)$ ⇒ the FRW line element $ds^2 = dt^2 - L(t)^2 d\xi^2$; comoving momenta fixed, frequencies redshift as the iso-energy curves of Ch. 2 — *the foundational observation of Part I was cosmology all along*. The sudden quenches of Part II are step-function scale-factor particle production; in 1+1D the equivalence is total (conformal flatness + massless conformal invariance; $m\Omega$ as the only dial — why everything organized into $mL$).

- Contents: the dictionary derivation; the 3D version $E = \frac{\pi^2}{2m}n_an_b\,g^{ab}$ and its coordinate invariance; what is and is not gauge: $(n,L)\to(\lambda n,\lambda L)$ = global Weyl rescaling of the spectrum (the bridge to Ch. 23); sudden/adiabatic as the standard cosmological mode-splitting (back-reference Ch. 5); the moving-mirror connection.
- **Figures:** `ch17_fig1_dictionary` (new schematic: same physics, two descriptions — wall moves vs metric grows; *signature diagram reused in Ch. 26*); `ch17_fig2_redshift_curves` (rebuild of iso-energy curves *re-labeled as comoving modes redshifting* — the deliberate "you've seen this before" moment).
- **Scripts:** `ch17_dictionary_check.py` (new, small: numerically verifies spectrum equality of moving-wall vs metric formulations).
- **Complexity:** ◇3. **Model: Opus 4.8 draft + Fable 5 rewrites the framing paragraphs** (the rhetorical placement of this chapter matters enormously). Length ≈ 7 pp.

---

### Ch. 18 — A lattice of cells: the triad, the missing shears, and the conformal obstruction

**Goal.** Promote one box to a field of cells. **Metric Emergence Theorem** ($g_{ij} = \sum_a e^i_a e^j_a$ from edge vectors, rank-1 updates from directional transitions, 6 components counted); recognize the edge vectors as a **vielbein/triad** with local $SO(3)$ redundancy; **Flatness Theorem** for single-cell transitions; then the structural diagnosis: directional transitions generate only the Cartan (diagonal) part of $GL(3)/SO(3)$ — no operator changes edge *angles*. Define the iso-energy **shear generators** $\hat S^{(ab)}$ and prove the **Conformal Obstruction Theorem**: dilation-gauging alone yields conformally flat metrics, zero Weyl tensor, hence no light deflection (Nordström's fate) and no gravitational waves — the configuration space itself, not just the action, fails to be generated. Dilations + shears + frame rotations close on $\mathfrak{gl}(3)$: 6 symmetric (metric) + 3 antisymmetric (gauge).

- **Pedagogy:** Toolbox "the Weyl tensor in three sentences"; Toolbox "vielbeins and why 9 = 6 + 3".
- **Figures:** `ch18_fig1_cell_lattice` (new schematic: tiled space, per-cell edge vectors); `ch18_fig2_gl3_decomposition` (new: the $1+5\,(\text{shape})$ vs $3+3$ generator decomposition diagram); `ch18_fig3_conformal_obstruction` (new schematic: conformally flat family vs the Weyl directions; light-bending/GW callouts).
- **Scripts:** `ch18_metric_updates.py` (port of metric-emergence plot logic, restyled).
- **Complexity:** ◇4. **Model: Fable 5.** Length ≈ 9 pp.

---

### Ch. 19 — Shape space: the Splitting Theorem

**Goal.** The quantum geometry of deformations. Define the fidelity (quantum) metric $\chi_n(\delta A)$ with its perturbation-theory equivalent; prove the **Shape-Space Splitting Theorem**: diagonal deformations are exactly overlap-flat but energy-carrying (hence iso-energy motion there must *jump* — the discrete ladder of Part I), while shears are iso-energy at rectangular points and metric-carrying (continuous flow). Read against the DeWitt supermetric: the model *derives* a clock/propagating split, with the honest caveat that the Lorentzian sign is dynamical (delivered in Ch. 22/25), not kinematic. Close with the identification-map principle (the Ch. 10 lesson recurring): fixed-coordinates vs embedding pictures are both legitimate and bound the physics.

- **Figures:** `ch19_fig1_shape_space` (new schematic: the 6-dim metric space, diagonal slice vs shear directions, level sets of $E$); `ch19_fig2_fidelity` (adapt existing: fidelity vs perturbation theory agreement; per-mode shear metric values).
- **Scripts:** `ch19_fidelity.py` (port of `dewitt_fidelity.py`).
- **Complexity:** ◇4–◇5. **Model: Fable 5.** Length ≈ 8 pp.

---

### Ch. 20 — Minisuperspace dynamics: the Kasner Selection Theorem

**Goal.** The first gravity payoff. Derive the continuum (minisuperspace) limit of the geometry ladder — Hamilton's equations for the wavepacket centroid on the tight-binding chain — and prove the **Kasner Selection Theorem**: the embedding-coupled ladder gives ballistic motion in $L$ (Milne-like class; isotropization with matter, anisotropy growth under contraction = the seed of BKL), while coupling sectors through the *dilation Noether generator* gives $t(n) = g(n+\tfrac12)$ (the Ch. 3 Lemma redeemed), conserved $n\cos\theta$, fixed points $\theta=\pm\pi/2$, and exact ballistic motion in $\alpha = \ln L$: the Kasner/Jacobs orbit class with $p_a = g_a/\sum g_b$. Then the quantum simulations: vacuum-Kasner point hit to 3–4 digits with one axis contracting, $\Sigma = 0.9997$ on the circle; equal couplings land on FRW $\Sigma = 1/3$; exponent self-averaging (quantum corrections $\sim \mathrm{Var}(n)/\langle n\rangle^2$, growing toward small $n$ — quantum-gravity effects exactly where expected).

- Contents: lattice dispersion $\omega(n,\theta) = 2t(n)\cos\theta$ and semiclassical equations; case (a) embedding law; relational exponents → $1/3$; case (b) generator law derivation (matrix-element integral in full); what neither law fixes ($\sum p_a^2$ = the constraint, deferred to Ch. 21); simulation methodology (genuine wavepackets, sparse propagation, no semiclassical shortcut) and results; animations referenced.
- **Figures:** `ch20_fig1_two_couplings` (new schematic: $L$-linear vs $\alpha$-linear orbits with the ladder diagram from Ch. 3); `ch20_fig2_relational_exponents` (adapt); `ch20_fig3_kasner_plane` (adapt: both dynamics on/inside the circle); `ch20_fig4_orbit_geometry` (adapt). Animations: `ch20_anim_kasner.gif`, `ch20_anim_isotropization.gif` (reuse, re-rendered with thesis styling).
- **Scripts:** `ch20_kasner_sim.py` (port of `kasner_sim.py`), `ch20_make_anims.py`.
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 11 pp.

---

### Ch. 21 — Closing the constraint: Friedmann universes and emergent attraction

**Goal.** Impose the one imported structural postulate of Part III — the Hamiltonian constraint $\mathcal H_{\text{grav}} + \mathcal H_{\text{matter}} = 0$ (its microscopic derivation is **[Open]**, Ch. 27) — and show the ladder becomes a Friedmann universe: coupling slaved to density, $g(\alpha) = g_0 e^{-\frac{3(1+w)}{2}(\alpha-\alpha_0)}$, exact $t^{2/3}$ and $t^{1/2}$ laws, quantum simulation agreeing to $|\Delta\alpha| \le 0.009$ over two e-folds. Then inhomogeneity at leading order: the separate-universe cell chain, the **Cell-Growth Theorem** ($\ddot\delta + 2H\dot\delta = 4\pi G\bar\rho\,\delta$ derived purely by comparing neighbouring cells), $\delta \propto a$ growth, turnaround at $\delta_{\text{lin}} = 1.062$, collapse at $1.686$, and the 41-cell simulation showing a potential well forming and walls drifting toward the overdensity: *matter attracting matter inside the model*.

- **Pedagogy:** Toolbox "the separate-universe approximation"; Toolbox "spherical collapse thresholds" (numbers derived, not quoted).
- **Figures:** `ch21_fig1_constraint_budget` (new schematic: geometry energy negative, matter positive, total zero); `ch21_fig2_friedmann` (adapt: quantum ladder vs analytic, exponents → $2/3, 1/2$); `ch21_fig3_collapse` (adapt: $\delta(a)$ with linear law and thresholds; the well profile). Animation: `ch21_anim_collapse.gif` (reuse).
- **Scripts:** `ch21_friedmann.py`, `ch21_collapse.py` (ports of `newton_sim.py` parts 2–3).
- **Complexity:** ◇4. **Model: Fable 5** (the constraint-import framing and growth derivation are load-bearing; the GR-side recaps lean on Ch. 16). Length ≈ 10 pp.

---

### Ch. 22 — The stiffness from matter loops: induced gravity and the Planck cell

**Goal.** Remove the last free dial. Top-down: Sakharov's induced gravity — heat-kernel derivation of $\Gamma \supset c_1\Lambda_{\text{UV}}^2 R$ per species, $G_{\text{ind}}^{-1} = N_{\text{eff}}\Lambda^2_{\text{UV}}/12\pi$ — with the cutoff *identified* as the cell scale ($\Lambda_{\text{UV}} = c/\bar L$: the heat kernel stops being continuum when the probe resolves cells), giving the **Planck-Cell Prediction** $\bar L \approx 1.6\,c\,\ell_P$. Bottom-up (1D, fully computable): the Casimir energy $E_C = -\pi/24L$ derived (zeta and heat-kernel, both shown); walls attract; occupied cells have equilibria, empty cells clump (**vacuum Jeans instability**); the chain-relaxation simulation; the measured induced action of the displacement field with $\kappa^{\text{vac}}_{\text{grad}} \approx -0.13 < 0$ vs matter $+3\pi^2n^2$ — the *sign* delivering, at once: gravity attracts, the conformal mode has wrong-sign kinetic energy (Euclidean conformal-factor problem as an *output*), matter stiffens geometry at short scales. The cosmological-constant term stated honestly every time it appears.

- **Pedagogy:** Toolbox "zeta-function regularization in one page" (cross-ref App. B); Toolbox "Sakharov's 1967 idea in plain words".
- **Figures:** `ch22_fig1_induced_logic` (new schematic: matter loops → geometry action, the "elasticity of the vacuum"); `ch22_fig2_cell_energetics` (adapt: $E(L)$ curves, equilibria, chain relaxation endpoint); `ch22_fig3_dispersion` (adapt: $\Pi(q)/\hat q^2$ measurement, negative vacuum band vs positive matter line); `ch22_fig4_planck_loop` (new schematic: $\bar L \leftrightarrow M_P$ self-consistency loop).
- **Scripts:** `ch22_casimir_chain.py` (port of `newton_sim.py` part 1), `ch22_dispersion.py` (port of `stiffness_dispersion.py`).
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 11 pp.

---

### Ch. 23 — Weyl gauging, the compensator identified, and the road to general relativity

**Goal.** The gauge-theoretic consolidation. From the spectrum's global Weyl rescaling (Ch. 17) through honest gauging: Einstein–Hilbert is *not* Weyl-invariant; the unique 2-derivative Weyl-invariant action needs a weight-$(-1)$ compensator $\sigma$; gauge-fixing $\sigma = \sigma_0$ yields EH + massive Weyl vector + cosmological constant with $M_P^2 = \sigma_0^2/6$ **[Standard, fully derived]**. Then the framework's two original contributions: the **Compensator Identification** $\sigma \sim 1/L_{\text{cell}}$ (exactly weight $-1$ under $(n,L)\to(\lambda n,\lambda L)$; "broken-Weyl vacuum" = "finite mean cell size"; consistency with Ch. 22's $M_P^2 \sim N_{\text{eff}}/\bar L^2$ fixes the $\mathcal O(1)$ constant), and the **Kibble–Sciama completion**: gauging the full frame group on the triad the model *derives* (Ch. 18) delivers Einstein–Cartan = GR for spinless sources; dilations survive as the determinant of $GL(3)$ — the Weyl sector is the trace of the larger gauge structure.

- Contents: gauging walkthrough (covariant derivative, $W_\mu$, why pure-$W$ invariants start at four derivatives/conformal gravity); compensator action term-by-term invariance check; gauge fixing and the mass spectrum; what is derived vs structural (table with P4 labels); the identification section; Kibble–Sciama in the model's variables; generator counting ($\mathfrak{gl}(3)$: nothing missing, nothing extra); the honest remainder (Lorentzian signature/lapse-shift; EH dominance over curvature-squared; $\Lambda$) cross-referenced to Ch. 27.
- **Figures:** `ch23_fig1_gauging_ladder` (new schematic: global Weyl → local Weyl → compensator → EH, with the two framework-supplied inputs highlighted); `ch23_fig2_compensator_id` (new schematic: $\sigma = c_\sigma/L_{\text{cell}}$ and the two independent routes to $M_P$ meeting).
- **Scripts:** none (structural chapter).
- **Complexity:** ◇4 (standard constructions + critical identifications). **Model: split — Opus 4.8 writes the standard Weyl/compensator construction to spec; Fable 5 writes the identification and Kibble–Sciama sections and the framing.** Length ≈ 10 pp.

---

### Ch. 24 — Newton's law and BKL chaos on the cell lattice

**Goal.** Spend the stiffness on the two most recognizable strong tests. Static weak field: linearize the constraint-closed cell dynamics → **Lattice Poisson Result** $\kappa\hat\nabla^2\delta\alpha = -c\bar\rho\,\delta$ with $\delta\alpha = -\Phi_{\text{Newton}}$ (cells *larger* near mass = the weak-field spatial metric), measured on a $128^3$ lattice: potential exponent $-1.008$, amplitude $96\%$ of $1/4\pi$, binding $-0.994$, force $-2.05$. Singularity approach: Bianchi IX with the $S^3$ curvature walls (functional form imported — stated honestly), Kasner flights caroming: the measured $u$-sequence $4.3 \to 3.3 \to 2.3 \to 1.3 \to 10/3$ matching the **BKL map** to $6\times10^{-4}$ *including era reversal*, constraint conserved to $10^{-13}$.

- **Figures:** `ch24_fig1_newton` (adapt: $1/r$ profile, $1/d$ binding, $1/d^2$ force); `ch24_fig2_well_schematic` (new: cells dilated around a mass = potential well); `ch24_fig3_mixmaster` (adapt: $\alpha_i(\tau)$ flights + Kasner-circle hopping). Animation: `ch24_anim_mixmaster.gif` (reuse).
- **Scripts:** `ch24_newton_lattice.py`, `ch24_mixmaster.py` (ports).
- **Complexity:** ◇3–◇4. **Model: split — Opus 4.8 writes the BKL exposition (leaning on Ch. 16) and simulation reportage to spec; Fable 5 writes the Poisson linearization derivation.** Length ≈ 9 pp.

---

### Ch. 25 — The graviton, measured

**Goal.** The decisive computation, told as the decisive computation. Metric perturbations split $1_{\text{conf}} + 2_{TT} + 3_{\text{diffeo}}$; GR's three structural claims (TT stiff and massless; conformal opposite sign; diffeos free). The obstacle as physics: a fixed lattice is a non-covariant cutoff — the diffeo zero-test *fails* ($\Pi \approx 0.42$) for the raw vacuum sum — and the cure is the covariant proper-time regulator $E(s) = \tfrac12\sum\omega_n e^{-s\omega_n}$ with its measured covariance window. Parameter-free prediction via Weyl counting: induced action $C(s)\int\sqrt g R^{(3)}$ with $C(s) = 1/48\pi^2 s^2$; quadratic forms worked for each sector (TT: $-\tfrac14\int(\partial h)^2$; conformal: opposite sign). **DeWitt-Signature Measurement**: $\kappa_{TT+}/C = -0.987$, $\kappa_{TT\times}/C = -1.008$ at $s=3$ (polarizations agreeing to 2% — earned, the lattice has no $45°$ rotation), conformal sign opposite (magnitude contaminated by $a_2$ terms, trending right — honesty scorecard), diffeo control bounding systematics, no induced mass. Static→dynamical bookkeeping: covariance supplies $\dot h^2$ with matched coefficient ⇒ $\omega = |k|$: a massless two-polarization wave. The model was given five ways to fail structurally and passed each.

- **Figures:** `ch25_fig1_sector_split` (new schematic: the six metric directions sorted into conf/TT/diffeo with predicted signs); `ch25_fig2_covariance_window` (new from data: diffeo $\Pi$ vs $s$ collapsing, physical signals persisting); `ch25_fig3_kappa_convergence` (adapt: $\kappa/C(s) \to -1$ for both polarizations); Animation: `ch25_anim_graviton.gif` (reuse: the TT wave crossing the cell lattice).
- **Scripts:** `ch25_graviton.py` (port of `graviton.py` + figure scripts; re-verify the table at $s = 2, 2.5, 3$).
- **Complexity:** ◇5. **Model: Fable 5.** Length ≈ 11 pp.

---

### Ch. 26 — Synthesis: three faces of one structure

**Goal.** Fuse the parts. In metric language, the Ch. 12 mechanism is **anomaly inflow on a moving boundary of the emergent spacetime**: charge is pumped when the boundary Dirac operator's spectral flow crosses zero as the scale factor evolves, the wall angles playing the Goldstone–Wilczek role. The expansion arrow (Part I), the baryon asymmetry (Part II), and gravity (Part III) are three faces of *matter on a dynamical geometry whose boundaries carry the CP physics*. Then the honest scorecard: the unbroken chain table (each link: result, strength, where proved/computed), what would falsify the program (e.g. TT polarizations disagreeing, conformal sign parallel, BKL map violated, crossing-density mechanism failing in 3+1D), and what the theory does *not* claim (cosmological constant; Lorentzian signature derivation; $\mathcal C_2^{\text{eff}}$ magnitude).

- **Figures:** `ch26_fig1_three_faces_closed` (the Ch. 1 schematic, now with every arrow labeled by a theorem name and chapter); `ch26_fig2_chain_table` (rendered table-figure of the verification chain with numerical strengths).
- **Scripts:** none.
- **Complexity:** ◇4. **Model: Fable 5.** Length ≈ 6 pp.

---

# PART IV — OUTLOOK

---

### Ch. 27 — Open problems, with attack plans

**Goal.** The research program, ordered by tractability, each item with: precise statement, why it matters, concrete first calculation, and falsification potential.

1. **Microscopic Hamiltonian constraint** — both sides of the budget are now computed objects; test conservation under ladder dynamics (finite calculation).
2. **$N_{\text{eff}}$ and the lattice constant $c$** — the spin-weighted species sum and the $s \leftrightarrow \bar L$ relation; converts induced coefficients into the physical Newton constant and sharpens the Planck-Cell Prediction.
3. **3+1D spherical chiral bag** — κ-resolved spectral asymmetry and crossing density with mismatched wall angles (Goldstone–Jaffe machinery applies); decides whether anisotropy controls magnitude.
4. **$\mathcal C_2^{\text{eff}}$ and boundary GIM evasion** — anchor the wall-angle mismatch $\Delta\theta$ in SM physics (boundary condensate vs bulk mass-matrix phase); the 10-order NDA gap.
5. **Lorentzian signature, lapse and shift** — the construction is spatial; time entered through transition dynamics; the measured negative compression stiffness is the dynamical seed — formalize it.
6. **Dynamical response $\Pi(\omega, q)$** — measure (not infer) the $\dot h^2$ coefficient and wave speed; also pins the conformal magnitude (larger lattices, smaller $qs$).
7. **Bianchi IX from intercell coupling** — replace the imported $S^3$ wall functional form with the model's own gradient coupling.
8. **Adiabatic dressing of the pump** — Landau–Zener corrections exactly at crossings, where the sudden approximation is weakest.
9. **The cosmological constant** — the $V/s^4$ term: stated, not solved; the one structure untouched.
10. **Lorentz-violation phenomenology** — if cell granularity is physical at $\bar L \sim \ell_P$, what dispersion signatures follow (kept honest and modest).

- **Figures:** `ch27_fig1_program_map` (new: open problems as a dependency graph with "one number: κ"-style annotations).
- **Complexity:** ◇4. **Model: Fable 5.** Length ≈ 7 pp.

---

# APPENDICES

- **App. A — Bogoliubov transformations in finite volume.** Scalar and Dirac; sector embeddings; α/β algebra; unitarity and its truncation breakdown (the technical backbone of Ch. 11). ◇2–◇3, **Opus 4.8** (Fable review). ≈ 6 pp.
- **App. B — Heat kernel, zeta function, Seeley–DeWitt.** Regulated sums, $\eta(t)$, Weyl counting $N(\omega)$, the $a_1$ coefficient used in Ch. 22/25, Richardson extrapolation. ◇3, **Opus 4.8** (Fable review). ≈ 6 pp.
- **App. C — Numerical methods and reproducibility.** Transfer-matrix solvers; complete-spectrum strategies; sparse `expm_multiply` wavepacket evolution; spectral derivatives and block-diagonalization on $63^3$ lattices; convergence/extrapolation standards; the `run_all.py` contract (every quoted number printed by exactly one script). ◇3, **Fable 5** (it owns all scripts). ≈ 5 pp.
- **App. D — Notation and glossary.** Global symbol table; the theorem-name registry (mirror of skeleton §2); index of figures. ◇1, **Opus 4.8**. ≈ 4 pp.

**Estimated totals:** 27 chapters + 4 appendices ≈ **230–250 pages**, ~75 figures (≈ 25 new schematics, ≈ 30 rebuilt/restyled, ≈ 14 adapted v14 plots, 5 animations), ≈ 35 scripts (≈ 22 ports, ≈ 13 new).

---

## 5. Coverage matrix (every legacy result → its v15 home)

Internal verification that nothing is lost and every superseded claim is correctly re-housed. (A fuller per-equation version will live in `v15/coverage_matrix.md`.)

| Legacy result / claim | Latest status | v15 home |
|---|---|---|
| Iso-energy curves, dilation homogeneity, selection rule | established | Ch. 2 |
| Hilbert space, $\hat T_\pm$, telescoping, geometric overlaps ($\sqrt{n/(n+m)}$, sinc) | established | Ch. 3 |
| Dilation generator $\hat D$, Ward identity, matrix elements | established (promoted: selects GR) | Ch. 3 → Ch. 20 |
| Tight-binding chain, band, Bessel dynamics, geometry superpositions, low-$n$ symmetry breaking | established | Ch. 4 |
| Sudden/adiabatic regimes, adiabatic parameter | established | Ch. 5 |
| Natural boundary at $n=1$; self-adjointness | established | Ch. 2 → Ch. 7 (E-breaking) |
| Norm theorems (expansion = 1; KE preservation; contraction deficit; kinetic divergence; resonant exceptions) | established | Ch. 6 |
| Multi-particle generator commutativity; path-dependence of physical amplitudes | established | Ch. 6 |
| Probability rule + factorization (paper 1 Eq. 13) | **upgraded to theorem** (golden rule in sector arena) | Ch. 6 |
| $N$-spectator exponential asymmetry, $N_{\text{crit}}$ | established (caveats kept) | Ch. 6 |
| $\hat E$, E-breaking (kinematic + dynamical), ECPT group, conjecture, two arrows | conjecture (framing tightened) | Ch. 7 |
| MIT bag spectrum, boundary-value suppression, massless/NR limits | standard | Ch. 8 |
| Second-quantised Dirac in a box, Fock space, Bogoliubov setup (incl. pair creation from occupied states) | standard/established | Ch. 8 + App. A |
| Relativistic iso-energy lattice (Dirac iso-energy condition; particle/antiparticle ladders; the gap) | established | Ch. 8 |
| Dirac spectator overlaps & contraction norms; CP phase does not weaken expansion preference; $f_{\text{neq}} \approx 1$ | established | Ch. 8 + Ch. 14 |
| Symmetrised charge operator, $Q_{\text{vac}} = -\eta/2$ | established | Ch. 8–9 |
| **Spectral-Flow Master Theorem** | proven | Ch. 9 |
| CP phase physical in a bag (chiral-symmetry breaking by walls) | standard/established | Ch. 10 |
| Slab eigenvalue theorem $\tan(pL) = -p/(m\cos\delta)$ (transverse momenta cancel) | established (vindicated) | Ch. 10 |
| Bulk-CP radial prescription; per-channel $\Delta Q_\kappa \neq 0$ | **trivialized** (dressing theorem; frame artifact) | Ch. 10 (diagnostic) |
| Spherical κ-cancellation theorem | true but contentless (secret CP conservation) | Ch. 10 (diagnostic) |
| v13.4 §46 opposite-sign radial system, CP-Coulomb $-2\kappa m_I/r$ | inconsistent system — **retired**; lesson kept as identification-map principle | Ch. 10 (lesson box); §53 gravity chain **not carried over** (superseded by Part III) |
| $\Delta Q \propto \sin 2\delta$ numerics, $A = 0.0764$, ECPT oddness check | reproducible but = polarization, not charge | Ch. 11 |
| **Polarization Identity**; empty-region non-definability | proven | Ch. 11 |
| Chiral two-wall bag: master equation, $\sin(\Delta/2)$ criterion, $-\Delta/2\pi$ law, crossing condition, quantized pump + control | proven/computed | Ch. 12 |
| Anisotropic escape (paper 3), rectangular factorisation conjecture, spin-counting collapse | re-assessed: anisotropy optional, crossings essential; conjecture stated | Ch. 14 |
| $\eta_B = 6.72\times10^{-10}\,\mathcal C_2^{\text{eff}}$ structural formula; GIM/NDA budget; E1–E3 | rebuilt on **Crossing-Density Amplitude** (new) | Ch. 14 |
| IMN: $\langle\phi^2\rangle$ finite, loop power counting; Casimir negative result; CCR breakdown; Lorentz violation | reinterpreted as no-continuum diagnostic → dynamical $L$ | Ch. 15 |
| Box = metric dictionary; $E = \frac{\pi^2}{2m}n_an_b g^{ab}$; FRW from $L(t)$; Parker/moving mirrors | established | Ch. 17 |
| Metric Emergence Theorem; commutativity/flatness; Weyl rescaling of $g_{ij}$ | established | Ch. 18 (+ Ch. 17) |
| Missing shears; $GL(3)$ closure; conformal-flatness obstruction (Nordström/GW argument) | established | Ch. 18 |
| Shape-space splitting theorem (fidelity metric) | proven/computed | Ch. 19 |
| Minisuperspace limit; coupling-selection theorem; Kasner sims ($\Sigma = 0.9997$, contracting axis); isotropization | proven/computed | Ch. 20 |
| Hamiltonian constraint closure → Friedmann (sub-percent); growth equation; turnaround/collapse sims | computed (constraint = imported postulate) | Ch. 21 |
| 1D Casimir stiffness; vacuum Jeans; chain relaxation; measured $\kappa^{\text{vac}}_{\text{grad}} < 0$; matter stiffness | derived/measured | Ch. 22 |
| Sakharov heat-kernel $G_{\text{ind}}^{-1} = N_{\text{eff}}\Lambda^2/12\pi$; cell-scale cutoff; $\bar L \approx 1.6\,c\,\ell_P$ | derived (finite sums open) | Ch. 22 |
| Weyl gauging corrected story; compensator action; EH + massive $W_\mu$ + $\Lambda$; $M_P^2 = \sigma_0^2/6$ | standard, fully derived | Ch. 23 |
| Compensator identified as $1/L_{\text{cell}}$; Kibble–Sciama completion | framework results | Ch. 23 |
| Lattice Poisson; $1/r$, $1/d$, $1/d^2$ measurements | computed | Ch. 24 |
| Bianchi IX / BKL map verification incl. era reversal | computed (wall form imported) | Ch. 24 |
| Diffeo-breaking lattice cutoff + covariant regulator; $C(s) = 1/48\pi^2s^2$; $\kappa_{TT} \to -C$; DeWitt signature; masslessness | proven/measured | Ch. 25 |
| "Three faces" unification (anomaly inflow on moving boundary) | synthesis | Ch. 26 |
| All open items (constraint origin, $N_{\text{eff}}$, 3+1D chiral bag, $\Lambda$, signature, $\Pi(\omega,q)$, …) | open, with attack plans | Ch. 27 |

**Deliberately dropped (with reason):** the v13.4 §53 CP-Coulomb→galaxy-rotation speculation (rested on the inconsistent radial system *and* had the wrong scaling — superseded by Part III's complete gravity program); v12-era claims already corrected in v13 (never resurface); Bessel-function overlap claims (corrected to elementary forms in Ch. 3).

---

## 6. Production plan

**Phase 0 — scaffolding (Fable 5).** Create `thesis/00_master.md` (title, abstract, reading guide, TOC), `appD` notation table, `scripts/style.mplstyle`, `run_all.py` skeleton. *Everything else depends on the notation file.*

**Phase 1 — Part I (Ch. 2–7).** Order: 2 → 3 → 4 → 5 → 6 → 7 (strict dependency chain). Scripts ported/validated *before* each chapter is written, so the text quotes verified numbers only.

**Phase 2 — Part II (Ch. 8–15).** Order: 8 → 9 → 10 → 11 → 12 → 13 → 14 → 15. The `bag1d` module (Ch. 8 script) is the shared engine for 10–12; build and test it first. Ch. 14's new crossing-density derivation is the riskiest single item in the thesis: prototype the script before writing the chapter.

**Phase 3 — Part III (Ch. 16–26).** Order: 16 → 17 → 18 → 19 → 20 → 21 → 22 → 23 → 24 → 25 → 26. Heavy simulation ports (20, 21, 22, 24, 25) re-run end-to-end once; cached outputs go to `v15/data/` so figure regeneration is cheap.

**Phase 4 — Ch. 1, Ch. 27, App. A–C, final pass.** The introduction is written *last* (it summarizes verified content). Fable 5 performs a full-thesis consistency pass: theorem-name registry enforced, no version references, every number traced to a script, every P1 announcement actually honored.

**Validation gates (per chapter):** (G1) all scripts run from `run_all.py` and print the quoted numbers; (G2) figures render and are referenced in text in order; (G3) honesty labels present on every major claim; (G4) Fable 5 review sign-off for Opus-authored chapters.

**Model workload summary:**

| Model | Chapters (author) | Share |
|---|---|---|
| **Fable 5** | 1, 3, 6, 7, 9, 10, 11, 12, 14, 18, 19, 20, 21, 22, 25, 26, 27, App. C + all scripts + all reviews + split-sections of 8, 15, 17, 23, 24 | ≈ 60% of pages, 100% of ◇4–◇5 |
| **Opus 4.8** | 2, 4, 5, 13, 16, App. A, B, D + main drafts of 8, 15, 17, 23 (spec'd sections), 24 (BKL exposition) | ≈ 40% of pages, ◇1–◇3 only |

---

*End of skeleton. Next step on approval: Phase 0 scaffolding, then Part I.*


