# Chapter 13 — Interlude: baryogenesis background

---

Chapter 12 built a pump; Chapter 14 must install it in the early universe. Between them sits a body of standard physics — what the asymmetry *is*, what any mechanism must supply (Sakharov), how produced charge becomes baryons (sphalerons), how it dilutes to today's number (entropy), and why the Standard Model's own CP violation is famously insufficient (GIM). This interlude derives each piece at the depth Ch. 14 consumes it, with the two numerical workhorses — $28/79$ and $45/2\pi^2 g_*$ — audited rather than quoted. Everything here is **[Standard]**; the only thesis-specific content is the final table saying which lines the framework rewrites.

## 13.1 The observation to explain

The asymmetry is conventionally quoted as the baryon-to-photon ratio, measured two independent ways. The cosmic microwave background's acoustic peaks fix the baryon density (the second-to-first peak ratio is the sensitive dial), giving $\eta_B = (6.143 \pm 0.019)\times10^{-10}$ (Planck 2018). Big-bang nucleosynthesis reads the same number from primordial abundances — deuterium above all, whose survival is exponentially sensitive to the baryon density at $T \sim 0.1$ MeV — and concords. Two epochs separated by 400{,}000 years of cosmic history, one number: $\eta_B$ is among the most secure quantities in cosmology.

Why is *nonzero but tiny* the puzzle, in both directions? Suppose the universe were exactly symmetric. Nucleons and antinucleons stay in equilibrium until annihilation $N\bar N \to$ pions outruns the expansion; freeze-out happens at $T \approx 20$ MeV, when the equilibrium abundance carries the Boltzmann penalty $e^{-m_N/T} \approx e^{-47}$. The relic that survives annihilation is set by (annihilation cross-section) $\times$ (Hubble rate), and the standard freeze-out estimate lands at

$$\frac{n_N}{n_\gamma}\bigg|_{\text{symmetric}} \sim 10^{-18}\text{–}10^{-19}$$

— *eight to nine orders of magnitude* too little matter, and equal antimatter alongside it (which the absence of annihilation gamma-ray signatures from galaxy clusters independently excludes at any interesting scale). A symmetric universe is not approximately our universe; it is empty. Conversely $\eta_B \sim 1$ would have rewritten nucleosynthesis and structure formation entirely. The number $6\times10^{-10}$ is a *target*, not an impression: one excess quark per ~1.6 billion quark–antiquark pairs at the primordial epoch.

## 13.2 Sakharov's three conditions

Any dynamical explanation must clear three hurdles (Sakharov, 1967). They are usually listed; they are more useful derived, because each derivation names the loophole a mechanism must thread.

**(1) Baryon-number violation.** Trivial direction: if $[\hat B, \hat H] = 0$, then $\dot B = 0$ and an initially symmetric state stays symmetric. The content is in the converse hunt: the violation need not be in the Lagrangian's perturbative vertices — the Standard Model's own $B$-violation (§13.3) is non-perturbative, and this framework's (Ch. 9, 12) is spectral flow under geometry change. *Loophole named:* "$B$-violation" means any dynamics whose charge bookkeeping is not invariant — boundaries count.

**(2) C and CP violation.** Suppose C were exact: for every history $h$ producing $\Delta B$, the conjugate history $Ch$ produces $-\Delta B$ *at the same rate* (the rates are matrix elements of a C-commuting evolution between C-conjugate states), so the ensemble average vanishes. Now suppose C is broken but CP exact: left-handed quarks may outproduce left-handed antiquarks, but CP maps the excess into the mirrored right-handed sector with opposite sign — summed over helicities, $\langle\Delta B\rangle = 0$ again. *Both* C and CP must break for the rate asymmetry to survive helicity summation. *Loophole named:* the breaking can sit in spatial structure (walls, interfaces) rather than in bulk couplings — Ch. 12's mechanism lives exactly there.

**(3) Departure from equilibrium.** The cleanest version is a three-line CPT argument. In thermal equilibrium, $\langle B\rangle = \mathrm{Tr}\,[e^{-\beta H}\hat B]/Z$. CPT invariance of $H$ supplies an antiunitary $\Theta$ with $\Theta H\Theta^{-1} = H$ and $\Theta\hat B\Theta^{-1} = -\hat B$ (baryon number is C-odd, PT-even). Then

$$\langle B\rangle = \frac{1}{Z}\mathrm{Tr}\big[\Theta\,e^{-\beta H}\hat B\,\Theta^{-1}\big] = \frac{1}{Z}\mathrm{Tr}\big[e^{-\beta H}(-\hat B)\big] = -\langle B\rangle = 0.$$

Equilibrium *erases* asymmetries as fast as anything writes them; generation must outrun thermalization. *Loophole named:* the departure can be a transient (a phase transition's bubble walls) or — the structural option this framework exercises — a *permanent* arrow in the dynamics (Ch. 6).

## 13.3 Sphalerons and the conversion factor $28/79$

The Standard Model violates $B$ on its own. The electroweak theory's chiral couplings make the $B{+}L$ current anomalous,

$$\partial_\mu j^\mu_{B+L} = \frac{3\,g^2}{16\pi^2}\,\mathrm{Tr}\,W_{\mu\nu}\tilde W^{\mu\nu} \quad (\text{while } B - L \text{ stays exact}),$$

so transitions between topologically distinct gauge vacua shift $B$ and $L$ by three units each (one per generation). At zero temperature such transitions tunnel through a barrier — the **sphaleron**, the saddle-point configuration, of height $E_{\text{sph}} \sim 8$–$10$ TeV — at the absurd rate $e^{-16\pi^2/g^2} \sim 10^{-160}$: protons are safe. Above the electroweak crossover, however, the barrier is *thermally stridden* rather than tunneled: the rate per volume $\sim \alpha_W^5 T^4$ is fast compared with the Hubble rate for all $T \gtrsim 130$ GeV. **Hot sphalerons are an equilibrium process**, with two consequences in opposite directions: they will gladly *erase* a pure $B{+}L$ asymmetry (a constraint on careless mechanisms), and they will faithfully *convert* any asymmetry they cannot erase — anything stored in $B{-}L$, or delivered as a generic charge excess — into baryons.

The conversion ratio is fixed by chemical equilibrium, and auditing it is a half-page of counting. Assign a chemical potential to each Standard-Model species; impose: sphalerons ($\sum_{\text{gen}}(3\mu_{q_L} + \mu_{\ell_L}) = 0$), Yukawas (relate left and right), gauge interactions (equalize within multiplets), and hypercharge neutrality of the plasma. Solving the linear system for three generations and one Higgs doublet expresses $B$ and $L$ in terms of the conserved $B - L$:

$$B \;=\; \frac{8 N_g + 4 N_H}{22 N_g + 13 N_H}\,(B - L) \;\;\xrightarrow{\;N_g = 3,\; N_H = 1\;}\;\; \boxed{\;C_{\text{sph}} = \frac{28}{79} \approx 0.354.\;}$$

The numerator and denominator are exactly the advertised audit: $8(3) + 4(1) = 28$; $22(3) + 13(1) = 79$. Freeze-out: sphalerons shut off when the Higgs condensate grows through the crossover ($T_* \approx 130$–$160$ GeV), after which $B$ is frozen. Chapter 14 books $C_{\text{sph}}$ as the price of converting the pump's charge excess into baryon number.

## 13.4 Entropy dilution: from produced charge to today's ratio

A mechanism computes an asymmetry *per comoving volume at production*; observation reports it *per photon today*. The exchange rate is entropy. In an adiabatically expanding universe the entropy density $s$ scales exactly as the comoving volume's inverse, so $n_B/s$ is frozen from production (post-sphaleron) to now. The bookkeeping:

$$s = \frac{2\pi^2}{45}\,g_*\,T^3, \qquad n_\gamma = \frac{2\zeta(3)}{\pi^2}\,T^3,$$

with $g_*$ the effective relativistic degrees of freedom. At the electroweak epoch the full Standard Model is relativistic; the audit, species by species (bosons count 1 per polarization, fermions $7/8$):

| sector | counting | $g_*$ |
|---|---|---|
| gluons | $8 \times 2$ | $16$ |
| electroweak gauge + Higgs | $W^\pm, Z, \gamma$ pre-breaking: $4\times2$; Higgs doublet: $4$ | $12$ |
| quarks | $6$ flavors $\times 3$ colors $\times 2$ spins $\times 2$ (q, $\bar q$) $\times \tfrac78$ | $63$ |
| charged leptons | $3 \times 2 \times 2 \times \tfrac78$ | $10.5$ |
| neutrinos | $3 \times 1 \times 2 \times \tfrac78$ | $5.25$ |
| **total** | | $\mathbf{106.75}$ |

Hence the dilution factor from "charge per entropy" to "baryons per photon today": $\eta_B = (n_B/s)\,\times\,s/n_\gamma|_{\text{today}}$, conventionally packaged (accounting also for $e^+e^-$ heating of photons after neutrino decoupling) as the prefactor Ch. 14 uses:

$$\frac{45}{2\pi^2 g_*} \;=\; 0.02136 \qquad (g_* = 106.75).$$

The structural moral: between production and observation the asymmetry is *diluted by everything else that was hot* — about one-and-a-half orders of magnitude — and this factor is not negotiable by model-building.

## 13.5 CP violation in the Standard Model: CKM, Jarlskog, and GIM

Now the heart of the difficulty, and of Ch. 14's budget. The Standard Model's only established CP violation (setting aside the strong-CP angle, empirically $\lesssim 10^{-10}$) lives in the CKM matrix: the misalignment between the unitary rotations diagonalizing up-type and down-type Yukawas leaves one physical phase in the charged-current mixing matrix $V_{\text{CKM}}$. Its reparametrization-invariant strength is the **Jarlskog invariant**

$$J_{CP} = \mathrm{Im}\big[V_{us}V_{cb}V^*_{ub}V^*_{cs}\big] \approx 3\times10^{-5},$$

— invariant because any phase redefinition of quark fields cancels in this quartet; *every* CP-odd observable built from CKM structure is proportional to it. Already a $3\times10^{-5}$ tax. But the deeper suppression is structural:

**The GIM mechanism, stated carefully.** $J_{CP}$ multiplies a quartet of CKM elements — which means any CP-odd *flavor-singlet* quantity (and a thermal expectation value in the early-universe plasma is flavor-summed, hence singlet) requires the interference of at least **four** charged-current insertions traversing the up- and down-type mass spectra. Unitarity of $V_{\text{CKM}}$ then enforces antisymmetrized sums over internal flavors: if any two same-type quark masses coincide, the relevant unitarity sum collapses and the observable vanishes identically. The leading CP-odd invariant is the famous Jarlskog determinant structure

$$\mathcal A_{CP} \;\propto\; J_{CP}\,(m_t^2 - m_c^2)(m_t^2 - m_u^2)(m_c^2 - m_u^2)\,(m_b^2 - m_s^2)(m_b^2 - m_d^2)(m_s^2 - m_d^2),$$

twelve powers of mass that must be made dimensionless by the only scale available in the hot plasma — the temperature. At $T \sim T_{EW} \sim 150$ GeV, even keeping the top factors at $\mathcal O(1)$, the light-mass-difference factors contribute at the scale of

$$\mathcal S_{\text{GIM}} \;\sim\; \frac{(m_b^2 - m_s^2)(m_c^2 - m_u^2)}{T^4} \;\sim\; \frac{(18\,\text{GeV}^2)(1.6\,\text{GeV}^2)}{(150\,\text{GeV})^4} \;\sim\; 10^{-7},$$

with the remaining factors and loop structure pushing the full classic estimate of effective SM CP violation in the symmetric plasma down to $\sim10^{-19}$-class. The honest caveats: the power counting is scheme- and observable-dependent (infrared-sensitive observables can soften it — one of Ch. 14's flagged enhancement routes), but *no* bulk, flavor-summed, perturbative observable escapes the quartet structure. Each loop also brings its measure, $1/16\pi^2 \approx 6\times10^{-3}$ per loop. These three numbers — $J_{CP} \sim 3\times10^{-5}$, $\mathcal S_{\text{GIM}} \sim 10^{-7}$, $(1/16\pi^2)^2 \sim 4\times10^{-5}$ — are the budget lines Ch. 14 inherits whenever wall physics is generated *perturbatively from the bulk*.

## 13.6 Why Standard-Model electroweak baryogenesis fails — twice

The classic scenario (electroweak baryogenesis, EWBG) deploys all the above at the electroweak transition: bubble walls of the broken phase sweep through the plasma; CP-violating reflection off the walls segregates charge; sphalerons outside the bubbles convert it; the broken phase inside freezes it. Two independent failures, and the distinction is strategic for this thesis:

**Failure 1 — equilibrium (kinematic).** The scenario needs a *first-order* transition: bubbles, walls, and out-of-equilibrium transport at the interface. Lattice computations settled the question: for the physical Higgs mass ($m_H = 125$ GeV $\gg$ the $\sim72$ GeV endpoint), the electroweak "transition" is a smooth **crossover** — no bubbles, no walls, no sustained departure from equilibrium. Sakharov's third condition simply is not met by Standard-Model cosmology.

**Failure 2 — CP magnitude (dynamic).** Even granting walls by hand, the transportable CP asymmetry is built from bulk CKM physics and pays the full §13.5 bill; estimates land near $\eta_B \sim 10^{-19}$, ten orders short.

The strategic point, stated exactly because Ch. 14 leans on it: **these failures are independent, and this framework structurally repairs the first while inheriting the second.** The arrow of expansion (Ch. 6) is a *permanent*, theorem-grade departure from equilibrium — no phase transition required, Sakharov-3 supplied at every epoch. The CP magnitude, by contrast, is repaired only if *boundary-localized* CP violation evades the bulk GIM structure — the wall breaks the chiral symmetry that organizes the quartet cancellation, so the question is genuinely open rather than settled-negative, but it is open, and the thesis budgets it as such (**[Open]**, Ch. 14.4, Ch. 27 item 4). A framework that claimed to fix both for free should not be believed; this one claims one theorem and one well-posed question.

## 13.7 The shopping list, filled

The interlude's content, organized as the checklist Ch. 14 walks through:

| Sakharov needs | Standard physics supplies | This framework supplies | Status |
|---|---|---|---|
| $B$ violation | sphaleron conversion of charge excess, $C_{\text{sph}} = 28/79$ | the excess itself: quantized spectral-flow pumping (Ch. 9, 12) | **[Standard]** + **[Theorem/Computed]** |
| C, CP violation | CKM phase, $J_{CP} \approx 3\times10^{-5}$; GIM-suppressed in bulk | wall-localized CP (mismatch $\Delta$); magnitude = $\varepsilon_{CP}$ budget | mechanism **[Theorem]**; magnitude **[Open]** |
| Non-equilibrium | — (crossover; SM has none) | the expansion arrow, permanent (Ch. 6) | **[Theorem]** |
| Bookkeeping to today | entropy dilution $45/2\pi^2 g_* = 0.0214$ | step/window counting $\to \mathcal P_\times$ (Ch. 14) | **[Standard]** + **[Theorem]** |

With the standard pieces audited and the division of labor explicit, the assembly can proceed.

---

**Validation.** `ch14_eta_assembly.py` prints the audited constants used here ($28/79 = 0.3544$; $45/2\pi^2g_* = 0.02136$ at $g_* = 106.75$). Figures: `ch13_fig1_sakharov_map.png`, `ch13_fig2_ewbg_failure_modes.png` (drawn schematics).
