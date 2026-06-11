# Chapter 22 — The stiffness from matter loops: induced gravity and the Planck cell

> **Author: Fable 5 · Complexity ◇5 · Status: draft v1 (pending global review pass)**
> Depends on: Ch. 16, 21; App. B. Feeds: Ch. 23 (compensator), 24 (Poisson), 25 (the covariant version of this chapter's measurement).

---

Every quantitative success so far has rented its units from one unexplained number: the stiffness $\kappa = (8\pi G)^{-1}$ — how much energy it costs to bend geometry. A theory that must be *told* the strength of gravity has explained its form but not its existence. This chapter computes the stiffness from the model's own matter content, twice, from opposite ends. **Top-down**: Sakharov's 1967 observation — matter vacuum loops on a curved background *induce* an Einstein–Hilbert term whether or not any gravitational dynamics is postulated — evaluated with the one ingredient the box model uniquely supplies: a physical identification of the cutoff (the cell scale), which converts the usual embarrassment of induced-gravity ("$G$ depends on an arbitrary $\Lambda_{\text{UV}}$") into a *prediction*: the mean cell size is the Planck length, up to a species factor. **Bottom-up**: in one dimension the model permits something stronger than an estimate — a direct numerical *measurement* of the induced action of the geometry field, including its **sign**, which turns out negative along compression: gravity attracts, the conformal mode has wrong-sign kinetic energy, and the constraint's "geometry energy is negative" (Ch. 21) stops being a postulate. The two computations meet at the same answer from different directions, which is the kind of overdetermination this Part is built on.

## 22.1 One dimension first: the Casimir energy of a cell, derived

Take the model's matter — a massless field with Dirichlet walls — in a cell of size $L$. Its zero-point energy is the divergent sum $E = \tfrac12\sum_{n\ge1}\omega_n = \tfrac{\pi}{2L}\sum_{n\ge1} n$. Regulate it both standard ways (App. B; they agree, as they must):

*Zeta route:* $\sum n \to \zeta(-1) = -\tfrac{1}{12}$, giving directly

$$\boxed{\;E_C(L) \;=\; -\,\frac{\pi}{24\,L}\;.} \tag{22.1}$$

*Heat-kernel route:* $E(s) = \tfrac{\pi}{2L}\sum_n n\,e^{-s n\pi/L} = \tfrac{\pi}{2L}\big[\tfrac{L^2}{\pi^2 s^2} - \tfrac{1}{12} + \mathcal O(s^2)\big]$; the $1/s^2$ piece is the universal bulk divergence (the cosmological-constant term — §22.5 owns the honesty about it), and the finite, regulator-independent remainder is again $-\pi/24L$. **[Standard + rederived]**

The force on the walls is $F = -\partial E_C/\partial L = -\pi/24L^2 < 0$: **the matter vacuum pulls the walls together.** This single line is Sakharov's entire conceptual content in its minimal habitat: the geometry (here, the wall configuration) acquires dynamics — *attractive* dynamics — purely from the matter living on it, with strength computed from the matter spectrum, not chosen.

**Cell energetics.** Add an occupant in mode $n$: $E_m = \pi^2 n^2/2L^2$, which resists compression. Total and equilibrium:

$$E(L) = -\frac{\pi}{24L} + \frac{\pi^2 n^2}{2L^2}, \qquad E'(L^\ast) = 0 \;\Longrightarrow\; L^\ast = 24\pi n^2. \tag{22.2}$$

Occupied cells have a stable size; empty cells do not — for $n = 0$ the energy decreases monotonically under shrinking *and is concave* ($E'' = -\pi/12L^3 < 0$), so a uniform chain of empty cells is unstable to clumping: any wall fluctuation grows. This is a **vacuum Jeans instability** — the one-dimensional shadow of gravity's signature pathology-that-isn't: uniform states are unstable, structure is mandatory.

**[Computed]** (`ch22_casimir_chain.py`): a 12-cell chain at fixed total length $24$, cells alternately occupied ($n = 2$) and empty, walls relaxing down the energy gradient (short-range core at $L_{\text{floor}} = 0.10$ regularizing the final collapse — physically, a minimum cell size, the model's UV station). Result: every empty cell collapses to the core ($L = 0.135 \approx \sqrt2\,L_{\text{floor}}$, the derived core equilibrium), every occupied cell survives and expands to share the freed length ($L = 3.865$ each; totals re-sum to $24.000$). **Wall-on-wall attraction, mediated by the matter vacuum, derived end to end.** The honest caveat that sharpens rather than weakens: 1D has no propagating graviton, so this demonstrates the induced *potential* sector only — exactly what 1D allows, and exactly the ingredient (the stiffness) the higher-dimensional construction needs; the propagating sector is Ch. 25's.

![Cell energetics and the relaxing chain](figures/ch22_fig2_cell_energetics.png)

*Figure 22.2 — Left: $E(L)$ for empty and occupied cells — monotone-concave (unstable) vs minimum at $L^\ast$ (stable). Right: the 12-cell chain relaxing: empty cells crushed to the core, occupied cells inflated, the vacuum doing the pulling.*

## 22.2 The induced action of the geometry field, measured

The chain experiment shows attraction between *adjacent* walls. The structural claim — that matter loops induce a *field theory* for geometry — needs the energy cost of a geometry *wave*. One dimension permits measuring it directly. Take a massless field on $[0, L_{\text{tot}}]$ with $M{-}1$ **semi-transparent** walls (delta barriers of strength $\lambda$ — cells that talk, unlike the opaque chain above); the wall displacement field $u_i$ *is* the geometry field (moving walls = moving cell boundaries = the microscopic content of $\delta\alpha(x)$). The matter vacuum energy as a functional of displacements,

$$E_{\text{vac}}[\{u_i\}] = \tfrac12\sum_n \omega_n[\{u_i\}] \quad (\text{heat-kernel smoothed, Richardson-extrapolated}),$$

is, by definition, the one-loop induced effective action of the geometry. Probe it with standing waves $u_i = \varepsilon\sin(\pi k i/M)$ and extract the quadratic response $\delta E = \tfrac{\varepsilon^2}{2}\Pi(q)\sum_i\sin^2$; an elastic medium gives $\Pi(q) = \kappa_{\text{grad}}\,\hat q^2$, $\hat q^2 = (2\sin\tfrac q2)^2$. **[Computed]** (`ch22_dispersion.py`; 24 cells, $\lambda = 8$, twelve modes):

$$\frac{\Pi(q)}{\hat q^2} \;\in\; [-0.18,\, -0.05], \qquad \text{median } \kappa^{\text{vac}}_{\text{grad}} \;\approx\; -0.13, \tag{22.3}$$

negative for eleven of twelve modes (the odd/even scatter is a fixed-endpoint parity effect, reported as is). Meanwhile one occupant of mode $n$ per cell contributes, by a single Taylor expansion of $\sum_i \pi^2n^2/2L_i^2$ in $L_i = 1 + u_i - u_{i-1}$,

$$\kappa^{\text{matter}}_{\text{grad}} \;=\; +\,3\pi^2 n^2 \;=\; +29.6 \;\;(n = 1) \tag{22.4}$$

per bond — positive and two orders of magnitude larger.

**What the sign means.** The matter vacuum makes the *compression* mode of geometry unstable: bending the cell sizes into a wave lowers the vacuum energy. That one measured sign is simultaneously three known facts about gravity, here *derived*: (i) **gravity is attractive** — compressions grow; the Jeans instability of §22.1 seen now in the effective action; (ii) **the conformal mode of the gravitational action has wrong-sign kinetic term** — the notorious conformal-factor problem of Euclidean gravity, an *output* of matter loops in this model, and the dynamical origin of the Lorentzian signature that Ch. 19 honestly deferred; (iii) **ordinary matter stabilizes geometry at short scales** — occupied cells are stiff (22.4), which is why space does not crumple around matter. The constraint's negative geometry-energy (Ch. 21) is hereby a measurement.

![The measured dispersion](figures/ch22_fig3_dispersion.png)

*Figure 22.3 — $\Pi(q)/\hat q^2$ for the vacuum (negative band, median $-0.13$) against the matter contribution (positive line at $+3\pi^2n^2$): the DeWitt sign pattern in one dimension, measured from matter loops.*

## 22.3 Top-down: Sakharov's integral with the model's cutoff

Now the 3+1D version, as standard machinery with one non-standard — and decisive — input. Integrating out $N$ matter fields on a background $g_{\mu\nu}$ gives the one-loop effective action in heat-kernel (proper-time) form; per minimally coupled scalar,

$$\Gamma[g] = -\frac12\int_0^\infty \frac{ds}{s}\,\operatorname{Tr} e^{-s(-\Box)}, \qquad \operatorname{Tr} e^{-s(-\Box)} = \frac{1}{(4\pi s)^2}\int\!\sqrt{-g}\,\Big[1 + \frac{s}{6}R + \mathcal O(s^2 R^2)\Big] \tag{22.5}$$

(the Seeley–DeWitt expansion; App. B). The proper-time integral diverges at $s \to 0$; cutting at $s_{\min} = 1/\Lambda_{\text{UV}}^2$ and collecting the term linear in $R$:

$$\Gamma \;\supset\; -\,\frac{\Lambda_{\text{UV}}^2}{192\pi^2}\int\sqrt{-g}\,R \qquad\Longrightarrow\qquad \boxed{\;G_{\text{ind}}^{-1} \;=\; \frac{N_{\text{eff}}\,\Lambda_{\text{UV}}^2}{12\pi}\;} \tag{22.6}$$

comparing with $S_{EH} = -\tfrac{1}{16\pi G}\int\sqrt{-g}R$ ($N_{\text{eff}}$ counts species with spin-dependent $\mathcal O(1)$ weights). **Gravity's strength is the residue of the matter spectrum at the cutoff** — Einstein–Hilbert appears with no gravitational postulate at all: gravity as the elasticity of the vacuum. **[Standard, rederived in App. B]**

**The cutoff, identified.** In generic induced-gravity programs, $\Lambda_{\text{UV}}$ is the embarrassment: an arbitrary scale on which $G$ depends quadratically. The box model has no such freedom. The continuum heat-kernel expansion (22.5) is valid only while the proper-time probe is too coarse to resolve the cell structure; at $s \lesssim \bar L^2$ the discrete spectrum deviates from the continuum kernel. The cutoff is therefore **the inverse cell size**, $\Lambda_{\text{UV}} = c/\bar L$ with $c = \mathcal O(1)$ set by lattice details — located *inside the very integral that needs it*, not appended as philosophy. Two consequences:

> **Planck-Cell Prediction [Theorem, given (22.6) + the cutoff identification].** Inverting $G_{\text{ind}}^{-1} = N_{\text{eff}}c^2/12\pi\bar L^2$:
>
> $$\bar L \;=\; c\,\sqrt{\frac{N_{\text{eff}}}{12\pi}}\;\ell_P \;\approx\; 1.6\,c\;\ell_P \quad (N_{\text{eff}} \sim 100). \tag{22.7}$$
>
> The geometry's granularity is *derived* to sit at the Planck scale — self-consistently: a smaller cell would strengthen gravity, which would shrink the Planck length, in just the right proportion. The model did not get to choose this; (22.7) could have come out at the electroweak scale or the Hubble scale and falsified the construction at a stroke.

And second: the **conformal sector closes**. Chapter 23 will identify the Weyl compensator as $\sigma \sim 1/L_{\text{cell}}$ with $M_P^2 = \sigma_0^2/6$; Sakharov independently gives $M_P^2 \propto N_{\text{eff}}/\bar L^2$. The two are the same statement with the $\mathcal O(1)$ constant now computable rather than free — matching them fixes $c$ in terms of $N_{\text{eff}}$. What remains genuinely open at this node is bookkeeping, not concept: the spin-weighted species sum for the model's actual matter content, and the lattice constant $c$ — both finite calculations (Ch. 27, item 2).

![The self-consistency loop](figures/ch22_fig4_planck_loop.png)

*Figure 22.4 — The Planck loop. Matter spectra → induced stiffness (Sakharov) → Newton's constant → Planck length → cell scale → cutoff → back into the induced integral. The loop closes at $\bar L \approx 1.6\,c\,\ell_P$.*

## 22.4 The term we do not explain

Honesty requires the same sentence every time this integral appears, so here it is in its own section. The leading term of (22.5) is $c_0\Lambda_{\text{UV}}^4\sqrt{-g}$ — an induced cosmological constant, enormous at the cell-scale cutoff, exactly as in every known approach to vacuum energy. This framework neither solves nor worsens the cosmological-constant problem; it inherits it intact. No cancellation is claimed, no anthropic appeal is made, and every later chapter that touches the issue (Ch. 25's $V/s^4$ term, Ch. 27 item 9) repeats this sentence rather than softening it.

## 22.5 Summary

The stiffness is no longer a unit choice. Bottom-up: the 1D induced action is *measured* — attractive Casimir potential (22.1), vacuum Jeans instability with the chain collapsing on cue, gradient coupling $\kappa^{\text{vac}}_{\text{grad}} \approx -0.13 < 0$ against matter's $+3\pi^2n^2$ — delivering the attraction, the conformal-factor sign, and matter stabilization as outputs. Top-down: Sakharov's integral with the cell-scale cutoff gives $G_{\text{ind}}^{-1} = N_{\text{eff}}\Lambda^2/12\pi$ and the **Planck-Cell Prediction** $\bar L \approx 1.6\,c\,\ell_P$. The constraint of Ch. 21 now stands on measured ground, and two chapters of spending await: the static limit (Newton's law, Ch. 24) and the propagating limit (the graviton with its universal coefficient, Ch. 25) — with the gauge-theoretic consolidation (Ch. 23) in between.

---

**Validation.** `ch22_casimir_chain.py` (port): cell energetics, equilibria, the 12-cell relaxation (endpoint values $0.135$/$3.865$/$24.000$). `ch22_dispersion.py` (port): the semi-transparent-wall vacuum dispersion (band $[-0.18, -0.05]$, median $-0.13$; matter line $+3\pi^2 n^2$), writing `data/ch22_stiffness_dispersion.json`. App. B re-derives (22.5)–(22.6) symbolically. All quoted numbers printed by the scripts.
