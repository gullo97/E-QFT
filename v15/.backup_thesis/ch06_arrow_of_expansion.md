# Chapter 6 — Many bodies and the arrow of expansion

> **Author: Fable 5 · Complexity ◇4–◇5 · Status: draft v1 (pending global review pass)**
> Depends on: Ch. 2–5. Feeds: Ch. 7 (E-breaking), Ch. 14 ($f_{\text{neq}}$), Ch. 20–21 (ladder dynamics with matter).

---

A single particle on the geometry ladder moves symmetrically: nothing in Chapters 3–4 distinguishes expansion from contraction. This chapter adds the one ingredient the real universe insists on — *other particles* — and proves that the symmetry shatters. When an active particle takes an iso-energy step, every other particle in the box (every **spectator**) must come along for the ride, and the price of the ride is direction-dependent: free for expansion, costly for contraction, and *multiplicatively* costly across spectators. The ratio of contraction to expansion probability falls exponentially with the number of particles. For a box holding more than a couple of dozen particles, contraction effectively never happens.

This is the framework's answer to one of cosmology's silent questions — *why does the universe expand rather than contract?* — and it is an answer of an unusual kind: not an initial condition, not a potential engineered to inflate, but a counting argument about overlap integrals, of the same logical species as the second law of thermodynamics.

The chapter has a second purpose, equally important for the thesis's defensibility. Two ingredients of this argument were, historically, *postulates*: the rule that an unrenormalized projection norm is a transition probability, and the hypothesis that spectator contributions factorize. Section 6.5 proves both — they are first-order perturbation theory in the sector arena of Ch. 3. The arrow of expansion thereby moves from "well-motivated model" to "theorem within the framework", and the assumptions that remain are exactly the two flagged in Ch. 3: the sector arena itself, and the form of the inter-sector coupling.

## 6.1 The cast: one active particle, $N$ spectators

Let $N+1$ non-interacting particles share the box $[0, L]$, in modes $n_1; n_{2,1}, \ldots, n_{2,N}$. The many-body state is the product

$$\Psi(x_0, x_1, \ldots, x_N) \;=\; \varphi_{n_1}(x_0; L)\prod_{\alpha=1}^N \varphi_{n_{2,\alpha}}(x_\alpha; L),$$

(statistics — symmetrization or antisymmetrization — does not affect the norms computed below as long as the occupied modes are distinct, and we comment where it would matter). The **active particle** is the one whose iso-energy step drives the geometry change:

$$\text{expansion: } n_1 \to n_1 + 1,\;\; L \to L' = \frac{n_1+1}{n_1}L, \qquad \text{contraction: } n_1 \to n_1 - 1,\;\; L \to L'' = \frac{n_1-1}{n_1}L .$$

The **spectators** have no say in the step; the box changes under them, suddenly (the regime of Ch. 5, $\eta_{\text{ad}} \ll 1$). Their state the instant after is their state the instant before — embedded ($\iota$) for expansion, restricted ($R$) for contraction, by the uniqueness argument of §3.3.

The question with dynamical content: *with what amplitude does the whole assembly arrive in the new sector?*

![Spectator projection](figures/ch06_fig1_spectator_projection.png)

*Figure 6.1 — A spectator rides a geometry change. The spectator's mode (blue) re-expanded over the new box's eigenbasis: dominant diagonal peak at $k = n_2$, sinc-tail leakage elsewhere. Under expansion (left) the weights sum to exactly 1; under contraction (right) a deficit remains — the weight of the wavefunction in the excised region.*

## 6.2 Expansion is free: the Expansion Norm and Kinetic-Energy Preservation Theorems

> **Expansion Norm Theorem [Theorem].** Under sudden expansion $L \to L' > L$, each spectator's total transfer weight is exactly unity:
>
> $$\sum_{k=1}^{\infty} \big|c_k(n_2, s)\big|^2 \;=\; \big\|\iota\,\varphi_{n_2}\big\|^2 \;=\; 1. \tag{6.1}$$
>
> *Proof.* Lemma 3.1(a): extension by zero is an isometry, and the $c_k$ are the components of $\iota\varphi_{n_2}$ in the complete new eigenbasis; Parseval does the rest. $\blacksquare$

The spectator is guaranteed to land *somewhere* in the new box — scattered over modes by the sinc tails of (3.8), but with no probability lost. And not only the probability survives; so does the energy:

> **Kinetic-Energy Preservation Theorem [Theorem].** Under sudden expansion the spectator's expected kinetic energy is unchanged:
>
> $$\sum_k |c_k|^2\, E_k(L') \;=\; E_{n_2}(L). \tag{6.2}$$
>
> *Proof.* The expectation of the kinetic operator in the embedded state is $\frac{1}{2m}\int_0^{L'} |\partial_x(\iota\varphi_{n_2})|^2 dx = \frac{1}{2m}\int_0^{L}|\varphi_{n_2}'|^2 dx = E_{n_2}(L)$, since the extension vanishes identically on $[L, L']$. (The embedded function has a kink at $x = L$ — $\varphi_{n_2}'(L^-) \ne 0$, derivative zero beyond — but a first-derivative jump is harmless: $|\partial_x\psi|^2$ remains integrable, and the left side of (6.2), computed spectrally, converges to it.) Expanding the same quadratic form in the new eigenbasis gives the sum rule. $\blacksquare$

So expansion asks nothing of the spectators: unit probability, energy intact. The "iso-energy" label of the transition is earned, for expansion, by the whole assembly and not just the active particle. **[Computed]** `ch06_norms.py` verifies (6.1) and (6.2) to the $10^{-6}$ truncation accuracy of 4000-mode sums across a grid of $(n_2, s)$.

## 6.3 Contraction costs: the Contraction Norm Theorem

Now run the step the other way: $L \to L'' = \frac{n_1 - 1}{n_1} L < L$. The spectator's state is *restricted* to the surviving interval, and Lemma 3.1(b) says restriction is a contraction — some norm is gone. How much, exactly?

> **Contraction Norm Theorem [Theorem].** For a spectator in mode $n_2$ under the iso-energy contraction driven by an active particle in mode $n_1$,
>
> $$\mathcal N_{\text{con}}(n_1, n_2) \;\equiv\; \big\|R\,\varphi_{n_2}\big\|^2 \;=\; \frac{n_1 - 1}{n_1} \;+\; \frac{\sin\!\big(2\pi n_2/n_1\big)}{2\pi n_2}\;<\;1. \tag{6.3}$$

*Derivation, in full.* The lost weight is the old wavefunction's probability in the excised region; the surviving weight is

$$\mathcal N_{\text{con}} = \int_0^{L''}\frac{2}{L}\sin^2\!\Big(\frac{n_2\pi x}{L}\Big)\,dx = \frac{2}{L}\left[\frac{x}{2} - \frac{L}{4n_2\pi}\sin\!\Big(\frac{2n_2\pi x}{L}\Big)\right]_0^{L''} = \frac{L''}{L} - \frac{\sin\!\big(2n_2\pi L''/L\big)}{2n_2\pi}.$$

Insert $L''/L = (n_1{-}1)/n_1$ and use $\sin\!\big(2\pi n_2 - \tfrac{2\pi n_2}{n_1}\big) = -\sin\tfrac{2\pi n_2}{n_1}$ (the $2\pi n_2$ winds away):

$$\mathcal N_{\text{con}} = \frac{n_1-1}{n_1} + \frac{\sin(2\pi n_2/n_1)}{2\pi n_2}. \qquad \blacksquare$$

Read the two terms. The first is purely geometric: the fraction of the box that survives. The second is the spectator's *self-interference* with the cut: it oscillates in sign with $n_2$, is bounded by $\frac{1}{2\pi n_2}$, and vanishes for high modes — a fast-oscillating spectator is uniformly spread, and loses exactly the geometric fraction. Representative values **[Computed]** (`ch06_norms.py`): $\mathcal N_{\text{con}}(5,3) = 0.76882$, $\mathcal N_{\text{con}}(8,5) = 0.85249$, $\mathcal N_{\text{con}}(12,7) = 0.90530$. The deficit $1 - \mathcal N_{\text{con}}$ is the probability that the spectator *would have been caught in the excised region* — an event with no final state. What that means dynamically is settled by the theorem of §6.5; first, one more structural fact about contraction.

![Norm deficit](figures/ch06_fig2_norm_deficit.png)

*Figure 6.2 — The contraction deficit $1 - \mathcal N_{\text{con}}(n_1, n_2)$ over the $(n_1, n_2)$ plane. The deficit approaches $1/n_1$ (the excised fraction) for fast spectators and oscillates around it for slow ones. There is no $(n_1, n_2)$ with zero deficit: contraction is never free.*

## 6.4 Contraction is violent: the Sudden-Contraction Divergence Theorem

The deficit (6.3) is finite, which might suggest contraction is merely *somewhat* suppressed. The energy bookkeeping says otherwise.

> **Sudden-Contraction Divergence Theorem [Theorem].** The restricted spectator state, re-expanded in the new box's eigenbasis with coefficients $\tilde c_k$, has divergent expected kinetic energy:
>
> $$\sum_k |\tilde c_k|^2 E_k(L'') \;=\; \infty, \tag{6.4}$$
>
> whenever $\varphi_{n_2}(L'') \neq 0$, i.e. for all but a measure-zero set of resonant spectator modes.

*Derivation.* The restricted function generally does *not* vanish at the new wall: $\varphi_{n_2}(L'') \ne 0$. A function with a boundary-value mismatch against the Dirichlet basis has sine-series coefficients with the slow universal tail obtained by two integrations by parts,

$$\tilde c_k = \int_0^{L''}\!\!\sqrt{\tfrac{2}{L''}}\sin\!\Big(\tfrac{k\pi x}{L''}\Big)\varphi_{n_2}(x)\,dx = \underbrace{(-1)^{k+1}\sqrt{\tfrac{2}{L''}}\,\frac{L''}{k\pi}\,\varphi_{n_2}(L'')}_{\text{boundary term}} + \mathcal O\!\Big(\tfrac{1}{k^2}\Big),$$

so $|\tilde c_k|^2 \sim \frac{2L''}{k^2\pi^2}\varphi_{n_2}(L'')^2$ while $E_k(L'') = \frac{k^2\pi^2}{2m L''^2}$ grows as $k^2$: the summand tends to the constant $\frac{\varphi_{n_2}(L'')^2}{m L''}$, and the sum diverges linearly in the mode cutoff. $\blacksquare$

Physically: slamming the wall onto a region where the spectator has finite amplitude creates a wavefunction discontinuity, and a discontinuity costs unbounded kinetic energy in the strict sudden limit. Two consequences matter downstream. First, *strict* sudden contraction is energetically self-inconsistent — honest contraction must be at least partially adiabatic, which further suppresses it (Ch. 5's crossover) and makes "iso-energy", as applied to the contraction branch, a label about the active particle only. We flag this asymmetry now and account for it in the interpretation of §6.7. Second, the *resonant exceptions* — spectator modes with $\varphi_{n_2}(L'') = 0$, i.e. $n_2$ a multiple of $n_1$ — are exactly the spectators sharing the active particle's iso-energy structure; they ride the contraction cheaply. Generic spectators do not.

## 6.5 From norms to rates: the Golden-Rule Embedding Theorem

Everything so far computes *norms*. To claim an arrow we need *probabilities of dynamical transitions*, and historically this step was bridged by two assumptions: (i) the unrenormalized projection norm is the transition probability weight; (ii) spectator weights multiply. Both are now theorems, and the proof is short enough to give completely.

> **Toolbox: first-order transitions in one page.** Let $H = H_0 + V$ with $V$ weak, and let $|i\rangle$ be an eigenstate of $H_0$ with energy $E_i$. To first order in $V$, the amplitude to reach eigenstate $|f\rangle \ne |i\rangle$ after time $t$ is $A_f(t) = -i\,\langle f|V|i\rangle \int_0^t e^{i(E_f - E_i)t'}dt' $. For times short against the inverse energy spread, the integral is $\approx t$ and $P_f(t) = t^2|\langle f|V|i\rangle|^2$; summed over a complete set of final states this is $t^2\,\|V|i\rangle\|^2$. (At longer times the familiar golden-rule linear-in-$t$ regime takes over; the short-time quadratic regime is the cleanest for the structural statement below, and the numerics confirm both regimes.)

> **Golden-Rule Embedding Theorem [Theorem].** In the sector arena $\mathcal H = \bigoplus_L \mathcal H(L)$, couple neighbouring sectors by the embedding,
>
> $$V \;=\; g\,\big(\iota + \iota^\dagger\big) \tag{6.5}$$
>
> ($\iota$ = extension-by-zero toward the larger sector, $\iota^\dagger = R$ toward the smaller; $g$ a coupling constant). Then for an initial product state $\Psi$ of $N{+}1$ particles in sector $L$, the total short-time probability of transfer to the neighbouring sector is
>
> $$P_{\text{transfer}}(t) \;=\; (gt)^2\,\big\|\,\iota_{(N+1)}\,\Psi\,\big\|^2 + \mathcal O(t^4), \qquad \iota_{(N+1)} = \iota^{\otimes(N+1)}, \tag{6.6}$$
>
> and the norm **factorizes over particles**:
>
> $$\big\|\iota^{\otimes(N+1)}\Psi\big\|^2 \;=\; \prod_{\beta=0}^{N}\big\|\iota\,\psi_\beta\big\|^2 . \tag{6.7}$$
>
> In particular: each spectator contributes factor $1$ (expansion) or $\mathcal N_{\text{con}}(n_1, n_{2,\alpha})$ (contraction), exactly the norms of §6.2–6.3.

*Proof.* Equation (6.6): apply the Toolbox with $V = g\iota$; summing $|\langle f|\iota|\Psi\rangle|^2$ over the complete eigenbasis $\{|f\rangle\}$ of the new sector is, by Parseval, the squared norm of $\iota\Psi$ in that sector. No renormalization of the embedded state is allowed to sneak in — the coupling (6.5) acts on the *unrenormalized* embedding, which is what first-order theory evaluates. Equation (6.7): on a product state, $\iota$ acts as the tensor product of single-particle embeddings (the box changes under every particle identically and independently — non-interaction is used exactly here), and norms of tensor products factorize. $\blacksquare$

*Unitarity bookkeeping.* For contraction the weights (6.7) multiply to less than one. Nothing is lost: the missing probability is, to this order, the amplitude that the system *remains in sector $L$* — the contraction attempt simply fails. "A spectator would have been caught in the excised region" is not an absorbed event; it is a non-event. The dynamics is unitary on $\bigoplus_L \mathcal H(L)$ throughout.

**[Computed]** (`ch06_golden_rule.py`, two-sector model, 220 modes per sector): Parseval sums reproduce $\mathcal N_{\text{con}}$ — $0.76749$ vs analytic $0.76882$ for $(n_1, n_2) = (5,3)$; $0.85112$ vs $0.85249$ for $(8,5)$; $0.90373$ vs $0.90530$ for $(12,7)$ — with residuals scaling as $1/n_{\max}$ exactly as the sinc tails dictate; expansion gives $1.00000$. Exact unitary evolution of the coupled two-sector Hamiltonian, no perturbation theory, gives $P(t)/(gt)^2 \to$ the Parseval value from below as $t \to 0$ (e.g. $0.965 \to 0.992$ toward $1$ for expansion at $n_2 = 3$), confirming both the theorem and its claimed regime of validity.

![Golden-rule convergence](figures/ch06_fig3_golden_rule_convergence.png)

*Figure 6.3 — The postulates become theorems. Left: Parseval sums vs the closed-form $\mathcal N_{\text{con}}$ for three $(n_1, n_2)$ pairs as the basis grows. Right: exact short-time evolution $P(t)/(gt)^2$ converging to the embedding norm as $t \to 0$.*

## 6.6 The arrow: the Expansion-Preference Theorem

Assemble the factors. The active particle's own weight is a Geometric Overlap of Ch. 3: its target mode has the *same physical wavenumber* in the new box, so its factor is the promoted/demoted-mode norm,

$$\text{expansion: } \big|\mathcal O_{n_1 \to n_1+1}\big|^2 = \frac{n_1}{n_1+1}, \qquad \text{contraction: } \big|\mathcal O_{n_1 \to n_1 - 1}\big|^2 = \frac{n_1 - 1}{n_1},$$

(the contraction value by the same same-wavenumber computation as (3.10), now with the integral cut at $L''$ where $\sin(2\pi(n_1{-}1)) = 0$ kills the boundary term exactly — the active particle is its own resonant exception, which is *why* the step is kinematically allowed at all). Each spectator contributes $1$ or $\mathcal N_{\text{con}}$. Dividing:

> **Expansion-Preference Theorem [Theorem].** For an active particle in mode $n_1$ and $N$ spectators in modes $\{n_{2,\alpha}\}$,
>
> $$\frac{P_{\text{con}}}{P_{\text{exp}}} \;=\; \frac{(n_1-1)/n_1}{\,n_1/(n_1+1)\,}\;\prod_{\alpha=1}^{N}\mathcal N_{\text{con}}(n_1, n_{2,\alpha}) \;=\; \Big(1 - \frac{1}{n_1^2}\Big)\prod_{\alpha=1}^{N}\mathcal N_{\text{con}}(n_1, n_{2,\alpha}). \tag{6.8}$$
>
> Since every factor $\mathcal N_{\text{con}} < 1$ (Theorem 6.3), the ratio decays exponentially in $N$:
>
> $$\frac{P_{\text{con}}}{P_{\text{exp}}} \;\sim\; e^{-\sigma N}, \qquad \sigma = -\big\langle \ln \mathcal N_{\text{con}}\big\rangle_{\text{spectators}} > 0 . \tag{6.9}$$

The rate is set by the active particle's step coarseness: since the deficit is dominated by the excised fraction $1/n_1$, $\sigma \approx -\ln(1 - 1/n_1) \approx 1/n_1$, giving $N_{\text{crit}}(10^{6}) \approx 6\ln 10/\sigma \approx 14\,n_1$. **[Computed]** (`ch06_norms.py`, spectators sampled over $n_2 \in [1, 4n_1]$): $\sigma = 0.216$ at $n_1 = 5$ ($N_{\text{crit}} \approx 64$ for $10^{6}$ suppression, $\approx 107$ for $10^{10}$); $\sigma = 0.082$ at $n_1 = 12$ ($N_{\text{crit}} \approx 168$ and $280$). Near the natural boundary ($n_1 = 2$–$3$, the deep quantum regime) $N_{\text{crit}}$ drops to a few tens. The constants matter less than the structure: the suppression is exponential with macroscopically enormous exponent — for any laboratory- or cosmology-scale content ($N \gtrsim 10^{20}$), the ratio is zero by any physical standard, at every $n_1$. *Expansion is not preferred by initial conditions; it is preferred by the measure of quantum mechanics on confined matter.*

![Asymmetry vs spectator number](figures/ch06_fig4_asymmetry_vs_N.png)

*Figure 6.4 — The arrow. $-\ln(P_{\text{con}}/P_{\text{exp}})$ grows linearly with spectator number (slope $\sigma$, set by $\approx 1/n_1$); the $10^{6}$-suppression threshold marked. Two active-mode cases shown ($n_1 = 5, 12$); ensemble-averaged over spectator modes.*

## 6.7 What precisely has been shown — and the three honest caveats

It is worth being exact, because this result carries weight downstream (Ch. 7 breaks a symmetry with it; Ch. 14 uses it as the departure-from-equilibrium factor).

**Established at theorem strength, given [Postulate: sector arena] + coupling (6.5):** the probability rule, the factorization, the per-spectator deficits, the exponential ratio (6.8)–(6.9), and energy bookkeeping (6.2), (6.4). The historic gap between "norm inequality" and "rate suppression" is closed.

**Caveat 1 — the coupling form.** The theorem used $V = g\,\iota$. Chapter 20 will find that the *dilation generator* coupling is the one that reproduces general relativity's orbits in the gravitational sector. The two couplings share the structural features the arrow needs (isometric toward larger boxes, contractive toward smaller — the generator's matrix elements inherit exactly the same embedding asymmetry), so the arrow is robust across the candidates; but a single coupling unifying Ch. 6 and Ch. 20 is the cleanest remaining structural question of Part I, and we flag it in Ch. 27.

**Caveat 2 — the strict sudden idealization.** Theorem 6.4 says strict sudden contraction is not merely suppressed but energetically inconsistent; real contractions are partially adiabatic, and adiabatic corrections *reduce* the deficit per spectator while *slowing* the step. The direction of the conclusion is unaffected (every regime suppresses contraction relative to expansion); the precise $\sigma$ is regime-dependent.

**Caveat 3 — small-$n_1$ examples vs the thermal regime.** Worked numbers at $n_1 = 5$ or $s = 2$ make vivid figures, but a thermal cosmological setting has $n_1 \gg 1$ and nearly infinitesimal steps. All theorems hold uniformly in $n_1$; only quantitative illustrations should be read at small $n_1$ with care. Chapter 14 uses the thermal regime.

**On statistics.** For fermions (Part II's case) the product-norm factorization (6.7) holds with Slater determinants of distinct modes, because the embedding is a single-particle map and determinants of isometries compose; Pauli blocking modifies which final states are *available*, an effect that further disfavors contraction (fewer escape channels for deficit amplitude). We use only the direction of this statement, never its magnitude.

## 6.8 Summary

Spectators turn the geometry ladder into a one-way street. The four norm theorems (expansion free and energy-preserving; contraction deficient by closed form, divergent in energy when strict), the **Golden-Rule Embedding Theorem** that converts norms into rates and proves factorization, and the **Expansion-Preference Theorem** with its exponential (6.9): together these constitute the framework's derivation of the arrow of expansion. The next chapter asks the question this arrow begs: if the universe's quantum mechanics breaks the symmetry between *growing* and *shrinking*, what discrete symmetry — if any — survives, and what must break along with it?

---

**Validation.** `ch06_norms.py`: closed forms (6.1)–(6.3), the deficit map (Fig. 6.2), sum rule (6.2), divergence rate of (6.4), the asymmetry plot (Fig. 6.4) with $\sigma$ and $N_{\text{crit}}$. `ch06_golden_rule.py`: two-sector Parseval vs analytic table, exact-evolution convergence (Fig. 6.3). Figure 6.1 (spectator projection coefficients) is produced by `ch06_norms.py --projection`.
