# Chapter 11 — The polarization diagnostic: what truncated pair-creation sums measure

> **Author: Fable 5 · Complexity ◇5 · Status: draft v1 (pending global review pass)**
> Depends on: Ch. 8–10, App. A. Feeds: Ch. 12 (the control quench), Ch. 14 (what the kinematic amplitude must *not* be).

---

Chapter 10 proved that a bulk CP phase produces exactly zero net charge under sudden expansion. And yet: compute the Bogoliubov coefficients for that very quench, sum the standard charge expression over a few dozen modes, and a clean, reproducible, nonzero number appears — odd in the phase as $\sin 2\delta$, stable under reasonable variations, passing every symmetry check one throws at it, including the ECPT-oddness identity of Ch. 7 to machine precision. At the parameters where such an analysis would naturally be run for the top quark at the electroweak crossover ($s = 2$, $mL = 1.085$, $\delta = 0.015$, $18$ modes per branch), the number is

$$A \;\equiv\; \frac{\Delta Q^{\text{trunc}}}{\sin 2\delta} \;\approx\; 0.076\,,$$

**[Computed]** ($0.0767$ from our independent solver; `ch11_polarization.py`). The Zero-Charge Theorem says this cannot be net charge. It is not a numerical error — the number converges. It is not noise — it is exactly odd in $\delta$. So the only scientifically acceptable response is to find out *what it is*. This chapter computes the infinite-basis limit of the truncated charge sums in closed form and identifies the answer: **a redistribution of vacuum polarization between the bag walls** — real physics, wrong observable. The diagnosis closes Part II's negative arc and hands Ch. 12 its two sharpest tools: the spectral-flow rule of practice, and the control-quench test.

## 11.1 The two truncated forms

Set up the quench precisely (machinery from Ch. 8/App. A). Old box $[0, L]$, new box $[0, L']$, old eigenmodes $\psi_n$ (energies $E_n$), new eigenmodes $\psi'_k$ (energies $E'_k$), overlaps

$$\alpha_{kn} = \int_0^{L} \psi_k'^{\dagger}(x)\,\psi_n(x)\,dx$$

(integration stops at $L$: old modes are extended by zero). Prepare the old vacuum; count, against the new vacuum. Two natural truncated expressions for "net charge produced", each summing $n_{\max}$ modes per branch:

$$\Delta Q_{\text{cross}} = \!\!\sum_{\substack{k:\,E'_k>0\\ n:\,E_n<0}}\!\! |\alpha_{kn}|^2 \;-\!\! \sum_{\substack{k:\,E'_k<0\\ n:\,E_n>0}}\!\! |\alpha_{kn}|^2, \qquad \Delta Q_{\text{same}} = \!\!\sum_{\substack{k:\,E'_k>0\\ n:\,E_n>0}}\!\! \big(\text{1} - \text{occupancy bookkeeping}\big)\ldots \tag{11.1}$$

— the *cross form* counts particles materializing from the old sea minus antiparticles from the old positive branch; the *same form* reaches the same formal object by tracking the conserved charge through the same-branch overlaps instead (the two are algebraically equal at infinite basis, by completeness). At $n_{\max} = 18$ at the representative point, they disagree by $20\%$ — and that discrepancy, far from being a blemish to average over, is **the smoking gun**: a unitary transformation truncated at finite rank is not unitary, and two expressions equal by unitarity must split by exactly the defect. A computation whose two honest forms disagree by 20% is telling you its target has not converged — or is not what you think it is. Both, here.

## 11.2 The Polarization Identity

Send the truncation to infinity and ask what the cross form actually converges to.

> **Polarization Identity [Theorem].** Let the old spectrum be $E \to -E$ symmetric (the case at hand, by the Spectral Mirror Theorem). Then
>
> $$\Delta Q_{\text{cross}}(\infty) \;=\; \frac12\,\operatorname{Tr}_{[0,L]}\operatorname{sgn}\big(H_{L'}\big) \;=\; \frac12\sum_k \operatorname{sgn}(E'_k)\;w_k\,, \qquad w_k \equiv \int_0^{L}\big|\psi'_k(x)\big|^2 dx \tag{11.2}$$
>
> — one half the *new* Hamiltonian's sign function, traced over the *old* region only.

*Proof.* Write the two branch sums as traces. With $\rho_\pm$ the projectors onto the old positive/negative branches (complete on $[0, L]$: $\rho_+ + \rho_- = \mathbb 1_{[0,L]}$) and $P'_\pm$ the new-branch projectors,

$$\Delta Q_{\text{cross}}(\infty) = \sum_{n:\,E_n<0}\langle\psi_n|P'_+|\psi_n\rangle \;-\; \sum_{n:\,E_n>0}\langle\psi_n|P'_-|\psi_n\rangle = \operatorname{Tr}_{[0,L]}\big[\rho_-\,P'_+ \;-\; \rho_+\,P'_-\big].$$

Substitute $\rho_\pm = \tfrac12(\mathbb 1 \pm S)$, $P'_\pm = \tfrac12(\mathbb 1 \pm S')$ with $S = \operatorname{sgn}(H_L)$, $S' = \operatorname{sgn}(H_{L'})$ and expand:

$$\rho_-P'_+ - \rho_+P'_- = \tfrac14\Big[(\mathbb 1 - S)(\mathbb 1 + S') - (\mathbb 1 + S)(\mathbb 1 - S')\Big] = \tfrac12\big(S' - S\big).$$

Now trace over $[0, L]$. The $S$ term gives $\tfrac12\operatorname{Tr}_{[0,L]} S = \tfrac12\sum_n \operatorname{sgn}(E_n) = \tfrac12\,\eta(L) = 0$ by the assumed mirror symmetry (old modes live entirely on $[0,L]$, so the regional trace is the full one). The $S'$ term gives $\tfrac12\sum_k \operatorname{sgn}(E'_k)\langle\psi'_k|\mathbb 1_{[0,L]}|\psi'_k\rangle$, which is (11.2). $\blacksquare$

**Reading the identity.** Recall (Ch. 8–9) that the new vacuum's charge *density* is $\rho'_{\text{vac}}(x) = -\tfrac12\sum_k \operatorname{sgn}(E'_k)\,|\psi'_k(x)|^2$ (regulated). So (11.2) says, with signs unwound:

$$\Delta Q_{\text{cross}}(\infty) \;=\; -\int_0^{L}\rho'_{\text{vac}}(x)\,dx \;=\; -\,\big(\text{new vacuum's charge residing in the old region}\big). \tag{11.3}$$

Integrated over the *whole* new box this would be $-Q_{\text{vac}}(L') = \tfrac12\eta(L') = 0$: the new vacuum is globally neutral, as Ch. 10 proved. But it is not *locally* neutral: at $\delta \ne 0$ the CP phase makes the two walls carry equal-and-opposite vacuum charge layers, $+q(\delta)$ piled against one wall, $-q(\delta)$ against the other, net zero. The sudden expansion strands the old wall-layer pattern away from the new walls, and the truncated Bogoliubov sums — patiently, convergently — measure the *mismatch between the two layouts over the old region*. A polarization. Odd in $\delta$ (the layers swap sign with the phase — hence the perfect $\sin 2\delta$ behaviour and the machine-precision ECPT oddness $\Delta Q(+\delta) = -\Delta Q(-\delta)$, both genuine properties *of this observable*). Finite. Convergent. And worth exactly zero baryons.

![Wall charge layers](figures/ch11_fig1_wall_charge_density.png)

*Figure 11.1 — The physics the sums were measuring. Regulated vacuum charge density $\rho_{\text{vac}}(x)$ in the bulk-phase bag at $\delta \neq 0$: equal-and-opposite layers hugging the two walls, zero net. After sudden expansion (dashed: new-vacuum layout) the old region $[0, L]$ contains an unbalanced slice of the new pattern — the Polarization Identity integrand. (New figure; `ch11_wall_density.py`.)*

## 11.3 The convergence study: every clause verified

The identity makes five falsifiable claims: the two truncated forms should (i) disagree at small basis, (ii) close on *each other* as the basis grows, (iii) converge to the independently computed polarization integral (11.2), while (iv) the whole-box spectral asymmetry stays identically zero, and (v) the limit is odd in $\delta$. All five hold. **[Computed]** At the benchmark point $\delta = \pi/4$, $s = 1.2$, $mL = 1$ (`ch11_polarization.py`, complete spectra, heat-kernel regulated, Richardson extrapolated):

| quantity | $n_{\max} = 18$ | $N = 85$ | limit |
|---|---|---|---|
| $\Delta Q_{\text{cross}}$ | $+0.0562$ | $+0.0640$ | $\to +0.0672$ |
| $\Delta Q_{\text{same}}$ | $+0.0689$ | $+0.0678$ | $\to +0.0672$ |
| polarization integral $\tfrac12\sum_k \operatorname{sgn}(E'_k) w_k$ | — | — | $+0.0672$ |
| whole-box $\eta(L'),\ \eta(L)$ | $0$ | $0$ | $0$ (machine) |

The historically quoted value at this point ($+0.069$) is recognizably the same-branch entry of the first column: a truncation snapshot of a polarization integral, taken on its way to $+0.0672$. The $20\%$ inter-form gap at $n_{\max} = 18$ is the non-unitarity defect of a rank-18 truncation chasing a slowly converging *boundary* integral — the $w_k$ weights decay like the sinc tails of Ch. 3, and Ch. 5's boxed principle applies with full force: an observable that needs the cutoff to mean anything is not the observable you wanted.

![Convergence anatomy](figures/ch11_fig2_convergence.png)

*Figure 11.2 — The two truncated forms versus basis size at the benchmark point, closing on each other and on the polarization integral (horizontal line), with the exact net-charge value (zero, by theorem) marked. Inset: $|{\Delta Q_{\text{cross}} - \Delta Q_{\text{same}}}|$ decaying as the unitarity defect.*

![Oddness in the phase](figures/ch11_fig3_ecpt_oddness.png)

*Figure 11.3 — $\Delta Q^{\text{trunc}}(\delta) + \Delta Q^{\text{trunc}}(-\delta)$ across the phase range: zero at machine precision. The ECPT-oddness identity of Ch. 7 is real — and is here a property of the polarization observable, illustrating that passing symmetry checks does not certify an observable's identity.*

## 11.4 The deeper instability: "empty fresh region" is not a state

One more layer down sits a conceptual landmine worth defusing publicly, because it threatens *any* naive expanding-box setup, not just this one. The initial condition "the old vacuum in $[0, L]$, and *no fermions* in the fresh region $[L, L']$" is, for a Dirac field, **not definable**. A literally empty interval — every level of the region unoccupied, positive and negative alike — differs from any filled-sea reference by one hole per negative level: its charge relative to every admissible vacuum diverges linearly with the mode cutoff. (The same statement in Ch. 8's language: "empty" is not a finite-energy, finite-charge state of a fermion field; seas are not optional.) The truncated sums of §11.1 produce finite answers only because the symmetric mode cutoff *implicitly half-fills* the fresh region — a regularization choice smuggled in by the truncation, which is precisely the kind of dependence Corollary 9.3 warned would disqualify an observable. The only quench data that survive this landmine are the cutoff-protected spectral quantities: $\eta$, its difference, and the level crossings that change it.

## 11.5 Rules of practice (final form)

The arc of Ch. 9–11, compressed into the three rules that govern every remaining computation in this thesis:

1. **Net charge = spectral flow. Always.** Compute complete spectra and take $-\tfrac12\Delta\eta$, or count zero crossings along the deformation. Never report a truncated branch-sum residue as charge — *even when it converges*: the limit exists and is the polarization (11.2), not the charge.
2. **Bogoliubov sums are for pair content** — spectra of created pairs, energy budgets, distributions — where they are honest and convergent.
3. **Run the control.** Any claimed charge-production mechanism must be re-run at parameters where the Master Theorem forces zero. If the pipeline still reports a nonzero value there, the pipeline is measuring polarization. (Ch. 12 institutionalizes this: its equal-angle control quench returns spectral flow $0$ while the truncated forms still read $\approx -0.06$ — the diagnostic working exactly as designed.)

The negative arc of Part II is complete: bulk CP phases are inert (Ch. 10), and the numbers that suggested otherwise are wall polarization (this chapter). What remains is the constructive question, and by now the theory itself has told us where to look — Ch. 10 closed every door except one: *the walls must differ*.

---

**Validation.** `ch11_polarization.py`: reproduces $A \approx 0.0767$ at the historical parameter point; the convergence table and Fig. 11.2; the whole-box $\eta$ zeros; `--oddness` produces Fig. 11.3 (machine-zero check across $\delta$). `ch11_wall_density.py` (new): regulated vacuum charge density profiles before/after expansion (Fig. 11.1), plus the regional integral matching the polarization value $+0.0672$ independently.
