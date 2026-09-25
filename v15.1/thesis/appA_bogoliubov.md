# Appendix A — Bogoliubov transformations in finite volume

---

This appendix is the machine room for every "re-expand the state in the new basis" step of the thesis. It develops the Bogoliubov apparatus in the finite-volume, two-sector setting the framework actually uses — scalar warm-up, Dirac case, the exact unitarity identities, and the quantitative anatomy of how truncation breaks them. The last section is the one to read before re-deriving anything in Chapter 11.

## A.1 Scalar warm-up

Let the box change $L \to L' > L$, with old real-scalar modes $u_n(x)$ (frequencies $\omega_n$) on $[0, L]$ and new modes $u'_k$ ($\omega'_k$) on $[0, L']$. The *same* field operator can be expanded in either family (old modes extended by zero, per Ch. 3):

$$\hat\phi = \sum_n \frac{1}{\sqrt{2\omega_n}}\big(\hat a_n u_n + \hat a_n^\dagger u_n^*\big) = \sum_k \frac{1}{\sqrt{2\omega'_k}}\big(\hat a'_k u'_k + \hat a'^\dagger_k u'^*_k\big),$$

and matching (together with the conjugate momentum) expresses new operators through old:

$$\hat a'_k \;=\; \sum_n \Big(\alpha_{kn}\,\hat a_n \;+\; \beta_{kn}\,\hat a^\dagger_n\Big), \tag{A.1}$$

with coefficients fixed by the mode overlaps and the frequency mismatch:

$$\alpha_{kn} = \tfrac12\Big(\sqrt{\tfrac{\omega'_k}{\omega_n}} + \sqrt{\tfrac{\omega_n}{\omega'_k}}\Big)\,\langle u'_k, u_n\rangle, \qquad \beta_{kn} = \tfrac12\Big(\sqrt{\tfrac{\omega'_k}{\omega_n}} - \sqrt{\tfrac{\omega_n}{\omega'_k}}\Big)\,\langle u'_k, u_n^*\rangle, \tag{A.2}$$

$\langle\cdot,\cdot\rangle$ the overlap on $[0, L']$. Structure to internalize: $\beta \neq 0$ exactly when frequencies mismatch between the bases — $\beta$ mixes creation into annihilation, so the old vacuum ($\hat a_n|0\rangle = 0$) is *not* annihilated by $\hat a'_k$:

$$\langle 0_{\text{old}}|\,\hat N'_k\,|0_{\text{old}}\rangle \;=\; \sum_n |\beta_{kn}|^2 \tag{A.3}$$

— the new observer counts quanta in the old vacuum: **pair production by geometry change**, the finite-volume skeleton of Parker's cosmological mechanism (Ch. 16.8). Formally, the old vacuum is a squeezed state over new pairs, $|0_{\text{old}}\rangle \propto \exp\big(\tfrac12\sum_{kk'}C_{kk'}\hat a'^\dagger_k \hat a'^\dagger_{k'}\big)|0_{\text{new}}\rangle$ with $C = \beta\,(\alpha^*)^{-1}$ (impose $\hat a_n|0_{\text{old}}\rangle = 0$ with the inverted transformation $\hat a_n = \sum_k(\alpha^*_{kn}\hat a'_k - \beta_{kn}\hat a'^\dagger_k)$; the condition is $C\alpha^* = \beta$ — single mode: $C = \beta/\alpha^* = e^{i\varphi}\tanh r$ for $\alpha = \cosh r$, $\beta = e^{i\varphi}\sinh r$). Pairs, never singles: a bilinear Hamiltonian cannot create an odd number of quanta, which is the scalar shadow of the charge bookkeeping that Part II's fermions sharpen into theorems.

## A.2 The Dirac case

For fermions the branch structure doubles the bookkeeping. Old eigenmodes $\psi_n$ split into positive- and negative-energy families; likewise new $\psi'_k$. Expanding the *same* field two ways (Ch. 8.6 conventions: $\hat b$ annihilate particles, $\hat d^\dagger$ create antiparticles = annihilate sea states),

$$\hat\psi = \sum_{n: E_n > 0}\hat b_n\psi_n + \sum_{n: E_n < 0}\hat d^\dagger_n\psi_n = \sum_{k: E'_k > 0}\hat b'_k\psi'_k + \sum_{k: E'_k < 0}\hat d'^\dagger_k\psi'_k,$$

orthonormality of the new basis projects out the transformation, whose coefficients are *pure overlaps* (no frequency dressing — first-order systems carry their normalization in the spinors):

$$\alpha_{kn} \;=\; \int_0^{L}\psi'^\dagger_k(x)\,\psi_n(x)\,dx \tag{A.4}$$

(integration stops at $L$: old modes are extended by zero — the embedding map of Ch. 3, which is the identification-map declaration this computation rests on; change the map and you change every number, the lesson of Ch. 10). Organizing by branches gives the four blocks

$$\hat b'_k = \sum_{n>0}\alpha^{++}_{kn}\hat b_n + \sum_{n<0}\alpha^{+-}_{kn}\hat d^\dagger_n, \qquad \hat d'^\dagger_k = \sum_{n>0}\alpha^{-+}_{kn}\hat b_n + \sum_{n<0}\alpha^{--}_{kn}\hat d^\dagger_n, \tag{A.5}$$

with the cross blocks $\alpha^{+-}, \alpha^{-+}$ playing $\beta$'s role: they populate the old vacuum with new particles and antiparticles,

$$\langle N'^{\,b}_k\rangle_{0_{\text{old}}} = \sum_{n<0}\big|\alpha^{+-}_{kn}\big|^2, \qquad \langle N'^{\,d}_k\rangle_{0_{\text{old}}} = \sum_{n>0}\big|\alpha^{-+}_{kn}\big|^2. \tag{A.6}$$

**Charge conjugation ties the blocks.** When both problems are C-covariant (Ch. 8.2: $\psi_c = \sigma_1\psi^*$ maps each eigenbasis onto itself with $E \to -E$), applying C inside (A.4) gives $\alpha^{-+}_{\bar k\bar n} = \big(\alpha^{+-}_{kn}\big)^*$ (bars: C-partners) — the two cross blocks are conjugate, their squared sums equal, and the *net* charge in (A.6) cancels pairwise. This is the operator-level reason mirror-symmetric situations produce pairs but never net charge — the Master Theorem's shadow at the level of coefficients.

## A.3 Unitarity and completeness identities

The exact identities, with their domains tracked carefully — the bookkeeping Ch. 11's Polarization Identity exploits.

**(i) Rows: old modes are complete targets.** Each old mode lives on $[0, L] \subset [0, L']$ and the new family is complete on $[0, L']$; Parseval gives

$$\sum_{k\ (\text{all branches})} \big|\alpha_{kn}\big|^2 \;=\; \|\psi_n\|^2 \;=\; 1 \qquad \text{for every } n. \tag{A.7}$$

**(ii) Columns: new modes are *not* complete targets.** The reverse sum runs the other way around the asymmetry:

$$\sum_{n\ (\text{all branches})} \big|\alpha_{kn}\big|^2 \;=\; \int_0^{L}\big|\psi'_k\big|^2\,dx \;\equiv\; w_k \;\le\; 1, \tag{A.8}$$

the *weight of the new mode inside the old region* — strictly less than 1 whenever $\psi'_k$ leaks into the fresh region $[L, L']$. The pair (A.7)/(A.8) is the precise content of "the old basis is complete on its domain; the old domain is not the whole new domain", and the $w_k$ of (A.8) are exactly the weights in Ch. 11's closed form $\Delta Q_{\text{cross}}(\infty) = \tfrac12\sum_k\operatorname{sgn}(E'_k)\,w_k$.

**(iii) Operator algebra — one direction only.** Inverting (A.5) expresses the *old* operators through the new ones, $\hat b_n = \sum_{k:\,E'_k>0}\alpha^*_{kn}\,\hat b'_k + \sum_{k:\,E'_k<0}\alpha^*_{kn}\,\hat d'^\dagger_k$ (and likewise $\hat d^\dagger_n$), and *this* map preserves the CAR exactly at infinite rank: $\{\hat b_n, \hat b^\dagger_{n'}\} = \sum_k \alpha^*_{kn}\,\alpha_{kn'} = \langle\psi_n, \psi_{n'}\rangle = \delta_{nn'}$ — the row identity (A.7), completeness of the new family on the whole new box. The reverse direction *cannot* close: computing $\{\hat b'_k, \hat b'^\dagger_{k'}\}$ from the old-mode expansion gives $\sum_n \alpha_{kn}\,\alpha^*_{k'n} = \langle\psi'_k|\,\mathbb 1_{[0,L]}\,|\psi'_{k'}\rangle$, a matrix whose diagonal is $w_k < 1$ by (A.8) — the old modes span only $L^2[0,L]$, and no algebra over them can manufacture operators for the fresh region. (This is §11.4's "an empty fresh region is not a state", seen from the operator side.) The change of basis is unitary on the *whole* new box; restricted to old-box data it is an isometry in one direction only, and every honest identity below is therefore written with old indices outside and new indices summed.

## A.4 Truncation: how unitarity fails at finite rank, and what it costs

Keep $N$ modes per branch — every numerical computation does — and define the defect operator on the *old* modes, in the direction where unitarity actually holds (§A.3, point iii):

$$\Delta^{(N)}_{nn'} \;\equiv\; \delta_{nn'} \;-\; \sum_{|k| \le N}\alpha_{kn}\,\alpha^*_{kn'}. \tag{A.9}$$

Its size is set by the discarded tails. From the Geometric Overlap Theorem (Ch. 3, eq. 3.8), bag-mode overlaps decay as the sinc envelope $|\alpha_{kn}| \sim c/|k - rn|$ at large index ($r$ the size ratio), so the discarded weight per old mode scales as $\sum_{|k|>N}|\alpha_{kn}|^2 \sim \mathcal O(1/N)$:

$$\big\|\Delta^{(N)}\big\| \;\sim\; \frac{C}{N} \tag{A.10}$$

— *slow*, harmonic-series-tail slow, which is the quantitative root of every Ch. 11 phenomenon. **[Computed]** At the Ch. 11 benchmark ($\delta = \pi/4$, $s = 1.2$, $mL = 1$): the maximal diagonal of $\Delta^{(N)}$ over the six lowest old modes falls $0.0130 \to 0.0070 \to 0.0032 \to 0.0022$ at $N = 10, 20, 40, 58$ — the $1/N$ law in action — while the wrong-way-round object $\delta_{kk'} - \sum_{n}\alpha_{kn}\alpha^*_{k'n}$ *saturates* at $1 - w_k \approx 0.09$–$0.20$, exactly as (A.8) demands (`ch11_polarization.py --defect`). The index direction is not pedantry; it is the difference between a vanishing error and a theorem-sized one. Two consequences of (A.10), both used as diagnostics in the main text:

**(a) Equivalent expressions split by the defect.** Any two formulas for the same observable that are equal *via* the row identity (A.7) — e.g. the cross-branch and same-branch charge forms of Ch. 11, whose infinite-basis difference is exactly $\eta(L)$ and hence zero *given a mirror-symmetric old spectrum* (Ch. 11.1) — differ, at rank $N$, by matrix elements of $\Delta^{(N)}$:

$$\big|\,\mathcal O_{\text{form 1}}^{(N)} - \mathcal O_{\text{form 2}}^{(N)}\big| \;\le\; \|\hat{\mathcal O}\|\cdot\big\|\Delta^{(N)}\big\| \;\sim\; \frac{C'}{N}. \tag{A.11}$$

Inverted, this is the **smoking-gun inequality**: a 20% split between equivalent forms at $N = 18$ is not noise to average — it is a *measurement of the defect*, announcing that the observable being chased converges no faster than the unitarity sum itself. Convergent-looking values under such a split are convergent to *something* (Ch. 11 computes what); they are not the unitarity-protected observable.

**(b) The $1/N$ residual law.** The same tails make truncated Parseval sums approach their limits as $1/N$ — precisely the scaling measured in the Golden-Rule validations (Ch. 6: residuals $0.00133, 0.00137, 0.00157$ at $N = 220$, consistent with $C/N$ at $C \sim 0.3$) and visible across the convergence studies of Ch. 11. When a finite-volume fermionic computation converges faster than $1/N$, something is being summed analytically; when slower, an identification map is inconsistent (Ch. 10's lesson). The $1/N$ law is the healthy baseline.

**Rules of practice** (the operational summary consumed by Ch. 9.5 and 11.5): compute *pair content* from truncated sums freely — (A.3)/(A.6) converge at the honest $1/N$ rate to unitarity-protected limits; compute *net charge* only spectrally — the difference of cross-block sums is a difference of $\mathcal O(1)$ quantities with $\mathcal O(1/N)$ errors *each*, chasing a limit that the Polarization Identity shows is not charge at all; and always run both equivalent forms — their split is a free, built-in error meter (A.11).

---

**Validation.** The identities (A.7)–(A.8) and the $1/N$ law are exercised numerically throughout: `ch06_golden_rule.py` (Parseval convergence with measured $1/N$ residuals), `ch11_polarization.py` (the two-form split and its closure; `--defect` prints the (A.9) scaling and the wrong-direction saturation quoted in A.4), `ch08_bag_spectrum.py` (basis completeness checks underlying (A.4)).
