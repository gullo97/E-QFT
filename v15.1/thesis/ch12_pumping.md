# Chapter 12 — The mechanism: wall-localized CP violation and quantized charge pumping

---

Three chapters of demolition have left exactly one door standing. Net charge requires spectral asymmetry (Ch. 9). A bulk CP phase cannot create it — it is a common shift of both wall angles in disguise, hence a mass renormalization (Ch. 10) — and the mode sums that pretended otherwise were measuring wall polarization (Ch. 11). The one parameter no bulk redefinition can reach is the **difference between the wall angles**. This chapter opens that door and finds behind it everything a baryogenesis mechanism could ask for: an exactly solvable spectrum governed by a single master equation; vacuum charge obeying a fractional law with a four-decade pedigree; level crossings at analytically computable box sizes; and — the payoff — *quantized*, anomaly-protected charge pumping under expansion, verified end to end, with the control experiment that Ch. 11's rules demand coming out exactly as the theorems require.

## 12.1 The two-wall chiral bag

The arena: a Dirac fermion of **real** mass $m$ (the chiral rotation of Ch. 10 has been spent; no bulk phase remains) on $[0, L]$, with independent chiral wall angles,

$$-\,i\,n\!\cdot\!\gamma\;e^{i\theta_0\gamma_5}\,\psi = \psi \;\;\text{at } x = 0, \qquad -\,i\,n\!\cdot\!\gamma\;e^{i\theta_L\gamma_5}\,\psi = \psi \;\;\text{at } x = L. \tag{12.1}$$

The two combinations that organize everything:

$$\Delta \;=\; \theta_L - \theta_0 \quad (\text{the mismatch}), \qquad \Sigma \;=\; \frac{\theta_0 + \theta_L}{2} \quad (\text{the mean}). \tag{12.2}$$

Chapter 10's verdict in this notation: a bulk phase shifts $\Sigma$ and leaves $\Delta$ untouched. Whatever physics lives in $\Delta$ is therefore beyond the reach of any bulk CP structure — it must be put on the walls by the microphysics (what *sets* $\Delta$ is Ch. 14's question; here we take $\theta_0 \ne \theta_L$ and solve).

**Boundary rays.** Each condition in (12.1) is a projector equation ($P\psi = \psi$ with $P^2 = \mathbb 1$), confining $\psi$ at the wall to a one-dimensional ray. Work them out once, completely. At $x = L$ ($n = +\hat x$, $\gamma^1 = i\sigma_2$, $\gamma_5 = \sigma_1$): $-i\gamma^1 e^{i\theta_L\gamma_5} = \sigma_2(\cos\theta_L + i\sigma_1\sin\theta_L) = \cos\theta_L\,\sigma_2 + \sin\theta_L\,\sigma_3$ (using $\sigma_2\sigma_1 = -i\sigma_3$). The $+1$ eigenvector of $\cos\theta_L\,\sigma_2 + \sin\theta_L\,\sigma_3$ is found by the half-angle rotation it suggests:

$$v_L \;\propto\; \begin{pmatrix}\cos b \\ i\,\sin b\end{pmatrix}, \qquad b \;=\; \frac{\pi}{4} - \frac{\theta_L}{2}; \tag{12.3}$$

*(check: $\sigma_2 v_L = (\sin b,\, i\cos b)^T$ and $\sigma_3 v_L = (\cos b,\, -i\sin b)^T$, so the operator gives $(\sin(b{+}\theta_L),\, i\cos(b{+}\theta_L))^T$, which equals $v_L$ iff $b + \theta_L = \tfrac\pi2 - b$.)* At $x = 0$ ($n = -\hat x$) the operator is $-(\cos\theta_0\,\sigma_2 + \sin\theta_0\,\sigma_3)$, with $+1$ eigenvector

$$v_0 \;\propto\; \begin{pmatrix}\cos a \\ -\,i\,\sin a\end{pmatrix}, \qquad a \;=\; \frac{\pi}{4} + \frac{\theta_0}{2}. \tag{12.4}$$

MIT check ($\theta = 0$): $a = b = \pi/4$ reproduces the rays $(1, -i)$ and $(1, +i)$ of §8.3. **[Computed]** ray formulas verified against direct eigendecomposition in `ch12_chiral_bag.py`.

## 12.2 The Chiral-Wall Master Equation

With real mass the bulk system (Ch. 10, eq. 10.5 at $m_I = 0$) is $\psi_1' = i(E + m)\psi_2$, $\psi_2' = i(E - m)\psi_1$, solved by

$$\psi_1 = A\cos px + B\sin px, \qquad \psi_2 = -\,\frac{i\,p}{E + m}\,\big(B\cos px - A\sin px\big), \qquad p^2 = E^2 - m^2 . \tag{12.5}$$

**Wall at $0$:** the ray (12.4) demands $\psi_2(0)/\psi_1(0) = -i\tan a$, giving $-\tfrac{ipB}{(E+m)A} = -i\tan a$, i.e.

$$B \;=\; \frac{(E + m)\tan a}{p}\,A. \tag{12.6}$$

**Wall at $L$:** the ray (12.3) demands $\psi_2(L)/\psi_1(L) = +i\tan b$:

$$-\,p\,\big(B\cos pL - A\sin pL\big) \;=\; (E + m)\tan b\,\big(A\cos pL + B\sin pL\big).$$

Insert (12.6), multiply through by $p\cos a\cos b/(E+m)$, and collect:

$$p\,(E+m)\,\sin(a + b)\,\cos pL \;=\; \sin pL\,\Big[\,p^2\cos a\cos b \;-\; (E+m)^2 \sin a \sin b\,\Big].$$

The right-hand bracket simplifies via $p^2 = (E+m)(E-m)$:

$$p^2\cos a\cos b - (E{+}m)^2\sin a\sin b = (E{+}m)\big[E\cos(a{+}b) - m\cos(a{-}b)\big],$$

and the common factor $(E + m)$ cancels. Now translate the angle combinations: from (12.3)–(12.4), $a + b = \tfrac\pi2 - \tfrac{\Delta}{2}$ and $a - b = \Sigma$, so $\sin(a{+}b) = \cos\tfrac\Delta2$ and $\cos(a{+}b) = \sin\tfrac\Delta2$:

> **Chiral-Wall Master Equation [Theorem].**
>
> $$\boxed{\;p\,\cos\!\Big(\frac{\Delta}{2}\Big)\,\cos(pL) \;=\; \sin(pL)\,\Big[\,E\,\sin\!\Big(\frac{\Delta}{2}\Big) \;-\; m\,\cos\Sigma\,\Big]\;,\qquad p^2 = E^2 - m^2,} \tag{12.7}$$
>
> with the below-gap branch ($|E| < m$) obtained by $p = iq$: $\;q\cos\frac\Delta2\,\cosh(qL) = \sinh(qL)\big[E\sin\frac\Delta2 - m\cos\Sigma\big]$.

**[Computed]** (12.7) verified against the independent transfer-matrix solver for five wall-angle pairs to $7\times10^{-8}$ (the solver's tolerance); the $(\theta_0, \theta_L) = (\delta, \delta)$ real-mass spectrum matches the bulk-complex-mass solver of Ch. 10 to $8\times10^{-8}$, closing the consistency loop promised in §10.6 (`ch12_chiral_bag.py`).

One transcendental equation now contains the entire phenomenology of Part II. Read it clause by clause.

**(a) The mismatch is the only charge knob.** The single term that knows the sign of $E$ carries $\sin(\Delta/2)$. Under $E \to -E$ (with $p \to p$), equation (12.7) maps to itself with $\Delta \to -\Delta$:

> **Wall-Mismatch Criterion [Theorem ("if"); Computed ("only if")].** If $\sin(\Delta/2) = 0$, the spectrum is $E \leftrightarrow -E$ symmetric — hence $\eta = 0$, hence $Q_{\text{vac}} = 0$, hence no net production. Conversely, $\sin(\Delta/2) \ne 0$ breaks the mirror *generically* — exactly in the massless limit ($\eta = \Delta/\pi$, Ch. 9) and at every sampled massive point — though the symmetry argument alone does not exclude an accidental mirror at isolated $(\Sigma, mL)$. Spectral asymmetry exists exactly when the walls disagree.

**(b) $\Delta = 0$ contains all of Chapter 10 in one line.** Set $\Delta = 0$: (12.7) collapses to $\tan(pL) = -\,p/(m\cos\Sigma)$. A common wall angle — and therefore also a bulk phase, which by (10.3) is the case $\Sigma = -\delta$, $\Delta = 0$ after the chiral rotation (only $\cos\Sigma$ enters, so the sign is spectrally invisible here — but it is the sign the Goldstone–Wilczek toolbox below lives or dies by) — merely renormalizes the mass, $m \to m\cos\Sigma$. The Dressing Theorem, the Spectral Mirror Theorem, and the Zero-Charge Theorem are the three shadows of this one specialization.

**(c) The gap can now be inhabited.** Unlike the $\Delta = 0$ case (where §10.4 proved the gap empty), the below-gap branch of (12.7) admits solutions once $\sin(\Delta/2) \ne 0$ — levels can *enter the gap and approach zero*. That is the geometric origin of everything in the next two sections.

## 12.3 Fractional vacuum charge: the $-\Delta/2\pi$ law

With $\Delta \neq 0$ the spectrum tilts, and the vacuum charges up. How much? Complete-spectrum numerics first, theory second.

**[Computed]** (`ch12_chiral_bag.py`): symmetric walls $(\theta_0, \theta_L) = (-\tfrac\Delta2, +\tfrac\Delta2)$ (so $\Sigma = 0$), $mL = 6$ (walls decoupled), all $\approx 1340$ levels below $E_{\max} = 350$, heat-kernel regulated and extrapolated:

$$Q_{\text{vac}} \;\approx\; -\,\frac{\Delta}{2\pi} \pmod 1 \tag{12.8}$$

across the full range tested — e.g. $Q_{\text{vac}} = -0.4455$ at $\Delta = 2.8$ against $-\Delta/2\pi = -0.4456$, residuals at the few-$10^{-4}$ level being the regulator extrapolation error and shrinking with it.

> **Fractional Charge Law [Computed + Standard].** Each chiral wall binds vacuum charge $\mp\theta/2\pi$ with the sign set by its orientation — $-\theta_L/2\pi$ at the right wall, $+\theta_0/2\pi$ at the left — so the two-wall bag carries $Q_{\text{vac}} = -\Delta/2\pi$ modulo integer rearrangements of the sea.

This is the 1+1D realization of the chiral-bag fractional fermion number — Goldstone and Jaffe's resolution of the baryon-number bookkeeping in hybrid bag models (Ch. 9's anchor, now carrying live weight). And there is a field-theory derivation of the coefficient worth having in full view:

> **Toolbox: the Goldstone–Wilczek current.** Let the *bulk mass angle* vary slowly in space, $m\,e^{i\delta(x)\gamma_5}$ (a wall is the sharp limit). Integrating out the heavy fermion induces, at leading order in gradients, the current
>
> $$\langle j^\mu\rangle \;=\; \frac{1}{2\pi}\,\epsilon^{\mu\nu}\,\partial_\nu\,\delta(x) \qquad \text{[Standard]}$$
>
> (the 1+1D Goldstone–Wilczek current; derivable by point-splitting, by bosonization — where $\delta$ shifts the boson and $j^\mu = \tfrac{1}{2\pi}\epsilon^{\mu\nu}\partial_\nu\phi$ makes it one line — or as the descent of the 2D anomaly). Now the one sign that matters: by (10.3), a bulk angle $\delta$ is spectrally equivalent to wall angles $\theta = -\delta$ — so the bulk profile whose sharp limit is the wall pair $(\theta_0, \theta_L)$ is $\delta(x)$ running from $-\theta_0$ to $-\theta_L$, and the charge between the walls is
>
> $$Q \;=\; \int_0^L \frac{\delta'(x)}{2\pi}\,dx \;=\; \frac{-\theta_L - (-\theta_0)}{2\pi} \;=\; -\,\frac{\Delta}{2\pi}$$
>
> up to the integers contributed by levels crossing zero — precisely the law (12.8), *including its sign*, found by brute-force spectral summation. (The relative minus between bulk and wall angles is not decoration; it is the same $\theta = -\delta$ map that clause (b) of §12.2 used, and consistent with it the bulk-phase polarization of Ch. 11 ($+0.067$) and the equal-wall control of §12.5 ($-0.06$) carry opposite signs.) The agreement between a regulated 1340-level sum and a one-line anomaly argument is the kind of overdetermination this thesis aims for everywhere.

![The fractional law](figures/ch12_fig1_qvac_law.png)

*Figure 12.1 — Vacuum charge vs wall mismatch: complete-spectrum $Q_{\text{vac}}$ (points) against $-\Delta/2\pi$ (line), $\Sigma = 0$, $mL = 6$. Residuals at the few-$10^{-4}$ level, attributable to the heat-kernel extrapolation.*

## 12.4 Level crossings: the Crossing Condition

The "$\bmod\ 1$" in (12.8) is where the integers live, and the integers are the baryons. Where exactly does the sea gain or lose a level? Set $E = 0$ in the below-gap branch of (12.7) ($E = 0 \Rightarrow q = m$):

> **Crossing Condition [Theorem].** A zero-energy level exists exactly when
>
> $$\boxed{\;\tanh(mL^*) \;=\; -\,\frac{\cos(\Delta/2)}{\cos\Sigma}\;.} \tag{12.9}$$
>
> Since $\tanh : (0,\infty) \to (0,1)$ is a bijection, the two-wall bag has **at most one** crossing size $L^*$, existing iff $\;0 < -\cos(\Delta/2)/\cos\Sigma < 1$, i.e. iff $\cos(\Delta/2)$ and $\cos\Sigma$ have opposite signs and $|\cos(\Delta/2)| < |\cos\Sigma|$.

### 12.4.1 Anatomy of the crossing: the zero mode built by hand

Equation (12.9) came out of the master equation in one substitution, which is efficient but opaque: it answers *where* the crossing is without showing *what* is crossing, or why a bag should single out one special size at all. Both questions have completely explicit answers, and they are worth having because "why is there a specific $L^*$?" is the question every reader should ask of a mechanism whose output is an integer per crossing.

**At $E = 0$ the Dirac system decouples.** Set $E = 0$ in the bulk equations (12.5-pre): $\psi_1' = im\,\psi_2$, $\psi_2' = -im\,\psi_1$. The combinations $u_\pm = \psi_1 \pm i\psi_2$ diagonalize the system,

$$u_\pm' \;=\; \pm\, m\,u_\pm \qquad\Longrightarrow\qquad u_\pm(x) \;=\; u_\pm(0)\,e^{\pm mx} : \tag{12.10}$$

a zero-energy solution is *always* a superposition of one exponentially growing and one exponentially decaying profile, with decay length $1/m$ — the fermion's Compton wavelength. (The $u_\pm$ are the eigenvectors of $\sigma_2$; this is the 1+1D avatar of the Jackiw–Rebbi zero-mode structure: at $E=0$, and only at $E=0$, the two chiralities propagate independently.) There is no oscillation and no free momentum at $E = 0$; the *only* length scale available is $1/m$. Whatever condition the walls impose can therefore only be satisfied at isolated values of $mL$ — the dimensional reason a specific size exists.

**The walls fix the mixture.** In the $u_\pm$ variables the boundary rays (12.3)–(12.4) read, after two lines of half-angle algebra,

$$\frac{u_+(0)}{u_-(0)}\bigg|_{\text{ray at }0} = -\,\cot\frac{\theta_0}{2}, \qquad \frac{u_+(L)}{u_-(L)}\bigg|_{\text{ray at }L} = +\,\tan\frac{\theta_L}{2}.$$

Propagating (12.10) from $0$ to $L$ multiplies the ratio by $e^{2mL}$, so a zero mode exists iff $-\cot(\theta_0/2)\,e^{2mL} = \tan(\theta_L/2)$, i.e.

$$\boxed{\;e^{2mL^*} \;=\; -\,\tan\!\frac{\theta_0}{2}\,\tan\!\frac{\theta_L}{2}\;} \tag{12.11}$$

— an equivalent, and more talkative, form of the Crossing Condition (one line of product-to-sum identities recovers (12.9) exactly). It says three things at a glance. *Existence:* the right side must exceed $1$ (since $mL^* > 0$), which requires the two wall tangents to have opposite signs and product beyond $-1$ — order-one angles, the threshold of Ch. 14 again. *Uniqueness:* the left side is strictly monotone in $L$, so at most one size works — the tanh-bijection argument of (12.9), now visible as the simple fact that an exponential crosses any level once. *Scale:* $mL^* = \tfrac12\ln\big[-\tan(\theta_0/2)\tan(\theta_L/2)\big]$ — the crossing size is the Compton wavelength times a logarithm of the wall geometry; for the demonstration point $(-2, +2)$, $\tan^2(1) = 2.4256$ and $\tfrac12\ln 2.4256 = 0.44302$, the number the spectrum scan found.

**What is crossing: two wall-bound states hybridizing.** Solve each wall alone on a half-line (the $L \to \infty$ decoupling limit). An evanescent state $\psi \propto e^{-qx}$ at the left wall must have $q = -(m + E)\tan(\tfrac\pi4 + \tfrac{\theta_0}{2}) > 0$, which (squaring, with $q^2 = m^2 - E^2$) pins its energy and existence:

$$E_{\text{left}} \;=\; -\,m\sin\theta_0 \;\;\big[\text{exists iff } \tan(\tfrac\pi4 + \tfrac{\theta_0}{2}) < 0\big], \qquad E_{\text{right}} \;=\; +\,m\sin\theta_L \;\;\big[\text{exists iff } \tan(\tfrac\pi4 - \tfrac{\theta_L}{2}) < 0\big] \tag{12.12}$$

— each chiral wall is a Jackiw–Rebbi defect that binds at most one sub-gap fermion state, at an energy set by its own angle. An ordinary MIT wall ($\theta = 0$) binds nothing — the microscopic content of the one-ordinary-wall no-go below. For the demonstration point both walls bind, degenerately at $E = +m\sin 2 = +0.909\,m$ (the $\Sigma = 0$ symmetry); one may check that (12.12) solves the $L\to\infty$ limit of the below-gap master equation, $q\cos\tfrac\Delta2 = E\sin\tfrac\Delta2 - m\cos\Sigma$, exactly. At finite $L$ the two evanescent tails overlap with strength $\sim e^{-2mL}$ and the levels hybridize and repel; as the box shrinks toward the Compton scale the interference grows strong enough to drag the lower branch across the entire gap. The crossing at $L^*$ *is* this dressed wall-bound state passing $E = 0$ — the point at which, by charge conjugation, the state's membership switches between "sea" and "sky", which is why exactly one unit of vacuum charge moves and why the location is fixed by wall data plus $1/m$ and nothing else.

**Why $E = 0$ and not some other energy?** Because $E = 0$ is the charge-conjugation-symmetric point of the spectrum: $Q_{\text{vac}} = -\tfrac12\sum\operatorname{sgn}(E_k)$ changes by $\pm1$ exactly when a level's sign flips, and nowhere else. Levels may wander arbitrarily inside the gap without consequence; only the zero matters, and the zero is available only to a state whose two chiral exponentials can *simultaneously* satisfy both walls — condition (12.11).

**The cosmological reading, previewed.** Through the dictionary of Ch. 17, $L^* \sim 1/m$ says: a cell pumps when the expansion carries its proper size through the fermion's Compton wavelength. This is why the epoch window of the baryogenesis estimate (Ch. 14) is an order-one window in $mL$ — at the electroweak crossover, $m_t/T \sim 1$, the top quark's Compton scale sits exactly where the relevant cells are — and it is why the crossing-density amplitude there is geometric rather than fine-tuned. The specific length is not a new scale postulated into the model; it is the one scale the fermion already owned, read out by the geometry.

Two immediate payoffs, one explanatory and one strategic.

*The explanatory payoff:* with one ordinary MIT wall ($\theta_0 = 0$, so $\Delta = \theta_L$, $\Sigma = \theta_L/2$), the right side of (12.9) is $-\cos(\theta_L/2)/\cos(\theta_L/2) = -1$, **never attained**: a bag with even one standard wall has *no crossing at any size*. A numerical hunt that scans only such configurations is guaranteed empty-handed — a trap this program fell into once, documented here so no reader repeats it.

*The strategic payoff:* the existence region in the $(\Sigma, \Delta)$ plane is sharply bounded and *requires order-one angles* — e.g. at $\Sigma = 0$ one needs $|\Delta| > \pi$. Wall CP violation must be strong to pump; this single fact reorganizes the cosmological estimate (Ch. 14) around the *population* of strong-CP-wall domains rather than around perturbative phase counting.

**The demonstration point.** Choose $\theta = (-2, +2)$: $\Sigma = 0$, $\Delta = 4$, $\cos(\Delta/2) = \cos 2 = -0.4161$, so (12.9) predicts a crossing at

$$mL^* = \operatorname{artanh}(0.41615) = 0.44302 .$$

**[Computed]** The complete-spectrum scan finds the level nearest zero diving through $E = 0$ at exactly this size, with $Q_{\text{vac}}(L)$ sitting on the plateau $+0.3640$ below $L^*$ (consistent with $-4/2\pi \equiv +0.3634 \bmod 1$), the plateau $-0.6361$ above, and a jump of $-1.00015$ — **one unit of charge, located at the analytically predicted size** — the jump residual ($1.5\times10^{-4}$) and the plateau offsets from the exact fractions ($6\times10^{-4}$) both being regulator-extrapolation bias, shrinking with tighter Richardson settings (Fig. 12.2).

![The crossing](figures/ch12_fig2_crossing.png)

*Figure 12.2 — The pump's escapement. Top: the level nearest zero vs box size, crossing at the predicted $L^* = 0.44302/m$ (vertical line). Bottom: $Q_{\text{vac}}(L)$ jumping by one unit between its two fractional plateaus, exactly at the crossing.*

![Spectral flow waterfall](figures/ch12_fig3_spectral_flow_waterfall.png)

*Figure 12.3 — The full spectrum vs $L$ ("waterfall"): mirror-symmetric pairs everywhere except the single gap-crossing level that carries the spectral flow. As the box grows through $L^*$ (at $\Delta > 0$), one state climbs out of the sea into the positive branch — the definition of "vacuum" changes by exactly one slot, $Q_{\text{vac}}$ jumps by $-1$, and the expansion pumps net charge $+1$. (New figure; `ch12_waterfall.py`; flow direction verified in `ch14_crossing_density.py`.)*

![The pump, live](figures/ch12_anim_pump.gif)

*Animation 12.A — The pump in motion. Left: the spectrum as the box expands, the gap level (orange) making its crossing. Right: the vacuum charge holding its fractional plateau, then jumping by exactly one unit at the analytically predicted $L^*$. (`make_anims_v15.py`; plateau values from the Fractional Charge Law, jump verified in `ch12_pump_control.py`.)*

## 12.5 The pump, run end to end — with its control

Everything is now assembled for the experiment that Part II has been building toward. Quench the box across the crossing and apply the Master Theorem's accounting; then re-run with the mismatch switched off, as Ch. 11's Rule 3 demands.

**The pump.** Expand $L = 0.30 \to 0.70$ ($m = 1$, walls $(-2, +2)$), straddling $L^* = 0.443$:

$$\Delta Q_{\text{net}} \;=\; Q_{\text{vac}}(0.30) - Q_{\text{vac}}(0.70) \;=\; +1.0000 \qquad \textbf{[Computed]}$$

— one quantized unit of net charge (`ch12_pump_control.py`, Richardson-extrapolated regulator; a coarser extrapolation gives $+0.996$, the difference being pure regulator error — the quantized value is the extrapolation-stable one). The expansion dragged one level through zero; the universe of this toy is one fermion richer, with no antifermion partner.

**The control.** Same geometry, same quench, equal walls $(0.7, 0.7)$ (so $\Delta = 0$; the Mirror Theorem applies): spectral flow gives **exactly $0$** ($-0.0000$ measured, same script). And the truncated Bogoliubov forms? In the control they converge to $\approx -0.06$ — nonzero, the wall polarization of Ch. 11, present even when net production vanishes identically. In the crossing quench they scatter ($+0.65$ and $-0.62$ at $N = 75$ for the two forms), nowhere near the true $+1$: confronted with genuine spectral flow, the truncated forms fail *in both directions at once*. **The truncated forms never measure net charge; the spectral flow always does** — Part II's methodological verdict, now demonstrated in a single figure.

![Pump vs control](figures/ch12_fig4_pump_vs_control.png)

*Figure 12.4 — The decisive experiment. Left: the crossing quench — spectral flow lands on $+1$ (quantized), truncated forms scatter. Right: the control quench — spectral flow exactly $0$, truncated forms converge to the polarization value $\approx -0.06$. One mechanism, one artifact, cleanly separated.*

## 12.6 What kind of mechanism this is: anomaly inflow

The physics deserves its proper name, because the name carries protections. Charge production here is **spectral flow of the Dirac operator under a deformation of its boundary data** — the same structure as Callan–Harvey anomaly inflow: a charge-violating-looking process on a defect (here, the moving wall) compensated by a topological current in the bulk, with the produced charge counting level crossings, an integer protected against smooth deformations, regulator choices, and truncation artifacts. Three structural consequences:

1. **Quantization is not approximate.** The yield per crossing is exactly one unit; all model-dependence is relegated to *whether and how often* crossings occur (the crossing density of Ch. 14) — a clean separation of the protected from the contingent.
2. **CP violation must be spatially structured.** A wall-angle difference is, after the inverse chiral rotation, a *spatially varying* phase $\delta(x)$ — a CP **domain wall** through the bag. Charge pumping requires CP violation *localized on boundaries or interfaces*, not a uniform bulk phase. This is structurally the electroweak-baryogenesis configuration — CP-violating bubble walls sweeping through the plasma — reached here from the opposite direction: not postulated as a phase-transition by-product, but *forced* by the spectral-flow theorems as the only configuration that can pump.
3. **Isotropy is no longer the enemy.** The historic worry that spherical geometry kills the asymmetry (the κ-cancellation) dissolves: that cancellation belonged to the inert bulk-phase mechanism (§10.5). For the genuine mechanism, the question in 3+1D is whether the spherical chiral bag's κ-channels have crossings and at what density — anisotropy may modulate the rate but is not required for its existence. The 3+1D execution is open and precisely specified (Ch. 27, item 3).

## 12.7 Summary

Mismatched chiral walls do everything bulk phases could not: the **Master Equation** (12.7) organizes the entire model class; the **Wall-Mismatch Criterion** localizes CP-effectiveness in $\sin(\Delta/2)$; the vacuum charges by the **Fractional Charge Law** $-\Delta/2\pi$ (anomaly-derived, numerically confirmed); the **Crossing Condition** (12.9) places at most one level crossing per bag at an analytically known size, requiring order-one wall angles; and the quench across it pumps **exactly one quantized unit**, while the control quench pumps exactly none. The mechanism is anomaly inflow on a moving boundary — a sentence that Ch. 26 will be able to read in two languages at once, because by then the "moving boundary" will be a feature of an emergent spacetime.

What stands between this chapter and a number for the universe is bookkeeping at cosmological scale: how often does an expanding universe of such bags cross? That is Chapter 14 — after a short interlude (Ch. 13) supplying the standard cosmology it needs.

---

**Validation.** `ch12_chiral_bag.py`: boundary rays vs eigendecomposition ($2\times10^{-16}$); master equation vs transfer matrix ($7\times10^{-8}$, five angle pairs); the $(\delta,\delta)$/bulk-phase consistency check ($8\times10^{-8}$); the fractional-law scan (Fig. 12.1); the crossing locator and plateau/jump values (Fig. 12.2); the truncated Bogoliubov forms for the pump and control quenches. `ch12_pump_control.py`: the spectral-flow pump ($+1.0000$) and control ($-0.0000$), Fig. 12.4. `ch12_waterfall.py`: the spectrum-vs-$L$ waterfall (Fig. 12.3). Every number quoted above is printed by these scripts.
