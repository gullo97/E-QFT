# Appendix B — Heat kernel, zeta function, Seeley–DeWitt

> **Author: Fable 5 (absorbing the former Opus slot) · Complexity ◇3 · Status: draft v1**
> Consumed by: Ch. 9 (η and its jumps), Ch. 22 (induced gravity), Ch. 25 (Weyl counting and $C(s)$). All content **[Standard]**, rederived at consumption depth.

---

Three chapters of this thesis stand on regulated spectral sums, and each leans on a different face of the same machine. This appendix builds the machine once: regulated sums and their regulator-(in)dependence (§B.1), the spectral asymmetry's properties (§B.2), the heat-kernel curvature expansion (§B.3), and the Weyl counting function with the universal coefficient $1/48\pi^2s^2$ derived in full (§B.4).

## B.1 Regulated sums: the toolbox

**The problem in miniature.** The sum $\sum_{n\ge1} n$ diverges; physics keeps producing it (zero-point energies, Ch. 22) and extracting finite answers. The two standard regulators, run side by side on this example so their relationship is visible:

*Heat kernel.* Damp and expand:

$$S(s) = \sum_{n\ge1} n\,e^{-sn} = \frac{e^{-s}}{(1 - e^{-s})^2} = \frac{1}{s^2} - \frac{1}{12} + \mathcal O(s^2). \tag{B.1}$$

The $1/s^2$ term is the *bulk divergence* — non-universal, regulator-shaped, scaling with the cutoff; in Casimir-type problems it is the piece that cancels between geometries or is absorbed into the cosmological term (Ch. 22.4 owns the honesty about it). The constant $-\tfrac{1}{12}$ is the *universal finite part*: regulator-independent (replace $e^{-sn}$ by any smooth $f(sn)$, $f(0) = 1$, decaying fast: Euler–Maclaurin gives $\int_0^\infty f(sx)x\,dx/s^0$-type bulk pieces plus the same $-B_2/2 = -\tfrac{1}{12}$).

*Zeta function.* Define $\zeta(z) = \sum n^{-z}$ where convergent, continue analytically, evaluate: $\zeta(-1) = -\tfrac{1}{12}$. The zeta route *silently discards* the bulk term and hands over the universal part directly — elegant, but the heat kernel's explicit bulk term is physically informative (it *is* the cosmological-constant problem), so the thesis quotes both routes wherever the result is load-bearing (Ch. 22.1).

**Richardson extrapolation**, the workhorse of the numerical suite (App. C): given $F(t) = F_0 + c_1 t + c_2 t^2 + \cdots$ sampled at $t, t/2$, the combination $2F(t/2) - F(t) = F_0 + \mathcal O(t^2)$ cancels the leading regulator error; iterated over a ladder of $t$ values it accelerates power-law convergence to machine-practical precision. Every "$t \to 0$ extrapolated" value in Parts II–III means this.

**When is the limit regulator-independent?** The general statement used throughout: if the regulated sum's small-$s$ expansion has the form (bulk powers of $1/s$) + (constant) + (vanishing terms), the constant is the same for every even, smooth, sufficiently decaying regulator — the bulk powers' *coefficients* differ between regulators, the constant does not (Mellin-transform bookkeeping: the constant is the residue at $z = 0$, a property of the spectral zeta function alone). All universal claims in this thesis sit in such constants.

## B.2 The spectral asymmetry η

For a Dirac-type operator with spectrum $\{E_k\}$ (both signs, symmetric growth), define

$$\eta(s) \;=\; \sum_k \operatorname{sgn}(E_k)\,e^{-s|E_k|}, \qquad \eta \;=\; \lim_{s\to0^+}\eta(s). \tag{B.2}$$

**Why η is cleaner than the sums it is built from.** Pair the spectrum at large $|E|$: by Weyl asymptotics the positive and negative branches have, level by level, matching densities up to $\mathcal O(1/|E|)$ corrections, so in $\eta(s)$ the dangerous high-energy tails cancel *pairwise inside the sum* before any limit is taken: $\eta(s)$ has **no bulk divergence** — its small-$s$ expansion is constant + vanishing. By §B.1's criterion the limit is then regulator-independent for every even damping. (Verified three ways on the toy spectra of Ch. 9: heat, Gaussian, smoothed-sharp regulators extrapolating to the same value — `ch09_eta_toy.py`.) $\eta$ is an *infrared* quantity wearing ultraviolet clothes: it counts the spectrum's imbalance, which is decided near $E = 0$.

**Behaviour at level crossings.** Deform the operator along a parameter $\lambda$. Away from zeros of the spectrum, each $\operatorname{sgn}(E_k(\lambda))$ is constant and $\eta(\lambda)$ varies only through the (uniformly convergent, hence continuous) damped sum — and its limit is locally constant: paired tails keep cancelling, small level motion is relabeling. When a single level crosses $E = 0$, its sign flips:

$$\Delta\eta \;=\; \pm2 \;\;\text{per crossing}, \tag{B.3}$$

the jump direction set by the crossing direction. This staircase behaviour — flat plateaus, integer jumps — is Corollary 9.2's engine and the quantization mechanism of Ch. 12. (The uniformity claim deferred from Ch. 9: on any compact parameter interval avoiding $E = 0$ degeneracies, the gap $\min_k|E_k(\lambda)|$ is bounded below, the damped sum converges uniformly in $(s, \lambda)$, and limit and deformation commute — the plateau is exact, not approximate.)

**Context [Standard].** $\eta$ is the finite-volume cousin of the Atiyah–Patodi–Singer eta-invariant — the boundary correction in the index theorem for manifolds with boundary — and the fractional fermion number results anchoring Ch. 9 (Jackiw–Rebbi; Goldstone–Jaffe) are its physics debut. The thesis uses only the properties derived above.

## B.3 The heat-kernel expansion

The curvature machinery behind Ch. 22. For a Laplace-type operator $-\Box + \xi R$ on a $d$-dimensional curved space, the trace of the heat operator has the small-$s$ **Seeley–DeWitt expansion**

$$\operatorname{Tr}\,e^{-s(-\Box + \xi R)} \;=\; \frac{1}{(4\pi s)^{d/2}}\int\sqrt{g}\;\Big[a_0 + s\,a_1 + s^2 a_2 + \cdots\Big], \qquad a_0 = 1, \quad a_1 = \Big(\tfrac16 - \xi\Big)R, \tag{B.4}$$

with $a_2$ the curvature-squared layer. Where the $\tfrac16$ comes from, at derivation depth (covariant perturbation theory, leading order): expand the metric in normal coordinates about a point, $g_{ij} = \delta_{ij} - \tfrac13 R_{ikjl}x^kx^l + \cdots$; the flat heat kernel $e^{-|x-y|^2/4s}/(4\pi s)^{d/2}$ gets corrected by the quadratic potential the curvature terms induce; Gaussian moments $\langle x^kx^l\rangle = 2s\,\delta^{kl}$ convert the curvature contraction into $+\tfrac{s}{6}R$ at coincidence — the universal $\tfrac16$ is, at bottom, the second moment of a Gaussian against the normal-coordinate expansion. Minimal coupling ($\xi = 0$), the thesis's matter, keeps $a_1 = R/6$.

**Spin weights** (the $N_{\text{eff}}$ inputs of Ch. 22/27): each field species contributes to the induced $\int\sqrt g\,R$ term in proportion to (number of field components) × (its $a_1$ coefficient), with the standard per-species tallies — scalar ($\xi = 0$): $+\tfrac16$ per component; Dirac fermion: the squared operator $-\slashed D^2 = -\Box + \tfrac R4$ gives $\big(\tfrac16 - \tfrac14\big) = -\tfrac{1}{12}$ per spinor component, times $(-1)$ for the fermion loop, net $+\tfrac{1}{12}\times$(4 components) $= +\tfrac13$ per Dirac field; gauge vectors similarly with ghost subtractions. The weighted sum over the model's matter content is the finite calculation flagged as Ch. 27 item 2; the thesis quotes only its $\mathcal O(10^2)$ scale.

**From (B.4) to induced gravity** (the two-line core of Ch. 22.3): the one-loop effective action $\Gamma = -\tfrac12\int_0^\infty\tfrac{ds}{s}\operatorname{Tr}e^{-s(\cdots)}$ picks up, from the $a_1$ term with the proper-time integral cut at $s_{\min} = 1/\Lambda^2$,

$$\Gamma \;\supset\; -\,\frac{(\tfrac16 - \xi)\,\Lambda^2}{32\pi^2}\int\sqrt{-g}\,R \qquad\Longrightarrow\qquad G^{-1}_{\text{ind}} = \frac{N_{\text{eff}}\Lambda^2}{12\pi} \tag{B.5}$$

after matching to $-\tfrac{1}{16\pi G}\int\sqrt{-g}R$ and summing species — Sakharov's induced Newton constant, the $a_0$ term meanwhile contributing the $\Lambda^4$ cosmological piece (§22.4's standing sentence).

## B.4 Weyl counting and the universal coefficient $C(s) = 1/48\pi^2 s^2$

Chapter 25's parameter-free prediction, derived in full. On a curved 3-space, the integrated density of Laplacian eigenfrequencies ($\omega = \sqrt{\text{eigenvalue}}$) has the **Weyl expansion** — the counting-function face of (B.4):

$$N(\omega) \;=\; \frac{V\,\omega^3}{6\pi^2} \;+\; \frac{\omega}{24\pi^2}\int\sqrt g\,R^{(3)} \;+\; \cdots \tag{B.6}$$

*(Derivation of the second coefficient: Mellin-transform (B.4) at $d = 3$. The trace $\operatorname{Tr}e^{-s\Delta} = \int_0^\infty e^{-s\omega^2}\,dN(\omega)$; inserting $N = A\omega^3 + B\omega$ and using $\int_0^\infty e^{-s\omega^2}\omega^2 d\omega = \tfrac{\sqrt\pi}{4}s^{-3/2}$, $\int_0^\infty e^{-s\omega^2}d\omega = \tfrac{\sqrt\pi}{2}s^{-1/2}$ matches $(4\pi s)^{-3/2}\int\sqrt g\,[1 + \tfrac s6 R]$ term by term: $A = V/6\pi^2$ and $3A\cdot\tfrac{\sqrt\pi}{4} = \tfrac{V}{(4\pi)^{3/2}}$ ✓, then $B\cdot\tfrac{\sqrt\pi}{2}s^{-1/2} = \tfrac{s}{(4\pi s)^{3/2}}\cdot\tfrac16\int\sqrt gR$ gives $B = \tfrac{1}{24\pi^2}\int\sqrt g R$.)*

Now evaluate the covariantly regulated vacuum energy of Ch. 25,

$$E(s) \;=\; \tfrac12\sum_n\omega_n e^{-s\omega_n} \;=\; \tfrac12\int_0^\infty \omega\,e^{-s\omega}\,dN(\omega),$$

term by term in (B.6). The volume term gives $\tfrac12 \cdot \tfrac{V}{2\pi^2}\int_0^\infty\omega^3 e^{-s\omega}d\omega = \tfrac{3V}{2\pi^2 s^4}$ — the bulk $V/s^4$ piece (cosmological term; standing sentence applies). The curvature term gives

$$E(s) \;\supset\; \frac{1}{2}\cdot\frac{1}{24\pi^2}\Big(\int_0^\infty \omega\,e^{-s\omega}\,d\omega\Big)\int\sqrt g\,R^{(3)} \;=\; \frac{1}{48\pi^2\,s^2}\int\sqrt g\,R^{(3)},$$

i.e. the induced spatial Einstein–Hilbert action with the **universal coefficient**

$$\boxed{\;C(s) \;=\; \frac{1}{48\pi^2\,s^2}\;} \tag{B.7}$$

— no model parameters, no fit freedom: one massless field, one regulator scale, one number. Chapter 25's measured TT stiffnesses land on $-C(s)$ within $1$–$2\%$, which is the only sentence of physics this appendix needs to motivate its existence.

---

**Validation.** §B.1–B.2: `ch09_eta_toy.py` (three-regulator agreement; plateau-and-jump behaviour). §B.3–B.4's algebra is closed-form; its consumers validate it wholesale — `ch22_casimir_chain.py` / `ch22_dispersion.py` (the $-\tfrac1{12}$ and the induced stiffness) and `ch25_graviton.py` (the $C(s)$ convergence table).
