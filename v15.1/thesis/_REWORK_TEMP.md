# v15 thesis — self-consistency review: portions needing rework

*Temporary working file, generated 2026-08-21. Five parallel adversarial review passes (Part I; ch08–12; ch13–15+App.A; ch16–21; ch22–27+App.B/C), each re-deriving key equations by hand and re-running every recomputable number. Delete this file once the items are cleared.*

**Status update (2026-08-21, fix round 1): C1 and C2 are RESOLVED** — see the annotations in place and the "Resolved this round" section at the end.

**Status update (2026-08-22, fix round 5): the MINORS are RESOLVED** (two accepted-as-is items, stated in the minors section). **The review is closed** — every finding in this file is fixed or explicitly accepted. This file can now be deleted, or kept as the audit trail of the round.

**Status update (2026-08-22, fix round 4): T1–T6 are RESOLVED** (plus minors F4 and m1 en route) — the Kasner theorem carries its sign vector and the honest reachable set, plateau exactness rests on the locality/GW argument instead of the false rescaling one, ch21 bills five imports including the fluid, ch22's hybrid is a declared postulate with its self-consistency check, Fig. 7.1's caption is group-theoretically possible, and ch25's error model is stated before its central values. **All CRITICAL and MAJOR items are now closed**; only minors remain.

**Status update (2026-08-22, fix round 3): B1–B4 are RESOLVED** (plus the ch14 minors F5–F8 and F15) — the $1/N_c$ factor is displayed, the budget is one assumption-labeled sixteen-order ladder with observation landing between the two GIM-evaded rungs, the spin null is stated as live, and Fig. 14.3 has no tuned bar. Remaining: T1–T6 and the non-ch14 minors.

**Status update (2026-08-22, fix round 2): S1–S4 are RESOLVED** — see the annotations on the sign/convention section. The headline: the ch22/ch25 contradiction dissolved *in the model's favor* — the quadratically divergent EH coefficient is regulator-scheme-dependent (sign included), and evaluating the same ch25 lattice spectra in the 4D-covariant proper-time scheme (now a second [Computed] table in `ch25_graviton.py` and §25.5) gives κ_TT > 0 (stable waves), κ_conf < 0 (attraction), diffeo at machine zero, and G_ind > 0 — the absolute signs GR requires, consistent with ch22's regulator-independent 1D anchor and with ch22's top-down Sakharov (which already used this scheme). Still open: the ch14 items (B1–B4), the theorem repairs (T1–T6), and the minors.

**Bottom line.** The mathematics is in unusually good shape: essentially every closed-form derivation re-derived by hand checked out, and every quoted number that could be recomputed reproduced to all printed digits (ch11 convergence table, crossing size to machine precision, contraction norms, Kasner/BKL/Sakharov arithmetic, ch13/14/15 chains, App. A/B algebra). All five previously catalogued bugs from the papers-critique round (resonance N_con=1, CP-Coulomb, spherical ansatz, N_spin=2, C₂^eff bound) are **dispatched, not inherited** — with the residues noted below. What needs rework falls into four clusters: **(1)** the Part III reproducibility apparatus is currently a fiction, including the make-or-break graviton table; **(2)** an unaudited sign/convention chain runs through ch22→ch24→ch25 (DeWitt signature) and a second one through ch01/ch09/ch12 (Δη and Goldstone–Wilczek); **(3)** ch14's η_B assembly has a silent factor 3, an incoherent CP-budget floor, and two residual over-claims; **(4)** two [Theorem]-labeled statements (Kasner Selection, plateau exactness) are supported by arguments that don't prove them as stated.

---

## CRITICAL — must fix before any external reading

### ~~C1~~ RESOLVED — ch25 graviton table now reproduced exactly by `scripts/ch25_graviton.py`

A new `ch25_graviton.py` implements the method App. C.4 describes (63³ lattice, 3969 blocks of 63, s = 2.0/2.5/3.0, k = 1..5, quartic fits, flat-space + uniform-shear unit tests with a hard abort, diffeo control incl. the raw s = 0 sum, one diagonalization per configuration). **The run reproduces the chapter's table to every printed digit** (−0.6895/−0.7786/+14.23/−3.46; −0.8932/−0.9402/+7.72/−1.341; −0.9866/−1.0077/+3.764/−0.4527) — the quoted table was from a genuine run whose script had never been committed. It writes `data/ch25_graviton_results.json` and regenerates Figs. 25.2–25.3 from live data. Chapter text updated where it disagreed with the script's own prints: uniform-shear unit test 3×10⁻⁸ → 5×10⁻⁹; raw diffeo Π 0.42 → 0.39 (three places); footer rewritten; App. C.4's "fit covariances propagated" replaced by the honest rms-fit-residual statement (a 5-point deterministic fit has no statistical covariance). T6 (percent-headline vs 45% systematic) is **still open** — it is a text/error-model question, not a provenance one.

### C1 (original finding, for the record). ch25 graviton table was not reproducible by any code in the repo
`ch25_graviton.md` ~59–67, 97; `scripts/ports/graviton.py`; `appC_numerics.md` ~30.
The chapter quotes a three-row convergence table at s = 2.0/2.5/3.0, 63³ lattice, k = 1..5, **quartic** fits, a uniform-shear unit test at 3×10⁻⁸, output `data/ch25_graviton_results.json`. The only graviton code in the repo (`ports/graviton.py`) runs s ∈ {2.0, 1.5}, 3D on **45³**, k = 1..4, **linear** fits, no shear unit test, never computes C(s) or κ/C, and writes into `ports/`. The s=2.5 and s=3.0 rows — including the headline −0.987/−1.008 — cannot be printed by extant code. By the thesis's own "one number, one script" contract the [Computed] labels on the flagship Part III result are unearned.
**Fix:** commit the script that actually produced the table, or rerun with the existing code and re-quote (and rewrite App. C.4's method description to match reality).

### ~~C2~~ RESOLVED — all cited scripts now exist and run; `--report` implemented; `data/` populated

Created as adaptations of the `ports/` originals (which stay as archive): `ch19_fidelity.py`, `ch20_kasner_sim.py`, `ch21_friedmann.py`, `ch21_collapse.py`, `ch22_casimir_chain.py`, `ch22_dispersion.py`, `ch24_newton_lattice.py`, `ch24_mixmaster.py`, plus the shared `ladder.py` module App. C.3 promised. All output to the contract paths (`data/*.npz|json`, `thesis/figures/chNN_*.png`), all run clean, and **every chapter headline number reproduces**: Kasner exponents/Σ_K 0.9997/0.3333, mixmaster 7.2e-18 → 8.9e-13 drift with the 1e-10 discard gate now a hard `sys.exit`, u-map to 6.2e-4, Newton A = 0.0766 / slope −1.008 / force −2.05, Friedmann 0.0090/0.0091, collapse 1.082/1168/1863, chain 0.135/3.865/24.000 with programmatic asserts, dispersion median −0.13 (k ≤ 8). `run_all.py --report` writes `data/provenance.json` (script → validates → status → captured printed numbers); missing scripts now count as failures instead of a silent "[skip]". Validation footers of ch19/ch20/ch21/ch22/ch24/ch25 updated (no more "pending wrapper" debts); App. B's script citations became valid automatically; `make_schematics.py` no longer hand-draws ch25_fig2 (it is a data figure owned by `ch25_graviton.py`).

### C2 (original finding, for the record). Part III reproducibility contract violated wholesale
`appC_numerics.md` ~3–10, 26, 34, 38; `appD_notation.md` D.7; `run_all.py` 46–65; `data/` (empty).
Chapters cite `ch19_fidelity.py`, `ch20_kasner_sim.py`, `ch21_friedmann.py`, `ch21_collapse.py`, `ch22_casimir_chain.py`, `ch22_dispersion.py`, `ch24_newton_lattice.py`, `ch24_mixmaster.py`, `ch25_graviton.py` — **none exist under those names**. The code lives in `scripts/ports/` under different names (`newton_sim.py`, `stiffness_dispersion.py`, `newton_lattice.py`, `mixmaster.py`, `graviton.py`, `kasner_sim.py`, `dewitt_fidelity.py`), which `run_all.py` never runs (it prints "[skip] (not yet written)"). `data/` is empty: no `provenance.json` (promised by App. C.6, App. D.7, 00_master), no result JSONs; `run_all.py` has no `--report` flag; App. C.3's shared module `ladder` does not exist; ch26's "run_all.py regenerates the complete numerical substrate in one command" is false today. `appB_heat_kernel.md` ~81 cites the nonexistent names too. Also App. C's "any run drifting beyond 10⁻¹⁰ is discarded by the script itself" — `mixmaster.py` only prints the drift, no gate.
**Fix:** port/rename the ports scripts to the cited chapter names, wire them into `run_all.py`, implement `--report`, populate `data/`, add the mixmaster drift assert — or rewrite App. C to describe what actually exists. (Note: the ports scripts' *contents* were verified to match the quoted setups where they exist, so this is bookkeeping/porting work, not missing science — except C1.)

---

## MAJOR — sign/convention chains (fix as one audit each)

### ~~S1~~ RESOLVED (2026-08-22) — two-scheme sign audit, computed and written into §25.5

The full derivation chain was audited from scratch. Invariant anchor (convention-proof, from the closed-FRW Friedmann curvature term): GR's static potential energy is E_pot = −(1/16πG)∫√g R³ — i.e. TT gradients cost energy, conformal gradients release it. The measured E_ind = +C(s)∫√gR³ has the opposite overall sign, and the naive matching to the covariant action crosses two hidden minuses (Γ_static = −ET, and √−g R⁴|_static = −√g R³ in (+,−,−,−)) which cancel — so the old §25.5 chain actually implied G_ind < 0. Resolution, verified numerically on the same lattice spectra: the quadratically divergent EH coefficient is regulator-scheme-dependent, sign included. In the 4D-covariant proper-time scheme (exact T=0 mode-sum reduction E_4D(s₀) = ½Σω erfc(√s₀ω) − Σe^{−s₀ω²}/(2√(πs₀)); Weyl coefficient C_4D = −1/(192π²s₀)) the same spectra give κ_TT/C_4D = −0.93/−0.92 (s₀=4; polarizations agree to 0.4%), κ_conf opposite, and the diffeo control at ~10⁻¹⁶ (machine zero — the 4D functional is covariant on the lattice to float precision). That is: κ_TT > 0, κ_conf < 0, G_ind = 1/(16π|C_4D|) > 0 — gravity's absolute signs. Implemented: second [Computed] table in `ch25_graviton.py` (+JSON `table_4dpt`); §25.5 rewritten as the sign audit (trap stated openly, both tables interpreted); §25.0 and §25.6 synced; ch22 §22.2 now flags the 1D number as the regulator-independent anchor (finite Casimir-type difference, no scheme ambiguity) and §22.5 links the top-down Sakharov (−Λ²/192π² per scalar — already the 4D scheme) to the audit; ch24 §24.1 reframed from energy-minimization (unbounded for the conformal mode) to constraint-elimination, with κ = |κ_grad| declared; `ports/graviton.py` docstring carries an archive correction note. M2 thereby resolves without demotion — "attraction derived" now has precise support (1D regulator-independent + covariant-scheme 3D) — and M3 is the ch24 reframe.

### S1 (original finding, for the record). The DeWitt-signature sign chain: ch22 vs ch24 vs ch25 contradict each other
- `ch22_stiffness.md` ~39–47 measures the conformal vacuum stiffness **negative** (κ ≈ −0.13) and reads *negative conformal = attraction = DeWitt*.
- `ch25_graviton.md` ~43–55, 79–85 predicts/measures κ_conf = **+C**, κ_TT = **−C**, and reads *that* as DeWitt. Under ch22's dictionary, ch25's sign would say compression is stable (no Jeans attraction). §25.5 silently conflates static **energy** with static **Lagrangian** (E = −L), dropping a minus in the E↔Γ↔H chain. Symptom: `ports/graviton.py`'s docstring wants κ_TT > 0, the opposite of eq. (25.6).
- `ch24_newton_bkl.md` ~9–17: the Poisson derivation uses E = Σ(κ/2)(δα_i−δα_j)² with κ > 0 for the conformal mode — the very mode ch22 measured negative; the chapter never says whether vacuum, matter, or net stiffness enters (24.1), and with the vacuum sign the functional is unbounded.
**Fix:** one dedicated sign-conventions box tracking energy → effective action → ADM Hamiltonian through ch22 §22.2, (24.1), and (25.6)/§25.5, then reconcile all three chapters and the graviton script docstring. Until done, demote ch22 §22.2/§22.5 and ch23 §23.4(3)'s "gravity is attractive — derived/measured" to "sign measured in 1D; 3+1D consistency pending" (the sign evidence is entirely 1D, and the Sakharov route inherits its sign from the S_EH convention — a matching, not a derivation).

### ~~S2~~ RESOLVED (2026-08-22) — convention pinned in Cor. 9.2: Δη ≡ η_after − η_before, ΔQ_net = +½Δη
Fixed at ch09 Cor. 9.2 (definition added), ch09 §9.6, the ch01 results table, ch11 Rule 1, and the ch26 chain table (the last two also carried the stale −½Δη). The pump's ΔQ = Q_vac(0.30) − Q_vac(0.70) = +1 line was checked against the new convention and is consistent.

### S2 (original finding, for the record). The Δη sign chain: master-theorem abbreviation contradicts its own corollary
`ch01_introduction.md` ~110; `ch09_master_theorem.md` ~121 vs ~45, 89, 93.
Boxed Thm 9.2 is ΔQ = −½[η_before − η_after], i.e. **+½Δη** under the standard Δ = after−before, which is exactly the convention Cor. 9.2 (9.3) uses (up-crossing → Δη = +2 → ΔQ = +1, confirmed by ch12's pump numerics). But the §9.6 summary and the ch01 results table print "ΔQ = **−**½Δη". Same symbol, opposite conventions, three paragraphs apart — and this program has been burned by exactly this sign typo before (v13 §46.2).
**Fix:** define Δη ≡ η_after − η_before once in ch09, write +½Δη in §9.6 and the ch01 table (or spell the bracket out in both places).

### ~~S3~~ RESOLVED (2026-08-22) — toolbox rewritten on the bulk angle δ(x) with the (10.3) map θ = −δ made explicit
The GW current now reads ⟨j^μ⟩ = ε^{μν}∂_ν δ(x)/2π with the bulk profile running −θ₀ → −θ_L, giving Q = −Δ/2π including the sign; the "sign conventions matched by ordering" hand-wave is gone. Also fixed: §12.2(b) "Σ = δ" → "Σ = −δ" (with the note that only cos Σ enters, so the sign is spectrally invisible there but load-bearing for the toolbox), and the Fractional Charge Law's per-wall clause now carries wall orientation (+θ₀/2π left, −θ_L/2π right — summing to −Δ/2π, not −Σ/π).

### S3 (original finding, for the record). Goldstone–Wilczek toolbox sign discrepancy, papered over
`ch12_pumping.md` ~85–93. The toolbox computes Q = ∫θ′/2π dx = **+Δ/2π** and calls it "precisely the law (12.8)", which is **−Δ/2π**; the mismatch is waved off as "sign conventions matched by the Q_vac = −η/2 ordering". The real resolution is already in the thesis: by (10.3) a bulk mass angle δ(x) corresponds to walls θ = −δ, so the bulk profile matching the wall configuration is δ(x) = −θ(x), giving −Δ/2π (independently confirmed numerically: bulk δ=+0.7 ↔ walls (0.7,0.7) give opposite-sign truncated charges). One-line fix, but as printed it's a sign error sitting in the flagship "overdetermination" paragraph. Related in the same chapter: ~67 "a bulk phase, which is the case Σ = δ" — by (10.3) it's Σ = −δ; and ~81 "each wall binds −θ/2π" summed gives −Σ/π, not −Δ/2π — the per-wall charge must carry wall orientation.

### ~~S4~~ RESOLVED (2026-08-22) — the two conventions now flagged at every crossing point
ch04 §4.2 already flagged the q → π − q relabeling; added the missing re-flag in §4.3's Bessel toolbox ("raw +J convention, the (−i)^d belongs to this one"), rewrote §4.2's "(its θ is this q)" gloss to name ch20's convention and the band-center coincidence of the two, and flagged (20.1) as the raw-+t convention with the pointer back to (4.3). No physics changed — the clash was purely notational, as suspected.

### S4 (original finding, for the record). Tight-binding dispersion sign clash ch04 vs ch20 (found independently by two reviewers)
`ch04_ladder_dynamics.md` (4.3) states E(q) = −2J cos q (and §4.2 adopts it, Fig. 4.1 agrees); ch20 (20.1) states ω = +2t cos θ for the same Hamiltonian; for +J hopping +2J cos q is correct. Also §4.2 vs §4.3 flip within ch04 itself with no re-flag; the Bessel propagator (−i)^d J_d is right only for the + form. Physics downstream unaffected (only swaps which fixed point expands), but two chapters currently contradict each other about the same operator.
**Fix:** pick one convention in ch04, flag it, and align ch20 §20.1's citation.

---

## MAJOR — ch14 η_B assembly

### ~~B1–B4~~ RESOLVED (2026-08-22, fix round 3)

**B1:** the $1/N_c$ quark→baryon factor is now a displayed factor in (14.1) (with a candid parenthetical about the earlier omission), in the §14.4 assembly ($N_{\text{dof}}/N_c = 1$), and in `ch14_eta_assembly.py`. New numbers, printed by the script: prefactor $5.33\times10^{-2}$ (was $1.6\times10^{-1}$); (14.9) required $\mathcal P_\times\varepsilon_{CP} = 1.15\times10^{-8}$ (was $3.8\times10^{-9}$); required $\varepsilon_{CP} = 1.2\times10^{-8}$–$1.2\times10^{-7}$ over $\mathcal P_\times = 1$–$0.1$.

**B2:** the four incompatible floors are replaced by one **assumption-labeled reference ladder** built from ch13.5's budget lines and printed by the script: bulk-perturbative $J\,\mathcal S_{\text{GIM}}(1/16\pi^2)^2 = 1.2\times10^{-16}$ | GIM-evaded two-loop $J(1/16\pi^2)^2 = 1.2\times10^{-9}$ | GIM-evaded one-loop ceiling $J/16\pi^2 = 1.9\times10^{-7}$ | maximal $\mathcal O(1)$ — a sixteen-order range replacing the underived $[10^{-10}, 1]$ "ten-order" bracket. "Nineteen orders" deleted; the honest placements (factor ~6 below the ceiling, ~25 above the two-loop rung, eight orders above bulk-perturbative at $\mathcal P_\times = 0.38$) are printed by the script. Loop-counting inconsistency (one vs two loops) resolved by making both rungs explicit. Synced: ch27 item 4 (stale $10^{-10}$ → the ladder), ch26's executioner ("nine" → "eight orders").

**B3:** §14.5 item 3 rewritten: spin factor is $[0, 2]$, the null is a live outcome that would set $\eta_B \equiv 0$ (ch27 item 3 cited), the factor explicitly *not* absorbed into the budget bars; "each bounded" → "the first two bounded, the third **not**"; §14.6 now carries the caveat too (per-channel unit is 1+1D; the channel sum is the open computation).

**B4:** the vacuous "consistency test it could easily have failed" is replaced by its inversion — window-membership carries no evidential weight (nothing below $5\times10^{-2}$ could fall outside), which is *why* $\varepsilon_{CP}$ is [Open]; the informative statement is that observation lands **between the two GIM-evaded rungs**. Fig. 14.3 regenerated with assumption-labeled rungs only (the tuned $\varepsilon = 10^{-8}$ "partial evasion" bar that landed 1% from observation is gone); caption updated to say none is tuned.

Also fixed in the same pass (ch14 minors from the review): F5 (footer "2000" → 400/390/101, matching the script); F6 (the $Q_{\text{vac}}$-jump [Computed] now attributed as slope-correlation + Cor. 9.2 bookkeeping, spectrally verified at the ch12 demo point); F7 ("exactly zero" properly scoped: exact for support inside $|\Delta| \le \pi$, doubly-exponentially small $\sim e^{-\pi^2/4\varsigma^2}$ for Gaussian tails, MC-resolution zero in the table); F8 ("quarter turn" → half turn of the $2\pi$-periodic angle); F15 (`ch14_eta_assembly.py` comment `2 + (7/8)*6*(4/11)` corrected).

### B1 (original finding, for the record). Silent factor 3: quark number equated with baryon number
`ch14_eta_b.md` ~13, 89–91. N_dof = N_c = 3 multiplies the per-bag pumped charge, but η_B counts **baryons** and n_B = n_quark/3; no 1/3 appears anywhere in (14.1)/(14.8). Honest prefactor: 7.04 × 0.02136 × (28/79) × 1 = 5.3×10⁻², so (14.9) becomes P×ε ≈ 1.15×10⁻⁸, a factor 3 off the quoted 3.84×10⁻⁹. In a chain carrying 28/79 to four digits, an unexamined factor 3 is inconsistent rigor.
**Fix:** insert the explicit 1/3 (or define pumped charge as baryon number per bag and drop N_dof=3); update (14.8), (14.9), the ε_CP ranges, `ch14_eta_assembly.py`, Fig. 14.3.

### B2. CP-budget floor stated four mutually incompatible ways; "nineteen orders" matches none
`ch14_eta_b.md` ~85–87, 97; `ch13_baryogenesis_interlude.md` ~90–92; ch27 item 4. Candidate floors in the text: O(10⁻¹⁰) (line 85, ch27, script scenarios) → gap 1.6–2.6 orders; ch13's inherited J·S_GIM·(1/16π²)² = 1.2×10⁻¹⁶ → 7.5–8.5 orders; "10⁻¹⁹-class" → 10.6–11.6 orders; nothing yields the printed "nineteen orders". Also: line 97 prices the electroweak scale as J/16π² (one loop) while ch13 line 92 hands over two loop factors; S_GIM as defined keeps only 4 of the Jarlskog determinant's 12 mass powers (ch13 is internally aware, ch14's composition isn't); and at the script's own "GIM-crushed" ε=10⁻¹⁰, η_B = 6.1×10⁻¹² — only two orders below observation, quietly contradicting the "GIM-crushed means dead" narrative.
**Fix:** derive one floor (decide which of J, S_GIM, and how many loop factors compose into ε_CP and why), propagate consistently through ch13 §13.5, ch14 lines 85/97, Fig. 14.3 labels, ch27 item 4; delete "nineteen orders".

### B3. Spin-cancellation risk understated: bracket [½,1] excludes the live null outcome
`ch14_eta_b.md` ~105, 109. The fatal N_spin=2 bug is gone (N_dof=3 is colors only — good), but item 3 asserts an unresolved spin factor "in [½,1]" with no derivation and line 105 calls all three uncertainty items "bounded". The falsified paper-3 result is that at k_⊥=0 the two spin polarisations pump exactly opposite dQ; until the per-channel 3+1D computation is done, exact cancellation (factor 0 ⇒ η_B ≡ 0) is a live outcome — ch27 item 3 says so itself ("a null across all channels would kill the baryogenesis uplift"), and ch10.A's twofold slab degeneracy is precisely the fact that made the paper-3 assembly go wrong.
**Fix:** bracket → [0,1]×spin, delete "each bounded", state explicitly that a null is possible and would zero the estimate, cite ch27 item 3.

### B4. Two residual over-claims against the "honest budget" promise
`ch14_eta_b.md` ~97 and Fig. 14.3. (i) "A consistency test it could easily have failed": the budget window spans 0 to 1.6×10⁻¹, so **no** observed η_B below 0.16 could fail it — window-membership is vacuous, not a passed test. Invert the sentence (the window is so wide observation cannot miss it; that is why ε_CP is [Open]). (ii) The "partial evasion" bar (P=0.38, ε=10⁻⁸) lands at 6.07×10⁻¹⁰ — 1% from observation — presented as a generic scenario; label it as the *required* value or choose scenario numbers not derived from the answer.

---

## MAJOR — theorem statements weaker than their labels

### ~~T1–T6~~ RESOLVED (2026-08-22, fix round 4)

**T1:** Thm 20.1(iii) restated with the sign vector $\sigma_a \in \{\pm1\}$ ($\theta_a = -\sigma_a\pi/2$), $p_a = \sigma_a g_a/\sum\sigma_b g_b$; the reachable set is stated as the *full plane* $\sum p_a = 1$ (generalized-Kasner family) with the all-positive case strictly inside the circle, the circle requiring mixed signs, $\Sigma_K > 1$ reachable (GR's negative-energy stiff matter), and circle-landing a codimension-one condition; the proof shows surjectivity ($\sigma_a = \operatorname{sgn}p_a$, $g_a \propto |p_a|$). §20.3's dichotomy extended to on/inside/**outside**, with excluding $\Sigma_K > 1$ named as part of ch21's constraint's job, and the $u=2$ run flagged as tuned-to-the-circle. §20.4 launch description now names $\sigma_a$.

**T2:** the false "pure rescaling / relabeling" justification for plateau exactness is replaced in Cor. 9.2, its explanatory parenthetical, and App. B.2 by the correct *locality* argument: the APS-local smooth part is an integral of a CP-odd local density, sourced only by bulk chiral gradients (absent for constant real mass) and wall data (held fixed) — equivalently the GW pinning of the fractional part at $-\Delta/2\pi$ at every $L$; the clause is labeled **[argument + Computed]** with the measured plateau flatness. Both texts now state explicitly that at $m \ne 0$ a size change is *not* a rescaling and the sub-gap level genuinely migrates.

**T3:** ch21 §21.1 now bills **five** ingredients: the old "form" bullet split into *sign measured* (Ch. 22 anchor, Ch. 25.5 audit) vs *exact balance* [Postulate], and the perfect fluid itself (dilution law, free $w$) added as an explicit [Postulate] with the note that the model owns a derivable equation of state never computed here (Ch. 27 item 1). "The graft is small" → "finite and itemized". ch26 row 14 label gains "+ fluid Postulates".

**T4:** ch22's cell energetics now declares the hybrid as a [Postulate] — massive nonrelativistic occupant ($E = \pi^2n^2/2mL^2$, $m=1$ units restored) + massless species' Casimir vacuum, two different fields — with the two failure modes of a single-field version stated (scale-free $1/L$ / exponentially suppressed Casimir) and the self-consistency check $n\pi/L^* = 1/24n \ll m$ printed in text.

**T5:** Fig. 7.1 caption rewritten: CPT exact *within sectors*, ECPT conjecturally exact on the full dynamics, and full-dynamics CPT broken *by group multiplication* $E = (ECPT)(CPT)$ — the impossible "both exact with E broken" claim is gone.

**T6:** ch25 §25.4 now states the error model first: $\kappa_{TT} = -C(s)$ within the $\pm45\%$ diffeo-control systematic; central values $1.3\%/0.8\%$, polarization agreement $2\%$, with the overestimate argument (diffeo has $R^{(1)} \ne 0$ structure TT lacks) given but "percent-level" no longer claimed as an error bar. Scorecard, §25.6, ch01 table, and the abstract synced ("central values … within stated systematics"). Also fixed en route: ch01/ch26 growth-collapse rows (F4: growth exponent $1.08$ vs thresholds $1.062/1.686$), ch26 row 15 induced-gravity-Postulate label (m1).

### T1 (original finding, for the record). Kasner Selection Theorem (ch20) cannot reach its own headline as stated
`ch20_kasner.md` ~41–57. Part (iii) states p_a = g_a/Σg_b with positive couplings — which puts every reachable point strictly *inside* the Kasner circle (all p_a > 0, Σp² < 1); the vacuum circle requires a negative exponent. The actual mechanism (used correctly in `kasner_sim.py`: `th = -pi/2*sign(p)`) is p_a = σ_a g_a/Σσ_b g_b with per-axis signs σ_a = ±1, which appears only in prose. With mixed signs the reachable set also includes Σ_K > 1 (e.g. g=(0.4,0.3,0.3), σ=(−1,+1,+1) → Σp² = 8.5) — negative-energy stiff matter in GR, not a Kasner/Jacobs orbit — so §20.3's dichotomy "on the circle or inside it" is incomplete and "precisely the Kasner/Jacobs orbit class" oversells: the theorem selects straight lines in α-space, a strictly larger family curtailed only by ch21's postulated constraint.
**Fix:** restate (iii) with σ_a; give the reachable set as all Σ_K ≥ 1/3 *and* Σ_K > 1; state that hitting the circle exactly is a codimension-one tuning of the g_a, and that excluding Σ_K > 1 is part of the ch21 constraint's job.

### T2. Plateau-exactness justification in Cor. 9.2 / App. B.2 is wrong for m ≠ 0
`ch09_master_theorem.md` ~95–99; `appB_heat_kernel.md` ~37. Exact quantization is justified by "pure rescalings move levels without tilting branch densities" — but changing L at fixed m and wall data is **not** a pure rescaling (mL is a shape parameter; the sub-gap level migrates across the whole gap, verified numerically: −0.46m → +0.36m between L=0.3 and 0.7). The rescaling argument covers only m=0. The conclusion is right, but its printed support inside a [Theorem]-labeled corollary is (i) numerics and (ii) the anomaly argument, neither of which is the argument printed.
**Fix:** prove dη/dL = 0 between crossings at fixed boundary data (APS variation route), or demote the plateau-exactness clause to [Computed] with the anomaly argument as theoretical support.

### T3. ch21's bill of imports omits the fluid itself
`ch21_friedmann.md` ~9–26. The four-ingredient [Postulate] accounting (constraint form, dynamics, κ, curvature offset) omits that the matter side is an externally supplied continuum perfect fluid with free w — nothing connects ρ ∝ e^{−3(1+w)Δα} to the model's actual matter (particle in a box, E ∝ 1/L², a specific derivable equation of state never computed). Also split the (21.1) bullet: "sign: measured (Ch. 22)" vs "exact balance: [Postulate]".
**Fix:** add the dilution law/e.o.s. to the imported list, or derive w for the model's own matter and use it.

### T4. ch22's cell energetics mixes incompatible matter treatments, and stability depends on the mix
`ch22_stiffness.md` ~19–25, 41–45. The Casimir term −π/24L is massless; the occupant term π²n²/2L² is nonrelativistic massive with m=1 unstated (as displayed the equation is dimensionally inconsistent, 1/L vs 1/L²). Massless occupant (E = nπ/L) → total ∝ 1/L, monotone, **no stable L\***, no occupied/empty dichotomy; massive field → Casimir not −π/24L. Same hybrid underlies κ_matter = 3π²n² in (22.4).
**Fix:** declare the hybrid (Part I nonrelativistic occupants + massless-limit vacuum) as an explicit modeling choice with regime of validity; write π²n²/2mL² with m=1.

### T5. Fig. 7.1 caption asserts a group-theoretic impossibility
`ch07_ecpt.md` ~47–49 vs 65–67. Caption: CPT exact and ECPT exact with E broken — but E = (ECPT)·(CPT), so both exact ⇒ E exact; the caption also contradicts itself ("all partial products containing E broken" — ECPT contains E) and drops the within-sector qualifier that §7.5 uses correctly.
**Fix:** caption → "CPT exact *within sectors* (fixed geometry); ECPT conjecturally exact on the full dynamics; on the full inter-sector dynamics CPT is then necessarily broken along with E."

### T6. ch25's percent-level headline outruns its own systematic bound
`ch25_graviton.md` ~67, 93. The chapter designates the diffeo control at s=3 as the conservative systematic bound (±45%), then headlines "converged at the percent level / within 1.3% and 0.8%". Central values 30× inside one's own systematic bound: either quote "consistent with −C(s) within ~45% systematics; central values 1–2% off" or justify a tighter error model. (Interacts with C1: re-quote from whatever code is committed.)

---

## ~~MINOR~~ RESOLVED (2026-08-22, fix round 5 — the wrap-up sweep)

Every item below is fixed, with two deliberate exceptions:
- **ch10 marginal points $E = \pm m$**: the Mirror Theorem's proof already handles $E = -m_R$ exactly and remarks on the $\delta = 0$ endpoint; the residual gap is cosmetic and left as is.
- **Fig. 3.2 (telescoping schematic)**: still the one carried-over drawing not regenerated by the suite — ch03 flags it and ch01's headline now carries the exception clause.

Highlights of the sweep: **appD's triad entry corrected** to the Gram form $g_{ab} = \sum_i e^i_a e^i_b$ with the ch16.9/ch18 index-convention swap spelled out; the ghost $\mathcal C_2^{\text{eff}}$ entry purged; D.5 gains the $\theta$-vs-$q$ and $c$-vs-$c_m$ rulings; D.6's Kasner digits corrected. **ch10's Zero-Charge "Full stop" scoped** to $|\delta| < \pi/2$ with the $(\pi/2, \pi]$ gap stated; **Figs. 10.2/10.3 renumbered** (files, scripts, captions, footer) to match display order. **ch12's Wall-Mismatch Criterion relabeled** [Theorem ("if"); Computed ("only if")], Fig. 12.1's caption corrected to few-$10^{-4}$, and the plateau/jump residuals separated ($6\times10^{-4}$ vs $1.5\times10^{-4}$). **ch15's boxed $\pi/8$ made dimensionless** (text, figure caption, validation, script y-label) and the divergence-degree sentence untangled. **ch05's Landau–Zener sentence** rewritten to match Cor. 9.3 (no transmission factor on the unit), and its stray $\hbar$ absorbed. **ch23's dangling formula** repaired and the Proca-sign question flagged as an explicit one-line [Open] audit instead of passing silently. **ch24's $E_{\text{int}}$ sign** fixed and the coupling renamed $c_m$. **ch16's mode equation** carries the conformal-coupling qualifier. **ch18's** $2\times10^{-12}$ harmonized and (18.5) marked as a vector-space (Cartan) decomposition. **ch22** gets the App. D.5 first-use $s$-flag. **Orphan figures**: 16.1 moved to §16.2, 17.1 / 18.3 / 19.1 now displayed with captions (numbering gaps closed). ch09's early $\Delta$ is glossed at first use; ch03's reachable-set clause covers both directions; ch01's $\sigma \approx 1/n_1$ carries its own [Computed] label; ch13's "quark pairs" line corrected to baryons-per-photon; appA's "(A.3 iii)" disambiguated.

Verified after the sweep: `run_all.py --fast --report` — 23 scripts, 0 failures.

**With this round, every item in this file — critical, major, and minor — is either fixed or explicitly accepted with a stated reason. The review is closed.**

## MINOR (original list, for the record) — by file

**ch01:** results-table "growth/collapse at 1.062/1.686" mislabels 1.062 (it is the *turnaround* threshold; growth is the exponent 1.08 — cf. ch21 ~64/71). σ ≈ 1/n₁ carried inside [Theorem] labels at ~31/108 (it's an ensemble-dependent estimate; ch06 keeps the distinction, ch01 compresses it). "Everything is regenerable" (~12) vs ch03 ~186's admitted non-regenerated Fig. 3.2 (and the C2 skips) — add "(with the exceptions the runner itself flags)".

**ch03:** ~88 justifies reachable set (3.2) by upward telescoping only; state downward reachability ($k < n_0$) explicitly.

**ch04:** internal §4.2 vs §4.3 sign flip (see S4). ~29's "(its θ is this q)" identification with ch20 is loose (q ∈ (0,π) here; ch20 uses both fixed points ±π/2).

**ch05:** ~39 keeps ħ explicit against the declared convention (only Ch. 2 is exempt). ~49's Landau–Zener sentence says a finite-speed crossing "rescales the probability that the unit is transferred" — the framework's own Cor. 9.3 and ch27 item 10 forbid exactly that transmission-factor reading (the unit passes with probability 1, or 0 if an interaction gap opens); reword to match Cor. 9.3.

**ch10:** Mirror Theorem scope |δ| < π/2 vs Zero-Charge Theorem's "at every order in δ … Full stop" — for |δ| ∈ (π/2, π] sub-gap solutions exist and the printed pairing proof doesn't cover them; scope the claim or extend the proof. Figure order: Fig. 10.3 appears before Fig. 10.2. Marginal points E = ±m (p=0) unhandled at δ=0 endpoint (cosmetic).

**ch12:** Wall-Mismatch Criterion "iff": only "if" is proven; converse (asymmetry for *every* (Σ, mL) when sin(Δ/2) ≠ 0) is asserted — weaken to "if + generically only if" or prove it. Fig. 12.1 caption "percent level" vs §12.3 text "few-10⁻⁴"; quoted plateaus (+0.3640/−0.6361) carry ~6×10⁻⁴ regulator bias vs exact ±(0.3634/0.6366) — tighten Richardson settings or quote honestly.

**ch09:** ~97 uses Δ and eq. (12.7) before ch12 defines them; add "(Δ ≡ θ_L − θ₀, Ch. 12)".

**ch13:** ~15 "one excess quark per ~1.6 billion quark–antiquark pairs" — 1.6×10⁹ is photons per excess *baryon*; per quark pair differs by ~an order. S_GIM definition (~90) keeps 4 of 12 mass powers — fine internally, but see B2.

**ch14:** validation says 2000 angle pairs, body and script say 400/390 (~119 vs ~35). ~63's [Computed] overstates: the script computes the slope correlation on ~10 surviving crossings, not the Q_vac jump spectrally. ~75/81 "exactly zero" for Gaussian ensembles is a finite-sample artifact (doubly-exponentially small, not zero; exact zero only for compactly supported P). ~37 "quarter turn" → half turn.

**ch15:** boxed 1+1D result π/(8L₀) has a dimensional slip — ⟨φ²⟩_mvn = π/8 is dimensionless (matches the chapter's own L₀^{1−d} scaling at d=1); fix box, validation paragraph, and `ch15_imn.py` Fig. 15.2 y-label. ~29 mixes actual (log) and superficial (linear) divergence degrees in one sentence.

**ch16:** ~121 the FRW mode equation ω² = k² + m²a² is the *conformally coupled* form; minimal coupling adds −a″/a. Add the qualifier (fine for the thesis's fermion, wrong as a generic [Standard] sentence). ~89 "three digits" vs ch20 "three–four" vs App. D "four digits" for the Kasner landing — settle on three (vacuum) / four (FRW).

**ch18:** §18.2 "10⁻¹²" vs validation "2×10⁻¹²". (18.5) gl(3) = so(3) ⊕ Sym(3) under a Lie-algebra heading: Sym(3) is not a subalgebra; add "(as a vector space; Cartan decomposition)".

**ch19/20/21:** chapter-named wrapper scripts self-flagged "pending" (part of C2). ch20 §20.5's "σ_n = 40 at ⟨n⟩ = 4100" not visibly produced by any run in `kasner_sim.py` — ground it. ch20 uses θ where App. D.5 rules q for ladder quasimomentum (and θ collides with Part II wall angles — either symbol needs a ruling).

**ch23:** ~52 literal dangling "…" mid-formula in the weight count. (23.7)'s −½σ₀²W² is the tachyonic Proca sign in (+,−,−,−) as written (inherited from the ghost-sign σ kinetic term) — one line of sign audit needed for the "Higgs-style mass" claim. (M_P² = σ₀²/6 normalization itself verified correct; the CCJ ratio checks.)

**ch24:** ~15 sign typo: E_int = −∫Φ₁ρ₂ with Φ₁ < 0 makes the displayed expression +1/d; drop the leading minus. Symbol c overloaded three ways (speed of light, lattice constant, matter coupling in (24.1)) with a collision *inside* ch24 (~13 vs ~48); rename the coupling or add a D.5 ruling.

**ch22:** s-collision not flagged at first use (App. D.5's own rule; ch25 does it right). ~25 "L = 0.135 ≈ √2·L_floor" is 5% off √2·0.10 — the script's force law explains it; say so. Script's κ median uses first 8 of 12 modes, text says "twelve modes".

**ch26:** row 15 upgrades L̄ ≈ 1.6c ℓ_P to "[Theorem + Computed/measured]", dropping the induced-gravity [Postulate] premise that ch22 §22.3 itself declares essential (row 14 does it honestly).

**appA:** ~65 "(A.3 iii)" reads as eq. (A.3); write "§A.3(iii)".

**appD:** **triad entry D.4 is wrong** — "g_ij = Σ_a e^i_a e^j_a" is EEᵀ, not ch18's Gram matrix g_ab = Σ_i e^i_a e^i_b (EᵀE); also ch16.9 and ch18 use opposite index conventions (frame vs edge) and D.4 landed on the incoherent hybrid — fix the entry and consider unifying the convention. C₂^eff entry (~54) is a dangling ghost of the retracted budget bound — no chapter uses it; purge it (also from 00_master if present). D.7's `run_all.py --report` claim is part of C2. D.5's q-ruling is stale vs ch20's θ (above). Kasner "four digits" (D.6) overstates.

**scripts:** `ch14_eta_assembly.py` line 18 comment `# 2 + 7/8*(4+6)*(4/11)` evaluates to 5.18, not 43/11 (code is right, comment lies). `ports/graviton.py` docstring contradicts (25.6) on the TT sign (part of S1). `ports/mixmaster.py` missing the drift gate App. C claims.

**figures:** `ch17_fig1_dictionary.png` and `ch18_fig3_conformal_obstruction.png` exist but are never referenced (ch17's first displayed figure is 17.2 — numbering gap will look like an error); Fig. 16.1 displayed at end of ch16 though it illustrates §16.2; ch19's Fig. 19.1 named only in the validation block.

---

## Verified clean (no action)

- Resonance bug fully repaired in ch06 (deficit/divergence roles correctly separated; "no (n₁,n₂) with zero deficit" now explicit). ch07 does not echo it.
- CP-Coulomb and spherical-ansatz bugs explicitly dispatched in ch10 §10.5's honesty note; the inertness proof (ch10) is complete at the level claimed — it is spectral/operator-level, so measure/anomaly subtleties genuinely don't enter.
- N_spin=2 bug not inherited in fatal form (but see B3); C₂^eff bound gone from chapters (but see appD ghost).
- η_B honesty largely delivered: no "9% prediction" anywhere; ch13 states the SM ten-order failure plainly; ch26.4 explicitly disclaims a derived η_B; perturbative regime honestly yields zero.
- ch21↔ch22 circularity audit clean (κ enters ch21 only as a unit; constraint labeled [Postulate]; ch22's two routes depend only on Part I matter).
- All Part I / Part II numbers recomputed reproduce exactly; ch11 table, pump/control pair, crossing size, massless tower η = Δ/π, contraction norms, N_crit, participation ratios, ch13 arithmetic, ch15 numerics, App. A algebra, Sakharov/Planck-cell/C(s)/Bianchi-IX/BKL arithmetic — all verified.
- All scripts and figures cited in Parts I–II exist; Part I–II cross-references all resolve.

## Resolved this round (2026-08-21, fix round 1) — beyond C1/C2 themselves

New findings surfaced while grounding the numbers, all fixed in place:

- **ch22 eq. (22.3) band was wrong as printed**: the actual dispersion band is [−0.178, **+0.060**] — the k = 2 mode is *positive*, which the band [−0.18, −0.05] silently excluded even though the prose already conceded "eleven of twelve". Also the "median −0.13" was the k ≤ 8 median (−0.131; all twelve give −0.142) while the text said "twelve modes". (22.3) and the footer now state both explicitly. *(Supersedes minor m7's second half.)*
- **App. C.3's update-interval claim was false**: "halving the update interval moves results at <10⁻³" — the check had never been run; `ch21_friedmann.py` now performs it and measures **4.8×10⁻³** (still half the quoted Friedmann deviation, so the qualitative claim survives). Text updated to the measured value.
- **ch20 §20.5's σ_n = 40 at ⟨n⟩ = 4100 is grounded**: it is the flat-hopping embedding run (n₀ = 2500 → n(T) = 4100, spread 40), printed by `ch20_kasner_sim.py`. *(Closes the open item in the ch19/20/21 minor.)*
- **ch19's machine-zero value was over-precise**: the quoted χ = 1.8×10⁻¹³ is build-dependent floating-point noise (this build prints −5.3×10⁻¹³ for the dedicated check, ~10⁻¹¹ for the Richardson-combined diagonal rows). Text now claims the honest order of magnitude.
- **ch25 Figs. 25.2/25.3 had no owning script** (the cited `ch25_make_figs` job did not exist; fig2 was hand-drawn in `make_schematics.py` from hard-coded table numbers). Both are now generated from live data by `ch25_graviton.py`.
- **ch21 collapse thresholds, extrapolated**: the new script also prints δ_lin at turnaround/collapse = 1.170/1.867 vs the spherical 1.062/1.686 (~10% high, consistent with the C₊-based prediction discussion in §21.4); quoted in the footer.
- Minor m8 (mixmaster drift gate) fixed — hard discard at 10⁻¹⁰. Minor m9 (App. B citing nonexistent script names) resolved automatically — the names exist now. Registry wording for ch21 updated (0.0091 is "~0.009", not "≤ 0.009").

Still open from the minors list: everything else (ch01/ch03/ch04/ch05/ch09/ch10/ch12/ch13/ch14/ch15/ch16/ch18/ch23/ch24/ch26 text minors, appD corrections, figure orphans, `ch14_eta_assembly.py` comment, `ports/graviton.py` docstring — the last is archive-only now, note it contradicts (25.6) if anyone reads it).

## Suggested fix order (remaining)

1. **S1–S4** (the four sign/convention audits — S1 is the deep one; S2/S3/S4 are one-liners once decided).
2. **B1–B4** (ch14 assembly; B1 changes printed numbers downstream, do before regenerating Fig. 14.3).
3. **T1–T6** (theorem-statement repairs; text-only — T6 now has the rms residuals from `ch25_graviton.py` to quote).
4. Minors, appD corrections, figure orphans.
