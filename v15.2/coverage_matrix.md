# v15 internal provenance and status (NOT part of the thesis)

Maps every v15 asset to its legacy source and records current build status.
The result-level coverage matrix (legacy claim → thesis chapter) lives in
`thesis_skeleton.md` §5 and is authoritative for content; this file tracks files.

## Chapter status

| File | Author | Status |
|---|---|---|
| **all chapters ch01–ch27, App A–D, 00_master, skeleton** | **Fable 5** | **draft v2 — adversarial verification round applied (see below)** |

(The original Fable/Opus split was dissolved by user decision; Fable 5 wrote the
former Opus slots too. The complexity ratings and section specs remain in the
skeleton as documentation of the design.)

## Verification round (draft v1 → v2)

Four independent adversarial referee passes (Part I, Part II, Part III ch16–21,
ch09+synthesis+appendices) re-derived every named theorem by hand, re-ran every
script, wrote independent numerical checks, and then applied their own fixes.
What changed, by chapter:

- **ch06**: false O(t⁴) error term → o(t²) with the measured ~t^0.6 approach;
  Expansion-Preference event definition made explicit; fermionic factorization
  corrected to the Hadamard inequality det G ≤ Π N_con (now printed by
  `ch06_norms.py`); golden-rule scope (short-time, unrenormalized embedding
  coupling) stated in the box and in ch01.
- **ch09**: (9.3) crossing-count sign fixed (up − down); Master Theorem restated
  for fixed domains + domain-growth ledger −½[η(L)+η_fresh−η(L′)]; the false
  "η is an integer staircase" corollary replaced by the APS decomposition
  dη = smooth local + 2·spectral-flow, with the m=0 chiral-tower counterexample
  (η = Δ/π) and fixed-angle plateau flatness added to `ch09_eta_toy.py`.
  The ch12 pump (L-only deformation) is explicitly unaffected.
- **ch10**: intro corrected to one proof (1+1D Mirror) + one diagnosis (the
  3+1D radial reduction is not a consistent projection) + one open problem;
  chiral-rotation sign −δ fixed; Slab Theorem now actually tested
  (`ch10_slab_check.py`, new: 2e-7 at random transverse momenta).
- **ch11**: `ch11_polarization.py` written (was cited but nonexistent);
  `ch11_wall_density.py` rewritten on analytic roots — regional integral now
  converges to −0.0672 (was −0.058 and drifting); ΔQ_same definition written out.
- **ch12**: new §12.4.1 "Anatomy of the crossing" (E=0 chirality decoupling,
  product form e^{2mL*} = −tan(θ0/2)tan(θL/2), isolated-wall bound states
  E = −m sin θ0 / +m sin θL, hybridization picture); `ch12_chiral_bag.py`
  written (was cited but nonexistent); plateaus/jump refreshed to
  +0.3640/−0.6361/−1.00015; fractional-law point −0.4455.
- **ch13/ch14**: the two entropy conversions un-conflated; the chain gains the
  today factor (s/n_γ)₀ = 7.04; assembled prefactor 2.3e-2 → **1.6e-1**;
  required P×·ε_CP 2.7e-8 → **3.8e-9**; n_bag ~ T³ normalization flagged
  [Postulate]; validation count corrected (390 pairs, not 2000); η_B
  observation value re-attributed with honest error bar.
- **ch15**: iso-energy motivation direction fixed (the ladder partner of mode n
  is L₀/n — the shortcut is an ad-hoc UV cure); commutator theorem restated
  (violation lives in [φ,π] and unequal-time spacelike [φ,φ]) and now computed
  by `ch15_imn.py`; 3+1D plateau value tail-corrected to 0.236.
- **ch16**: conformal-time mode frequency fixed (ω² = k² + m²a²); 3D
  Weyl/Cotton fine print added; two oversold framings softened.
- **ch17**: the dynamical dictionary re-scoped — (17.4) is a *definition*; the
  −(L̇/L)D̂ discrepancy operator between moving-wall and metric pictures
  displayed; Parker identification qualified; `ch17_dictionary_check.py` made a
  real two-matrix check (4e-6 at the discretization floor; was tautological).
- **ch18**: (18.3) restricted to rectangular cells (the 26–29% sheared-cell
  failure is now a printed negative control); Metric Emergence restated on the
  Gram matrix g_ab = e_a·e_b; Conformal Obstruction Theorem rewritten (Cotton
  tensor in 3D; loss of × polarization and SO(3) covariance; false
  light-bending claim removed — Schwarzschild's conformally flat slices stated
  in print); `ch18_metric_updates.py` written (was cited but nonexistent).
- **ch20**: "Noether" dressing replaced by similarity/homogeneity language;
  derived-vs-dialed scorecard added; contracting-axis bounce caveat
  (n_min = |C|) added; simulated hopping form disclosed.
- **ch21**: per-cell curvature offset declared as a fourth [Postulate] (flat
  closure forbids turnaround — stated); growing-mode amplitude C₊ = 9.6×10⁻⁴
  quoted so the 1168 turnaround is checkable; 1.062/1.686 derived in text;
  collapse sim honestly described as classical per-cell FRW.
- **ch22**: cutoff identification *derived* (exact lattice heat kernel
  e^{−z}I₀(z): continuum Weyl form above s ~ L̄², saturation at one mode/cell
  below); Planck-Cell Prediction audited as a consistency inversion (three
  premises; falsifiable at its joints; the "could have come out at the
  electroweak scale" overclaim removed); voxel caveats added.
- **ch23**: obstruction gloss synced to corrected ch18 (× polarization, not
  light bending).
- **ch24**: new §24.2.1 — the −1.008/−2.05/96% figures identified as simulation
  residuals with their three mundane sources; the model's genuine falsifiable
  deviations listed with scales (Planck-suppressed (ℓ_P/r)² Newton corrections;
  GW birefringence at the cell scale; dispersion must start at n = 2 or the
  model is dead today; Weyl-vector constraint; boundary-GIM kill test).
- **ch25**: opener synced to corrected ch18 (TT sector completed by shears).
- **ch26/ch27**: microcausality/no-interacting-extension added to not-claimed
  list and as open problem 12; item 11 (declared state of a suddenly enlarged
  box) added; Landau–Zener item corrected (net charge protocol-independent —
  only pair content is dressed); executioner-item mapping fixed; pump number
  synced to +1.0000; ch27 item 1 now names the curvature-offset import.
- **appA**: unitarity/defect index direction fixed (old-in-new closes the CAR;
  new-from-old saturates at 1−w_k); squeezed-state coefficient C = β(α*)⁻¹.
- **appB**: "plateau is exact" lemma replaced by the APS decomposition;
  N_eff normalization (scalar ≡ 1) stated; boundary-term caveat added.
- **appC/appD**: provenance numbers synced to what scripts actually print;
  Mixmaster residual normalization defined; q-collision ruling added.
- Cross-file: ch08 equation renumbering (8.14)–(8.16) propagated to ch09;
  numpy < 2.0 compat aliases (np.trapezoid) in ch08 + ports; Fig. 10.1
  schematic sign fixed; ch14_crossing_density.py tautological counter fixed.

## Fix round (2026-08-21): reproducibility criticals (C1/C2 of `thesis/_REWORK_TEMP.md`)

The Part III provenance gap is closed. New chapter-named scripts, adapted from the
`ports/` archive with contract-conformant output paths, all run and reproduce every
chapter headline number:

| new script | source | key prints (verified vs chapter) |
|---|---|---|
| ch19_fidelity.py | ports/dewitt_fidelity.py | shear χ 0.049/5.32/5.32/3.73 = PT to 4 digits; diagonal machine-zero; → data/ch19_dewitt_fidelity.json, Fig 19.2 |
| ch20_kasner_sim.py | ports/kasner_sim.py | p = (−0.2855, 0.4284, 0.8571); Σ_K 0.9997 / 0.3333; σ_n = 40 @ ⟨n⟩ = 4100 (flat run); → data/ch20_kasner_traj.npz, Figs 20.2–20.4 |
| ch21_friedmann.py | ports/newton_sim.py §B | deviations 0.0090/0.0091; exponents 0.664/0.498; update-interval check 4.8e-3 (new); → Fig 21.2 |
| ch21_collapse.py | ports/newton_sim.py §C | growth 1.082; turnaround 1168; collapse 1863; δ_lin 1.170/1.867; → data/ch21_collapse_traj.npz, Fig 21.3 |
| ch22_casimir_chain.py | ports/newton_sim.py §A | endpoints 0.135/3.865/24.000 now asserted; core ratio 0.956; → Fig 22.2 |
| ch22_dispersion.py | ports/stiffness_dispersion.py | band [−0.178, +0.060] (k=2 outlier now stated in (22.3)); medians −0.131 (k≤8) / −0.142 (all); → data/ch22_stiffness_dispersion.json, Fig 22.3 |
| ch24_newton_lattice.py | ports/newton_lattice.py | A = 0.0766; slope −1.008; force −2.053; → Fig 24.1 |
| ch24_mixmaster.py | ports/mixmaster.py | drift 8.9e-13 with hard 1e-10 discard gate (new); u-map to 6.2e-4; → data/ch24_mixmaster_traj.npz, Fig 24.3 |
| **ch25_graviton.py** | **written fresh** (method of App. C.4; ports/graviton.py was 45³/linear-fit and could not print the table) | **reproduces the ch25 table to every digit** (−0.987/−1.008 at s=3 etc.); unit tests flat 1.2e-15 / shear 5.3e-9; raw diffeo Π +0.39; → data/ch25_graviton_results.json, Figs 25.2–25.3 |
| ladder.py | new shared module | App. C.3's promised module (used by ch20/ch21) |

`run_all.py`: `--report` implemented (writes data/provenance.json with captured
per-script output); missing scripts are now failures, not silent skips.
`make_schematics.py` no longer hand-draws ch25_fig2 (data figure, owned by
ch25_graviton.py). Text synced where scripts contradicted it: ch25 (5e-9 unit
test, 0.39 raw diffeo, footer), (22.3) band + medians, App. C.3 (4.8e-3 measured),
C.4 (rms residuals, not fit covariances), C.5 (endpoint asserts as implemented),
C.6/D.7 (provenance format), ch19/ch20/ch21/ch22/ch24 validation footers.

## Script status

| v15 script | Provenance | Status |
|---|---|---|
| ch08_bag_spectrum.py (`bag1d`) | new (physics from v14 files1) | **written + self-tests pass** (A: 3.6e-8; B: mirror 7e-9, η ~ 1e-10; C: 3.8e-7; D: 8.7e-8) |
| ch14_crossing_density.py | **new derivation** | **written + validated** (390/390 existence+uniqueness; locations 1e-6; demo L* = 0.44302; flow dir consistent; P_x/eps_CP tables) |
| ch11_wall_density.py | new | **written + run** (nets = 0; regional −0.058 @ loose regulator → +0.067 benchmark) |
| ch12_waterfall.py | new | **written + run** (crossing observed 0.44314 vs 0.44302) |
| ch09_eta_toy.py | new | **written + run** (3 regulators → 2; mirror → 0) |
| ch03_overlaps.py | new (formulas re-derived) | **written + run** (closed form vs quadrature 8e-16; completeness 1e-8; Lemma values exact) |
| ch06_norms.py | new | **written + run** (N_con exact; KE preservation 6e-5; σ(n1) measured) |
| ch06_golden_rule.py | port-rewrite of golden_rule_toy | **written + run** (0.76749/0.85112/0.90373 — reproduces v14 exactly) |
| ch02_iso_energy.py | restyle of v13 plot | **written + run** |
| ch14_eta_assembly.py | new | **written + run** (prefactor 2.27e-2; required P_x·ε = 2.71e-8) |
| ch17_dictionary_check.py | new | **written + run** (spectra equal to 2e-16) |
| ch04_ladder_dynamics.py | new | **written + run** (Bessel propagator 4e-16; R(1) = 2.52) |
| ch05_regimes.py | new | **written + run** (regime diagram) |
| ch10_mirror_spectrum.py | new | **written + run** (mirror spectra vs δ; η ~ 1e-9 at solver tol) |
| ch10_radial_similarity.py | **new independent validation** | **written + run** (dressed radial ODE vs m cos δ bag: 4.5e-11; ratio = e^{m_I r} to 4.4e-11) |
| ch12_pump_control.py | new | **written + run** (**pump +1.0000, control −0.0000** — quantization exact under Richardson) |
| ch16_kasner_circle.py | new | **written + run** (circle + BKL cobweb; orbit 4.3→3.3→2.3→1.3→3.3333) |
| make_schematics.py | new | **written + run** (all 28 concept figures, reproducible) |
| regen_part1_figs.py | new (replaces adapted v13 figs) | **written + run** (ch03_fig3, ch06_fig1 w/ deficit check 8e-5, ch08_fig1) |
| ch15_imn.py | new (replaces port) | **written + run** (π/8 closed form; 3+1D plateau; Casimir negative result) |
| make_anims_v15.py | new | **written + run** (ch04_anim_spreading.gif, ch12_anim_pump.gif) |
| ports/*.py (kasner_sim, newton_sim, newton_lattice, stiffness_dispersion, dewitt_fidelity, mixmaster, graviton, quench1d, chiral1d2, chiral_exact, check1_radial_similarity, imn_uv_regularization, anim makers) | v14/v13 verbatim | **archive only** since the 2026-08-21 fix round — chapter-named adaptations above are the scripts of record (note: ports/graviton.py docstring wants κ_TT > 0, contradicting (25.6); do not cite it) |
| ch19_fidelity, ch20_kasner_sim, ch21_friedmann, ch21_collapse, ch22_casimir_chain, ch22_dispersion, ch24_newton_lattice, ch24_mixmaster, ch25_graviton, ladder (module) | fix round 2026-08-21 (see table above) | **written + run + all chapter numbers reproduced** |

## Figure status (thesis/figures/) — COMPLETE: every referenced figure exists (77 files)

- **Generated fresh by v15 scripts (data):** ch02_fig1/2, ch03_fig1, ch04_fig1/2/3, ch05_fig1, ch06_fig2/3/4, ch09_fig1/2, ch10_fig2/3, ch11_fig1, ch12_fig3/4, ch14_fig2/3, ch16_fig2/3, ch17_fig2, ch25_fig2.
- **Generated by `make_schematics.py` (concept drawings):** ch01_fig1/2, ch07_fig1/2, ch08_fig2/3, ch10_fig1, ch13_fig1/2, ch14_fig1, ch15_fig3, ch16_fig1, ch17_fig1, ch18_fig2/3, ch19_fig1, ch20_fig1, ch21_fig1, ch22_fig1/4, ch23_fig1/2, ch24_fig2, ch25_fig1, ch26_fig1/2, ch27_fig1.
- **Adapted (copied from validated v13/v14 runs; restyle optional):** ch03_fig2, ch11_fig2/3, ch12_fig1/2, ch18_fig1, anims (ch20×2, ch21, ch24, ch25). (ch03_fig3, ch06_fig1, ch08_fig1, ch15_fig1/2 regenerated fresh earlier; ch19_fig2, ch20_fig2/3/4, ch21_fig2/3, ch22_fig2/3, ch24_fig1/3, ch25_fig2/3 regenerated fresh by the 2026-08-21 chapter scripts.)
- **Animations (7 total):** new — ch04_anim_spreading, ch12_anim_pump (`make_anims_v15.py`); adapted from validated archive — ch20×2, ch21, ch24, ch25. All embedded in their dynamics sections.

## Known deviations from legacy text (deliberate, validated)

1. **Flow direction at a crossing**: legacy prose said the gap level "dives into the sea" at θ=(−2,+2); the legacy *numbers* (Q_vac plateaus) and the new measurement agree the level **rises** (sea → positive branch) for sin(Δ/2) > 0, Q_vac jumps by −1, pumped ΔQ = +1. v15 text follows the numbers.
2. **N_crit**: legacy "15–20 spectators" holds only at n1 = 2–3; v15 quotes the measured σ(n1) ≈ 1/n1 law instead.
3. **Spectator-coefficient closed form**: the v13 pedagogical summary's c_{n2} formula missed a 1/√s factor; v15 (3.8) is re-derived and validated against quadrature (8e-16).
4. **K-matrix in legacy report Part II §4**: displayed matrix had the m_I terms misplaced (typo); the v15 form (10.5) reproduces ψ″ = −p²ψ and all downstream results.
5. **B-coefficient (10.7)**: legacy appendix relation B = pA/(E+m_R−m_I) is superseded by B = A(E+m_R−m_I)/p, which reproduces tan(pL) = −p/(m cos δ) exactly; validated numerically.

## Cloud-sync cautions (two distinct traps, both encountered in this build)

1. **Appended null bytes**: OneDrive intermittently appends `\x00` runs to freshly
   edited files (crashes Python imports loudly). Strip them — but see trap 2 first.
2. **Stale truncated reads**: the sandbox mount can serve a *partial* (truncated)
   copy of a recently written file for many minutes. Therefore: (a) never judge file
   completeness from the mount; (b) **never write a file back through the mount based
   on a mount read** (a null-strip sweep on a truncated copy would truncate the real
   file). Clean nulls only on files verified complete, or do the cleanup on the
   Windows side. Word counts / reference checks run on the mount undercount fresh
   files; the editor-side copies are authoritative.
