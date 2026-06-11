# v15 — The Expanding Quantum Box (thesis build)

The complete, self-contained, pedagogical statement of the theory: the iso-energy
framework, the arrow of expansion, spectral-flow baryogenesis, and emergent gravity.
No prior material (papers/, v13/, v14/) is required or referenced by the thesis text.

## Layout

- **`thesis_skeleton.md`** — the construction blueprint: chapter specs, theorem-name
  registry, figure/script plans, complexity ratings, model assignments. Read first.
- **`thesis/00_master.md`** — title, abstract, reading guide, linked TOC, status board.
- **`thesis/chNN_*.md`** — the chapters. Fable 5 chapters are full text (draft v1);
  Opus 4.8 slots are binding section-level specs marked `OPUS 4.8 WRITING SLOT`.
- **`thesis/figures/`** — generated + adapted figures (status in `coverage_matrix.md`).
- **`scripts/`** — validation suite (`run_all.py`; every quoted number printed by a
  named script). `scripts/ports/` holds the validated legacy scripts pending rename.
- **`coverage_matrix.md`** — internal provenance/status; **not** part of the thesis.
- **`data/`** — cached simulation outputs.

## Build state

- **All 27 chapters + 4 appendices written** (draft v1, single author: Fable 5) — the
  document is complete end to end, ~75k words of text plus ~65 figures.
- **18 validation scripts written and run, all checks passing**, including the thesis's
  new derivations re-validated independently in this build: the One-Crossing Theorem
  (390/390), the Dressing Theorem by direct radial ODE (4.5e-11), the quantized pump
  (+1.0000) with its control (0.0000), and the Bessel propagator (4e-16).
- All schematic/concept figures generated programmatically (`make_schematics.py`) —
  even the drawings are reproducible.
- Legacy heavy simulations: validated figures adapted; ports staged in `scripts/ports/`
  (mixmaster re-run verified end-to-end; remaining re-runs are cosmetic restyling).

## Next steps (polish, not construction)

1. Global review pass: gates G1–G4 in `thesis_skeleton.md` §6 (cross-references,
   theorem-name registry, number-to-script tracing).
2. Re-run the heavy port tier (`kasner`, `graviton`, `newton_sim`) under v15 styling.
3. Optional: export to PDF/docx for the supervisor.
