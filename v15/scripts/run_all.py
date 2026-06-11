#!/usr/bin/env python3
"""
v15 master driver: regenerates every thesis figure and prints every quoted number.

Contract (App. C):
  * every number quoted in the thesis is printed by exactly one script;
  * every script is standalone (numpy/scipy/matplotlib only) and writes its
    figures to ../thesis/figures/ with the chNN_figM_name.png convention;
  * scripts are grouped by cost; `python run_all.py --fast` runs the cheap tier
    only, `--all` includes the heavy simulation tier (minutes).

Each entry: (script, tier, what it validates).
"""
import argparse
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent

REGISTRY = [
    # (script, tier, validates)
    ("ch02_iso_energy.py",        "fast",  "iso-energy curves figure"),
    ("ch03_overlaps.py",          "fast",  "overlap closed forms; dilation matrix elements 10.4762/40.4938/160.4983"),
    ("ch04_ladder_dynamics.py",   "fast",  "band structure; Bessel propagator 4e-16; participation R(1)=2.52"),
    ("ch05_regimes.py",           "fast",  "sudden/adiabatic regime diagram"),
    ("ch06_norms.py",             "fast",  "norm theorems; N_con table; asymmetry vs N"),
    ("ch06_golden_rule.py",       "fast",  "golden-rule embedding theorem (Parseval vs N_con vs exact evolution)"),
    ("ch08_bag_spectrum.py",      "fast",  "MIT/chiral bag spectra; tan(pL)=-p/m checks (bag1d module self-test)"),
    ("ch09_eta_toy.py",           "fast",  "heat-kernel eta regulator on solvable example"),
    ("ch10_radial_similarity.py", "fast",  "bulk-phase dressing theorem (5e-10 eigenvalues; e^{m_I r} ratio)"),
    ("ch10_mirror_spectrum.py",   "fast",  "spectral mirror theorem; eta = O(1e-16) checks"),
    ("ch11_polarization.py",      "slow",  "convergence of truncated forms to polarization integral +0.0672"),
    ("ch11_wall_density.py",      "fast",  "vacuum charge density wall layers (new)"),
    ("ch12_chiral_bag.py",        "slow",  "master equation 5e-8; Q_vac law; crossing at L*=0.44302; pump +0.996 / control 0"),
    ("ch12_waterfall.py",         "fast",  "spectral-flow waterfall figure (new)"),
    ("ch12_pump_control.py",      "slow",  "pump +1.0000 vs control 0.0000 (spectral side, new)"),
    ("ch14_crossing_density.py",  "fast",  "crossing-existence criterion; nu maps; brute-force validation (new)"),
    ("ch14_eta_assembly.py",      "fast",  "eta_B assembly and budget tables (new)"),
    ("ch15_imn.py",               "fast",  "IMN phi^2 plateaus; per-mode decay; Casimir negative result"),
    ("ch16_kasner_circle.py",     "fast",  "Kasner circle / u-map figures"),
    ("ch17_dictionary_check.py",  "fast",  "moving-wall vs metric spectrum equality (new)"),
    ("ch18_metric_updates.py",    "fast",  "rank-1 metric updates figure"),
    ("ch19_fidelity.py",          "slow",  "splitting theorem: chi values; diagonal flatness 1e-13"),
    ("ch20_kasner_sim.py",        "heavy", "Kasner exponents to 3-4 digits; Sigma_K = 0.9997 / 0.3333"),
    ("ch21_friedmann.py",         "slow",  "Friedmann ladder |Delta alpha| <= 0.009; exponents 2/3, 1/2"),
    ("ch21_collapse.py",          "slow",  "growth exponent 1.08; turnaround/collapse thresholds"),
    ("ch22_casimir_chain.py",     "fast",  "cell energetics; chain relaxation endpoints (0.135 / 3.865 / 24.000)"),
    ("ch22_dispersion.py",        "slow",  "kappa_grad^vac ~ -0.13 vs matter +3 pi^2 n^2"),
    ("ch24_newton_lattice.py",    "slow",  "1/r exponent -1.008; amplitude 0.0766; force -2.05"),
    ("ch24_mixmaster.py",         "fast",  "BKL u-map sequence 4.3 -> 3.3 -> 2.3 -> 1.3 -> 10/3"),
    ("ch25_graviton.py",          "heavy", "kappa_TT/C -> -0.987/-1.008 at s=3; conformal sign; diffeo control"),
    ("make_schematics.py",        "fast",  "all drawn/concept figures (28), reproducibly"),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="cheap tier only")
    ap.add_argument("--all", action="store_true", help="include heavy tier")
    ap.add_argument("--only", type=str, default=None, help="comma-separated script substrings")
    args = ap.parse_args()

    tiers = {"fast"} if args.fast else ({"fast", "slow", "heavy"} if args.all else {"fast", "slow"})
    selected = [(s, t, w) for (s, t, w) in REGISTRY if t in tiers]
    if args.only:
        keys = [k.strip() for k in args.only.split(",")]
        selected = [(s, t, w) for (s, t, w) in selected if any(k in s for k in keys)]

    failures = []
    for script, tier, what in selected:
        path = HERE / script
        if not path.exists():
            print(f"[skip ] {script:30s} (not yet written) — {what}")
            continue
        print(f"[run  ] {script:30s} [{tier}] — {what}")
        t0 = time.time()
        r = subprocess.run([sys.executable, str(path)], cwd=HERE)
        dt = time.time() - t0
        status = "ok" if r.returncode == 0 else "FAIL"
        print(f"[{status:5s}] {script:30s} ({dt:.1f}s)")
        if r.returncode != 0:
            failures.append(script)

    print("\n=== run_all summary ===")
    print(f"selected: {len(selected)}, failures: {len(failures)}")
    for f in failures:
        print(f"  FAILED: {f}")
    sys.exit(1 if failures else 0)

if __name__ == "__main__":
    main()
