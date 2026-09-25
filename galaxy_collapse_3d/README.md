# Galaxy collapse in 3D — a cosmological simulation on a discrete cell lattice

A faithful 3-D extension of the gravitational-collapse simulation in **The Expanding
Quantum Box** (v15 thesis, Part III). It evolves a proto-galaxy over-density on a
discrete grid of cells using the model's *own* derived dynamics and renders the
collapse as a 3-D animation. **Nothing in `v15/` is modified or imported** — the
physics is transcribed from the thesis and first reproduces its validated numbers.

```
python galaxy_collapse_3d.py            # full run: validation + summary + GIF
python galaxy_collapse_3d.py --quick    # fast, low-res test render
python galaxy_collapse_3d.py --no-anim  # numbers + static summary only
```

Outputs (written next to the script):
- `galaxy_collapse_3d.gif` — the 3-D collapse animation
- `galaxy_collapse_summary.png` — static multi-panel summary
- `galaxy_collapse_traj.npz` — the trajectory archive

## The model, exactly as used

Each discrete cell of space **is** a Friedmann patch — this is the thesis's central
identity, not an analogy:

| Thesis | Equation | Role here |
|---|---|---|
| **Ch. 17** dictionary | a box of size `L` ≡ a unit cell with metric `g = L²`; `a(t)=L(t)` | every cell carries its own scale factor |
| **Ch. 18** cell lattice | `g_ij = Σ_a e_a^i e_a^j`; gravitational variable `δα = ln a` | the discrete grid of cells |
| **Ch. 21** constraint closure | `ä = −(H₀²/2)·Ω/a²`, `Ω = 1+δ₀(x)`, `ȧ₀ = H₀√(Ω/a₀+(1−Ω))` | separate-universe collapse dynamics |
| **Ch. 21** mass conservation | `δ(x,t) = (1+δ₀)(a_bg/a)³ − 1` | density contrast |
| **Ch. 24** lattice Poisson | `κ ∇²δα = −c ρ δ`, `δα = −Φ` (FFT solve) | the Newtonian potential well |

An **over-dense cell (Ω>1) is a slightly closed universe**: it expands ever slower
than the background, reaches **turnaround**, and **recollapses** — gravitational
instability with no force law postulated. The 1-D version is the thesis's
Animation 21.A; this script lifts it to a genuine 3-D grid and renders the galaxy.

### Exact solver

The cell equation `ä = −(H₀²/2)Ω/a²` with constraint-consistent initial data is
solved **exactly**: for an over-dense (closed) cell the solution is the recollapsing
cycloid `a=A(1−cos η)`, `t=B(η−sin η)`, `A=Ω/2(Ω−1)`, `B=Ω/2H₀(Ω−1)^{3/2}`. This is
the model solved in closed form (no stiffness from integrating the recollapse cell
by cell). The run cross-checks it against a direct ODE integration of the central
cell — agreement to `4.8×10⁻⁵`.

### Faithful 3-D matter map

For the spherically-symmetric proto-galaxy the tracer-matter displacement is the
exact 3-D radial analogue of the thesis's 1-D "walls falling toward the
over-density" (the comoving cumulative-width map of `make_collapse_anim.py`):

```
r̃(R,t) = (1/a_bg) ∫₀ᴿ a(r,t) dr      # background-normalized physical radius
```

Inner shells whose cells collapse (`a → a_min`) fall toward the centre while the
outskirts keep riding the Hubble flow — matter attracting matter, inside the model.

**Virial core (visualization only).** The separate-universe model has no
shell-crossing, so spherical collapse runs to a point. For display, collapsed cells
are shown as a finite virialized halo of radius `R_VIR` (a standard structure-
formation visualization) rather than a singular point — this affects *only* the
rendering, never the dynamics or the validation numbers.

## Validation against the thesis (Ch. 21)

The run prints, and reproduces, the thesis numbers for the central cell:

| quantity | this 3-D sim | thesis (Ch. 21) |
|---|---|---|
| linear growth exponent `d ln δ / d ln a` | **1.082** | 1.08 (theory 1) |
| central turnaround `a_bg/a₀` | **1172** | 1168 |
| central collapse `a_bg/a₀` | **1899** | 1863 |
| cycloid vs direct-ODE cross-check | rel. diff `4.8×10⁻⁵` | — |

with the spherical-collapse thresholds `δ_lin = 1.062` (turnaround) and `1.686`
(collapse) framing the arc.

## Adding angular momentum — the stable-orbit velocity curve

`galaxy_rotation_3d.py` answers: *what tangential-velocity profile puts every shell
on a stable orbit so nothing collapses at large t?*

```
python galaxy_rotation_3d.py            # rotation curve + 3-D cold-vs-rotating animation
python galaxy_rotation_3d.py --no-anim  # numbers + the v_c(r) figure only
```

The model's gravity is **exactly Newtonian** (Ch. 24: measured `1/r` potential,
`1/r²` force), so a shell at radius `r` feels `dΦ/dr = G·M(<r)/r²`. With specific
angular momentum `ℓ`, the effective potential `V_eff = Φ + ℓ²/2r²` has a circular
(non-collapsing) orbit at `V_eff'=0`, giving the **equilibrium / no-collapse
condition**

```
v_t(r) = v_c(r) = sqrt( r · dΦ/dr ) = sqrt( G·M(<r)/r )
```

stable where `d(r·v_c)/dr > 0` (Rayleigh). For the proto-galaxy's Gaussian mass
`ρ ∝ exp(−(r/R₀)²)`, `M(<r)/M_tot = erf(x) − (2/√π)x·e^{−x²}` (`x=r/R₀`), the curve:

- **rises as solid-body** `v∝r` in the core (`ρ≈const`, `M∝r³`) — measured slope **+0.952**,
- **peaks** near `r ≈ 1.5 R₀`,
- **declines Keplerian** `v∝1/√r` in the outskirts — measured slope **−0.500**,
- is **not flat**: a flat curve needs `M(<r)∝r` ⇔ `ρ∝1/r²` (isothermal halo). Since
  E-QFT reproduces Newton exactly, an isolated luminous mass gives a *declining*
  curve — the dark-matter problem, in the model's own language.

The N-body run confirms it: cold tracers (`v=0`) collapse (`r₁/₂` 4.31→2.7), while
tracers given `v=v_c(r)` stay on circular orbits (`r₁/₂` 4.31→4.32, orbit drift
0.3%) — a **stable rotating galaxy that never collapses**. Outputs:
`galaxy_rotation_curve.png`, `galaxy_rotation_3d.gif`, `galaxy_rotation_traj.npz`.

## Conventions

`ħ = c = 1`. Collapse run: `H₀ = 1`, `a₀ = 0.02`, `a_min = 0.25 a₀`,
`δ₀(r) = 0.08·exp(−(r/4)²)` — the thesis's exact initial data as a function of 3-D
radius. Rotation run: `G = 1`, `R₀ = 4`, Newtonian gravity (Ch. 24). The dynamical
parameters are fixed by the model; only rendering parameters are free.
