# Chapter 18 — A lattice of cells: the triad, the missing shears, and the conformal obstruction

> **Author: Fable 5 · Complexity ◇4 · Status: draft v1 (pending global review pass)**
> Depends on: Ch. 3, 16, 17. Feeds: Ch. 19–20 (shape space, dynamics), Ch. 23 (frame gauging), Ch. 25 (the TT sector lives in what this chapter builds).

---

The dictionary of Ch. 17 converts one box into one metric — a *homogeneous* geometry, a single dial per direction. General relativity is a theory of geometry that varies from place to place. The minimal upgrade is a **lattice of cells**: tile space, give each cell its own shape, and let the iso-energy machinery of Part I act cell by cell. This chapter builds the kinematics of that upgrade and discovers two things, one gratifying and one humbling. Gratifying: the cell description automatically produces the *triad* (vielbein) — the variable on which every modern gauge-theoretic formulation of gravity is built, and the one that standard treatments must postulate. Humbling: the transition operators inherited from Part I generate only *half* the metric — and the missing half is not a detail, because a theorem stands behind it: geometry built from dilations alone is conformally flat, and conformally flat gravity is experimentally dead (no light bending, no gravitational waves). Identifying the missing generators — the iso-energy **shears** — and verifying that the completed set closes on $\mathfrak{gl}(3)$ is this chapter's constructive payoff.

## 18.1 The cell lattice and its variables

Tile space with cells; the cell at lattice site $\mathbf x$ is a parallelepiped spanned by three **edge vectors** $\mathbf e_1(\mathbf x), \mathbf e_2(\mathbf x), \mathbf e_3(\mathbf x)$ (not necessarily orthogonal, not necessarily equal in length). The matter content of each cell is the box quantum mechanics of Part I, with quantum numbers $\mathbf n(\mathbf x)$. Two questions immediately: what is the coordinate-free description of a cell's shape? and what do the transition operators of Ch. 3 do to it?

## 18.2 The Metric Emergence Theorem

A cell's shape is exactly the information needed to measure lengths within it: the extent of the cell in direction $\hat u$ is $\ell(\hat u)^2 = \sum_a (\mathbf e_a \cdot \hat u)^2$, a quadratic form. Quadratic forms on directions *are* metrics:

> **Metric Emergence Theorem [Theorem].** The cell shape is completely and coordinate-freely characterized by
>
> $$g_{ij}(\mathbf x) \;=\; \sum_{a=1}^{3} e_a^i(\mathbf x)\,e_a^j(\mathbf x) \tag{18.1}$$
>
> — symmetric and positive-definite by construction (a sum of outer products of three independent vectors) — and the directional iso-energy transitions act on it as computable rank-one updates:
>
> $$\hat T^{(a)}_+:\quad g_{ij} \;\longrightarrow\; g_{ij} + \Big[\Big(\tfrac{n_a + 1}{n_a}\Big)^2 - 1\Big]\,e^i_a\,e^j_a \;\;\approx\;\; g_{ij} + \frac{2}{n_a}\,e^i_a e^j_a \quad (n_a \gg 1). \tag{18.2}$$
>
> The particle energy is the coordinate-invariant pairing
>
> $$E \;=\; \frac{\pi^2}{2m}\,n_a\,n_b\;g^{ab}, \tag{18.3}$$
>
> with $g^{ab}$ the inverse metric — the 3D dictionary of Ch. 17 for arbitrary cell shapes.
>
> *Proof.* (18.1): symmetry is manifest; positivity because $g_{ij}v^iv^j = \sum_a(\mathbf e_a\cdot\mathbf v)^2 > 0$ for $\mathbf v \ne 0$ when the $\mathbf e_a$ span. Components: $g_{ii} = |\mathbf e_i|^2$ (squared edge lengths), $g_{ij} = \mathbf e_i\cdot\mathbf e_j$ (edge angles); $6$ independent entries, matching the spatial metric of GR. (18.2): the transition dilates one edge, $\mathbf e_a \to \tfrac{n_a+1}{n_a}\mathbf e_a$; insert into (18.1) and subtract. (18.3): for the rectangular case it is Ch. 17's formula; for general parallelepipeds, transform coordinates to the cell frame — $n_a$ transforms as a covector and $g^{ab}$ as an inverse metric, leaving $E$ invariant. $\blacksquare$

The metric is not postulated; it is the unique coordinate-free bookkeeping of cell shape, and it *moves* when transitions fire. **[Computed]** (18.2) exercised over random transition sequences in `ch18_metric_updates.py` (Fig. 18.1).

![Metric updates on the lattice](figures/ch18_fig1_cell_lattice.png)

*Figure 18.1 — Cells as geometry. A patch of the lattice with per-cell edge vectors (top); a sequence of directional transitions deforming one cell's shape ellipsoid and its $g_{ij}$ entries (bottom).*

## 18.3 The triad recognized — and its gauge redundancy

Equation (18.1) deserves a stronger name than "bookkeeping". It is exactly the **vielbein (triad) decomposition** $g_{ij} = e^a_i e^a_j$ that underlies every gauge-theoretic formulation of gravity — Cartan's, Kibble–Sciama's (Ch. 16.9), and the Ashtekar variables of loop quantum gravity. The accounting that matters:

$$\underbrace{9}_{\text{triad } e^a_i} \;=\; \underbrace{6}_{\text{metric } g_{ij}} \;+\; \underbrace{3}_{\text{local } SO(3) \text{ frame rotations}},$$

since $\mathbf e_a \to R_{ab}\mathbf e_b$ with $R \in SO(3)$ leaves (18.1) untouched. In the standard treatments the triad is *introduced* so that spinors can couple to gravity, and its rotation redundancy is decreed. Here both arrive unbidden: the edges are physical objects of the cell model, and the freedom to relabel them is *visibly* a redundancy of description — the first conceptual step of Ch. 23's frame-group gauging, supplied by the microscopic theory for free. (Part II quietly needed this too: spinor boundary conditions live on frames, not metrics.)

## 18.4 Flatness of the single-cell kinematics

> **Flatness (Commutativity) Theorem [Theorem].** Directional transitions on a single cell commute:
>
> $$\big[\hat T^{(a)}_\pm,\, \hat T^{(b)}_\pm\big] = 0 \qquad \forall\, a, b.$$
>
> *Proof.* $\hat T^{(a)}$ acts only on the pair $(n_a, \mathbf e_a)$; different $a$ touch disjoint variables; both orderings land on the same final configuration (the explicit two-path check is one line per component). $\blacksquare$

The honest reading, calibrated carefully: commutativity is the discrete analogue of a *flat connection* — paths between the same endpoints in shape space compose identically — which is exactly right for a *single* cell, whose geometry is homogeneous (one point's worth of metric; nothing to curve). Curvature, in the gauge-theoretic sense Ch. 23 will need, must come from the *lattice*: transport around a closed loop of cells whose geometries disagree. The theorem is thus a consistency check with a built-in forecast: when inter-cell coupling arrives (Ch. 21–22), the deformed non-commutativity of transport is where the Riemann tensor will live. (And recall Ch. 3's slogan — kinematics telescopes, dynamics remembers: physical *amplitudes* along different paths already differ through spectator dressing; it is the bare kinematics that is flat.)

## 18.5 What the inherited generators span — and what they miss

Now the audit. Decompose the metric's six dimensions at a point, GR-style (Ch. 16.4):

$$g_{ij} = e^{2\phi}\,\hat g_{ij}, \quad \det\hat g = 1: \qquad 6 \;=\; \underbrace{1}_{\text{conformal }\phi} + \underbrace{5}_{\text{unimodular shape }\hat g}.$$

Group-theoretically, metrics live on the coset $GL(3)/SO(3)$, whose tangent space is the six-dimensional space of symmetric matrices. The directional transitions (18.2) add multiples of $e^i_ae^j_a$ — for a rectangular cell, *diagonal* matrices. Three independent directions: the **Cartan subalgebra** of the symmetric space. Composing them reaches any diagonal metric — the conformal factor plus two of the five shape directions (the anisotropic *stretches*; precisely the Bianchi-I configurations Ch. 20 will evolve). What no composition of them reaches: the three **off-diagonal** components — the *angles* between edges. No operator of Part I changes the angle between two edges; the rotational transitions one might hope to repurpose act on $\mathbf n$, not on the cell, and do not fill the gap.

> **Definition (iso-energy shear generators) [Definition].** The missing generators are
>
> $$\hat S^{(ab)}:\quad \mathbf e_a \;\to\; \mathbf e_a + \varepsilon\,\mathbf e_b \qquad (a \ne b), \tag{18.4}$$
>
> tilting the cell at fixed (to $\mathcal O(\varepsilon)$) volume. Iso-energy versions exist because the invariant energy (18.3) defines level sets in the *full six-parameter* metric space, not only along its diagonal: the selection-rule logic of Ch. 2, applied to $g^{ab}$ rather than to $L_a$. In fact the shears are *better* than iso-energy-by-selection-rule — Ch. 19 proves their energy cost vanishes at first order identically, no quantum-number jump required.

## 18.6 Why the gap is fatal if unfilled: the Conformal Obstruction Theorem

One might hope the diagonal sector is "enough gravity to start with". It is not, and the obstruction is a theorem about experiment, not aesthetics.

> **Conformal Obstruction Theorem [Theorem; Standard ingredients].** A theory whose geometry sector is generated by local dilations alone produces only metrics of the form $g_{ij}(\mathbf x) = \Omega^2(\mathbf x)\,\delta_{ij}$ — conformally flat. Such metrics have identically vanishing Weyl tensor. Consequently the theory possesses:
>
> 1. **no light deflection of the GR type** — the bending of starlight around the sun is a Weyl-curvature effect; the unique conformally-flat competitor theory (Nordström's, 1913) predicts zero deflection and was refuted by the 1919 eclipse expedition;
> 2. **no gravitational waves** — the two LIGO polarizations are transverse-traceless Weyl modes; a conformal factor cannot wave.
>
> *Proof structure.* Dilations exponentiate to position-dependent rescalings of a fixed reference metric, which is the definition of conformal flatness. The Riemann tensor splits into the Ricci part (algebraically tied to local sources by Einstein's equations) and the Weyl part (the free, propagating remainder); for $g = \Omega^2\delta$, direct computation gives $C^{i}{}_{jkl} \equiv 0$ (the Weyl tensor is conformally invariant and vanishes for flat space). Light bending and GW polarizations live in the Weyl part **[Standard]**. $\blacksquare$

The obstruction is *kinematic*, which is what makes it sharp: it is not that dilation-only dynamics has the wrong action (Ch. 23 will add that, too — pure-Weyl invariants start at four derivatives), but that the **configuration space itself** is too small. No clever Lagrangian on a conformally-flat configuration space produces a Weyl tensor. Either the shears (18.4) exist as physical generators of the cell dynamics, or the program stops here. The thesis takes the first branch, with consequences immediately testable in this Part: the shear sector must carry the gravitational waves — a prediction Ch. 25 confirms with the measured coefficient.

![The obstruction and its cure](figures/ch18_fig2_gl3_decomposition.png)

*Figure 18.2 — The generator audit. The six metric directions (symmetric matrices) split into the diagonal/Cartan sector reached by Part I's transitions (conformal + 2 stretches) and the off-diagonal shear sector (3 angles) reached only by $\hat S^{(ab)}$; the three frame rotations are gauge. Callouts: which experimental physics lives in which sector.*

## 18.7 Closure: the full kinematic algebra

Assemble the inventory: three directional dilations (diagonal, symmetric), three shears $\hat S^{(ab)}$ (off-diagonal, symmetric, by symmetrizing (18.4)), three frame rotations (antisymmetric, gauge). Their commutators close on the Lie algebra

$$\mathfrak{gl}(3) \;=\; \underbrace{\mathfrak{so}(3)}_{\text{3: gauge rotations}} \;\oplus\; \underbrace{\text{Sym}(3)}_{\text{6: metric directions}}, \tag{18.5}$$

with the symmetric part acting transitively on the space of metrics $GL(3)/SO(3)$ and the trace of $\mathfrak{gl}(3)$ — the overall dilation — being exactly the Weyl direction of Ch. 17.5. *Nothing missing, nothing extra*: the completed cell kinematics carries precisely the local degrees of freedom of a spatial geometry plus the frame gauge that spinors require. **[Computed]** the closure relations and the transitivity statement are verified symbolically/numerically in `ch18_metric_updates.py --algebra`.

## 18.8 Summary

A lattice of cells upgrades the dictionary to a *field* of metrics, with the triad derived (18.1), transitions acting as rank-one updates (18.2), single-cell kinematics flat (as it must be), and the energy pairing (18.3) coordinate-invariant. The audit found Part I's generators spanning only the diagonal sector; the **Conformal Obstruction Theorem** shows that stopping there forfeits light bending and gravitational waves; the iso-energy **shears** (18.4) complete the algebra to $\mathfrak{gl}(3)$ (18.5). Two questions now stand, and they organize the next two chapters: *what is the quantum geometry of the completed shape space?* (Ch. 19 — where the diagonal and shear sectors turn out to play structurally different roles, in exact parallel to the DeWitt split), and *what dynamics does the sector traffic of Part I induce on it?* (Ch. 20 — where the coupling law selects, or fails to select, general relativity's orbits).

---

**Validation.** `ch18_metric_updates.py`: rank-one update sequences and shape-ellipsoid evolution (Fig. 18.1); `--algebra`: commutator closure of {dilations, shears, rotations} on $\mathfrak{gl}(3)$ and transitivity on sampled metrics. The Conformal Obstruction Theorem's computational content ($C \equiv 0$ for $\Omega^2\delta$) is checked symbolically in the same script.
