# Appendix D — Notation and glossary

> **Author: Fable 5 · Status: complete (symbol table + glossary). No chapter may silently deviate from this table.**

## D.1 Units and signature

| Convention | Value |
|---|---|
| Units | $\hbar = c = 1$ (Ch. 2 keeps $\hbar$ explicit for pedagogy) |
| Metric signature | $(+,-,-,-)$ |
| Reduced Planck mass | $M_P^2 = 1/8\pi G$ |
| 1+1D gamma matrices | $\gamma^0 = \sigma_3,\ \gamma^1 = i\sigma_2,\ \gamma_5 = \gamma^0\gamma^1 = \sigma_1$ |

## D.2 Foundations (Part I)

| Symbol | Meaning | Introduced |
|---|---|---|
| $L$, $L_a$ | box size; directional side lengths | Ch. 2 |
| $n$, $n_a$, $\mathbf n$ | mode quantum number(s) | Ch. 2 |
| $E_n(L) = n^2\pi^2\hbar^2/2mL^2$ | box spectrum | Ch. 2 |
| $s$ | expansion ratio $L'/L$ (e.g. $s = (n_1{+}1)/n_1$) | Ch. 2 |
| $\varphi_n(x) = \sqrt{2/L}\sin(n\pi x/L)$ | box mode functions | Ch. 2 |
| $\mathcal H = \bigoplus_L \mathcal H(L)$ | geometry-sector decomposition | Ch. 3 |
| $\iota$ | embedding map between sectors (extension-by-zero / restriction) | Ch. 3 |
| $\hat T_\pm^{(a)}$ | iso-energy transition operators (directional) | Ch. 3 |
| $\hat D = \tfrac12(\hat x\hat p + \hat p\hat x)$ | dilation generator | Ch. 3 |
| $\mathcal O_{n\to n+m} = \sqrt{n/(n+m)}$ | promoted-mode geometric overlap | Ch. 3 |
| $c_k(n,s)$ | spectator overlap coefficients | Ch. 6 |
| $\mathcal N_{\text{con}}(n_1,n_2)$ | contraction norm $= \frac{n_1-1}{n_1} + \frac{\sin(2\pi n_2/n_1)}{2\pi n_2}$ | Ch. 6 |
| $N$, $N_{\text{crit}}$ | spectator number; critical value for fixed suppression ($N_{\text{crit}}(10^6) \approx 14\,n_1$) | Ch. 6 |
| $g$ | sector coupling constant (golden rule); later slaved to matter (Ch. 21) | Ch. 6 |
| $\hat E$ | expansion-reversal operator | Ch. 7 |
| $J$, $q$ | tight-binding hopping; ladder quasimomentum | Ch. 4 |
| $\eta_{\text{ad}}$ | adiabatic parameter (regime diagnostic) | Ch. 5 |

## D.3 Charge from geometry (Part II)

| Symbol | Meaning | Introduced |
|---|---|---|
| $m_R = m\cos\delta$, $m_I = m\sin\delta$ | real/imaginary parts of CP mass $m\,e^{i\delta\gamma_5}$ | Ch. 8/10 |
| $p$, $q$ | momentum, $p^2 = E^2 - m^2$; below-gap $p = iq$ | Ch. 8 |
| $\theta_0,\ \theta_L$ | chiral wall angles at $x=0$, $x=L$ | Ch. 8 |
| $\Delta = \theta_L - \theta_0$ | wall mismatch (the charge knob) | Ch. 12 |
| $\Sigma = (\theta_0+\theta_L)/2$ | wall mean (a mass-renormalization knob) | Ch. 12 |
| $\hat Q$ | symmetrically ordered charge operator | Ch. 8 |
| $\eta = \sum_k \operatorname{sgn} E_k$ | spectral asymmetry (heat-kernel regulated) | Ch. 9 |
| $Q_{\text{vac}} = -\eta/2$ | vacuum charge | Ch. 9 |
| $\alpha_{kn}^{ss'}$ | Bogoliubov coefficients between sectors | Ch. 8 / App. A |
| $w_k = \int_0^L \lvert\psi_k'\rvert^2 dx$ | weight of new mode $k$ in old region | Ch. 11 |
| $\Delta Q_{\text{cross}},\ \Delta Q_{\text{same}}$ | truncated Bogoliubov charge forms | Ch. 11 |
| $L^*$ | crossing size: $\tanh(mL^*) = -\cos(\Delta/2)/\cos\Sigma$ | Ch. 12 |
| $\nu$ | crossing density per expansion step | Ch. 14 |
| $\varepsilon_{CP}$ | sign-bias of the wall-angle distribution | Ch. 14 |
| $\eta_B$ | baryon-to-photon ratio | Ch. 13/14 |
| $C_{\text{sph}} = 28/79$, $g_*$ | sphaleron conversion; relativistic d.o.f. | Ch. 13 |
| $J_{CP}$, $\mathcal S_{\text{GIM}}$, $\mathcal C_2^{\text{eff}}$ | Jarlskog invariant; GIM suppression; CP-coefficient budget | Ch. 13/14 |

## D.4 Emergent gravity (Part III)

| Symbol | Meaning | Introduced |
|---|---|---|
| $\xi$ | fixed (comoving) coordinate on the unit box | Ch. 17 |
| $g_{ij} = \sum_a e^i_a e^j_a$ | cell metric from edge vectors (triad $e_a$) | Ch. 18 |
| $\hat S^{(ab)}$ | iso-energy shear generators | Ch. 18 |
| $\chi_n(\delta A)$ | fidelity (quantum) metric of eigenmode $n$ | Ch. 19 |
| $\alpha_a = \ln a_a = \ln L_a$ | log scale factors | Ch. 16/20 |
| $p_a = \dot\alpha_a/\sum_b\dot\alpha_b$ | relational Kasner exponents | Ch. 16/20 |
| $\Sigma_K = \sum_a p_a^2$ | anisotropy invariant (1 = vacuum Kasner, 1/3 = FRW) | Ch. 16/20 |
| $u$ | Kasner parameter; BKL map $u\to u-1$, $u\to 1/(u-1)$ | Ch. 16/24 |
| $t(n)$ | ladder hopping law; generator coupling $t(n) = g(n+\tfrac12)$ | Ch. 20 |
| $\kappa = (8\pi G)^{-1}$ | stiffness (constraint normalization) | Ch. 21/22 |
| $w$ | equation-of-state parameter ($\rho \propto a^{-3(1+w)}$) | Ch. 21 |
| $\delta$ (density) | density contrast $\delta\rho/\bar\rho$ — context disambiguates from CP phase | Ch. 21 |
| $E_C(L) = -\pi/24L$ | 1D Casimir energy | Ch. 22 |
| $\Lambda_{\text{UV}} = c/\bar L$ | cell-scale cutoff ($\bar L$ mean cell size; $c$ lattice constant) | Ch. 22 |
| $N_{\text{eff}}$ | spin-weighted species count | Ch. 22 |
| $\sigma$, $\sigma_0$ | Weyl compensator ($\sim 1/L_{\text{cell}}$); its VEV | Ch. 23 |
| $W_\mu$, $F^{(W)}_{\mu\nu}$ | Weyl vector and field strength | Ch. 23 |
| $s$ (proper time) | covariant regulator in $E(s) = \tfrac12\sum\omega_n e^{-s\omega_n}$ — context disambiguates from expansion ratio | Ch. 25 |
| $C(s) = 1/48\pi^2 s^2$ | universal one-loop EH coefficient | Ch. 25 |
| $\kappa_{TT},\ \kappa_{\text{conf}},\ \kappa_{\text{diffeo}}$ | induced gradient couplings by sector | Ch. 25 |
| $h_{ij}$, $P_{TT+} = \mathrm{diag}(0,1,-1)$, $P_{TT\times}$ | metric perturbation; TT polarizations | Ch. 25 |
| $G_{ijkl}$ | DeWitt supermetric | Ch. 16/19 |

## D.5 Reserved-name collisions (explicit rulings)

- $\delta$: CP phase in Part II; density contrast in Ch. 21 (always written with "density contrast" nearby). Never both in one equation.
- $\Sigma$: wall-angle mean (Part II) vs Kasner invariant — the latter is always $\Sigma_K$.
- $s$: expansion ratio (Parts I–II) vs proper-time regulator (Ch. 22/25); flagged at first use in each chapter.
- $\kappa$: Dirac angular quantum number (Ch. 8/10 spherical sections) vs stiffness (Part III); never both in one chapter without explicit note.
- $\eta$: spectral asymmetry (Part II) vs $\eta_B$ (always subscripted) vs adiabatic parameter (always $\eta_{\text{ad}}$).

## D.6 Glossary

Plain-language entries, one breath each, with the home chapter.

**Bogoliubov transformation** — the linear dictionary between two observers' notions of "particle": when the Hamiltonian (here: the geometry) changes, old creation/annihilation operators become mixtures of new ones; the mixing coefficients are mode overlaps, and the cross terms mean the old vacuum contains new pairs. (Ch. 8, App. A.)

**Dirac sea** — the bookkeeping that stabilizes a fermion spectrum unbounded below: declare all negative-energy levels filled and call that the vacuum; particles are filled positive levels, antiparticles are holes in the sea. Equivalent to normal ordering — except for one physical residue, the vacuum charge. (Ch. 8.6–8.7.)

**Spectral asymmetry (η)** — the regulated count of how lopsided a spectrum is between positive and negative energies, $\eta = \sum_k\operatorname{sgn}E_k$. The vacuum's charge is $-\eta/2$; η is flat under smooth deformation and jumps by ±2 when a level crosses zero. (Ch. 9, App. B.2.)

**Fractional fermion number** — the established fact that a vacuum can carry non-integer charge (½ on a Jackiw–Rebbi kink; $-\theta/2\pi$ per chiral wall), with no particle anywhere: the sharp physical reading of $Q_{\text{vac}} = -\eta/2$. (Ch. 9.2, 12.3.)

**Spectral flow** — the net number of eigenvalues crossing zero (with sign) as a parameter is deformed; an integer, robust against everything smooth. In this thesis: net charge production *is* spectral flow. (Ch. 9, Cor. 9.2.)

**Anomaly inflow (Callan–Harvey)** — the mechanism by which charge apparently created on a defect (wall, string) is supplied by a topological current in the bulk; the defect's Dirac operator's spectral flow does the counting. The pump of Ch. 12 is its moving-boundary instance. (Ch. 12.6, 26.1.)

**Chiral bag** — a fermion confined by a wall that reflects with a chiral phase $\theta$ ($-i\,n\!\cdot\!\gamma\,e^{i\theta\gamma_5}\psi = \psi$); MIT bag = $\theta = 0$. The wall angle is physical (chiral symmetry is broken by confinement), and its mismatch between walls is the only charge knob. (Ch. 8.4, 10, 12.)

**Sphaleron** — the saddle-point of electroweak field space between topologically distinct vacua; hot-plasma transitions over it violate $B{+}L$ rapidly above the crossover, converting any charge excess into baryons at the fixed exchange rate $28/79$. (Ch. 13.3.)

**Jarlskog invariant ($J_{CP}$)** — the unique reparametrization-invariant measure of CKM CP violation, $\approx 3\times10^{-5}$; every CP-odd CKM observable is proportional to it. (Ch. 13.5.)

**GIM mechanism** — the unitarity-enforced cancellation that makes flavor-summed CP-odd observables vanish with quark mass *differences*; the structural reason Standard-Model baryogenesis starves, and the budget line this framework must evade at its walls or fail. (Ch. 13.5, 14.4.)

**Vielbein / triad** — the factorization of a metric into frame vectors, $g_{ij} = e^a_ie^a_j$; the variable spinors require, postulated in standard gravity, *derived* here as cell edge vectors. (Ch. 16.9, 18.)

**Weyl rescaling** — multiplying the metric by a position-dependent scale factor, $g \to e^{2\epsilon}g$; the box spectrum's homogeneity is its global version, and gauging it is Ch. 23's road. (Ch. 17.5, 23.)

**Compensator** — a scalar of Weyl weight −1 whose constant gauge-fixed value breaks scale symmetry and sets the Planck mass ($M_P^2 = \sigma_0^2/6$); identified here with the inverse cell size. (Ch. 23.3–23.4.)

**ADM split / lapse / shift** — the slicing of spacetime into space-at-an-instant: the spatial metric is dynamical; the lapse (proper time per slice step) and shift (coordinate slide between slices) are gauge choices of the slicing. (Ch. 16.2.)

**Hamiltonian constraint** — the lapse's equation of motion: total energy of geometry plus matter vanishes, geometry entering negatively along the volume direction; for FRW it *is* the Friedmann equation. (Ch. 16.3, 21.)

**DeWitt supermetric** — the natural inner product on metric velocities, with indefinite signature: the conformal (volume) direction is timelike (a clock), the five shape directions spacelike (they propagate). Derived kinematically in Ch. 19, with its sign measured in Ch. 22/25. (Ch. 16.4.)

**Kasner circle** — the locus $\sum p_a = \sum p_a^2 = 1$ of vacuum Bianchi-I exponents; every point but one has a contracting axis. The generator-coupled ladder lands on it to four digits. (Ch. 16.6, 20.)

**BKL map** — the chaotic recursion $u \to u - 1$, $u \to 1/(u-1)$ governing the sequence of Kasner epochs in a collapsing anisotropic universe; conjugate to the continued-fraction map. Reproduced by the cell model to $6\times10^{-4}$. (Ch. 16.7, 24.)

**Separate-universe approximation** — treating each long-wavelength region as its own homogeneous universe with its own parameters; exact at leading order in gradients, native to a cell model, and the derivation route of the growth equation. (Ch. 21.3.)

**Sakharov induced gravity** — the 1967 observation that integrating out matter on a curved background *generates* an Einstein–Hilbert term ($G^{-1} \sim N\Lambda^2$) whether or not gravity was postulated: gravity as vacuum elasticity. The cell scale supplies the cutoff here, predicting $\bar L \sim \ell_P$. (Ch. 22.)

**Heat kernel / zeta regularization** — the two standard ways to assign finite values to divergent spectral sums (damp-and-expand vs analytically continue); they agree on the universal constants and differ only in how visibly they carry the bulk (cosmological) divergence. (App. B.1.)

**Generalized free field** — a field whose correlators are Gaussian but whose two-point function fails to come from a local canonical theory; the mode-volume normalization of Ch. 15 produces one, with smooth equal-time commutators violating **microcausality** — the requirement that spacelike-separated observables commute, the licence for any local interacting extension. (Ch. 15.4.)

## D.7 Figure index

Generated programmatically: `run_all.py --report` writes the figure-to-script-to-chapter table into `data/provenance.json`; the file list with status lives in `coverage_matrix.md`.
