"""
Q1c - if each iso-energy step moves the wall by dL = lambda/2 in t0 = dL/c (massless matter),
a cell's wall can move at most at c, so the steps of the N occupants must not overlap:
    total step rate  sum_a Gamma_a <= 1/t0  (for equal energies)
Monte-Carlo: Poisson step proposals at rate Gamma per occupant; a proposal arriving while the
wall is already moving is either (i) blocked (causal rule) or (ii) allowed to superpose
(velocities add, as a naive rate model would). Record mean wall speed dL/dt / c.
"""
import numpy as np
rng = np.random.default_rng(4)
def run(N, Gamma, t0=1.0, T=4000.0, rule="block"):
    t, busy_until, moved, active = 0.0, 0.0, 0.0, []
    R = N*Gamma
    while t < T:
        t += rng.exponential(1/R)
        if rule == "block":
            if t >= busy_until:
                busy_until = t + t0; moved += 1.0
        else:
            moved += 1.0
    return moved*1.0/T   # each step moves dL = c t0 = 1
for N in (1, 4, 16, 64):
    for G in (0.01, 0.05, 0.2, 1.0):
        b = run(N, G); f = run(N, G, rule="add")
        x = N*G
        print(f"N={N:3d} Gamma*t0={G:5}: wall speed/c causal {b:.3f} (theory x/(1+x)={x/(1+x):.3f}), naive {f:.3f} (=N Gamma t0 {x:.3f})")
