"""
Sim 4 - dynamics after a local step in a 1D chain with shared, massive walls.
Cells i=0..M-1 each hold one occupant (mode n_i), adiabatic pressure p_i = pi^2 n_i^2/(m L_i^3).
Wall j (between cell j and j+1): M_w x_j'' = p_j - p_{j+1}.  Outer walls fixed (far away).
P1: occupant of centre cell promoted n->n+1 at fixed walls, then released.
P2: thesis-literal iso-energy step: centre cell jumps to L(n+1)/n (walls moved outward by delta/2
    each; neighbours contracted), occupant n->n+1; then released.
"""
import numpy as np
m, Mw, n0, M = 1.0, 200.0, 10, 601
c = M//2
def run(protocol, T=60.0, dt=2e-3, nsave=300):
    n = np.full(M, n0, float)
    x = np.arange(1, M, dtype=float)              # interior walls; cell i spans [x_{i-1}, x_i], x_{-1}=0, x_{M-1}=M
    v = np.zeros_like(x)
    L = lambda x: np.diff(np.concatenate([[0.0], x, [float(M)]]))
    E0 = (np.pi**2*n**2/(2*m*L(x)**2)).sum()
    n[c] += 1
    if protocol == "P2":
        d = 1.0/n0
        x[c-1] -= d/2; x[c] += d/2
    Ecell = lambda x: (np.pi**2*n**2/(2*m*L(x)**2))
    E_inj = Ecell(x).sum() - E0
    p = lambda x: np.pi**2*n**2/(m*L(x)**3)
    acc = lambda x: (p(x)[:-1] - p(x)[1:])/Mw
    a = acc(x); hist = []; ts = []
    steps = int(T/dt); every = steps//nsave
    for k in range(steps+1):
        if k % every == 0:
            hist.append(L(x).copy()); ts.append(k*dt)
            if k == 0: Etot0 = Ecell(x).sum() + 0.5*Mw*(v**2).sum()
        x += v*dt + 0.5*a*dt*dt; an = acc(x); v += 0.5*(a+an)*dt; a = an
    Etot1 = Ecell(x).sum() + 0.5*Mw*(v**2).sum()
    return np.array(ts), np.array(hist), E_inj, Etot0, Etot1, np.abs(v).max()

vp = np.pi*n0/1.0/m; cs = np.sqrt(3*m/Mw)*vp
print(f"particle speed v_p={vp:.3f}; predicted sound speed c_s = sqrt(3m/Mw) v_p = {cs:.4f} cells/time")
out = {}
for P in ("P1", "P2"):
    ts, H, Einj, E0, E1, vmax = run(P)
    dL = H - 1.0
    # front: outermost cell with |dL| > 1e-4 on the right
    front = np.array([np.max(np.where(np.abs(h[c:]) > 2e-4)[0]) if (np.abs(h[c:]) > 2e-4).any() else 0 for h in dL])
    sel = (ts > 5) & (ts < 55)
    vf = np.polyfit(ts[sel], front[sel], 1)[0]
    late = ts > 40
    Lc_late = H[late, c].mean()
    print(f"{P}: energy injected by the step {Einj:+.4f} (cell energy {np.pi**2*n0**2/2:.3f}); "
          f"energy drift {E1-E0:+.2e}; max wall speed/v_p {vmax/vp:.3e}")
    print(f"    front speed {vf:.4f} (c_s={cs:.4f});  centre cell late-time mean L = {Lc_late:.5f}, "
          f"predicted ((n+1)/n)^(2/3) = {((n0+1)/n0)**(2/3):.5f}; iso-energy rule (n+1)/n = {(n0+1)/n0:.5f}")
    far = (np.arange(M) > c+5) & (np.arange(M) < c+100)
    print(f"    after the front passes, cells at X+dX return to L: mean |L-1| in cells c+5..c+100 at t=55: {np.abs(H[-1][far]-1).mean():.2e}")
    out[P] = (ts, H)
np.savez("figs/sim4.npz", ts=out["P1"][0], H1=out["P1"][1], H2=out["P2"][1], c=c, cs=cs, n0=n0)
