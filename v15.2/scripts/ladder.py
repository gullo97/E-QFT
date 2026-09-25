"""
ladder.py — shared module for Part III quantum ladder simulations (App. C.3).

Sparse tridiagonal hopping Hamiltonians on the geometry ladder; wavepacket
evolution by Krylov `expm_multiply`; observables <n>, Var(n), alpha = ln<n>;
relational exponents p_a by centered differences against the common time grid
(clock-free by construction: any reparametrization cancels in p_a).

Used by `ch20_kasner_sim.py` and `ch21_friedmann.py` / `ch21_collapse.py`
(the only permitted inter-script imports besides `bag1d`, contract C.1.2).
"""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply


def hop_hamiltonian(N, tfun):
    """H = sum_n t(n) (|n+1><n| + h.c.) on rungs n = 1..N (bond n -> n+1)."""
    n = np.arange(1, N + 1)
    t_n = tfun(n[:-1])
    return diags([t_n, t_n], [-1, 1], format="csc"), n


def gaussian_packet(n, n0, sigma, theta):
    """Gaussian in n with quasimomentum theta imprinted as a phase ramp."""
    psi = np.exp(-(n - n0) ** 2 / (4 * sigma**2) + 1j * theta * n)
    return psi / np.linalg.norm(psi)


def evolve_chain(N, tfun, n0, sigma, theta, times, snap_stride=4):
    """Evolve a Gaussian wavepacket; return <n>(t), Var(t), snapshots, n-grid."""
    H, n = hop_hamiltonian(N, tfun)
    psi = gaussian_packet(n, n0, sigma, theta)
    means, vars_, snaps = [], [], []
    tprev = 0.0
    for tt in times:
        if tt > tprev:
            psi = expm_multiply((-1j * (tt - tprev)) * H, psi)
            tprev = tt
        P = np.abs(psi) ** 2
        mu = float(P @ n)
        means.append(mu)
        vars_.append(float(P @ n**2) - mu**2)
        snaps.append(P[::snap_stride].astype(np.float32))
    return np.array(means), np.array(vars_), np.array(snaps), n[::snap_stride]


def relational(times, means_list):
    """alpha_a = ln<n_a>, relational exponents p_a = alphadot_a / sum_b, Sum p^2."""
    al = np.array([np.log(m) for m in means_list])
    dal = np.array([np.gradient(a, times) for a in al])
    tot = dal.sum(axis=0)
    p = dal / tot
    return al, p, (p**2).sum(axis=0)
