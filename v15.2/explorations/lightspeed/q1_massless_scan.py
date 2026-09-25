import numpy as np
from q1_massless import analyse
us = np.concatenate([np.linspace(0.05, 0.9, 18), [0.95, 0.98, 0.995, 1.0]])
out = {}
for n in (3, 10, 30):
    out[n] = np.array([(u,)+analyse(n, u, ds=1e-5 if n == 30 else 2e-5)[:2] for u in us])
    print(n, "u=1:", np.round(out[n][-1], 5), " u=0.5:", np.round(out[n][np.argmin(abs(us-0.5))], 5))
np.save("figs/q1_massless.npy", out, allow_pickle=True)
