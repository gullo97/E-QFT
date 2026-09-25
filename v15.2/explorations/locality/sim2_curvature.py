"""
Sim 2 - curvature carried by a single-cell step when cells are glued abstractly
(edge lengths fundamental, no ambient embedding).  Exact 2D Gaussian curvature.
  directional step: g = diag(e^{2a}, 1),  K = -e^{-a} d_y^2 e^{a}
  isotropic step:   g = e^{2a} delta,     K = -e^{-2a} lap(a)
a(x,y) = ln(s) * chi_w(x,y), chi_w a smoothed indicator of the unit cell.
"""
import numpy as np
from vstyle import *

def chi(X, Y, w):
    f = lambda u: 0.25*(1+np.tanh((u+0.5)/w))*(1-np.tanh((u-0.5)/w))
    return f(X)*f(Y)

def fields(s, w, N=801, R=2.5):
    x = np.linspace(-R, R, N); h = x[1]-x[0]
    X, Y = np.meshgrid(x, x, indexing="xy")          # rows = y, cols = x
    a = np.log(s)*chi(X, Y, w)
    ea = np.exp(a)
    d2y = lambda f: np.gradient(np.gradient(f, h, axis=0), h, axis=0)
    d2x = lambda f: np.gradient(np.gradient(f, h, axis=1), h, axis=1)
    Kdir = -d2y(ea)/ea
    Kiso = -(d2x(a)+d2y(a))*np.exp(-2*a)
    return x, h, a, Kdir, Kiso

s = 4/3                                              # n = 3 -> 4 step, (n+1)/n
for w in (0.02, 0.04, 0.08):
    x, h, a, Kd, Ki = fields(s, w)
    ea = np.exp(a)
    tot_d = (Kd*ea).sum()*h*h                        # sqrt(g) = e^a
    tot_i = (Ki*np.exp(2*a)).sum()*h*h
    j0 = np.argmin(abs(x)); up = x > 0
    dip = np.sum((x[up]-0.5)*Kd[up, j0]*ea[up, j0])*h   # across top face, per unit length
    faceK = abs(Kd[:, j0]).max()
    # isotropic: curvature on x-faces too
    print(f"w={w}: total curvature dir {tot_d:+.2e}, iso {tot_i:+.2e}; "
          f"dipole across top face {dip:+.5f} (predicted -(s-1) = {-(s-1):+.5f}); max|K| {faceK:.1f}"
          f"; max|K| on x-face (dir) {abs(Kd[j0, :]).max():.2e}")
np.savez("figs/sim2.npz", s=s)
