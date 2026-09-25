import numpy as np
from cellmedium import Cell, hessian_F
from itertools import product
for K in (44, 56):
    c = Cell(K=K); h = 2e-3; I = np.eye(2)
    E = lambda i, j: np.outer(np.eye(2)[i], np.eye(2)[j])
    Wc = lambda F: c.W(F)
    C1212 = (Wc(I+2*h*E(0,1)) - 2*Wc(I) + Wc(I-2*h*E(0,1)))/(4*h*h)
    print(f"K={K}: C1212 = {C1212:.3f}")
