import numpy as np
from concurrent.futures import ProcessPoolExecutor
import q2_production as P
from q2_analyse import speeds
def one(args):
    kind, M, T, seed, kick = args
    P.kick = kick
    return P.one((kind, M, T, seed))
if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=30) as ex:
        for kind, mu in ((1, 1.0), (0, 2.0)):
            M = mu*16 if kind == 0 else mu*32
            pred = np.sqrt(1/(1+mu)) if kind else np.sqrt(3/(1+mu))
            T = 110/pred
            for kick in (1.0, 0.3):
                R = 240 if kick == 1.0 else 720
                d = np.array(list(ex.map(one, [(kind, M, T, 5000+s, kick) for s in range(R)])))
                c0 = P.Nc//2; right = d[:, :, c0+1:]; left = -d[:, :, c0::-1][:, :, :right.shape[2]]
                sym = np.concatenate([right, left]); r = dict(T=T, mean=sym.mean(0), sem=sym.std(0)/np.sqrt(len(sym)), arrive=np.full(right.shape[2], np.nan))
                vs = speeds(r)[0]
                print(f"kind={kind} mu={mu} kick={kick}: pulse speed {vs:.3f}  prediction {pred:.3f}  ratio {vs/pred:.3f}", flush=True)
