"""fw1atomo_x1.py - exploration AFTER the pre-registration (declared): at theta = 0,5 and h = 0,6 / 0,8 the Ritz of P1c
found no exact level (no residual < 1e-4) near the 1s.  Is the 1s still bound?  Sublattice 0, spin-up continuum 1s,
1600 periods: survival lines and the weight left inside r < 6 lambda_C (leak into the box), plus h = 0,4 as control.
Also the free-rule gaps of FW1 at theta = 0,5 (lowest energy away from Gamma, distance from the pi edge) for each h.
Run: OMP_NUM_THREADS=1 timeout 590 python lab/fw1atomo_x1.py <h>; python lab/fw1atomo_x1.py sum"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_core as fc

RES = os.path.join(fc.SCR, "x1")
os.makedirs(RES, exist_ok=True)


def free_gaps(h, theta=0.5, n=24):
    import senzadoppi_core as s, senzadoppi_r5 as r5
    k = (np.arange(n) + 0.5) / n * 2 * np.pi - np.pi
    K = np.stack([q.ravel() for q in np.meshgrid(k, k, k, indexing="ij")], -1)
    q = s.quasi(r5.FW1(h).bloch(K, theta), 2)
    G = np.zeros(len(K), bool)
    for Gc in [(0, 0, 0)] + [tuple(x) for x in np.array(np.meshgrid(*[[0, np.pi]] * 3)).reshape(3, -1).T]:
        d = np.sqrt(np.sum(((K - np.array(Gc) + np.pi) % (2 * np.pi) - np.pi)**2, 1))
        G |= d < 0.45
    E0 = float(np.min(np.abs(q[~G])))
    Epi = float(np.min(np.pi / 2 - np.abs(q)))       # period T = 2 ticks: edge at pi/T = pi/2 per tick
    return E0, Epi


def run(h):
    L = fc.Lattice(40, 0.5, "angolo")
    R = L.rule("FW1", h)
    s = L.dirac1s() * (L.cls == 0)[None]
    psi0 = s / np.linalg.norm(s)
    psi = psi0.copy()
    c = np.zeros(1601, complex); c[0] = 1
    inner = []
    m = L.r < 6
    for p in range(1, 1601):
        psi = R.period(psi)
        c[p] = np.vdot(psi0, psi)
        if p % 50 == 0:
            inner.append(float(np.sum(np.abs(psi[:, m])**2)))
    E0, Epi = free_gaps(h)
    out = {"h": h, "lines": sorted(fc.lines(c, 1.0, thr=0.01), key=lambda x: -x[1])[:6], "inner_every_50": inner,
           "free_E0_per_tick": E0, "free_Epi_per_tick": Epi}
    json.dump(out, open(os.path.join(RES, f"h{h}.json"), "w"))
    print(json.dumps(out))


if __name__ == "__main__":
    if sys.argv[1] == "sum":
        res = {f[:-5]: json.load(open(os.path.join(RES, f))) for f in sorted(os.listdir(RES))}
        json.dump(res, open(os.path.join(os.path.dirname(__file__), "out", "fw1atomo_x1.json"), "w"), indent=1)
        for k, v in res.items():
            print(k, v["lines"][:3], v["inner_every_50"][::8], v["free_E0_per_tick"], v["free_Epi_per_tick"])
    else:
        run(float(sys.argv[1]))
