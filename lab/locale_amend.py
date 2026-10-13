"""locale_amend.py - REGOLA LOCALE: declared AMENDMENTS after the pre-registration (written after locale_r1 ran).
Nothing in lab/out/locale_prereg.json is changed; the FAIL of T1b stays a FAIL in the table.
 A1 (tool control, added after r1): the degree integral of locale_core.su2_degree on a NON-local map of known degree,
    a(k) = (m + sum cos k_j, sin k_x, sin k_y, sin k_z)/|.|: degree 1 (m = -2), -2 (m = 0), 1 (m = 2), 0 (m = 4).
    Without this control the 14 zeros of T2 would not prove that the integral can see a non-zero degree.
 A2 (T1b redone): r1 used a 32^3 periodic box; after 32 ticks the packet (rms width ~16 cells) wraps around the box
    and the centroid of a periodic box is ill-defined (the FAIL measures the box, not the copies). Same packet, same
    32 ticks, same 1% criterion, box 96^3 (packet far from the edges).
Writes lab/out/locale_amend.json.  Run: OMP_NUM_THREADS=1 timeout 590 python lab/locale_amend.py"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from locale_core import *
from locale_cands import fermion

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()
res, verdict = {}, []


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


class Ctrl:
    def __init__(self, m):
        self.m = m

    def bloch(self, K):
        a0 = self.m + np.cos(K).sum(1); a = np.sin(K)
        n = np.sqrt(a0 ** 2 + (a ** 2).sum(1))
        return (a0 / n)[:, None, None] * I2[None] + 1j * np.einsum("ni,ijk->njk", a / n[:, None], PAULI)


res["A1"] = {}
for m, exp in ((-2.0, 1), (0.0, -2), (2.0, 1), (4.0, 0)):
    d = su2_degree(Ctrl(m), 40)
    res["A1"][str(m)] = d
    check(f"A1 control m={m}", abs(d - exp) <= 0.02, f"degree {d:+.6f} (expected {exp})")

# A2
n = 96
rng = np.random.default_rng(SEED)        # same seed and draw order as locale_r1 is not needed: new spinor declared
w = fermion("F0", 3, 0.1)
R = RealWalk(w, n, 3)
idx = np.indices((n,) * 3).astype(float)
c0 = n / 2
g = np.exp(-((idx[0] - c0) ** 2 + (idx[1] - c0) ** 2 + (idx[2] - c0) ** 2) / (4 * 3.0 ** 2) + 1j * 0.3 * idx[0])
s = rng.normal(size=4) + 1j * rng.normal(size=4); s /= np.linalg.norm(s)
full = s[:, None, None, None] * g[None]
mask = ((idx[0] % 2 == 0) & (idx[1] % 2 == 0) & (idx[2] % 2 == 0))
one = full * mask[None]


def moments(psi):
    rho = np.sum(np.abs(psi) ** 2, axis=0)
    m = rho.sum()
    c = np.array([np.sum(rho * idx[a]) for a in range(3)]) / m
    wd = np.sqrt(np.sum(rho * ((idx[0] - c[0]) ** 2 + (idx[1] - c[1]) ** 2 + (idx[2] - c[2]) ** 2)) / m)
    return c, wd


traj = []
for t in range(16):
    full = R.period(full); one = R.period(one)
    cf, wf = moments(full); co, wo = moments(one)
    traj.append([float(cf[0]), float(co[0]), float(wf), float(wo)])
outside = float(np.max(np.abs(one * (~mask)[None])))
dc = max(abs(a - b) for a, b, _, _ in traj)
dw = max(abs(c - d) for _, _, c, d in traj)
wref = traj[-1][2]
edge = float(np.sum(np.abs(full[:, :4]) ** 2 + np.abs(full[:, -4:]) ** 2) / np.sum(np.abs(full) ** 2))
res["A2"] = dict(outside_max=outside, centroid_max_diff=dc, width_max_diff=dw, width=wref, traj=traj, edge_fraction=edge)
check("A2 one coset stays on its coset (96^3)", outside == 0.0, f"max |psi| outside = {outside:.1e}")
check("A2 one coset = same physics (96^3)", max(dc, dw) <= 0.01 * wref,
      f"centroid diff {dc:.2e}, width diff {dw:.2e} cells (width {wref:.2f}; 1% = {0.01 * wref:.3f}); "
      f"probability within 4 cells of the box edge {edge:.1e}")
res["verdict"] = verdict
res["seconds"] = time.time() - t0
with open(os.path.join(OUT, "locale_amend.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
print(f"=== {sum(v for _, v in verdict)} PASS / {sum(not v for _, v in verdict)} FAIL ({time.time() - t0:.0f} s)")
