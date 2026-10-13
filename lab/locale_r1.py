"""locale_r1.py - REGOLA LOCALE, production run 1: general exams G1-G4 for every candidate, T1 (sheets) and
T2 (3-winding / chiralities). Pre-registration: lab/out/locale_prereg.json. Writes lab/out/locale_r1.json.
Run: OMP_NUM_THREADS=1 timeout 590 python lab/locale_r1.py"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from locale_core import *
from locale_cands import fermion, light, LIGHTS, FERMIONS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()
rng = np.random.default_rng(SEED)
res, verdict = {}, []


def check(name, ok, msg):
    verdict.append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def radius(w):
    r_f = max(np.max(np.abs(iv)) if iv is not None else 0 for iv in w.ivec)
    # period radius: largest |displacement| a component can accumulate = sum of |v| over factors (upper bound)
    r_p = sum(np.linalg.norm(w.B @ iv) * (2 if f[0] == "link" else 1) for f, iv in zip(w.pre, w.ivec) if iv is not None
              and f[0] == "shift")
    return float(r_f), float(r_p)


# ---------------------------------------------------------------- G1-G4
print("=== G1-G4 unitarity, locality, charge, Bloch == running walk ===", flush=True)
cands = [(n, fermion(n, 3)) for n in FERMIONS] + [(n + "_2D", fermion(n, 2)) for n in ("F0", "F2")] + \
        [(n, light(n)[0]) for n in LIGHTS]
res["G"] = {}
for name, w in cands:
    dim = w.dim
    K = rng.uniform(-PI, PI, (2000, 3)); K[:, dim:] = 0
    W = w.bloch(K)
    u = float(np.max(np.abs(np.einsum("nba,nbc->nac", W.conj(), W) - np.eye(w.n)[None])))
    n = 12 if dim == 3 else 32
    V = rng.uniform(0, 2 * PI, (n,) * dim)
    R = RealWalk(w, n, dim, V=V)
    psi = rng.normal(size=(w.n,) + (n,) * dim) + 1j * rng.normal(size=(w.n,) + (n,) * dim)
    N0 = np.sum(np.abs(psi) ** 2)
    for _ in range(10):
        psi = R.period(psi)
    dN = float(abs(np.sum(np.abs(psi) ** 2) - N0) / N0)
    # plane wave
    R0 = RealWalk(w, n, dim)
    klat = 2 * PI * rng.integers(0, n, dim) / n
    idx = np.indices((n,) * dim)
    s0 = rng.normal(size=w.n) + 1j * rng.normal(size=w.n)
    pw = s0.reshape((-1,) + (1,) * dim) * np.exp(1j * np.tensordot(klat, idx, axes=(0, 0)))[None]
    out = R0.period(pw)
    Kr = np.linalg.solve(w.B.T, np.concatenate([klat, np.zeros(3 - dim)]))
    pred = w.bloch(Kr[None])[0] @ s0
    b = float(np.max(np.abs(out[(slice(None),) + (0,) * dim] - pred)))
    rf, rp = radius(w)
    res["G"][name] = dict(unitarity=u, charge=dN, bloch=b, radius_factor=rf, radius_period_bound=rp, T=w.T,
                          n_factors=len(w.pre), n_shift=sum(f[0] == "shift" for f in w.pre),
                          n_link=sum(f[0] == "link" for f in w.pre), n_coin=sum(f[0] == "coin" for f in w.pre))
    check(f"G1 {name}", u <= 1e-12, f"||W^dag W - 1|| = {u:.1e}")
    check(f"G3 {name}", dN <= 1e-12, f"|dN|/N after 10 periods with random on-site phases = {dN:.1e}")
    check(f"G4 {name}", b <= 1e-12, f"running walk vs Bloch = {b:.1e}")
    print(f"     G2 {name}: factor radius {rf:.0f} cells, period T = {w.T:.3f}, factors {len(w.pre)} "
          f"(shift {res['G'][name]['n_shift']}, link {res['G'][name]['n_link']}, coin {res['G'][name]['n_coin']})", flush=True)

# ---------------------------------------------------------------- T1 sheets
print("=== T1 sheets: the frozen rule is 2^d independent copies on the lattice 2Z^d ===", flush=True)
res["T1"] = {}
for dim in (3, 2):
    w = fermion("F0", dim, 0.1)
    K = rng.uniform(-PI, PI, (2000, 3)); K[:, dim:] = 0
    d = 0.0
    for j in range(dim):
        e = np.zeros(3); e[j] = PI
        d = max(d, float(np.max(np.abs(w.bloch(K + e) - w.bloch(K)))))
    res["T1"][f"periodicity_{dim}D"] = d
    check(f"T1a {dim}-D W(k + pi e_j) = W(k)", d <= 1e-13, f"max difference {d:.1e} -> {2 ** dim} copies")
# real space
n = 32
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
for t in range(16):                       # 16 periods = 32 ticks
    full = R.period(full)
    one = R.period(one)
    cf, wf = moments(full); co, wo = moments(one)
    traj.append([float(cf[0]), float(co[0]), wf, wo])
outside = float(np.max(np.abs(one * (~mask)[None])))
dc = max(abs(a - b) for a, b, _, _ in traj)
dw = max(abs(c - d) for _, _, c, d in traj)
wref = traj[-1][2]
res["T1"]["realspace"] = dict(outside_max=outside, centroid_max_diff=dc, width_max_diff=dw, width=wref, traj=traj)
check("T1b one coset stays on its coset", outside == 0.0, f"max |psi| outside the coset after 32 ticks = {outside:.1e}")
check("T1b one coset = same physics", max(dc, dw) <= 0.01 * wref,
      f"centroid diff {dc:.2e}, width diff {dw:.2e} cells (width {wref:.2f}; threshold 1% = {0.01 * wref:.3f})")

# ---------------------------------------------------------------- T2 degree
print("=== T2 3-winding (degree) of Weyl walks: predicted 0 for every split-step walk ===", flush=True)
res["T2"] = {}


def random_weyl(seed):
    r = np.random.default_rng(seed)
    vecs = [np.array(v, float) for v in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (0, 1, 1),
                                          (1, 0, 1), (1, 1, 1), (1, -1, 1), (2, 0, 0), (0, 2, 1)]]
    fac = []
    for _ in range(r.integers(4, 9)):
        t = r.integers(0, 3)
        if t == 0:
            v = vecs[r.integers(len(vecs))] * r.choice([1, -1])
            n = v / np.linalg.norm(v)
            if r.random() < 0.5:   # spin axis not along v (still integer spectrum)
                n = r.normal(size=3); n /= np.linalg.norm(n)
            fac.append(("shift", v, sdot(n, 0.5)))
        elif t == 1:
            a = r.normal(size=3)
            fac.append(("coin", np.cos(np.linalg.norm(a)) * I2 - 1j * np.sin(np.linalg.norm(a)) * sdot(a / np.linalg.norm(a), 0.5)))
        else:
            a = r.normal(size=3)
            C = np.cos(np.linalg.norm(a)) * I2 - 1j * np.sin(np.linalg.norm(a)) * sdot(a / np.linalg.norm(a), 0.5)
            fac.append(("link", vecs[r.integers(len(vecs))], np.diag([0, 1]).astype(complex), C))
    return SplitWalk(fac, 2, T_phys=1.0, name=f"random{seed}")


weyls = [("U_xyz", make_walk("weyl", [CUBE_AXES])), ("U_zyx", make_walk("weyl", [CUBE_AXES[::-1]])),
         ("W_frozen", make_walk("weyl", zigzag_ticks(3))), ("S4_weyl", composed_walk("weyl", S4_SEQ))] + \
        [(f"random{s}", random_weyl(s)) for s in range(10)]
for name, w in weyls:
    ngrid = 56 if name == "S4_weyl" else 40
    dg = su2_degree(w, ngrid)
    res["T2"][name] = dg
    check(f"T2 degree {name}", abs(dg) <= 0.02, f"3-winding = {dg:+.4f}")
    if time.time() - t0 > 500:
        print("time budget: stopping T2 early", flush=True)
        break
nodes = find_nodes(make_walk("weyl", zigzag_ticks(3)), 16)
s_plus = sum(c for k, s, c, _ in nodes if s > 0)
s_minus = sum(c for k, s, c, _ in nodes if s < 0)
res["T2"]["frozen_W_nodes"] = [dict(k_over_pi=(k / PI).round(4).tolist(), U=s, chirality=c) for k, s, c, _ in nodes]
check("T2 chirality sums of the frozen W", s_plus == 0 and s_minus == 0,
      f"{sum(1 for _, s, _, _ in nodes if s > 0)} nodes at quasi-energy 0 (sum of chiralities {s_plus}), "
      f"{sum(1 for _, s, _, _ in nodes if s < 0)} at pi (sum {s_minus})")

res["verdict"] = verdict
res["seconds"] = time.time() - t0
with open(os.path.join(OUT, "locale_r1.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
npass = sum(v for _, v in verdict)
print(f"=== {npass} PASS / {len(verdict) - npass} FAIL ({time.time() - t0:.0f} s)", flush=True)
