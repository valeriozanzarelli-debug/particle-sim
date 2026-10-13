"""locale_s4_run.py - PHASE 2 production run of REGOLA LOCALE, exactly the exams of lab/out/locale_s4_prereg.json
(P3 hydrogen, P4 Schwinger, P5 folding, P6 causality / precursors, P7 variants).  Checkpointed: every finished item is
written to lab/out/locale_s4_run.json; the script stops after ~500 s with exit code 3 ("run me again").
Run:  OMP_NUM_THREADS=1 timeout 590 python lab/locale_s4_run.py   (repeat until exit code 0)."""
import os, sys, json, time, itertools
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.linalg import expm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import (composed_walk, make_walk, zigzag_ticks, RealWalk, S4_SEQ, n_for_tol, dirs_for, fib_sphere)
from locale_s4_core import SEQ_S4, SEQ_FROZEN, seq_c, cone, dirac1_bloch, RadialS4
from zz_core import dirac_formula, ALPHA

PI = np.pi
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_s4_run.json")
T0 = time.time()
BUDGET = 500.0
R = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    with open(OUT + ".tmp", "w") as f:
        json.dump(R, f, indent=1)
    os.replace(OUT + ".tmp", OUT)


def todo(key):
    if key in R:
        return False
    if time.time() - T0 > BUDGET:
        save(); print("budget reached, rerun"); sys.exit(3)
    print(f"[{time.time() - T0:6.1f}s] {key}", flush=True)
    return True


SEQS = {"S4": SEQ_S4, "frozen": SEQ_FROZEN}

# ------------------------------------------------------------------ P3 hydrogen (radial)
LEVELS = [("1s", 1, -1), ("2s", 2, -1), ("2p1/2", 2, 1), ("2p3/2", 2, -2)]
for sname, seq in SEQS.items():
    for th in (0.2, 0.1, 0.05):
        for lab, n, kap in LEVELS:
            key = f"P3|{sname}|{th}|{lab}"
            if not todo(key):
                continue
            E = dirac_formula(n, kap)
            rs = RadialS4(kap, th, 30 / ALPHA, seq)
            om, x, res = rs.eigen(E, iters=12)
            R[key] = dict(E_lat=om, E_exact=float(E), rel_binding_err=float(((1 - om) - (1 - E)) / (1 - E)),
                          resid=res, sites=int(rs.N))
            save()

# ------------------------------------------------------------------ P4 Schwinger (Landau-Zener in a constant field)
TH4 = 0.1


def bands(W):
    w, v = np.linalg.eig(W)
    eps = -np.angle(w)
    i = np.argsort(eps)
    return eps[i], v[:, i]                       # [negative, positive] (both within the window |eps| < pi)


def lz_lattice(seq, EEc, kmax):
    E = EEc * TH4 ** 2
    T = float(sum(seq))
    nper = int(round(2 * kmax / (E * T)))
    k = -kmax
    W0, _ = dirac1_bloch(k, TH4, seq)
    _, v0 = bands(W0)
    psi = v0[:, 0].copy()
    A = lambda t: -E * t
    for p in range(nper):
        W, _ = dirac1_bloch(k, TH4, seq, A=A, t0=p * T)
        psi = W @ psi
    kend = k + E * nper * T
    W1, _ = dirac1_bloch(kend, TH4, seq)
    eps1, v1 = bands(W1)
    P = float(abs(np.vdot(v1[:, 1], psi)) ** 2)
    return P, float(kend), nper, eps1.tolist()


def lz_continuum(EEc, kmax, kend):
    """continuum 2-level Dirac i dpsi/dt = (k(t) sz + m sx) psi, same finite sweep (information only)"""
    E = EEc * TH4 ** 2
    tend = (kend + kmax) / E
    n = int(np.ceil(tend / 0.01))
    dt = tend / n
    sz = np.diag([1.0, -1.0]); sx = np.array([[0, 1.0], [1.0, 0]])
    H = lambda k: k * sz + TH4 * sx
    _, v = np.linalg.eigh(H(-kmax))
    psi = v[:, 0].astype(complex)
    for j in range(n):
        kk = -kmax + E * (j + 0.5) * dt
        # 2x2 exact exponential of the midpoint Hamiltonian (2nd order in dt, dt = 0.01 -> error ~1e-6)
        e = np.sqrt(kk ** 2 + TH4 ** 2)
        U = np.cos(e * dt) * np.eye(2) - 1j * np.sin(e * dt) * H(kk) / e
        psi = U @ psi
    _, v = np.linalg.eigh(H(kend))
    return float(abs(np.vdot(v[:, 1], psi)) ** 2)


for sname, seq in SEQS.items():
    for kmax in (0.4, 0.8):
        for EEc in (0.5, 1, 2, 4):
            key = f"P4|{sname}|{kmax}|{EEc}"
            if not todo(key):
                continue
            P, kend, nper, eps1 = lz_lattice(seq, EEc, kmax)
            Ps = float(np.exp(-PI / EEc))
            R[key] = dict(P=P, P_schwinger=Ps, ratio_m1=P / Ps - 1, k_end=kend, periods=nper,
                          P_continuum_same_sweep=lz_continuum(EEc, kmax, kend))
            save()


# ------------------------------------------------------------------ P5 folding
def min_sep(W):
    eps = -np.angle(np.linalg.eigvals(W))
    d = []
    for i in range(len(eps)):
        for j in range(i + 1, len(eps)):
            d.append(abs(np.angle(np.exp(1j * (eps[i] - eps[j])))))
    return min(d)


def ray_crossings(walk, n, kmin, kmax, steps=600):
    """interior local minima of the smallest eigenphase separation along the ray k n, refined by ternary search"""
    ks = np.linspace(kmin, kmax, steps)
    K = ks[:, None] * n[None]
    Ws = walk.bloch(K)
    d = np.array([min_sep(W) for W in Ws])
    out = []
    for i in range(1, steps - 1):
        if d[i] <= d[i - 1] and d[i] <= d[i + 1]:
            a, b = ks[i - 1], ks[i + 1]
            f = lambda k: min_sep(walk.bloch((k * n)[None])[0])
            for _ in range(60):
                m1, m2 = a + (b - a) / 3, b - (b - a) / 3
                if f(m1) < f(m2):
                    b = m2
                else:
                    a = m1
            km = 0.5 * (a + b)
            out.append((float(km), float(f(km))))
    return out


def fold_scan(walk, dirs, Nmin=13.0, kmax_extra=1.25):
    kmax = 2 * PI / Nmin * kmax_extra
    res = []
    for nn in dirs:
        for km, g in ray_crossings(walk, nn, 0.02, kmax):
            res.append(dict(dir=[float(x) for x in nn], k=km, N=2 * PI / km, gap=g))
    bad = [r for r in res if r["N"] >= Nmin and r["gap"] > 1e-6]
    first = min([r for r in res if r["gap"] > 1e-6], key=lambda r: -r["N"], default=None)
    return dict(crossings=res, n_bad=len(bad), bad=bad[:10], first_avoided=first,
                max_gap_protected=max([r["gap"] for r in res if r["gap"] <= 1e-6], default=None))


if todo("P5a|2D"):
    R["P5a|2D"] = fold_scan(composed_walk("light", S4_SEQ, dim=2), dirs_for(2)[::3])
    save()
if todo("P5a|3D"):
    R["P5a|3D"] = fold_scan(composed_walk("light", S4_SEQ, dim=3), fib_sphere(40, half=True))
    save()
for sname, seq in SEQS.items():
    key = f"P5b|{sname}"
    if not todo(key):
        continue
    T = float(sum(seq))
    # first fold: the two bands meet near eps = +-pi (i.e. E T = pi). separation measured mod 2 pi.
    kp = np.sqrt(max((PI / T) ** 2 - TH4 ** 2, 0.0)) if PI / T < PI else PI
    ks = np.linspace(0.7 * kp, min(1.3 * kp, PI), 4001) if T > 1 else np.linspace(0.5 * PI, PI, 4001)
    seps = np.array([min_sep(dirac1_bloch(k, TH4, seq)[0]) for k in ks])
    i = int(np.argmin(seps))
    a, b = ks[max(i - 1, 0)], ks[min(i + 1, len(ks) - 1)]
    f = lambda k: min_sep(dirac1_bloch(k, TH4, seq)[0])
    for _ in range(60):
        m1, m2 = a + (b - a) / 3, b - (b - a) / 3
        if f(m1) < f(m2):
            b = m2
        else:
            a = m1
    kf = 0.5 * (a + b)
    R[key] = dict(T=T, k_fold=float(kf), gap_rad_per_period=float(f(kf)), window_E_over_m=float(PI / T / TH4),
                  E_fold_over_m=float(np.sqrt(kf ** 2 + TH4 ** 2) / TH4))
    save()


# ------------------------------------------------------------------ P6 causality (2-D light, real space)
NB = 384
C0 = NB // 2
x = np.arange(NB) - C0
X, Y = np.meshgrid(x, x, indexing="ij")
RR = np.sqrt(X ** 2 + Y ** 2)
RECV = {"axis": (C0 + 100, C0), "diag": (C0 + int(round(100 / np.sqrt(2))), C0 + int(round(100 / np.sqrt(2))))}
EPS = (1e-3, 1e-6, 1e-9, 1e-12)
SOURCES = [("smooth", 3), ("smooth", 6), ("smooth", 12), ("smooth", 24), ("disk", 12)]


def source(kind, Rs):
    if kind == "smooth":
        f = np.where(RR < Rs, (1 - (RR / Rs) ** 2) ** 4, 0.0)
    else:
        f = np.where(RR < Rs, 1.0, 0.0)
    psi = np.zeros((3, NB, NB), complex)
    psi[2] = f
    return psi / np.sqrt(np.sum(np.abs(psi) ** 2))


def run_causal(sname, kind, Rs):
    w = composed_walk("light", S4_SEQ, dim=2) if sname == "S4" else make_walk("light", zigzag_ticks(2), dim=2)
    rw = RealWalk(w, NB, 2)
    T = w.T
    psi = source(kind, Rs)
    nper = int(round(120 / T))
    rec = {k: [] for k in RECV}
    times, outside = [], {}
    axis_line = None
    for p in range(1, nper + 1):
        psi = rw.period(psi)
        t = p * T
        rho = np.sum(np.abs(psi) ** 2, axis=0)
        times.append(t)
        for k, (i, j) in RECV.items():
            rec[k].append(float(rho[i, j]))
        if any(abs(t - tt) < 1e-9 for tt in (30, 60, 90, 120)):
            outside[str(int(round(t)))] = float(np.sum(rho[RR > t + Rs + 1]) / np.sum(rho))
        if abs(t - 120) < 1e-9:
            axis_line = rho[C0:, C0].tolist()
            diag_line = [float(rho[C0 + q, C0 + q]) for q in range(C0)]
            # where is the outside energy? radius of the outermost cell with rho >= 1e-30 and the 'front radius'
            # r_eps(120) = largest r (axis and diagonal) with rho >= eps * max(rho on the same line)
            front = {}
            for name, line, step in (("axis", np.array(axis_line), 1.0), ("diag", np.array(diag_line), np.sqrt(2))):
                mx = line.max()
                for e in EPS:
                    idx = np.nonzero(line >= e * mx)[0]
                    front[f"{name}|{e:g}"] = float(idx.max() * step) if len(idx) else None
    norm = float(np.sum(np.abs(psi) ** 2))
    arr = {}
    for k in RECV:
        s = np.array(rec[k]); tt = np.array(times)
        pk = s.max()
        dist = 100.0 if k == "axis" else float(np.hypot(*(np.array(RECV[k]) - C0)))
        for e in EPS:
            thr = e * pk
            j = np.nonzero(s >= thr)[0]
            if len(j) == 0:
                arr[f"{k}|{e:g}"] = None; continue
            j = j[0]
            if j == 0:
                ta = tt[0]
            else:              # log-linear interpolation between the two samples (time resolution = one period T)
                l0, l1 = np.log(max(s[j - 1], 1e-300)), np.log(s[j])
                ta = tt[j - 1] + (np.log(thr) - l0) / (l1 - l0) * (tt[j] - tt[j - 1])
            arr[f"{k}|{e:g}"] = dict(t=float(ta), v=float((dist - Rs) / ta), dist=dist)
    return dict(T=T, norm_end=norm, outside=outside, arrival=arr, front_t120=front,
                axis_line_t120=[float(v) for v in axis_line[:200:1]], diag_line_t120=diag_line[:180],
                recv_series={k: rec[k] for k in rec}, times=times)


for sname in ("frozen", "S4"):
    for kind, Rs in SOURCES:
        key = f"P6|{sname}|{kind}|{Rs}"
        if not todo(key):
            continue
        R[key] = run_causal(sname, kind, Rs)
        save()

# ------------------------------------------------------------------ P7 variants
if todo("P7|enum"):
    found = {}
    for c, L in ((2, 13), (3, 13), (4, 11)):
        vals = [m for m in range(-c, c + 1) if m != 0]
        best = None
        for Ln in range(3, L + 1):
            h = (Ln + 1) // 2
            for half in itertools.product(vals, repeat=h):
                seq = half + half[::-1][Ln % 2:]
                if max(abs(m) for m in seq) != c:
                    continue
                s1 = sum(seq)
                if s1 <= 0 or sum(m ** 3 for m in seq) != 0:
                    continue
                cn = cone(seq)
                if best is None or cn < best[0] - 1e-12:
                    best = (cn, seq)
        found[str(c)] = dict(best_cone=best[0] if best else None, best_seq=list(best[1]) if best else None,
                             bound=(c ** 2 + 1) / (c ** 2 - 1),
                             ok=(best is None) or best[0] >= (c ** 2 + 1) / (c ** 2 - 1) - 1e-12)
    R["P7|enum"] = found
    save()
if todo("P7|c4_N1"):
    w4 = composed_walk("light", seq_c(4), dim=2)
    R["P7|c4_N1"] = dict(seq_len=len(seq_c(4)), T=w4.T, cone=cone(seq_c(4)),
                         N1=float(n_for_tol(w4, dirs_for(2)[::3], hi=600.0)),
                         cost_blocks_per_time=sum(abs(m) for m in seq_c(4)) / sum(seq_c(4)))
    save()
if todo("P7|c4_fold"):
    w4 = composed_walk("light", seq_c(4), dim=2)
    # first avoided crossing of the light bands (expected near k T = 2 pi, i.e. N ~ T)
    res = []
    for nn in dirs_for(2)[::9]:
        for km, g in ray_crossings(w4, nn, 0.005, 2 * PI / 13, steps=1500):
            res.append(dict(k=km, N=2 * PI / km, gap=g))
    av = [r for r in res if r["gap"] > 1e-6]
    R["P7|c4_fold"] = dict(n_avoided=len(av), largest_N_avoided=max([r["N"] for r in av], default=None),
                           max_gap=max([r["gap"] for r in av], default=None))
    p = [min_sep(dirac1_bloch(k, TH4, seq_c(4))[0]) for k in np.linspace(0, 1.0, 2001)]
    R["P7|c4_fold"]["dirac1_min_sep_k<1"] = float(min(p[1:]))
    R["P7|c4_fold"]["dirac1_window_E_over_m"] = float(PI / sum(seq_c(4)) / TH4)
    save()

save()
print("ALL DONE", time.time() - T0)
