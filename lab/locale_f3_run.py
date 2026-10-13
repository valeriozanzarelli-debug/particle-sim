"""locale_f3_run.py - PHASE 3 of REGOLA LOCALE (pre-registered in lab/out/locale_f3_prereg.json):
 (1) Schwinger / vacuum pairs per valley, frozen F0 vs twisted Wilson F2; (1b) h from a principle (max doubler gap);
 (2) S4 radial hydrogen with the true finite proton + softened-Coulomb diagnostic; (3) light precursors in 3-D.
Checkpointed into lab/out/locale_f3_run.json; exit code 3 = run again."""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.optimize import minimize_scalar
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import make_walk, composed_walk, zigzag_ticks, RealWalk, S4_SEQ
from locale_cands import fermion
from locale_s4_core import SEQ_S4, SEQ_FROZEN, RadialS4
from zz_core import dirac_formula, ALPHA

PI = np.pi
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f3_run.json")
T0 = time.time()
BUDGET = 480.0
R = json.load(open(OUT)) if os.path.exists(OUT) else {}


def save():
    with open(OUT + ".tmp", "w") as f:
        json.dump(R, f, indent=1)
    os.replace(OUT + ".tmp", OUT)


def todo(key, need=0.0):
    if key in R:
        return False
    if time.time() - T0 > BUDGET - need:
        save(); print("budget reached, rerun", flush=True); sys.exit(3)
    print(f"[{time.time() - T0:6.1f}s] {key}", flush=True)
    return True


TH = 0.1
VALLEYS = {"G": ((0, 0, 0), 0.0), "N_xy": ((.5, .5, 0), 0.0), "N_xz": ((.5, 0, .5), 0.0), "N_yz": ((0, .5, .5), 0.0),
           "X_x": ((.5, 0, 0), PI), "X_y": ((0, .5, 0), PI), "X_z": ((0, 0, .5), PI), "R_xyz": ((.5, .5, .5), PI)}


def rel_phase(eps, c):
    return np.angle(np.exp(-1j * (eps - c)))      # phase measured from the valley energy c, in (-pi, pi]


def subspace(W, c, sign):
    """orthonormal basis (N, 4, 2) of the 2 eigenvectors of W with sign(rel phase) = sign"""
    w, v = np.linalg.eig(W)
    d = rel_phase(-np.angle(w), c)
    idx = np.argsort(sign * d, axis=1)[:, ::-1][:, :2] if sign > 0 else np.argsort(d, axis=1)[:, :2]
    B = np.take_along_axis(v, idx[:, None, :], axis=2)
    Q, _ = np.linalg.qr(B)
    return Q


# ------------------------------------------------------------------ (1b) h principle (first: needed for h*)
def zone_gap(w, n=24, kc=0.4):
    g = (np.arange(n) + 0.5) / n * PI - PI / 2
    K = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    K = K[np.linalg.norm(K, axis=1) >= kc]
    eps, _ = w.phases(K)
    g0 = np.min(np.abs(eps)) / w.T
    gpi = np.min(PI - np.abs(eps)) / w.T
    vg = {}
    for name, (kv, c) in VALLEYS.items():
        if name == "G":
            continue
        e, _ = w.phases(PI * np.array(kv, float)[None])
        vg[name] = float(np.min(np.abs(rel_phase(e[0], c))) / w.T)
    return float(min(g0, gpi, min(vg.values()))), float(g0), float(gpi), vg


def E_low(w, K):
    eps, _ = w.phases(np.atleast_2d(K))
    return np.min(np.abs(eps), axis=1) / w.T


def d_exams(w, th):
    dirs = [np.array([1, 0, 0.]), np.array([1, 1, 0]) / 2 ** .5, np.array([1, 1, 1]) / 3 ** .5]
    errs, mk = [], []
    for d in dirs:
        P = np.linspace(0.1, 1, 10)[:, None] * th * d[None]
        errs.append(np.max(np.abs(E_low(w, P) / np.hypot(np.linalg.norm(P, axis=1), th) - 1)))
        ps = np.linspace(0.02, 0.1, 9) * th
        E = E_low(w, ps[:, None] * d[None]); E0 = E_low(w, np.zeros((1, 3)))[0]
        c = np.linalg.lstsq(np.stack([ps ** 2, ps ** 4], 1), E ** 2 - E0 ** 2, rcond=None)[0]
        mk.append(1.0 / c[0] - 1)
    return float(max(errs)), [float(x) for x in mk]


HGRID = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, PI / 2]
for th in (0.1, 0.05, 0.02):
    key = f"H|scan|{th}"
    if not todo(key):
        continue
    rows = []
    for h in HGRID:
        G, g0, gpi, vg = zone_gap(fermion("F2", 3, th, h=h))
        rows.append(dict(h=h, G=G, g0=g0, gpi=gpi, valleys=vg))
        print(f"   th={th} h={h:.3f} G={G:.4f} (0: {g0:.4f}, pi: {gpi:.4f})", flush=True)
    Gs = np.array([r["G"] for r in rows]); i = int(np.argmax(Gs))
    interior = 0 < i < len(HGRID) - 1
    hstar = HGRID[i]
    if interior:
        r = minimize_scalar(lambda h: -zone_gap(fermion("F2", 3, th, h=h))[0], bracket=(HGRID[i - 1], HGRID[i], HGRID[i + 1]),
                            tol=1e-3)
        hstar = float(r.x)
    Gstar = zone_gap(fermion("F2", 3, th, h=hstar))
    disp, mk = d_exams(fermion("F2", 3, th, h=hstar), th)
    disp0, mk0 = d_exams(fermion("F0", 3, th), th)
    R[key] = dict(rows=rows, interior_max=bool(interior), h_star=hstar, G_star=Gstar[0], G_star_parts=Gstar[1:3],
                  disp_hstar=disp, mkin_hstar=mk, disp_frozen=disp0, mkin_frozen=mk0, G_over_m=Gstar[0] / th)
    save()

HSTAR = R["H|scan|0.1"]["h_star"]


# ------------------------------------------------------------------ (1) Schwinger per valley
def schwinger_valley(w, kv, c, EEc, half=0.25, step=0.0125):
    E = EEc * TH ** 2
    T = w.T
    g = np.arange(-half, half + 1e-12, step)
    KX, KY = np.meshgrid(g, g, indexing="ij")
    kperp = np.stack([KX.ravel(), KY.ravel()], 1)
    kv = PI * np.array(kv, float)
    N = len(kperp)
    K0 = np.zeros((N, 3)); K0[:, 0] = kv[0] + kperp[:, 0]; K0[:, 1] = kv[1] + kperp[:, 1]
    nper = int(round(0.8 / (E * T)))

    def Kz(kz):
        K = K0.copy(); K[:, 2] = kz; return K
    psi = subspace(w.bloch(Kz(kv[2] - 0.4)), c, -1)                      # (N, 4, 2)
    for p in range(nper):
        psi = w.bloch(Kz(kv[2] - 0.4 + E * (p + 0.5) * T)) @ psi
    Qp = subspace(w.bloch(Kz(kv[2] - 0.4 + E * nper * T)), c, +1)
    P = np.sum(np.abs(np.einsum("nai,naj->nij", Qp.conj(), psi)) ** 2, axis=(1, 2))   # pairs per mode (2 spins)
    Pc = 2 * np.exp(-PI * (TH ** 2 + np.sum(kperp ** 2, 1)) / E)
    e, _ = w.phases(kv[None])
    gap = float(np.min(np.abs(rel_phase(e[0], c))) / T)
    return dict(R=float(P.sum() / Pc.sum()), P_center=float(P[N // 2]), half_gap=gap, periods=nper)


WALKS = {"F0": lambda: fermion("F0", 3, TH), "F2_h0.1": lambda: fermion("F2", 3, TH, h=0.1),
         "F2_hstar": lambda: fermion("F2", 3, TH, h=HSTAR)}
for wn, mk in WALKS.items():
    w = mk()
    for EEc in (1, 2, 4):
        for vn, (kv, c) in VALLEYS.items():
            key = f"S|{wn}|{EEc}|{vn}"
            if not todo(key):
                continue
            R[key] = schwinger_valley(w, kv, c, EEc)
            print(f"   {key}: R = {R[key]['R']:.4e} (half gap {R[key]['half_gap']:.4f})", flush=True)
            save()

# ------------------------------------------------------------------ (2) finite nucleus and softened Coulomb
LAMBDA_C_FM = 386.159
RP = 0.8414 / LAMBDA_C_FM          # proton radius in 1/m


def radial(kap, th, seq, Vfun=None, rmax=30 / ALPHA):
    rs = RadialS4(kap, th, rmax, seq)
    if Vfun is not None:
        V, hx, hz, nrm = rs._loc
        rs._loc = (Vfun(rs.rw.r), hx, hz, nrm)
    return rs


def v_proton(r):
    return np.where(r < RP, -ALPHA * (3 - (r / RP) ** 2) / (2 * RP), -ALPHA / r)


LEV = (("1s", 1, -1), ("2s", 2, -1), ("2p1/2", 2, 1))
for sname, seq in (("S4", SEQ_S4), ("frozen", SEQ_FROZEN)):
    for th in (0.1, 0.05, 0.025):
        for lab, n, kap in LEV:
            key = f"N1|{sname}|{th}|{lab}"
            if not todo(key):
                continue
            E = dirac_formula(n, kap)
            rp = radial(kap, th, seq, v_proton)
            ra = radial(kap, th, seq)
            n_inside = int(np.sum(rp.rw.r < RP))
            ep, _, resp = rp.eigen(E, iters=12)
            ea, _, _ = ra.eigen(E, iters=12)
            R[key] = dict(rel_err_finite=float(((1 - ep) - (1 - E)) / (1 - E)), rel_err_point=float(((1 - ea) - (1 - E)) / (1 - E)),
                          finite_minus_point_rel=float((ea - ep) / (1 - E)), cells_inside_proton=n_inside, resid=resp)
            save()

for b_a0 in (0.02, 0.1):
    b = b_a0 / ALPHA
    Vs = lambda r, b=b: -ALPHA / np.sqrt(r ** 2 + b ** 2)
    for lab, n, kap in LEV:
        key = f"N2|ref|{b_a0}|{lab}"
        if todo(key):
            guess = dirac_formula(n, kap)
            es = []
            for th in (0.05, 0.025, 0.0125):
                e, _, _ = radial(kap, th, SEQ_FROZEN, Vs).eigen(guess, iters=12)
                es.append(e); guess = e
            # E(th) = E0 + a th^2 + b th^4 at th = h, h/2, h/4
            A = np.array([[1, t ** 2, t ** 4] for t in (0.05, 0.025, 0.0125)])
            E0 = float(np.linalg.solve(A, np.array(es))[0])
            R[key] = dict(E_ref=E0, frozen=[float(x) for x in es])
            save()
        for th in (0.1, 0.07, 0.05, 0.04, 0.03, 0.025):
            key = f"N2|S4|{b_a0}|{th}|{lab}"
            if not todo(key):
                continue
            E0 = R[f"N2|ref|{b_a0}|{lab}"]["E_ref"]
            e, _, _ = radial(kap, th, SEQ_S4, Vs).eigen(E0, iters=12)
            R[key] = dict(rel_err=float(((1 - e) - (1 - E0)) / (1 - E0)))
            save()

# ------------------------------------------------------------------ (3) precursors in 3-D
NB = 192
C0 = NB // 2
EPS = (1e-3, 1e-6, 1e-9, 1e-12)


def run3d(sname, Rs):
    w = composed_walk("light", S4_SEQ, dim=3) if sname == "S4" else make_walk("light", zigzag_ticks(3), dim=3)
    rw = RealWalk(w, NB, 3)
    x = np.arange(NB) - C0
    r2 = x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2
    RR = np.sqrt(r2.astype(float)); del r2
    psi = np.zeros((3, NB, NB, NB), complex)
    psi[2] = np.where(RR < Rs, (1 - (RR / Rs) ** 2) ** 4, 0.0)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2))
    q = int(round(50 / np.sqrt(3)))
    recv = {"axis": ((C0 + 50, C0, C0), 50.0), "body": ((C0 + q, C0 + q, C0 + q), float(q * np.sqrt(3)))}
    rec = {k: [] for k in recv}; times = []
    nper = int(round(60 / w.T))
    out = None
    for p in range(1, nper + 1):
        psi = rw.period(psi)
        t = p * w.T; times.append(t)
        rho = np.sum(np.abs(psi) ** 2, axis=0)
        for k, (ijk, _) in recv.items():
            rec[k].append(float(rho[ijk]))
        if p == nper:
            out = float(np.sum(rho[RR > t + Rs + 1]) / np.sum(rho))
            axis_line = rho[C0:, C0, C0].tolist()
            body_line = [float(rho[C0 + s, C0 + s, C0 + s]) for s in range(C0)]
    arr = {}
    for k, (_, dist) in recv.items():
        s = np.array(rec[k]); tt = np.array(times); pk = s.max()
        for e in EPS:
            j = np.nonzero(s >= e * pk)[0]
            if len(j) == 0:
                arr[f"{k}|{e:g}"] = None; continue
            j = j[0]
            if j == 0:
                ta = tt[0]
            else:
                l0, l1 = np.log(max(s[j - 1], 1e-300)), np.log(s[j])
                ta = tt[j - 1] + (np.log(e * pk) - l0) / (l1 - l0) * (tt[j] - tt[j - 1])
            arr[f"{k}|{e:g}"] = dict(t=float(ta), v=float((dist - Rs) / ta))
    return dict(T=w.T, outside_t60=out, arrival=arr, axis_line=axis_line, body_line=body_line)


for Rs in (12, 6):
    for sname in ("frozen", "S4"):
        key = f"C3|{sname}|{Rs}"
        if not todo(key, need=250):
            continue
        R[key] = run3d(sname, Rs)
        print(f"   {key}: outside {R[key]['outside_t60']:.2e}", flush=True)
        save()

save()
print("ALL DONE", time.time() - T0)
