"""motore_f3_run.py - phase 3 runs (pre-registration lab/out/motore_f3_prereg.json, committed before these runs).
Usage: python lab/motore_f3_run.py lot                 part 1 (bell selection, G3a/b/c bis)
       python lab/motore_f3_run.py sweep               part 3 (Intel box sweep + automatic choice)
       python lab/motore_f3_run.py ultra               part 3 (Ultra at the chosen box; checkpointed: re-run until DONE)
       python lab/motore_f3_run.py card                part 2 states at Intel +-27 a0 (checkpointed)
       python lab/motore_f3_run.py truth <scenario>    part 2 grid vs cartolina (k0.05 k0.2 k0.5 k1.0 S1 S2)
       python lab/motore_f3_run.py judge
Writes lab/out/motore_f3_<part>.json; state vectors lab/out/motore_f3_vec_*.npy (NOT committed).
Choice made while writing this script (before the runs, declared): the cartolina is renormalized to the norm of the
grid state (charge conservation), both for the distance (a) and for the dipole velocity (b)."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys, json, time
import numpy as np
import scipy.fft as sfft
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motore_f2_core as f2
import motore_f3_core as f3
import motore_core as mc

PRE = json.load(open(os.path.join(f3.OUT, "motore_f3_prereg.json")))
CF = PRE["configs"]
T0 = time.time()
J = lambda p: os.path.join(f3.OUT, "motore_f3_%s.json" % p)
VEC = lambda name: os.path.join(f3.OUT, "motore_f3_vec_%s.npy" % name)
A0 = f2.A0


def load(p):
    return json.load(open(J(p))) if os.path.exists(J(p)) else {}


def save(p, d):
    json.dump(d, open(J(p), "w"), indent=1, default=float)


# =============================== part 1 ===============================
def part_lot():
    import reticolo_core as rc
    out = dict(calibration=[], choice=None)
    cands = [(L, k) for L in (32, 64, 128) for k in ("sin2", "sin3", "cinf")]
    for L, kind in cands:
        r = f3.electron_roundtrip(CF["calib_k0_a0"], L, L // 2, kind, CF["halvings"])
        out["calibration"].append(r)
        print("calib", r, flush=True)
        if abs(r["dN_roundtrip"]) <= 1e-7 and abs(r["dEkin_roundtrip"]) <= 1e-7 and abs(r["dN_coarse"]) <= 1e-7:
            out["choice"] = dict(L=L, e=L // 2, bell=kind)
            break
    L, e, kind = out["choice"]["L"], out["choice"]["e"], out["choice"]["bell"]
    out["electron"] = [f3.electron_roundtrip(k0, L, e, kind, CF["halvings"]) for k0 in CF["test_k0_a0"]]
    print("G3c_bis", out["electron"], flush=True)
    out["light"] = []
    for Nc in (32, 16, 12, 8):
        Nf = 2 * Nc
        nf = int(np.ceil(27 * Nf / L) * L)
        Ff = rc.packet((nf, nf), 2, Nf, 30.0, [nf / 2, nf / 2], 3 * Nf)
        Ef = mc.energy(Ff)
        t0 = time.time()
        Fc = f3.lot_restrict(Ff, L, e, (0, 1), kind)
        Fb = f3.lot_prolong(Fc, L, e, (0, 1), kind)
        dt = time.time() - t0
        row = dict(N_coarse=Nc, dE_restrict=4 * mc.energy(Fc) / Ef - 1, dE_roundtrip=mc.energy(Fb) / Ef - 1,
                   fidelity=float(abs(mc.overlap(Ff, Fb)) ** 2), sec_per_fine_cell=dt / Ff.size * 3 / 3)
        out["light"].append(row)
        print("G3a_bis", row, flush=True)
    # locality with the declared support
    n = 8 * L
    rng = np.random.default_rng(1)
    base = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    pert = base.copy()
    b0 = 4 * L
    pert[b0:b0 + L, b0:b0 + L] += rng.normal(size=(L, L))
    dC = np.abs(f3.lot_restrict(pert, L, e, (0, 1), kind) - f3.lot_restrict(base, L, e, (0, 1), kind))
    lo, hi = (b0 - L) // 2 - e // 2, (b0 + 2 * L) // 2 + e // 2
    xs = np.arange(n // 2)
    ins = (xs >= lo) & (xs < hi)
    outside = max(float(dC[~ins, :].max()), float(dC[:, ~ins].max()))
    rows_changed = np.where(dC.max(1) > 1e-14)[0]
    out["locality"] = dict(block_fine=[b0, b0 + L], declared_support_coarse=[lo, hi], max_change_outside=outside,
                           max_change_inside=float(dC.max()), measured_support=[int(rows_changed.min()), int(rows_changed.max()) + 1],
                           support_width_fine_cells=2 * (hi - lo))
    print("G3b_bis", out["locality"], flush=True)
    el = out["electron"]
    out["G3c_bis"] = "PASS" if all(abs(r["dN_coarse"]) <= 1e-6 and abs(r["dN_roundtrip"]) <= 1e-6 and abs(r["dEkin_roundtrip"]) <= 1e-6 for r in el) else "FAIL"
    out["G3a_bis"] = "PASS" if all(abs(r["dE_restrict"]) <= 1e-6 and abs(r["dE_roundtrip"]) <= 1e-6 for r in out["light"] if r["N_coarse"] >= 12) else "FAIL"
    out["G3b_bis"] = "PASS" if outside <= 1e-14 else "FAIL"
    print(out["choice"], out["G3a_bis"], out["G3b_bis"], out["G3c_bis"])
    save("lot", out)


# =============================== lean octant atom + single-vector LOBPCG ===============================
class Lean(f2.OctAtom):
    """as OctAtom but without the cached exponentials (memory): W computes them on the fly"""
    def __init__(self, lev, n, sector):
        self.n, self.lev, self.sector = n, lev, sector
        self.h, self.tau = lev["h"], lev["tau"]
        w = lev["w"]
        j = np.arange(-w, w + 1)
        kc = np.pi * np.arange(n) / n
        ks = np.pi * (np.arange(n) + 1) / n
        Cx = np.cos(np.outer(kc, j))
        Cz = np.cos(np.outer(ks if sector == "odd" else kc, j))
        self.G = np.einsum("abc,pa,qb,rc->pqr", lev["g"], Cx, Cx, Cz, optimize=True)
        inv, self.x = f2.cell_avg_inv_r_octant(n, self.h)
        self.V = -f2.zz.ALPHA * inv
        del inv
        self.pre = 1.0 / (self.G / self.tau + 0.5 * f2.zz.ALPHA ** 2)
        self.napply = 0

    def W(self, f):
        D = np.exp(-0.5j * self.V * self.tau)
        g = self.fwd(D * f)
        g *= np.exp(-1j * self.G)
        return D * self.bwd(g)

    def Hv(self, v):
        return self.H(v.reshape((self.n,) * 3)).ravel()

    def Pv(self, v):
        return self.bwd(self.pre * self.fwd(v.reshape((self.n,) * 3))).ravel()


def lobpcg1(at, x, Y=None, maxiter=80, tol=1e-10, ckpt=None, tmax=480):
    """single-vector LOBPCG (Rayleigh-Ritz on x, P r, p), orthogonal to the columns of Y; checkpoint the vector"""
    def orth(v):
        if Y is not None:
            for y in Y:
                v = v - y * np.dot(y, v)
        return v
    x = orth(x); x /= np.linalg.norm(x)
    Ax = at.Hv(x)
    p = Ap = None
    it, rn = 0, np.inf
    for it in range(1, maxiter + 1):
        lam = np.dot(x, Ax)
        r = Ax - lam * x
        rn = float(np.linalg.norm(r))
        if rn <= tol:
            break
        w = orth(at.Pv(r)); w /= np.linalg.norm(w)
        Aw = at.Hv(w)
        S = [x, w] + ([p] if p is not None else [])
        AS = [Ax, Aw] + ([Ap] if p is not None else [])
        Hs = np.array([[np.dot(a, b) for b in AS] for a in S])
        Ms = np.array([[np.dot(a, b) for b in S] for a in S])
        Hs = (Hs + Hs.T) / 2; Ms = (Ms + Ms.T) / 2
        import scipy.linalg as sl
        ev, C = sl.eigh(Hs, Ms)
        c = C[:, 0]
        xn = sum(ci * si for ci, si in zip(c, S))
        Axn = sum(ci * si for ci, si in zip(c, AS))
        pn = sum(ci * si for ci, si in zip(c[1:], S[1:]))
        Apn = sum(ci * si for ci, si in zip(c[1:], AS[1:]))
        nx = np.linalg.norm(xn); np_ = np.linalg.norm(pn)
        x, Ax = xn / nx, Axn / nx
        p, Ap = pn / np_, Apn / np_
        del S, AS, w, Aw, xn, Axn, pn, Apn
        if ckpt and it % 10 == 0:
            np.save(ckpt, x)
        if time.time() - T0 > tmax:
            if ckpt:
                np.save(ckpt, x)
            return x, rn, it, False
    return x, rn, it, True


def upsample_from(v, n_from, n_to, sector):
    """spectral (DCT/DST) zero-padding of an octant vector to a finer grid on the same box (initial guess only)"""
    F = Lean.__new__(Lean); F.sector = sector
    c = F.fwd(v.reshape((n_from,) * 3))
    C = np.zeros((n_to,) * 3)
    C[:n_from, :n_from, :n_from] = c
    T = Lean.__new__(Lean); T.sector = sector
    out = T.bwd(C).ravel()
    return out / np.linalg.norm(out)


def kick_numbers(at_e, at_o, v1s, v2p, q):
    n = at_e.n
    X, Y_, Z = at_e.grid()
    psi1 = v1s.reshape((n,) * 3)
    ev = np.cos(q * Z / A0) * psi1
    c1 = float(np.vdot(v1s, ev.ravel()))
    del ev
    od = np.sin(q * Z / A0) * psi1
    c2p = float(np.vdot(v2p, od.ravel()))
    return c1 ** 2, c2p ** 2


# =============================== part 3: sweep ===============================
def part_sweep():
    lev = f2.load_level3(512)
    res = load("sweep")
    rows = res.setdefault("rows", {})
    q = CF["q_a0"]
    for b in (12, 15, 18, 21, 24, 27):
        if str(b) in rows:
            continue
        if time.time() - T0 > 420:
            save("sweep", res); print("again"); return
        n = 8 * b
        ae = Lean(lev, n, "even")
        v1, r1, i1, _ = lobpcg1(ae, ae.guess("1s"))
        E1 = ae.energy(v1)[0] * f2.MC2
        ao = Lean(lev, n, "odd")
        v2, r2, i2, _ = lobpcg1(ao, ao.guess("2p"))
        E2 = ao.energy(v2)[0] * f2.MC2
        P1s, P2p = kick_numbers(ae, ao, v1, v2, q)
        rows[str(b)] = dict(box_a0=b, octant=n, E1s=E1, E2p=E2, line=E2 - E1, P1s=P1s, P2p=P2p, resid=[r1, r2], iters=[i1, i2])
        if b == 18:
            np.save(VEC("intel18_1s"), v1); np.save(VEC("intel18_2p"), v2)
        print(rows[str(b)], flush=True)
        save("sweep", res)
    bs = sorted(int(k) for k in rows)
    ch = None
    for b in bs[:-1]:
        a, c = rows[str(b)], rows[str(b + 3)]
        d = dict(P2p=abs(a["P2p"] / c["P2p"] - 1), line=abs(a["line"] / c["line"] - 1), P1s=abs(a["P1s"] / c["P1s"] - 1))
        rows[str(b)]["change_to_next"] = d
        if ch is None and max(d.values()) <= CF["tol_ultra"] / 4:
            ch = b
    res["chosen_box_a0"] = ch
    res["box_error_estimate_P2p"] = abs(rows[str(ch)]["P2p"] / rows[str(bs[-1])]["P2p"] - 1) if ch else None
    print("chosen", ch, res["box_error_estimate_P2p"])
    save("sweep", res)


# =============================== part 3: Ultra at the chosen box ===============================
def part_ultra():
    sw = load("sweep")
    b = sw["chosen_box_a0"]
    n = 16 * b
    lev = f2.load_level3(256)
    res = load("ultra")
    res.update(box_a0=b, octant=n)
    q = CF["q_a0"]
    for kind, sec in (("1s", "even"), ("2p", "odd")):
        if kind in res:
            continue
        at = Lean(lev, n, sec)
        ck = VEC("ultra_%s_ckpt" % kind)
        if os.path.exists(ck):
            x0 = np.load(ck)
        elif b == 18:
            x0 = upsample_from(np.load(VEC("intel18_" + kind)), 8 * b, n, sec)
        else:
            x0 = at.guess(kind)
        t0 = time.time()
        v, rn, it, done = lobpcg1(at, x0, ckpt=ck, tmax=500)
        if not done and rn > 1e-10:
            save("ultra", res); print("again", kind, rn, it); return
        np.save(VEC("ultra_" + kind), v)
        E, rr, lam = at.energy(v)
        res[kind] = dict(E_eV=E * f2.MC2, rule_resid=rr, lobpcg_resid=rn, iters=it, seconds=time.time() - t0)
        print(kind, res[kind], flush=True)
        save("ultra", res)
        del at
        if time.time() - T0 > 300:
            print("again"); return
    if "P2p" not in res:
        ae = Lean.__new__(Lean); ae.n = n; ae.h = lev["h"]; ae.x = (np.arange(n) + 0.5) * lev["h"]
        P1s, P2p = kick_numbers(ae, None, np.load(VEC("ultra_1s")), np.load(VEC("ultra_2p")), q)
        res.update(P1s=P1s, P2p=P2p, line_eV=res["2p"]["E_eV"] - res["1s"]["E_eV"])
        save("ultra", res)
    print(res); print("DONE")


# =============================== part 2: cartolina states at Intel +-27 ===============================
def guess_n3(at, kind):
    X, Y, Z = at.grid()
    x, y, z = X / A0, Y / A0, Z / A0
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    g = {"3s": (27 - 18 * r + 2 * r ** 2) * np.exp(-r / 3), "3d0": (3 * z ** 2 - r ** 2) * np.exp(-r / 3),
         "3dxy": (x ** 2 - y ** 2) * np.exp(-r / 3), "3p": z * (6 - r) * np.exp(-r / 3)}[kind]
    return (g / np.linalg.norm(g)).ravel()


def part_card():
    lev = f2.load_level3(512)
    n = 216
    res = load("card")
    res.update(level="Intel a0/8", octant=n, box_a0=n / 8)
    # odd sector first (smaller set)
    if "3p" not in res:
        ao = Lean(lev, n, "odd")
        for kind, Yk in (("2p", None), ("3p", ["2p"])):
            if kind in res:
                continue
            Y = [np.load(VEC("card_" + k)) for k in Yk] if Yk else None
            x0 = ao.guess("2p") if kind == "2p" else guess_n3(ao, "3p")
            v, rn, it, done = lobpcg1(ao, x0, Y=Y, maxiter=150, tmax=520)
            if not done and rn > 1e-10:
                print("again", kind, rn); save("card", res); return
            np.save(VEC("card_" + kind), v)
            E, rr, lam = ao.energy(v)
            res[kind] = dict(E_eV=E * f2.MC2, rule_resid=rr, lobpcg_resid=rn, iters=it, lam=[lam.real, lam.imag])
            print(kind, res[kind], flush=True); save("card", res)
        del ao
        print("again"); return
    ae = Lean(lev, n, "even")
    for kind, Yk in (("1s", None), ("2s", ["1s"])):
        if kind in res:
            continue
        Y = [np.load(VEC("card_" + k)) for k in Yk] if Yk else None
        v, rn, it, done = lobpcg1(ae, ae.guess(kind), Y=Y, maxiter=150, tmax=520)
        if not done and rn > 1e-10:
            print("again", kind, rn); save("card", res); return
        np.save(VEC("card_" + kind), v)
        E, rr, lam = ae.energy(v)
        res[kind] = dict(E_eV=E * f2.MC2, rule_resid=rr, lobpcg_resid=rn, iters=it, lam=[lam.real, lam.imag])
        print(kind, res[kind], flush=True); save("card", res)
        if time.time() - T0 > 300:
            print("again"); return
    if "n3" not in res:
        N3 = n ** 3
        Y = np.stack([np.load(VEC("card_1s")), np.load(VEC("card_2s"))], 1)
        X0 = np.stack([guess_n3(ae, k) for k in ("3s", "3d0", "3dxy")], 1)
        X0 = X0 - Y @ (Y.T @ X0)
        A = spl.LinearOperator((N3, N3), matvec=ae.Hv, matmat=lambda M: np.stack([ae.Hv(c) for c in M.T], 1), dtype=float)
        P = spl.LinearOperator((N3, N3), matvec=ae.Pv, matmat=lambda M: np.stack([ae.Pv(c) for c in M.T], 1), dtype=float)
        w, V, hist = spl.lobpcg(A, X0, M=P, Y=Y, tol=1e-10, maxiter=60, largest=False, retResidualNormsHistory=True)
        rows = []
        g3s = guess_n3(ae, "3s")
        for i in range(3):
            v = V[:, i] / np.linalg.norm(V[:, i])
            E, rr, lam = ae.energy(v)
            s_like = float(np.dot(v, g3s) ** 2)
            np.save(VEC("card_n3_%d" % i), v)
            rows.append(dict(E_eV=E * f2.MC2, rule_resid=rr, s_overlap=s_like, lam=[lam.real, lam.imag]))
        res["n3"] = dict(states=rows, lobpcg_resid=float(np.max(hist[-1])), iters=len(hist))
        print("n3", res["n3"], flush=True); save("card", res)
    print("DONE")


def basis(n):
    ev = [("1s", VEC("card_1s")), ("2s", VEC("card_2s"))] + [("n3_%d" % i, VEC("card_n3_%d" % i)) for i in range(3)]
    od = [("2p", VEC("card_2p")), ("3p", VEC("card_3p"))]
    return ev, od


def part_truth(scn):
    lev = f2.load_level3(512)
    C = load("card")
    n = C["octant"]
    ae, ao = Lean(lev, n, "even"), Lean(lev, n, "odd")
    X, Y_, Z = ae.grid()
    del X, Y_
    evb, odb = basis(n)
    lam = {k: complex(*C[k]["lam"]) for k in ("1s", "2s", "2p", "3p")}
    for i in range(3):
        lam["n3_%d" % i] = complex(*C["n3"]["states"][i]["lam"])
    v1 = np.load(VEC("card_1s")).reshape((n,) * 3)
    if scn.startswith("k"):
        q = float(scn[1:])
        fe = (np.cos(q * Z / A0) * v1).astype(complex)
        fo = (1j * np.sin(q * Z / A0) * v1).astype(complex)
    else:
        v2p = np.load(VEC("card_2p")).reshape((n,) * 3)
        if scn == "S1":
            fe = v1.astype(complex); fo = 0.3j * v2p
        else:
            fe = v1 + 0.2 * np.load(VEC("card_2s")).reshape((n,) * 3)
            # 3s = the n3 state with the largest s overlap
            i3s = int(np.argmax([s["s_overlap"] for s in C["n3"]["states"]]))
            fe = (fe + 0.1 * np.load(VEC("card_n3_%d" % i3s)).reshape((n,) * 3)).astype(complex)
            fo = 0.2j * np.load(VEC("card_3p")).reshape((n,) * 3)
    del v1
    nrm = np.sqrt(np.sum(np.abs(fe) ** 2) + np.sum(np.abs(fo) ** 2))
    fe /= nrm; fo /= nrm
    # cartolina coefficients and dipole matrix elements
    ce, co, zeo = {}, {}, {}
    for k, fn in evb:
        ce[k] = complex(np.vdot(np.load(fn), fe.ravel()))
    for k, fn in odb:
        co[k] = complex(np.vdot(np.load(fn), fo.ravel()))
    for ke, fne in evb:
        ve = np.load(fne).reshape((n,) * 3)
        for ko, fno in odb:
            zeo[(ke, ko)] = float(np.sum(ve * Z * np.load(fno).reshape((n,) * 3)))
        del ve
    w_in = sum(abs(c) ** 2 for c in ce.values()) + sum(abs(c) ** 2 for c in co.values())
    w_out = 1 - w_in
    Pgs = abs(ce["1s"]) ** 2
    rho = w_out / max(1 - Pgs, 1e-300)
    tol = CF["tol_intel"]
    admitted = (w_out <= tol ** 2) and (rho <= tol)

    def z_cart(s):
        acc = 0j
        for ke, cev in ce.items():
            for ko, cov in co.items():
                acc += np.conj(cev * lam[ke] ** s) * cov * lam[ko] ** s * zeo[(ke, ko)]
        return 2 * acc.real / w_in

    def z_grid(fe, fo):
        return float(2 * np.real(np.sum(np.conj(fe) * Z * fo)) / (np.sum(np.abs(fe) ** 2) + np.sum(np.abs(fo) ** 2)))
    K = CF["v_fit_steps"]
    M = CF["M_truth"]
    tau = lev["tau"]
    zg = [z_grid(fe, fo)]
    t0 = time.time()
    for s in range(1, M + 1):
        fe = ae.W(fe); fo = ao.W(fo)
        if s <= K:
            zg.append(z_grid(fe, fo))
        if time.time() - T0 > 560:
            print("TIMEOUT at step", s); return
    sec_step = (time.time() - t0) / M
    t = np.arange(K + 1) * tau
    zc = [z_cart(s) for s in range(K + 1)]
    vg = float(np.polyfit(t, zg, 2)[1])
    vc = float(np.polyfit(t, zc, 2)[1])
    # distance after M steps (cartolina renormalized)
    num = np.sum(np.abs(fe) ** 2) + np.sum(np.abs(fo) ** 2)
    sc = 1 / np.sqrt(w_in)
    pe = sum(ce[k] * lam[k] ** M * np.load(fn) for k, fn in evb).reshape((n,) * 3) * sc
    de = np.sum(np.abs(fe - pe) ** 2); del pe
    po = sum(co[k] * lam[k] ** M * np.load(fn) for k, fn in odb).reshape((n,) * 3) * sc
    do = np.sum(np.abs(fo - po) ** 2); del po
    dist = float(np.sqrt((de + do) / num))
    # subspace check (as G2): grid projected vs cartolina phases
    sub = max([abs(complex(np.vdot(np.load(fn), fe.ravel())) - ce[k] * lam[k] ** M) for k, fn in evb] +
              [abs(complex(np.vdot(np.load(fn), fo.ravel())) - co[k] * lam[k] ** M) for k, fn in odb])
    truth = dist <= tol and abs(vc / vg - 1) <= tol
    row = dict(scenario=scn, P={k: abs(c) ** 2 for k, c in {**ce, **co}.items()}, w_out=w_out, rho=rho, admitted=admitted,
               dist_M=dist, v_grid=vg, v_cart=vc, v_ratio=vc / vg, truth=truth, subspace_diff=sub,
               v_expected_hbar_q_over_m=(float(scn[1:]) / A0 if scn.startswith("k") else None),
               sec_per_step_grid=sec_step, z_grid=zg, z_cart=zc)
    print(json.dumps(row, default=float, indent=1))
    res = load("truth")
    res[scn] = row
    save("truth", res)


# =============================== judge ===============================
def judge():
    Lo, Sw, U, C, T = load("lot"), load("sweep"), load("ultra"), load("card"), load("truth")
    P = PRE["part2"]["predictions"]["exact_P"]
    out = dict(part1=dict(choice=Lo["choice"], G3a_bis=Lo["G3a_bis"], G3b_bis=Lo["G3b_bis"], G3c_bis=Lo["G3c_bis"]))
    # part 3
    ex2p, ex1s, exl = P["2p"], P["1s"], 0.5 * f2.zz.ALPHA ** 2 * f2.MC2 * 0.75
    b = Sw["chosen_box_a0"]
    B1 = dict(P2p_rel=U["P2p"] / ex2p - 1, P1s_rel=U["P1s"] / ex1s - 1, line_rel=U["line_eV"] / exl - 1)
    B1["PASS"] = all(abs(v) <= CF["tol_ultra"] for v in B1.values())
    A2 = json.load(open(os.path.join(f3.OUT, "motore_f2_atom.json")))
    dU = A2["ultra"]["kick"]["P2p"] / U["P2p"] - 1
    dI = Sw["rows"]["12"]["P2p"] / Sw["rows"][str(b)]["P2p"] - 1
    out["part3"] = dict(chosen_box=b, box_error_estimate=Sw["box_error_estimate_P2p"], B1=B1,
                        B2=dict(delta_ultra=dU, delta_intel=dI, PASS=abs(dU - dI) <= 0.005))
    # part 2
    k5 = T.get("k0.5", {})
    lat = k5.get("P", {})
    st = C["n3"]["states"]
    i3s = int(np.argmax([s["s_overlap"] for s in st]))
    p3s = lat.get("n3_%d" % i3s, np.nan)
    p3d = sum(v for k, v in lat.items() if k.startswith("n3_") and k != "n3_%d" % i3s)
    c1 = {k: lat.get(k, np.nan) / P[k] - 1 for k in ("1s", "2s", "2p", "3p")}
    c1.update({"3s": p3s / P["3s"] - 1, "3d": p3d / P["3d"] - 1})
    C1 = all(abs(c1[k]) <= CF["tol_intel"] for k in ("1s", "2s", "2p", "3p"))
    resid = max([C[k]["rule_resid"] for k in ("1s", "2s", "2p", "3p")] + [s["rule_resid"] for s in st])
    C2 = all(r["subspace_diff"] <= CF["M_truth"] * resid * 2 + 1e-12 for r in T.values())
    C3 = {}
    for scn, r in T.items():
        if scn.startswith("k"):
            e = f3.exact_kick_probs(float(scn[1:]))
            C3[scn] = dict(lattice=r["w_out"], exact=1 - sum(e.values()), rho=r["rho"])
    C4 = all(r["truth"] for r in T.values() if r["admitted"])
    C5 = all(T[s]["admitted"] for s in ("S1", "S2") if s in T)
    C6 = {s: r["v_grid"] / r["v_expected_hbar_q_over_m"] - 1 for s, r in T.items() if s.startswith("k")}
    out["part2"] = dict(C1=dict(rel=c1, PASS=C1), C2=dict(max_resid=resid, PASS=C2), C3=C3, C4=C4, C5=C5,
                        C6=dict(rel=C6, PASS=all(abs(v) <= CF["tol_intel"] for v in C6.values())),
                        table={s: dict(w_out=r["w_out"], rho=r["rho"], admitted=r["admitted"], dist=r["dist_M"],
                                       v_ratio=r["v_ratio"], truth=r["truth"]) for s, r in T.items()})
    print(json.dumps(out, indent=1, default=float))
    save("judge", out)


if __name__ == "__main__":
    p = sys.argv[1]
    if p == "truth":
        part_truth(sys.argv[2])
    else:
        {"lot": part_lot, "sweep": part_sweep, "ultra": part_ultra, "card": part_card, "judge": judge}[p]()
