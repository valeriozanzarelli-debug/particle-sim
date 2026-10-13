"""nucleo1s_3d_run.py - production runs of NUCLEO_1S phase 2 (3-D), as pre-registered in lab/out/nucleo1s_3d_prereg.json.
Usage: timeout 590 python lab/nucleo1s_3d_run.py <part> [args]     (each call <= 10 min; re-run until 'done')
parts: lab | learn | scan <a> <cfg> | excited <a> <cfg> | glue | judge
Results -> lab/out/nucleo1s_3d_results.json (checkpoints of long listens -> lab/out/nucleo1s_3d_ck_*.npz, NOT committed)."""
import os, sys, json, time, fcntl
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucleo1s_3d import Colour3D, radial_state, spinor_field, pencil, _dot  # noqa: E402
from zz_core import ALPHA, dirac_formula  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RES = os.path.join(OUT, "nucleo1s_3d_results.json")
PRE = json.load(open(os.path.join(OUT, "nucleo1s_3d_prereg.json")))
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def load():
    return json.load(open(RES)) if os.path.exists(RES) else {}


def save(key, val):
    lk = open(RES + ".lock", "w")
    fcntl.flock(lk, fcntl.LOCK_EX)
    d = load()
    d[key] = val
    tmp = RES + ".tmp"
    json.dump(d, open(tmp, "w"), indent=1, default=float)
    os.replace(tmp, RES)
    fcntl.flock(lk, fcntl.LOCK_UN)


# ------------------------------------------------------------------ the lattice with an (optional) agglomerated nucleus
def make_colour(n, a, Za, s=(1.0, 1.0), q=(0, 0, 0)):
    """s = (s0, s1): factors on the Coulomb phase of the 2 nearest sites (r0 = sqrt3/2 a) and of the 6 next ones
    (r1 = sqrt11/2 a) of the colour. (1, 1) = the bare rule."""
    r0, r1 = np.sqrt(3) / 2 * a, np.sqrt(11) / 2 * a

    def Vf(rr):
        V = -Za / rr
        V = np.where(np.abs(rr - r0) < 1e-6 * a, s[0] * V, V)
        return np.where(np.abs(rr - r1) < 1e-6 * a, s[1] * V, V)
    return Colour3D(n, a, Za, Vfun=Vf, q=q)


def listen_ck(col, starts, nt, tag, refs=None, chunk_s=480):
    """listen from the start state `starts` (E, O) for nt ticks, recording <ref_i|psi(t)> for every ref (default: start).
    Checkpointed: returns None if the time budget of this call ran out (call again)."""
    ck = os.path.join(OUT, f"nucleo1s_3d_ck_{tag}.npz")
    refs = refs if refs is not None else [starts]
    if os.path.exists(ck):
        z = np.load(ck)
        E, O, C, t0 = z["E"], z["O"], z["C"], int(z["t"])
        col.ticks = t0
    else:
        E, O = starts[0].copy(), starts[1].copy()
        C = np.zeros((len(refs), nt + 1), complex)
        t0 = 0
    Eb, Ob = np.empty_like(E), np.empty_like(O)
    tstart = time.time()
    for t in range(t0, nt + 1):
        for i, (re, ro) in enumerate(refs):
            C[i, t] = _dot(re, E) + _dot(ro, O)
        if t == nt:
            break
        Eb, Ob = col.tick(E, O, Eb, Ob)
        E, Eb, O, Ob = Eb, E, Ob, O
        if time.time() - tstart > chunk_s:
            np.savez(ck, E=E, O=O, C=C, t=t + 1)
            log(f"  checkpoint {tag} at tick {t + 1}/{nt}")
            return None
    if os.path.exists(ck):
        os.remove(ck)
    return C


def lines(C, a, wlo, whi, wmin=0.02):
    out = [dict(omega=o, weight=w) for o, w, _, _ in pencil(C, a, wlo, whi) if w >= wmin]
    out2 = [dict(omega=o, weight=w) for o, w, _, _ in pencil(C[:int(0.8 * len(C))], a, wlo, whi) if w >= wmin]
    for x in out:     # stability between two signal lengths
        x["stab"] = min([abs(x["omega"] - y["omega"]) for y in out2] or [np.inf])
    return out


def quartet(L):
    """weighted centroid and weighted std of the lines of the 1s multiplet"""
    w = np.array([x["weight"] for x in L]); o = np.array([x["omega"] for x in L])
    c = float((w * o).sum() / w.sum())
    return dict(centroid=c, std=float(np.sqrt((w * (o - c)**2).sum() / w.sum())), wsum=float(w.sum()),
                spread=float(o.max() - o.min()), n=len(L))


Za = PRE["config"]["Z"] * ALPHA
D1S = dirac_formula(1, -1, Za)


def start_1s(col, a):
    om_r, r, g, f = radial_state(-1, 1, a, Za)
    return col.field(spinor_field(-1, r, g, f, a))


def learned_s():
    d = load().get("learn")
    if not d or "s" not in d:
        raise SystemExit("learn first")
    return tuple(d["s"])


# ------------------------------------------------------------------ parts
def part_lab():
    cfg = PRE["config"]["lab"]
    n, a, nt = cfg["n"], cfg["a"], cfg["nt"]
    out = load().get("lab", {})
    for q in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)):
        if str(q) in out:
            continue
        col = make_colour(n, a, Za, q=q)
        st = start_1s(col, a)
        C = listen_ck(col, st, nt, f"lab_{q[0]}{q[1]}{q[2]}")
        if C is None:
            save("lab", out); return
        L = lines(C[0], a, 0.75, 1.0)
        out[str(q)] = dict(lines=L, C_full_re=C[0].real[::2].tolist(), C_full_im=C[0].imag[::2].tolist())
        log("colour", q, [(round(x["omega"], 6), round(x["weight"], 4)) for x in L])
        save("lab", out)
    # the Laboratorio's own p_exc (read only)
    import csv
    p = os.path.join(PRE["config"]["lab_data"], "fotogrammi.csv")
    rows = list(csv.DictReader(open(p)))
    t = np.array([float(r["t"]) for r in rows]); pe = np.array([float(r["p_exc"]) for r in rows])
    x = (pe - pe.mean()) * np.hanning(len(pe))
    om = np.linspace(0.002, 0.12, 4000)
    S = np.abs(np.exp(1j * np.outer(om, t)) @ x)
    pk = [i for i in range(1, len(om) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1]]
    pk = sorted(pk, key=lambda i: -S[i])[:2]
    out["lab_pexc_peaks"] = [dict(omega=float(om[i]), period_ticks=float(2 * np.pi / om[i] / a), amp=float(S[i])) for i in pk]
    L0 = out["(0, 0, 0)"]["lines"]
    strong = [x for x in L0 if x["weight"] >= 0.1 and 0.80 <= x["omega"] <= 0.92]
    diffs = sorted({abs(x["omega"] - y["omega"]) for x in strong for y in strong if abs(x["omega"] - y["omega"]) > 1e-3})
    match = []
    for pkd in out["lab_pexc_peaks"]:
        best = min(diffs, key=lambda d: abs(d / pkd["omega"] - 1)) if diffs else None
        match.append(dict(peak=pkd["omega"], line_diff=best, rel=None if best is None else best / pkd["omega"] - 1))
    out["B1"] = dict(n_strong=len(strong), match=match,
                     PASS=bool(len(strong) >= 3 and all(m["rel"] is not None and abs(m["rel"]) <= PRE["exams"]["B1"]["tol"] for m in match)))
    oms = [[x["omega"] for x in out[str(q)]["lines"] if x["weight"] >= 0.05] for q in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))]
    ref = sorted(set(round(o, 7) for o in oms[0]))
    dmax = 0.0
    for o in oms[1:]:
        for x in o:
            dmax = max(dmax, min(abs(x - y) for y in ref))
    out["B2"] = dict(max_line_diff=dmax, PASS=bool(dmax <= PRE["exams"]["B2"]["tol"]))
    log("B1", out["B1"], "B2", out["B2"])
    save("lab", out)


def part_learn():
    """the agglomerator learns the 2 numbers of the nucleus at the Laboratorio lattice: criterion = the 1s multiplet of
    the Laboratorio start state becomes ONE level (weighted std of its lines -> min). Nelder-Mead, checkpointed."""
    cfg = PRE["config"]["learn"]
    n, a, nt = cfg["n"], cfg["a"], cfg["nt"]
    d = load().get("learn", {"evals": []})
    cache = {tuple(e["s"]): e["std"] for e in d["evals"]}

    def f(s):
        s = (round(float(s[0]), 6), round(float(s[1]), 6))
        if s in cache:
            return cache[s]
        if time.time() - T0 > 520:
            raise TimeoutError
        if not (0 <= s[0] <= 2.5 and 0 <= s[1] <= 2.5):
            return 1.0
        col = make_colour(n, a, Za, s=s)
        C = listen_ck(col, start_1s(col, a), nt, "learn", chunk_s=1e9)
        L = [x for x in lines(C[0], a, 0.80, 0.95) if x["weight"] >= 0.03]
        qq = quartet(L)
        cache[s] = qq["std"]
        d["evals"].append(dict(s=list(s), std=qq["std"], centroid=qq["centroid"], lines=L))
        save("learn", d)
        log("eval", s, qq)
        return qq["std"]
    from scipy.optimize import minimize
    x0 = cfg["start"]
    try:
        r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-3, fatol=2e-5, maxfev=cfg["maxfev"],
                                                               initial_simplex=np.array(cfg["simplex"])))
    except TimeoutError:
        log("time budget: run again"); return
    best = min(d["evals"], key=lambda e: e["std"])
    d["s"] = best["s"]; d["std"] = best["std"]
    bare = load().get("lab", {}).get("(0, 0, 0)")
    if bare:
        qb = quartet([x for x in bare["lines"] if 0.80 <= x["omega"] <= 0.95 and x["weight"] >= 0.03])
        d["bare_std"] = qb["std"]
        d["N1"] = dict(ratio=best["std"] / qb["std"], PASS=bool(best["std"] <= PRE["exams"]["N1"]["max_ratio"] * qb["std"]))
    log("learned", d.get("s"), d.get("N1"))
    save("learn", d)


def part_scan(a, cfgname):
    """1s multiplet (Laboratorio-like start: radial 1s of the same a) on the box of side 24, bare or learned"""
    key = f"scan_{cfgname}_{a}"
    if key in load():
        log("have", key); return
    box = PRE["config"]["scan"]["box"]
    n = int(round(box / a / 2)) * 2
    s = (1.0, 1.0) if cfgname == "bare" else learned_s()
    nt = int(round(PRE["config"]["scan"]["T"] / a))
    col = make_colour(n, a, Za, s=s)
    C = listen_ck(col, start_1s(col, a), nt, f"scan_{cfgname}_{a}")
    if C is None:
        return
    L = lines(C[0], a, 0.80, 0.95)
    qq = quartet([x for x in L if x["weight"] >= 0.03])
    rec = dict(a=a, n=n, s=s, nt=nt, lines=L, quartet=qq, C_re=C[0].real.tolist(), C_im=C[0].imag.tolist(), binding=1 - qq["centroid"], dirac_binding=1 - D1S,
               rel=(1 - qq["centroid"]) / (1 - D1S) - 1)
    log(key, rec["quartet"], f"rel {100 * rec['rel']:+.3f} %")
    save(key, rec)


def start_excited(col, a):
    """large components of radial 2s, 2p1/2 (m=1/2), 2p3/2 m=3/2 and m=1/2 (all along z) of the same a"""
    y = {}
    out = []
    for name, kappa, nl in (("2s", -1, 2), ("2p1/2", 1, 2), ("2p3/2", -2, 2)):
        om, r, g, f = radial_state(kappa, nl, a, Za, rmax=60.0)
        y[name] = (om, r, g)
    c11 = -np.sqrt(3 / (8 * np.pi)); c10 = np.sqrt(3 / (4 * np.pi))

    def radial(name, X, Y, Z):
        om, r, g = y[name]
        rr = np.maximum(np.sqrt(X**2 + Y**2 + Z**2), 1e-9)
        x = rr / a - 0.5
        i = np.clip(np.floor(np.maximum(x, 0)).astype(int), 0, len(g) - 2)
        t = np.clip(x - i, 0, 1)
        return (g[i] * (1 - t) + g[i + 1] * t) / rr, rr

    def mk(kind):
        def fun(X, Y, Z):
            if kind == "2s":
                G, rr = radial("2s", X, Y, Z)
                up = np.stack([G / np.sqrt(4 * np.pi), 0 * G], -1)
            elif kind == "2p1/2":
                G, rr = radial("2p1/2", X, Y, Z)
                y0 = 0.5 / np.sqrt(np.pi)
                up = np.stack([-y0 * Z / rr * G, -y0 * (X + 1j * Y) / rr * G], -1)
            elif kind == "2p3/2,3/2":
                G, rr = radial("2p3/2", X, Y, Z)
                up = np.stack([c11 * (X + 1j * Y) / rr * G, 0 * G], -1)
            else:
                G, rr = radial("2p3/2", X, Y, Z)
                up = np.stack([np.sqrt(2 / 3) * c10 * Z / rr * G, np.sqrt(1 / 3) * c11 * (X + 1j * Y) / rr * G], -1)
            s = 1 / np.sqrt(2)
            return np.concatenate([up * s, up * s], -1)     # chi = 0: L = R = phi / sqrt2
        return fun
    kinds = ["2s", "2p1/2", "2p3/2,3/2", "2p3/2,1/2"]
    for k in kinds:
        out.append(col.field(mk(k)))
    return kinds, out, {k: y[k][0] for k in y}


def part_excited(a, cfgname):
    key = f"exc_{cfgname}_{a}"
    if key in load():
        log("have", key); return
    box = PRE["config"]["excited"]["box"]
    n = int(round(box / a / 2)) * 2
    s = (1.0, 1.0) if cfgname == "bare" else learned_s()
    nt = int(round(PRE["config"]["excited"]["T"] / a))
    col = make_colour(n, a, Za, s=s)
    kinds, refs, omr = start_excited(col, a)
    ph = np.exp(2j * np.pi * np.array([0.0, 0.31, 0.57, 0.83]))
    E = sum(p * r[0] for p, r in zip(ph, refs)); O = sum(p * r[1] for p, r in zip(ph, refs))
    nr = np.sqrt(np.sum(np.abs(E)**2) + np.sum(np.abs(O)**2)); E /= nr; O /= nr
    C = listen_ck(col, (E, O), nt, f"exc_{cfgname}_{a}", refs=refs)
    if C is None:
        return
    rec = dict(a=a, n=n, s=s, nt=nt, radial_same_a=omr, C_re=C.real.tolist(), C_im=C.imag.tolist())
    for i, k in enumerate(kinds):
        L = lines(C[i], a, 0.955, 0.995, wmin=0.01)
        rec[k] = L
        log(key, k, [(round(x["omega"], 6), round(x["weight"], 3)) for x in L])
    save(key, rec)


def part_centre(a):
    """nucleus on a site (centre of a cell), amendment nucleo1s_3d_amend1.json: all 4 colours, bare rule"""
    AM = json.load(open(os.path.join(OUT, "nucleo1s_3d_amend1.json")))["centre"]
    key = f"centre_{a}"
    out = load().get(key, {})
    n = int(round(PRE["config"]["scan"]["box"] / a / 2)) * 2
    nt = int(round(AM["nt_T"] / a))
    cube = (3 * np.log(2 + np.sqrt(3)) - np.pi / 2) / a

    def Vf(rr):
        return np.where(rr < 1e-9 * a, -Za * cube, -Za / np.maximum(rr, 1e-9 * a))
    for q in ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)):
        if str(q) in out:
            continue
        col = Colour3D(n, a, Za, Vfun=Vf, q=q, origin=(a / 2, a / 2, a / 2))
        C = listen_ck(col, start_1s(col, a), nt, f"centre_{a}_{q[0]}{q[1]}{q[2]}")
        if C is None:
            save(key, out); return
        L = lines(C[0], a, 0.55, 0.97)
        strong = [x for x in L if x["weight"] >= 0.05 and 0.80 <= x["omega"] <= 0.95]
        L1s = [x for x in L if x["weight"] >= 0.03 and 0.80 <= x["omega"] <= 0.95]
        out[str(q)] = dict(lines=L, quartet=quartet(L1s) if L1s else None, C_re=C[0].real.tolist(), C_im=C[0].imag.tolist(),
                           split=float(max(x["omega"] for x in strong) - min(x["omega"] for x in strong)) if strong else None,
                           r_min=float(min(col.rE.min(), col.rO.min())))
        log(key, q, out[str(q)]["split"], [(round(x["omega"], 6), round(x["weight"], 3)) for x in L])
        save(key, out)


def part_glue():
    from nucleo1s_glue import glue_scan
    out = glue_scan(PRE["config"]["glue"], log)
    save("glue", out)


if __name__ == "__main__":
    p = sys.argv[1]
    if p == "lab":
        part_lab()
    elif p == "learn":
        part_learn()
    elif p == "scan":
        part_scan(float(sys.argv[2]), sys.argv[3])
    elif p == "excited":
        part_excited(float(sys.argv[2]), sys.argv[3])
    elif p == "centre":
        part_centre(float(sys.argv[2]))
    elif p == "glue":
        part_glue()
    log("done", p)
