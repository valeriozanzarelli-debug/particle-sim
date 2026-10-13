"""nucleo1s_f5_run.py - phase 5 runs (lab/out/nucleo1s_f5_prereg.json). Z-aware versions of the phase-3/4 parts.
Usage: timeout 590 python lab/nucleo1s_f5_run.py bare <Z> <a> | learn <Z> | use <Z> <a> <pred>
Results -> lab/out/nucleo1s_f5_results.json (checkpoints lab/out/nucleo1s_3d_ck_f5_*.npz, NOT committed)."""
import os, sys, json, time, fcntl
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nucleo1s_3d_run as B  # noqa: E402  (phase-2 helpers: make_colour, listen_ck, lines, quartet)
from nucleo1s_3d import radial_state, spinor_field  # noqa: E402
from zz_core import ALPHA, dirac_formula  # noqa: E402

OUT = B.OUT
RES = os.path.join(OUT, "nucleo1s_f5_results.json")
P5 = json.load(open(os.path.join(OUT, "nucleo1s_f5_prereg.json")))
_AM = os.path.join(OUT, "nucleo1s_f5_amend1.json")
AM = json.load(open(_AM)) if os.path.exists(_AM) else {}
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def load():
    return json.load(open(RES)) if os.path.exists(RES) else {}


def save(key, val):
    lk = open(RES + ".lock", "w")
    fcntl.flock(lk, fcntl.LOCK_EX)
    d = load(); d[key] = val
    json.dump(d, open(RES + ".tmp", "w"), indent=1, default=float); os.replace(RES + ".tmp", RES)
    fcntl.flock(lk, fcntl.LOCK_UN)


def window(Z):
    """listening window of the 1s: [1 - 2 B, 1 - 0.5 B], B = Dirac binding (only WHERE to listen; for Z = 60 it is
    [0.798, 0.950], the phase-2 window [0.80, 0.95])"""
    Bd = 1 - dirac_formula(1, -1, Z * ALPHA)
    return 1 - 2 * Bd, 1 - 0.5 * Bd


def box(Z):
    return P5["config"]["box"][str(Z)]


def start(col, a, Za):
    om_r, r, g, f = radial_state(-1, 1, a, Za, rmax=max(40.0, box(int(round(Za / ALPHA)))))
    return col.field(spinor_field(-1, r, g, f, a))


def U_of(L):
    s = sorted(L, key=lambda x: -x["omega"]); w = o = 0.0
    for x in s:
        w += x["weight"]; o += x["weight"] * x["omega"]
        if w >= 0.4:
            break
    return o / w


def scan(Z, a, s, key, T):
    if key in load():
        log("have", key); return
    Za = Z * ALPHA
    n = int(round(box(Z) / a / 2)) * 2
    nt = int(round(T / a))
    col = B.make_colour(n, a, Za, s=s)
    C = B.listen_ck(col, start(col, a, Za), nt, f"f5_{key}")
    if C is None:
        return
    lo, hi = window(Z)
    L = B.lines(C[0], a, lo, hi)
    qq = B.quartet([x for x in L if x["weight"] >= 0.03])
    D = dirac_formula(1, -1, Za)
    rec = dict(Z=Z, a=a, n=n, nt=nt, s=list(s), lines=L, quartet=qq, U=U_of(L), dirac=D,
               rel=(1 - qq["centroid"]) / (1 - D) - 1, rel_U=(1 - U_of(L)) / (1 - D) - 1)
    log(key, qq, f"rel {100 * rec['rel']:+.3f} %  U {100 * rec['rel_U']:+.3f} %")
    save(key, rec)


def part_learn(Z):
    cfg = P5["config"]["learn"]
    a, nt = cfg["a"], int(round(cfg["T"] / cfg["a"]))
    Za = Z * ALPHA
    n = int(round(box(Z) / a / 2)) * 2
    T = AM.get("bare_T", {}).get(str(Z), P5["config"]["bare_T"])
    U = load()[f"bare_Z{Z}_{a}" if T == P5["config"]["bare_T"] else f"bare_Z{Z}_{a}_T{int(T)}"]["U"]
    key = f"learn_Z{Z}"
    d = load().get(key, {"evals": [], "U_bare": U})
    if d.get("done"):
        log("have", key); return
    cache = {tuple(e["s"]): e["obj"] for e in d["evals"]}
    lo, hi = window(Z)

    def f(s):
        s = (round(float(s[0]), 6), round(float(s[1]), 6))
        if s in cache:
            return cache[s]
        if time.time() - T0 > 500:
            raise TimeoutError
        if not (0 <= s[0] <= 2.5 and 0 <= s[1] <= 2.5):
            return 1.0
        col = B.make_colour(n, a, Za, s=s)
        C = B.listen_ck(col, start(col, a, Za), nt, f"f5_{key}", chunk_s=1e9)
        L = [x for x in B.lines(C[0], a, lo, hi) if x["weight"] >= 0.03]
        qq = B.quartet(L)
        obj = qq["std"]**2 + (qq["centroid"] - U)**2
        cache[s] = obj
        d["evals"].append(dict(s=list(s), obj=obj, std=qq["std"], centroid=qq["centroid"], lines=L))
        save(key, d)
        log("eval", Z, s, qq, "obj", obj)
        return obj
    from scipy.optimize import minimize
    try:
        minimize(f, cfg["start"], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-12, maxfev=cfg["maxfev"],
                                                                     initial_simplex=np.array(cfg["simplex"])))
    except TimeoutError:
        log("time budget: run again"); return
    best = min(d["evals"], key=lambda e: e["obj"])
    d.update(s=best["s"], std=best["std"], centroid=best["centroid"], dU=best["centroid"] - U, done=True)
    log("learned", Z, d["s"], d["std"], d["dU"])
    save(key, d)


if __name__ == "__main__":
    p = sys.argv[1]
    if p == "bare":
        Z, a = int(sys.argv[2]), float(sys.argv[3])
        T = AM.get("bare_T", {}).get(str(Z), P5["config"]["bare_T"])
        key = f"bare_Z{Z}_{a}" if T == P5["config"]["bare_T"] else f"bare_Z{Z}_{a}_T{int(T)}"
        scan(Z, a, (1.0, 1.0), key, T)
    elif p == "learn":
        part_learn(int(sys.argv[2]))
    elif p == "use":
        Z, a, pred = int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
        s = P5["predictions"][pred][f"Z{Z}_a{a}"]
        scan(Z, a, tuple(s), f"use_{pred}_Z{Z}_{a}", P5["config"]["use_T"])
    log("done", p)
