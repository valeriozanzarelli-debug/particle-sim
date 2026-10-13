"""fw1atomo_p3.py - P3 of the pre-registration: the 3-D atom (equivalent of lab/zz_atom3d.py, ZIGZAG A1) with FW1 and
with the frozen rule.  Z = 60, nucleus on the corner of 8 cells, box ~16 lambda_C, a = 0,4 / 0,3 / 0,24 / 0,2.
Steps:  ritz_<rule>_<start>_<a>   Ritz in [0,78; 0,96] (dw 0,01, 600 lambda_C/c) from start 'dirac' (continuum
                                   Dirac 1s, spin up) or 'cieca' (random spinors on the 4^3 sites around the nucleus)
        v0                         validation: F0 with the zz_atom3d protocol at a = 0,4, box 24 (projected blob,
                                   Gaussian lock-in std 27, Rayleigh quotient) -> 0,130043
Run: OMP_NUM_THREADS=1 timeout 590 python lab/fw1atomo_p3.py <step>  (exit 3 = checkpoint); list; sum"""
import os
import sys
import json
import time
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import fw1atomo_core as fc
import zz_core as zz

T0 = time.time()
OUT = os.path.join(os.path.dirname(__file__), "out", "fw1atomo_p3.json")
RES = os.path.join(fc.SCR, "p3")
os.makedirs(RES, exist_ok=True)
AS = {0.4: 40, 0.3: 54, 0.24: 66, 0.2: 80}
ZA = fc.ZA60
B_DIRAC = 1 - zz.dirac_formula(1, -1, ZA)
STEPS = ["v0"] + [f"ritz_{r}_{s}_{a}_c{c}" for a in AS for r in ("FW1", "F0") for s in ("dirac", "cieca") for c in ((0, 1) if s == "dirac" else (0,))]


def run_ritz(step):
    _, rule, start, a, c = step.split("_")
    a, c = float(a), int(c[1:])
    L = fc.Lattice(AS[a], a, "angolo")
    R = L.rule(rule, 0.4)
    S = L.dirac1s((1, 0)) if start == "dirac" else L.blind()
    S = S * (L.cls == c)[None]                    # amendment A3: one period-sublattice (U2 is block diagonal)
    S = [S / np.linalg.norm(S)]
    tp = time.time()
    for _ in range(3):
        R.period(S[0])
    t_per = (time.time() - tp) / 3
    lv = fc.ritz(R, S, a, 0.78, 0.96, 0.01, int(round(600 / (2 * a))), ck=os.path.join(fc.SCR, "p3_" + step), tmax=540, t0=T0)
    out = []
    for o in lv:
        if o["resid"] < 1e-4 and o["w_start"][0] >= 1e-4:
            d = L.diag(o["v"])
            d.update(omega=o["omega"], resid=o["resid"], w=o["w_start"][0])
            out.append(d)
    return {"levels": out, "s_per_period": t_per, "n": AS[a], "a": a}


def run_v0():
    """zz_atom3d.run(0.4) with the rules of fw1atomo_core (F0): same blob, projector, window, Rayleigh"""
    a, half, S_T = 0.4, 12.0, 27.0
    n = int(round(2 * half / a))
    L = fc.Lattice(n, a, "angolo")
    R = L.rule("F0")
    r = L.r
    psi = np.exp(-r * ZA)[None] * (np.array([1, 0, 1, 0], complex) / np.sqrt(2))[:, None, None, None]
    pk = np.fft.fftn(psi, axes=(1, 2, 3)).reshape(4, -1)
    kk = 2 * np.pi * np.fft.fftfreq(n)
    K = np.stack([q.ravel() for q in np.meshgrid(kk, kk, kk, indexing="ij")], -1)
    out = np.zeros_like(pk)
    for b0 in range(0, K.shape[0], 200000):
        sl = slice(b0, b0 + 200000)
        w, v = np.linalg.eig(zz.bloch_batch(K[sl], a, 0.5))
        om = -np.angle(w)
        o = np.argsort(-om, axis=1)[:, :2]
        Q, _ = np.linalg.qr(np.take_along_axis(v, o[:, None, :], axis=2))
        c = np.einsum("nab,na->nb", Q.conj(), pk[:, sl].T)
        out[:, sl] = np.einsum("nab,nb->na", Q, c).T
    psi = np.fft.ifftn(out.reshape(4, n, n, n), axes=(1, 2, 3))
    psi /= np.sqrt(np.sum(np.abs(psi)**2))

    def rayleigh(p):
        q = R.tick(R.tick(p.copy(), 0), 1)
        return -np.angle(np.vdot(p, q) / np.vdot(p, p)) / (2 * a)
    om0 = rayleigh(psi)
    nt = int(round(4.5 * S_T / a))
    nt += nt % 2
    c0, s = nt / 2, S_T / a
    acc = np.zeros_like(psi)
    for t in range(nt + 1):
        acc += np.exp(-0.5 * ((t - c0) / s)**2) * np.exp(1j * om0 * a * t) * psi
        psi = R.tick(psi, t)
    acc /= np.linalg.norm(acc)
    om1 = rayleigh(acc)
    return {"omega_start": float(om0), "omega": float(om1), "binding": float(1 - om1), "ref_zz_atom3d": 0.13004250114978633}


def run(step):
    fn = os.path.join(RES, step + ".json")
    if os.path.exists(fn):
        print("done already", step)
        return
    res = run_v0() if step == "v0" else run_ritz(step)
    json.dump(res, open(fn, "w"), indent=1)
    if step != "v0":
        for d in res["levels"]:
            if d["resid"] < 1e-6 and d["w"] > 2e-3:
                print(f"{d['omega']:.5f} res {d['resid']:.1e} w {d['w']:.4f} G {d['val_Gamma']:.3f} N {d['val_N']:.3f} "
                      f"pi {d['val_pi']:.3f} r {d['r_mean']:.2f} S {np.round(d['spin'], 2)}", flush=True)
    else:
        print(res)
    print("DONE", step, time.time() - T0)


def groups_of(lv, wmin_rel=0.05, lo=0.78, hi=0.95, tol=2e-3):
    ex = [l for l in lv if l["resid"] < 1e-6 and lo <= l["omega"] <= hi]
    tot = sum(l["w"] for l in ex)
    ex = sorted([l for l in ex if l["w"] >= wmin_rel * tot], key=lambda l: l["omega"])
    gr = []
    for l in ex:
        if gr and l["omega"] - gr[-1][-1]["omega"] <= tol:
            gr[-1].append(l)
        else:
            gr.append([l])
    return [{"omega": [round(l["omega"], 5) for l in g], "w_rel": round(sum(l["w"] for l in g) / tot, 4),
             "Gamma": round(float(np.mean([l["val_Gamma"] for l in g])), 3), "r": round(float(np.mean([l["r_mean"] for l in g])), 2)} for g in gr]


def summary():
    res = {"B_dirac": B_DIRAC, "per_a": {}}
    fn = os.path.join(RES, "v0.json")
    if os.path.exists(fn):
        v0 = json.load(open(fn))
        res["V0"] = dict(v0, PASS=bool(abs(v0["binding"] - v0["ref_zz_atom3d"]) <= 2e-4))
    for a in AS:
        W = zz.RadialWalk(-1, a, 40.0, ZA)
        om_rad = float(W.eigen(1 - 0.5 * ZA**2 * 1.2)[0])
        e = {"n": AS[a], "radial_omega": om_rad, "radial_B": 1 - om_rad}
        for rule in ("FW1", "F0"):
            for st in ("dirac_0.%s_c0", "dirac_0.%s_c1", "cieca_0.%s_c0"):
                key = f"ritz_{rule}_{st.split('_')[0]}_{a}_{st.split('_')[-1]}"
                f = os.path.join(RES, key + ".json")
                if not os.path.exists(f):
                    continue
                d = json.load(open(f))
                lv = [l for l in d["levels"] if l["resid"] < 1e-6]
                g = groups_of(lv)
                r = {"groups": g, "s_per_period": d["s_per_period"]}
                if lv and st.startswith("dirac"):
                    m = max(lv, key=lambda l: l["w"])
                    r.update(main_omega=m["omega"], main_B=1 - m["omega"], main_w=m["w"], main_Gamma=m["val_Gamma"], main_r=m["r_mean"],
                             main_spin=m["spin"])
                    near = [l for l in lv if abs(l["omega"] - m["omega"]) < 0.03 and l is not m and l["w"] > 0.01]
                    # spin partner (pre-registration: "coppia di spin mediata"): nearest level with opposite spin
                    part = [l for l in near if np.dot(l["spin"], m["spin"]) < 0]
                    if part:
                        q = min(part, key=lambda l: abs(l["omega"] - m["omega"]))
                        r.update(partner_omega=q["omega"], split=abs(q["omega"] - m["omega"]),
                                 pair_omega=0.5 * (q["omega"] + m["omega"]))
                        r["main_B"] = 1 - r["pair_omega"]
                    r["others_near"] = [(round(l["omega"], 5), round(l["w"], 3), round(l["val_Gamma"], 3)) for l in near]
                e[key] = r
        res["per_a"][str(a)] = e
    # extrapolation of FW1 binding (sublattice 0) and of the radial rule
    aa = np.array([a for a in AS if f"ritz_FW1_dirac_{a}_c0" in res["per_a"][str(a)]])
    if len(aa) >= 3:
        B = np.array([res["per_a"][str(a)][f"ritz_FW1_dirac_{a}_c0"]["main_B"] for a in aa])
        gam = np.sqrt(1 - ZA**2)
        fin = np.argsort(aa)[:3]
        for nm, p in (("a2gamma", 2 * gam), ("a2", 2.0)):
            X = np.vstack([np.ones(3), aa[fin]**p]).T
            c = np.linalg.lstsq(X, B[fin], rcond=None)[0]
            res[f"extrap_FW1_{nm}"] = {"B0": float(c[0]), "rel": float(c[0] / B_DIRAC - 1)}
        res["A1"] = {"PASS": abs(res["extrap_FW1_a2gamma"]["rel"]) <= 0.01}
        res["A1b"] = {"PASS": abs(res["extrap_FW1_a2gamma"]["rel"]) <= 0.03}
        res["A1c"] = {str(a): float(res["per_a"][str(a)][f"ritz_FW1_dirac_{a}_c0"]["main_B"] / res["per_a"][str(a)]["radial_B"] - 1) for a in aa}
    A3 = {}
    for a in AS:
        e = res["per_a"][str(a)]
        for rule in ("FW1", "F0"):
            k = f"ritz_{rule}_cieca_{a}_c0"
            if k in e:
                A3[f"{rule}_{a}"] = len(e[k]["groups"])
    res["A3_groups"] = A3
    json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "per_a"}, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "list":
        print(" ".join(STEPS))
    elif sys.argv[1] == "sum":
        summary()
    else:
        run(sys.argv[1])

