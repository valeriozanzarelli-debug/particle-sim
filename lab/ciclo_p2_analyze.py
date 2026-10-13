"""ciclo_p2_analyze.py - verdicts of phase 2 (exams A1-A8 of lab/out/ciclo_prereg.json, fase2) -> lab/out/ciclo_p2.json"""
import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def J(n):
    return json.load(open(os.path.join(OUT, n)))


def loginterp(Bs, ys, B, signed=False):
    """value at B by linear interpolation in log B (log |y| unless signed); None if B is outside the range"""
    Bs = np.array(Bs, float); ys = np.array(ys, float)
    if B < Bs.min() - 1e-9 or B > Bs.max() + 1e-9:
        return None
    if signed:
        return float(np.interp(np.log(B), np.log(Bs), ys))
    return float(np.exp(np.interp(np.log(B), np.log(Bs), np.log(np.abs(ys)))))


def taps_at(levels, B, key="B"):
    """taps of the level nearest (in log B) to B"""
    i = int(np.argmin([abs(np.log(l[key] / B)) for l in levels]))
    return levels[i]


def verdict(scores):
    """phi wins / loses / tie against the best of {2, 3}"""
    out = []
    for s in scores:
        p, b2, b3 = s["phi"], s["2"], s["3"]
        best = min(x for x in (b2, b3) if x is not None)
        if p is None:
            out.append("LOSES (does not reach B)")
        elif p <= 0.8 * best:
            out.append("WINS")
        elif p >= 1.2 * best:
            out.append("LOSES")
        else:
            out.append("TIE")
    if all(o == "WINS" for o in out):
        return "WINS", out
    if any(o.startswith("LOSES") for o in out):
        return "LOSES", out
    return "TIE", out


def main():
    res = dict(exams={})
    d1 = {(m, s): J(f"ciclo_p2_1d_{m}_{s}.json") for m in ("free", "coulomb") for s in ("2", "3", "fib")}
    # A1
    worst, Cs = 1.0, set()
    for (m, s), d in d1.items():
        for l in d["levels"]:
            if "flag" in l:
                continue
            worst = min(worst, min(l["hand_overlap"].values()))
            Cs.add(l["C"])
    res["exams"]["A1"] = dict(value=dict(min_overlap=worst, C_used=sorted(Cs)), ok=worst >= 0.99 and Cs == {1})
    # A2
    a2 = J("ciclo_p2_3d_A2.json")
    res["exams"]["A2"] = dict(value=a2["overlap"], eigs=a2["eigs"], ok=a2["overlap"] >= 0.999)
    # A3 / A4
    def chain_ok(d):
        adm = [l for l in d["levels"] if "flag" not in l]
        allok = all("flag" not in l for l in d["levels"])
        return dict(last_B=adm[-1]["B"] if adm else 0, acc_last=adm[-1].get("acc_err") if adm else None,
                    all_admitted=allok, ok=bool(allok and adm and adm[-1]["B"] >= 27 and adm[-1]["acc_err"] <= 1e-2))
    a3 = {s: chain_ok(d1[("free", s)]) for s in ("2", "3")}
    res["exams"]["A3"] = dict(value=a3, ok=all(v["ok"] for v in a3.values()))
    a4 = chain_ok(d1[("free", "fib")])
    res["exams"]["A4"] = dict(value=a4, ok=a4["ok"])
    # A5
    scores = []
    for Bc in (8.0, 16.0):
        sc = {}
        for s in ("2", "3", "fib"):
            adm = [l for l in d1[("free", s)]["levels"] if "flag" not in l]
            v = loginterp([l["B"] for l in adm], [l["acc_err"] for l in adm], Bc) if adm else None
            if v is None:
                sc["phi" if s == "fib" else s] = None
                continue
            l = taps_at(adm, Bc)
            sc["phi" if s == "fib" else s] = v * (2 * l["w"] + 1) * l["C"] ** 2
        sc["B"] = Bc
        scores.append(sc)
    v5, per = verdict(scores)
    res["exams"]["A5"] = dict(value=scores, per_B=per, outcome=v5, prediction="LOSES", ok=(v5 == "LOSES"))
    # A6
    a6 = {}
    for s in ("2", "3", "fib"):
        adm = [l for l in d1[("coulomb", s)]["levels"] if "flag" not in l]
        if not adm or max(l["B"] for l in adm) < 4 / 1.5:
            a6[s] = dict(reached=False)
            continue
        l = taps_at(adm, 4.0)
        a6[s] = dict(reached=True, B=l["B"], rel=l["bound_rel_err"], ok=max(l["bound_rel_err"]) <= 1e-2,
                     score=max(l["bound_rel_err"]) * (2 * l["w"] + 1))
    res["exams"]["A6"] = dict(value=a6, ok=all(v["ok"] for v in a6.values() if v.get("reached")))
    # A7
    d3 = {s: J(f"ciclo_p2_3d_{s}.json") for s in ("2", "3", "phi")}
    a7, sc7, sc7acc = {}, {}, {}
    for s, d in d3.items():
        adm = [l for l in d["levels"] if "flag" not in l]
        reach = max(l["B"] for l in adm) >= 400
        dE = loginterp([l["B"] for l in adm], [l["dE1s_rel"] for l in adm], 400.0, signed=True)
        acc = loginterp([l["B"] for l in adm], [l["acc_rel"] for l in adm], 400.0)
        l = taps_at(adm, 400.0)
        a7[s] = dict(reached=reach, dE1s_rel_400=dE, acc_rel_400=acc, taps=l["taps"], levels_to_400=sum(1 for x in adm if 16 < x["B"] <= 470),
                     flags=[x.get("flag") for x in d["levels"] if "flag" in x], ok=bool(reach and dE is not None and abs(dE) <= 1e-3))
        sc7[s] = None if dE is None else abs(dE) * l["taps"]
        sc7acc[s] = None if acc is None else acc * l["taps"]
    v7, per7 = verdict([sc7])
    v7a, _ = verdict([sc7acc])
    res["exams"]["A7"] = dict(value=a7, score=sc7, outcome=v7, score_acc_rel_reported=sc7acc, outcome_acc_rel=v7a,
                             ok=all(v["ok"] for v in a7.values()))
    # A8
    an = J("ciclo_p2_1d_anchor.json")
    res["exams"]["A8"] = dict(value=an["rel"], ok=max(an["rel"]) <= 1e-2)
    # per-decade numbers (reported)
    dec = {}
    for s, d in d3.items():
        adm = [l for l in d["levels"] if "flag" not in l]
        B = np.array([l["B"] for l in adm]); a = np.array([l["acc_rel"] for l in adm])
        dec[s] = dict(levels_per_decade=1 / np.log10(d["b"]), acc_growth_per_decade=float(np.polyfit(np.log10(B[1:]), np.log10(a[1:]), 1)[0]))
    res["per_decade_3d"] = dec
    json.dump(res, open(os.path.join(OUT, "ciclo_p2.json"), "w"), indent=1, default=float)
    for k, v in res["exams"].items():
        print(f"[{'PASS' if v['ok'] else 'FAIL'}] {k}: {json.dumps(v.get('value'), default=float)[:300]} {v.get('outcome', '')}")
    print(json.dumps(dec))


if __name__ == "__main__":
    main()
