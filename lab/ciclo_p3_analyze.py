"""ciclo_p3_analyze.py - verdicts of phase 3 (exams B1-B6 of lab/out/ciclo_prereg.json, fase3) -> lab/out/ciclo_p3.json"""
import os, json
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
VMAX, CONE0, B = 0.029, 36.0, 9


def J(n):
    return json.load(open(os.path.join(OUT, f"ciclo_p3_{n}.json")))


def nref(d):
    return int(np.sum(d["refined"]))


def cone(d):
    worst = -1e9
    for ev in d["events"]:
        for I in ev["cells"]:
            xc = (I + 0.5) * B
            worst = max(worst, abs(xc - d["xd"]) - (CONE0 + VMAX * ev["t"]))
    return worst                     # <= 0: inside the cone


def conserv(d):
    nn = [h["norm"] for h in d["hist"]]
    E = [h["E_own"] for h in d["hist"] if h["t"] > d["Tp"]]
    return dict(norm_drift=float(max(abs(np.array(nn) - 1))), norm_loss_refine=d["norm_loss"],
                E_rel_after=float((max(E) - min(E)) / abs(np.mean(E))))


def main():
    res = dict(exams={})
    q, f = J("Q_mix"), J("F_mix")
    res["exams"]["B1"] = dict(value=dict(refined=nref(q), fid_min=q["fid_min"], floor=max(h["emax"] for h in q["hist"])),
                              ok=nref(q) == 0 and q["fid_min"] >= 0.999)
    b2, b3, b6 = {}, {}, {}
    for s in ("P", "P2"):
        m, c, u = J(f"{s}_mix"), J(f"{s}_coarse"), J(f"{s}_b3")
        tr = J(f"{s}_truth")
        b2[s] = dict(cone_worst=cone(m), frac=nref(m) / 243, ok=cone(m) <= 0 and nref(m) / 243 <= 0.5)
        po, pt, pc = m["P_out"][-1], m["P_out_truth"][-1], c["P_out"][-1]
        b3[s] = dict(P_out=po, P_out_truth=pt, rel=abs(po - pt) / pt, fid_min=m["fid_min"], coarse_rel=abs(pc - pt) / pt,
                     coarse_fid=c["fid_min"], b3_rel=abs(u["P_out"][-1] - pt) / pt, b3_fid=u["fid_min"],
                     ok=abs(po - pt) / pt <= 0.05 and m["fid_min"] >= 0.99 and abs(pc - pt) / pt > 0.3)
        unif = [("B=9", c), ("B=3", u), ("base", dict(fid_min=1.0, cost=tr["cost"]))]
        ok_u = [(n, d["cost"]) for n, d in unif if d["fid_min"] >= m["fid_min"]]
        nm, cu = min(ok_u, key=lambda z: z[1])
        b6[s] = dict(cheapest_uniform=nm, cost_uniform=cu, cost_mix=m["cost"], cost_mix_steps=m["cost_steps"],
                     cost_mix_check=m["cost_check"], saving=cu / m["cost"], ok=cu / m["cost"] >= 1.5)
    res["exams"]["B2"] = dict(value=b2, ok=all(v["ok"] for v in b2.values()))
    res["exams"]["B3"] = dict(value=b3, ok=all(v["ok"] for v in b3.values()))
    b4 = {s: conserv(J(f"{s}_mix")) for s in ("Q", "P", "P2", "F")}
    res["exams"]["B4"] = dict(value=b4, ok=all(v["norm_drift"] <= 1e-8 and v["E_rel_after"] <= 1e-6 for v in b4.values()))
    res["exams"]["B5"] = dict(value=dict(refined=nref(f), fid_min=f["fid_min"]), ok=nref(f) == 0 and f["fid_min"] >= 0.999)
    res["exams"]["B6"] = dict(value=b6, ok=all(v["ok"] for v in b6.values()))
    json.dump(res, open(os.path.join(OUT, "ciclo_p3.json"), "w"), indent=1, default=float)
    for k, v in res["exams"].items():
        print(f"[{'PASS' if v['ok'] else 'FAIL'}] {k}: {json.dumps(v['value'], default=float)[:400]}")


if __name__ == "__main__":
    main()
