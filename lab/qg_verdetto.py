"""qg_verdetto.py - grades of the pre-registered exams (lab/out/qg_prereg.json) -> lab/out/qg_verdetto.json.

    python3 lab/qg_verdetto.py [parte1] [parte2]
Part 1 reads <SCR>/qg_prod_analisi.json (lab/qg_analisi1.py prod); part 2 reads the run diagnostics and fields.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qg_run import SCR  # noqa: E402

LAB = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(LAB, "out", "qg_verdetto.json")


def parte1():
    an = json.load(open(os.path.join(SCR, "qg_prod_analisi.json")))
    v = {"n_misure": an["n_meas"]}
    p, pe = an["plaq"]
    v["A_placchetta"] = {"noi": [p, pe], "letteratura": 0.5492, "scarto": p - 0.5492,
                         "voto": "PASS" if abs(p - 0.5492) <= 0.002 else "FAIL",
                         "nota": f"{p:.4f}({pe * 1e4:.0f}) contro 0,5492"}
    a, ae = an["cornell"]["a_sqrt_sigma"]
    ok = abs(a - 0.40) <= 0.04 + 2 * ae
    v["B_tensione"] = {"noi": [a, ae], "letteratura": 0.40, "scarto_rel": (a - 0.40) / 0.40,
                       "cornell": an["cornell"], "voto": "PASS" if ok else "FAIL",
                       "nota": f"a√σ = {a:.3f}({ae * 1e3:.0f}) contro 0,40"}
    fl = an["flux"]["4_4"]
    r, re = fl["R_Y"]
    fc, fce = fl["f"]["centroide"]
    if r - 1 > 2 * re and fc > 3 * fce:
        vt = "PASS (Y)"
    elif 1 - r > 2 * re:
        vt = "FAIL (Delta)"
    else:
        vt = "FAIL (non concluso)"
    v["C_Y_contro_Delta"] = {"R_Y": [r, re], "f_centroide": [fc, fce], "f_medi_lati": fl["f"]["medi_lati"],
                             "f_origine": fl["f"]["origine"], "voto": vt,
                             "nota": f"R_Y = {r:.2f} ± {re:.2f} (d=4, T=4)",
                             "riportati": {k: {"R_Y": an["flux"][k]["R_Y"], "f": an["flux"][k]["f"],
                                               "eps": an["flux"][k]["eps"]} for k in an["flux"]}}
    E = {}
    for d, e in an["E3Q"].items():
        i = e["T"].index(3)
        E[d] = [e["E"][i], e["err"][i]]
    v["D_energia_3Q"] = {"E3Q_T3": E, "voto": "riportato"}
    ds = sorted(E, key=int)
    if len(ds) >= 2:
        d0, d1 = ds[0], ds[-1]
        slope = (E[d1][0] - E[d0][0]) / ((int(d1) - int(d0)) * np.sqrt(6))
        v["D_energia_3Q"]["pendenza_per_lunghezza_Y"] = slope
        v["D_energia_3Q"]["sigma_qqbar"] = an["cornell"]["sigma"]
    return v


def diag(name):
    return json.load(open(os.path.join(SCR, f"qg_ym_{name}_diag.json")))


def fine(name):
    z = np.load(os.path.join(SCR, f"qg_ym_{name}_fine.npz"))
    return (z["eE"] + z["eB"]).astype(float)


def D(name, ref, amp, amp_ref=0.005):
    e = fine(name)
    r = fine(ref) * (amp / amp_ref) ** 2
    return float(np.abs(e - r).sum() / e.sum())


def parte2():
    pr = json.load(open(os.path.join(LAB, "out", "qg_prereg.json")))["parte2"]["corse"]
    v = {}
    g1, g2 = {}, {}
    for name in pr:
        try:
            dg = diag(name)
        except FileNotFoundError:
            continue
        H = np.array(dg["H"])
        ga = np.array(dg["gauss_abs"])
        g1[name] = {"iniziale_rel": dg["gauss_rel"][0], "max_rel": max(dg["gauss_rel"]),
                    "deriva_rel": float(np.abs(ga - ga[0]).max() / max(ga[0], 1e-300))}
        g2[name] = float(np.abs(H - H[0]).max() / H[0])
    okg = all(x["iniziale_rel"] <= 1e-3 and (x["deriva_rel"] <= 1e-3 or x["max_rel"] < 1e-12) for x in g1.values())
    v["G1_gauss"] = {"corse": g1, "voto": "PASS" if okg else "FAIL",
                     "nota": "per il colore uguale il residuo e' ~1e-15 (zero di macchina): li' la deriva relativa "
                             "non ha senso e conta il valore assoluto"}
    v["G2_energia"] = {"corse": g2, "voto": "PASS" if max(g2.values()) <= 1e-3 else "FAIL"}
    amps = {"a002": 0.02, "a005": 0.05, "a01": 0.1, "a03": 0.3, "a1": 1.0}
    Ds = {k: D(k, "rif", a) for k, a in amps.items() if os.path.exists(os.path.join(SCR, f"qg_ym_{k}_fine.npz"))}
    v["D"] = Ds
    v["G3_lineare"] = {"D_0.02": Ds.get("a002"), "voto": "PASS" if Ds.get("a002", 1) <= 0.01 else "FAIL"}
    v["G4_non_lineare"] = {"D_1.0": Ds.get("a1"), "voto": "PASS" if Ds.get("a1", 0) >= 0.10 else "FAIL"}
    ds = D("a1_stesso", "rif_stesso", 1.0) if os.path.exists(os.path.join(SCR, "qg_ym_a1_stesso_fine.npz")) else None
    v["G5_colore"] = {"D_stesso_1.0": ds, "D_1.0": Ds.get("a1"),
                      "voto": "PASS" if ds is not None and ds <= Ds.get("a1", 0) / 3 else "FAIL"}
    if all(k in Ds for k in ("a002", "a005", "a01")):
        x = np.log([0.02, 0.05, 0.1])
        y = np.log([Ds["a002"], Ds["a005"], Ds["a01"]])
        sl = float(np.polyfit(x, y, 1)[0])
        v["G6_scala"] = {"pendenza": sl, "previsione": 2, "voto": "PASS" if 1.6 <= sl <= 2.4 else "FAIL"}
    fid = {}
    for name in ("video_piccolo", "video_grande"):
        try:
            dg = diag(name)
            fid[name] = {"min_r": min(dg["fid_r"]), "n_fotogrammi": len(dg["fid_r"])}
        except (FileNotFoundError, KeyError, ValueError):
            pass
    if fid:
        v["G7_palline"] = {"corse": fid, "voto": "PASS" if len(fid) == 2 and all(f["min_r"] >= 0.9 for f in fid.values()) else "FAIL"}
    return v


if __name__ == "__main__":
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    if "parte1" in sys.argv:
        out["parte1"] = parte1()
    if "parte2" in sys.argv:
        out["parte2"] = parte2()
    json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
    for k, part in out.items():
        for e, r in part.items():
            if isinstance(r, dict) and "voto" in r:
                print(k, e, r["voto"], r.get("nota", ""))
