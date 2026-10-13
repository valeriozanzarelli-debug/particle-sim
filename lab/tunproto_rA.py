"""tunproto_rA.py - corse A e B: la misura con l'ancella ridà tau_m (regola tick per tick e Dirac continuo),
bordi lisci, delta finito. Uso: python3 tunproto_rA.py regola | dirac | bordi   (un pezzo per chiamata, checkpoint JSON)."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunproto_core import OUT, rule_run, spectral_rule, dirac_run

PRE = json.load(open(os.path.join(OUT, "tunproto_prereg.json")))


def piece_rule():
    A = PRE["A_misura"]
    res = {}
    for a in A["passi"]:
        for V0 in A["V0"]:
            t0 = time.time()
            r = rule_run(a, A["E0"], V0, A["L"], A["sigma"], A["delta"], "mV")
            sp = {c: spectral_rule(a, A["E0"], V0, A["L"], A["sigma"], A["delta"], c) for c in ("m", "V")}
            res["a%g_V%.1f" % (a, V0)] = dict(run=r, spett=sp, sec=time.time() - t0)
            print(a, V0, "m", r["m"]["tz"], sp["m"]["tz"], r["m"]["ty"], sp["m"]["ty"], "%.1fs" % (time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(OUT, "tunproto_rA_regola.json"), "w"), indent=1)


def piece_dirac():
    A = PRE["A_misura"]
    res = {}
    for V0 in A["V0"]:
        t0 = time.time()
        r = dirac_run(A["E0"], V0, A["L"], A["sigma"], A["delta"], "m", dx=A["dirac_split"]["dx"], dt=A["dirac_split"]["dt"])
        res["V%.1f" % V0] = r
        print(V0, r, "%.1fs" % (time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(OUT, "tunproto_rA_dirac.json"), "w"), indent=1)


def piece_bordi():
    B = PRE["B_sistematici"]["bordi_lisci"]
    a = B["a"]
    res = {}
    for w in B["w"][1:]:
        for prof in (False, True):
            key = "w%.2f_%s" % (w, "profilo" if prof else "rettangolo")
            t0 = time.time()
            r = rule_run(a, 2.0, 2.0, 4.0, B["sigma"], B["delta"], "mV", w=w, clock_profile=prof)
            res[key] = r
            print(key, r["m"]["tz"], B["val"][key]["spettrale_L4"]["m"]["tz"], "%.1fs" % (time.time() - t0), flush=True)
    # delta finito sulla regola tick per tick (a = 0,025), V0 = 2, L = 4
    for d in (0.2, 0.8):
        r = rule_run(a, 2.0, 2.0, 4.0, 10.0, d, "m")
        sp = spectral_rule(a, 2.0, 2.0, 4.0, 10.0, d, "m")
        sp0 = spectral_rule(a, 2.0, 2.0, 4.0, 10.0, 1e-3, "m")
        res["delta%.1f" % d] = dict(run=r, spett=sp, spett_debole=sp0)
        print("delta", d, r["m"]["tz"], sp["tz"], sp0["tz"], flush=True)
    json.dump(res, open(os.path.join(OUT, "tunproto_rA_bordi.json"), "w"), indent=1)


if __name__ == "__main__":
    {"regola": piece_rule, "dirac": piece_dirac, "bordi": piece_bordi}[sys.argv[1]]()
