"""tavg_f5_judge.py - FASE 5: verdetti come in lab/out/tavg_prereg_f5.json (C1 NON primo cieco: dichiarato).
Uso: python tavg_f5_judge.py  -> lab/out/tavg_f5_verdetti.json"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_blocks as tb
import tavg_corr as tc

OUT = os.path.join(HERE, "out"); HA_EV = 27.211386245988
R = json.load(open(os.path.join(OUT, "tavg_f5_runs.json")))
XC = (4, 5) if "--em1" in sys.argv else (3, 4)       # emendamento 1: core da X = 4, 5
FULL = {"O": (8, "1s2 2s2 2p4"), "O+": (8, "1s2 2s2 2p3"), "N": (7, "1s2 2s2 2p3"), "N+": (7, "1s2 2s2 2p2")}


def ec(key):
    return R[key]["Ec_val"]          # E_tot - E_ROHF (fine level), all electrons in the reference


def cbs(a, b, Ea, Eb):
    return (b ** 3 * Eb - a ** 3 * Ea) / (b ** 3 - a ** 3)


def main():
    V = {"nota": "seconda prova dopo il C1 della fase 4: verdetto onesto ma NON primo cieco"}
    Ehe = tc._ecg("he1S_Z8")[0]; ex = Ehe - tb.get_ae(8, "1s2").E
    a_, b_ = XC
    k = cbs(a_, b_, ec(f"He_Z8_AEc_X{a_}"), ec(f"He_Z8_AEc_X{b_}"))
    V["XC"] = XC
    V["K5"] = dict(Ec_Xa=ec(f"He_Z8_AEc_X{a_}"), Ec_Xb=ec(f"He_Z8_AEc_X{b_}"), Ec_cbs=k, esatta=ex, diff=k - ex,
                   esito="PASS" if abs(k - ex) <= 5e-4 else "FAIL")
    pr = {}
    for sp, (Z, cfg) in FULL.items():
        val = cbs(4, 5, ec(f"{sp}_FCv_X4"), ec(f"{sp}_FCv_X5"))
        d3 = ec(f"{sp}_AEc_X{a_}") - ec(f"{sp}_FCc_X{a_}"); d4 = ec(f"{sp}_AEc_X{b_}") - ec(f"{sp}_FCc_X{b_}")
        core = cbs(a_, b_, d3, d4)
        e1s = tc._ecg(f"he1S_Z{Z}")[0] - tb.get_ae(Z, "1s2").E
        pr[sp] = dict(E_HF=tb.get_ae(Z, cfg).E, Ec_val=val, Ec_core_tot=core, Ec_cv_calcolato=core - e1s,
                      Ec=val + core, d_core_X=[d3, d4])
    V["previsioni"] = pr
    ip = {a: (pr[b]["E_HF"] - pr[a]["E_HF"] + pr[b]["Ec"] - pr[a]["Ec"]) * HA_EV for a, b in (("O", "O+"), ("N", "N+"))}
    c1 = (pr["O"]["Ec"] + 0.257902) / -0.257902
    V["C1"] = dict(valore=pr["O"]["Ec"], rel=c1, esito="PASS" if abs(c1) <= 0.05 else "FAIL", nota="non primo cieco")
    c2 = (ip["O"] - 13.618) / 13.618
    V["C2"] = dict(valore=ip["O"], rel=c2, esito="PASS" if abs(c2) <= 0.03 else "FAIL")
    V["C3"] = dict(valore=[ip["N"], ip["O"]], esito="PASS" if ip["N"] > ip["O"] else "FAIL")
    V["C4"] = dict(valore=ip["N"] - ip["O"], esito="PASS" if abs(ip["N"] - ip["O"] - 0.916) <= 0.3 else "FAIL")
    c5 = (ip["N"] - 14.534) / 14.534
    V["C5"] = dict(valore=ip["N"], rel=c5, esito="PASS" if abs(c5) <= 0.03 else "FAIL")
    V["R1"] = dict(valore=pr["N"]["Ec"], chiave=-0.1888, rel=(pr["N"]["Ec"] + 0.1888) / -0.1888)
    V["R2"] = dict(cv_O_calcolato=pr["O"]["Ec_cv_calcolato"], cv_O_stima_fase4=-0.04227)
    V["export"] = all(V[x]["esito"] == "PASS" for x in ("K5", "C1", "C2", "C3", "C4", "C5"))
    json.dump(V, open(os.path.join(OUT, "tavg_f5_verdetti" + ("_em1" if XC == (4, 5) else "") + ".json"), "w"), indent=1)
    for x in ("K5", "C1", "C2", "C3", "C4", "C5", "R1", "R2"):
        print(x, {k: v for k, v in V[x].items() if k not in ("Ec_Xa", "Ec_Xb")})
    for sp in pr:
        print(sp, {k: round(v * 1e3, 2) for k, v in pr[sp].items() if isinstance(v, float) and k != "E_HF"})
    print("export", V["export"])


if __name__ == "__main__":
    main()
