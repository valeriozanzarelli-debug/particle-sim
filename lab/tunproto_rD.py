"""tunproto_rD.py - corse D: il protocollo sulle guide d'onda binarie (modi accoppiati), orologio = polarizzazione.
Uso: python3 tunproto_rD.py      (checkpoint: un JSON)."""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunproto_core import OUT, wg_run

PRE = json.load(open(os.path.join(OUT, "tunproto_prereg.json")))
D = PRE["D_guide"]


def go(V0, Lg, **kw):
    kw.setdefault("delta", D["delta_cm"])
    return wg_run(D["kappa_cm"], D["sigma_massa_cm"], V0, Lg, D["sx_guide"], D["E0_cm"], Z=D["Z_cm"],
                  margin=D["margin"], **kw)


def main():
    t0 = time.time()
    res = {"pulito": {}}
    for V0 in D["V0_cm"]:
        for Lg in D["L_guide"]:
            r = {c: go(V0, Lg, clock=c) for c in ("m", "V")}
            res["pulito"]["V%.1f_L%d" % (V0, Lg)] = r
            p = D["previsioni_dirac"]["V%.1f_L%d" % (V0, Lg)]
            print(V0, Lg, "wg tz_m %.3f ty_m %.3f | dirac %.3f %.3f | V wg %.3f dirac %.3f" % (
                r["m"]["tz_m"], r["m"]["ty_m"], p["m"]["tz"], p["m"]["ty"], r["V"]["tz_m"], p["V"]["tz"]), flush=True)
    # contaminazione: accoppiamento dipendente dalla polarizzazione (eta), con e senza orologio
    eta = D["eta_kappa"]
    on = go(3.0, 18, clock="m", eta_kappa=eta)
    off = go(3.0, 18, clock="m", eta_kappa=eta, delta=0.0)
    clean = res["pulito"]["V3.0_L18"]["m"]
    sig = D["sigma_massa_cm"]
    tz_corr = -(np.arctanh(on["sz"]) - np.arctanh(off["sz"])) / D["delta_cm"] * sig
    res["contaminazione"] = dict(eta=eta, tz_m_on=on["tz_m"], sz_off=off["sz"], tz_m_corretto=float(tz_corr),
                                 tz_m_pulito=clean["tz_m"])
    print("contaminazione", res["contaminazione"], flush=True)
    # disordine di fabbricazione: un campione = una realizzazione (lettura per campione)
    vals = []
    for s in range(D["realizzazioni_disordine"]):
        r = go(3.0, 18, clock="m", disorder=D["disordine_cm"], seed=100 + s)
        vals.append(r["tz_m"])
    vals = np.array(vals)
    res["disordine"] = dict(std_cm=D["disordine_cm"], media=float(vals.mean()), dispersione=float(vals.std()),
                            valori=vals.tolist(), pulito=clean["tz_m"])
    print("disordine", vals.mean(), vals.std(), clean["tz_m"], "%.0fs" % (time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(OUT, "tunproto_rD.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
