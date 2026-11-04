# secondo tentativo dichiarato: scansione 1D in R a CBS(2,3), angoli al minimo X=2
import json, os, numpy as np
import acqua_dimero_plan as P
TD, TA = 4.8, -55.3; RS = [2.88, 2.92, 2.97]
def cbs(R):
    g = (R, TD, TA); hf3 = P.eint(3, g, "hf"); c2, c3 = P.eint(2, g, "corr"), P.eint(3, g, "corr")
    return hf3 + (27 * c3 - 8 * c2) / 19
E = np.array([cbs(R) for R in RS]); a, b, c = np.polyfit(RS, E, 2)
Rm = -b / (2 * a); Em = np.polyval([a, b, c], Rm); De = -Em * P.EH_EV
out = dict(R=RS, Eint_cbs=E.tolist(), R_min=Rm, D_e_eV=De, verdetti=dict(
    D2b_ROO=dict(valore=Rm, rif=2.91, rel=(Rm - 2.91) / 2.91, esito="PASS" if abs(Rm - 2.91) / 2.91 <= 0.02 else "FAIL"),
    D1b_De=dict(valore=De, rif=0.217, rel=(De - 0.217) / 0.217, esito="PASS" if abs(De - 0.217) / 0.217 <= 0.10 else "FAIL")))
json.dump(out, open("out/acqua_dimero2_verdetti.json", "w"), indent=1); print(json.dumps(out))
