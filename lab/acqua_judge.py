"""acqua_judge.py - verdetti H2O contro lab/out/acqua_prereg.json -> lab/out/acqua_verdetti.json"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_plan as PL, acqua_run as AR, acqua_basis as AB

A2B = AR.A2B; EH_EV = 27.211386245988; AU_D = 2.541746473; CM = 219474.6313632
ME_U = 1822.888486209; MH = 1836.15267343; MO = 15.99491461957 * ME_U - 8


def ec(p):
    return p["Ec_ccsd"] + p["Ec_T"]


def Ecbs(p2, p3):
    return p3["E_hf"] + (27 * ec(p3) - 8 * ec(p2)) / 19


def internal(x):
    O, H1, H2 = x[:3], x[3:6], x[6:]
    a, b = H1 - O, H2 - O
    return np.array([np.linalg.norm(a), np.linalg.norm(b), math.acos(a @ b / np.linalg.norm(a) / np.linalg.norm(b))])


def freqs(r, th, F):
    x = np.array(AR.geom(r, r, th)).ravel(); Bm = np.zeros((3, 9)); h = 1e-5
    for k in range(9):
        e = np.zeros(9); e[k] = h; Bm[:, k] = (internal(x + e) - internal(x - e)) / (2 * h)
    Hc = Bm.T @ F @ Bm; m = np.repeat([MO, MH, MH], 3); Hm = Hc / np.sqrt(np.outer(m, m))
    w = np.sort(np.linalg.eigvalsh(Hm))[-3:]
    return list(np.sqrt(w) * CM)


def main():
    s = PL.st_load(); P = json.load(open(os.path.join(AB.OUT, "acqua_prereg.json"))); ex = P["esami"]
    c = s["c3"]; S = PL.surface(c, "cbs"); S3 = PL.surface(c, 3); S2 = PL.surface(s["c2"], 2)
    r, th = s["min_cbs"]
    out = dict(centro_griglia=c, minimo_cbs=S["min"], minimo_X3=S3["min"], minimo_X2=S2["min"], hf2_min=s["hf2_min"])
    # force constants (bohr, rad)
    H = np.array(S["hess"]); cr = 1 / A2B; ct = 180 / math.pi
    Ssym = H[0, 0] * cr * cr; frt = 0.5 * H[0, 1] * cr * ct; ftt = H[1, 1] * ct * ct
    E0 = Ecbs(PL.pt(2, r, r, th), PL.pt(3, r, r, th))
    dE = [Ecbs(PL.pt(2, r + d, r - d, th), PL.pt(3, r + d, r - d, th)) - E0 for d in PL.DELTAS]
    d = PL.DELTAS[0] * A2B
    Aa = (16 * dE[0] - dE[1]) / (12 * d * d)
    frr = 0.5 * (Ssym / 2 + Aa); frr2 = 0.5 * (Ssym / 2 - Aa)
    F = np.array([[frr, frr2, frt], [frr2, frr, frt], [frt, frt, ftt]])
    w = freqs(S["min"][0], S["min"][1], F)
    out["costanti_forza_au"] = dict(f_rr=frr, f_rr1=frr2, f_rt=frt, f_tt=ftt)
    out["frequenze_cm-1"] = dict(piegamento=w[0], stiramento_sim=w[1], stiramento_asim=w[2])

    def dip(getE):
        h = 0.002
        E = {F_: getE(F_) for F_ in PL.FIELDS}
        dEdF = (8 * (E[h] - E[-h]) - (E[2 * h] - E[-2 * h])) / (12 * h)
        return abs(dEdF) * AU_D
    mu3 = dip(lambda F_: PL.pt(3, r, r, th, F_)["E_tot"])
    mucbs = dip(lambda F_: Ecbs(PL.pt(2, r, r, th, F_), PL.pt(3, r, r, th, F_)))
    out["dipolo_D"] = dict(X3=mu3, cbs=mucbs, hf_X3=abs(PL.pt(3, r, r, th)["mu_hf_z"]) * AU_D)
    at = {(el, X): AR.run(["atom", el, str(X)]) for el in "OH" for X in (2, 3)}
    EO = Ecbs(at[("O", 2)], at[("O", 3)]); EH = at[("H", 3)]["E_hf"]
    De = EO + 2 * EH - E0
    De3 = at[("O", 3)]["E_tot"] + 2 * EH - PL.pt(3, r, r, th)["E_tot"]
    out["D_e"] = dict(cbs_Eh=De, cbs_eV=De * EH_EV, X3_Eh=De3, HF_X3_Eh=at[("O", 3)]["E_hf"] + 2 * EH - PL.pt(3, r, r, th)["E_hf"])
    v = {}
    rr, tt = S["min"]
    v["A1_r_OH"] = dict(valore=rr, rif=ex["A1_r_OH"]["rif_A"], rel=(rr - 0.9572) / 0.9572, esito="PASS" if abs(rr - 0.9572) / 0.9572 <= 0.005 else "FAIL")
    v["A2_angolo"] = dict(valore=tt, rif=104.52, diff=tt - 104.52, esito="PASS" if abs(tt - 104.52) <= 1.0 else "FAIL")
    v["A3_D_e"] = dict(valore=De, rif=0.3705, rel=(De - 0.3705) / 0.3705, esito="PASS" if abs(De - 0.3705) / 0.3705 <= 0.03 else "FAIL")
    v["A4_dipolo"] = dict(valore=mu3, rif=1.855, rel=(mu3 - 1.855) / 1.855, esito="PASS" if abs(mu3 - 1.855) / 1.855 <= 0.05 else "FAIL")
    ref = [1649, 3832, 3943]; rel = [(a - b) / b for a, b in zip(w, ref)]
    v["A5_frequenze"] = dict(valore=w, rif=ref, rel=rel, esito="PASS" if max(abs(x) for x in rel) <= 0.03 else "FAIL")
    out["verdetti"] = v
    json.dump(out, open(os.path.join(AB.OUT, "acqua_verdetti.json"), "w"), indent=1)
    print(json.dumps({k: (round(x["valore"], 4) if not isinstance(x["valore"], list) else [round(y) for y in x["valore"]], x["esito"]) for k, x in v.items()}))


if __name__ == "__main__":
    main()
