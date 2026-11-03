"""acqua_checks.py - controlli di macchina del motore molecolare (studio ACQUA).
  C-INT: integrali contro quadratura a forza bruta (acqua_quad.py, soglia 1e-8)      -> letto da out/acqua_quad.out
  C-CC : CCSD e (T) a guscio chiuso contro il codice a spin-orbitali di tavg_cc (stessi integrali, soglia 1e-9)
  C-H  : atomo H, HF nella base X=4 aug contro -0.5 esatto (soglia 0.2 mEh)
  C-H2 : H2, CCSD (= full CI per 2 elettroni), CBS(3,4): R_e 1.4009 entro 0.003 bohr, D_e 4.7477 eV entro 0.02 eV
  C-O  : atomo O 3P, CCSD(T) a tutti gli elettroni CBS(3,4) aug: E_c = E_tot - E_HF(livello fine, -74.809398)
         entro il 2 % di E_c tavg fase 5 (-0.256225 Eh); riportata anche E_c rispetto al nostro ROHF
Soglie fissate in questo file PRIMA di eseguire i controlli.
Uso: python acqua_checks.py <pezzo>   pezzo in {cc, h2pts, opts, judge}"""
import os, sys, json, math, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_run as AR, acqua_basis as AB

EH_EV = 27.211386245988
RH2 = [1.30, 1.35, 1.40, 1.45, 1.50]
F = os.path.join(AB.OUT, "acqua_checks.json")


def cbs(e2, e3, X2, X3):
    return (X3 ** 3 * e3 - X2 ** 3 * e2) / (X3 ** 3 - X2 ** 3)


def check_cc():
    import acqua_int as I, acqua_qc as Q
    r = 1.81; th = math.radians(104.5)
    cen = [[0, 0, 0], [r * math.sin(th / 2), 0, r * math.cos(th / 2)], [-r * math.sin(th / 2), 0, r * math.cos(th / 2)]]
    cen, Z, sh = AB.shells_for(["O", "H", "H"], cen, 2, aug=False)
    B = I.Basis(cen, Z, sh); ints = I.integrals(B); res = Q.scf(B, ints, 5)
    h, G = Q.mo_ints(B, ints, res); G2 = G.copy()
    E, t1, t2, c = Q.rccsd(G, res["eps"], 5); ET = Q.rccsd_t(G, res["eps"], 5, t1, t2)
    o = Q.so_ccsd_t(h, G2, 5, 0)
    return dict(M=B.M, rccsd=E, so_ccsd=o["Ec_ccsd"], rT=ET, so_T=o["Ec_T"],
                diff=max(abs(E - o["Ec_ccsd"]), abs(ET - o["Ec_T"])), esito="PASS" if max(abs(E - o["Ec_ccsd"]), abs(ET - o["Ec_T"])) < 1e-9 else "FAIL")


def judge():
    out = json.load(open(F)) if os.path.exists(F) else {}
    q = open(os.path.join(AB.OUT, "acqua_quad.out")).read().strip().splitlines()[-1]
    w = eval(q); out["C-INT"] = dict(errori=w, esito="PASS" if max(w.values()) < 1e-8 else "FAIL")
    H = {X: AR.run(["atom", "H", str(X)]) for X in (2, 3, 4)}
    out["C-H"] = dict(E_hf={X: H[X]["E_hf"] for X in H}, esito="PASS" if abs(H[4]["E_hf"] + 0.5) < 2e-4 else "FAIL")
    res = {}
    for X in (2, 3, 4):
        res[X] = [AR.run(["h2", str(X), f"{R:.6f}"]) for R in RH2]
    curves = {}
    for X in (2, 3, 4):
        curves[X] = np.array([p["E_tot"] for p in res[X]])
    hf4 = np.array([p["E_hf"] for p in res[4]])
    ec = {X: np.array([p["E_tot"] - p["E_hf"] for p in res[X]]) for X in (3, 4)}
    curves["cbs"] = hf4 + cbs(ec[3], ec[4], 3, 4)
    Rg = np.linspace(1.3, 1.5, 2001); tab = {}
    for k, E in curves.items():
        c = np.polyfit(RH2, E, 4); Ev = np.polyval(c, Rg); i = int(np.argmin(Ev))
        eh = H[4]["E_hf"] if k == "cbs" else H[k]["E_hf"]
        tab[str(k)] = dict(R_e=float(Rg[i]), D_e_eV=float((2 * eh - Ev[i]) * EH_EV))
    t = tab["cbs"]
    out["C-H2"] = dict(tabella=tab, esito="PASS" if abs(t["R_e"] - 1.4009) <= 0.003 and abs(t["D_e_eV"] - 4.7477) <= 0.02 else "FAIL")
    O = {X: AR.run(["atom", "O", str(X)]) for X in (2, 3, 4)}
    ecO = {X: O[X]["E_tot"] - O[X]["E_hf"] for X in O}
    Ecb = cbs(ecO[3], ecO[4], 3, 4); Etot = O[4]["E_hf"] + Ecb
    Ec_fine = Etot - (-74.80939846910752); ref = -0.25622540436465935
    out["C-O"] = dict(E_hf={X: O[X]["E_hf"] for X in O}, Ec_proprio={X: ecO[X] for X in O}, Ec_cbs_proprio=Ecb,
                      E_tot_cbs=Etot, Ec_vs_HF_fine=Ec_fine, rif=ref, rel=(Ec_fine - ref) / ref,
                      esito="PASS" if abs((Ec_fine - ref) / ref) <= 0.02 else "FAIL")
    json.dump(out, open(F, "w"), indent=1)
    print(json.dumps({k: v.get("esito") for k, v in out.items()}))


if __name__ == "__main__":
    p = sys.argv[1]
    if p == "cc":
        out = json.load(open(F)) if os.path.exists(F) else {}
        out["C-CC"] = check_cc(); json.dump(out, open(F, "w"), indent=1); print(out["C-CC"])
    elif p == "judge":
        judge()
