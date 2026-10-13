"""tav_analyze.py - giudizio automatico degli esami pre-registrati (lab/out/tav_prereg.json) sulle corse
(lab/out/tav_runs.json + basi in lab/out/tav_<job>.npz). Scrive lab/out/tav_results.json.
Uso: python lab/tav_analyze.py li | be"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tav_ecg as T

PRE = json.load(open(os.path.join(T.OUT, "tav_prereg.json")))
REF = PRE["riferimenti"]
RJ = os.path.join(T.OUT, "tav_runs.json")
OUTJ = os.path.join(T.OUT, "tav_results.json")
ALPHA = 7.2973525693e-3
AU_S = 2.4188843265857e-17
HC_EVNM = 1239.841984
NUC = {3: 12786.393, 4: 16424.20}


def load_svm(job):
    j = dict(PRE["svm_common"]); j.update(PRE["jobs"][job])
    z = np.load(os.path.join(T.OUT, f"tav_{job}.npz"))
    s = T.System(j["Z"], j["spin"], L=j["L"], m_kin=T.M_T_THETA_A)
    sv = T.SVM(s, np.random.default_rng(0), j["Z"])
    sv.A = z["A"].copy(); sv.u = z["u"].copy()
    sv.rebuild()
    return sv


def dipole(svS, svP, chunk=20000):
    E1, c1 = svS.state(); E2, c2 = svP.state()
    KS, KP = len(c1), len(c2)
    I, J = np.meshgrid(np.arange(KS), np.arange(KP), indexing="ij")
    I = I.ravel(); J = J.ravel()
    Dl = np.zeros(len(I)); Dv = np.zeros(len(I))
    for i0 in range(0, len(I), chunk):
        sl = slice(i0, i0 + chunk)
        Dl[sl], Dv[sl] = T.dipole_me(svS.s, svP.s, svS.A[I[sl]], svP.A[J[sl]], svP.u[J[sl]])
    zl = float(c1 @ Dl.reshape(KS, KP) @ c2); zv = float(c1 @ Dv.reshape(KS, KP) @ c2)
    return zl, zv


def element(el):
    R = json.load(open(RJ))
    res = json.load(open(OUTJ)) if os.path.exists(OUTJ) else {}
    sc = T.M0 / T.M_T_THETA_A                       # theta -> 0 (scala esatta)
    if el == "li":
        Z, g, i, p, key = 3, "li2S", "liP1S", "li2P", "Li"
    else:
        Z, g, i, p, key = 4, "be1S", "beP2S", "be1P", "Be"
    M = NUC[Z]; mu = M / (M + 1)
    rg, ri, rp = R[g], R[i], R[p]
    Eg, Ei, Ep = rg["E_mT"] * sc, ri["E_mT"] * sc, rp["E_mT"] * sc
    out = dict(E_theta0=dict(ground=Eg, ion=Ei, P=Ep), E_mT=dict(ground=rg["E_mT"], ion=ri["E_mT"], P=rp["E_mT"]),
               barre=dict(ground=rg["barra"], ion=ri["barra"], P=rp["barra"]), K=dict(ground=rg["K"], ion=ri["K"], P=rp["K"]))
    IP = (Ei - Eg) * mu * T.EH_EV
    dE = (Ep - Eg)
    line = dE * mu * T.EH_EV
    out["IP_eV"] = IP; out["line_eV"] = line; out["line_nm_vac"] = HC_EVNM / line
    # dipolo sulle basi finali (alla massa m_T; f invariante per la scala)
    svS = load_svm(g); svP = load_svm(p)
    zl, zv = dipole(svS, svP)
    dEm = svP.E - svS.E
    fL = 2 * dEm * zl ** 2 * T.M_T_THETA_A                  # f con la massa dell'operatore (regola TRK: somma = N)
    fV = 2 * zv ** 2 / (T.M_T_THETA_A * dEm)
    out.update(z_len=zl, z_vel=zv, f_L=fL, f_V=fV)
    # vita del P (massa ridotta nel dE, f di lunghezza)
    dEmu = dE * mu
    A = 2 * ALPHA ** 3 * dEmu ** 2 * fL / 3
    out["tau_P_ns"] = 1 / A * AU_S * 1e9
    ex = rg["expect"]
    lsc = T.M_T_THETA_A / T.M0                                  # lunghezze: r(m0) = r(m_T) m_T/m0
    out["r2"] = ex["r2"] * lsc ** 2; out["rho0"] = ex["delta"] / lsc ** 3; out["virial"] = ex["virial_-V/T"]
    out["virial_ion"] = ri["expect"]["virial_-V/T"]
    ok = {}
    if el == "li":
        ok["LI1"] = abs(Eg - REF["Li_E"]) <= 1.0e-4
        ok["LI2"] = abs(Ei - REF["LiP_E"]) <= 2.0e-5
        ok["LI3"] = abs(Ep - REF["Li2P_E"]) <= 1.0e-4
        ok["LI4a"] = abs(IP / REF["Li_IP_exp_eV"] - 1) <= 0.01
        ok["LI4b"] = abs(IP - REF["Li_IP_exp_eV"]) <= 0.003
        ok["LI5a"] = abs(line / REF["Li_line_eV"] - 1) <= 0.01 and abs(out["line_nm_vac"] / REF["Li_line_nm_vac"] - 1) <= 0.01
        ok["LI5b"] = abs(line - REF["Li_line_eV"]) <= 0.003
        ok["LI6"] = abs(fL / REF["Li_f"] - 1) <= 0.01 and abs(fV / REF["Li_f"] - 1) <= 0.01 and abs(fL - fV) / fL <= 0.01
        ok["LI7"] = abs(out["tau_P_ns"] / REF["Li_tau2p_ns"] - 1) <= 0.02
        ok["LI8"] = abs(out["r2"] / REF["Li_r2"] - 1) <= 0.01 and abs(out["rho0"] / REF["Li_rho0"] - 1) <= 0.01
        ok["LI9"] = abs(out["virial"] - 2) <= 1e-4
        out["diff"] = dict(E_mEh=(Eg - REF["Li_E"]) * 1e3, Eion_mEh=(Ei - REF["LiP_E"]) * 1e3, EP_mEh=(Ep - REF["Li2P_E"]) * 1e3,
                           IP_meV=(IP - REF["Li_IP_exp_eV"]) * 1e3, line_meV=(line - REF["Li_line_eV"]) * 1e3,
                           f_L_rel=fL / REF["Li_f"] - 1, f_V_rel=fV / REF["Li_f"] - 1,
                           tau_rel=out["tau_P_ns"] / REF["Li_tau2p_ns"] - 1, r2_rel=out["r2"] / REF["Li_r2"] - 1,
                           rho0_rel=out["rho0"] / REF["Li_rho0"] - 1)
        out["passa_tavola"] = bool(ok["LI1"] and ok["LI4a"] and ok["LI5a"] and ok["LI6"])
    else:
        ok["BE1"] = abs(Eg - REF["Be_E"]) <= 1.5e-3
        ok["BE2"] = abs(Ei - REF["BeP_E"]) <= 1.0e-4
        ok["BE3a"] = abs(IP / REF["Be_IP_exp_eV"] - 1) <= 0.01
        ok["BE3b"] = abs(IP - REF["Be_IP_exp_eV"]) <= 0.030
        ok["BE4a"] = abs(line / REF["Be_line_eV"] - 1) <= 0.01 and abs(out["line_nm_vac"] / REF["Be_line_nm_vac"] - 1) <= 0.01
        ok["BE4b"] = abs(line - REF["Be_line_eV"]) <= 0.030
        ok["BE5"] = abs(fL / REF["Be_f"] - 1) <= 0.03 and abs(fV / REF["Be_f"] - 1) <= 0.03 and abs(fL - fV) / fL <= 0.03
        ok["BE7"] = abs(out["virial"] - 2) <= 3e-4
        ok["X1"] = Ep > Eg + 0.15
        out["diff"] = dict(E_mEh=(Eg - REF["Be_E"]) * 1e3, Eion_mEh=(Ei - REF["BeP_E"]) * 1e3,
                           EP_mEh_vs_memoria=(Ep - REF["Be_1P_E_memoria"]) * 1e3,
                           IP_meV=(IP - REF["Be_IP_exp_eV"]) * 1e3, line_meV=(line - REF["Be_line_eV"]) * 1e3,
                           f_L_rel=fL / REF["Be_f"] - 1, f_V_rel=fV / REF["Be_f"] - 1,
                           r2_vs_memoria=out["r2"] / REF["Be_r2_memoria"] - 1, rho0_vs_memoria=out["rho0"] / REF["Be_rho0_memoria"] - 1)
        out["passa_tavola"] = bool(ok["BE1"] and ok["BE3a"] and ok["BE4a"] and ok["BE5"])
    out["esami"] = {k: bool(v) for k, v in ok.items()}
    res[key] = out
    T.save_json(OUTJ, res)
    print(json.dumps(out, indent=1, default=float))
    return out


if __name__ == "__main__":
    element(sys.argv[1])
