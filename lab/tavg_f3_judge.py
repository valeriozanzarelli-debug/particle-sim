"""tavg_f3_judge.py - FASE 3: leggi di coppia, verifica a 3 elettroni, previsione alla cieca di O e N, verdetti.
Tutto come in lab/out/tavg_prereg_f3.json (nessuna scelta fatta qui dentro dopo aver visto i bersagli).
Uso: python tavg_f3_judge.py      -> lab/out/tavg_f3_verdetti.json, lab/out/tavg_f3_figura.png"""
import os, sys, json, math
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_pecg as P
import tavg_blocks as tb
import tavg_corr as tc

OUT = os.path.join(HERE, "out")
HA_EV = 27.211386245988
TYPES = ["s2_1S", "sp_3P", "sp_1P", "pp_3P", "pp_1D"]
PRE = json.load(open(os.path.join(OUT, "tavg_prereg_f3.json")))
EX = PRE["esami"]


def J(g):
    return json.load(open(os.path.join(OUT, f"tavg_f3_{g}.json")))


def fit_law(xs, ys):
    xs, ys = np.asarray(xs), np.asarray(ys)
    out = {}
    for form, phi in (("lin", lambda x: x), ("inv", lambda x: 1 / x)):
        X = np.stack([np.ones_like(xs), phi(xs)], 1)
        c = np.linalg.lstsq(X, ys, rcond=None)[0]
        loo = []
        for k in range(len(xs)):
            m = np.arange(len(xs)) != k
            ck = np.linalg.lstsq(X[m], ys[m], rcond=None)[0]
            loo.append(float(ys[k] - X[k] @ ck))
        out[form] = dict(c=c.tolist(), loo_rms=float(np.sqrt(np.mean(np.square(loo)))),
                         res=(ys - X @ c).tolist())
    best = min(out, key=lambda f: out[f]["loo_rms"])
    return best, out


def law(L, t, x):
    f = L[t]["forma"]; a, b = L[t]["fit"][f]["c"]
    return a + b * (x if f == "lin" else 1 / x)


def mean_r(a, i):
    return float(a.g.integ(a.P[i] ** 2 * a.g.r))


def orb_index(a, n, l):
    for i, (nn, ll, occ) in enumerate(a.config):
        if nn == n and ll == l:
            return i
    return None


def valence_pairs(L, rs, rp, n_sp_same, n_sp_opp, n_pp_same, n_pp_opp):
    xs2 = 1 / rs; xsp = 1 / math.sqrt(rs * rp); xpp = 1 / rp
    e = dict(s2=law(L, "s2_1S", xs2), sp3=law(L, "sp_3P", xsp), sp1=law(L, "sp_1P", xsp),
             pp3=law(L, "pp_3P", xpp), pp1D=law(L, "pp_1D", xpp))
    e["pp1S"] = e["pp1D"]                                       # debt, pre-registered
    parts = dict(s2=e["s2"], sp_stesso=n_sp_same * e["sp3"], sp_opposto=n_sp_opp * 0.5 * (e["sp1"] + e["sp3"]),
                 pp_stesso=n_pp_same * e["pp3"], pp_opposto=n_pp_opp * (5 * e["pp1D"] + e["pp1S"] + 3 * e["pp3"]) / 9)
    return sum(parts.values()), parts, dict(x_s2=xs2, x_sp=xsp, x_pp=xpp)


COUNTS = {"O": (8, "1s2 2s2 2p4", 4, 4, 3, 3), "O+": (8, "1s2 2s2 2p3", 3, 3, 3, 0),
          "N": (7, "1s2 2s2 2p3", 3, 3, 3, 0), "N+": (7, "1s2 2s2 2p2", 2, 2, 1, 0)}
NVAL = {"O": (2, 4), "O+": (2, 3), "N": (2, 3), "N+": (2, 2)}


def main():
    V = {}
    # ---------------- M1
    m = json.load(open(os.path.join(OUT, "tavg_pecg2_macchina.json")))
    ok = all(m[k]["PASS"] for k in ("O_noee_2s2", "He_nudo"))
    V["M1"] = dict(valore={k: dict(K=m[k]["K"], dE=m[k]["dE"], monotona=m[k]["monotona"]) for k in ("O_noee_2s2", "He_nudo")},
                   esito="PASS" if ok else "FAIL")
    # ---------------- pair laws
    g2 = J("g2")
    L = {}
    for t in TYPES:
        rows = [g2[f"P_Z{Z}_{t}"] for Z in (5, 6, 7, 8)]
        xs = [1 / math.sqrt(r["r_medio"][0] * r["r_medio"][1]) for r in rows]
        ys = [r["Ec"] for r in rows]
        best, fits = fit_law(xs, ys)
        L[t] = dict(forma=best, fit=fits, x=xs, Ec=ys, barra=[r["barra"] for r in rows])
    V["leggi"] = L
    V["pp_1S_riportato"] = {Z: g2[f"P_Z{Z}_pp_1S"]["Ec"] for Z in (5, 6, 7, 8)}
    # ---------------- V1
    g1 = J("g1")
    d = {st: g1[f"V_O_{st}"]["E"] - g2[f"P_Z8_{st}"]["E_inf"] for st in ("s2_1S", "sp_3P", "sp_1P")}
    V["V1"] = dict(valore=d, esito="PASS" if max(abs(v) for v in d.values()) <= EX["V1"]["soglia_Eh"] else "FAIL")
    # ---------------- A1
    g3 = J("g3"); a1 = {}
    for Z in (6, 8):
        f = P.BlockField(Z)
        a = f.radial_atom([(1, 0, 2), (2, 1, 1)]); a.scf()
        rs, rp = mean_r(a, 0), mean_r(a, 1)
        ec_ex = g3[f"A3_Z{Z}_2s2_2p"]["E"] - a.E
        ec_pr, parts, xs = valence_pairs(L, rs, rp, 1, 1, 0, 0)
        a1[Z] = dict(E_ECG=g3[f"A3_Z{Z}_2s2_2p"]["E"], K=g3[f"A3_Z{Z}_2s2_2p"]["K"], E_ROHF=a.E, Ec_esatta=ec_ex,
                     Ec_coppie=ec_pr, rel=(ec_pr - ec_ex) / ec_ex, parti=parts, **xs)
    V["A1"] = dict(valore=a1, esito="PASS" if all(abs(v["rel"]) <= EX["A1"]["soglia_rel"] for v in a1.values()) else "FAIL")
    # ---------------- core and core-valence
    core, cv = {}, {}
    for Z in (7, 8):
        Ehe, _ = tc._ecg(f"he1S_Z{Z}")
        core[Z] = Ehe - tb.get_ae(Z, "1s2").E
        cv[Z] = {}
        for lab, cfg in (("2S", "1s2 2s1"), ("2P", "1s2 2p1")):
            cv[Z][lab] = (g2[f"Li_Z{Z}_{lab}"]["E"] - tb.get_ae(Z, cfg).E) - core[Z]
    V["core_1s2"] = core; V["core_valenza_per_elettrone"] = cv
    # ---------------- blind
    pr = {}
    for name, (Z, cfg, nss, nso, nps, npo) in COUNTS.items():
        a = tb.get_ae(Z, cfg)
        rs = mean_r(a, orb_index(a, 2, 0)); rp = mean_r(a, orb_index(a, 2, 1))
        ev, parts, xs = valence_pairs(L, rs, rp, nss, nso, nps, npo)
        ns, npp = NVAL[name]
        ecv = ns * cv[Z]["2S"] + npp * cv[Z]["2P"]
        ec = core[Z] + ecv + ev
        pr[name] = dict(E_HF=a.E, r_2s=rs, r_2p=rp, Ec=ec, Ec_core=core[Z], Ec_cv=ecv, Ec_val=ev, parti=parts, **xs)
    V["previsioni"] = pr
    ipO = (pr["O+"]["E_HF"] - pr["O"]["E_HF"] + pr["O+"]["Ec"] - pr["O"]["Ec"]) * HA_EV
    ipN = (pr["N+"]["E_HF"] - pr["N"]["E_HF"] + pr["N+"]["Ec"] - pr["N"]["Ec"]) * HA_EV
    c1 = (pr["O"]["Ec"] - EX["C1"]["chiave_Eh"]) / EX["C1"]["chiave_Eh"]
    V["C1"] = dict(valore=pr["O"]["Ec"], rel=c1, esito="PASS" if abs(c1) <= EX["C1"]["soglia_rel"] else "FAIL")
    c2 = (ipO - EX["C2"]["chiave_eV"]) / EX["C2"]["chiave_eV"]
    V["C2"] = dict(valore=ipO, rel=c2, esito="PASS" if abs(c2) <= EX["C2"]["soglia_rel"] else "FAIL")
    V["C3"] = dict(valore=[ipN, ipO], esito="PASS" if ipN > ipO else "FAIL")
    V["C4"] = dict(valore=ipN - ipO, esito="PASS" if abs(ipN - ipO - EX["C4"]["chiave_eV"]) <= EX["C4"]["soglia_eV"] else "FAIL")
    c5 = (ipN - EX["C5"]["chiave_eV"]) / EX["C5"]["chiave_eV"]
    V["C5"] = dict(valore=ipN, rel=c5, esito="PASS" if abs(c5) <= EX["C5"]["soglia_rel"] else "FAIL")
    V["R1"] = dict(valore=pr["N"]["Ec"], chiave=-0.1888, nota="riportato")
    V["R2"] = dict(valore=(pr["O+"]["Ec"] - pr["O"]["Ec"]) * 1e3, chiave_mEh=63.7, nota="riportato")
    json.dump(V, open(os.path.join(OUT, "tavg_f3_verdetti.json"), "w"), indent=1)
    for k in ("M1", "V1", "A1", "C1", "C2", "C3", "C4", "C5"):
        print(k, V[k]["esito"], json.dumps(V[k].get("rel", V[k]["valore"]))[:200])
    print("O:", {k: (round(v * 1e3, 2) if isinstance(v, float) else v) for k, v in pr["O"].items() if k.startswith("Ec")})
    figure(V)


def figure(V):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    L = V["leggi"]; pr = V["previsioni"]
    for t, c in zip(TYPES, ["C0", "C1", "C2", "C3", "C4"]):
        xs = np.array(L[t]["x"]); ax[0].plot(xs, np.array(L[t]["Ec"]) * 1e3, "o", color=c, label=f"{t} ({L[t]['forma']})")
        xx = np.linspace(0.6, max(xs) * 1.05, 100)
        ax[0].plot(xx, [law(L, t, x) * 1e3 for x in xx], "-", color=c, lw=1)
    for k, mk in (("x_s2", "s"), ("x_sp", "^"), ("x_pp", "v")):
        ax[0].axvline(pr["O"][k], color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("x = 1/sqrt(<r>a <r>b)  [1/a0]   (punteggiate: O neutro)"); ax[0].set_ylabel("E_c coppia [mEh]")
    ax[0].set_title("a) coppie a 2 e- nei blocchi B, C, N, O e leggi"); ax[0].legend(fontsize=7)
    names = ["O", "O+", "N", "N+"]
    b0 = np.array([pr[n]["Ec_core"] for n in names]) * 1e3
    b1 = np.array([pr[n]["Ec_cv"] for n in names]) * 1e3
    b2 = np.array([pr[n]["Ec_val"] for n in names]) * 1e3
    ax[1].bar(names, b0, label="1s2"); ax[1].bar(names, b1, bottom=b0, label="core-valenza")
    ax[1].bar(names, b2, bottom=b0 + b1, label="valenza (coppie)")
    ax[1].plot(["O"], [-257.9], "k*", ms=12, label="O chiave -257.9"); ax[1].plot(["N"], [-188.8], "kx", ms=10, label="N riportato -188.8")
    ax[1].set_ylabel("E_c [mEh]"); ax[1].legend(fontsize=7)
    ax[1].set_title(f"b) cieco: E_c(O) {pr['O']['Ec'] * 1e3:.1f} mEh, IP(O) {V['C2']['valore']:.3f} eV")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "tavg_f3_figura.png"), dpi=110)


if __name__ == "__main__":
    main()
