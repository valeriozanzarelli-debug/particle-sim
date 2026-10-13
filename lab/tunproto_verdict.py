"""tunproto_verdict.py - giudizio delle corse contro lab/out/tunproto_prereg.json -> lab/out/tunproto_verdict.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunproto_core import OUT

PRE = json.load(open(os.path.join(OUT, "tunproto_prereg.json")))
ld = lambda n: json.load(open(os.path.join(OUT, n)))
V = {}


def put(k, ok, **info):
    V[k] = dict(esito="PASS" if ok else "FAIL", **info)
    print(k, V[k]["esito"], {kk: vv for kk, vv in info.items() if not isinstance(vv, (list, dict))})


# ---------------------------------------------------------------- A
RA = ld("tunproto_rA_regola.json")
A = PRE["A_misura"]
d1 = 0
for k, v in RA.items():
    for c in ("m", "V"):
        for q in ("tz", "ty"):
            d1 = max(d1, abs(v["run"][c][q] - v["spett"][c][q]))
put("A1", d1 <= 1e-5, max_scarto=d1)
d2, d2q, rows = 0, 0, {}
for V0 in A["V0"]:
    pr = A["previsioni_dirac_spettrale"][str(V0)]
    for c in ("m", "V"):
        for q in ("tz", "ty"):
            t1 = RA["a0.025_V%.1f" % V0]["run"][c][q]
            t2 = RA["a0.0125_V%.1f" % V0]["run"][c][q]
            t0 = RA["a0.05_V%.1f" % V0]["run"][c][q]
            lin = 2 * t2 - t1
            quad = (4 * t2 - t1) / 3
            d2 = max(d2, abs(lin - pr[c][q]))
            d2q = max(d2q, abs(quad - pr[c][q]))
            rows["V%.1f_%s_%s" % (V0, c, q)] = dict(a005=t0, a0025=t1, a00125=t2, lin=lin, quad=quad, dirac=pr[c][q])
put("A2", d2 <= 1e-3, max_scarto_lineare=d2, max_scarto_quadratico_posthoc=d2q, righe=rows)
RD = ld("tunproto_rA_dirac.json")
d3, d3p = 0, 0
for V0 in A["V0"]:
    pr = A["previsioni_dirac_spettrale"][str(V0)]["m"]
    r = RD["V%.1f" % V0]
    d3 = max(d3, abs(r["tz"] - pr["tz"]), abs(r["ty"] - pr["ty"]))
    d3p = max(d3p, abs(r["PT"] / pr["PT"] - 1))
put("A3", d3 <= 5e-3 and d3p <= 1e-3, max_scarto_tau=d3, max_scarto_PT_rel=d3p)
ok4, info4 = True, {}
for V0 in A["V0"]:
    st = A["stazionario_dirac_E0"][str(V0)]
    W = A["limite_debole_stazionario"][str(V0)]
    for c, (ky, kz) in (("m", ("tau_my", "tau_mz")), ("V", ("tau_y", "tau_z"))):
        dev = [max(abs(W[s][c]["ty"] - st[ky]), abs(W[s][c]["tz"] - st[kz])) for s in ("10", "20", "40", "80")]
        info4["V%.1f_%s" % (V0, c)] = dev
        ok4 &= dev[-1] <= 1e-2 and all(np.diff(dev) <= 1e-12)
put("A4", ok4, scarti_sigma_10_20_40_80=info4)
r5 = {V0: RA["a0.0125_V%.1f" % V0]["run"] for V0 in A["V0"]}
L = A["L"]
c5 = dict(segno_V24=(r5[2.4]["m"]["ty"] < 0 and r5[2.4]["V"]["ty"] > 0),
          tmz_su_L=r5[2.0]["m"]["tz"] / L, tz_su_L=r5[2.0]["V"]["tz"] / L,
          rapporto_16=r5[1.6]["m"]["ty"] / r5[1.6]["V"]["ty"], rapporto_24=r5[2.4]["m"]["ty"] / r5[2.4]["V"]["ty"])
ok5 = (c5["segno_V24"] and 0.97 <= c5["tmz_su_L"] <= 1.03 and c5["tz_su_L"] < 0.15
       and abs(c5["rapporto_16"] - 0.4) <= 0.05 and abs(c5["rapporto_24"] + 0.4) <= 0.05)
put("A5", ok5, **c5)

# ---------------------------------------------------------------- B
RB = ld("tunproto_rA_bordi.json")
B = PRE["B_sistematici"]
d = 0
for k, r in RB.items():
    if k.startswith("w"):
        sp = B["bordi_lisci"]["val"][k]["spettrale_L4"]
        for c in ("m", "V"):
            for q in ("tz", "ty"):
                d = max(d, abs(r[c][q] - sp[c][q]))
put("B1", d <= 1e-4, max_scarto=d)
sl = {k: v["slope_m"] for k, v in B["bordi_lisci"]["val"].items()}
okr = all(abs(sl["w%.2f_rettangolo" % w] - 1) <= 0.05 for w in (0.0, 0.25, 0.5))
put("B1r", okr, pendenze=sl)
bd = B["delta_finito"]["val"]
base = bd["0.001"]["tz"]
okb, inf = True, {}
for k, v in bd.items():
    rel = abs(v["tz"] / base - 1)
    inf[k] = dict(delta_tau=float(k) * base, scarto_rel=rel)
    if float(k) * base <= 0.8:
        okb &= rel <= 0.01
dd = max(abs(RB["delta%.1f" % x]["run"]["m"]["tz"] - RB["delta%.1f" % x]["spett"]["tz"]) for x in (0.2, 0.8))
put("B2", okb and dd <= 1e-4, scarti=inf, regola_vs_spettrale=dd)

# ---------------------------------------------------------------- C
RC = ld("tunproto_rC_pulito.json")
C = PRE["C_fibra"]
th = C["theta"]
dc, worst = 0, ""
for V0 in C["V0"]:
    for N in C["N_siti"]:
        sp = C["val"]["V%.1f_N%d" % (V0, N)]["spett"]
        r = RC["d0.15_V%.1f_N%d" % (V0, N)]
        for c in ("m", "V"):
            for q in ("tz", "ty"):
                tol = max(0.03 * abs(sp[c][q]), 0.05)
                ex = abs(r[c][q] - sp[c][q]) / tol
                if ex > dc:
                    dc, worst = ex, "V%.1f_N%d_%s_%s" % (V0, N, c, q)
put("C1", dc <= 1, max_scarto_su_tolleranza=dc, peggiore=worst)
Ls = np.array([4, 6, 8]) * th
slm = np.polyfit(Ls, [RC["d0.15_V2.0_N%d" % N]["m"]["tz"] for N in (4, 6, 8)], 1)[0]
slv = np.polyfit(Ls, [RC["d0.15_V2.0_N%d" % N]["V"]["tz"] for N in (4, 6, 8)], 1)[0]
put("C2", abs(slm * np.cos(th) - 1) <= 0.10 and slv < 0.3, pendenza_m=slm, pendenza_V=slv, attesa=1 / np.cos(th))
MC = ld("tunproto_rC_mc.json")
put("C3", MC["delta0.15"]["z_D1_mediana"] >= 4.5, **MC["delta0.15"])
put("C4", MC["delta0.15"]["z_D2_mediana"] >= 4.5, z=MC["delta0.15"]["z_D2_mediana"])
put("C3b", MC["delta0.60"]["z_D1_mediana"] >= 4.5, z=MC["delta0.60"]["z_D1_mediana"])
put("C4b", MC["delta0.60"]["z_D2_mediana"] >= 4.5, z=MC["delta0.60"]["z_D2_mediana"])
RN = ld("tunproto_rC_rumore.json")
for name, rr in RN.items():
    ok = (abs(rr["pendenza_m"] / slm - 1) <= 0.10 and rr["V2.4_N8"]["m"]["ty"] < 0 and rr["V2.4_N8"]["m"]["r"] >= 0.5)
    put("C5_" + name, ok, pendenza_m=rr["pendenza_m"], pendenza_V=rr["pendenza_V"], tmy_V24=rr["V2.4_N8"]["m"]["ty"],
        r=rr["V2.4_N8"]["m"]["r"], rel=rr["pendenza_m"] / slm - 1)

# ---------------------------------------------------------------- D
RW = ld("tunproto_rD.json")
Dp = PRE["D_guide"]
dmax = 0
for k, r in RW["pulito"].items():
    p = Dp["previsioni_dirac"][k]["m"]["tz"]
    dmax = max(dmax, abs(r["m"]["tz_m"] / p - 1))
Lc = np.array(Dp["L_guide"]) / 6.0
slw = np.polyfit(Lc, [RW["pulito"]["V3.0_L%d" % Lg]["m"]["tz_m"] for Lg in Dp["L_guide"]], 1)[0]
put("D1", dmax <= 0.10 and abs(slw - 1) <= 0.15, max_scarto_rel=dmax, pendenza=slw)
co = RW["contaminazione"]
put("D2", abs(co["tz_m_on"] / co["tz_m_pulito"] - 1) <= 0.10 and abs(co["tz_m_corretto"] / co["tz_m_pulito"] - 1) <= 0.01,
    spostamento=co["tz_m_on"] / co["tz_m_pulito"] - 1, residuo=co["tz_m_corretto"] / co["tz_m_pulito"] - 1)
di = RW["disordine"]
put("D3", abs(di["media"] / di["pulito"] - 1) <= 0.05 and di["dispersione"] / di["pulito"] <= 0.10,
    media_rel=di["media"] / di["pulito"] - 1, dispersione_rel=di["dispersione"] / di["pulito"])

n_pass = sum(v["esito"] == "PASS" for v in V.values())
print("PASS %d / %d" % (n_pass, len(V)))
json.dump(dict(esami=V, pass_=n_pass, totale=len(V)), open(os.path.join(OUT, "tunproto_verdict.json"), "w"), indent=1, default=float)
