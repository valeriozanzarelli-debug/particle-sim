"""Idrogeno di Bohr, fase 4: esami C (cucitura), V (convergenza nel blocco), R (energia N=256), M (massa del maestro)
contro lab/out/idrogeno_bohr_fase4_prereg.json. Legge i log delle corse Rust (cartella data come argomento: C_pav.log,
C_ibr.log, V_<ni>.log, R.log, m3_2p_000.log, m3_2p_111.log) e i json della regola; G e' in idrogeno_bohr_fase4_grosso.py.
Scrive lab/out/idrogeno_bohr_fase4_esame.json.   Uso: python3 lab/idrogeno_bohr_fase4_esame.py <cartella log>"""
import os, re, sys, json
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
D = sys.argv[1]
out = {"esami": {}, "numeri": {}}
num, ex = out["numeri"], out["esami"]


def rd(fn):
    p = os.path.join(D, fn)
    return open(p).read() if os.path.exists(p) else ""


def chk(nome, ok, det):
    ex[nome] = dict(passed=bool(ok), dettaglio=det)
    print(f"{'PASS' if ok else 'FAIL'}  {nome}: {det}")


def risultati(t):
    r = {}
    m = re.search(r"hybrid/coarse-1 (\S+)", t); r["ibrido_su_grosso"] = float(m.group(1)) if m else None
    m = re.search(r"RISULTATO ibrido .*?E_hybrid (\S+) eV \((\S+) %\)", t); r["E_eV"], r["E_rel_pct"] = (float(m.group(1)), float(m.group(2))) if m else (None, None)
    m = re.search(r"coarse density max\|drho\|/max rho (\S+), fine block (\S+)", t); r["staz_grosso"], r["staz_fine"] = (float(m.group(1)), float(m.group(2))) if m else (None, None)
    m = re.search(r"non coperte (\S+)", t); r["staz_fuori"] = float(m.group(1)) if m else None
    js = re.findall(r"DIAG (\d+) .*j_coperte (\S+) j_fuori (\S+)", t)
    r["j"] = [(int(a), float(b), float(c)) for a, b, c in js]
    pp = [(int(a), float(b), float(c)) for a, b, c in re.findall(r"DIAG (\d+) perdita (\S+) sfasamento (\S+)", t)]
    r["perdita_mediana_2a_meta"] = float(np.median([abs(p) for k, p, s in pp if k > max(k for k, _, _ in pp) / 2])) if pp else None
    r["sfasamento_medio_2a_meta"] = float(np.mean([s for k, p, s in pp if k > max(k for k, _, _ in pp) / 2])) if pp else None
    m = re.search(r"finestra \[(\d+),(\d+)\]: E_hybrid (\S+) eV \((\S+) %\)\n(?!RISULTATO finestra)", t)
    fin = re.findall(r"finestra \[(\d+),(\d+)\]: E_hybrid (\S+) eV \((\S+) %\)", t)
    r["finestre"] = [(int(a), int(b), float(c), float(d)) for a, b, c, d in fin]
    r["partenza"] = re.findall(r"partenza veff iter (\d+): V_eff nucleo (\S+) eV.*?E\(G\+V\+V_eff\) (\S+) eV \((\S+) %\), 1-\|<prima\|dopo>\|\^2 = (\S+)", t)
    return r


# ------------------------------------------------------------------ C: cucitura
pav, ibr = risultati(rd("C_pav.log")), risultati(rd("C_ibr.log"))
if pav["j"] and ibr["j"]:
    num["C_pavimento"], num["C_ibrido"] = pav, ibr
    chk("C1_perdita", ibr["perdita_mediana_2a_meta"] <= 2e-7, f"mediana |perdita| passi 121..240 = {ibr['perdita_mediana_2a_meta']:.2e} per passo (soglia 2e-7; fase 3 ~5e-6)")
    chk("C2_stazionarieta_fuori", ibr["staz_fuori"] <= 10 * pav["staz_fuori"],
        f"{ibr['staz_fuori']:.3e} contro 10 x pavimento = {10 * pav['staz_fuori']:.3e} ({ibr['staz_fuori'] / pav['staz_fuori']:.1f} x)")
    chk("C2b_stazionarieta_fase3", ibr["staz_grosso"] <= 1e-4, f"stazionarieta grosso {ibr['staz_grosso']:.3e} (soglia 1e-4; fine {ibr['staz_fine']:.3e})")
    jf, jp = ibr["j"][-1][2], pav["j"][-1][2]
    chk("C3_corrente_fuori", jf <= 10 * jp, f"j_fuori al passo 240 = {jf:.3e} contro 10 x pavimento = {10 * jp:.3e} ({jf / jp:.0f} x) [previsione: FAIL]")
    num["C3b_j_coperte"] = ibr["j"][-1][1]

# ------------------------------------------------------------------ V: convergenza nel blocco
V = {}
for ni in (3, 5, 7, 9, 11):
    r = risultati(rd(f"V_{ni}.log"))
    if r["finestre"]:
        w = [f for f in r["finestre"] if f[0] == 72]
        V[ni] = dict(E_eV=w[0][2] if w else r["finestre"][-1][2], ibrido_su_grosso=r["ibrido_su_grosso"], staz_fuori=r["staz_fuori"])
num["V"] = V
if len(V) == 5:
    E = {k: v["E_eV"] for k, v in V.items()}
    ks = sorted(E)
    passi = [(E[b] - E[a]) / abs(E[a]) for a, b in zip(ks, ks[1:])]       # < 0: more bound
    chk("V1_monotona", all(p <= 1e-5 for p in passi), "variazioni relative 3->5->7->9->11: " + ", ".join(f"{p:+.2e}" for p in passi) + " (negativo = piu' legato)")
    d119, d97 = abs(E[11] - E[9]) / abs(E[9]), abs(E[9] - E[7]) / abs(E[7])
    chk("V2_plateau", d119 < 3e-5 and d97 < 1e-4, f"|E11-E9|/|E| = {d119:.2e} (soglia 3e-5), |E9-E7|/|E| = {d97:.2e} (soglia 1e-4)")
    pk = abs(E[7] - 0.5 * (E[5] + E[9])) / abs(E[7])
    chk("V3_picco", pk < 1e-4, f"|E7-(E5+E9)/2|/|E| = {pk:.2e}")

# ------------------------------------------------------------------ R: energia N = 256
R = risultati(rd("R.log"))
if R["finestre"]:
    num["R"] = R
    chk("C4_R_energia", abs(R["E_rel_pct"]) <= 0.10, f"E1 ibrido N=256 blocco 9 = {R['E_eV']:.6f} eV ({R['E_rel_pct']:+.4f} % da Dirac x mu; soglia 0.10 %)")

# ------------------------------------------------------------------ M: massa del maestro
rg = json.load(open(os.path.join(OUT, "idrogeno_bohr_fase4_regola.json")))
ms = rg["masse_per_j"]
num["M_masse_per_j"] = ms
js = sorted(int(j) for j in ms)
dm = {j: abs(ms[str(j)] - ms[str(j - 1)]) for j in js[1:]}
chk("M1_convergenza", rg["maestro_j"] <= 4 and dm[rg["maestro_j"]] < 1e-5, f"maestro scelto j = {rg['maestro_j']}; |dm| per j: " + ", ".join(f"{j}:{d:.2e}" for j, d in dm.items()))
lvs = rg["chain_A_3d"]["levels"]
mm = [l["controlli"]["m_assoluta"] - 1 for l in lvs]
chk("M2_massa", max(mm) <= 2e-5, "m - 1 (a0/16, a0/8) = " + ", ".join(f"{x:.3e}" for x in mm) + " (prima 3.73e-4)")
okk = all(l["controlli"]["kcheck"]["packet"] <= 1e-3 and l["controlli"]["kcheck"]["energy_rel"] <= 3e-4 and l["controlli"]["kcheck"]["rest_rel"] <= 3e-4
          and l["controlli"]["kcheck"]["bz_rel90"] <= 0.2 and l["controlli"]["kcheck"]["neg_gen"] <= 0.01 for l in lvs)
chk("M4_kcheck", okk, "; ".join(f"B={l['B_fine_cells']}: w {l['w']}, packet {l['controlli']['kcheck']['packet']:.1e}, energia {l['controlli']['kcheck']['energy_rel']:.1e}, bz90 {l['controlli']['kcheck']['bz_rel90']:.2f}" for l in lvs))
p2 = []
for t in ("000", "111"):
    m = re.search(r"2p: E_rule = (\S+) eV vs Dirac×mu centroid (\S+): (\S+) %", rd(f"m3_2p_{t}.log"))
    if m:
        p2.append((float(m.group(1)), float(m.group(2)), float(m.group(3))))
if len(p2) == 2:
    Ea = 0.5 * (p2[0][0] + p2[1][0])
    r = 100 * (Ea / p2[0][1] - 1)
    num["M3_2p"] = dict(E_000=p2[0][0], E_111=p2[1][0], rel_000_pct=p2[0][2], rel_111_pct=p2[1][2], media_eV=Ea, media_rel_pct=r)
    chk("M3_2p_rust", abs(r) <= 0.005, f"2p (Rust, a0/8, N=256, media bordi 000/111) {Ea:.6f} eV, {r:+.5f} % da Dirac x mu (000 {p2[0][2]:+.4f} %, 111 {p2[1][2]:+.4f} %; prima +0.0338 %)")

G = os.path.join(OUT, "idrogeno_bohr_fase4_grosso_esame.json")
if os.path.exists(G):
    g = json.load(open(G))
    for k, v in g["esami"].items():
        if "passed" in v:
            ex[k] = dict(passed=v["passed"], dettaglio={kk: vv for kk, vv in v.items() if kk != "passed"})
    num["G_livelli"] = g["livelli"]
json.dump(out, open(os.path.join(OUT, "idrogeno_bohr_fase4_esame.json"), "w"), indent=1, default=float)
print("PASS:", [k for k, v in ex.items() if v["passed"]])
print("FAIL:", [k for k, v in ex.items() if not v["passed"]])
