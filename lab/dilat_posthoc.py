"""dilat_posthoc.py - CONTROLLI POST-HOC (dichiarati, dopo aver visto i risultati delle corse pre-registrate).
Non cambiano i verdetti PASS/FAIL della pre-registrazione: li spiegano. Scrive lab/out/dilat_emendamenti.json.

A1  E1a, E2b, E6 falliscono la soglia di precisione a a = 1/32, ma l'errore va esattamente come a^2 (rapporti 4,00):
    estrapolazione di Richardson (a = 1/16, 1/32) -> confronto col continuo.
A2  E2a, phi = pi/2: il residuo 3e-12 e' il pavimento numerico della fase (delta = 1e-6). Rifatto con la media
    spettrale esatta <r(k)> su |g|^2 (stessa regola, nessun arg).
A3  E3: la previsione usava il continuo SCALARE (H = sqrt(p^2 + M^2), come Smith-Ahmadi/PPK). Per una particella di
    Dirac lo spinore di energia positiva dipende dalla massa: lo stato preparato per la massa m non e' di energia
    positiva per m +- delta/2 -> transitorio di Zitterbewegung nella fase dell'orologio. Si confronta (i) la regola
    con il continuo di DIRAC con la stessa preparazione, (ii) il ritmo asintotico (fit della fase su t in [50, 250])
    con lo scalare e con Paige-Plato-Kim.
"""
import json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dilat_core import (two_packets, u_plus_lattice, u_plus_dirac, spinor_state, overlap_spec, mean_rate)

HERE = os.path.dirname(os.path.abspath(__file__))
R = lambda n: json.load(open(os.path.join(HERE, "out", "dilat_r_%s.json" % n)))
PRE = json.load(open(os.path.join(HERE, "out", "dilat_prereg.json")))["esami"]
m = 1.0
out = {}

# ---------------------------------------------------------------- A1
e1 = R("E1")
rows = []
for r in e1["righe"]:
    a16, a32 = r["R"][str(1 / 16)], r["R"][str(1 / 32)]
    ex = (4 * a32 - a16) / 3
    rows.append(dict(pbar=r["pbar"], R_rich=ex, R_cont=r["R_cont"], rel=abs(ex - r["R_cont"]) / r["R_cont"]))
e2 = R("E2")
rows2 = []
for r in e2["E2b"]["righe"]:
    ex = (4 * r["dq"][str(1 / 32)] - r["dq"][str(1 / 16)]) / 3
    rows2.append(dict(s=r["s"], alpha=r["alpha"], phi=r["phi"], dq_rich=ex, dq_cont=r["dq_cont"],
                      rel=abs(ex - r["dq_cont"]) / abs(r["dq_cont"]), dq_SA_NR=r["dq_SA_NR"]))
e6 = R("E6")
t6 = np.array(e6["tau"][str(1 / 32)])
t6b = np.array(e6["tau"][str(1 / 16)])
prev6 = np.array(PRE["E6"]["previsione"]["tau_cont"])
rich6 = (4 * t6 - t6b) / 3
out["A1"] = dict(E1=rows, E1_max_rel=max(x["rel"] for x in rows), E2b=rows2, E2b_max_rel=max(x["rel"] for x in rows2),
                 E6_max_diff_rich=float(np.max(np.abs(rich6 - prev6))))
print("A1 E1 max rel %.2e  E2b max rel %.2e  E6 max diff %.2e" % (out["A1"]["E1_max_rel"], out["A1"]["E2b_max_rel"],
                                                                 out["A1"]["E6_max_diff_rich"]))

# ---------------------------------------------------------------- A2
a = 1 / 32
k = 2 * np.pi * np.fft.fftfreq(int(round(1024 / a)), d=a)
res = []
for pr in PRE["E2a"]["previsione"]:
    if abs(pr["phi"] - np.pi / 2) > 1e-9:
        continue
    g, g1, g2 = two_packets(k, np.cos(pr["alpha"]), np.sin(pr["alpha"]), pr["phi"], pr["p1"], pr["p2"], 0.02)
    dq = mean_rate(g, k, a) - (np.cos(pr["alpha"]) ** 2 * mean_rate(g1, k, a) + np.sin(pr["alpha"]) ** 2 * mean_rate(g2, k, a))
    res.append(abs(dq))
out["A2"] = dict(max_abs_dq_phi_pi2=max(res))
print("A2 phi=pi/2 max |dq| (media spettrale) %.2e" % max(res))

# ---------------------------------------------------------------- A3
pg = 2 * np.pi * np.fft.fftfreq(2 ** 16, d=4096 / 2 ** 16)
rowsA3 = []
for pb, dl in [(0.05, 0.2), (0.05, 0.5), (0.1, 0.2), (0.1, 0.5), (0.3, 0.2), (0.3, 0.5)]:
    gc, _, _ = two_packets(pg, 1, 0, 0, pb, pb, 0.01)
    w = np.abs(gc) ** 2 / np.sum(np.abs(gc) ** 2)
    Mp, Mm = m + dl / 2, m - dl / 2
    D_scal = 1 - float(np.sum(w * (np.sqrt(pg ** 2 + Mp ** 2) - np.sqrt(pg ** 2 + Mm ** 2))) / dl)
    D_ppk = float(np.sum(w * pg ** 2 / (2 * Mp * Mm)))
    D_mean = float(np.sum(w * pg ** 2 / 2))
    Od = overlap_spec(spinor_state(gc, u_plus_dirac(pg, m)), pg, [2.0], m, dl, "dirac")[0]
    D_dir2 = 1 + np.angle(Od) / (dl * 2.0)
    row = dict(pbar=pb, delta=dl, D_scal=D_scal, D_PPK=D_ppk, D_media=D_mean, D_dirac_t2=D_dir2, regola_t2={}, fit={})
    ts = np.arange(0, 250.0 + 1e-9, 0.25)
    for a in [1 / 8, 1 / 16, 1 / 32]:
        k = 2 * np.pi * np.fft.fftfreq(int(round(1024 / a)), d=a)
        g, _, _ = two_packets(k, 1, 0, 0, pb, pb, 0.01)
        psi = spinor_state(g, u_plus_lattice(k, a, m))
        ns = [int(round(t / a)) for t in ts]
        O = overlap_spec(psi, k, ns, m, dl, "lattice", a)
        ph = np.unwrap(np.angle(O))
        row["regola_t2"][str(a)] = 1 + ph[8] / (dl * 2.0)
        sel = ts >= 50
        slope = np.polyfit(ts[sel], ph[sel], 1)[0]
        row["fit"][str(a)] = 1 + slope / dl
        if a == 1 / 32:
            row["fase_transitoria"] = (ph + dl * (1 - D_scal) * ts)[:200:8].tolist()   # fase - parte lineare scalare
    f16, f32 = row["fit"][str(1 / 16)], row["fit"][str(1 / 32)]
    row["fit_rich"] = (4 * f32 - f16) / 3
    row["rel_rich_scal"] = abs(row["fit_rich"] - D_scal) / D_scal
    row["rel_t2_dirac"] = abs(row["regola_t2"][str(1 / 32)] - D_dir2) / D_dir2
    rowsA3.append(row)
    print("A3 p=%.2f d=%.1f  t=2: regola %.6e Dirac %.6e | fit asintotico (Rich) %.7e scalare %.7e PPK %.7e media %.7e" % (
        pb, dl, row["regola_t2"][str(1 / 32)], D_dir2, row["fit_rich"], D_scal, D_ppk, D_mean))
sel = [r for r in rowsA3 if r["pbar"] == 0.05 and r["delta"] == 0.5][0]
out["A3"] = dict(righe=rowsA3, test_PPK_asintotico=dict(dPPK=abs(sel["fit_rich"] - sel["D_PPK"]),
                                                        gap=abs(sel["D_PPK"] - sel["D_media"])))
print("A3 PPK asintotico: |D - PPK| = %.2e   |PPK - media| = %.2e" % (out["A3"]["test_PPK_asintotico"]["dPPK"],
                                                                      out["A3"]["test_PPK_asintotico"]["gap"]))
with open(os.path.join(HERE, "out", "dilat_emendamenti.json"), "w") as f:
    json.dump(out, f, indent=1, default=float)
