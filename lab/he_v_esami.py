"""he_v_esami.py - exams E1-E5 of VIDEO ELIO (pre-registration lab/out/he_v_prereg.json) -> lab/out/he_v_esami.json"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import he_v_render as RV

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
D = RV.dati()
pairs = D["pairs"]; d = D["dens"]
es = {}
# E1: radial distribution of all balls vs the exact density
r = np.linalg.norm(pairs.reshape(-1, 3), axis=1)
edges = np.arange(0, 6.0001, 0.05)
h, _ = np.histogram(r, edges)
p_ball = h / len(r)
rr, rho = d["rho1_esatta_r"], d["rho1_esatta"]
fine = np.linspace(0, 6, 6001)
f4 = 4 * np.pi * fine ** 2 * np.interp(fine, rr, rho)
cum = np.concatenate([[0], np.cumsum(0.5 * (f4[1:] + f4[:-1]) * np.diff(fine))])
p_ex = np.diff(np.interp(edges, fine, cum))
L1 = float(np.abs(p_ball - p_ex).sum())
es["E1"] = dict(L1=L1, n_balls=int(len(r)), frac_beyond_6=float(np.mean(r > 6)), PASS=bool(L1 <= 0.03))
# E2
r12 = np.linalg.norm(pairs[:, 0] - pairs[:, 1], axis=1)
mr, mr12 = float(r.mean()), float(r12.mean())
es["E2"] = dict(mean_r=mr, ref_r=0.929472, rel_r=mr / 0.929472 - 1, mean_r12=mr12, ref_r12=1.42207, rel_r12=mr12 / 1.42207 - 1,
                E_Eh_quoted=RV.MAN["sorgente"]["E_Eh"], E_exact=RV.MAN["sorgente"]["E_esatta_Eh"],
                PASS=bool(abs(mr / 0.929472 - 1) <= 5e-3 and abs(mr12 / 1.42207 - 1) <= 5e-3))
# E3: Coulomb hole vs shuffled partners
rng = np.random.default_rng(5)
sh = np.linalg.norm(pairs[:, 0] - pairs[rng.permutation(len(pairs)), 1], axis=1)
Pr, Ps = float(np.mean(r12 < 0.5)), float(np.mean(sh < 0.5))
es["E3"] = dict(P_real=Pr, P_shuffled=Ps, ratio=Pr / Ps, mean_r12_shuffled=float(sh.mean()), PASS=bool(Pr / Ps <= 0.85))
# E4: g from the pairs themselves (rotated so that r1 -> +z) vs the map used for the light
both = np.concatenate([pairs, pairs[:, ::-1]])
n1 = np.linalg.norm(both[:, 0], axis=1)
sel = both[(n1 >= 0.9) & (n1 <= 1.1)]
u = sel[:, 0] / np.linalg.norm(sel[:, 0], axis=1)[:, None]
z = np.sum(sel[:, 1] * u, 1)
rho_c = np.sqrt(np.maximum(np.sum(sel[:, 1] ** 2, 1) - z ** 2, 0))
eb = np.arange(-3, 3.0001, 0.25); er = np.arange(0, 3.0001, 0.25)
H, _, _ = np.histogram2d(z, rho_c, [eb, er])
zc, rc = 0.5 * (eb[1:] + eb[:-1]), 0.5 * (er[1:] + er[:-1])
ZZ, RR = np.meshgrid(zc, rc, indexing="ij")
vol = 2 * np.pi * RR * 0.25 * 0.25
cond = H / len(sel) / vol
rho1 = np.interp(np.sqrt(ZZ ** 2 + RR ** 2), rr, rho)
g_pairs = cond / rho1
g_map = D["g_ip"](np.stack([ZZ.ravel(), RR.ravel()], 1)).reshape(ZZ.shape)
ok = H >= 200
pear = float(np.corrcoef(g_pairs[ok], g_map[ok])[0, 1])
med = float(np.median(np.abs(g_pairs[ok] - g_map[ok])))
es["E4"] = dict(n_selected=int(len(sel)), n_cells=int(ok.sum()), pearson=pear, median_abs_diff=med,
                g_map_min=float(g_map.min()), PASS=bool(pear >= 0.9 and med <= 0.05))
np.save(os.path.join(RV.SCR, "he_v_E4.npy"), np.stack([g_pairs, g_map, ok.astype(float)]))
# E5: mean of g over the balls in every B frame
gm = []
for f in range(RV.NA, RV.NF):
    dd = RV.carica(f)
    gm.append(float(RV.g_of(np.array(dd["meta"]["ref"]), dd["P"]).mean()))
gm = np.array(gm)
es["E5"] = dict(min=float(gm.min()), max=float(gm.max()), max_dev=float(np.max(np.abs(gm - 1))), PASS=bool(np.max(np.abs(gm - 1)) <= 0.03))
json.dump(es, open(os.path.join(OUT, "he_v_esami.json"), "w"), indent=1)
for k, v in es.items():
    print(k, v)
