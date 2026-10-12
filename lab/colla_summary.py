"""colla_summary.py - evaluates the pre-registered COLLA exams (colla_prereg.json) from the result files and writes
lab/out/colla_results.json (+ prints the table). Shape symmetry scores are taken in the best-fit orientation
(the noise of the start can leave the knot slightly rotated: rotations are zero modes): for each candidate group
we minimise the summed asym over a global rotation R (Nelder-Mead on the rotation vector) and report the asyms there."""
import json, os, math, glob, numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
import nuclei_core as nc
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
H = 0.2
ld = lambda f: json.load(open(os.path.join(OUT, f))) if os.path.exists(os.path.join(OUT, f)) else None
units = ld("colla_units.json"); AN = units["AN"]; Eu = AN["E_unit"]
EXP = {"d": 2.2246, "3He": 7.718, "4He": 28.296}
NAME = {2: "d", 3: "3He", 4: "4He"}
GHS_E1 = 938.27


def run(job):
    return ld(f"colla_{job}.json")


GROUPS = {
    "tetra": dict(C2x=nc.rot([1, 0, 0], math.pi), C2y=nc.rot([0, 1, 0], math.pi), C2z=nc.rot([0, 0, 1], math.pi),
                  C3=nc.rot([1, 1, 1], 2 * math.pi / 3)),
    "cube": dict(Cz90=nc.rot([0, 0, 1], math.pi / 2), C3=nc.rot([1, 1, 1], 2 * math.pi / 3), inv=-np.eye(3)),
}


def sym_fit(job, group):
    d = np.load(os.path.join(OUT, f"colla_fields_{job}.npz")); w = np.clip(d["b"].astype(float), 0, None)
    _, c, _, _ = nc.moments(w, H)
    G = GROUPS[group]

    def cost(rv):
        R = Rotation.from_rotvec(rv).as_matrix()
        return sum(nc.asym(w, H, c, R @ g @ R.T) for g in G.values())
    best = None
    for s in ([0, 0, 0], [0.2, -0.2, 0.3]):
        r = minimize(cost, s, method="Nelder-Mead", options=dict(xatol=1e-3, fatol=1e-5, maxiter=150))
        if best is None or r.fun < best.fun:
            best = r
    R = Rotation.from_rotvec(best.x).as_matrix()
    out = {k: nc.asym(w, H, c, R @ g @ R.T) for k, g in G.items()}
    out["inv"] = nc.asym(w, H, c, -np.eye(3))
    out["rot_deg"] = float(np.degrees(np.linalg.norm(best.x)))
    return out


def E(job):
    r = run(job); return None if r is None else r["final"]["E"]


def best(jobs):
    rs = [(E(j), j) for j in jobs if E(j) is not None]
    return min(rs) if rs else (None, None)


res = dict(exams={}, numbers={})
X = res["exams"]; N = res["numbers"]


def exam(k, ok, detail):
    X[k] = dict(verdict=("PASS" if ok else "FAIL") if ok is not None else "N/A", detail=detail)
    print(f"{k:8s} {X[k]['verdict']:5s} {detail}")


# ------------------------------------------------------------------ C1
exam("C1a", abs(AN["e"] / 4.84 - 1) < 0.02 and abs(AN["F_pi"] / 108 - 1) < 0.02,
     f"e = {AN['e']:.4f} (4.84), F_pi = {AN['F_pi']:.2f} MeV (108), m = {AN['m']:.4f}")
E1l, E1d = E("an_B1"), units["hedgehog"]["an"]["E"]
exam("C1b", abs(E1l / E1d - 1) < 0.01, f"lattice E1 = {E1l:.4f} vs 1D {E1d:.4f} ({100*(E1l/E1d-1):+.3f} %)")
sh = {}
r2 = run("an_B2")["final"]
sh[2] = dict(Cax=[r2["asym_axis"][f"Cax{d}"] for d in (30, 45, 90)], ring_center_over_max=r2["ring_center_over_max"])
sh[3] = sym_fit("an_B3", "tetra"); sh[4] = sym_fit("an_B4", "cube")
ok2 = max(sh[2]["Cax"]) < 0.02 and sh[2]["ring_center_over_max"] < 0.1
ok3 = max(sh[3][k] for k in ("C2x", "C2y", "C2z", "C3")) < 0.02 and sh[3]["inv"] > 0.3
ok4 = max(sh[4][k] for k in ("Cz90", "C3", "inv")) < 0.02
exam("C1c", ok2 and ok3 and ok4, f"B2 Cax30/45/90 = {np.round(sh[2]['Cax'], 4).tolist()}, centre/max {sh[2]['ring_center_over_max']:.3f}; "
     f"B3 {({k: round(v, 4) for k, v in sh[3].items()})}; B4 {({k: round(v, 4) for k, v in sh[4].items()})}")
N["shapes_AN"] = sh
cl = {B: (B * E1l - E(f"an_B{B}")) * Eu for B in (2, 3, 4)}
N["AN_classical_MeV"] = {NAME[B]: cl[B] for B in cl}
N["AN_E_over_12pi2B"] = {B: E(f"an_B{B}") / (12 * math.pi**2 * B) for B in (1, 2, 3, 4)}
exam("C1d", cl[2] > 10 * EXP["d"] and cl[4] > 3 * EXP["4He"],
     f"classical binding (AN units): d {cl[2]:.1f}, 3He {cl[3]:.1f}, 4He {cl[4]:.1f} MeV")
e0 = E("ghs0_B4")
exam("C1e", None if e0 is None else abs(e0 / (12 * math.pi**2 * 4) / 1.307 - 1) < 0.01,
     "n/a" if e0 is None else f"m=1 cube E/12pi^2B = {e0/(48*math.pi**2):.4f} (BMS 2007: 1.307)")

# ------------------------------------------------------------------ C2a
q = ld("colla_quant.json")
if q:
    b1 = q["B"]["1"]; Lam1 = AN["Lam"]
    dev = max(abs(np.diag(np.array(b1[t])) * s / Lam1 - 1).max() for t, s in (("U", 1), ("V", 1), ("W", -1)))
    exam("C2a-i", dev < 0.03, f"lattice U,V,-W diag vs 1D Lambda = {Lam1:.3f}: max dev {100*dev:.2f} %  (Lambda_lat = {b1['Lam']:.3f})")
    c = q["B"]["2"]["checks"]
    ok = c["U11_U22"] < 0.01 and c["V11_V22"] < 0.01 and abs(c["V33_over_U33"] / 4 - 1) < 0.03 and abs(c["W33_over_U33"] / 2 - 1) < 0.03 and c["W11_over_U11"] < 0.02
    exam("C2a-ii", ok, json.dumps({k: round(v, 4) for k, v in c.items()}))
    bq = q["binding"]
    exam("C2a-iii", bq["d"]["quantised"] > bq["d"]["classical"] and bq["4He"]["quantised"] > bq["4He"]["classical"],
         f"d {bq['d']['classical']:.1f} -> {bq['d']['quantised']:.1f}; 3He {bq['3He']['classical']:.1f} -> {bq['3He']['quantised']:.1f}; "
         f"4He {bq['4He']['classical']:.1f} -> {bq['4He']['quantised']:.1f} MeV (M_N = {bq['M_N_MeV']:.1f}, M_Delta = {bq['M_Delta_MeV']:.1f})")
    exam("C2a-iv", 79 < bq["d"]["quantised"] < 316, f"d rigid rotor {bq['d']['quantised']:.1f} MeV vs Braaten-Carson 158 MeV")
    N["C2a"] = bq

# ------------------------------------------------------------------ C2b
m1, m2 = ld("colla_modes_an_B1.json"), ld("colla_modes_an_B2.json")
if m1 and m2 and q:
    hb = AN["hbar"]; m = AN["m"]
    w1 = np.array(m1["omega2_over_m2"]); w2 = np.array(m2["omega2_over_m2"])
    z1, z2 = int((w1 < 0.05).sum()), int((w2 < 0.05).sum())
    exam("C2b-i", z1 == 6 and z2 == 8, f"zero modes B1 {z1}, B2 {z2}; omega^2/m^2 B1 {np.round(w1, 3).tolist()}; B2 {np.round(w2, 3).tolist()}")
    b1v = w1[(w1 >= 0.05) & (w1 < 1)]; b2v = w2[(w2 >= 0.05) & (w2 < 1)]
    zp1 = 0.5 * hb * m * np.sqrt(b1v).sum() * Eu; zp2 = 0.5 * hb * m * np.sqrt(b2v).sum() * Eu
    dzp = zp2 - 2 * zp1
    reached = bool(w1.max() >= 1 and w2.max() >= 1)
    d_after = q["binding"]["d"]["quantised"] - dzp
    exam("C2b-ii", dzp > 0 and d_after > 10 * EXP["d"],
         f"bound modes B1 {len(b1v)}, B2 {len(b2v)} (spectrum reached threshold: {reached}); zero-point B1 {zp1:.1f}, B2 {zp2:.1f} MeV; "
         f"Delta = {dzp:.1f} MeV -> d binding after C2a+C2b = {d_after:.1f} MeV")
    N["C2b"] = dict(zp1=zp1, zp2=zp2, delta=dzp, d_after=d_after, bound_B1=b1v.tolist(), bound_B2=b2v.tolist(), reached=reached)

# ------------------------------------------------------------------ C3a
alphas = {0.0: "ghs0", 0.5: "ghs50", 0.8: "ghs80", 0.95: "ghs95"}
scan = {}
for a, mo in alphas.items():
    e1 = E(f"{mo}_B1")
    if e1 is None:
        continue
    row = dict(E1=e1)
    for B in (2, 3, 4):
        eb, jb = best([f"{mo}_B{B}", f"{mo}_B{B}pair", f"{mo}_B{B}tet"])
        if eb is None:
            continue
        frac = (B * e1 - eb) / (B * e1)
        row[NAME[B]] = dict(job=jb, frac=frac, MeV=frac * B * GHS_E1, per_nucleon=frac * GHS_E1)
        f = run(jb)["final"]
        row[NAME[B]].update(n_peaks=f["n_peaks"], ring_center_over_max=f["ring_center_over_max"], q=f["q"])
    scan[a] = row
N["C3a_scan"] = scan
if 0.95 in scan:
    s = scan[0.95]
    e81 = units["hedgehog"]["ghs95"]["E_over_8pi2"]
    ok = abs(e81 / 1.1082 - 1) < 0.01
    det = f"1D E1/8pi^2 = {e81:.4f} (1.1082)"
    for nm, lit in (("d", 1.84), ("4He", 4.56)):
        if nm in s:
            ok = ok and abs(s[nm]["per_nucleon"] / lit - 1) < 0.3
            det += f"; {nm} BEPN = {s[nm]['per_nucleon']:.2f} MeV ({lit})"
    exam("C3a-i", ok, det)
fr = {nm: [scan[a][nm]["frac"] for a in sorted(scan) if nm in scan[a]] for nm in ("d", "4He")}
exam("C3a-ii", all(np.all(np.diff(v) < 0) for v in fr.values()) and all(len(v) == 4 for v in fr.values()),
     f"binding fractions vs alpha {sorted(scan)}: " + "; ".join(f"{k} {np.round(np.array(v)*100, 3).tolist()} %" for k, v in fr.items()))
shp = {}
for a, mo in alphas.items():
    for B in (2, 4):
        eb, jb = best([f"{mo}_B{B}", f"{mo}_B{B}pair", f"{mo}_B{B}tet"])
        if jb is None:
            continue
        f = run(jb)["final"]
        dd = dict(job=jb, n_peaks=f["n_peaks"], ring_center_over_max=f["ring_center_over_max"], Cax=[f["asym_axis"][f"Cax{x}"] for x in (30, 45, 90)])
        if B == 4:
            dd.update(sym_fit(jb, "cube"))
        shp[f"{a}_B{B}"] = dd
N["C3a_shapes"] = shp
try:
    lo = shp["0.5_B2"]; lo4 = shp["0.5_B4"]; hi = shp["0.95_B2"]; hi4 = shp["0.95_B4"]
    kept = max(lo["Cax"]) < 0.02 and lo["ring_center_over_max"] < 0.1 and max(lo4[k] for k in ("Cz90", "C3", "inv")) < 0.02
    split = (hi["ring_center_over_max"] > 0.3 or hi["n_peaks"] == 2) and (hi4["inv"] > 0.1 or hi4["n_peaks"] == 4)
    exam("C3a-iii", kept and split, f"alpha=0.5 kept: {kept} (B2 Cax {np.round(lo['Cax'],3).tolist()}, c/max {lo['ring_center_over_max']:.3f}; "
         f"B4 Cz90/C3/inv {lo4['Cz90']:.3f}/{lo4['C3']:.3f}/{lo4['inv']:.3f}); alpha=0.95 split: {split} (B2 c/max {hi['ring_center_over_max']:.3f}, "
         f"peaks {hi['n_peaks']}; B4 inv {hi4['inv']:.3f}, peaks {hi4['n_peaks']})")
except KeyError as ex:
    exam("C3a-iii", None, f"missing {ex}")
# alpha*: log-linear interpolation of the 4He binding (MeV) in alpha
al = sorted(a for a in scan if "4He" in scan[a])
star = None
for a0, a1 in zip(al[:-1], al[1:]):
    y0, y1 = scan[a0]["4He"]["MeV"], scan[a1]["4He"]["MeV"]
    if (y0 - EXP["4He"]) * (y1 - EXP["4He"]) <= 0 and y0 > 0 and y1 > 0:
        t = (math.log(EXP["4He"]) - math.log(y0)) / (math.log(y1) - math.log(y0))
        star = dict(alpha=a0 + t * (a1 - a0), bracket=[a0, a1], t=t)
        for nm in ("d", "3He"):
            if nm in scan[a0] and nm in scan[a1] and scan[a0][nm]["MeV"] > 0 and scan[a1][nm]["MeV"] > 0:
                star[nm] = math.exp(math.log(scan[a0][nm]["MeV"]) + t * (math.log(scan[a1][nm]["MeV"]) - math.log(scan[a0][nm]["MeV"])))
        break
N["C3a_alpha_star"] = star
if star:
    okd = "d" in star and 0.5 < star["d"] / EXP["d"] < 2
    okh = "3He" in star and 0.5 < star["3He"] / EXP["3He"] < 2
    exam("C3a-iv", okd and okh, f"alpha* = {star['alpha']:.4f} (FIT to 4He = 28.3 MeV); predicted d = {star.get('d', float('nan')):.2f} MeV (2.22), "
         f"3He = {star.get('3He', float('nan')):.2f} MeV (7.72)")
else:
    ys = {a: round(scan[a]["4He"]["MeV"], 2) for a in al}
    exam("C3a-iv", None, f"4He binding not bracketed by the scan: {ys} - alpha* outside [0, 0.95]")

# post-hoc (not an exam): direct runs at alpha* (model ghsS in colla_units.json) to check the interpolation
e1s = E("ghsS_B1")
if e1s is not None:
    ds = dict(alpha=units["models"]["ghsS"]["alpha"])
    for B in (2, 3, 4):
        eb, jb = best([f"ghsS_B{B}"])
        if eb is not None:
            fr_ = (B * e1s - eb) / (B * e1s)
            f = run(jb)["final"]
            ds[NAME[B]] = dict(frac=fr_, MeV=fr_ * B * GHS_E1, n_peaks=f["n_peaks"], q=f["q"])
    N["C3a_direct_alpha_star"] = ds
    print("post-hoc alpha* direct:", {k: (round(v["MeV"], 2) if isinstance(v, dict) else v) for k, v in ds.items()})

# ------------------------------------------------------------------ C3b
sx = {}
for c6, mo in ((0, "an"), (50, "sx50"), (200, "sx200")):
    e1 = E(f"{mo}_B1")
    if e1 is None:
        continue
    sx[c6] = {NAME[B]: dict(frac=(B * e1 - E(f"{mo}_B{B}")) / (B * e1), MeV=(B * e1 - E(f"{mo}_B{B}")) * Eu)
              for B in (2, 4) if E(f"{mo}_B{B}") is not None}
N["C3b"] = sx
if 200 in sx and 50 in sx:
    mono = all(sx[0][k]["frac"] > sx[50][k]["frac"] > sx[200][k]["frac"] for k in ("d", "4He"))
    exam("C3b-i", mono and sx[200]["4He"]["MeV"] > 3 * EXP["4He"],
         "; ".join(f"c6={c}: d {v['d']['MeV']:.1f} MeV ({100*v['d']['frac']:.2f} %), 4He {v['4He']['MeV']:.1f} MeV ({100*v['4He']['frac']:.2f} %)" for c, v in sx.items()))
json.dump(res, open(os.path.join(OUT, "colla_results.json"), "w"), indent=1, default=float)
