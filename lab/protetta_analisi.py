"""protetta_analisi.py - exams D1-D3 from lab/out/protetta_vestito.json, gap table and figure.
-> lab/out/protetta_esami.json, lab/out/protetta_figura.png.   Run: python3 lab/protetta_analisi.py
Debiti numbers (bare light walker, origin/master 35afcc5, lab/out/vuoto_f3_exams.json / VUOTO_POLARIZZA.md par. 14.2) are
quoted, not recomputed: isotropic dv(a->0) = -0.0668 (R, K12), +0.0455 (R, K12/3), +0.0444 (C, K12/3), +0.0531 (continuum)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from zz_core import ALPHA

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
D = json.load(open(os.path.join(OUT, "protetta_vestito.json")))


def get(rule, cut, th, pd):
    return D.get("%s|%s|%g|%g" % (rule, cut, th, pd))


ex = {"D1": {}, "D2": {}, "D3": {}}
curves = {}
for cut in ("K12", "K12/3"):
    for ax, an in ((0, "x"), (1, "y")):
        for pd in (1, 3, 10):
            ths, ds = [], []
            for key, r in D.items():
                rule, c, th, p = key.split("|")
                if rule != "R" or c != cut or float(p) != pd:
                    continue
                th = float(th)
                if th < 5e-4:
                    continue
                r0 = get("R", cut, th, 0.05)
                if r0 is None:
                    continue
                ths.append(th); ds.append(r["dv_over_alpha"][ax] - r0["dv_over_alpha"][ax])
            o = np.argsort(ths); ths = np.array(ths)[o]; ds = np.array(ds)[o]
            if len(ths) >= 3:
                c2 = np.polyfit(ths, ds, 2)
                d0 = float(c2[-1]); ok = abs(d0) < 1e-3
            else:
                d0 = float(ds[0]) if len(ds) else None; ok = None
            curves[(cut, an, pd)] = (ths, ds)
            ex["D1"]["%s|%s|pd=%g" % (cut, an, pd)] = {"theta": ths.tolist(), "D_over_alpha": ds.tolist(), "d0": d0,
                                                        "n_points": int(len(ths)), "PASS": ok}
ex["D1"]["PASS"] = bool(all(v["PASS"] for v in ex["D1"].values() if isinstance(v, dict) and v["PASS"] is not None))
ex["D1"]["max_abs_D_all_points"] = float(max(np.max(np.abs(v["D_over_alpha"])) for v in ex["D1"].values()
                                               if isinstance(v, dict) and v.get("D_over_alpha")))
# D2 cyclic
rows = {}
for th in (0.003125, 0.0015625, 0.00078125):
    a, b = get("C", "K12/3", th, 10), get("C", "K12/3", th, 0.05)
    if a and b:
        rows["%g" % th] = {"dv_p<<m": b["dv_over_alpha"][0], "D": a["dv_over_alpha"][0] - b["dv_over_alpha"][0]}
ex["D2"] = {"rows": rows, "PASS": bool(rows and all(abs(v["D"]) < 1e-3 for v in rows.values()))}
# D3 massless limit
rows = {}
th = 3.90625e-4
for cut in ("K12", "K12/3"):
    b = get("R", cut, th, 0.05)
    for pd in (25.6, 51.2, 76.8):
        a = get("R", cut, th, pd)
        if a and b:
            for ax, an in ((0, "x"), (1, "y")):
                rows["%s|%s|p=%.3g" % (cut, an, pd * th)] = {"dv_p<<m": b["dv_over_alpha"][ax],
                                                            "D": a["dv_over_alpha"][ax] - b["dv_over_alpha"][ax]}
ex["D3"] = {"rows": rows, "PASS": bool(rows and all(abs(v["D"]) < 2e-3 for v in rows.values()))}

# gap table (c_gamma - c_e)/alpha, isotropic part (mean over x, y, z with z = x)
iso = lambda x, y: (2 * x + y) / 3
DEB = {"R K12": -0.0668, "R K12/3": 0.0455, "C K12/3": 0.0444, "continuo": 0.0531}
gap = {}
for k, v in DEB.items():
    gap["luce nuda | " + k] = -v
gap["composto neutro | R K12/3"] = -DEB["R K12/3"]
for cut in ("K12", "K12/3"):
    dx = ex["D1"]["%s|x|pd=10" % cut]["d0"]; dy = ex["D1"]["%s|y|pd=10" % cut]["d0"]
    gap["composto carico, solo autoenergia (<=, griglia) | R " + cut] = iso(dx, dy)
# D3 point as the massless constituent (p = 0.03)
for cut in ("K12", "K12/3"):
    r = ex["D3"]["rows"].get("%s|x|p=0.03" % cut), ex["D3"]["rows"].get("%s|y|p=0.03" % cut)
    if all(r):
        gap["composto carico, theta->0, p=77 m (<=) | R " + cut] = iso(r[0]["D"], r[1]["D"])
ex["gap_iso_over_alpha"] = gap
G = json.load(open(os.path.join(OUT, "protetta_griglia.json")))
em1 = {}
for cut in ("K12/3", "K12"):
    for th, pd in ((0.00078125, 10), (0.000390625, 76.8)):
        a, b = G["R|%s|%g|%g" % (cut, th, pd)], G["R|%s|%g|0.05" % (cut, th)]
        pa, pb = get("R", cut, th, pd), get("R", cut, th, 0.05)
        for ax, an in ((0, "x"), (1, "y")):
            em1["%s|%s|theta=%g|pd=%g" % (cut, an, th, pd)] = {
                "D_fine": a["dv_over_alpha"][ax] - b["dv_over_alpha"][ax],
                "D_production": pa["dv_over_alpha"][ax] - pb["dv_over_alpha"][ax]}
ex["EM1_grid"] = em1
ex["experimental_requirement_over_alpha"] = {"Crab electrons (indicative)": 2.5e-20 / ALPHA, "photon decay": 6e-21 / ALPHA,
                                             "LEP synchrotron": 5e-15 / ALPHA}
json.dump(ex, open(os.path.join(OUT, "protetta_esami.json"), "w"), indent=1)

# ---------------- figure
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
ax = axs[0]
names = list(gap)
vals = np.array([gap[n] for n in names])
y = np.arange(len(names))
cols = ["#888888"] * 4 + ["#c77"] + ["#3a7"] * (len(names) - 5)
ax.barh(y, np.abs(vals), color=cols)
for i, v in enumerate(vals):
    ax.text(max(abs(v), 1e-19) * 1.3, i, "%+.2g" % v, va="center", fontsize=8)
ax.axvline(2.5e-20 / ALPHA, color="k", ls="--", lw=1)
ax.text(2.5e-20 / ALPHA * 1.5, 4.0, "limite sperimentale\n(Granchio, 2.5e-20 / alpha)", fontsize=8)
ax.axvline(5e-15 / ALPHA, color="k", ls=":", lw=1)
ax.text(5e-15 / ALPHA * 1.5, -0.6, "LEP", fontsize=8)
ax.set_xscale("log"); ax.set_xlim(1e-19, 1)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8); ax.invert_yaxis()
ax.set_xlabel("|c_luce - c_elettrone| / alpha   (a -> 0)")
ax.set_title("(a) gap isotropo c_luce - c_e per variante (segno accanto)\nverde = fotone composto: solo autoenergie, termine di scambio NON calcolato", fontsize=9)
ax = axs[1]
mk = {1: "o", 3: "s", 10: "^"}
for (cut, an, pd), (ths, ds) in curves.items():
    if len(ths) == 0:
        continue
    ax.loglog(ths, np.abs(ds) + 1e-12, marker=mk[pd], ls="-" if cut == "K12" else "--",
              color={"x": "C0", "y": "C3"}[an], label="%s %s p=%g m" % (cut, an, pd), ms=4)
tt = np.array([7e-4, 1.3e-2])
ax.loglog(tt, 2e-3 * (tt / 1e-2) ** 2, "k:", label="~ theta^2")
ax.set_xlabel("theta = m a  (passo del reticolo a massa fissa)")
ax.set_ylabel("|D| / alpha = |dv(p) - dv(p<<m)| / alpha")
ax.set_title("(b) c vestito ad alto impulso = c vestito a riposo?", fontsize=10)
ax.legend(fontsize=6, ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "protetta_figura.png"), dpi=130)
print(json.dumps({k: v for k, v in ex.items() if k != "D1"}, indent=1))
print("D1 PASS", ex["D1"]["PASS"], "max|D|", ex["D1"]["max_abs_D_all_points"])
for k, v in ex["D1"].items():
    if isinstance(v, dict):
        print(k, v["n_points"], v["d0"], v["PASS"])
