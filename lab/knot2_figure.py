"""knot2_figure.py -> lab/out/knot2_figure.png (S1 dispersion, S2 boosted knot, S3 Lambda, S4 annihilation)"""
import os, json, glob, math, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
J = lambda n: json.load(open(os.path.join(OUT, n)))
ex = lambda n: os.path.exists(os.path.join(OUT, n))
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

# (a) dispersion
a = ax[0, 0]
kh = np.linspace(1e-3, np.pi, 400)
w2 = 2 * np.sin(kh / 2); w4 = (27 * np.sin(kh / 2) - np.sin(3 * kh / 2)) / 12
a.plot(kh, w2 / kh, "C0", label="2° ordine (vecchio, c2): fase")
a.plot(kh, np.cos(kh / 2), "C0--", label="2° ordine: gruppo")
a.plot(kh, w4 / kh, "C3", label="c4: fase")
a.plot(kh, (27 * np.cos(kh / 2) - 3 * np.cos(3 * kh / 2)) / 24, "C3--", label="c4: gruppo")
if ex("knot2_s1.json"):
    s = J("knot2_s1.json")
    for key, c in (("S1a_old_axis", "C0"), ("S1a_c4_axis", "C3")):
        if key in s:
            a.plot(s[key]["kh"], s[key]["v_phase"], "o", color=c, ms=3)
for v in (0.6, 0.8):
    a.axhline(v, color="gray", lw=0.7, ls=":"); a.text(0.1, v + 0.01, f"v = {v}", fontsize=8, color="gray")
a.set_xlabel("k h"); a.set_ylabel("velocità dell'onda / c"); a.set_ylim(0, 1.05)
a.set_title("S1: onde di pioni piccole sul reticolo\n(punti = misurate con la dinamica)", fontsize=10); a.legend(fontsize=7, loc="lower left")

# (b) boosted knot: lag x_c(t) - x0 - v t
a = ax[0, 1]
sty = {"old": ":", "c2": "--", "c4": "-"}
for v, col in ((0.3, "C2"), (0.6, "C1"), (0.8, "C3")):
    for f in sorted(glob.glob(os.path.join(OUT, f"knot2_boost_*_v{v:g}.json"))):
        d = json.load(open(f)); m = d["meta"]; r = d["records"]
        t = np.array([q["t"] for q in r]); x = np.array([q["b_cx"] for q in r])
        a.plot(t, x - x[0] - v * t, ls=sty[m["scheme"]], color=col, lw=1.8 if m["h"] == 0.15 else 1.0,
               label=f"{m['scheme']} h={m['h']:g} v={v:g}")
    for hh in (0.2, 0.15):
        n = f"knot_k2_h{hh:g}_v{v:g}.json"
        if ex(n) and v > 0.3:
            d = J(n); r = d["records"]; t = np.array([q["t"] for q in r]); x = np.array([q["b_cx"] for q in r])
            a.plot(t, x - x[0] - v * t, ls=":", color=col, lw=1.8 if hh == 0.15 else 1.0, label=f"old h={hh:g} v={v:g}")
a.axhline(0, color="k", lw=0.6)
a.set_xlabel("t"); a.set_ylabel("ritardo del centro  x(t) − x(0) − v t")
a.set_title("S2: il nodo che corre rimane indietro?\n(0 = corre esattamente a v)", fontsize=10); a.legend(fontsize=6, ncol=2)

# (c) Lambda
a = ax[0, 2]
if ex("knot2_s3.json"):
    s = J("knot2_s3.json")
    for sch, c in (("old", "C0"), ("c4", "C3")):
        for hh, mk in ((0.2, "o"), (0.15, "s")):
            rows = [v for k, v in s["rows"].items() if k.startswith(sch) and v["h"] == hh]
            if not rows: continue
            R = [1 / q["R_eff"] for q in rows]
            a.plot(R, [q["Lambda_lat"] for q in rows], mk, mfc="none", color=c, label=f"{sch} h={hh:g} grezzo")
            a.plot(R, [q["Lambda_corr"] for q in rows], mk, color=c, label=f"{sch} h={hh:g} + coda")
    a.axhline(106.77, color="k", lw=1); a.text(0.15, 107.5, "ODE 106.77", fontsize=8)
    a.axhspan(106.77 * 0.98, 106.77 * 1.02, color="gray", alpha=0.15)
a.set_xlabel("1 / R_eff (scatola)"); a.set_ylabel("Λ (momento d'inerzia)"); a.set_xlim(0, 0.32)
a.set_title("S3: Λ del reticolo con la coda 1/R aggiunta", fontsize=10); a.legend(fontsize=7)

# (d) energy out of r < 4.5
a = ax[1, 0]
for v, col in ((0.3, "C2"), (0.6, "C1")):
    n = f"knot2_s4_v{v:g}.json"
    if ex(n):
        d = J(n); r = d["records"]; t = np.array([q["t"] for q in r]); E = np.array([q["E_in"] for q in r])
        a.plot(t, 1 - E / E[0], color=col, lw=2, label=f"c4, scatola 80³, v={v:g}")
    n = f"knot_k3_SA_v{v:g}.json"
    if ex(n):
        d = J(n); r = d["records"]; t = np.array([q["t"] for q in r]); E = np.array([q["E_in"] for q in r])
        a.plot(t, 1 - E / E[0], color=col, ls=":", label=f"K3 (vecchio, 64³), v={v:g}")
a.axhline(0.7, color="k", ls="--", lw=0.8); a.axvline(24, color="gray", lw=0.6)
a.set_xlabel("t"); a.set_ylabel("frazione di energia uscita da r < 4.5")
a.set_title("S4: nodo + antinodo → onde (tempo lungo)", fontsize=10); a.legend(fontsize=7)

# (e) r-t map of the energy
a = ax[1, 1]
n = "knot2_s4_v0.6.json"
if ex(n):
    d = J(n); r = d["records"]; rb = np.array(d["rbins"])
    M = np.array([q["e_r"][1:len(rb)] for q in r]); t = np.array([q["t"] for q in r])
    rc = rb[:-1] + 0.1
    dEdr = M / 0.2                      # energy per unit radius: constant along r = t + const for free outgoing waves
    im = a.pcolormesh(rc, t, np.log10(dEdr + 1e-6), shading="auto", cmap="magma", vmin=-1, vmax=2.5)
    for t0 in (6, 10, 14):
        a.plot([1, 7.9], [t0, t0 + 6.9], "c--", lw=0.8)
    a.axvline(4.5, color="w", lw=0.6); a.axvline(6.5, color="w", lw=0.6, ls=":")
    plt.colorbar(im, ax=a, label="log10 energia per unità di raggio")
a.set_xlabel("r"); a.set_ylabel("t"); a.set_title("S4 (v = 0.6): l'energia esce a velocità c\n(linee tratteggiate: r = t + cost.)", fontsize=10)

# (f) spectrum
a = ax[1, 2]
for v, col in ((0.6, "C1"), (0.3, "C2")):
    n = f"knot2_s4_v{v:g}.json"
    if ex(n) and "spectral" in J(n):
        sp = J(n)["spectral"]; k = np.array(sp["k_center"]); w = np.array(sp["omega_est"]); Pd = np.array(sp["Pd"])
        use = np.isin(np.round(k, 6), np.round(sp["k_bins_used"], 6))
        a.scatter(k[use], w[use], s=8 + 200 * Pd[use] / Pd.max(), color=col, alpha=0.7, label=f"v={v:g}: ω stimata (90 % energia)")
        a.plot(k[k > 0], np.array(sp["omega_lat"])[k > 0], color=col, lw=0.6, alpha=0.6)
kk = np.linspace(0, 12, 50); a.plot(kk, kk, "k--", lw=1, label="ω = c k (senza massa)")
a.plot(kk, np.sqrt(kk**2 + 1), "k:", lw=0.8, label="ω = √(k²+1) (con massa 1)")
a.set_xlabel("k"); a.set_ylabel("ω"); a.set_xlim(0, 12); a.set_ylim(0, 12)
a.set_title("S4: spettro delle onde uscenti (guscio 4.5-6.3)\nesame pre-registrato FAIL (stimatore rovinato dal vuoto ruotato)", fontsize=10); a.legend(fontsize=7)
fig.suptitle("Nodo 3D, seconda parte: stencil del reticolo, Lorentz, Λ con la coda, annichilazione in onde", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "knot2_figure.png"), dpi=110)
print("ok")
