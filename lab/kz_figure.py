"""kz_figure.py - lab/out/kz_figure.png (K1 alpha, sensitivity, K2 nu, K3 defect motion)."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out"); sys.path.insert(0, HERE)
import kz_k2

k1 = json.load(open(os.path.join(OUT, "kz_k1.json")))
k2 = json.load(open(os.path.join(OUT, "kz_k2.json")))
k3 = json.load(open(os.path.join(OUT, "kz_k3.json")))
tr = np.load(os.path.join(OUT, "kz_k3_tracks.npz"))

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
# 1. density vs tauQ
a = ax[0, 0]
for key, lab, mk in (("main_L48_a2_b0.5", "L=48 (principale)", "o"), ("L32", "L=32", "s"),
                     ("L64_window_tau>=256", "L=64", "^"), ("meanfield_protocol_L48", "conteggio a t~tauQ^0.5", "x")):
    if key in k1:
        d = k1[key]
        a.errorbar(d["tauQ"], d["mean_density"], yerr=d["sem_density"], fmt=mk, ms=6, capsize=2,
                   label="%s: alpha=%.3f" % (lab, d["alpha"]), mfc="none" if mk != "o" else None)
t = np.array([64, 4096.])
m = k1["main_L48_a2_b0.5"]; n0 = m["mean_density"][0]
a.plot(t, n0 * (t / 64) ** -0.5, "k--", lw=1, label="campo medio 0.5")
a.plot(t, n0 * (t / 64) ** -0.5733, "r-", lw=1, label="3D XY modello A 0.573")
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("tempo di raffreddamento tauQ (sweep)")
a.set_ylabel("lunghezza di vortice per sito"); a.set_title("K1: difetti congelati vs velocita di quench")
a.legend(fontsize=7.5)
# 2. alpha sensitivity
a = ax[0, 1]
items = [(k.replace("_", " "), v[0], v[1]) for k, v in k1["sensitivity_a_b_L48"].items()]
for key, lab in (("window_tau>=256_L48", "L48 tau>=256"), ("L32", "L32"), ("L64_window_tau>=256", "L64 tau>=256"),
                 ("meanfield_protocol_L48", "t~tauQ^0.5")):
    if key in k1:
        items.append((lab, k1[key]["alpha"], k1[key]["ci95"]))
y = np.arange(len(items))
for i, (lab, al, ci) in enumerate(items):
    a.errorbar(al, i, xerr=[[al - ci[0]], [ci[1] - al]], fmt="o", color="C3" if lab == "a2 b0.5" else "C0", capsize=3)
a.set_yticks(y); a.set_yticklabels([it[0] for it in items], fontsize=8)
a.axvline(0.5, color="k", ls="--", lw=1, label="campo medio 0.5")
a.axvline(0.5733, color="r", lw=1, label="3D XY 0.573")
a.axvline(0.668, color="g", ls=":", lw=1, label="elio (modello F) 0.668")
a.set_xlabel("alpha (IC 95% bootstrap)"); a.set_title("K1: robustezza (a = quando conto, b = quanto liscio)")
a.legend(fontsize=7.5)
# 3. K2
a = ax[0, 2]
Ls = k2["L"]
for name, idx, mk in (("d ln<|m|>/d beta", "dlnm_dbeta", "o"), ("d ln<m^2>/d beta", "dlnm2_dbeta", "s")):
    v = np.array([k2["L%d" % L][idx][0] for L in Ls]); e = np.array([k2["L%d" % L][idx][1] for L in Ls])
    a.errorbar(Ls, v, yerr=e, fmt=mk, mfc="none", label=name)
    p = np.polyfit(np.log(Ls), np.log(v), 1)
    a.plot(Ls, np.exp(np.polyval(p, np.log(Ls))), "-", lw=1, label="pendenza 1/nu -> nu=%.3f" % (1 / p[0]))
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("L"); a.set_ylabel("derivata a beta_c")
a.set_title("K2: nu = %.3f +- %.3f (atteso 0.672)" % tuple(k2["nu"]))
a.legend(fontsize=7.5)
ins = a.inset_axes([0.58, 0.08, 0.38, 0.33])
Ts = np.linspace(2.17, 2.235, 60)
for L in Ls:
    d = np.load(os.path.join(OUT, "kz_k2_L%d.npz" % L))
    ins.plot(Ts, [kz_k2.obs(d["E"], d["M"], beta=1 / T)[0] for T in Ts], lw=1, label=str(L))
ins.axvline(2.20184, color="k", lw=0.6); ins.set_title("Binder U4(T)", fontsize=7); ins.tick_params(labelsize=6)
# 4. K3a pair
a = ax[1, 0]
ta = tr["a_t"]; sep = tr["a_xm"] - tr["a_xp"]
a.plot(ta, sep, "C0", label="d=24 (pre-registrata): bloccata dal reticolo")
a.plot(tr["a8_t"], tr["a8_xm"] - tr["a8_xp"], "C2", lw=2, label="d=8 (a posteriori): si annichila")
a.set_xlim(0, 120)
a.set_xlabel("tempo"); a.set_ylabel("distanza (siti)")
a2 = a.twinx(); a2.plot(tr["a8_t"], np.array(tr["a8_Ekin"]) / np.array(tr["a8_E"]), "C1", lw=1, label="d=8: E_cin/E_tot")
a2.set_ylabel("d=8: E_cin / E_tot", color="C1")
ka = k3["a_pair_d8_posthoc"]
a.set_title("K3a: coppia vortice-antivortice (dinamica inerziale)")
a.legend(loc="lower left", fontsize=8)
# 5. K3a wave front
a = ax[1, 1]
sh = tr["a8_shells"]
a.plot(sh[:, 0], sh[:, 2], "C3", label="raggio 90%% energia: v=%.2f" % ka["front_speed_r90"])
a.plot(sh[:, 0], sh[:, 1], "C2", label="raggio 50%% energia: v=%.2f" % ka["front_speed_r50"])
a.plot(sh[:, 0], sh[:, 0] + sh[0, 2], "k:", lw=1, label="velocita c = 1")
a.set_xlabel("tempo dopo l'annichilazione"); a.set_ylabel("raggio (siti)")
a.set_title("K3a (d=8): dopo l'annichilazione l'energia parte come anello"); a.legend(fontsize=8)
# 6. K3c/d
a = ax[1, 2]
c = tr["c"]
a.plot(c[:, 0], c[:, 1] / c[0, 1], "C0o-", ms=3, label="2D: vortici (n0=%d)" % c[0, 1])
d = tr["d"]
if d.size:
    a.plot(d[:, 0], d[:, 1] / d[0, 1], "C3s-", ms=3, label="3D: lunghezza linee")
a.set_xlabel("tempo dinamica inerziale (gamma=0.02)"); a.set_ylabel("frazione rimasta")
a.set_title("K3c/d: i difetti nati si muovono e si annichilano"); a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "kz_figure.png"), dpi=110)
print("saved")
