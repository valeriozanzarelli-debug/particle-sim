"""tavg_figure.py - figura dello studio tavg_ (ossigeno): lab/out/tavg_figure.png"""
import os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
EV = 27.211386245988
c = json.load(open(os.path.join(OUT, "tavg_corr.json")))
g = json.load(open(os.path.join(OUT, "tavg_grid.json")))
s = json.load(open(os.path.join(OUT, "tavg_simm.json")))
m = c["models"][c["scelto"]]

fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
# (a) correlazione: allenamento e cieco
names = [t["name"] for t in c["training"]]
ec = [-1e3 * t["Ec"] for t in c["training"]]
fitp = [-1e3 * (l["Ec"] + l["err_mEh"] / 1e3) for l in m["loo"]]
x = range(len(names))
ax[0].plot(x, ec, "o", color="#2b6cb0", label="esatto della regola (allenamento)")
ax[0].plot(x, fitp, "x", color="#c05621", label="previsto lascia-uno-fuori")
keys = {"N": 188.8, "O": 257.9}
for i, k in enumerate(("N", "O")):
    xi = len(names) + 1 + i
    ax[0].plot(xi, -1e3 * c["cieco"][k]["Ec_pred"], "s", color="#c05621", ms=8)
    ax[0].plot(xi, keys[k], "*", color="k", ms=12)
ax[0].plot([], [], "s", color="#c05621", label="previsto alla cieca")
ax[0].plot([], [], "*", color="k", label="valore noto (solo esame)")
ax[0].set_xticks(list(x) + [len(names) + 1, len(names) + 2])
ax[0].set_xticklabels(names + ["N", "O"], rotation=70, fontsize=7)
ax[0].set_ylabel("-E_c (mEh)")
ax[0].set_title("a) correlazione imparata: regge sulle coppie 1s,\nnon sulle coppie 2s/2p dell'ossigeno")
ax[0].legend(fontsize=7)

# (b) IP
ip_hf = {"O": (g["atom_Z8_q1"]["E"] - g["atom_Z8_q0"]["E"]) * EV,
         "N": (g["atom_Z7_q1"]["E"] - g["atom_Z7_q0"]["E"]) * EV}
dc = {k: (c["cieco"][k + "+"]["Ec_pred"] - c["cieco"][k]["Ec_pred"]) * EV for k in ("N", "O")}
exp = {"N": 14.534, "O": 13.618}
for i, k in enumerate(("N", "O")):
    ax[1].bar(i - 0.25, ip_hf[k], 0.25, color="#a0aec0", label="HF griglia" if i == 0 else None)
    ax[1].bar(i, ip_hf[k] + dc[k], 0.25, color="#c05621", label="HF + correlazione imparata" if i == 0 else None)
    ax[1].bar(i + 0.25, exp[k], 0.25, color="k", label="sperimentale (esame)" if i == 0 else None)
    ax[1].fill_between([i - 0.4, i + 0.4], exp[k] * 0.97, exp[k] * 1.03, color="#38a169", alpha=0.2)
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["N", "O"])
ax[1].set_ylim(11, 15.2); ax[1].set_ylabel("energia di ionizzazione (eV)")
ax[1].set_title("b) IP: banda verde = 3 %. N passa, O no (-3,0 %)")
ax[1].legend(fontsize=7)

# (c) rottura di simmetria
lab, val = [], []
for key, lbl in (("Z8_h0.2500_L16_fermo", "a0/4 fermo"), ("Z8_h0.2500_L16_rilassato", "a0/4 rilassato"),
                 ("Z8_h0.1667_L12_fermo", "a0/6 fermo"), ("Z8_h0.1667_L12_rilassato", "a0/6 rilassato")):
    lab.append(lbl); val.append(s[key]["dE_vs_radiale_mEh"])
ax[2].bar(range(4), val, color=["#a0aec0", "#2b6cb0"] * 2)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].set_xticks(range(4)); ax[2].set_xticklabels(lab, fontsize=8)
ax[2].set_ylabel("E griglia - E radiale, stesso blocco (mEh)")
ax[2].set_title("c) ossigeno: orbitali p uguali -> errore piccolo;\nliberi -> -3,1 mEh a ogni cella (fisica HF)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "tavg_figure.png"), dpi=110)
print("ok")
