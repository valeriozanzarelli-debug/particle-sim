"""trap_plots.py - figure di docs/LUCE_INTRAPPOLATA.md dai risultati in lab/out/trap_*.json|npz.
uso: python lab/trap_plots.py   (dopo trap_null.py 2d|torus|hopf, trap_kerr.py, trap_derrick.py, trap_dirac.py)
"""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(__file__), "out")
J = lambda n: json.load(open(os.path.join(OUT, n)))

# ------------------------------------------------------------- 1) luce chiusa nel vuoto lineare
d2, dt, dh = J("trap_null_2d.json"), J("trap_null_torus.json"), J("trap_null_hopf.json")["hopfion"]
z = np.load(os.path.join(OUT, "trap_null_ring2d.npz"))
fig = plt.figure(figsize=(14, 7.5))
L = float(z["L"])
for i, (key, title) in enumerate([("e0", "t = 0"), ("e_mid", "t = 4 λ/c"), ("e_end", "t = 10 λ/c")]):
    ax = fig.add_subplot(2, 4, i + 1)
    img = z[key]
    n = img.shape[0]
    ext = L / 2
    ax.imshow(np.sqrt(img.T), origin="lower", extent=[-ext, ext, -ext, ext], cmap="magma")
    ax.set_xlim(-14, 14); ax.set_ylim(-14, 14)
    ax.set_title(f"anello R=4λ, w=λ: {title}", fontsize=9)
    ax.set_xlabel("x [λ]", fontsize=8)
ax = fig.add_subplot(2, 4, 4)
cur = d2["2d_ring_curve"]
ax.plot(cur["t"], cur["annulus"], label="nel tubo |r-R|<2w")
ax.plot(cur["t"], cur["disk"], label="nel disco r<2R")
ax.set_xlabel("t [λ/c]"); ax.set_ylabel("frazione di energia"); ax.legend(fontsize=8)
ax.set_title("2D, R=4λ: l'energia esce", fontsize=9)
ax = fig.add_subplot(2, 4, 5)
for rec in d2["2d_ring"]:
    ax.plot(rec["ray_pred_annulus"], rec["t_half_annulus"], "o", color="C0")
    ax.plot(rec["ray_pred_disk"], rec["t_half_disk"], "s", color="C1")
for rec in dt["3d_torus"]:
    ax.plot(rec["ray_pred_tube"], rec["t_half_tube"], "^", color="C2")
    ax.plot(rec["ray_pred_ball"], rec["t_half_ball"], "v", color="C3")
xx = np.linspace(0, 30, 2)
ax.plot(xx, xx, "k--", lw=0.8)
ax.set_xlabel("tempo di fuga dei raggi rettilinei [λ/c]")
ax.set_ylabel("t½ misurato [λ/c]")
ax.set_title("tempo di fuga = geometria dei raggi\n(o,s: 2D anello; ^,v: 3D toro)", fontsize=9)
ax = fig.add_subplot(2, 4, 6)
ax.plot(dh["t"], dh["f_r1"], label="r<1"); ax.plot(dh["t"], dh["f_r2"], label="r<2"); ax.plot(dh["t"], dh["f_r4"], label="r<4")
ax.set_xlabel("t [a/c]"); ax.set_ylabel("frazione di energia"); ax.legend(fontsize=8)
ax.set_title("hopfione (luce annodata): si apre", fontsize=9)
ax = fig.add_subplot(2, 4, 7)
H = np.array(dh["helicity_m"]) + np.array(dh["helicity_e"])
ax.plot(dh["t"], H / H[0], "C4", label="elicità (nodo) / iniziale")
ax.plot(dh["t"], np.array(dh["energy"]) / dh["energy"][0], "k--", label="energia / iniziale")
ax.plot(dh["t"], dh["rpeak"], "C5", label="raggio del guscio [a]")
ax.set_xlabel("t [a/c]"); ax.legend(fontsize=8)
ax.set_title("il nodo si conserva, ma vola via a c", fontsize=9)
ax = fig.add_subplot(2, 4, 8)
ax.axis("off")
bc = d2["ball_check"]
ax.text(0, 0.9, "controllo 'a palline' (somma di ondine\nritardate, formula di Poisson 2D)\n"
        f"vs propagatore spettrale, t = {bc['t']}:\nerrore max {bc['max_err_rel_to_max_u0']:.1e}\n"
        "(interpolazione bilineare del campo iniziale)", fontsize=9, va="top")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "trap_null.png"), dpi=90)
plt.close(fig)

# ------------------------------------------------------------- 2) regole non lineari / topologiche
dd = J("trap_derrick.json")
kerr_ring = [J(os.path.basename(p)) for p in sorted(glob.glob(os.path.join(OUT, "trap_kerr_ring*.json")))]
kb = J("trap_kerr_beam.json") if os.path.exists(os.path.join(OUT, "trap_kerr_beam.json")) else None
fig, axs = plt.subplots(2, 3, figsize=(14, 8))
ax = axs[0, 0]
if kb:
    z2 = np.arange(len(kb["nls_2d_kerr"]["0.95"]["peak"])) * 0.5
    for fac in ("0.95", "1.0", "1.05"):
        pk = kb["nls_2d_kerr"][fac]["peak"]
        ax.semilogy(np.arange(len(pk)) * 0.5, pk, label=f"Townes × {fac}")
    pk = kb["nls_2d_sat"]["peak"]
    ax.semilogy(np.arange(len(pk)) * kb["nls_2d_sat"]["dz"] / 10, pk, "k--", label="saturabile (+10%, ellisse) — z/10")
    ax.set_xlabel("distanza di propagazione z"); ax.set_ylabel("picco |ψ|"); ax.legend(fontsize=8)
    ax.set_title("fascio di Kerr 2+1: collassa o diffrange;\nsaturabile: stabile (ma corre a c lungo z)", fontsize=9)
ax = axs[0, 1]
for rec in [x["ring_fullwave"][0] for x in kerr_ring]:
    lab = f"{rec['kind']} p={rec['p']} A={rec['A']}"
    if rec.get("blowup_t") is not None:
        lab += f" (esplode t={rec['blowup_t']:.0f})"
    ax.plot(rec["t"], rec["frac"], label=lab)
ax.set_xlabel("t [λ/c]"); ax.set_ylabel("energia entro r<2R / iniziale"); ax.legend(fontsize=7)
ax.set_title("luce FERMA in un mezzo di Kerr (Maxwell 2D completo)", fontsize=9)
ax = axs[0, 2]
der = dd["derrick_3d"]
ax.plot(der["mu"], der["E"], "o-", label="scalare 3D senza topologia (max: instabile)")
sk = dd["skyrmion"]
ax.plot(sk["scale_mu"], np.array(sk["scale_E"]) / sk["M_coef"] * der["E"][3], "s-",
        label="skyrmione (riscalato; min: stabile)")
ax.set_xlabel("fattore di dilatazione μ"); ax.set_ylabel("energia")
ax.legend(fontsize=8); ax.set_title("teorema di Derrick", fontsize=9)
ax = axs[1, 0]
ax.plot(sk["profile_x"], sk["profile_F"], "C3")
ax.set_xlabel("x = e F_π r"); ax.set_ylabel("F(r)")
ax.set_title(f"skyrmione B = {sk['B']:.4f}: profilo a riccio", fontsize=9)
ax = axs[1, 1]
zb = np.load(os.path.join(OUT, "trap_baby.npz"))
q = zb["q"]
Lb = float(zb["L"])
ax.imshow(np.abs(q).T, origin="lower", extent=[-Lb / 2, Lb / 2, -Lb / 2, Lb / 2], cmap="viridis")
bb = dd["baby_skyrmion"]
ax.set_title(f"baby-skyrmione 2D: densità topologica\nQ = {bb['Q_end']:.3f}, E = {bb['E']:.3f} (radiale {bb['E_radial']:.3f})", fontsize=9)
ax = axs[1, 2]
h0 = bb["no_skyrme_history"]
ax.plot([h["it"] for h in h0], [h["size"] for h in h0], "o-", label="raggio (n3<0)")
ax.plot([h["it"] for h in h0], [abs(h["Q"]) for h in h0], "s-", label="|Q|")
ax.set_xlabel("iterazioni del rilassamento"); ax.legend(fontsize=8)
ax.set_title("senza termine di Skyrme: si restringe e 'cade'\ntra i nodi del reticolo (collasso di Derrick)", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "trap_solitons.png"), dpi=90)
plt.close(fig)

# ------------------------------------------------------------- 3) Dirac: Zitterbewegung e Bohm
dz = J("trap_dirac.json")
fig, axs = plt.subplots(1, 2, figsize=(12, 4))
for key, st in (("sigma5.0_misto", "-"), ("sigma1.0_misto", "-"), ("sigma5.0_solo_E_pos", "--")):
    r = dz["zb_1d"][key]
    axs[0].plot(r["t"], np.array(r["x"]) - r["x"][0], st, label=key.replace("sigma", "σ=").replace("_", " "))
axs[0].set_xlabel("t [ħ/mc²]"); axs[0].set_ylabel("<x>(t) - <x>(0)  [ħ/mc]"); axs[0].legend(fontsize=8)
axs[0].set_title("Zitterbewegung 1D: frequenza 2mc²/ħ, ampiezza ħ/2mc", fontsize=9)
sig = [dz["bohm_2d"][k]["sigma"] for k in dz["bohm_2d"]]
axs[1].plot(sig, [dz["bohm_2d"][k]["vphi_at_sigma"] for k in dz["bohm_2d"]], "o-", label="v_φ di Bohm-Dirac a r = σ")
axs[1].plot(sig, [dz["bohm_2d"][k]["max_speed_field"] for k in dz["bohm_2d"]], "s-", label="max |v| nel pacchetto")
axs[1].axhline(1, color="k", lw=0.8)
axs[1].set_xscale("log"); axs[1].set_xlabel("larghezza σ [ħ/mc]"); axs[1].set_ylabel("velocità [c]")
axs[1].legend(fontsize=8); axs[1].set_title("2D, energia positiva, spin su: la corrente gira, sempre ≤ c", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "trap_dirac.png"), dpi=90)
print("scritte lab/out/trap_null.png, trap_solitons.png, trap_dirac.png")
