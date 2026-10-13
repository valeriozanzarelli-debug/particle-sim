"""tunnel_figure.py - figura dello studio: lab/out/tunnel_figura.png"""
import json, os, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
R1 = json.load(open(os.path.join(OUT, "tunnel_r1.json")))
R3 = json.load(open(os.path.join(OUT, "tunnel_r3.json")))
PK1 = json.load(open(os.path.join(OUT, "tunnel_r2_PK1_0.0125.json")))
PRE = json.load(open(os.path.join(OUT, "tunnel_prereg.json")))
E = 2.0
V0S = [1.2, 1.6, 2.0, 2.4, 2.8]
LS = [1, 2, 4, 8, 12, 16]
cols = plt.cm.viridis(np.linspace(0, 0.9, 5))

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))

# (a) pacchetti PK1
a_ = ax[0, 0]
s = PK1["snap"]
x = np.array(s["x"])
P = PRE["pacchetti"]["PK1"]
a_.semilogy(x, np.maximum(s["rho0"], 1e-30), color="0.6", lw=1, label="pacchetto in arrivo (t = 0)")
a_.semilogy(x, np.maximum(s["rhoF"], 1e-30), color="C0", lw=1.2, label="stesso pacchetto senza barriera (t finale)")
a_.semilogy(x, np.maximum(s["rhoT"], 1e-30), color="C3", lw=1.2, label="con la barriera (t finale)")
pk = PK1["peak"]
a_.axvline(pk["x0"] + pk["t"], color="k", ls="--", lw=1, label="raggio di luce dal picco iniziale")
a_.axvspan(0, P["L"], color="C1", alpha=0.3, label="barriera (V0 = E = 5 mc²)")
a_.set_xlim(-320, 260)
a_.set_ylim(1e-24, 1e-1)
a_.set_xlabel("x  [ħ/mc]")
a_.set_ylabel("densità per cella")
a_.set_title("PK1: il picco trasmesso esce 4,05 davanti alla luce\n(ma il fronte, a destra, non supera mai c)", fontsize=10)
a_.legend(fontsize=7, loc="lower left")

# (b) Hartman
b = ax[0, 1]
for c, V0 in zip(cols, V0S):
    tl = [R1["V0"][str(V0)][str(L)]["lattice"]["tphi"][3] for L in LS]
    tc = [R1["V0"][str(V0)][str(L)]["continuum"]["tphi"] for L in LS]
    b.plot(LS, tc, "-", color=c, lw=1)
    b.plot(LS, tl, "o", color=c, ms=4, label=f"V0 = {V0} (ε = {E - V0:+.1f})")
b.plot([0, 16], [0, 16], "k--", lw=1, label="L/c (luce)")
b.set_ylim(0, 4)
b.set_xlabel("spessore L  [ħ/mc]")
b.set_ylabel("tempo di fase τ_φ  [ħ/mc²]")
b.set_title("Hartman nella regola: τ_φ satura (punti: regola a = 0,0125; linee: Dirac)", fontsize=10)
b.legend(fontsize=7)

# (c) orologio delle rotazioni: rapporto e pendenze (barriera spessa)
c_ = ax[0, 2]
eps = np.array([E - V for V in V0S])
rat = [R1["verdict"]["E11b_rapporto"]["ratio"][str(V)] for V in V0S]
smz = [R1["verdict"]["E11c_pendenze"]["slopes"][str(V)]["slope_mz"] for V in V0S]
sz = [R1["verdict"]["E11c_pendenze"]["slopes"][str(V)]["slope_z"] for V in V0S]
ee = np.linspace(-0.95, 0.95, 200)
c_.plot(ee, ee, "C0-", lw=1)
c_.plot(eps, rat, "C0o", label="τ_my/τ_y (rotazioni / Larmor, nel piano)  → ε/m")
c_.plot(ee, 1 / np.sqrt(1 - ee ** 2), "C3-", lw=1)
c_.plot(eps, smz, "C3s", label="dτ_mz/dL (rotazioni, inclinazione)  → m/κ ≥ 1")
c_.plot(ee, ee / np.sqrt(1 - ee ** 2), "C2-", lw=1)
c_.plot(eps, sz, "C2^", label="dτ_z/dL (Büttiker–Landauer)  → ε/κ")
c_.axhline(1, color="k", ls="--", lw=1)
c_.text(-0.9, 1.05, "velocità della luce", fontsize=7)
c_.set_ylim(-2, 2.5)
c_.set_xlabel("ε = (E − V0)/mc²  (posizione nel gap)")
c_.set_title("Orologio che conta le rotazioni (barriera spessa)\npunti: regola a → 0; linee: formule", fontsize=10)
c_.legend(fontsize=7, loc="lower right")

# (d) limite NR, forma Ramos
d = ax[1, 0]
xs = [0.3, 0.5, 0.7, 0.9]
for kw, cc in (("1.0", "C0"), ("2.0", "C3")):
    yr = [R3["E10"]["data"][f"{kw}_{x}"]["rule_V0_001"][0] for x in xs]
    ys = [R3["E10"]["data"][f"{kw}_{x}"]["schrod"][0] for x in xs]
    d.plot(xs, ys, "-", color=cc, lw=1)
    d.plot(xs, yr, "o", color=cc, label=f"barriera gaussiana κ0 w = {kw}")
for xb in ("0.2", "0.4", "0.6", "0.8"):
    e9 = R3["E9"]["data"][xb]
    d.plot(float(xb), e9["NR_limit"][0], "k^", ms=5)
    d.plot(float(xb), e9["buttiker"][0], "k_", ms=12)
d.plot([], [], "k^", label="rettangolo κ0 L = 3 (regola, V0/m → 0) e Büttiker (—)")
d.set_xlabel("E / V0")
d.set_ylabel("τ_y · V0/ħ  (Larmor nel piano)")
d.set_title("Limite non relativistico: regola = Büttiker / Schrödinger\n(come Ramos 2020: energia più bassa → tempo più corto)", fontsize=10)
d.legend(fontsize=7)

# (e) surfisti (Bohm) PK1
e_ = ax[1, 1]
H = np.array(PK1["bohm"]["Xq_hist"], float)
t = H[:, 0]
for j in range(1, H.shape[1]):
    e_.plot(H[:, j], t, color="C3", lw=0.8)
e_.axvspan(0, 10, color="C1", alpha=0.3)
e_.plot([-10 + 0, 260], [0, 270], "k--", lw=1, label="luce dal bordo del pacchetto iniziale")
e_.legend(fontsize=7, loc="upper left")
e_.set_xlim(-30, 120)
e_.set_xlabel("x  [ħ/mc]")
e_.set_ylabel("t  [ħ/mc²]")
e_.set_title("Surfisti trasmessi (Bohm = quantili 1-D), PK1\npartono davanti, aspettano ~115 nella barriera, mai oltre c", fontsize=10)

# (f) tabella dei tempi PK1/PK2
f = ax[1, 2]
f.axis("off")
V = json.load(open(os.path.join(OUT, "tunnel_r2_verdict.json")))
PK2 = json.load(open(os.path.join(OUT, "tunnel_r2_PK2_0.0125.json")))
rows = [["", "PK1 (E=V0=5, L=10)", "PK2 (E=V0=2, L=8)"],
        ["L/c (luce)", "10", "8"],
        ["fase τ_φ (Wigner)", f"{PK1['stationary']['tphi_k0']:.3f}", f"{PK2['stationary']['tphi_k0']:.3f}"],
        ["Larmor τ_y (= dwell)", f"{PK1['tau_y']:.3f}", f"{PK2['tau_y']:.3f}"],
        ["Larmor τ_z", f"{PK1['tau_z']:.3f}", f"{PK2['tau_z']:.3f}"],
        ["rotazioni τ_my", f"{PK1['tau_my']:.4f}", f"{PK2['tau_my']:.4f}"],
        ["rotazioni τ_mz", f"{PK1['tau_mz']:.3f}", f"{PK2['tau_mz']:.3f}"],
        ["Bohm (surfisti) trasmessi", f"{PK1['bohm']['tauB_T_traj']:.1f}", f"{PK2['bohm']['tauB_T_traj']:.1f}"],
        ["T (trasmissione)", f"{PK1['PT']:.2e}", f"{PK2['PT']:.2e}"]]
tb = f.table(cellText=rows, loc="center", cellLoc="center")
tb.auto_set_font_size(False)
tb.set_fontsize(9)
tb.scale(1, 1.6)
f.set_title("Un tunnel, sette orologi (regola, a = 0,0125; unità ħ/mc²)", fontsize=10)

fig.suptitle("Tempo di tunnel sulla regola zig-zag (1+1 D): corri a c, poi ruota", fontsize=13)
fig.tight_layout()
fn = os.path.join(OUT, "tunnel_figura.png")
fig.savefig(fn, dpi=110)
os.makedirs("/mnt/project-files/problemi-aperti", exist_ok=True)
shutil.copy(fn, "/mnt/project-files/problemi-aperti/tunnel_figura.png")
print("ok", fn)
