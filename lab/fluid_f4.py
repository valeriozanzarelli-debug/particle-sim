"""fluid_f4.py - F4 'rendering in negative' + the summary figure out/fluid_figure.png.
Draws only the defects (density-deficit isoline, vortex core isoline, kink markers), counts primitives vs balls.
Also evaluates the F1/F2/F3 exams from the run outputs into out/fluid_verdict.json.
POST-HOC overlay (declared, not pre-registered): parameter-free inviscid COMPRESSIBLE 2D cylinder mass
M(s)/(rho pi R^2) = -K1(q R)/(q R K1'(q R)), q = s/c, drawn on panel (a)."""
import sys, os, json
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.special import k0, k1
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
J = lambda n: json.load(open(os.path.join(OUT, n)))

# ---------------------------------------------------------------- F4 on the F1 snapshot (R_w = 10)
sn = np.load(os.path.join(OUT, "fluid_f1_snapshot.npz"))
pos, X, L, Rw = sn["pos"], sn["X"], float(sn["L"]), float(sn["Rw"])
r10 = J("fluid_f1_R10.json"); rho_b = r10["rho_b"]
g = 0.5; nb = int(round(L / g))
H, xe, ye = np.histogram2d(pos[:, 0] % L, pos[:, 1] % L, bins=nb, range=[[0, L], [0, L]])
rho_cg = gaussian_filter(H / g ** 2, 1.5 / g, mode="wrap")
defi = rho_b - rho_cg
xc = 0.5 * (xe[1:] + xe[:-1])
XX, YY = np.meshgrid(xc, xc, indexing="ij")
fig0 = plt.figure(); cs = plt.contour(XX, YY, rho_cg, levels=[rho_b / 2]); plt.close(fig0)
segs = cs.allsegs[0]
n_prim = int(sum(max(len(s) - 1, 1) for s in segs))
n_loops = len(segs)
dx = (XX - X[0] + L / 2) % L - L / 2; dy = (YY - X[1] + L / 2) % L - L / 2
mask = defi > rho_b / 2
cx = float((dx[mask] * defi[mask]).sum() / defi[mask].sum()); cy = float((dy[mask] * defi[mask]).sum() / defi[mask].sum())
off = float(np.hypot(cx, cy))
win = np.hypot(dx - cx, dy - cy) < Rw + 6
N_field = float(defi[win].sum() * g ** 2)
f4 = dict(balls=int(len(pos)), primitives=n_prim, loops=n_loops, compression=n_prim / len(pos),
          centroid_offset=off, N_miss_field=N_field, N_miss_profile=r10["N_miss"],
          N_rel=N_field / r10["N_miss"] - 1,
          exams=dict(E4a=bool(n_prim / len(pos) < 0.01), E4b=bool(off < 1.0 and abs(N_field / r10["N_miss"] - 1) < 0.2)))
# F2 vortex: |psi|^2 on a 2D grid (what a full GP renderer would draw) vs one isoline
f2n = np.load(os.path.join(OUT, "fluid_f2.npz")); rr, ff = f2n["r"], f2n["f"]
gx = np.arange(-25, 25, 0.25); GX, GY = np.meshgrid(gx, gx, indexing="ij"); RR = np.hypot(GX, GY)
dens2 = np.interp(RR, rr, ff) ** 2
fig0 = plt.figure(); c2 = plt.contour(GX, GY, dens2, levels=[0.5]); plt.close(fig0)
v_prim = int(sum(len(s) - 1 for s in c2.allsegs[0]))
f4["vortex"] = dict(grid_cells=int(dens2.size), atoms_in_box=float(50 * 50), primitives=v_prim)
f4["kink"] = dict(balls=1600, primitives=1)

# ---------------------------------------------------------------- verdicts
F1 = {t: J("fluid_f1_%s.json" % t) for t in ["R3", "R6", "R10", "R6b"]}
main = [F1["R3"], F1["R6"], F1["R10"]]
Nm = np.array([d["N_miss"] for d in main]); Ms = np.array([d["M_star"] for d in main])
slope = float(np.polyfit(np.log(Nm), np.log(Ms), 1)[0])
v1 = dict(
    E1=bool(all(d["M_hi"] <= d["N_miss"] <= d["M_lo"] for d in main)),
    E2=bool(all(0.75 <= d["ratio_star"] <= 1.25 for d in main)), ratios=[d["ratio_star"] for d in main],
    E3=bool(0.85 <= slope <= 1.15), slope=slope,
    E4=bool(abs(F1["R6b"]["M_star"] - F1["R6"]["M_star"]) / F1["R6"]["M_star"] < 0.10),
    E4_rel=abs(F1["R6b"]["M_star"] - F1["R6"]["M_star"]) / F1["R6"]["M_star"])
F2 = J("fluid_f2.json"); F3 = J("fluid_f3.json")
q = 1.0
v1["posthoc_inviscid_compressible_at_sR_over_c_1"] = float(k1(q) / (q * (k0(q) + k1(q) / q)))
verdict = dict(F1=v1, F2=F2["exams"], F3=F3["exams"], F4=f4["exams"], f4=f4)
json.dump(verdict, open(os.path.join(OUT, "fluid_verdict.json"), "w"), indent=1)
print(json.dumps(verdict, indent=1))

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(2, 4, figsize=(21, 9.5))
cols = {"R3": "#1f77b4", "R6": "#ff7f0e", "R10": "#2ca02c", "R6b": "#9467bd"}
a = ax[0, 0]
for t, d in F1.items():
    z = np.load(os.path.join(OUT, "fluid_f1_%s.npz" % t))
    x = z["sgrid"] * d["R_eff"] / d["c_S"]
    lab = "R_w=%g (N_miss=%.0f)" % (d["Rw"], d["N_miss"]) + (", M0/4" if t == "R6b" else "")
    a.plot(x, z["Ms"] / d["N_miss"], color=cols[t], ls="--" if t == "R6b" else "-", label=lab)
xs = np.geomspace(0.05, 10, 200)
a.plot(xs, k1(xs) / (xs * (k0(xs) + k1(xs) / xs)), "k:", label="fluido ideale comprimibile (post-hoc, 0 parametri)")
a.axhline(1, color="gray", lw=0.8); a.axvspan(0.5, 2, color="0.92")
a.axvline(1, color="gray", lw=0.8, ls="--")
a.set_xscale("log"); a.set_ylim(0, 3); a.set_xlim(0.05, 10)
a.set_xlabel("s R_eff / c_s  (frequenza di Laplace adimensionale)"); a.set_ylabel("massa del buco / palline mancanti")
a.set_title("F1  inerzia del buco ~ palline mancanti?"); a.legend(fontsize=7)
a = ax[0, 1]
a.loglog(Nm, Ms, "o-", label="M_fl(s*), slope %.2f" % slope)
a.loglog(Nm, Nm, "k--", label="M = m N_miss")
a.set_xlabel("palline mancanti N_miss"); a.set_ylabel("massa idrodinamica (s* = c/R)"); a.legend(fontsize=8)
a.set_title("F1  scala con il numero di palline mancanti")
a = ax[0, 2]
R = np.array(F2["R"])
a.semilogx(R, F2["E"], "o-", label="energia vortice E")
a.semilogx(R, F2["E_swirl"], "s-", label="energia del vortice che gira")
a.semilogx(R, F2["N_miss"], "^-", label="palline mancanti x c^2")
a.semilogx(R, np.pi * np.log(F2["A"] * R), "k:", label="pi ln(%.3f R/xi)" % F2["A"])
a.set_xlabel("raggio del sistema R / xi"); a.set_title("F2  vortice GP: E ~ ln R, = (mancanti) c^2 + cost")
a.legend(fontsize=8)
a = ax[0, 3]
for h, mk in [("h0.25", "o"), ("h0.5", "s"), ("h1", "^")]:
    rows = F3["single"][h]
    v = [r["v"] for r in rows]
    a.plot(v, [1 + r["wg_rel"] for r in rows], mk + "-", label="larghezza*gamma, " + h)
    a.plot(v, [1 + r["Eg_rel"] for r in rows], mk + "--", label="E/(gamma E0), " + h)
a.axhspan(0.95, 1.05, color="0.9"); a.set_ylim(0.85, 1.3)
a.set_xlabel("velocita del kink / c_s"); a.set_title("F3  contrazione di Lorentz e E = gamma m c^2")
a.legend(fontsize=7)
a = ax[1, 0]
Fs = [d["F"] for d in F3["driven"]]; post = J("fluid_f3_post.json")
FF = np.geomspace(0.03, 1, 100)
a.semilogx(FF, 1 / np.sqrt(1 + (4 * 0.05 / (np.pi * FF)) ** 2), "k-", label="McLaughlin-Scott (continuo)")
a.semilogx(Fs, [d["v"] for d in F3["driven"]], "o", label="catena h=0.25, T=150 (pre-reg)")
a.semilogx(Fs, [d["v"] for d in post], "x", ms=9, label="T=600 (post-hoc)")
a.axhline(1, color="r", ls="--", label="c_s")
a.set_xlabel("forza F"); a.set_ylabel("velocita terminale"); a.set_title("F3  la velocita limite e' c_s")
a.legend(fontsize=8)
a = ax[1, 1]
z = np.load(os.path.join(OUT, "fluid_f3.npz"))
for k in ["h1_v0.05", "h1_v0.2", "h1_v0.4", "h0.25_v0.2", "h1_v0.8"]:
    hh = z["hist_" + k]; a.plot(hh[:, 0], hh[:, 1], label=k)
a.set_xlabel("tempo"); a.set_ylabel("energia in |x|<10"); a.set_title("F3  kink + antikink: niente annichilazione")
a.legend(fontsize=8)
a = ax[1, 2]
d = ((pos - X + L / 2) % L) - L / 2
a.scatter(d[:, 0], d[:, 1], s=1, c="0.75")
a.contour(dx, dy, rho_cg, levels=[rho_b / 2], colors="r")
a.set_xlim(-25, 25); a.set_ylim(-25, 25); a.set_aspect("equal")
a.set_title("F4  in negativo: %d palline (grigio) -> %d segmenti (rosso)" % (len(pos), n_prim), fontsize=10)
a = ax[1, 3]
names = ["F1 buco", "F2 vortice", "F3 kink"]
nb_ = [len(pos), dens2.size, 1600]; npr = [n_prim, v_prim, 1]
xx = np.arange(3)
a.bar(xx - 0.2, nb_, 0.4, label="palline / celle"); a.bar(xx + 0.2, npr, 0.4, label="primitive disegnate")
a.set_yscale("log"); a.set_xticks(xx); a.set_xticklabels(names); a.legend(fontsize=8)
for i in range(3): a.text(xx[i], max(nb_[i], 1) * 1.3, "x%.0f" % (nb_[i] / npr[i]), ha="center")
a.set_title("F4  compressione")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fluid_figure.png"), dpi=110)
print("figure written")
