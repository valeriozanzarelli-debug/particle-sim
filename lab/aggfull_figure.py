"""aggfull_figure.py - figure of AGGLOMERATORE COMPLETO -> lab/out/aggfull_figure.png"""
import json, os, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import aggfull_core as ac
from zz_core import dirac_formula

OUT = ac.OUT
R = json.load(open(os.path.join(OUT, "aggfull_atom.json")))
VD = json.load(open(os.path.join(OUT, "aggfull_verdict.json")))
GEO = json.load(open(os.path.join(OUT, "aggfull_geo_bz.json")))
ST = pickle.load(open(os.path.join(OUT, "aggfull_learn.pkl"), "rb"))
MC2, HEV = ac.MC2, 4.135667696e-15
fig, ax = plt.subplots(3, 3, figsize=(17, 14.5))
L1 = json.load(open(os.path.join(OUT, "aggfull_l1.json")))
L1V = json.load(open(os.path.join(OUT, "aggfull_l1_verdict.json")))

# (a) the spin-orbit part of the coupling: fine-walk data vs the learned stencil vs Dirac sigma.(k' x k)/4
a = ax[0, 0]
B = 256
h = B * ac.TH_A
rng = np.random.default_rng(5)
k = rng.normal(size=(300, 3)) * 0.02
kp = rng.normal(size=(300, 3)) * 0.02
Gd = ac.coarse_from_fine(k, kp, ac.TH_A, B)[3]
s = [x for x in ST if x["B"] == B][0]
Gl = ac.eval_coupling(s["C"], s["offs"], k * h, kp * h)
cr = np.cross(kp, k)
sod = np.einsum("nab,ba->n", Gd, ac.PAULI[2]).imag / 2
sol = np.einsum("nab,ba->n", Gl, ac.PAULI[2]).imag / 2
a.plot(cr[:, 2] / 4, sod, ".", ms=4, label="fine zig-zag walk (data)")
a.plot(cr[:, 2] / 4, sol, "x", ms=3, label="learned 7-site stencil")
xx = np.linspace(cr[:, 2].min(), cr[:, 2].max(), 2) / 4
a.plot(xx, xx, "k-", lw=0.8, label="Dirac: spin-orbit 1/4m$^2$")
a.set_xlabel(r"$(k'\times k)_z/4$  [$1/\lambda_C^2$]")
a.set_ylabel(r"Im tr$(\Gamma\sigma_z)/2$")
a.set_title("(a) how the coarse cell feels the field:\nspin-orbit comes out of the data (B = 256)")
a.legend(fontsize=8)

# (b) level errors: old hand-written restriction vs learned (a0/8)
a = ax[0, 1]
lev = VD.get("levels", {})
names, old, new = [], [], []
pairs = [("1s", "ca8_1s_N192", "e8_1s", 1), ("2s", "ca8_2s_N192", "e8_2s", 2), ("2p", None, "e8_2p", 2)]
for nm, o, n_, n in pairs:
    names.append(nm)
    old.append(abs(R[o]["rel_to_schroedinger"]) if o and o in R else np.nan)
    new.append(abs(R[n_]["rel_to_schroedinger"]) if n_ in R else np.nan)
oldA = {"1s": 0.00426, "2s": 0.01003, "2p": 0.00289}      # AGGLOMERATO_ZIGZAG A2 (old box +-12 a0, cell average)
x = np.arange(3)
a.bar(x - 0.27, [oldA[n] for n in names], 0.27, label="old study (cell average, box $\\pm$12 a$_0$)")
a.bar(x, old, 0.27, label="cell average, this run")
a.bar(x + 0.27, new, 0.27, label="learned coupling + Galerkin field (box $\\pm$16 a$_0$)")
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(names)
a.axhline(0.01, color="r", ls="--", lw=0.8); a.axhline(0.001, color="g", ls="--", lw=0.8)
a.set_ylabel("|E/E$_{Schr}$ - 1|  (cell a$_0$/8)")
a.set_title("(b) levels with cells a$_0$/8")
a.legend(fontsize=7)

# (c) the 2p sextet from the learned coupling vs Dirac
a = ax[0, 2]
if "q8_2p" in R:
    ev = np.array(R["q8_2p"]["dV_eigen_eV"])
    ev = (ev - ev.mean()) / HEV / 1e9
    a.plot(np.zeros(6) + np.arange(6) * 0.05, ev, "o", label="6x6 first order, learned rule")
    d12 = -(1 - dirac_formula(2, 1)) * MC2
    d32 = -(1 - dirac_formula(2, -2)) * MC2
    m = (2 * d12 + 4 * d32) / 6
    for d, nm in ((d12, "2p$_{1/2}$ Dirac"), (d32, "2p$_{3/2}$ Dirac")):
        a.axhline((d - m) / HEV / 1e9, color="k", lw=0.8)
        a.text(0.3, (d - m) / HEV / 1e9 + 0.2, nm)
    a.set_xlim(-0.2, 0.8)
    a.set_ylabel("GHz (from the centroid)")
a.set_title("(c) fine structure of 2p: 2 + 4 states")
a.set_xticks([])

# (d) recoil
a = ax[1, 0]
want = -1 / (1 + ac.MP)
lab, val = [], []
for e, p in (("e16_1s", "p16_1s"), ("e8_1s", "p8_1s"), ("e8_2s", "p8_2s"), ("e8_2p", "p8_2p")):
    if e in R and p in R:
        lab.append(e[1:].replace("_", " ").replace("16", "a0/16").replace("8 ", "a0/8 "))
        val.append((R[p]["E_eV"] / R[e]["E_eV"] - 1) / want)
a.bar(range(len(val)), val)
a.axhline(1, color="k", lw=0.8); a.axhspan(0.99, 1.01, color="g", alpha=0.15)
a.set_xticks(range(len(val))); a.set_xticklabels(lab, fontsize=8)
a.set_ylim(0.97, 1.03)
a.set_title("(d) recoil from the proton as a heavy ball:\nshift / reduced-mass prediction")

# (e) geometry of the Brillouin zone
a = ax[1, 1]
fr = [GEO[k]["mode_fraction"] for k in GEO]
er = [abs(GEO[k]["rel"]) for k in GEO]
a.plot([GEO[k]["mode_fraction"] for k in GEO if k != "cube"], [abs(GEO[k]["rel"]) for k in GEO if k != "cube"], "o-", label="ball of modes")
a.plot([1.0], [abs(GEO["cube"]["rel"])], "s", ms=9, label="full cube")
a.set_xlabel("fraction of the cube's modes kept"); a.set_ylabel("|1s error| (exact kinetic, a$_0$/8)")
a.set_yscale("log"); a.legend(fontsize=8)
a.set_title("(e) a rounder zone (bcc/fcc) would save only ~10 %")

# (f) level 1 for the helium study: mass, k^4, parity - old learner vs Taylor-pinned
a = ax[1, 2]
lab = ["m/m$_{fine}$-1\n(x1e5)", "k$^4$ axis", "k$^4$ face", "k$^4$ body", "odd part\n(kh=1, x10)"]
def row(k):
    v = L1[k]
    return [1e5 * v["m_over_mfine_minus_1"], v["k4frac"]["axis"], v["k4frac"]["face"], v["k4frac"]["body"], 10 * v["parity_odd_rel"]["1"]]
x = np.arange(5)
a.bar(x - 0.2, row("A_Z1_w2_U_s1"), 0.4, label="old learner (level 1 as used by helium)")
a.bar(x + 0.2, row("A_Z10_w4_T_s1"), 0.4, label="Taylor-pinned + parity (T)")
a.axhline(1, color="k", lw=0.6); a.axhline(0, color="k", lw=0.6)
a.set_xticks(x); a.set_xticklabels(lab, fontsize=8); a.legend(fontsize=8)
a.set_title("(f) level 1 (cell a$_0$/256): one cure, three symptoms")

# (g) walker kernel: negative weight vs beta
a = ax[2, 0]
for k, nm in (("A_Z1_w2_U_s1", "old, 125 taps"), ("A_Z10_w4_T_s1", "T, 729 taps"), ("A_Z2_w4_M_s1", "Markov M, 729 taps")):
    kn = L1[k]["kernel_neg"]
    b = np.array(sorted(map(float, kn))); v = np.array([max(kn[f"{x:g}"], 1e-16) for x in b])
    a.loglog(b, v, "o-", label=nm)
a.axvline(16 * ac.TH_A * ac.zz.ALPHA ** 2, color="gray", ls=":", lw=0.8); a.text(3e-5, 1e-3, "1 tick", fontsize=8, color="gray")
a.set_xlabel(r"$\beta$ (a.u.)"); a.set_ylabel(r"negative weight of K$_\beta$")
a.legend(fontsize=8); a.set_title("(g) the walkers' kernel: a probability at every $\\beta$ only for M")

# (h) the hop-length bound
a = ax[2, 1]
ws = [2, 3, 4]
dv = [max(L1[f"A_Z2_w{w}_M_s1"]["max_c2_rel_dev"], L1[f"A_Z2_w{w}_M_s1"]["max_c4_rel_dev"]) for w in ws]
a.semilogy(ws, np.maximum(dv, 1e-12), "s-", ms=9)
a.axvline(np.sqrt(3) / (16 * ac.TH_A), color="r", ls="--", label=r"bound $w \geq \sqrt{3}\,\lambda_C/h$ = 3.24")
a.set_xlabel("stencil half-width w (cells)"); a.set_ylabel("max error of mass / k$^4$ (positive stencil)")
a.legend(fontsize=8); a.set_title("(h) positive + relativistic needs jumps of ~1.7 $\\lambda_C$")

# (i) verdict table
a = ax[2, 2]
a.axis("off")
allv = VD["verdict"] + L1V["verdict"]
txt = "\n".join(f"{'PASS' if ok else 'FAIL'}  {nm[:58]}" for nm, ok, _ in allv)
a.text(0, 1, txt, va="top", family="monospace", fontsize=8)
a.set_title("(i) pre-registered exams (phases 2, amendments, 3)")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "aggfull_figure.png"), dpi=110)
print("saved")
