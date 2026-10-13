"""Figura e numeri dello studio BATTIMENTO_1S (docs/BATTIMENTO_1S.md).

Legge le serie scritte dai test di src/atomo/test_battimento.rs (BATT_DIR=<dir>):
  accoppiato_a600.csv, libero_1s_a600.csv, dopo_libero_a600.csv, esame_p.csv, e i livelli liv_*.bin, radiale.bin
Uso: python3 lab/battimento_figura.py <dir>  ->  lab/out/battimento_figure.png, lab/out/battimento_results.json
"""
import glob, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = sys.argv[1]
OUT = os.path.join(os.path.dirname(__file__), "out")
N = 48


def load(n):
    return np.genfromtxt(os.path.join(D, n), delimiter=",", names=True)


def spettro(t, c, lo, hi, dw):
    T = t[-1] - t[0]
    h = np.sin(np.pi * (t - t[0]) / T) ** 2
    w = np.arange(lo, hi, dw)
    return w, np.abs(np.exp(1j * np.outer(w, t)) @ (c * h)) / h.sum()


# --- free two-tick walk: the 4 valleys ---------------------------------------------------------------
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1, -1]).astype(complex)
I2 = np.eye(2); SIG = [sx, sy, sz]; TH = 0.5


def S(j, k):
    M = np.zeros((4, 4), complex)
    M[:2, :2] = np.cos(k) * I2 + 1j * np.sin(k) * SIG[j]
    M[2:, 2:] = np.cos(k) * I2 - 1j * np.sin(k) * SIG[j]
    return M


ROT = np.zeros((4, 4), complex)
ROT[:2, :2] = ROT[2:, 2:] = np.cos(TH) * I2
ROT[:2, 2:] = ROT[2:, :2] = -1j * np.sin(TH) * I2


def U2(k):
    return (ROT @ S(0, k[0]) @ S(1, k[1]) @ S(2, k[2])) @ (ROT @ S(2, k[2]) @ S(1, k[1]) @ S(0, k[0]))


def omega_min(k):
    return np.abs(np.angle(np.linalg.eigvals(U2(k)))).min()


rng = np.random.default_rng(1)
equiv = max(np.max(np.abs(np.sort(np.angle(np.linalg.eigvals(U2(d)))) - np.sort(np.angle(np.linalg.eigvals(U2(np.array(q) + d))))))
            for q in [(np.pi / 2, np.pi / 2, 0), (np.pi / 2, 0, np.pi / 2), (0, np.pi / 2, np.pi / 2)]
            for d in rng.normal(size=(100, 3)) * 0.3)
g = np.linspace(-np.pi, np.pi, 97)
mappa = np.array([[omega_min((a, b, 0.0)) for a in g] for b in g])

# --- valley content of the exact levels -----------------------------------------------------------------
k1 = np.fft.fftfreq(N) * 2 * np.pi
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
VALLI = [(0, 0, 0), (np.pi / 2, np.pi / 2, 0), (np.pi / 2, 0, np.pi / 2), (0, np.pi / 2, np.pi / 2)]


def carica(f):
    a = np.fromfile(f, dtype="<f8")
    return (a[0::2] + 1j * a[1::2]).reshape(N, N, N, 4)


def valli(v):
    P = np.sum(np.abs(np.fft.fftn(v, axes=(0, 1, 2))) ** 2, axis=3); P /= P.sum()
    fk = [np.mod(K + np.pi / 4, np.pi) - np.pi / 4 for K in (KX, KY, KZ)]
    d = [sum(np.minimum(np.abs(fk[a] - c[a]), np.pi - np.abs(fk[a] - c[a])) ** 2 for a in range(3)) for c in VALLI]
    lab = np.argmin(np.array(d), axis=0)
    return [float(P[lab == i].sum()) for i in range(4)]


rad = carica(os.path.join(D, "radiale.bin"))
livelli = []
for f in sorted(glob.glob(os.path.join(D, "liv_*.bin"))):
    v = carica(f)
    om = float(os.path.basename(f)[4:-4])
    livelli.append({"omega": om, "peso_radiale": float(abs(np.vdot(v, rad)) ** 2), "valli": valli(v)})
valli_rad = valli(rad)
valli_lab = valli(carica(os.path.join(D, "phi1s_a600.bin")))

# --- time series -------------------------------------------------------------------------------------------
acc = load("accoppiato_a600.csv"); lib = load("libero_1s_a600.csv"); es = load("esame_p.csv")
p_acc = 1 - (acc["cre"] ** 2 + acc["cim"] ** 2); p_lib = 1 - (lib["cre"] ** 2 + lib["cim"] ** 2)
ev_a = slice(0, None, 8); m = min(len(acc), len(lib))
diff = float(np.max(np.abs(p_acc[:m:2] - p_lib[:m:2])))
dopo = load("dopo_libero_a600.csv")
td = dopo["t"][::2]; c1 = (dopo["c1sre"] + 1j * dopo["c1sim"])[::2]; cd = (dopo["cdre"] + 1j * dopo["cdim"])[::2]
w, Sd = spettro(td, cd, 0.80, 1.02, 0.0002)
w2, Sp = spettro(td, 1 - np.abs(c1) ** 2 - np.mean(1 - np.abs(c1) ** 2), 0.002, 0.1, 0.0002)
riga = float(w2[np.argmax(Sp)])

fig, ax = plt.subplots(2, 2, figsize=(12, 9))
a = ax[0, 0]
a.plot(acc["t"][ev_a], p_acc[ev_a], color="#c0392b", lw=1.4, label="accoppiata, 1s ascoltato (prima)")
a.plot(lib["t"][ev_a], p_lib[ev_a], color="#7f8c8d", lw=1.0, ls="--", label="stesso 1s, SENZA luce")
a.plot(es["t"], es["p_exc"], color="#2471a3", lw=1.8, label="accoppiata, 1s esatto (cura)")
a.axvline(450, color="k", lw=0.6, ls=":"); a.text(455, 0.93, "luce uscita", fontsize=8)
a.set_xlim(0, 1600); a.set_ylim(0, 1); a.set_xlabel("t (λ_C/c)"); a.set_ylabel("p_exc = 1 − |⟨1s|ψ⟩|²")
a.set_title("a. il battimento c'è anche senza luce; con l'1s esatto sparisce"); a.legend(fontsize=8, loc="center right")

a = ax[0, 1]
a.semilogy(w, Sd, color="#2c3e50", lw=1)
for om, nome in [(0.84037, "0,8404"), (0.88627, "0,8863"), (0.8995, "1s 0,8993/0,8996"), (0.9786, "2p3/2")]:
    a.axvline(om, color="#e67e22", lw=0.6, ls="--"); a.text(om + 0.001, 0.4, nome, rotation=90, fontsize=8, va="top")
a.set_ylim(1e-4, 1); a.set_xlabel("ω (mc²/ħ)"); a.set_ylabel("peso")
a.set_title(f"b. righe dello stato dopo l'impulso (t = 500)\nbattimento p_exc: {riga:.4f} → periodo {2*np.pi/riga:.0f} λ_C/c")

a = ax[1, 0]
im = a.imshow(mappa, extent=[-np.pi, np.pi, -np.pi, np.pi], origin="lower", cmap="magma_r", vmin=1, vmax=2.2)
for q in [(0, 0), (np.pi / 2, np.pi / 2), (-np.pi / 2, np.pi / 2), (np.pi / 2, -np.pi / 2), (-np.pi / 2, -np.pi / 2), (np.pi, 0), (0, np.pi), (np.pi, np.pi)]:
    a.plot(*q, "o", mfc="none", mec="#1abc9c", ms=10)
a.set_xlabel("k_x (per cella)"); a.set_ylabel("k_y"); plt.colorbar(im, ax=a, label="ω più bassa del cammino libero a due tick")
a.set_title(f"c. le valli del cammino libero (k_z = 0): copie esatte a {equiv:.0e}")

a = ax[1, 1]
x = np.arange(len(livelli)); b0 = np.zeros(len(livelli))
col = ["#2471a3", "#e67e22", "#27ae60", "#8e44ad"]
nomi = ["k = 0", "(π/2,π/2,0)", "(π/2,0,π/2)", "(0,π/2,π/2)"]
for i in range(4):
    h = np.array([l["peso_radiale"] * l["valli"][i] for l in livelli])
    a.bar(x, h, bottom=b0, color=col[i], label=nomi[i]); b0 += h
a.set_xticks(x, [f"{l['omega']:.4f}" for l in livelli]); a.set_ylabel("peso nell'1s radiale liscio")
a.set_title("d. l'1s liscio (98 % valle k = 0) si divide su 4 livelli esatti,\nognuno 25/25/25/25 sulle 4 valli")
a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "battimento_figure.png"), dpi=110)

res = {"livelli_esatti_multipletto_1s": livelli, "valli_1s_radiale": valli_rad, "valli_1s_ascoltato_600": valli_lab,
       "equivalenza_valli_cammino_libero": float(equiv), "max_diff_p_exc_accoppiata_vs_libera": diff,
       "riga_battimento": riga, "periodo_battimento": float(2 * np.pi / riga),
       "p_exc_libera_1s_ascoltato_max": float(p_lib[::2].max())}
json.dump(res, open(os.path.join(OUT, "battimento_results.json"), "w"), indent=1, ensure_ascii=False)
print(json.dumps(res, indent=1, ensure_ascii=False))
