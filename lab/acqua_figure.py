"""acqua_figure.py - densita' elettronica di H2O (RHF X=3 aug al minimo CCSD(T) CBS) come palline luminose in una
fetta sottile nel piano della molecola; numeri dei verdetti sopra. Anche l'export 'approssimato' (lab/out/acqua_dati/).
Uso: python acqua_figure.py [export]"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_int as I, acqua_qc as Q, acqua_basis as AB, acqua_run as AR, acqua_plan as PL


def ao_values(B, P):
    """values of all spherical AOs at points P (n,3)"""
    out = np.zeros((len(P), B.M))
    for k in range(len(B.sL)):
        L, a, C, N = B.sL[k], B.sa[k], B.sC[k], B.sN[k]
        d = P - C; g = np.exp(-a * np.sum(d * d, 1)) * N
        mono = np.stack([d[:, 0] ** x * d[:, 1] ** y * d[:, 2] ** z for x, y, z in I.cart_list(L)], 1)
        out[:, B.off[k]:B.off[k + 1]] = (mono @ I.SPH[L, :2 * L + 1, :mono.shape[1]].T) * g[:, None]
    return out


def density():
    s = PL.st_load(); r, th = s["min_cbs"]
    cen, Z, sh = AB.shells_for(["O", "H", "H"], AR.geom(r, r, th), 3, aug=True)
    B = I.Basis(cen, Z, sh); ints = I.integrals(B); res = Q.scf(B, ints, 5)
    C = res["C"][:, :5]
    return B, C, res, r, th


def sample(B, C, n, rng, slab=None):
    """electron positions drawn from rho (rejection on a box); slab: |y| < slab"""
    pts = []; box = 4.5; rho_max = None
    while sum(len(p) for p in pts) < n:
        P = rng.uniform(-box, box, size=(200000, 3)); P[:, 2] += 1.0
        if slab:
            P[:, 1] = rng.uniform(-slab, slab, len(P))
        phi = ao_values(B, P) @ C; rho = 2 * np.sum(phi ** 2, 1)
        w = np.minimum(rho, 2.0)                       # cap the 1s peak so the valence shows (drawing only)
        if rho_max is None:
            rho_max = 2.0
        keep = rng.random(len(P)) < w / rho_max
        pts.append(np.c_[P[keep], rho[keep]])
    return np.vstack(pts)[:n]


def figure(path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    B, C, res, r, th = density(); rng = np.random.default_rng(1)
    S = sample(B, C, 60000, rng, slab=0.35)
    v = json.load(open(os.path.join(AB.OUT, "acqua_verdetti.json")))["verdetti"]
    fig = plt.figure(figsize=(9, 9), facecolor="black"); ax = fig.add_axes([0, 0, 1, 1]); ax.set_facecolor("black")
    lr = np.log10(S[:, 3]); t = np.clip((lr + 2.5) / 3.0, 0, 1)
    col = np.c_[0.35 + 0.65 * t, 0.55 + 0.35 * t, 1.0 - 0.5 * t]
    for size, alpha in ((60, 0.012), (14, 0.04), (2, 0.25)):
        ax.scatter(S[:, 0], S[:, 2], s=size, c=col, alpha=alpha, linewidths=0)
    for x in B.R:
        ax.scatter([x[0]], [x[2]], s=30 if x[2] == 0 else 12, color="white")
    ax.set_xlim(-3.6, 3.6); ax.set_ylim(-2.4, 4.8); ax.set_aspect("equal"); ax.axis("off")
    T = lambda y, s_, c="white", fs=12: ax.text(0.03, y, s_, color=c, fontsize=fs, transform=ax.transAxes, va="top")
    T(0.97, "Acqua H₂O dalla regola: CCSD(T) a tutti gli elettroni, basi nostre, CBS(2,3)", fs=14)
    T(0.93, f"r(OH) {v['A1_r_OH']['valore']:.4f} Å (exp 0,9572)   angolo {v['A2_angolo']['valore']:.2f}° (104,52)", "#bbbbbb", 11)
    T(0.90, f"D_e {v['A3_D_e']['valore']:.4f} Eh (0,3705)   dipolo {v['A4_dipolo']['valore']:.3f} D (1,855)", "#bbbbbb", 11)
    w = v["A5_frequenze"]["valore"]
    T(0.87, f"frequenze armoniche {w[0]:.0f}, {w[1]:.0f}, {w[2]:.0f} cm⁻¹ (1649, 3832, 3943)", "#bbbbbb", 11)
    ax.text(0.03, 0.04, "APPROSSIMATO · 5 esami su 5 PASS · zero manopole fittate · densità da RHF (correlazione solo nei numeri)",
            color="#ffb347", fontsize=10, transform=ax.transAxes)
    ax.text(0.97, 0.08, "fetta sottile nel piano della molecola · 1 a₀ ≈ 0,53 Å", color="#888888", fontsize=10,
            transform=ax.transAxes, ha="right")
    fig.savefig(path, dpi=110, facecolor="black")


def export():
    D = os.path.join(AB.OUT, "acqua_dati"); os.makedirs(D, exist_ok=True)
    B, C, res, r, th = density(); rng = np.random.default_rng(2)
    v = json.load(open(os.path.join(AB.OUT, "acqua_verdetti.json")))
    geo = dict(unita="bohr", atomi=["O", "H", "H"], Z=[8, 1, 1], posizioni=B.R.tolist(), r_OH_A=v["minimo_cbs"][0],
               angolo_gradi=v["minimo_cbs"][1], nota="asse C2 = z, molecola nel piano xz, O nell'origine")
    json.dump(geo, open(os.path.join(D, "acqua_geometria.json"), "w"), indent=1)
    liv = dict(D_e=v["D_e"], dipolo_D=v["dipolo_D"], frequenze_cm1=v["frequenze_cm-1"], costanti_forza_au=v["costanti_forza_au"],
               E_HF_X3=res["E"], esiti={k: x["esito"] for k, x in v["verdetti"].items()})
    json.dump(liv, open(os.path.join(D, "acqua_livelli.json"), "w"), indent=1)
    S = sample(B, C, 10 * 20000, rng)[:, :3].astype(np.float32).reshape(-1, 10, 3)
    np.save(os.path.join(D, "acqua_campioni.npy"), S)
    man = {"cosa": "Acqua (H2O), tier 'approssimato'", "tier": "approssimato",
           "sorgente": {"elettronica": "RHF + CCSD(T) a tutti gli elettroni, motore molecolare lab/acqua_*.py, basi even-tempered nostre X=2,3 aug, CBS X^-3",
                        "verdetti": "lab/out/acqua_verdetti.json (pre-registrazione lab/out/acqua_prereg.json)", "documento": "docs/ACQUA.md"},
           "unita": "bohr, Hartree; nuclei fermi (Born-Oppenheimer), niente relativita', niente ZPE",
           "file": {"acqua_geometria.json": "geometria di equilibrio CBS", "acqua_livelli.json": "energie, dipolo, frequenze, esiti",
                    "acqua_campioni.npy": "float32 (n, 10, 3): posizioni dei 10 elettroni estratte INDIPENDENTI dalla densita' RHF X=3 (solo per disegnare; NON tracciato)"},
           "rigenera": "python lab/acqua_figure.py export",
           "debiti": ["densita' e campioni da RHF (correlazione solo nelle energie)", "triple perturbative", "base finita X=2,3 + X^-3",
                      "nuclei fermi: niente moto di punto zero"]}
    json.dump(man, open(os.path.join(D, "manifest.json"), "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "export":
        export()
    else:
        figure(os.path.join(AB.OUT, "acqua_figure.png"))
