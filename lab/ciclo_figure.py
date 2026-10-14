"""ciclo_figure.py - figure of CICLO_LIVELLI (lab/out/ciclo_figure.png). Panels are drawn only if their data exist.
Run: python lab/ciclo_figure.py"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
HBARC_EVM = 197.3269804e-9        # eV m
MEC2 = 0.51099895e6               # eV
ALPHA = 7.2973525693e-3
LC = HBARC_EVM / MEC2
A0 = LC / ALPHA


def load(name):
    p = os.path.join(OUT, name)
    return json.load(open(p)) if os.path.exists(p) else None


def panel_shapes(ax):
    z = np.load(os.path.join(OUT, "ciclo_p1_shapes.npz"))
    mp = load("ciclo_p1_map.json")
    cols = plt.cm.viridis(np.linspace(0, 0.9, 5))
    for i, r in enumerate(mp["chains"]["1D_auto"]):
        k = z[f"1D_auto_B{r['B']}_k"]; S = z[f"1D_auto_B{r['B']}_S_axis"]
        ax.plot(k, S, color=cols[i], lw=1.6, label=f"B={r['B']} (cella {r['h']:.2g} $\\lambda_C$)")
    k = np.linspace(0, np.pi / 2, 100)
    ax.plot(k, k * k, "k--", lw=1.2, label="punto fisso $\\kappa^2$ (Schrödinger)")
    ax.set_xlabel("$\\kappa$ = k · cella"); ax.set_ylabel("forma normalizzata $S_n(\\kappa)=G_n/c_2$")
    ax.set_title("a) la regola di ogni livello, normalizzata (1-D, catena automatica)", fontsize=9)
    ax.legend(fontsize=6.5, loc="upper left")


def panel_c4(ax):
    mp = load("ciclo_p1_map.json")
    sty = {"1D_ladder": ("o", "C0", "1-D, b=2"), "1D_auto": ("s", "C1", "1-D, automatica"),
           "3D_A": ("^", "C3", "3-D libera (catena A)"), "proton": ("v", "C4", "protone 3-D libero"),
           "T": ("D", "C2", "3-D con Taylor fissato (T)")}
    for c, (m, col, lab) in sty.items():
        L = mp["chains"][c]
        h = np.array([r["h"] for r in L]); y = np.abs(np.array([r["c4_over_c2"] for r in L]))
        ax.loglog(h, y, m, color=col, ms=6, mfc="none" if c in ("3D_A", "proton") else col, label=lab)
    hh = np.logspace(-1.3, 4.6, 50)
    ax.loglog(hh, 1 / (4 * hh * hh), "k-", lw=1, label="Dirac: $1/(4h^2)$, y = −2")
    ax.set_xlabel("cella h ($\\lambda_C$ della particella)"); ax.set_ylabel("$|c_4/c_2|$ (unità di cella)")
    ax.set_title("b) il termine relativistico per livello: irrilevante (y=−2)", fontsize=9)
    ax.legend(fontsize=6.5, loc="lower left")


def panel_ledger(ax):
    led = load("ciclo_p1_ledger.json")
    L = np.logspace(-14, 8, 600)
    kin = (HBARC_EVM / L) ** 2 / (np.sqrt((HBARC_EVM / L) ** 2 + MEC2 ** 2) + MEC2)
    ax.loglog(L, np.full_like(L, MEC2), color="0.5", lw=1, label="massa $mc^2$")
    ax.loglog(L, kin, "C0", lw=2, label="cinetico (confinamento)")
    rel = (HBARC_EVM / L) ** 4 / (8 * MEC2 ** 3)
    m = L > LC
    ax.loglog(L[m], rel[m], "C0:", lw=1.5, label="relativistico $k^4$ (y=−2)")
    coul = ALPHA * HBARC_EVM / L
    m = L <= 2 * A0
    ax.loglog(L[m], coul[m], "C3", lw=2, label="Coulomb (y=+1 rispetto al cinetico)")
    dar = ALPHA * HBARC_EVM * LC ** 2 / (8 * L ** 3)
    m = (L > LC) & (L <= 2 * A0)
    ax.loglog(L[m], dar[m], "C3:", lw=1.5, label="Darwin / spin-orbita (y=−2)")
    if led:
        g = led["gravity"]; sol = g["h2_solid_engine"]
        ax.loglog([3 * A0, 1e8], [4.75, 4.75], "C1", lw=1.5, label="legame H$_2$ (4,75 eV)")
        ax.loglog([1e-9, 1e8], [sol["eps_eV"]] * 2, "C2", lw=1.5, label="coesione H$_2$ solido (motore, statica)")
        rho = sol["rho_kg_m3"]; m2 = 2 * 1.67262192369e-27
        Eg = 4 * np.pi / 3 * 6.6743e-11 * rho * m2 * L ** 2 / 1.602176634e-19
        mm = L > 1e-9
        ax.loglog(L[mm], Eg[mm], "k", lw=2, label="gravità per molecola (y=+2)")
        ax.axvline(g["R_cross_h2_engine_m"], color="k", ls="--", lw=0.8)
        ax.text(g["R_cross_h2_engine_m"] * 0.6, 1e-9, f"{g['R_cross_h2_engine_m'] / 1e3:.0f} km", rotation=90, fontsize=7)
    for x, t in ((LC, "$\\lambda_C$"), (A0, "$a_0$")):
        ax.axvline(x, color="0.6", ls="--", lw=0.8); ax.text(x * 1.2, 1e-10, t, fontsize=8)
    for B in (16, 256, 512, 1024):
        ax.plot([B * A0 / 4096], [3e-12], "|", color="C5", ms=10)  # existing levels (cells of chain A)
    ax.set_ylim(1e-12, 1e7); ax.set_xlim(1e-14, 1e8)
    ax.set_xlabel("scala L (m)"); ax.set_ylabel("energia per particella (eV)")
    ax.set_title("c) registro delle influenze: chi comanda a ogni scala", fontsize=9)
    ax.legend(fontsize=5.8, loc="upper right", ncol=2)


def panel_gravity(ax):
    led = load("ciclo_p1_ledger.json"); g = led["gravity"]
    R = np.logspace(3, 8, 200)
    G = 6.6743e-11; MP = 1.67262192369e-27; EV = 1.602176634e-19
    sol = g["h2_solid_engine"]
    for lab, m, v, eps, col in (("H$_2$ solido (motore)", 2 * MP, sol["v_m3"], sol["eps_eV"], "C2"),
                                ("ghiaccio (solo esame)", 18.015 * 1.66053907e-27, 18.015 * 1.66053907e-27 / 917, 0.49, "C0"),
                                ("atomi H (1 Rydberg)", MP, 4 * np.pi / 3 * (5.29177e-11) ** 3, 13.6057, "C3")):
        Eg = 4 * np.pi / 3 * G * (m / v) * m * R ** 2 / EV
        ax.loglog(R / 1e3, Eg / eps, color=col, lw=2, label=lab)
    ax.axhline(1, color="k", lw=0.8)
    fy = g["f_y_implied_ice"]
    ax.axhline(fy, color="C0", ls=":", lw=1)
    ax.text(1.2, fy * 1.3, f"f_y ≈ {fy:.0e} (cedimento/coesione: non derivabile oggi)", fontsize=6.5, color="C0")
    ax.axvspan(150, 300, color="C0", alpha=0.15)
    ax.text(160, 3e-6, "corpi ghiacciati\nrotondi (esame)", fontsize=6.5)
    ax.set_xlabel("raggio del corpo (km)"); ax.set_ylabel("peso di una colonna R / coesione (per molecola)")
    ax.set_title("d) quando la gravità vince la chimica (y=+2)", fontsize=9)
    ax.set_ylim(1e-7, 1e3); ax.legend(fontsize=6.5, loc="upper left")


def main():
    have2 = os.path.exists(os.path.join(OUT, "ciclo_p2.json"))
    have3 = os.path.exists(os.path.join(OUT, "ciclo_p3.json"))
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
    panel_shapes(axs[0, 0]); panel_c4(axs[0, 1]); panel_ledger(axs[0, 2]); panel_gravity(axs[1, 0])
    if have2:
        import ciclo_p2_figure as f2
        f2.panel(axs[1, 1])
    else:
        axs[1, 1].axis("off")
    if have3:
        import ciclo_p3_figure as f3
        f3.panel(axs[1, 2])
    else:
        axs[1, 2].axis("off")
    fig.suptitle("Il ciclo tra i livelli: quanto cambia la regola salendo, chi comanda a ogni scala", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "ciclo_figure.png"), dpi=110)
    print("wrote ciclo_figure.png")


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
