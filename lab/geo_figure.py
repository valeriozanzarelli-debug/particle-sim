"""geo_figure.py - figures for docs/FORMA_DALLE_PALLINE.md
  lab/out/geo_forma.png    drawn (hidden) shapes vs shapes rebuilt from the balls' mutual information
  lab/out/geo_inverso.png  Part B: masses -> shape (train/test), isospectral drums and their couplings
Reads the json/npz written by geo_space.py, geo_inverse.py, geo_drums.py. The colours use the drawn
coordinates ONLY to let the eye compare; the reconstruction never saw them. The hidden-circle MDS pair is
picked blindly: the first pair of (near) equal MDS eigenvalues whose angle is NOT a harmonic of the
angle of an earlier pair (new circle, not an overtone of the big loop).
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geo_core as g

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
load = lambda tag: (json.load(open(os.path.join(OUT, tag + ".json"), encoding="utf-8")),
                    np.load(os.path.join(OUT, tag + ".npz")))


def circle_pairs(X, ev, kmax=30, tol=0.05):
    """blind: list of MDS pairs (k, k+1) with equal eigenvalues that are NEW circles"""
    pairs, angles = [], []
    k = 0
    while k < kmax - 1:
        if abs(ev[k] - ev[k + 1]) < tol * ev[k]:
            phi = np.arctan2(X[:, k + 1], X[:, k])
            harmonic = any(abs(np.mean(np.exp(1j * (phi - m * a)))) > 0.3
                           for a in angles for m in range(-12, 13) if m != 0)
            if not harmonic:
                pairs.append(k); angles.append(phi)
            k += 2
        else:
            k += 1
    return pairs


def fig_space():
    fig = plt.figure(figsize=(17, 12.5))
    # ---------------- row 1: the drawn shapes
    ax = fig.add_subplot(3, 4, 1)
    t = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    ax.scatter(np.cos(t), np.sin(t), c=np.arange(48), cmap="viridis", s=18)
    ax.set_title("disegnato: anello di 48 siti"); ax.set_aspect("equal"); ax.axis("off")

    ax = fig.add_subplot(3, 4, 2, projection="3d")
    R = 18 / (2 * np.pi)
    xs, ws = np.meshgrid(np.arange(36), np.arange(18), indexing="ij")
    ax.scatter(xs, R * np.cos(2 * np.pi * ws / 18), R * np.sin(2 * np.pi * ws / 18), c=xs.ravel(),
               cmap="viridis", s=4)
    ax.set_box_aspect((35, 2 * R, 2 * R)); ax.set_title("disegnato: linea 36 x cerchio nascosto 18")
    ax.set_axis_off()

    S = g.stadium_with_hole(30, 16)
    ax = fig.add_subplot(3, 4, 3)
    ax.scatter(S.xy[:, 0], S.xy[:, 1], c=S.xy[:, 0], cmap="viridis", s=8)
    ax.set_title("disegnato: stadio con un buco"); ax.set_aspect("equal"); ax.axis("off")

    ax = fig.add_subplot(3, 4, 4)
    ax.axis("off")
    ax.text(0, 0.95, "Come si legge", fontsize=12, weight="bold", va="top")
    ax.text(0, 0.82, "Riga 1: la forma che ho usato per muovere\nle palline (la ricostruzione NON la vede).\n\n"
            "Riga 2: la forma ricostruita SOLO dall'informazione\nmutua del vuoto (MDS sulle geodetiche).\n"
            "Colore = coordinata disegnata, solo per l'occhio.\n\n"
            "Riga 3: il cerchio nascosto che compare da solo\n(coppia di componenti MDS scelta alla cieca),\n"
            "la dimensione che cambia con la scala,\ne i buchi contati (omologia persistente).",
            fontsize=10, va="top")

    # ---------------- row 2: reconstructions (vacuum)
    J, z = load("geo_space_vac_ring48")
    X, ev = g.mds(z["G0"], 4)
    ax = fig.add_subplot(3, 4, 5)
    ax.scatter(X[:, 0], X[:, 1], c=z["xy"][z["sites0"], 0], cmap="viridis", s=18)
    p = J["pieces"][0]
    ax.set_title(f"ricostruito: b0={p['betti'][0]}, b1={p['betti'][1]}, D={p['D_small']:.1f}\n"
                 f"R/estensione err {100*J['grade']['radius_err'][0]:.3f} %")
    ax.set_aspect("equal"); ax.axis("off")

    J, z = load("geo_space_vac_cyl36x18")
    X, ev = g.mds(z["G0"], 30)
    Xcyl = X
    xy = z["xy"][z["sites0"]]
    pr = circle_pairs(X, ev)
    big, small = pr[0], pr[1]
    a1 = np.arctan2(X[:, big + 1], X[:, big]); a2 = np.arctan2(X[:, small + 1], X[:, small])
    A1, A2 = np.sqrt(ev[big]), np.sqrt(ev[small])
    ax = fig.add_subplot(3, 4, 6, projection="3d")
    ax.scatter((A1 + A2 * np.cos(a2)) * np.cos(a1), (A1 + A2 * np.cos(a2)) * np.sin(a1), A2 * np.sin(a2),
               c=xy[:, 0], cmap="viridis", s=4)
    ax.set_box_aspect((1, 1, A2 / A1 * 1.2)); ax.set_axis_off()
    p = J["pieces"][0]
    ax.set_title(f"ricostruito: un anello-tubo (b1={p['betti'][1]})\nle due copie della linea si chiudono "
                 f"agli estremi;\nla sezione del tubo e' il cerchio nascosto", fontsize=10)

    J, z = load("geo_space_vac_stadium30x16")
    X, ev = g.mds(z["G0"], 3)
    ax = fig.add_subplot(3, 4, 7)
    ax.scatter(X[:, 0], X[:, 1], c=z["xy"][z["sites0"], 0], cmap="viridis", s=8)
    p = J["pieces"][0]
    ax.set_title(f"ricostruito: b0={p['betti'][0]}, b1={p['betti'][1]} (il buco), D={p['D_small']:.1f}")
    ax.set_aspect("equal"); ax.axis("off")

    J, z = load("geo_space_vac_torus16")
    ax = fig.add_subplot(3, 4, 8)
    p = J["pieces"][0]
    ax.axis("off")
    ax.text(0, 0.95, "toro 16x16 (non disegnato in 3D)", weight="bold", va="top")
    ax.text(0, 0.8, f"pezzi trovati: {J['pieces_found']} (due 'fogli' intrecciati:\n"
            f"x+y pari / dispari - il raddoppio della regola)\n"
            f"ogni foglio: b = {p['betti']} (toro: 1,2,1)\nD = {p['D_small']:.2f}\n"
            f"Spearman vs grafo causale {p['spearman_causal']:.3f}", va="top", fontsize=10)

    # ---------------- row 3: the hidden circle, dimension vs scale, barcode
    J, z = load("geo_space_vac_cyl36x18")
    ax = fig.add_subplot(3, 4, 9)
    sc = ax.scatter(Xcyl[:, small], Xcyl[:, small + 1], c=xy[:, 1], cmap="twilight", s=10)
    ax.set_title(f"il cerchio nascosto, da solo:\ncomponenti MDS {small+1}-{small+2} (colore = w disegnato)")
    ax.set_aspect("equal"); ax.axis("off")
    plt.colorbar(sc, ax=ax, fraction=0.04, label="posizione sul cerchio nascosto (0..17)")

    ax = fig.add_subplot(3, 4, 10)
    for tag, lab in (("geo_space_vac_cyl36x18", "linea x cerchio"), ("geo_space_vac_ring48", "anello"),
                     ("geo_space_vac_torus16", "toro (un foglio)")):
        _, zz = load(tag)
        rs, Nr = g.dimension_curve(zz["G0"], 40)
        Dloc = np.gradient(np.log(Nr), np.log(rs))
        h = np.median(np.where(np.eye(len(zz["G0"])) > 0, np.inf, zz["G0"]).min(axis=1))
        ax.plot(rs / h, Dloc, label=lab)
    ax.axhline(1, color="gray", lw=0.5); ax.axhline(2, color="gray", lw=0.5)
    ax.set_xscale("log"); ax.set_ylim(0, 3); ax.set_xlabel("scala r (in passi)")
    ax.set_ylabel("dimensione locale d logN / d log r")
    ax.set_title("da vicino 2-D, da lontano 1-D\n(il cerchio nascosto sparisce a grande scala)")
    ax.legend(fontsize=8)

    ax = fig.add_subplot(3, 4, 11)
    D = np.array(J["pieces"][0]["top_H1_bars"])
    full = load("geo_space_vac_cyl36x18")[1]
    from ripser import ripser
    dg = ripser(full["G0"], distance_matrix=True, maxdim=1)["dgms"][1]
    pers = dg[:, 1] - dg[:, 0]
    o = np.argsort(-pers)[:12]
    for i, k in enumerate(o):
        ax.plot([dg[k, 0], dg[k, 1]], [i, i], lw=3, color="C3" if pers[k] > 1.5 * J["pieces"][0]["h"] else "0.6")
    ax.set_yticks([]); ax.set_xlabel("scala (distanza informazionale)")
    ax.set_title("buchi (H1) della linea x cerchio:\nrosso = buchi veri; grigio = rumore del reticolo")
    ax.annotate("cerchio nascosto", (dg[o[1], 1], 1), xytext=(15, 3), arrowprops=dict(arrowstyle="->"))
    ax.annotate("linea raddoppiata\n(chiusa agli estremi)", (dg[o[0], 1], 0), xytext=(22, 6),
                arrowprops=dict(arrowstyle="->"))

    Jq, zq = load("geo_space_quench_cyl36x18")
    ax = fig.add_subplot(3, 4, 12)
    ax.axis("off")
    pv, pq = J["pieces"][0], Jq["pieces"][0]
    ax.text(0, 0.95, "Raggio nascosto / estensione visibile", weight="bold", va="top")
    ax.text(0, 0.82, f"vero (disegnato): {(18/(2*np.pi))/35:.4f}\n"
            f"dal vuoto: {pv['R_over_diam']:.4f} (err {100*J['grade']['radius_err'][0]:+.1f} %)\n"
            f"dal quench (palline lasciate andare da\nuno stato a caso): {pq['R_over_diam']:.4f} "
            f"(err {100*Jq['grade']['radius_err'][0]:+.1f} %)\n\n"
            f"dimensione vuoto: {pv['D_small']:.2f} -> {pv['D_large']:.2f}\n"
            f"dimensione quench: {pq['D_small']:.2f} -> {pq['D_large']:.2f}\n\n"
            f"MA: la linea aperta compare chiusa ad anello\n(b1 = 2 invece di 1): le due copie della\n"
            f"regola si saldano agli estremi.", va="top", fontsize=10)
    fig.suptitle("La forma letta dalle palline: nessuno dice alla ricostruzione chi e' vicino a chi",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "geo_forma.png"), dpi=110)
    print("pairs (blind circle picker):", pr[:4])


def fig_inverse():
    fig = plt.figure(figsize=(16, 5.2))
    J = json.load(open(os.path.join(OUT, "geo_inverse.json"), encoding="utf-8"))
    ax = fig.add_subplot(1, 3, 1)
    for name, rows in J["families"].items():
        ax.plot([r["k"] for r in rows], [r["test_median"] for r in rows], "o-",
                label=f"{name} ({rows[0]['p']} param.)")
    for P in ("3", "6"):
        rows = J["poly"][P]
        ax.plot([r["k"] for r in rows], [r["test_median"] for r in rows], "x--", label=f"F4 polinomio ({P} param.)")
    ax.axhline(1e-3, color="gray", lw=0.6, ls=":"); ax.text(14, 1.2e-3, "rumore dei dati", fontsize=8)
    ax.set_yscale("log"); ax.set_xlabel("k = masse usate per imparare (le altre fino a 40 = esame)")
    ax.set_ylabel("errore relativo mediano sulle masse d'esame")
    ax.set_title("dalle masse alla forma: famiglia giusta\ncon pochi parametri = predice")
    ax.legend(fontsize=7)

    D = json.load(open(os.path.join(OUT, "geo_drums.json"), encoding="utf-8"))
    z = np.load(os.path.join(OUT, "geo_drums.npz"))
    for i, (key, Xk, Ek, Uk) in enumerate((("outline1", "X1", "E1", "U1"), ("outline2", "X2", "E2", "U2"))):
        ax = fig.add_subplot(1, 3, 2 + i)
        X, E, U = z[Xk], z[Ek], z[Uk]
        ax.tripcolor(X[:, 0], X[:, 1], E, U[:, 0], shading="gouraud", cmap="magma")
        for s in D[key]:
            s = np.array(s); ax.plot(s[:, 0], s[:, 1], "w-", lw=2)
        c = D["n"]["32"]["coup1" if i == 0 else "coup2"]
        ev = D["n"]["32"]["eig1" if i == 0 else "eig2"]
        ax.set_title(f"tamburo {i+1}: stesse frequenze {np.round(ev[:3], 3).tolist()}...\n"
                     f"accoppiamenti diversi: g3={c['g3']:.3f}, g4={c['g4']:.3f}, |g123|={c['g123']:.3f}",
                     fontsize=9)
        ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "geo_inverso.png"), dpi=110)


if __name__ == "__main__":
    fig_space()
    fig_inverse()
