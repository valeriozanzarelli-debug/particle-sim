"""panel (f) of lab/out/ciclo_figure.png: phase 3, the electron knocked out of the atom by a local field; where the
check refined (shaded) and the density of the adaptive run against the base everywhere and the coarse level alone"""
import os, json
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def panel(ax):
    m = json.load(open(os.path.join(OUT, "ciclo_p3_P_mix.json")))
    c = json.load(open(os.path.join(OUT, "ciclo_p3_P_coarse.json")))
    r = json.load(open(os.path.join(OUT, "ciclo_p3.json")))["exams"]
    x = (np.arange(len(m["dens_final"])) * 3 + 1.5 - 1093.5) * 1.6       # lambda_C from the atom
    ax.semilogy(x, np.maximum(m["dens_final_truth"], 1e-12), color="0.2", lw=2.2, label="base ovunque (verità)")
    ax.semilogy(x, np.maximum(m["dens_final"], 1e-12), color="C1", lw=1.2, ls="--", label="adattivo (grosso + toppe)")
    ax.semilogy(x, np.maximum(c["dens_final"], 1e-12), color="C0", lw=1, ls=":", label="solo livello grosso")
    ref = np.array(m["refined"], bool)
    for I in np.where(ref)[0]:
        ax.axvspan(((I * 9) - 1093.5) * 1.6, ((I * 9 + 9) - 1093.5) * 1.6, color="C1", alpha=0.12, lw=0)
    b6 = r["B6"]["value"]["P"]
    ax.set_ylim(1e-9, 1); ax.set_xlim(-1200, 1200)
    ax.set_xlabel("x dall'atomo (λ_C)"); ax.set_ylabel("densità a fine corsa (per 3 celle base)")
    ax.set_title(f"f) controllo al volo: raffina solo dove l'elettrone esce\n(arancio, {ref.sum()}/243 celle); "
                 f"risparmio {b6['saving']:.1f}x", fontsize=8.5)
    ax.legend(fontsize=6.5, loc="lower center")
