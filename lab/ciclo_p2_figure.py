"""panel (e) of lab/out/ciclo_figure.png: phase 2, accumulated error against the blocking scheme"""
import os, json
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def panel(ax):
    cols = {"2": "C0", "3": "C1", "phi": "C3", "fib": "C3"}
    for s in ("2", "3", "phi"):
        d = json.load(open(os.path.join(OUT, f"ciclo_p2_3d_{s}.json")))
        adm = [l for l in d["levels"] if "flag" not in l]
        ax.loglog([l["B"] for l in adm], [l["acc_rel"] for l in adm], "o-", color=cols[s], lw=1.8,
                  label=f"3-D, b = {'φ' if s == 'phi' else s} ({len(adm) - 1} livelli)")
        for l in d["levels"]:
            if "flag" in l:
                ax.plot([l["B"]], [2e-3], "x", color=cols[s], ms=9)
    for s in ("2", "3", "fib"):
        d = json.load(open(os.path.join(OUT, f"ciclo_p2_1d_free_{s}.json")))
        adm = [l for l in d["levels"] if "flag" not in l]
        B = [l["B"] * 32 for l in adm]          # in fine cells (base = 32 fine cells)
        ax.loglog(B, [l["acc_err"] for l in adm], "s--", color=cols[s], mfc="none", lw=1,
                  label=f"1-D, {'Fibonacci' if s == 'fib' else 'b = ' + s}")
        for l in d["levels"]:
            if "flag" in l:
                ax.plot([l["B"] * 32], [2e-2], "X", color=cols[s], ms=10)
    ax.axvline(400, color="0.6", ls=":", lw=0.8)
    ax.set_xlabel("cella del livello (celle fini)"); ax.set_ylabel("errore accumulato contro la regola fine")
    ax.set_title("e) fattore di blocco 2 / 3 / φ: vince il 3, φ perde (× = livello non ammesso)", fontsize=9)
    ax.legend(fontsize=6.5, loc="upper left")
