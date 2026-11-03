"""tavg_f5_export.py - export dell'ossigeno, tier "approssimato" (fase 5 PASS: K5, C1-C5 con l'emendamento 1).
Scrive lab/out/tavg_dati_ossigeno/: manifest.json (tracciato), ossigeno_livelli.json (tracciato), e i campioni
NON tracciati ossigeno_densita.npz e ossigeno_campioni.npy (rigenerare: python lab/tavg_f5_export.py)."""
import os, sys, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_blocks as tb

OUT = os.path.join(HERE, "out"); D = os.path.join(OUT, "tavg_dati_ossigeno"); HA_EV = 27.211386245988


def main(n=20000, seed=8):
    os.makedirs(D, exist_ok=True)
    V = json.load(open(os.path.join(OUT, "tavg_f5_verdetti_em1.json")))
    assert V["export"], "export non permesso dai verdetti"
    pr = V["previsioni"]
    a = tb.get_ae(8, "1s2 2s2 2p4"); r = a.g.r
    # spherically averaged HF density (Hund determinant, shells summed)
    rho = sum(c[2] * a.P[i] ** 2 for i, c in enumerate(a.config)) / (4 * math.pi * r * r)
    np.savez(os.path.join(D, "ossigeno_densita.npz"), r=r, rho=rho)
    # independent electron positions from the HF radial density (uncorrelated: declared)
    rng = np.random.default_rng(seed)
    w = 4 * math.pi * r * r * rho * r * a.g.h; cdf = np.cumsum(w) / w.sum()
    rr = np.interp(rng.random((n, 8)), cdf, r)
    u = rng.normal(size=(n, 8, 3)); u /= np.linalg.norm(u, axis=-1, keepdims=True)
    np.save(os.path.join(D, "ossigeno_campioni.npy"), (u * rr[..., None]).astype(np.float32))
    lev = dict(E_HF=pr["O"]["E_HF"], E_c=pr["O"]["Ec"], E_totale=pr["O"]["E_HF"] + pr["O"]["Ec"],
               E_HF_catione=pr["O+"]["E_HF"], E_c_catione=pr["O+"]["Ec"], IP_eV=V["C2"]["valore"],
               esami={k: V[k]["esito"] for k in ("K5", "C1", "C2", "C3", "C4", "C5")},
               chiavi=dict(E_c=-0.257902, IP_eV=13.618))
    json.dump(lev, open(os.path.join(D, "ossigeno_livelli.json"), "w"), indent=1)
    man = {"cosa": "Ossigeno (O 3P), tier 'approssimato'",
           "tier": "approssimato",
           "sorgente": {"HF": "ROHF del livello fine a tutti gli elettroni (lab/tavg_radial.py via tavg_blocks.get_ae)",
                        "correlazione": "CCSD(T) a tutti gli elettroni, nucleo nudo, lab/tavg_cc.py (fase 5 + emendamento 1): valenza X=4,5 e core X=4,5, estrapolazione X^-3, zero manopole",
                        "blocco": "il blocco-nucleo dell'agglomeratore NON e' usato in questo export",
                        "verdetti": "lab/out/tavg_f5_verdetti_em1.json (seconda prova dopo il C1 della fase 4: non primo cieco)",
                        "documento": "docs/TAVOLA_GRIGLIA.md, fasi 4-5"},
           "unita": "bohr, Hartree; nucleo nell'origine, fermo (massa infinita, niente relativita')",
           "file": {"ossigeno_livelli.json": "energie (HF, correlazione, totale, catione, IP) ed esiti degli esami",
                    "ossigeno_densita.npz": "r, rho(r) media sferica della densita' HF (NON tracciato)",
                    "ossigeno_campioni.npy": "float32 (n, 8, 3) posizioni degli 8 elettroni estratte INDIPENDENTI dalla densita' HF radiale (senza correlazione, senza forma p: solo per disegnare la nuvola; NON tracciato)"},
           "rigenera": "python lab/tavg_f5_export.py",
           "debiti": ["densita' e campioni da HF (la correlazione e' solo nelle energie)", "triple perturbative", "base finita + X^-3"]}
    json.dump(man, open(os.path.join(D, "manifest.json"), "w"), indent=1, ensure_ascii=False)
    print("export fatto:", D)


if __name__ == "__main__":
    main()
