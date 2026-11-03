"""acqua_basis.py - basi gaussiane even-tempered costruite e ottimizzate QUI sui nostri atomi (studio ACQUA).

Regola (scelta dichiarata, non un fit a proprieta' dell'acqua):
  - s e p di O, s di H: serie even-tempered a_k = a0 * beta^k, (a0, beta) per ogni l ottimizzati sull'energia
    ROHF dell'atomo (O 3P; H 2S, esatta -0.5).
  - polarizzazione (l >= 1 per H, l >= 2 per O): serie even-tempered centrata, a_k = zeta_l * bp^(k-(n_l-1)/2),
    zeta_l per ogni l e un bp comune ottimizzati sull'energia CCSD a tutti gli elettroni dell'atomo (O) o di H2 a
    R = 1.4 bohr (H; H2 e' nostro: AGGLOMERATO_DUE_ATOMI).
  - diffuse ("aug"): un termine in piu' per ogni l, a = a_min(l) / beta_s(atomo) (la serie s continua di un passo).
  - numero di funzioni per X (crescente): vedi NFUN.
I parametri ottimizzati stanno in lab/out/acqua_basi.json."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out"); FBAS = os.path.join(OUT, "acqua_basi.json")

# X -> {elemento: [n_s, n_p, n_d, n_f, n_g]}; le prime due (O) / la prima (H) sono la parte "HF"
NFUN = {2: {"O": [9, 5, 1], "H": [4, 1]},
        3: {"O": [11, 6, 2, 1], "H": [5, 2, 1]},
        4: {"O": [13, 8, 3, 2, 1], "H": [6, 3, 2, 1]}}
NHF = {"O": 2, "H": 1}
ZNUC = {"O": 8, "H": 1}


def load():
    return json.load(open(FBAS)) if os.path.exists(FBAS) else {}


def save(d):
    json.dump(d, open(FBAS, "w"), indent=1)


def exps(el, X, par, aug=False):
    """list of (L, exponent) for element el at size X; par = {'hf': [[a0,beta],...], 'pol': [[zeta...], bp]}"""
    n = NFUN[X][el]; out = []
    for l in range(NHF[el]):
        a0, b = par["hf"][l]
        out += [(l, a0 * b ** k) for k in range(n[l])]
    if "pol" in par:
        zs, bp = par["pol"]
        for l in range(NHF[el], len(n)):
            z = zs[l - NHF[el]]; m = n[l]
            out += [(l, z * bp ** (k - (m - 1) / 2)) for k in range(m)]
    if aug:
        bs = par["hf"][0][1]
        for l in range(len(n)):
            amin = min(a for L, a in out if L == l)
            out.append((l, amin / bs))
    return out


def shells_for(atoms, coords, X, aug=True, bas=None):
    """atoms: list of 'O'/'H'; returns (centers, Z, shells) for acqua_int.Basis"""
    bas = bas or load()
    sh = []
    for i, el in enumerate(atoms):
        par = bas[f"{el}{X}"]
        sh += [(i, L, a) for L, a in exps(el, X, par, aug)]
    return np.asarray(coords, float), [ZNUC[e] for e in atoms], sh
