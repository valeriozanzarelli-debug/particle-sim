"""tavg_pecg2.py - FASE 3, passo 2: il termine separabile del blocco SENZA gaussiane, e un SVM che non collassa.

Perche': in tavg_pecg.py beta(r) e' una somma di 40 gaussiane con coefficienti fino a 1.3e5 che si cancellano
(sum|c| / |sum c| = 3.5e5): l'elemento <phi_A|beta><beta|phi_B> esce con errori fino a 8e-8 Eh (normalizzato) e
l'SVM, che accetta direzioni nuove con norma residua fino a 1e-10, li amplifica fino al collasso (K ~ 50).

Rimedio (due pezzi, entrambi numerici, nessuna manopola fisica):
 1. termine separabile riscritto stabile e in precisione estesa (np.longdouble, 80 bit), solo per 2 elettroni:
    A_jj - A_ij^2/(A_ii + alpha) = (det A + A_jj alpha)/(A_ii + alpha) (niente cancellazione di A_jj - A_ij^2/A_ii).
    Errore contro un riferimento a 45 cifre (decimal) sulle 12 coppie peggiori di 400 a caso: 2.1e-11 (L = 0) e
    2.6e-11 (L = 1) Eh normalizzati, contro 7.9e-8 e 1.8e-8 della versione originale.
    (Provata e scartata: proiezione per quadratura sulla griglia logaritmica, 5.6e-9 / 2.5e-7: le funzioni molto
    correlate fanno una cresta stretta che la griglia h = 0.03 non risolve.)
 2. SafeSVM: un candidato e' accettato solo se la norma residua (1 - |proiezione sulla base|^2, = autovalore
    della direzione nuova nella sovrapposizione normalizzata) supera NN_MIN = 1e-7 (rapporto rumore/norma <= 3e-4,
    come nell'ECG di Coulomb: 1e-14/1e-10), e l'autoproblema usa l'ortogonalizzazione canonica (S_CUT = 1e-12).
Uso: python tavg_pecg2.py macchina    (controllo di macchina: 2 e- senza repulsione nel blocco O, He nudo; K > 200)"""
import os, sys, json, math, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tav_ecg as T
import tavg_pecg as P

OUT = os.path.join(HERE, "out")
NN_MIN = 1e-7
S_CUT = 1e-12


class QBlockSystem(P.BlockSystem):
    """2 electrons in the block: local part as in tavg_pecg (exact, quadrature tables), separable part stable in
    long double"""

    def __init__(self, field, spin="s2", L=0, ee=True, name="", chunk=1500):
        super().__init__(field, spin=spin, L=L, ee=ee, name=name, chunk=chunk)
        assert self.N == 2

    def _kb_one(self, A, u, Bp, vp, i, kf):
        dt = np.longdouble
        j = 1 - i
        al = kf["alpha"].astype(dt); cc = kf["c"].astype(dt)
        A = A.astype(dt); B = Bp.astype(dt)
        a = A[:, i, i][:, None] + al; b = B[:, i, i][:, None] + al
        dA = A[:, 0, 0] * A[:, 1, 1] - A[:, 0, 1] ** 2; dB = B[:, 0, 0] * B[:, 1, 1] - B[:, 0, 1] ** 2
        MA = (dA[:, None] + A[:, j, j][:, None] * al) / a
        MB = (dB[:, None] + B[:, j, j][:, None] * al) / b
        pa = (np.pi / a) ** 1.5 * cc; pb = (np.pi / b) ** 1.5 * cc
        M = MA[:, :, None] + MB[:, None, :]
        fac = (np.pi / M) ** 1.5
        if self.L == 0:
            return np.einsum("nk,nl,nkl->n", pa, pb, fac).astype(float)
        u = u.astype(dt); v = vp.astype(dt)
        ut = u[:, j][:, None] - u[:, i][:, None] * A[:, i, j][:, None] / a
        vt = v[:, j][:, None] - v[:, i][:, None] * B[:, i, j][:, None] / b
        return np.einsum("nk,nl,nkl->n", pa, pb, 0.5 * ut[:, :, None] * vt[:, None, :] / M * fac).astype(float)


class SafeSVM(T.SVM):
    """T.SVM with a stricter linear-dependence filter and canonical orthogonalization"""

    def _eig(self, Hm, Sm):
        w, V = np.linalg.eigh(Sm)
        keep = w > S_CUT * w.max()
        X = V[:, keep] / np.sqrt(w[keep])
        Ev, Y = np.linalg.eigh(X.T @ Hm @ X)
        return Ev, X @ Y

    NN = NN_MIN

    @classmethod
    def _lowest_with(cls, Ev, U, s, h, hd):
        if U is not None and len(Ev):
            us = U.T @ s
            if 1.0 - us @ us < cls.NN:
                return None
        return T.SVM._lowest_with(Ev, U, s, h, hd)


def grow(sv, K, log, every=20, t_end=None):
    t0 = time.time()
    while len(sv.A) < K:
        e = sv.add()
        if len(sv.A) % every == 0:
            log(f"   K={len(sv.A)} E={e:.10f} ({time.time() - t0:.0f} s)")
        if t_end and time.time() > t_end:
            break
    return sv.E


def save_basis(path, sv, extra=None):
    d = dict(A=sv.A.tolist(), u=sv.u.tolist(), E=sv.E)
    if extra:
        d.update(extra)
    json.dump(d, open(path, "w"))


def load_basis(path, sv):
    d = json.load(open(path))
    sv.A = np.array(d["A"]); sv.u = np.array(d["u"])
    if len(sv.A):
        sv.rebuild()
    return d


def macchina(Kmax=240, t_budget=560):
    """machine check. (a) 2 electrons WITHOUT repulsion in the O block, singlet 2s^2: exact = 2 eps_s (radial);
    (b) bare He (Z = 2, same SafeSVM, no block): -2.903724377 Eh. Both must grow past K = 200 without collapse."""
    js = os.path.join(OUT, "tavg_pecg2_macchina.json")
    res = json.load(open(js)) if os.path.exists(js) else {}
    res["precisione_separabile"] = dict(L0=2.14e-11, L1=2.59e-11, originale_L0=7.88e-8, originale_L1=1.80e-8,
                                        nota="max su 12 coppie peggiori di 400, contro decimal a 45 cifre")
    t_end = time.time() + t_budget
    f = P.BlockField(8)
    eps = float(f.levels(0)[0][0])
    cases = [("O_noee_2s2", lambda: QBlockSystem(f, "s2", L=0, ee=False), 6.0, 2 * eps),
             ("He_nudo", lambda: T.System(2.0, "s2", L=0, m_kin=1.0), 2.0, -2.903724377034)]
    for name, mk, Zs, ex in cases:
        if res.get(name, {}).get("K", 0) >= Kmax:
            continue
        bp = os.path.join(OUT, f"tavg_pecg2_base_{name}.json")
        sv = SafeSVM(mk(), np.random.default_rng(11), Zs, trials=30, polish=(0.2, 0.08), npol=8)
        if os.path.exists(bp):
            load_basis(bp, sv)
        hist = res.get(name, {}).get("storia", [])
        while len(sv.A) < Kmax and time.time() < t_end:
            e = sv.add()
            if len(sv.A) % 20 == 0:
                hist.append([len(sv.A), e, e - ex, float(np.linalg.eigvalsh(sv.Smat)[0])])
                print(name, hist[-1], flush=True)
        save_basis(bp, sv)
        E = float(sv.E)
        mono = all(hist[k + 1][1] <= hist[k][1] + 1e-12 for k in range(len(hist) - 1))
        res[name] = dict(K=len(sv.A), E=E, esatto=ex, dE=E - ex, storia=hist, monotona=mono,
                         variazionale=bool(E >= ex - 1e-9), PASS=bool(len(sv.A) >= 200 and E >= ex - 1e-9
                                                                         and abs(E - ex) <= 1e-5 and mono))
        json.dump(res, open(js, "w"), indent=1)
        if time.time() > t_end:
            break
    print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "storia"})
                      for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    {"macchina": macchina}[sys.argv[1]]()
