"""acqua_qc.py - SCF e coupled cluster per molecole (studio ACQUA), sopra gli integrali di acqua_int.

- scf(): ROHF ad alto spin (nd orbitali doppi, ns singoli alpha; ns = 0 -> RHF), Fock efficace di tipo Guest-Saunders
  (blocchi diagonali = (Fa+Fb)/2, core-aperti = Fb, aperti-virtuali = Fa), DIIS. Campo elettrico finito opzionale
  lungo z (h -> h + F z per l'elettrone, E_nuc -> E_nuc - F sum Z z).
- mo_inplace(): trasformazione AO->MO dei (pq|rs) sul posto (memoria: un solo array M^4).
- rccsd()/rccsd_t(): CCSD(T) a guscio chiuso in orbitali spaziali (equazioni a intermedi spin-adattate, forma
  standard di letteratura per riferimento RHF), controllato contro il CCSD(T) a spin-orbitali di tavg_cc.
- so_ccsd_t(): riferimento ROHF (atomi aperti, es. O 3P): usa SpinOrb/ccsd/triples di tavg_cc (sola lettura).
"""
import os, sys, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import acqua_int as I


def _JK(G, D):
    M = D.shape[0]
    J = (G.reshape(M * M, M * M) @ D.ravel()).reshape(M, M)
    K = np.empty((M, M))
    for p in range(M):
        K[p] = np.tensordot(G[p], D, axes=([0, 2], [0, 1]))
    return J, K


def sph_sym(B):
    """for a one-centre s,p basis: O_h average of an AO matrix (s-p blocks -> 0, p-p blocks -> trace/3 * identity):
    keeps orbitals spherical (same radial function for px, py, pz), as in the Hund determinant of the fine level"""
    assert set(B.sL.tolist()) <= {0, 1}
    pidx = [B.off[k] for k in range(len(B.sL)) if B.sL[k] == 1]; sidx = [B.off[k] for k in range(len(B.sL)) if B.sL[k] == 0]

    def f(A):
        out = np.zeros_like(A)
        out[np.ix_(sidx, sidx)] = A[np.ix_(sidx, sidx)]
        for i in pidx:
            for j in pidx:
                v = (A[i, j] + A[i + 1, j + 1] + A[i + 2, j + 2]) / 3
                for m in range(3):
                    out[i + m, j + m] = v
        return out
    return f


def scf(B, ints, nd, ns=0, field=0.0, tol=1e-10, maxit=200, C0=None, log=None, sym=None):
    S, h, G = ints["S"], ints["T"] + ints["V"] + field * ints["D"][2], ints["G"]
    Enuc = B.enuc() - field * float(np.sum(B.Z * B.R[:, 2]))
    w, U = np.linalg.eigh(S); keep = w > 1e-8
    X = U[:, keep] / np.sqrt(w[keep])
    if C0 is None:
        e, V = np.linalg.eigh(X.T @ h @ X); C = X @ V
    else:
        C = C0
    no = nd + ns; hist = []; Eold = 0.0
    for it in range(maxit):
        Da = C[:, :no] @ C[:, :no].T; Db = C[:, :nd] @ C[:, :nd].T
        J, Ka = _JK(G, Da)
        Jb, Kb = _JK(G, Db) if ns else (J, Ka)
        if not ns:
            Kb = Ka; Jt = 2 * J
        else:
            Jt = J + Jb
        Fa = h + Jt - Ka; Fb = h + Jt - Kb
        E = 0.5 * float(np.sum((Da + Db) * h) + np.sum(Da * Fa) + np.sum(Db * Fb)) + Enuc
        Fc = 0.5 * (Fa + Fb)
        R = C.T @ Fc @ C; Ra = C.T @ Fa @ C; Rb = C.T @ Fb @ C
        c = slice(0, nd); o = slice(nd, no); v = slice(no, C.shape[1])
        R[c, o] = Rb[c, o]; R[o, c] = Rb[o, c]; R[o, v] = Ra[o, v]; R[v, o] = Ra[v, o]
        Emo = np.zeros_like(R)
        for a_, b_ in ((c, o), (c, v), (o, v)):
            Emo[a_, b_] = R[a_, b_]; Emo[b_, a_] = -R[b_, a_]
        err = np.abs(Emo).max()
        SC = S @ C
        Rao = SC @ R @ SC.T; eao = SC @ Emo @ SC.T
        if sym is not None:
            Rao = sym(Rao); eao = sym(eao); err = np.abs(eao).max()
        hist.append((Rao, eao)); hist = hist[-8:]
        if len(hist) > 1:
            n = len(hist); Bm = -np.ones((n + 1, n + 1)); Bm[-1, -1] = 0
            for i in range(n):
                for j in range(n):
                    Bm[i, j] = np.sum(hist[i][1] * hist[j][1])
            r = np.zeros(n + 1); r[-1] = -1
            try:
                cc = np.linalg.solve(Bm, r)[:n]; Rao = sum(x * hh[0] for x, hh in zip(cc, hist))
            except np.linalg.LinAlgError:
                pass
        if log:
            log(f"  scf {it} E={E:.12f} err={err:.1e}")
        if err < tol and abs(E - Eold) < 1e-11:
            break
        Eold = E
        e, V = np.linalg.eigh(X.T @ Rao @ X); C = X @ V
    # final canonicalisation inside blocks (semicanonical for ROHF)
    for blk in (slice(0, nd), slice(nd, no), slice(no, C.shape[1])):
        if blk.stop > blk.start:
            ee, UU = np.linalg.eigh(R[blk, blk]); R[blk, blk] = np.diag(ee); C[:, blk] = C[:, blk] @ UU
    eps = np.diag(C.T @ Fc @ C).copy()
    return dict(E=E, C=C, eps=eps, h=h, Enuc=Enuc, conv=err < tol, it=it, Fa=Fa, Fb=Fb)


def mo_inplace(G, C):
    """G (M,M,M,M) AO -> MO with square C (M,M); returns G (same memory)"""
    M = G.shape[0]; assert C.shape == (M, M)
    for m in range(M):
        t = np.tensordot(G[m], C, axes=([2], [0]))          # (n, l, s)
        t = np.tensordot(t, C, axes=([1], [0]))             # (n, s, r)
        G[m] = t.transpose(0, 2, 1)
    for r in range(M):
        X = G[:, :, r, :]                                   # (m, n, s)
        X = np.tensordot(C, X, axes=([0], [0]))             # (p, n, s)
        X = np.tensordot(C, X, axes=([0], [1]))             # (q, p, s)
        G[:, :, r, :] = X.transpose(1, 0, 2)
    return G


def mo_ints(B, ints, res):
    """MO h and MO (pq|rs); handles a non-square C (removed near-dependencies) with a copy"""
    C = res["C"]; G = ints["G"]; M = G.shape[0]
    hmo = C.T @ res["h"] @ C
    if C.shape[1] == M:
        Gm = mo_inplace(G, C); ints["G"] = None
    else:
        Cp = np.zeros((M, M)); Cp[:, :C.shape[1]] = C
        Gm = mo_inplace(G, Cp)[:C.shape[1], :C.shape[1], :C.shape[1], :C.shape[1]].copy(); ints["G"] = None
    return hmo, Gm


# ------------------------------------------------------------------ closed-shell CCSD(T)
def rccsd(G, eps, nocc, nfz=0, maxit=100, tol=1e-10, log=None):
    """G: MO (pq|rs) full; eps: canonical RHF energies; nfz frozen core orbitals; returns E_ccsd, t1, t2"""
    o = slice(nfz, nocc); v = slice(nocc, G.shape[0]); no = nocc - nfz; nv = G.shape[0] - nocc
    eo, ev = eps[o], eps[v]
    ovov = G[o, v, o, v]; oovv = G[o, o, v, v]; ovvo = G[o, v, v, o]; ovoo = G[o, v, o, o]; oooo = G[o, o, o, o]
    ovvv = np.ascontiguousarray(G[o, v, v, v])
    eia = eo[:, None] - ev[None, :]; D2 = eia[:, None, :, None] + eia[None, :, None, :]
    t1 = np.zeros((no, nv)); t2 = (ovov.transpose(0, 2, 1, 3) / D2).copy()
    E = lambda t1, t2: float(2 * np.einsum("ijab,iajb", t2 + np.einsum("ia,jb->ijab", t1, t1), ovov)
                             - np.einsum("ijab,ibja", t2 + np.einsum("ia,jb->ijab", t1, t1), ovov))
    es = lambda *a: np.einsum(*a, optimize=True)
    hist = []; Eold = E(t1, t2)
    for it in range(maxit):
        Foo = 2 * es("kcld,ilcd->ki", ovov, t2) - es("kdlc,ilcd->ki", ovov, t2) \
            + 2 * es("kcld,ic,ld->ki", ovov, t1, t1) - es("kdlc,ic,ld->ki", ovov, t1, t1)
        Fvv = -2 * es("kcld,klad->ac", ovov, t2) + es("kdlc,klad->ac", ovov, t2) \
            - 2 * es("kcld,ka,ld->ac", ovov, t1, t1) + es("kdlc,ka,ld->ac", ovov, t1, t1)
        Fov = 2 * es("kcld,ld->kc", ovov, t1) - es("kdlc,ld->kc", ovov, t1)
        # T1
        r1 = es("ac,ic->ia", Fvv, t1) - es("ki,ka->ia", Foo, t1) + 2 * es("kc,kica->ia", Fov, t2) \
            - es("kc,ikca->ia", Fov, t2) + es("kc,ic,ka->ia", Fov, t1, t1) \
            + 2 * es("kcai,kc->ia", ovvo, t1) - es("kiac,kc->ia", oovv, t1) \
            + 2 * es("kdac,ikcd->ia", ovvv, t2) - es("kcad,ikcd->ia", ovvv, t2) \
            + 2 * es("kdac,kd,ic->ia", ovvv, t1, t1) - es("kcad,kd,ic->ia", ovvv, t1, t1) \
            - 2 * es("lcki,klac->ia", ovoo, t2) + es("kcli,klac->ia", ovoo, t2) \
            - 2 * es("lcki,lc,ka->ia", ovoo, t1, t1) + es("kcli,lc,ka->ia", ovoo, t1, t1)
        # T2
        tmp2 = es("kibc,ka->abic", oovv, -t1) + ovvv.transpose(1, 3, 0, 2)
        tmp = es("abic,jc->ijab", tmp2, t1); r2 = tmp + tmp.transpose(1, 0, 3, 2)
        tmp2 = es("kcai,jc->akij", ovvo, t1) + ovoo.transpose(1, 3, 0, 2)
        tmp = es("akij,kb->ijab", tmp2, t1); r2 -= tmp + tmp.transpose(1, 0, 3, 2)
        r2 += ovov.transpose(0, 2, 1, 3)
        Loo = Foo + 2 * es("lcki,lc->ki", ovoo, t1) - es("kcli,lc->ki", ovoo, t1)
        Lvv = Fvv + 2 * es("kdac,kd->ac", ovvv, t1) - es("kcad,kd->ac", ovvv, t1)
        Woooo = es("lcki,jc->klij", ovoo, t1) + es("kclj,ic->klij", ovoo, t1) + es("kcld,ijcd->klij", ovov, t2) \
            + es("kcld,ic,jd->klij", ovov, t1, t1) + oooo.transpose(0, 2, 1, 3)
        Wvoov = es("kcad,id->akic", ovvv, t1) - es("kcli,la->akic", ovoo, t1) + ovvo.transpose(2, 0, 3, 1) \
            - 0.5 * es("ldkc,ilda->akic", ovov, t2) - 0.5 * es("lckd,ilad->akic", ovov, t2) \
            - es("ldkc,id,la->akic", ovov, t1, t1) + es("ldkc,ilad->akic", ovov, t2)
        Wvovo = es("kdac,id->akci", ovvv, t1) - es("lcki,la->akci", ovoo, t1) + oovv.transpose(2, 0, 3, 1) \
            - 0.5 * es("lckd,ilda->akci", ovov, t2) - es("lckd,id,la->akci", ovov, t1, t1)
        tau = t2 + np.einsum("ia,jb->ijab", t1, t1)
        r2 += es("klij,klab->ijab", Woooo, tau)
        # Wvvvv . tau: (ac|bd) part by slices over a, t1 parts through ovvv
        for a in range(nv):
            Ga = G[nocc + a, nocc:, nocc:, nocc:]             # [c, b, d] = (ac|bd)
            r2[:, :, a, :] += np.tensordot(tau, Ga, axes=([2, 3], [0, 2]))
        X2 = es("kcbd,ijcd->kbij", ovvv, tau)
        r2 -= es("kb,kaji->ijab", t1, X2) + es("ka,kbij->ijab", t1, X2)
        tmp = es("ac,ijcb->ijab", Lvv, t2); r2 += tmp + tmp.transpose(1, 0, 3, 2)
        tmp = es("ki,kjab->ijab", Loo, t2); r2 -= tmp + tmp.transpose(1, 0, 3, 2)
        tmp = 2 * es("akic,kjcb->ijab", Wvoov, t2) - es("akci,kjcb->ijab", Wvovo, t2); r2 += tmp + tmp.transpose(1, 0, 3, 2)
        tmp = es("akic,kjbc->ijab", Wvoov, t2); r2 -= tmp + tmp.transpose(1, 0, 3, 2)
        tmp = es("bkci,kjac->ijab", Wvovo, t2); r2 -= tmp + tmp.transpose(1, 0, 3, 2)
        # move the diagonal Fock parts (Foo, Fvv above exclude f: canonical orbitals, f diagonal)
        t1n = r1 / eia; t2n = r2 / D2
        vec = np.concatenate([t1n.ravel(), t2n.ravel()]); er = vec - np.concatenate([t1.ravel(), t2.ravel()])
        hist.append((vec, er)); hist = hist[-8:]
        if len(hist) > 2:
            n = len(hist); Bm = -np.ones((n + 1, n + 1)); Bm[-1, -1] = 0
            for x in range(n):
                for y in range(n):
                    Bm[x, y] = hist[x][1] @ hist[y][1]
            rr = np.zeros(n + 1); rr[-1] = -1
            try:
                vec = sum(c_ * hh[0] for c_, hh in zip(np.linalg.solve(Bm, rr)[:n], hist))
            except np.linalg.LinAlgError:
                pass
        t1 = vec[:t1.size].reshape(no, nv); t2 = vec[t1.size:].reshape(no, no, nv, nv)
        En = E(t1, t2); d = np.abs(er).max()
        if log:
            log(f"  ccsd {it} E={En:.12f} d={d:.1e}")
        if abs(En - Eold) < tol and d < 1e-8:
            return En, t1, t2, True
        Eold = En
    return En, t1, t2, False


# ------------------------------------------------------------------ open shell through tavg_cc (spin orbitals)
class _BMO:
    def __init__(self, h, G):
        self.M = h.shape[0]; self.h = h; self.G = G

    def Gslice(self, p):
        return self.G[p]


def so_ccsd_t(hmo, Gmo, nd, ns, nfz=0, triples=True, ck=None, t_budget=None):
    """ROHF-CCSD(T) a spin-orbitali (tavg_cc). ck: file .npz di checkpoint (ampiezze, DIIS, triple parziali);
    con t_budget [s] si ferma e salva; ritorna conv=False e 'parziale'=True se non ha finito."""
    import tavg_cc as tc
    t_end = time.time() + t_budget if t_budget else None
    B = _BMO(hmo, Gmo)
    occ = [(p, s) for p in range(nd) for s in (0, 1)] + [(p, 0) for p in range(nd, nd + ns)]
    fz = [(p, s) for p in range(nfz) for s in (0, 1)]
    S = tc.SpinOrb(B, occ, frozen=fz)
    st = None; hist = []; done = False; stT = None
    if ck and os.path.exists(ck):
        z = np.load(ck); st = (z["t1"], z["t2"]); done = bool(z["conv"])
        hist = [(z[f"v{i}"], z[f"e{i}"]) for i in range(int(z["n"]))]
        if "TE" in z:
            stT = (float(z["TE"]), int(z["Ti"]))
    if done:
        t1, t2 = st; conv = True
        o = np.arange(S.no); v = np.arange(S.no, S.no + S.nv)
        oovv = S.anti(o, o, v, v); fov = S.f[np.ix_(o, v)]
        Ec = float(np.einsum("ia,ia->", fov, t1) + 0.25 * np.einsum("ijab,ijab->", oovv, t2)
                   + 0.5 * np.einsum("ijab,ia,jb->", oovv, t1, t1))
    else:
        Ec, (t1, t2), conv = tc.ccsd(S, log=lambda s: None, state=st, hist=hist, t_end=t_end)

    def save(extra=None):
        if ck:
            dd = dict(t1=t1, t2=t2, conv=conv, n=len(hist))
            for i, (a_, b_) in enumerate(hist):
                dd[f"v{i}"] = a_; dd[f"e{i}"] = b_
            if extra:
                dd.update(extra)
            np.savez(ck, **dd)
    if not conv:
        save(); return dict(E_ref=S.E_ref, Ec_ccsd=Ec, Ec_T=float("nan"), conv=False, parziale=True)
    ET = 0.0
    if triples:
        ET = tc.triples(S, t1, t2, t_end=t_end, start=stT)
        if isinstance(ET, tuple):
            save(dict(TE=ET[0], Ti=ET[1])); return dict(E_ref=S.E_ref, Ec_ccsd=Ec, Ec_T=float("nan"), conv=False, parziale=True)
    save()
    return dict(E_ref=S.E_ref, Ec_ccsd=Ec, Ec_T=ET, conv=True)


def rccsd_t(G, eps, nocc, t1, t2, nfz=0):
    """closed-shell (T), loop over i,j,k (Lee-Rendell form): W = P[sum_d (bd|ck) t_ij^ad - sum_l (ck|jl) t_il^ab],
    V = W + (bj|ck) t_i^a + (ai|ck) t_j^b + (ai|bj) t_k^c, E = sum (4W_abc + W_bca + W_cab)(V_abc - V_cba)/(3D)"""
    o = slice(nfz, nocc); v = slice(nocc, G.shape[0]); no = nocc - nfz
    eo, ev = eps[o], eps[v]
    vvvo = np.ascontiguousarray(G[v, v, v, o])        # (bd|ck) -> [b,d,c,k]
    vooo = np.ascontiguousarray(G[v, o, o, o])        # (ck|jl) -> [c,k,j,l]
    vovo = np.ascontiguousarray(G[v, o, v, o])        # (ai|bj) -> [a,i,b,j]
    Dv = -ev[:, None, None] - ev[None, :, None] - ev[None, None, :]
    E = 0.0
    for i in range(no):
        for j in range(no):
            for k in range(no):
                def w(i, j, k):        # returns [a,b,c] for (ia)(jb)(kc)
                    x = np.einsum("bdc,ad->abc", vvvo[:, :, :, k], t2[i, j])
                    return x - np.einsum("cl,lab->abc", vooo[:, k, j, :], t2[i])
                W = (w(i, j, k) + w(i, k, j).transpose(0, 2, 1) + w(j, i, k).transpose(1, 0, 2)
                     + w(j, k, i).transpose(2, 0, 1) + w(k, i, j).transpose(1, 2, 0) + w(k, j, i).transpose(2, 1, 0))
                V = W + np.einsum("bc,a->abc", vovo[:, j, :, k], t1[i]) + np.einsum("ac,b->abc", vovo[:, i, :, k], t1[j]) \
                    + np.einsum("ab,c->abc", vovo[:, i, :, j], t1[k])
                D = Dv + eo[i] + eo[j] + eo[k]
                E += np.sum((4 * W + W.transpose(1, 2, 0) + W.transpose(2, 0, 1)) * (V - V.transpose(2, 1, 0)) / (3 * D))
    return float(E)
