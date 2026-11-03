"""tavg_cc.py - FASE 4 dello studio tavg_: correlazione del guscio di valenza calcolata DIRETTAMENTE, con tutti gli
elettroni di valenza insieme nel blocco congelato (niente somma di coppie, niente manopole fittate).

Risolutore: coupled cluster a singole e doppie (CCSD, equazioni a spin-orbitali di Stanton-Gauss 1991, Fock non
diagonale ammessa: riferimento ROHF) + triple perturbative (T) in orbitali semicanonici con il termine f_ia di
Watts-Gauss-Bartlett 1993 (ROHF-CCSD(T)). Con 2 elettroni CCSD = full CI; il codice ha anche una full CI a
determinanti (piccola) per i controlli.

Orbitali: radiali per l = 0..X (virtuali ruotati a diagonalizzare h_l + (N-1)/N V_H: solo una rotazione, i
risultati CCSD non cambiano; servono denominatori buoni; iterazioni con DIIS a 8 vettori): gli orbitali ROHF del blocco (tr.Atom nel campo del blocco, determinante di Hund) +
X+1-l funzioni di Laguerre x^(l+1) e^(-x/2) L_k^(2l+2)(x), x = 2 lam_l r, con lam_l = (2l+3)/(2 <r>_v) (la prima
funzione ha lo stesso <r> dell'orbitale di valenza piu' esterno: regola fissa, non un fit), ortonormalizzate contro gli
occupati. Angolare: Y_lm complesse. Integrali: R^k sulla griglia del livello fine (mappa Lk di tavg_radial),
coefficienti di Gaunt c^k(lm, l'm') = (-1)^m sqrt((2l+1)(2l'+1)) (l k l'; 0 0 0)(l k l'; -m m-m' m').
Uno-elettrone: h_l del blocco (cinetico + Vloc + KB per l = 0) sulla stessa griglia.
Base crescente X = 2..5 (l_max = X), estrapolazione E_c(X) = E_c(inf) + A X^-3 da X = 4, 5.
Uso: python tavg_cc.py <Z> "<config>" <X> [nome]     es. 8 "2s2 2p4" 4"""
import os, sys, json, math, time
from fractions import Fraction
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.special import eval_genlaguerre
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_pecg as P
import tavg_radial as tr

OUT = os.path.join(HERE, "out")
LBL = "spdfghi"


# ------------------------------------------------------------------ angular
def _f(n):
    return math.factorial(n)


def w3j(j1, j2, j3, m1, m2, m3):
    if m1 + m2 + m3 != 0 or j3 < abs(j1 - j2) or j3 > j1 + j2 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    tri = Fraction(_f(j1 + j2 - j3) * _f(j1 - j2 + j3) * _f(-j1 + j2 + j3), _f(j1 + j2 + j3 + 1))
    pre = tri * _f(j1 + m1) * _f(j1 - m1) * _f(j2 + m2) * _f(j2 - m2) * _f(j3 + m3) * _f(j3 - m3)
    s = Fraction(0)
    for t in range(0, j1 + j2 + j3 + 1):
        args = (t, j3 - j2 + t + m1, j3 - j1 + t - m2, j1 + j2 - j3 - t, j1 - t - m1, j2 - t + m2)
        if min(args) < 0:
            continue
        den = 1
        for x in args:
            den *= _f(x)
        s += Fraction((-1) ** t, den)
    return (1 if (j1 - j2 - m3) % 2 == 0 else -1) * math.sqrt(pre) * float(s)


from functools import lru_cache


@lru_cache(maxsize=None)
def gaunt(l1, m1, k, l2, m2):
    return (1 if m1 % 2 == 0 else -1) * math.sqrt((2 * l1 + 1) * (2 * l2 + 1)) * w3j(l1, k, l2, 0, 0, 0) * w3j(l1, k, l2, -m1, m1 - m2, m2)


# ------------------------------------------------------------------ basis and integrals
class Basis:
    def __init__(self, Z, config, X, log=print, ae=False, core=False):
        """config: valence configuration in the block, e.g. '2s2 2p4' (block numbering: first s = 2s, first p = 2p).
        ae=True (fase 5): all electrons, bare nucleus Z, config like '1s2 2s2 2p4' (fine-level ROHF, tavg_blocks.get_ae);
        core=True adds tight Laguerre sets with the extent of the 1s: n = X-1 (l = 0), X-l (1 <= l <= X-1)"""
        t0 = time.time()
        if ae:
            import tavg_blocks as tb
            cfg = tr.config_from_string(config)
            a = tb.get_ae(Z, config)
            self.f = None
        else:
            self.f = f = P.BlockField(Z)
            cfg = []
            for tok in config.split():
                l = LBL.index(tok[1]); cfg.append((l + 1, l, int(tok[2:])))
            a = f.radial_atom(cfg); a.scf()
        self.cfg = cfg
        self.atom, self.E_rohf = a, float(a.E)
        g = a.g; r = g.r; self.g = g; w = r * g.h; self.w = w
        rv = max(float(g.integ(a.P[i] ** 2 * r)) for i in range(len(cfg)))
        self.r_v = rv
        self.r_c = min(float(g.integ(a.P[i] ** 2 * r)) for i in range(len(cfg)))
        Nel = sum(c[2] for c in cfg)
        rho = sum(c[2] * a.P[i] ** 2 for i, c in enumerate(cfg))
        self.VH = (Nel - 1) / Nel * (g.Lk(0) @ rho)
        rad, rl = [], []                      # radial functions (on grid), their l
        self.occ_rad = {}                     # (shell index) -> radial index
        H = {}
        for l in range(X + 1):
            Pl = [a.P[i] for i, c in enumerate(cfg) if c[1] == l]
            idx_occ = [i for i, c in enumerate(cfg) if c[1] == l]
            lam = (2 * l + 3) / (2 * rv)
            x = 2 * lam * r
            Lg = [x ** (l + 1) * np.exp(-x / 2) * eval_genlaguerre(k, 2 * l + 2, x) for k in range(X + 1 - l)]
            if core:
                nc = (X - 1) if l == 0 else (X - l)
                xc = 2 * (2 * l + 3) / (2 * self.r_c) * r
                Lg += [xc ** (l + 1) * np.exp(-xc / 2) * eval_genlaguerre(k, 2 * l + 2, xc) for k in range(max(nc, 0))]
            V = np.array(Pl + Lg)             # (n, M)
            # Gram-Schmidt for occupied first, then canonical orthonormalization of the rest
            Q = []
            for v in V[:len(Pl)]:
                Q.append(v / math.sqrt(np.sum(v * v * w)))
            Rst = V[len(Pl):].copy()
            for q in Q:
                Rst -= np.outer(Rst @ (q * w), q)
            S = (Rst * w) @ Rst.T
            sw, sv = np.linalg.eigh(S)
            keep = sw > 1e-10 * sw.max()
            Rn = ((sv[:, keep] / np.sqrt(sw[keep])).T @ Rst)
            Hl = a.h_matrix(l)
            Hl = 0.5 * (Hl + Hl.T)
            # virtuals rotated to diagonalize h_l + (N-1)/N V_Hartree (only a rotation inside the virtual span:
            # nothing changes for CCSD, it just gives good denominators)
            if len(Rn):
                Uv = Rn / np.sqrt(r)[None, :]
                Fv = g.h * (Uv @ (Hl + np.diag(r * r * self.VH)) @ Uv.T)
                ew, ev = np.linalg.eigh(0.5 * (Fv + Fv.T))
                Rn = ev.T @ Rn
            allf = np.array(Q + list(Rn))
            U = allf / np.sqrt(r)[None, :]
            H[l] = g.h * (U @ Hl @ U.T)
            for j, v in enumerate(allf):
                if j < len(idx_occ):
                    self.occ_rad[idx_occ[j]] = len(rad)
                rad.append(v); rl.append(l)
        self.rad = np.array(rad); self.rl = np.array(rl); Nr = len(rad)
        # spatial orbitals (radial index, l, m)
        orbs = []
        for ir in range(Nr):
            l = rl[ir]
            for m in range(l, -l - 1, -1):
                orbs.append((ir, l, m))
        self.orbs = orbs; M = len(orbs); self.M = M
        oi = np.array([o[0] for o in orbs]); ol = np.array([o[1] for o in orbs]); om = np.array([o[2] for o in orbs])
        self.oi, self.ol, self.om = oi, ol, om
        # one-electron
        h = np.zeros((M, M))
        start = {}
        for ir in range(Nr):
            start.setdefault(rl[ir], []).append(ir)
        for p in range(M):
            for q in range(M):
                if ol[p] == ol[q] and om[p] == om[q]:
                    l = ol[p]; lst = start[l]
                    h[p, q] = H[l][lst.index(oi[p]), lst.index(oi[q])]
        self.h = h
        # radial integrals R^k[a,b,c,d] = int int P_a P_b (1) Lk P_c P_d (2)
        D = (self.rad[:, None, :] * self.rad[None, :, :]).reshape(Nr * Nr, -1)
        kmax = 2 * X
        Rk = {}
        for k in range(kmax + 1):
            Lk = g.Lk(k)
            Rk[k] = ((D * w) @ (Lk @ D.T)).reshape(Nr, Nr, Nr, Nr)
            Rk[k] = 0.5 * (Rk[k] + Rk[k].transpose(2, 3, 0, 1))
        # Gaunt tables c^k[p, q]
        ck = {}
        for k in range(kmax + 1):
            c = np.zeros((M, M))
            for p in range(M):
                for q in range(M):
                    if (ol[p] + ol[q] + k) % 2 == 0 and abs(ol[p] - ol[q]) <= k <= ol[p] + ol[q] and abs(om[p] - om[q]) <= k:
                        c[p, q] = gaunt(int(ol[p]), int(om[p]), k, int(ol[q]), int(om[q]))
            ck[k] = c
        # chemist (pq|rs) = d(m_p - m_q, m_s - m_r) sum_k c^k[p,q] c^k[s,r] R^k[p q, r s]
        self._Rk, self._ck, self.kmax = Rk, ck, kmax
        self._cache = {}
        if M <= LOWMEM_M:
            G = np.zeros((M, M, M, M))
            for p in range(M):
                G[p] = self._slice(p)
            self.G = G
        else:
            self.G = LazyG(self)                           # integrals in pieces (one first index at a time)
        self.secondi = time.time() - t0
        log(f"   base X={X}: {Nr} radiali, {M} orbitali spaziali, integrali in {self.secondi:.0f} s")

    def _slice(self, p):
        """(p q|r s) for all q, r, s"""
        M, oi, om, ck = self.M, self.oi, self.om, self._ck
        mask = ((om[p] - om)[:, None, None] == (om[None, None, :] - om[None, :, None]))
        acc = np.zeros((M, M, M))
        for k in range(self.kmax + 1):
            if not np.any(ck[k][p]):
                continue
            acc += ck[k][p][:, None, None] * ck[k].T[None, :, :] * self._Rk[k][oi[p]][oi][:, oi][:, :, oi]
        return acc * mask

    def Gslice(self, p):
        if isinstance(self.G, np.ndarray):
            return self.G[p]
        v = self._cache.get(p)
        if v is None:
            v = self._slice(p)
            if len(self._cache) >= 4:
                self._cache.pop(next(iter(self._cache)))
            self._cache[p] = v
        return v

    def reference(self):
        """occupied spatial orbitals of the Hund determinant: (spatial index, spin) with spin 0 = up, 1 = down"""
        occ = []
        for (si, ms, sp) in tr.hund_determinant(self.cfg):
            ir = self.occ_rad[si]; l = self.cfg[si][1]
            p = [i for i, o in enumerate(self.orbs) if o[0] == ir and o[2] == ms][0]
            occ.append((p, 0 if sp > 0 else 1))
        return occ

    def core_spinorbs(self):
        """spin-orbitals of the 1s shell (shell 0 with n = 1) in the all-electron reference"""
        out = []
        for (si, ms, sp) in tr.hund_determinant(self.cfg):
            if self.cfg[si][0] == 1 and self.cfg[si][1] == 0:
                ir = self.occ_rad[si]
                p = [i for i, o in enumerate(self.orbs) if o[0] == ir and o[2] == ms][0]
                out.append((p, 0 if sp > 0 else 1))
        return out


LOWMEM_M = 100


class LazyG:
    """G[p] or G[np.ix_(a, b, c, d)] computed slice by slice (memory: a few M^3 arrays)"""

    def __init__(self, B):
        self.B = B

    def __getitem__(self, key):
        if isinstance(key, (int, np.integer)):
            return self.B.Gslice(int(key))
        k0, k1, k2, k3 = [np.asarray(x).ravel() for x in key]
        out = np.empty((len(k0), len(k1), len(k2), len(k3)))
        for n, p in enumerate(k0):
            out[n] = self.B.Gslice(int(p))[np.ix_(k1, k2, k3)]
        return out


# ------------------------------------------------------------------ spin-orbital machinery
class SpinOrb:
    def __init__(self, B, occ, frozen=()):
        """occ: all occupied (spatial, spin); frozen: subset kept out of the correlation (still in the Fock/reference)"""
        M = B.M
        so = [(p, s) for s in (0, 1) for p in range(M)]
        occset = set(occ)
        o_list = [x for x in so if x in occset]; v_list = [x for x in so if x not in occset]
        self.so = o_list + v_list; self.no = len(o_list); self.nv = len(v_list)
        self.sp = np.array([x[0] for x in self.so]); self.ss = np.array([x[1] for x in self.so])
        self.B = B
        n = len(self.so)
        hs = B.h[self.sp][:, self.sp] * (self.ss[:, None] == self.ss[None, :])
        o = np.arange(self.no)
        # Fock: f_pq = h_pq + sum_m <pm||qm>
        allp = np.arange(n)
        J = np.zeros((n, n))
        for m in o:                                      # <pm||qm> = <mp||mq>: first index occupied
            J += self.anti([m], allp, [m], allp)[0, :, 0, :]
        self.f = hs + J
        self.E_ref = float(np.trace(hs[:self.no, :self.no]) + 0.5 * np.einsum("ijij->", self.anti(o, o, o, o)))
        if frozen:
            fz = set(frozen)
            keep = [i for i, x in enumerate(self.so) if x not in fz]
            self.so = [self.so[i] for i in keep]; self.no -= len(fz)
            self.sp = self.sp[keep]; self.ss = self.ss[keep]; self.f = self.f[np.ix_(keep, keep)]

    def _phys(self, P1, P2, P3, P4, diag_mm=False):
        """<pq||rs> for index sets (p in P1, q in P2, r in P3, s in P4); if diag_mm returns (P1, P2, P3, P2) with q=s"""
        if diag_mm:
            out = np.zeros((len(P1), len(P2), len(P3), len(P2)))
            for j, m in enumerate(P2):
                out[:, j, :, j] = self.anti(P1, [m], P3, [m])[:, 0, :, 0]
            return out
        return self.anti(P1, P2, P3, P4)

    def anti(self, P1, P2, P3, P4):
        """<pq||rs> = <pq|rs> - <pq|sr>, <pq|rs> = (pr|qs) d(sp,sr) d(sq,ss)"""
        G = self.B.G; sp, ss = self.sp, self.ss
        P1, P2, P3, P4 = map(np.asarray, (P1, P2, P3, P4))
        a = G[np.ix_(sp[P1], sp[P3], sp[P2], sp[P4])].transpose(0, 2, 1, 3)       # (pr|qs) -> [p,q,r,s]
        a = a * ((ss[P1][:, None, None, None] == ss[P3][None, None, :, None]) & (ss[P2][None, :, None, None] == ss[P4][None, None, None, :]))
        b = G[np.ix_(sp[P1], sp[P4], sp[P2], sp[P3])].transpose(0, 2, 3, 1)       # (ps|qr) -> [p,q,r,s]
        b = b * ((ss[P1][:, None, None, None] == ss[P4][None, None, None, :]) & (ss[P2][None, :, None, None] == ss[P3][None, None, :, None]))
        return a - b


def ccsd(S, maxit=200, tol=1e-9, t_end=None, state=None, log=print, hist=None):
    no, nv = S.no, S.nv
    o = np.arange(no); v = np.arange(no, no + nv)
    f = S.f
    foo, fvv, fov = f[np.ix_(o, o)], f[np.ix_(v, v)], f[np.ix_(o, v)]
    oooo = S.anti(o, o, o, o); ooov = S.anti(o, o, o, v); oovv = S.anti(o, o, v, v)
    ovov = S.anti(o, v, o, v)
    ovvv = np.empty((no, nv, nv, nv), dtype=np.float32)   # stored in single precision, used in double per m
    for m in range(no):
        ovvv[m] = S.anti([m], v, v, v)[0]
    ov = lambda m: ovvv[m].astype(float)
    ovvo = -ovov.transpose(0, 1, 3, 2)                         # <mb||ej> = -<mb||je>
    Do = np.diag(foo); Dv = np.diag(fvv)
    D1 = Do[:, None] - Dv[None, :]
    D2 = Do[:, None, None, None] + Do[None, :, None, None] - Dv[None, None, :, None] - Dv[None, None, None, :]
    if state is not None:
        t1, t2 = state
    else:
        t1 = fov / D1; t2 = oovv / D2
    sp, ss = S.sp, S.ss
    vsp, vss = sp[v], ss[v]

    def ladder(tau):
        """0.5 sum_ef tau_ij^ef <ab||ef>, batched over a with spatial integrals"""
        out = np.zeros((no, no, nv, nv))
        for ia in np.argsort(vsp, kind="stable"):
            pa, sa = vsp[ia], vss[ia]
            # <ab||ef> = (ae|bf) d(sa,se) d(sb,sf) - (af|be) d(sa,sf) d(sb,se)
            A = S.B.Gslice(pa)[np.ix_(vsp, vsp, vsp)]           # [e, b, f] = (a e|b f)
            m1 = (vss == sa)[:, None, None] & (vss[None, :, None] == vss[None, None, :])
            A *= m1
            W = A.transpose(1, 0, 2)                             # [b, e, f] (view)
            W = W - W.transpose(0, 2, 1); del A
            out[:, :, ia, :] = 0.5 * np.einsum("ijef,bef->ijb", tau, W, optimize=True)
        return out

    E_old = 0.0
    hist = [] if hist is None else hist
    for it in range(maxit):
        tau = t2 + np.einsum("ia,jb->ijab", t1, t1) - np.einsum("ib,ja->ijab", t1, t1)
        taut = t2 + 0.5 * (np.einsum("ia,jb->ijab", t1, t1) - np.einsum("ib,ja->ijab", t1, t1))
        Fae = fvv - np.diag(Dv) - 0.5 * np.einsum("me,ma->ae", fov, t1) - 0.5 * np.einsum("mnaf,mnef->ae", taut, oovv)
        for m in range(no):
            Fae += np.einsum("f,afe->ae", t1[m], ov(m))
        Fmi = foo - np.diag(Do) + 0.5 * np.einsum("ie,me->mi", t1, fov) + np.einsum("ne,mnie->mi", t1, ooov) \
            + 0.5 * np.einsum("inef,mnef->mi", taut, oovv)
        Fme = fov + np.einsum("nf,mnef->me", t1, oovv)
        Wmnij = oooo + np.einsum("je,mnie->mnij", t1, ooov) - np.einsum("ie,mnje->mnij", t1, ooov) \
            + 0.25 * np.einsum("ijef,mnef->mnij", tau, oovv)
        Wmbej = ovvo + np.stack([np.einsum("jf,bef->bej", t1, ov(m)) for m in range(no)]) + np.einsum("nb,mnje->mbej", t1, ooov) \
            - np.einsum("jnfb,mnef->mbej", 0.5 * t2 + np.einsum("jf,nb->jnfb", t1, t1), oovv)
        # T1
        r1 = fov.copy() + np.einsum("ie,ae->ia", t1, Fae) - np.einsum("ma,mi->ia", t1, Fmi) \
            + np.einsum("imae,me->ia", t2, Fme) - np.einsum("nf,naif->ia", t1, ovov) \
            + 0.5 * np.einsum("mnae,nmie->ia", t2, ooov)
        for m in range(no):
            r1 -= 0.5 * np.einsum("ief,aef->ia", t2[:, m], ov(m))
        # T2
        r2 = oovv.copy()
        Xe = Fae - 0.5 * np.einsum("mb,me->be", t1, Fme)
        tmp = np.einsum("ijae,be->ijab", t2, Xe); r2 += tmp - tmp.transpose(0, 1, 3, 2)
        Xm = Fmi + 0.5 * np.einsum("je,me->mj", t1, Fme)
        tmp = np.einsum("imab,mj->ijab", t2, Xm); r2 -= tmp - tmp.transpose(1, 0, 2, 3)
        r2 += 0.5 * np.einsum("mnab,mnij->ijab", tau, Wmnij)
        # Wabef contribution (see module doc): ladder + t1 piece + oooo piece
        r2 += ladder(tau)
        Zx = np.stack([0.5 * np.einsum("ijef,aef->ija", tau, ov(m)) for m in range(no)], -1)   # <am||ef> = -<ma||ef>
        tmp = np.einsum("ijam,mb->ijab", -Zx, t1); r2 -= tmp - tmp.transpose(0, 1, 3, 2)
        Y = np.einsum("mnef,ijef->mnij", oovv, tau)
        r2 += 0.125 * np.einsum("mnab,mnij->ijab", tau, Y)
        tmp = np.einsum("imae,mbej->ijab", t2, Wmbej) - np.einsum("ie,ma,mbej->ijab", t1, t1, ovvo)
        r2 += tmp - tmp.transpose(1, 0, 2, 3) - tmp.transpose(0, 1, 3, 2) + tmp.transpose(1, 0, 3, 2)
        # P(ij) t_i^e <ab||ej>, <ab||ej> = -<je||ab>... = -ovvv[j,e,a,b] -> careful: <ab||ej> = <ej||ab> = -<je||ab>
        tmp = np.stack([-np.einsum("ie,eab->iab", t1, ov(j)) for j in range(no)], 1); r2 += tmp - tmp.transpose(1, 0, 2, 3)
        # -P(ab) t_m^a <mb||ij>, <mb||ij> = <ij||mb> = ooov[i,j,m,b]
        tmp = np.einsum("ma,ijmb->ijab", t1, ooov); r2 -= tmp - tmp.transpose(0, 1, 3, 2)
        t1n = r1 / D1; t2n = r2 / D2
        dt = max(np.abs(t1n - t1).max(), np.abs(t2n - t2).max())
        vec = np.concatenate([t1n.ravel(), t2n.ravel()]); err = vec - np.concatenate([t1.ravel(), t2.ravel()])
        hist.append((vec.astype(np.float32), err.astype(np.float32))); hist[:] = hist[-6:]
        if len(hist) >= 3:
            nh = len(hist)
            Bm = -np.ones((nh + 1, nh + 1)); Bm[-1, -1] = 0
            for x in range(nh):
                for y in range(nh):
                    Bm[x, y] = float(hist[x][1].astype(float) @ hist[y][1].astype(float))
            rhs = np.zeros(nh + 1); rhs[-1] = -1
            try:
                cc = np.linalg.solve(Bm, rhs)[:nh]
                vec = sum(c * h[0].astype(float) for c, h in zip(cc, hist))
            except np.linalg.LinAlgError:
                pass
        t1 = vec[:t1.size].reshape(t1.shape); t2 = vec[t1.size:].reshape(t2.shape)
        E = float(np.einsum("ia,ia->", fov, t1) + 0.25 * np.einsum("ijab,ijab->", oovv, t2)
                  + 0.5 * np.einsum("ijab,ia,jb->", oovv, t1, t1))
        log(f"      CCSD it {it}: E_c = {E:.10f}  dt = {dt:.1e}")
        if abs(E - E_old) < tol and dt < 1e-7:
            return E, (t1, t2), True
        E_old = E
        if t_end and time.time() > t_end:
            return E, (t1, t2), False
    return E, (t1, t2), False


def triples(S, t1, t2, log=print, t_end=None, start=None):
    """returns E(T), or (partial E, next triple index) if t_end is reached (resumable with start=(E, idx))"""
    """ROHF-CCSD(T) in semicanonical orbitals (Watts, Gauss, Bartlett 1993): E[T] + E[ST] with the f_ia term"""
    no, nv = S.no, S.nv
    o = np.arange(no); v = np.arange(no, no + nv)
    f = S.f
    wo, Uo = np.linalg.eigh(f[np.ix_(o, o)]); wv, Uv = np.linalg.eigh(f[np.ix_(v, v)])
    # transform amplitudes and integrals to the semicanonical basis (real orthogonal rotations)
    t1 = Uo.T @ t1 @ Uv
    t2 = np.einsum("ijab,iI,jJ,aA,bB->IJAB", t2, Uo, Uo, Uv, Uv, optimize=True)
    fov = Uo.T @ f[np.ix_(o, v)] @ Uv
    oovv = np.einsum("ijab,iI,jJ,aA,bB->IJAB", S.anti(o, o, v, v), Uo, Uo, Uv, Uv, optimize=True)
    ovvv = np.empty((no, nv, nv, nv), dtype=np.float32)
    for m in range(no):
        x = S.anti([m], v, v, v)[0]
        x = np.einsum("abc,aA->Abc", x, Uv, optimize=True)
        x = np.einsum("Abc,bB->ABc", x, Uv, optimize=True)
        ovvv[m] = np.einsum("ABc,cC->ABC", x, Uv, optimize=True)
    flat = ovvv.reshape(no, -1)
    for c0 in range(0, flat.shape[1], 1 << 20):
        sl = slice(c0, c0 + (1 << 20))
        flat[:, sl] = (Uo.T @ flat[:, sl].astype(float)).astype(np.float32)
    ooov = np.einsum("ijka,iI,jJ,kK,aA->IJKA", S.anti(o, o, o, v), Uo, Uo, Uo, Uv, optimize=True)
    # single precision for the triples arrays (memory); energy accumulated in double
    f32 = np.float32
    t1, t2, fov, oovv, ooov = (x.astype(f32) for x in (t1, t2, fov, oovv, ooov))
    wo32, wv32 = wo.astype(f32), wv.astype(f32)
    E, i0 = (0.0, 0) if start is None else (float(start[0]), int(start[1]))
    trip = [(i, j, k) for i in range(no) for j in range(i + 1, no) for k in range(j + 1, no)]
    for it_, (i, j, k) in enumerate(trip):
        if it_ < i0:
            continue
        if t_end and time.time() > t_end:
            return (E, it_)
        if True:
            if True:
                def conn(i, j, k):
                    # sum_e t_jk^ae <ei||bc> - sum_m t_im^bc <ma||jk>;  <ei||bc> = -<ie||bc>; <ma||jk> = <jk||ma>
                    return -np.einsum("ae,ebc->abc", t2[j, k], ovvv[i]) - np.einsum("mbc,ma->abc", t2[i], ooov[j, k])

                def disc(i, j, k):
                    return np.einsum("a,bc->abc", t1[i], oovv[j, k]) + np.einsum("a,bc->abc", fov[i], t2[j, k])

                def Pabc(x):
                    return x - x.transpose(1, 0, 2) - x.transpose(2, 1, 0)

                c = Pabc(conn(i, j, k) - conn(j, i, k) - conn(k, j, i))
                d = Pabc(disc(i, j, k) - disc(j, i, k) - disc(k, j, i))
                d += c
                c /= (wo32[i] + wo32[j] + wo32[k] - wv32[:, None, None] - wv32[None, :, None] - wv32[None, None, :])
                E += float(np.sum(c * d, dtype=np.float64)) / 6.0       # (c/D) (c + d)
                del c, d           # i<j<k restricted: 1/36 * 6 ... sum over a,b,c full -> /6
    return float(E)


def fci(S, nexc=None):
    """small determinant full CI with the same spin-orbitals (fixed M_S, M_L of the reference); for checks only"""
    from itertools import combinations
    no = S.no; n = len(S.so)
    ml = np.array([S.B.om[p] for p in S.sp]); ms = S.ss
    ML = int(sum(ml[:no])); MS = int(sum(ms[:no]))
    dets = [c for c in combinations(range(n), no) if sum(ml[list(c)]) == ML and sum(ms[list(c)]) == MS]
    allp = np.arange(n)
    hs = S.B.h[S.sp][:, S.sp] * (S.ss[:, None] == S.ss[None, :])
    A = S.anti(allp, allp, allp, allp)
    N = len(dets); Hm = np.zeros((N, N)); sets = [set(d) for d in dets]
    for I in range(N):
        dI = dets[I]
        for Jx in range(I, N):
            dJ = dets[Jx]
            diff = sets[I] ^ sets[Jx]
            if len(diff) > 4:
                continue
            # bring to maximal coincidence via permutation sign
            if len(diff) == 0:
                e = sum(hs[p, p] for p in dI) + 0.5 * sum(A[p, q, p, q] for p in dI for q in dI)
            else:
                pa = [p for p in dI if p not in sets[Jx]]; qa = [q for q in dJ if q not in sets[I]]
                sI = list(dI); sJ = list(dJ)
                sign = 1
                for p in pa:
                    pos = sI.index(p); sign *= (-1) ** pos; sI.pop(pos); sI.insert(0, p)
                for q in qa:
                    pos = sJ.index(q); sign *= (-1) ** pos; sJ.pop(pos); sJ.insert(0, q)
                # after moving, the moved ones are in reverse order at the front: same reversal on both sides
                if len(pa) == 1:
                    p, q = pa[0], qa[0]
                    e = hs[p, q] + sum(A[p, m, q, m] for m in dI if m != p)
                else:
                    e = A[sI[0], sI[1], sJ[0], sJ[1]]
                e *= sign
            Hm[I, Jx] = Hm[Jx, I] = e
    return float(np.linalg.eigvalsh(Hm)[0]), N


def run(Z, config, X, name=None, t_budget=540, log=print, ae=False, core=False, frozen_core=False):
    """one species, one basis size: E_ROHF(block), E_ref, CCSD, (T); checkpoint of amplitudes (npz, untracked)"""
    t_end = time.time() + t_budget
    B = Basis(Z, config, X, log=log, ae=ae, core=core)
    S = SpinOrb(B, B.reference(), frozen=B.core_spinorbs() if frozen_core else ())
    ck = os.path.join(OUT, "tavg_cc_amp", f"{name or ('Z%d_' % Z + config.replace(' ', '_'))}_X{X}.npz")
    os.makedirs(os.path.dirname(ck), exist_ok=True)
    st = None; done = False
    if os.path.exists(ck):
        z = np.load(ck); st = (z["t1"], z["t2"]); done = bool(z["conv"]) if "conv" in z else False
    if done:
        t1, t2 = st; conv = True
        o = np.arange(S.no); v = np.arange(S.no, S.no + S.nv)
        oovv = S.anti(o, o, v, v); fov = S.f[np.ix_(o, v)]
        Ec = float(np.einsum("ia,ia->", fov, t1) + 0.25 * np.einsum("ijab,ijab->", oovv, t2)
                   + 0.5 * np.einsum("ijab,ia,jb->", oovv, t1, t1))
    else:
        ckD = ck.replace(".npz", "_diis.npz"); hist = []
        if st is not None and os.path.exists(ckD):
            zd = np.load(ckD); hist = [(zd[f"v{i}"], zd[f"e{i}"]) for i in range(int(zd["n"]))]
        Ec, (t1, t2), conv = ccsd(S, t_end=t_end, state=st, log=lambda s: None, hist=hist)
        np.savez(ck, t1=t1, t2=t2, conv=conv)
        dd = {"n": len(hist)}
        for i, (a_, b_) in enumerate(hist):
            dd[f"v{i}"] = a_; dd[f"e{i}"] = b_
        np.savez(ckD, **dd)
    ckT = ck.replace(".npz", "_T.npz")
    out = dict(Z=Z, config=config, X=X, M=B.M, no=S.no, nv=S.nv, E_rohf_blocco=B.E_rohf, E_ref=S.E_ref,
               dE_ref=S.E_ref - B.E_rohf, Ec_ccsd=Ec, convergiuto=conv)
    if conv:
        stT = None
        if os.path.exists(ckT):
            z = np.load(ckT); stT = (float(z["E"]), int(z["i"]))
        rT = triples(S, t1, t2, t_end=t_end, start=stT)
        if isinstance(rT, tuple):
            np.savez(ckT, E=rT[0], i=rT[1]); out["convergiuto"] = False; out["T_parziale"] = list(rT)
            return out
        out["Ec_T"] = rT
        out["E_tot"] = S.E_ref + Ec + out["Ec_T"]
        out["Ec_val"] = out["E_tot"] - B.E_rohf
    return out


if __name__ == "__main__":
    print(json.dumps(run(int(sys.argv[1]), sys.argv[2], int(sys.argv[3])), indent=1))
