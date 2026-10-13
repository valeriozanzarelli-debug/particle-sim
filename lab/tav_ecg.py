"""tav_ecg.py - studio TAVOLA (litio e berillio, docs/TAVOLA_LITIO_BERILLIO.md): livello ESATTO a pochi corpi per
3 e 4 elettroni con l'operatore della regola (stessa strada dell'elio del thread principale, docs/MOLTI_ELETTRONI.md,
e di H2 in docs/AGGLOMERATO_DUE_ATOMI.md par. 12).

Operatore (unita' atomiche, nucleo fermo nell'origine, DEBT-3):
   H = sum_i p_i^2 / (2 m_kin) - Z sum_i 1/r_i + sum_{i<j} 1/r_ij
   m_kin = massa cinetica imparata della regola (livello T6, multi_e_v5.json: m_T(theta_A) = 1,000373269), poi
   theta -> 0 con la scala esatta E(m) = m E(1) (m0 = 1,0000000046).

Base: gaussiane correlate esplicite (ECG), scritte qui da zero (stesse formule di lab/multi_e_ecg.py del thread
principale, piu' le funzioni P):
   phi(x) = (u^T x_z)^l exp(-x^T A x),  x = (N, 3) posizioni degli elettroni, l = 0 (stati S) o 1 (stati P, M = 0).
   (u^T x_z) exp(-x^T A x) e' esattamente L = 1, M = 0 (vettore per funzione scalare), parita' dispari.
Antisimmetria con lo spin giusto (DEBT-2): funzione di spin accoppiata chi (doppietto per 3 e, singoletto per 4 e),
somma su TUTTE le N! permutazioni con il peso c_P = sign(P) <chi | P chi>. Cosi' Be 2s2p 1P e' davvero singoletto
(con la sola funzione primitiva a M_S = 0 il minimo sarebbe 3P, piu' basso).

Elementi di matrice (C = A + B', Sigma = C^-1 / 2 covarianza per componente, S0 = (pi^N / det C)^{3/2}):
   l = 0:  S = S0; T = S0 (3/m) tr[A C^-1 B']; <1/|w.x|> = S0 sqrt(2/pi) / sigma_w, sigma_w^2 = w^T Sigma w
   l = 1:  S = S0 u^T Sigma v
           T = S0/(2m) [u.v - 2 v^T Sigma B' u - 2 u^T Sigma A v + 12 (u^T Sigma v) tr(A B' Sigma) + 4 u^T Sigma (A B' + B' A) Sigma v]
           <(u.x_z)(v.x_z)/|y|> = S0 sqrt(2/pi)/sigma [u^T Sigma v - g_u g_v / (3 sigma^2)],  g = Sigma w
   dipolo  <0| sum z_i |1> = S0 1^T Sigma v;   <0| sum d/dz_i |1> = S0 [1.v - 2 v^T Sigma B' 1]
(derivate a mano con il teorema di Isserlis; controllate in tav_dev.py contro quadrature Monte Carlo e contro sistemi
non interagenti esatti.)"""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time, itertools
import numpy as np
import scipy.linalg as sla
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SQ2PI = np.sqrt(2 / np.pi)
EH_EV = 27.211386245988
# masse della regola (letture del thread principale, multi_e_v5.json, checks R3v5 e tables T6): ingressi imparati
M_T_THETA_A = 1.0003732690576126
M0 = 1.0000000045603115


# =============================================================================================== spin
def spin_vec(kind):
    """funzione di spin accoppiata come tensore (2,)*N.
    'd3': doppietto  [ab - ba]/sqrt2 x a   (Li, Be+)
    's2': singoletto [ab - ba]/sqrt2        (Li+, He-like)
    's4': singoletto x singoletto           (Be)
    'up1': una pallina su                    (idrogenoidi)
    'p3', 'p4': funzioni primitive a M_S fissato (aab, aabb): danno il minimo del settore M_S, che per Li 2S, Li 2P e
                Be 1S e' lo stato voluto (gli stati di spin piu' alto stanno molto piu' su); per Be 1P serve 's4'"""
    up, dn = np.array([1., 0.]), np.array([0., 1.])
    th = (np.kron(up, dn) - np.kron(dn, up)) / np.sqrt(2)
    if kind == "up1":
        v = up
    elif kind == "s2":
        v = th
    elif kind == "d3":
        v = np.kron(th, up)
    elif kind == "s4":
        v = np.kron(th, th)
    elif kind == "p3":                                  # primitiva a,a,b (M_S = 1/2): 2 permutazioni
        v = np.kron(np.kron(up, up), dn)
    elif kind == "p4":                                  # primitiva a,a,b,b (M_S = 0): 4 permutazioni
        v = np.kron(np.kron(up, up), np.kron(dn, dn))
    elif kind == "t3":                                  # quartetto aaa (solo per i controlli)
        v = np.kron(np.kron(up, up), up)
    elif kind == "b3":                                  # 'bosoni' (controllo D3: senza antisimmetria)
        v = np.kron(th, up)
    else:
        raise ValueError(kind)
    N = int(round(np.log2(len(v))))
    return v.reshape((2,) * N), N


def perm_weights(kind, bosonic=False):
    """[(c_P, Pm)] con Pm = eye[p] ((Pm x)_i = x_p(i)); c_P = sign(P) <chi|P chi> (sign = 1 se bosonic)"""
    chi, N = spin_vec(kind)
    out = []
    for p in itertools.permutations(range(N)):
        Pm = np.eye(N)[list(p)]
        if bosonic:                                     # simmetrizzatore totale (nessuna antisimmetria, controllo)
            out.append((1.0, Pm)); continue
        sg = float(round(np.linalg.det(Pm)))
        c = sg * float(np.sum(chi * np.transpose(chi, np.argsort(p))))   # (P chi)(s) = chi(s o p), come lo spazio
        if abs(c) > 1e-12:
            out.append((c, Pm))
    return out


# =============================================================================================== operatore
class System:
    def __init__(self, Z, spin="d3", L=0, m_kin=M_T_THETA_A, ee=True, bosonic=False, omega=0.0, name=""):
        self.Z, self.L, self.m, self.ee, self.omega, self.name = float(Z), int(L), m_kin, ee, omega, name
        self.spin = spin
        self.perms = perm_weights(spin, bosonic)
        self.N = self.perms[0][1].shape[0]
        self.pairs = [(i, j) for i in range(self.N) for j in range(i + 1, self.N)]
        W = [np.eye(self.N)[i] for i in range(self.N)]
        self.wn = np.array(W)                                     # nucleo-elettrone
        self.wp = np.array([W[i] - W[j] for i, j in self.pairs]) if self.pairs else np.zeros((0, self.N))

    def me(self, A, u, B, v, parts=False):
        """A, B (n,N,N), u, v (n,N) (ignorati se L = 0) -> S, H (n,) sommati sulle permutazioni pesate"""
        n, N = A.shape[0], self.N
        S_t = np.zeros(n); H_t = np.zeros(n)
        if parts:
            P = {k: np.zeros(n) for k in ("T", "Ven", "Vee", "r2", "delta", "P4")}
        for c, Pm in self.perms:
            Bp = np.einsum("ji,njk,kl->nil", Pm, B, Pm)
            C = A + Bp
            Ci = np.linalg.inv(C)
            S0 = (np.pi ** N / np.linalg.det(C)) ** 1.5
            Sig = 0.5 * Ci
            if self.L == 0:
                S = S0
                T = S0 * (3 / self.m) * np.einsum("nij,nji->n", A @ Ci, Bp)
                sn2 = np.einsum("ki,nij,kj->nk", self.wn, Sig, self.wn)
                Ven = -self.Z * S0 * SQ2PI * np.sum(sn2 ** -0.5, axis=1)
                Vee = S0 * SQ2PI * np.sum(np.einsum("ki,nij,kj->nk", self.wp, Sig, self.wp) ** -0.5, axis=1) if (self.ee and self.pairs) else 0 * S0
                Vh = 0.5 * self.omega ** 2 * S0 * 3 * np.trace(Sig, axis1=1, axis2=2) if self.omega else 0
            else:
                vp = np.einsum("ji,nj->ni", Pm, v)
                Sv = np.einsum("nij,nj->ni", Sig, vp)              # Sigma v
                Su = np.einsum("nij,nj->ni", Sig, u)               # Sigma u
                uSv = np.einsum("ni,ni->n", u, Sv)
                S = S0 * uSv
                AB = A @ Bp
                trABS = np.einsum("nij,nji->n", AB, Sig)
                T = S0 / (2 * self.m) * (np.einsum("ni,ni->n", u, vp)
                                         - 2 * np.einsum("ni,nij,nj->n", Sv, Bp, u)
                                         - 2 * np.einsum("ni,nij,nj->n", Su, A, vp)
                                         + 12 * uSv * trABS
                                         + 4 * np.einsum("ni,nij,nj->n", Su, AB + Bp @ A, Sv))
                def pot(W):
                    s2 = np.einsum("ki,nij,kj->nk", W, Sig, W)            # (n, k)
                    gu = np.einsum("ki,ni->nk", W, Su)                   # u^T Sigma w
                    gv = np.einsum("ki,ni->nk", W, Sv)
                    return S0 * SQ2PI * np.sum(s2 ** -0.5 * (uSv[:, None] - gu * gv / (3 * s2)), axis=1)
                Ven = -self.Z * pot(self.wn)
                Vee = pot(self.wp) if (self.ee and self.pairs) else 0 * S0
                Vh = 0.5 * self.omega ** 2 * S0 * (3 * uSv * np.trace(Sig, axis1=1, axis2=2) + 2 * np.einsum("ni,ni->n", Su, Sv)) if self.omega else 0
            S_t += c * S; H_t += c * (T + Ven + Vee + Vh)
            if parts:
                P["T"] += c * T; P["Ven"] += c * Ven; P["Vee"] += c * Vee
                if self.L == 0:
                    dg = np.diagonal(Sig, axis1=1, axis2=2)
                    P["r2"] += c * S0 * 3 * np.sum(dg, axis=1)
                    P["delta"] += c * S0 * np.sum((2 * np.pi * dg) ** -1.5, axis=1)
                    s4 = 2 * np.diagonal(A @ Ci @ Bp, axis1=1, axis2=2)    # varianza d'impulso per componente
                    P["P4"] += c * S0 * np.sum(15 * s4 ** 2, axis=1)     # solo approssimato (vedi expect_parts)
        if parts:
            return S_t, H_t, P
        return S_t, H_t


def cross_weights(kind_bra, kind_ket):
    """c_P = sign(P) <chi_bra | P chi_ket> (stati con funzioni di spin diverse, stesso S e M_S: es. Be 1S aabb e 1P s4)"""
    cb, N = spin_vec(kind_bra)
    ck, _ = spin_vec(kind_ket)
    out = []
    for p in itertools.permutations(range(N)):
        Pm = np.eye(N)[list(p)]
        c = float(round(np.linalg.det(Pm))) * float(np.sum(cb * np.transpose(ck, np.argsort(p))))
        if abs(c) > 1e-12:
            out.append((c, Pm))
    return out


def dipole_me(sysS, sysP, A, B, v):
    """<phi_A (S, l=0) | sum_i z_i | P phi_B (l=1)> e <..| sum_i d/dz_i |..> sommati con i pesi di spin incrociati
    c_P = sign(P) <chi_S | P chi_P>. A (n,N,N), B (n,N,N), v (n,N)"""
    n, N = A.shape[0], sysS.N
    one = np.ones(N)
    D_len = np.zeros(n); D_vel = np.zeros(n)
    for c, Pm in cross_weights(sysS.spin, sysP.spin):
        Bp = np.einsum("ji,njk,kl->nil", Pm, B, Pm)
        vp = np.einsum("ji,nj->ni", Pm, v)
        C = A + Bp
        Ci = np.linalg.inv(C)
        S0 = (np.pi ** N / np.linalg.det(C)) ** 1.5
        Sv = 0.5 * np.einsum("nij,nj->ni", Ci, vp)
        D_len += c * S0 * (Sv @ one)
        D_vel += c * S0 * (vp @ one - 2 * np.einsum("ni,nij,j->n", Sv, Bp, one))
    return D_len, D_vel


# =============================================================================================== funzioni a caso
def build_A(N, a_en, a_ee):
    A = np.diag(a_en).astype(float)
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            A[i, i] += a_ee[k]; A[j, j] += a_ee[k]; A[i, j] -= a_ee[k]; A[j, i] -= a_ee[k]
            k += 1
    return A


def random_fn(rng, N, Zs, Lmin=0.01, Lmax=30.0, L=0, p_pair=0.8):
    """stessa famiglia del thread principale (random_A di multi_e_ecg.py): lunghezze log-uniformi in [Lmin, Lmax]/Zs
    per elettrone-nucleo, coppie presenti con prob. 0,8 (segno + 4 volte su 5); per L = 1 un vettore u a caso
    (meta' delle volte un solo elettrone, meta' gaussiano)"""
    while True:
        Ls = np.exp(rng.uniform(np.log(Lmin / Zs), np.log(Lmax / Zs), N))
        a_en = 1 / Ls ** 2
        npair = N * (N - 1) // 2
        a_ee = np.zeros(npair)
        for k in range(npair):
            if rng.random() > p_pair:
                continue
            L12 = np.exp(rng.uniform(np.log(0.05 / Zs), np.log(Lmax / Zs)))
            a_ee[k] = (1 if rng.random() < 0.8 else -1) / L12 ** 2
        A = build_A(N, a_en, a_ee)
        if np.linalg.eigvalsh(A).min() > 1e-12 * np.abs(A).max():
            break
    if L == 0:
        return A, np.zeros(N)
    if rng.random() < 0.5:
        u = np.zeros(N); u[rng.integers(N)] = 1.0
    else:
        u = rng.normal(size=N); u /= np.linalg.norm(u)
    return A, u


def mutate(rng, A, u, Zs, Lmax, L, sig=0.4, p_pair=0.3, chol=False):
    """candidato 'figlio' di una funzione gia' nella base: A' = D A D (D diagonale, scale log-normali: resta definita
    positiva) e, con probabilita' p_pair, un termine di coppia nuovo; u' = u + rumore. (Aggiunta di questo studio:
    con le sole funzioni a caso il litio converge male, vedi tav_dev / documento par. 2.)"""
    N = A.shape[0]
    if chol:                                            # mossa generale: rumore sul fattore di Cholesky
        Lc = np.linalg.cholesky(A)
        sc = np.sqrt(np.abs(np.diag(Lc)))[:, None] * np.sqrt(np.abs(np.diag(Lc)))[None, :]
        Lc2 = np.tril(Lc * np.exp(sig * rng.normal(size=(N, N))) + 0.5 * sig * sc * np.tril(rng.normal(size=(N, N)), -1))
        A2 = Lc2 @ Lc2.T
        if np.linalg.eigvalsh(A2).min() <= 1e-12 * np.abs(A2).max():
            A2 = A.copy()
    else:
        d = np.exp(sig * rng.normal(size=N))
        A2 = A * np.outer(d, d)
    if rng.random() < p_pair and N > 1:
        i, j = rng.choice(N, 2, replace=False)
        L12 = np.exp(rng.uniform(np.log(0.05 / Zs), np.log(Lmax / Zs)))
        a = 1 / L12 ** 2
        e = np.zeros(N); e[i] = 1; e[j] = -1
        A2 = A2 + a * np.outer(e, e)
    if L == 1:
        u2 = u + 0.3 * rng.normal(size=N); u2 /= np.linalg.norm(u2)
    else:
        u2 = u
    return A2, u2


# =============================================================================================== SVM
class SVM:
    """crescita stocastica (una funzione alla volta, la migliore fra 'trials' a caso, aggiornamento di rango 1 con
    l'equazione secolare) + giri di raffinamento; stessa logica di multi_e_ecg.SVM, riscritta per (A, u)."""

    def __init__(self, sys_, rng, Zs, Lmax=30.0, trials=30, mut=0.5, p_pair=0.8, polish=(), npol=10, chol=False, nm=0):
        self.s, self.rng, self.Zs, self.Lmax, self.trials, self.mut = sys_, rng, Zs, Lmax, trials, mut
        self.p_pair, self.polish, self.npol, self.chol, self.nm = p_pair, tuple(polish), npol, chol, nm
        N = sys_.N
        self.A = np.zeros((0, N, N)); self.u = np.zeros((0, N))
        self.nset = np.zeros(0)
        self.Smat = np.zeros((0, 0)); self.Hmat = np.zeros((0, 0))
        self.E = None

    def _eig(self, Hm, Sm):
        try:
            return sla.eigh(Hm, Sm)
        except (np.linalg.LinAlgError, ValueError):
            w, V = np.linalg.eigh(Sm)
            keep = w > 1e-13 * w.max()
            X = V[:, keep] / np.sqrt(w[keep])
            Ev, Y = np.linalg.eigh(X.T @ Hm @ X)
            return Ev, X @ Y

    @staticmethod
    def _lowest_with(Ev, U, s, h, hd):
        if U is None or len(Ev) == 0:
            return hd
        us = U.T @ s; uh = U.T @ h
        nn = 1.0 - us @ us
        if nn < 1e-10:
            return None
        c = (uh - Ev * us) / np.sqrt(nn)
        d = (hd - 2 * us @ uh + np.sum(Ev * us ** 2)) / nn
        f = lambda lam: d - lam - np.sum(c ** 2 / (Ev - lam))
        hi = Ev[0] - 1e-14 * max(1, abs(Ev[0]))
        if f(hi) > 0:
            return Ev[0]
        lo = min(d, Ev[0]) - np.sqrt(np.sum(c ** 2)) - 1.0
        while f(lo) < 0:
            lo -= 2 * (abs(lo) + 1)
        return brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=300)

    def _rows(self, Ac, uc, Aset, uset):
        m, K = len(Ac), len(Aset)
        if K:
            S, H = self.s.me(np.repeat(Ac, K, 0), np.repeat(uc, K, 0), np.tile(Aset, (m, 1, 1)), np.tile(uset, (m, 1)))
            S = S.reshape(m, K); H = H.reshape(m, K)
        else:
            S = H = np.zeros((m, 0))
        Sd, Hd = self.s.me(Ac, uc, Ac, uc)
        return S, H, Sd, Hd

    def _best(self, Aset, uset, Smat, Hmat, nset, extra=None):
        if len(Aset):
            Ev, U = self._eig(Hmat, Smat)
        else:
            Ev, U = np.zeros(0), None
        nmut = int(self.mut * self.trials) if len(Aset) else 0
        cs = [random_fn(self.rng, self.s.N, self.Zs, Lmax=self.Lmax, L=self.s.L, p_pair=self.p_pair) for _ in range(self.trials - nmut)]
        cs += [mutate(self.rng, Aset[k], uset[k], self.Zs, self.Lmax, self.s.L)
               for k in self.rng.integers(len(Aset), size=nmut)] if nmut else []
        Ac = np.array([c[0] for c in cs]); uc = np.array([c[1] for c in cs])
        if extra is not None:
            Ac = np.concatenate([Ac, extra[0][None]]); uc = np.concatenate([uc, extra[1][None]])
        best = self._score(Ac, uc, Aset, uset, nset, Ev, U, extra is not None)
        # lucidatura: qualche giro di figli del vincitore (ricerca locale a caso, sigma che cala)
        for sg in self.polish:
            if best is None:
                break
            ch = [mutate(self.rng, best[1], best[2], self.Zs, self.Lmax, self.s.L, sig=sg, p_pair=0.0,
                         chol=self.chol and self.rng.random() < 0.5)
                  for _ in range(self.npol)]
            b2 = self._score(np.array([c[0] for c in ch]), np.array([c[1] for c in ch]), Aset, uset, nset, Ev, U, False)
            if b2 is not None and b2[0] < best[0]:
                best = b2
        if self.nm and best is not None:
            best = self._nelder_mead(best, Aset, uset, nset, Ev, U)
        return best

    def _nelder_mead(self, best, Aset, uset, nset, Ev, U):
        """ottimizzazione locale (Nelder-Mead, maxfev = self.nm) dei parametri non lineari del vincitore:
        A = L L^T, L triangolare inferiore con diagonale exp(.), piu' u (per L = 1)"""
        from scipy.optimize import minimize
        N = self.s.N
        Lc = np.linalg.cholesky(best[1])
        il = np.tril_indices(N)
        d = np.sqrt(np.diag(best[1]))
        x0 = []
        for i, j in zip(*il):
            x0.append(np.log(Lc[i, j]) if i == j else Lc[i, j] / d[i])
        if self.s.L == 1:
            x0 += list(best[2])
        x0 = np.array(x0)
        npar = len(il[0])

        def unpack(x):
            L = np.zeros((N, N))
            for k, (i, j) in enumerate(zip(*il)):
                L[i, j] = np.exp(x[k]) if i == j else x[k] * d[i]
            u = x[npar:] if self.s.L == 1 else np.zeros(N)
            return L @ L.T, u
        cache = {}

        def f(x):
            A, u = unpack(x)
            if not np.all(np.isfinite(A)) or np.linalg.eigvalsh(A).min() <= 1e-12 * np.abs(A).max():
                return 1e10
            if self.s.L == 1 and np.linalg.norm(u) < 1e-6:
                return 1e10
            b = self._score(A[None], u[None], Aset, uset, nset, Ev, U, False)
            if b is None:
                return 1e10
            cache[b[0]] = b
            return b[0]
        sim = [x0] + [x0 + 0.25 * np.eye(len(x0))[k] for k in range(len(x0))]
        r = minimize(f, x0, method="Nelder-Mead", options=dict(maxfev=self.nm, initial_simplex=np.array(sim),
                                                                xatol=1e-4, fatol=1e-10))
        bb = min(cache) if cache else None
        if bb is not None and bb < best[0]:
            return cache[bb]
        return best

    def _score(self, Ac, uc, Aset, uset, nset, Ev, U, last_is_old):
        S, H, Sd, Hd = self._rows(Ac, uc, Aset, uset)
        best = None
        for i in range(len(Ac)):
            if not (Sd[i] > 0):
                continue
            nr = 1 / np.sqrt(Sd[i])
            s = S[i] * nr * nset; h = H[i] * nr * nset
            e = self._lowest_with(Ev, U, s, h, Hd[i] * nr ** 2)
            if e is not None and (best is None or e < best[0]):
                best = (e, Ac[i], uc[i], s, h, Hd[i] * nr ** 2, last_is_old and i == len(Ac) - 1, nr)
        return best

    def add(self):
        K = len(self.A)
        b = self._best(self.A, self.u, self.Smat, self.Hmat, self.nset)
        e, A, u, s, h, hd, _, nr = b
        Sm = np.zeros((K + 1, K + 1)); Hm = np.zeros((K + 1, K + 1))
        Sm[:K, :K] = self.Smat; Hm[:K, :K] = self.Hmat
        Sm[K, :K] = Sm[:K, K] = s; Hm[K, :K] = Hm[:K, K] = h
        Sm[K, K] = 1.0; Hm[K, K] = hd
        self.A = np.concatenate([self.A, A[None]]); self.u = np.concatenate([self.u, u[None]])
        self.nset = np.concatenate([self.nset, [nr]])
        self.Smat, self.Hmat = Sm, Hm
        self.E = float(self._eig(Hm, Sm)[0][0])
        return self.E

    def refine_one(self, k):
        K = len(self.A)
        keep = np.array([i for i in range(K) if i != k])
        Sm = self.Smat[np.ix_(keep, keep)]; Hm = self.Hmat[np.ix_(keep, keep)]
        b = self._best(self.A[keep], self.u[keep], Sm, Hm, self.nset[keep], extra=(self.A[k], self.u[k]))
        e, A, u, s, h, hd, is_old, nr = b
        if is_old:
            return
        self.A[k] = A; self.u[k] = u; self.nset[k] = nr
        self.Smat[k, keep] = self.Smat[keep, k] = s; self.Hmat[k, keep] = self.Hmat[keep, k] = h
        self.Smat[k, k] = 1.0; self.Hmat[k, k] = hd

    def rebuild(self, sys_new=None, chunk=20000):
        if sys_new is not None:
            self.s = sys_new
        K = len(self.A)
        I, J = np.triu_indices(K)
        Sv = np.zeros(len(I)); Hv = np.zeros(len(I))
        for i0 in range(0, len(I), chunk):
            sl = slice(i0, i0 + chunk)
            Sv[sl], Hv[sl] = self.s.me(self.A[I[sl]], self.u[I[sl]], self.A[J[sl]], self.u[J[sl]])
        Sd = np.zeros((K, K)); Hd = np.zeros((K, K))
        Sd[I, J] = Sv; Sd[J, I] = Sv; Hd[I, J] = Hv; Hd[J, I] = Hv
        nr = 1 / np.sqrt(np.diag(Sd))
        self.nset = nr
        self.Smat = Sd * np.outer(nr, nr); self.Hmat = Hd * np.outer(nr, nr)
        self.E = float(self._eig(self.Hmat, self.Smat)[0][0])
        return self.E

    def energy_subset(self, K):
        return float(self._eig(self.Hmat[:K, :K], self.Smat[:K, :K])[0][0])

    def state(self):
        """energia e coefficienti delle funzioni NON normalizzate (A, u) dello stato piu' basso, norma <psi|psi> = 1"""
        Ev, U = self._eig(self.Hmat, self.Smat)
        c = U[:, 0] / np.sqrt(U[:, 0] @ self.Smat @ U[:, 0])
        return float(Ev[0]), c * self.nset

    def expect(self, chunk=20000):
        """<T>, <Ven>, <Vee>, <sum r_i^2>, <sum delta(r_i)> sullo stato piu' basso (L = 0 per r2/delta)"""
        K = len(self.A)
        E, c = self.state()
        I, J = np.triu_indices(K)
        w = c[I] * c[J] * np.where(I == J, 1.0, 2.0)
        acc = {}
        nrm = 0.0; Etot = 0.0
        for i0 in range(0, len(I), chunk):
            sl = slice(i0, i0 + chunk)
            S, H, P = self.s.me(self.A[I[sl]], self.u[I[sl]], self.A[J[sl]], self.u[J[sl]], parts=True)
            nrm += np.sum(w[sl] * S); Etot += np.sum(w[sl] * H)
            for k, v in P.items():
                acc[k] = acc.get(k, 0.0) + float(np.sum(w[sl] * v))
        out = {k: v / nrm for k, v in acc.items()}
        out["E"] = Etot / nrm; out["norm"] = nrm
        out["virial_-V/T"] = -(out["Ven"] + out["Vee"]) / out["T"]
        return out


# =============================================================================================== funzione d'onda
class Wave:
    """psi(x) = sum_k c_k sum_P c_P (u_k^T (P x)_z)^l exp(-(Px)^T A_k (Px)) (spazio; con il peso di spin gia' dentro).
    Per il campionamento di |psi|^2 si usa la densita' spaziale dello stato con spin: rho(x) = sum_P c_P psi_0(x) psi_0(Px)
    non e' positiva punto per punto; quindi si campiona la parte spaziale per un'etichetta di spin fissata
    (componente up,up,down / up,down,up,down), che e' una funzione d'onda antisimmetrica ordinaria nelle coordinate
    con lo stesso spin: Psi_sigma(x) = sum_P sign(P) chi(P sigma) phi(P x) ecc. (vedi spatial_for_spin)."""

    def __init__(self, sys_, A, u, c):
        self.s, self.A, self.u, self.c = sys_, A, u, c
        chi, N = spin_vec(sys_.spin)
        self.chi = chi
        self.N = N

    def spatial_for_spin(self, sigma):
        """lista di (peso, permutazione) per Psi(x, sigma) = sum_P sign(P) chi(sigma o P) phi(P x)"""
        out = []
        for p in itertools.permutations(range(self.N)):
            Pm = np.eye(self.N)[list(p)]
            sg = float(round(np.linalg.det(Pm)))
            # (P Psi)(x, sigma) = Psi(x_p, sigma_p): chi valutata sulle etichette permutate
            w = sg * self.chi[tuple(sigma[i] for i in p)]
            if abs(w) > 1e-12:
                out.append((w, np.array(p)))
        return out

    def value(self, x, terms):
        """x (n, N, 3) -> Psi (n,)"""
        tot = np.zeros(len(x))
        for w, p in terms:
            xp = x[:, p, :]
            q = np.einsum("nix,kij,njx->nk", xp, self.A, xp)
            f = np.exp(-q)
            if self.s.L == 1:
                f = f * np.einsum("ki,ni->nk", self.u, xp[:, :, 2])
            tot += w * (f @ self.c)
        return tot


def metropolis(wave, sigma, n_out, rng, nw=2000, step=None, burn=300, every=10, scale=1.0, log=print):
    """campioni (n_out, N, 3) float32 da |Psi(x, sigma)|^2; passo adattato all'accettazione ~0,5 per elettrone"""
    terms = wave.spatial_for_spin(sigma)
    N = wave.N
    x = rng.normal(0, scale, (nw, N, 3))
    p = wave.value(x, terms) ** 2
    st = np.full(N, 0.3 * scale) if step is None else np.full(N, step)
    out = []; acc = np.zeros(N); tot = np.zeros(N)
    it = 0
    while sum(len(o) for o in out) < n_out:
        for i in range(N):                                      # una pallina alla volta
            xn = x.copy(); xn[:, i] += rng.normal(0, 1, (nw, 3)) * st[i]
            pn = wave.value(xn, terms) ** 2
            a = rng.random(nw) < pn / np.maximum(p, 1e-300)
            x[a] = xn[a]; p[a] = pn[a]; acc[i] += a.sum(); tot[i] += nw
        it += 1
        if it < burn and it % 20 == 0:
            r = acc / np.maximum(tot, 1)
            st *= np.exp(r - 0.5)
            acc[:] = 0; tot[:] = 0
        if it >= burn and it % every == 0:
            out.append(x.copy())
    S = np.concatenate(out)[:n_out]
    return S.astype(np.float32), dict(accept=(acc / np.maximum(tot, 1)).tolist(), step=st.tolist(), iters=it)


def save_json(path, obj):
    tmp = path + ".tmp"
    json.dump(obj, open(tmp, "w"), indent=1, default=float, ensure_ascii=False)
    os.replace(tmp, path)


def load_json(path):
    return json.load(open(path)) if os.path.exists(path) else {}
