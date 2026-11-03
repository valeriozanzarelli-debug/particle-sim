"""tavg_ci2.py - FASE 3: due elettroni nel blocco con interazione di configurazioni a onde parziali (tutti i tipi di
coppia: 2s2 1S, 2s2p 3P/1P, 2p2 3P/1D/1S), piu' l'Hartree-Fock della stessa configurazione (stesso termine).

Base radiale: per ogni l, n_rad funzioni di Laguerre x^(l+1) e^(-x/2) L_k^(2l+2)(x), x = 2 lam r, lam = sqrt(-2 eps_s), sulla
griglia del livello fine (r <= R_box), ortogonalizzate (canonica) e diagonalizzate con h_l del blocco (cinetico + Vloc
+ l(l+1)/2r^2 + proiettore KB per l = 0 con il chi VERO del blocco): la parte a un elettrone e' diagonale. Repulsione: sviluppo in multipoli,
  <l1 l2 L|C^k.C^k|l3 l4 L> = (-1)^(l2+l3+L) {l1 l2 L; l4 l3 k} <l1||C^k||l3><l2||C^k||l4>,
integrali radiali R^k con la mappa Lk della griglia (tavg_radial). Spin: spazio (anti)simmetrico (1 + s P12),
P12 |a l1, b l2; L> = (-1)^(l1+l2-L) |b l2, a l1; L>.
HF: una sola configurazione (n l)(n' l') col termine dato, orbitali radiali ottimizzati nella stessa base (iterazione
di Fock col tensore della CI, smorzata). Correlazione di coppia = E_CI(l -> inf) - E_HF.
Coda in l: incrementi dE_l ~ (l + 1/2)^-p con p dal rapporto degli ultimi due incrementi (limitato a [3, 8]),
somma della coda analitica (regola fissata nella pre-registrazione).
Uso: python tavg_ci2.py <Z> <stato> [lmax n_rad R_box]   stati: s2_1S sp_3P sp_1P pp_3P pp_1D pp_1S"""
import os, sys, json, math, time
from fractions import Fraction
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import scipy.linalg as sla
from scipy.special import eval_genlaguerre
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_pecg as P
import tavg_radial as tr

OUT = os.path.join(HERE, "out")
# state: (L, S, parity, reference configuration (l1, l2), root index of the reference in its symmetry)
LAM = float(os.environ.get("TAVG_LAM", 1.0))   # Laguerre scale / sqrt(-2 eps_s): numerical, convergence tested
STATES = {"s2_1S": (0, 0, 0, (0, 0)), "sp_3P": (1, 1, 1, (0, 1)), "sp_1P": (1, 0, 1, (0, 1)),
          "pp_3P": (1, 1, 0, (1, 1)), "pp_1D": (2, 0, 0, (1, 1)), "pp_1S": (0, 0, 0, (1, 1))}


# ------------------------------------------------------------------ angular algebra (exact rationals)
def _f(n):
    return math.factorial(n)


def _delta(a, b, c):
    return Fraction(_f(a + b - c) * _f(a - b + c) * _f(-a + b + c), _f(a + b + c + 1))


def w3j0(a, b, c):
    """(a b c; 0 0 0)"""
    J = a + b + c
    if J % 2 or c < abs(a - b) or c > a + b:
        return 0.0
    g = J // 2
    v = (-1) ** g * math.sqrt(_delta(a, b, c)) * _f(g) / (_f(g - a) * _f(g - b) * _f(g - c))
    return v


def w6j(a, b, c, d, e, f):
    """{a b c; d e f} (integers), Racah formula"""
    for (x, y, z) in ((a, b, c), (a, e, f), (d, b, f), (d, e, c)):
        if z < abs(x - y) or z > x + y:
            return 0.0
    pre = math.sqrt(_delta(a, b, c) * _delta(a, e, f) * _delta(d, b, f) * _delta(d, e, c))
    s = 0
    for t in range(max(a + b + c, a + e + f, d + b + f, d + e + c), min(a + b + d + e, a + c + d + f, b + c + e + f) + 1):
        s += Fraction((-1) ** t * _f(t + 1), _f(t - a - b - c) * _f(t - a - e - f) * _f(t - d - b - f) * _f(t - d - e - c)
                      * _f(a + b + d + e - t) * _f(a + c + d + f - t) * _f(b + c + e + f - t))
    return pre * float(s)


def redC(l, k, lp):
    return (-1) ** l * math.sqrt((2 * l + 1) * (2 * lp + 1)) * w3j0(l, k, lp)


# ------------------------------------------------------------------ radial basis
class PairCI:
    def __init__(self, field, lmax=6, n_rad=24, R_box=12.0, log=print):
        self.f, self.lmax, self.n, self.Rb, self.log = field, lmax, n_rad, R_box, log
        g = field.g; r = g.r
        self.g = g
        self.m = int(np.searchsorted(r, R_box))
        self.w = r * g.h                                # int f dr = sum f w
        atom = field.radial_atom([(1, 0, 1)])
        zz = float(getattr(field, "Zv", getattr(field, "Z", 1.0)))
        self.eps, self.Pb = {}, {}
        m = self.m
        rr = r[:m]
        H0 = 0.5 * (atom.h_matrix(0)[:m, :m] + atom.h_matrix(0)[:m, :m].T)
        e0 = float(tr.geig_low(H0, r[:m] ** 2, 1, -0.6 * zz ** 2 - 20)[0][0])
        self.lam = LAM * math.sqrt(-2 * e0)               # Laguerre scale from the lowest s level of the block
        for l in range(lmax + 1):
            H = atom.h_matrix(l)[:m, :m]
            H = 0.5 * (H + H.T)
            x = 2 * self.lam * rr
            Pst = np.stack([x ** (l + 1) * np.exp(-x / 2) * eval_genlaguerre(k, 2 * l + 2, x) for k in range(n_rad)], 1)
            Pst /= np.sqrt(np.sum(Pst ** 2 * self.w[:m, None], 0))
            Ust = Pst / np.sqrt(rr)[:, None]
            Hb = Ust.T @ H @ Ust; Sb = Ust.T @ (Ust * (rr ** 2)[:, None])
            sw, sv = np.linalg.eigh(Sb)
            keep = sw > 1e-13 * sw.max()
            X = sv[:, keep] / np.sqrt(sw[keep])
            w, Y = np.linalg.eigh(X.T @ Hb @ X)
            Pm = Pst @ (X @ Y)
            nr = np.sqrt(np.sum(Pm ** 2 * self.w[:m, None], 0))
            self.eps[l] = w; self.Pb[l] = (Pm / nr).T    # (n, m), n = kept functions
        self.n = min(len(self.eps[l]) for l in self.eps)
        for l in self.eps:
            self.eps[l] = self.eps[l][:self.n]; self.Pb[l] = self.Pb[l][:self.n]
        self.Lk = {}
        self._blk = {}

    def _L(self, k):
        if k not in self.Lk:
            self.Lk[k] = self.g.Lk(k)[:self.m, :self.m]
        return self.Lk[k]

    def Rk(self, k, l1, l3, l2, l4):
        """R^k[a,b,c,d] = int int P_a P_c (1) r<^k/r>^(k+1) P_b P_d (2), a in l1, c in l3, b in l2, d in l4"""
        n, w = self.n, self.w[:self.m]
        D13 = (self.Pb[l1][:, None, :] * self.Pb[l3][None, :, :]).reshape(n * n, -1)
        D24 = (self.Pb[l2][:, None, :] * self.Pb[l4][None, :, :]).reshape(n * n, -1)
        V = D24 @ self._L(k).T                           # (bd, i): (Lk @ D_bd)_i
        R = (D13 * w) @ V.T                              # (ac, bd)
        return R.reshape(n, n, n, n).transpose(0, 2, 1, 3)   # a b c d

    def direct(self, l1, l2, l3, l4, L):
        """<a l1, b l2; L| H |c l3, d l4; L> (no symmetrization), as (n, n, n, n) [a, b, c, d]"""
        key = (l1, l2, l3, l4, L)
        if key in self._blk:
            return self._blk[key]
        n = self.n
        M = np.zeros((n, n, n, n))
        if (l1, l2) == (l3, l4):
            e = self.eps[l1][:, None] + self.eps[l2][None, :]
            idx = np.arange(n)
            M[idx[:, None], idx[None, :], idx[:, None], idx[None, :]] = e
        for k in range(0, 2 * self.lmax + 1):
            c1, c2 = redC(l1, k, l3), redC(l2, k, l4)
            if c1 == 0 or c2 == 0:
                continue
            ang = (-1) ** (l2 + l3 + L) * w6j(l1, l2, L, l4, l3, k) * c1 * c2
            if ang != 0:
                M += ang * self.Rk(k, l1, l3, l2, l4)
        self._blk[key] = M
        return M

    def channels(self, L, par):
        return [(l1, l2) for l1 in range(self.lmax + 1) for l2 in range(l1, self.lmax + 1)
                if abs(l1 - l2) <= L <= l1 + l2 and (l1 + l2) % 2 == par]

    def T(self, l1, l2, l3, l4, L, s):
        """<a l1 b l2| H (1 + s P12) |c l3 d l4>"""
        X = self.direct(l1, l2, l4, l3, L).transpose(0, 1, 3, 2)     # <ab|H|d c>_(l4 l3)
        return self.direct(l1, l2, l3, l4, L) + s * (-1) ** (l3 + l4 - L) * X

    def solve(self, L, S, par, lmax=None, nroot=1):
        s = 1 if S == 0 else -1
        n = self.n
        chs = [c for c in self.channels(L, par) if lmax is None or c[1] <= lmax]
        basis = []                                        # per channel: list of (a, b), norms
        for (l1, l2) in chs:
            if l1 == l2:
                prs = [(a, b) for a in range(n) for b in range(a, n)]
            else:
                prs = [(a, b) for a in range(n) for b in range(n)]
            nr = []
            for a, b in prs:
                pp = (-1) ** (l1 + l2 - L) if (l1 == l2 and a == b) else 0
                nr.append(2 * (1 + s * pp))
            nr = np.array(nr, float)
            keep = nr > 1e-12
            basis.append((l1, l2, np.array(prs)[keep], nr[keep]))
        dims = [len(b[2]) for b in basis]
        off = np.concatenate([[0], np.cumsum(dims)])
        H = np.zeros((off[-1], off[-1]))
        for I, (l1, l2, p1, n1) in enumerate(basis):
            for J in range(I, len(basis)):
                l3, l4, p2, n2 = basis[J]
                Tm = self.T(l1, l2, l3, l4, L, s)
                blk = 2 * Tm[p1[:, 0][:, None], p1[:, 1][:, None], p2[:, 0][None, :], p2[:, 1][None, :]]
                blk = blk / np.sqrt(np.outer(n1, n2))
                H[off[I]:off[I + 1], off[J]:off[J + 1]] = blk
                H[off[J]:off[J + 1], off[I]:off[I + 1]] = blk.T
        H = 0.5 * (H + H.T)
        w = sla.eigh(H, eigvals_only=True, subset_by_index=[0, nroot - 1])
        return w, off[-1]

    def hf(self, L, S, l1, l2, maxit=400, tol=1e-12):
        """single configuration (l1)(l2) of term (L, S), radial orbitals optimized in the basis"""
        s = 1 if S == 0 else -1
        Tm = self.T(l1, l2, l1, l2, L, s)
        n = self.n
        c = np.zeros(n); c[0] = 1.0
        d = np.zeros(n); d[0] = 1.0
        same = (l1 == l2)
        E_old = None
        for it in range(maxit):
            pp = (-1) ** (l1 + l2 - L) if same else 0
            if same:
                F = np.einsum("abcd,b,d->ac", Tm, c, c)
                F = 0.5 * (F + F.T)
                w, U = np.linalg.eigh(F)
                cn = U[:, 0] * np.sign(U[0, 0] if abs(U[0, 0]) > 1e-8 else 1)
                c = 0.5 * c + 0.5 * cn; c /= np.linalg.norm(c); d = c
            else:
                F = np.einsum("abcd,b,d->ac", Tm, d, d); F = 0.5 * (F + F.T)
                w, U = np.linalg.eigh(F); c = U[:, 0]
                G = np.einsum("abcd,a,c->bd", Tm, c, c); G = 0.5 * (G + G.T)
                w, U = np.linalg.eigh(G); d = U[:, 0]
            E = float(np.einsum("abcd,a,b,c,d->", Tm, c, d, c, d) / (1 + s * pp * (c @ d) ** 2))
            if E_old is not None and abs(E - E_old) < tol:
                break
            E_old = E
        rr = self.g.r[:self.m]
        Pc = c @ self.Pb[l1]; Pd = d @ self.Pb[l2]
        r1 = float(np.sum(Pc ** 2 * rr * self.w[:self.m])); r2 = float(np.sum(Pd ** 2 * rr * self.w[:self.m]))
        return E, it, r1, r2


def tail(Es, ls):
    """E(l -> inf) from the last increments: dE_l ~ (l + 1/2)^-p"""
    d = np.diff(Es)
    if len(d) < 2 or abs(d[-1]) < 1e-12:
        return Es[-1], None
    ratio = d[-1] / d[-2]
    if ratio <= 0 or ratio >= 1:
        p = 4.0
    else:
        p = math.log(ratio) / math.log((ls[-2] + 0.5) / (ls[-1] + 0.5))
    p = min(max(p, 3.0), 8.0)
    lm = ls[-1]
    s = sum(((l + 0.5) / (lm + 0.5)) ** (-p) for l in range(lm + 1, lm + 2000))
    return Es[-1] + d[-1] * s, p


def run_state(Z, state, lmax=6, n_rad=24, R_box=12.0, field=None, log=print):
    L, S, par, (l1, l2) = STATES[state]
    f = field or P.BlockField(Z)
    t0 = time.time()
    ci = PairCI(f, lmax=lmax, n_rad=n_rad, R_box=R_box, log=log)
    Ehf, it, r1, r2 = ci.hf(L, S, l1, l2)
    nroot = 2 if state == "pp_1S" else 1
    Es, dims, ls = [], [], []
    lmin = max(l1, l2) + (1 if state in ("pp_3P", "pp_1D") else 0)
    for lm in range(max(1, l2), lmax + 1):
        w, dim = ci.solve(L, S, par, lmax=lm, nroot=nroot)
        Es.append(float(w[nroot - 1])); dims.append(int(dim)); ls.append(lm)
        log(f"  Z={Z} {state} lmax={lm} dim={dim} E={Es[-1]:.9f} ({time.time() - t0:.0f} s)")
    Einf, p = tail(Es, ls)
    return dict(Z=Z, stato=state, lmax=lmax, n_rad=n_rad, R_box=R_box, E_HF=Ehf, HF_iter=it, r_medio=[r1, r2],
                E_l=Es, l=ls, dim=dims, E_inf=Einf, p_coda=p, Ec=Einf - Ehf, Ec_lmax=Es[-1] - Ehf,
                secondi=time.time() - t0)


if __name__ == "__main__":
    Z = int(sys.argv[1]); st = sys.argv[2]
    a = [float(x) for x in sys.argv[3:]]
    kw = {}
    if a:
        kw = dict(lmax=int(a[0]), n_rad=int(a[1]), R_box=a[2])
    print(json.dumps(run_state(Z, st, **kw), indent=1))
