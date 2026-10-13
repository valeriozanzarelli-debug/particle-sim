"""locale_f5_kramers.py - phase-5 amendment (lab/out/locale_f5_amend_prereg.json): time reversal / Kramers pairing of
the walks (K1), and the time-reversal-respecting 'spin-link' twisted Wilson (K2). Writes lab/out/locale_f5_kramers.json."""
import os, sys, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import SplitWalk, sdot, zigzag_ticks, XCHI, YCHI, QR, CUBE_AXES, RealWalk, PAULI
from locale_cands import fermion
from locale_f4_core import walk as f4walk, measure

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "locale_f5_kramers.json")
SY = np.array([[0, -1j], [1j, 0]])


def spinlink_factors(h, vecs):
    """links L(x) <-> R(x+v), generator A_v = (sigma_y on chirality) (x) sigma_j for v = +-2 e_j; palindromic; the
    on-site part is undone EXACTLY by the inverse of the link product at k = 0"""
    allv = []
    for v in vecs:
        allv += [v, -v]
    seq = allv + allv[::-1]
    out, P0 = [], np.eye(4, dtype=complex)
    for v in seq:
        j = int(np.argmax(np.abs(v)))
        A = np.kron(SY, PAULI[j])
        C = np.cos(h / 2) * np.eye(4) - 1j * np.sin(h / 2) * A
        out.append(("link", v, QR, C))
        P0 = C @ P0
    out.append(("coin", np.linalg.inv(P0)))
    return out


def sl_walk(theta, h):
    fac = []
    for it, tk in enumerate(zigzag_ticks(3)):
        for d in tk:
            M = np.zeros((4, 4), complex); M[:2, :2] = -sdot(d, 0.5); M[2:, 2:] = sdot(d, 0.5)
            fac.append(("shift", d, M))
        vecs = [2 * e for e in CUBE_AXES] if it == 0 else [-2 * e for e in CUBE_AXES]
        fac += spinlink_factors(h / 2, vecs)
        fac.append(("coin", np.cos(theta) * np.eye(4) - 1j * np.sin(theta) * XCHI))
        fac += spinlink_factors(h / 2, vecs)
    return SplitWalk(fac, 4, np.eye(3), dim=3)


def kramers(w, n=6, theta=0.37, seed=7):
    rng = np.random.default_rng(seed)
    V = rng.normal(size=(n, n, n))
    rw = RealWalk(w, n, 3, V=V)
    D = 4 * n ** 3
    M = np.zeros((D, D), complex)
    for k in range(D):
        e = np.zeros(D, complex); e[k] = 1
        M[:, k] = rw.period(e.reshape(4, n, n, n)).ravel()
    ph = np.sort(np.angle(np.linalg.eigvals(M)))
    pair = float(np.max(np.abs(ph[0::2] - ph[1::2])))
    gap = float(np.min(ph[2::2] - ph[1:-1:2]))
    return dict(max_pair_split=pair, min_gap_between_pairs=gap, unitarity=float(np.abs(M.conj().T @ M - np.eye(D)).max()))


R = {}
th = 0.37
for name, w in (("F0", fermion("F0", 3, th)), ("F2_Y_h0.1", fermion("F2", 3, th, h=0.1)),
                ("F2b_X_h0.1", fermion("F2b", 3, th, h=0.1)), ("SL_h0.1", sl_walk(th, 0.1)),
                ("SL_h0.6", sl_walk(th, 0.6)), ("breath_mu0.3", f4walk(th, mu_x=0.3)),
                ("breath_mu0.3_SL0.1", None)):
    if w is None:
        continue
    R[f"K1|{name}"] = kramers(w)
    print(name, R[f"K1|{name}"], flush=True)
ref = fermion("F0", 3, 0.1)
for h in (0.1, 0.3, 0.6):
    r = measure(sl_walk(0.1, h), 0.1, ref)
    R[f"K2|SL|{h}"] = r
    print("SL", h, f"N {r['N_xy']:.3f} zone0 {r['zone_g0']:.3f} pi {r['X_x']:.3f} E0 {r['E0_rel']:.0e} drift {r['drift']:.0e} "
          f"mk {['%.1e' % x for x in r['mkin_vs_frozen']]} slope {r['dE_vs_frozen']['slope']}", flush=True)
json.dump(R, open(OUT, "w"), indent=1)


# ------------------------------------------------------------------ declared correction of the K1 test (after the run):
# the potential was applied once at the END of the period (RealWalk V). For the frozen rule V commutes with the mass coin
# and the period stays palindromic, but for walks with LINK factors V(x) does not commute with the links, so that test
# breaks the palindrome by itself. Fair test: V/2 applied together with each MASS coin (the coin not preceded by a link).
class RealWalkVcoin(RealWalk):
    def __init__(self, walk, n, V):
        super().__init__(walk, n, 3, V=None)
        self.Vc = np.exp(-0.5j * V)
        self.mass = [f[0] == "coin" and (i == 0 or walk.pre[i - 1][0] != "link") for i, f in enumerate(walk.pre)]

    def period(self, psi):
        for i, (f, iv) in enumerate(zip(self.w.pre, self.w.ivec)):
            if f[0] == "shift":
                _, v, lam, P = f
                c = np.tensordot(P.conj().T, psi, axes=(1, 0))
                c = np.array([self._roll(c[q], lam[q] * iv) if lam[q] else c[q] for q in range(len(lam))])
                psi = np.tensordot(P, c, axes=(1, 0))
            elif f[0] == "coin":
                psi = np.tensordot(f[1], psi, axes=(1, 0))
                if self.mass[i]:
                    psi = psi * self.Vc[None]
            else:
                _, v, Q, C = f
                qi = np.real(np.diag(Q)) > 0.5
                psi = np.array([self._roll(psi[q], iv) if qi[q] else psi[q] for q in range(len(qi))])
                psi = np.tensordot(C, psi, axes=(1, 0))
                psi = np.array([self._roll(psi[q], -iv) if qi[q] else psi[q] for q in range(len(qi))])
        return psi


def kramers_fair(w, n=6, seed=7):
    rng = np.random.default_rng(seed)
    V = rng.normal(size=(n, n, n))
    rw = RealWalkVcoin(w, n, V)
    D = 4 * n ** 3
    M = np.zeros((D, D), complex)
    for k in range(D):
        e = np.zeros(D, complex); e[k] = 1
        M[:, k] = rw.period(e.reshape(4, n, n, n)).ravel()
    ph = np.sort(np.angle(np.linalg.eigvals(M)))
    return dict(max_pair_split=float(np.max(np.abs(ph[0::2] - ph[1::2]))), n_mass_coins=int(sum(rw.mass)))


for name, w in (("F0", fermion("F0", 3, th)), ("F2_Y_h0.1", fermion("F2", 3, th, h=0.1)),
                ("F2b_X_h0.1", fermion("F2b", 3, th, h=0.1)), ("SL_h0.1", sl_walk(th, 0.1)),
                ("breath_mu0.3", f4walk(th, mu_x=0.3)), ("breath_mu0.3_Y0.1", f4walk(th, mu_x=0.3, hu=0.1))):
    R[f"K1fair|{name}"] = kramers_fair(w)
    print("fair", name, R[f"K1fair|{name}"], flush=True)
json.dump(R, open(OUT, "w"), indent=1)
