"""multi_e_ecg.py - livello ESATTO a pochi corpi (N <= 4, L = 0): gaussiane correlate esplicite (ECG) con crescita
stocastica della base (SVM), con l'operatore della regola (docs/MOLTI_ELETTRONI_PIANO.md par. 2 L4, 3.4).

  phi_A(r) = exp(-r^T (A x I3) r), r = posizioni degli elettroni rispetto al nucleo (massa infinita, DEBT-3).
  Antisimmetria (DEBT-2): funzione di spin primitiva up..up down..down, somma con segno sulle permutazioni degli
  elettroni con lo stesso spin (S_Nup x S_Ndown), B -> P^T B P.
  S = (pi^N/det C)^{3/2}, C = A + B';  <1/|w^T r|> = S (2/sqrt pi) (w^T C^-1 w)^{-1/2};
  cinetica NR: S (3/m_kin) tr[A C^-1 B'];  cinetica piena della regola: S sum_i F(s_i), s_i = 2 [A C^-1 B']_ii
  (F = media della banda imparata su una gaussiana isotropa di varianza s per componente, tabella di multi_e_disp.py).

Uso:  python lab/multi_e_ecg.py g1          gate G1 (H, Hooke, non interagenti Li/Be, Monte Carlo)
      python lab/multi_e_ecg.py run <job>  una corsa (vedi JOBS), checkpoint in lab/out/multi_e_ecg.json
      python lab/multi_e_ecg.py all        tutte le corse NR di M1 in fila (per nohup)"""
import os, sys, time, json, itertools
from multi_e_core import OUT, ALPHA, load_json, save_json
import numpy as np
import scipy.linalg as sla
from scipy.optimize import brentq

SQPI = np.sqrt(np.pi)
JS = "multi_e_ecg.json"


# =============================================================================================== operator
class System:
    def __init__(self, N, n_up, Z, m_kin=1.0, kin="NR", Ftab=None, omega=0.0, ee=True, bosonic=False, name=""):
        self.N, self.n_up, self.Z, self.m, self.kin = N, n_up, Z, m_kin, kin
        self.Ftab, self.omega, self.ee, self.name = Ftab, omega, ee, name
        ups, dns = list(range(n_up)), list(range(n_up, N))
        self.perms = []
        for pu in itertools.permutations(ups):
            for pd in itertools.permutations(dns):
                p = list(pu) + list(pd)
                sg = 1 if bosonic else perm_sign(p)
                P = np.zeros((N, N)); P[np.arange(N), p] = 1
                self.perms.append((sg, P))
        self.pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]

    def F(self, s):
        sv, R = self.Ftab
        ls = np.log(np.clip(s, sv[0], sv[-1]))
        return 1.5 * s * np.interp(ls, np.log(sv), R)

    def me(self, A, B, parts=False):
        """A, B (n, N, N) -> S, H (n,) summed over signed permutations of the ket"""
        n = A.shape[0]
        S_t = np.zeros(n); H_t = np.zeros(n); T_t = np.zeros(n); Ve_t = np.zeros(n); Vee_t = np.zeros(n); P4_t = np.zeros(n)
        for sg, P in self.perms:
            Bp = np.einsum("ji,njk,kl->nil", P, B, P)
            C = A + Bp
            Ci = np.linalg.inv(C)
            dC = np.linalg.det(C)
            S = (np.pi ** self.N / dC) ** 1.5
            ACB = A @ Ci @ Bp
            if self.kin == "NR":
                T = S * (3 / self.m) * np.trace(ACB, axis1=1, axis2=2)
            else:
                s = 2 * np.diagonal(ACB, axis1=1, axis2=2)
                T = S * np.sum(self.F(s), axis=1)
            dg = np.diagonal(Ci, axis1=1, axis2=2)
            Ve = -self.Z * S * (2 / SQPI) * np.sum(dg ** -0.5, axis=1) if self.Z else 0 * S
            Vee = 0 * S
            if self.ee:
                for i, j in self.pairs:
                    Vee = Vee + S * (2 / SQPI) * (Ci[:, i, i] + Ci[:, j, j] - 2 * Ci[:, i, j]) ** -0.5
            Vh = 0.5 * self.omega ** 2 * S * 1.5 * np.sum(dg, axis=1) if self.omega else 0
            S_t += sg * S; H_t += sg * (T + Ve + Vee + Vh)
            if parts:
                s4 = 2 * np.diagonal(ACB, axis1=1, axis2=2)
                T_t += sg * T; Ve_t += sg * Ve; Vee_t += sg * Vee; P4_t += sg * S * np.sum(15 * s4 ** 2, axis=1)
        if parts:
            return S_t, H_t, dict(T=T_t, Ven=Ve_t, Vee=Vee_t, P4=P4_t)
        return S_t, H_t


def perm_sign(p):
    p = list(p); s = 1
    for i in range(len(p)):
        while p[i] != i:
            j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
    return s


# =============================================================================================== random functions
def build_A(N, a_en, a_ee):
    """A = sum_i a_i e_i e_i^T + sum_{i<j} a_ij (e_i - e_j)(e_i - e_j)^T"""
    A = np.diag(a_en).astype(float)
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            A[i, i] += a_ee[k]; A[j, j] += a_ee[k]; A[i, j] -= a_ee[k]; A[j, i] -= a_ee[k]
            k += 1
    return A


def random_A(rng, N, Zs, Lmax_en=20.0):
    while True:
        L = np.exp(rng.uniform(np.log(0.01 / Zs), np.log(Lmax_en / Zs), N))
        a_en = 1 / L ** 2
        npair = N * (N - 1) // 2
        a_ee = np.zeros(npair)
        for k in range(npair):
            u = rng.random()
            if u < 0.2:
                continue
            L12 = np.exp(rng.uniform(np.log(0.05 / Zs), np.log(30 / Zs)))
            a_ee[k] = (1 if u < 0.84 else -1) / L12 ** 2       # + con prob. 0.64/0.8 = 0.8 fra i presenti, - con 0.2
        A = build_A(N, a_en, a_ee)
        if np.linalg.eigvalsh(A).min() > 1e-12 * np.abs(A).max():
            return A


# =============================================================================================== SVM
class SVM:
    def __init__(self, sys_, rng, Zs, Lmax=20.0, trials=30):
        self.s, self.rng, self.Zs, self.Lmax, self.trials = sys_, rng, Zs, Lmax, trials
        self.A = np.zeros((0, sys_.N, sys_.N))
        self.nset = np.zeros(0)                       # 1/sqrt(<A|A>) delle funzioni della base
        self.Smat = np.zeros((0, 0)); self.Hmat = np.zeros((0, 0))
        self.E = None; self.U = None; self.Ev = None
        self.hist = []

    def norm_A(self, A):
        S, _ = self.s.me(A[None], A[None])
        return S[0]

    def _cand_rows(self, Ac, Aset):
        """matrix elements of candidates Ac (m,N,N) with the set Aset (K,N,N) -> (m,K), plus diagonals"""
        m, K = Ac.shape[0], Aset.shape[0]
        if K:
            S, H = self.s.me(np.repeat(Ac, K, 0), np.tile(Aset, (m, 1, 1)))
            S = S.reshape(m, K); H = H.reshape(m, K)
        else:
            S = H = np.zeros((m, 0))
        Sd, Hd = self.s.me(Ac, Ac)
        return S, H, Sd, Hd

    def _eig(self, Hm, Sm):
        Ev, U = sla.eigh(Hm, Sm)
        return Ev, U

    def _lowest_with(self, Ev, U, s, h, sd, hd):
        """lowest eigenvalue after adding one function (rank-1 border, secular equation); None if near-dependent"""
        if U is None or len(Ev) == 0:
            return hd / sd
        us = U.T @ s; uh = U.T @ h
        nn = sd - us @ us
        if nn / sd < 1e-10:
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

    def _best_of(self, Aset, Smat, Hmat, nset, extra=None):
        """try self.trials random candidates (+ extra) against the set; return (E, A, row S, row H, Sd, Hd)"""
        if len(Aset):
            Ev, U = self._eig(Hmat, Smat)
        else:
            Ev, U = np.zeros(0), None
        cands = np.array([random_A(self.rng, self.s.N, self.Zs, self.Lmax) for _ in range(self.trials)])
        if extra is not None:
            cands = np.concatenate([cands, extra[None]])
        # normalizzo
        Sd0, _ = self.s.me(cands, cands)
        ok = Sd0 > 0
        cands = cands[ok]
        S, H, Sd, Hd = self._cand_rows(cands, Aset)
        nrm = 1 / np.sqrt(Sd)
        best = None
        for i in range(len(cands)):
            s = S[i] * nrm[i] * nset; h = H[i] * nrm[i] * nset
            e = self._lowest_with(Ev, U, s, h, 1.0, Hd[i] * nrm[i] ** 2)
            if e is not None and (best is None or e < best[0]):
                best = (e, cands[i], s, h, 1.0, Hd[i] * nrm[i] ** 2, i == len(cands) - 1 and extra is not None, nrm[i])
        return best

    def add(self):
        K = len(self.A)
        b = self._best_of(self.A, self.Smat, self.Hmat, self.nset)
        e, A, s, h, sd, hd, _, nn = b
        self.nset = np.concatenate([self.nset, [nn]])
        # salvo la funzione normalizzata: memorizzo A e lo scalare di norma nelle matrici
        Sm = np.zeros((K + 1, K + 1)); Hm = np.zeros((K + 1, K + 1))
        Sm[:K, :K] = self.Smat; Hm[:K, :K] = self.Hmat
        Sm[K, :K] = Sm[:K, K] = s; Hm[K, :K] = Hm[:K, K] = h
        Sm[K, K] = sd; Hm[K, K] = hd
        self.A = np.concatenate([self.A, A[None]])
        self.Smat, self.Hmat = Sm, Hm
        self.E = float(sla.eigh(Hm, Sm, eigvals_only=True, subset_by_index=[0, 0])[0])
        return self.E

    def refine_sweep(self):
        K = len(self.A)
        for k in range(K):
            keep = np.array([i for i in range(K) if i != k])
            Sm = self.Smat[np.ix_(keep, keep)]; Hm = self.Hmat[np.ix_(keep, keep)]
            b = self._best_of(self.A[keep], Sm, Hm, self.nset[keep], extra=self.A[k])
            e, A, s, h, sd, hd, is_old, nn = b
            if is_old:
                continue
            self.A[k] = A; self.nset[k] = nn
            # riscrivo riga/colonna k (le righe s, h sono rispetto a keep)
            self.Smat[k, keep] = self.Smat[keep, k] = s; self.Hmat[k, keep] = self.Hmat[keep, k] = h
            self.Smat[k, k] = sd; self.Hmat[k, k] = hd
        self.E = float(sla.eigh(self.Hmat, self.Smat, eigvals_only=True, subset_by_index=[0, 0])[0])
        return self.E

    def energy_subset(self, K):
        return float(sla.eigh(self.Hmat[:K, :K], self.Smat[:K, :K], eigvals_only=True, subset_by_index=[0, 0])[0])

    def rebuild(self, sys_new=None):
        """recompute S, H for the current basis (optionally with another operator, e.g. the full band)"""
        if sys_new is not None:
            self.s = sys_new
        K = len(self.A)
        I, J = np.triu_indices(K)
        S, H = self.s.me(self.A[I], self.A[J])
        Sd = np.zeros((K, K)); Hd = np.zeros((K, K))
        Sd[I, J] = S; Sd[J, I] = S; Hd[I, J] = H; Hd[J, I] = H
        nr = 1 / np.sqrt(np.diag(Sd))
        self.nset = nr
        self.Smat = Sd * np.outer(nr, nr); self.Hmat = Hd * np.outer(nr, nr)
        self.E = float(sla.eigh(self.Hmat, self.Smat, eigvals_only=True, subset_by_index=[0, 0])[0])
        return self.E

    def state(self):
        """coefficients (for normalized functions) of the lowest state"""
        Ev, U = sla.eigh(self.Hmat, self.Smat, subset_by_index=[0, 0])
        return Ev[0], U[:, 0]

    def expect_parts(self):
        """<T>, <Ven>, <Vee>, <sum p^4>, s-values weights (per G4) for the lowest state"""
        K = len(self.A)
        E, c = self.state()
        I, J = np.triu_indices(K)
        S, H, P = self.s.me(self.A[I], self.A[J], parts=True)
        nrm = self.nset
        w = c[I] * c[J] * nrm[I] * nrm[J] * np.where(I == J, 1.0, 2.0)
        out = {k: float(np.sum(w * v)) for k, v in P.items()}
        out["norm"] = float(np.sum(w * S)); out["E"] = float(np.sum(w * H))
        for k in ("T", "Ven", "Vee", "P4"):
            out[k] /= out["norm"]
        out["E"] /= out["norm"]
        return out


def run_svm(sys_, Kmax, Zs, seed, Lmax=20.0, trials=30, sweeps=3, log=print, tmax=None):
    rng = np.random.default_rng(seed)
    sv = SVM(sys_, rng, Zs, Lmax, trials)
    t0 = time.time()
    EK = []
    while len(sv.A) < Kmax:
        e = sv.add()
        if len(sv.A) % 10 == 0:
            EK.append((len(sv.A), e))
            log(f"   K={len(sv.A)} E={e:.12f} ({time.time() - t0:.0f} s)")
    Es = []
    for sw in range(sweeps):
        e = sv.refine_sweep()
        Es.append(e)
        log(f"   sweep {sw + 1}: E={e:.12f} ({time.time() - t0:.0f} s)")
    return sv, EK, Es


# =============================================================================================== G1
def g1():
    res = load_json(JS)
    out = res.get("G1", {})
    t0 = time.time()
    # (a) idrogeno
    if "H" not in out:
        s = System(1, 1, 1.0, 1.0, name="H")
        sv, EK, Es = run_svm(s, 40, 1.0, 11, trials=30, sweeps=3)
        out["H"] = dict(E=sv.E, err=sv.E + 0.5, K=40, PASS=bool(abs(sv.E + 0.5) < 1e-8), EK=EK, sweeps=Es)
        print("G1 H:", sv.E, sv.E + 0.5, flush=True)
        res["G1"] = out; save_json(JS, res)
    # (b) Hooke omega = 1/2
    if "Hooke" not in out:
        s = System(2, 1, 0.0, 1.0, omega=0.5, name="Hooke")
        sv, EK, Es = run_svm(s, 200, 1.0, 12, trials=30, sweeps=3)
        out["Hooke"] = dict(E=sv.E, err=sv.E - 2.0, K=200, PASS=bool(abs(sv.E - 2.0) < 1e-8), EK=EK, sweeps=Es)
        print("G1 Hooke:", sv.E, sv.E - 2, flush=True)
        res["G1"] = out; save_json(JS, res)
    # (c) non interagenti, base prodotto (DV2)
    for name, Z, n, nup, N in (("Li_nonint", 3.0, 14, 2, 3), ("Be_nonint", 4.0, 10, 2, 4)):
        if name in out:
            continue
        a = 0.0075 * Z ** 2 * 2.9 ** np.arange(n)
        # una particella, codice indipendente
        AA, BB = np.meshgrid(a, a, indexing="ij")
        S1 = (np.pi / (AA + BB)) ** 1.5
        T1 = 3 * AA * BB / (AA + BB) * S1
        V1 = -Z * 2 / SQPI * np.sqrt(AA + BB) * S1
        e1 = sla.eigh(T1 + V1, S1, eigvals_only=True)
        target = (2 * e1[0] + e1[1]) if N == 3 else (2 * e1[0] + 2 * e1[1])
        ups = list(itertools.combinations(range(n), 2))
        dns = list(range(n)) if N == 3 else list(itertools.combinations(range(n), 2))
        funcs = []
        for u in ups:
            for d in dns:
                idx = list(u) + ([d] if N == 3 else list(d))
                funcs.append(np.diag(a[idx]))
        Fa = np.array(funcs)
        K = len(Fa)
        r = {}
        for bos in (False, True):
            s = System(N, nup, Z, 1.0, ee=False, bosonic=bos)
            I, J = np.triu_indices(K)
            Sv = np.zeros(len(I)); Hv = np.zeros(len(I))
            ch = 200000
            for i0 in range(0, len(I), ch):
                Sv[i0:i0 + ch], Hv[i0:i0 + ch] = s.me(Fa[I[i0:i0 + ch]], Fa[J[i0:i0 + ch]])
            Sm = np.zeros((K, K)); Hm = np.zeros((K, K))
            Sm[I, J] = Sv; Sm[J, I] = Sv; Hm[I, J] = Hv; Hm[J, I] = Hv
            nr = 1 / np.sqrt(np.diag(Sm))
            Sm *= np.outer(nr, nr); Hm *= np.outer(nr, nr)
            if bos:
                # nel caso bosonico alcune funzioni (a = b) non ci sono; tolgo le direzioni quasi dipendenti
                w_, V_ = np.linalg.eigh(Sm)
                keep = w_ > 1e-10
                X = V_[:, keep] / np.sqrt(w_[keep])
                Eb = float(np.linalg.eigvalsh(X.T @ Hm @ X)[0])
                r["E_bosonic"] = Eb
            else:
                E = float(sla.eigh(Hm, Sm, eigvals_only=True, subset_by_index=[0, 0])[0])
                r["E_fermionic"] = E
        exact = -Z ** 2 * (2 + 0.25) / 2 if N == 3 else -Z ** 2 * (2 + 0.5) / 2
        out[name] = dict(K=K, n1=n, E=r["E_fermionic"], target_same_basis=float(target), diff=r["E_fermionic"] - target,
                         E_bosonic=r["E_bosonic"], three_1s_same_basis=float(N * e1[0]), exact_limit=exact,
                         basis_error_report_only=r["E_fermionic"] - exact,
                         PASS=bool(abs(r["E_fermionic"] - target) < 1e-8 and r["E_bosonic"] < r["E_fermionic"] - 1e-3))
        print("G1", name, out[name], f"{time.time() - t0:.0f}s", flush=True)
        res["G1"] = out; save_json(JS, res)
    # (d) Monte Carlo di <1/r_ij> e <1/r_i>
    if "MC" not in out:
        rng = np.random.default_rng(5)
        N = 3
        recs = []
        for t in range(3):
            A = random_A(rng, N, 2.0); B = random_A(rng, N, 2.0)
            C = A + B
            Ci = np.linalg.inv(C)
            Lc = np.linalg.cholesky(Ci / 2)                  # r ~ N(0, C^-1/2) per componente
            acc = np.zeros(2); nsamp = 0
            for ch in range(40):
                z = rng.normal(size=(10 ** 7, N, 3)) if False else rng.normal(size=(N, 3, 10 ** 7 // 4))
                r = np.einsum("ij,jkn->ikn", Lc, z)
                acc[0] += np.sum(1 / np.linalg.norm(r[0] - r[1], axis=0))
                acc[1] += np.sum(1 / np.linalg.norm(r[2], axis=0))
                nsamp += z.shape[-1]
            mc = acc / nsamp
            an = np.array([(2 / SQPI) * (Ci[0, 0] + Ci[1, 1] - 2 * Ci[0, 1]) ** -0.5, (2 / SQPI) * Ci[2, 2] ** -0.5])
            recs.append(dict(analytic=an.tolist(), mc=mc.tolist(), rel=(mc / an - 1).tolist(), n=nsamp))
            print("MC", recs[-1], flush=True)
        worst = max(abs(x) for r_ in recs for x in r_["rel"])
        out["MC"] = dict(recs=recs, worst_rel=worst, PASS=bool(worst < 1e-4),
                         note="campioni 1e8 per elemento (prereg: 4e6, che da' sigma ~4e-4 > soglia: cambio dichiarato)")
        res["G1"] = out; save_json(JS, res)
    out["PASS"] = bool(all(out[k]["PASS"] for k in ("H", "Hooke", "Li_nonint", "Be_nonint", "MC")))
    res["G1"] = out; save_json(JS, res)
    print("G1:", {k: (v["PASS"] if isinstance(v, dict) else v) for k, v in out.items()})


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "g1":
    g1()


# =============================================================================================== corse di produzione (M1)
# ATTENZIONE: G1 e' FALLITO (Hooke 2.7e-8 contro 1e-8, Monte Carlo 1.34e-4 contro 1e-4): per la regola del gate queste
# corse NON si lanciano finche' Valerio non decide (docs/MOLTI_ELETTRONI.md, "Dove Valerio puo' aiutare"). Il codice e'
# pronto; 'run' rifiuta di partire senza la variabile MULTI_E_G1_DECISION che rimanda alla decisione scritta.
def job_list():
    J = {}
    for Zi in range(1, 11):
        J[f"he1S_Z{Zi}"] = dict(N=2, nup=1, Z=float(Zi), K=200, Lmax=30.0 if Zi == 1 else 20.0)
    for Zi in range(2, 11):
        J[f"he3S_Z{Zi}"] = dict(N=2, nup=2, Z=float(Zi), K=150, Lmax=20.0)
    for Zi in (1, 2, 10):
        J[f"one_Z{Zi}"] = dict(N=1, nup=1, Z=float(Zi), K=40, Lmax=20.0)
    return J


def run_job(name, seed=2026):
    if not os.environ.get("MULTI_E_G1_DECISION"):
        raise SystemExit("G1 FALLITO: produzione ECG ferma finche' Valerio non decide (vedi docs/MOLTI_ELETTRONI.md)")
    res = load_json(JS)
    if name in res.get("runs", {}):
        return
    m_A = load_json("multi_e_kernel.json")["R3"]["m_kin_theta_A_used"]
    j = job_list()[name]
    s = System(j["N"], j["nup"], j["Z"], m_A, name=name)
    sv, EK, Es = run_svm(s, j["K"], j["Z"], seed, Lmax=j["Lmax"])
    sv.rebuild()
    parts = sv.expect_parts()
    rec = dict(E=sv.E, E_Kminus50=sv.energy_subset(j["K"] - 50) if j["K"] > 50 else None, EK=EK, sweeps=Es, parts=parts,
               m_kin=m_A, K=j["K"], decision=os.environ.get("MULTI_E_G1_DECISION"))
    np.savez(os.path.join(OUT, f"multi_e_ecg_{name}.npz"), A=sv.A)
    res = load_json(JS)
    res.setdefault("runs", {})[name] = rec
    save_json(JS, res)


if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[1] == "run":
    run_job(sys.argv[2])
