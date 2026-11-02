"""multi_e_meanfield.py - campo medio "ascoltato" dall'operatore della regola (DEBT-5): ogni orbitale e' lo stato
stazionario dell'operatore non relativistico della regola (cinetica k^2/2m_kin, nucleo puntiforme) nel campo del nucleo
e degli altri elettroni (Hartree dai gusci di Gauss, scambio dall'antisimmetria con Y^k e fattori angolari s/p esatti).
In linguaggio standard: Hartree-Fock radiale (ROHF a media sferica per i gusci aperti). Gate G3 (unit test con valori
della letteratura, m_kin = 1) e G5 (rapporto 2s/1s monotono per He 2^3S).

Griglia: x = ln r uniforme, P(r) = r^{1/2} f(x); -(1/2m) f'' + (l+1/2)^2/(2m) f + r^2 V f + [non locale] = eps r^2 f;
derivata seconda alle differenze finite del 4o ordine; integrali con la regola del trapezio in x (dr = r dx).
Run: python lab/multi_e_meanfield.py g3 | g5 | he"""
import os, sys, time
from multi_e_core import load_json, save_json
import numpy as np
import scipy.linalg as sla

JS = "multi_e_mf.json"
W3 = {(0, 0): {0: 1.0}, (0, 1): {1: 1 / 3}, (1, 0): {1: 1 / 3}, (1, 1): {0: 1 / 3, 2: 2 / 15}}   # (l k l';000)^2


class Grid:
    def __init__(self, Z, n, rmax=60.0, rmin_fac=1e-9):
        self.x = np.linspace(np.log(rmin_fac / Z), np.log(rmax), n)
        self.h = self.x[1] - self.x[0]
        self.r = np.exp(self.x)
        self.n = n
        e = np.ones(n)
        D2 = (np.diag(-30 * e) + np.diag(16 * e[:-1], 1) + np.diag(16 * e[:-1], -1) - np.diag(e[:-2], 2) - np.diag(e[:-2], -2))
        self.D2 = D2 / (12 * self.h ** 2)
        self.R32 = self.r ** 1.5
        rr = self.r
        rmin = np.minimum.outer(rr, rr); rmax = np.maximum.outer(rr, rr)
        self.Kk = {k: rmin ** k / rmax ** (k + 1) for k in (0, 1, 2)}      # r_<^k / r_>^{k+1}, una volta per griglia
        del rmin, rmax

    def P(self, f):
        return np.sqrt(self.r) * f

    def integ(self, g):                       # int g dr
        return np.sum(g * self.r) * self.h


class Atom:
    """shells: list of (name, l, q) ; A[a][b], B[(a,b,k)] coefficients of E = sum q_a h_a + 1/2 sum (A F0 - B G^k)"""

    def __init__(self, Z, shells, A, B, m_kin=1.0, n=3000, rmax=60.0):
        self.Z, self.shells, self.A, self.B, self.m = Z, shells, A, B, m_kin
        self.g = Grid(Z, n, rmax)
        g = self.g
        self.ns = len(shells)
        self.T = {l: -g.D2 / (2 * self.m) + np.eye(g.n) * (l + 0.5) ** 2 / (2 * self.m) for l in (0, 1)}
        self.Mdiag = g.r ** 2

    # --- potentials
    def Yk_over_r(self, Pa, Pb, k):
        """Y^k[Pa Pb](r)/r on the grid: int Pa Pb(r') r_<^k / r_>^{k+1} dr'"""
        g = self.g
        return g.Kk[k] @ (Pa * Pb * g.r) * g.h

    def fock(self, a, orbs, l):
        """matrix (in f, generalized form) of the Fock operator of shell a, given the current orbitals (P on grid)"""
        g = self.g
        q = self.shells[a][2]
        V = -self.Z / g.r
        loc = np.zeros(g.n)
        for b in range(self.ns):
            if self.A[a][b]:
                loc += self.A[a][b] * self.Yk_over_r(orbs[b], orbs[b], 0) / q
        F = self.T[l] + np.diag(g.r ** 2 * (V + loc))
        for (aa, b, k), cf in self.B.items():
            if aa != a or not cf:
                continue
            Pb = orbs[b]
            K = g.Kk[k]
            # (K_b P)(r) = Pb(r) int Pb(r') P(r') r_<^k/r_>^{k+1} dr'  ->  r^{3/2} [..] r'^{3/2} h  in f-space
            Kf = (g.R32 * Pb)[:, None] * K * (Pb * g.R32)[None, :] * g.h
            F = F - (cf / q) * Kf
        return F

    def energy(self, orbs, fs):
        g = self.g
        E = 0.0
        parts = {}
        for a, (nm, l, q) in enumerate(self.shells):
            f = fs[a]
            ha = g.h * f @ ((self.T[l]) @ f) + g.integ(orbs[a] ** 2 * (-self.Z / g.r))
            E += q * ha
            parts[f"h_{nm}"] = float(ha)
        for a in range(self.ns):
            for b in range(self.ns):
                if self.A[a][b]:
                    E += 0.5 * self.A[a][b] * g.integ(orbs[a] ** 2 * self.Yk_over_r(orbs[b], orbs[b], 0))
        for (a, b, k), cf in self.B.items():
            if cf:
                E -= 0.5 * cf * g.integ(orbs[a] * orbs[b] * self.Yk_over_r(orbs[a], orbs[b], k))
        return float(E), parts

    def solve_l(self, F, l, exclude=None, nroot=1, sigma=None, iters=200, tol=1e-14, guess=None):
        """lowest roots of F f = eps M f (M = diag r^2) by SHIFT-INVERT subspace iteration + Rayleigh-Ritz in the M metric
        (no 1/r scaling: the direct symmetric form has entries ~1/r_min^2 ~ 1e18 and loses all precision),
        optionally in the M-orthogonal complement of the vectors in exclude"""
        g = self.g
        M = g.r ** 2
        if sigma is None:
            sigma = -0.6 * self.Z ** 2 - 1.0
        # constraint x ._M e = 0 for e in exclude: bordered shift-invert system (exact constrained resolvent)
        E_ = np.array(exclude).T if exclude else np.zeros((g.n, 0))
        k = E_.shape[1]
        Big = np.zeros((g.n + k, g.n + k))
        Big[:g.n, :g.n] = F - sigma * np.diag(M)
        Big[:g.n, g.n:] = M[:, None] * E_
        Big[g.n:, :g.n] = (M[:, None] * E_).T
        lu = sla.lu_factor(Big)

        def proj(X):
            return X

        def solve(RHS):
            Z_ = np.zeros((g.n + k, RHS.shape[1]))
            Z_[:g.n] = RHS
            return sla.lu_solve(lu, Z_)[:g.n]
        m = nroot + 3
        rng = np.random.default_rng(0)
        X = rng.normal(size=(g.n, m)) * (g.r ** (l + 0.5) * np.exp(-g.r))[:, None]
        if guess is not None:
            X[:, :len(guess)] = np.array(guess).T
        if k:
            X = solve(M[:, None] * X)
        wold = None
        for it in range(iters):
            Y = solve(M[:, None] * X)
            Fs = Y.T @ F @ Y; Ms = Y.T @ (M[:, None] * Y)
            Fs = 0.5 * (Fs + Fs.T); Ms = 0.5 * (Ms + Ms.T)
            w, C = sla.eigh(Fs, Ms)
            X = Y @ C
            X = X / np.sqrt(np.sum(X * X * M[:, None], 0) * g.h)[None, :] if False else X / np.sqrt(np.sum(X * X * M[:, None], 0))[None, :]
            if wold is not None and np.max(np.abs(w[:nroot] - wold[:nroot])) < tol * max(1, abs(w[0])):
                break
            wold = w
        fs = [X[:, i] / np.sqrt(g.h * np.sum(X[:, i] ** 2 * M)) for i in range(nroot)]
        fs = [f * np.sign(f[np.argmax(np.abs(f))]) for f in fs]
        return w[:nroot], fs

    def start(self, kind, rng=None):
        """starting orbitals: 'hydrogen' (bare-nucleus eigenvectors of the rule's operator) or 'random' (band-limited)"""
        g = self.g
        fs = [None] * self.ns
        by_l = {}
        for a, (nm, l, q) in enumerate(self.shells):
            by_l.setdefault(l, []).append(a)
        for l, idx in by_l.items():
            if kind == "hydrogen":
                F = self.T[l] + np.diag(g.r ** 2 * (-self.Z / g.r))
                _, ff = self.solve_l(F, l, nroot=len(idx))
            else:
                ff = []
                for i in idx:
                    t = (g.x - g.x[0]) / (g.x[-1] - g.x[0])
                    c = rng.normal(size=12)
                    f = sum(c[j] * np.sin((j + 1) * np.pi * t) for j in range(12)) * np.exp(-g.r / 4) * g.r ** (l + 0.5)
                    ff.append(f)
                # Gram-Schmidt in the metric r^2
                out = []
                for f in ff:
                    for e in out:
                        f = f - (g.h * np.sum(f * e * g.r ** 2)) * e
                    out.append(f / np.sqrt(g.h * np.sum(f * f * g.r ** 2)))
                ff = out
            for a, f in zip(idx, ff):
                fs[a] = f
        return fs

    def scf(self, fs, maxit=400, tol=1e-12, mix=0.5, log=None):
        g = self.g
        Fold = None
        Eold = None
        for it in range(maxit):
            orbs = [g.P(f) for f in fs]
            Fm = [self.fock(a, orbs, self.shells[a][1]) for a in range(self.ns)]
            if Fold is not None:
                Fm = [mix * F + (1 - mix) * Fo for F, Fo in zip(Fm, Fold)]
            Fold = Fm
            new = list(fs)
            groups = {}
            for a, (nm, l, q) in enumerate(self.shells):
                groups.setdefault(l, []).append(a)
            for l, idx in groups.items():
                # shells with the same Fock operator (closed shells: identical) -> lowest roots of that operator
                same = all(np.allclose(Fm[idx[0]], Fm[b]) for b in idx)
                if same:
                    _, ff = self.solve_l(Fm[idx[0]], l, nroot=len(idx), guess=[fs[a] for a in idx])
                    for a, f in zip(idx, ff):
                        new[a] = f
                else:
                    # open shells: each orbital from its own Fock operator, orthogonal to the others of the same l;
                    # then the optimal rotation between them (direct minimisation of E)
                    # Gauss-Seidel by blocks: each orbital = lowest root of ITS OWN Fock operator (fresh, unmixed),
                    # in the complement of the other orbitals of the same l (exact constrained minimum when E is
                    # quadratic in that orbital: monotone descent)
                    for a in idx:
                        cur = [g.P(f) for f in new]
                        Fa = self.fock(a, cur, l)
                        others = [new[b] for b in idx if b != a]
                        _, ff = self.solve_l(Fa, l, exclude=others, nroot=1, guess=[new[a]])
                        new[a] = ff[0]
                    if len(idx) == 2:
                        new = self.best_rotation(new, idx[0], idx[1])
            fs = new
            orbs = [g.P(f) for f in fs]
            E, _ = self.energy(orbs, fs)
            if log:
                log(it, E)
            if Eold is not None and abs(E - Eold) < tol and it > 5:
                break
            Eold = E
        return E, fs, it

    def best_rotation(self, fs, a, b):
        from scipy.optimize import minimize_scalar
        g = self.g

        def rot(t):
            c, s = np.cos(t), np.sin(t)
            out = list(fs)
            out[a] = c * fs[a] + s * fs[b]; out[b] = -s * fs[a] + c * fs[b]
            return out

        def E(t):
            o = rot(t)
            return self.energy([g.P(f) for f in o], o)[0]
        r = minimize_scalar(E, bounds=(-0.3, 0.3), method="bounded", options={"xatol": 1e-12})
        return rot(r.x) if E(r.x) < E(0.0) else fs


# ---------------------------------------------------------------- atoms
def closed(Z, shells, m_kin=1.0, n=3000, rmax=60.0):
    """all shells closed: q = 2(2l+1); A = q_a q_b; B_abk = q_a q_b/2 (l_a k l_b;000)^2"""
    A = [[sa[2] * sb[2] for sb in shells] for sa in shells]
    B = {}
    for a, sa in enumerate(shells):
        for b, sb in enumerate(shells):
            for k, c in W3[(sa[1], sb[1])].items():
                B[(a, b, k)] = sa[2] * sb[2] / 2 * c
    return Atom(Z, shells, A, B, m_kin, n, rmax)


def make(name, m_kin=1.0, n=3000):
    if name == "He":
        return closed(2.0, [("1s", 0, 2)], m_kin, n)
    if name == "H-":
        return closed(1.0, [("1s", 0, 2)], m_kin, n, rmax=80.0)
    if name == "Be":
        return closed(4.0, [("1s", 0, 2), ("2s", 0, 2)], m_kin, n)
    if name == "Ne":
        return closed(10.0, [("1s", 0, 2), ("2s", 0, 2), ("2p", 1, 6)], m_kin, n)
    if name == "Li":            # 1s^2 2s ROHF: E = 2h1 + h2 + F0(11) + 2F0(12) - G0(12)
        A = [[4, 2], [2, 0]]
        B = {(0, 0, 0): 2, (0, 1, 0): 1, (1, 0, 0): 1}
        return Atom(3.0, [("1s", 0, 2), ("2s", 0, 1)], A, B, m_kin, n)
    if name == "He3S":          # 1s 2s triplet: E = h1 + h2 + F0(12) - G0(12)
        A = [[0, 1], [1, 0]]
        B = {(0, 1, 0): 1, (1, 0, 0): 1}
        return Atom(2.0, [("1s", 0, 1), ("2s", 0, 1)], A, B, m_kin, n)
    raise ValueError(name)


TARGET_G3 = {"He": -2.8616800, "Li": -7.4327269, "Be": -14.5730232, "Ne": -128.5470981, "H-": -0.4879297}

if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "g3"
    res = load_json(JS)
    if what == "g3":
        out = res.get("G3", {})
        names = sys.argv[2].split(",") if len(sys.argv) > 2 else list(TARGET_G3)
        for nm in names:
            rec = out.get(nm, {})
            t0 = time.time()
            Es = {}
            for n in (1500, 3000, 6000):
                if f"E_n{n}" in rec:
                    Es[n] = rec[f"E_n{n}"]; continue
                at = make(nm, 1.0, n)
                E, fs, it = at.scf(at.start("hydrogen"))
                Es[n] = E; rec[f"E_n{n}"] = E; rec[f"iters_n{n}"] = it
                print(f"{nm} n={n}: E = {E:.10f} ({it} it, {time.time() - t0:.0f} s)", flush=True)
                out[nm] = rec; res["G3"] = out; save_json(JS, res)
            e1, e2, e3 = Es[1500], Es[3000], Es[6000]
            p = np.log2(abs(e1 - e2) / abs(e2 - e3)) if abs(e2 - e3) > 0 else float("nan")
            ER = e3 + (e3 - e2) / (2 ** p - 1) if np.isfinite(p) and p > 0.5 else e3
            rel = ER / TARGET_G3[nm] - 1
            rec.update(order_p=float(p), E_richardson=float(ER), target=TARGET_G3[nm], rel_err=float(rel),
                       PASS=bool(abs(rel) < 2e-5))
            print(f"{nm}: Richardson {ER:.10f} (p = {p:.2f}) vs {TARGET_G3[nm]} rel {rel:.2e}", flush=True)
            out[nm] = rec; res["G3"] = out; save_json(JS, res)
    if what == "starts":
        # indipendenza dalla partenza: 5 campi casuali a banda limitata + idrogenoide (stessa griglia n = 1500)
        out = res.get("starts", {})
        for nm in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["He", "He3S", "Li", "H-"]):
            if nm in out:
                continue
            at = make(nm, 1.0, 1500)
            Eh_, _, _ = at.scf(at.start("hydrogen"), tol=1e-13)
            Er = []
            for sd in range(5):
                rng = np.random.default_rng(100 + sd)
                E, _, it = at.scf(at.start("random", rng), tol=1e-13, maxit=600)
                Er.append(E)
            dmax = float(max(abs(e - Eh_) for e in Er))
            out[nm] = dict(E_hydrogen_start=Eh_, E_random=Er, max_diff=dmax, PASS=bool(dmax < 1e-9))
            print(nm, out[nm], flush=True)
            res["starts"] = out; save_json(JS, res)
    if what == "carrier":
        # portatore dei nodi per la DMC dell'elio: He 1S e He 2^3S a m_kin(theta_A), n = 6000; G5; orbitali salvati
        m_kin = load_json("multi_e_kernel.json")["R3"]["m_kin_theta_A_used"]
        out = res.get("carrier", {})
        for nm in ("He", "He3S", "H-"):
            if nm in out:
                continue
            at = make(nm, m_kin, 6000)
            E, fs, it = at.scf(at.start("hydrogen"), tol=1e-13)
            g = at.g
            P = [g.P(f) for f in fs]
            rec = dict(E=E, m_kin=m_kin, iters=it, E_over_m=E / m_kin)
            if nm == "He3S":
                # orbitali canonici: 1s senza nodi; rapporto 2s/1s sui punti con |phi_1s| > 1e-300
                ok = np.abs(P[0]) > 1e-300
                ratio = P[1][ok] / P[0][ok]
                d = np.diff(ratio)
                mono = bool(np.all(d > 0) or np.all(d < 0))
                nz1 = int(np.sum(np.diff(np.sign(P[0][ok])) != 0))
                rec.update(G5_monotone=mono, nodes_1s=nz1, ratio_first=float(ratio[0]), ratio_last=float(ratio[-1]),
                           n_points=int(ok.sum()), worst_step=float(d.max() if d.mean() < 0 else -d.min()))
                print("G5:", mono, flush=True)
            np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", f"multi_e_orb_{nm}.npz"),
                     x=g.x, r=g.r, P=np.array(P), m_kin=m_kin, E=E)
            out[nm] = rec
            print(nm, rec, flush=True)
            res["carrier"] = out; save_json(JS, res)
