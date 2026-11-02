"""tavg_pecg.py - FASE 3 dello studio tavg_: il BLOCCO-NUCLEO dentro l'ECG (gaussiane correlate esplicite).

Gli elettroni di valenza (1, 2, 3) si muovono nel campo del blocco imparato (lab/out/tavg_blocchi/Z??.json) invece che
nel campo -Z/r del nucleo nudo. Il blocco ha due pezzi:
  * locale  Vloc(r) (canale p; coda -Zv/r oltre r_tail): elemento di matrice ESATTO per quadratura radiale. Per una
    funzione ECG la distribuzione marginale dell'elettrone i e' una gaussiana isotropa di varianza s2 = Sigma_ii per
    componente, quindi  <V(r_i)> = S0 G0(s2)  (stati S),  e per stati P (u.z)(v.z):
        S0 [ (u'Sv - g_u g_v / s2) G0(s2) + g_u g_v / s2^2 * G2(s2) / 3 ],  G0 = E[V(|y|)], G2 = E[|y|^2 V(|y|)]
    (stessa algebra della formula di Coulomb di tav_ecg). G0, G2 = parte di Coulomb analitica + parte corta
    dV = Vloc + Zv/r per quadratura sulla griglia logaritmica del livello fine, tabulata in log s2 (spline cubica).
  * separabile  |beta><beta| / E_kb, solo l = 0 (forma di Kleinman-Bylander del blocco). beta(r) = chi(r)/(r sqrt(4 pi))
    viene scritto come somma di gaussiane sum_k c_k exp(-alpha_k r^2) (minimi quadrati in L2, 36 termini fra 0.05 e
    500; regola automatica: il numero di termini piu' piccolo che ridà i 3 livelli s del blocco entro 2e-7 Eh sul
    livello radiale). Allora l'integrale sull'elettrone i e' analitico:
        int exp(-alpha r_i^2) exp(-x'Ax) dr_i = (pi/(a+alpha))^{3/2} exp(-R'(A_R - t t'/(a+alpha))R)   (t = A[i, R])
    e per le funzioni P il fattore lineare diventa (u_R - u_i t/(a+alpha)).z_R. Resta un integrale gaussiano su R.
  Il chi e il Vloc del blocco hanno rumore numerico per r < 1e-5 a0 (dall'inversione): chi sostituito da una parabola
  in r fittata su [1e-4, 1e-3], Vloc dal suo valore a 1e-4 (effetto sui livelli ~1e-12 Eh).
Spin e antisimmetria, cinetico, repulsione: tav_ecg.py (studio litio-berillio, commit b3a4ad7), importato senza
modifiche; m_kin = 1 (il blocco e' costruito a m_kin = 1).
Uso: python tavg_pecg.py check <Z>    controllo di macchina (1 elettrone, livelli s e p del blocco)"""
import os, sys, json, math, time
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
from scipy.interpolate import CubicSpline
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tav_ecg as T
import tavg_radial as tr

OUT = os.path.join(HERE, "out")
SQ2PI = math.sqrt(2 / math.pi)
KB_SETS = [(24, 0.05, 300.0), (30, 0.05, 400.0), (36, 0.05, 500.0), (40, 0.05, 800.0)]
KB_TOL = 2e-7


class BlockField:
    """the learned block of element Z on the fine radial grid, plus its ECG representation"""

    def __init__(self, Z, kb_set=None):
        b = json.load(open(os.path.join(OUT, "tavg_blocchi", f"Z{Z:02d}.json")))
        self.Z, self.Zv, self.rec = Z, float(b["Zv"]), b
        g = tr.Grid(Z); r = g.r; n = len(b["r"])
        assert np.allclose(g.r[:n], b["r"], rtol=1e-12)
        self.g = g
        V = np.empty(g.M); V[:n] = b["Vloc"]; V[n:] = -self.Zv / r[n:]
        # numerical noise of the inversion at r < 1e-5 a0 (|V| up to 1e12 at 1e-12): replaced by the value at 1e-4
        V[r < 1e-5] = V[np.searchsorted(r, 1e-4)]
        self.Vloc = V
        self.kb = []
        for kb in b["kb"]:
            chi = np.zeros(g.M); chi[:n] = kb["chi"]
            i0, i1 = np.searchsorted(r, 1e-4), np.searchsorted(r, 1e-3)
            p = np.polyfit(r[i0:i1], chi[i0:i1] / r[i0:i1], 2)
            m = r < 1e-5
            chi[m] = np.polyval(p, r[m]) * r[m]
            self.kb.append(dict(l=kb["l"], chi=chi, E=kb["E"]))
        assert all(k["l"] == 0 for k in self.kb), "solo proiettori s previsti"
        self._tables()
        self._fit_kb(kb_set)

    # ---------------------------------------------------------------- radial reference (the block itself)
    def radial_atom(self, config, det=None, chis=None):
        ext = {0: [(k["chi"], k["E"]) for k in self.kb]} if chis is None else {0: chis}
        return tr.Atom(self.Zv, config, grid=self.g, nuc=self.Vloc, extNL=ext, det=det)

    def levels(self, l, n=3, chis=None):
        a = self.radial_atom([(l + 1, l, 1)], chis=chis)
        w, U = tr.geig_low(a.h_matrix(l), self.g.r ** 2, n, -0.6 * self.Zv ** 2 - 20)
        P = [a._norm(U[:, i] * np.sqrt(self.g.r)) for i in range(n)]
        return w, P

    # ---------------------------------------------------------------- local part
    def _tables(self):
        r, h = self.g.r, self.g.h
        dV = self.Vloc + self.Zv / r
        s2 = np.exp(np.linspace(math.log(1e-10), math.log(1e6), 4000))
        sg = np.sqrt(s2)
        Q0 = np.empty(len(s2)); Q2 = np.empty(len(s2))
        for i0 in range(0, len(s2), 500):
            ss = s2[i0:i0 + 500, None]
            rho = 4 * math.pi * r ** 2 * (2 * math.pi * ss) ** -1.5 * np.exp(-r ** 2 / (2 * ss))
            Q0[i0:i0 + 500] = np.sum(rho * dV * r * h, axis=1)
            Q2[i0:i0 + 500] = np.sum(rho * dV * r ** 3 * h, axis=1)
        self._q0 = CubicSpline(np.log(s2), Q0 * sg)
        self._q2 = CubicSpline(np.log(s2), Q2 / sg)

    def G0G2(self, s2):
        sg = np.sqrt(s2); ls = np.log(s2)
        G0 = -self.Zv * SQ2PI / sg + self._q0(ls) / sg
        G2 = -self.Zv * 2 * sg * SQ2PI + self._q2(ls) * sg
        return G0, G2

    # ---------------------------------------------------------------- separable part as gaussians
    def _fit_kb(self, kb_set):
        r, h = self.g.r, self.g.h
        w = np.sqrt(r ** 3 * h)
        ref = self.levels(0)[0]
        self.kbfit = []
        for kb in self.kb:
            f = kb["chi"] / (r * math.sqrt(4 * math.pi))
            sets = [kb_set] if kb_set else KB_SETS
            for (nal, a0, a1) in sets:
                al = np.geomspace(a0, a1, nal)
                G = np.exp(-np.outer(r ** 2, al))
                c, *_ = np.linalg.lstsq(G * w[:, None], f * w, rcond=1e-16)
                chf = r * math.sqrt(4 * math.pi) * (G @ c)
                d = self.levels(0, chis=[(chf, kb["E"])])[0] - ref
                if np.max(np.abs(d)) <= KB_TOL or kb_set:
                    break
            self.kbfit.append(dict(alpha=al, c=c, E=kb["E"], dlev=d.tolist(), n=nal, range=(a0, a1),
                                   chi_fit=chf))


class BlockSystem(T.System):
    """valence electrons in the block field (no bare nucleus). spin: 'up1', 's2', 't2' (two up), 'd3', 'p3'"""

    def __init__(self, field, spin="s2", L=0, ee=True, name="", chunk=1500):
        if spin == "t2":
            super().__init__(0.0, "s2", L=L, m_kin=1.0, ee=ee, name=name)
            self.spin = "t2"
            self.perms = [(1.0, np.eye(2)), (-1.0, np.eye(2)[[1, 0]])]
        else:
            super().__init__(0.0, spin, L=L, m_kin=1.0, ee=ee, name=name)
        self.f = field
        self.chunk = chunk

    def me(self, A, u, B, v, parts=False):
        n = A.shape[0]
        if n > self.chunk:
            outs = [self.me(A[i:i + self.chunk], u[i:i + self.chunk], B[i:i + self.chunk], v[i:i + self.chunk], parts)
                    for i in range(0, n, self.chunk)]
            S = np.concatenate([o[0] for o in outs]); H = np.concatenate([o[1] for o in outs])
            if parts:
                return S, H, {k: np.concatenate([o[2][k] for o in outs]) for k in outs[0][2]}
            return S, H
        res = super().me(A, u, B, v, parts=parts)
        S, H = res[0], res[1]
        Vb = self.block_me(A, u, B, v)
        H = H + Vb
        if parts:
            P = res[2]; P["Ven"] = P["Ven"] + Vb
            return S, H, P
        return S, H

    def block_me(self, A, u, B, v):
        n, N = A.shape[0], self.N
        out = np.zeros(n)
        for c, Pm in self.perms:
            Bp = np.einsum("ji,njk,kl->nil", Pm, B, Pm)
            C = A + Bp
            Ci = np.linalg.inv(C)
            S0 = (np.pi ** N / np.linalg.det(C)) ** 1.5
            Sig = 0.5 * Ci
            vp = np.einsum("ji,nj->ni", Pm, v) if self.L else None
            loc = np.zeros(n)
            for i in range(N):
                s2 = Sig[:, i, i]
                G0, G2 = self.f.G0G2(s2)
                if self.L == 0:
                    loc += S0 * G0
                else:
                    Su = np.einsum("nij,nj->ni", Sig, u); Sv = np.einsum("nij,nj->ni", Sig, vp)
                    uSv = np.einsum("ni,ni->n", u, Sv)
                    gu = Su[:, i]; gv = Sv[:, i]
                    loc += S0 * ((uSv - gu * gv / s2) * G0 + gu * gv / s2 ** 2 * G2 / 3)
            nl = np.zeros(n)
            for kf in self.f.kbfit:
                for i in range(N):
                    nl += self._kb_one(A, u, Bp, vp, i, kf) / kf["E"]
            out += c * (loc + nl)
        return out

    def _kb_one(self, A, u, Bp, vp, i, kf):
        al, cc = kf["alpha"], kf["c"]
        N = self.N
        R = [j for j in range(N) if j != i]
        a = A[:, i, i][:, None] + al[None, :]               # (n, k)
        b = Bp[:, i, i][:, None] + al[None, :]
        pa = (np.pi / a) ** 1.5 * cc[None, :]
        pb = (np.pi / b) ** 1.5 * cc[None, :]
        if N == 1:
            if self.L:
                return np.zeros(A.shape[0])
            return pa.sum(1) * pb.sum(1)
        tA = A[:, i, R]; tB = Bp[:, i, R]                    # (n, m)
        AR = A[:, R][:, :, R]; BR = Bp[:, R][:, :, R]       # (n, m, m)
        MA = AR[:, None] - tA[:, None, :, None] * tA[:, None, None, :] / a[:, :, None, None]   # (n,k,m,m)
        MB = BR[:, None] - tB[:, None, :, None] * tB[:, None, None, :] / b[:, :, None, None]
        M = MA[:, :, None] + MB[:, None, :]                  # (n,k,k',m,m)
        m = len(R)
        if m == 1:
            det = M[..., 0, 0]
        else:
            det = M[..., 0, 0] * M[..., 1, 1] - M[..., 0, 1] * M[..., 1, 0]
        fac = (np.pi ** m / det) ** 1.5
        if self.L == 0:
            return np.einsum("nk,nl,nkl->n", pa, pb, fac)
        ut = u[:, R][:, None, :] - u[:, i][:, None, None] * tA[:, None, :] / a[:, :, None]     # (n,k,m)
        vt = vp[:, R][:, None, :] - vp[:, i][:, None, None] * tB[:, None, :] / b[:, :, None]
        if m == 1:
            q = ut[:, :, None, 0] * vt[:, None, :, 0] / M[..., 0, 0]
        else:
            inv = np.stack([np.stack([M[..., 1, 1], -M[..., 0, 1]], -1),
                            np.stack([-M[..., 1, 0], M[..., 0, 0]], -1)], -2) / det[..., None, None]
            q = np.einsum("nkx,nklxy,nly->nkl", ut, inv, vt)
        return np.einsum("nk,nl,nkl->n", pa, pb, 0.5 * q * fac)


def run_svm(sys_, K, Zs, seed=2026, sweeps=2, trials=30, Lmax=30.0, log=print, polish=(0.2, 0.08), nm=0):
    rng = np.random.default_rng(seed)
    sv = T.SVM(sys_, rng, Zs, Lmax=Lmax, trials=trials, polish=polish, npol=8, nm=nm)
    t0 = time.time()
    while len(sv.A) < K:
        e = sv.add()
        if len(sv.A) % 20 == 0:
            log(f"   K={len(sv.A)} E={e:.10f} ({time.time() - t0:.0f} s)")
    for s in range(sweeps):
        for k in range(len(sv.A)):
            sv.refine_one(k)
        e = sv.rebuild()
        log(f"   giro {s + 1}: E={e:.10f} ({time.time() - t0:.0f} s)")
    return sv


def canon_eig(Hm, Sm, thr=1e-13):
    w, V = np.linalg.eigh(Sm)
    keep = w > thr * w.max()
    X = V[:, keep] / np.sqrt(w[keep])
    Ev, Y = np.linalg.eigh(X.T @ Hm @ X)
    return Ev, X @ Y


def check(Z):
    t0 = time.time()
    f = BlockField(Z)
    out = dict(Z=Z, kbfit=[dict(n=k["n"], range=k["range"], dlev_radiale=k["dlev"]) for k in f.kbfit])
    ws, Ps = f.levels(0); wp, Pp = f.levels(1)
    out["radiale"] = dict(s=ws.tolist(), p=wp.tolist())
    for lab, L, ref in (("s", 0, ws[0]), ("p", 1, wp[0])):
        sysx = BlockSystem(f, "up1", L=L, name=f"Z{Z} 1e {lab}")
        # one electron: the ECG basis is a 1-parameter family; deterministic even-tempered sets instead of SVM
        Es = []
        for K, a0, a1 in ((40, 0.005, 3e4), (50, 0.005, 1e5), (60, 0.003, 1e5), (70, 0.002, 2e5), (80, 0.001, 3e5)):
            al = np.geomspace(a0, a1, K)
            sv = T.SVM(sysx, np.random.default_rng(0), 1.0)
            sv.A = al[:, None, None].copy(); sv.u = np.ones((K, 1))
            sv.rebuild()
            sv._eig = canon_eig                       # even-tempered sets are near-dependent: canonical orthog.
            Es.append(float(canon_eig(sv.Hmat, sv.Smat)[0][0]))
        out[lab] = dict(E_ecg=Es, E_radiale=float(ref), dE=[e - float(ref) for e in Es],
                        PASS=bool(abs(Es[-1] - ref) <= 1e-6))
        # orbital: ECG one-electron function on the radial grid vs the block orbital
        E, c = sv.state()
        r = f.g.r
        al = np.array([A[0, 0] for A in sv.A])
        if L == 0:
            psi = np.exp(-np.outer(r ** 2, al)) @ c                      # psi(r) (s)
            P = psi * r * math.sqrt(4 * math.pi)
        else:
            # (u z) e^{-a r^2}: radial part r e^{-a r^2} times z/r -> P(r) = r^2 e^{-a r^2} * sqrt(4 pi/3)
            P = (np.exp(-np.outer(r ** 2, al)) * (r ** 2)[:, None]) @ (c * sv.u[:, 0]) * math.sqrt(4 * math.pi / 3)
        P = P / math.sqrt(f.g.integ(P * P))
        Pref = Ps[0] if L == 0 else Pp[0]
        ov = f.g.integ(P * Pref)
        out[lab]["1-|overlap|"] = float(1 - abs(ov))
        out[lab]["max|dP|"] = float(np.max(np.abs(np.sign(ov) * P - Pref)))
    out["secondi"] = time.time() - t0
    js = os.path.join(OUT, "tavg_pecg_check.json")
    allr = json.load(open(js)) if os.path.exists(js) else {}
    allr[str(Z)] = out
    json.dump(allr, open(js, "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "check":
        for z in sys.argv[2:]:
            check(int(z))
