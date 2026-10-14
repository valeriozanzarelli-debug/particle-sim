"""tav_dev.py - controlli di MACCHINA (sviluppo, dichiarati: nessun sistema d'esame) per tav_ecg.py.
  mc      : elementi di matrice S, T, V_en, V_ee per funzioni P (l = 1) e dipolo, contro Monte Carlo (gaussiana come
            peso d'importanza), + simmetria bra/ket
  hyd     : idrogeno 1s (-0,5), 2p (-0,125), dipolo 1s-2p <z> = 128 sqrt2/243, f = 0,4162 (forme lunghezza/velocita')
  nonint  : 3 e 4 elettroni senza repulsione (Z = 3, 4): 1s^2 2s e 1s^2 2p (-10,125), 1s^2 2s^2 e 1s^2 2s 2p 1P (-20);
            versione 'bosonica' piu' bassa (l'antisimmetria lavora)
  he      : elio 1S (gia' validato dal thread principale: qui solo per vedere che la macchina a pesi di spin lo rifa')
  timing  : costo per elemento di matrice a 3 e 4 elettroni
Uscita: lab/out/tav_dev.json"""
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import tav_ecg as T

JS = os.path.join(T.OUT, "tav_dev.json")


def mc_check():
    rng = np.random.default_rng(7)
    recs = []
    for N, spin in ((3, "d3"), (4, "s4")):
        s = T.System(3.0, spin, L=1, m_kin=1.0)
        s.perms = [(1.0, np.eye(N))]                      # solo identita': controllo delle formule
        A, u = T.random_fn(rng, N, 3.0, Lmax=6.0, L=1)
        B, v = T.random_fn(rng, N, 3.0, Lmax=6.0, L=1)
        S, H, P = s.me(A[None], u[None], B[None], v[None], parts=True)
        S2, H2 = s.me(B[None], v[None], A[None], u[None])
        C = A + B; Ci = np.linalg.inv(C); S0 = (np.pi ** N / np.linalg.det(C)) ** 1.5
        L = np.linalg.cholesky(Ci / 2)
        acc = np.zeros(5); acc2 = np.zeros(5); n = 0
        for ch in range(20):
            z = rng.normal(size=(200000, N, 3))
            x = np.einsum("ij,njk->nik", L, z)
            fu = x[:, :, 2] @ u; gv = x[:, :, 2] @ v
            Ax = np.einsum("ij,njk->nik", A, x); Bx = np.einsum("ij,njk->nik", B, x)
            # grad f / e^-xAx = u_k e_z - 2 (u.x_z)(Ax)_k
            gf = -2 * fu[:, None, None] * Ax; gf[:, :, 2] += u
            gg = -2 * gv[:, None, None] * Bx; gg[:, :, 2] += v
            tk = np.sum(gf * gg, axis=(1, 2)) / 2.0
            ven = -3.0 * np.sum(1 / np.linalg.norm(x, axis=2), axis=1) * fu * gv
            vee = sum(1 / np.linalg.norm(x[:, i] - x[:, j], axis=1) for i in range(N) for j in range(i + 1, N)) * fu * gv
            dz = np.sum(x[:, :, 2], axis=1) * gv                # <f0|sum z|g>, f0 = exp(-xAx)
            arr = np.stack([fu * gv, tk, ven, vee, dz], 1)
            acc += arr.sum(0); acc2 += (arr ** 2).sum(0); n += len(arr)
        mc = S0 * acc / n
        err = S0 * np.sqrt(acc2 / n - (acc / n) ** 2) / np.sqrt(n)
        sS = T.System(3.0, spin, L=0, m_kin=1.0); sS.perms = [(1.0, np.eye(N))]
        Dl, Dv = T.dipole_me(sS, s, A[None], B[None], v[None])
        an = np.array([S[0], P["T"][0], P["Ven"][0], P["Vee"][0], Dl[0]])
        z = (an - mc) / err
        recs.append(dict(N=N, analytic=an.tolist(), mc=mc.tolist(), z=z.tolist(), herm_S=float(S[0] - S2[0]),
                         herm_H=float(H[0] - H2[0])))
        print("MC", N, np.round(z, 2), S[0] - S2[0], H[0] - H2[0], flush=True)
    ok = all(max(abs(np.array(r["z"]))) < 4 and abs(r["herm_H"]) < 1e-10 * max(1, abs(r["analytic"][1])) for r in recs)
    return dict(recs=recs, PASS=bool(ok))


def hydrogen():
    out = {}
    rng = np.random.default_rng(3)
    sS = T.System(1.0, "up1", L=0, m_kin=1.0)
    s1, EK, _ = run(sS, 30, 1.0, rng, sweeps=2)
    sP = T.System(1.0, "up1", L=1, m_kin=1.0)
    s2, EK2, _ = run(sP, 30, 1.0, rng, sweeps=2)
    E1, c1 = s1.state(); E2, c2 = s2.state()
    I, J = np.meshgrid(np.arange(len(c1)), np.arange(len(c2)), indexing="ij")
    Dl, Dv = T.dipole_me(sS, sP, s1.A[I.ravel()], s2.A[J.ravel()], s2.u[J.ravel()])
    zl = float(c1 @ Dl.reshape(I.shape) @ c2); zv = float(c1 @ Dv.reshape(I.shape) @ c2)
    dE = E2 - E1
    out = dict(E1s=E1, E2p=E2, z=abs(zl), z_exact=128 * np.sqrt(2) / 243, f_len=2 * dE * zl ** 2, f_vel=2 * zv ** 2 / dE,
               f_exact=0.41620)
    out["PASS"] = bool(abs(E1 + 0.5) < 1e-6 and abs(E2 + 0.125) < 1e-6 and abs(out["f_len"] / 0.41620 - 1) < 1e-3
                       and abs(out["f_vel"] / 0.41620 - 1) < 1e-3)
    print("H", out, flush=True)
    return out


def run(sys_, K, Zs, rng, sweeps=1, Lmax=30.0, trials=30, log=None):
    sv = T.SVM(sys_, rng, Zs, Lmax, trials)
    EK = []
    while len(sv.A) < K:
        e = sv.add()
        if len(sv.A) % 10 == 0:
            EK.append((len(sv.A), e))
            if log:
                log(f"  K={len(sv.A)} E={e:.10f}")
    for sw in range(sweeps):
        for k in range(len(sv.A)):
            sv.refine_one(k)
        sv.E = float(sv._eig(sv.Hmat, sv.Smat)[0][0])
    return sv, EK, sv.E


def nonint():
    out = {}
    rng = np.random.default_rng(5)
    for name, Z, spin, L, ex in (("Li_1s2_2s", 3.0, "d3", 0, -10.125), ("Li_1s2_2p", 3.0, "d3", 1, -10.125),
                                 ("Be_1s2_2s2", 4.0, "s4", 0, -20.0), ("Be_1s2_2s2p_1P", 4.0, "s4", 1, -20.0)):
        t0 = time.time()
        s = T.System(Z, spin, L=L, m_kin=1.0, ee=False)
        sv, EK, E = run(s, 60, Z, rng, sweeps=1)
        sb = T.System(Z, spin, L=L, m_kin=1.0, ee=False, bosonic=True)
        svb, _, Eb = run(sb, 30, Z, rng, sweeps=0)
        out[name] = dict(E=E, exact=ex, err=E - ex, E_bosonic=Eb, sec=time.time() - t0,
                         PASS=bool(0 <= E - ex < 2e-3 and Eb < ex - 0.5))
        print(name, out[name], flush=True)
    return out


def helium():
    rng = np.random.default_rng(11)
    s = T.System(2.0, "s2", L=0, m_kin=1.0)
    t0 = time.time()
    sv, EK, E = run(s, 100, 2.0, rng, sweeps=1)
    out = dict(E=E, ref=-2.903724377, err=E + 2.903724377, sec=time.time() - t0, PASS=bool(abs(E + 2.903724377) < 2e-5))
    print("He", out, flush=True)
    return out


def timing():
    rng = np.random.default_rng(1)
    out = {}
    for N, spin, L in ((3, "d3", 0), (3, "d3", 1), (4, "s4", 0), (4, "s4", 1)):
        s = T.System(3.0, spin, L=L)
        fs = [T.random_fn(rng, N, 3.0, L=L) for _ in range(2000)]
        A = np.array([f[0] for f in fs]); u = np.array([f[1] for f in fs])
        t0 = time.time()
        s.me(A, u, A[::-1].copy(), u[::-1].copy())
        dt = (time.time() - t0) / 2000
        out[f"N{N}_L{L}"] = dict(us_per_me=dt * 1e6, nperm=len(s.perms))
        print(out, flush=True)
    return out


if __name__ == "__main__" and sys.argv[1:] != ["cross"]:
    res = T.load_json(JS)
    what = sys.argv[1:] or ["timing", "mc", "hyd", "nonint", "he"]
    for w in what:
        f = dict(mc=mc_check, hyd=hydrogen, nonint=nonint, he=helium, timing=timing)[w]
        res[w] = f()
        res["nota"] = "controlli di sviluppo della macchina (nessun sistema d'esame), dichiarati nella pre-registrazione"
        T.save_json(JS, res)


def fval(svS, svP):
    """f_L, f_V tra gli stati piu' bassi di due SVM (S e P)"""
    E1, c1 = svS.state(); E2, c2 = svP.state()
    KS, KP = len(c1), len(c2)
    I, J = np.meshgrid(np.arange(KS), np.arange(KP), indexing="ij")
    Dl, Dv = T.dipole_me(svS.s, svP.s, svS.A[I.ravel()], svP.A[J.ravel()], svP.u[J.ravel()])
    zl = float(c1 @ Dl.reshape(I.shape) @ c2); zv = float(c1 @ Dv.reshape(I.shape) @ c2)
    dE = E2 - E1
    return dict(E_S=E1, E_P=E2, dE=dE, f_L=2 * dE * zl ** 2, f_V=2 * zv ** 2 / (svS.s.m ** 2 * dE))


def cross_spin():
    """pesi di spin incrociati (non esame): ione Li-like Z = 5 con 2S in p3 e in d3, 2P in d3 e in p3;
    ione Be-like Z = 6 con 1S in p4 e s4, 1P in s4"""
    out = {}
    kw = dict(Lmax=20.0, trials=20, mut=0.5, p_pair=0.4, polish=(0.3, 0.1), npol=10, nm=60)
    def grow(Z, spin, L, K, seed):
        sv = T.SVM(T.System(Z, spin, L=L, m_kin=1.0), np.random.default_rng(seed), Z, **kw)
        while len(sv.A) < K:
            sv.add()
        return sv
    t0 = time.time()
    S3p = grow(5.0, "p3", 0, 100, 1); S3d = grow(5.0, "d3", 0, 100, 2); P3d = grow(5.0, "d3", 1, 100, 3); P3p = grow(5.0, "p3", 1, 100, 4)
    out["Z5"] = {k: fval(a, b) for k, (a, b) in dict(pp=(S3p, P3p), pd=(S3p, P3d), dd=(S3d, P3d), dp=(S3d, P3p)).items()}
    print(out["Z5"], time.time() - t0, flush=True)
    S4p = grow(6.0, "p4", 0, 80, 5); S4s = grow(6.0, "s4", 0, 80, 6); P4s = grow(6.0, "s4", 1, 60, 7)
    out["Z6"] = {k: fval(a, b) for k, (a, b) in dict(ps=(S4p, P4s), ss=(S4s, P4s)).items()}
    print(out["Z6"], time.time() - t0, flush=True)
    return out


if __name__ == "__main__" and sys.argv[1:] == ["cross"]:
    res = T.load_json(JS); res["cross"] = cross_spin(); T.save_json(JS, res)
