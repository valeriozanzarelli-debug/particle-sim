"""aggfull_run.py - production jobs of AGGLOMERATORE COMPLETO (configurations fixed in lab/out/aggfull_prereg.json).
One job per call (each <= ~9 minutes; a job that does not converge saves its vector and continues at the next call).
Usage: timeout 590 python lab/aggfull_run.py JOB      (python lab/aggfull_run.py list  -> the jobs and their state)
Results: lab/out/aggfull_atom.json (vectors lab/out/aggfull_vec_<job>.npy, Galerkin fields cached as .npy: not committed)."""
import os, sys, time, json, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import aggfull_atom as aa

DEV = aa.DEV
OUT = ac.OUT
PFX = "dev_" if DEV else ""
# amendment runs (aggfull_amend.json, aggfull_amend2.json): jobs whose name has "K" in second place use the re-learned
# kinetic (aggfull_relearn.pkl); they write to their own file aggfull_atomK.json
# phase 4: second letter "T" = the Taylor-pinned + parity kinetic (aggfull_relearnT.pkl), file aggfull_atomT.json
KMODE = len(sys.argv) > 1 and all(len(a) > 1 and a[1] in "KT" for a in sys.argv[1:])
KL = sys.argv[1][1] if KMODE else ""
JS = os.path.join(OUT, PFX + {"K": "aggfull_atomK.json", "T": "aggfull_atomT.json", "": "aggfull_atom.json"}[KL])
res = json.load(open(JS)) if os.path.exists(JS) else {}
TH = ac.TH_A
ALPHA = ac.zz.ALPHA

# job: (level B of chain A, N, state, recoil, depends-on (orthogonal to), restriction)
JOBS = {
    "e16_1s": (256, 256, "1s", False, None, "learned"),
    "p16_1s": (256, 256, "1s", True, None, "learned"),
    "e8_1s": (512, 256, "1s", False, None, "learned"),
    "e8_2s": (512, 256, "2s", False, "e8_1s", "learned"),
    "e8_2p": (512, 256, "2p", False, None, "learned"),
    "p8_1s": (512, 256, "1s", True, None, "learned"),
    "p8_2s": (512, 256, "2s", True, "p8_1s", "learned"),
    "p8_2p": (512, 256, "2p", True, None, "learned"),
    "e8_1s_N128": (512, 128, "1s", False, None, "learned"),            # box check (+- 8 a0 vs +- 16 a0)
    "e8_1s_N192": (512, 192, "1s", False, None, "learned"),
    "e8_2s_N192": (512, 192, "2s", False, "e8_1s_N192", "learned"),     # box check for 2s (+- 12 a0, the old box)
    "ca8_1s_N192": (512, 192, "1s", False, None, "cellavg"),            # the OLD hand-written restriction, old box
    "ca8_2s_N192": (512, 192, "2s", False, "ca8_1s_N192", "cellavg"),
    "ca8_1s": (512, 256, "1s", False, None, "cellavg"),                 # the old restriction in the large box
    "ca8_2s": (512, 256, "2s", False, "ca8_1s", "cellavg"),
    "e8_2p_N192": (512, 192, "2p", False, None, "learned"),           # phase 4: float64 field (N <= 192), box +-12 a0
    "e8_2p_c": (512, 256, "2p", False, None, "learned"),               # phase 5: 'clean' solver (parity penalty + projection)
}
if DEV:      # development: same jobs on small grids (pipeline test only)
    DN = int(os.environ.get("AGGFULL_DEVN", "96"))
    JOBS = {k: (v[0], min(v[1], DN)) + v[2:] for k, v in JOBS.items()}
PERT = {"q16_1s": "e16_1s", "q8_1s": "e8_1s", "q8_2s": "e8_2s", "q8_2p": "e8_2p", "q8_2p_N192": "e8_2p_N192", "q8_2p_c": "e8_2p_c"}
for _L in "KT":
    for _j in ("e16_1s", "p16_1s", "e8_1s", "e8_2s", "e8_2p", "p8_1s", "p8_2s", "p8_2p", "e8_2p_N192", "e8_2p_c"):
        _v = JOBS[_j]
        JOBS[_j[0] + _L + _j[1:]] = _v[:4] + ((_v[4][0] + _L + _v[4][1:]) if _v[4] else None,) + _v[5:]
    for _q, _s in [x for x in list(PERT.items()) if x[0][1] not in "KT"]:
        PERT[_q[0] + _L + _q[1:]] = _s[0] + _L + _s[1:]
MOVE = {"m8_1s": (512, 128)}


def level(B):
    if KMODE:
        fn = {"K": "aggfull_relearn.pkl", "T": "aggfull_relearnT.pkl"}[KL]
        return [l for l in pickle.load(open(os.path.join(OUT, fn), "rb")) if l["B"] == B][0]
    return [l for l in ac.load_chain_A() if l["B"] == B][0]


def proton_level(B):
    st = pickle.load(open(os.path.join(OUT, PFX + "aggfull_proton.pkl"), "rb"))
    Bp = {256: 2 ** 19, 512: 2 ** 20}[B]
    return [l for l in st["levels"] if l["B"] == Bp][0]


def coupling(B):
    st = pickle.load(open(os.path.join(OUT, PFX + "aggfull_learn.pkl"), "rb"))
    s = [l for l in st if l["B"] == B][0]
    C, offs = s["C"], s["offs"]
    g0 = complex(np.trace(ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]) / 2)
    return C, offs, g0


def field(N, h, kind="learned"):
    if kind == "cellavg":
        V, x = ac.coulomb_cellavg(N, h)

        class F:
            pass
        f = F()
        f.x = x
        f.apply = lambda v: V * v
        f.apply_k = lambda Fv: ac.fftn(V * ac.ifftn(Fv))
        return f
    fn = os.path.join(OUT, f"aggfull_field_N{N}_h{h:.6f}.npz")
    if os.path.exists(fn):
        d = np.load(fn)
        f = ac.GalerkinField.__new__(ac.GalerkinField)
        f.N, f.h, f.x = N, h, d["x"]
        f.subs = [(sx, sy, sz) for sx in (0, 1) for sy in (0, 1) for sz in (0, 1)]
        f.V = {s: d[f"V{s[0]}{s[1]}{s[2]}"] for s in f.subs}
        f.k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
        f.ph1 = np.exp(1j * f.k1 * h / 2)
        return f
    t0 = time.time()
    f = ac.GalerkinField(N, h)
    np.savez(fn, x=f.x, **{f"V{s[0]}{s[1]}{s[2]}": f.V[s] for s in f.subs})
    print(f"  Galerkin field N={N} built in {time.time() - t0:.0f} s", flush=True)
    return f


def save():
    json.dump(res, open(JS, "w"), indent=1, default=float)


def vecfile(job):
    return os.path.join(OUT, f"{PFX}aggfull_vec_{job}.npy")


def run_job(job):
    B, N, kind, recoil, dep, restr = JOBS[job]
    if res.get(job, {}).get("converged"):
        print(f"{job}: done already", flush=True)
        return
    t0 = time.time()
    h = B * TH
    rec = level(B)
    T, odd, spin = aa.kinetic(rec, N, h)
    info = dict(B=B, N=N, h=h, cell_over_a0=h / ac.A0, box_a0=N * h / ac.A0, state=kind, recoil=recoil, restriction=restr,
                odd_part_e=odd, spin_part_e=spin)
    if recoil:
        prec = proton_level(B)
        Tp, oddp, spinp = aa.kinetic(prec, N, h)
        T = T + Tp                                  # T_p(-k) = T_p(k) (parity-even part)
        info.update(odd_part_p=oddp, spin_part_p=spinp)
        del Tp
    C, offs, g0 = coupling(B)
    gam = g0.real if restr == "learned" else 1.0
    info["gamma0"] = [g0.real, g0.imag]
    f = field(N, h, restr)
    at = aa.Atom(T, f, gamma0=gam, real=True, sym="oddz" if kind == "2p" else "even", clean=job.endswith("_c"))
    Y = None
    if dep is not None:
        if not res.get(dep, {}).get("converged"):
            print(f"{job}: needs {dep} first", flush=True)
            return
        Y = np.load(vecfile(dep))[:, None]
    warm = vecfile(job.replace(KL, "", 1)) if KMODE else None      # K jobs start from the converged original vector
    if job.endswith("_c") and not os.path.exists(vecfile(job)):    # phase 5: start from the same job without '_c' (projected)
        w2 = vecfile(job[:-2])
        warm = w2 if os.path.exists(w2) else warm
    x0 = (np.load(vecfile(job)) if os.path.exists(vecfile(job)) else
          np.load(warm) if warm and os.path.exists(warm) else aa.guess(f.x, kind))
    budget = 560 - (time.time() - t0)
    it_total = res.get(job, {}).get("iterations", 0)
    conv = False
    while time.time() - t0 < budget - 150:
        v, lres, it = aa.solve(at, x0, Y=Y, maxiter=5)
        x0 = v
        it_total += it
        E, r = at.energy(v)
        np.save(vecfile(job), v)
        print(f"  {job}: it {it_total}: E = {E * ac.MC2:.8f} eV, residual {r:.2e} ({time.time() - t0:.0f} s)", flush=True)
        if r < 2e-8:      # (was 2e-9; float32 field: residual floor ~1e-8, energy error ~r^2/gap < 1e-12 eV)
            conv = True
            break
    n = 1 if kind == "1s" else 2
    E, r = at.energy(x0)
    info.update(E_eV=E * ac.MC2, residual=r, iterations=it_total, converged=conv, schroedinger_eV=-ac.RY / n ** 2,
                rel_to_schroedinger=E * ac.MC2 / (-ac.RY / n ** 2) - 1, seconds=res.get(job, {}).get("seconds", 0) + time.time() - t0,
                applications=at.napply)
    res[job] = info
    save()
    print(f"{job}: E = {E * ac.MC2:.8f} eV ({100 * info['rel_to_schroedinger']:+.5f} % vs Schroedinger), residual {r:.1e}, "
          f"converged {conv}", flush=True)


def run_pert(job):
    src = PERT[job]
    if job in res:
        print(f"{job}: done already")
        return
    if not res.get(src, {}).get("converged"):
        print(f"{job}: needs {src}")
        return
    t0 = time.time()
    B, N, kind, recoil, dep, restr = JOBS[src]
    h = B * TH
    rec = level(B)
    T, odd, spin = aa.kinetic(rec, N, h)
    mk, mk_all = aa.kinetic_mass(rec, h)
    C, offs, g0 = coupling(B)
    f = field(N, h)
    v = np.load(vecfile(src)).reshape((N,) * 3)
    v = v / np.linalg.norm(v)
    if kind == "2p":
        pz = v
        px = np.transpose(pz, (2, 1, 0))
        py = np.transpose(pz, (0, 2, 1))
        phis = [px, py, pz]
    else:
        phis = [v]
    S = aa.coupling_blocks(f, phis, offs)
    M = aa.first_order_matrix(C, offs, S, g0.real)
    ev, evec = np.linalg.eigh(M)
    Rk = aa.kinetic_remainder(T, h, v, mk)
    out = dict(source=src, m_kin_e=mk, m_kin_e_dirs=mk_all, kinetic_remainder_eV=Rk * ac.MC2,
               dV_eigen_eV=(ev * ac.MC2).tolist(), dV_matrix_re=M.real.tolist(), dV_matrix_im=M.imag.tolist(),
               seconds=time.time() - t0)
    if kind == "2p":
        # spin-orbit content: decompose on L.S (j = 1/2: 2 states, j = 3/2: 4 states)
        out["splitting_eV"] = float((ev[2:].mean() - ev[:2].mean()) * ac.MC2)
        out["lower_pair_spread_eV"] = float((ev[1] - ev[0]) * ac.MC2)
        out["upper_quartet_spread_eV"] = float((ev[5] - ev[2]) * ac.MC2)
        out["interval_ratio"] = float((ev[2:].mean() - ev.mean()) / (ev.mean() - ev[:2].mean()))
    res[job] = out
    save()
    print(f"{job}: m_kin {mk:.6f}; <T - k^2/2m> = {Rk * ac.MC2:.4e} eV; dV eigenvalues (eV) " +
          " ".join(f"{e * ac.MC2:.5e}" for e in ev) + f"  ({time.time() - t0:.0f} s)", flush=True)


def run_move(job):
    """the atom moving as a whole: ground state of H_P(k) = T_e(k + P/2) + T_p(-k + P/2) + V (split 1/2: no masses used)
    for P = n x 4 pi/L (so that P/2 is a grid vector: exact relabelling, no twist), n = 0, 1, 2"""
    if job in res and len(res[job].get("E_eV", [])) == 3:
        print(f"{job}: done already")
        return
    B, N = MOVE[job]
    t0 = time.time()
    h = B * TH
    rec, prec = level(B), proton_level(B)
    Te, _, _ = aa.kinetic(rec, N, h)
    Tp, _, _ = aa.kinetic(prec, N, h)
    C, offs, g0 = coupling(B)
    f = field(N, h)
    out = res.get(job, dict(P=[], E_eV=[], residual=[]))
    for n in range(len(out["E_eV"]), 3):
        sh = (n, 0, 0)                                    # P/2 = n * 2 pi / L along x
        # T_e(k + P/2) + T_p(-k + P/2) = T_e(k + P/2) + T_p(k - P/2)  (T_p parity-even)
        T = np.roll(Te, shift=-n, axis=0) + np.roll(Tp, shift=n, axis=0)
        at = aa.Atom(T, f, gamma0=g0.real, real=False)
        v, lres, it = aa.solve(at, aa.guess(f.x, "1s").astype(complex), maxiter=40)
        E, r = at.energy(v)
        P = 2 * n * 2 * np.pi / (N * h)
        out["P"].append(P)
        out["E_eV"].append(E * ac.MC2)
        out["residual"].append(r)
        res[job] = out
        save()
        print(f"  {job}: P = {P:.6e} (1/lambda_C) = {P * ac.A0:.4f}/a0: E = {E * ac.MC2:.10f} eV, residual {r:.1e} "
              f"({time.time() - t0:.0f} s)", flush=True)
        if time.time() - t0 > 400:
            break
    if len(out["E_eV"]) == 3:
        P = np.array(out["P"]); E = np.array(out["E_eV"]) / ac.MC2
        c = np.polyfit(P ** 2, E, 1)[0]
        out["M_atom"] = float(1 / (2 * c))
        out["M_over_me_plus_mp"] = float(1 / (2 * c) / (1 + ac.MP))
        res[job] = out
        save()
        print(f"{job}: M_atom/(m_e + m_p) = {out['M_over_me_plus_mp']:.6f}")


if __name__ == "__main__":
    for job in sys.argv[1:]:
        if job == "list":
            for j in list(JOBS) + list(PERT) + list(MOVE):
                print(j, "converged" if res.get(j, {}).get("converged") or (j in PERT and j in res) else "-")
        elif job in JOBS:
            run_job(job)
        elif job in PERT:
            run_pert(job)
        elif job in MOVE:
            run_move(job)
