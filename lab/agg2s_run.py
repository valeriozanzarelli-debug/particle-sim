"""agg2s_run.py - production jobs of "2s all'agglomeratore" (configurations fixed in lab/out/agg2s_prereg.json).
The coarse hydrogen atom of AGGLOMERATORE COMPLETO (learned kinetic K = aggfull_relearn.pkl or T = aggfull_relearnT.pkl,
learned coupling gamma0 = 1, alias-free Galerkin field of -alpha/r) in the box L = 32 a0 (N = 256, cells a0/8), with
TWISTED boundary conditions (agg2s_core.TwistAtom): tw = 000 periodic (the aggfull_ runs), 111 anti-periodic, 100, 110.
Kinetic 'S' = the exact Schroedinger k^2/2 (m = 1): reference for the cell (basis) error.
Usage: timeout 590 python lab/agg2s_run.py JOB [JOB ...]   (list -> jobs and state). Each call <= ~9 min; a job that
has not converged saves its vector and continues at the next call. Results: lab/out/agg2s_atom.json; vectors
lab/out/agg2s_vec_<job>.npy (not committed)."""
import os, sys, time, json, pickle
import numpy as np
import agg2s_core as c

DEV = os.environ.get("AGG2S_DEV") == "1"
PFX = "dev_" if DEV else ""
JS = os.path.join(c.OUT, PFX + "agg2s_atom.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
TOL = 2e-8                       # residual criterion (aggfull_run: float32 field floor ~1e-8)

# job -> (kinetic, B (cell a0/16: 256, a0/8: 512), N, state, twist, deflation job, warm-start file)
def _jobs():
    J = {}
    agg = lambda nm: os.path.join(c.OUT, f"aggfull_vec_{nm}.npy")
    for kin in "KTS":
        for tw in ("000", "111", "100", "110"):
            if kin != "K" and tw in ("100", "110"):
                continue
            J[f"{kin}1s_{tw}"] = (kin, 512, 256, "1s", tw, None, agg("eK8_1s"))
            J[f"{kin}2s_{tw}"] = (kin, 512, 256, "2s", tw, f"{kin}1s_{tw}", agg("eK8_2s"))
            if tw in ("000", "111"):
                J[f"{kin}2p_{tw}"] = (kin, 512, 256, "2p", tw, None, agg("eK8_2p"))
    J["T1s_000_a16"] = ("T", 256, 256, "1s", "000", None, agg("eK16_1s"))
    J["K1s_000_a16"] = ("K", 256, 256, "1s", "000", None, agg("eK16_1s"))
    return J


JOBS = _jobs()
if DEV:
    DN = int(os.environ.get("AGG2S_DEVN", "64"))
    JOBS = {k: (v[0], v[1], DN) + v[3:6] + (None,) for k, v in JOBS.items()}


def level(kin, B):
    fn = {"K": "aggfull_relearn.pkl", "T": "aggfull_relearnT.pkl"}[kin]
    return [l for l in pickle.load(open(os.path.join(c.OUT, fn), "rb")) if l["B"] == B][0]


def taylor_iso(rec, h):
    """isotropic Taylor coefficients of the learned spin-blind band T(k) = k^2/2m + c4 k^4 + ...:
    1/2m = -(h/6) sum_j Re s_j |j|^2 / h^0 ... (cell units -> physical), c4 = (h^3/120) sum_j Re s_j |j|^4"""
    K, offs = rec["K"], np.asarray(rec["offs"], float)
    s = (0.5 * (K[0, 0] + K[1, 1])).real
    j2 = np.sum(offs ** 2, 1)
    inv2m = -(h / 6) * np.sum(s * j2)
    c4 = (h ** 3 / 120) * np.sum(s * j2 ** 2)
    return float(1 / (2 * inv2m)), float(c4)


def vecfile(job):
    return os.path.join(c.OUT, f"{PFX}agg2s_vec_{job}.npy")


def run_job(job):
    kin, B, N, kind, tws, dep, warm = JOBS[job]
    if res.get(job, {}).get("converged"):
        print(f"{job}: done already", flush=True)
        return
    t0 = time.time()
    tw = tuple(int(x) for x in tws)
    h = B * c.TH
    if kin == "S":
        T = c.kin_schroedinger(N, h, tw)
        m_kin, c4 = 1.0, 0.0
    else:
        rec = level(kin, B)
        T, spin = c.kin_learned(rec, N, h, tw)
        m_kin, c4 = taylor_iso(rec, h)
    f = c.galerkin_field(N, h)
    at = c.TwistAtom(T, f, tw, gamma0=1.0, sym="oddz" if kind == "2p" else "even")
    Y = None
    if dep is not None:
        if not res.get(dep, {}).get("converged"):
            print(f"{job}: needs {dep} first", flush=True)
            return
        Y = np.load(vecfile(dep))[:, None]
    if os.path.exists(vecfile(job)):
        x0 = np.load(vecfile(job))
    elif warm and os.path.exists(warm):
        x0 = at.project(np.load(warm).reshape((N,) * 3)).ravel()
        x0 /= np.linalg.norm(x0)
    else:
        x0 = c.guess(at, kind)
    prev = res.get(job, {})
    budget = 540 - (time.time() - t0)
    log = lambda s: print(f"  {job}{s}", flush=True)
    v, E, r, it, ok = c.lobpcg1(at, x0, Y=Y, tol_res=TOL, budget=max(budget - 120, 60), log=log)
    np.save(vecfile(job), v)
    Tk, Vp = at.parts(v)
    n = int(kind[0])
    vv = v.reshape((N,) * 3)
    h_a0 = h / c.A0
    nrm = float(np.sum(vv ** 2)) * h_a0 ** 3
    rho = float(vv[N // 2, N // 2, N // 2] ** 2 / nrm) if kind != "2p" else 0.0
    F = c.sfft.fftn(vv, workers=1)
    w = np.empty(F.shape)
    np.multiply(F.real, F.real, out=w)
    w += F.imag ** 2
    del F
    ka = c.kgrid_axes(N, h_a0, tw)
    kyz2 = ka[1][:, None] ** 2 + ka[2][None, :] ** 2
    sw = sk4 = srem = 0.0
    for i in range(N):                                   # slice by slice (memory)
        k2 = ka[0][i] ** 2 + kyz2
        wi = w[i]
        sw += wi.sum()
        sk4 += np.sum(wi * k2 ** 2)
        if kin != "S":
            srem += np.sum(wi * (T[i] - k2 / c.A0 ** 2 / (2 * m_kin)))
    del w
    k4 = float(sk4 / sw)
    krem = float(srem / sw) * c.MC2      # <T - k^2/2m> (relativistic part of the learned band), eV
    info = dict(kin=kin, B=B, N=N, h_a0=h_a0, box_a0=N * h_a0, state=kind, twist=tws, E_eV=E * c.MC2, residual=r,
                iterations=prev.get("iterations", 0) + it, converged=bool(ok), T_eV=Tk * c.MC2, V_eV=Vp * c.MC2,
                rho_c_a0=rho, k4_c_a0=k4, m_kin=m_kin, c4=c4, kin_remainder_eV=krem,
                seconds=prev.get("seconds", 0) + time.time() - t0)
    res[job] = info
    json.dump(res, open(JS, "w"), indent=1)
    print(f"{job}: E = {E * c.MC2:.9f} eV, residual {r:.2e}, converged {ok}, rho_c {rho:.6g}, <k^4> {k4:.6g}, "
          f"kin.rem {krem:.4e} eV ({time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    for job in sys.argv[1:]:
        if job == "list":
            for j in JOBS:
                print(j, "converged" if res.get(j, {}).get("converged") else ("started" if j in res else "-"))
        else:
            run_job(job)
