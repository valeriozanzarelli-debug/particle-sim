"""aggfull_p6_atom.py - phase 6, coarse-level alpha^4 exam in 3-D: the T level (aggfull_relearnT.pkl: Taylor-pinned +
parity kinetic) + the coupling learned in the LARGE-COMPONENT language (aggfull_p6_learn.pkl), first order on the
converged T-level vectors (Hermitian H0 = T + g0 V_G: left = right eigenvector, so the first-order shift of the
non-Hermitian LC level is <v|dV_LC|v> exactly; for the 2p the 6x6 block is diagonalised as it is, not symmetrised).
For comparison on the SAME vectors: the phase-2 coupling (FW language, aggfull_learn.pkl). Language independence
(similarity) predicts the two agree.
Usage: timeout 590 python lab/aggfull_p6_atom.py JOB [JOB...]  with JOB in eT8_1s eT8_2s eT16_1s eT8_2p_c
-> lab/out/aggfull_p6_atom.json
"""
import os, sys, json, pickle, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggfull_core as ac
import aggfull_atom as aa
from zz_core import dirac_formula

DEV = aa.DEV
PFX = "dev_" if DEV else ""
OUT = ac.OUT
TH = ac.TH_A
JS = os.path.join(OUT, PFX + "aggfull_p6_atom.json")
JOBS = {"eT8_1s": (512, 256, "1s", 1, -1), "eT8_2s": (512, 256, "2s", 2, -1), "eT16_1s": (256, 256, "1s", 1, -1),
        "eT8_2p_c": (512, 256, "2p", 2, None)}
if DEV:
    JOBS = {"e8_2p": (512, 96, "2p", 2, None), "e8_1s": (512, 96, "1s", 1, -1)}


def first_order_full(C, offs, S, gamma0):
    dC = C.copy()
    i0 = [i for i, o in enumerate(offs) if not np.any(o)][0]
    dC[:, :, i0, i0] -= gamma0 * np.eye(2)
    P = S.shape[0]
    return np.einsum("abjl,qjpl->qapb", dC, S).reshape(2 * P, 2 * P)


def g0_of(C, offs):
    return complex(np.trace(ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]) / 2).real


def dirac_shift_eV(n, kap):
    return -(1 - dirac_formula(n, kap)) * ac.MC2 + ac.RY / n ** 2


def run(job):
    res = json.load(open(JS)) if os.path.exists(JS) else {}
    if job in res:
        print(job, "done already"); return
    t0 = time.time()
    B, N, kind, n, kap = JOBS[job]
    h = B * TH
    rec = [l for l in pickle.load(open(os.path.join(OUT, "aggfull_relearnT.pkl"), "rb")) if l["B"] == B][0]
    T, _, _ = aa.kinetic(rec, N, h)
    mk, _ = aa.kinetic_mass(rec, h)
    lc = [s for s in pickle.load(open(os.path.join(OUT, PFX + "aggfull_p6_learn.pkl"), "rb")) if s["B"] == B][0]
    fw = [s for s in pickle.load(open(os.path.join(OUT, PFX + "aggfull_learn.pkl"), "rb")) if s["B"] == B][0]
    assert np.array_equal(np.asarray(lc["offs"]), np.asarray(fw["offs"]))
    vf = os.path.join(OUT, f"{PFX}aggfull_vec_{job}.npy")
    v = np.load(vf).reshape((N,) * 3)
    v = v / np.linalg.norm(v)
    phis = [np.transpose(v, (2, 1, 0)), np.transpose(v, (0, 2, 1)), v] if kind == "2p" else [v]
    f = ac.GalerkinField(N, h) if DEV else None
    if f is None:
        import importlib
        sys.argv = [sys.argv[0], "eT8_1s"]                      # aggfull_run reads argv at import (T mode file names)
        ar = importlib.import_module("aggfull_run")
        f = ar.field(N, h)
    S = aa.coupling_blocks(f, phis, lc["offs"])
    Rk = aa.kinetic_remainder(T, h, v, mk)
    out = dict(B=B, N=N, h=h, state=kind, m_kin=mk, kinetic_remainder_eV=Rk * ac.MC2)
    for name, s in (("LC", lc), ("FW", fw)):
        g0 = g0_of(s["C"], s["offs"])
        M = first_order_full(s["C"], s["offs"], S, g0)
        ev = np.linalg.eigvals(M)
        ev = ev[np.argsort(ev.real)]
        out[name] = dict(g0=g0, dV_eV=[[e.real * ac.MC2, e.imag * ac.MC2] for e in ev],
                         antiherm_rel=float(np.linalg.norm(M - M.conj().T) / max(np.linalg.norm(M), 1e-300)))
        if kind == "2p":
            d12, d32 = dirac_shift_eV(2, 1), dirac_shift_eV(2, -2)
            e12 = Rk * ac.MC2 + ev[:2].real.mean() * ac.MC2
            e32 = Rk * ac.MC2 + ev[2:].real.mean() * ac.MC2
            out[name].update(shift12_eV=e12, shift32_eV=e32, err12=e12 / d12 - 1, err32=e32 / d32 - 1)
        else:
            d = dirac_shift_eV(n, kap)
            e = Rk * ac.MC2 + ev.real.mean() * ac.MC2
            out[name].update(shift_eV=e, dirac_eV=d, err=e / d - 1, darwin_part_eV=ev.real.mean() * ac.MC2,
                             darwin_dirac_eV=0.5 * ac.zz.ALPHA ** 4 * ac.MC2 / n ** 3,
                             kinetic_dirac_eV=-(ac.zz.ALPHA ** 4 * ac.MC2 / (2 * n ** 4)) * (4 * n / (0 + 0.5) - 3) / 4)
        print(job, name, {k: v for k, v in out[name].items() if k not in ("dV_eV",)}, flush=True)
    out["seconds"] = time.time() - t0
    res = json.load(open(JS)) if os.path.exists(JS) else {}
    res[job] = out
    json.dump(res, open(JS, "w"), indent=1, default=float)
    print(job, f"R {Rk * ac.MC2:.4e} eV ({out['seconds']:.0f} s)", flush=True)


if __name__ == "__main__":
    for j in sys.argv[1:]:
        run(j)
