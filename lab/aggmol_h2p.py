"""aggmol_h2p.py - task (3) of study aggmol_: COLLAUDO A DUE CENTRI. One electron, two protons held fixed (H2+ with clamped
nuclei: a single electron, declared as a test of the two-centre machine; NOT the H2 bond, which needs two electrons and the
many-electron rule). The coarse rule: learned kinetic (T) + learned coupling gamma0 x the Galerkin field of the two nuclei
(input: -alpha/r of each proton and its position, nothing else). Exam key: lab/aggmol_exact.py (exact two-centre solution).

Jobs (one call = one job, each <= ~9 min; results appended to lab/out/aggmol_h2p.json):
  atom_<f>            1s of one H atom at cell a0/f, box +-8 a0 (for D_e at the same level and box)
  curve_<f>_R<R>      1 sigma_g at separation R (a0), molecule along z, box +-8 a0, parity sector (+,+,+)
  exc_<f>_<state>     R = 2: 1su (+,+,-), 1pu (-,+,+), 2sg (second state of (+,+,+), node test), box +-12 a0
  shift_<f>_<i>       R = 2, both protons shifted by a fraction of a cell (translation covariance), no parity imposed
                      except where exact (half-cell shift: reflections about the half-cell planes)
  rot_<f>             R = 2, molecule along the body diagonal (1,1,1): only the inversion sector
Usage: timeout 590 python lab/aggmol_h2p.py JOB [JOB ...] | list"""
import json, os, sys, time
import numpy as np
import scipy.fft as sfft
from scipy.ndimage import map_coordinates
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_core as mc
import aggmol_exact as ex

PFX = "dev_" if mc.DEV else ""
JS = os.path.join(mc.OUT, PFX + "aggmol_h2p.json")
res = json.load(open(JS)) if os.path.exists(JS) else {}
FRAC_B = {16: 256, 8: 512, 4: 1024, 2: 2048}
RS = (1.0, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 3.0, 4.0)
SHIFTS = [(-0.5, -0.5, -0.5), (0.31, 0.17, 0.43), (0.5, 0.0, 0.25)]
EH = mc.ALPHA ** 2                       # Hartree in units of m c^2


def save():
    json.dump(res, open(JS, "w"), indent=1, default=float)


def setup(f, box):
    B = FRAC_B[f]
    h = B * mc.TH
    N = int(round(2 * box * mc.A0 / h))
    rec = mc.level_rec(B)
    C, offs, g0 = mc.coupling_rec(B)
    T = mc.kinetic_grid(rec, N, h)
    return B, h, N, T, g0.real


def run(mol, X0, tag, budget=520, rtol=2e-8, Y=None):
    X, E, r, conv, its = mc.solve_until(mol, X0, budget, rtol=rtol, Y=Y, tag=tag)
    return X, [e / EH for e in E], r, conv, its


def job_atom(f):
    B, h, N, T, g0 = setup(f, 8)
    fld = mc.MultiField(N, h, [(0, 0, 0)], dtype=np.float32 if N > 192 else np.float64)
    m = mc.Mol(T, fld, gamma0=g0, sym=(1, 1, 1))
    X, E, r, conv, its = run(m, mc.guess_mol(N, h, [(0, 0, 0)], "1s")[:, None], f"atom a0/{f}")
    return dict(f=f, N=N, E_Eh=E[0], residual=r[0], converged=conv, its=its, rel_vs_exact=E[0] / -0.5 - 1)


def centres(R, h, axis=(0, 0, 1), shift=(0, 0, 0)):
    u = np.array(axis, float) / np.linalg.norm(axis)
    c = R * mc.A0 / 2 / h * u
    return [tuple(c + np.array(shift)), tuple(-c + np.array(shift))]


def job_curve(f, R):
    B, h, N, T, g0 = setup(f, 8)
    cs = centres(R, h)
    fld = mc.MultiField(N, h, cs, dtype=np.float32 if N > 192 else np.float64)
    m = mc.Mol(T, fld, gamma0=g0, sym=(1, 1, 1))
    X, E, r, conv, its = run(m, mc.guess_mol(N, h, cs, "1sg")[:, None], f"1sg a0/{f} R={R}")
    Eex = ex.energy(R)[0]
    return dict(f=f, R=R, N=N, E_el_Eh=E[0], E_tot_Eh=E[0] + 1 / R, exact_el_Eh=Eex, err_mEh=1e3 * (E[0] - Eex),
                rel=E[0] / Eex - 1, residual=r[0], converged=conv, its=its)


def ray_nodes(v, N, h, R, ndir=24):
    """zero crossings of the state along rays from the midpoint (x-z half plane, angles from the axis 0..90 deg and two
    azimuths), by band-limited (spectral, zero-padded x2) + cubic interpolation; returns xi = (r_a + r_b)/R and r"""
    F = sfft.fftn(v.reshape((N,) * 3))
    M = 2 * N
    G = np.zeros((M, M, M), complex)
    k = np.fft.fftfreq(N) * N
    idx = np.rint(k).astype(int) % M
    G[np.ix_(idx, idx, idx)] = F
    up = np.real(sfft.ifftn(G)) * (M / N) ** 3
    del G, F
    out = []
    Rl = R * mc.A0
    for az in (0.0, np.pi / 4):
        for th in np.linspace(0, np.pi / 2, ndir):
            u = np.array([np.sin(th) * np.cos(az), np.sin(th) * np.sin(az), np.cos(th)])
            rr = np.linspace(0.3, 5.0, 1200) * mc.A0
            P = rr[:, None] * u[None]
            ii = (P / (h / 2) + M // 2).T          # index on the 2N grid (origin at index M//2)
            val = map_coordinates(up, ii, order=3, mode="grid-wrap")
            sg = np.where(np.sign(val[:-1]) != np.sign(val[1:]))[0]
            if len(sg) == 0:
                continue
            j = sg[0]
            r0 = rr[j] - val[j] * (rr[j + 1] - rr[j]) / (val[j + 1] - val[j])
            p0 = r0 * u
            ra = np.linalg.norm(p0 - np.array([0, 0, Rl / 2])); rb = np.linalg.norm(p0 + np.array([0, 0, Rl / 2]))
            out.append(dict(theta=float(th), az=float(az), r=float(r0 / mc.A0), xi=float((ra + rb) / Rl)))
    return out


def job_exc(f, state):
    B, h, N, T, g0 = setup(f, 12)
    R = 2.0
    cs = centres(R, h)
    fld = mc.MultiField(N, h, cs)
    sym = {"1su": (1, 1, -1), "1pu": (-1, 1, 1), "2sg": (1, 1, 1)}[state]
    m = mc.Mol(T, fld, gamma0=g0, sym=sym)
    if state == "2sg":
        X0 = np.stack([mc.guess_mol(N, h, cs, "1sg"), mc.guess_mol(N, h, cs, "2sg")], 1)
    else:
        X0 = mc.guess_mol(N, h, cs, state)[:, None]
    X, E, r, conv, its = run(m, X0, f"{state} a0/{f}", budget=500)
    q = dict(st=(0, 1, 0), pu=(0, 0, 1))
    Eex = ex.energy(R, *ex.STATES[state])[0]
    e = E[-1]
    out = dict(f=f, R=R, N=N, state=state, E_el_Eh=e, exact_el_Eh=Eex, rel=e / Eex - 1, err_mEh=1e3 * (e - Eex),
               residual=r[-1], converged=conv, its=its, all_E=E)
    if state == "2sg":
        out["exact_1sg_el_Eh"] = ex.energy(R)[0]
        nodes = ray_nodes(X[:, 1], N, h, R)
        xi = np.array([n["xi"] for n in nodes]); rr = np.array([n["r"] for n in nodes])
        xi0_ex = ex.radial_node(R, 1)[0][0]
        out.update(nodes=nodes, xi_mean=float(xi.mean()), xi_spread_rel=float(xi.std() / xi.mean()),
                   r_spread_rel=float(rr.std() / rr.mean()), xi0_exact=xi0_ex,
                   xi0_err_rel=float((xi.mean() - xi0_ex) / (xi0_ex - 1)))
        np.save(os.path.join(mc.OUT, PFX + f"aggmol_vec_2sg_a0_{f}.npy"), X[:, 1].astype(np.float32))
        np.save(os.path.join(mc.OUT, PFX + f"aggmol_vec_1sg_a0_{f}.npy"), X[:, 0].astype(np.float32))
    return out


def job_shift(f, i):
    B, h, N, T, g0 = setup(f, 8)
    R = 2.0
    sh = SHIFTS[i]
    cs = centres(R, h, shift=sh)
    fld = mc.MultiField(N, h, cs)
    sym = (2, 2, 2) if all(abs(s + 0.5) < 1e-12 for s in sh) else (0, 0, 0)
    m = mc.Mol(T, fld, gamma0=g0, sym=sym)
    X, E, r, conv, its = run(m, mc.guess_mol(N, h, cs, "1sg")[:, None], f"shift{i} a0/{f}")
    return dict(f=f, R=R, shift_cells=sh, sym=sym, E_el_Eh=E[0], residual=r[0], converged=conv, its=its)


def job_rot(f):
    B, h, N, T, g0 = setup(f, 8)
    R = 2.0
    cs = centres(R, h, axis=(1, 1, 1))
    fld = mc.MultiField(N, h, cs)
    m = mc.Mol(T, fld, gamma0=g0, sym=(0, 0, 0), inv=1)
    X, E, r, conv, its = run(m, mc.guess_mol(N, h, cs, "1sg")[:, None], f"rot a0/{f}")
    return dict(f=f, R=R, axis="111", E_el_Eh=E[0], residual=r[0], converged=conv, its=its)


def dispatch(job):
    p = job.split("_")
    if p[0] == "atom":
        return job_atom(int(p[1]))
    if p[0] == "curve":
        return job_curve(int(p[1]), float(p[2][1:]))
    if p[0] == "exc":
        return job_exc(int(p[1]), p[2])
    if p[0] == "shift":
        return job_shift(int(p[1]), int(p[2]))
    if p[0] == "rot":
        return job_rot(int(p[1]))
    raise ValueError(job)


if __name__ == "__main__":
    for job in sys.argv[1:]:
        if job in res and res[job].get("converged"):
            print(f"{job}: done", flush=True)
            continue
        t0 = time.time()
        out = dispatch(job)
        out["seconds"] = time.time() - t0
        res[job] = out
        save()
        print(f"{job}: " + ", ".join(f"{k}={v:.8g}" if isinstance(v, float) else f"{k}={v}" for k, v in out.items()
                                     if k not in ("nodes", "all_E")), flush=True)
