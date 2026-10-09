"""fluid_f1.py - F1 production. Usage:
  python fluid_f1.py cs                         -> sound speed (writes out/fluid_f1_cs.json)
  python fluid_f1.py run R_w L reps seed M0f dV tag  -> kicks (writes out/fluid_f1_<tag>.npz/.json)
Protocol in out/fluid_prereg.json (F1)."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fluid_core import FluidMD, wca_f_over_r, wca_u, RC

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
RHO0, T0 = 0.75, 1.0


def thermo(rho, T, seed):
    md = FluidMD(32, rho, -10.0, M0=1.0, T0=T, seed=seed)  # wall never active
    for k in range(1500):
        md.step()
        if k % 10 == 0: md.rescale(T)
    Ps, Es = [], []
    for k in range(4000):
        md.step()
        if k % 10 == 0: md.rescale(T)
        if k % 20 == 0:
            d = md.mic(md.pos[md.i] - md.pos[md.j]); r2 = (d * d).sum(1); m = r2 < RC * RC
            Ps.append((md.N * T + 0.5 * (wca_f_over_r(r2[m]) * r2[m]).sum()) / md.L ** 2)
            Es.append(wca_u(r2[m]).sum() / md.N + T)
    return float(np.mean(Ps)), float(np.mean(Es)), md.N / md.L ** 2


def sound():
    r1 = thermo(RHO0 * 0.97, T0, 11); r2 = thermo(RHO0 * 1.03, T0, 12)
    r3 = thermo(RHO0, T0 * 0.95, 13); r4 = thermo(RHO0, T0 * 1.05, 14)
    cT2 = (r2[0] - r1[0]) / (r2[2] - r1[2])
    dPdT = (r4[0] - r3[0]) / (0.1 * T0); cv = (r4[1] - r3[1]) / (0.1 * T0); rr = 0.5 * (r3[2] + r4[2])
    cS2 = cT2 + T0 * dPdT ** 2 / (rr ** 2 * cv)
    res = dict(runs=[r1, r2, r3, r4], c_T=float(np.sqrt(cT2)), c_S=float(np.sqrt(cS2)), c_v=float(cv), rho=float(rr))
    json.dump(res, open(os.path.join(OUT, "fluid_f1_cs.json"), "w"), indent=1)
    print(res)


def run(Rw, L, reps, seed, M0f, dV, tag):
    cS = json.load(open(os.path.join(OUT, "fluid_f1_cs.json")))["c_S"]
    M0 = M0f * RHO0 * np.pi * Rw ** 2
    md = FluidMD(L, RHO0, Rw, M0=M0, T0=T0, seed=seed)
    t0 = time.time()
    for k in range(2000):
        md.step()
        if k % 10 == 0: md.rescale(T0)
    nst = int(round(L / cS / md.dt))
    edges = np.arange(0, L / 2 + 1e-9, 0.25); H = np.zeros(len(edges) - 1); nH = 0
    P = M0 * dV; R = []
    for rep in range(reps):
        for k in range(400):
            md.step()
            if k % 10 == 0: md.rescale(T0)
            if k % 50 == 0: H += md.radial(edges); nH += 1
        st = md.state(); res = []
        for sg in (1, -1):
            md.set_state(st); md.V = md.V + sg * np.array([dV, 0.0]); vs = np.empty(nst)
            for k in range(nst):
                md.step(); vs[k] = md.V[0]
            res.append(vs)
        md.set_state(st)
        R.append(0.5 * (res[0] - res[1]))
        if rep == 0 and tag == "R10":
            np.savez_compressed(os.path.join(OUT, "fluid_f1_snapshot.npz"), pos=st[0], X=st[2], L=L, Rw=Rw)
    R = np.array(R)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2); dens = H / nH / area
    lo = edges[:-1]
    rb = float(dens[(lo >= Rw + 8) & (edges[1:] <= L / 2 - 1)].mean())
    rcut = Rw + 6
    cnt = (H / nH)[edges[1:] <= rcut + 1e-9].sum()
    Nmiss = float(rb * np.pi * rcut ** 2 - cnt)
    Reff = float(np.sqrt(Nmiss / (np.pi * rb)))
    t = np.arange(1, nst + 1) * md.dt
    def Mfl(s, Rr):
        m = Rr.mean(0)
        f = m * np.exp(-s * t)
        Vh = md.dt * (0.5 * dV + f[:-1].sum() + 0.5 * f[-1])
        return float(P / (s * Vh) - M0)
    sgrid = np.geomspace(0.05, 10, 60)
    Ms = np.array([Mfl(s, R) for s in sgrid])
    # half-split error
    Mh = np.array([[Mfl(s, R[: reps // 2]) for s in sgrid], [Mfl(s, R[reps // 2:]) for s in sgrid]])
    err = np.abs(Mh[0] - Mh[1]) / 2
    sstar = cS / Reff
    out = dict(tag=tag, Rw=Rw, L=L, reps=reps, seed=seed, M0=M0, dV=dV, N=md.N, rho_b=rb, N_miss=Nmiss, R_eff=Reff,
               c_S=cS, s_star=sstar, M_star=Mfl(sstar, R), M_lo=Mfl(cS / (2 * Reff), R), M_hi=Mfl(2 * cS / Reff, R),
               M_star_err=float(abs(Mfl(sstar, R[: reps // 2]) - Mfl(sstar, R[reps // 2:])) / 2),
               ratio_star=Mfl(sstar, R) / Nmiss, tail_V_end=float(R.mean(0)[-1]), seconds=time.time() - t0)
    np.savez_compressed(os.path.join(OUT, "fluid_f1_%s.npz" % tag), t=t, Vmean=R.mean(0), Vse=R.std(0) / np.sqrt(reps),
                        sgrid=sgrid, Ms=Ms, Merr=err, edges=edges, dens=dens)
    json.dump(out, open(os.path.join(OUT, "fluid_f1_%s.json" % tag), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "cs":
        sound()
    else:
        a = sys.argv[2:]
        run(float(a[0]), float(a[1]), int(a[2]), int(a[3]), float(a[4]), float(a[5]), a[6])
