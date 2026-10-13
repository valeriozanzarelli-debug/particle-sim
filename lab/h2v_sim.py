"""h2v_sim.py - VIDEO H2: the balls. Every electron ball is a sample of the computed two-electron state (Metropolis on
|Psi|^2; its motion is the sampler, NOT a trajectory: declared). Protons: two glows at +-R/2.

    python3 lab/h2v_sim.py exams        # V1-V4b, V7 (analytic side + sampler side) -> lab/out/h2v_esami.json
    python3 lab/h2v_sim.py frames       # writes <SCR>/h2v_frames/f_XXXX.npz, resumable (run again until FRAMES DONE)

Timeline (24 fps): A 0-240 still molecule (R = 1.40), B 240-528 correlation, C 528-912 vibration.
"""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h2v_core as C

OUT = C.OUT
SCR = os.environ.get("H2V_SCR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "h2v_scr"))
FR = os.path.join(SCR, "h2v_frames")
os.makedirs(FR, exist_ok=True)
FPS = 24
NA, NB, NC = 240, 288, 384
NF = NA + NB + NC
NW = 12000                 # pair walkers -> 24000 electron balls
NO = 8000                  # orange walkers (other electron given the reference one)
STEP = 0.07                # a0 per coordinate per Metropolis move during the video (shimmer)
SLOW_PERIOD_S = 4.0        # one vibration period shown in 4 s
JE = os.path.join(OUT, "h2v_esami.json")
esami = json.load(open(JE)) if os.path.exists(JE) else {}


def save_esami():
    json.dump(esami, open(JE, "w"), indent=1, default=float)


class Target:
    """pi(X) = sum_R w_R |Psi_R(X)|^2"""
    def __init__(self):
        self.wfs = [None] * len(C.RGRID)

    def wf(self, k):
        if self.wfs[k] is None:
            self.wfs[k] = C.WF(C.RGRID[k])          # loaded when first needed (bases still being rebuilt)
        return self.wfs[k]

    def __call__(self, X, w):
        p = np.zeros(len(X))
        for k in range(len(self.wfs)):
            if w[k] > 1e-4:
                p += w[k] * self.wf(k).psi(X) ** 2
        return p


def metropolis(X, logp_fn, step, rng, nstep, p=None):
    if p is None:
        p = logp_fn(X)
    acc = 0
    for _ in range(nstep):
        Y = X + rng.normal(0, step, X.shape)
        q = logp_fn(Y)
        ok = rng.random(len(X)) * p < q
        X[ok] = Y[ok]; p[ok] = q[ok]
        acc += ok.mean()
    return X, p, acc / max(nstep, 1)


def init_walkers(n, R, rng):
    X = rng.normal(0, 0.8, (n, 2, 3))
    X[:, :, 2] += rng.choice([-R / 2, R / 2], (n, 2))
    return X


def ref_path(u):
    """the reference electron in part B: a CHOSEN slow loop (declared), u in [0, 1]; returns (x, y, z) in a0"""
    a = 2 * np.pi * u
    return np.array([0.55 * np.sin(2 * a), 0.0, -1.25 * np.cos(a)])


# ---------------------------------------------------------------- exams
def job_exams():
    rng = np.random.default_rng(7)
    wf = C.WF(1.4)
    t0 = time.time()
    if "V1" not in esami:
        ex = wf.expect()
        V = C.vib()
        dEdR = float(V["spline"](1.4, 1))
        Etot = ex["E_el"] + 1 / 1.4
        Tpred = -Etot - 1.4 * dEdR
        esami["expect_R1.40"] = ex
        esami["V1"] = dict(T=ex["T"], T_virial=Tpred, diff=ex["T"] - Tpred, ratio_T_over_minusE=ex["T"] / -Etot,
                           norm_minus_1=ex["norm_vs_solver"] - 1, E_el_diff=ex["E_el"] - ex["E_el_solver"],
                           PASS=bool(abs(ex["T"] - Tpred) <= 1e-3 and abs(ex["norm_vs_solver"] - 1) <= 1e-6))
        esami["V2"] = dict(inv_r12=ex["inv_r12"], ref=0.5873, rel=ex["inv_r12"] / 0.5873 - 1,
                           PASS=bool(abs(ex["inv_r12"] / 0.5873 - 1) <= 0.01),
                           note="riferimento citato a memoria, da verificare")
        save_esami(); print(esami["V1"], esami["V2"], flush=True)
    # sampler: 1.2e6 balls from 50 well separated snapshots
    edges = np.linspace(-4, 4, 41)
    tgt = lambda X: wf.psi(X) ** 2
    X = init_walkers(NW, 1.4, rng)
    X, p, acc0 = metropolis(X, tgt, 0.35, rng, 150)
    H = np.zeros((40, 40)); n = 0
    opp = 0; nn = 0; cl = 0; cr = 0; az = 0.0
    for snap in range(50):
        X, p, acc = metropolis(X, tgt, 0.35, rng, 12, p)
        for e in range(2):
            h, _, _ = np.histogram2d(X[:, e, 0], X[:, e, 2], [edges, edges])
            H += h
        n += 2 * len(X)
        z1, z2 = X[:, 0, 2], X[:, 1, 2]
        opp += np.sum(z1 * z2 < 0); nn += len(X)
        for a, b in ((z1, z2), (z2, z1)):
            m = a < -0.35
            cl += m.sum(); cr += np.sum(m & (b > 0))
        az += np.sum(np.abs(X[:, :, 2]))
    ana = wf.bins_xz(edges)
    Hn = H / n * 2               # electrons per cell (of 2)
    inside = ana.sum()
    pear = float(np.corrcoef(Hn.ravel(), ana.ravel())[0, 1])
    L1 = float(np.abs(Hn - ana).sum() / ana.sum())
    esami["V3"] = dict(n_balls=n, pearson=pear, L1=L1, frac_inside_box_analytic=float(inside / 2),
                       acc=float(acc), PASS=bool(pear >= 0.995 and L1 <= 0.05))
    Popp_a = wf.prob_z(C.images(-np.inf, 0, 0, np.inf)) * 2      # (z1<0, z2>0) or (z1>0, z2<0)
    Popp_s = opp / nn
    esami["V4"] = dict(P_opp_analytic=Popp_a, P_opp_balls=float(Popp_s), HF_value=0.5,
                       PASS=bool(abs(Popp_a - Popp_s) <= 0.01 and Popp_a >= 0.52))
    # V4b: P(z2 > 0 | z1 < -0.35)
    num = wf.prob_z(C.images(-np.inf, -0.35, 0, np.inf))
    den = wf.prob_z(C.images(-np.inf, -0.35, -np.inf, np.inf))
    esami["V4b"] = dict(P_analytic=num / den, P_balls=float(cr / cl), PASS=bool(abs(num / den - cr / cl) <= 0.02))
    esami["abs_z_R1.40"] = dict(analytic=wf.abs_z(), balls=float(az / (2 * nn)))
    esami["sec_exams"] = time.time() - t0
    save_esami()
    np.save(os.path.join(SCR, "h2v_V3_hist.npy"), np.stack([Hn, ana]))
    for k in ("V3", "V4", "V4b", "abs_z_R1.40"):
        print(k, esami[k], flush=True)


# ---------------------------------------------------------------- frames
def job_frames():
    t0 = time.time()
    st_file = os.path.join(SCR, "h2v_state.npz")
    i14 = C.RGRID.index(1.4)
    V = C.vib()
    period_au = 2 * np.pi / (V["E"][1] - V["E"][0])
    rng = np.random.default_rng(11)
    tgt = Target()
    w14 = np.zeros(len(C.RGRID)); w14[i14] = 1
    if os.path.exists(st_file):
        z = np.load(st_file, allow_pickle=True)
        X, Y, f0 = z["X"], z["Y"], int(z["f"])
        rng.bit_generator.state = json.loads(str(z["rng"]))
    else:
        X = init_walkers(NW, 1.4, rng)
        X, _, _ = metropolis(X, lambda Q: tgt(Q, w14), 0.35, rng, 200)
        Y = None; f0 = 0
    xR = V["x"]
    for f in range(f0, NF):
        need_wait = f == NA + NB and not all(os.path.exists(os.path.join(OUT, f"h2v_basis_R{R:.2f}.npz")) for R in C.RGRID)
        if time.time() - t0 > 420 or f % 120 == 0 or need_wait:
            np.savez(st_file, X=X, Y=Y if Y is not None else np.zeros(0), f=f, rng=json.dumps(rng.bit_generator.state))
            if time.time() - t0 > 420 or need_wait:
                print(f"checkpoint at frame {f}: run again" + (" (waiting for the bases)" if need_wait else ""), flush=True); return
        meta = dict(f=f)
        if f < NA + NB:
            w = w14; R = 1.4
            X, _, acc = metropolis(X, lambda Q: tgt(Q, w), STEP, rng, 2)
            balls = [X.reshape(-1, 3)]; cls = [np.zeros(2 * NW, np.uint8)]
            meta.update(seg="A" if f < NA else "B", R=R, Rmean=R)
            if f >= NA:
                u = (f - NA) / NB
                r1 = ref_path(u)
                cond = lambda Q: tgt.wf(i14).psi(np.concatenate([np.broadcast_to(r1, (len(Q), 1, 3)), Q[:, None, :]], 1)) ** 2
                if Y is None or len(Y) == 0:
                    Y = rng.normal(0, 0.8, (NO, 3)); Y[:, 2] += rng.choice([-0.7, 0.7], NO)
                    Y, _, _ = metropolis(Y, cond, 0.35, rng, 150)
                Y, _, acc2 = metropolis(Y, cond, 0.09, rng, 3)
                balls += [Y.copy(), r1[None]]; cls += [np.ones(NO, np.uint8), np.array([2], np.uint8)]
                meta.update(ref=r1.tolist(), acc_cond=float(acc2))
        else:
            tC = (f - NA - NB) / FPS
            t_au = tC / SLOW_PERIOD_S * period_au
            dens = C.packet(V, t_au)
            w = C.grid_weights(V, dens)
            Rm = float(np.sum(xR * dens) * V["dx"])
            X, _, acc = metropolis(X, lambda Q: tgt(Q, w), STEP, rng, 2)
            balls = [X.reshape(-1, 3)]; cls = [np.zeros(2 * NW, np.uint8)]
            meta.update(seg="C", R=Rm, Rmean=Rm, t_fs=t_au * C.AU_FS, w=w.tolist())
        meta["acc"] = float(acc)
        prot = np.array([[0, 0, meta["Rmean"] / 2], [0, 0, -meta["Rmean"] / 2]])
        np.savez(os.path.join(FR, f"f_{f:04d}.npz"), p=np.concatenate(balls).astype(np.float32),
                 c=np.concatenate(cls), prot=prot, meta=json.dumps(meta))
        if f % 24 == 0:
            print(f"frame {f} {meta.get('seg')} R {meta['Rmean']:.3f} acc {acc:.2f} ({time.time() - t0:.0f} s)", flush=True)
    print("FRAMES DONE", flush=True)


# ---------------------------------------------------------------- part C, version 2 (after the V7 FAIL, declared)
NWC = 8000      # pair walkers in C v2
STEPC = 0.32    # a0, Metropolis step in C v2
NSTC = 6        # moves per frame in C v2 (v1: 2 moves of 0.07 a0: the cloud lagged behind the protons, V7 FAIL)
WMIN = 1e-2     # grid weights below this are dropped and the rest renormalised (speed; declared)


def job_frames_c2():
    t0 = time.time()
    st_file = os.path.join(SCR, "h2v_state_c2.npz")
    V = C.vib()
    period_au = 2 * np.pi / (V["E"][1] - V["E"][0])
    rng = np.random.default_rng(21)
    tgt = Target()
    xR = V["x"]

    def weights(f):
        t_au = (f - NA - NB) / FPS / SLOW_PERIOD_S * period_au
        dens = C.packet(V, t_au)
        w = C.grid_weights(V, dens); w[w < WMIN] = 0; w /= w.sum()
        return w, t_au, float(np.sum(xR * dens) * V["dx"])
    if os.path.exists(st_file):
        z = np.load(st_file, allow_pickle=True)
        X, f0 = z["X"], int(z["f"])
        rng.bit_generator.state = json.loads(str(z["rng"]))
    else:
        w, _, _ = weights(NA + NB)
        X = init_walkers(NWC, 1.5, rng)
        X, _, _ = metropolis(X, lambda Q: tgt(Q, w), 0.35, rng, 200)
        f0 = NA + NB
    for f in range(f0, NF):
        if time.time() - t0 > 400:
            np.savez(st_file, X=X, f=f, rng=json.dumps(rng.bit_generator.state))
            print(f"checkpoint at frame {f}: run again", flush=True); return
        w, t_au, Rm = weights(f)
        X, _, acc = metropolis(X, lambda Q: tgt(Q, w), STEPC, rng, NSTC)
        meta = dict(f=f, seg="C", R=Rm, Rmean=Rm, t_fs=t_au * C.AU_FS, w=w.tolist(), acc=float(acc), version="C2")
        prot = np.array([[0, 0, Rm / 2], [0, 0, -Rm / 2]])
        np.savez(os.path.join(FR, f"f_{f:04d}.npz"), p=X.reshape(-1, 3).astype(np.float32),
                 c=np.zeros(2 * NWC, np.uint8), prot=prot, meta=json.dumps(meta))
        if f % 12 == 0:
            print(f"frame {f} C2 R {Rm:.3f} acc {acc:.2f} ({time.time() - t0:.0f} s)", flush=True)
    print("C2 DONE", flush=True)


if __name__ == "__main__":
    dict(exams=job_exams, frames=job_frames, framesC2=job_frames_c2)[sys.argv[1]]()
