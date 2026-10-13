"""aggfull_learn.py - P1/P2: the agglomerator learns HOW THE COARSE CELL FEELS THE FIELD (coupling Gamma), level by level.

Chain A (same levels, cells and blocks as lab/agg_zz_3d.py: B = 16 -> 256 -> 512 fine cells, theta = a0/4096):
  level 1  data = the fine 3-D zig-zag walk, 16 ticks, first-order response to a weak potential plane wave (aggfull_core)
  level 2  data = the LEARNED level 1 (its generator and its learned coupling), 16 steps, first order (Duhamel)
  level 3  data = the LEARNED level 2, 2 steps
For each level: wave-vector pairs (k, k') drawn like the observable state (momentum density of a 1s-like state of size
a0 up to the a0/16 cell, 2 a0 for the a0/8 cell; 75 %) plus uniform pairs in the Brillouin cube (25 %, weight 0.03),
exactly the recipe of the free learner. Candidates: stencils 'point' (1 site: a purely local field) and 'star' (7 sites);
CHOSEN AUTOMATICALLY: the smallest stencil that passes the check (pre-registered budgets, aggfull_prereg.py):
  (c1) absolute error |Gamma_fit - Gamma_data| on held-out state pairs, mean, <= EPS_ABS (relative to Gamma(0,0) ~ 1)
  (c2) relative error of the k-dependent part (Gamma - Gamma(0,0)) on held-out state pairs, mean, <= EPS_REL
  (c3) tameness: |Gamma_fit| <= 1.5 on uniform pairs of the cube.
Then also: every level against the fine walk DIRECTLY (same B, held-out pairs), and the physics read off the learned
stencil (aggfull_core.physics_of_coupling): g0, Darwin (gk, gx) and spin-orbit (gso) against Dirac (Foldy-Wouthuysen).
Writes lab/out/aggfull_learn.json and lab/out/aggfull_learn.pkl (stencils).  Run: timeout 590 python lab/aggfull_learn.py
"""
import json, os, pickle, sys, time
import numpy as np
import aggfull_core as ac
import agg_zz_core as azc

DEV = os.environ.get("AGGFULL_DEV") == "1"
PRE = json.load(open(os.path.join(ac.OUT, "aggfull_prereg.json"))) if not DEV else None
EPS_ABS = PRE["learn"]["eps_abs"] if PRE else 1e-4
EPS_REL = PRE["learn"]["eps_rel"] if PRE else 2e-2
NPAIR = 3000
TAG = "dev_aggfull_learn" if DEV else "aggfull_learn"
t_start = time.time()
levels = ac.load_chain_A()
TH = ac.TH_A
rng = np.random.default_rng(20261013)


def sample_pairs(n, h, k_s):
    """pairs (k, k') in PHYSICAL units: 75 % like the state (size 1/k_s), 25 % uniform in the cube of cell h"""
    ng = int(0.75 * n)
    kc = azc.sample_state_k(rng, ng, 3, k_s * h)
    kpc = azc.sample_state_k(rng, ng, 3, k_s * h)
    ok = np.all(np.abs(kc) < np.pi, 1) & np.all(np.abs(kpc) < np.pi, 1)
    kc, kpc = kc[ok], kpc[ok]
    nu = n - len(kc)
    ku = np.pi * (2 * rng.random((nu, 3)) - 1)
    kpu = np.pi * (2 * rng.random((nu, 3)) - 1)
    w = np.concatenate([np.ones(len(kc)), 0.03 * np.ones(nu)])
    return np.concatenate([kc, ku]) / h, np.concatenate([kpc, kpu]) / h, w, len(kc)


class Lev:
    """a learned level: free generator (chain A stencil) + learned coupling stencil"""

    def __init__(self, rec, C=None, offs=None):
        self.rec, self.B, self.h = rec, rec["B"], rec["B"] * TH
        self.rule = ac.genrule(rec)
        self.phi0 = rec["phi0"]
        self.C, self.offs = C, offs

    def gen(self, kphys):
        G = self.rule.gen_at(kphys * self.h)
        return 0.5 * (G + np.conj(np.swapaxes(G, -1, -2)))

    def coup(self, k, kp):
        return ac.eval_coupling(self.C, self.offs, k * self.h, kp * self.h)

    def response(self, k, kp, nsteps):
        """first-order change of nsteps steps (unit potential plane wave), and the nsteps-step free maps"""
        Gk, Gp = self.gen(k), self.gen(kp)
        tau = self.h
        X = tau * self.coup(k, kp)
        # one-step first-order change from Duhamel: dW = U' [ -i X~ e^{-i(l'+l)/2} sinc ] U^dag e^{i phi0}
        ek, Uk = np.linalg.eigh(Gk)
        ep, Up = np.linalg.eigh(Gp)
        Xt = np.conj(np.swapaxes(Up, -1, -2)) @ X @ Uk
        diff = ep[:, :, None] - ek[:, None, :]
        fac = -1j * np.exp(-0.5j * (ep[:, :, None] + ek[:, None, :])) * np.sinc(diff / (2 * np.pi))
        dW = Up @ (Xt * fac) @ np.conj(np.swapaxes(Uk, -1, -2)) * np.exp(1j * self.phi0)
        Wk = Uk @ (np.exp(-1j * ek)[..., None] * np.conj(np.swapaxes(Uk, -1, -2))) * np.exp(1j * self.phi0)
        Wp = Up @ (np.exp(-1j * ep)[..., None] * np.conj(np.swapaxes(Up, -1, -2))) * np.exp(1j * self.phi0)
        D = np.zeros_like(dW)
        for t in range(nsteps):
            D += np.linalg.matrix_power(Wp, nsteps - 1 - t) @ dW @ np.linalg.matrix_power(Wk, t)
        return np.linalg.matrix_power(Wk, nsteps), np.linalg.matrix_power(Wp, nsteps), D


def data_from_below(below, k, kp, B_new):
    """coupling data Gamma(kp,k) of the new level (cell B_new) from the level below (fine walk if below is None)"""
    if below is None:
        Gk, Gp, phi0, Gam, D = ac.coarse_from_fine(k, kp, TH, B_new)
        return Gam
    b = B_new // below.B
    Tk, Tp, D = below.response(k, kp, b)
    phi0 = float(np.angle(np.mean(np.linalg.eigvals(below.response(np.zeros((1, 3)), np.zeros((1, 3)), b)[0][0]))))
    Gk = ac.gen_from_step(Tk, phi0)
    Gp = ac.gen_from_step(Tp, phi0)
    X = ac.duhamel_invert(D, Gk, Gp, phi0)
    return X / (B_new * TH)


def check(C, offs, h, k, kp, Gam, nstate, ku, kpu):
    Gf = ac.eval_coupling(C, offs, k[:nstate] * h, kp[:nstate] * h)
    G00 = ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]
    D00 = Gam[:nstate]
    err_abs = float(np.mean(np.linalg.norm(Gf - D00, axis=(1, 2))) / np.linalg.norm(G00) * np.sqrt(2))
    dep = D00 - G00[None]
    err_rel = float(np.mean(np.linalg.norm(Gf - D00, axis=(1, 2))) / np.mean(np.linalg.norm(dep, axis=(1, 2))))
    Gu = ac.eval_coupling(C, offs, ku * h, kpu * h)
    tame = float(np.max(np.linalg.norm(Gu, axis=(1, 2)) / np.sqrt(2)))
    return dict(err_abs=err_abs, err_rel=err_rel, tame_max=tame,
                ok=bool(err_abs <= EPS_ABS and err_rel <= EPS_REL and tame <= 1.5))


out = dict(theta=TH, eps_abs=EPS_ABS, eps_rel=EPS_REL, levels=[])
stencils = []
below = None
for rec in levels:
    B, h = rec["B"], rec["B"] * TH
    k_s = 1 / ac.A0 if B <= 256 else 1 / (2 * ac.A0)
    t0 = time.time()
    k, kp, w, ns = sample_pairs(NPAIR, h, k_s)
    Gam = data_from_below(below, k, kp, B)
    kh, kph, wh, nsh = sample_pairs(1500, h, k_s)              # held out
    Gh = data_from_below(below, kh, kph, B)
    ku = np.pi * (2 * rng.random((2000, 3)) - 1) / h
    kpu = np.pi * (2 * rng.random((2000, 3)) - 1) / h
    # relative least squares (as the free learner): every pair's error counted against its own k-dependent part
    G0d = data_from_below(below, np.zeros((1, 3)), np.zeros((1, 3)), B)[0]
    dep = np.linalg.norm(Gam - G0d[None], axis=(1, 2))
    wr = w / (dep + 1e-3 * np.median(dep[:ns]))
    tried, chosen = [], None
    for name in ("point", "star", "cube"):
        offs = ac.stencil(name)
        C = ac.fit_coupling(k * h, kp * h, Gam, offs, wts=wr * np.sqrt(np.median(dep[:ns])))
        ck = check(C, offs, h, kh, kph, Gh, nsh, ku, kpu)
        tried.append(dict(stencil=name, **ck))
        print(f"B={B:4d} h={h:.4f} lambda_C: candidate {name:5s}: abs {ck['err_abs']:.2e} rel(k-dep) {ck['err_rel']:.2e} "
              f"tame {ck['tame_max']:.4f} {'OK' if ck['ok'] else ''}", flush=True)
        if ck["ok"]:
            chosen = (name, C, offs)
            break
    if chosen is None:
        # the automatic rule: no admissible coupling at this level -> flag; this level's coupling is not made and the
        # next level learns from the last admissible one (here: the fine walk itself)
        print(f"  no admissible stencil at B={B}: FLAG, the next level learns from the level below this one", flush=True)
        out["levels"].append(dict(B=B, cell_lambdaC=h, chosen=None, candidates=tried, flag="no admissible coupling"))
        continue
    name, C, offs = chosen
    lev = Lev(rec, C, offs)
    # direct comparison with the fine walk at the same B (held-out state pairs)
    Gfine = ac.coarse_from_fine(kh[:nsh], kph[:nsh], TH, B)[3]
    Gl = lev.coup(kh[:nsh], kph[:nsh])
    G00 = lev.coup(np.zeros((1, 3)), np.zeros((1, 3)))[0]
    dep = Gfine - np.eye(2)[None] * Gfine[:, :1, :1].mean()
    vs_fine_abs = float(np.mean(np.linalg.norm(Gl - Gfine, axis=(1, 2))) / np.sqrt(2))
    vs_fine_rel = float(np.mean(np.linalg.norm(Gl - Gfine, axis=(1, 2))) / np.mean(np.linalg.norm(Gfine - G00[None], axis=(1, 2))))
    phys = ac.physics_of_coupling(C, offs, h)
    rec_out = dict(B=B, cell_lambdaC=h, cell_over_a0=h / ac.A0, k_s_over_m=k_s, chosen=name, n_sites=len(offs),
                   n_params_real=int(4 * len(offs) ** 2), candidates=tried, vs_fine_abs=vs_fine_abs,
                   vs_fine_rel=vs_fine_rel, physics={kk: [v.real, v.imag] if isinstance(v, complex) else v for kk, v in phys.items()},
                   seconds=time.time() - t0)
    print(f"  CHOSEN {name} ({len(offs)} sites); vs fine walk directly: abs {vs_fine_abs:.2e}, rel(k-dep) {vs_fine_rel:.2e}; "
          f"g0 {phys['g0'].real:.7f}, gk {phys['gk'].real:+.5f} (Dirac -0.125), gx {phys['gx'].real:+.5f} (+0.25), "
          f"gso {phys['gso'].real:+.5f} (+0.25), SO shape resid {phys['so_shape_resid']:.1e}  ({time.time() - t0:.0f} s)", flush=True)
    out["levels"].append(rec_out)
    stencils.append(dict(B=B, h=h, C=C, offs=offs, name=name))
    below = lev
out["seconds"] = time.time() - t_start
json.dump(out, open(os.path.join(ac.OUT, TAG + ".json"), "w"), indent=1, default=float)
pickle.dump(stencils, open(os.path.join(ac.OUT, TAG + ".pkl"), "wb"))
print(f"done ({out['seconds']:.0f} s)")
