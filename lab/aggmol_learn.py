"""aggmol_learn.py - task (2) of study aggmol_: the agglomerator learns the next COARSER levels needed for the two-atom
scale (a0/4 and a0/2), with the same machine and the same exact constraints as phase 3-4 of AGGLOMERATO_COMPLETO.

Starting point (frozen, learned by aggfull_): level a0/8 = B 512 of the Taylor-pinned chain (lab/out/aggfull_relearnT.pkl)
and its learned coupling (lab/out/aggfull_learn.pkl).

For every candidate block b (largest first: b = 4 -> a0/2, then b = 2 -> a0/4; then REPEAT from the chosen level):
  LEARN (free)    fit mode "T" of aggfull_l1 (frozen copy aggmol_fz_l1): spin-blind band with c2, c4, c6 pinned EXACTLY to
                  the level below along 13 directions (read at kappa <= 1.5 per new cell), parity G(-k) = G(k), traceless
                  spin part fitted; data = plane waves through b steps of the level below (75 % like a 1s state of size
                  a0, 25 % uniform with weight 3e-4); width w = the smallest of (4, 5) passing the check.
  POSITIVITY      the moment bound: a positive (Markov) rule needs |c4|/c2 >= 1/12 (cell units) while the level below
                  asks for ~1/(4 h^2) << 1/12 -> Markov must be infeasible: fit mode "M" (w = 4) and report the constraint
                  deviation (predicted >> 1).
  CHECK (free)    agg_zz_core.kcheck in 'state' mode (observable: 1s-like state of size a0): packet <= 1e-3, energy and
                  rest <= 3e-4, bz_rel90 <= 0.2, no negative motion energy (<= 0.01).
  LEARN (field)   the cell-field coupling Gamma(k', k), sandwich stencil point/star/cube chosen automatically, from the
                  first-order response of the LEARNED level below (Duhamel), as aggfull_learn.py; check abs <= 1e-4,
                  rel(k-dep) <= 2e-2, tameness <= 1.5; Darwin and spin-orbit read off afterwards.
  CHECK (atom)    the 1s with the new level (kinetic + gamma0 Galerkin field, box +-8 a0) against the 1s of the level below
                  (same box): |dE| <= 1.594 mEh (chemical accuracy, 1 kcal/mol). Pass -> the level is ADMISSIBLE for atoms
                  and molecules; fail -> flagged 'preview only' (kept for the Intel tier, never used for exams of accuracy).
Writes lab/out/aggmol_levels.pkl (not committed) and lab/out/aggmol_learn.json.
Usage: timeout 590 python lab/aggmol_learn.py [dev]"""
import json, os, sys, pickle, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import aggmol_core as mc
import aggmol_fz_core as ac
import aggmol_fz_l1 as L
import agg_zz_core as azc
from agg_zz_core import Level

DEV = len(sys.argv) > 1 and sys.argv[1] == "dev"
TAG = "dev_" if DEV else ""
JS = os.path.join(mc.OUT, TAG + "aggmol_learn.json")
PK = os.path.join(mc.OUT, TAG + "aggmol_levels.pkl")
TH = mc.TH
A0 = mc.A0
D = L.dirs13()
SEED = 777 if DEV else 20261013
CHEM = 1.594e-3                       # 1 kcal/mol in Eh
EPS_ABS, EPS_REL, TAME = 1e-4, 2e-2, 1.5
res = json.load(open(JS)) if os.path.exists(JS) else {}
st = pickle.load(open(PK, "rb")) if os.path.exists(PK) else dict(free=[], coupling=[])


def save():
    json.dump(res, open(JS, "w"), indent=1, default=float)
    pickle.dump(st, open(PK, "wb"))


# ------------------------------------------------------------------------------------------- coupling (as aggfull_learn)
class Lev:
    def __init__(self, rec, C, offs):
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
        Gk, Gp = self.gen(k), self.gen(kp)
        X = self.h * self.coup(k, kp)
        ek, Uk = np.linalg.eigh(Gk)
        ep, Up = np.linalg.eigh(Gp)
        Xt = np.conj(np.swapaxes(Up, -1, -2)) @ X @ Uk
        diff = ep[:, :, None] - ek[:, None, :]
        fac = -1j * np.exp(-0.5j * (ep[:, :, None] + ek[:, None, :])) * np.sinc(diff / (2 * np.pi))
        dW = Up @ (Xt * fac) @ np.conj(np.swapaxes(Uk, -1, -2)) * np.exp(1j * self.phi0)
        Wk = Uk @ (np.exp(-1j * ek)[..., None] * np.conj(np.swapaxes(Uk, -1, -2))) * np.exp(1j * self.phi0)
        Wp = Up @ (np.exp(-1j * ep)[..., None] * np.conj(np.swapaxes(Up, -1, -2))) * np.exp(1j * self.phi0)
        Dd = np.zeros_like(dW)
        for t in range(nsteps):
            Dd += np.linalg.matrix_power(Wp, nsteps - 1 - t) @ dW @ np.linalg.matrix_power(Wk, t)
        return np.linalg.matrix_power(Wk, nsteps), np.linalg.matrix_power(Wp, nsteps), Dd


def data_from_below(below, k, kp, B_new):
    b = B_new // below.B
    Tk, Tp, Dd = below.response(k, kp, b)
    z = np.zeros((1, 3))
    phi0 = float(np.angle(np.mean(np.linalg.eigvals(below.response(z, z, b)[0][0]))))
    Gk = ac.gen_from_step(Tk, phi0)
    Gp = ac.gen_from_step(Tp, phi0)
    X = ac.duhamel_invert(Dd, Gk, Gp, phi0)
    return X / (B_new * TH)


def sample_pairs(rng, n, h, k_s):
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


def learn_coupling(below, B, rng):
    h = B * TH
    k_s = 1 / A0
    k, kp, w, ns = sample_pairs(rng, 3000, h, k_s)
    Gam = data_from_below(below, k, kp, B)
    kh, kph, wh, nsh = sample_pairs(rng, 1500, h, k_s)
    Gh = data_from_below(below, kh, kph, B)
    ku = np.pi * (2 * rng.random((2000, 3)) - 1) / h
    kpu = np.pi * (2 * rng.random((2000, 3)) - 1) / h
    G0d = data_from_below(below, np.zeros((1, 3)), np.zeros((1, 3)), B)[0]
    dep = np.linalg.norm(Gam - G0d[None], axis=(1, 2))
    wr = w / (dep + 1e-3 * np.median(dep[:ns]))
    tried = []
    for name in ("point", "star", "cube"):
        offs = ac.stencil(name)
        C = ac.fit_coupling(k * h, kp * h, Gam, offs, wts=wr * np.sqrt(np.median(dep[:ns])))
        Gf = ac.eval_coupling(C, offs, kh[:nsh] * h, kph[:nsh] * h)
        G00 = ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]
        D00 = Gh[:nsh]
        err_abs = float(np.mean(np.linalg.norm(Gf - D00, axis=(1, 2))) / np.linalg.norm(G00) * np.sqrt(2))
        err_rel = float(np.mean(np.linalg.norm(Gf - D00, axis=(1, 2))) / np.mean(np.linalg.norm(D00 - G00[None], axis=(1, 2))))
        tame = float(np.max(np.linalg.norm(ac.eval_coupling(C, offs, ku * h, kpu * h), axis=(1, 2)) / np.sqrt(2)))
        ok = err_abs <= EPS_ABS and err_rel <= EPS_REL and tame <= TAME
        tried.append(dict(stencil=name, err_abs=err_abs, err_rel=err_rel, tame_max=tame, ok=bool(ok)))
        print(f"   coupling B={B} {name}: abs {err_abs:.2e} rel {err_rel:.2e} tame {tame:.4f} {'OK' if ok else ''}", flush=True)
        if ok:
            phys = ac.physics_of_coupling(C, offs, h)
            return dict(B=B, h=h, C=C, offs=offs, name=name), dict(
                chosen=name, candidates=tried, physics={kk: [v.real, v.imag] if isinstance(v, complex) else v for kk, v in phys.items()})
    return None, dict(chosen=None, candidates=tried)


# ------------------------------------------------------------------------------------------- free level
def learn_free(prev, b, rng_seed):
    B = prev.B * b
    cons = L.taylor6(prev, b, D, kmax=1.5)
    rng = np.random.default_rng(rng_seed)
    kfac = 1 / A0
    kap, Gs, sw, phi0 = L.samples(prev, b, kfac * TH * B, rng, 8000, 0.25, 3e-4)
    rec = dict(B=B, b=b, tried=[], c6_over_dirac=float(np.mean(cons[:, 2] / (1 / (16 * (B * TH) ** 5)))),
               c4_over_dirac=float(np.mean(cons[:, 1] / (-1 / (8 * (B * TH) ** 3)))),
               m_kin_from_c2=float(np.mean(1 / (2 * cons[:, 0] * B * TH))))
    chosen = None
    for w in (4, 5):
        r = L.fit(kap, Gs, sw, phi0, w, "T", cons=cons, D=D)
        cand = Level("learned", 3, "pauli2", 2, B, TH, rule=r, b=b, w=w)
        ck = azc.kcheck(prev, cand, b, np.random.default_rng(B + w + rng_seed), 1000, kfac * TH, sig=max(1.0, 2.0 / (kfac * TH * B)),
                        k0_cells=0.5 * kfac * TH * B, energy_mode="state")
        ok = ck["packet"] <= 1e-3 and ck["energy_rel"] <= 3e-4 and ck["rest_rel"] <= 3e-4 and ck["bz_rel90"] <= 0.2 and ck["neg_gen"] <= 0.01
        g, o = L.band_of(r)
        c2 = np.array([-0.5 * np.sum(g * (o @ d) ** 2) for d in D]); c4 = np.array([np.sum(g * (o @ d) ** 4) / 24 for d in D])
        rec["tried"].append(dict(w=w, check=ck, ok=bool(ok), max_c2_dev=float(np.max(np.abs(c2 / cons[:, 0] - 1))),
                                 max_c4_dev=float(np.max(np.abs(c4 / cons[:, 1] - 1)))))
        print(f"   free B={B} T w={w}: " + ", ".join(f"{k} {v:.2e}" for k, v in ck.items() if isinstance(v, float)) + f" ok={ok}", flush=True)
        if ok:
            chosen = (w, r, cand)
            break
    rM = L.fit(kap, Gs, sw, phi0, 4, "M", cons=cons[:, :2], D=D)
    g, o = L.band_of(rM)
    c2 = np.array([-0.5 * np.sum(g * (o @ d) ** 2) for d in D]); c4 = np.array([np.sum(g * (o @ d) ** 4) / 24 for d in D])
    rec["markov_w4"] = dict(max_c2_dev=float(np.max(np.abs(c2 / cons[:, 0] - 1))), max_c4_dev=float(np.max(np.abs(c4 / cons[:, 1] - 1))),
                            ratio_needed=float(np.mean(np.abs(cons[:, 1]) / cons[:, 0])), ratio_min_positive=1 / 12)
    print(f"   Markov B={B} w=4: c2 dev {rec['markov_w4']['max_c2_dev']:.2e} c4 dev {rec['markov_w4']['max_c4_dev']:.2e} "
          f"(|c4|/c2 needed {rec['markov_w4']['ratio_needed']:.2e}, positive >= 1/12)", flush=True)
    return chosen, rec


# ------------------------------------------------------------------------------------------- atom check
def atom_1s(rec_free, coup, B, budget=400):
    h = B * TH
    N = int(round(16 * A0 / h))
    T = mc.kinetic_grid(rec_free, N, h)
    f = mc.MultiField(N, h, [(0, 0, 0)])
    C, offs = coup["C"], coup["offs"]
    g0 = complex(np.trace(ac.eval_coupling(C, offs, np.zeros((1, 3)), np.zeros((1, 3)))[0]) / 2).real
    m = mc.Mol(T, f, gamma0=g0, sym=(1, 1, 1))
    X, E, r, conv, its = mc.solve_until(m, mc.guess_mol(N, h, [(0, 0, 0)], "1s")[:, None], budget, tag=f"1s B={B}")
    return dict(N=N, box_a0=N * h / A0, E_Eh=E[0] / mc.ALPHA ** 2, residual=r[0], converged=conv, its=its, gamma0=g0)


if __name__ == "__main__":
    t0 = time.time()
    rec512 = mc.level_rec(512)
    c512 = [l for l in pickle.load(open(os.path.join(mc.OUT, "aggfull_learn.pkl"), "rb")) if l["B"] == 512][0]
    levels = {512: (rec512, dict(C=c512["C"], offs=c512["offs"]))}
    if "atom_512" not in res:
        res["atom_512"] = atom_1s(rec512, levels[512][1], 512)
        save()
    for r in st["free"]:
        levels[r["B"]] = (r, [c for c in st["coupling"] if c["B"] == r["B"]][0])

    def try_level(Bprev, b):
        B = Bprev * b
        key = f"B{B}_from{Bprev}"
        if key in res and res[key].get("done"):
            return res[key]
        print(f"== candidate B={B} (a0/{int(round(A0 / (B * TH)))}... cell {B * TH / A0:.4f} a0) from B={Bprev}, b={b}", flush=True)
        rp, cp = levels[Bprev]
        prev = Level("learned", 3, "pauli2", 2, Bprev, TH, rule=ac.genrule(rp), b=None, w=rp.get("w"))
        chosen, rec = learn_free(prev, b, SEED + B)
        out = dict(free=rec)
        if chosen is None:
            out.update(done=True, admissible=False, reason="no free candidate passes")
            res[key] = out; save(); return out
        w, rule, cand = chosen
        frec = dict(B=B, b=b, w=w, desc="pauli2", C=2, K=rule.K, offs=rule.offs, phi0=rule.phi0, from_B=Bprev)
        below = Lev(rp, cp["C"], cp["offs"])
        crec, cinfo = learn_coupling(below, B, np.random.default_rng(SEED + 2 * B))
        out["coupling"] = cinfo
        if crec is None:
            out.update(done=True, admissible=False, reason="no admissible coupling")
            res[key] = out; save(); return out
        at = atom_1s(frec, crec, B)
        ref = res[f"atom_{Bprev}"]["E_Eh"]
        at["dE_vs_below_mEh"] = 1e3 * (at["E_Eh"] - ref)
        at["dE_vs_exact_mEh"] = 1e3 * (at["E_Eh"] + 0.5)
        adm = abs(at["E_Eh"] - ref) <= CHEM
        out.update(atom=at, admissible=bool(adm), done=True)
        res[f"atom_{B}"] = at
        print(f"   atom check B={B}: 1s {at['E_Eh']:.8f} Eh, vs below {at['dE_vs_below_mEh']:+.4f} mEh, vs -1/2 "
              f"{at['dE_vs_exact_mEh']:+.4f} mEh -> {'ADMISSIBLE' if adm else 'preview only'}", flush=True)
        st["free"] = [x for x in st["free"] if x["B"] != B] + [dict(frec, admissible=bool(adm))]
        st["coupling"] = [x for x in st["coupling"] if x["B"] != B] + [crec]
        levels[B] = (frec, crec)
        res[key] = out
        save()
        return out

    # the automatic choice: largest block first, from a0/8
    chosen_B = None
    for b in (4, 2):
        o = try_level(512, b)
        if o.get("admissible"):
            chosen_B = 512 * b
            break
        if time.time() - t0 > 420:
            print("checkpoint: run again", flush=True); sys.exit(0)
    res["chosen_two_atom_level_B"] = chosen_B
    save()
    # repeat from the chosen level (b = 2), to see whether a coarser level is admissible; else flagged 'preview'
    if chosen_B is not None and not res.get(f"B{2 * chosen_B}_from{chosen_B}", {}).get("done"):
        if time.time() - t0 > 400:
            print("checkpoint: run again", flush=True); sys.exit(0)
        try_level(chosen_B, 2)
    res["seconds_last_call"] = time.time() - t0
    save()
    print("LEARN DONE", flush=True)
