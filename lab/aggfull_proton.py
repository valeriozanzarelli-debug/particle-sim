"""aggfull_proton.py - P4: the PROTON as a ball of the same rule, agglomerated by the same machine.

The proton is the 3-D zig-zag walk (lab/zz_core.py, frozen) with its own flip theta_p = m_p c a_p / hbar on its own fine
lattice a_p (input: m_p/m_e = 1836.15267343, a mass of a nucleus). Its fine cell is chosen so that the coarse cells of
the electron's chain A (a0/16 = 256 x 0.0334561 lambda_C and a0/8) are exact powers of two of it:
   a_p = (a0/16) / 2^19   ->  theta_p = m_p a_p = 0.029995   (kinetic-mass lattice error theta_p^2/3 = 3.0e-4, as the
                                                               electron's fine walk: 3.7e-4)
The SAME automatic agglomerator of lab/agg_zz_core.py (auto_level: learn the generator from the level below, check
packets + state-weighted energy + tameness, choose the coarsest, repeat from the LEARNED level) builds the hierarchy up to
the cell a0/16 and then a0/8. Only change: blocks b in {16, 8, 4, 2} (powers of two, so that the proton's coarse cells
land exactly on the electron's; declared). Observable: the atom's relative momentum (state of size a0, then 2 a0), the
same k_obs as chain A. Checkpoint lab/out/aggfull_proton.pkl, results lab/out/aggfull_proton.json.
Run (repeat until 'hierarchy complete'): timeout 590 python lab/aggfull_proton.py
"""
import os, sys, json, pickle, time
import numpy as np
import aggfull_core as ac
import agg_zz_core as azc
from agg_zz_core import Level, auto_level, band_phase, conv_k
import zz_core as zz

OUT = ac.OUT
A0 = ac.A0
H16 = ac.TH_A * 256                       # electron coarse cell a0/16 (lambda_C)
BP16 = 2 ** 19
AP = H16 / BP16                           # proton fine cell
THP = ac.MP * AP                          # proton flip angle
DEV = os.environ.get("AGGFULL_DEV") == "1"
TAG = "dev_aggfull_proton" if DEV else "aggfull_proton"
PK = os.path.join(OUT, TAG + ".pkl")
TARGETS = [(BP16, zz.ALPHA * ac.MP), (2 * BP16, 0.5 * zz.ALPHA * ac.MP)]     # (max cell, k_obs = alpha m_p (/n))
KCAP = 1 / 16
_orig_ssk = azc.sample_state_k
azc.sample_state_k = lambda rng, n, d, k_s: _orig_ssk(rng, n, d, min(k_s, KCAP))
_orig_kcheck = azc.kcheck


def _kcheck_capped(prev, cand, b, rng, T_obs, k_obs, sig=4.0, nk=3000, k0_cells=0.5, energy_mode="window"):
    """same check, with the observable capped at 0.5 per cell of the candidate (packets and energy)"""
    ko = min(k_obs, KCAP / cand.B)
    return _orig_kcheck(prev, cand, b, rng, T_obs, ko, sig=max(1.0, 2.0 / (ko * cand.B)), nk=nk, k0_cells=0.5 * ko * cand.B,
                        energy_mode=energy_mode)


azc.kcheck = _kcheck_capped
t_start = time.time()
BUDGET = 480
st = pickle.load(open(PK, "rb")) if os.path.exists(PK) else dict(levels=[], log=[], done=False, seed=20261013)
rng = np.random.default_rng(st["seed"] + 101 * len(st["levels"]))


def log(s):
    print(s, flush=True)
    st["log"].append(s)


def rebuild(r):
    return Level("learned", 3, r["desc"], r["C"], r["B"], THP, rule=azc.GenRule(r["K"], r["offs"], 3, r["C"], r["phi0"]),
                 b=r["b"], w=r["w"], check=r["check"])


def kinetic_mass(lv):
    """kinetic mass of a learned level from the curvature of its generator along x and the diagonal (physical units)"""
    h = lv.B * AP
    out = {}
    for nm, u in (("x", np.array([1.0, 0, 0])), ("diag", np.ones(3) / np.sqrt(3))):
        kc = np.array([1 / 64, 1 / 32, 3 / 64, 1 / 16])[:, None] * u[None]   # cell units: the observable's scale (<= KCAP)
        G = lv.rule.gen_at(kc)
        g = np.real(0.5 * (G[:, 0, 0] + G[:, 1, 1]))
        g0 = np.real(0.5 * np.trace(lv.rule.gen_at(np.zeros((1, 3)))[0]))
        e = (g - g0) / h                                                 # energy per unit time (mc^2)
        kp = np.linalg.norm(kc, axis=1) / h
        c2 = np.polyfit(kp ** 2, e, 2)[-2]                               # e = k^2/(2 m_kin) + ...
        out[nm] = float(1 / (2 * c2))
    return out


if not st["done"]:
    prev = Level("fine3", 3, "dirac4", 4, 1, THP, period=2)
    for r in st["levels"]:
        prev = rebuild(r)
    while time.time() - t_start < BUDGET:
        tgt = [t for t in TARGETS if prev.B < t[0]]
        if not tgt:
            st["done"] = True
            break
        Bmax, kfac = tgt[0]
        ksr = lambda B: -kfac * AP * B
        log(f"level {len(st['levels']) + 1}: from B = {prev.B} ({prev.desc}), target {Bmax}, theta_p = {THP:.6f}")
        # observable of the proton chain: a proton moving at v = alpha c (k = alpha m_p), i.e. the SAME velocity scale as the
        # electron's chain A (k = alpha m_e): in fine cells k a_p = alpha theta_p = 2.2e-4, as chain A's 2.4e-4. The atom's
        # own proton momenta (k ~ 1/a0 = 1e-7 per proton fine cell) are below what a least-squares fit resolves in double
        # precision at the fine levels; they are an exam of the final level (kinetic mass at k << 1/h). When the
        # observable's k exceeds half the Brillouin cube of a level (coarse levels: the whole cube is non-relativistic),
        # it is capped at KCAP = 1/16 per cell (as chain A: the state of size a0 seen by cells a0/16; wrappers below).
        lv, tried = auto_level(prev, rng, 20, eps_packet=1e-3, eps_energy=3e-4, T_obs=1000, k_obs=kfac * AP,
                               bs=(16, 8, 4, 2), ws=(1, 2, 3, 4), bmax_total=Bmax, log=log, ks_rule=ksr,
                               learner="auto", packets="phys", energy_mode="state")
        if lv is None:
            log(f"  no admissible candidate: the proton hierarchy stops at B = {prev.B}")
            st["done"] = True
            st["stopped_no_candidate"] = prev.B
            break
        mk = kinetic_mass(lv) if lv.desc == "pauli2" else {}
        rec = dict(B=lv.B, b=lv.b, w=lv.w, desc=lv.desc, C=lv.C, K=lv.rule.K, offs=lv.rule.offs, check=lv.check,
                   phi0=getattr(lv.rule, "phi0", None), cell_lambdaC=lv.B * AP, m_kin=mk, candidates=tried,
                   seconds=time.time() - t_start)
        st["levels"].append(rec)
        log(f"  CHOSEN B={lv.B} (cell {lv.B * AP:.4e} lambda_C) {lv.desc} b={lv.b} w={lv.w}: packet {lv.check['packet']:.1e} "
            f"energy {lv.check['energy_rel']:.1e}; m_kin/m_p = " + ", ".join(f"{k} {v / ac.MP:.6f}" for k, v in mk.items()))
        prev = lv
        pickle.dump(st, open(PK, "wb"))
    pickle.dump(st, open(PK, "wb"))

if st["done"]:
    js = dict(theta_p=THP, a_p=AP, m_p=ac.MP, levels=[{k: v for k, v in r.items() if k not in ("K", "offs")} for r in st["levels"]],
              log=st["log"], stopped_no_candidate=st.get("stopped_no_candidate"))
    json.dump(js, open(os.path.join(OUT, TAG + ".json"), "w"), indent=1, default=float)
    log("hierarchy complete")
else:
    log("budget reached: run again")
