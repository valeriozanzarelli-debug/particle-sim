"""agg_zz_3d.py - the automatic agglomerator on the 3-D zig-zag walk (chain A: theta = a0/4096, or chain B: a0/512).
Levels are learned one after the other (each from the LEARNED level below), checkpointed in lab/out/agg_zz_3d_<chain>.pkl
(re-run the script to continue: every call does at most ~8 minutes of work).
Observable: hydrogen energies -> momentum density of a 1s-like state of size a0 (k_s = 1/a0) up to the 1s target cell
a0/16, then of size 2 a0 (n = 2) for the level with cell <= a0/8. Learner: plane-wave least squares (learn_h for the
Pauli generator, learn_k for 4-component kernels); check in k-space (packets + state-weighted energy error).
Run: python lab/agg_zz_3d.py A   (or B)"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import json, pickle, sys, time
import numpy as np
from agg_zz_core import *
from agg_zz_core import GenRule

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
A0 = 1 / zz.ALPHA
chain = sys.argv[1] if len(sys.argv) > 1 else "A"
SUB = {"A": 256, "B": 32}[chain]
TH = A0 / 16 / SUB
EPS_E = float(sys.argv[2]) if len(sys.argv) > 2 else 3e-4
TARGETS = [(SUB, 1 / A0), (2 * SUB, 1 / (2 * A0))] if chain == "A" else [(SUB, 1 / A0)]   # (max cell, k_s = 1/(n a0))
NC = 20
PK = os.path.join(OUT, f"agg_zz_3d_{chain}.pkl")
t_start = time.time()
BUDGET = 480

if os.path.exists(PK):
    st = pickle.load(open(PK, "rb"))
else:
    st = dict(levels=[], log=[], done=False, seed=20261012)
rng = np.random.default_rng(st["seed"] + 101 * len(st["levels"]))


def log(s):
    print(s, flush=True)
    st["log"].append(s)


def rebuild(lvrec):
    if lvrec.get("phi0") is not None:
        r = GenRule(lvrec["K"], lvrec["offs"], 3, lvrec["C"], lvrec["phi0"])
    else:
        r = Rule(lvrec["K"], lvrec["offs"], 3, lvrec["C"])
    return Level("learned", 3, lvrec["desc"], lvrec["C"], lvrec["B"], TH, rule=r, b=lvrec["b"], w=lvrec["w"],
                 check=lvrec["check"])


def disp_vs_fine(lv):
    """positive-band motion energy of the learned level vs the fine walk's own (Bloch, B/2 two-tick steps), along x and
    the diagonal, kappa in [0, pi/2]"""
    out = {}
    for name, u in (("x", np.array([1.0, 0, 0])), ("diag", np.ones(3) / np.sqrt(3))):
        kap = np.linspace(0, np.pi / 2, 33)[:, None] * u[None, :]
        pc = band_phase(lv.rule.mult_at(kap), lv.desc)
        Wf = np.linalg.matrix_power(zz.bloch_batch(kap / lv.B, TH), lv.B // 2)
        pf = band_phase(conv_k(Wf, "dirac4", lv.desc), lv.desc)
        ec = -np.angle(np.exp(1j * (pc - pc[0]))) / lv.B
        ef = -np.angle(np.exp(1j * (pf - pf[0]))) / lv.B
        out[name] = dict(kappa=np.linalg.norm(kap, axis=1).tolist(), eps_learned=ec.tolist(), eps_fine=ef.tolist())
    out["rest_phase_learned"] = float(band_phase(lv.rule.mult_at(np.zeros((1, 3))), lv.desc)[0])
    out["rest_phase_fine_theta_B"] = float(np.angle(np.exp(-1j * TH * lv.B)))
    if lv.desc == "pauli2":
        K = lv.rule.K
        tr = 0.5 * (K[0, 0] + K[1, 1])
        out["spin_part_rel"] = float(np.max(np.abs(K - tr[None, None] * np.eye(2)[:, :, None])) / np.max(np.abs(K)))
    return out


if not st["done"]:
    prev = Level("fine3", 3, "dirac4", 4, 1, TH, period=2)
    for rec in st["levels"]:
        prev = rebuild(rec)
    while time.time() - t_start < BUDGET:
        tgt = [t for t in TARGETS if prev.B < t[0]]
        if not tgt:
            st["done"] = True
            break
        Bmax, kfac = tgt[0]
        ksr = lambda B: -kfac * TH * B          # learn on wave vectors drawn like the observable state
        log(f"level {len(st['levels']) + 1}: from B = {prev.B} ({prev.desc}), cell target <= {Bmax}, k_obs = {kfac} m")
        lv, tried = auto_level(prev, rng, NC, eps_packet=1e-3, eps_energy=EPS_E, T_obs=1000, k_obs=kfac * TH,
                               bs=(16, 8, 4, 3, 2), ws=(1, 2, 3, 4), bmax_total=Bmax, log=log, ks_rule=ksr,
                               learner="auto", packets="phys", energy_mode="state")
        if lv is None:
            log(f"  no admissible candidate: hierarchy stops at B = {prev.B} for this observable")
            if len(tgt) > 1:
                # the next observable (larger orbit) may still allow a coarser level
                TARGETS.remove(tgt[0])
                continue
            st["done"] = True
            st["stopped_no_candidate"] = prev.B
            break
        d = disp_vs_fine(lv)
        rec = dict(B=lv.B, b=lv.b, w=lv.w, desc=lv.desc, C=lv.C, K=lv.rule.K, offs=lv.rule.offs, check=lv.check,
                   phi0=getattr(lv.rule, "phi0", None), rule_type=type(lv.rule).__name__,
                   k_obs_over_m=kfac, ks=tried[-1]["ks"], candidates=tried, disp=d, cell_lambdaC=lv.B * TH,
                   cell_over_a0=lv.B * TH / A0, seconds=time.time() - t_start)
        st["levels"].append(rec)
        dev = np.max(np.abs(np.array(d["x"]["eps_learned"]) - np.array(d["x"]["eps_fine"])))
        log(f"  CHOSEN B={lv.B} ({lv.B * TH:.3f} lambda_C = a0/{A0 / (lv.B * TH):.1f}) {lv.desc} b={lv.b} w={lv.w}: "
            f"packet {lv.check['packet']:.1e} energy {lv.check['energy_rel']:.1e} rest {lv.check['rest_rel']:.1e}; "
            f"rest phase {d['rest_phase_learned']:+.5f} vs fine {d['rest_phase_fine_theta_B']:+.5f}; max |eps-eps_fine| (x) "
            f"{dev:.2e}" + (f"; spin part {d['spin_part_rel']:.1e}" if "spin_part_rel" in d else ""))
        prev = lv
        pickle.dump(st, open(PK, "wb"))
    pickle.dump(st, open(PK, "wb"))

if st["done"]:
    js = dict(chain=chain, theta=TH, levels=[{k: v for k, v in r.items() if k not in ("K", "offs")} for r in st["levels"]],
              log=st["log"], stopped_no_candidate=st.get("stopped_no_candidate"))
    json.dump(js, open(os.path.join(OUT, f"agg_zz_3d_{chain}.json"), "w"), indent=1, default=float)
    log("hierarchy complete")
else:
    log("budget reached: run again to continue")
