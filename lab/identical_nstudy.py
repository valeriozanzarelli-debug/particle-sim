"""TASK 3 — how many identical balls make a photon?

  python lab/identical_nstudy.py free  RULES Ns [seeds]      e.g.  free A,C 50,100,200 q,1,2,3
  python lab/identical_nstudy.py slit  RULES Ns [seeds]      e.g.  slit A 50,100,200,500 q

free: width error at zR, 2zR, 5zR.  'q' = regular (quantile) start: error vs theory.  Integer seed = random
      start (balls drawn from |ψ|²): error vs theory AND vs the start sample's own width (the part the rule is
      responsible for; the rest is the sampling noise of the start, ~1/sqrt(2N)).
      A = HDW (implicit midpoint for q, explicit Verlet for random starts), C = ψ carried (m = 8, P = 2).
slit: coarse correlation and fringe-height error vs the exact RS flux, balls in dark zones, order swaps.
Runs in parallel on 4 processes.
"""
import sys, time, numpy as np
from multiprocessing import Pool


def free_job(args):
    rule, N, seed = args
    import identical_free as F
    r = "Ai" if (rule == "A" and seed is None) else rule
    kw = dict(m=8) if rule == "C" else {}
    if r == "Ai":                                                     # implicit step: smaller when the balls are dense
        kw["dzmax"] = min(0.2, 0.2 * (2000 / N) ** 2)
    t = time.time()
    res, info = F.run(r, N, seed, tmax=300, **kw)
    return rule, N, seed, res, info, time.time() - t


def slit_job(args):
    rule, N, seed = args
    import identical_slit as S
    t = time.time()
    xl, info = S.run_slit(rule, N, seed, 1.0, tmax=540)
    if xl is None or info.get("fail") or not np.all(np.isfinite(xl)):
        return rule, N, seed, None, info, time.time() - t
    vals, line, base = S.evaluate(xl, 1.0, N, seed, info, "")
    return rule, N, seed, (vals, base), info, time.time() - t


if __name__ == "__main__":
    mode, rules = sys.argv[1], sys.argv[2].split(",")
    Ns = [int(a) for a in sys.argv[3].split(",")]
    seeds = [None if a == "q" else int(a) for a in (sys.argv[4].split(",") if len(sys.argv) > 4 else ["q"])]
    jobs = [(r, N, s) for r in rules for N in Ns for s in seeds]
    jobs.sort(key=lambda j: -j[1])                                   # big ones first
    with Pool(4) as pool:
        for out in pool.imap_unordered(free_job if mode == "free" else slit_job, jobs):
            rule, N, seed, res, info, dt = out
            st = "quantile" if seed is None else f"seed {seed}"
            if mode == "free":
                txt = "  ".join(("FAIL" if not np.isfinite(e) else f"{100*e:+.3f}% (own {100*e2:+.3f}%)")
                                for z, w, wt, e, e2 in res)
                print(f"free {rule} N={N:6d} {st:9s} start width/w0-1 {100*(info['w_start_ratio']-1):+.3f}%  "
                      f"zR, 2zR, 5zR: {txt}  swaps {info['swaps']}  [{dt:.0f}s]", flush=True)
            else:
                if res is None:
                    print(f"slit {rule} N={N:6d} {st:9s} FAILED at z={info['z']:.1f} swaps {info['swaps']} [{dt:.0f}s]", flush=True)
                else:
                    (c, r, mx, cp, df), base = res
                    print(f"slit {rule} N={N:6d} {st:9s} corr {c:.4f}  fringe rms {100*r:.1f}% max {100*mx:.1f}%  "
                          f"dark {100*df:.1f}% (th 2.0%)  swaps {info['swaps']}  [{dt:.0f}s]\n    {base.strip()}", flush=True)
