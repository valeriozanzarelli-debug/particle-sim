"""zoomknot_run.py CONFIG [bare|learned]  - relax one pre-registered configuration (lab/out/zoomknot_prereg.json)
with checkpoints. Each invocation runs <= ~8.5 min; call again until the result json exists.
Stop rule (pre-registered): rel. E change over 5 outer steps < 2e-8, or 45 outer steps, or 30 min total.
Writes lab/out/zoomknot_ckpt/<tag>.pt and, when finished, lab/out/zoomknot_<tag>.json."""
import sys, os, json, math, time, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc
import zoomknot_core as zc

torch.set_num_threads(1)
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
PRE = json.load(open(os.path.join(OUT, "zoomknot_prereg.json")))
name = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else "bare"
tag = name + ("_learned" if mode == "learned" else "")
levels = PRE["configs"][name]
k2 = k4 = None
if mode == "learned":
    K = json.load(open(os.path.join(OUT, "zoomknot_learn.json")))
    k2, k4 = K["k2"], K["k4"]
H = zc.Hier(levels, k2, k4)
os.makedirs(os.path.join(OUT, "zoomknot_ckpt"), exist_ok=True)
ck = os.path.join(OUT, "zoomknot_ckpt", tag + ".pt")
res_path = os.path.join(OUT, f"zoomknot_{tag}.json")
if os.path.exists(res_path):
    print("done already", res_path); sys.exit(0)
t0 = time.time()
if os.path.exists(ck):
    S = torch.load(ck)
    us, hist, wall_used = S["us"], S["hist"], S["wall"]
else:
    us, hist, wall_used = H.init_u(kc.radial_bvp(60.0)), [], 0.0
    E0 = float(H.energy(us)); print(f"{tag}: sites {H.sites}  start E/12pi^2 {E0/(12*math.pi**2):.6f}", flush=True)
done = False
while True:
    left = min(510 - (time.time() - t0), 1800 - wall_used - (time.time() - t0))
    if left < 30:
        break
    us, hist, done = zc.relax(H, us, maxit=45, tol=4e-9, wall=left, hist=hist)
    print(f"{tag}: outer {len(hist)}  E/12pi^2 {hist[-1]/(12*math.pi**2):.9f}  t {time.time()-t0:.0f}s", flush=True)
    if done or len(hist) >= 45:
        done = True; break
wall_used += time.time() - t0
torch.save({"us": us, "hist": hist, "wall": wall_used}, ck)
if done or wall_used >= 1800 - 30:
    m = zc.measure(H, us, os.path.join(OUT, "press_profiles.npz"))
    ex = PRE["exact_reference"]
    m.update(config=name, mode=mode, outer_steps=len(hist), wall_s=wall_used,
             stop=("converged" if (len(hist) >= 11 and abs(hist[-6] - hist[-1]) / hist[-1] < 2e-8) else
                   ("max_steps" if len(hist) >= 45 else "wall")),
             E_over_12pi2=m["E"] / (12 * math.pi**2), hist_tail=[h / (12 * math.pi**2) for h in hist[-6:]])
    m["err_E"] = m["E_over_12pi2"] / ex["E_over_12pi2"] - 1
    m["err_p0"] = m["p0_central"] / ex["p0"] - 1
    m["err_p0_fb"] = m["p0_fb"] / ex["p0"] - 1
    m["err_Lam"] = m["Lam"] / ex["Lambda"] - 1
    json.dump(m, open(res_path, "w"), indent=1)
    print(json.dumps({k: v for k, v in m.items() if k not in ("profile", "levels")}, indent=1))
else:
    print(f"{tag}: checkpoint saved, wall used {wall_used:.0f}s, call again")
