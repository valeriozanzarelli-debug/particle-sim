"""locale_r3.py - REGOLA LOCALE, production run 3: the light exams M1-M4 and T3 (pre-registration
lab/out/locale_prereg.json). Writes lab/out/locale_r3.json (checkpoint after every candidate: rerun to continue).
Run: OMP_NUM_THREADS=1 timeout 590 python lab/locale_r3.py   (exit code 3 = run me again)"""
import os, sys, json, time, itertools
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from locale_core import *
from locale_cands import light, LIGHTS

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
CK = os.path.join(OUT, "locale_r3.json")
t0 = time.time()
BUDGET = 520
res = json.load(open(CK)) if os.path.exists(CK) else {"M": {}, "verdict": []}


def check(name, ok, msg):
    res["verdict"].append((name, bool(ok)))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {msg}", flush=True)


def save():
    with open(CK, "w") as f:
        json.dump(res, f, indent=1, default=float)


def dirset(dim):
    return dirs_for(2, 46) if dim == 2 else fib_sphere(320, half=True)


ref_cost = {}
for name in LIGHTS:
    if name in res["M"]:
        continue
    if time.time() - t0 > BUDGET:
        save(); print("budget reached: rerun", flush=True); sys.exit(3)
    w, dens, dim = light(name)
    D = dirset(dim)
    r = {"T": w.T, "density": dens, "dim": dim}
    # M1
    k = 2 * PI / 1000
    om, _ = light_band(w, k * D)
    e1 = float(np.max(np.abs(om / k - 1)))
    r["M1"] = e1
    check(f"M1 {name}", e1 <= 1e-5, f"|omega/ck - 1| at 1000 cells/lambda = {e1:.1e}")
    # M2 static mode
    k16 = 2 * PI / 16
    W = w.bloch(k16 * D)
    lam, V = np.linalg.eig(W)
    j = np.argmin(np.abs(lam - 1), axis=1)
    dev = float(np.max(np.abs(lam[np.arange(len(j)), j] - 1)))
    v0 = V[np.arange(len(j)), :, j]
    v0 = v0 / np.linalg.norm(v0, axis=1)[:, None]
    # the static mode of Maxwell is F parallel to k (longitudinal); in the F = (Fx, Fy, Fz) basis the overlap with k_hat
    ov = float(np.max(1 - np.abs(np.einsum("ni,ni->n", v0.conj(), D / np.linalg.norm(D, axis=1)[:, None])) ** 2))
    r["M2"] = dict(eig_dev=dev, overlap_defect_N16=ov)
    check(f"M2 {name}", dev <= 1e-12, f"static eigenvalue |lambda - 1| = {dev:.1e}; longitudinal overlap defect at 16 cells = {ov:.1e}")
    # M3 N_1%
    n1 = {key: n_for_tol(w, D, 0.01, key, lo=3.0, hi=400.0) for key in ("phase", "group", "diffraction")}
    Nw = max(n1.values())
    r["N1"] = n1; r["N1_worst"] = Nw
    nshift = sum(f[0] == "shift" for f in w.pre)
    r["factors_per_time"] = nshift / w.T
    r["sites_per_lambda_d"] = dens * Nw ** dim
    r["cost"] = r["sites_per_lambda_d"] * r["factors_per_time"]
    # errors on a fixed grid of N (for the figure)
    r["grid"] = {}
    for N in (8, 12, 16, 24, 32, 48, 64):
        r["grid"][str(N)] = worst(light_metrics(w, N, D))
    # superluminal check of the two moving bands for |k| <= 2 pi / N_1%
    kmax = 2 * PI / Nw if np.isfinite(Nw) else 2 * PI / 400
    vmax = 0.0
    for sgn in (+1, -1):
        for s in np.linspace(0.1, 1, 10):
            K = s * kmax * D
            w0, v0b = light_band(w, K, sign=sgn)
            g = np.zeros((len(K), dim))
            hh = 1e-6
            for a in range(dim):
                e = np.zeros(3); e[a] = hh
                g[:, a] = (light_band(w, K + e, ref=(w0, v0b))[0] - light_band(w, K - e, ref=(w0, v0b))[0]) / (2 * hh)
            vmax = max(vmax, float(np.max(np.linalg.norm(g, axis=1))))
    r["vmax_moving_resolved"] = vmax
    # whole zone (information)
    rng = np.random.default_rng(SEED + 3)
    kr = 4.3 if name.startswith("L1") else PI          # hexagonal zone reaches |k| = 4 pi / 3
    KZ = rng.uniform(-kr, kr, (2000, 3)); KZ[:, dim:] = 0
    vz, gz = group_velocities(w, KZ)
    okz = gz > 1e-3
    supz = (vz > 1 + 1e-6) & okz
    kk = np.repeat(np.linalg.norm(KZ, axis=1)[:, None], vz.shape[1], 1)
    r["zone"] = dict(vmax=float(vz[okz].max()), frac=float(supz.mean()),
                     kmin=float(kk[supz].min()) if supz.any() else None)
    res["M"][name] = r
    print(f"     M3 {name}: N_1% phase {n1['phase']:.1f}, group {n1['group']:.1f}, diffraction {n1['diffraction']:.1f} "
          f"-> worst {Nw:.1f} cells/lambda; sites/lambda^{dim} {r['sites_per_lambda_d']:.3g}; "
          f"shift factors per unit time {r['factors_per_time']:.2f}; cost {r['cost']:.3g}", flush=True)
    print(f"     (info) {name}: max |v_g| moving bands for |k| <= 2pi/N_1%: {vmax:.6f} c; whole zone max {r['zone']['vmax']:.4f} c, "
          f"fraction above c {r['zone']['frac']:.2e}, smallest |k| above c {r['zone']['kmin']}", flush=True)
    save()

# verdict of M3 (debt #4 reduced?) after all candidates
if "M3_done" not in res:
    for name in LIGHTS:
        r = res["M"][name]
        ref = res["M"]["L0_3D" if r["dim"] == 3 else "L0_2D"]
        ok = (r["N1_worst"] <= 20) and (r["cost"] <= ref["cost"]) and (r["vmax_moving_resolved"] <= 1 + 1e-6)
        check(f"M3 {name}", ok, f"worst N_1% = {r['N1_worst']:.1f} (<= 20?), cost {r['cost'] / ref['cost']:.3f} x frozen "
              f"(<= 1?), max v_g resolved {r['vmax_moving_resolved']:.6f} (<= 1 + 1e-6?)")
    res["M3_done"] = True
    save()

# T3: order
if "T3" not in res:
    if time.time() - t0 > BUDGET - 120:
        save(); print("budget reached: rerun", flush=True); sys.exit(3)
    res["T3"] = {}
    for name, lo, hi in (("L0_3D", 3, 5), ("L0_2D", 3, 5), ("L2_3D", 12, 20), ("L2_2D", 12, 20)):
        g = res["M"][name]["grid"]
        ratio = g["16"]["phase"] / g["32"]["phase"]
        res["T3"][name] = ratio
        check(f"T3 order {name}", lo <= ratio <= hi, f"phase error ratio 16 -> 32 cells = {ratio:.2f}")
    # positive unit-step words (2-D)
    X, Y = CUBE_AXES[0], CUBE_AXES[1]
    D2 = dirs_for(2, 46)
    best = []
    for L in (4, 6, 8, 10, 12):
        for half in itertools.product((0, 1), repeat=L // 2):
            word = list(half) + list(half[::-1])
            if sum(word) != L // 2 or word[0] == 1:
                continue
            w = make_walk("light", [[X if c == 0 else Y for c in word]], dim=2)
            best.append((n_for_tol(w, D2, 0.01, "diffraction"), "".join("xy"[c] for c in word)))
    best.sort()
    frozen = res["M"]["L0_2D"]["N1"]["diffraction"]
    res["T3"]["positive_words_best"] = best[:10]
    check("T3 positive words", best[0][0] >= 0.98 * frozen,
          f"best positive palindromic word '{best[0][1]}' diffraction N_1% = {best[0][0]:.2f} vs frozen {frozen:.2f} "
          f"({len(best)} words tried)")
    save()

npass = sum(v for _, v in res["verdict"])
print(f"=== {npass} PASS / {len(res['verdict']) - npass} FAIL ({time.time() - t0:.0f} s this piece)", flush=True)
