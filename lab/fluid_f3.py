"""fluid_f3.py - F3 production: Frenkel-Kontorova kinks (Lorentz contraction, E = gamma E0, limiting speed,
kink-antikink collisions). Writes out/fluid_f3.json, out/fluid_f3.npz. Protocol in out/fluid_prereg.json (F3)."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fluid_core import fk_acc, fk_energy_density, kink_profile, fit_kink

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
t0 = time.time()


def single(h, v, T=60.0, Lx=400.0):
    n = int(Lx / h); x = (np.arange(n) - n / 2) * h; x0 = -0.4 * Lx / 2 if v > 0 else 0.0
    g = 1 / np.sqrt(1 - v * v)
    phi = kink_profile(x, x0, g); dphi = -v * g * 2 / np.cosh(g * (x - x0))
    dt = 0.2 * h; a = fk_acc(phi, h); ts, Xs, ws = [], [], []; t = 0.0
    every = int(round(1 / dt))
    for s in range(int(round(T / dt))):
        dphi += 0.5 * dt * a; phi += dt * dphi; a = fk_acc(phi, h); dphi += 0.5 * dt * a; t += dt
        if (s + 1) % every == 0 and t > 10:
            X, w = fit_kink(x, phi); ts.append(t); Xs.append(X); ws.append(w)
    vm = float(np.polyfit(ts, Xs, 1)[0]); gm = 1 / np.sqrt(1 - vm * vm)
    e = fk_energy_density(phi, dphi, h); Ek = float(e[np.abs(x - Xs[-1]) < 15].sum())
    snap = (x[np.abs(x - Xs[-1]) < 8] - Xs[-1], phi[np.abs(x - Xs[-1]) < 8])
    return dict(v0=v, v=vm, gamma=float(gm), w=float(np.mean(ws)), w_std=float(np.std(ws)), E=Ek,
                E_total=float(e.sum())), snap


def driven(h, F, alpha=0.05, T=150.0, Lx=600.0):
    n = int(Lx / h); x = (np.arange(n) - n / 2) * h; x0 = 0.0  # BUGFIX after first run: was -0.4 Lx, kink moves left and hit the free end
    phi = kink_profile(x, x0, 1.0); dphi = np.zeros_like(phi)
    # background tilt: uniform state phi = asin F far away; kink between asin F and 2 pi + asin F
    phi = phi + np.arcsin(F)
    dt = 0.2 * h; a = fk_acc(phi, h, F, alpha, dphi); ts, Xs = [], []; t = 0.0
    for s in range(int(round(T / dt))):
        dphi += 0.5 * dt * a; phi += dt * dphi; a = fk_acc(phi, h, F, alpha, dphi); dphi += 0.5 * dt * a; t += dt
        if s % int(round(1 / dt)) == 0 and t > T - 50:
            c = np.where(np.diff(np.floor((phi - np.pi - np.arcsin(F)) / (2 * np.pi))) != 0)[0]  # BUGFIX: unique crossing
            ts.append(t); Xs.append(x[c[0]])
    v = float(abs(np.polyfit(ts, Xs, 1)[0]))  # kink moves towards -x for F > 0
    vms = 1 / np.sqrt(1 + (4 * alpha / (np.pi * F)) ** 2)
    return dict(F=F, v=v, v_MS=float(vms), rel=float(v / vms - 1))


def kak(h, v, Tafter=1000.0, d=10.0):
    T = d / v + Tafter; Lx = 2 * (T + 50); n = int(Lx / h); x = (np.arange(n) - n / 2) * h; g = 1 / np.sqrt(1 - v * v)
    phi = 4 * np.arctan(np.exp(g * (x + d))) + 4 * np.arctan(np.exp(-g * (x - d))) - 2 * np.pi
    dphi = -v * g * 2 / np.cosh(g * (x + d)) - v * g * 2 / np.cosh(g * (x - d))
    dt = 0.2 * h; a = fk_acc(phi, h); hist = []; every = int(round(5 / dt)); t = 0.0
    for s in range(int(round(T / dt))):
        dphi += 0.5 * dt * a; phi += dt * dphi; a = fk_acc(phi, h); dphi += 0.5 * dt * a; t += dt
        if s % every == 0:
            e = fk_energy_density(phi, dphi, h); hist.append((t, float(e[np.abs(x) < 10].sum())))
    e = fk_energy_density(phi, dphi, h); E = float(e.sum())
    w = int(round(10 / h)); cs = np.convolve(e, np.ones(2 * w + 1), "same")
    lumps = []; c2 = cs.copy()
    while c2.max() > 2.0:
        k = int(np.argmax(c2)); lumps.append((float(x[k]), float(cs[k]))); c2[max(0, k - 2 * w):k + 2 * w + 1] = 0
    El = sum(l[1] for l in lumps); frac_rad = 1 - El / E
    if frac_rad > 0.5: cls = "annihilated into phonons"
    elif len(lumps) == 1 and abs(lumps[0][0]) < 10: cls = "bound (breather)"
    elif len(lumps) == 2: cls = "pass-through"
    else: cls = "other"
    hist = np.array(hist)
    rate = None
    if cls == "bound (breather)":
        late = hist[:, 0] > T - 600
        sl = np.polyfit(hist[late, 0], np.log(hist[late, 1]), 1)[0]
        rate = float(-sl)
    return dict(h=h, v=v, E=E, frac_rad=float(frac_rad), lumps=lumps, cls=cls, breather_decay_rate=rate,
                E_centre_end=float(hist[-1, 1])), hist


out = dict(single={}, driven=[], kak=[])
snaps = {}
for h in [0.25, 0.5, 1.0]:
    rows = []
    for v in [0, 0.2, 0.4, 0.6, 0.8, 0.9]:
        r, sn = single(h, v); rows.append(r); snaps["h%g_v%g" % (h, v)] = np.array(sn)
    w0, E0 = rows[0]["w"], rows[0]["E"]
    for r in rows:
        r["wg_rel"] = r["w"] * r["gamma"] / w0 - 1; r["Eg_rel"] = r["E"] / (r["gamma"] * E0) - 1
    out["single"]["h%g" % h] = rows
for F in [0.05, 0.1, 0.2, 0.4, 0.8]:
    out["driven"].append(driven(0.25, F))
hists = {}
for h in [0.25, 1.0]:
    for v in [0.05, 0.1, 0.2, 0.4, 0.8]:
        r, hi = kak(h, v); out["kak"].append(r); hists["h%g_v%g" % (h, v)] = hi
m = out["single"]["h0.25"]
ok = [r for r in m if r["v0"] <= 0.8]
out["exams"] = dict(
    E3a=bool(all(abs(r["wg_rel"]) < 0.05 for r in ok)),
    E3b=bool(all(abs(r["Eg_rel"]) < 0.05 for r in ok)),
    E3c=bool(all(d["v"] < 1 for d in out["driven"]) and out["driven"][-1]["v"] > 0.95
             and all(abs(d["rel"]) < 0.02 for d in out["driven"])),
    E3d_prediction_held=bool(all(k["cls"] == "pass-through" for k in out["kak"] if k["h"] == 0.25)
                             and all(k["cls"] == "bound (breather)" for k in out["kak"] if k["h"] == 1.0 and k["v"] <= 0.4)
                             and not any(k["cls"] == "annihilated into phonons" for k in out["kak"])))
out["seconds"] = time.time() - t0
json.dump(out, open(os.path.join(OUT, "fluid_f3.json"), "w"), indent=1)
np.savez_compressed(os.path.join(OUT, "fluid_f3.npz"), **snaps, **{"hist_" + k: v for k, v in hists.items()})
print(json.dumps(out, indent=1))
