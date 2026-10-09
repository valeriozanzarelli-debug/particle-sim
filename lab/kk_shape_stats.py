"""kk_shape_stats.py - TEST 2: level statistics (integrable vs chaotic hidden shapes) and one scarred mode.
Criteria: lab/out/kk_shape_prereg.json (T2). Writes lab/out/kk_shape_stats.json and kk_shape_stats.npz."""
import json, time, sys
import numpy as np
from scipy.sparse.linalg import eigsh
import kk_shape_core as S

t0 = time.time()
PART = sys.argv[1] if len(sys.argv) > 1 else "all"      # lap | walk | walkdisk | scar | all (each < 5 min)
try:
    out = json.load(open("out/kk_shape_stats.json"))
    save = dict(np.load("out/kk_shape_stats.npz"))
except (OSError, ValueError):
    out, save = {}, {}
out["reference"] = {"Poisson": 0.3863, "GOE": 0.5307, "GUE": 0.5996}

# ---------------- 2a. Laplacian (Dirichlet), desymmetrized quarter shapes ----------------
A = 12000
if PART not in ("lap", "all"):
    shapes_on = False
else:
    shapes_on = True
gold = (1 + np.sqrt(5)) / 2
shapes = {}
shapes["quarter disk"] = S.shape_mask("qdisk", A)[0]
shapes["quarter stadium"] = S.shape_mask("qstadium", A, 1.0)[0]
b = int(round(np.sqrt(A / gold))); a = int(round(gold * b))
m = np.zeros((a + 4, b + 4), bool); m[2:a + 2, 2:b + 2] = True
shapes["golden rectangle"] = m
for nm, mask in (shapes.items() if shapes_on else []):
    L, _ = S.laplacian(mask)
    lam = np.sort(eigsh(L.tocsc(), k=460, sigma=0, which="LM")[0])
    k = np.sqrt(lam)
    sel = lam[50:][k[50:] < 0.6]
    r, _ = S.ratio_stat(sel)
    sp_ = S.unfold(np.sqrt(sel))
    out[f"laplacian {nm}"] = dict(levels=len(sel), r_mean=r, sites=int(mask.sum()))
    save[f"lap_{nm.replace(' ', '_')}"] = sp_
    print(f"[Laplacian] {nm:16s} {len(sel)} levels  <r~> = {r:.3f}  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
if shapes_on:
  pq = out["laplacian quarter disk"]["r_mean"] < 0.45 and out["laplacian golden rectangle"]["r_mean"] < 0.45
  ps = 0.49 <= out["laplacian quarter stadium"]["r_mean"] <= 0.57
  out["laplacian PASS"] = bool(pq and ps)
  print("  PASS (integrable < 0.45, stadium in [0.49, 0.57]):", out["laplacian PASS"])

# ---------------- 2b. the walk (Dirac billiard with mass wall) on the quarter stadium ----------------
for w, label in [(w_, l_) for w_, l_ in ((0, "sharp wall (pre-registered)"), (8, "soft wall w=8 (post-hoc)"))
                 if PART in ("walk", "all", f"walk{w_}")]:
    mask, info = S.shape_mask("qstadium", 6000 if w == 0 else 4000, 1.0, margin=4 + w)   # soft: 4000 sites (time limit), same as the disk control
    th = S.soft_wall(mask, w)
    E = S.rest_masses(th, Emax=0.8)
    E = E[E > 0.05]
    r_raw, _ = S.ratio_stat(E)
    mid, spl, worst = S.pair_up(E)
    r_mid, _ = S.ratio_stat(mid)
    out[f"walk qstadium {label}"] = dict(levels=len(E), r_raw=r_raw, r_pairs_merged=r_mid,
                                         median_split_over_spacing=float(np.median(spl[:-1] / np.diff(mid))),
                                         worst_split_over_gap=worst)
    save[f"walk_w{w}"] = S.unfold(mid)
    print(f"[walk] quarter stadium, {label}: {len(E)} levels, <r~> raw = {r_raw:.3f}, doubler pairs merged = {r_mid:.3f}"
          f"  (median pair splitting / spacing {out[f'walk qstadium {label}']['median_split_over_spacing']:.3f})"
          f"  ({time.time() - t0:.0f}s)")
    sys.stdout.flush()
# integrable-ish control for the walk: quarter disk, soft wall (smaller: 4000 sites)
if PART in ("walkdisk", "all"):
  mask, info = S.shape_mask("qdisk", 4000, margin=12)
  E = S.rest_masses(S.soft_wall(mask, 8), Emax=0.8); E = E[E > 0.05]
  mid, spl, worst = S.pair_up(E)
  out["walk qdisk soft wall w=8 (post-hoc)"] = dict(levels=len(E), r_raw=S.ratio_stat(E)[0], r_pairs_merged=S.ratio_stat(mid)[0])
  save["walkdisk_w8"] = S.unfold(mid)
  print(f"[walk] quarter disk, soft wall: {len(E)} levels, <r~> raw = {S.ratio_stat(E)[0]:.3f}, pairs merged = {S.ratio_stat(mid)[0]:.3f}")

# ---------------- 2c. scar: full stadium, Dirichlet Laplacian ----------------
if PART in ('scar', 'all'):
  mask, info = S.shape_mask("stadium", 8000, 1.0)
  L, pts = S.laplacian(mask)
  lam, V = eigsh(L.tocsc(), k=400, sigma=0, which="LM")
  o = np.argsort(lam); lam, V = lam[o], V[:, o]
  x = pts[:, 0] - info["c1"]; y = pts[:, 1] - info["c2"]
  a = info["a"]
  best = None
  enh = []
  for i in range(len(lam)):
      k = np.sqrt(lam[i]); hw = np.pi / (2 * k)
      tube = (np.abs(y) < hw) & (np.abs(x) > a)
      frac = tube.mean()
      wgt = np.sum(V[tube, i]**2)
      e = wgt / frac
      enh.append(e)
      if best is None or e > best[1]:
          best = (i, e, k)
  enh = np.array(enh)
  i, e, k = best
  img = np.full(mask.shape, np.nan); img[pts[:, 0], pts[:, 1]] = V[:, i]**2
  save["scar_img"] = img; save["scar_mask"] = mask
  # a typical (median-enhancement) mode for comparison
  j = int(np.argsort(enh)[len(enh) // 2])
  img2 = np.full(mask.shape, np.nan); img2[pts[:, 0], pts[:, 1]] = V[:, j]**2
  save["typical_img"] = img2
  out["scar"] = dict(mode_index=int(i), k=float(k), enhancement_caps_tube=float(e), median_enhancement=float(np.median(enh)),
                     n_modes_above_3=int(np.sum(enh >= 3)), PASS=bool(e >= 3), stadium=dict(a=info["a"], r=info["r"]))
  print(f"[scar] full stadium, 400 modes: best mode #{i} (k = {k:.3f}) has {e:.2f}x the uniform weight in the "
        f"caps part of the tube around the unstable horizontal orbit; median {np.median(enh):.2f}; modes >= 3x: "
        f"{int(np.sum(enh >= 3))}  PASS {e >= 3}")
# control for the pair-merging step (post-hoc): what merging does to levels with NO repulsion
rg = np.random.default_rng(S.SEED)
Ep = np.cumsum(rg.exponential(size=20000))
E0 = np.cumsum(rg.exponential(size=10000)); Dd = np.sort(np.concatenate([E0, E0 + 0.16 * rg.random(10000)]))
out["merge_artifact"] = dict(poisson=S.ratio_stat(Ep)[0], poisson_pair_merged=S.ratio_stat(S.pair_up(Ep)[0])[0],
                             true_doublets_merged=S.ratio_stat(S.pair_up(Dd)[0])[0])
print("[control] pair-merging a Poisson sequence gives <r~> =", round(out["merge_artifact"]["poisson_pair_merged"], 3),
      "; true doublets of a Poisson sequence, merged:", round(out["merge_artifact"]["true_doublets_merged"], 3))
out["seconds"] = time.time() - t0
with open("out/kk_shape_stats.json", "w") as fh:
    json.dump(out, fh, indent=1, default=float)
np.savez_compressed("out/kk_shape_stats.npz", **save)
print(f"done {time.time() - t0:.0f}s")
