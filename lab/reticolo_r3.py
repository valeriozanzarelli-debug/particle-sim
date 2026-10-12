"""reticolo_r3.py - R3: ONE world with a 'gamma' pulse at level 0 and a 'radio' pulse at level K, the levels in
between existing only as a nested chimney around the gamma.  Usage:
   python reticolo_r3.py ladder     # K = 20: build, cells/memory per level, short run (2048 gamma ticks), frames
   python reticolo_r3.py both       # K = 10: run 8 radio ticks = 8192 gamma ticks, both pulses checked
Writes lab/out/reticolo_r3_<mode>.json, lab/out/reticolo_r3_<mode>.npz (frames data)."""
import sys, os, json, time, pickle
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reticolo_core as rc
from reticolo_core import Level, World, packet, band_vectors, kgrid, OUT, single_level_world

B = 16
LAM = 24
GSIG = (24.0, 72.0)          # gamma envelope (along x, across) in h0
RSIG = 36.0                  # radio envelope in h_K
RANG = 60.0
TOL_DETAIL = 1e-5
TOL_TOP = 1e-4
EVERY = 8
BUDGET = 470.0               # s per invocation, then checkpoint (each script <= 10 min)
PAR_MARGIN = 2               # coarse blocks around the parents of the finer level's blocks


def dil(c, r, d=2):
    offs = np.array(list(np.ndindex(*([2 * r + 1] * d)))) - r
    return np.unique((c[:, None, :] + offs[None]).reshape(-1, d), axis=0)


def box_blocks(g0, shape):
    lo = np.floor_divide(np.asarray(g0), B)
    hi = np.floor_divide(np.asarray(g0) + np.asarray(shape) - 1, B)
    return np.stack(np.meshgrid(*[np.arange(a, b + 1) for a, b in zip(lo, hi)], indexing="ij"), -1).reshape(-1, 2)


def gamma_dense(sig=GSIG, shp=(384, 1024)):
    F = packet(shp, 2, LAM, 0.0, [shp[0] / 2, shp[1] / 2], sig, amp=1.0)
    return np.array([-shp[0] // 2, -shp[1] // 2]), F


def radio_dense(center_cells):
    shp = (512, 512)
    F = packet(shp, 2, LAM, RANG, [shp[0] / 2, shp[1] / 2], RSIG, amp=0.5)
    g0 = np.array(center_cells) - 256
    return g0, F



def build(K, radio_center, gsig=GSIG, gshape=(384, 1024)):
    g0g, Fg = gamma_dense(gsig, gshape)
    g0r, Fr = radio_dense(radio_center)
    coords = [box_blocks(g0g, Fg.shape[:2])]
    for k in range(1, K):
        coords.append(dil(np.floor_divide(coords[-1], 2), PAR_MARGIN))
    top = np.concatenate([box_blocks(g0r, Fr.shape[:2]), dil(np.floor_divide(coords[-1], 2), PAR_MARGIN)])
    coords.append(top)
    levels = [Level(k, 2, B, coords[k]) for k in range(K + 1)]
    w = World(levels)
    levels[K].set_dense(g0r, Fr)
    for i in range(K - 1, -1, -1):
        w.pairs[i].prolong_all()
    levels[0].set_dense(g0g, Fg, add=True)
    w.sync_all_up()
    w.tol = [TOL_DETAIL] * K + [TOL_TOP]
    w.halo = 1
    w.par_margin = PAR_MARGIN
    for i in range(K + 1):
        w.regrid(i)
    for i in range(K - 1, -1, -1):     # second pass so the chimney is consistent top-down
        w.regrid(i)
    w.regrid_every = {i: EVERY for i in range(K + 1)}
    return w, (g0g, Fg), (g0r, Fr)


def level_table(w):
    rows = []
    for L in w.levels:
        rows.append(dict(k=L.k, h=L.h, blocks=int(L.nb), cells=int(L.nb * L.Bd),
                         bytes=int(L.F.nbytes + L.tmp.nbytes)))
    return rows


def gamma_exact(Fg, g0g, T, box=(1024, 2048)):
    """exact uniform-lattice evolution of the gamma (spectral, periodic box large enough), returned as a
    function of global level-0 cell indices"""
    Lx, Ly = box
    big = np.zeros((Lx, Ly, 3), complex)
    ox, oy = (Lx - Fg.shape[0]) // 2, (Ly - Fg.shape[1]) // 2
    big[ox:ox + Fg.shape[0], oy:oy + Fg.shape[1]] = Fg
    bg0 = np.array([g0g[0] - ox, g0g[1] - oy])
    Fk = np.fft.fftn(big, axes=(0, 1)).reshape(-1, 3)
    KA = kgrid((Lx, Ly), 2).reshape(-1, 3)
    ch = 1 << 16
    for a in range(0, len(KA), ch):
        vs, om = band_vectors(KA[a:a + ch])
        v = vs[0]
        c = np.einsum("ni,ni->n", v.conj(), Fk[a:a + ch]) * np.exp(-1j * om * T)
        Fk[a:a + ch] = c[:, None] * v
    Ft = np.fft.ifftn(Fk.reshape(Lx, Ly, 3), axes=(0, 1))

    def at(G):
        ix = np.mod(G[:, 0] - bg0[0], Lx)
        iy = G[:, 1] - bg0[1]
        ok = (iy >= 0) & (iy < Ly)
        out = np.zeros((len(G), 3), complex)
        out[ok] = Ft[ix[ok], iy[ok]]
        return out
    return at, float((np.abs(Ft) ** 2).sum())


def composite_dense(w, lo, shp):
    """level-0 resolution picture of the composite field: each cell from the finest level holding it"""
    grids = np.stack(np.meshgrid(*[np.arange(n) + o for n, o in zip(shp, lo)], indexing="ij"), -1).reshape(-1, 2)
    out = np.zeros((len(grids), 3), complex)
    todo = np.ones(len(grids), bool)
    for L in w.levels:
        if not todo.any():
            break
        g = np.floor_divide(grids[todo], 1 << L.k)
        fl, ok = L.lookup(g)
        idx = np.nonzero(todo)[0]
        out[idx[ok]] = L.flat()[fl[ok]]
        todo[idx[ok]] = False
    return out.reshape(tuple(shp) + (3,))


def gamma_check(w, Fg, g0g, T, highpass, box=(1024, 2048)):
    L0 = w.levels[0]
    G = L0.gcoords().reshape(-1, 2)
    m = L0.count_mask().reshape(-1)
    F = L0.flat()[:L0.nb * L0.Bd].copy()
    at, Eex = gamma_exact(Fg, g0g, T, box)
    Ex = at(G)
    if highpass:
        # remove the (locally smooth) radio field: dense box around the gamma, smooth window, high-pass |k| < 0.1
        xc = T
        lo = np.array([xc - 192, -1024])
        shp = (384, 2048)
        D = composite_dense(w, lo, shp)
        X = np.arange(shp[0]) + lo[0] + 0.5 - xc
        Y = np.arange(shp[1]) + lo[1] + 0.5
        win = (np.exp(-(X / 120.0) ** 8)[:, None] * np.exp(-(Y / 900.0) ** 8)[None, :])[..., None]
        Dk = np.fft.fftn(D * win, axes=(0, 1))
        K = kgrid(shp, 2)
        Dk[np.linalg.norm(K, axis=-1) < 0.1] = 0
        Dg = np.fft.ifftn(Dk, axes=(0, 1))
        grids = np.stack(np.meshgrid(np.arange(shp[0]) + lo[0], np.arange(shp[1]) + lo[1], indexing="ij"), -1)
        Ed = at(grids.reshape(-1, 2)).reshape(shp + (3,)) * win
        A, Bv = Dg.reshape(-1, 3), Ed.reshape(-1, 3)
    else:
        A, Bv = F[m], Ex[m]
    ov = np.vdot(Bv, A)
    fid = abs(ov) / np.sqrt(np.vdot(A, A).real * np.vdot(Bv, Bv).real)
    e_lvl0 = float((np.abs(F[m]) ** 2).sum())
    xcent = float((np.abs(F[m]) ** 2).sum(-1) @ (G[m, 0] + 0.5) / max(e_lvl0, 1e-300))
    return dict(fidelity=float(fid), phase=float(np.angle(ov)), amp_ratio=float(np.sqrt(np.vdot(A, A).real / np.vdot(Bv, Bv).real)),
                energy_level0_over_exact=e_lvl0 / Eex, x_centroid=xcent)


def render(w, center, half, npix=256, kmin=0):
    """picture of what the world holds at resolution >= level kmin: each pixel takes the finest level >= kmin
    holding it (nearest cell). Returns complex Fz and the level map"""
    xs = center[0] + (np.arange(npix) + 0.5) / npix * 2 * half - half
    ys = center[1] + (np.arange(npix) + 0.5) / npix * 2 * half - half
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    P = np.stack([X.reshape(-1), Y.reshape(-1)], 1)
    val = np.full(len(P), np.nan + 0j)
    lev = np.full(len(P), -1)
    for L in reversed(w.levels):
        if L.k < kmin:
            break
        g = np.floor(P / L.h).astype(np.int64)
        fl, ok = L.lookup(g)
        Fz = L.flat()[fl, 2]
        val[ok] = Fz[ok]
        lev[ok] = L.k
    return val.reshape(npix, npix), lev.reshape(npix, npix)


def block_rects(L, center, half):
    out = []
    for c in L.coords:
        x0, y0 = c[0] * B * L.h, c[1] * B * L.h
        s = B * L.h
        if x0 + s < center[0] - half or x0 > center[0] + half or y0 + s < center[1] - half or y0 > center[1] + half:
            continue
        out.append((x0, y0, s))
    return np.array(out) if out else np.zeros((0, 3))


def snapshot_frames(w, center, tag, wide=None):
    frames = {}
    if wide is not None:                      # (k, centre, half): one extra wide view
        k, c, half = wide
        v, lev = render(w, c, half, npix=512, kmin=k)
        frames["wide_v"], frames["wide_lev"], frames["wide_c"], frames["wide_half"] = v, lev, np.array(c), np.array(half)
        frames["wide_k"] = np.array(k)
        for M in w.levels:
            frames["wide_rect_%02d" % M.k] = block_rects(M, c, half)
    for L in w.levels:
        half = 80 * L.h
        v, lev = render(w, center, half, kmin=L.k)
        frames["v%02d" % L.k] = v
        frames["lev%02d" % L.k] = lev
        frames["half%02d" % L.k] = np.array(half)
        for M in w.levels:
            if M.k >= L.k - 1 and M.k <= L.k + 1:
                frames["rect%02d_%02d" % (L.k, M.k)] = block_rects(M, center, half)
    frames["center"] = np.array(center)
    np.savez_compressed(os.path.join(OUT, "reticolo_r3_%s.npz" % tag), **frames)


def main(mode):
    t0 = time.time()
    res = {"mode": mode, "sync_every": rc.SYNC_EVERY, "W_band": rc.W_BAND}
    if mode == "ladder":
        K, NT = 20, 2048
        rcen = (-40, -30)          # radio centre in level-20 cells; the gamma (origin) sits ~1.4 sigma from it
    else:
        K, NT = 10, 8 * 1024
        rcen = (-260, 0)           # radio far from the gamma (7 sigma): checked separately
    state = os.path.join(OUT, "reticolo_r3_%s_state.pkl" % mode)
    if os.path.exists(state):
        with open(state, "rb") as fh:
            st = pickle.load(fh)
        w, alone, hist, maxcells, maxbytes, res, Etot0, g0g, Fg, g0r, Fr, run_acc = st
        print("resumed at t =", w.t, flush=True)
    else:
        if mode == "ladder":
            w, (g0g, Fg), (g0r, Fr) = build(K, rcen)
        else:
            w, (g0g, Fg), (g0r, Fr) = build(K, rcen, gsig=(24.0, 144.0), gshape=(384, 2048))
        res["K"] = K
        res["gamma_sigma"] = [24.0, 72.0] if mode == "ladder" else [24.0, 144.0]
        res["build_s"] = time.time() - t0
        res["levels_t0"] = level_table(w)
        res["cells_t0"] = w.cells()
        res["bytes_t0"] = sum(r["bytes"] for r in res["levels_t0"])
        print("built", K, "cells", w.cells(), "MB %.1f" % (res["bytes_t0"] / 1e6), "%.0fs" % res["build_s"],
              [(r["k"], r["blocks"]) for r in res["levels_t0"]], flush=True)
        if mode == "ladder":
            snapshot_frames(w, (0.0, 0.0), "ladder_t0")
        Etot0 = w.energy()
        alone = single_level_world(K, 2, B, g0r, Fr, TOL_TOP, every=EVERY) if mode == "both" else None
        hist = []
        maxcells, maxbytes = w.cells(), res["bytes_t0"]
        run_acc = 0.0
    t1 = time.time()
    chunk = 256
    while w.t < NT:
        w.step_fine(chunk)
        cells = w.cells()
        by = sum(L.F.nbytes for L in w.levels)
        maxcells, maxbytes = max(maxcells, cells), max(maxbytes, by)
        hist.append((w.t, cells, [int(L.nb) for L in w.levels]))
        if (w.t // chunk) % 8 == 0:
            print("t", w.t, "cells", cells, "lvl0 blocks", w.levels[0].nb, "%.0fs" % (time.time() - t1), flush=True)
        if time.time() - t0 > BUDGET and w.t < NT:
            with open(state, "wb") as fh:
                pickle.dump((w, alone, hist, maxcells, maxbytes, res, Etot0, g0g, Fg, g0r, Fr,
                             run_acc + time.time() - t1), fh)
            print("checkpoint at t =", w.t, "- run again to continue", flush=True)
            return
    run_acc += time.time() - t1
    if os.path.exists(state):
        os.remove(state)
    res["run_s"] = run_acc
    res["ticks_level0"] = w.t
    res["ticks_per_level"] = [int(L.ticks) for L in w.levels]
    res["s_per_level0_tick"] = res["run_s"] / NT
    res["max_cells"] = int(maxcells)
    res["max_bytes"] = int(maxbytes)
    res["levels_end"] = level_table(w)
    res["hist"] = hist
    res["rebuilds"] = w.nrebuild
    full_sync = (w.t % (1 << (K + 1)) == 0)
    res["all_levels_synchronised"] = bool(full_sync)
    res["energy_total_ratio"] = (w.energy() / Etot0) if full_sync else None
    res["gamma"] = gamma_check(w, Fg, g0g, w.t, highpass=(mode == "ladder"),
                               box=(1024, 2048) if mode == "ladder" else (1024, 4096))
    res["gamma"]["x_expected"] = float(w.t)
    if mode == "ladder":
        res["note_time"] = ("Berger-Oliger: at t=0 every level ticks once and waits ahead in time; levels >= 12 "
                            "are at their next tick time (2^k h0/c), level 20 at 2^20: the radio has made 1 tick")
        snapshot_frames(w, (float(w.t), 0.0), "ladder_tend")
    else:
        alone.step_top(w.levels[-1].ticks)
        A = alone.levels[0]
        lo = np.array(rcen) - 200
        Da = A.get_dense(lo, (400, 400))
        Dc = w.levels[-1].get_dense(lo, (400, 400))
        res["radio_rel_L2_vs_alone"] = float(np.sqrt((np.abs(Dc - Da) ** 2).sum() / (np.abs(Da) ** 2).sum()))
        e = (np.abs(Da) ** 2).sum(-1)
        X = np.stack(np.meshgrid(np.arange(400) + lo[0] + 0.5, np.arange(400) + lo[1] + 0.5, indexing="ij"), -1)
        cen = (e[..., None] * X).sum((0, 1)) / e.sum()
        e0 = (np.abs(Fr) ** 2).sum(-1)
        X0 = np.stack(np.meshgrid(np.arange(512) + g0r[0] + 0.5, np.arange(512) + g0r[1] + 0.5, indexing="ij"), -1)
        cen0 = (e0[..., None] * X0).sum((0, 1)) / e0.sum()
        dvec = np.array([np.cos(np.radians(RANG)), np.sin(np.radians(RANG))])
        res["radio_ticks"] = int(w.levels[-1].ticks)
        res["radio_centroid_move_cells"] = float((cen - cen0) @ dvec)
        res["radio_speed_over_c"] = res["radio_centroid_move_cells"] / w.levels[-1].ticks
        ec = (np.abs(Dc) ** 2).sum()
        res["radio_energy_ratio"] = float(ec / (np.abs(Fr) ** 2).sum())
        h = w.levels[-1].h
        snapshot_frames(w, (float(w.t), 0.0), "both_tend",
                        wide=(K, (0.5 * (rcen[0] * h + w.t), 0.0), 330 * h))
    res["wall_s"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "reticolo_r3_%s.json" % mode), "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("hist", "levels_t0", "levels_end")}, indent=1))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "frames":       # re-draw the t = 0 ladder frames only
        w, _, _ = build(20, (-40, -30))
        snapshot_frames(w, (0.0, 0.0), "ladder_t0", wide=(20, (-40 * 2.0 ** 20, -30 * 2.0 ** 20), 200 * 2.0 ** 20))
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else "ladder")
