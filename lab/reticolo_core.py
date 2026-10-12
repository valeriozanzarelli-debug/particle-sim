"""reticolo_core.py - SPARSE HIERARCHICAL lattice for the zig-zag light rule ("il reticolo dell'universo").

Design (frozen with the pre-registration lab/reticolo_prereg.py, before any validation run):

* LEVELS.  Level k has cell size h_k = 2^k h0 and tick dt_k = h_k / c.  The massless spin-1 walker
  (F = E + icB, (S_k)_ij = -i eps_kij, src/zigzag/walk.rs) has no length of its own, so the SAME rule runs
  on every level: one tick = for each axis j (x,y[,z] on even ticks, reversed on odd ticks) split F on the
  eigenvectors of S_j, move the +1 part one cell along +j, the -1 part one cell along -j, keep the 0 part.
* SPARSE BLOCKS.  Each level stores only ACTIVE blocks of B^d cells (B = 16 in 2-D, 8 in 3-D), keyed by
  integer block coordinates; missing blocks are exactly zero (a shared zero "sentinel" block).
  Block (b) of level k+1 covers exactly 2^d blocks of level k: the octree / quadtree.
* NESTING.  Every active block of level k lies inside level k+1 with >= 1 coarse block of margin
  (needed by the interpolation stencil).  The coarsest level is open (zero outside).
* TIME.  Berger-Oliger recursion: advance(k) = tick(k); advance(k-1) twice; sync(k-1 -> k).
  A level-k tick lasts 2^k level-0 ticks.  Each level keeps its own tick parity (palindromic pair).
* INTERFACE ("fascia di cucitura", as docs/ZOOM_NODO.md):  at every sync (fine and coarse at the same time)
  (1) RESTRICTION: every coarse cell whose 4^d-cell cubic stencil lies in the fine INTERIOR is overwritten by
      the 4th-order midpoint interpolation of the fine field (weights (-1, 9, 9, -1)/16 per axis);
  (2) PROLONGATION: the fine BAND (the outer W = 4 fine cells of the fine region, = 2 coarse cells) is
      overwritten by tensor-cubic Lagrange interpolation of the coarse field (weights at +-1/4 of a coarse cell).
  Between syncs the band evolves with the fine rule (no time interpolation needed): zeros from outside the
  fine region can corrupt at most 2 band cells in 2 fine ticks, never the interior.
  Every point of space is counted ONCE in the composite energy: fine interior cells, plus coarse cells
  not entirely made of fine-interior children.
* REGRID (dynamic).  A level-k block is needed where the field has detail the coarser level cannot carry:
  max |F_k - P(F_{k+1})| > tol (coarsest level: max |F| > tol), plus the parents of the finer level's
  blocks; flagged blocks are dilated by `halo` blocks.  New blocks are filled by prolongation.
Units: h0 = 1, c = 1.  Fields are point samples (not amplitudes per cell): energy = sum |F|^2 h^d.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys
import numpy as np
import numba as nb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zz_core import bloch_batch  # noqa: E402  (two-tick Bloch matrix of the frozen zig-zag rule)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
S2 = 1.0 / np.sqrt(2.0)
W_BAND = 4
WP_EVEN = np.array([-0.0390625, 0.2734375, 0.8203125, -0.0546875])   # fine g = 2m : coarse nodes m-2..m+1
WP_ODD = np.array([-0.0546875, 0.8203125, 0.2734375, -0.0390625])    # fine g = 2m+1: coarse nodes m-1..m+2
WR = np.array([-1.0, 9.0, 9.0, -1.0]) / 16.0                        # coarse gc <- fine 2gc-1..2gc+2
OFF = 1 << 20
BITS = 21


# ------------------------------------------------------------------ kernels (single thread)
@nb.njit(cache=True)
def _shift(F, nbr, j, B, tmp):
    """one axis of the light rule on all active blocks. F (nb+1, B^d, 3); last block = zero sentinel."""
    nbk = nbr.shape[0]
    Bd = F.shape[1]
    j1 = (j + 1) % 3
    j2 = (j + 2) % 3
    st = 1
    for _ in range(j):
        st *= B
    s2 = 0.7071067811865476
    for b in range(nbk + 1):
        for s in range(Bd):
            a = F[b, s, j1]
            c = F[b, s, j2]
            tmp[b, s, 0] = (a - 1j * c) * s2
            tmp[b, s, 1] = (a + 1j * c) * s2
    lo = nbr[:, 2 * j]
    hi = nbr[:, 2 * j + 1]
    for b in range(nbk):
        bl = lo[b]
        bh = hi[b]
        for s in range(Bd):
            cj = (s // st) % B
            if cj > 0:
                p = tmp[b, s - st, 0]
            else:
                p = tmp[bl, s + (B - 1) * st, 0]
            if cj < B - 1:
                m = tmp[b, s + st, 1]
            else:
                m = tmp[bh, s - (B - 1) * st, 1]
            F[b, s, j1] = (p + m) * s2
            F[b, s, j2] = 1j * (p - m) * s2


@nb.njit(cache=True)
def _gather(dst, tgt, src_arr, src, w):
    """dst[tgt[i]] = sum_q w[i,q] src_arr[src[i,q]]  (flat (N,3) views)"""
    for i in range(tgt.shape[0]):
        a0 = 0j
        a1 = 0j
        a2 = 0j
        for q in range(src.shape[1]):
            ww = w[i, q]
            r = src[i, q]
            a0 += ww * src_arr[r, 0]
            a1 += ww * src_arr[r, 1]
            a2 += ww * src_arr[r, 2]
        dst[tgt[i], 0] = a0
        dst[tgt[i], 1] = a1
        dst[tgt[i], 2] = a2


@nb.njit(cache=True)
def _gather_fixed(dst, tgt, src_arr, src, w):
    """same as _gather with the same weights w[q] for every target (restriction)"""
    for i in range(tgt.shape[0]):
        a0 = 0j
        a1 = 0j
        a2 = 0j
        for q in range(src.shape[1]):
            ww = w[q]
            r = src[i, q]
            a0 += ww * src_arr[r, 0]
            a1 += ww * src_arr[r, 1]
            a2 += ww * src_arr[r, 2]
        dst[tgt[i], 0] = a0
        dst[tgt[i], 1] = a1
        dst[tgt[i], 2] = a2


@nb.njit(cache=True)
def _detail_max(F, src_arr, src, w, nbk, Bd):
    """per block max |F - P(coarse)| over all cells"""
    out = np.zeros(nbk)
    for i in range(nbk * Bd):
        a0 = 0j
        a1 = 0j
        a2 = 0j
        for q in range(src.shape[1]):
            ww = w[i, q]
            r = src[i, q]
            a0 += ww * src_arr[r, 0]
            a1 += ww * src_arr[r, 1]
            a2 += ww * src_arr[r, 2]
        d = abs(F[i, 0] - a0) ** 2 + abs(F[i, 1] - a1) ** 2 + abs(F[i, 2] - a2) ** 2
        b = i // Bd
        if d > out[b]:
            out[b] = d
    return np.sqrt(out)


# ------------------------------------------------------------------ helpers
def enc(c):
    c = np.asarray(c, np.int64) + OFF
    if np.any(c < 0) or np.any(c >= (1 << BITS)):
        raise ValueError("block coordinate out of the prototype key range (+-2^20 blocks per level)")
    key = c[:, 0].copy()
    for j in range(1, c.shape[1]):
        key = (key << BITS) | c[:, j]
    return key


def local_coords(B, d):
    s = np.arange(B ** d)
    return np.stack([(s // B ** j) % B for j in range(d)], axis=1)       # (Bd, d), x fastest


class Level:
    def __init__(self, k, d, B, coords, period=None):
        self.k, self.d, self.B = k, d, B
        self.Bd = B ** d
        self.h = float(2.0 ** k)
        self.period = np.zeros(d, np.int64) if period is None else np.asarray(period, np.int64)
        self.loc = local_coords(B, d)
        self.ticks = 0
        self.has_parent = False
        self.owned = None        # (nb, Bd) bool: cell value is a copy of finer interior (not counted)
        self.set_blocks(coords)

    # --- structure
    def _wrap(self, bc):
        bc = np.array(bc, np.int64, copy=True)
        for j in range(self.d):
            if self.period[j] > 0:
                bc[:, j] %= self.period[j]
        return bc

    def set_blocks(self, coords, init=None):
        coords = self._wrap(np.asarray(coords, np.int64).reshape(-1, self.d))
        keys = enc(coords)
        keys, idx = np.unique(keys, return_index=True)
        self.coords = coords[idx]
        self.keys = keys
        self.nb = len(keys)
        self.owned = None
        self.F = np.zeros((self.nb + 1, self.Bd, 3), complex)
        self.tmp = np.zeros((self.nb + 1, self.Bd, 2), complex)
        nbr = np.empty((self.nb, 2 * self.d), np.int64)
        for j in range(self.d):
            for sgn, col in ((-1, 2 * j), (1, 2 * j + 1)):
                c2 = self.coords.copy()
                c2[:, j] += sgn
                nbr[:, col] = self.find(c2)
        self.nbr = nbr
        self._band()

    def find(self, bc):
        bc = self._wrap(np.asarray(bc, np.int64).reshape(-1, self.d))
        if self.nb == 0:
            return np.zeros(len(bc), np.int64)
        key = enc(bc)
        i = np.searchsorted(self.keys, key)
        i = np.minimum(i, max(self.nb - 1, 0))
        ok = (self.nb > 0) & (self.keys[i] == key) if self.nb > 0 else np.zeros(len(key), bool)
        return np.where(ok, i, self.nb)

    def _band(self):
        """band mask: fine cells within W of any missing neighbour block (all 3^d directions)"""
        B, d, W = self.B, self.d, W_BAND
        band = np.zeros((self.nb, self.Bd), bool)
        for off in np.ndindex(*([3] * d)):
            o = np.array(off) - 1
            if not o.any():
                continue
            c2 = self.coords + o[None, :]
            miss = self.find(c2) == self.nb
            if not miss.any():
                continue
            cond = np.ones(self.Bd, bool)
            for j in range(d):
                if o[j] == -1:
                    cond &= self.loc[:, j] < W
                elif o[j] == 1:
                    cond &= self.loc[:, j] >= B - W
            band[miss] |= cond[None, :]
        self.band = band

    def gcoords(self):
        """global cell indices (nb, Bd, d)"""
        return self.coords[:, None, :] * self.B + self.loc[None, :, :]

    def centers(self):
        """physical cell centres in units of h0 (nb, Bd, d)"""
        return (self.gcoords() + 0.5) * self.h

    def lookup(self, g):
        """global cell indices (..., d) -> flat index into F.reshape(-1,3); missing -> sentinel cell"""
        sh = g.shape[:-1]
        g = g.reshape(-1, self.d)
        bc = np.floor_divide(g, self.B)
        lc = g - bc * self.B
        bi = self.find(bc)
        s = np.zeros(len(g), np.int64)
        for j in range(self.d):
            s += lc[:, j] * self.B ** j
        flat = np.where(bi < self.nb, bi * self.Bd + s, self.nb * self.Bd)
        return flat.reshape(sh), (bi < self.nb).reshape(sh)

    def flat(self):
        return self.F.reshape(-1, 3)

    # --- dynamics
    def tick(self):
        if self.nb == 0:
            self.ticks += 1
            return
        order = range(self.d) if self.ticks % 2 == 0 else range(self.d - 1, -1, -1)
        for j in order:
            _shift(self.F, self.nbr, j, self.B, self.tmp)
        self.ticks += 1

    def energy_density(self):
        return (np.abs(self.F[:self.nb]) ** 2).sum(-1)              # (nb, Bd)

    def count_mask(self):
        m = np.ones((self.nb, self.Bd), bool)
        if self.has_parent:
            m &= ~self.band
        if self.owned is not None:
            m &= ~self.owned
        return m

    # --- dense in/out
    def set_dense(self, g0, arr, add=False):
        """write dense array arr[x, y, (z), 3] whose cell [0,0] has global index g0, into ACTIVE blocks"""
        g0 = np.asarray(g0, np.int64)
        shp = np.array(arr.shape[:self.d])
        G = self.gcoords().reshape(-1, self.d) - g0
        inside = np.all((G >= 0) & (G < shp), axis=1)
        vals = arr[tuple(G[inside].T)]
        Ff = self.flat()
        idx = np.nonzero(inside)[0]
        if add:
            Ff[idx] += vals
        else:
            Ff[idx] = vals

    def get_dense(self, g0, shape):
        g0 = np.asarray(g0, np.int64)
        grids = np.meshgrid(*[np.arange(n) for n in shape], indexing="ij")
        G = np.stack(grids, -1) + g0
        fl, ok = self.lookup(G)
        out = self.flat()[fl]
        out[~ok] = 0
        return out


# ------------------------------------------------------------------ inter-level tables
def stencil_prolong(G, d):
    """for fine global indices G (N,d): coarse node indices (N, 4^d, d) and weights (N, 4^d)"""
    N = len(G)
    par = np.mod(G, 2)
    m = np.floor_divide(G, 2)
    start = np.where(par == 0, m - 2, m - 1)
    wax = np.where(par[..., None] == 0, WP_EVEN[None, None, :], WP_ODD[None, None, :])  # (N,d,4)
    offs = np.array(list(np.ndindex(*([4] * d))))            # (4^d, d)
    nodes = start[:, None, :] + offs[None, :, :]
    w = np.ones((N, len(offs)))
    for j in range(d):
        w *= wax[:, j, offs[:, j]]
    return nodes, w


def stencil_restrict(GC, d):
    offs = np.array(list(np.ndindex(*([4] * d))))
    nodes = 2 * GC[:, None, :] - 1 + offs[None, :, :]
    w = np.ones(len(offs))
    for j in range(d):
        w *= WR[offs[:, j]]
    return nodes, w


class Pair:
    """interface tables between fine level f and coarse level c = f+1"""

    def __init__(self, f, c):
        self.f, self.c = f, c
        d = f.d
        f.has_parent = True
        f._band()
        # (2) prolongation into the band
        G = f.gcoords().reshape(-1, d)
        bflat = np.nonzero(f.band.reshape(-1))[0]
        nodes, w = stencil_prolong(G[bflat], d)
        src, ok = c.lookup(nodes)
        self.p_tgt, self.p_src, self.p_w = bflat.astype(np.int64), src.astype(np.int64), w
        self.p_missing = int((~ok).sum())
        # (1) restriction from the fine interior
        interior = (~f.band).reshape(-1)
        GI = G[interior]
        if len(GI):
            GC = np.unique(np.floor_divide(GI, 2), axis=0)
            nodes, wr = stencil_restrict(GC, d)
            ffl, fok = f.lookup(nodes)
            iok = fok & interior[np.minimum(ffl, len(interior) - 1)]
            good = iok.all(axis=1)
            cfl, cok = c.lookup(GC)
            good &= cok
            self.r_tgt, self.r_src, self.r_w = cfl[good].astype(np.int64), ffl[good].astype(np.int64), wr
            # ownership: coarse cells whose 2^d children are all fine interior
            offs2 = np.array(list(np.ndindex(*([2] * d))))
            ch = 2 * GC[:, None, :] + offs2[None]
            chf, chok = f.lookup(ch)
            own = (chok & interior[np.minimum(chf, len(interior) - 1)]).all(axis=1) & cok
            owned = np.zeros(c.nb * c.Bd, bool)
            owned[cfl[own]] = True
        else:
            self.r_tgt = np.zeros(0, np.int64)
            self.r_src = np.zeros((0, 4 ** d), np.int64)
            self.r_w = np.ones(4 ** d) / 4 ** d
            owned = np.zeros(c.nb * c.Bd, bool)
        c.owned = owned.reshape(c.nb, c.Bd)

    def sync(self):
        _gather_fixed(self.c.flat(), self.r_tgt, self.f.flat(), self.r_src, self.r_w)
        _gather(self.f.flat(), self.p_tgt, self.c.flat(), self.p_src, self.p_w)

    def restrict_only(self):
        _gather_fixed(self.c.flat(), self.r_tgt, self.f.flat(), self.r_src, self.r_w)

    def prolong_all(self, mask_flat=None):
        """fill fine cells (all, or mask) by interpolation of the coarse level"""
        f, c = self.f, self.c
        G = f.gcoords().reshape(-1, f.d)
        idx = np.arange(len(G)) if mask_flat is None else np.nonzero(mask_flat)[0]
        if len(idx) == 0:
            return
        nodes, w = stencil_prolong(G[idx], f.d)
        src, _ = c.lookup(nodes)
        _gather(f.flat(), idx.astype(np.int64), c.flat(), src.astype(np.int64), w)

    def detail(self):
        f, c = self.f, self.c
        G = f.gcoords().reshape(-1, f.d)
        nodes, w = stencil_prolong(G, f.d)
        src, _ = c.lookup(nodes)
        return _detail_max(f.flat(), c.flat(), src.astype(np.int64), w, f.nb, f.Bd)


# ------------------------------------------------------------------ the world
class World:
    """levels[i] has k = kmin + i (contiguous).  levels[0] finest."""

    def __init__(self, levels):
        self.levels = levels
        for a, b in zip(levels[:-1], levels[1:]):
            assert b.k == a.k + 1
        self.pairs = [Pair(levels[i], levels[i + 1]) for i in range(len(levels) - 1)]
        self.regrid_every = None      # dict i -> interval in own ticks (None = static)
        self.tol = None
        self.halo = 1
        self.nrebuild = 0
        self.max_cells = 0

    def rebuild_pairs(self, i):
        for p in (i - 1, i):
            if 0 <= p < len(self.pairs):
                self.pairs[p] = Pair(self.levels[p], self.levels[p + 1])

    def sync_all_up(self):
        """bottom-up restriction then top-down prolongation of bands (used after initialisation)"""
        for p in self.pairs:
            p.restrict_only()
        for p in reversed(self.pairs):
            _gather(p.f.flat(), p.p_tgt, p.c.flat(), p.p_src, p.p_w)

    def advance(self, i):
        L = self.levels[i]
        L.tick()
        if i > 0:
            self.advance(i - 1)
            self.advance(i - 1)
            self.pairs[i - 1].sync()
            if self.regrid_every is not None:
                ev = self.regrid_every[i - 1]
                if ev and self.levels[i - 1].ticks % ev == 0:
                    self.regrid(i - 1)

    def step_top(self, n=1):
        top = len(self.levels) - 1
        for _ in range(n):
            self.advance(top)
            if self.regrid_every is not None:
                ev = self.regrid_every[top]
                if ev and self.levels[top].ticks % ev == 0:
                    self.regrid(top)

    # --- dynamic refinement
    def _dilate(self, coords, r, d):
        if len(coords) == 0:
            return coords
        offs = np.array(list(np.ndindex(*([2 * r + 1] * d)))) - r
        return (coords[:, None, :] + offs[None]).reshape(-1, d)

    def needed_blocks(self, i):
        L = self.levels[i]
        d = L.d
        if i == len(self.levels) - 1:
            mx = np.sqrt(L.energy_density().max(axis=1)) if L.nb else np.zeros(0)
        else:
            mx = self.pairs[i].detail() if L.nb else np.zeros(0)
        flag = L.coords[mx > self.tol[i]]
        need = [self._dilate(flag, self.halo, d)]
        if i > 0:
            ch = self.levels[i - 1].coords
            par = np.floor_divide(ch, 2)
            need.append(self._dilate(np.unique(par, axis=0), 1, d))
        need = np.concatenate(need, axis=0) if need else np.zeros((0, d), np.int64)
        if len(need) == 0:
            return need
        need = np.unique(need, axis=0)
        if i < len(self.levels) - 1:
            # keep only blocks whose parent neighbourhood exists at the coarser level
            P = self.levels[i + 1]
            par = np.floor_divide(need, 2)
            ok = np.ones(len(need), bool)
            for off in np.ndindex(*([3] * d)):
                ok &= P.find(par + np.array(off) - 1) < P.nb
            need = need[ok]
        return need

    def regrid(self, i):
        L = self.levels[i]
        need = self.needed_blocks(i)
        newkeys = np.unique(enc(L._wrap(need))) if len(need) else np.zeros(0, np.int64)
        if len(newkeys) == len(L.keys) and np.array_equal(newkeys, L.keys):
            return False
        oldF, oldkeys = L.F, L.keys
        L.set_blocks(need)
        # copy surviving blocks
        if len(oldkeys) and L.nb:
            pos = np.searchsorted(oldkeys, L.keys)
            pos = np.minimum(pos, len(oldkeys) - 1)
            keep = oldkeys[pos] == L.keys
            L.F[:L.nb][keep] = oldF[pos[keep]]
            fresh = ~keep
        else:
            fresh = np.ones(L.nb, bool)
        self.rebuild_pairs(i)
        if i < len(self.levels) - 1 and fresh.any():
            m = np.zeros((L.nb, L.Bd), bool)
            m[fresh] = True
            self.pairs[i].prolong_all(m.reshape(-1))
        self.nrebuild += 1
        self.max_cells = max(self.max_cells, self.cells())
        return True

    # --- measures
    def cells(self):
        return sum(L.nb * L.Bd for L in self.levels)

    def blocks(self):
        return [L.nb for L in self.levels]

    def energy(self, region=None):
        """composite energy (each point counted once); region(x (N,d) in h0 units) -> bool mask"""
        E = 0.0
        for L in self.levels:
            if L.nb == 0:
                continue
            m = L.count_mask()
            e = L.energy_density() * L.h ** L.d
            if region is not None:
                m = m & region(L.centers().reshape(-1, L.d)).reshape(m.shape)
            E += float(e[m].sum())
        return E

    def moments(self, region=None):
        """composite energy and centroid"""
        E, M = 0.0, 0.0
        for L in self.levels:
            if L.nb == 0:
                continue
            m = L.count_mask()
            X = L.centers()
            if region is not None:
                m = m & region(X.reshape(-1, L.d)).reshape(m.shape)
            e = L.energy_density() * L.h ** L.d
            E += float(e[m].sum())
            M = M + (e[m][:, None] * X[m]).sum(0)
        return E, M / max(E, 1e-300)


# ------------------------------------------------------------------ spectral tools (dense, periodic)
def kgrid(shape, d):
    ks = np.meshgrid(*[2 * np.pi * np.fft.fftfreq(n) for n in shape], indexing="ij")
    KA = np.zeros(shape + (3,))
    for j in range(d):
        KA[..., j] = ks[j]
    return KA


def band_projectors(shape, d):
    """lattice eigen-projectors of the two-tick light walk: P+, P-, P0, and omega+ (all on the k grid)"""
    KA = kgrid(shape, d).reshape(-1, 3)
    W = bloch_batch(KA, 0.0, spin=1)
    lam, V = np.linalg.eig(W)
    om = -np.angle(lam) / 2.0
    o = np.argsort(om, axis=1)
    N = len(KA)
    P = []
    for col in (2, 0, 1):        # max omega = +helicity, min = -helicity, middle = static
        v = V[np.arange(N), :, o[:, col]]
        v = v / np.linalg.norm(v, axis=1, keepdims=True)
        P.append(np.einsum("ni,nj->nij", v, v.conj()))
    omp = om[np.arange(N), o[:, 2]]
    return [p.reshape(shape + (3, 3)) for p in P], omp.reshape(shape)


def pol_plus(khat, d):
    """continuum helicity +1 eigenvector of k.S for unit k (in the x-y plane for d=2)"""
    S = [np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]]), np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]]),
         np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]])]
    M = sum(khat[j] * S[j] for j in range(3))
    w, v = np.linalg.eigh(M)
    return v[:, np.argmax(w)]


def packet(shape, d, lam, theta_deg, center, sigma, amp=1.0, helicity=+1, periodic_axes=(), phi_deg=0.0):
    """dense lattice-eigen packet: Gaussian envelope x carrier, projected on the +(or -) band of the walk.
    lam in cells; direction (cos th cos ph, sin th cos ph, sin ph); sigma (scalar or per axis) in cells."""
    th, ph = np.radians(theta_deg), np.radians(phi_deg)
    khat = np.array([np.cos(th) * np.cos(ph), np.sin(th) * np.cos(ph), np.sin(ph)])
    k0 = 2 * np.pi / lam * khat
    for ax in periodic_axes:      # quantise the carrier along periodic axes
        n = shape[ax]
        k0[ax] = 2 * np.pi * np.round(k0[ax] * n / (2 * np.pi)) / n
    grids = np.meshgrid(*[np.arange(n) + 0.5 for n in shape], indexing="ij")
    sig = np.broadcast_to(np.asarray(sigma, float), (d,))
    env = np.ones(shape)
    phase = np.zeros(shape)
    for j in range(d):
        x = grids[j] - center[j]
        if j in periodic_axes:
            env = env
        else:
            env = env * np.exp(-x ** 2 / (2 * sig[j] ** 2))
        phase = phase + k0[j] * grids[j]
    e = pol_plus(khat / np.linalg.norm(khat), d) if helicity > 0 else pol_plus(-khat, d)
    F = (amp * env * np.exp(1j * phase))[..., None] * e
    Fk = np.fft.fftn(F, axes=tuple(range(d)))
    P, _ = band_projectors(shape, d)
    Pr = P[0] if helicity > 0 else P[1]
    Fk = np.einsum("...ij,...j->...i", Pr, Fk)
    return np.fft.ifftn(Fk, axes=tuple(range(d)))


def helicity_fractions(F, d):
    shape = F.shape[:d]
    Fk = np.fft.fftn(F, axes=tuple(range(d)))
    P, _ = band_projectors(shape, d)
    out = []
    for Pr in P:
        G = np.einsum("...ij,...j->...i", Pr, Fk)
        out.append(float((np.abs(G) ** 2).sum()))
    tot = sum(out)
    return [o / tot for o in out]


def restrict_dense(F, d, periodic_axes=()):
    """cubic midpoint restriction of a dense fine array (even sizes) -> coarse array"""
    out = F
    for ax in range(d):
        n = out.shape[ax]
        if ax in periodic_axes:
            r = lambda s: np.roll(out, -s, axis=ax)
            a = WR[0] * r(-1) + WR[1] * r(0) + WR[2] * r(1) + WR[3] * r(2)
        else:
            pad = [(0, 0)] * out.ndim
            pad[ax] = (1, 2)
            P = np.pad(out, pad)
            sl = lambda s: np.take(P, np.arange(n) + s + 1, axis=ax)
            a = WR[0] * sl(-1) + WR[1] * sl(0) + WR[2] * sl(1) + WR[3] * sl(2)
        out = np.take(a, np.arange(0, n, 2), axis=ax)
    return out


def single_level_world(k, d, B, g0, arr, tol, period=None, halo=1, every=8):
    """one-level world holding a dense array, with dynamic regrid on |F| > tol"""
    shp = np.array(arr.shape[:d])
    lo = np.floor_divide(g0, B)
    hi = np.floor_divide(np.asarray(g0) + shp - 1, B)
    coords = np.stack(np.meshgrid(*[np.arange(a, b + 1) for a, b in zip(lo, hi)], indexing="ij"), -1).reshape(-1, d)
    L = Level(k, d, B, coords, period)
    L.set_dense(g0, arr)
    w = World([L])
    w.tol = [tol]
    w.halo = halo
    w.regrid_every = {0: every}
    w.regrid(0)
    return w
