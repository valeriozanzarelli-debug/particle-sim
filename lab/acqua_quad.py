"""acqua_quad.py - controllo a forza bruta degli integrali di acqua_int (quadratura numerica, niente Boys, niente
ricorrenze di Hermite). 1/r = (2/sqrt(pi)) int_0^inf exp(-t^2 r^2) dt: per ogni t l'integrale si fattorizza in x, y, z;
in ogni direzione griglie uniformi (trapezi, esatti a livello spettrale per gaussiane), per l'ERI coordinate
(x+, x-) con passo in x- che si stringe come 1/t. Quadratura di Gauss-Legendre in s, t = s/(1-s).
Uso: python acqua_quad.py  -> stampa l'errore massimo per tipo; esce con codice 1 se > 1e-8."""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import acqua_int as I


def f1(x, A, i, a):
    y = x - A
    return y ** i * np.exp(-a * y * y)


def tgrid(n):
    s, w = np.polynomial.legendre.leggauss(n); s = 0.5 * (s + 1); w = 0.5 * w
    return s / (1 - s), w / (1 - s) ** 2


def cart_to_sph(blk, Ls, Ns):
    out = blk
    for k, L in enumerate(Ls):
        C = I.SPH[L, :2 * L + 1, :(L + 1) * (L + 2) // 2]
        out = np.moveaxis(np.tensordot(C, out, axes=([1], [k])), 0, k)
    return out * np.prod(Ns)


def one_e_ref(shA, shB, CA, CB, Zs, Rs):
    (la, a), (lb, b) = shA, shB
    x = np.linspace(-14, 14, 5601)
    S1 = np.zeros((3, la + 1, lb + 1)); T1 = np.zeros((3, la + 1, lb + 1)); X1 = np.zeros((3, la + 1, lb + 1))
    for d in range(3):
        for i in range(la + 1):
            for j in range(lb + 1):
                fa = f1(x, CA[d], i, a); fb = f1(x, CB[d], j, b); y = x - CB[d]
                d2 = (j * (j - 1) * y ** (j - 2 if j >= 2 else 0) * (j >= 2) - 2 * b * (2 * j + 1) * y ** j
                      + 4 * b * b * y ** (j + 2)) * np.exp(-b * y * y)
                h = x[1] - x[0]
                S1[d, i, j] = np.sum(fa * fb) * h; T1[d, i, j] = -0.5 * np.sum(fa * d2) * h
                X1[d, i, j] = np.sum(fa * fb * x) * h
    ca, cb = I.cart_list(la), I.cart_list(lb)
    S = np.zeros((len(ca), len(cb))); T = S.copy(); D = np.zeros((3,) + S.shape); V = S.copy()
    for p, A_ in enumerate(ca):
        for q, B_ in enumerate(cb):
            s = [S1[d, A_[d], B_[d]] for d in range(3)]
            S[p, q] = np.prod(s)
            T[p, q] = sum(T1[d, A_[d], B_[d]] * np.prod([s[e] for e in range(3) if e != d]) for d in range(3))
            for d in range(3):
                D[d, p, q] = X1[d, A_[d], B_[d]] * np.prod([s[e] for e in range(3) if e != d])
    # nuclear attraction
    ts, ws = tgrid(300)
    for Z, C in zip(Zs, Rs):
        acc = np.zeros(S.shape)
        for t, w in zip(ts, ws):
            span = min(14.0, 12.0 / max(t, 1e-9)); hh = min(0.005, 0.1 / max(t, 1e-9))
            K = np.zeros((3, la + 1, lb + 1))
            for d in range(3):
                xx = C[d] + np.linspace(-span, span, 2 * int(span / hh) + 1)
                h = xx[1] - xx[0]; ker = np.exp(-t * t * (xx - C[d]) ** 2)
                Fa = np.stack([f1(xx, CA[d], i, a) for i in range(la + 1)]) * ker
                Fb = np.stack([f1(xx, CB[d], j, b) for j in range(lb + 1)])
                K[d] = Fa @ Fb.T * h
            for p, A_ in enumerate(ca):
                for q, B_ in enumerate(cb):
                    acc[p, q] += w * K[0, A_[0], B_[0]] * K[1, A_[1], B_[1]] * K[2, A_[2], B_[2]]
        V += -Z * 2 / math.sqrt(math.pi) * acc
    N = [I.norm_prim(la, a), I.norm_prim(lb, b)]
    return [cart_to_sph(X, (la, lb), N) for X in (S, T, V)] + [np.stack([cart_to_sph(D[d], (la, lb), N) for d in range(3)])]


def eri_ref(shs, Cs, nt=240):
    (la, a), (lb, b), (lc, c), (ld, dd) = shs
    ts, ws = tgrid(nt)
    lst = [I.cart_list(L) for L in (la, lb, lc, ld)]
    acc = np.zeros(tuple(len(l) for l in lst))
    xp = np.linspace(-9, 9, 721); hp = xp[1] - xp[0]
    for t, w in zip(ts, ws):
        ext = min(14.0, 9.0 / t); hm = min(0.025, 0.12 / t)
        xm = np.linspace(-ext, ext, 2 * int(ext / hm) + 1); h2 = (xm[1] - xm[0]) * hp
        X1 = xp[:, None] + 0.5 * xm[None, :]; X2 = xp[:, None] - 0.5 * xm[None, :]
        ker = np.exp(-t * t * xm * xm)[None, :]
        K = []
        for d in range(3):
            Fab = np.stack([f1(X1, Cs[0][d], i, a) * f1(X1, Cs[1][d], j, b) for i in range(la + 1) for j in range(lb + 1)])
            Fcd = np.stack([f1(X2, Cs[2][d], k, c) * f1(X2, Cs[3][d], l, dd) for k in range(lc + 1) for l in range(ld + 1)])
            K.append(np.einsum("pxy,qxy->pq", Fab * ker, Fcd, optimize=True).reshape(la + 1, lb + 1, lc + 1, ld + 1) * h2)
        for idx in itertools.product(*[range(len(l)) for l in lst]):
            A_, B_, C_, D_ = [lst[k][idx[k]] for k in range(4)]
            acc[idx] += w * np.prod([K[d][A_[d], B_[d], C_[d], D_[d]] for d in range(3)])
    acc *= 2 / math.sqrt(math.pi)
    return cart_to_sph(acc, (la, lb, lc, ld), [I.norm_prim(L, e) for L, e in shs])


if __name__ == "__main__":
    cen = np.array([[0.0, 0.0, 0.0], [0.4, -0.3, 1.2], [-1.1, 0.5, -0.4]])
    shells = [(0, 0, 1.7), (0, 1, 0.9), (1, 2, 0.75), (2, 3, 0.6), (0, 4, 1.1), (1, 0, 0.35), (2, 1, 2.5)]
    B = I.Basis(cen, [8, 1, 1], shells)
    o = I.integrals(B)
    worst = {}
    only = [int(x) for x in sys.argv[1:]]          # optional: only these ERI quartets (one-electron part skipped)
    for A, Bq in ([] if only else [(0, 0), (1, 2), (2, 3), (3, 4), (4, 4), (5, 6), (6, 1), (3, 3), (4, 0), (2, 6)]):
        if True:
            ref = one_e_ref(shells[A][1:], shells[Bq][1:], cen[shells[A][0]], cen[shells[Bq][0]], B.Z, cen)
            sa = slice(B.off[A], B.off[A + 1]); sb = slice(B.off[Bq], B.off[Bq + 1])
            for nm, X, Rf in (("S", o["S"], ref[0]), ("T", o["T"], ref[1]), ("V", o["V"], ref[2])):
                worst[nm] = max(worst.get(nm, 0), np.abs(X[sa, sb] - Rf).max())
            worst["D"] = max(worst.get("D", 0), np.abs(o["D"][:, sa, sb] - ref[3]).max())
            print("pair", A, Bq, {k: float(v) for k, v in worst.items()}, flush=True)
    quartets = [(0, 0, 0, 0), (1, 2, 5, 3), (3, 2, 1, 0), (4, 5, 2, 1), (6, 3, 6, 2), (3, 3, 3, 3)]
    for q in ([quartets[i] for i in only] if only else quartets):
        ref = eri_ref([shells[k][1:] for k in q], [cen[shells[k][0]] for k in q])
        sl = tuple(slice(B.off[k], B.off[k + 1]) for k in q)
        err = np.abs(o["G"][sl] - ref).max()
        worst["ERI"] = max(worst.get("ERI", 0), err); print("quartet", q, "err", err, "max", np.abs(ref).max(), flush=True)
    print({k: float(v) for k, v in worst.items()})
    sys.exit(0 if max(worst.values()) < 1e-8 else 1)
