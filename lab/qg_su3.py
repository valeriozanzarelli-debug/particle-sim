"""qg_su3.py - SU(3) matrices on a lattice, numpy, layout (3, 3, *lattice) (colour axes first, lattice axes last).

Shared by both parts of the QUARK_GLUONI study:
  part 1 (imaginary time, 3+1 D Monte Carlo): Cabibbo-Marinari heat bath + over-relaxation, Wilson action;
  part 2 (real time, 3-D classical Yang-Mills): exp of su(3), forces, energy density, Poynting vector.

The only number of the rule (part 1) is beta = 6/g^2 (Wilson action S = beta sum_P (1 - 1/3 Re tr U_P)).
No mass, no colour charge put by hand anywhere.
"""
import numpy as np

EIN = "ij...,jk...->ik..."


def mm(a, b):
    return np.einsum(EIN, a, b, optimize=False)


def mm3(a, b, c):
    return mm(mm(a, b), c)


def dag(a):
    return np.conj(np.swapaxes(a, 0, 1))


def tr(a):
    return a[0, 0] + a[1, 1] + a[2, 2]


def rtr(a):
    return (a[0, 0] + a[1, 1] + a[2, 2]).real


def det3(a):
    return (a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1])
            - a[0, 1] * (a[1, 0] * a[2, 2] - a[1, 2] * a[2, 0])
            + a[0, 2] * (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]))


def ident(shape, dtype=np.complex128):
    u = np.zeros((3, 3) + tuple(shape), dtype)
    for i in range(3):
        u[i, i] = 1.0
    return u


def sh(a, d, n=1):
    """field at x + n e_d (d = lattice axis index, 0-based on the lattice axes)."""
    return np.roll(a, -n, axis=2 + d)


def reunit(u):
    """Gram-Schmidt on the rows, third row = conj(row0 x row1): exact SU(3)."""
    r0 = u[0]
    r0 = r0 / np.sqrt((np.abs(r0) ** 2).sum(0))
    r1 = u[1]
    r1 = r1 - (np.conj(r0) * r1).sum(0) * r0
    r1 = r1 / np.sqrt((np.abs(r1) ** 2).sum(0))
    r2 = np.conj(np.stack([r0[1] * r1[2] - r0[2] * r1[1],
                           r0[2] * r1[0] - r0[0] * r1[2],
                           r0[0] * r1[1] - r0[1] * r1[0]]))
    return np.stack([r0, r1, r2])


def eps_eps(A, B, C):
    """(1/6) eps_abc eps_a'b'c' A_aa' B_bb' C_cc'  (= det A when A = B = C; = 1 for identities)."""
    tA, tB, tC = tr(A), tr(B), tr(C)
    AB, AC, BC = mm(A, B), mm(A, C), mm(B, C)
    return (tA * tB * tC - tA * tr(BC) - tB * tr(AC) - tC * tr(AB) + tr(mm(AB, C)) + tr(mm(AC, B))) / 6.0


# ---------------------------------------------------------------- SU(2) subgroup tools
SUB = ((0, 1), (0, 2), (1, 2))


def sub_quat(W, i, j):
    """real quaternion (r0, r1, r2, r3) of the 2x2 block (i, j) of W: block ~ r0 + i r.sigma (projection)."""
    w00, w01, w10, w11 = W[i, i], W[i, j], W[j, i], W[j, j]
    return np.stack([(w00 + w11).real / 2, (w01 + w10).imag / 2, (w01 - w10).real / 2, (w00 - w11).imag / 2])


def quat_mat(q):
    """2x2 complex entries of r0 + i r.sigma."""
    a00 = q[0] + 1j * q[3]
    a01 = q[2] + 1j * q[1]
    a10 = -q[2] + 1j * q[1]
    a11 = q[0] - 1j * q[3]
    return a00, a01, a10, a11


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                     a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
                     a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
                     a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1)])


def qdag(a):
    return np.stack([a[0], -a[1], -a[2], -a[3]])


def left_sub(M, q, i, j):
    """M <- a M where a = quaternion q embedded in rows/cols (i, j). In place."""
    a00, a01, a10, a11 = quat_mat(q)
    ri, rj = M[i].copy(), M[j].copy()
    M[i] = a00 * ri + a01 * rj
    M[j] = a10 * ri + a11 * rj


def kp(alpha, rng):
    """Kennedy-Pendleton: x0 with density sqrt(1 - x0^2) exp(alpha x0)."""
    n = alpha.size
    x0 = np.empty(n)
    todo = np.arange(n)
    while todo.size:
        aa = alpha[todo]
        r1 = 1.0 - rng.random(todo.size)
        r2 = rng.random(todo.size)
        r3 = 1.0 - rng.random(todo.size)
        r4 = rng.random(todo.size)
        lam2 = -(np.log(r1) + np.cos(2 * np.pi * r2) ** 2 * np.log(r3)) / (2 * aa)
        ok = r4 * r4 <= 1.0 - lam2
        x0[todo[ok]] = 1.0 - 2.0 * lam2[ok]
        todo = todo[~ok]
    return x0


def rand_su2(x0, rng):
    r = np.sqrt(np.clip(1.0 - x0 * x0, 0, None))
    ct = 2 * rng.random(x0.size) - 1
    ph = 2 * np.pi * rng.random(x0.size)
    st = np.sqrt(1 - ct * ct)
    return np.stack([x0, r * st * np.cos(ph), r * st * np.sin(ph), r * ct])


# ---------------------------------------------------------------- su(3) generators (Gell-Mann / 2)
def gellmann():
    l = np.zeros((8, 3, 3), complex)
    l[0][0, 1] = l[0][1, 0] = 1
    l[1][0, 1] = -1j; l[1][1, 0] = 1j
    l[2][0, 0] = 1; l[2][1, 1] = -1
    l[3][0, 2] = l[3][2, 0] = 1
    l[4][0, 2] = -1j; l[4][2, 0] = 1j
    l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2] = -1j; l[6][2, 1] = 1j
    l[7][0, 0] = l[7][1, 1] = 1 / np.sqrt(3); l[7][2, 2] = -2 / np.sqrt(3)
    return l / 2          # T^a, tr(T^a T^b) = delta/2


TA = gellmann()


def from_comp(c):
    """c: (8, *lat) real -> hermitian traceless (3,3,*lat) = c^a T^a."""
    return np.einsum("a...,aij->ij...", c, TA)


def to_comp(H):
    """hermitian traceless -> components c^a = 2 tr(T^a H)."""
    return 2 * np.einsum("aji,ij...->a...", TA, H).real


def herm_tl(X):
    """traceless hermitian part of (X - X^+)/(2i)."""
    Y = (X - dag(X)) / 2j
    t = tr(Y) / 3
    for i in range(3):
        Y[i, i] = Y[i, i] - t
    return Y


def expi(H, order=12):
    """exp(i H) for hermitian H with small norm (Taylor + scaling and squaring, then exact SU(3) re-projection)."""
    nrm = np.sqrt((np.abs(H) ** 2).sum((0, 1))).max()
    s = 0
    while nrm / 2 ** s > 0.25:
        s += 1
    X = 1j * H / 2 ** s
    out = ident(H.shape[2:])
    term = ident(H.shape[2:])
    for k in range(1, order + 1):
        term = mm(term, X) / k
        out = out + term
    for _ in range(s):
        out = mm(out, out)
    return out
