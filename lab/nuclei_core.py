"""nuclei_core.py - multi-Skyrmions ('nuclei made of knots') on the lattice of balls of knot_core.py (read-only reuse).

 - rational-map initial data (Houghton-Manton-Sutcliffe 1998): U = cos f + i sin f n_R . tau,
   n_R = (2 Re R, 2 Im R, 1 - |R|^2) / (1 + |R|^2),  z = tan(theta/2) e^{i phi} = (x + i y) / (r + z_coord)
 - angular integrals B = (1/4pi) int J dOmega, I = (1/4pi) int J^2 dOmega, J = [(1+|z|^2)/(1+|R|^2) |R'|]^2
 - rational-map profile f(r): minimise 4 pi int [r^2 f'^2 + 2B (f'^2 + 1) sin^2 f + I sin^4 f / r^2] dr (Skyrme units)
 - product-ansatz clusters of B = 1 hedgehogs (with isorotations) for the 'aggregation' runs
 - shape measures: rotation-invariant multipole spectrum q_l of the baryon density, peaks ('vertices'), inertia
   tensors, symmetry scores asym(g) = ||b - g b|| / ||b|| (cubic interpolation about the centroid).
Units: Skyrme units of knot_core (energy F_pi/(4e), length 2/(e F_pi)); Bogomolny bound 12 pi^2 |B|.
"""
import math, os, sys, numpy as np, torch
from scipy import ndimage
from scipy.special import sph_harm_y
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import knot_core as kc

S3 = math.sqrt(3.0)
MAPS = {  # (numerator coeffs, denominator coeffs), highest power first (numpy.polyval convention)
    1: ([1, 0], [1]),                                             # R = z  (hedgehog)
    2: ([1, 0, 0], [1]),                                          # R = z^2  (torus)
    3: ([S3 * 1j, 0, -1], [1, 0, -S3 * 1j, 0]),                   # R = (sqrt3 i z^2 - 1)/(z^3 - sqrt3 i z)  (tetrahedron)
    4: ([1, 0, 2 * S3 * 1j, 0, 1], [1, 0, -2 * S3 * 1j, 0, 1]),   # R = (z^4 + 2 sqrt3 i z^2 + 1)/(z^4 - 2 sqrt3 i z^2 + 1) (cube)
}
LIT = {  # Battye & Sutcliffe, Rev. Math. Phys. 14 (2002) 29 [hep-th/0103026], table 14: massless, E/(12 pi^2 B), +-0.001.
    "E_over_B": {1: 1.2322, 2: 1.1791, 3: 1.1462, 4: 1.1201},
    # older values (Battye & Sutcliffe, PRL 79 (1997) 363, table 2, 100^3 grid, E_dis/B_dis) - the ones in the brief
    "E_over_B_BS1997": {1: 1.232, 2: 1.171, 3: 1.143, 4: 1.116},
    # rational-map ansatz (HMS 1998): minimal I and E/(12 pi^2 B) of the ansatz (not the true minimum)
    "I_min": {1: 1.0, 2: 5.81, 3: 13.58, 4: 20.65},
    "E_over_B_ansatz": {1: 1.232, 2: 1.208, 3: 1.184, 4: 1.137},
    "symmetry": {1: "O(3)", 2: "D_inf_h (torus)", 3: "T_d (tetrahedron)", 4: "O_h (cube)"},
    # rms baryon radius Delta r in Skyrme units (BS 1997, table 2)
    "Delta_r_BS1997": {1: 1.034, 2: 1.416, 3: 1.636, 4: 1.860},
}


def Rmap(B, z):
    p, q = MAPS[B]
    return np.polyval(p, z) / np.polyval(q, z)


def dRmap(B, z):
    p, q = MAPS[B]
    P, Q = np.polyval(p, z), np.polyval(q, z)
    return (np.polyval(np.polyder(p), z) * Q - P * np.polyval(np.polyder(q), z)) / Q**2


def angular_integrals(B, n=1200):
    """B and I by quadrature on the sphere (Gauss-Legendre in cos theta, uniform in phi)."""
    x, w = np.polynomial.legendre.leggauss(n)
    th = np.arccos(x); ph = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    z = np.tan(TH / 2) * np.exp(1j * PH)
    J = ((1 + abs(z)**2) / (1 + abs(Rmap(B, z))**2) * abs(dRmap(B, z)))**2
    W = w[:, None] * (2 * np.pi / (2 * n))
    return float((J * W).sum() / (4 * np.pi)), float((J**2 * W).sum() / (4 * np.pi))


def rm_profile(B, I, R=12.0, n=6000):
    """Rational-map profile by direct 1D minimisation (L-BFGS), f(0) = pi, f(R) = 0. Returns r, f, E."""
    r = torch.linspace(0, R, n + 1, dtype=torch.float64); h = R / n
    rn = r.numpy()
    f0 = np.pi * np.exp(-(rn / (0.9 * B**0.5))**1.5) * (1 - (rn / R)**4)
    fi = torch.tensor(f0[1:-1], requires_grad=True)
    opt = torch.optim.LBFGS([fi], lr=1, max_iter=5000, tolerance_grad=1e-12, tolerance_change=1e-15,
                            history_size=50, line_search_fn="strong_wolfe")
    rm = 0.5 * (r[1:] + r[:-1])
    pi_ = torch.tensor([np.pi], dtype=torch.float64); z_ = torch.tensor([0.0], dtype=torch.float64)

    def E():
        f = torch.cat([pi_, fi, z_]); fp = (f[1:] - f[:-1]) / h
        s2 = torch.sin(0.5 * (f[1:] + f[:-1]))**2
        return 4 * np.pi * h * torch.sum(rm**2 * fp**2 + 2 * B * (fp**2 + 1) * s2 + I * s2**2 / rm**2)

    def clo():
        opt.zero_grad(); e = E(); e.backward(); return e
    for _ in range(4):
        opt.step(clo)
    with torch.no_grad():
        f = torch.cat([pi_, fi, z_]).numpy()
        return rn, f, float(E())


def rm_field(B, X, Y, Z, r_prof, f_prof):
    x, y, zc = X.numpy(), Y.numpy(), Z.numpy()
    r = np.sqrt(x * x + y * y + zc * zc)
    z = (x + 1j * y) / (r + zc)          # grids have even L: no point on the negative z axis
    Rv = Rmap(B, z)
    a2 = abs(Rv)**2
    n = np.stack([2 * Rv.real, 2 * Rv.imag, 1 - a2]) / (1 + a2)
    f = np.interp(r, r_prof, f_prof, right=0.0)
    phi = np.concatenate([np.cos(f)[None], np.sin(f)[None] * n])
    return torch.tensor(phi, dtype=X.dtype)


def rot_pi(axis):
    """isorotation by pi about a coordinate axis (3x3 on (pi1, pi2, pi3))."""
    d = -np.ones(3); d["xyz".index(axis)] = 1
    return np.diag(d)


def cluster(X, Y, Z, prof, centers, isos):
    """product ansatz U = U_1 U_2 ... of hedgehogs at 'centers' with isorotations 'isos' (3x3 or None)."""
    phi = None
    for c, R in zip(centers, isos):
        u = kc.hedgehog(X, Y, Z, prof, center=c, iso=R)
        phi = u if phi is None else kc.qmul(phi, u)
    return kc.normalize(phi)


# ----------------------------------------------------------------------------------------------------- measures
def axes(shape, h):
    return [(np.arange(n) - (n - 1) / 2) * h for n in shape]


def moments(w, h):
    """total, centroid, second-moment tensor about centroid (normalised), inertia tensor (normalised)."""
    ax = axes(w.shape, h)
    X, Y, Z = np.meshgrid(*ax, indexing="ij")
    tot = w.sum()
    c = np.array([(w * A).sum() / tot for A in (X, Y, Z)])
    D = [X - c[0], Y - c[1], Z - c[2]]
    M = np.array([[(w * D[i] * D[j]).sum() / tot for j in range(3)] for i in range(3)])
    I = np.trace(M) * np.eye(3) - M
    return float(tot * h**3), c, M, I


def multipoles(w, h, lmax=6):
    """rotation-invariant q_l = sqrt(4pi/(2l+1) sum_m |Q_lm|^2) / int w r^l (about the centroid); q_l in [0, 1]."""
    _, c, _, _ = moments(w, h)
    ax = axes(w.shape, h)
    X, Y, Z = np.meshgrid(*ax, indexing="ij")
    x, y, z = X - c[0], Y - c[1], Z - c[2]
    r = np.sqrt(x * x + y * y + z * z) + 1e-12
    th = np.arccos(np.clip(z / r, -1, 1)); ph = np.arctan2(y, x)
    q = {}
    for l in range(1, lmax + 1):
        s = 0.0
        for m in range(-l, l + 1):
            Q = (w * r**l * np.conj(sph_harm_y(l, m, th, ph))).sum()
            s += abs(Q)**2
        q[l] = float(math.sqrt(4 * math.pi / (2 * l + 1) * s) / (w * r**l).sum())
    return q


def peaks(w, h, frac=0.5):
    """local maxima of w above frac * max (3x3x3 neighbourhood), refined by a weighted centroid of the neighbourhood."""
    mx = ndimage.maximum_filter(w, size=3, mode="constant")
    idx = np.argwhere((w == mx) & (w > frac * w.max()))
    ax = axes(w.shape, h)
    out = []
    for i, j, k in idx:
        sl = w[i - 1:i + 2, j - 1:j + 2, k - 1:k + 2]
        g = np.meshgrid(ax[0][i - 1:i + 2], ax[1][j - 1:j + 2], ax[2][k - 1:k + 2], indexing="ij")
        out.append([float((sl * G).sum() / sl.sum()) for G in g] + [float(w[i, j, k])])
    return np.array(out).reshape(-1, 4)


def distance_classes(P):
    """sorted pairwise distances of peak positions (N x 3)."""
    d = [np.linalg.norm(P[i] - P[j]) for i in range(len(P)) for j in range(i + 1, len(P))]
    return np.sort(np.array(d))


def transform(w, h, c, Rm):
    """(g w)(x) = w(c + Rm^T (x - c)), cubic interpolation, zero outside."""
    n = np.array(w.shape)
    ax = axes(w.shape, h)
    X, Y, Z = np.meshgrid(*ax, indexing="ij")
    P = np.stack([X - c[0], Y - c[1], Z - c[2]]).reshape(3, -1)
    Q = Rm.T @ P + c[:, None]
    idx = Q / h + ((n - 1) / 2)[:, None]
    return ndimage.map_coordinates(w, idx, order=3, mode="constant").reshape(w.shape)


def asym(w, h, c, Rm):
    return float(np.linalg.norm(w - transform(w, h, c, Rm)) / np.linalg.norm(w))


def rot(axis, ang):
    a = np.asarray(axis, float); a = a / np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K
