"""mol_acqua_geom.py - study mol_: geometries of water dimers and trimers built from a RIGID monomer whose shape
(r_OH, angle HOH) comes from the energy minimisation of the monomer (never from a table). Pure geometry, no energies.

Dimer coordinates (6 for rigid monomers): R = O-O distance, and the orientation of each monomer (Euler angles).
  scan_dimer(...) -> list of geometries: a quasi-random (Halton) set of orientations at each R, plus the 'linear
  H-bond' family (donor H on the O-O axis, acceptor plane tilted by theta_a), the minimum's neighbourhood.
Trimer: cyclic (C1-like) arrangement with three equal O-O distances (equilateral triangle) and each monomer donating
one H along a side: geometry for the 3-body term E3 = E(ABC) - sum E(AB) + sum E(A) ..."""
import numpy as np

ANG = 1.8897261246                  # bohr per angstrom (constant)


def monomer(r_oh, ang_deg):
    """O at origin, molecule in the x-z plane, bisector along +z (H's at +z)"""
    a = np.deg2rad(ang_deg) / 2
    return np.array([[0, 0, 0], [r_oh * np.sin(a), 0, r_oh * np.cos(a)], [-r_oh * np.sin(a), 0, r_oh * np.cos(a)]], float)


def rot(axis, ang):
    axis = np.asarray(axis, float); axis /= np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K


def euler(a, b, c):
    return rot((0, 0, 1), a) @ rot((0, 1, 0), b) @ rot((0, 0, 1), c)


def halton(n, base):
    out = np.zeros(n)
    for i in range(n):
        f, r, k = 1.0, 0.0, i + 1
        while k:
            f /= base; r += f * (k % base); k //= base
        out[i] = r
    return out


def random_rotations(n):
    """quasi-uniform rotations (Halton -> unit quaternions, Shoemake)"""
    u1, u2, u3 = halton(n, 2), halton(n, 3), halton(n, 5)
    q = np.stack([np.sqrt(1 - u1) * np.sin(2 * np.pi * u2), np.sqrt(1 - u1) * np.cos(2 * np.pi * u2),
                  np.sqrt(u1) * np.sin(2 * np.pi * u3), np.sqrt(u1) * np.cos(2 * np.pi * u3)], 1)
    Rs = []
    for x, y, z, w in q:
        Rs.append(np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]))
    return Rs


def place(mono, Rm, centre):
    return mono @ Rm.T + np.asarray(centre, float)


def linear_hbond(mono, R_oo, theta_a_deg=57.0, beta_deg=0.0):
    """donor A at origin with one O-H along +z towards B; acceptor B at (0, 0, R_oo) with its bisector tilted by
    theta_a from the -z axis (lone pairs towards the donor H). Returns (A, B) coordinates, nuclei order O, H, H."""
    r_oh = np.linalg.norm(mono[1]); half = np.arccos(mono[1][2] / r_oh)
    # rotate donor so that H1 lies on +z: rotate about y by -half
    A = place(mono, rot((0, 1, 0), -half) @ rot((0, 0, 1), np.deg2rad(beta_deg)), (0, 0, 0))
    # acceptor: bisector (initially +z) -> pointing away from donor (+z), tilted by theta_a about y
    B = place(mono, rot((0, 1, 0), np.deg2rad(theta_a_deg)), (0, 0, R_oo))
    return A, B


def scan_dimer(mono, R_list, n_orient):
    """returns list of dict(kind, R, A, B): linear H-bond family and quasi-random orientations at each R"""
    out = []
    Rs = random_rotations(2 * n_orient)
    for R in R_list:
        for th in (0.0, 30.0, 57.0, 90.0):
            A, B = linear_hbond(mono, R, th)
            out.append(dict(kind=f"hbond_theta{th:.0f}", R=R, A=A, B=B))
        for k in range(n_orient):
            A = place(mono, Rs[2 * k], (0, 0, 0)); B = place(mono, Rs[2 * k + 1], (0, 0, R))
            out.append(dict(kind=f"rand{k}", R=R, A=A, B=B))
    return out


def cyclic_trimer(mono, R_oo):
    """three O on an equilateral triangle (side R_oo) in the x-y plane; each monomer donates one H along the side
    to the next O (sequential donor-acceptor ring); free H's alternate up/down"""
    r_oh = np.linalg.norm(mono[1]); half = np.arccos(mono[1][2] / r_oh)
    Rc = R_oo / np.sqrt(3)
    Os = [Rc * np.array([np.cos(t), np.sin(t), 0]) for t in np.deg2rad([90, 210, 330])]
    mols = []
    for i in range(3):
        O, On = Os[i], Os[(i + 1) % 3]
        u = (On - O) / np.linalg.norm(On - O)
        # local frame: H1 along u; H2 at angle 2*half from u, out of plane alternately
        ang = 2 * half
        w = np.array([0, 0, 1.0]) * (1 if i % 2 == 0 else -1)
        v = np.cross(w, u); v /= np.linalg.norm(v)
        dirH2 = np.cos(ang) * u + np.sin(ang) * (0.5 * v + 0.866 * w)
        dirH2 /= np.linalg.norm(dirH2)
        mols.append(np.array([O, O + r_oh * u, O + r_oh * dirH2]))
    return mols
