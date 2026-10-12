"""nuclei_summary.py - applies the pre-registered exams (lab/out/nuclei_prereg.json) to the run outputs nuclei_<job>.json
and fields nuclei_fields_<job>.npz; N3 naive sum, N3c binding in ANW fit units, N4 isosurface Euler characteristics.
Writes lab/out/nuclei_results.json."""
import json, os, math, numpy as np
from scipy import ndimage
import nuclei_core as nc, knot_core as kc
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out")
H = 0.2
LIT = {1: 1.2322, 2: 1.1791, 3: 1.1462, 4: 1.1201}
LIT97 = {1: 1.232, 2: 1.171, 3: 1.143, 4: 1.116}
RMS97 = {1: 1.034, 2: 1.416, 3: 1.636, 4: 1.860}
UNIT = 129.0 / (4 * 5.44)          # MeV per Skyrme energy unit (ANW fit, docs/NODO_3D.md)
U12 = 12 * math.pi**2 * UNIT
NUC = {2: ("deuterone", 2.224), 3: ("3He / 3H", 7.72), 4: ("4He", 28.30)}


def load(job):
    p = os.path.join(OUT, f"nuclei_{job}.json")
    if not os.path.exists(p): return None, None
    d = json.load(open(p)); f = np.load(os.path.join(OUT, f"nuclei_fields_{job}.npz"))
    return d, f


def euler(b, frac=0.25):
    T = kc.isosurface(b, frac * b.max(), H, -(np.array(b.shape) - 1) / 2 * H)
    V = np.round(T.reshape(-1, 3) / H * 1e5).astype(np.int64)
    _, inv = np.unique(V, axis=0, return_inverse=True); inv = inv.reshape(-1, 3)
    good = (inv[:, 0] != inv[:, 1]) & (inv[:, 1] != inv[:, 2]) & (inv[:, 0] != inv[:, 2])
    inv = inv[good]
    E = set()
    for a, c in ((0, 1), (1, 2), (0, 2)):
        for u, v in zip(inv[:, a], inv[:, c]): E.add((min(u, v), max(u, v)))
    nV = len(np.unique(inv)); nE = len(E); nF = len(inv)
    lab, ncomp = ndimage.label(b > frac * b.max())
    return dict(triangles=int(len(T)), V=nV, E=nE, F=nF, chi=nV - nE + nF, components=int(ncomp))


def ok(x): return "PASS" if x else "FAIL"


def main():
    R = {"N1": {}, "N2": {}, "N3": {}, "N4": {}}
    rm = {B: load(f"rm{B}") for B in (1, 2, 3, 4)}
    fin = {B: (rm[B][0]["final"] if rm[B][0] else None) for B in rm}
    # ---------------- N1
    ex = {}
    for B in (1, 2, 3, 4):
        m = fin[B]
        if m is None: continue
        a = m["asym_grid"]; ev = sorted(m["M_eig"]); q = {int(k): v for k, v in m["q"].items()}
        e = dict(E_over=m["E_over"], lit=LIT[B], dev=m["E_over"] / LIT[B] - 1, dev97=m["E_over"] / LIT97[B] - 1,
                 N1b=ok(abs(m["E_over"] / LIT[B] - 1) < 0.02), B_meas=m["B_meas"], N1d=ok(abs(m["B_meas"] / B - 1) < 0.01),
                 rms_B=m["rms_B"], rms_lit97=RMS97[B], converged=rm[B][0]["converged"], its=rm[B][0]["outer_iterations"],
                 grad_rms=rm[B][0]["grad_rms"], E_start_over=rm[B][0]["start"]["E_over"])
        if B == 2:
            c = [a["Cz30"] < .03, a["Cz45"] < .03, a["Cz90"] < .03, a["Pz"] < .03, a["C2x"] < .03, a["Cx90"] > .3,
                 m["center_over_max"] < .3, abs(ev[2] - ev[1]) / ev[2] < .02, ev[0] / ev[2] < .6, q[2] > .1]
        elif B == 3:
            c = [a["C2x"] < .03, a["C2y"] < .03, a["C2z"] < .03, a["C3"] < .03, a["inv"] > .3, a["Cz90"] > .3, q[2] < .01, q[3] > .1]
        elif B == 4:
            c = [a["Cz90"] < .03, a["Cx90"] < .03, a["C3"] < .03, a["inv"] < .03, a["Cz45"] > .3, q[2] < .01, q[3] < .01, q[4] > .08]
        else:
            c = None
        if c is not None: e["N1a"] = ok(all(c)); e["N1a_checks"] = [bool(x) for x in c]
        e["asym"] = a; e["q"] = q; e["M_eig"] = ev
        ex[B] = e
    if 1 in ex:
        for B in ex: ex[B]["ratio_to_B1_lattice"] = ex[B]["E_over"] / ex[1]["E_over"]; ex[B]["ratio_lit"] = LIT[B] / LIT[1]
    R["N1"]["per_B"] = ex
    if len(ex) == 4:
        R["N1"]["N1c"] = ok(all(ex[B]["E_over"] > ex[B + 1]["E_over"] for B in (1, 2, 3)))
    d64, _ = load("rm4_L64")
    if d64:
        m = d64["final"]
        R["N1"]["box_check_rm4_L64"] = dict(E_over=m["E_over"], dev=m["E_over"] / LIT[4] - 1, B_meas=m["B_meas"], rms_B=m["rms_B"],
                                           q=m["q"], converged=d64["converged"], change_vs_L48=(m["E_over"] / fin[4]["E_over"] - 1) if fin[4] else None)
    # ---------------- N2
    for job in ("pair25", "pair30", "tet4", "sq4"):
        d, f = load(job)
        if d is None: continue
        m = d["final"]; a = m["asym_axis"]; ev = sorted(m["M_eig"]); q = {int(k): v for k, v in m["q"].items()}
        r = dict(E_over=m["E_over"], start_E_over=d["start"]["E_over"], B_meas=m["B_meas"], converged=d["converged"],
                 its=d["outer_iterations"], grad_rms=d["grad_rms"], q=q, M_eig=ev, axis=m["axis"], asym_axis=a,
                 asym_grid=m["asym_grid"], n_peaks=m["n_peaks"], ring_center_over_max=m["ring_center_over_max"],
                 ring_var=m["ring_var"], rms_B=m["rms_B"])
        if job.startswith("pair") and fin[2]:
            c = [abs(m["E_over"] / fin[2]["E_over"] - 1) < .005, a["Cax30"] < .05, a["Cax45"] < .05, a["Cax90"] < .05,
                 a["Pax"] < .05, m["ring_center_over_max"] < .3, m["ring_var"] < .1,
                 abs(ev[2] - ev[1]) / ev[2] < .02, ev[0] / ev[2] < .6]
            r["torus_test"] = ok(all(c)); r["checks"] = [bool(x) for x in c]
            r["axis_y_prediction"] = ok(abs(m["axis"][1]) > .95)
            r["dE_vs_rm2"] = m["E_over"] / fin[2]["E_over"] - 1
        if job in ("tet4", "sq4") and fin[4]:
            q4r = {int(k): v for k, v in fin[4]["q"].items()}[4]
            c = [abs(m["E_over"] / fin[4]["E_over"] - 1) < .005, q[2] < .01, q[3] < .01, abs(q[4] / q4r - 1) < .15]
            r["cube_test"] = ok(all(c)); r["checks"] = [bool(x) for x in c]; r["dE_vs_rm4"] = m["E_over"] / fin[4]["E_over"] - 1
        R["N2"][job] = r
    # ---------------- N3
    if fin[1] and fin[2]:
        b1 = rm[1][1]["b"].astype(np.float64); b1 = np.clip(b1, 0, None)
        r1, r2 = fin[1]["rms_B"], fin[2]["rms_B"]
        d = 2 * math.sqrt(max(r2**2 - r1**2, 0))
        s = d / 2 / H
        nv = ndimage.shift(b1, (s, 0, 0), order=3, mode="constant") + ndimage.shift(b1, (-s, 0, 0), order=3, mode="constant")
        nv = np.clip(nv, 0, None)
        tot, c, M, In = nc.moments(nv, H); ev = sorted(np.linalg.eigvalsh(M))
        pk = nc.peaks(nv, H, 0.5)
        e2 = sorted(fin[2]["M_eig"])
        N3 = dict(d_naive=d, naive=dict(rms=float(math.sqrt(sum(ev))), M_eig=ev, n_peaks=int(len(pk)), q=nc.multipoles(nv, H, 4),
                                        center_over_max=float(nv[23:25, 23:25, 23:25].mean() / nv.max())),
                  B2=dict(rms=r2, M_eig=e2, n_peaks=fin[2]["n_peaks"], q=fin[2]["q"], ring_center_over_max=fin[2]["ring_center_over_max"],
                          ring_var=fin[2]["ring_var"]))
        N3["N3a"] = ok(abs(e2[2] - e2[1]) / e2[2] < .03 and e2[0] < .8 * e2[2] and abs(ev[1] - ev[0]) / ev[1] < .03 and ev[2] > 1.2 * ev[1])
        N3["N3b"] = ok(fin[2]["ring_center_over_max"] < .3 and fin[2]["ring_var"] < .1 and len(pk) == 2)
        R["N3"]["naive_vs_B2"] = N3
    per = {}
    for B in (1, 2, 3, 4):
        if fin[B] is None: continue
        per[B] = dict(rms_B=fin[B]["rms_B"], rms_E=fin[B]["rms_E"], M_eig=sorted(fin[B]["M_eig"]), inertia_eig=sorted(fin[B]["inertia_eig"]),
                      n_peaks=fin[B]["n_peaks"], rms_naive_sqrtB=fin[1]["rms_B"] * B ** (1 / 3) if fin[1] else None)
    R["N3"]["per_B"] = per
    if fin[1]:
        E1 = fin[1]["E_over"]
        bind = {}
        for B in (2, 3, 4):
            if fin[B] is None: continue
            eb_lat = B * (E1 - fin[B]["E_over"]) * U12
            eb_lit = B * (LIT[1] - LIT[B]) * U12
            name, exp = NUC[B]
            bind[B] = dict(nucleus=name, E_bind_lattice_MeV=eb_lat, E_bind_lit_MeV=eb_lit, exp_MeV=exp, ratio_lattice=eb_lat / exp,
                           ratio_lit=eb_lit / exp, fraction_of_mass_lattice=eb_lat / (B * E1 * U12), fraction_exp=exp / (B * 938.9))
        R["N3"]["binding"] = bind
        R["N3"]["units"] = dict(MeV_per_unit=UNIT, MeV_per_12pi2=U12, M1_classical_MeV=E1 * U12, note="ANW fit units (2 inputs)")
        if 2 in bind:
            R["N3"]["N3c"] = ok(0.5 < bind[2]["ratio_lattice"] < 2)
    # ---------------- N4
    for job in ("rm1", "rm2", "rm3", "rm4", "pair25", "pair30", "tet4", "sq4"):
        d, f = load(job)
        if d is None: continue
        R["N4"][job] = euler(np.clip(f["b"].astype(np.float64), 0, None))
    if "rm1" in R["N4"] and "rm2" in R["N4"]:
        R["N4"]["exam"] = ok(R["N4"]["rm1"]["chi"] == 2 and R["N4"]["rm2"]["chi"] == 0)
    json.dump(R, open(os.path.join(OUT, "nuclei_results.json"), "w"), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
