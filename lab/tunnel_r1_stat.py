"""tunnel_r1_stat.py - corsa R1: tempi stazionari ESATTI della regola (matrice di trasferimento sito per sito).
Esami E2, E3, E4, E6 (parte stazionaria), E8, E11a-f della pre-registrazione. Uscita lab/out/tunnel_r1.json"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tunnel_core import times_lattice, times_dirac, flips_direct

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
PRE = json.load(open(os.path.join(OUT, "tunnel_prereg.json")))
A = PRE["passi_reticolari"]
E = 2.0
V0S = [1.2, 1.6, 2.0, 2.4, 2.8]
LS = [1, 2, 4, 8, 12, 16]
t0 = time.time()


def extrap(vals):
    """fit c0 + c1 a + c2 a^2 sui 4 passi; ritorna c0 e il coefficiente lineare (per vedere se c'e' un errore O(a))"""
    a = np.array(A)
    M = np.stack([np.ones_like(a), a, a ** 2], 1)
    c, *_ = np.linalg.lstsq(M, np.array(vals), rcond=None)
    return float(c[0]), float(c[1])


KEYS = ["tphi", "tau_y", "tau_z", "tau_my", "tau_mz", "tau_yR", "dwell", "T", "rot_time"]
res = dict(V0={}, free={})
for V0 in V0S:
    res["V0"][str(V0)] = {}
    for L in LS:
        per_a = {k: [] for k in KEYS}
        for a in A:
            N = int(round(L / a))
            th = a
            r = times_lattice(E, a, th, np.full(N, th), np.full(N, V0))
            per_a["tphi"].append(r["tphi"])
            per_a["tau_y"].append(r["TV"].real)
            per_a["tau_z"].append(-r["TV"].imag)
            per_a["tau_my"].append(r["Tm"].real)
            per_a["tau_mz"].append(-r["Tm"].imag)
            per_a["tau_yR"].append(r["TVr"].real)
            per_a["dwell"].append(r["dwell"])
            per_a["T"].append(r["T"])
            # orologio delle rotazioni letto dal contatore: tau_rot = i f_w / m (f_w contato con il peso z)
            if a == A[1]:
                fd = flips_direct(E, a, th, np.full(N, th), np.full(N, V0))
                per_a.setdefault("fw_identity_err", abs(fd - r["fw"]) / abs(r["fw"]))
            per_a["rot_time"].append(complex(1j * r["fw"]))   # m = 1
        cont = times_dirac(E, 1.0, V0, L)
        entry = dict(lattice={k: (per_a[k] if k != "rot_time" else [[z.real, z.imag] for z in per_a[k]]) for k in KEYS},
                     fw_identity_err=per_a.get("fw_identity_err"))
        ex = {k: extrap(per_a[k]) for k in KEYS if k != "rot_time"}
        rr = np.array(per_a["rot_time"])
        ex["rot_time_re"] = extrap(rr.real)
        ex["rot_time_im"] = extrap(rr.imag)
        entry["extrap"] = ex
        entry["continuum"] = dict(tphi=cont["tphi"], tau_y=cont["TV"].real, tau_z=-cont["TV"].imag,
                                  tau_my=cont["Tm"].real, tau_mz=-cont["Tm"].imag, tau_yR=cont["TVr"].real, T=cont["T"])
        res["V0"][str(V0)][str(L)] = entry
    print(f"V0={V0} fatto, {time.time() - t0:.1f} s", flush=True)

# zona libera: orologio delle rotazioni = tempo proprio
for L in (4, 8):
    vals = []
    for a in A:
        N = int(round(L / a))
        r = times_lattice(E, a, a, np.full(N, a), np.zeros(N))
        vals.append(r["Tm"].real)
    res["free"][str(L)] = dict(lattice=vals, extrap=extrap(vals), proper_time=L / np.sqrt(E ** 2 - 1))

# ---------------- giudizio
ex = lambda V0, L, k: res["V0"][str(V0)][str(L)]["extrap"][k][0]
co = lambda V0, L, k: res["V0"][str(V0)][str(L)]["continuum"][k]
lat = lambda V0, L, k, ia: res["V0"][str(V0)][str(L)]["lattice"][k][ia]
verdict = {}
e2 = max(abs(ex(V, L, "tphi") / co(V, L, "tphi") - 1) for V in V0S for L in LS)
verdict["E2_hartman_continuo"] = dict(max_rel_err=e2, soglia=1e-3, PASS=bool(e2 < 1e-3))
ia = A.index(0.025)
e3 = max(abs(lat(V, 16, "tphi", ia) - lat(V, 8, "tphi", ia)) / lat(V, 16, "tphi", ia) for V in V0S)
verdict["E3_saturazione"] = dict(max_rel_change_8_16=e3, soglia=1e-3, PASS=bool(e3 < 1e-3))
sp = {str(V): 4 / lat(V, 4, "tphi", ia) for V in V0S}
verdict["E4_velocita_apparente"] = dict(L_over_tphi_L4=sp, PASS=bool(min(sp.values()) > 1))
e6 = max(max(abs(ex(V, L, "tau_y") / ex(V, L, "dwell") - 1), abs(ex(V, L, "tau_yR") / ex(V, L, "dwell") - 1)) for V in V0S for L in LS)
e6fin = max(abs(lat(V, L, "tau_y", 3) / lat(V, L, "dwell", 3) - 1) for V in V0S for L in LS)
verdict["E6_larmor_dwell_stazionario"] = dict(max_rel_err_extrap=e6, max_rel_err_a00125=e6fin, soglia=1e-3, PASS=bool(e6 < 1e-3))
e8 = 0.0
for V in V0S:
    for L in (1, 4, 8):
        for k in ("tau_y", "tau_z"):
            e8 = max(e8, abs(ex(V, L, k) - co(V, L, k)) / max(abs(co(V, L, k)), 1e-3 * abs(co(V, L, "tau_y")) + 1e-12))
verdict["E8_larmor_continuo"] = dict(max_rel_err=e8, soglia=1e-3, PASS=bool(e8 < 1e-3))
e11a = max(abs(res["free"][L]["extrap"][0] / res["free"][L]["proper_time"] - 1) for L in ("4", "8"))
idn = max(res["V0"][str(V)][str(L)]["fw_identity_err"] for V in V0S for L in LS)
verdict["E11a_rotazioni_libere"] = dict(rel_err_proper_time=e11a, identity_err=idn, PASS=bool(e11a < 1e-3 and idn < 1e-6))
rat = {str(V): ex(V, 16, "tau_my") / ex(V, 16, "tau_y") for V in V0S}
e11b = max(abs(rat[str(V)] - (E - V)) for V in V0S)
verdict["E11b_rapporto"] = dict(ratio=rat, eps={str(V): E - V for V in V0S}, max_abs_err=e11b, soglia="1% (assoluto su eps/m)", PASS=bool(e11b < 0.01))
sl = {}
ok = True
for V in V0S:
    eps = E - V
    kap = np.sqrt(1 - eps ** 2)
    smz = (ex(V, 16, "tau_mz") - ex(V, 8, "tau_mz")) / 8
    sz = (ex(V, 16, "tau_z") - ex(V, 8, "tau_z")) / 8
    sl[str(V)] = dict(slope_mz=smz, m_over_kappa=1 / kap, slope_z=sz, eps_over_kappa=eps / kap)
    ok &= abs(smz * kap - 1) < 0.01 and abs(sz - eps / kap) < 0.01 * max(1 / kap, 1e-9)
verdict["E11c_pendenze"] = dict(slopes=sl, PASS=bool(ok))
d = {str(L): ex(2.0, L, "tau_mz") / L for L in LS}
my = max(abs(ex(2.0, L, "tau_my")) / ex(2.0, L, "tau_y") for L in (4, 8, 12, 16))
verdict["E11d_centro_gap"] = dict(tau_mz_over_L=d, max_my_over_y_Lge4=my,
                                  PASS=bool(all(abs(d[str(L)] - 1) < 0.01 for L in (4, 8, 12, 16)) and my < 0.01))
sg = all(np.sign(ex(V, L, "tau_my")) == np.sign(E - V) for V in (1.2, 1.6, 2.4, 2.8) for L in (8, 12, 16))
verdict["E11e_segno"] = dict(PASS=bool(sg))
th = {str(L): ex(2.0, L, "tau_mz") / L for L in (1, 2)}
verdict["E11f_non_causale_sottile"] = dict(tau_mz_over_L=th, continuum_L1=co(2.0, 1, "tau_mz"), PASS=bool(th["1"] < 1 and abs(ex(2.0, 1, "tau_mz") - co(2.0, 1, "tau_mz")) < 1e-3))
res["verdict"] = verdict
res["time_s"] = time.time() - t0
json.dump(res, open(os.path.join(OUT, "tunnel_r1.json"), "w"), indent=1, default=float)
for k, v in verdict.items():
    print(k, "PASS" if v["PASS"] else "FAIL", {kk: vv for kk, vv in v.items() if kk != "PASS"})
print("tempo %.1f s" % (time.time() - t0))
