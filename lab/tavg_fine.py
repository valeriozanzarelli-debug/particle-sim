"""tavg_fine.py - livello fine (zoom radiale, tutti gli elettroni) per gli elementi Z = 5..36 e i loro cationi.

Per ogni elemento: (1) aufbau: si provano le configurazioni candidate (generatore meccanico, uguale per tutti) e si tiene
quella di energia HF piu' bassa (la configurazione "emerge"); (2) energie ROHF del termine di Hund del neutro e del
catione; (3) struttura a gusci: massimi di r^2 rho, cuspide di Kato, rho(0), <r> per guscio; (4) stima relativistica al
primo ordine (massa-velocita' + Darwin del nucleo) sul totale e sul potenziale di ionizzazione (solo una STIMA per
segnalare dove la regola non relativistica non basta).
Esiti in lab/out/tavg_fine.json (checkpoint per elemento). Uso: python tavg_fine.py [Zmin Zmax]  (exit 3 = rilanciami)
"""
import os, sys, json, time, math
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_radial as tr
from tavg_blocks import NOBLE, SYMB, NAMES

OUT = os.path.join(HERE, "out")
FN = os.path.join(OUT, "tavg_fine.json")
C_LIGHT = 137.035999084
EH_EV = 27.211386245988


def core_noble(Z):
    for zc in (36, 18, 10, 2):
        if Z > zc:
            return zc
    return 0


def candidates(Z, charge=0):
    """mechanical candidate generator (the same for every element): noble core + valence electrons spread over the
    outer subshells. Rows 2-3: ns^a np^b. K..Kr: 4s^a 3d^c 4p^b (for Ga..Kr also 3d^10 kept, 4s/4p moved)."""
    nel = Z - charge
    zc = core_noble(Z)
    nv = nel - zc
    core = NOBLE.get(zc, "")
    out = []
    if Z > 30:
        core = NOBLE[18] + " 3d10"; nv = nel - 28
        for a in range(0, 3):
            b = nv - a
            if 0 <= b <= 6:
                out.append(f"{core} 4s{a} 4p{b}")
    elif zc in (2, 10) or (zc == 0):
        n = {0: 1, 2: 2, 10: 3}[zc]
        for a in range(0, 3):
            b = nv - a
            if 0 <= b <= 6 and (zc > 0 or b == 0):
                out.append(f"{core} {n}s{a} {n}p{b}")
    elif zc == 18:
        for a in range(0, 3):
            for b in range(0, 2):
                c = nv - a - b
                if 0 <= c <= 10:
                    out.append(f"{core} 3d{c} 4s{a} 4p{b}")
    cl = []
    for s in out:
        toks = [t for t in s.split() if not t.endswith(("s0", "p0", "d0"))]
        cl.append(" ".join(toks))
    return sorted(set(c for c in cl if c))


def nshells_maxima(a):
    r = a.g.r
    D = 4 * math.pi * r * r * a.rho()
    sel = (r > 1e-3 / a.Z) & (r < 20)
    d = D[sel]
    mx = (d[1:-1] > d[:-2]) & (d[1:-1] > d[2:]) & (d[1:-1] > 1e-3 * d.max())
    return int(mx.sum())


def rel_first_order(a):
    """first-order Breit-Pauli one-electron terms on the HF orbitals: mass-velocity -(1/8c^2)<p^4> and nuclear Darwin
    (pi Z / 2c^2) rho(0). Returns (E_mv, E_darwin) in Eh."""
    g = a.g; r = g.r
    Emv = 0.0
    for i, (n, l, q) in enumerate(a.config):
        u = a.P[i] / np.sqrt(r)
        Ppp = (g.D2 @ u - 0.25 * u) / r ** 1.5
        f = Ppp - l * (l + 1) * a.P[i] / r ** 2
        Emv += q * g.integ(f * f)
    Emv *= -1.0 / (8 * C_LIGHT ** 2)
    Ed = math.pi * a.Z / (2 * C_LIGHT ** 2) * a.rho0()
    return Emv, Ed


def run_atom(Z, cfg, grid=None):
    a = tr.Atom(Z, tr.config_from_string(cfg), grid=grid)
    a.scf()
    return a


def describe(a):
    shells = {}
    for i, (n, l, q) in enumerate(a.config):
        shells[f"{n}{tr.L_LETTER[l]}"] = dict(occ=q, eps=a.eps.get(i), r_mean=a.mean_r(i, 1), r2_mean=a.mean_r(i, 2))
    Emv, Ed = rel_first_order(a)
    return dict(E=a.E, config=" ".join(f"{n}{tr.L_LETTER[l]}{q}" for n, l, q in a.config), shells=shells,
                rho0=a.rho0(), kato=a.kato(), maxima=nshells_maxima(a), E_mv=Emv, E_darwin=Ed, it=len(a.hist))


def main():
    Zs = list(range(int(sys.argv[1]), int(sys.argv[2]) + 1)) if len(sys.argv) > 2 else list(range(5, 37))
    res = json.load(open(FN)) if os.path.exists(FN) else {}
    t0 = time.time()
    for Z in Zs:
        key = str(Z)
        rec = res.get(key, {})
        if rec.get("done"):
            continue
        for charge in (0, 1):
            ck = f"q{charge}"
            if ck in rec and rec[ck].get("done"):
                continue
            cands = candidates(Z, charge)
            cr = rec.get(ck, {}).get("candidati", {})
            grid = tr.Grid(Z)
            for c in cands:
                if c in cr:
                    continue
                if time.time() - t0 > 520:
                    rec.setdefault(ck, {})["candidati"] = cr; res[key] = rec
                    json.dump(res, open(FN, "w"), indent=1); print("checkpoint: rilanciami", flush=True); sys.exit(3)
                ts = time.time()
                try:
                    a = run_atom(Z, c, grid)
                    cr[c] = describe(a); cr[c]["secondi"] = time.time() - ts
                    print(f"Z={Z} q={charge} {c}: E={a.E:.8f} ({time.time() - ts:.0f} s)", flush=True)
                except Exception as e:                     # declared: a candidate that does not converge is recorded
                    cr[c] = dict(E=None, errore=repr(e))
                    print(f"Z={Z} q={charge} {c}: ERRORE {e!r}", flush=True)
            best = min((c for c in cr if cr[c].get("E") is not None), key=lambda c: cr[c]["E"])
            rec[ck] = dict(candidati=cr, migliore=best, E=cr[best]["E"], done=True)
            res[key] = rec
            json.dump(res, open(FN, "w"), indent=1)
        rec["done"] = True
        rec["IP_HF_eV"] = (rec["q1"]["E"] - rec["q0"]["E"]) * EH_EV
        r0 = rec["q0"]["candidati"][rec["q0"]["migliore"]]; r1 = rec["q1"]["candidati"][rec["q1"]["migliore"]]
        rec["dIP_rel_eV"] = ((r1["E_mv"] + r1["E_darwin"]) - (r0["E_mv"] + r0["E_darwin"])) * EH_EV
        rec["simbolo"] = SYMB[Z]; rec["nome"] = NAMES[Z]
        res[key] = rec
        json.dump(res, open(FN, "w"), indent=1)
        print(f"== Z={Z} {SYMB[Z]}: {rec['q0']['migliore']} E={rec['q0']['E']:.8f}  IP_HF={rec['IP_HF_eV']:.4f} eV  "
              f"dIP_rel={rec['dIP_rel_eV']:.4f} eV", flush=True)
    print("FINE DONE", flush=True)


if __name__ == "__main__":
    main()
