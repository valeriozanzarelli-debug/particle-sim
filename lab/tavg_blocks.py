"""tavg_blocks.py - impara i blocchi-nucleo (tavg_core.Block) degli elementi e li scrive in lab/out/tavg_blocchi/Z??.json
(file piccoli, JSON: si possono committare). Uso:  python tavg_blocks.py 5 6 7 8 9 10   (numeri atomici)

Regola automatica per r_c (dell'agglomeratore): per ogni canale l, r_c e' il PIU' PICCOLO raggio (blocco piu' piccolo,
piu' trasferibile) tale che il livello grosso riesce a portare l'orbitale liscio: energia cinetica oltre
KCUT = 0.8 * pi / h (h = cella della griglia) <= TAIL_BUDGET. Manopole: TAIL_BUDGET e il fattore 0.8 (dichiarate).
"""
import os, sys, json, time, math
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tavg_radial as tr
import tavg_core as tc

OUT = os.path.join(HERE, "out")
BLOCKDIR = os.path.join(OUT, "tavg_blocchi")

NOBLE = {2: "1s2", 10: "1s2 2s2 2p6", 18: "1s2 2s2 2p6 3s2 3p6", 36: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6"}
NAMES = {1: "idrogeno", 2: "elio", 3: "litio", 4: "berillio", 5: "boro", 6: "carbonio", 7: "azoto", 8: "ossigeno",
         9: "fluoro", 10: "neon", 11: "sodio", 12: "magnesio", 13: "alluminio", 14: "silicio", 15: "fosforo",
         16: "zolfo", 17: "cloro", 18: "argon", 19: "potassio", 20: "calcio", 21: "scandio", 22: "titanio",
         23: "vanadio", 24: "cromo", 25: "manganese", 26: "ferro", 27: "cobalto", 28: "nichel", 29: "rame",
         30: "zinco", 31: "gallio", 32: "germanio", 33: "arsenico", 34: "selenio", 35: "bromo", 36: "kripton"}
SYMB = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg",
        13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V",
        24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se",
        35: "Br", 36: "Kr"}


def core_of(Z):
    """noble-gas core string and number of core shells; Ga..Kr also put the full 3d in the block"""
    if Z <= 2:
        return "", 0
    if Z <= 10:
        return NOBLE[2], 1
    if Z <= 18:
        return NOBLE[10], 3
    if Z <= 30:
        return NOBLE[18], 5
    return NOBLE[18] + " 3d10", 6


def valence_shells(Z):
    if Z <= 2:
        return [(1, 0)], None
    if Z <= 10:
        return [(2, 0), (2, 1)], None
    if Z <= 18:
        return [(3, 0), (3, 1)], None
    if Z <= 30:
        return [(4, 0), (3, 2)], (4, 1)
    return [(4, 0), (4, 1)], None


def ground_config(Z, charge=0):
    """HF candidate list is in tavg_aufbau.py; here: the reference (experimental-ground-like) configuration used as the
    block's reference. Only used as a START; the aufbau study decides which is lowest."""
    nel = Z - charge
    core, nc = core_of(Z)
    ncore = sum(int(t[2:]) for t in core.split()) if core else 0
    nv = nel - ncore
    if Z <= 2:
        return f"1s{nv}"
    if Z <= 18 or Z >= 31:
        n = valence_shells(Z)[0][0][0]
        s = min(2, nv); p = nv - s
        return (core + f" {n}s{s}" + (f" {n}p{p}" if p else "")).strip() if s else core
    # K..Zn
    special = {24: (5, 1), 29: (10, 1)}
    if charge == 0:
        if Z in special:
            d, s = special[Z]
        else:
            s = min(2, nv); d = nv - s
    else:
        d, s = TM_CATION.get((Z, charge), (max(nv - 1, 0), min(1, nv)))
    parts = [core]
    if d:
        parts.append(f"3d{d}")
    if s:
        parts.append(f"4s{s}")
    return " ".join(parts)


# cations of K..Zn: experimental ground configurations (used as the START of the aufbau search only)
TM_CATION = {(19, 1): (0, 0), (20, 1): (0, 1), (21, 1): (1, 1), (22, 1): (2, 1), (23, 1): (4, 0), (24, 1): (5, 0),
             (25, 1): (5, 1), (26, 1): (6, 1), (27, 1): (8, 0), (28, 1): (9, 0), (29, 1): (10, 0), (30, 1): (10, 1)}


CACHE = os.path.join(OUT, "tavg_cache")


def get_ae(Z, cfg, mkin=1.0, grid=None):
    """all-electron atom on the fine level, cached on disk (lab/out/tavg_cache/*.npz: never committed)"""
    os.makedirs(CACHE, exist_ok=True)
    key = f"ae_Z{Z:02d}_{cfg.replace(' ', '_')}_m{mkin:.9f}.npz"
    fn = os.path.join(CACHE, key)
    a = tr.Atom(Z, tr.config_from_string(cfg), mkin=mkin, grid=grid)
    if os.path.exists(fn):
        z = np.load(fn)
        if len(z["r"]) == a.g.M and np.allclose(z["r"], a.g.r):
            a.P = {i: z[f"P{i}"] for i in range(len(a.config))}
            a.E = float(z["E"]); a.eps = {i: float(e) for i, e in enumerate(z["eps"])}
            a.hist = list(z["hist"])
            return a
    a.scf()
    np.savez(fn, r=a.g.r, E=a.E, eps=np.array([a.eps.get(i, np.nan) for i in range(len(a.config))]),
             hist=np.array(a.hist), **{f"P{i}": a.P[i] for i in range(len(a.config))})
    return a


def outer_rmax(P, r):
    """radius of the outermost maximum of |P|"""
    a = np.abs(P)
    i = int(np.argmax(a * (r < 30)))
    # outermost local maximum: scan from large r inwards
    for j in range(len(r) - 2, 1, -1):
        if r[j] < 30 and a[j] >= a[j - 1] and a[j] >= a[j + 1] and a[j] > 0.05 * a[i]:
            return float(r[j])
    return float(r[i])


def block_to_json(b, kappa, extra):
    r = b.g.r
    keep = r < 40
    d = dict(Z=b.Z, Zv=float(b.Zv), simbolo=SYMB.get(b.Z), nome=NAMES.get(b.Z), regola_rc=kappa,
             rc={str(l): v for l, v in b.rc.items()}, lloc=b.lloc, r_tail=b.r_tail,
             E_core=float(b.E_core), E_AE_ref=float(b.ae.E), E_ps_ref=float(b.ps_ref["E"]),
             eps={str(l): v for l, v in b.eps.items()},
             config_ref=" ".join(f"{n}{tr.L_LETTER[l]}{q}" for n, l, q in b.ae.config),
             r=r[keep].tolist(), Vloc=b.Vloc[keep].tolist(),
             V={str(l): v[keep].tolist() for l, v in b.V.items()},
             phi={str(l): v[keep].tolist() for l, v in b.phi.items()},
             kb=[dict(l=l, chi=chi[keep].tolist(), E=float(E)) for l, (chi, E) in b.kb.items()],
             mkin=b.mkin, extra=extra)
    return d


TAIL_BUDGET = 2e-4          # Eh per canale
KFRAC = 0.8


def learn(Z, h_grid=0.25, mkin=1.0, rc_override=None, write=True, ae=None, log=print, budget=TAIL_BUDGET):
    kcut = KFRAC * math.pi / h_grid
    kappa = dict(h_grid=h_grid, kcut=kcut, budget=budget)
    t0 = time.time()
    cfg = ground_config(Z)
    core, nc = core_of(Z)
    if ae is None:
        ae = get_ae(Z, cfg, mkin)
    r = ae.g.r
    vs, extra_l = valence_shells(Z)
    rc = {}; extra_refs = {}
    for (n, l) in vs:
        idx = [i for i, c in enumerate(ae.config) if c[0] == n and c[1] == l]
        if idx:
            rc[l] = tc.choose_rc(ae.g, ae.P[idx[0]], l, kcut, budget)[0]
    # channels not occupied in the reference: from an excited/ionic configuration of the same element
    missing = [(n, l) for (n, l) in vs if l not in rc]
    if extra_l is not None:
        missing.append(extra_l)
    nel_v = Z - (sum(int(t[2:]) for t in core.split()) if core else 0)
    for (n, l) in missing:
        # put one valence electron in that shell (the rest stays where it is)
        base = [(nn, ll, qq) for (nn, ll, qq) in tr.config_from_string(cfg)]
        # remove one electron from the outermost occupied shell
        last = max(range(len(base)), key=lambda i: (base[i][0], -base[i][1]) if i >= nc else (-1, 0))
        nn, ll, qq = base[last]
        base[last] = (nn, ll, qq - 1)
        base = [x for x in base if x[2] > 0] + [(n, l, 1)]
        s = " ".join(f"{a}{tr.L_LETTER[b_]}{c}" for a, b_, c in base)
        extra_refs[l] = s
        aex = get_ae(Z, s, mkin, grid=ae.g)
        j = [i for i, c in enumerate(aex.config) if c[0] == n and c[1] == l][0]
        rc[l] = tc.choose_rc(ae.g, aex.P[j], l, kcut, budget)[0]
    if rc_override:
        rc.update(rc_override)
    lloc = 1 if 1 in rc else max(rc)
    b = tc.Block(Z, cfg, nc, rc, mkin=mkin, lloc=lloc, extra_refs=extra_refs, ae=ae)
    d = block_to_json(b, kappa, dict(extra_refs=extra_refs, secondi=time.time() - t0))
    if write:
        os.makedirs(BLOCKDIR, exist_ok=True)
        json.dump(d, open(os.path.join(BLOCKDIR, f"Z{Z:02d}.json"), "w"))
    log(f"Z={Z} blocco: Zv={b.Zv} rc={ {l: round(v, 3) for l, v in rc.items()} } E_core={b.E_core:.6f} "
        f"({time.time() - t0:.0f} s)")
    return b, d


if __name__ == "__main__":
    for z in sys.argv[1:]:
        learn(int(z))
