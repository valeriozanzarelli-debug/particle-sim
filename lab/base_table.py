"""Collects lab/out/base_<exp>_<ver>.json into the markdown results table of docs/REGOLA_BASE.md.
  python lab/base_table.py v1
"""
import json
try:
    sys_ = __import__("sys"); sys_.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ["free2d", "free3d", "gouy", "double", "slit3d", "single", "grating", "airy", "arago", "talbot", "edge",
         "which", "malus", "three", "mz_phase_0_90", "mz_phase_180_270", "bomb", "wheeler"]


def fmt_meta(m):
    if isinstance(m, list):
        ph = sum(x.get("photons", 0) for x in m)
        cl = sum(x.get("clicks", x.get("counts", {}).get("click", 0)) for x in m)
        rt = sum(x.get("runtime", 0) for x in m)
        b = m[0].get("balls", {})
    elif isinstance(m, dict) and "runs" in m:
        return fmt_meta(m["runs"])[:3] + (m.get("runtime", 0),)
    elif isinstance(m, dict) and all(isinstance(v, dict) and "photons" in v for v in m.values()):
        return fmt_meta(list(m.values()))
    else:
        ph = m.get("photons", 0)
        cl = m.get("clicks", m.get("counts", {}).get("click", 0) if "counts" in m else 0)
        rt = m.get("runtime", 0)
        b = m.get("balls", {})
    return ph, cl, b, rt


def main(ver):
    print("| # | esperimento / grandezza | misurato | teoria | tolleranza | esito | fotoni / clic | palline-onda per fotone (livello 0 · tutti i livelli) | tempo |")
    print("|---|---|---|---|---|---|---|---|---|")
    k = 0
    for name in ORDER:
        p = os.path.join(HERE, "out", f"base_{name}_{ver}.json")
        if not os.path.exists(p):
            print(f"| | {name} | (not run) | | | | | | |")
            continue
        d = json.load(open(p, encoding="utf-8"))
        ph, cl, b, rt = fmt_meta(d["meta"])
        balls = f"{b.get('balls_level0', 0):,} · {b.get('balls_all_levels', 0):.2g}".replace(",", ".") if b else "-"
        for i, r in enumerate(d["rows"]):
            k += 1
            extra = f" ({r['extra']})" if r.get("extra") else ""
            print(f"| {k} | {r['name']} | {r['measured']} | {r['theory']} | {r['tol']} | {'**PASS**' if r['ok'] else '**FAIL**'} | "
                  f"{ph} / {cl}{'' if i == 0 else ' (idem)'} | {balls if i == 0 else 'idem'} | {rt:.0f} s{extra if False else ''} |")
    print()
    for name in ORDER:
        p = os.path.join(HERE, "out", f"base_{name}_{ver}.json")
        if os.path.exists(p):
            d = json.load(open(p, encoding="utf-8"))
            for r in d["rows"]:
                if r.get("extra"):
                    print(f"* {r['name']}: {r['extra']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "v1")
