"""box_amend.py - declared amendment to lab/out/box_prereg.json, written after the first engine runs and BEFORE the
graded runs. Writes lab/out/box_amend.json."""
import json, os, subprocess, datetime

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
head = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
A = {
    "meta": {"written": datetime.datetime.now().isoformat(timespec="seconds"), "git_head_before": head},
    "what_happened": "first A1a trial runs with the pre-registered wall (reflection into the opposite y-mover, phase i): "
        "P_lat/P_pred = 0.66-0.89 for d = 16, 64 (and 1.2-6.8 at d = 4). Cause found: that wall carries a CHIRAL EDGE "
        "STATE (for kx > 0 one band crosses the gap, half-gap 0.036 at kx = 1.2, d = 16; same for phases -i and mixed "
        "i/-i). The walk is a Floquet topological phase and this reflecting wall is a domain wall: the 'vacuum' of the "
        "strip is then not gapped and the pair count is contaminated. With the old wall A1a is a FAIL (kept, reported).",
    "amendment": "new wall = 'mass wall': the y direction is a ring of d strip sites (coin th = 0.1) + 12 wall sites "
        "where the coin angle is th_w = 1.5 (a large mass of the SAME sign as th: the continuum 'infinite-mass' "
        "boundary, which has no edge states). Checked: the gap minimum is at kx = 0 and no band crosses the gap "
        "(d = 8, 16; th_w = 0.5, 1.0, 1.5); with th_w < 0 (opposite sign) Jackiw-Rebbi bound states appear (rejected). "
        "Pairs: initial sea = states with w in (-1.4, 0) at kx = -1.0; pairs = probability ending in (0, 1.4) at "
        "kx = +1.0 (the wall bands live at |w| >= 1.5). m_j = positive quasi-energies below 1.4 at kx = 0.",
    "engine_fact_changes": "with the mass wall the lowest transverse momentum is NOT pi/(2d) any more: d = 8 -> "
        "m_eff = 0.242, d = 16 -> 0.164 (pi/(2d) would give 0.220, 0.140; pi/d gives 0.405, 0.220). Added, not graded: "
        "the continuum infinite-mass Dirac strip formula p d = pi - arctan(p/m) (m_eff = sqrt(m^2 + p^2)).",
    "criteria": "unchanged (A1a 10% for every (d, x), A1b vs sqrt(th^2 + (pi/d)^2) 10%, A1c, A2, A3, A4).",
    "trial_numbers_already_seen_with_new_wall": "d = 4, 16, 64, x = 1, 2, 4: ratios 0.84, 0.55, 1.02 / 0.89, 0.95, 0.95 / "
        "0.99, 1.00, 0.99 (so A1a is expected to FAIL at d = 4 x = 2 and d = 16 x = 1: lattice limits when m_eff is "
        "not small, declared now).",
}
json.dump(A, open(os.path.join(OUT, "box_amend.json"), "w"), indent=1)
print("written")
