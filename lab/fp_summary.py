"""fp_summary.py - judges the pre-registered exams (out/fp_prereg.json + fp_amend.json) -> out/fp_summary.json"""
import json, os
import numpy as np
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
res = {}
for pol, target in (('circ+', 1.0), ('circ-', -1.0), ('lin_x', 0.0)):
    f = os.path.join(OUT, f'fp_run_{pol}.json')
    if not os.path.exists(f):
        continue
    r = json.load(open(f))
    if 'F1_end' not in r:
        continue
    c0, tv0 = r['F1_t0']; c1, tv1 = r['F1_end']
    F1 = c1 >= 0.95 and tv1 <= tv0 + 0.05
    sb = float(np.mean(r['spin_balls'])); sf = float(np.mean(r['spin_field']))
    F2 = abs(sb - target) <= 0.05 if target else abs(sb) <= 0.05
    zc = np.unwrap(np.array(r['zc_field']) * 2 * np.pi / r['shape'][2]) * r['shape'][2] / (2 * np.pi)
    vfield = float((zc[-1] - zc[0]) / (len(zc) - 1))
    vb = float(np.mean(r['vz_balls']))
    F3 = abs(vb / vfield - 1) <= 0.01
    res[pol] = {
        'F1': {'corr_t0': c0, 'TV_t0': tv0, 'corr_t96': c1, 'TV_t96': tv1, 'pass': bool(F1)},
        'F2': {'spin_palline_medio': sb, 'spin_campo_medio': sf, 'atteso': target * 1.0012,
               'spin_palline_min_max': [float(np.min(r['spin_balls'])), float(np.max(r['spin_balls']))],
               'spin_campo_min_max': [float(np.min(r['spin_field'])), float(np.max(r['spin_field']))], 'pass': bool(F2)},
        'F3': {'vz_palline': vb, 'v_baricentro_campo': vfield, 'scarto': vb / vfield - 1, 'pass': bool(F3)},
    }
    print(pol, 'F1', F1, round(c1, 4), round(tv1, 4), '| F2', F2, round(sb, 4), '(campo', round(sf, 4), ') | F3', F3, round(vb / vfield - 1, 5))
json.dump(res, open(os.path.join(OUT, 'fp_summary.json'), 'w'), indent=1)
