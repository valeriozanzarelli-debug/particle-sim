"""mg_gpu_g5.py - giudice dell'esame G5 (lab/out/mg_gpu_prereg.json + mg_gpu_amend.json) e figura dello studio
motore_grafico_gpu. Legge lab/out/mg_gpu_g5_{bassa,alta}_{circp,circm,lin}.json (scritti da src/atomo/gpu_esami.rs)
e lab/out/mg_gpu_g123.json; scrive lab/out/mg_gpu_g5.json e docs/fig/motore_gpu.png.
Soglie (scritte prima): (a) cos(P, d) >= 0.99; (b) circolari: segno giusto, |S| entro 5 % di 1 ad alta, 10 % a bassa,
|S_alta - S_bassa| <= 0.10; (c) lineare |S| <= 0.05 alta, <= 0.10 bassa; (d) elicita' media: segno giusto, entro 0.05
di +-1 / 0; (e) apertura theta entro 10 % di lambda/(pi w0) (circolare +). Medie sui tick 24..48."""
import json, os
import numpy as np

QUI = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(QUI, 'out')
FIG = os.path.join(QUI, '..', 'docs', 'fig')
ATTESO = {'circp': 1.0, 'circm': -1.0, 'lin': 0.0}


def carica(preset, pol):
    p = os.path.join(OUT, f'mg_gpu_g5_{preset}_{pol}.json')
    return json.load(open(p)) if os.path.exists(p) else None


def medie(d):
    r = [x['m'] for x in d['righe'] if 24 <= x['tick'] <= 48]
    dirv = np.array(d['dir'])
    cos = [np.dot(x['poynting'], dirv) / np.linalg.norm(x['poynting']) for x in r]
    return {'cos_flusso_min': float(np.min(cos)), 'spin': float(np.mean([x['spin'] for x in r])),
            'elicita': float(np.mean([x['elicita'] for x in r])), 'theta': d['theta'], 'theta_maxwell': d['theta_maxwell'],
            'cammino': d['cammino'], 'celle_per_lambda': d['celle_per_lambda'], 'n': d['n'], 'palline': d['palline'],
            'ms_tick': d['ms_tick']}


res, ok, mancano = {}, True, []
for pol in ('circp', 'circm', 'lin'):
    for preset in ('bassa', 'alta'):
        d = carica(preset, pol)
        if d is None:
            mancano.append(f'{preset}_{pol}')
            continue
        m = medie(d)
        att = ATTESO[pol]
        tol_s = 0.05 if preset == 'alta' else 0.10
        ca = m['cos_flusso_min'] >= 0.99
        if pol == 'lin':
            cb = abs(m['spin']) <= tol_s
        else:
            cb = np.sign(m['spin']) == np.sign(att) and abs(abs(m['spin']) - 1.0) <= tol_s
        cd = abs(m['elicita'] - att) <= 0.05
        ce = None
        if pol == 'circp':
            ce = m['theta'] is not None and abs(m['theta'] / m['theta_maxwell'] - 1.0) <= 0.10
        m.update({'a_flusso': bool(ca), 'b_c_spin': bool(cb), 'd_elicita': bool(cd), 'e_diffrazione': ce})
        ok &= ca and cb and cd and (ce is not False)
        res[f'{preset}_{pol}'] = m
for pol in ('circp', 'circm'):
    a, b = res.get(f'alta_{pol}'), res.get(f'bassa_{pol}')
    if a and b:
        dd = abs(a['spin'] - b['spin'])
        res[f'confronto_{pol}'] = {'spin_alta_meno_bassa': dd, 'passa': bool(dd <= 0.10)}
        ok &= dd <= 0.10
res['mancano'] = mancano
res['esito'] = 'passa' if ok and not mancano else ('incompleto' if mancano else 'fallisce')
json.dump(res, open(os.path.join(OUT, 'mg_gpu_g5.json'), 'w'), indent=1)
for k, v in res.items():
    print(k, v)

# ── figura: spin dalle palline tick per tick, bassa (CPU) e alta (GPU); G3: CPU contro GPU ad alta ──
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
col = {'circp': '#1f77b4', 'circm': '#d62728', 'lin': '#7f7f7f'}
for pol in ('circp', 'circm', 'lin'):
    for preset, ls in (('bassa', '--'), ('alta', '-')):
        d = carica(preset, pol)
        if d is None:
            continue
        t = [x['tick'] for x in d['righe']]
        s = [x['m']['spin'] for x in d['righe']]
        ax[0].plot(t, s, ls, color=col[pol], lw=1.4, label=f'{pol} {preset} ({d["cammino"]}, {d["celle_per_lambda"]:.0f} c/λ)')
ax[0].axvspan(24, 48, color='#eeeeee', zorder=0)
ax[0].axhline(1, color='k', lw=0.6, ls=':'); ax[0].axhline(-1, color='k', lw=0.6, ls=':')
ax[0].set_xlabel('tick'); ax[0].set_ylabel('spin dalle palline (ħ per fotone)')
ax[0].set_title('G5: stessi fenomeni a bassa (CPU) e alta (GPU)')
ax[0].legend(fontsize=7, ncol=2)
p = os.path.join(OUT, 'mg_gpu_g123.json')
if os.path.exists(p):
    g = json.load(open(p))
    t = [r['tick'] for r in g['righe']]
    ax[1].semilogy(t, [abs(r['gpu']['energia'] / r['cpu']['energia'] - 1) + 1e-12 for r in g['righe']], label='energia')
    ax[1].semilogy(t, [np.linalg.norm(np.subtract(r['gpu']['poynting'], r['cpu']['poynting'])) / np.linalg.norm(r['cpu']['poynting']) + 1e-12 for r in g['righe']], label='Poynting')
    ax[1].semilogy(t, [abs(r['gpu']['spin'] / r['cpu']['spin'] - 1) + 1e-12 for r in g['righe']], label='spin dalle palline')
    ax[1].axhline(1e-4, color='k', lw=0.8, ls=':'); ax[1].axhline(5e-3, color='k', lw=0.8, ls='--')
    ax[1].set_xlabel('tick'); ax[1].set_ylabel('|GPU − CPU| relativo')
    ax[1].set_title(f'G3 ad alta (192³, {g["palline"]} palline)')
    ax[1].legend(fontsize=8)
fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
fig.savefig(os.path.join(FIG, 'motore_gpu.png'), dpi=120)
print('figura scritta')
