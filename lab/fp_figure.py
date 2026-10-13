"""fp_figure.py - figure of FOTONE A PALLINE: the photon seen from behind (circ+, linear, circ-) and the spin read
from the balls vs the field, tick by tick -> out/fp_figure.png"""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
fig = plt.figure(figsize=(16, 7.2), facecolor='#05060a')
names = {'circ+': 'circolare +', 'lin_x': 'lineare', 'circ-': 'circolare −'}
for i, pol in enumerate(('circ+', 'lin_x', 'circ-')):
    ax = fig.add_axes([0.005 + i * 0.333, 0.43, 0.323, 0.53])
    ax.imshow(Image.open(os.path.join(OUT, f'fp_{pol}_dentro_t60.png')))
    ax.set_axis_off()
    ax.set_title(f'{names[pol]} — visto da dietro, tick 60', color='w', fontsize=12)
ax = fig.add_axes([0.07, 0.09, 0.88, 0.30], facecolor='#05060a')
col = {'circ+': '#ffb347', 'lin_x': '#9ad0ff', 'circ-': '#d58cff'}
for pol in ('circ+', 'lin_x', 'circ-'):
    r = json.load(open(os.path.join(OUT, f'fp_run_{pol}.json')))
    ax.plot(r['t'], r['spin_balls'], '-', color=col[pol], lw=2, label=f'{names[pol]}: dalle palline')
    ax.plot(r['t'], r['spin_field'], ':', color=col[pol], lw=1.2, label=f'{names[pol]}: dal campo')
for y in (1.0012, -1.0012):
    ax.axhline(y, color='w', lw=0.6, alpha=0.5)
ax.set_xlabel('tick (24 tick = una lunghezza d\'onda)', color='w')
ax.set_ylabel('spin per fotone (ħ)', color='w')
ax.tick_params(colors='w')
for s in ax.spines.values():
    s.set_color('#666')
ax.legend(ncol=3, fontsize=8, facecolor='#111', labelcolor='w', loc='center right')
fig.savefig(os.path.join(OUT, 'fp_figure.png'), dpi=110, facecolor=fig.get_facecolor())
