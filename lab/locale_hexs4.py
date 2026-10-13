"""locale_hexs4.py - EXPLORATION after the pre-registration (no exam): S4 composition on the hexagonal 2-D lattice (Valerio's request: geometric ideas). Output printed; numbers copied in docs/REGOLA_LOCALE.md."""
import numpy as np, sys, warnings, itertools
warnings.filterwarnings("ignore")
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from locale_core import *
from locale_cands import HEXD
B=LATTICES['hex'][0]
def hexS4(orders):
    ticks=[]
    for o in orders:
        for m in S4_SEQ:
            ticks.append([(HEXD[i],m) for i in o]); ticks.append([(HEXD[i],m) for i in o[::-1]])
    return make_walk('light',ticks,basis=B,dim=2)
D=dirs_for(2,46)
for name,orders in (("hexS4 one order",[(0,1,2)]),("hexS4 3 cyclic",[(0,1,2),(1,2,0),(2,0,1)])):
    w=hexS4(orders)
    r={k:round(n_for_tol(w,D,0.01,k),1) for k in ('phase','group','diffraction')}
    nf=sum(f[0]=='shift' for f in w.pre)
    N=max(r.values()); cost=2/np.sqrt(3)*N**2*nf/w.T
    print(name,'T',w.T,r,'cost',round(cost,1),'(frozen 2D 2640, S4 square 620, hex cyc 1090)',flush=True)
